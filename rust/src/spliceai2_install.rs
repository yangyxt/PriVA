//! `priva spliceai2-install {check,register}`: the two checks around
//! `priva spliceai2-prepare prepare` that `SpliceAI2_install` ran as inline
//! Python blocks before these commands.
//!
//! - `check --config C --manifest M --output-dir O`, before preparation: the
//!   configuration names hg38/GRCh38; the manifest and the output folder are
//!   absolute paths without symlinks, `..` or a `.partial` component; the
//!   manifest is an owner-confirmed `downloads_complete` GRCh38 v2.0 manifest
//!   at the pinned source revision, naming no `.partial` path; an existing
//!   output folder holds no `.building-*` entry and, unless it is empty, its
//!   `adapter.ready.json`. Prints the SHA-256 of the configuration and of the
//!   manifest, separated by a space.
//! - `register --config C --manifest M --output-dir O --receipt R
//!   --config-sha S --manifest-sha T`, after preparation: under the folder's
//!   `.install.lock` (taken without waiting), the receipt must describe a
//!   complete, not yet qualified preparation of that manifest, and both
//!   outputs (SNV and indel, each a BGZF VCF and its `.tbi`) must be distinct
//!   non-empty regular files under the folder with the recorded sizes and
//!   SHA-256 sums. The configuration must still have the SHA-256 `check`
//!   printed. Then `spliceai_snv_prescore` and `spliceai_indel_prescore` are
//!   written into it: a line `KEY: ...` at the start of a line is replaced
//!   with `KEY: "PATH"` (JSON quoting), a missing key is added at the end, and
//!   the result must read back as the same settings apart from those two keys.
//!   The file is replaced in one rename, keeping its mode.
//!
//! A failure prints the reason and exits with 1, as the Python blocks did.
//! Compared with the Python blocks by `rust/tests/install_checks_compare.py`.
use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use crate::sha256::{sha256_file, sha256_hex};
use crate::yaml_config::{self, ConfigDocument, ConfigValue, Parsed, Scalar};
use std::collections::BTreeSet;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::os::unix::fs::{MetadataExt, OpenOptionsExt, PermissionsExt};
use std::os::unix::io::AsRawFd;
use std::path::{Component, Path, PathBuf};

const USAGE: &str = "usage: priva spliceai2-install check --config CONFIG --manifest MANIFEST --output-dir OUTPUT_DIR\n\
       priva spliceai2-install register --config CONFIG --manifest MANIFEST --output-dir OUTPUT_DIR --receipt RECEIPT --config-sha SHA256 --manifest-sha SHA256\n\
Check a SpliceAI2 GRCh38 preparation before it runs (check) and record its two outputs in the configuration after it ran (register).";
const REVISION: &str = "cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9";
type Checked<T> = std::result::Result<T, String>;

/// A path as `pathlib.PurePosixPath` writes it: no empty or `.` components
/// and no trailing slash. `..` stays.
fn pure(path: &Path) -> PathBuf {
    let mut out = PathBuf::new();
    for component in path.components() {
        match component {
            Component::CurDir => {}
            other => out.push(other.as_os_str()),
        }
    }
    out
}

fn text(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

fn has_partial_part(path: &Path) -> bool {
    path.components().any(|part| {
        part.as_os_str()
            .to_string_lossy()
            .to_lowercase()
            .contains(".partial")
    })
}

fn is_symlink(path: &Path) -> bool {
    fs::symlink_metadata(path).is_ok_and(|m| m.file_type().is_symlink())
}

/// `path.resolve() == path` for an absolute path: no `..` and no symlink
/// among the parts that exist.
fn resolves_to_itself(path: &Path) -> bool {
    if path.components().any(|part| part == Component::ParentDir) {
        return false;
    }
    let mut prefix = PathBuf::new();
    for component in path.components() {
        prefix.push(component.as_os_str());
        match fs::symlink_metadata(&prefix) {
            Ok(metadata) if metadata.file_type().is_symlink() => return false,
            Ok(_) => {}
            Err(_) => return true,
        }
    }
    true
}

/// `(st_dev, st_ino, st_size, st_mtime_ns, st_ctime_ns)` of `os.stat`.
fn identity(path: &Path) -> Checked<(u64, u64, u64, i64, i64)> {
    let m = fs::metadata(path).map_err(|error| format!("{error}: {}", text(path)))?;
    Ok((
        m.dev(),
        m.ino(),
        m.size(),
        m.mtime() * 1_000_000_000 + m.mtime_nsec(),
        m.ctime() * 1_000_000_000 + m.ctime_nsec(),
    ))
}

fn read_config(bytes: &[u8]) -> Checked<ConfigDocument> {
    let source = std::str::from_utf8(bytes)
        .map_err(|error| format!("The PriVA config is not UTF-8 text: {error}"))?;
    match yaml_config::parse_config(source) {
        Ok(Parsed::Mapping(document)) => Ok(document),
        Ok(_) => Err("The PriVA config is not a YAML mapping".to_owned()),
        Err(error) => Err(format!("Cannot read the PriVA config: {error}")),
    }
}

fn read_json(path: &Path) -> Checked<Json> {
    let bytes = fs::read(path).map_err(|error| format!("{error}: {}", text(path)))?;
    let source = std::str::from_utf8(&bytes).map_err(|error| format!("{}: {error}", text(path)))?;
    json_text::parse(source, false).map_err(|error| format!("{}: {error}", text(path)))
}

/// Whether any value of a `path` key, at any depth, is not a string or names
/// a `.partial` component.
fn names_partial(value: &Json) -> bool {
    match value {
        Json::Object(items) => items.iter().any(|(key, item)| {
            (key == "path"
                && item
                    .as_str()
                    .is_none_or(|path| has_partial_part(Path::new(path))))
                || names_partial(item)
        }),
        Json::Array(items) => items.iter().any(names_partial),
        _ => false,
    }
}

pub(crate) fn check(config: &Path, manifest: &Path, output: &Path) -> Checked<(String, String)> {
    let content = fs::read(config).map_err(|error| format!("{error}: {}", text(config)))?;
    // `yaml.safe_load`: anything but a mapping naming hg38 or GRCh38 stops here.
    let source = std::str::from_utf8(&content)
        .map_err(|error| format!("The PriVA config is not UTF-8 text: {error}"))?;
    let hg38 = match yaml_config::parse_config(source) {
        Ok(Parsed::Mapping(settings)) => {
            matches!(settings.scalar("assembly"), Some(Scalar::Str(name)) if name == "hg38" || name == "GRCh38")
        }
        Ok(_) => false,
        Err(error) => return Err(format!("Cannot read the PriVA config: {error}")),
    };
    if !hg38 {
        return Err(
            "SpliceAI2 preparation requires hg38/GRCh38; no fallback is performed.".to_owned(),
        );
    }
    if !manifest.is_absolute() || !output.is_absolute() {
        return Err("SpliceAI2 manifest and output directory must be absolute paths.".to_owned());
    }
    let (manifest, output) = (pure(manifest), pure(output));
    if has_partial_part(&manifest)
        || has_partial_part(&output)
        || !resolves_to_itself(&manifest)
        || !resolves_to_itself(&output)
        || !manifest.is_file()
    {
        return Err(
            "SpliceAI2 requires a complete regular manifest and a non-partial output path."
                .to_owned(),
        );
    }
    let raw = fs::read(&manifest).map_err(|error| format!("{error}: {}", text(&manifest)))?;
    let data = read_json(&manifest)?;
    let field = |key: &str| data.get(key).and_then(Json::as_str);
    let confirmed = field("confirmed_by").is_some_and(|who| !who.trim().is_empty());
    if !matches!(data, Json::Object(_))
        || field("status") != Some("downloads_complete")
        || !confirmed
        || field("assembly") != Some("GRCh38")
        || field("version") != Some("2.0")
        || field("source_revision") != Some(REVISION)
    {
        return Err("SpliceAI2 requires an owner-confirmed complete GRCh38 v2.0 source manifest at the pinned revision.".to_owned());
    }
    if names_partial(&data) {
        return Err("SpliceAI2 does not consume .partial source paths.".to_owned());
    }
    if output.exists() {
        let entries: Vec<String> = match fs::read_dir(&output) {
            Ok(entries) => entries
                .filter_map(|entry| entry.ok())
                .map(|entry| entry.file_name().to_string_lossy().into_owned())
                .collect(),
            Err(_) => {
                return Err(
                    "SpliceAI2 output is incomplete or still being prepared; refusing reuse."
                        .to_owned(),
                );
            }
        };
        if entries.iter().any(|name| name.starts_with(".building-")) {
            return Err(
                "SpliceAI2 output is incomplete or still being prepared; refusing reuse."
                    .to_owned(),
            );
        }
        if !output.join("adapter.ready.json").exists() && !entries.is_empty() {
            return Err(
                "SpliceAI2 output has no completed receipt; refusing duplicate conversion."
                    .to_owned(),
            );
        }
    }
    Ok((sha256_hex(&content), sha256_hex(&raw)))
}

/// `no_links_or_partial`: absolute, already normal, and neither the path nor
/// a folder above it is a symlink or has `.partial` in its name.
fn plain_path(path: &Path) -> Checked<()> {
    let normal = path.is_absolute()
        && !path.components().any(|part| part == Component::ParentDir)
        && pure(path) == path;
    if !normal {
        return Err("SpliceAI2 receipt paths must be absolute and canonical.".to_owned());
    }
    for part in path.ancestors() {
        let name = part
            .file_name()
            .map(|name| name.to_string_lossy().to_lowercase())
            .unwrap_or_default();
        if name.contains(".partial") || is_symlink(part) {
            return Err(format!(
                "SpliceAI2 refuses symlinks or .partial paths: {}",
                text(path)
            ));
        }
    }
    Ok(())
}

fn verified_sha(path: &Path) -> Checked<String> {
    println!("SpliceAI2: verify SHA256 {}", text(path));
    let checksum = sha256_file(path).map_err(|error| format!("{error}: {}", text(path)))?;
    #[cfg(test)]
    tests::AFTER_HASH.with(|hook| {
        if let Some(hook) = hook.borrow_mut().as_mut() {
            hook(path)
        }
    });
    Ok(checksum)
}

fn is_sha(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}

/// The two keys written into the text, as the Python block wrote them.
fn with_keys(original: &str, updates: &[(&str, String)]) -> Checked<String> {
    let mut text = original.to_owned();
    for (key, value) in updates {
        let head = format!("{key}:");
        let starts: Vec<usize> = std::iter::once(0)
            .chain(text.match_indices('\n').map(|(at, _)| at + 1))
            .filter(|&at| text[at..].starts_with(&head))
            .collect();
        if starts.len() > 1 {
            return Err(format!("Duplicate SpliceAI2 config key: {key}"));
        }
        let replacement = format!("{key}: {}", json_text::quote(value, true));
        match starts.first() {
            Some(&start) => {
                let end = text[start..]
                    .find(['\r', '\n'])
                    .map_or(text.len(), |at| start + at);
                text.replace_range(start..end, &replacement);
            }
            None => {
                if !text.is_empty() && !text.ends_with('\n') {
                    text.push('\n');
                }
                text.push_str(&replacement);
                text.push('\n');
            }
        }
    }
    Ok(text)
}

/// Whether `after` holds the settings of `before` with `updates` applied.
fn same_settings(
    before: &ConfigDocument,
    after: &ConfigDocument,
    updates: &[(&str, String)],
) -> bool {
    let keys = |document: &ConfigDocument| -> BTreeSet<String> {
        document
            .entries
            .iter()
            .map(|entry| entry.key.clone())
            .collect()
    };
    let mut expected = keys(before);
    expected.extend(updates.iter().map(|(key, _)| key.to_string()));
    expected == keys(after)
        && expected
            .iter()
            .all(|key| match updates.iter().find(|(name, _)| name == key) {
                Some((_, value)) => {
                    after.get(key) == Some(&ConfigValue::Scalar(Scalar::Str(value.clone())))
                }
                None => after.get(key) == before.get(key),
            })
}

pub(crate) fn register(
    config: &Path,
    manifest: &Path,
    output: &Path,
    receipt: &Path,
    config_sha: &str,
    manifest_sha: &str,
) -> Checked<()> {
    plain_path(output)?;
    plain_path(receipt)?;
    let regular = fs::symlink_metadata(receipt).is_ok_and(|m| m.is_file() && m.len() > 0);
    if receipt.parent() != Some(output) || !regular {
        return Err(
            "SpliceAI2 requires a nonempty regular adapter receipt under OUTPUT_DIR.".to_owned(),
        );
    }
    let lock_path = output.join(".install.lock");
    plain_path(&lock_path)?;
    let lock = OpenOptions::new()
        .append(true)
        .create(true)
        .mode(0o666)
        .open(&lock_path)
        .map_err(|error| format!("{error}: {}", text(&lock_path)))?;
    // SAFETY: a valid descriptor owned by `lock`, held until the process ends.
    if unsafe { libc::flock(lock.as_raw_fd(), libc::LOCK_EX | libc::LOCK_NB) } != 0 {
        return Err("SpliceAI2 preparation is active; refusing receipt reuse.".to_owned());
    }
    let building = fs::read_dir(output)
        .map_err(|error| format!("{error}: {}", text(output)))?
        .filter_map(|entry| entry.ok())
        .any(|entry| {
            entry
                .file_name()
                .to_string_lossy()
                .starts_with(".building-")
        });
    if building {
        return Err(
            "SpliceAI2 output is incomplete or still being prepared; refusing reuse.".to_owned(),
        );
    }
    let mut tracked = vec![
        (receipt.to_path_buf(), identity(receipt)?),
        (manifest.to_path_buf(), identity(manifest)?),
        (config.to_path_buf(), identity(config)?),
    ];
    let data = read_json(receipt)?;
    if !matches!(data, Json::Object(_)) {
        return Err("SpliceAI2 adapter receipt must be an object.".to_owned());
    }
    let text_field = |value: &str| Json::Str(value.to_owned());
    let expected = [
        ("schema_version", Json::Int(1)),
        ("resource", text_field("SpliceAI2")),
        ("version", text_field("2.0")),
        ("assembly", text_field("GRCh38")),
        ("status", text_field("prepared_annotation_only")),
        ("scientific_qualified", Json::Bool(false)),
        ("coverage", text_field("complete_published")),
        ("source_revision", text_field(REVISION)),
        ("source_manifest_sha256", text_field(manifest_sha)),
    ];
    for (key, value) in &expected {
        let same = match (data.get(key), value) {
            (Some(Json::Int(a)), Json::Int(b)) => a == b,
            (Some(Json::Bool(a)), Json::Bool(b)) => a == b,
            (Some(Json::Str(a)), Json::Str(b)) => a == b,
            _ => false,
        };
        if !same {
            return Err(format!(
                "SpliceAI2 receipt does not qualify full cache preparation: {key}"
            ));
        }
    }
    if verified_sha(manifest)? != manifest_sha {
        return Err("SpliceAI2 source manifest changed during preparation.".to_owned());
    }
    let original = fs::read(config).map_err(|error| format!("{error}: {}", text(config)))?;
    if sha256_hex(&original) != config_sha {
        return Err(
            "PriVA config changed during SpliceAI2 preparation; refusing to overwrite it."
                .to_owned(),
        );
    }
    let outputs = match data.get("outputs") {
        Some(Json::Object(items)) => items,
        _ => return Err("SpliceAI2 receipt requires both SNV and indel outputs.".to_owned()),
    };
    let roles: BTreeSet<&str> = outputs.iter().map(|(role, _)| role.as_str()).collect();
    if roles != BTreeSet::from(["snv", "indel"]) {
        return Err("SpliceAI2 receipt requires both SNV and indel outputs.".to_owned());
    }
    let mut paths = Vec::new();
    let mut seen: Vec<PathBuf> = Vec::new();
    for role in ["snv", "indel"] {
        let Some(item @ Json::Object(_)) = data.get("outputs").and_then(|o| o.get(role)) else {
            return Err(format!("SpliceAI2 output entry must be an object: {role}"));
        };
        for (key, checksum_key) in [("path", "sha256"), ("index_path", "index_sha256")] {
            let (Some(value), Some(checksum)) = (
                item.get(key).and_then(Json::as_str),
                item.get(checksum_key).and_then(Json::as_str),
            ) else {
                return Err(format!(
                    "SpliceAI2 output path/checksum is missing or invalid: {role}"
                ));
            };
            if !is_sha(checksum) {
                return Err(format!(
                    "SpliceAI2 output path/checksum is missing or invalid: {role}"
                ));
            }
            let path = pure(Path::new(value));
            plain_path(&path)?;
            let under = path.starts_with(output) && path != output;
            let regular = fs::metadata(&path).is_ok_and(|m| m.is_file() && m.len() > 0);
            if !under || !regular || seen.contains(&path) {
                return Err(format!(
                    "SpliceAI2 output must be a distinct nonempty regular file under OUTPUT_DIR: {}",
                    text(&path)
                ));
            }
            seen.push(path.clone());
            tracked.push((path.clone(), identity(&path)?));
            if key == "path" {
                let size = fs::metadata(&path).map(|m| m.len() as i128).unwrap_or(-1);
                if !matches!(item.get("bytes"), Some(Json::Int(bytes)) if *bytes == size) {
                    return Err(format!("SpliceAI2 output size mismatch: {role}"));
                }
                let name = text(&path);
                if !name.ends_with(".vcf.gz") && !name.ends_with(".vcf.bgz") {
                    return Err(format!(
                        "SpliceAI2 output must have a BGZF VCF filename: {role}"
                    ));
                }
                paths.push(name);
            } else if Some(value.to_owned())
                != item
                    .get("path")
                    .and_then(Json::as_str)
                    .map(|path| format!("{path}.tbi"))
            {
                return Err(format!(
                    "SpliceAI2 tabix index must match its VCF path: {role}"
                ));
            }
            if verified_sha(&path)? != checksum {
                return Err(format!("SpliceAI2 output SHA256 mismatch: {}", text(&path)));
            }
        }
    }
    let updates = [
        ("spliceai_snv_prescore", paths[0].clone()),
        ("spliceai_indel_prescore", paths[1].clone()),
    ];
    let settings = read_config(&original)?;
    let source = String::from_utf8(original)
        .map_err(|error| format!("The PriVA config is not UTF-8 text: {error}"))?;
    let updated = with_keys(&source, &updates)?;
    if !same_settings(&settings, &read_config(updated.as_bytes())?, &updates) {
        return Err(
            "SpliceAI2 config update would alter unrelated settings; refusing it.".to_owned(),
        );
    }
    for (path, before) in &tracked {
        if identity(path)? != *before {
            return Err(format!(
                "SpliceAI2 verification input changed: {}",
                text(path)
            ));
        }
    }
    // Keep every other key and comment; provenance stays in the receipt.
    let mode = fs::metadata(config)
        .map_err(|error| error.to_string())?
        .mode()
        & 0o7777;
    let folder = config.parent().unwrap_or(Path::new("."));
    let name = config
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_default();
    let temporary = folder.join(format!(".{name}.spliceai2-{}", std::process::id()));
    let written = (|| -> std::io::Result<()> {
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .mode(0o600)
            .open(&temporary)?;
        file.write_all(updated.as_bytes())?;
        file.sync_all()?;
        fs::set_permissions(&temporary, fs::Permissions::from_mode(mode))?;
        fs::rename(&temporary, config)
    })();
    if let Err(error) = written {
        let _ = fs::remove_file(&temporary);
        return Err(format!(
            "Cannot write the PriVA config {}: {error}",
            text(config)
        ));
    }
    println!("SpliceAI2 preparation recorded; scientific activation remains pending.");
    drop(lock);
    Ok(())
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args
        .first()
        .is_none_or(|arg| arg == "--help" || arg == "-h")
    {
        println!("{USAGE}");
        return Ok(());
    }
    let mode = args[0].as_str();
    let names: &[&str] = match mode {
        "check" => &["config", "manifest", "output-dir"],
        "register" => &[
            "config",
            "manifest",
            "output-dir",
            "receipt",
            "config-sha",
            "manifest-sha",
        ],
        _ => return Err(invalid(USAGE)),
    };
    let mut values: Vec<Option<String>> = vec![None; names.len()];
    let mut rest = args[1..].iter();
    while let Some(arg) = rest.next() {
        let Some(at) = arg
            .strip_prefix("--")
            .and_then(|name| names.iter().position(|n| *n == name))
        else {
            return Err(invalid(format!("unrecognized argument {arg}\n{USAGE}")));
        };
        values[at] = Some(
            rest.next()
                .ok_or_else(|| invalid(format!("{arg} needs a value")))?
                .clone(),
        );
    }
    let Some(values) = values.into_iter().collect::<Option<Vec<String>>>() else {
        return Err(invalid(USAGE));
    };
    let path = |at: usize| Path::new(&values[at]);
    let outcome = if mode == "check" {
        check(path(0), path(1), path(2)).map(|(config, manifest)| println!("{config} {manifest}"))
    } else {
        register(path(0), path(1), path(2), path(3), &values[4], &values[5])
    };
    if let Err(message) = outcome {
        eprintln!("priva: {message}");
        std::process::exit(1);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;

    type Hook = Box<dyn FnMut(&Path)>;

    thread_local! {
        /// Runs after each checksum, to change a file between its hash and
        /// the final identity check.
        pub(super) static AFTER_HASH: RefCell<Option<Hook>> = const { RefCell::new(None) };
    }

    #[test]
    fn output_changed_after_its_checksum_is_refused_before_the_config_write() {
        let root = std::env::temp_dir().join(format!("priva-spliceai2-{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        let output = root.join("out");
        fs::create_dir_all(&output).unwrap();
        let config = root.join("config.yaml");
        let original = b"assembly: hg38\nspliceai_snv_prescore: old\n".to_vec();
        fs::write(&config, &original).unwrap();
        let manifest = root.join("manifest.json");
        fs::write(&manifest, b"{}\n").unwrap();
        let mut outputs = Vec::new();
        for role in ["snv", "indel"] {
            let vcf = output.join(format!("{role}.vcf.gz"));
            let index = output.join(format!("{role}.vcf.gz.tbi"));
            fs::write(&vcf, format!("{role} records")).unwrap();
            fs::write(&index, format!("{role} index")).unwrap();
            outputs.push(format!(
                r#""{role}": {{"path": "{}", "sha256": "{}", "index_path": "{}", "index_sha256": "{}", "bytes": {}}}"#,
                text(&vcf),
                sha256_file(&vcf).unwrap(),
                text(&index),
                sha256_file(&index).unwrap(),
                fs::metadata(&vcf).unwrap().len()
            ));
        }
        let manifest_sha = sha256_file(&manifest).unwrap();
        let receipt = output.join("adapter.ready.json");
        fs::write(
            &receipt,
            format!(
                r#"{{"schema_version": 1, "resource": "SpliceAI2", "version": "2.0", "assembly": "GRCh38", "status": "prepared_annotation_only", "scientific_qualified": false, "coverage": "complete_published", "source_revision": "{REVISION}", "source_manifest_sha256": "{manifest_sha}", "outputs": {{{}}}}}"#,
                outputs.join(", ")
            ),
        )
        .unwrap();
        let target = output.join("snv.vcf.gz");
        let changed = std::rc::Rc::new(std::cell::Cell::new(false));
        let seen = changed.clone();
        AFTER_HASH.with(|hook| {
            *hook.borrow_mut() = Some(Box::new(move |path: &Path| {
                if path == target && !seen.get() {
                    let mut file = OpenOptions::new().append(true).open(path).unwrap();
                    file.write_all(b"changed after checksum").unwrap();
                    seen.set(true);
                }
            }))
        });
        let outcome = register(
            &config,
            &manifest,
            &output,
            &receipt,
            &sha256_hex(&original),
            &manifest_sha,
        );
        AFTER_HASH.with(|hook| *hook.borrow_mut() = None);
        let config_after = fs::read(&config).unwrap();
        fs::remove_dir_all(&root).unwrap();
        assert!(changed.get());
        assert_eq!(
            outcome,
            Err(format!(
                "SpliceAI2 verification input changed: {}",
                text(&root.join("out/snv.vcf.gz"))
            ))
        );
        assert_eq!(config_after, original);
    }

    #[test]
    fn keys_written_as_the_python_block_wrote_them() {
        let updates = [
            ("spliceai_snv_prescore", "/a/snv.vcf.gz".to_owned()),
            ("b", "x\u{e9}".to_owned()),
        ];
        assert_eq!(
            with_keys("a: 1\nspliceai_snv_prescore: old # note\r\nc: 2", &updates).unwrap(),
            "a: 1\nspliceai_snv_prescore: \"/a/snv.vcf.gz\"\r\nc: 2\nb: \"x\\u00e9\"\n"
        );
        assert!(
            with_keys("b: 1\nb: 2\n", &updates)
                .unwrap_err()
                .contains("Duplicate")
        );
        // A key inside a nested block is not at the start of a line.
        assert_eq!(
            with_keys("", &updates[..1]).unwrap(),
            "spliceai_snv_prescore: \"/a/snv.vcf.gz\"\n"
        );
    }

    #[test]
    fn normal_paths() {
        assert_eq!(pure(Path::new("/a//b/./c/")), PathBuf::from("/a/b/c"));
        assert!(!resolves_to_itself(Path::new("/a/../b")));
        assert!(resolves_to_itself(Path::new("/no/such/folder/here")));
        assert!(has_partial_part(Path::new("/x/Y.PARTIAL/z")));
    }
}
