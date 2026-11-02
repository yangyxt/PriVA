//! Installation boundary between legacy cache builders and native runtime readers.
//!
//! Builders still produce/read their original pickles. `prepare` gives them a
//! private config; `publish` converts completed files and changes the user's
//! existing runtime keys only after every selected cache succeeds. Per-archive
//! receipts retain the exact builder inputs for the next installation. Nothing
//! here changes a scientific producer, its order, or its cache refresh policy.
use crate::cache_convert as caches;
use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use crate::sha256::{sha256_file, sha256_hex};
use crate::yaml_config::{self, ConfigDocument, Parsed, Scalar};
use std::collections::BTreeMap;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

const KEYS: [&str; 10] = [
    "clinvar_aa_stat",
    "clinvar_splice_stat",
    "alphamissense_tranx_domain_map",
    "clinvar_patho_exon_af_stat",
    "clinvar_gene_stat",
    "clinvar_patho_af_stat",
    "pm1_regions_pkl",
    "all_intolerant_domains",
    "interpro_mapping_pickle",
    "clingen_map",
];
const HELP: &str = "usage: priva cache-install {prepare|publish} --config CONFIG --build-config PRIVATE_CONFIG --source-root REPO\n\
prepare creates a private builder config; publish converts/reuses completed caches and atomically updates CONFIG.\n\
Keep the private config and its adjacent receipt until publish succeeds. No annotation is rerun by this command.";
static NEXT: AtomicU64 = AtomicU64::new(0);

fn text(path: &Path) -> Result<&str> {
    path.to_str()
        .ok_or_else(|| invalid("cache-install requires UTF-8 paths"))
}

fn nonempty(path: &Path) -> Result<()> {
    let info = fs::metadata(path)?;
    if !info.is_file() || info.len() == 0 {
        return Err(invalid(format!(
            "{}: required nonempty file is missing",
            path.display()
        )));
    }
    Ok(())
}

fn document(contents: &str) -> Result<ConfigDocument> {
    match yaml_config::parse_config(contents).map_err(|e| invalid(e.to_string()))? {
        Parsed::Mapping(doc) if !doc.has_duplicate_keys() => Ok(doc),
        _ => Err(invalid(
            "cache-install needs one YAML mapping without duplicate keys",
        )),
    }
}

fn setting<'a>(doc: &'a ConfigDocument, key: &str) -> Result<Option<&'a str>> {
    match doc.scalar(key) {
        None if doc.get(key).is_none() => Ok(None),
        Some(Scalar::Null) => Ok(None),
        Some(Scalar::Str(value)) if value.is_empty() => Ok(None),
        Some(Scalar::Str(value)) => Ok(Some(value)),
        _ => Err(invalid(format!("{key} must be a path string or null"))),
    }
}

fn required<'a>(doc: &'a ConfigDocument, key: &str) -> Result<&'a str> {
    setting(doc, key)?.ok_or_else(|| invalid(format!("configure {key} first")))
}

fn profile(doc: &ConfigDocument) -> Result<String> {
    let assembly = match required(doc, "assembly")? {
        "hg19" | "GRCh37" => "GRCh37",
        "hg38" | "GRCh38" => "GRCh38",
        _ => {
            return Err(invalid(
                "unsupported assembly; expected hg19/GRCh37 or hg38/GRCh38",
            ));
        }
    };
    let release = match doc.scalar("vep_installed_version") {
        Some(Scalar::Str(value) | Scalar::Int(value))
            if !value.is_empty() && value.bytes().all(|b| b.is_ascii_digit()) =>
        {
            value
        }
        _ => return Err(invalid("vep_installed_version must be a numeric release")),
    };
    Ok(format!("{assembly}/VEP{release}"))
}

fn resolve(value: &str, repo: &Path) -> PathBuf {
    let path = Path::new(value);
    if path.is_absolute() {
        path.to_owned()
    } else {
        repo.join(path)
    }
}

fn rewrite(contents: &str, updates: &BTreeMap<String, String>) -> Result<String> {
    let doc = document(contents)?;
    let mut spans = Vec::new();
    for (key, value) in updates {
        let entry = doc
            .entry(key)
            .ok_or_else(|| invalid(format!("missing original key {key}")))?;
        spans.push((
            entry.value_start,
            entry.value_end,
            yaml_config::quoted_for_yaml(value),
        ));
    }
    spans.sort_by_key(|x| std::cmp::Reverse(x.0));
    let mut result = contents.to_owned();
    for (start, end, value) in spans {
        result.replace_range(start..end, &value);
    }
    Ok(result)
}

fn sidecar(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(suffix);
    PathBuf::from(name)
}

fn write_new(path: &Path, data: &[u8]) -> Result<()> {
    let mut file = OpenOptions::new().write(true).create_new(true).open(path)?;
    if let Err(error) = file.write_all(data).and_then(|_| file.sync_all()) {
        let _ = fs::remove_file(path);
        return Err(error.into());
    }
    Ok(())
}

/// Replace by rename, never by truncating an archive/config a running reader uses.
fn atomic_write(path: &Path, data: &[u8]) -> Result<()> {
    let temporary = sidecar(
        path,
        &format!(
            ".{}.{}.tmp",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ),
    );
    write_new(&temporary, data)?;
    let result = (|| -> Result<()> {
        if path.exists() {
            fs::set_permissions(&temporary, fs::metadata(path)?.permissions())?;
        }
        fs::rename(&temporary, path)?;
        if let Some(parent) = path.parent() {
            fs::File::open(parent)?.sync_all()?;
        }
        Ok(())
    })();
    if result.is_err() {
        let _ = fs::remove_file(&temporary);
    }
    result
}

fn json(path: &Path) -> Result<Json> {
    json_text::parse(&fs::read_to_string(path)?, true).map_err(invalid)
}

fn field<'a>(value: &'a Json, key: &str) -> Result<&'a str> {
    value
        .get(key)
        .and_then(Json::as_str)
        .ok_or_else(|| invalid(format!("receipt lacks string {key}")))
}

fn receipt(fields: &[(&str, &str)]) -> Result<Vec<u8>> {
    let value = Json::Object(
        fields
            .iter()
            .map(|(k, v)| ((*k).to_owned(), Json::Str((*v).to_owned())))
            .collect(),
    );
    Ok((json_text::dumps_indented(&value).map_err(invalid)? + "\n").into_bytes())
}

fn bundle(key: &str, repo: &Path, selected_profile: &str) -> Option<(PathBuf, PathBuf)> {
    let name = match key {
        "interpro_mapping_pickle" => "data/InterPro/Interpro_entry_mapping".to_owned(),
        "clingen_map" => format!(
            "data/clingen/clingen_map.{}",
            if selected_profile.starts_with("GRCh37/") {
                "hg19"
            } else {
                "hg38"
            }
        ),
        // The full-value ClinVar gene archive exceeds the distribution limit.
        // Keep its small pickle in Git and convert under cache_root at install.
        _ => return None,
    };
    Some((
        repo.join(format!("{name}.pkl.gz")),
        repo.join(format!("{name}.v1.rkyv")),
    ))
}

fn same_file(left: &Path, right: &Path) -> bool {
    matches!((fs::canonicalize(left), fs::canonicalize(right)), (Ok(a), Ok(b)) if a == b)
}

fn prepare(config: &Path, private: &Path, repo: &Path) -> Result<()> {
    let contents = fs::read_to_string(config)?;
    let doc = document(&contents)?;
    let selected_profile = profile(&doc)?;
    let original = fs::canonicalize(config)?;
    let original_hash = sha256_hex(contents.as_bytes());
    let context_path = sidecar(private, ".cache-install.json");
    let mut replace_completed = false;
    if private.exists() || context_path.exists() {
        nonempty(private)?;
        let context = json(&context_path)?;
        if field(&context, "config")? != text(&original)?
            || field(&context, "source_root")? != text(repo)?
        {
            return Err(invalid(
                "private config belongs to another source root/config; retained unchanged",
            ));
        }
        if field(&context, "config_sha256")? == original_hash {
            if profile(&document(&fs::read_to_string(private)?)?)? != selected_profile {
                return Err(invalid(
                    "pending builder profile changed; private config retained unchanged",
                ));
            }
            println!(
                "Resuming the existing private builder config; its completed updates are retained."
            );
            return Ok(());
        }
        // The planned runtime digest is written before config publication, so
        // a crash immediately after the rename is recoverable. Once a prior
        // publication is complete, explicit subsequent user edits are inputs
        // for a fresh builder view; pending edits are never guessed away.
        replace_completed = context.get("published_sha256").and_then(Json::as_str)
            == Some(&original_hash)
            || context.get("state").and_then(Json::as_str) == Some("complete");
        if !replace_completed {
            return Err(invalid(
                "original config changed during a pending installation; preserve/reconcile its private config first",
            ));
        }
    }
    let mut updates = BTreeMap::new();
    for key in KEYS {
        let Some(value) = setting(&doc, key)? else {
            continue;
        };
        let native = resolve(value, repo);
        if native.extension().and_then(|s| s.to_str()) != Some("rkyv") {
            continue;
        }
        nonempty(&native)?;
        let source = if let Some((source, archive)) =
            bundle(key, repo, &selected_profile).filter(|(_, archive)| same_file(archive, &native))
        {
            // Repository bundles need existence/nonzero only. Actual mapped
            // reads still enforce the format's memory-safety contract.
            nonempty(&archive)?;
            source
        } else {
            let record = json(&sidecar(&native, ".source.json"))?;
            if field(&record, "key")? != key
                || field(&record, "profile")? != selected_profile
                || !same_file(Path::new(field(&record, "archive")?), &native)
            {
                return Err(invalid(format!(
                    "{key}: native cache receipt does not match this path/profile"
                )));
            }
            PathBuf::from(field(&record, "source")?)
        };
        nonempty(&source)?;
        updates.insert(key.to_owned(), text(&source)?.to_owned());
    }
    let context = receipt(&[
        ("config", text(&original)?),
        ("config_sha256", &original_hash),
        ("source_root", text(repo)?),
        ("profile", &selected_profile),
    ])?;
    let builder = rewrite(&contents, &updates)?;
    if replace_completed {
        atomic_write(private, builder.as_bytes())?;
        atomic_write(&context_path, &context)?;
    } else {
        write_new(private, builder.as_bytes())?;
        if let Err(error) = write_new(&context_path, &context) {
            let _ = fs::remove_file(private);
            return Err(error);
        }
    }
    println!("Builder config prepared; runtime config remains unchanged.");
    Ok(())
}

fn check_context(config: &Path, private: &Path, repo: &Path) -> Result<Json> {
    let context = json(&sidecar(private, ".cache-install.json"))?;
    if field(&context, "config")? != text(&fs::canonicalize(config)?)?
        || field(&context, "source_root")? != text(repo)?
        || field(&context, "config_sha256")? != sha256_file(config)?
    {
        return Err(invalid(
            "original config/source root changed after prepare; refusing to overwrite it",
        ));
    }
    Ok(context)
}

fn directory(root: &Path, requested: &Path) -> Result<PathBuf> {
    let relative = requested
        .strip_prefix(root)
        .map_err(|_| invalid("native cache destination must be inside cache_root"))?;
    let mut path = root.to_owned();
    for part in relative.components() {
        let std::path::Component::Normal(name) = part else {
            return Err(invalid("cache destination contains a non-normal component"));
        };
        path.push(name);
        match fs::symlink_metadata(&path) {
            Ok(info) if info.is_dir() && !info.file_type().is_symlink() => (),
            Ok(_) => {
                return Err(invalid(format!(
                    "{}: managed cache directory is not a real directory",
                    path.display()
                )));
            }
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => fs::create_dir(&path)?,
            Err(e) => return Err(e.into()),
        }
    }
    Ok(path)
}

fn target(key: &str, source_hash: &str, doc: &ConfigDocument, repo: &Path) -> Result<PathBuf> {
    let root_input = Path::new(required(doc, "cache_root")?);
    if !root_input.is_absolute() {
        return Err(invalid("cache_root must be an existing absolute directory"));
    }
    let root = fs::canonicalize(root_input)?;
    if !root.is_dir() || root.parent().is_none() {
        return Err(invalid("cache_root must be a dedicated existing directory"));
    }
    let family = if key.starts_with("clinvar_") {
        setting(doc, "clinvar_vcf_dir")?
            .map(|v| resolve(v, repo))
            .unwrap_or_else(|| root_input.join("ClinVar"))
    } else if [
        "alphamissense_tranx_domain_map",
        "pm1_regions_pkl",
        "all_intolerant_domains",
    ]
    .contains(&key)
    {
        setting(doc, "vep_plugins_cachedir")?
            .map(|v| resolve(v, repo))
            .unwrap_or_else(|| root_input.join("VEP_plugin_cache"))
            .join("AlphaMissense")
    } else {
        root_input.join(if key == "clingen_map" {
            "ClinGen"
        } else {
            "InterPro"
        })
    };
    // A user may select a symlink for cache_root itself, but not redirect a
    // managed child. Translate the one approved root alias before checking.
    let family = if let Ok(relative) = family.strip_prefix(root_input) {
        root.join(relative)
    } else {
        family
    };
    let family = directory(&root, &family)?;
    // Source-specific names preserve the previous runtime profile during a
    // refresh. The full digest is verified in the adjacent receipt.
    let selected_profile = profile(doc)?.replace('/', ".");
    Ok(family.join(format!(
        "{key}.{selected_profile}.{}.v1.rkyv",
        &source_hash[..16]
    )))
}

fn validate(key: &str, path: &Path) -> Result<()> {
    let path = text(path)?;
    match key {
        "clinvar_aa_stat" => {
            caches::ClinvarAaCache::open(path)?;
        }
        "clinvar_splice_stat" => {
            caches::clinvar_splice_stat::ClinvarSpliceCache::open(path)?;
        }
        "alphamissense_tranx_domain_map" => {
            caches::transcript_domains::TranscriptDomainCache::open(path)?;
        }
        "clinvar_patho_exon_af_stat" => {
            caches::clinvar_exon_af::ClinvarExonAfCache::open(path)?;
        }
        "clinvar_gene_stat" => {
            caches::clinvar_gene_stat::ClinvarGeneCache::open(path)?;
        }
        "clinvar_patho_af_stat" => {
            caches::clinvar_gene_af::ClinvarGeneAfCache::open(path)?;
        }
        "pm1_regions_pkl" => {
            caches::pm1_regions::Pm1RegionCache::open(path)?;
        }
        "all_intolerant_domains" => {
            caches::intolerant_domains::IntolerantDomainCache::open(path)?;
        }
        "interpro_mapping_pickle" => {
            caches::interpro_mapping::InterproCache::open(path)?;
        }
        "clingen_map" => {
            caches::clingen_map::ClingenCache::open(path)?;
        }
        _ => return Err(invalid("unknown native cache")),
    }
    Ok(())
}

fn install_one(
    key: &str,
    source: &Path,
    doc: &ConfigDocument,
    repo: &Path,
    selected_profile: &str,
) -> Result<PathBuf> {
    nonempty(source)?;
    if let Some((pickle, native)) =
        bundle(key, repo, selected_profile).filter(|(pickle, _)| same_file(pickle, source))
    {
        nonempty(&pickle)?;
        nonempty(&native)?;
        println!("{key}: reuse repository-bundled native cache");
        return Ok(native);
    }
    // This fingerprint identifies the bytes converted; it does not check a
    // repository-bundled pickle against an expected publisher checksum.
    let source_hash = sha256_file(source)?;
    let output = target(key, &source_hash, doc, repo)?;
    let record_path = sidecar(&output, ".source.json");
    if output.exists() || record_path.exists() {
        // Do not overwrite a supplied/unregistered file or silently adopt a
        // partial archive. Every reusable managed archive has a closed receipt.
        let record = json(&record_path)?;
        if field(&record, "key")? != key
            || field(&record, "format")? != "PRIVARKV-v1"
            || field(&record, "source_sha256")? != source_hash
            || field(&record, "profile")? != selected_profile
        {
            return Err(invalid(format!(
                "{}: conflicting cache receipt; original retained",
                output.display()
            )));
        }
        nonempty(&output)?;
        if sha256_file(&output)? != field(&record, "archive_sha256")? {
            return Err(invalid(format!(
                "{}: archive checksum differs from completed conversion; original retained",
                output.display()
            )));
        }
        validate(key, &output)?;
        let source_path = fs::canonicalize(source)?;
        if field(&record, "source")? != text(&source_path)? {
            let mut record = record;
            record.set("source", Json::Str(text(&source_path)?.to_owned()));
            atomic_write(
                &record_path,
                (json_text::dumps_indented(&record).map_err(invalid)? + "\n").as_bytes(),
            )?;
        }
        println!("{key}: reuse completed native conversion");
        return Ok(output);
    }
    caches::dispatch(&[
        key.to_owned(),
        "--pickle".to_owned(),
        text(source)?.to_owned(),
        "--output".to_owned(),
        text(&output)?.to_owned(),
    ])?;
    if sha256_file(source)? != source_hash {
        return Err(invalid(format!(
            "{} changed while converting; runtime config was not published",
            source.display()
        )));
    }
    let archive_hash = sha256_file(&output)?;
    let canonical_source = fs::canonicalize(source)?;
    let source_bytes = fs::metadata(source)?.len().to_string();
    let archive_bytes = fs::metadata(&output)?.len().to_string();
    atomic_write(
        &record_path,
        &receipt(&[
            ("format", "PRIVARKV-v1"),
            ("key", key),
            ("profile", selected_profile),
            ("source", text(&canonical_source)?),
            ("source_sha256", &source_hash),
            ("source_bytes", &source_bytes),
            ("archive", text(&output)?),
            ("archive_sha256", &archive_hash),
            ("archive_bytes", &archive_bytes),
        ])?,
    )?;
    Ok(output)
}

fn publish(config: &Path, private: &Path, repo: &Path) -> Result<()> {
    let context = check_context(config, private, repo)?;
    let contents = fs::read_to_string(private)?;
    let doc = document(&contents)?;
    let selected_profile = profile(&doc)?;
    if selected_profile != field(&context, "profile")? {
        return Err(invalid("assembly/VEP profile changed during installation"));
    }
    let root = fs::canonicalize(required(&doc, "cache_root")?)?;
    if !root.is_dir() || root.parent().is_none() {
        return Err(invalid("cache_root must be a dedicated directory"));
    }
    let lock_path = root.join(".priva-native-caches.lock");
    if fs::symlink_metadata(&lock_path).is_ok_and(|m| !m.is_file() || m.file_type().is_symlink()) {
        return Err(invalid("native cache lock is not a regular file"));
    }
    let lock = OpenOptions::new()
        .create(true)
        .truncate(false)
        .read(true)
        .write(true)
        .open(&lock_path)?;
    lock.try_lock().map_err(|e| {
        invalid(format!(
            "another native-cache publisher holds cache_root: {e}"
        ))
    })?;
    let mut updates = BTreeMap::new();
    for key in KEYS {
        let Some(value) = setting(&doc, key)? else {
            continue;
        };
        let source = resolve(value, repo);
        if source.extension().and_then(|s| s.to_str()) == Some("rkyv") {
            return Err(invalid(format!(
                "{key}: builder config still selects a native archive"
            )));
        }
        let output = install_one(key, &source, &doc, repo, &selected_profile)?;
        updates.insert(key.to_owned(), text(&output)?.to_owned());
    }
    let runtime = rewrite(&contents, &updates)?;
    // Builders may take hours. A concurrent user edit must not be discarded.
    check_context(config, private, repo)?;
    if fs::read_to_string(private)? != contents {
        return Err(invalid("builder config changed during conversion"));
    }
    let mut context = context;
    context.set(
        "published_sha256",
        Json::Str(sha256_hex(runtime.as_bytes())),
    );
    let context_path = sidecar(private, ".cache-install.json");
    atomic_write(
        &context_path,
        (json_text::dumps_indented(&context).map_err(invalid)? + "\n").as_bytes(),
    )?;
    atomic_write(config, runtime.as_bytes())?;
    context.set("state", Json::Str("complete".to_owned()));
    atomic_write(
        &context_path,
        (json_text::dumps_indented(&context).map_err(invalid)? + "\n").as_bytes(),
    )?;
    println!(
        "Runtime config published with {} native cache paths; original builder sources retained.",
        updates.len()
    );
    Ok(())
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args == ["--help"] || args == ["-h"] {
        println!("{HELP}");
        return Ok(());
    }
    if args.len() != 7 || !["prepare", "publish"].contains(&args[0].as_str()) {
        return Err(invalid(HELP));
    }
    let mut options = BTreeMap::new();
    for pair in args[1..].chunks_exact(2) {
        if !["--config", "--build-config", "--source-root"].contains(&pair[0].as_str())
            || options
                .insert(pair[0].as_str(), PathBuf::from(&pair[1]))
                .is_some()
        {
            return Err(invalid(HELP));
        }
    }
    let config = options.get("--config").ok_or_else(|| invalid(HELP))?;
    let private = options.get("--build-config").ok_or_else(|| invalid(HELP))?;
    let repo = fs::canonicalize(options.get("--source-root").ok_or_else(|| invalid(HELP))?)?;
    crate::rust_cache::ensure_distinct_paths(config, private)?;
    match args[0].as_str() {
        "prepare" => prepare(config, private, &repo),
        _ => publish(config, private, &repo),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    struct Temp(PathBuf);
    impl Temp {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "priva-cache-install-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir(&path).unwrap();
            Self(path)
        }
        fn config(&self, extra: &str) -> (PathBuf, PathBuf) {
            let config = self.0.join("config.yaml");
            fs::write(&config, format!("# retain comments\nassembly: hg38\nvep_installed_version: 116\ncache_root: {}\n{extra}", self.0.display())).unwrap();
            (config, self.0.join("private.yaml"))
        }
    }
    impl Drop for Temp {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    #[test]
    fn install_roundtrip_restores_builder_inputs_and_reuses_immutable_archive() {
        let tmp = Temp::new();
        let pickle = tmp.0.join("gene.pkl");
        fs::write(&pickle, b"\x80\x04}\x94.").unwrap();
        let (config, private) = tmp.config(&format!(
            "clinvar_patho_af_stat: '{}' # keep\nother: [a, b]\n",
            pickle.display()
        ));
        prepare(&config, &private, &tmp.0).unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        let runtime = fs::read_to_string(&config).unwrap();
        assert!(runtime.contains("# keep\nother: [a, b]"));
        let doc = document(&runtime).unwrap();
        let archive = PathBuf::from(required(&doc, "clinvar_patho_af_stat").unwrap());
        assert!(archive.starts_with(tmp.0.join("ClinVar")));
        let before = fs::metadata(&archive).unwrap().modified().unwrap();
        let private2 = private.clone();
        prepare(&config, &private2, &tmp.0).unwrap();
        assert_eq!(
            setting(
                &document(&fs::read_to_string(&private2).unwrap()).unwrap(),
                "clinvar_patho_af_stat"
            )
            .unwrap(),
            pickle.to_str()
        );
        publish(&config, &private2, &tmp.0).unwrap();
        assert_eq!(before, fs::metadata(&archive).unwrap().modified().unwrap());
        // A new completed builder result gets a separate generation. The old
        // archive remains usable until the final runtime config rename.
        prepare(&config, &private2, &tmp.0).unwrap();
        fs::write(&pickle, b"\x80\x04}\x94\x8c\x04GENE\x94Ns.").unwrap();
        publish(&config, &private2, &tmp.0).unwrap();
        let updated = fs::read_to_string(&config).unwrap();
        assert_ne!(runtime, updated);
        assert_eq!(before, fs::metadata(&archive).unwrap().modified().unwrap());
    }

    #[test]
    fn refused_publish_keeps_user_edits_and_previous_config() {
        let tmp = Temp::new();
        let (config, private) = tmp.config("");
        prepare(&config, &private, &tmp.0).unwrap();
        fs::write(&config, "user edit\n").unwrap();
        assert!(publish(&config, &private, &tmp.0).is_err());
        assert_eq!(fs::read_to_string(&config).unwrap(), "user edit\n");
        let (config, private2) = tmp.config("clinvar_gene_stat: missing.pkl\n");
        let private2 = private2.with_file_name("missing.yaml");
        let before = fs::read(&config).unwrap();
        prepare(&config, &private2, &tmp.0).unwrap();
        assert!(publish(&config, &private2, &tmp.0).is_err());
        assert_eq!(fs::read(&config).unwrap(), before);
    }

    #[test]
    fn pending_build_config_is_preserved_but_user_profile_changes_are_refused() {
        let tmp = Temp::new();
        let (config, private) = tmp.config("");
        prepare(&config, &private, &tmp.0).unwrap();
        let changed = fs::read_to_string(&private).unwrap() + "new_builder_output: new.vcf.gz\n";
        fs::write(&private, &changed).unwrap();
        prepare(&config, &private, &tmp.0).unwrap();
        assert_eq!(fs::read_to_string(&private).unwrap(), changed);
        let original = fs::read_to_string(&config).unwrap().replace("hg38", "hg19");
        fs::write(&config, original).unwrap();
        assert!(prepare(&config, &private, &tmp.0).is_err());
        assert_eq!(fs::read_to_string(&private).unwrap(), changed);
    }

    #[test]
    fn damaged_managed_archive_is_not_reused_or_overwritten() {
        let tmp = Temp::new();
        let pickle = tmp.0.join("gene.pkl");
        fs::write(&pickle, b"\x80\x04}\x94.").unwrap();
        let (config, private) =
            tmp.config(&format!("clinvar_patho_af_stat: '{}'\n", pickle.display()));
        prepare(&config, &private, &tmp.0).unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        let runtime = fs::read_to_string(&config).unwrap();
        let doc = document(&runtime).unwrap();
        let archive = PathBuf::from(required(&doc, "clinvar_patho_af_stat").unwrap());
        fs::write(&archive, "damaged").unwrap();
        prepare(&config, &private, &tmp.0).unwrap();
        assert!(publish(&config, &private, &tmp.0).is_err());
        assert_eq!(fs::read_to_string(&config).unwrap(), runtime);
        assert_eq!(fs::read_to_string(&archive).unwrap(), "damaged");
    }

    #[test]
    fn small_repository_bundle_uses_only_existence_and_nonzero_checks() {
        let tmp = Temp::new();
        fs::create_dir_all(tmp.0.join("data/InterPro")).unwrap();
        let pickle = tmp.0.join("data/InterPro/Interpro_entry_mapping.pkl.gz");
        let archive = tmp.0.join("data/InterPro/Interpro_entry_mapping.v1.rkyv");
        fs::write(&pickle, "fixture").unwrap();
        fs::write(&archive, "fixture").unwrap();
        let (config, private) =
            tmp.config("interpro_mapping_pickle: data/InterPro/Interpro_entry_mapping.pkl.gz\n");
        prepare(&config, &private, &tmp.0).unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        assert_eq!(
            setting(
                &document(&fs::read_to_string(&config).unwrap()).unwrap(),
                "interpro_mapping_pickle"
            )
            .unwrap(),
            archive.to_str()
        );
        let second = tmp.0.join("second.yaml");
        prepare(&config, &second, &tmp.0).unwrap();
        assert_eq!(
            setting(
                &document(&fs::read_to_string(&second).unwrap()).unwrap(),
                "interpro_mapping_pickle"
            )
            .unwrap(),
            pickle.to_str()
        );
    }

    #[test]
    fn reject_unknown_native_inputs_duplicate_keys_and_redirected_cache_dirs() {
        let tmp = Temp::new();
        fs::write(tmp.0.join("unregistered.rkyv"), "x").unwrap();
        let (config, private) = tmp.config("clinvar_gene_stat: unregistered.rkyv\n");
        assert!(prepare(&config, &private, &tmp.0).is_err());
        assert!(!private.exists());
        assert!(document("assembly: hg38\nassembly: hg19\n").is_err());
        #[cfg(unix)]
        {
            std::os::unix::fs::symlink(&tmp.0, tmp.0.join("redirect")).unwrap();
            assert!(directory(&tmp.0, &tmp.0.join("redirect/out")).is_err());
        }
    }
}
