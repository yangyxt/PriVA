//! Native `scripts/vep_plugins_release.py`: validate and record the official
//! VEP plugin release in a plugin directory, carry a recorded release switch,
//! and reproduce the approved UTRAnnotator 116 reporting profile.
//!
//! `check`, `install`, `transition`, `clear-transition` and
//! `ensure-utr-reporting` keep the script's checks, receipts (bytes, modes and
//! temporary-file names' prefixes), backups, printed line and exit statuses.
//! `check` exits 1 silently when the receipt does not match. A failure that
//! main() printed as `ERROR: message` exits 1 here as
//! `priva: CATEGORY: message`, CATEGORY being the Python exception class.
//! Differences:
//! - The mode comes first; argparse also accepted it after the options.
//! - The UTRAnnotator 116 reporting module is read from the `data` folder
//!   that `PRIVA_DATA_DIR` names (the script read its own checkout's).
//! - Temporary files carry six random characters after their prefix, not
//!   eight; they never remain after a completed run.
mod difflib;
mod utr;

use crate::builder_error;
use crate::cli::{Result, invalid};
use crate::installer_config_paths::files::{
    copy2, exists, is_dir, is_file, is_link, metadata, read_text, rename, resolve, temporary,
    unlink_missing,
};
use crate::installer_prerequisites::paths::{join, pathlib};
use crate::json_text::{self, DumpsOptions, Json, py_equal};
use crate::python_args::{parse_long_option_occurrences, require_options};
use crate::python_path::repr;
use crate::python_text::is_python_space;
use crate::sha256::Sha256;
use std::collections::HashMap;
use std::error::Error;
use std::fs;
use std::io::{Read, Write};
use std::path::Path;

const HELP: &str = "usage: priva vep-plugins {check,install,transition,clear-transition,ensure-utr-reporting}\n\
       --release RELEASE [--directory DIRECTORY] [--plugins PLUGINS] [--stage STAGE]\n\
       [--transition-file TRANSITION_FILE] [--from-release FROM_RELEASE] [--config-json CONFIG_JSON]\n\
Validate and record official VEP plugin releases without touching other files; the checks,\n\
receipts, printed line and exit statuses of scripts/vep_plugins_release.py. The mode comes first.";
const MODES: [&str; 5] = [
    "check",
    "install",
    "transition",
    "clear-transition",
    "ensure-utr-reporting",
];
const NAMES: [&str; 7] = [
    "directory",
    "release",
    "plugins",
    "stage",
    "transition-file",
    "from-release",
    "config-json",
];
const REPOSITORY: &str = "https://raw.githubusercontent.com/Ensembl/VEP_plugins";
const RECEIPT: &str = ".priva-vep-plugins.json";

fn error(category: &'static str, message: impl Into<String>) -> Box<dyn Error> {
    builder_error::category(category, invalid(message))
}
fn value_error(message: impl Into<String>) -> Box<dyn Error> {
    error("ValueError", message)
}
fn require(condition: bool, message: impl Into<String>) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(value_error(message))
    }
}

fn source(release: &str) -> String {
    format!("{REPOSITORY}/release/{release}")
}

/// `re.fullmatch(r"[1-9][0-9]*", text)`
fn numeric(text: &str) -> bool {
    let bytes = text.as_bytes();
    !bytes.is_empty() && (b'1'..=b'9').contains(&bytes[0]) && bytes.iter().all(u8::is_ascii_digit)
}
/// `re.fullmatch(r"[0-9a-f]{64}", text)`
fn hex64(text: &str) -> bool {
    text.len() == 64
        && text
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
/// `re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", text)`
fn module_name(text: &str) -> bool {
    let bytes = text.as_bytes();
    !bytes.is_empty()
        && bytes[0].is_ascii_alphabetic()
        && bytes
            .iter()
            .all(|b| b.is_ascii_alphanumeric() || *b == b'_')
}
/// `re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*\.pm", text)`
fn plugin_file(text: &str) -> bool {
    text.strip_suffix(".pm").is_some_and(module_name)
}
/// `str(value)` matched by a pattern that only a string or a decimal integer
/// can satisfy (`[1-9][0-9]*` or `[0-9a-f]{64}`).
fn text_of(value: Option<&Json>) -> String {
    match value {
        Some(Json::Str(text)) => text.clone(),
        Some(Json::Int(number)) => number.to_string(),
        _ => String::new(),
    }
}
fn is(value: Option<&Json>, text: &str) -> bool {
    value.is_some_and(|value| py_equal(value, &Json::Str(text.to_owned())))
}
fn is_one(value: Option<&Json>) -> bool {
    value.is_some_and(|value| py_equal(value, &Json::Int(1)))
}
/// Python truth value.
fn truthy(value: Option<&Json>) -> bool {
    match value {
        None | Some(Json::Null) => false,
        Some(Json::Bool(flag)) => *flag,
        Some(Json::Int(number)) => *number != 0,
        Some(Json::Float(number)) => *number != 0.0,
        Some(Json::Str(text)) => !text.is_empty(),
        Some(Json::Array(items)) => !items.is_empty(),
        Some(Json::Object(items)) => !items.is_empty(),
    }
}
fn python_type(value: &Json) -> &'static str {
    match value {
        Json::Null => "NoneType",
        Json::Bool(_) => "bool",
        Json::Int(_) => "int",
        Json::Float(_) => "float",
        Json::Str(_) => "str",
        Json::Array(_) => "list",
        Json::Object(_) => "dict",
    }
}
fn object(items: Vec<(&str, Json)>) -> Json {
    Json::Object(
        items
            .into_iter()
            .map(|(key, value)| (key.to_owned(), value))
            .collect(),
    )
}
fn text(value: &str) -> Json {
    Json::Str(value.to_owned())
}
/// `value.setdefault(key, {})` on an object.
fn setdefault<'a>(value: &'a mut Json, key: &str) -> &'a mut Json {
    let Json::Object(items) = value else {
        unreachable!("setdefault on a non-object")
    };
    let index = match items.iter().position(|(name, _)| name == key) {
        Some(index) => index,
        None => {
            items.push((key.to_owned(), Json::Object(Vec::new())));
            items.len() - 1
        }
    };
    &mut items[index].1
}

/// `str(Path(path).parent)`
fn parent(path: &str) -> String {
    match path.rfind('/') {
        None => ".".to_owned(),
        Some(_) => crate::installer_prerequisites::paths::parent(path),
    }
}

/// `digest(path)`: the SHA-256 of a non-empty regular file that is not a link.
fn digest(path: &str) -> Result<String> {
    let at = Path::new(path);
    if is_link(at)? || !is_file(at)? || metadata(at)?.len() == 0 {
        return Err(value_error(format!(
            "Missing or empty regular plugin file: {path}"
        )));
    }
    Ok(sha256_hex(&read_bytes(path)?))
}
fn read_bytes(path: &str) -> Result<Vec<u8>> {
    let at = Path::new(path);
    let mut file = fs::File::open(crate::installer_config_paths::files::native(at)?)
        .map_err(|cause| crate::installer_config_paths::files::os_error(cause, at, None))?;
    let mut data = Vec::new();
    file.read_to_end(&mut data)
        .map_err(|cause| crate::installer_config_paths::files::os_error(cause, at, None))?;
    Ok(data)
}
fn sha256_hex(data: &[u8]) -> String {
    let mut hash = Sha256::new();
    hash.update(data);
    hash.finish_hex()
}

/// `json.loads(path.read_text())`, with the decode error's class.
fn load_json(path: &str) -> Result<Json> {
    let content = read_text(Path::new(path))?;
    json_text::parse(&content, false).map_err(|message| error("JSONDecodeError", message))
}

/// `load_receipt(directory)`: the receipt object, or `{}` when it cannot be read.
fn load_receipt(directory: &str) -> Json {
    match load_json(&join(directory, RECEIPT)) {
        Ok(value @ Json::Object(_)) => value,
        _ => Json::Object(Vec::new()),
    }
}

/// `write_receipt(path, value)`: indented, sorted JSON published by rename.
fn write_receipt(path: &str, value: &Json) -> Result<()> {
    let body = json_text::dumps_with_options(
        value,
        DumpsOptions {
            indent: Some(2),
            sort_keys: true,
            ..DumpsOptions::default()
        },
    )
    .map_err(value_error)?;
    let (mut file, temporary) = temporary(Path::new(&parent(path)), ".priva-plugin-receipt-")?;
    let written = file
        .write_all(format!("{body}\n").as_bytes())
        .map_err(|cause| crate::installer_config_paths::files::os_error(cause, &temporary, None));
    drop(file);
    let published = written.and_then(|()| rename(&temporary, Path::new(path)));
    let removed = unlink_missing(&temporary);
    published.and(removed)
}

/// `local_patch_matches(directory, receipt, filename)`: any failure is false.
fn local_patch_matches(directory: &str, receipt: &Json, filename: &str) -> bool {
    // A non-mapping receipt field or patch fails .get(); a missing key fails [].
    if !matches!(receipt.get("local_patches"), None | Some(Json::Object(_))) {
        return false;
    }
    let empty = Json::Object(Vec::new());
    let patch = receipt
        .get("local_patches")
        .and_then(|patches| patches.get(filename))
        .unwrap_or(&empty);
    let Json::Object(_) = patch else {
        return false;
    };
    let Some(official) = receipt.get("plugins").and_then(|plugins| match plugins {
        Json::Object(_) => plugins.get(filename),
        _ => None,
    }) else {
        return false;
    };
    let upstream = patch.get("upstream_sha256");
    if !py_equal(upstream.unwrap_or(&Json::Null), official) {
        return false;
    }
    // re.fullmatch on a non-string value raises TypeError.
    let field = |key: &str| patch.get(key).map_or(Some(""), Json::as_str);
    let (Some(installed), Some(recorded_patch)) =
        (field("installed_sha256"), field("patch_sha256"))
    else {
        return false;
    };
    if !hex64(installed) || !hex64(recorded_patch) {
        return false;
    }
    if digest(&join(directory, filename)).ok().as_deref() != Some(installed) {
        return false;
    }
    let file_digest = |key: &str| match patch.get(key) {
        Some(Json::Str(path)) => digest(&pathlib(path)).ok(),
        _ => None,
    };
    let Some(upstream) = upstream else {
        return false;
    };
    match file_digest("original_file") {
        Some(value) if py_equal(upstream, &Json::Str(value.clone())) => {}
        _ => return false,
    }
    if file_digest("patch_file").as_deref() != Some(recorded_patch) {
        return false;
    }
    truthy(patch.get("profile"))
}

/// `filename in receipt.get("local_patches", {})` for any JSON value.
fn listed_as_patched(receipt: &Json, filename: &str) -> Result<bool> {
    match receipt.get("local_patches") {
        None => Ok(false),
        Some(Json::Object(items)) => Ok(items.iter().any(|(name, _)| name == filename)),
        Some(Json::Array(items)) => Ok(items.iter().any(|item| is(Some(item), filename))),
        Some(Json::Str(text)) => Ok(text.contains(filename)),
        Some(other) => Err(error(
            "TypeError",
            format!(
                "argument of type {} is not iterable",
                repr(python_type(other))
            ),
        )),
    }
}

/// `matches(directory, receipt, release, names)`
fn matches(directory: &str, receipt: &Json, release: &str, names: &[String]) -> Result<bool> {
    let Some(Json::Object(_)) = receipt.get("plugins") else {
        return Ok(false);
    };
    if !is_one(receipt.get("schema_version"))
        || !is(receipt.get("release"), release)
        || !is(receipt.get("source"), &source(release))
        || !hex64(&text_of(Some(
            receipt
                .get("plugin_config_sha256")
                .unwrap_or(&Json::Str(String::new())),
        )))
    {
        return Ok(false);
    }
    for name in names {
        let filename = format!("{name}.pm");
        let matched = if listed_as_patched(receipt, &filename)? {
            local_patch_matches(directory, receipt, &filename)
        } else {
            let recorded = receipt
                .get("plugins")
                .and_then(|plugins| plugins.get(&filename));
            let observed = digest(&join(directory, &filename))?;
            is(recorded, &observed)
        };
        if !matched {
            return Ok(false);
        }
    }
    Ok(true)
}

/// `capture_transition(path, config, previous, release)`
fn capture_transition(
    path: &str,
    config: &Json,
    previous: Option<&str>,
    release: &str,
) -> Result<()> {
    let previous = previous.unwrap_or_default();
    require(
        numeric(previous) && previous != release,
        "A release transition needs distinct numeric source and target releases",
    )?;
    require(
        matches!(config, Json::Object(_)),
        "The installer configuration must be a mapping",
    )?;
    let root = config.get("cache_root");
    let directory = if root.is_some_and(|root| !matches!(root, Json::Null) && !is(Some(root), "")) {
        let root = match root {
            Some(Json::Str(root)) if absolute_plain(root) => root,
            _ => {
                return Err(value_error(
                    "cache_root must be an absolute path without spaces or shell characters",
                ));
            }
        };
        let root = crate::python_path::tempfile_parent(Path::new(root))?
            .display()
            .to_string();
        require(
            !root.trim_matches('/').is_empty(),
            "cache_root must name a dedicated cache directory",
        )?;
        let empty = Json::Object(Vec::new());
        let folders = config
            .get("cache_folders")
            .filter(|value| truthy(Some(value)))
            .unwrap_or(&empty);
        require(
            matches!(folders, Json::Object(_)),
            "cache_folders must be a mapping",
        )?;
        let folder = match folders.get("vep_plugin") {
            None => "VEP_plugin",
            Some(Json::Str(folder)) => folder.as_str(),
            Some(_) => "",
        };
        require(
            !folder.is_empty()
                && folder
                    .bytes()
                    .all(|b| b.is_ascii_alphanumeric() || b"_.-".contains(&b))
                && folder != "."
                && folder != "..",
            "cache_folders.vep_plugin must be a single folder name",
        )?;
        join(&pathlib(&root), folder)
    } else {
        match config.get("vep_plugins_dir") {
            Some(value) if truthy(Some(value)) => match value {
                Json::Str(directory) => pathlib(directory),
                other => {
                    return Err(error(
                        "TypeError",
                        format!(
                            "argument should be a str or an os.PathLike object where __fspath__ \
                             returns a str, not {}",
                            repr(python_type(other))
                        ),
                    ));
                }
            },
            _ => ".".to_owned(),
        }
    };
    require(
        directory.starts_with('/'),
        "vep_plugins_dir must resolve to an absolute path before a release switch",
    )?;
    let directory = resolve(Path::new(&directory))?.display().to_string();
    let mut prior = Vec::new();
    if let Ok(entries) = fs::read_dir(&directory) {
        for entry in entries.flatten() {
            let name = entry.file_name().to_string_lossy().into_owned();
            let candidate = join(&directory, &name);
            if plugin_file(&name) && !is_link(Path::new(&candidate))? {
                prior.push((name, text(&digest(&candidate)?)));
            }
        }
    }
    let environment = resolve(Path::new(&parent(path)))?.display().to_string();
    write_receipt(
        path,
        &object(vec![
            ("schema_version", Json::Int(1)),
            ("from_release", text(previous)),
            ("to_release", text(release)),
            ("environment_prefix", text(&environment)),
            ("plugins_directory", text(&directory)),
            ("plugins", Json::Object(prior)),
        ]),
    )
}

/// `re.fullmatch(r"/[A-Za-z0-9_./-]+", root)`
fn absolute_plain(root: &str) -> bool {
    root.len() > 1
        && root.starts_with('/')
        && root
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"_./-".contains(&b))
}

/// `matching_transition(path, directory, release)`: the recorded switch to
/// this release for this environment and plugin directory, else None.
fn matching_transition(path: Option<&str>, directory: &str, release: &str) -> Option<Json> {
    let path = path?;
    let value = load_json(path).ok()?;
    let Json::Object(_) = value else {
        return None;
    };
    let from = value.get("from_release");
    let accepted = is_one(value.get("schema_version"))
        && is(value.get("to_release"), release)
        && !is(from, release)
        && numeric(&text_of(Some(from.unwrap_or(&Json::Str(String::new())))))
        && is(
            value.get("environment_prefix"),
            &resolve(Path::new(&parent(path)))
                .ok()?
                .display()
                .to_string(),
        )
        && is(
            value.get("plugins_directory"),
            &resolve(Path::new(directory)).ok()?.display().to_string(),
        )
        && matches!(value.get("plugins"), Some(Json::Object(_)));
    accepted.then_some(value)
}

/// `re.search(r"(?m)^\s*package\s+NAME\s*;", text)`
fn declares_package(content: &str, name: &str) -> bool {
    let chars: Vec<char> = content.chars().collect();
    let name: Vec<char> = name.chars().collect();
    let keyword: Vec<char> = "package".chars().collect();
    let skip = |mut at: usize| {
        while at < chars.len() && is_python_space(chars[at]) {
            at += 1;
        }
        at
    };
    let literal = |at: usize, word: &[char]| chars.get(at..at + word.len()) == Some(word);
    (0..=chars.len())
        .filter(|&start| start == 0 || chars[start - 1] == '\n')
        .any(|start| {
            let at = skip(start);
            if !literal(at, &keyword) {
                return false;
            }
            let after = at + keyword.len();
            let at = skip(after);
            if at == after || !literal(at, &name) {
                return false;
            }
            let at = skip(at + name.len());
            chars.get(at) == Some(&';')
        })
}

/// `install(directory, release, names, stage, transition_file)`
fn install(
    directory: &str,
    release: &str,
    names: &[String],
    stage: &str,
    transition_file: Option<&str>,
) -> Result<()> {
    let transition = matching_transition(transition_file, directory, release);
    let config_sha256 = digest(&join(stage, "plugin_config.txt"))?;
    let mut hashes: Vec<(String, String)> = Vec::new();
    for name in names {
        let filename = format!("{name}.pm");
        let value = digest(&join(stage, &filename))?;
        hashes.push((filename, value));
    }
    for name in names {
        let path = join(stage, &format!("{name}.pm"));
        if !declares_package(&read_text(Path::new(&path))?, name) {
            return Err(value_error(format!(
                "Downloaded plugin does not declare package {name}: {path}"
            )));
        }
    }
    let previous = load_receipt(directory);
    let empty = Json::Object(Vec::new());
    let no_release = Json::Str(String::new());
    let old_release = previous.get("release").unwrap_or(&no_release);
    let old_plugins = previous.get("plugins").unwrap_or(&empty);
    let local_patches = previous.get("local_patches").unwrap_or(&empty);
    let Json::Object(patches) = local_patches else {
        return Err(value_error(
            "Malformed local plugin patch receipt; installed files were preserved",
        ));
    };
    let downloaded = |filename: &str| {
        hashes
            .iter()
            .find(|(name, _)| name == filename)
            .map(|(_, value)| value.as_str())
    };
    for (filename, patch) in patches {
        let reviewed = local_patch_matches(directory, &previous, filename)
            && is(Some(old_release), release)
            && is(previous.get("source"), &source(release))
            && downloaded(filename).is_some_and(|value| is(patch.get("upstream_sha256"), value));
        if !reviewed {
            return Err(value_error(format!(
                "Recorded local plugin patch needs review before an upstream change: {filename}. \
                 Installed files were preserved; restore the recorded original before changing releases."
            )));
        }
    }
    let mut allow_release_change = false;
    if let (Json::Str(old), Json::Object(old_files)) = (old_release, old_plugins)
        && numeric(old)
        && old != release
    {
        let old_names: Vec<&str> = old_files
            .iter()
            .filter(|(name, _)| plugin_file(name))
            .map(|(name, _)| &name[..name.len() - 3])
            .collect();
        // A failed publication can leave selected files at the new release.
        // Retry only if old files match their recorded hash or this download.
        let metadata_valid = is_one(previous.get("schema_version"))
            && is(previous.get("source"), &source(old))
            && hex64(&text_of(Some(
                previous
                    .get("plugin_config_sha256")
                    .unwrap_or(&Json::Str(String::new())),
            )));
        let mut known_files = Some(true);
        for name in &old_names {
            let filename = format!("{name}.pm");
            let Ok(observed) = digest(&join(directory, &filename)) else {
                known_files = None;
                break;
            };
            let recorded = old_plugins.get(&filename).unwrap_or(&Json::Null);
            if let Json::Array(_) | Json::Object(_) = recorded {
                return Err(error(
                    "TypeError",
                    format!("unhashable type: {}", repr(python_type(recorded))),
                ));
            }
            if !is(Some(recorded), &observed) && downloaded(&filename) != Some(observed.as_str()) {
                known_files = Some(false);
                break;
            }
        }
        if let Some(known_files) = known_files {
            allow_release_change = metadata_valid && !old_names.is_empty() && known_files;
        }
    }
    let carried = |filename: &str| patches.iter().any(|(name, _)| name == filename);
    for (filename, expected) in &hashes {
        let destination = join(directory, filename);
        let at = Path::new(&destination);
        if is_link(at)? {
            return Err(value_error(format!(
                "Refusing to replace a symlinked VEP plugin: {destination}"
            )));
        }
        if !carried(filename)
            && exists(at)?
            && digest(&destination)? != *expected
            && !allow_release_change
            && !is(
                transition
                    .as_ref()
                    .and_then(|value| value.get("plugins"))
                    .and_then(|plugins| plugins.get(filename)),
                &digest(&destination)?,
            )
        {
            return Err(value_error(format!(
                "Existing VEP plugin differs from official release {release}: {destination}. \
                 Choose an empty vep_plugins_dir for this release or explicitly switch the \
                 environment's major release; existing customized or unrecorded plugins were preserved."
            )));
        }
    }
    let mut receipt = object(vec![
        ("schema_version", Json::Int(1)),
        ("release", text(release)),
        ("source", text(&source(release))),
        ("plugin_config_sha256", text(&config_sha256)),
        (
            "plugins",
            Json::Object(
                hashes
                    .iter()
                    .map(|(name, value)| (name.clone(), text(value)))
                    .collect(),
            ),
        ),
    ]);
    if !patches.is_empty() {
        receipt.set("local_patches", local_patches.clone());
    }
    for (filename, expected) in &hashes {
        let destination = join(directory, filename);
        if carried(filename)
            || (exists(Path::new(&destination))? && digest(&destination)? == *expected)
        {
            continue;
        }
        let (file, staged) = temporary(Path::new(directory), ".priva-plugin-")?;
        drop(file);
        let published = copy2(Path::new(&join(stage, filename)), &staged)
            .and_then(|()| rename(&staged, Path::new(&destination)));
        let removed = unlink_missing(&staged);
        published.and(removed)?;
    }
    write_receipt(&join(directory, RECEIPT), &receipt)
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(mode) = args.first() else {
        return Err(invalid(format!(
            "the following arguments are required: mode, --release\n{HELP}"
        )));
    };
    if matches!(mode.as_str(), "-h" | "--help") {
        println!("{HELP}");
        return Ok(());
    }
    if !MODES.contains(&mode.as_str()) {
        return Err(invalid(format!(
            "argument mode: invalid choice: {} (choose from {})",
            repr(mode),
            MODES.map(repr).join(", ")
        )));
    }
    let Some(occurrences) =
        parse_long_option_occurrences(&args[1..], &NAMES, &[]).map_err(invalid)?
    else {
        println!("{HELP}");
        return Ok(());
    };
    let values: HashMap<String, String> = occurrences.into_iter().collect();
    require_options(&values, &["release"]).map_err(invalid)?;
    let option = |name: &str| values.get(name).map(String::as_str);
    let path_option = |name: &str| option(name).map(pathlib);
    let outcome = builder_error::finish(run(
        mode,
        option("release").unwrap_or_default(),
        path_option("directory").as_deref(),
        option("plugins").unwrap_or_default(),
        path_option("stage").as_deref(),
        path_option("transition-file").as_deref(),
        option("from-release"),
        option("config-json"),
    ))?;
    if !outcome {
        std::process::exit(1);
    }
    Ok(())
}

/// main() after argument parsing; false is the silent `check` mismatch.
#[allow(clippy::too_many_arguments)]
fn run(
    mode: &str,
    release: &str,
    directory: Option<&str>,
    plugins: &str,
    stage: Option<&str>,
    transition_file: Option<&str>,
    from_release: Option<&str>,
    config_json: Option<&str>,
) -> Result<bool> {
    require(
        numeric(release),
        "VEP release must be a numeric major release",
    )?;
    if mode == "transition" {
        let Some(path) = transition_file else {
            return Err(value_error(
                "A release transition requires an environment transition file",
            ));
        };
        let config = json_text::parse(
            config_json.filter(|v| !v.is_empty()).unwrap_or("null"),
            false,
        )
        .map_err(|message| error("JSONDecodeError", message))?;
        capture_transition(path, &config, from_release, release)?;
        return Ok(true);
    }
    let Some(directory) = directory else {
        return Err(value_error("A plugin directory is required"));
    };
    if mode == "clear-transition" {
        if let Some(path) = transition_file
            && matching_transition(Some(path), directory, release).is_some_and(|v| truthy(Some(&v)))
        {
            fs::remove_file(path).map_err(|cause| {
                crate::installer_config_paths::files::os_error(cause, Path::new(path), None)
            })?;
        }
        return Ok(true);
    }
    if mode == "ensure-utr-reporting" {
        require(
            is_dir(Path::new(directory))?,
            format!("VEP plugin directory does not exist: {directory}"),
        )?;
        let result = utr::ensure(directory, release)?;
        println!("{}", json_text::dumps(&result, true).map_err(value_error)?);
        return Ok(true);
    }
    let mut names: Vec<String> = Vec::new();
    for name in plugins.split(',') {
        if !names.iter().any(|seen| seen == name) {
            names.push(name.to_owned());
        }
    }
    require(
        names.iter().all(|name| module_name(name)),
        "VEP plugin names must be simple module names",
    )?;
    require(
        is_dir(Path::new(directory))?,
        format!("VEP plugin directory does not exist: {directory}"),
    )?;
    if mode == "check" {
        return matches(directory, &load_receipt(directory), release, &names);
    }
    let Some(stage) = stage else {
        return Err(value_error(
            "Installing plugin files requires a completed staging directory",
        ));
    };
    install(directory, release, &names, stage, transition_file)?;
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn package_declarations_follow_the_multiline_pattern() {
        assert!(declares_package("package SpliceAI;\n", "SpliceAI"));
        assert!(declares_package(
            "# x\n  package\n\tSpliceAI \n;",
            "SpliceAI"
        ));
        assert!(declares_package("\u{a0}package SpliceAI;", "SpliceAI"));
        assert!(!declares_package("package SpliceAI2;", "SpliceAI"));
        assert!(!declares_package("packageSpliceAI;", "SpliceAI"));
        assert!(!declares_package("x package SpliceAI;", "SpliceAI"));
        assert!(!declares_package("x\rpackage SpliceAI;", "SpliceAI"));
    }

    #[test]
    fn patterns_match_like_their_regular_expressions() {
        assert!(numeric("116") && !numeric("0") && !numeric("016") && !numeric(""));
        assert!(plugin_file("UTRAnnotator.pm") && !plugin_file("_x.pm") && !plugin_file(".pm"));
        assert!(absolute_plain("/a/b") && !absolute_plain("/") && !absolute_plain("/a b"));
        assert_eq!(text_of(Some(&Json::Int(115))), "115");
        assert_eq!(text_of(Some(&Json::Bool(true))), "");
    }

    #[test]
    fn admission_errors_exit_two_and_checks_exit_one() {
        let args = |items: &[&str]| items.iter().map(|v| (*v).to_owned()).collect::<Vec<_>>();
        for bad in [
            args(&[]),
            args(&["fetch", "--release", "116"]),
            args(&["check"]),
        ] {
            let failure = dispatch(&bad).unwrap_err();
            assert!(!failure.is::<builder_error::RuntimeError>(), "{failure}");
        }
        let failure = dispatch(&args(&["check", "--release", "v116"])).unwrap_err();
        assert_eq!(
            failure.to_string(),
            "ValueError: VEP release must be a numeric major release"
        );
        let failure = dispatch(&args(&["install", "--release", "116"])).unwrap_err();
        assert_eq!(
            failure.to_string(),
            "ValueError: A plugin directory is required"
        );
    }
}
