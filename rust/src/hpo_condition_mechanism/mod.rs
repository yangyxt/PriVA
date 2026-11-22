//! Native `scripts/build_hpo_condition_mechanism_cache.py`: the HPO-framed
//! condition and pathogenic-mechanism cache builder with the script's options.
//!
//! The HPO gene-condition frame and its MONDO merge, the G2P/Orphadata
//! condition mechanisms and penetrance, ClinGen haploinsufficiency, the
//! GoFCards exact alleles with their ClinVar condition links, the deduced
//! recessive LOF, the final scope gate, the structural validation, the
//! provenance check and the atomic publication follow the script. Printed
//! counts, the JSON (compact or `--pretty`; gzip through the same zlib, so
//! the compressed bytes match too) and exit statuses follow the script; an
//! exception the script let escape
//! exits 1 here as `priva: CATEGORY: message`, and its `SystemExit` messages
//! do the same with the category `SystemExit`.
//!
//! Differences:
//! - `--source-root CHECKOUT` is required. The default ClinGen dosage table
//!   comes from it, and the cache records the checkout's Python builder as
//!   its `builder_script` source, so a cache built by either one passes the
//!   other's provenance check.
//! - An input table is decoded completely before its rows are read, so
//!   invalid UTF-8 is reported before an error in an earlier row.
mod build;
pub(crate) mod tsv;
mod validate;

use crate::builder_error;
use crate::cli::{Result, invalid};
use crate::gene_nonlof_mechanism as shared;
use crate::json_text::{self, DumpsOptions, Json};
use crate::python_args::{parse_long_options_converted, require_options};
use crate::python_path::{os_error, path_value};
use std::collections::HashMap;
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};
use std::os::unix::fs::{OpenOptionsExt, PermissionsExt};
use std::path::{Path, PathBuf};

pub(super) use shared::{equal, error, get, py_str, repr, truthy, type_name};

const HELP: &str = "usage: priva hpo-condition-mechanism --source-root CHECKOUT\n\
       [--hpo-assertions HPO_ASSERTIONS] [--mechanism-evidence MECHANISM_EVIDENCE]\n\
       [--clingen-dosage CLINGEN_DOSAGE] [--gofcards-variants GOFCARDS_VARIANTS]\n\
       [--hpo-release HPO_RELEASE] [--mondo-release MONDO_RELEASE] [--output OUTPUT]\n\
       [--validate-only CACHE_JSON] [--pretty]\n\
Build PriVA's HPO-framed condition and pathogenic-mechanism JSON cache; the options,\n\
outputs and exit statuses of scripts/build_hpo_condition_mechanism_cache.py.";

const VALUES: &[&str] = &[
    "source-root",
    "hpo-assertions",
    "mechanism-evidence",
    "clingen-dosage",
    "gofcards-variants",
    "hpo-release",
    "mondo-release",
    "output",
    "validate-only",
];
const FLAGS: &[&str] = &["pretty"];

const BUILDER: &str = "scripts/build_hpo_condition_mechanism_cache.py";
const DEFAULT_CLINGEN_DOSAGE: &str = "data/clingen/gene_dosage_sensitivity.hg38.tsv";

/// The script's `_clean` for a text value.
pub(super) fn clean(value: &str) -> String {
    let text = crate::hpo_penetrance::py_strip(value);
    if matches!(text.to_lowercase().as_str(), "nan" | "none" | "null") {
        String::new()
    } else {
        text.to_owned()
    }
}

/// The script's `_clean` for a parsed JSON value (`str(value)` first).
pub(super) fn clean_json(value: &Json) -> String {
    match value {
        Json::Null => String::new(),
        other => clean(&py_str(other)),
    }
}

/// `Path(path).name.endswith(".gz")`.
fn file_name_ends_gz(path: &Path) -> bool {
    path.file_name()
        .is_some_and(|name| name.to_string_lossy().ends_with(".gz"))
}

/// The bytes of a file, through gzip when asked, with Python's exceptions.
fn read_maybe_gzip(path: &Path, gz: bool) -> Result<Vec<u8>> {
    let file = File::open(path).map_err(|e| os_error(e, path, None))?;
    let mut bytes = Vec::new();
    if gz {
        crate::python_io::PythonGzipReader::new(file)
            .read_to_end(&mut bytes)
            .map_err(|e| match crate::python_io::GzipFailure::of(&e) {
                Some(crate::python_io::GzipFailure::EndOfFile) => error("EOFError", e.to_string()),
                Some(crate::python_io::GzipFailure::Bad(_)) => error("BadGzipFile", e.to_string()),
                Some(crate::python_io::GzipFailure::Zlib(_)) => error("error", e.to_string()),
                None => os_error(e, path, None),
            })?;
    } else {
        let mut file = file;
        file.read_to_end(&mut bytes)
            .map_err(|e| os_error(e, path, None))?;
    }
    Ok(bytes)
}

/// The script's `_source_metadata`.
pub(super) fn source_metadata(path: &Path) -> Result<Json> {
    let source = shared::resolve(path)?;
    let meta = fs::metadata(&source).map_err(|e| os_error(e, &source, None))?;
    Ok(Json::Object(vec![
        ("path".into(), Json::Str(source.display().to_string())),
        ("size_bytes".into(), Json::Int(i128::from(meta.len()))),
        ("sha256".into(), Json::Str(shared::sha256_file(&source)?)),
    ]))
}

/// `json.dump(payload, ensure_ascii=True, indent=..., separators=...,
/// sort_keys=True)` followed by a newline.
fn render(payload: &Json, pretty: bool) -> Result<String> {
    let body = json_text::dumps_with_options(
        payload,
        DumpsOptions {
            indent: pretty.then_some(2),
            sort_keys: true,
            ensure_ascii: true,
            compact: !pretty,
            allow_nan: true,
        },
    )
    .map_err(|e| error("ValueError", e))?;
    Ok(body + "\n")
}

/// `tempfile.mkstemp(prefix=..., suffix=..., dir=...)`: a new file, mode 0600.
fn temporary_file(directory: &Path, prefix: &str, suffix: &str) -> Result<(PathBuf, File)> {
    const LETTERS: &[u8] = b"abcdefghijklmnopqrstuvwxyz0123456789_";
    let mut random =
        File::open("/dev/urandom").map_err(|e| os_error(e, Path::new("/dev/urandom"), None))?;
    for _ in 0..100 {
        let mut bytes = [0u8; 8];
        random
            .read_exact(&mut bytes)
            .map_err(|e| os_error(e, Path::new("/dev/urandom"), None))?;
        let name: String = bytes
            .iter()
            .map(|b| char::from(LETTERS[usize::from(*b) % LETTERS.len()]))
            .collect();
        let path = directory.join(format!("{prefix}{name}{suffix}"));
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .mode(0o600)
            .open(&path)
        {
            Ok(file) => return Ok((path, file)),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(os_error(e, &path, None)),
        }
    }
    Err(error(
        "FileExistsError",
        "[Errno 17] No usable temporary file name found",
    ))
}

/// The script's `write_json_atomic`.
fn write_json_atomic(payload: &Json, output: &Path, pretty: bool) -> Result<()> {
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."))
        .to_owned();
    fs::create_dir_all(&parent).map_err(|e| os_error(e, &parent, None))?;
    let gz = output.to_string_lossy().ends_with(".gz");
    let name = output
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let (temporary, file) = temporary_file(
        &parent,
        &format!(".{name}."),
        if gz { ".tmp.gz" } else { ".tmp" },
    )?;
    let result = (|| -> Result<()> {
        let body = render(payload, pretty)?;
        let bytes = if gz {
            // open_text writes through a TextIOWrapper, which flushes the GzipFile on close.
            crate::python_gzip::compress(body.as_bytes(), 9, true)?
        } else {
            body.into_bytes()
        };
        let mut file = file;
        file.write_all(&bytes)
            .map_err(|e| os_error(e, &temporary, None))?;
        validate::load_and_validate(&temporary)?;
        fs::set_permissions(&temporary, fs::Permissions::from_mode(0o644))
            .map_err(|e| os_error(e, &temporary, None))?;
        fs::rename(&temporary, output).map_err(|e| os_error(e, &temporary, Some(output)))?;
        // Some network filesystems refuse a directory fsync; the rename above
        // still keeps readers from seeing partial content.
        if let Ok(directory) = File::open(&parent) {
            let _ = directory.sync_all();
        }
        Ok(())
    })();
    let _ = fs::remove_file(&temporary);
    result
}

struct Options {
    source_root: PathBuf,
    values: HashMap<String, String>,
    pretty: bool,
}

fn parse_options(args: &[String]) -> Result<Option<Options>> {
    let converted =
        parse_long_options_converted(args, VALUES, FLAGS, |_, _| Ok(())).map_err(invalid)?;
    let Some(occurrences) = converted else {
        println!("{HELP}");
        return Ok(None);
    };
    let mut values: HashMap<String, String> = occurrences.into_iter().collect();
    require_options(&values, &["source-root"]).map_err(invalid)?;
    let pretty = values.remove("pretty").is_some();
    Ok(Some(Options {
        source_root: path_value(&values["source-root"]),
        values,
        pretty,
    }))
}

fn system_exit(message: &str) -> Box<dyn std::error::Error> {
    error("SystemExit", message.to_owned())
}

fn run(options: Options) -> Result<()> {
    let value = |name: &str| options.values.get(name).cloned();
    let builder = options.source_root.join(BUILDER);
    if let Some(cache) = value("validate-only") {
        let provenance = [
            value("hpo-assertions"),
            value("mechanism-evidence"),
            value("clingen-dosage"),
            value("gofcards-variants"),
            value("hpo-release"),
            value("mondo-release"),
        ];
        let counts = if provenance.iter().any(Option::is_some) {
            if provenance.iter().any(Option::is_none) {
                return Err(system_exit(
                    "Provenance validation requires all build inputs and releases",
                ));
            }
            let [
                hpo,
                mechanism,
                clingen,
                gofcards,
                hpo_release,
                mondo_release,
            ] = provenance.map(Option::unwrap);
            validate::validate_provenance(
                &path_value(&cache),
                &[
                    ("builder_script", builder),
                    ("hpo_assertions", path_value(&hpo)),
                    ("mechanism_evidence", path_value(&mechanism)),
                    ("clingen_dosage", path_value(&clingen)),
                    ("gofcards_variants", path_value(&gofcards)),
                ],
                &hpo_release,
                &mondo_release,
            )?
        } else {
            validate::load_and_validate(&path_value(&cache))?
        };
        println!(
            "{}",
            json_text::dumps(&counts, true).map_err(|e| error("ValueError", e))?
        );
        return Ok(());
    }
    let clingen = value("clingen-dosage")
        .filter(|v| !v.is_empty())
        .unwrap_or_else(|| {
            options
                .source_root
                .join(DEFAULT_CLINGEN_DOSAGE)
                .display()
                .to_string()
        });
    let required = [
        ("--hpo-assertions", value("hpo-assertions")),
        ("--mechanism-evidence", value("mechanism-evidence")),
        ("--clingen-dosage", Some(clingen.clone())),
        ("--gofcards-variants", value("gofcards-variants")),
        ("--output", value("output")),
    ];
    let missing: Vec<&str> = required
        .iter()
        .filter(|(_, v)| v.as_deref().unwrap_or("").is_empty())
        .map(|(name, _)| *name)
        .collect();
    if !missing.is_empty() {
        return Err(system_exit(&format!(
            "Missing required build arguments: {}",
            missing.join(", ")
        )));
    }
    let payload = build::build_cache_payload(&build::Inputs {
        builder,
        hpo_assertions: path_value(&value("hpo-assertions").unwrap()),
        mechanism_evidence: path_value(&value("mechanism-evidence").unwrap()),
        clingen_dosage: path_value(&clingen),
        gofcards_variants: path_value(&value("gofcards-variants").unwrap()),
        hpo_release: value("hpo-release").unwrap_or_default(),
        mondo_release: value("mondo-release").unwrap_or_default(),
    })?;
    write_json_atomic(
        &payload,
        &path_value(&value("output").unwrap()),
        options.pretty,
    )?;
    let counts = payload
        .get("_meta")
        .and_then(|meta| meta.get("counts"))
        .cloned()
        .unwrap_or(Json::Null);
    println!(
        "{}",
        json_text::dumps(&counts, true).map_err(|e| error("ValueError", e))?
    );
    Ok(())
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(options) = parse_options(args)? else {
        return Ok(());
    };
    builder_error::finish(run(options))
}
