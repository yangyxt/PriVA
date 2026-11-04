//! Optional native procedures for the explicitly pinned BioMart crosswalks.
//! Source assets identify the helper and its repository bundles; they are never
//! executed. Existing Python installation/scientific callers are unchanged.
mod cli;
mod ftp;
mod identity;
mod receipt;
mod table;
mod verify;
use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use crate::python_path::{metadata, os_error, repr};
use crate::sha256::Sha256;
pub(crate) use cli::dispatch;
pub(crate) use identity::{
    artifact_route, atomic_json, changed_target, check_receipt, managed_artifact_path,
    mapping_identity,
};
pub(crate) use receipt::artifact_receipt;
use std::fs::{self, File};
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

const ORIGINAL: &[u8] = include_bytes!("../../../scripts/biomart_nm_to_enst.py");
const FTP_ORIGINAL: &[u8] = include_bytes!("../../../scripts/biomart_ftp_nm_to_enst.py");
const GRCH37: &str = "ensembl-grch37-mart116";
const FTP: &str = "ensembl-grch38-mart116-ftp";
const DATASET: &str = "hsapiens_gene_ensembl";

/// Explicit module location, independent of the native executable's location.
/// Construction is the source-admission step, before invoking any procedure.
pub(crate) struct Context {
    helper_root: PathBuf,
}
impl Context {
    pub(crate) fn new(helper_root: &Path) -> Result<Self> {
        for (name, expected) in [
            ("biomart_nm_to_enst.py", ORIGINAL),
            ("biomart_ftp_nm_to_enst.py", FTP_ORIGINAL),
        ] {
            let path = helper_root.join("scripts").join(name);
            let actual = fs::read(&path).map_err(|e| os_error(e, &path, None))?;
            require(
                actual == expected,
                format!(
                    "BioMart helper source asset differs from the native implementation: {}",
                    path.display()
                ),
            )?;
        }
        Ok(Self {
            helper_root: helper_root.to_path_buf(),
        })
    }

    fn source_root(&self) -> Result<PathBuf> {
        // Python derives the default bundle root from __file__.resolve(),
        // including a symlinked helper directory or individual source asset.
        // Do this at the procedure's original use point, not during admission.
        let source = crate::clinvar_contract::resolve_non_strict(
            &self.helper_root.join("scripts/biomart_nm_to_enst.py"),
        )?;
        source
            .parent()
            .and_then(Path::parent)
            .map(Path::to_path_buf)
            .ok_or_else(|| error("IndexError", "1"))
    }
}
fn error(category: &'static str, message: impl Into<String>) -> Box<dyn std::error::Error> {
    crate::builder_error::category(category, invalid(message))
}
fn require(ok: bool, message: impl Into<String>) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(error("ValueError", message))
    }
}
fn object<const N: usize>(items: [(&str, Json); N]) -> Json {
    Json::Object(items.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}
fn string(value: impl Into<String>) -> Json {
    Json::Str(value.into())
}
fn get<'a>(value: &'a Json, key: &str) -> &'a Json {
    value.get(key).unwrap_or(&Json::Null)
}
fn text(value: &Json) -> &str {
    value.as_str().expect("fixed BioMart source string")
}
fn array(value: &Json) -> &[Json] {
    if let Json::Array(items) = value {
        items
    } else {
        &[]
    }
}
fn equal(a: &Json, b: &Json) -> bool {
    json_text::py_equal(a, b)
}
fn facts() -> &'static Json {
    static FACTS: OnceLock<Json> = OnceLock::new();
    FACTS.get_or_init(|| {
        json_text::parse(include_str!("facts.json"), false).expect("reviewed BioMart source facts")
    })
}
fn source(source_id: &Json) -> Result<&'static Json> {
    array(get(facts(), "sources"))
        .iter()
        .find(|s| equal(get(s, "source_id"), source_id))
        .ok_or_else(|| {
            error(
                "ValueError",
                format!(
                    "Unknown explicitly pinned BioMart source: {}",
                    py_str(source_id)
                ),
            )
        })
}
fn source_for_basename(name: &str) -> Option<&'static Json> {
    array(get(facts(), "sources"))
        .iter()
        .find(|s| get(s, "basename").as_str() == Some(name))
}
/// The reviewed facts of one pinned mapping source that the installer
/// prerequisite checks read (`source_for_assembly`/`source_for_basename`).
pub(crate) struct PinnedSource {
    pub(crate) source_id: String,
    pub(crate) relative_path: String,
    pub(crate) aliases: Vec<String>,
}
impl PinnedSource {
    fn of(selected: &Json) -> Self {
        Self {
            source_id: text(get(selected, "source_id")).to_owned(),
            relative_path: text(get(selected, "relative_path")).to_owned(),
            aliases: array(get(selected, "aliases"))
                .iter()
                .map(|alias| text(alias).to_owned())
                .collect(),
        }
    }
}
/// `source_for_assembly(assembly, vep_release)`; `assembly` is `str()` of the
/// configured value, which no non-string configuration value can turn into an alias.
pub(crate) fn pinned_source_for_assembly(
    assembly: &str,
    release: Option<&Json>,
) -> Result<PinnedSource> {
    cli::source_for_assembly(assembly, release).map(PinnedSource::of)
}
/// The pinned source of a `PINNED_BASENAMES` file name, or `None` for any other name.
pub(crate) fn pinned_source_for_basename(basename: &str) -> Option<PinnedSource> {
    source_for_basename(basename).map(PinnedSource::of)
}
/// `verify_pinned(mapping, assembly)`, discarding the returned record.
pub(crate) fn verify_pinned_mapping(mapping: &Path, assembly: Option<&str>) -> Result<()> {
    verify::pinned(mapping, assembly).map(|_| ())
}
fn name(path: &Path) -> String {
    use std::path::Component;
    match path.components().next_back() {
        Some(Component::Normal(n)) => n.to_string_lossy().into_owned(),
        Some(Component::ParentDir) => "..".into(),
        _ => String::new(),
    }
}
fn provenance_path(path: &Path) -> PathBuf {
    PathBuf::from(format!("{}.provenance.json", path.display()))
}
fn sha256(path: &Path) -> Result<String> {
    let mut input = File::open(path).map_err(|e| os_error(e, path, None))?;
    let mut digest = Sha256::new();
    let mut buffer = vec![0; 1024 * 1024];
    loop {
        let n = input
            .read(&mut buffer)
            .map_err(|e| os_error(e, path, None))?;
        if n == 0 {
            break;
        }
        digest.update(&buffer[..n]);
    }
    Ok(digest.finish_hex())
}
fn hash(bytes: &[u8]) -> String {
    let mut h = Sha256::new();
    h.update(bytes);
    h.finish_hex()
}
fn parse_json(raw: &str) -> Result<Json> {
    json_text::parse(raw, false).map_err(json_failure)
}
fn json_failure(e: String) -> Box<dyn std::error::Error> {
    error(
        if json_text::is_refusal(&e) {
            "RuntimeError"
        } else if e.starts_with("Exceeds the limit (4300 digits)") {
            "ValueError"
        } else {
            "JSONDecodeError"
        },
        e,
    )
}
fn read_json_text(path: &Path) -> Result<String> {
    // TextIOWrapper.read() decodes the complete document before JSON admission.
    let raw = fs::read(path).map_err(|e| os_error(e, path, None))?;
    let text = crate::python_io::decode_utf8(&raw).map_err(|e| error("UnicodeDecodeError", e))?;
    Ok(text.replace("\r\n", "\n").replace('\r', "\n"))
}
fn read_json(path: &Path) -> Result<Json> {
    parse_json(&read_json_text(path)?)
}
fn exists(path: &Path) -> Result<bool> {
    match fs::metadata(path) {
        Ok(_) => Ok(true),
        Err(e)
            if matches!(
                e.raw_os_error(),
                Some(libc::ENOENT | libc::ENOTDIR | libc::EBADF | libc::ELOOP)
            ) =>
        {
            Ok(false)
        }
        Err(e) => Err(os_error(e, path, None)),
    }
}
fn is_link(path: &Path) -> Result<bool> {
    match fs::symlink_metadata(path) {
        Ok(s) => Ok(s.file_type().is_symlink()),
        Err(e)
            if matches!(
                e.raw_os_error(),
                Some(libc::ENOENT | libc::ENOTDIR | libc::EBADF | libc::ELOOP)
            ) =>
        {
            Ok(false)
        }
        Err(e) => Err(os_error(e, path, None)),
    }
}
fn truth(value: &Json) -> bool {
    match value {
        Json::Null => false,
        Json::Bool(v) => *v,
        Json::Int(v) => *v != 0,
        Json::Float(v) => *v != 0.0,
        Json::Str(v) => !v.is_empty(),
        Json::Array(v) => !v.is_empty(),
        Json::Object(v) => !v.is_empty(),
    }
}
fn py_repr(value: &Json) -> String {
    match value {
        Json::Null => "None".into(),
        Json::Bool(v) => (if *v { "True" } else { "False" }).into(),
        Json::Int(v) => v.to_string(),
        Json::Float(v) => json_text::python_float_repr(*v),
        Json::Str(v) => repr(v),
        Json::Array(v) => format!("[{}]", v.iter().map(py_repr).collect::<Vec<_>>().join(", ")),
        Json::Object(v) => format!(
            "{{{}}}",
            v.iter()
                .map(|(k, v)| format!("{}: {}", repr(k), py_repr(v)))
                .collect::<Vec<_>>()
                .join(", ")
        ),
    }
}
fn py_str(value: &Json) -> String {
    value
        .as_str()
        .map(str::to_owned)
        .unwrap_or_else(|| py_repr(value))
}
fn hex64(value: &Json) -> bool {
    value.as_str().is_some_and(|s| {
        s.len() == 64
            && s.bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    })
}
fn exact_keys(value: &Json, keys: &[&str]) -> bool {
    value
        .as_object()
        .is_some_and(|v| v.len() == keys.len() && v.iter().all(|(k, _)| keys.contains(&k.as_str())))
}
fn snapshot_identity(path: &Path, source_id: &Json) -> Result<Json> {
    let selected = source(source_id)?;
    let observed = object([
        ("bytes", Json::Int(metadata(path)?.len() as i128)),
        ("sha256", string(sha256(path)?)),
    ]);
    let expected = object([
        ("bytes", get(selected, "expected_bytes").clone()),
        ("sha256", get(selected, "expected_sha256").clone()),
    ]);
    require(
        equal(&observed, &expected),
        format!(
            "Pinned BioMart default snapshot size/SHA256 mismatch: {}; retain the source response and recover explicitly",
            path.display()
        ),
    )?;
    Ok(observed)
}
