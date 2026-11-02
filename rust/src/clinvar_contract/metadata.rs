//! Source-bound AA/splice completion receipts. Filesystem metadata checks are
//! ordered exactly as the Python procedure; payload bytes are never hashed.
use super::{error, object, required};
use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use crate::python_path::{absolute, is_file, os_error, path_value, tempfile_parent};
use std::collections::{HashMap, HashSet};
use std::ffi::{CString, OsStr};
use std::fs::{self, File};
use std::io::Write;
use std::os::fd::FromRawFd;
use std::os::unix::ffi::OsStrExt;
use std::os::unix::fs::MetadataExt;
use std::path::{Path, PathBuf};

fn path(value: &str) -> Result<PathBuf> {
    if value.contains('\0') {
        return Err(error("ValueError", "embedded null byte"));
    }
    Ok(path_value(value))
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
fn is_symlink(path: &Path) -> Result<bool> {
    match fs::symlink_metadata(path) {
        Ok(value) => Ok(value.file_type().is_symlink()),
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

// Path.resolve(strict=False) resolves existing symlink components while keeping
// missing tails. canonicalize() alone would reject those permitted output paths.
fn resolve_parts(mut base: PathBuf, tail: &Path, active: &mut HashSet<PathBuf>) -> Result<PathBuf> {
    if tail.is_absolute() {
        base = PathBuf::from("/");
    }
    for part in tail.as_os_str().as_bytes().split(|b| *b == b'/') {
        if part.is_empty() || part == b"." {
            continue;
        }
        if part == b".." {
            base.pop();
            if base.as_os_str().is_empty() {
                base.push("/");
            }
            continue;
        }
        let candidate = base.join(OsStr::from_bytes(part));
        match fs::symlink_metadata(&candidate) {
            Ok(stat) if stat.file_type().is_symlink() => {
                if !active.insert(candidate.clone()) {
                    return Err(error(
                        "RuntimeError",
                        format!(
                            "Symlink loop from {}",
                            crate::python_path::repr(&candidate.display().to_string())
                        ),
                    ));
                }
                let target =
                    fs::read_link(&candidate).map_err(|e| os_error(e, &candidate, None))?;
                base = resolve_parts(base, &target, active)?;
                active.remove(&candidate);
            }
            // Python's non-strict realpath retains this component after an OS
            // lookup failure; the later exists/stat stage owns that failure.
            _ => base = candidate,
        }
    }
    Ok(base)
}
fn resolve(path: &Path) -> Result<PathBuf> {
    resolve_parts(std::env::current_dir()?, path, &mut HashSet::new())
}

pub(super) fn distinct(protected: &[String], metadata: &[String]) -> Result<()> {
    let mut protected_paths = Vec::new();
    let mut resolved = Vec::new();
    for value in protected {
        let p = path(value)?;
        resolved.push(resolve(&p)?);
        protected_paths.push(p);
    }
    if resolved.iter().collect::<HashSet<_>>().len() != resolved.len() {
        return Err(error(
            "ValueError",
            "ClinVar pair input/output paths must be distinct",
        ));
    }
    let mut metadata_paths = Vec::new();
    for value in metadata {
        let p = path(value)?;
        if is_symlink(&p)? || resolved.contains(&resolve(&p)?) {
            return Err(error(
                "ValueError",
                "ClinVar pair metadata must not alias an input or output",
            ));
        }
        resolved.push(resolve(&p)?);
        metadata_paths.push(p);
    }
    let mut existing = Vec::new();
    for p in protected_paths.iter().chain(&metadata_paths) {
        if exists(p)? {
            existing.push(p);
        }
    }
    // Preserve the source's pairwise samefile/stat sequence. Do not collapse
    // this into a hash of cached metadata, which changes error precedence.
    for (index, p) in existing.iter().enumerate() {
        for other in &existing[index + 1..] {
            let a = fs::metadata(p).map_err(|e| os_error(e, p, None))?;
            let b = fs::metadata(other).map_err(|e| os_error(e, other, None))?;
            if a.dev() == b.dev() && a.ino() == b.ino() {
                return Err(error(
                    "ValueError",
                    "ClinVar pair paths must not share an existing file",
                ));
            }
        }
    }
    Ok(())
}

pub(super) fn file(value: &str) -> Result<Json> {
    let p = path(value)?;
    let absolute = absolute(&p)?;
    let p = Path::new(&absolute);
    let stat = fs::metadata(p).map_err(|e| os_error(e, p, None))?;
    if !is_file(p)? || stat.len() == 0 {
        return Err(error(
            "ValueError",
            format!("Missing/non-file/empty ClinVar pair dependency: {absolute}"),
        ));
    }
    Ok(object([
        ("path", Json::Str(absolute.clone())),
        ("target", Json::Str(resolve(p)?.display().to_string())),
        ("inode", Json::Int(stat.ino() as i128)),
        ("size", Json::Int(stat.len() as i128)),
        (
            "mtime_ns",
            Json::Int(stat.mtime() as i128 * 1_000_000_000 + stat.mtime_nsec() as i128),
        ),
        (
            "ctime_ns",
            Json::Int(stat.ctime() as i128 * 1_000_000_000 + stat.ctime_nsec() as i128),
        ),
    ]))
}
fn source_hash(value: &str) -> Result<Json> {
    let p = path(value)?;
    Ok(Json::Str(
        crate::sha256::sha256_file(&p).map_err(|e| os_error(e, &p, None))?,
    ))
}
pub(super) fn inputs(vcf: &str, producer: &str, contract: &str) -> Result<Json> {
    let mut indexes = Vec::new();
    for suffix in [".csi", ".tbi"] {
        let p = path(&format!("{vcf}{suffix}"))?;
        if exists(&p)? {
            indexes.push(p.display().to_string());
        }
    }
    if indexes.is_empty() {
        return Err(error(
            "ValueError",
            "ClinVar AA/splice input requires an index",
        ));
    }
    let vcf = file(vcf)?;
    let indexes = indexes
        .iter()
        .map(|p| file(p))
        .collect::<Result<Vec<_>>>()?;
    let producer = source_hash(producer)?;
    let contract = source_hash(contract)?;
    Ok(object([
        ("vcf", vcf),
        ("indexes", Json::Array(indexes)),
        ("producer_sha256", producer),
        ("contract_sha256", contract),
    ]))
}

fn read_json(value: &str, target: &Json) -> Result<Json> {
    let p = path(value)?;
    let bytes = fs::read(&p).map_err(|e| os_error(e, &p, None))?;
    let text = crate::python_io::decode_utf8(&bytes).map_err(|e| error("UnicodeDecodeError", e))?;
    json_text::parse_metadata_comparison(&text.replace("\r\n", "\n").replace('\r', "\n"), target)
        .map_err(|e| {
            if e == "JSON nested too deeply" {
                error(
                    "RecursionError",
                    "maximum recursion depth exceeded while decoding a JSON value",
                )
            } else if e.starts_with("Exceeds the limit (4300 digits)") {
                error("ValueError", e)
            } else {
                error("JSONDecodeError", e)
            }
        })
}
fn publish(value: &str, data: &Json) -> Result<()> {
    let target = path(value)?;
    let parent = target
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let name = target.file_name().unwrap_or_default().to_string_lossy();
    let template = tempfile_parent(parent)?.join(format!(".{name}.XXXXXX"));
    let mut template = CString::new(template.as_os_str().as_bytes())
        .map_err(|_| error("ValueError", "embedded null byte"))?
        .into_bytes_with_nul();
    // Exclusive 0600 temporary creation and replacement reproduce mkstemp and
    // os.replace. Only the unpredictable temporary basename is implementation-specific.
    let fd = unsafe { libc::mkostemp(template.as_mut_ptr().cast(), libc::O_CLOEXEC) };
    let temporary = PathBuf::from(OsStr::from_bytes(&template[..template.len() - 1]));
    if fd < 0 {
        return Err(os_error(std::io::Error::last_os_error(), &temporary, None));
    }
    let result = (|| {
        let mut stream = unsafe { File::from_raw_fd(fd) };
        let text = json_text::dumps(data, true).map_err(invalid)? + "\n";
        stream
            .write_all(text.as_bytes())
            .map_err(|e| os_error(e, &temporary, None))?;
        drop(stream);
        fs::rename(&temporary, &target).map_err(|e| os_error(e, &temporary, Some(&target)))?;
        Ok(())
    })();
    if exists(&temporary)? {
        fs::remove_file(&temporary).map_err(|e| os_error(e, &temporary, None))?;
    }
    result
}

pub(super) fn pair(command: &str, options: &HashMap<String, String>) -> Result<()> {
    let vcf = required(options, "vcf")?;
    let producer = required(options, "producer")?;
    let contract = required(options, "contract-source")?;
    let aa = required(options, "aa")?;
    let splice = required(options, "splice")?;
    let receipt = required(options, "receipt")?;
    let captured = options.get("captured").filter(|s| !s.is_empty());
    let mut protected = vec![
        vcf.to_owned(),
        format!("{vcf}.tbi"),
        format!("{vcf}.csi"),
        producer.to_owned(),
        contract.to_owned(),
        aa.to_owned(),
        splice.to_owned(),
    ];
    if let Some(config) = options.get("config").filter(|s| !s.is_empty()) {
        protected.push(config.clone());
    }
    let mut metadata = vec![receipt.to_owned()];
    if let Some(captured) = captured {
        metadata.push(captured.clone());
    }
    distinct(&protected, &metadata)?;
    let inputs = inputs(vcf, producer, contract)?;
    if command == "pair-capture" {
        let captured =
            captured.ok_or_else(|| error("ValueError", "pair-capture requires --captured"))?;
        return publish(captured, &inputs);
    }
    let mut state = object([
        ("schema", Json::Int(1)),
        (
            "identity_policy",
            Json::Str("path-target-inode-size-mtime-ctime; source-content".into()),
        ),
        ("inputs", inputs.clone()),
        ("aa", file(aa)?),
        ("splice", file(splice)?),
    ]);
    if command == "pair-check" {
        let previous = read_json(receipt, &state)?;
        if !matches!(previous, Json::Object(_)) {
            let name = match previous {
                Json::Null => "NoneType",
                Json::Bool(_) => "bool",
                Json::Int(_) => "int",
                Json::Float(_) => "float",
                Json::Str(_) => "str",
                Json::Array(_) => "list",
                Json::Object(_) => unreachable!(),
            };
            return Err(error(
                "AttributeError",
                format!("'{name}' object has no attribute 'get'"),
            ));
        }
        if state.as_object().unwrap().iter().any(|(k, v)| {
            previous
                .get(k)
                .is_none_or(|old| !json_text::py_equal(old, v))
        }) {
            return Err(error(
                "ValueError",
                "ClinVar AA/splice completion association is stale",
            ));
        }
        Ok(())
    } else {
        if captured.is_none()
            || !json_text::py_equal(&read_json(captured.unwrap(), &inputs)?, &inputs)
        {
            return Err(error(
                "ValueError",
                "ClinVar AA/splice inputs changed during generation",
            ));
        }
        state.set("producer_path", Json::Str(absolute(&path(producer)?)?));
        publish(receipt, &state)
    }
}
