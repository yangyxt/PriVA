//! Local filesystem and JSON helpers shared by the AlphaMissense and ClinVar
//! annotation checks and the mapped-artifact publication.
use super::error;
use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use crate::python_path::{os_error, path_value, repr, tempfile_parent};
use std::collections::HashSet;
use std::ffi::{CString, OsStr};
use std::fs::{self, OpenOptions};
use std::os::unix::ffi::OsStrExt;
use std::os::unix::fs::OpenOptionsExt;
use std::path::{Path, PathBuf};

pub(crate) fn string(value: impl Into<String>) -> Json {
    Json::Str(value.into())
}
pub(crate) fn pstr(path: &Path) -> Json {
    string(path.display().to_string())
}
pub(crate) fn path(value: &Json) -> Result<PathBuf> {
    let Some(value) = value.as_str() else {
        return Err(error(
            "TypeError",
            format!(
                "argument should be a str or an os.PathLike object where __fspath__ returns a str, not '{}'",
                type_name(value)
            ),
        ));
    };
    if value.contains('\0') {
        return Err(error("ValueError", "embedded null byte"));
    }
    Ok(path_value(value))
}
pub(crate) fn type_name(value: &Json) -> &'static str {
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
pub(crate) fn truth(value: &Json) -> bool {
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
pub(crate) fn get<'a>(value: &'a Json, key: &str) -> Result<&'a Json> {
    match value {
        Json::Object(_) => value.get(key).ok_or_else(|| error("KeyError", repr(key))),
        _ => Err(error(
            "TypeError",
            format!("'{}' object is not subscriptable", type_name(value)),
        )),
    }
}
pub(super) fn text(value: &Json) -> String {
    match value {
        Json::Str(s) => s.clone(),
        Json::Null => "None".into(),
        Json::Bool(true) => "True".into(),
        Json::Bool(false) => "False".into(),
        Json::Int(n) => n.to_string(),
        Json::Float(n) => json_text::python_float_repr(*n),
        _ => json_text::dumps(value, false).unwrap_or_default(),
    }
}
pub(crate) fn equal(a: &Json, b: &Json) -> bool {
    json_text::py_equal(a, b)
}
pub(crate) fn suffix(path: &Path, extra: &str) -> PathBuf {
    PathBuf::from(format!("{}{extra}", path.display()))
}
pub(crate) fn lexists(path: &Path) -> Result<bool> {
    probe(path, true).map(|m| m.is_some())
}
pub(crate) fn symlink(path: &Path) -> Result<bool> {
    probe(path, true).map(|m| m.is_some_and(|m| m.file_type().is_symlink()))
}
fn probe(path: &Path, link: bool) -> Result<Option<fs::Metadata>> {
    match if link {
        fs::symlink_metadata(path)
    } else {
        fs::metadata(path)
    } {
        Ok(m) => Ok(Some(m)),
        Err(e)
            if matches!(
                e.raw_os_error(),
                Some(libc::ENOENT | libc::ENOTDIR | libc::EBADF | libc::ELOOP)
            ) =>
        {
            Ok(None)
        }
        Err(e) => Err(os_error(e, path, None)),
    }
}
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
            Ok(m) if m.file_type().is_symlink() => {
                if !active.insert(candidate.clone()) {
                    return Err(error(
                        "RuntimeError",
                        format!(
                            "Symlink loop from {}",
                            repr(&candidate.display().to_string())
                        ),
                    ));
                }
                let target =
                    fs::read_link(&candidate).map_err(|e| os_error(e, &candidate, None))?;
                base = resolve_parts(base, &target, active)?;
                active.remove(&candidate);
            }
            _ => base = candidate,
        }
    }
    Ok(base)
}
pub(crate) fn resolve(path: &Path) -> Result<PathBuf> {
    resolve_parts(std::env::current_dir()?, path, &mut HashSet::new())
}
pub(crate) fn time(seconds: i64, nanos: i64) -> Json {
    Json::Int(seconds as i128 * 1_000_000_000 + nanos as i128)
}
pub(crate) fn digest(path: &Path) -> Result<String> {
    crate::sha256::sha256_file(path).map_err(|e| os_error(e, path, None))
}
pub(crate) fn formatted_json(value: &Json) -> Result<String> {
    json_text::dumps_with_options(
        value,
        json_text::DumpsOptions {
            sort_keys: true,
            indent: Some(2),
            ..Default::default()
        },
    )
    .map(|value| value + "\n")
    .map_err(invalid)
}
pub(crate) fn read_text(path: &Path) -> Result<String> {
    let bytes = fs::read(path).map_err(|e| os_error(e, path, None))?;
    Ok(crate::python_io::decode_utf8(&bytes)
        .map_err(|e| error("UnicodeDecodeError", e))?
        .replace("\r\n", "\n")
        .replace('\r', "\n"))
}
pub(crate) fn mkdir(path: &Path) -> Result<()> {
    fs::create_dir_all(path).map_err(|e| os_error(e, path, None))
}
pub(crate) fn rename(source: &Path, target: &Path) -> Result<()> {
    fs::rename(source, target).map_err(|e| os_error(e, source, Some(target)))
}
pub(crate) fn unlink(path: &Path) -> Result<()> {
    fs::remove_file(path).map_err(|e| os_error(e, path, None))
}
pub(crate) fn stage(parent: &Path, prefix: &str) -> Result<PathBuf> {
    let mut name = CString::new(
        tempfile_parent(parent)?
            .join(format!("{prefix}XXXXXX"))
            .as_os_str()
            .as_bytes(),
    )
    .map_err(|_| error("ValueError", "embedded null byte"))?
    .into_bytes_with_nul();
    let result = unsafe { libc::mkdtemp(name.as_mut_ptr().cast()) };
    let path = PathBuf::from(OsStr::from_bytes(&name[..name.len() - 1]));
    if result.is_null() {
        Err(os_error(std::io::Error::last_os_error(), &path, None))
    } else {
        Ok(path)
    }
}
pub(crate) fn lock(directory: &Path) -> Result<crate::file_lock::ExclusiveLock> {
    let file = OpenOptions::new()
        .read(true)
        .custom_flags(libc::O_DIRECTORY)
        .open(directory)
        .map_err(|e| os_error(e, directory, None))?;
    crate::file_lock::ExclusiveLock::acquire(file, crate::file_lock::Wait::Nonblocking).map_err(
        |e| {
            if e.kind() == std::io::ErrorKind::WouldBlock {
                error(
                    "BlockingIOError",
                    "[Errno 11] Resource temporarily unavailable",
                )
            } else {
                os_error(e, directory, None)
            }
        },
    )
}
