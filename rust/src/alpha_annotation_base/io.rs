//! Local filesystem and JSON operations for the annotation-base protocol.
//! Absolute spellings, resolved targets and rename-stable identities are distinct.
use super::{error, object, require};
use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use crate::python_path::{absolute, metadata, os_error, path_value, repr, tempfile_parent};
use std::collections::HashSet;
use std::ffi::{CString, OsStr};
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};
use std::os::fd::FromRawFd;
use std::os::unix::ffi::OsStrExt;
use std::os::unix::fs::{MetadataExt, OpenOptionsExt};
use std::os::unix::process::ExitStatusExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

pub(super) fn string(value: impl Into<String>) -> Json {
    Json::Str(value.into())
}
pub(super) fn pstr(path: &Path) -> Json {
    string(path.display().to_string())
}
pub(super) fn path(value: &Json) -> Result<PathBuf> {
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
pub(super) fn type_name(value: &Json) -> &'static str {
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
pub(super) fn truth(value: &Json) -> bool {
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
pub(super) fn get<'a>(value: &'a Json, key: &str) -> Result<&'a Json> {
    match value {
        Json::Object(_) => value.get(key).ok_or_else(|| error("KeyError", repr(key))),
        _ => Err(error(
            "TypeError",
            format!("'{}' object is not subscriptable", type_name(value)),
        )),
    }
}
pub(super) fn optional<'a>(value: &'a Json, key: &str) -> Result<&'a Json> {
    match value {
        Json::Object(_) => Ok(value.get(key).unwrap_or(&Json::Null)),
        _ => Err(error(
            "AttributeError",
            format!("'{}' object has no attribute 'get'", type_name(value)),
        )),
    }
}
pub(super) fn list(value: &Json) -> Result<&Vec<Json>> {
    match value {
        Json::Array(v) => Ok(v),
        _ => Err(error(
            "TypeError",
            format!("'{}' object is not iterable", type_name(value)),
        )),
    }
}
pub(super) fn members(value: &Json) -> Result<&Vec<(String, Json)>> {
    value.as_object().ok_or_else(|| {
        error(
            "AttributeError",
            format!("'{}' object has no attribute 'items'", type_name(value)),
        )
    })
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
pub(super) fn equal(a: &Json, b: &Json) -> bool {
    json_text::py_equal(a, b)
}
pub(super) fn suffix(path: &Path, extra: &str) -> PathBuf {
    PathBuf::from(format!("{}{extra}", path.display()))
}
pub(super) fn exists(path: &Path) -> Result<bool> {
    probe(path, false).map(|m| m.is_some())
}
pub(super) fn lexists(path: &Path) -> Result<bool> {
    probe(path, true).map(|m| m.is_some())
}
pub(super) fn symlink(path: &Path) -> Result<bool> {
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
pub(super) fn resolve(path: &Path) -> Result<PathBuf> {
    resolve_parts(std::env::current_dir()?, path, &mut HashSet::new())
}
pub(super) fn time(seconds: i64, nanos: i64) -> Json {
    Json::Int(seconds as i128 * 1_000_000_000 + nanos as i128)
}
pub(super) fn identity(path: &Path, missing: bool) -> Result<Json> {
    let path = PathBuf::from(absolute(path)?);
    if missing && !lexists(&path)? {
        return Ok(Json::Null);
    }
    let m = metadata(&path)?;
    require(
        m.is_file() && (missing || m.len() > 0),
        format!("Missing/empty/nonregular dependency: {}", path.display()),
    )?;
    Ok(object([
        ("path", pstr(&path)),
        ("target", pstr(&resolve(&path)?)),
        ("inode", Json::Int(m.ino() as i128)),
        ("size", Json::Int(m.len() as i128)),
        ("mtime_ns", time(m.mtime(), m.mtime_nsec())),
        ("ctime_ns", time(m.ctime(), m.ctime_nsec())),
    ]))
}
pub(super) fn digest(path: &Path) -> Result<String> {
    crate::sha256::sha256_file(path).map_err(|e| os_error(e, path, None))
}
pub(super) fn json_digest(value: &Json) -> Result<String> {
    let out = json_text::dumps_with_options(
        value,
        json_text::DumpsOptions {
            sort_keys: true,
            compact: true,
            ..Default::default()
        },
    )
    .map_err(invalid)?;
    Ok(crate::sha256::sha256_hex(out.as_bytes()))
}
pub(super) fn formatted_json(value: &Json) -> Result<String> {
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
pub(super) fn read_text(path: &Path) -> Result<String> {
    let bytes = fs::read(path).map_err(|e| os_error(e, path, None))?;
    Ok(crate::python_io::decode_utf8(&bytes)
        .map_err(|e| error("UnicodeDecodeError", e))?
        .replace("\r\n", "\n")
        .replace('\r', "\n"))
}
pub(super) fn parse(text: &str) -> Result<Json> {
    json_text::parse(text, false).map_err(|e| {
        error(
            if json_text::is_refusal(&e) {
                "UnsupportedMetadata"
            } else {
                "JSONDecodeError"
            },
            e,
        )
    })
}
pub(super) fn read_json(path: &Path) -> Result<Json> {
    parse(&read_text(path)?)
}
pub(super) fn compare_metadata(text: &str, target: &Json) -> Result<Json> {
    json_text::parse_metadata_comparison(text, target).map_err(|e| {
        error(
            if json_text::is_refusal(&e) {
                "UnsupportedMetadata"
            } else {
                "JSONDecodeError"
            },
            e,
        )
    })
}
pub(super) fn mkdir(path: &Path) -> Result<()> {
    fs::create_dir_all(path).map_err(|e| os_error(e, path, None))
}
pub(super) fn rename(source: &Path, target: &Path) -> Result<()> {
    fs::rename(source, target).map_err(|e| os_error(e, source, Some(target)))
}
pub(super) fn unlink(path: &Path) -> Result<()> {
    fs::remove_file(path).map_err(|e| os_error(e, path, None))
}
pub(super) fn stage(parent: &Path, prefix: &str) -> Result<PathBuf> {
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
pub(super) fn atomic_json(path: &Path, value: &Json) -> Result<()> {
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let prefix = format!(
        ".{}.XXXXXX",
        path.file_name().unwrap_or_default().to_string_lossy()
    );
    let mut name = CString::new(tempfile_parent(parent)?.join(prefix).as_os_str().as_bytes())
        .map_err(|_| error("ValueError", "embedded null byte"))?
        .into_bytes_with_nul();
    let fd = unsafe { libc::mkostemp(name.as_mut_ptr().cast(), libc::O_CLOEXEC) };
    let temporary = PathBuf::from(OsStr::from_bytes(&name[..name.len() - 1]));
    if fd < 0 {
        return Err(os_error(std::io::Error::last_os_error(), &temporary, None));
    }
    let result = (|| {
        let mut file = unsafe { File::from_raw_fd(fd) };
        let data = formatted_json(value)?;
        file.write_all(data.as_bytes())
            .map_err(|e| os_error(e, &temporary, None))?;
        file.sync_all().map_err(|e| os_error(e, &temporary, None))?;
        drop(file);
        rename(&temporary, path)
    })();
    if exists(&temporary)? {
        unlink(&temporary)?;
    }
    result
}
pub(super) fn lock(directory: &Path) -> Result<crate::file_lock::ExclusiveLock> {
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
pub(super) fn copy_verified(
    source: &Path,
    target: &Path,
    expected: Option<&Json>,
) -> Result<String> {
    let before = identity(source, false)?;
    let mut input = File::open(source).map_err(|e| os_error(e, source, None))?;
    let mut output = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(target)
        .map_err(|e| os_error(e, target, None))?;
    let mut hash = crate::sha256::Sha256::new();
    let mut buffer = vec![0; 8 * 1024 * 1024];
    loop {
        let n = input
            .read(&mut buffer)
            .map_err(|e| os_error(e, source, None))?;
        if n == 0 {
            break;
        }
        output
            .write_all(&buffer[..n])
            .map_err(|e| os_error(e, target, None))?;
        hash.update(&buffer[..n]);
    }
    output.sync_all().map_err(|e| os_error(e, target, None))?;
    drop(output);
    drop(input);
    require(
        equal(&identity(source, false)?, &before),
        "Annotation source changed while copying",
    )?;
    let hash = hash.finish_hex();
    require(
        expected.is_none_or(|v| matches!(v, Json::Null) || equal(v, &string(&hash))),
        "Verified native source SHA256 differs",
    )?;
    Ok(hash)
}
pub(super) fn which(name: &str) -> Option<PathBuf> {
    let choices = if name.contains('/') {
        vec![PathBuf::from(name)]
    } else {
        std::env::split_paths(&std::env::var_os("PATH").unwrap_or_else(|| "/bin:/usr/bin".into()))
            .map(|p| p.join(name))
            .collect()
    };
    choices.into_iter().find(|p| {
        fs::metadata(p).is_ok_and(|m| !m.is_dir())
            && CString::new(p.as_os_str().as_bytes())
                .is_ok_and(|p| unsafe { libc::access(p.as_ptr(), libc::X_OK) == 0 })
    })
}
pub(super) fn run(argv: &[String], merge_stderr: bool) -> Result<String> {
    let mut command = Command::new(&argv[0]);
    command.args(&argv[1..]).stdout(Stdio::piped());
    let (status, data) = if merge_stderr {
        // One actual pipe preserves subprocess.STDOUT interleaving. Appending
        // independently buffered stderr can select the wrong first version.
        let mut descriptors = [0; 2];
        if unsafe { libc::pipe2(descriptors.as_mut_ptr(), libc::O_CLOEXEC) } != 0 {
            return Err(error(
                "OSError",
                std::io::Error::last_os_error().to_string(),
            ));
        }
        let mut read = unsafe { File::from_raw_fd(descriptors[0]) };
        let write = unsafe { File::from_raw_fd(descriptors[1]) };
        command
            .stdout(Stdio::from(write.try_clone()?))
            .stderr(Stdio::from(write));
        let mut child = command
            .spawn()
            .map_err(|e| os_error(e, Path::new(&argv[0]), None))?;
        drop(command);
        let mut bytes = Vec::new();
        let read_result = read.read_to_end(&mut bytes);
        let status = child.wait()?;
        read_result?;
        (status, bytes)
    } else {
        command.stderr(Stdio::inherit());
        let result = command
            .output()
            .map_err(|e| os_error(e, Path::new(&argv[0]), None))?;
        (result.status, result.stdout)
    };
    if !status.success() {
        let list = format!(
            "[{}]",
            argv.iter().map(|s| repr(s)).collect::<Vec<_>>().join(", ")
        );
        let message = if let Some(code) = status.code() {
            format!("Command '{list}' returned non-zero exit status {code}.")
        } else {
            format!(
                "Command '{list}' died with {}.",
                crate::python_subprocess::signal_repr(status.signal().unwrap_or(1))
            )
        };
        return Err(error("CalledProcessError", message));
    }
    Ok(crate::python_io::decode_utf8(&data)
        .map_err(|e| error("UnicodeDecodeError", e))?
        .replace("\r\n", "\n")
        .replace('\r', "\n"))
}
