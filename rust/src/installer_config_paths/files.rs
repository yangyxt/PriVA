use super::*;
use std::ffi::{CString, OsStr};
use std::fs::{self, File};
use std::io::{Read, Write};
use std::os::fd::FromRawFd;
use std::os::unix::ffi::OsStrExt;
use std::os::unix::fs::{MetadataExt, PermissionsExt};

pub(crate) fn native(path: &Path) -> Result<PathBuf> {
    crate::python_unicode::native(path).map_err(|e| error("UnicodeEncodeError", e))
}
pub(crate) fn os_error(
    cause: std::io::Error,
    path: &Path,
    dest: Option<&Path>,
) -> Box<dyn std::error::Error> {
    let kind = match cause.kind() {
        std::io::ErrorKind::NotFound => "FileNotFoundError",
        std::io::ErrorKind::PermissionDenied => "PermissionError",
        std::io::ErrorKind::AlreadyExists => "FileExistsError",
        std::io::ErrorKind::IsADirectory => "IsADirectoryError",
        std::io::ErrorKind::NotADirectory => "NotADirectoryError",
        _ => "OSError",
    };
    let text = crate::python_path::os_error(cause, path, dest).to_string();
    error(
        kind,
        text.strip_prefix(&format!("{kind}: ")).unwrap_or(&text),
    )
}
pub(crate) fn metadata(path: &Path) -> Result<fs::Metadata> {
    fs::metadata(native(path)?).map_err(|e| os_error(e, path, None))
}
pub(crate) fn resolve(path: &Path) -> Result<PathBuf> {
    Ok(crate::python_unicode::from_native(
        &crate::clinvar_contract::resolve_non_strict(&native(path)?)?,
    ))
}
pub(crate) fn is_file(path: &Path) -> Result<bool> {
    Ok(probe(path, false)?.is_some_and(|m| m.is_file()))
}

fn probe(path: &Path, link: bool) -> Result<Option<fs::Metadata>> {
    // Pathlib's exists/is_file/is_dir/is_symlink suppress ValueError, which
    // includes UnicodeEncodeError. stat and actual writes do not suppress it.
    let Ok(actual) = crate::python_unicode::native(path) else {
        return Ok(None);
    };
    match if link {
        fs::symlink_metadata(&actual)
    } else {
        fs::metadata(&actual)
    } {
        Ok(v) => Ok(Some(v)),
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
pub(crate) fn exists(path: &Path) -> Result<bool> {
    Ok(probe(path, false)?.is_some())
}
pub(crate) fn is_dir(path: &Path) -> Result<bool> {
    Ok(probe(path, false)?.is_some_and(|m| m.is_dir()))
}
pub(crate) fn is_link(path: &Path) -> Result<bool> {
    Ok(probe(path, true)?.is_some_and(|m| m.file_type().is_symlink()))
}
pub(super) fn suffix(path: &Path, value: &str) -> PathBuf {
    let mut text = path.as_os_str().to_owned();
    text.push(value);
    text.into()
}
pub(crate) fn cache_identity(path: &Path) -> Result<Json> {
    let m = metadata(path)?;
    Ok(Json::Array(vec![
        Json::Int(m.ino() as i128),
        Json::Int(m.len() as i128),
        Json::Int(m.mtime() as i128 * 1_000_000_000 + m.mtime_nsec() as i128),
        Json::Int(m.ctime() as i128 * 1_000_000_000 + m.ctime_nsec() as i128),
    ]))
}
pub(crate) fn sha256_file(path: &Path) -> Result<String> {
    let mut input = File::open(native(path)?).map_err(|e| os_error(e, path, None))?;
    let mut h = crate::sha256::Sha256::new();
    let mut block = vec![0; 8 * 1024 * 1024];
    loop {
        let n = input
            .read(&mut block)
            .map_err(|e| os_error(e, path, None))?;
        if n == 0 {
            break;
        }
        h.update(&block[..n]);
    }
    Ok(h.finish_hex())
}
pub(crate) fn read_text(path: &Path) -> Result<String> {
    let bytes = fs::read(native(path)?).map_err(|e| os_error(e, path, None))?;
    Ok(crate::python_io::decode_utf8(&bytes)
        .map_err(|e| error("UnicodeDecodeError", e))?
        .replace("\r\n", "\n")
        .replace('\r', "\n"))
}
pub(crate) fn mkdir(path: &Path) -> Result<()> {
    fs::create_dir_all(native(path)?).map_err(|e| os_error(e, path, None))
}
pub(crate) fn rename(source: &Path, target: &Path) -> Result<()> {
    fs::rename(native(source)?, native(target)?).map_err(|e| os_error(e, source, Some(target)))
}
pub(crate) fn unlink_missing(path: &Path) -> Result<()> {
    match fs::remove_file(native(path)?) {
        Ok(()) => Ok(()),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(e) => Err(os_error(e, path, None)),
    }
}
pub(crate) fn temporary(parent: &Path, prefix: &str) -> Result<(File, PathBuf)> {
    let mut raw = CString::new(
        native(&tempfile_parent(parent)?)?
            .join(format!("{prefix}XXXXXX"))
            .as_os_str()
            .as_bytes(),
    )
    .map_err(|_| error("ValueError", "embedded null byte"))?
    .into_bytes_with_nul();
    let fd = unsafe { libc::mkostemp(raw.as_mut_ptr().cast(), libc::O_CLOEXEC) };
    let path =
        crate::python_unicode::from_native(Path::new(OsStr::from_bytes(&raw[..raw.len() - 1])));
    if fd < 0 {
        Err(os_error(std::io::Error::last_os_error(), &path, None))
    } else {
        Ok((unsafe { File::from_raw_fd(fd) }, path))
    }
}
pub(super) fn copystat(source: &Path, target: &Path) -> Result<()> {
    let m = metadata(source)?;
    let src = CString::new(native(source)?.as_os_str().as_bytes())
        .map_err(|_| error("ValueError", "embedded null byte"))?;
    let dst = CString::new(native(target)?.as_os_str().as_bytes())
        .map_err(|_| error("ValueError", "embedded null byte"))?;
    let times = [
        libc::timespec {
            tv_sec: m.atime(),
            tv_nsec: m.atime_nsec(),
        },
        libc::timespec {
            tv_sec: m.mtime(),
            tv_nsec: m.mtime_nsec(),
        },
    ];
    if unsafe { libc::utimensat(libc::AT_FDCWD, dst.as_ptr(), times.as_ptr(), 0) } != 0 {
        return Err(os_error(std::io::Error::last_os_error(), target, None));
    }
    // shutil.copystat copies Linux extended attributes before chmod. Follow
    // symlinks and preserve its narrowly ignored errno set, not all failures.
    let ignored = |e: &std::io::Error| {
        matches!(
            e.raw_os_error(),
            Some(libc::ENOTSUP | libc::ENODATA | libc::EINVAL)
        )
    };
    let len = unsafe { libc::listxattr(src.as_ptr(), std::ptr::null_mut(), 0) };
    if len < 0 {
        let e = std::io::Error::last_os_error();
        if !ignored(&e) {
            return Err(os_error(e, source, None));
        }
    } else if len > 0 {
        let mut names = vec![0u8; len as usize];
        let n = unsafe { libc::listxattr(src.as_ptr(), names.as_mut_ptr().cast(), names.len()) };
        if n < 0 {
            return Err(os_error(std::io::Error::last_os_error(), source, None));
        }
        for key in names[..n as usize]
            .split(|b| *b == 0)
            .filter(|v| !v.is_empty())
        {
            let key = CString::new(key).unwrap();
            let len =
                unsafe { libc::getxattr(src.as_ptr(), key.as_ptr(), std::ptr::null_mut(), 0) };
            let copied = (|| -> std::io::Result<()> {
                if len < 0 {
                    return Err(std::io::Error::last_os_error());
                }
                let mut data = vec![0u8; len as usize];
                let n = unsafe {
                    libc::getxattr(
                        src.as_ptr(),
                        key.as_ptr(),
                        data.as_mut_ptr().cast(),
                        data.len(),
                    )
                };
                if n < 0 {
                    return Err(std::io::Error::last_os_error());
                }
                if unsafe {
                    libc::setxattr(
                        dst.as_ptr(),
                        key.as_ptr(),
                        data.as_ptr().cast(),
                        n as usize,
                        0,
                    )
                } < 0
                {
                    return Err(std::io::Error::last_os_error());
                }
                Ok(())
            })();
            if let Err(e) = copied {
                if !ignored(&e) && e.raw_os_error() != Some(libc::EPERM) {
                    return Err(os_error(e, target, None));
                }
            }
        }
    }
    fs::set_permissions(
        native(target)?,
        fs::Permissions::from_mode(m.mode() & 0o7777),
    )
    .map_err(|e| os_error(e, target, None))
}
pub(crate) fn copy2(source: &Path, target: &Path) -> Result<()> {
    let m = metadata(source)?;
    if m.mode() & libc::S_IFMT == libc::S_IFIFO {
        return Err(error(
            "SpecialFileError",
            format!("`{}` is a named pipe", source.display()),
        ));
    }
    let mut input = File::open(native(source)?).map_err(|e| os_error(e, source, None))?;
    let mut output = File::create(native(target)?).map_err(|e| os_error(e, target, None))?;
    let mut block = vec![0; 8 * 1024 * 1024];
    loop {
        let n = input
            .read(&mut block)
            .map_err(|e| os_error(e, source, None))?;
        if n == 0 {
            break;
        }
        output
            .write_all(&block[..n])
            .map_err(|e| os_error(e, target, None))?;
    }
    drop(output);
    copystat(source, target)
}
pub(crate) fn stamp() -> Result<String> {
    let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)?;
    let seconds = now.as_secs() as libc::time_t;
    let mut tm = unsafe { std::mem::zeroed::<libc::tm>() };
    unsafe { libc::gmtime_r(&seconds, &mut tm) };
    Ok(format!(
        "{:04}{:02}{:02}T{:02}{:02}{:02}.{:06}Z",
        tm.tm_year + 1900,
        tm.tm_mon + 1,
        tm.tm_mday,
        tm.tm_hour,
        tm.tm_min,
        tm.tm_sec,
        now.subsec_micros()
    ))
}
