use super::*;
use std::collections::HashMap;
use std::ffi::{CString, OsStr};
use std::fs::File;
use std::io::Write;
use std::os::fd::FromRawFd;
use std::os::unix::ffi::OsStrExt;
#[derive(Debug)]
pub(super) struct PathLoop(String);
impl std::fmt::Display for PathLoop {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for PathLoop {}
pub(super) fn loop_error(path: &Path) -> Box<dyn std::error::Error> {
    Box::new(PathLoop(format!(
        "Symlink loop from {}",
        repr(&path.display().to_string())
    )))
}
pub(super) fn read_text(path: &Path) -> Result<String> {
    let bytes = fs::read(path).map_err(|e| os_error(e, path, None))?;
    text(&bytes)
}
pub(super) fn exists(path: &Path) -> Result<bool> {
    kind(path, false).map(|v| v.is_some())
}
pub(super) fn kind(path: &Path, link: bool) -> Result<Option<fs::Metadata>> {
    match if link {
        fs::symlink_metadata(path)
    } else {
        fs::metadata(path)
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
pub(super) fn is_link(path: &Path) -> Result<bool> {
    Ok(kind(path, true)?.is_some_and(|m| m.file_type().is_symlink()))
}
pub(super) fn is_dir(path: &Path) -> Result<bool> {
    Ok(kind(path, false)?.is_some_and(|m| m.is_dir()))
}

// Non-strict posixpath.realpath for CADD's Snakemake identity. Unlike a lexical
// abspath, '..' after a symbolic link is evaluated at the resolved target.
pub(super) fn realpath(path: &Path, pathlib: bool) -> Result<PathBuf> {
    fn walk(
        mut path: PathBuf,
        rest: &Path,
        seen: &mut HashMap<PathBuf, Option<PathBuf>>,
    ) -> Result<PathBuf> {
        if rest.is_absolute() {
            path = PathBuf::from("/");
        }
        let spelling = rest.as_os_str().as_bytes();
        for part in spelling.split(|b| *b == b'/') {
            if part.is_empty() || part == b"." {
                continue;
            }
            if part == b".." {
                path.pop();
                continue;
            }
            let next = path.join(OsStr::from_bytes(part));
            let linked = fs::symlink_metadata(&next).is_ok_and(|m| m.file_type().is_symlink());
            if !linked {
                path = next;
                continue;
            }
            if let Some(found) = seen.get(&next) {
                path = found.clone().unwrap_or(next);
                continue;
            }
            seen.insert(next.clone(), None);
            let target = fs::read_link(&next).map_err(|e| os_error(e, &next, None))?;
            path = walk(path, &target, seen)?;
            seen.insert(next, Some(path.clone()));
        }
        Ok(path)
    }
    let out = walk(std::env::current_dir()?, path, &mut HashMap::new())?;
    if pathlib && fs::metadata(&out).is_err_and(|e| e.raw_os_error() == Some(libc::ELOOP)) {
        return Err(loop_error(&out));
    }
    Ok(out)
}
pub(super) fn publish(path: &Path, payload: &str) -> Result<()> {
    let parent = path.parent().unwrap();
    let name = path.file_name().unwrap().to_string_lossy();
    let template = crate::python_path::tempfile_parent(parent)?.join(format!(".{name}.XXXXXX"));
    let mut template = CString::new(template.as_os_str().as_bytes())
        .map_err(|_| error("ValueError", "embedded null byte"))?
        .into_bytes_with_nul();
    let fd = unsafe { libc::mkostemp(template.as_mut_ptr().cast(), libc::O_CLOEXEC) };
    let temporary = PathBuf::from(OsStr::from_bytes(&template[..template.len() - 1]));
    if fd < 0 {
        return Err(os_error(std::io::Error::last_os_error(), &temporary, None));
    }
    let result = (|| {
        let mut file = unsafe { File::from_raw_fd(fd) };
        // NamedTemporaryFile publishes its name before write/flush/fsync. Close
        // completes before rename, and finally unlink can supersede an error.
        file.write_all(payload.as_bytes())
            .map_err(|e| os_error(e, &temporary, None))?;
        file.sync_all().map_err(|e| os_error(e, &temporary, None))?;
        drop(file);
        fs::rename(&temporary, path).map_err(|e| os_error(e, &temporary, Some(path)))?;
        Ok(())
    })();
    match fs::remove_file(&temporary) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
        Err(e) => return Err(os_error(e, &temporary, None)),
    }
    result
}
