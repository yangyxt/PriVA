//! tempfile.TemporaryDirectory cleanup, including its permission retry.
use super::*;
use std::os::unix::fs::PermissionsExt;

pub(super) struct TemporaryDirectory(pub(super) PathBuf, bool);
impl TemporaryDirectory {
    /// Python's context exit runs cleanup even after an operation failed, and
    /// a cleanup exception replaces that operation's exception. Drop alone
    /// cannot express that behavior, so every normal exit calls this method.
    pub(super) fn finish<T>(mut self, outcome: Result<T>) -> Result<T> {
        self.1 = true;
        remove_tree(&self.0, &self.0, false)?;
        outcome
    }
    pub(super) fn new(parent: &Path, prefix: &str) -> Result<Self> {
        let parent = crate::python_path::tempfile_parent(parent)?;
        let mut template = std::ffi::CString::new(
            parent
                .join(format!("{prefix}XXXXXX"))
                .as_os_str()
                .as_encoded_bytes(),
        )
        .map_err(|_| error("ValueError", "embedded null byte"))?
        .into_bytes_with_nul();
        let result = unsafe { libc::mkdtemp(template.as_mut_ptr().cast()) };
        if result.is_null() {
            return Err(os_error(io::Error::last_os_error(), &parent, None));
        }
        template.pop();
        use std::os::unix::ffi::OsStringExt;
        Ok(Self(
            PathBuf::from(std::ffi::OsString::from_vec(template)),
            false,
        ))
    }
}

impl Drop for TemporaryDirectory {
    fn drop(&mut self) {
        // Only unwinding can reach this fallback. Ordinary success and error
        // paths explicitly preserve cleanup failures through finish().
        if !self.1 {
            let _ = remove_tree(&self.0, &self.0, false);
        }
    }
}

fn reset_permissions(path: &Path) -> Result<()> {
    let info = fs::symlink_metadata(path).map_err(|e| os_error(e, path, None))?;
    // Linux Python's _dont_follow_symlinks skips chmod on symlinks. The staged
    // selected publisher table is a symlink and must never be chmod'ed.
    if !info.file_type().is_symlink() {
        fs::set_permissions(path, fs::Permissions::from_mode(0o700))
            .map_err(|e| os_error(e, path, None))?;
    }
    Ok(())
}

fn cleanup_error(path: &Path, root: &Path, repeated: bool, cause: io::Error) -> Result<()> {
    if cause.kind() == io::ErrorKind::NotFound {
        return Ok(());
    }
    if cause.kind() != io::ErrorKind::PermissionDenied || (repeated && path == root) {
        return Err(os_error(cause, path, None));
    }
    // tempfile._rmtree's permission handler resets the affected parent and
    // entry, then retries unlink before deciding whether to recurse. It does
    // not reset permissions on every file in an otherwise healthy directory.
    let retried = (|| {
        if path != root {
            reset_permissions(path.parent().unwrap())?;
        }
        reset_permissions(path)?;
        match fs::remove_file(path) {
            Ok(()) => Ok(()),
            Err(e) if e.kind() == io::ErrorKind::IsADirectory => remove_tree(path, path, false),
            Err(e) if e.kind() == io::ErrorKind::PermissionDenied && path.is_dir() => {
                remove_tree(path, path, path == root)
            }
            Err(e) => Err(os_error(e, path, None)),
        }
    })();
    // A path disappearing during cleanup is ignored by tempfile's handler.
    match retried {
        Err(cause) if cause.to_string().starts_with("FileNotFoundError: ") => Ok(()),
        result => result,
    }
}

fn remove_tree(path: &Path, root: &Path, repeated: bool) -> Result<()> {
    let info = match fs::symlink_metadata(path) {
        Ok(info) => info,
        Err(cause) => return cleanup_error(path, root, repeated, cause),
    };
    if info.file_type().is_symlink() {
        return Err(error("OSError", "Cannot call rmtree on a symbolic link"));
    }
    let entries = match fs::read_dir(path) {
        Ok(entries) => entries,
        Err(cause) => return cleanup_error(path, root, repeated, cause),
    };
    for entry in entries {
        let entry = entry.map_err(|e| os_error(e, path, None))?;
        let child = entry.path();
        let kind = match entry.file_type() {
            Ok(kind) => kind,
            Err(cause) => {
                cleanup_error(&child, root, repeated, cause)?;
                continue;
            }
        };
        if kind.is_dir() {
            remove_tree(&child, root, repeated)?;
        } else if let Err(cause) = fs::remove_file(&child) {
            cleanup_error(&child, root, repeated, cause)?;
        }
    }
    if let Err(cause) = fs::remove_dir(path) {
        cleanup_error(path, root, repeated, cause)?;
    }
    Ok(())
}
