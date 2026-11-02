//! Shared POSIX Pathlib spellings and Python-style filesystem failures.
//! Extracted without changing the already qualified Conservation helper bodies.
use crate::cli::{Result, invalid};
use std::path::{Path, PathBuf};
use std::{fs, io};

fn error(category: &'static str, message: impl Into<String>) -> Box<dyn std::error::Error> {
    crate::builder_error::category(category, invalid(message))
}

/// Pathlib collapses `.` and redundant separators but keeps `..` and two
/// leading slashes. In particular absolute() must not resolve a symlink.
pub(crate) fn path_value(value: &str) -> PathBuf {
    let prefix = if value.starts_with("//") && !value.starts_with("///") {
        "//"
    } else if value.starts_with('/') {
        "/"
    } else {
        ""
    };
    let parts = value
        .split('/')
        .filter(|part| !part.is_empty() && *part != ".")
        .collect::<Vec<_>>();
    let text = format!("{prefix}{}", parts.join("/"));
    PathBuf::from(if text.is_empty() {
        ".".to_owned()
    } else {
        text
    })
}
pub(crate) fn absolute(path: &Path) -> Result<String> {
    if path.is_absolute() {
        return Ok(path.display().to_string());
    }
    Ok(path_value(&format!(
        "{}/{}",
        std::env::current_dir()?.display(),
        path.display()
    ))
    .display()
    .to_string())
}
pub(crate) fn repr(text: &str) -> String {
    let quote = if text.contains('\'') && !text.contains('"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::from(quote);
    for c in text.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c == quote => {
                out.push('\\');
                out.push(c);
            }
            c if c.is_control() && (c as u32) <= 255 => {
                out.push_str(&format!("\\x{:02x}", c as u32))
            }
            c if c.is_control() || matches!(c, '\u{2028}' | '\u{2029}') => {
                out.push_str(&format!("\\u{:04x}", c as u32))
            }
            c => out.push(c),
        }
    }
    out.push(quote);
    out
}
pub(crate) fn os_error(
    cause: io::Error,
    path: &Path,
    destination: Option<&Path>,
) -> Box<dyn std::error::Error> {
    let category = match cause.kind() {
        io::ErrorKind::NotFound => "FileNotFoundError",
        io::ErrorKind::PermissionDenied => "PermissionError",
        io::ErrorKind::AlreadyExists => "FileExistsError",
        io::ErrorKind::IsADirectory => "IsADirectoryError",
        io::ErrorKind::NotADirectory => "NotADirectoryError",
        _ => "OSError",
    };
    let message = if let Some(errno) = cause.raw_os_error() {
        let description = cause.to_string();
        let description = description
            .strip_suffix(&format!(" (os error {errno})"))
            .unwrap_or(&description);
        format!(
            "[Errno {errno}] {description}: {}{}",
            repr(&path.display().to_string()),
            destination.map_or(String::new(), |p| format!(
                " -> {}",
                repr(&p.display().to_string())
            ))
        )
    } else {
        cause.to_string()
    };
    error(category, message)
}
pub(crate) fn is_file(path: &Path) -> Result<bool> {
    match fs::metadata(path) {
        Ok(m) => Ok(m.is_file()),
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
pub(crate) fn metadata(path: &Path) -> Result<fs::Metadata> {
    fs::metadata(path).map_err(|e| os_error(e, path, None))
}
