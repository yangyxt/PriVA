//! Pathlib and os.path spellings, Python 3.12 file tests, text reading and the
//! optional inspection-root policy (`checked_path`/`readable_file`).
use super::{Result, Value, error, require, value_error};
use crate::python_io::GzipFailure;
use crate::python_lines::Lines;
use crate::python_path::{metadata, os_error, path_value, tempfile_parent};
use crate::python_text::is_python_space;
use std::error::Error;
use std::ffi::CString;
use std::fs::{self, File};
use std::io::{self, Read};
use std::path::Path;

/// `str(Path(value))`
pub(super) fn pathlib(value: &str) -> String {
    path_value(value).display().to_string()
}

/// `str(Path(base) / relative)` for a pathlib-normalized `base`.
pub(super) fn join(base: &str, relative: &str) -> String {
    if relative.starts_with('/') {
        pathlib(relative)
    } else if base.ends_with('/') {
        pathlib(&format!("{base}{relative}"))
    } else {
        pathlib(&format!("{base}/{relative}"))
    }
}

/// `os.path.abspath(value)`
pub(super) fn abspath(value: &str) -> Result<String> {
    Ok(tempfile_parent(Path::new(value))?.display().to_string())
}

/// `Path(path).resolve()`: symbolic links followed where they exist.
pub(super) fn resolve(path: &str) -> Result<String> {
    Ok(
        crate::clinvar_contract::resolve_non_strict(Path::new(path))?
            .display()
            .to_string(),
    )
}

/// `Path(path).name`
pub(super) fn name(path: &str) -> &str {
    path.rsplit('/').next().unwrap_or("")
}

/// `Path(path).suffix`
pub(super) fn suffix(path: &str) -> &str {
    let name = name(path);
    match name.rfind('.') {
        Some(index) if index > 0 && index < name.len() - 1 => &name[index..],
        _ => "",
    }
}

/// `Path(path).parts` of a normalized path.
pub(super) fn parts(path: &str) -> Vec<&str> {
    let (anchor, rest) = if path.starts_with("//") && !path.starts_with("///") {
        ("//", &path[2..])
    } else if path.starts_with('/') {
        ("/", path.trim_start_matches('/'))
    } else {
        ("", path)
    };
    let mut parts = Vec::new();
    if !anchor.is_empty() {
        parts.push(anchor);
    }
    for part in rest.split('/') {
        if !part.is_empty() && part != "." {
            parts.push(part);
        }
    }
    parts
}

fn from_parts(parts: &[&str]) -> String {
    match parts.first() {
        Some(&anchor) if anchor.starts_with('/') => format!("{anchor}{}", parts[1..].join("/")),
        Some(_) => parts.join("/"),
        None => ".".to_owned(),
    }
}

/// `Path(path).parent`
fn parent(path: &str) -> String {
    let all = parts(path);
    if all.len() <= 1 {
        return path.to_owned();
    }
    from_parts(&all[..all.len() - 1])
}

/// `Path(path) == Path(base) or Path(base) in Path(path).parents`
pub(super) fn is_relative_to(path: &str, base: &str) -> bool {
    let (path, base) = (parts(path), parts(base));
    path.len() >= base.len() && path[..base.len()] == base[..]
}

/// The errors pathlib's file tests treat as "no such file".
fn ignored(failure: &io::Error) -> bool {
    matches!(
        failure.raw_os_error(),
        Some(libc::ENOENT | libc::ENOTDIR | libc::EBADF | libc::ELOOP)
    )
}

fn tested(
    path: &str,
    result: io::Result<fs::Metadata>,
    test: fn(&fs::Metadata) -> bool,
) -> Result<bool> {
    match result {
        Ok(metadata) => Ok(test(&metadata)),
        Err(failure) if ignored(&failure) => Ok(false),
        Err(failure) => Err(os_error(failure, Path::new(path), None)),
    }
}

/// `Path(path).exists()`
pub(super) fn exists(path: &str) -> Result<bool> {
    tested(path, fs::metadata(path), |_| true)
}

/// `Path(path).is_file()`
pub(super) fn is_file(path: &str) -> Result<bool> {
    tested(path, fs::metadata(path), fs::Metadata::is_file)
}

/// `Path(path).is_dir()`
pub(super) fn is_dir(path: &str) -> Result<bool> {
    tested(path, fs::metadata(path), fs::Metadata::is_dir)
}

/// `Path(path).is_symlink()`
pub(super) fn is_symlink(path: &str) -> Result<bool> {
    tested(path, fs::symlink_metadata(path), |metadata| {
        metadata.file_type().is_symlink()
    })
}

/// `os.path.lexists(path)`: every failure reads as absent.
pub(super) fn lexists(path: &str) -> bool {
    fs::symlink_metadata(path).is_ok()
}

/// `os.access(path, os.R_OK)`
fn accessible(path: &str) -> bool {
    CString::new(path.as_bytes())
        .is_ok_and(|text| unsafe { libc::access(text.as_ptr(), libc::R_OK) } == 0)
}

/// `path.is_file() and os.access(path, os.R_OK) and path.stat().st_size`
pub(super) fn nonempty_readable(path: &str) -> Result<bool> {
    if !is_file(path)? || !accessible(path) {
        return Ok(false);
    }
    Ok(metadata(Path::new(path))?.len() != 0)
}

/// `checked_path(value, key)` for a configured value.
pub(super) fn checked(value: &Value, key: &str) -> Result<String> {
    match value {
        Value::Str(text) => checked_text(text, key, 0),
        _ => Err(unset(key)),
    }
}

fn unset(key: &str) -> Box<dyn Error> {
    value_error(format!(
        "Set {key} to the original, explicitly supplied input path"
    ))
}

/// `checked_path(value, key, _depth)`: the absolute path, or, when
/// PRIVA_INSTALL_INSPECT_ROOTS names roots, the path inside them with each link
/// followed only to a target that stays inside them.
pub(super) fn checked_text(value: &str, key: &str, depth: usize) -> Result<String> {
    if value.chars().all(is_python_space) || value == "null" {
        return Err(unset(key));
    }
    require(
        !value.starts_with("/paedyl01/"),
        format!("Unresolved old-host path {key}: {value}; supply its actual local path"),
    )?;
    let path = abspath(value)?;
    let setting = std::env::var_os("PRIVA_INSTALL_INSPECT_ROOTS").unwrap_or_default();
    let roots = setting
        .to_string_lossy()
        .split(':')
        .filter(|root| !root.is_empty())
        .map(abspath)
        .collect::<Result<Vec<_>>>()?;
    if roots.is_empty() {
        return Ok(path);
    }
    if depth > 40 {
        return Err(value_error(format!("Too many symlinks for {key}: {value}")));
    }
    // max() keeps the first of equally deep roots.
    let mut selected: Option<&String> = None;
    for root in roots.iter().filter(|root| is_relative_to(&path, root)) {
        if selected.is_none_or(|best| parts(root).len() > parts(best).len()) {
            selected = Some(root);
        }
    }
    let root = selected.ok_or_else(|| {
        value_error(format!(
            "Path outside permitted inspection roots for {key}: {value}"
        ))
    })?;
    let mut cursor = root.clone();
    for part in &parts(&path)[parts(root).len()..] {
        cursor = join(&cursor, part);
        if is_symlink(&cursor)? {
            // Do not follow a link into data outside the selected policy.
            let link = fs::read_link(&cursor)
                .map_err(|failure| os_error(failure, Path::new(&cursor), None))?;
            let target = abspath(&join(&parent(&cursor), &link.to_string_lossy()))?;
            require(
                roots.iter().any(|allowed| is_relative_to(&target, allowed)),
                format!("Symlink leaves permitted inspection roots for {key}: {value}"),
            )?;
            cursor = checked_text(&target, key, depth + 1)?;
        }
    }
    Ok(cursor)
}

/// `readable_file(value, key)` for a configured value.
pub(super) fn readable(value: &Value, key: &str) -> Result<String> {
    let path = checked(value, key)?;
    let shown = value.text().unwrap_or_default();
    require(
        nonempty_readable(&path)?,
        format!("Missing, empty, or unreadable {key}: {shown}"),
    )?;
    Ok(path)
}

/// `readable_file(value, key)` for a path text.
pub(super) fn readable_text(value: &str, key: &str) -> Result<String> {
    let path = checked_text(value, key, 0)?;
    require(
        nonempty_readable(&path)?,
        format!("Missing, empty, or unreadable {key}: {value}"),
    )?;
    Ok(path)
}

/// A failure while reading text, as the Python reader reports it.
pub(super) fn read_failure(failure: Box<dyn Error>, path: &str) -> Box<dyn Error> {
    if let Some(cause) = failure.downcast_ref::<io::Error>() {
        if let Some(gzip) = GzipFailure::of(cause) {
            let category = match gzip {
                GzipFailure::EndOfFile => "EOFError",
                GzipFailure::Bad(_) => "BadGzipFile",
                GzipFailure::Zlib(_) => "error",
            };
            return error(category, gzip.to_string());
        }
        if let Some(code) = cause.raw_os_error() {
            return os_error(io::Error::from_raw_os_error(code), Path::new(path), None);
        }
    }
    failure
}

/// The lines of a UTF-8 text stream, read as Python's text mode reads them:
/// 8 KiB decode blocks and universal newlines translated to `\n`.
pub(super) struct TextLines<R: Read> {
    lines: Lines<R>,
    path: String,
}

impl<R: Read> TextLines<R> {
    pub(super) fn new(source: R, path: &str) -> Self {
        Self {
            lines: Lines::new(source),
            path: path.to_owned(),
        }
    }
    /// The next line with its `\n`, or `None` at the end.
    pub(super) fn readline(&mut self) -> Result<Option<String>> {
        let line = self
            .lines
            .readline(usize::MAX)
            .map_err(|failure| read_failure(failure, &self.path))?;
        Ok((!line.is_empty()).then_some(line))
    }
}

/// `open(path)`
pub(super) fn open(path: &str) -> Result<File> {
    File::open(path).map_err(|failure| os_error(failure, Path::new(path), None))
}

/// `open(path, encoding="utf-8").read()`
pub(super) fn read_text(path: &str) -> Result<String> {
    let mut lines = TextLines::new(open(path)?, path);
    let mut text = String::new();
    while let Some(line) = lines.readline()? {
        text.push_str(&line);
    }
    Ok(text)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pathlib_spellings() {
        assert_eq!(join("/a/b", "c/./d"), "/a/b/c/d");
        assert_eq!(join("/", "x"), "/x");
        assert_eq!(join("//", "x"), "//x");
        assert_eq!(join(".", "x"), "x");
        assert_eq!(join("/a", "/abs//p"), "/abs/p");
        assert_eq!(join("/a", "../b"), "/a/../b");
        assert_eq!(abspath("/a/../b/./c").unwrap(), "/b/c");
        assert_eq!(name("/a/b.tar.gz"), "b.tar.gz");
        assert_eq!(suffix("/a/b.tar.gz"), ".gz");
        assert_eq!(suffix("/a/.hidden"), "");
        assert_eq!(suffix("/a/trailing."), "");
        assert_eq!(parts("/a/b"), ["/", "a", "b"]);
        assert_eq!(parts("//a"), ["//", "a"]);
        assert_eq!(parent("/a/b"), "/a");
        assert_eq!(parent("/a"), "/");
        assert!(is_relative_to("/a/b", "/a"));
        assert!(is_relative_to("/a", "/a"));
        assert!(!is_relative_to("/ab", "/a"));
        assert!(!is_relative_to("/a", "/a/b"));
    }

    #[test]
    fn unset_and_old_host_values_fail_before_any_file_test() {
        for value in ["", "  ", "\u{1c}", "null"] {
            let failure = checked_text(value, "key", 0).err().unwrap();
            assert_eq!(
                super::super::exception(failure.as_ref()).unwrap().1,
                "Set key to the original, explicitly supplied input path"
            );
        }
        let failure = checked(&Value::Int("1".into()), "key").err().unwrap();
        assert!(failure.to_string().contains("Set key"));
        let failure = checked_text("/paedyl01/x", "key", 0).err().unwrap();
        assert!(failure.to_string().contains("Unresolved old-host path key"));
    }
}
