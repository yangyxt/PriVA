//! `priva loftee-resources validate --ancestor FASTA_GZ --conservation SQLITE
//! --assembly ASSEMBLY`: the checks of LOFTEE's ancestral FASTA and
//! conservation database that `validate_loftee_resources` ran as an inline
//! Python block before this command. The shell function still checks that the
//! files exist and runs `gzip -t`; this command checks that:
//!
//! - the ancestor's first line is a FASTA header;
//! - every line of its `.fai` index has five fields, the second to fifth
//!   decimal numbers and a length of at least 1, and no chromosome (with or
//!   without `chr`) has two lengths;
//! - chromosome 1 has the length of the configured assembly (249,250,621 for
//!   hg19/GRCh37, 248,956,422 otherwise);
//! - the conservation file, opened read-only through SQLite, passes
//!   `PRAGMA quick_check` and has at least one table.
//!
//! A failure prints `ERROR: Invalid LoFtee resource: REASON` and exits with
//! 1, as the Python block did. SQLite is the copy vendored for the crate
//! (`third_party/sqlite3`). Compared with the Python block by
//! `rust/tests/install_checks_compare.py`.
use crate::cli::{Result, invalid};
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::ffi::{CStr, CString, c_char, c_int, c_void};
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::os::unix::ffi::OsStrExt;
use std::path::Path;
use std::ptr;

const USAGE: &str = "usage: priva loftee-resources validate --ancestor FASTA_GZ --conservation SQLITE --assembly ASSEMBLY\n\
Check the LOFTEE ancestral FASTA, its .fai index and the conservation SQLite database.";

type Handle = *mut c_void;
const SQLITE_OPEN_READONLY: c_int = 0x1;
const SQLITE_ROW: c_int = 100;
const SQLITE_DONE: c_int = 101;

unsafe extern "C" {
    fn sqlite3_open_v2(
        path: *const c_char,
        db: *mut Handle,
        flags: c_int,
        vfs: *const c_char,
    ) -> c_int;
    fn sqlite3_close_v2(db: Handle) -> c_int;
    fn sqlite3_errmsg(db: Handle) -> *const c_char;
    fn sqlite3_prepare_v2(
        db: Handle,
        sql: *const c_char,
        bytes: c_int,
        stmt: *mut Handle,
        tail: *mut *const c_char,
    ) -> c_int;
    fn sqlite3_step(stmt: Handle) -> c_int;
    fn sqlite3_finalize(stmt: Handle) -> c_int;
    fn sqlite3_column_count(stmt: Handle) -> c_int;
    fn sqlite3_column_text(stmt: Handle, index: c_int) -> *const u8;
}

/// A read-only connection, closed when dropped.
struct Database(Handle);

impl Drop for Database {
    fn drop(&mut self) {
        // SAFETY: the handle came from sqlite3_open_v2 and is closed once.
        unsafe { sqlite3_close_v2(self.0) };
    }
}

impl Database {
    fn open(path: &Path) -> std::result::Result<Self, String> {
        let name = CString::new(path.as_os_str().as_bytes())
            .map_err(|_| "embedded null byte".to_owned())?;
        let mut handle = ptr::null_mut();
        // SAFETY: valid NUL-terminated path and out pointer.
        let code = unsafe {
            sqlite3_open_v2(
                name.as_ptr(),
                &mut handle,
                SQLITE_OPEN_READONLY,
                ptr::null(),
            )
        };
        if handle.is_null() {
            return Err("out of memory".to_owned());
        }
        let database = Self(handle);
        if code != 0 {
            return Err(database.message());
        }
        Ok(database)
    }

    fn message(&self) -> String {
        // SAFETY: SQLite returns a NUL-terminated message owned by the handle.
        unsafe { CStr::from_ptr(sqlite3_errmsg(self.0)) }
            .to_string_lossy()
            .into_owned()
    }

    /// Every row of a query, each column as text (`None` for NULL).
    fn rows(&self, sql: &str) -> std::result::Result<Vec<Vec<Option<String>>>, String> {
        let text = CString::new(sql).map_err(|_| "embedded null byte".to_owned())?;
        let mut statement = ptr::null_mut();
        // SAFETY: live handle, NUL-terminated SQL, out pointer.
        let code = unsafe {
            sqlite3_prepare_v2(self.0, text.as_ptr(), -1, &mut statement, ptr::null_mut())
        };
        if code != 0 {
            return Err(self.message());
        }
        let mut rows = Vec::new();
        let outcome = loop {
            // SAFETY: the statement was prepared above and is finalized below.
            match unsafe { sqlite3_step(statement) } {
                SQLITE_ROW => {
                    let columns = unsafe { sqlite3_column_count(statement) };
                    let row = (0..columns)
                        .map(|column| {
                            let value = unsafe { sqlite3_column_text(statement, column) };
                            (!value.is_null()).then(|| {
                                unsafe { CStr::from_ptr(value.cast()) }
                                    .to_string_lossy()
                                    .into_owned()
                            })
                        })
                        .collect();
                    rows.push(row);
                }
                SQLITE_DONE => break Ok(rows),
                _ => break Err(self.message()),
            }
        };
        // SAFETY: finalized once.
        unsafe { sqlite3_finalize(statement) };
        outcome
    }
}

fn ascii_number(text: &str) -> Option<u64> {
    (!text.is_empty() && text.bytes().all(|b| b.is_ascii_digit()))
        .then(|| text.parse().ok())
        .flatten()
}

fn check_ancestor(ancestor: &str) -> std::result::Result<HashMap<String, u64>, String> {
    let file = File::open(ancestor).map_err(|error| error.to_string())?;
    let mut first = Vec::new();
    BufReader::new(MultiGzDecoder::new(BufReader::new(file)))
        .read_until(b'\n', &mut first)
        .map_err(|error| error.to_string())?;
    if !first.starts_with(b">") {
        return Err("LoFtee ancestor is not FASTA".to_owned());
    }
    let index = fs::read(format!("{ancestor}.fai")).map_err(|error| error.to_string())?;
    let index = String::from_utf8(index).map_err(|error| error.to_string())?;
    let mut lengths = HashMap::new();
    // Python iterates the lines of a text file: `\n`, `\r\n` or `\r` ends one.
    let mut text = index.as_str();
    while !text.is_empty() {
        let end = text.find(['\n', '\r']).unwrap_or(text.len());
        let line = &text[..end];
        let mut next = end;
        if text[end..].starts_with("\r\n") {
            next += 2;
        } else if end < text.len() {
            next += 1;
        }
        text = &text[next..];
        let fields: Vec<&str> = line.split('\t').collect();
        let numbers: Vec<Option<u64>> = fields
            .iter()
            .skip(1)
            .take(4)
            .map(|f| ascii_number(f))
            .collect();
        if fields.len() < 5 || numbers.iter().any(Option::is_none) || numbers[0] < Some(1) {
            return Err("Malformed LoFtee ancestor FASTA index".to_owned());
        }
        let length = numbers[0].unwrap_or_default();
        let name = fields[0]
            .strip_prefix("chr")
            .unwrap_or(fields[0])
            .to_owned();
        if lengths.get(&name).is_some_and(|known| *known != length) {
            return Err("Conflicting LoFtee ancestor chromosome lengths".to_owned());
        }
        lengths.insert(name, length);
    }
    Ok(lengths)
}

fn check_conservation(conservation: &str) -> std::result::Result<(), String> {
    // Python opens Path(conservation).resolve() as a read-only URI.
    let path =
        fs::canonicalize(conservation).unwrap_or_else(|_| Path::new(conservation).to_path_buf());
    let database = Database::open(&path)?;
    if database.rows("PRAGMA quick_check")? != vec![vec![Some("ok".to_owned())]] {
        return Err("LoFtee conservation SQLite integrity check failed".to_owned());
    }
    if database
        .rows("SELECT name FROM sqlite_master WHERE type='table' LIMIT 1")?
        .is_empty()
    {
        return Err("LoFtee conservation SQLite has no tables".to_owned());
    }
    Ok(())
}

pub(crate) fn validate(
    ancestor: &str,
    conservation: &str,
    assembly: &str,
) -> std::result::Result<(), String> {
    let lengths = check_ancestor(ancestor)?;
    let expected = if matches!(assembly, "hg19" | "GRCh37") {
        249_250_621
    } else {
        248_956_422
    };
    if lengths.get("1") != Some(&expected) {
        return Err("LoFtee ancestor FASTA index has the wrong assembly".to_owned());
    }
    check_conservation(conservation)
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args
        .first()
        .is_some_and(|arg| arg == "--help" || arg == "-h")
    {
        println!("{USAGE}");
        return Ok(());
    }
    if args.first().map(String::as_str) != Some("validate") {
        return Err(invalid(USAGE));
    }
    let (mut ancestor, mut conservation, mut assembly) = (None, None, None);
    let mut rest = args[1..].iter();
    while let Some(arg) = rest.next() {
        let slot = match arg.as_str() {
            "--ancestor" => &mut ancestor,
            "--conservation" => &mut conservation,
            "--assembly" => &mut assembly,
            "--help" | "-h" => {
                println!("{USAGE}");
                return Ok(());
            }
            _ => return Err(invalid(format!("unrecognized argument {arg}\n{USAGE}"))),
        };
        *slot = Some(
            rest.next()
                .ok_or_else(|| invalid(format!("{arg} needs a value")))?
                .clone(),
        );
    }
    let (Some(ancestor), Some(conservation), Some(assembly)) = (ancestor, conservation, assembly)
    else {
        return Err(invalid(USAGE));
    };
    if let Err(message) = validate(&ancestor, &conservation, &assembly) {
        eprintln!("ERROR: Invalid LoFtee resource: {message}");
        std::process::exit(1);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn ancestor_index_and_database() {
        let folder = std::env::temp_dir().join(format!("priva-loftee-{}", std::process::id()));
        fs::create_dir_all(&folder).unwrap();
        let ancestor = folder.join("ancestor.fa.gz");
        let mut gz = flate2::write::GzEncoder::new(
            File::create(&ancestor).unwrap(),
            flate2::Compression::default(),
        );
        gz.write_all(b">1\nACGT\n").unwrap();
        gz.finish().unwrap();
        let ancestor = ancestor.to_string_lossy().into_owned();
        fs::write(
            format!("{ancestor}.fai"),
            "1\t248956422\t3\t60\t61\nchr2\t5\t9\t60\t61\n",
        )
        .unwrap();
        let missing = folder.join("absent.sqlite").to_string_lossy().into_owned();
        assert_eq!(
            validate(&ancestor, &missing, "hg38"),
            Err("unable to open database file".to_owned())
        );
        assert_eq!(
            validate(&ancestor, &missing, "hg19"),
            Err("LoFtee ancestor FASTA index has the wrong assembly".to_owned())
        );
        let text = folder.join("text.sqlite");
        fs::write(&text, b"not a database at all, just some bytes of text").unwrap();
        assert_eq!(
            validate(&ancestor, &text.to_string_lossy(), "hg38"),
            Err("file is not a database".to_owned())
        );
        fs::write(
            format!("{ancestor}.fai"),
            "1\t248956422\t3\t60\t61\nchr1\t5\t9\t60\t61\n",
        )
        .unwrap();
        assert_eq!(
            validate(&ancestor, &missing, "hg38"),
            Err("Conflicting LoFtee ancestor chromosome lengths".to_owned())
        );
        fs::write(format!("{ancestor}.fai"), "1\t248956422\t3\t60\n").unwrap();
        assert_eq!(
            validate(&ancestor, &missing, "hg38"),
            Err("Malformed LoFtee ancestor FASTA index".to_owned())
        );
    }
}
