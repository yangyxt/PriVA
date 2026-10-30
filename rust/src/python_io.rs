//! File and text behaviour of Python's standard library that the SpliceAI
//! scripts show in their error messages:
//! - `decode_utf8`: the text of a `UnicodeDecodeError` for a UTF-8 decode;
//! - `bytes_repr`: `repr()` of a `bytes` value, as in `Not a gzipped file (b'##')`;
//! - `resolve_strict`: `Path.resolve(strict=True)`, which names the first path
//!   component that cannot be examined (`/a/b` when `/a` exists and `/a/b`
//!   does not, even for `/a/b/c`), not the whole path;
//! - `PythonGzipReader`: reads gzip files the way `gzip.GzipFile` does, member
//!   by member, with its checks and its error texts. A damaged file fails with
//!   `EOFError` (a missing end), `BadGzipFile` (wrong magic bytes, unknown
//!   method, failed CRC or length check) or `zlib.error`. The last of these
//!   carries zlib's own reason in Python ("invalid distance too far back"); the
//!   deflate decoder used here does not report one, so that text is generic.
use flate2::{Crc, Decompress, FlushDecompress, Status};
use std::collections::HashMap;
use std::fmt;
use std::fs;
use std::io::{self, Read};
use std::path::{Path, PathBuf};

// ---------------------------------------------------------------- text

/// `bytes.decode("utf-8")`, or the message of the `UnicodeDecodeError` Python raises.
pub(crate) fn decode_utf8(bytes: &[u8]) -> Result<&str, String> {
    match std::str::from_utf8(bytes) {
        Ok(text) => Ok(text),
        Err(error) => {
            let start = error.valid_up_to();
            let (end, reason) = match error.error_len() {
                None => (bytes.len(), "unexpected end of data"),
                Some(length) => {
                    let lead = bytes[start];
                    let reason = if matches!(lead, 0x80..=0xc1 | 0xf5..=0xff) {
                        "invalid start byte"
                    } else {
                        "invalid continuation byte"
                    };
                    (start + length, reason)
                }
            };
            Err(if end == start + 1 {
                format!(
                    "'utf-8' codec can't decode byte 0x{:02x} in position {start}: {reason}",
                    bytes[start]
                )
            } else {
                format!(
                    "'utf-8' codec can't decode bytes in position {start}-{}: {reason}",
                    end - 1
                )
            })
        }
    }
}

/// `repr()` of a `bytes` value.
pub(crate) fn bytes_repr(bytes: &[u8]) -> String {
    let quote = if bytes.contains(&b'\'') && !bytes.contains(&b'"') {
        b'"'
    } else {
        b'\''
    };
    let mut out = String::from("b");
    out.push(quote as char);
    for &byte in bytes {
        match byte {
            b'\\' => out.push_str("\\\\"),
            b'\t' => out.push_str("\\t"),
            b'\n' => out.push_str("\\n"),
            b'\r' => out.push_str("\\r"),
            byte if byte == quote => {
                out.push('\\');
                out.push(byte as char);
            }
            0x20..=0x7e => out.push(byte as char),
            byte => out.push_str(&format!("\\x{byte:02x}")),
        }
    }
    out.push(quote as char);
    out
}

// ---------------------------------------------------------------- paths

/// `posixpath.split`: everything up to the last slash (without trailing
/// slashes unless it is all slashes) and the rest.
fn split_path(path: &str) -> (String, String) {
    let cut = path.rfind('/').map_or(0, |at| at + 1);
    let (head, tail) = path.split_at(cut);
    let head = if !head.is_empty() && head.bytes().any(|byte| byte != b'/') {
        head.trim_end_matches('/')
    } else {
        head
    };
    (head.to_owned(), tail.to_owned())
}

/// `posixpath.join` of a path and one name.
fn join_name(path: &str, name: &str) -> String {
    if path.is_empty() {
        name.to_owned()
    } else if path.ends_with('/') {
        format!("{path}{name}")
    } else {
        format!("{path}/{name}")
    }
}

type Walk = Result<String, (io::Error, String)>;

/// `posixpath._joinrealpath` in strict mode: the resolved path, or the error
/// of the first component that cannot be examined with the path it names.
fn join_real_path(start: String, rest: &str, seen: &mut HashMap<String, Option<String>>) -> Walk {
    let mut path = start;
    let mut remaining = rest;
    if let Some(stripped) = remaining.strip_prefix('/') {
        remaining = stripped;
        path = "/".to_owned();
    }
    while !remaining.is_empty() {
        let (name, tail) = remaining.split_once('/').unwrap_or((remaining, ""));
        remaining = tail;
        if name.is_empty() || name == "." {
            continue;
        }
        if name == ".." {
            if path.is_empty() {
                path = "..".to_owned();
            } else {
                let (head, last) = split_path(&path);
                path = if last == ".." {
                    format!("{}/../..", path.trim_end_matches('/'))
                } else {
                    head
                };
            }
            continue;
        }
        let new_path = join_name(&path, name);
        let metadata =
            fs::symlink_metadata(&new_path).map_err(|error| (error, new_path.clone()))?;
        if !metadata.file_type().is_symlink() {
            path = new_path;
            continue;
        }
        match seen.get(&new_path) {
            Some(Some(resolved)) => {
                path = resolved.clone();
                continue;
            }
            Some(None) => {
                // Too many levels of symbolic links.
                return Err((io::Error::from_raw_os_error(40), new_path));
            }
            None => {}
        }
        seen.insert(new_path.clone(), None);
        let target = fs::read_link(&new_path).map_err(|error| (error, new_path.clone()))?;
        path = join_real_path(path, &target.to_string_lossy(), seen)?;
        seen.insert(new_path, Some(path.clone()));
    }
    Ok(path)
}

/// `Path(path).resolve(strict=True)` for an absolute path. A failure carries
/// the path Python's error message shows.
pub(crate) fn resolve_strict(path: &Path) -> Result<PathBuf, (io::Error, PathBuf)> {
    match fs::canonicalize(path) {
        Ok(real) => Ok(real),
        Err(system_error) => {
            let text = path.to_string_lossy().into_owned();
            match join_real_path(String::new(), &text, &mut HashMap::new()) {
                Err((error, named)) => Err((error, PathBuf::from(named))),
                // The walk found nothing wrong where the system did: keep the system's error.
                Ok(_) => Err((system_error, path.to_path_buf())),
            }
        }
    }
}

// ---------------------------------------------------------------- gzip

/// How `gzip.GzipFile` fails on a damaged file.
#[derive(Debug)]
pub(crate) enum GzipFailure {
    /// `EOFError`: the file ends before the end of a member.
    EndOfFile,
    /// `gzip.BadGzipFile`, an `OSError`.
    Bad(String),
    /// `zlib.error`
    Zlib(String),
}

pub(crate) const COMPRESSED_FILE_ENDED: &str =
    "Compressed file ended before the end-of-stream marker was reached";

impl fmt::Display for GzipFailure {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GzipFailure::EndOfFile => write!(formatter, "{COMPRESSED_FILE_ENDED}"),
            GzipFailure::Bad(text) | GzipFailure::Zlib(text) => write!(formatter, "{text}"),
        }
    }
}

impl std::error::Error for GzipFailure {}

impl GzipFailure {
    /// The failure inside an error this reader returned, if it is one.
    pub(crate) fn of(error: &io::Error) -> Option<&GzipFailure> {
        error
            .get_ref()
            .and_then(|inner| inner.downcast_ref::<GzipFailure>())
    }

    fn into_error(self) -> io::Error {
        let kind = match self {
            GzipFailure::EndOfFile => io::ErrorKind::UnexpectedEof,
            _ => io::ErrorKind::InvalidData,
        };
        io::Error::new(kind, self)
    }
}

/// `gzip.READ_BUFFER_SIZE`
const READ_BLOCK: usize = 128 * 1024;

#[derive(PartialEq, Eq)]
enum State {
    /// Before the header of the next member.
    NewMember,
    /// Inside a member's deflate data.
    InMember,
    /// The deflate data has ended; the trailer is still to be read.
    MemberEnded,
    Finished,
}

/// Concatenated gzip members read like `gzip.GzipFile`, with Python's checks:
/// the trailer of a member is checked when the next read after its last data
/// comes (so earlier lines can be used first), zero padding between members is
/// skipped, and anything else after a member must be another member.
pub(crate) struct PythonGzipReader<R: Read> {
    source: R,
    block: Vec<u8>,
    start: usize,
    end: usize,
    state: State,
    inflater: Decompress,
    crc: Crc,
    produced: u64,
}

impl<R: Read> PythonGzipReader<R> {
    pub(crate) fn new(source: R) -> PythonGzipReader<R> {
        PythonGzipReader {
            source,
            block: vec![0; READ_BLOCK],
            start: 0,
            end: 0,
            state: State::NewMember,
            inflater: Decompress::new(false),
            crc: Crc::new(),
            produced: 0,
        }
    }

    /// Makes at least one unread byte available; false at the end of the file.
    fn fill(&mut self) -> io::Result<bool> {
        if self.start < self.end {
            return Ok(true);
        }
        loop {
            match self.source.read(&mut self.block) {
                Ok(0) => return Ok(false),
                Ok(length) => {
                    self.start = 0;
                    self.end = length;
                    return Ok(true);
                }
                Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
                Err(error) => return Err(error),
            }
        }
    }

    /// `fp.read(1)`
    fn next_byte(&mut self) -> io::Result<Option<u8>> {
        if !self.fill()? {
            return Ok(None);
        }
        let byte = self.block[self.start];
        self.start += 1;
        Ok(Some(byte))
    }

    /// `_read_exact`: exactly `count` bytes, or `EOFError`.
    fn read_exact_bytes(&mut self, count: usize) -> io::Result<Vec<u8>> {
        let mut bytes = Vec::with_capacity(count);
        while bytes.len() < count {
            match self.next_byte()? {
                Some(byte) => bytes.push(byte),
                None => return Err(GzipFailure::EndOfFile.into_error()),
            }
        }
        Ok(bytes)
    }

    /// `_read_gzip_header`: false when the file ends cleanly before a member.
    fn read_header(&mut self) -> io::Result<bool> {
        let mut magic = Vec::with_capacity(2);
        while magic.len() < 2 {
            match self.next_byte()? {
                Some(byte) => magic.push(byte),
                None => break,
            }
        }
        if magic.is_empty() {
            return Ok(false);
        }
        if magic != [0x1f, 0x8b] {
            return Err(
                GzipFailure::Bad(format!("Not a gzipped file ({})", bytes_repr(&magic)))
                    .into_error(),
            );
        }
        let fixed = self.read_exact_bytes(8)?;
        let (method, flags) = (fixed[0], fixed[1]);
        if method != 8 {
            return Err(GzipFailure::Bad("Unknown compression method".to_owned()).into_error());
        }
        if flags & 4 != 0 {
            let length = self.read_exact_bytes(2)?;
            let length = u16::from_le_bytes([length[0], length[1]]) as usize;
            self.read_exact_bytes(length)?;
        }
        for flag in [8, 16] {
            if flags & flag != 0 {
                // A file name or comment ends with a zero byte or with the file.
                loop {
                    match self.next_byte()? {
                        None | Some(0) => break,
                        Some(_) => {}
                    }
                }
            }
        }
        if flags & 2 != 0 {
            self.read_exact_bytes(2)?;
        }
        Ok(true)
    }

    /// `_read_eof`: the trailer of the member that just ended, then the padding.
    fn read_trailer(&mut self) -> io::Result<()> {
        let trailer = self.read_exact_bytes(8)?;
        let expected_crc = u32::from_le_bytes([trailer[0], trailer[1], trailer[2], trailer[3]]);
        let expected_size = u32::from_le_bytes([trailer[4], trailer[5], trailer[6], trailer[7]]);
        if expected_crc != self.crc.sum() {
            return Err(GzipFailure::Bad(format!(
                "CRC check failed {:#x} != {:#x}",
                expected_crc,
                self.crc.sum()
            ))
            .into_error());
        }
        if expected_size != (self.produced & 0xffff_ffff) as u32 {
            return Err(
                GzipFailure::Bad("Incorrect length of data produced".to_owned()).into_error(),
            );
        }
        // Zero bytes may pad a gzip file; the first other byte starts the next member.
        loop {
            match self.next_byte()? {
                Some(0) => {}
                Some(_) => {
                    self.start -= 1;
                    return Ok(());
                }
                None => return Ok(()),
            }
        }
    }
}

impl<R: Read> Read for PythonGzipReader<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        loop {
            match self.state {
                State::Finished => return Ok(0),
                State::MemberEnded => {
                    self.read_trailer()?;
                    self.state = State::NewMember;
                }
                State::NewMember => {
                    if !self.read_header()? {
                        self.state = State::Finished;
                        return Ok(0);
                    }
                    self.inflater = Decompress::new(false);
                    self.crc = Crc::new();
                    self.produced = 0;
                    self.state = State::InMember;
                }
                State::InMember => {
                    if !self.fill()? {
                        return Err(GzipFailure::EndOfFile.into_error());
                    }
                    let in_before = self.inflater.total_in();
                    let out_before = self.inflater.total_out();
                    let status = self
                        .inflater
                        .decompress(
                            &self.block[self.start..self.end],
                            out,
                            FlushDecompress::None,
                        )
                        .map_err(|_| {
                            GzipFailure::Zlib(
                                "Error -3 while decompressing data: invalid deflate data"
                                    .to_owned(),
                            )
                            .into_error()
                        })?;
                    let used = (self.inflater.total_in() - in_before) as usize;
                    let made = (self.inflater.total_out() - out_before) as usize;
                    self.start += used;
                    self.crc.update(&out[..made]);
                    self.produced += made as u64;
                    if status == Status::StreamEnd {
                        self.state = State::MemberEnded;
                    }
                    if made > 0 {
                        return Ok(made);
                    }
                    if status != Status::StreamEnd && used == 0 {
                        // No input used and nothing produced: the decoder is stuck.
                        return Err(GzipFailure::Zlib(
                            "Error -5 while decompressing data: incomplete or truncated stream"
                                .to_owned(),
                        )
                        .into_error());
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;

    fn gzip(data: &[u8]) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(data).unwrap();
        encoder.finish().unwrap()
    }

    fn read_all(bytes: &[u8]) -> io::Result<Vec<u8>> {
        let mut out = Vec::new();
        PythonGzipReader::new(bytes).read_to_end(&mut out)?;
        Ok(out)
    }

    fn failure_text(error: &io::Error) -> String {
        GzipFailure::of(error)
            .map(|f| f.to_string())
            .unwrap_or_else(|| error.to_string())
    }

    #[test]
    fn utf8_errors_read_like_python() {
        assert_eq!(decode_utf8(b"abc"), Ok("abc"));
        assert_eq!(
            decode_utf8(b"##source=\xff"),
            Err("'utf-8' codec can't decode byte 0xff in position 9: invalid start byte".into())
        );
        assert_eq!(
            decode_utf8(b"a\xe2\x82"),
            Err("'utf-8' codec can't decode bytes in position 1-2: unexpected end of data".into())
        );
        assert_eq!(
            decode_utf8(b"\xe2\x82A"),
            Err(
                "'utf-8' codec can't decode bytes in position 0-1: invalid continuation byte"
                    .into()
            )
        );
        assert_eq!(
            decode_utf8(b"x\xc3("),
            Err(
                "'utf-8' codec can't decode byte 0xc3 in position 1: invalid continuation byte"
                    .into()
            )
        );
        assert_eq!(
            decode_utf8(b"x\xe2"),
            Err(
                "'utf-8' codec can't decode byte 0xe2 in position 1: unexpected end of data".into()
            )
        );
    }

    #[test]
    fn bytes_repr_follows_python() {
        assert_eq!(bytes_repr(b"##"), "b'##'");
        assert_eq!(bytes_repr(b"ga"), "b'ga'");
        assert_eq!(bytes_repr(b"\x1f\x00"), "b'\\x1f\\x00'");
        assert_eq!(bytes_repr(b"a'b"), "b\"a'b\"");
        assert_eq!(bytes_repr(b"a'\"b\n"), "b'a\\'\"b\\n'");
    }

    #[test]
    fn a_missing_path_names_its_first_missing_component() {
        let base = std::env::temp_dir().join(format!("priva_python_io_{}", std::process::id()));
        let _ = fs::remove_dir_all(&base);
        fs::create_dir_all(&base).unwrap();
        let (error, named) = resolve_strict(&base.join("nothing/deeper/file")).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::NotFound);
        assert_eq!(named, base.join("nothing"));
        let (_, named) =
            resolve_strict(&base.join("..").join("priva_python_io_missing_x")).unwrap_err();
        assert!(named.ends_with("priva_python_io_missing_x"));
        fs::write(base.join("plain"), b"x").unwrap();
        let (error, named) = resolve_strict(&base.join("plain/inside")).unwrap_err();
        assert_eq!(error.raw_os_error(), Some(20));
        assert_eq!(named, base.join("plain/inside"));
        std::os::unix::fs::symlink(base.join("elsewhere"), base.join("link")).unwrap();
        let (_, named) = resolve_strict(&base.join("link/inner")).unwrap_err();
        assert_eq!(named, base.join("elsewhere"));
        assert!(resolve_strict(&base.join("plain")).is_ok());
        let _ = fs::remove_dir_all(&base);
    }

    #[test]
    fn members_are_joined_and_padding_is_skipped() {
        let mut bytes = gzip(b"first\n");
        bytes.extend_from_slice(&[0, 0, 0]);
        bytes.extend_from_slice(&gzip(b"second\n"));
        assert_eq!(read_all(&bytes).unwrap(), b"first\nsecond\n");
        assert_eq!(read_all(&gzip(b"")).unwrap(), b"");
    }

    #[test]
    fn damaged_files_fail_with_pythons_texts() {
        let good = gzip(b"hello\n");
        // The trailer is cut off.
        let error = read_all(&good[..good.len() - 4]).unwrap_err();
        assert_eq!(failure_text(&error), COMPRESSED_FILE_ENDED);
        // Data is cut off.
        let error = read_all(&good[..good.len() - 12]).unwrap_err();
        assert_eq!(failure_text(&error), COMPRESSED_FILE_ENDED);
        // The header is cut off.
        let error = read_all(&good[..6]).unwrap_err();
        assert_eq!(failure_text(&error), COMPRESSED_FILE_ENDED);
        // Garbage after the member.
        let mut garbage = good.clone();
        garbage.extend_from_slice(b"garbage");
        assert_eq!(
            failure_text(&read_all(&garbage).unwrap_err()),
            "Not a gzipped file (b'ga')"
        );
        // Unknown method.
        let mut method = good.clone();
        method[2] = 7;
        assert_eq!(
            failure_text(&read_all(&method).unwrap_err()),
            "Unknown compression method"
        );
        // Not gzip at all.
        assert_eq!(
            failure_text(&read_all(b"##x").unwrap_err()),
            "Not a gzipped file (b'##')"
        );
        // A wrong checksum.
        let mut checksum = good.clone();
        let at = checksum.len() - 8;
        checksum[at] ^= 0xff;
        assert!(failure_text(&read_all(&checksum).unwrap_err()).starts_with("CRC check failed 0x"));
        // A wrong length.
        let mut size = good.clone();
        let at = size.len() - 1;
        size[at] ^= 0x01;
        assert_eq!(
            failure_text(&read_all(&size).unwrap_err()),
            "Incorrect length of data produced"
        );
    }

    #[test]
    fn the_trailer_is_checked_after_the_data_is_handed_out() {
        let mut truncated = gzip(b"line one\nline two\n");
        truncated.truncate(truncated.len() - 3);
        let mut reader = PythonGzipReader::new(&truncated[..]);
        let mut buffer = [0u8; 64];
        let first = reader.read(&mut buffer).unwrap();
        assert_eq!(&buffer[..first], b"line one\nline two\n");
        assert!(reader.read(&mut buffer).is_err());
    }
}
