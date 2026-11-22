//! The Python file behaviour these steps show: `Path.read_text`, text-mode
//! line iteration (UTF-8 decoded in 8 KiB blocks with universal newlines),
//! `os.stat` identities, streamed digests, gzip prefixes, `Path.mkdir(
//! parents=True, exist_ok=True)` and link-based publication.
use super::{Result, error, require};
use crate::installer_prerequisites::paths::{is_dir, open, parent, read_failure};
use crate::json_text::Json;
use crate::md5::Md5;
use crate::python_io::{PythonGzipReader, decode_utf8};
use crate::python_path::{metadata, os_error};
use crate::sha256::Sha256;
use std::fs;
use std::io::{self, Read};
use std::os::unix::fs::MetadataExt;
use std::path::Path;

pub(super) const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
/// The script's digest and decompression block.
pub(super) const BLOCK: usize = 4 * 1024 * 1024;
/// `TextIOWrapper._CHUNK_SIZE`
const TEXT_BLOCK: usize = 8192;
const NANOSECONDS: i128 = 1_000_000_000;

/// Universal newlines: CRLF and lone CR become LF.
fn translate(text: &str) -> String {
    text.replace("\r\n", "\n").replace('\r', "\n")
}

fn decode_error(bytes: &[u8]) -> Box<dyn std::error::Error> {
    error(
        "UnicodeDecodeError",
        decode_utf8(bytes).err().unwrap_or_default(),
    )
}

/// `Path(path).read_text()`: the whole file decoded at once.
pub(super) fn read_text(path: &str) -> Result<String> {
    let bytes = fs::read(path).map_err(|failure| os_error(failure, Path::new(path), None))?;
    match std::str::from_utf8(&bytes) {
        Ok(text) => Ok(translate(text)),
        Err(_) => Err(decode_error(&bytes)),
    }
}

/// One read call, retried when interrupted (`read1`).
fn read_once(source: &mut impl Read, block: &mut [u8]) -> io::Result<usize> {
    loop {
        match source.read(block) {
            Err(failure) if failure.kind() == io::ErrorKind::Interrupted => continue,
            other => return other,
        }
    }
}

/// Up to `block.len()` bytes, fewer only at the end (`read(n)`).
pub(super) fn read_block(source: &mut impl Read, block: &mut [u8]) -> io::Result<usize> {
    let mut filled = 0;
    while filled < block.len() {
        let count = read_once(source, &mut block[filled..])?;
        if count == 0 {
            break;
        }
        filled += count;
    }
    Ok(filled)
}

/// Text-mode iteration over a stream, as TextIOWrapper decodes it: each block
/// is decoded with the bytes an earlier block left incomplete, a decoding
/// failure names its position in that buffer, and a CR at the end of a block
/// waits for the next one before newlines are translated.
pub(super) struct TextLines<R: Read> {
    source: R,
    path: String,
    pending: Vec<u8>,
    carriage: bool,
    text: String,
    at: usize,
    finished: bool,
}

impl<R: Read> TextLines<R> {
    pub(super) fn new(source: R, path: &str) -> Self {
        Self {
            source,
            path: path.to_owned(),
            pending: Vec::new(),
            carriage: false,
            text: String::new(),
            at: 0,
            finished: false,
        }
    }

    fn decode_block(&mut self) -> Result<()> {
        let mut block = vec![0; TEXT_BLOCK];
        let count = read_once(&mut self.source, &mut block)
            .map_err(|failure| read_failure(Box::new(failure), &self.path))?;
        let last = count == 0;
        let mut data = std::mem::take(&mut self.pending);
        data.extend_from_slice(&block[..count]);
        let mut output = match std::str::from_utf8(&data) {
            Ok(text) => text.to_owned(),
            Err(failure) if failure.error_len().is_none() && !last => {
                let valid = failure.valid_up_to();
                self.pending = data[valid..].to_vec();
                String::from_utf8_lossy(&data[..valid]).into_owned()
            }
            Err(_) => return Err(decode_error(&data)),
        };
        if self.carriage && (!output.is_empty() || last) {
            output.insert(0, '\r');
            self.carriage = false;
        }
        if output.ends_with('\r') && !last {
            output.pop();
            self.carriage = true;
        }
        self.text.replace_range(..self.at, "");
        self.at = 0;
        self.text.push_str(&translate(&output));
        self.finished = last;
        Ok(())
    }

    /// The next line with its LF, or `None` at the end.
    pub(super) fn readline(&mut self) -> Result<Option<String>> {
        loop {
            if let Some(offset) = self.text[self.at..].find('\n') {
                let line = self.text[self.at..self.at + offset + 1].to_owned();
                self.at += offset + 1;
                return Ok(Some(line));
            }
            if self.finished {
                if self.at == self.text.len() {
                    return Ok(None);
                }
                let line = self.text[self.at..].to_owned();
                self.at = self.text.len();
                return Ok(Some(line));
            }
            self.decode_block()?;
        }
    }
}

/// `[st_dev, st_ino, st_size, st_mtime_ns, st_ctime_ns]` of `path.stat()`.
pub(super) type Identity = [i128; 5];

pub(super) fn identity(path: &str) -> Result<Identity> {
    let stat = metadata(Path::new(path))?;
    Ok([
        i128::from(stat.dev()),
        i128::from(stat.ino()),
        i128::from(stat.size()),
        i128::from(stat.mtime()) * NANOSECONDS + i128::from(stat.mtime_nsec()),
        i128::from(stat.ctime()) * NANOSECONDS + i128::from(stat.ctime_nsec()),
    ])
}

pub(super) fn identity_json(identity: &Identity) -> Json {
    Json::Array(identity.iter().map(|value| Json::Int(*value)).collect())
}

/// `path.stat().st_mtime_ns`
pub(super) fn mtime_ns(path: &str) -> Result<i128> {
    Ok(identity(path)?[3])
}

/// `path.stat().st_size`
pub(super) fn file_size(path: &str) -> Result<u64> {
    Ok(metadata(Path::new(path))?.len())
}

#[derive(Clone, Copy)]
pub(super) enum Algorithm {
    Sha256,
    Md5,
}

/// `digest(path, algorithm)`: the hexadecimal digest, refused when the file
/// changed while it was read.
pub(super) fn digest(path: &str, algorithm: Algorithm) -> Result<String> {
    let before = identity(path)?;
    let mut file = open(path)?;
    let mut sha256 = Sha256::new();
    let mut md5 = Md5::new();
    let mut block = vec![0; BLOCK];
    loop {
        let count = read_block(&mut file, &mut block)
            .map_err(|failure| read_failure(Box::new(failure), path))?;
        if count == 0 {
            break;
        }
        match algorithm {
            Algorithm::Sha256 => sha256.update(&block[..count]),
            Algorithm::Md5 => md5.update(&block[..count]),
        }
    }
    require(
        identity(path)? == before,
        format!("Reference resource changed during verification: {path}; wait for its writer"),
    )?;
    Ok(match algorithm {
        Algorithm::Sha256 => sha256.finish_hex(),
        Algorithm::Md5 => md5.finish_hex(),
    })
}

/// `path.open("rb").read(count)`
pub(super) fn prefix(path: &str, count: usize) -> Result<Vec<u8>> {
    let mut block = vec![0; count];
    let read = read_block(&mut open(path)?, &mut block)
        .map_err(|failure| read_failure(Box::new(failure), path))?;
    block.truncate(read);
    Ok(block)
}

/// `gzip.open(path, "rb").read(count)`
pub(super) fn gzip_prefix(path: &str, count: usize) -> Result<Vec<u8>> {
    let mut block = vec![0; count];
    let mut reader = PythonGzipReader::new(open(path)?);
    let read = read_block(&mut reader, &mut block)
        .map_err(|failure| read_failure(Box::new(failure), path))?;
    block.truncate(read);
    Ok(block)
}

/// Whether the file starts with the gzip magic bytes.
pub(super) fn compressed(path: &str) -> Result<bool> {
    Ok(prefix(path, 2)? == GZIP_MAGIC)
}

/// `Path(path).mkdir(parents=True, exist_ok=True)`
pub(super) fn make_directories(path: &str) -> Result<()> {
    make_directory(path, true)
}

fn make_directory(path: &str, parents: bool) -> Result<()> {
    match fs::create_dir(path) {
        Ok(()) => Ok(()),
        Err(failure) if failure.kind() == io::ErrorKind::NotFound => {
            let above = parent(path);
            if !parents || above == path {
                return Err(os_error(failure, Path::new(path), None));
            }
            make_directory(&above, true)?;
            make_directory(path, false)
        }
        Err(failure) => {
            if is_dir(path)? {
                Ok(())
            } else {
                Err(os_error(failure, Path::new(path), None))
            }
        }
    }
}

/// `Path(path).unlink(missing_ok=True)`
pub(super) fn unlink_missing_ok(path: &str) -> Result<()> {
    match fs::remove_file(path) {
        Err(failure) if failure.kind() != io::ErrorKind::NotFound => {
            Err(os_error(failure, Path::new(path), None))
        }
        _ => Ok(()),
    }
}

/// `publish_new(partial, destination)`: the link fails if a concurrent
/// process supplied the final file, which is never replaced.
pub(super) fn publish_new(partial: &str, destination: &str) -> Result<()> {
    fs::hard_link(partial, destination)
        .map_err(|failure| os_error(failure, Path::new(partial), Some(Path::new(destination))))?;
    fs::remove_file(partial).map_err(|failure| os_error(failure, Path::new(partial), None))
}

/// `os.replace(source, destination)`
pub(super) fn replace(source: &str, destination: &str) -> Result<()> {
    fs::rename(source, destination)
        .map_err(|failure| os_error(failure, Path::new(source), Some(Path::new(destination))))
}

/// `Path(path).write_text(text)`
pub(super) fn write_text(path: &str, text: &str) -> Result<()> {
    fs::write(path, text).map_err(|failure| os_error(failure, Path::new(path), None))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lines(bytes: &[u8]) -> Result<Vec<String>> {
        let mut reader = TextLines::new(bytes, "fixture");
        let mut found = Vec::new();
        while let Some(line) = reader.readline()? {
            found.push(line);
        }
        Ok(found)
    }

    fn failure_text(failure: Box<dyn std::error::Error>) -> String {
        failure.to_string()
    }

    #[test]
    fn text_lines_translate_newlines_across_blocks() {
        let mut bytes = vec![b'a'; TEXT_BLOCK - 1];
        bytes.extend_from_slice(b"\r\nb\rc\n\xc3\xa9");
        let found = lines(&bytes).unwrap();
        assert_eq!(found.len(), 4);
        assert_eq!(found[0].len(), TEXT_BLOCK);
        assert_eq!(found[1..], ["b\n", "c\n", "\u{e9}"]);
        assert_eq!(lines(b"x\r").unwrap(), ["x\n"]);
        assert!(lines(b"").unwrap().is_empty());
    }

    #[test]
    fn decoding_failures_name_their_position_in_the_block() {
        let failure = lines(b"ok\n\xff\n").unwrap_err();
        assert_eq!(
            failure_text(failure),
            "UnicodeDecodeError: 'utf-8' codec can't decode byte 0xff in position 3: \
             invalid start byte"
        );
        let failure = lines(b"ok\n\xe2\x82").unwrap_err();
        assert_eq!(
            failure_text(failure),
            "UnicodeDecodeError: 'utf-8' codec can't decode bytes in position 0-1: \
             unexpected end of data"
        );
        // A failure in a later block stays unread while earlier lines suffice.
        let mut late = vec![b'a'; TEXT_BLOCK - 1];
        late.extend_from_slice(b"\n\xff");
        let mut reader = TextLines::new(late.as_slice(), "fixture");
        assert_eq!(reader.readline().unwrap().unwrap().len(), TEXT_BLOCK);
        assert!(reader.readline().is_err());
    }

    #[test]
    fn whole_file_text_reports_absolute_positions() {
        let name = format!("priva-reference-text-{}", std::process::id());
        let path = std::env::temp_dir().join(name);
        fs::write(&path, b"a\r\nb\xff").unwrap();
        let shown = path.display().to_string();
        let failure = read_text(&shown).unwrap_err();
        assert!(failure_text(failure).ends_with("byte 0xff in position 4: invalid start byte"));
        fs::write(&path, b"a\r\nb\rc").unwrap();
        assert_eq!(read_text(&shown).unwrap(), "a\nb\nc");
        fs::remove_file(path).unwrap();
    }
}
