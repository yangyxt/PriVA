//! The evidence builder's pandas text profile and deterministic sibling writer.
//! Existing generic readers retain their own options and refusal boundaries.
use super::{error, exists};
use crate::cli::Result;
use crate::pandas_table::{self, Column, Table};
use crate::python_path::{metadata, os_error};
use std::fs::{self, File};
use std::io::{BufWriter, Read, Seek, Write};
use std::path::Path;

pub(super) fn read(path: &Path, delimiter: u8, skip_bad: bool) -> Result<Table> {
    let file = File::open(path).map_err(|e| os_error(e, path, None))?;
    if file.metadata()?.is_dir() {
        return Err(os_error(
            std::io::Error::from_raw_os_error(libc::EISDIR),
            path,
            None,
        ));
    }
    let name = path.to_string_lossy().to_lowercase();
    if name.ends_with(".zip") {
        return decode_and_parse(&zip_bytes(file, path)?, delimiter, skip_bad);
    }
    if [".tar", ".tar.gz", ".tar.bz2", ".tar.xz"]
        .iter()
        .any(|suffix| name.ends_with(suffix))
    {
        return decode_and_parse(&tar_bytes(file, path)?, delimiter, skip_bad);
    }
    let mut source: Box<dyn Read> = if name.ends_with(".gz") && !name.ends_with(".tar.gz") {
        Box::new(crate::python_io::PythonGzipReader::new(file))
    } else if name.ends_with(".bz2") && !name.ends_with(".tar.bz2") {
        Box::new(bzip2::read::MultiBzDecoder::new(file))
    } else if name.ends_with(".xz") && !name.ends_with(".tar.xz") {
        Box::new(xz2::read::XzDecoder::new_multi_decoder(file))
    } else if name.ends_with(".zst") {
        Box::new(zstd::stream::read::Decoder::new(file)?)
    } else {
        Box::new(file)
    };
    let mut bytes = Vec::new();
    source
        .read_to_end(&mut bytes)
        .map_err(|e| match crate::python_io::GzipFailure::of(&e) {
            Some(crate::python_io::GzipFailure::EndOfFile) => error("EOFError", e.to_string()),
            Some(crate::python_io::GzipFailure::Bad(_)) => error("BadGzipFile", e.to_string()),
            Some(crate::python_io::GzipFailure::Zlib(_)) => error("error", e.to_string()),
            None => error("OSError", e.to_string()),
        })?;
    decode_and_parse(&bytes, delimiter, skip_bad)
}

fn member_names(names: &[String]) -> String {
    format!(
        "[{}]",
        names
            .iter()
            .map(|name| crate::python_path::repr(name))
            .collect::<Vec<_>>()
            .join(", ")
    )
}
fn zip_error(cause: zip::result::ZipError) -> Box<dyn std::error::Error> {
    match cause {
        zip::result::ZipError::InvalidArchive("Could not find central directory end") => {
            error("BadZipFile", "File is not a zip file")
        }
        cause => error("BadZipFile", cause.to_string()),
    }
}
fn zip_bytes(file: File, path: &Path) -> Result<Vec<u8>> {
    let mut archive = zip::ZipArchive::new(file).map_err(zip_error)?;
    let names = (0..archive.len())
        .map(|i| {
            archive
                .by_index(i)
                .map(|entry| entry.name().to_owned())
                .map_err(zip_error)
        })
        .collect::<Result<Vec<_>>>()?;
    if names.is_empty() {
        return Err(error(
            "ValueError",
            format!("Zero files found in ZIP file {}", path.display()),
        ));
    }
    if names.len() != 1 {
        return Err(error(
            "ValueError",
            format!(
                "Multiple files found in ZIP file. Only one file per ZIP: {}",
                member_names(&names)
            ),
        ));
    }
    let mut bytes = Vec::new();
    archive
        .by_index(0)
        .map_err(zip_error)?
        .read_to_end(&mut bytes)?;
    Ok(bytes)
}
fn tar_bytes(mut file: File, path: &Path) -> Result<Vec<u8>> {
    let mut magic = [0u8; 6];
    let used = file.read(&mut magic)?;
    file.rewind()?;
    let reader: Box<dyn Read> = if magic[..used].starts_with(&[0x1f, 0x8b]) {
        Box::new(crate::python_io::PythonGzipReader::new(file))
    } else if magic[..used].starts_with(b"BZh") {
        Box::new(bzip2::read::MultiBzDecoder::new(file))
    } else if magic[..used].starts_with(&[0xfd, b'7', b'z', b'X', b'Z', 0]) {
        Box::new(xz2::read::XzDecoder::new_multi_decoder(file))
    } else {
        Box::new(file)
    };
    let mut archive = tar::Archive::new(reader);
    let mut names = Vec::new();
    let mut first = Vec::new();
    let mut regular = false;
    for entry in archive
        .entries()
        .map_err(|e| error("ReadError", e.to_string()))?
    {
        let mut entry = entry.map_err(|e| error("ReadError", e.to_string()))?;
        names.push(entry.path()?.to_string_lossy().into_owned());
        if names.len() == 1 {
            regular = entry.header().entry_type().is_file()
                || entry.header().entry_type().is_contiguous();
            entry.read_to_end(&mut first)?;
        }
    }
    if names.is_empty() {
        return Err(error(
            "ValueError",
            format!("Zero files found in TAR archive {}", path.display()),
        ));
    }
    if names.len() != 1 {
        return Err(error(
            "ValueError",
            format!(
                "Multiple files found in TAR archive. Only one file per TAR archive: {}",
                member_names(&names)
            ),
        ));
    }
    if !regular {
        return Err(error("AssertionError", ""));
    }
    Ok(first)
}

/// TextIOWrapper's bounded `read(256 Ki characters)` requests used by pandas.
/// Its incremental UTF-8 decoder reports offsets within the current decode
/// block; an incomplete terminal codepoint is retried alone with final=True.
struct DecodeFailure {
    delivered: String,
    cause: Box<dyn std::error::Error>,
}
fn decode_and_parse(bytes: &[u8], delimiter: u8, skip_bad: bool) -> Result<Table> {
    match decode_text(bytes) {
        Ok(text) => Ok(parse(&text, delimiter, skip_bad, true)?.unwrap()),
        Err(failed) => {
            // pandas tokenizes completed TextIOWrapper.read requests before
            // asking for more text. A bad later decode cannot hide an earlier
            // width error, but a partial read that itself failed was not seen.
            parse(&failed.delivered, delimiter, skip_bad, false)?;
            Err(failed.cause)
        }
    }
}
fn decode_text(bytes: &[u8]) -> std::result::Result<String, DecodeFailure> {
    let mut output = String::new();
    let mut at = 0;
    let mut pending = Vec::new();
    let mut need = 256 * 1024;
    let mut ratio = 1.0_f64;
    let mut initial = true;
    loop {
        let request = ((need as f64 * ratio).ceil() as usize).max(8192);
        let count = (bytes.len() - at).min(request);
        pending.extend_from_slice(&bytes[at..at + count]);
        at += count;
        if initial {
            initial = false;
            if pending.starts_with(&[0xef, 0xbb, 0xbf]) {
                pending.drain(..3);
            }
        }
        let complete = match std::str::from_utf8(&pending) {
            Ok(_) => pending.len(),
            Err(e) if e.error_len().is_none() && count > 0 => e.valid_up_to(),
            Err(_) => {
                let characters = output.chars().count();
                let delivered = characters / (256 * 1024) * (256 * 1024);
                let end = output
                    .char_indices()
                    .nth(delivered)
                    .map_or(output.len(), |(at, _)| at);
                output.truncate(end);
                return Err(DecodeFailure {
                    delivered: output,
                    cause: error(
                        "UnicodeDecodeError",
                        crate::python_io::decode_utf8(&pending).unwrap_err(),
                    ),
                });
            }
        };
        let text = std::str::from_utf8(&pending[..complete]).unwrap();
        let characters = text.chars().count();
        if characters > 0 {
            ratio = count as f64 / characters as f64;
        }
        output.push_str(text);
        pending.drain(..complete);
        need = if characters >= need {
            256 * 1024 - (characters - need) % (256 * 1024)
        } else {
            need - characters
        };
        if count == 0 {
            break;
        }
    }
    Ok(output)
}

fn parse(text: &str, delimiter: u8, skip_bad: bool, final_input: bool) -> Result<Option<Table>> {
    let mut names = None;
    let mut columns: Vec<Vec<Option<String>>> = Vec::new();
    let mut width = None;
    let mut count = 0;
    pandas_table::tokenize_builder_separated(text, delimiter, final_input, |record, line| {
        // pandas' C string conversion stops at embedded NUL, while framing
        // still consumes the complete field (including quotes/delimiters).
        let record = record.into_iter().map(|v| v.split('\0').next().unwrap_or("").to_owned()).collect::<Vec<_>>();
        if names.is_none() {
            columns = vec![Vec::new(); record.len()];
            names = Some(pandas_table::column_names(record));
            return Ok(());
        }
        let expected = *width.get_or_insert_with(|| columns.len().max(record.len()));
        if record.len() > expected {
            if skip_bad { return Ok(()); }
            return Err(error("ParserError", format!("Error tokenizing data. C error: Expected {expected} fields in line {line}, saw {}\n", record.len())));
        }
        let offset = expected - columns.len();
        for (index, column) in columns.iter_mut().enumerate() {
            let value = record.get(index + offset).map_or("", String::as_str);
            column.push(Some(if pandas_table::is_na(value) { String::new() } else { value.to_owned() }));
        }
        count += 1;
        Ok(())
    }).map_err(|e| crate::builder_error::category("ParserError", e))?;
    if !final_input {
        return Ok(None);
    }
    Ok(Some(Table {
        names: names.ok_or_else(|| error("EmptyDataError", "No columns to parse from file"))?,
        columns: columns.into_iter().map(Column::Text).collect(),
        rows: count,
    }))
}

pub(super) fn publish(path: &Path, data: &Table) -> Result<usize> {
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    fs::create_dir_all(parent).map_err(|e| os_error(e, parent, None))?;
    let name = path.file_name().ok_or_else(|| {
        error(
            "ValueError",
            format!(
                "PosixPath({}) has an empty name",
                crate::python_path::repr(&path.display().to_string())
            ),
        )
    })?;
    let mut temporary = name.to_os_string();
    temporary.push(".tmp");
    let temporary = path.with_file_name(temporary);
    let file = File::create(&temporary).map_err(|e| os_error(e, &temporary, None))?;
    let mut stream = BufWriter::new(file);
    pandas_table::write_table_to(&mut stream, data)?;
    stream.flush()?;
    drop(stream);
    if exists(path)?
        && metadata(path)?.len() == metadata(&temporary)?.len()
        && crate::sha256::sha256_file(path).map_err(|e| os_error(e, path, None))?
            == crate::sha256::sha256_file(&temporary).map_err(|e| os_error(e, &temporary, None))?
    {
        fs::remove_file(&temporary).map_err(|e| os_error(e, &temporary, None))?;
    } else {
        fs::rename(&temporary, path).map_err(|e| os_error(e, &temporary, Some(path)))?;
    }
    Ok(data.rows)
}
