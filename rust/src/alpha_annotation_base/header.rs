//! Bounded alphamissense_cache.check_annotation semantics, independent of
//! pysam record readers. Plain VCF and ordinary gzip are source-supported.
use super::io::*;
use super::profile::{assembly, major};
use super::{error, object, require};
use crate::cli::Result;
use crate::json_text::Json;
use crate::python_io::{GzipFailure, PythonGzipReader};
use crate::python_path::{absolute, is_file, metadata, os_error};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Read};
use std::os::unix::fs::MetadataExt;
use std::path::{Path, PathBuf};

type State = (u64, u64, u32, u64, u64, i64, i64, i64, i64);
fn fields(m: fs::Metadata) -> State {
    (
        m.dev(),
        m.ino(),
        m.mode(),
        m.nlink(),
        m.len(),
        m.mtime(),
        m.mtime_nsec(),
        m.ctime(),
        m.ctime_nsec(),
    )
}
fn state(path: &Path) -> Result<(State, State)> {
    Ok((
        fields(fs::symlink_metadata(path).map_err(|e| os_error(e, path, None))?),
        fields(metadata(path)?),
    ))
}
fn gzip_error(cause: io::Error) -> Box<dyn std::error::Error> {
    let name = match GzipFailure::of(&cause) {
        Some(GzipFailure::EndOfFile) => "EOFError",
        Some(GzipFailure::Bad(_)) => "BadGzipFile",
        Some(GzipFailure::Zlib(_)) => "error",
        None => "OSError",
    };
    error(name, cause.to_string())
}
fn split(value: &str) -> Result<Vec<String>> {
    let mut chars = value.chars().peekable();
    let mut words = Vec::new();
    let mut word = String::new();
    let mut active = false;
    let mut quote = None;
    while let Some(c) = chars.next() {
        match quote {
            Some(q) if c == q => quote = None,
            Some('"') if c == '\\' => {
                let Some(next) = chars.next() else {
                    return Err(error("ValueError", "No escaped character"));
                };
                if next != '"' && next != '\\' {
                    word.push('\\');
                }
                word.push(next);
            }
            Some(_) => word.push(c),
            None if c == '\'' || c == '"' => {
                quote = Some(c);
                active = true;
            }
            None if c == '\\' => {
                word.push(
                    chars
                        .next()
                        .ok_or_else(|| error("ValueError", "No escaped character"))?,
                );
                active = true;
            }
            None if matches!(c, ' ' | '\t' | '\r' | '\n') => {
                if active {
                    words.push(std::mem::take(&mut word));
                    active = false;
                }
            }
            None => {
                word.push(c);
                active = true;
            }
        }
    }
    if quote.is_some() {
        return Err(error("ValueError", "No closing quotation"));
    }
    if active {
        words.push(word);
    }
    Ok(words)
}
fn check_header(
    path: &Path,
    selected: &str,
    requested: &str,
    raw_assembly: &Json,
    resource: &str,
) -> Result<Json> {
    let mut first = File::open(path).map_err(|e| os_error(e, path, None))?;
    let mut magic = [0; 2];
    let n = first
        .read(&mut magic)
        .map_err(|e| os_error(e, path, None))?;
    drop(first);
    let file = File::open(path).map_err(|e| os_error(e, path, None))?;
    let reader: Box<dyn Read> = if n == 2 && magic == [31, 139] {
        Box::new(PythonGzipReader::new(file))
    } else {
        Box::new(file)
    };
    let mut reader = BufReader::new(reader);
    let mut headers = Vec::new();
    let mut remaining = 4 * 1024 * 1024usize;
    let mut complete = false;
    for _ in 0..50000 {
        let mut line = Vec::new();
        let mut limited = (&mut reader).take((remaining + 1) as u64);
        limited.read_until(b'\n', &mut line).map_err(gzip_error)?;
        require(
            !line.is_empty() && line.len() <= remaining && line.ends_with(b"\n"),
            format!("{resource} annotation has no complete bounded VCF header"),
        )?;
        remaining -= line.len();
        if line.starts_with(b"##VEP=") {
            while line.last().is_some_and(|c| matches!(c, b'\r' | b'\n')) {
                line.pop();
            }
            headers.push(
                crate::python_io::decode_utf8(&line)
                    .map_err(|e| error("UnicodeDecodeError", e))?
                    .to_owned(),
            );
        }
        if line.starts_with(b"#CHROM\t") {
            complete = true;
            break;
        }
        require(
            line.starts_with(b"##"),
            format!("{resource} annotation has no complete VCF column header"),
        )?;
    }
    require(
        complete,
        format!("{resource} annotation VCF header exceeds the bounded line limit"),
    )?;
    require(
        headers.len() == 1,
        format!("{resource} annotation requires exactly one unambiguous ##VEP header"),
    )?;
    let mut versions = Vec::new();
    let mut assemblies = Vec::new();
    for token in split(&headers[0])? {
        if let Some((k, v)) = token.split_once('=') {
            if k == "##VEP" {
                versions.push(v.to_owned());
            }
            if k == "assembly" {
                assemblies.push(v.to_owned());
            }
        }
    }
    require(
        versions.len() == 1 && assemblies.len() == 1,
        format!("{resource} annotation requires unambiguous VEP release and assembly metadata"),
    )?;
    let version = &versions[0];
    let valid = version.strip_prefix('v').is_some_and(|v| {
        let mut parts = v.split('.');
        let first = parts.next().unwrap();
        major(first)
            && first == requested
            && parts.all(|p| !p.is_empty() && p.bytes().all(|c| c.is_ascii_digit()))
    });
    require(
        valid,
        format!(
            "{resource} annotation VEP release {version} differs from requested major {requested}"
        ),
    )?;
    let observed = &assemblies[0];
    let normalized = observed
        .split_once(".p")
        .filter(|(base, patch)| matches!(*base, "GRCh37" | "GRCh38") && major(patch))
        .map_or(observed.as_str(), |(base, _)| base);
    require(
        assembly(&string(normalized))? == selected,
        format!(
            "{resource} annotation assembly {observed} differs from requested {}",
            text(raw_assembly)
        ),
    )?;
    Ok(object([
        ("path", pstr(path)),
        ("assembly", string(selected)),
        ("vep_major", string(requested)),
        ("vep_version", string(&version[1..])),
    ]))
}
pub(crate) fn check(path: &Path, raw_assembly: &Json, vep_major: &Json) -> Result<Json> {
    check_resource(path, raw_assembly, vep_major, "AlphaMissense")
}
/// `check_vep_annotation(..., resource=resource)`: the same bounded check for
/// another derived annotation, such as the managed ClinVar VCF.
pub(crate) fn check_resource(
    path: &Path,
    raw_assembly: &Json,
    vep_major: &Json,
    resource: &str,
) -> Result<Json> {
    let selected = assembly(raw_assembly)?;
    let requested = text(vep_major);
    require(
        major(&requested),
        "Requested VEP major must be a positive integer",
    )?;
    let path = PathBuf::from(absolute(path)?);
    require(
        is_file(&path)? && metadata(&path)?.len() > 0,
        format!("Missing/empty {resource} annotation: {}", path.display()),
    )?;
    let before = state(&path)?;
    let value = check_header(&path, selected, &requested, raw_assembly, resource);
    if state(&path).ok() != Some(before) {
        return Err(error(
            "AnnotationHeaderChangedError",
            format!("{resource} annotation changed during its header check"),
        ));
    }
    value
}
