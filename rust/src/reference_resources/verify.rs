//! Read-only checks: FASTA indexes, the forward chain, pinned-source receipts
//! and `preflight`, which inspects every planned resource before any lock,
//! directory, index or download.
use super::files::{
    Algorithm, BLOCK, GZIP_MAGIC, TextLines, compressed, digest, file_size, gzip_prefix,
    identity, identity_json, mtime_ns, prefix, read_block, read_text,
};
use super::plan::{Kind, Resource, checked_path};
use super::{Result, error, require, source};
use crate::conservation_resource::python_isdigit;
use crate::installer_prerequisites::checks::{python_int, splitlines};
use crate::installer_prerequisites::paths::{exists, is_file, lexists, open, read_failure};
use crate::json_text::{self, Json, py_equal};
use crate::python_io::PythonGzipReader;
use crate::python_text::is_python_space;
use crate::sha256::Sha256;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};

/// `CHR1_LENGTH[assembly]` as Python prints it.
pub(super) fn chr1_length(assembly: &str) -> &'static str {
    if assembly == "hg19" {
        "249250621"
    } else {
        "248956422"
    }
}

/// `valid_index(path, assembly)`: a current `.fai` whose rows are complete and
/// whose chromosome 1 has the assembly's length. Unreadable or malformed
/// indexes are not valid, rather than failures.
pub(super) fn valid_index(path: &str, assembly: &str) -> Result<bool> {
    let index = checked_path(&format!("{path}.fai"))?;
    if !is_file(&index)? || mtime_ns(&index)? < mtime_ns(path)? {
        return Ok(false);
    }
    Ok(index_identifies(&index, assembly))
}

/// `line.rstrip().split("\t")`
fn fields(line: &str) -> Vec<&str> {
    let kept = line.trim_end_matches(is_python_space);
    kept.split('\t').collect()
}

/// A `.fai` row: five fields, the last four made of digits.
fn complete(row: &[&str]) -> bool {
    row.len() == 5 && row[1..].iter().all(|part| python_isdigit(part))
}

fn index_identifies(index: &str, assembly: &str) -> bool {
    let Ok(text) = read_text(index) else {
        return false;
    };
    let rows: Vec<Vec<&str>> = splitlines(&text).into_iter().map(fields).collect();
    if rows.is_empty() || !rows.iter().all(|row| complete(row)) {
        return false;
    }
    let mut names = HashSet::new();
    if !rows.iter().all(|row| names.insert(row[0])) {
        return false;
    }
    let mut lengths = HashMap::new();
    for row in &rows {
        match python_int(row[1]) {
            Ok(length) => lengths.insert(row[0], length),
            Err(_) => return false,
        };
    }
    let observed = lengths.get("chr1").or_else(|| lengths.get("1"));
    observed.is_some_and(|length| length == chr1_length(assembly))
}

/// `valid_index_at(index, assembly)`: a newly written index; an integer that
/// does not convert is a failure here.
pub(super) fn valid_index_at(index: &str, assembly: &str) -> Result<bool> {
    let mut lines = TextLines::new(open(index)?, index);
    let mut lengths = HashMap::new();
    while let Some(line) = lines.readline()? {
        let row = fields(&line);
        if row.len() == 5 {
            lengths.insert(row[0].to_owned(), python_int(row[1])?);
        }
    }
    let observed = lengths.get("chr1").or_else(|| lengths.get("1"));
    Ok(observed.is_some_and(|length| length == chr1_length(assembly)))
}

/// `bytes.split()` and `bytes.strip()` white space.
fn ascii_space(byte: &u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
}

/// `fasta_chr1(path)`: the bases actually counted under the first `1` or
/// `chr1` record, or `None` without one. A name is never inferred.
fn fasta_chr1(path: &str) -> Result<Option<u64>> {
    let source: Box<dyn Read> = if compressed(path)? {
        Box::new(PythonGzipReader::new(open(path)?))
    } else {
        Box::new(open(path)?)
    };
    let mut reader = BufReader::new(source);
    let mut found = false;
    let mut length = 0;
    let mut line = Vec::new();
    loop {
        line.clear();
        let count = reader
            .read_until(b'\n', &mut line)
            .map_err(|failure| read_failure(Box::new(failure), path))?;
        if count == 0 {
            break;
        }
        if line.starts_with(b">") {
            if found {
                return Ok(Some(length));
            }
            let Some(name) = line[1..].split(ascii_space).find(|word| !word.is_empty()) else {
                return Err(error("IndexError", "list index out of range"));
            };
            found = name == b"1" || name == b"chr1";
        } else if found {
            let start = line.iter().position(|byte| !ascii_space(byte));
            let end = line.iter().rposition(|byte| !ascii_space(byte));
            if let (Some(start), Some(end)) = (start, end) {
                length += (end + 1 - start) as u64;
            }
        }
    }
    Ok(found.then_some(length))
}

/// `validate_fasta(path, assembly)`
fn validate_fasta(path: &str, assembly: &str) -> Result<()> {
    let mut first = prefix(path, 2)?;
    if first == GZIP_MAGIC {
        first = gzip_prefix(path, 1)?;
    }
    require(
        first.starts_with(b">"),
        format!("Not a FASTA reference: {path}; existing bytes retained"),
    )?;
    if !valid_index(path, assembly)? {
        let expected = chr1_length(assembly).parse::<u64>().ok();
        require(
            fasta_chr1(path)? == expected,
            format!(
                "Reference does not identify the required {assembly} chromosome 1 length: \
                 {path}; existing bytes retained"
            ),
        )?;
    }
    Ok(())
}

/// The words of the first line that starts with "chain ", or none.
fn chain_header<R: Read>(mut lines: TextLines<R>) -> Result<Vec<String>> {
    while let Some(line) = lines.readline()? {
        if line.starts_with("chain ") {
            let words = line.split(is_python_space).filter(|word| !word.is_empty());
            return Ok(words.map(str::to_owned).collect());
        }
    }
    Ok(Vec::new())
}

/// `validate_chain(path)`: UCSC's forward hg19-to-hg38 chromosome 1 chain.
pub(super) fn validate_chain(path: &str) -> Result<()> {
    let header = if compressed(path)? {
        chain_header(TextLines::new(PythonGzipReader::new(open(path)?), path))?
    } else {
        chain_header(TextLines::new(open(path)?, path))?
    };
    let unprefixed = |index: usize| {
        let word = header[index].as_str();
        word.strip_prefix("chr").unwrap_or(word)
    };
    require(
        header.len() == 13
            && unprefixed(2) == "1"
            && unprefixed(7) == "1"
            && header[3] == chr1_length("hg19")
            && header[8] == chr1_length("hg38"),
        format!("Expected the forward UCSC hg19-to-hg38 chain: {path}; existing bytes retained"),
    )
}

/// `receipt_path(path)`
pub(super) fn receipt_path(path: &str) -> Result<String> {
    checked_path(&format!("{path}.source.json"))
}

/// `type(value).__name__` of a decoded JSON value.
fn json_type(value: &Json) -> &'static str {
    match value {
        Json::Null => "NoneType",
        Json::Bool(_) => "bool",
        Json::Int(_) => "int",
        Json::Float(_) => "float",
        Json::Str(_) => "str",
        Json::Array(_) => "list",
        Json::Object(_) => "dict",
    }
}

fn lower_hex(byte: u8) -> bool {
    matches!(byte, b'0'..=b'9' | b'a'..=b'f')
}

/// `re.fullmatch(r"[0-9a-f]{64}", text)`
fn hex64(text: &str) -> bool {
    text.len() == 64 && text.bytes().all(lower_hex)
}

/// `verify_managed(resource)`: the SHA-256 of a managed default, proved by its
/// receipt or, without one, by the exact pinned compressed source beside it.
pub(super) fn verify_managed(resource: &Resource) -> Result<String> {
    let path = resource.path.as_str();
    let source = source(resource.assembly);
    let receipt = receipt_path(path)?;
    if exists(&receipt)? {
        let metadata = json_text::parse(&read_text(&receipt)?, false).map_err(|message| {
            let category = if json_text::is_refusal(&message) {
                "NotImplementedError"
            } else {
                "JSONDecodeError"
            };
            error(category, message)
        })?;
        if !matches!(metadata, Json::Object(_)) {
            return Err(error(
                "AttributeError",
                format!("'{}' object has no attribute 'get'", json_type(&metadata)),
            ));
        }
        let named = |key: &str, expected: &Json| {
            metadata
                .get(key)
                .is_some_and(|value| py_equal(value, expected))
        };
        require(
            named("source", &source.json()) && named("output", &Json::Str(path.to_owned())),
            format!("Pinned UCSC source identity mismatch: {receipt}; existing bytes retained"),
        )?;
        let sha256 = match metadata.get("sha256") {
            Some(Json::Str(text)) => text.clone(),
            _ => String::new(),
        };
        require(
            hex64(&sha256),
            format!("Invalid UCSC source receipt: {receipt}"),
        )?;
        let current = identity_json(&identity(path)?);
        if !named("identity", &current) {
            require(
                sha256 == digest(path, Algorithm::Sha256)?,
                format!("Pinned UCSC reference checksum mismatch: {path}; existing bytes retained"),
            )?;
        }
        return Ok(sha256);
    }
    // Recover a publication interrupted before its receipt, or an earlier
    // installation whose exact pinned compressed source is still available.
    let archive = checked_path(&format!("{path}.gz"))?;
    require(
        is_file(&archive)?
            && file_size(&archive)? == source.bytes
            && digest(&archive, Algorithm::Md5)? == source.md5,
        format!(
            "Existing default reference lacks verified pinned UCSC provenance: {path}. Select an \
             original reference at its own supplied path or move this file aside explicitly; it \
             will not be replaced."
        ),
    )?;
    let expected = source_digest(&archive, None)?;
    require(
        digest(path, Algorithm::Sha256)? == expected,
        format!(
            "Existing default differs from its pinned UCSC source: {path}; existing bytes retained"
        ),
    )?;
    Ok(expected)
}

/// `source_digest(archive, output)`: the SHA-256 of the decompressed archive,
/// written to `output` as it is read.
pub(super) fn source_digest(
    archive: &str,
    mut output: Option<(&mut File, &str)>,
) -> Result<String> {
    let before = identity(archive)?;
    let mut reader = PythonGzipReader::new(open(archive)?);
    let mut sha256 = Sha256::new();
    let mut block = vec![0; BLOCK];
    loop {
        let count = read_block(&mut reader, &mut block)
            .map_err(|failure| read_failure(Box::new(failure), archive))?;
        if count == 0 {
            break;
        }
        sha256.update(&block[..count]);
        if let Some((file, name)) = output.as_mut() {
            file.write_all(&block[..count])
                .map_err(|failure| read_failure(Box::new(failure), name))?;
        }
    }
    require(
        identity(archive)? == before,
        format!("UCSC source archive changed during decompression: {archive}"),
    )?;
    Ok(sha256.finish_hex())
}

/// `preflight(resources)`: read-only; checks every supplied input before any
/// lock, directory or download.
pub(crate) fn preflight(resources: &[Resource]) -> Result<()> {
    for resource in resources {
        let path = resource.path.as_str();
        if lexists(path) {
            require(
                is_file(path)? && file_size(path)? > 0,
                format!(
                    "Missing, empty, or non-file reference resource: {path}; existing bytes \
                     retained"
                ),
            )?;
            if resource.managed {
                verify_managed(resource)?;
            }
            match resource.kind {
                Kind::Fasta => validate_fasta(path, resource.assembly)?,
                Kind::Chain => validate_chain(path)?,
            }
        } else {
            require(
                resource.managed,
                format!(
                    "Explicit reference resource is missing ({}): {path}; restore the selected \
                     original input. A UCSC download will not replace it.",
                    resource.key
                ),
            )?;
            if resource.kind == Kind::Fasta {
                let index = checked_path(&format!("{path}.fai"))?;
                require(
                    !lexists(&index),
                    format!(
                        "Orphan FASTA index exists without its reference: {index}; resolve it \
                         explicitly"
                    ),
                )?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder_error::RuntimeError;
    use std::fs;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT: AtomicU64 = AtomicU64::new(0);

    fn folder() -> String {
        let path = std::env::temp_dir().join(format!(
            "priva-reference-verify-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&path).unwrap();
        crate::installer_prerequisites::paths::abspath(&path.display().to_string()).unwrap()
    }

    fn failure(result: Result<()>) -> (&'static str, String) {
        let failure = result.unwrap_err();
        let known = failure.downcast_ref::<RuntimeError>().unwrap();
        (known.category_name(), known.message())
    }

    fn resource(path: &str, assembly: &'static str, managed: bool, kind: Kind) -> Resource {
        Resource {
            path: path.to_owned(),
            assembly,
            key: "ref_genome".to_owned(),
            managed,
            kind,
        }
    }

    #[test]
    fn indexes_follow_the_python_row_rules() {
        let root = folder();
        let fasta = format!("{root}/x.fa");
        fs::write(&fasta, ">chr1\nACGT\n").unwrap();
        let index = format!("{fasta}.fai");
        for (text, valid) in [
            ("chr1\t248956422\t6\t60\t61\n", true),
            ("1\t248956422\t6\t60\t61\r\n", true),
            ("chr1\t\u{ff12}48956422\t6\t\u{b2}\t61\n", true),
            ("chr1\t248956422\t6\t60\n", false),
            ("chr1\t\u{b2}48956422\t6\t60\t61\n", false),
            ("chr1\t248956422\t6\t60\t61\nchr1\t1\t1\t1\t1\n", false),
            ("", false),
        ] {
            fs::write(&index, text).unwrap();
            assert_eq!(valid_index(&fasta, "hg38").unwrap(), valid, "{text:?}");
        }
        fs::write(&index, b"chr1\t248956422\t6\t60\t61\n\xff").unwrap();
        assert!(!valid_index(&fasta, "hg38").unwrap());
        fs::write(&index, "chr1\tlong\t6\t60\t61\n").unwrap();
        let message = valid_index_at(&index, "hg38").unwrap_err().to_string();
        assert!(message.contains("invalid literal for int()"), "{message}");
    }

    #[test]
    fn fasta_lengths_are_counted_only_without_a_valid_index() {
        let root = folder();
        let fasta = format!("{root}/x.fa");
        fs::write(&fasta, ">chr2\nAA\n>1 description\nAC GT \r\nAA\n>chr3\nA\n").unwrap();
        assert_eq!(fasta_chr1(&fasta).unwrap(), Some(7));
        fs::write(&fasta, ">\nA\n").unwrap();
        let message = fasta_chr1(&fasta).unwrap_err().to_string();
        assert_eq!(message, "IndexError: list index out of range");
        fs::write(&fasta, "hello\n").unwrap();
        let (category, message) = failure(validate_fasta(&fasta, "hg19"));
        assert_eq!(category, "ValueError");
        assert!(message.starts_with("Not a FASTA reference"));
        fs::write(&fasta, ">chr1\nACGT\n").unwrap();
        let (_, message) = failure(validate_fasta(&fasta, "hg19"));
        assert!(message.starts_with("Reference does not identify the required hg19"));
    }

    #[test]
    fn only_the_forward_chain_is_accepted() {
        let root = folder();
        let chain = format!("{root}/forward.chain");
        let forward = "chain 1 chr1 249250621 + 0 1 1 248956422 + 0 1 1\n1\n";
        fs::write(&chain, format!("# note\n{forward}")).unwrap();
        validate_chain(&chain).unwrap();
        let reverse = "chain 1 chr1 248956422 + 0 1 chr1 249250621 + 0 1 1\n";
        fs::write(&chain, reverse).unwrap();
        let (_, message) = failure(validate_chain(&chain));
        assert!(message.starts_with("Expected the forward UCSC hg19-to-hg38 chain"));
        fs::write(&chain, "chain 1 chr1 0249250621 + 0 1 1 248956422 + 0 1 1\n").unwrap();
        assert!(validate_chain(&chain).is_err());
    }

    #[test]
    fn preflight_refuses_missing_originals_and_orphan_indexes() {
        let root = folder();
        let original = format!("{root}/original.fa");
        let explicit = resource(&original, "hg38", false, Kind::Fasta);
        let (_, message) = failure(preflight(&[explicit]));
        assert!(message.starts_with("Explicit reference resource is missing (ref_genome)"));
        let managed = resource(&original, "hg38", true, Kind::Fasta);
        preflight(std::slice::from_ref(&managed)).unwrap();
        fs::write(format!("{original}.fai"), "").unwrap();
        let (_, message) = failure(preflight(&[managed]));
        assert!(message.starts_with("Orphan FASTA index exists without its reference"));
        fs::write(&original, "").unwrap();
        let empty = resource(&original, "hg38", false, Kind::Fasta);
        let (_, message) = failure(preflight(&[empty]));
        assert!(message.starts_with("Missing, empty, or non-file reference resource"));
    }

    #[test]
    fn receipts_must_name_the_pinned_source_and_output() {
        let root = folder();
        let path = format!("{root}/ucsc.hg38.fasta");
        fs::write(&path, ">chr1\nA\n").unwrap();
        let managed = resource(&path, "hg38", true, Kind::Fasta);
        let receipt = format!("{path}.source.json");
        let (_, message) = failure(verify_managed(&managed).map(|_| ()));
        assert!(message.starts_with("Existing default reference lacks verified pinned UCSC"));
        fs::write(&receipt, "[]").unwrap();
        let (category, message) = failure(verify_managed(&managed).map(|_| ()));
        assert_eq!(category, "AttributeError");
        assert_eq!(message, "'list' object has no attribute 'get'");
        let sha256 = digest(&path, Algorithm::Sha256).unwrap();
        let mut metadata = Json::Object(vec![
            ("source".to_owned(), source("hg38").json()),
            ("output".to_owned(), Json::Str(path.clone())),
            ("sha256".to_owned(), Json::Str(sha256.clone())),
            ("identity".to_owned(), Json::Array(vec![])),
        ]);
        fs::write(&receipt, json_text::dumps(&metadata, true).unwrap()).unwrap();
        assert_eq!(verify_managed(&managed).unwrap(), sha256);
        metadata.set("sha256", Json::Str("0".repeat(64)));
        fs::write(&receipt, json_text::dumps(&metadata, true).unwrap()).unwrap();
        let (_, message) = failure(verify_managed(&managed).map(|_| ()));
        assert!(message.starts_with("Pinned UCSC reference checksum mismatch"));
        metadata.set("output", Json::Str(format!("{path}x")));
        fs::write(&receipt, json_text::dumps(&metadata, true).unwrap()).unwrap();
        let (_, message) = failure(verify_managed(&managed).map(|_| ()));
        assert!(message.starts_with("Pinned UCSC source identity mismatch"));
    }
}
