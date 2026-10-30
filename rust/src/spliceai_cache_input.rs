//! Rust version of `scripts/spliceai_cache_input.py`.
//!
//! The script inspects the SpliceAI score files an installation is given,
//! registers compatible ones in the configuration, and adopts a completed
//! SpliceAI2 cache into the managed cache folder. Two commands keep the
//! Python script's options:
//! - `install --config C --cache-dir D --snv S --indel I`: reads only the
//!   bounded header of each score file (gzip is found by its magic bytes, not
//!   its name); a native SpliceAI2 table is handed to `install_utils.sh
//!   SpliceAI2_install` (found in the folder named by the environment variable
//!   `PRIVA_SCRIPT_DIR`, because Rust has no `__file__`; no option is added); a
//!   SpliceAI2 VCF must have its completed adapter receipt, a SpliceAI v1 VCF
//!   its tabix index, which must return the file's first record; then the two
//!   `spliceai_*_prescore` values of the YAML configuration are replaced in
//!   place, keeping every comment, quote and line ending around them;
//! - `adopt-prepared --source-dir A --output-dir B --cache-root C`: copies a
//!   completed SpliceAI2 cache, member by member with SHA-256 checks and a
//!   resumable partial copy, writes the relocated receipt and a completion
//!   record, and publishes the folder with one rename.
//! The report is printed as sorted JSON on one line. Exit codes follow the
//! script: 2 for a rejected input (a ValueError or OSError in Python), 1 for
//! an error the script does not handle (a traceback in Python), the helper
//! script's own exit code when `install_utils.sh` fails, 0 on success.
//!
//! The command line is the script's `argparse` parser reproduced (see the
//! command line section below): same options, abbreviations, usage lines, help
//! texts and error messages, with the program name `priva spliceai-cache-input`.
//!
//! Library functions for later conversions: `inspect_header`,
//! `canonical_assembly`, `Header`.
//!
//! Differences from the Python script (checked by `spliceai_cache_parity.py`):
//! - YAML that PyYAML reads but this tool does not (anchors, aliases, tags, block
//!   scalars, multi-line scalars, keys with spaces or other punctuation, nesting
//!   deeper than one block) is refused with exit code 1 and a message saying so,
//!   where the script goes on (or fails in its own way). The real `config.yaml`
//!   files hold none of them. A byte order mark, a leading `---` line, quoted
//!   keys and tabs inside comments and quoted text are read as PyYAML reads them;
//! - for text PyYAML itself refuses, the exit code agrees (2 here, 1 for the
//!   cutoffs script) but PyYAML's multi-line diagnostic is replaced by one line;
//! - `install_utils.sh` is looked up in the folder named by the environment
//!   variable `PRIVA_SCRIPT_DIR` (the script finds it beside itself); without
//!   it a native-table install stops with exit code 1;
//! - for a corrupt deflate stream the script prints zlib's reason ("invalid
//!   distance too far back"); here the reason is generic. A damaged gzip header
//!   or trailer, a missing end and trailing bytes read as Python reads them;
//! - the first record is decoded on its own, while Python decodes 8 KiB chunks,
//!   so a bad byte after the first line is not noticed and the position in a
//!   decode error counts from the start of the record;
//! - the tabix lookup is done by this crate instead of pysam; every index it
//!   cannot read gets htslib's "Invalid argument" line and pysam's message;
//! - the copied files do not get the extended attributes that
//!   `shutil.copystat` copies (permissions and times are copied);
//! - JSON integers beyond 128 bits, lone surrogates in JSON strings and nesting
//!   deeper than 400 levels are refused with exit code 1.
use crate::cli::Result;
use crate::json_text::{self, Json};
use crate::python_io::{self, COMPRESSED_FILE_ENDED, GzipFailure, PythonGzipReader};
use crate::python_text::{is_python_space, python_int, python_strip};
use crate::sha256::{Sha256, sha256_file, sha256_hex};
use crate::tabix_index;
use crate::yaml_config::{
    self, ConfigDocument, ConfigValue, Parsed as YamlParsed, Scalar, YamlError, python_str_repr,
};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufRead, BufReader, Read, Write};
use std::os::unix::fs::{MetadataExt, OpenOptionsExt, PermissionsExt};
use std::path::{Component, Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

const MAX_LINE: usize = 64 * 1024;
const MAX_HEADER: usize = 1024 * 1024;
const MAX_LINES: usize = 4096;
const REVISION: &str = "cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9";
const MANIFEST_NAME: &str = "SpliceAI2_GRCh38.downloads-complete.json";
const RELOCATION_KIND: &str = "priva_spliceai2_prepared_relocation_v1";
const COPY_BLOCK: usize = 8 * 1024 * 1024;
const CACHE_KEYS: [&str; 2] = ["spliceai_snv_prescore", "spliceai_indel_prescore"];
const CHROMS: [&str; 24] = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
    "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
    "chr22", "chrX", "chrY",
];
const VCF_COLUMNS: [&str; 8] = [
    "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
];

/// Why a command stopped.
#[derive(Debug)]
pub(crate) enum Failure {
    /// A ValueError or OSError in Python: printed as "rejected", exit code 2.
    Rejected(String),
    /// An exception the Python script does not catch: exit code 1.
    Crash(String),
    /// The helper script failed with this exit code.
    Helper(i32),
}

pub(crate) type Done<T> = std::result::Result<T, Failure>;

fn os_error_text(error: &io::Error) -> String {
    let text = error.to_string();
    let text = match text.rfind(" (os error ") {
        Some(at) => text[..at].to_owned(),
        None => text,
    };
    match error.raw_os_error() {
        Some(code) => format!("[Errno {code}] {text}"),
        None => text,
    }
}

/// `bytes.decode("utf-8")` with the text of Python's `UnicodeDecodeError` as the error.
pub(crate) fn python_utf8_decode(bytes: &[u8]) -> std::result::Result<String, String> {
    python_io::decode_utf8(bytes).map(str::to_owned)
}

impl Failure {
    /// An operating-system error about a path, worded as Python words it.
    pub(crate) fn os(error: &io::Error, path: &Path) -> Failure {
        Failure::Rejected(format!(
            "{}: {}",
            os_error_text(error),
            python_str_repr(&path.to_string_lossy())
        ))
    }
}

impl From<io::Error> for Failure {
    fn from(error: io::Error) -> Failure {
        Failure::Rejected(os_error_text(&error))
    }
}

fn reject<T>(message: impl Into<String>) -> Done<T> {
    Err(Failure::Rejected(message.into()))
}

fn crash<T>(message: impl Into<String>) -> Done<T> {
    Err(Failure::Crash(message.into()))
}

fn require(condition: bool, message: &str) -> Done<()> {
    if condition { Ok(()) } else { reject(message) }
}

// ---------------------------------------------------------------- paths

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct Identity {
    device: u64,
    inode: u64,
    size: u64,
    modified_ns: i128,
    changed_ns: i128,
}

fn identity_of(metadata: &fs::Metadata) -> Identity {
    Identity {
        device: metadata.dev(),
        inode: metadata.ino(),
        size: metadata.size(),
        modified_ns: metadata.mtime() as i128 * 1_000_000_000 + metadata.mtime_nsec() as i128,
        changed_ns: metadata.ctime() as i128 * 1_000_000_000 + metadata.ctime_nsec() as i128,
    }
}

fn identity(path: &Path) -> Done<Identity> {
    let metadata = fs::metadata(path).map_err(|error| Failure::os(&error, path))?;
    Ok(identity_of(&metadata))
}

fn identity_json(value: &Identity) -> Json {
    Json::Array(vec![
        Json::Int(value.device as i128),
        Json::Int(value.inode as i128),
        Json::Int(value.size as i128),
        Json::Int(value.modified_ns),
        Json::Int(value.changed_ns),
    ])
}

/// `Path.absolute()`: the working directory is put in front of a relative
/// path, repeated separators and `.` are dropped, `..` is kept.
fn absolute_path(path: &Path) -> PathBuf {
    let joined = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir().unwrap_or_default().join(path)
    };
    joined.components().collect()
}

/// `os.path.normpath` for an absolute path: `..` and `.` are resolved by text alone.
fn lexical_normalize(path: &Path) -> PathBuf {
    let mut parts: Vec<std::ffi::OsString> = Vec::new();
    for component in path.components() {
        match component {
            Component::Normal(name) => parts.push(name.to_owned()),
            Component::ParentDir => {
                parts.pop();
            }
            _ => {}
        }
    }
    let mut result = PathBuf::from("/");
    for part in parts {
        result.push(part);
    }
    result
}

/// `Path.resolve()` without `strict`: links are followed as far as the path exists.
fn resolve_non_strict(path: &Path) -> PathBuf {
    let mut existing = path.to_path_buf();
    let mut tail: Vec<std::ffi::OsString> = Vec::new();
    loop {
        if let Ok(real) = fs::canonicalize(&existing) {
            let mut result = real;
            for part in tail.iter().rev() {
                result.push(part);
            }
            return result;
        }
        let name = existing.file_name().map(|name| name.to_owned());
        let parent = existing.parent().map(|parent| parent.to_path_buf());
        match (name, parent) {
            (Some(name), Some(parent)) => {
                tail.push(name);
                existing = parent;
            }
            _ => return path.to_path_buf(),
        }
    }
}

fn is_unfinished(part: &str) -> bool {
    let lower = part.to_lowercase();
    lower.contains(".partial") || lower == ".building" || lower.starts_with(".building-")
}

fn has_unfinished_part(path: &Path) -> bool {
    path.components()
        .any(|part| is_unfinished(&part.as_os_str().to_string_lossy()))
}

/// `Path(text).name`
fn path_name(text: &str) -> String {
    match Path::new(text).components().next_back() {
        Some(Component::Normal(name)) => name.to_string_lossy().into_owned(),
        Some(Component::ParentDir) => "..".to_owned(),
        _ => String::new(),
    }
}

/// `directory / name`, which is the directory itself for an empty name.
fn join_name(directory: &Path, name: &str) -> PathBuf {
    if name.is_empty() {
        directory.to_path_buf()
    } else {
        directory.join(name)
    }
}

fn is_symlink(path: &Path) -> bool {
    fs::symlink_metadata(path)
        .map(|metadata| metadata.file_type().is_symlink())
        .unwrap_or(false)
}

fn path_text(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

fn file_name_text(path: &Path) -> String {
    path.file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_default()
}

fn is_lower_hex_64(text: &str) -> bool {
    text.len() == 64
        && text
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

/// `Path.resolve(strict=True)`: an error names the first component that cannot be examined.
fn resolve_strict(path: &Path) -> Done<PathBuf> {
    python_io::resolve_strict(path).map_err(|(error, named)| Failure::os(&error, &named))
}

/// The resolved path of a complete, non-empty regular file; never a staging path.
fn complete_path(path: &Path) -> Done<PathBuf> {
    let absolute = absolute_path(path);
    require(
        !has_unfinished_part(&absolute),
        &format!("Incomplete/staging input path: {}", absolute.display()),
    )?;
    let resolved = resolve_strict(&absolute)?;
    require(
        !has_unfinished_part(&resolved),
        &format!("Incomplete/staging input target: {}", absolute.display()),
    )?;
    let metadata = fs::metadata(&resolved).map_err(|error| Failure::os(&error, &resolved))?;
    require(
        metadata.is_file() && metadata.len() > 0,
        &format!(
            "Input is not a nonempty regular file: {}",
            absolute.display()
        ),
    )?;
    Ok(resolved)
}

// ---------------------------------------------------------------- headers

/// What a score file's header says it is.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Header {
    /// `spliceai2_tsv`, `vcf_spliceai2` or `vcf_v1`.
    pub(crate) format: &'static str,
    pub(crate) assembly: Option<&'static str>,
}

fn native_columns() -> Vec<String> {
    let mut columns: Vec<String> = ["chrom", "pos", "ref", "alt", "strand"]
        .iter()
        .map(|name| name.to_string())
        .collect();
    for effect in [
        "donor_gain",
        "donor_loss",
        "acceptor_gain",
        "acceptor_loss",
        "jxn_gain",
        "jxn_loss",
    ] {
        let mut fields = vec!["delta_score", "ref_score", "alt_score"];
        if effect.starts_with("jxn") {
            fields.extend(["donor_dist", "acceptor_dist"]);
        } else {
            fields.push("dist");
        }
        for field in &fields {
            for index in 0..3 {
                columns.push(format!("{effect}_{field}_{index}"));
            }
        }
    }
    columns.push("spliceai2_summary_score".to_owned());
    columns
}

/// The assembly name behind an alias, as the script accepts them.
pub(crate) fn canonical_assembly(assembly: &str) -> Done<&'static str> {
    match assembly.to_lowercase().as_str() {
        "hg19" | "grch37" => Ok("GRCh37"),
        "hg38" | "grch38" => Ok("GRCh38"),
        _ => reject(format!(
            "Unsupported assembly: {}",
            python_str_repr(assembly)
        )),
    }
}

/// `canonical_assembly` for a value read from a YAML configuration.
fn canonical_assembly_of_config(value: Option<&ConfigValue>) -> Done<&'static str> {
    match value {
        Some(ConfigValue::Scalar(Scalar::Str(text))) => canonical_assembly(text),
        Some(ConfigValue::Scalar(other)) => {
            reject(format!("Unsupported assembly: {}", other.python_repr()))
        }
        Some(ConfigValue::Opaque(text)) => reject(format!("Unsupported assembly: {text}")),
        None => reject("Unsupported assembly: None"),
    }
}

/// The assemblies named in a text: whole words `GRCh37`, `GRCh38`, `hg19`, `hg38`.
fn assemblies_in(value: &str) -> Vec<&'static str> {
    let bytes = value.as_bytes();
    let lower: Vec<u8> = bytes.iter().map(|byte| byte.to_ascii_lowercase()).collect();
    let mut found: Vec<&'static str> = Vec::new();
    for (token, canonical) in [
        ("grch37", "GRCh37"),
        ("grch38", "GRCh38"),
        ("hg19", "GRCh37"),
        ("hg38", "GRCh38"),
    ] {
        let wanted = token.as_bytes();
        let mut start = 0;
        while start + wanted.len() <= lower.len() {
            if &lower[start..start + wanted.len()] == wanted {
                let before_ok = start == 0 || !bytes[start - 1].is_ascii_alphanumeric();
                let after = start + wanted.len();
                let after_ok = after == bytes.len() || !bytes[after].is_ascii_alphanumeric();
                if before_ok && after_ok && !found.contains(&canonical) {
                    found.push(canonical);
                }
            }
            start += 1;
        }
    }
    found
}

/// The values of `assembly=` in a `##contig=<...>` line's inside, per the
/// pattern `(?:^|,)assembly=("[^"\n]*"|[^,>]+)` read case-insensitively.
fn contig_assembly_values(inner: &str) -> Vec<String> {
    let bytes = inner.as_bytes();
    let lower: Vec<u8> = bytes.iter().map(|byte| byte.to_ascii_lowercase()).collect();
    let key = b"assembly=";
    let mut values = Vec::new();
    let mut at = 0;
    while at < bytes.len() {
        let key_start = if at == 0 && lower[at..].starts_with(key) {
            Some(at)
        } else if bytes[at] == b',' && lower[at + 1..].starts_with(key) {
            Some(at + 1)
        } else {
            None
        };
        let Some(key_start) = key_start else {
            at += 1;
            continue;
        };
        let group = key_start + key.len();
        let mut end = None;
        if bytes.get(group) == Some(&b'"') {
            let mut scan = group + 1;
            while scan < bytes.len() && bytes[scan] != b'"' && bytes[scan] != b'\n' {
                scan += 1;
            }
            if scan < bytes.len() && bytes[scan] == b'"' {
                end = Some(scan + 1);
            }
        }
        if end.is_none() {
            let mut scan = group;
            while scan < bytes.len() && bytes[scan] != b',' && bytes[scan] != b'>' {
                scan += 1;
            }
            if scan > group {
                end = Some(scan);
            }
        }
        match end {
            Some(end) => {
                values.push(inner[group..end].to_owned());
                at = end;
            }
            None => at += 1,
        }
    }
    values
}

/// `re.search(r"spliceai(?:2|[ _-]+v?2)(?:\b|_)", value, re.I)`
fn mentions_spliceai2(value: &str) -> bool {
    let chars: Vec<char> = value.chars().collect();
    let needle: Vec<char> = "spliceai".chars().collect();
    let boundary_after = |index: usize| index >= chars.len() || !chars[index].is_alphanumeric();
    for start in 0..chars.len() {
        if start + needle.len() > chars.len() {
            break;
        }
        let matches_name = (0..needle.len())
            .all(|offset| chars[start + offset].to_ascii_lowercase() == needle[offset]);
        if !matches_name {
            continue;
        }
        let at = start + needle.len();
        if at < chars.len() && chars[at] == '2' && boundary_after(at + 1) {
            return true;
        }
        let mut run = at;
        while run < chars.len() && matches!(chars[run], ' ' | '_' | '-') {
            run += 1;
        }
        if run > at {
            let mut next = run;
            if next < chars.len() && chars[next].to_ascii_lowercase() == 'v' {
                next += 1;
            }
            if next < chars.len() && chars[next] == '2' && boundary_after(next + 1) {
                return true;
            }
        }
    }
    false
}

/// `^2(\.\d+)*$` when `first_digits` is "2", `^[12](\.\d+)*$` when it is "12".
fn is_version_text(value: &str, first_digits: &str) -> bool {
    let mut parts = value.split('.');
    let head = parts.next().unwrap_or("");
    if head.chars().count() != 1 || !first_digits.contains(head) {
        return false;
    }
    parts.all(|part| !part.is_empty() && part.bytes().all(|byte| byte.is_ascii_digit()))
}

/// `_metadata`: the key and value of a `##key=value` line, with structured
/// values (`<...>`) checked for closed brackets and quotes.
fn metadata(line: &str) -> Done<(String, String)> {
    let body = line.strip_prefix("##").unwrap_or(line);
    let key_length = body
        .bytes()
        .take_while(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'.' | b'-'))
        .count();
    let valid = key_length > 0
        && body.as_bytes().get(key_length) == Some(&b'=')
        && body.len() > key_length + 1;
    require(valid, "Invalid VCF metadata line")?;
    let key = body[..key_length].to_owned();
    let value = body[key_length + 1..].to_owned();
    if value.starts_with('<') {
        require(value.ends_with('>'), "Truncated structured VCF metadata")?;
        let characters: Vec<char> = value.chars().collect();
        let inside = if characters.len() >= 2 {
            &characters[1..characters.len() - 1]
        } else {
            &characters[0..0]
        };
        let (mut quoted, mut escaped) = (false, false);
        for character in inside {
            if escaped {
                escaped = false;
            } else if *character == '\\' && quoted {
                escaped = true;
            } else if *character == '"' {
                quoted = !quoted;
            } else if (*character == '<' || *character == '>') && !quoted {
                return reject("Invalid structured VCF metadata brackets");
            }
        }
        require(!quoted && !escaped, "Unclosed VCF metadata quotation")?;
    }
    Ok((key, value))
}

struct Limits {
    line: usize,
    header: usize,
    lines: usize,
}

const LIMITS: Limits = Limits {
    line: MAX_LINE,
    header: MAX_HEADER,
    lines: MAX_LINES,
};

/// `readline(limit)` on bytes: up to `limit` bytes, through the first newline.
fn read_bounded_line(reader: &mut dyn BufRead, limit: usize) -> io::Result<Vec<u8>> {
    let mut line: Vec<u8> = Vec::new();
    loop {
        let (used, done) = {
            let available = reader.fill_buf()?;
            if available.is_empty() {
                (0, true)
            } else {
                let room = limit - line.len();
                let window = available.len().min(room);
                match available[..window].iter().position(|byte| *byte == b'\n') {
                    Some(found) => {
                        line.extend_from_slice(&available[..=found]);
                        (found + 1, true)
                    }
                    None => {
                        line.extend_from_slice(&available[..window]);
                        (window, line.len() >= limit)
                    }
                }
            }
        };
        reader.consume(used);
        if done {
            break;
        }
    }
    Ok(line)
}

fn starts_with_gzip_magic(path: &Path) -> Done<bool> {
    let mut file = File::open(path).map_err(|error| Failure::os(&error, path))?;
    let mut magic = [0u8; 2];
    let mut filled = 0;
    while filled < 2 {
        let read = file
            .read(&mut magic[filled..])
            .map_err(|error| Failure::os(&error, path))?;
        if read == 0 {
            break;
        }
        filled += read;
    }
    Ok(filled == 2 && magic == [0x1f, 0x8b])
}

/// Reads the bounded header of a score file and says what it is.
pub(crate) fn inspect_header(path: &Path) -> Done<Header> {
    inspect_header_limited(path, &LIMITS)
}

fn inspect_header_limited(path: &Path, limits: &Limits) -> Done<Header> {
    let resolved = complete_path(path)?;
    let before = identity(&resolved)?;
    let compressed = starts_with_gzip_magic(&resolved)?;
    let file = File::open(&resolved).map_err(|error| Failure::os(&error, &resolved))?;
    let mut stream: Box<dyn BufRead> = if compressed {
        Box::new(BufReader::new(PythonGzipReader::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let corrupt = |detail: String| -> Failure {
        Failure::Rejected(format!(
            "Corrupt/truncated cache header: {}: {detail}",
            resolved.display()
        ))
    };
    // The script reports EOFError, UnicodeError and zlib.error as a corrupt
    // header; BadGzipFile and operating-system errors are OSErrors it lets through.
    let read_failure = |error: io::Error| -> Failure {
        match GzipFailure::of(&error) {
            Some(GzipFailure::EndOfFile) => corrupt(COMPRESSED_FILE_ENDED.to_owned()),
            Some(GzipFailure::Zlib(text)) => corrupt(text.clone()),
            Some(GzipFailure::Bad(text)) => Failure::Rejected(text.clone()),
            None => Failure::Rejected(os_error_text(&error)),
        }
    };
    let native_header = native_columns().join("\t");
    let mut assemblies: Vec<&'static str> = Vec::new();
    let mut versions: Vec<u8> = Vec::new();
    let mut total = 0usize;
    let mut result: Option<Header> = None;
    for number in 1..=limits.lines {
        let raw = read_bounded_line(stream.as_mut(), limits.line + 1).map_err(&read_failure)?;
        require(!raw.is_empty(), "Missing complete VCF/SpliceAI2 header")?;
        total += raw.len();
        require(
            raw.len() <= limits.line && total <= limits.header,
            "Cache header exceeds bounded size limit",
        )?;
        require(
            raw.last() == Some(&b'\n'),
            "Truncated cache header line (missing newline)",
        )?;
        let mut bytes = &raw[..raw.len() - 1];
        if bytes.last() == Some(&b'\r') {
            bytes = &bytes[..bytes.len() - 1];
        }
        let line = python_io::decode_utf8(bytes).map_err(&corrupt)?;
        require(
            !line
                .chars()
                .any(|c| ((c as u32) < 32 && c != '\t') || c as u32 == 127),
            "Control byte in cache header",
        )?;
        if number == 1 {
            if line == native_header {
                result = Some(Header {
                    format: "spliceai2_tsv",
                    assembly: Some("GRCh38"),
                });
                break;
            }
            let supported = line.len() == 20
                && line.starts_with("##fileformat=VCFv4.")
                && matches!(line.as_bytes()[19], b'0'..=b'5');
            require(
                supported,
                "Unsupported cache header: expected VCF or exact 84-column SpliceAI2 TSV",
            )?;
            continue;
        }
        if line.starts_with("##") {
            let (key, value) = metadata(line)?;
            let lower = key.to_lowercase();
            require(lower != "fileformat", "Duplicate VCF fileformat line")?;
            if lower == "reference" || lower == "assembly" {
                for assembly in assemblies_in(&value) {
                    if !assemblies.contains(&assembly) {
                        assemblies.push(assembly);
                    }
                }
            } else if lower == "contig" {
                let characters: Vec<char> = value.chars().collect();
                let inner: String = if characters.len() >= 2 {
                    characters[1..characters.len() - 1].iter().collect()
                } else {
                    String::new()
                };
                for group in contig_assembly_values(&inner) {
                    for assembly in assemblies_in(&group) {
                        if !assemblies.contains(&assembly) {
                            assemblies.push(assembly);
                        }
                    }
                }
            }
            if (lower.starts_with("spliceai2") || (lower == "source" && mentions_spliceai2(&value)))
                && !versions.contains(&2)
            {
                versions.push(2);
            }
            if lower == "spliceai2version" {
                require(
                    is_version_text(&value, "2"),
                    "Unsupported explicit SpliceAI2 version",
                )?;
            }
            if lower == "spliceaiversion" || lower == "spliceai_version" {
                require(
                    is_version_text(&value, "12"),
                    "Unsupported explicit SpliceAI version",
                )?;
                let version = if value.starts_with('1') { 1 } else { 2 };
                if !versions.contains(&version) {
                    versions.push(version);
                }
            }
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        require(
            columns.len() >= 8 && columns[..8] == VCF_COLUMNS,
            "Missing valid complete #CHROM VCF header",
        )?;
        let samples = &columns[columns.len().min(9)..];
        let unique_samples: HashSet<&&str> = samples.iter().collect();
        let sample_columns_ok = columns.len() == 8
            || (columns.len() >= 10
                && columns[8] == "FORMAT"
                && samples
                    .iter()
                    .all(|name| !name.is_empty() && !name.chars().any(is_python_space))
                && unique_samples.len() == samples.len());
        require(sample_columns_ok, "Invalid VCF sample columns")?;
        require(
            assemblies.len() <= 1,
            "Conflicting assembly metadata in VCF header",
        )?;
        require(versions.len() <= 1, "Conflicting SpliceAI version metadata")?;
        let version_two = versions.contains(&2);
        require(
            !version_two || assemblies != ["GRCh37"],
            "SpliceAI2 VCF conflicts with GRCh37 assembly",
        )?;
        result = Some(Header {
            format: if version_two {
                "vcf_spliceai2"
            } else {
                "vcf_v1"
            },
            assembly: assemblies.first().copied().or(if version_two {
                Some("GRCh38")
            } else {
                None
            }),
        });
        break;
    }
    let header = match result {
        Some(header) => header,
        None => return reject("Cache header exceeds bounded line limit"),
    };
    // Also finds a truncated gzip footer for header-only files; it checks no record.
    let mut one = [0u8; 1];
    stream.read(&mut one).map_err(&read_failure)?;
    let again = resolve_strict(path)?;
    require(
        identity(&resolved)? == before && again == resolved,
        &format!("Input changed during header inspection: {}", path.display()),
    )?;
    Ok(header)
}

struct PairReport {
    format: &'static str,
    assembly: &'static str,
    snv: Header,
    indel: Header,
}

/// Inspects both inputs and enforces matching format and assembly.
fn validate_pair(assembly: &'static str, snv: &Path, indel: &Path) -> Done<PairReport> {
    let snv_header = inspect_header(snv)?;
    let indel_header = inspect_header(indel)?;
    require(
        snv_header.format == indel_header.format,
        "Mixed SpliceAI cache formats/versions are not supported",
    )?;
    for (kind, header) in [("snv", &snv_header), ("indel", &indel_header)] {
        require(
            header.assembly.is_none() || header.assembly == Some(assembly),
            &format!("{kind} cache assembly conflicts with {assembly}"),
        )?;
    }
    Ok(PairReport {
        format: snv_header.format,
        assembly,
        snv: snv_header,
        indel: indel_header,
    })
}

fn header_json(header: &Header) -> Json {
    Json::Object(vec![
        ("format".to_owned(), Json::Str(header.format.to_owned())),
        (
            "assembly".to_owned(),
            match header.assembly {
                Some(name) => Json::Str(name.to_owned()),
                None => Json::Null,
            },
        ),
    ])
}

fn report_json(report: &PairReport) -> Json {
    Json::Object(vec![
        ("format".to_owned(), Json::Str(report.format.to_owned())),
        ("assembly".to_owned(), Json::Str(report.assembly.to_owned())),
        ("snv".to_owned(), header_json(&report.snv)),
        ("indel".to_owned(), header_json(&report.indel)),
    ])
}

// ---------------------------------------------------------------- JSON files

/// A `JSONDecodeError` is a rejection (exit 2); what this tool refuses to read is not.
fn json_error(error: String) -> Failure {
    if json_text::is_refusal(&error) {
        Failure::Crash(error)
    } else {
        Failure::Rejected(error)
    }
}

/// `json.loads(path.read_text())`: the file is decoded as UTF-8 (a byte order
/// mark is an error) and every line break is read as `\n`.
fn read_json(path: &Path, reject_duplicates: bool) -> Done<Json> {
    let bytes = fs::read(path).map_err(|error| Failure::os(&error, path))?;
    json_from_file_bytes(&bytes, reject_duplicates)
}

/// `json.loads(text)` for the bytes of a file read with `read_text()`.
fn json_from_file_bytes(bytes: &[u8], reject_duplicates: bool) -> Done<Json> {
    let text = python_io::decode_utf8(bytes)
        .map_err(Failure::Rejected)?
        .replace("\r\n", "\n")
        .replace('\r', "\n");
    json_text::parse(&text, reject_duplicates).map_err(json_error)
}

/// `json.loads(data)` for bytes: Python decodes them as UTF-8 and drops a leading byte order mark.
fn json_from_bytes(bytes: &[u8], reject_duplicates: bool) -> Done<Json> {
    let body = bytes.strip_prefix(b"\xef\xbb\xbf").unwrap_or(bytes);
    let text = python_io::decode_utf8(body).map_err(Failure::Rejected)?;
    json_text::parse(text, reject_duplicates).map_err(json_error)
}

/// `type(value).__name__` of the Python value a JSON document gives.
fn json_type_name(value: &Json) -> &'static str {
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

/// `value.get(key)` where a missing key gives `None`.
fn dict_get<'a>(value: &'a Json, key: &str) -> Done<Option<&'a Json>> {
    match value {
        Json::Object(_) => Ok(value.get(key)),
        other => crash(format!(
            "AttributeError: '{}' object has no attribute 'get'",
            json_type_name(other)
        )),
    }
}

/// `value[key]`
fn dict_index<'a>(value: &'a Json, key: &str) -> Done<&'a Json> {
    match value {
        Json::Object(_) => match value.get(key) {
            Some(found) => Ok(found),
            None => crash(format!("KeyError: {}", python_str_repr(key))),
        },
        Json::Array(_) => crash("TypeError: list indices must be integers or slices, not str"),
        Json::Str(_) => crash("TypeError: string indices must be integers, not 'str'"),
        other => crash(format!(
            "TypeError: '{}' object is not subscriptable",
            json_type_name(other)
        )),
    }
}

fn is_string(value: Option<&Json>) -> bool {
    matches!(value, Some(Json::Str(_)))
}

/// The text of a value the script hands to `Path(...)`.
fn string_of(value: &Json) -> Done<&str> {
    match value {
        Json::Str(text) => Ok(text),
        other => crash(format!(
            "TypeError: argument should be a str or an os.PathLike object where __fspath__ returns a str, not '{}'",
            json_type_name(other)
        )),
    }
}

fn dumps_indented_line(value: &Json) -> Done<Vec<u8>> {
    let mut text = json_text::dumps_indented(value).map_err(json_error)?;
    text.push('\n');
    Ok(text.into_bytes())
}

// ---------------------------------------------------------------- manifests and receipts

fn is_published(kind: &str, chrom: &str) -> bool {
    (kind == "snv" || kind == "indel")
        && CHROMS.contains(&chrom)
        && !(kind == "indel" && chrom == "chrY")
}

/// The adjacent downloads-complete manifest of two native tables, after a cheap check of all 47 shards.
fn discover_manifest(snv: &Path, indel: &Path) -> Done<PathBuf> {
    let selected = [
        ("snv", complete_path(snv)?),
        ("indel", complete_path(indel)?),
    ];
    let mut candidates: Vec<PathBuf> = Vec::new();
    for (_, path) in &selected {
        let candidate = path.parent().unwrap_or(Path::new("/")).join(MANIFEST_NAME);
        if !candidates.contains(&candidate) {
            candidates.push(candidate);
        }
    }
    require(
        candidates.len() == 1,
        "Selected native tables must share an adjacent downloads-complete manifest",
    )?;
    let path = complete_path(&candidates[0])?;
    let before = identity(&path)?;
    require(
        before.size <= MAX_HEADER as u64,
        "Source manifest exceeds bounded size limit",
    )?;
    let manifest = read_json(&path, true)?;
    require(
        matches!(manifest, Json::Object(_)),
        "Source manifest must be an object",
    )?;
    let confirmed = manifest
        .get("confirmed_by")
        .and_then(Json::as_str)
        .map(|text| !python_strip(text).is_empty())
        .unwrap_or(false);
    require(
        matches!(manifest.get("schema_version"), Some(Json::Int(1)))
            && manifest.get("status").and_then(Json::as_str) == Some("downloads_complete")
            && confirmed,
        "A confirmed downloads_complete manifest is required",
    )?;
    require(
        manifest.get("resource").and_then(Json::as_str) == Some("SpliceAI2")
            && manifest.get("version").and_then(Json::as_str) == Some("2.0")
            && manifest.get("assembly").and_then(Json::as_str) == Some("GRCh38")
            && manifest.get("source_revision").and_then(Json::as_str) == Some(REVISION),
        "Wrong SpliceAI2 source version/assembly/revision",
    )?;
    let files = match manifest.get("files") {
        Some(Json::Array(files)) if files.len() == 47 => files,
        _ => return reject("Full manifest requires all 47 published shards"),
    };
    let mut roles: HashSet<(String, String)> = HashSet::new();
    let mut paths: HashMap<PathBuf, (String, String)> = HashMap::new();
    for item in files {
        require(
            matches!(item, Json::Object(_)),
            "Invalid source shard entry",
        )?;
        let (kind, chrom) = match (
            item.get("kind").and_then(Json::as_str),
            item.get("chrom").and_then(Json::as_str),
        ) {
            (Some(kind), Some(chrom)) => (kind.to_owned(), chrom.to_owned()),
            _ => return reject("Invalid source shard role"),
        };
        let role = (kind.clone(), chrom.clone());
        require(
            is_published(&kind, &chrom) && !roles.contains(&role),
            "Duplicate/unsupported source shard role",
        )?;
        roles.insert(role.clone());
        let name = match item.get("path").and_then(Json::as_str) {
            Some(name) if Path::new(name).is_absolute() => name,
            _ => return reject("Manifest source paths must be absolute"),
        };
        require(
            path_name(name) == format!("SpliceAI2_GRCh38_{kind}_{chrom}.tsv.gz"),
            "Source filename conflicts with manifest role",
        )?;
        let target = complete_path(Path::new(name))?;
        require(
            !paths.contains_key(&target),
            "One source path cannot represent multiple shard roles",
        )?;
        let size_matches = match item.get("bytes") {
            Some(Json::Int(bytes)) => {
                *bytes > 0
                    && fs::metadata(&target)
                        .map_err(|error| Failure::os(&error, &target))?
                        .len() as i128
                        == *bytes
            }
            _ => false,
        };
        require(
            size_matches,
            "Source shard size does not match complete manifest",
        )?;
        let pinned = item
            .get("sha256")
            .and_then(Json::as_str)
            .map(is_lower_hex_64)
            .unwrap_or(false);
        require(pinned, "Source shard needs a pinned SHA256")?;
        paths.insert(target, role);
    }
    require(
        roles.len() == 47,
        "Full manifest requires all 47 published shards",
    )?;
    for (kind, target) in &selected {
        let listed = paths
            .get(target)
            .map(|role| role.0 == *kind)
            .unwrap_or(false);
        require(
            listed,
            &format!("Selected {kind} input is absent from its manifest role"),
        )?;
    }
    require(
        identity(&path)? == before,
        "Source manifest changed during inspection",
    )?;
    Ok(path)
}

/// An input and what it looked like when it was validated.
#[derive(Clone)]
struct Tracked {
    requested: String,
    resolved: PathBuf,
    before: Identity,
}

fn unchanged_inputs(inputs: &[Tracked]) -> Done<()> {
    for item in inputs {
        let requested = Path::new(&item.requested);
        let again = resolve_strict(requested)?;
        let now = identity(&item.resolved)?;
        require(
            again == item.resolved && now == item.before,
            &format!("Cache input changed after validation: {}", item.requested),
        )?;
    }
    Ok(())
}

/// A completed adapter publication, held under the driver's own lock for as long as it lives.
struct Projection {
    _lock: File,
    tracked: Vec<Tracked>,
}

fn has_expected_receipt_fields(data: &Json) -> bool {
    let texts = [
        ("resource", "SpliceAI2"),
        ("version", "2.0"),
        ("assembly", "GRCh38"),
        ("status", "prepared_annotation_only"),
        ("coverage", "complete_published"),
        ("source_revision", REVISION),
    ];
    matches!(data.get("schema_version"), Some(Json::Int(1)))
        && matches!(data.get("scientific_qualified"), Some(Json::Bool(false)))
        && texts
            .iter()
            .all(|(key, wanted)| data.get(key).and_then(Json::as_str) == Some(*wanted))
}

/// Cheap publication check of a SpliceAI2 VCF pair: the completed receipt of
/// its folder, no pending build, outputs that match the receipt.
fn completed_projection(snv: &Path, indel: &Path, require_existing_lock: bool) -> Done<Projection> {
    let selected = [
        ("snv", complete_path(snv)?),
        ("indel", complete_path(indel)?),
    ];
    let folder = selected[0]
        .1
        .parent()
        .unwrap_or(Path::new("/"))
        .to_path_buf();
    let other = selected[1].1.parent().unwrap_or(Path::new("/"));
    require(
        folder == other,
        "SpliceAI2 VCF pair must share one completed adapter receipt",
    )?;
    let receipt_path = folder.join("adapter.ready.json");
    let receipt = complete_path(&receipt_path)?;
    require(
        receipt == receipt_path,
        "Adapter receipt must be a regular adjacent file, not a symlink",
    )?;
    let lock_path = folder.join(".install.lock");
    require(
        !is_symlink(&lock_path),
        "Adapter installation lock must not be a symlink",
    )?;
    if require_existing_lock {
        require(
            fs::metadata(&lock_path)
                .map(|m| m.is_file())
                .unwrap_or(false),
            "Original completed adapter must already contain its installation lock",
        )?;
    }
    let lock = if require_existing_lock {
        File::open(&lock_path)
    } else {
        OpenOptions::new()
            .append(true)
            .create(true)
            .open(&lock_path)
    }
    .map_err(|error| Failure::os(&error, &lock_path))?;
    match lock.try_lock() {
        Ok(()) => {}
        Err(fs::TryLockError::WouldBlock) => {
            return reject("SpliceAI2 preparation is active; refusing VCF registration");
        }
        Err(fs::TryLockError::Error(error)) => return Err(Failure::os(&error, &lock_path)),
    }
    let entries = fs::read_dir(&folder).map_err(|error| Failure::os(&error, &folder))?;
    for entry in entries {
        let entry = entry.map_err(|error| Failure::os(&error, &folder))?;
        let name = entry.file_name().to_string_lossy().into_owned();
        require(
            !(name == ".building" || name.starts_with(".building-")),
            "SpliceAI2 output has a pending .building directory",
        )?;
    }
    let before = identity(&receipt)?;
    require(
        before.size <= MAX_HEADER as u64,
        "Adapter receipt exceeds bounded size limit",
    )?;
    let data = read_json(&receipt, true)?;
    require(
        matches!(data, Json::Object(_)),
        "Adapter receipt must be an object",
    )?;
    require(
        has_expected_receipt_fields(&data),
        "SpliceAI2 VCF requires a completed full GRCh38 adapter receipt",
    )?;
    require(
        data.get("source_manifest_sha256")
            .and_then(Json::as_str)
            .map(is_lower_hex_64)
            .unwrap_or(false),
        "Adapter receipt lacks source manifest identity",
    )?;
    let outputs = match data.get("outputs") {
        Some(found) => found,
        None => return reject("Adapter receipt needs both output roles"),
    };
    let has_both_roles = match outputs {
        Json::Object(items) => {
            items.len() == 2 && outputs.get("snv").is_some() && outputs.get("indel").is_some()
        }
        _ => false,
    };
    require(has_both_roles, "Adapter receipt needs both output roles")?;
    let mut tracked = vec![Tracked {
        requested: path_text(&receipt),
        resolved: receipt.clone(),
        before,
    }];
    let mut seen: Vec<PathBuf> = Vec::new();
    for (role, target) in &selected {
        let item = dict_index(outputs, role)?;
        require(
            matches!(item, Json::Object(_)),
            "Invalid adapter output entry",
        )?;
        for (key, checksum) in [("path", "sha256"), ("index_path", "index_sha256")] {
            let name = match item.get(key).and_then(Json::as_str) {
                Some(name) if Path::new(name).is_absolute() => name,
                _ => return reject("Adapter output paths must be absolute"),
            };
            let path = complete_path(Path::new(name))?;
            require(
                path.parent() == Some(folder.as_path()) && !seen.contains(&path),
                "Adapter outputs must be distinct files beside their receipt",
            )?;
            require(
                item.get(checksum)
                    .and_then(Json::as_str)
                    .map(is_lower_hex_64)
                    .unwrap_or(false),
                "Adapter output lacks recorded SHA256",
            )?;
            if key == "path" {
                require(
                    path == *target,
                    &format!("Selected {role} VCF is not its completed receipt output"),
                )?;
                let size_matches = match item.get("bytes") {
                    Some(Json::Int(bytes)) => {
                        fs::metadata(&path)
                            .map_err(|error| Failure::os(&error, &path))?
                            .len() as i128
                            == *bytes
                    }
                    _ => false,
                };
                require(size_matches, "Adapter VCF size does not match receipt")?;
            } else {
                let vcf_name = item.get("path").and_then(Json::as_str).unwrap_or("");
                require(
                    name == format!("{vcf_name}.tbi"),
                    "Adapter index does not match its VCF",
                )?;
            }
            seen.push(path.clone());
            let before = identity(&path)?;
            tracked.push(Tracked {
                requested: name.to_owned(),
                resolved: path,
                before,
            });
        }
    }
    if data.get("relocation").is_some() {
        tracked.extend(relocation_provenance(&folder, &receipt, &data)?);
    }
    unchanged_inputs(&tracked)?;
    Ok(Projection {
        _lock: lock,
        tracked,
    })
}

/// The receipt of a copied cache: the original's, with every cache location
/// moved into `destination` and the relocation record added.
fn relocated_receipt(original: &Json, destination: &Path, relocation: &Json) -> Done<Json> {
    let mut value = original.clone();
    let outputs = match dict_index(&value, "outputs")? {
        Json::Object(items) => items.clone(),
        other => {
            return crash(format!(
                "AttributeError: '{}' object has no attribute 'values'",
                json_type_name(other)
            ));
        }
    };
    let mut moved: Vec<(String, Json)> = Vec::new();
    for (role, item) in outputs {
        let mut item = item;
        let path = path_name(string_of(dict_index(&item, "path")?)?);
        item.set("path", Json::Str(path_text(&join_name(destination, &path))));
        let index_path = path_name(string_of(dict_index(&item, "index_path")?)?);
        item.set(
            "index_path",
            Json::Str(path_text(&join_name(destination, &index_path))),
        );
        moved.push((role, item));
    }
    value.set("outputs", Json::Object(moved));
    let context = dict_index(&value, "gene_context")?.clone();
    let context_name = path_name(string_of(dict_index(&context, "path")?)?);
    let mut context = context;
    context.set(
        "path",
        Json::Str(path_text(&join_name(destination, &context_name))),
    );
    value.set("gene_context", context);
    value.set("relocation", relocation.clone());
    Ok(value)
}

fn relocation_provenance(folder: &Path, receipt: &Path, data: &Json) -> Done<Vec<Tracked>> {
    let relocation = dict_index(data, "relocation")?;
    require(
        matches!(relocation, Json::Object(_))
            && relocation.get("kind").and_then(Json::as_str) == Some(RELOCATION_KIND),
        "Unknown prepared SpliceAI2 relocation provenance",
    )?;
    let mut tracked: Vec<Tracked> = Vec::new();
    let mut related: HashMap<&str, PathBuf> = HashMap::new();
    for (key, basename) in [
        ("original_receipt", "adapter.original.ready.json"),
        ("record", "relocation.ready.json"),
    ] {
        let expected = path_text(&folder.join(basename));
        require(
            relocation.get(key).and_then(Json::as_str) == Some(expected.as_str()),
            "Relocation metadata must be adjacent to the copied cache",
        )?;
        let path = complete_path(Path::new(&expected))?;
        let size = fs::metadata(&path)
            .map_err(|error| Failure::os(&error, &path))?
            .len();
        require(
            path == folder.join(basename) && size <= MAX_HEADER as u64,
            "Relocation metadata must be a bounded regular adjacent file",
        )?;
        let before = identity(&path)?;
        tracked.push(Tracked {
            requested: path_text(&path),
            resolved: path.clone(),
            before,
        });
        related.insert(key, path);
    }
    let original_path = &related["original_receipt"];
    let original = fs::read(original_path).map_err(|error| Failure::os(&error, original_path))?;
    require(
        relocation
            .get("original_receipt_sha256")
            .and_then(Json::as_str)
            == Some(sha256_hex(&original).as_str()),
        "Copied original adapter receipt identity changed",
    )?;
    let original_data = json_from_bytes(&original, true)?;
    let rebuilt = relocated_receipt(&original_data, folder, relocation)?;
    require(
        json_text::py_equal(&rebuilt, data),
        "Relocated receipt changes more than cache locations and relocation provenance",
    )?;
    let record = read_json(&related["record"], true)?;
    let receipt_sha = sha256_file(receipt).map_err(|error| Failure::os(&error, receipt))?;
    let binds_both = dict_get(&record, "kind")?.and_then(Json::as_str) == Some(RELOCATION_KIND)
        && dict_get(&record, "status")?.and_then(Json::as_str) == Some("complete")
        && matches!(
            dict_get(&record, "scientific_qualified")?,
            Some(Json::Bool(false))
        )
        && dict_get(&record, "original_receipt_sha256")?.and_then(Json::as_str)
            == relocation
                .get("original_receipt_sha256")
                .and_then(Json::as_str)
        && dict_get(&record, "relocated_receipt_sha256")?.and_then(Json::as_str)
            == Some(receipt_sha.as_str());
    require(
        binds_both,
        "Relocation completion record does not bind both adapter receipts",
    )?;
    Ok(tracked)
}

// ---------------------------------------------------------------- copying a cache

static TEMPORARY_COUNTER: AtomicU64 = AtomicU64::new(0);

/// A new file `prefix + random text` in `directory`, readable only by the owner.
fn create_temporary(directory: &Path, prefix: &str) -> Done<(File, PathBuf)> {
    for _ in 0..100 {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|elapsed| elapsed.as_nanos() as u64)
            .unwrap_or(0);
        let count = TEMPORARY_COUNTER.fetch_add(1, Ordering::Relaxed);
        let token = nanos ^ (u64::from(std::process::id()) << 20) ^ count.wrapping_mul(0x9e37_79b9);
        let path = directory.join(format!("{prefix}{token:x}"));
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .mode(0o600)
            .open(&path)
        {
            Ok(file) => return Ok((file, path)),
            Err(error) if error.kind() == io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(Failure::os(&error, &path)),
        }
    }
    reject("could not create a temporary file")
}

fn write_adoption_metadata(path: &Path, value: &[u8]) -> Done<()> {
    let linked = is_symlink(path);
    let links = fs::metadata(path).ok().map(|metadata| metadata.nlink());
    require(
        !linked && links.map(|count| count == 1).unwrap_or(true),
        "Refusing linked adoption metadata",
    )?;
    let directory = path.parent().unwrap_or(Path::new("."));
    let (mut file, temporary) = create_temporary(directory, ".adoption-metadata-")?;
    let outcome = (|| -> io::Result<()> {
        file.write_all(value)?;
        file.flush()?;
        file.sync_all()?;
        fs::rename(&temporary, path)
    })();
    if outcome.is_err() {
        let _ = fs::remove_file(&temporary);
    }
    outcome.map_err(|error| Failure::os(&error, path))
}

fn read_full(reader: &mut dyn Read, buffer: &mut [u8]) -> io::Result<usize> {
    let mut filled = 0;
    while filled < buffer.len() {
        let read = reader.read(&mut buffer[filled..])?;
        if read == 0 {
            break;
        }
        filled += read;
    }
    Ok(filled)
}

/// `shutil.copystat` without extended attributes: access and modification times, then permissions.
fn copy_times_and_mode(source: &Path, destination: &Path) -> Done<()> {
    let metadata = fs::metadata(source).map_err(|error| Failure::os(&error, source))?;
    let file = OpenOptions::new()
        .write(true)
        .open(destination)
        .map_err(|error| Failure::os(&error, destination))?;
    let mut times = fs::FileTimes::new();
    if let Ok(accessed) = metadata.accessed() {
        times = times.set_accessed(accessed);
    }
    if let Ok(modified) = metadata.modified() {
        times = times.set_modified(modified);
    }
    file.set_times(times)
        .map_err(|error| Failure::os(&error, destination))?;
    fs::set_permissions(
        destination,
        fs::Permissions::from_mode(metadata.mode() & 0o7777),
    )
    .map_err(|error| Failure::os(&error, destination))
}

/// Copies one member, resuming an owned partial copy, and checks it against the original SHA-256.
fn copy_verified_member(source: &Path, target: &Path, expected_sha: &str) -> Done<Json> {
    require(
        is_lower_hex_64(expected_sha),
        "Prepared member lacks its original SHA256",
    )?;
    let before = identity(source)?;
    require(
        !is_symlink(target),
        "Refusing symlinked adoption staging member",
    )?;
    let target_metadata = fs::metadata(target).ok();
    require(
        target_metadata
            .as_ref()
            .map(|metadata| metadata.nlink() == 1)
            .unwrap_or(true),
        "Refusing hardlinked adoption staging member",
    )?;
    let partial = target.with_file_name(format!("{}.copying", file_name_text(target)));
    require(
        !is_symlink(&partial),
        "Refusing symlinked adoption partial member",
    )?;
    let partial_metadata = fs::metadata(&partial).ok();
    require(
        partial_metadata
            .as_ref()
            .map(|metadata| metadata.nlink() == 1)
            .unwrap_or(true),
        "Refusing hardlinked adoption partial member",
    )?;
    if let Some(metadata) = &target_metadata {
        let same_bytes = metadata.is_file()
            && metadata.len() == before.size
            && sha256_file(source).map_err(|error| Failure::os(&error, source))? == expected_sha
            && sha256_file(target).map_err(|error| Failure::os(&error, target))? == expected_sha;
        require(
            same_bytes,
            "Completed adoption staging bytes differ from original receipt",
        )?;
    } else {
        require(
            partial_metadata
                .as_ref()
                .map(|metadata| metadata.is_file())
                .unwrap_or(true),
            "Invalid adoption partial member",
        )?;
        let mut hasher = Sha256::new();
        let mut original = File::open(source).map_err(|error| Failure::os(&error, source))?;
        let mut options = OpenOptions::new();
        options.read(true).write(true);
        if partial_metadata.is_none() {
            options.create_new(true);
        }
        let mut copied = options
            .open(&partial)
            .map_err(|error| Failure::os(&error, &partial))?;
        let mut remaining = copied.metadata()?.len();
        require(
            remaining <= before.size,
            "Adoption partial is longer than its original input",
        )?;
        let mut block = vec![0u8; COPY_BLOCK];
        let mut check = vec![0u8; COPY_BLOCK];
        while remaining > 0 {
            let wanted = remaining.min(COPY_BLOCK as u64) as usize;
            let got = read_full(&mut original, &mut block[..wanted])?;
            let mut same = got > 0;
            if same {
                let copied_got = read_full(&mut copied, &mut check[..got])?;
                same = copied_got == got && block[..got] == check[..got];
            }
            require(
                same,
                "Retained adoption partial prefix differs from the original input",
            )?;
            hasher.update(&block[..got]);
            remaining -= got as u64;
        }
        loop {
            let got = original.read(&mut block)?;
            if got == 0 {
                break;
            }
            copied.write_all(&block[..got])?;
            hasher.update(&block[..got]);
        }
        copied.flush()?;
        copied.sync_all()?;
        drop(copied);
        require(
            hasher.finish_hex() == expected_sha,
            "Original prepared member SHA256 differs from its receipt",
        )?;
        let partial_size = fs::metadata(&partial)
            .map_err(|error| Failure::os(&error, &partial))?
            .len();
        require(
            partial_size == before.size
                && sha256_file(&partial).map_err(|error| Failure::os(&error, &partial))?
                    == expected_sha,
            "Copied prepared member SHA256 differs from its receipt",
        )?;
        copy_times_and_mode(source, &partial)?;
        fs::rename(&partial, target).map_err(|error| Failure::os(&error, target))?;
    }
    require(
        identity(source)? == before,
        "Original prepared member changed during adoption",
    )?;
    Ok(Json::Object(vec![
        ("source".to_owned(), Json::Str(path_text(source))),
        ("destination".to_owned(), Json::Str(path_text(target))),
        ("bytes".to_owned(), Json::Int(before.size as i128)),
        ("sha256".to_owned(), Json::Str(expected_sha.to_owned())),
        ("source_identity".to_owned(), identity_json(&before)),
        (
            "copied_identity".to_owned(),
            identity_json(&identity(target)?),
        ),
    ]))
}

fn declares_both_roles(data: &Json) -> bool {
    let outputs = match data.get("outputs") {
        Some(found) => found,
        None => return false,
    };
    let two_roles = match outputs {
        Json::Object(items) => items.len() == 2,
        _ => false,
    };
    if !two_roles {
        return false;
    }
    ["snv", "indel"].iter().all(|role| match outputs.get(role) {
        Some(item) => matches!(item, Json::Object(_)) && is_string(item.get("path")),
        None => false,
    })
}

/// Copies a completed SpliceAI2 cache into the managed cache folder and returns its completion record.
pub(crate) fn adopt_completed_projection(
    source_dir: &Path,
    cache_dir: &Path,
    cache_root: &Path,
) -> Done<Json> {
    require(
        source_dir.is_absolute() && cache_root.is_absolute() && cache_dir.is_absolute(),
        "Prepared cache adoption paths must be absolute",
    )?;
    let source = resolve_strict(source_dir)?;
    let canonical_names = lexical_normalize(cache_root) == cache_root
        && lexical_normalize(cache_dir) == cache_dir
        && !cache_dir.components().any(|part| {
            let text = part.as_os_str().to_string_lossy().into_owned();
            text.to_lowercase().contains(".partial")
                || text == ".building"
                || text.starts_with(".building-")
        });
    require(
        canonical_names,
        "Prepared cache adoption needs complete canonical destination paths",
    )?;
    require(
        cache_dir.starts_with(cache_root) && cache_dir != cache_root,
        "Prepared cache adoption requires a dedicated destination inside cache_root",
    )?;
    require(
        !is_symlink(cache_root) || cache_root.is_dir(),
        "Selected cache_root alias must name an existing directory",
    )?;
    // A deliberately selected root alias is supported by the installer. Resolve
    // that boundary once; never follow links created beneath the managed root.
    let root = resolve_non_strict(cache_root);
    require(
        root != Path::new("/"),
        "Prepared cache adoption requires a dedicated cache_root directory",
    )?;
    let relative = cache_dir.strip_prefix(cache_root).unwrap_or(Path::new(""));
    let destination = root.join(relative);
    let traverses_link = destination
        .ancestors()
        .filter(|part| part.starts_with(&root))
        .any(is_symlink);
    require(
        !traverses_link,
        "Prepared cache adoption must not traverse destination symlinks",
    )?;
    require(
        !destination.starts_with(&source) && !source.starts_with(&destination),
        "Prepared cache source and destination must be separate directories",
    )?;
    let receipt_path = source.join("adapter.ready.json");
    complete_path(&receipt_path)?;
    let receipt_size = fs::metadata(&receipt_path)
        .map_err(|error| Failure::os(&error, &receipt_path))?
        .len();
    require(
        receipt_size <= MAX_HEADER as u64,
        "Original adapter receipt exceeds bounded size limit",
    )?;
    let original_before = identity(&receipt_path)?;
    let original = fs::read(&receipt_path).map_err(|error| Failure::os(&error, &receipt_path))?;
    let data = json_from_bytes(&original, true)?;
    require(
        declares_both_roles(&data),
        "Prepared source receipt must declare both VCF output roles",
    )?;
    let outputs = dict_index(&data, "outputs")?;
    let snv = string_of(dict_index(dict_index(outputs, "snv")?, "path")?)?.to_owned();
    let indel = string_of(dict_index(dict_index(outputs, "indel")?, "path")?)?.to_owned();
    let pair = validate_pair("GRCh38", Path::new(&snv), Path::new(&indel))?;
    require(
        pair.format == "vcf_spliceai2",
        "Prepared adoption requires completed SpliceAI2 VCFs",
    )?;
    let projection = completed_projection(Path::new(&snv), Path::new(&indel), true)?;
    let Projection {
        _lock: _source_lock,
        mut tracked,
    } = projection;
    require(
        identity(&receipt_path)? == original_before,
        "Original adapter receipt changed before adoption",
    )?;
    let snv_parent = fs::canonicalize(&snv).map(|path| path.parent().map(Path::to_path_buf));
    let indel_parent = fs::canonicalize(&indel).map(|path| path.parent().map(Path::to_path_buf));
    let same_source = match (snv_parent, indel_parent) {
        (Ok(Some(a)), Ok(Some(b))) => a == source && b == source,
        _ => false,
    };
    require(
        same_source,
        "Original prepared receipt belongs to a different source directory",
    )?;
    let mut members: Vec<(PathBuf, String)> = Vec::new();
    for role in ["snv", "indel"] {
        let item = dict_index(outputs, role)?;
        for (key, checksum) in [("path", "sha256"), ("index_path", "index_sha256")] {
            let path = complete_path(Path::new(string_of(dict_index(item, key)?)?))?;
            members.push((path, string_of(dict_index(item, checksum)?)?.to_owned()));
        }
    }
    let context = match data.get("gene_context") {
        Some(found) => found.clone(),
        None => Json::Object(Vec::new()),
    };
    require(
        matches!(context, Json::Object(_))
            && is_string(context.get("path"))
            && is_string(context.get("sha256")),
        "Prepared source receipt lacks gene-context provenance",
    )?;
    let context_path = complete_path(Path::new(string_of(dict_index(&context, "path")?)?))?;
    require(
        context_path.parent() == Some(source.as_path()),
        "Prepared gene context must be beside its original receipt",
    )?;
    members.push((
        context_path.clone(),
        string_of(dict_index(&context, "sha256")?)?.to_owned(),
    ));
    let distinct_names: HashSet<String> = members
        .iter()
        .map(|(path, _)| file_name_text(path))
        .collect();
    require(
        distinct_names.len() == 5,
        "Prepared cache members must have distinct filenames",
    )?;
    let context_before = identity(&context_path)?;
    tracked.push(Tracked {
        requested: path_text(&context_path),
        resolved: context_path.clone(),
        before: context_before,
    });
    let original_sha = sha256_hex(&original);
    let parent = destination.parent().unwrap_or(Path::new("/")).to_path_buf();
    fs::create_dir_all(&parent).map_err(|error| Failure::os(&error, &parent))?;
    let adoption_lock_path =
        parent.join(format!(".{}.adoption.lock", file_name_text(&destination)));
    let adoption_lock = OpenOptions::new()
        .append(true)
        .create(true)
        .open(&adoption_lock_path)
        .map_err(|error| Failure::os(&error, &adoption_lock_path))?;
    match adoption_lock.try_lock() {
        Ok(()) => {}
        Err(fs::TryLockError::WouldBlock) => {
            return Err(Failure::Rejected(
                "[Errno 11] Resource temporarily unavailable".to_owned(),
            ));
        }
        Err(fs::TryLockError::Error(error)) => {
            return Err(Failure::os(&error, &adoption_lock_path));
        }
    }
    if destination.exists() {
        require(
            destination.is_dir(),
            "Prepared adoption destination is occupied",
        )?;
        let copied_path = destination.join("adapter.ready.json");
        let copied_text =
            fs::read(&copied_path).map_err(|error| Failure::os(&error, &copied_path))?;
        let copied = json_from_file_bytes(&copied_text, false)?;
        let recorded = match dict_get(&copied, "relocation")? {
            Some(found) => dict_get(found, "original_receipt_sha256")?,
            None => None,
        };
        require(
            recorded.and_then(Json::as_str) == Some(original_sha.as_str()),
            "Prepared adoption destination belongs to a different original receipt",
        )?;
        let copied_outputs = dict_index(&copied, "outputs")?;
        let copied_snv =
            string_of(dict_index(dict_index(copied_outputs, "snv")?, "path")?)?.to_owned();
        let copied_indel =
            string_of(dict_index(dict_index(copied_outputs, "indel")?, "path")?)?.to_owned();
        {
            let _destination_lock =
                completed_projection(Path::new(&copied_snv), Path::new(&copied_indel), false)?;
            for (path, checksum) in &members {
                let target = destination.join(file_name_text(path));
                copy_verified_member(path, &target, checksum)?;
            }
        }
        unchanged_inputs(&tracked)?;
        return Ok(Json::Object(vec![
            ("kind".to_owned(), Json::Str(RELOCATION_KIND.to_owned())),
            ("status".to_owned(), Json::Str("reused".to_owned())),
            ("destination".to_owned(), Json::Str(path_text(&destination))),
            (
                "original_receipt_sha256".to_owned(),
                Json::Str(original_sha),
            ),
        ]));
    }
    let stage = parent.join(format!(
        ".{}.adoption-{}",
        file_name_text(&destination),
        &original_sha[..16]
    ));
    require(!is_symlink(&stage), "Refusing symlinked adoption workspace")?;
    match fs::create_dir(&stage) {
        Ok(()) => {}
        Err(error) if error.kind() == io::ErrorKind::AlreadyExists && stage.is_dir() => {}
        Err(error) => return Err(Failure::os(&error, &stage)),
    }
    let owner = stage.join("adoption.plan.json");
    let planned = Json::Object(vec![
        ("kind".to_owned(), Json::Str(RELOCATION_KIND.to_owned())),
        ("source".to_owned(), Json::Str(path_text(&source))),
        ("destination".to_owned(), Json::Str(path_text(&destination))),
        (
            "original_receipt_sha256".to_owned(),
            Json::Str(original_sha.clone()),
        ),
    ]);
    if owner.exists() {
        let existing = read_json(&owner, false)?;
        require(
            json_text::py_equal(&existing, &planned),
            "Occupied adoption workspace belongs to another source",
        )?;
    } else {
        let occupied = fs::read_dir(&stage)
            .map_err(|error| Failure::os(&error, &stage))?
            .next()
            .is_some();
        require(!occupied, "Occupied adoption workspace has no owned plan")?;
        write_adoption_metadata(&owner, &dumps_indented_line(&planned)?)?;
    }
    let mut records: Vec<Json> = Vec::new();
    for (path, checksum) in &members {
        let target = stage.join(file_name_text(path));
        let mut record = copy_verified_member(path, &target, checksum)?;
        let final_name = destination.join(file_name_text(&target));
        record.set("destination", Json::Str(path_text(&final_name)));
        records.push(record);
    }
    write_adoption_metadata(&stage.join("adapter.original.ready.json"), &original)?;
    let relocation = Json::Object(vec![
        ("kind".to_owned(), Json::Str(RELOCATION_KIND.to_owned())),
        (
            "original_receipt_sha256".to_owned(),
            Json::Str(original_sha.clone()),
        ),
        (
            "original_receipt".to_owned(),
            Json::Str(path_text(&destination.join("adapter.original.ready.json"))),
        ),
        (
            "record".to_owned(),
            Json::Str(path_text(&destination.join("relocation.ready.json"))),
        ),
    ]);
    let relocated = relocated_receipt(&data, &destination, &relocation)?;
    let ready_bytes = dumps_indented_line(&relocated)?;
    let mut record_items = planned.as_object().cloned().unwrap_or_default();
    record_items.push(("status".to_owned(), Json::Str("complete".to_owned())));
    record_items.push(("scientific_qualified".to_owned(), Json::Bool(false)));
    record_items.push(("members".to_owned(), Json::Array(records)));
    record_items.push((
        "relocated_receipt_sha256".to_owned(),
        Json::Str(sha256_hex(&ready_bytes)),
    ));
    record_items.push(("original_receipt_preserved".to_owned(), Json::Bool(true)));
    record_items.push(("converter_executed".to_owned(), Json::Bool(false)));
    let record = Json::Object(record_items);
    write_adoption_metadata(
        &stage.join("relocation.ready.json"),
        &dumps_indented_line(&record)?,
    )?;
    write_adoption_metadata(&stage.join("adapter.ready.json"), &ready_bytes)?;
    unchanged_inputs(&tracked)?;
    require(
        !destination.exists() && !is_symlink(&destination),
        "Prepared adoption destination became occupied",
    )?;
    fs::rename(&stage, &destination).map_err(|error| Failure::os(&error, &destination))?;
    let moved_outputs = dict_index(&relocated, "outputs")?;
    let moved_snv = string_of(dict_index(dict_index(moved_outputs, "snv")?, "path")?)?.to_owned();
    let moved_indel =
        string_of(dict_index(dict_index(moved_outputs, "indel")?, "path")?)?.to_owned();
    let _published = completed_projection(Path::new(&moved_snv), Path::new(&moved_indel), false)?;
    Ok(record)
}

// ---------------------------------------------------------------- configuration

struct ConfigRead {
    path: PathBuf,
    before: Identity,
    text: String,
    document: ConfigDocument,
}

fn read_config(config: &Path) -> Done<ConfigRead> {
    let path = complete_path(config)?;
    let before = identity(&path)?;
    let bytes = fs::read(&path).map_err(|error| Failure::os(&error, &path))?;
    // `path.open(newline="")` decodes as UTF-8 and keeps the line breaks as written.
    let text = python_utf8_decode(&bytes).map_err(Failure::Rejected)?;
    let document = match yaml_config::parse_config(&text) {
        Ok(YamlParsed::Mapping(document)) => document,
        Ok(_) => return reject("Expected a block-mapping YAML config"),
        Err(error @ YamlError::Invalid(_)) => {
            return reject(format!("Invalid YAML config: {error}"));
        }
        // YAML that PyYAML reads but this reader does not (anchors, aliases, tags,
        // block scalars, multi-line scalars, keys with spaces, deeper nesting) is
        // refused with exit code 1, so that it is not mistaken for the exit code 2
        // the script gives text that PyYAML itself refuses.
        Err(error @ YamlError::Unsupported(_)) => {
            return crash(format!("YAML error: {error}"));
        }
    };
    require(
        !document.has_duplicate_keys(),
        "Duplicate/unsupported YAML config key",
    )?;
    require(identity(&path)? == before, "Config changed during reading")?;
    Ok(ConfigRead {
        path,
        before,
        text,
        document,
    })
}

/// Replaces just the two scalar values, keeping the rest of the file as written.
fn register_vcfs(
    config: &Path,
    snv: &Path,
    indel: &Path,
    expected_identity: Identity,
    inputs: &[Tracked],
) -> Done<()> {
    let read = read_config(config)?;
    require(
        read.before == expected_identity,
        "Config changed after cache validation",
    )?;
    let selected: Vec<(&str, String)> = vec![
        (CACHE_KEYS[0], path_text(&absolute_path(snv))),
        (CACHE_KEYS[1], path_text(&absolute_path(indel))),
    ];
    unchanged_inputs(inputs)?;
    let already = selected.iter().all(|(key, value)| {
        matches!(read.document.scalar(key), Some(Scalar::Str(found)) if found == value)
    });
    if already {
        // Existing VCF reuse must leave even the quoting and modification time untouched.
        return Ok(());
    }
    let mut text = read.text.clone();
    let mut changes: Vec<(usize, usize, String)> = Vec::new();
    let mut additions = String::new();
    for (key, value) in &selected {
        let quoted = yaml_config::quoted_for_yaml(value);
        match read.document.entry(key) {
            None => additions.push_str(&format!("{key}: {quoted}\n")),
            Some(entry) => {
                require(
                    matches!(entry.value, ConfigValue::Scalar(_)),
                    &format!("Expected a single-line scalar for {key}"),
                )?;
                // An empty value has no space before it yet.
                let spaced = entry.value_start == 0
                    || text[..entry.value_start]
                        .chars()
                        .next_back()
                        .map(is_python_space)
                        .unwrap_or(false);
                let prefix = if spaced { "" } else { " " };
                changes.push((
                    entry.value_start,
                    entry.value_end,
                    format!("{prefix}{quoted}"),
                ));
            }
        }
    }
    changes.sort();
    changes.reverse();
    for (start, end, value) in changes {
        text = format!("{}{}{}", &text[..start], value, &text[end..]);
    }
    if !additions.is_empty() {
        if !text.ends_with('\n') {
            text.push('\n');
        }
        text.push_str(&additions);
    }
    // Validate the generated YAML before publishing either key.
    let updated = match yaml_config::parse_config(&text) {
        Ok(YamlParsed::Mapping(document)) => document,
        _ => return crash("the updated configuration is not a YAML mapping"),
    };
    let mut keys: Vec<&str> = read
        .document
        .entries
        .iter()
        .map(|e| e.key.as_str())
        .collect();
    for (key, _) in &selected {
        if !keys.contains(key) {
            keys.push(*key);
        }
    }
    let mut updated_keys: Vec<&str> = updated.entries.iter().map(|e| e.key.as_str()).collect();
    updated_keys.sort();
    updated_keys.dedup();
    let mut expected_keys = keys.clone();
    expected_keys.sort();
    let mut same = updated_keys == expected_keys;
    for key in &keys {
        let wanted = match selected.iter().find(|(name, _)| name == key) {
            Some((_, value)) => Some(ConfigValue::Scalar(Scalar::Str(value.clone()))),
            None => read.document.get(key).cloned(),
        };
        if wanted.as_ref() != updated.get(key) {
            same = false;
        }
    }
    require(same, "Config update would modify unrelated YAML values")?;
    let mode = fs::metadata(&read.path)
        .map_err(|error| Failure::os(&error, &read.path))?
        .mode()
        & 0o7777;
    let directory = read.path.parent().unwrap_or(Path::new("."));
    let (mut file, temporary) = create_temporary(directory, ".spliceai-config-")?;
    let outcome = (|| -> Done<()> {
        file.set_permissions(fs::Permissions::from_mode(mode))?;
        file.write_all(text.as_bytes())?;
        file.flush()?;
        file.sync_all()?;
        require(
            identity(&read.path)? == read.before,
            "Config changed before atomic registration",
        )?;
        unchanged_inputs(inputs)?;
        fs::rename(&temporary, &read.path).map_err(|error| Failure::os(&error, &read.path))
    })();
    if outcome.is_err() {
        let _ = fs::remove_file(&temporary);
    }
    outcome
}

// ---------------------------------------------------------------- SpliceAI v1 indexes

fn check_first_record(path: &Path, index_path: &Path) -> Done<()> {
    let wrap = |detail: String| -> Failure {
        Failure::Rejected(format!(
            "SpliceAI VCF/index validation failed for {}: {detail}",
            path.display()
        ))
    };
    let unwrap_failure = |failure: Failure| -> Failure {
        match failure {
            Failure::Rejected(message) => wrap(message),
            other => other,
        }
    };
    let file = File::open(path).map_err(|error| Failure::os(&error, path))?;
    let mut stream = BufReader::new(PythonGzipReader::new(file));
    // The script lets zlib.error escape as a traceback (exit 1); EOFError,
    // BadGzipFile and operating-system errors are reported as rejections.
    let read_failure = |error: io::Error| -> Failure {
        match GzipFailure::of(&error) {
            Some(GzipFailure::Zlib(text)) => Failure::Crash(format!("zlib.error: {text}")),
            Some(other) => wrap(other.to_string()),
            None => wrap(os_error_text(&error)),
        }
    };
    let mut line_bytes: Vec<u8> = Vec::new();
    let mut found_record = false;
    for _ in 0..=MAX_LINES {
        let line = read_bounded_line(&mut stream, MAX_LINE + 1).map_err(&read_failure)?;
        require(
            line.len() <= MAX_LINE,
            &format!(
                "SpliceAI record exceeds validation bound: {}",
                path.display()
            ),
        )
        .map_err(&unwrap_failure)?;
        if !line.is_empty() && line[0] != b'#' {
            line_bytes = line;
            found_record = true;
            break;
        }
        require(
            !line.is_empty(),
            &format!(
                "SpliceAI VCF contains no record for index validation: {}",
                path.display()
            ),
        )
        .map_err(&unwrap_failure)?;
    }
    if !found_record {
        return Err(wrap(format!(
            "SpliceAI header exceeds validation bound: {}",
            path.display()
        )));
    }
    // The script decodes in 8 KiB chunks, so a bad byte beyond the first line
    // can fail it earlier, and a position counts from the start of its chunk;
    // here only the first record is decoded and the position counts from its start.
    let text = python_io::decode_utf8(&line_bytes).map_err(&wrap)?;
    let line = text.trim_end_matches(['\r', '\n']);
    let fields: Vec<&str> = line.split('\t').collect();
    require(
        fields.len() >= 8,
        &format!("Malformed first SpliceAI VCF record: {}", path.display()),
    )
    .map_err(&unwrap_failure)?;
    let position = match python_int(fields[1]) {
        Some(position) => position,
        None => {
            return Err(wrap(format!(
                "invalid literal for int() with base 10: {}",
                python_str_repr(fields[1])
            )));
        }
    };
    require(
        position > 0,
        &format!("Invalid first SpliceAI VCF position: {}", path.display()),
    )
    .map_err(&unwrap_failure)?;
    // pysam reports every unreadable index with this one message, naming the data
    // file, after htslib has written its own line about the index to stderr
    // ("Invalid argument" is the reason htslib gives for an index it cannot parse).
    let index = match tabix_index::read_index(index_path) {
        Ok(index) => index,
        Err(_) => {
            eprintln!(
                "[E::hts_idx_load3] Could not load local index file '{}' : Invalid argument",
                index_path.display()
            );
            return Err(wrap(format!(
                "could not open index for `{}`",
                path.display()
            )));
        }
    };
    require(
        index.names.iter().any(|name| name == fields[0]),
        &format!(
            "SpliceAI index lacks first-record contig: {}",
            index_path.display()
        ),
    )
    .map_err(&unwrap_failure)?;
    let mut found = false;
    let mut bound_error: Option<String> = None;
    let begin = (position - 1) as i64;
    let outcome = tabix_index::fetch_lines(
        path,
        &index,
        fields[0],
        begin,
        begin + 1,
        &mut |number, row| {
            if number >= 10000 {
                bound_error = Some(format!(
                    "SpliceAI index query exceeds validation bound: {}",
                    path.display()
                ));
                return Err("query bound".to_owned());
            }
            if !row.is_ascii() {
                return Err("an indexed row is not ASCII text".to_owned());
            }
            if row == line.as_bytes() {
                found = true;
                return Ok(true);
            }
            Ok(false)
        },
    );
    if let Some(message) = bound_error {
        return Err(wrap(message));
    }
    outcome.map_err(&wrap)?;
    require(
        found,
        &format!(
            "SpliceAI index does not retrieve the original first record: {}",
            index_path.display()
        ),
    )
    .map_err(&unwrap_failure)
}

/// Bounded index checks for SpliceAI v1 VCFs; the original files stay as they are.
fn validate_v1_indexes(snv: &str, indel: &str) -> Done<Vec<Tracked>> {
    let mut inputs: Vec<Tracked> = Vec::new();
    for name in [snv, indel] {
        let path = complete_path(Path::new(name))?;
        let index_name = format!("{name}.tbi");
        let index = match complete_path(Path::new(&index_name)) {
            Ok(index) => index,
            Err(Failure::Rejected(_)) => {
                return reject(format!(
                    "Missing/incomplete SpliceAI tabix index: {index_name}; supply the matching .tbi or run tabix -p vcf on the completed original VCF"
                ));
            }
            Err(other) => return Err(other),
        };
        let path_before = identity(&path)?;
        let index_before = identity(&index)?;
        inputs.push(Tracked {
            requested: name.to_owned(),
            resolved: path.clone(),
            before: path_before,
        });
        inputs.push(Tracked {
            requested: index_name.clone(),
            resolved: index.clone(),
            before: index_before,
        });
        check_first_record(&path, &index)?;
    }
    unchanged_inputs(&inputs)?;
    Ok(inputs)
}

// ---------------------------------------------------------------- install

fn install(config: &str, cache_dir: &str, snv: &str, indel: &str) -> Done<PairReport> {
    let read = read_config(Path::new(config))?;
    let mut inputs: Vec<Tracked> = Vec::new();
    for name in [snv, indel] {
        let resolved = complete_path(Path::new(name))?;
        let before = identity(&resolved)?;
        inputs.push(Tracked {
            requested: name.to_owned(),
            resolved,
            before,
        });
    }
    let assembly = canonical_assembly_of_config(read.document.get("assembly"))?;
    let report = validate_pair(assembly, Path::new(snv), Path::new(indel))?;
    if report.format == "spliceai2_tsv" {
        let manifest = discover_manifest(Path::new(snv), Path::new(indel))?;
        // The script finds install_utils.sh beside itself. A compiled program has no
        // such place, so the folder comes from the environment; without it the run
        // stops with exit code 1 (a refusal, not the script's exit code 2).
        let script_dir = match std::env::var("PRIVA_SCRIPT_DIR") {
            Ok(folder) if !folder.is_empty() => folder,
            _ => {
                return crash(
                    "preparing a native SpliceAI2 table needs the environment variable PRIVA_SCRIPT_DIR set to the folder holding install_utils.sh",
                );
            }
        };
        let helper = Path::new(&script_dir).join("install_utils.sh");
        let output = absolute_path(Path::new(cache_dir))
            .join("SpliceAI2")
            .join("GRCh38");
        let config_absolute = absolute_path(Path::new(config));
        let again = resolve_strict(Path::new(config))?;
        require(
            again == read.path && identity(&read.path)? == read.before,
            "Config changed after cache validation",
        )?;
        unchanged_inputs(&inputs)?;
        let status = Command::new("bash")
            .arg(&helper)
            .arg("SpliceAI2_install")
            .arg(&config_absolute)
            .arg(&manifest)
            .arg(&output)
            .status()
            .map_err(|error| Failure::os(&error, Path::new("bash")))?;
        if !status.success() {
            return Err(Failure::Helper(match status.code() {
                Some(code) if code > 0 => code,
                _ => 1,
            }));
        }
    } else {
        let again = resolve_strict(Path::new(config))?;
        require(
            again == read.path,
            "Config target changed after cache validation",
        )?;
        if report.format == "vcf_spliceai2" {
            let projection = completed_projection(Path::new(snv), Path::new(indel), false)?;
            inputs.extend(projection.tracked.iter().cloned());
            register_vcfs(
                Path::new(config),
                Path::new(snv),
                Path::new(indel),
                read.before,
                &inputs,
            )?;
        } else {
            inputs.extend(validate_v1_indexes(snv, indel)?);
            register_vcfs(
                Path::new(config),
                Path::new(snv),
                Path::new(indel),
                read.before,
                &inputs,
            )?;
        }
    }
    Ok(report)
}

// ---------------------------------------------------------------- command line
//
// The script's `argparse` parser, reproduced with the program name this
// command has (`priva spliceai-cache-input`): the same options and abbreviations,
// the same order of checks, the same usage lines and help texts (worded and
// wrapped as argparse does at 80 columns), and the same error messages. The
// script's `install` takes exactly `--config --cache-dir --snv --indel`; the
// folder holding `install_utils.sh`, which Python finds beside the script, is
// taken from the environment variable PRIVA_SCRIPT_DIR instead of an option,
// so that the options (and what an abbreviation such as `--s` stands for) are
// the script's own.

const PROGRAM: &str = "priva spliceai-cache-input";
const MAIN_USAGE: &str = "usage: priva spliceai-cache-input [-h] {install,adopt-prepared} ...";
const INSTALL_USAGE: &str = "usage: priva spliceai-cache-input install [-h] --config CONFIG --cache-dir\n                                          CACHE_DIR --snv SNV --indel INDEL";
const ADOPT_USAGE: &str = "usage: priva spliceai-cache-input adopt-prepared [-h] --source-dir SOURCE_DIR\n                                                 --output-dir OUTPUT_DIR\n                                                 --cache-root CACHE_ROOT";
const MAIN_HELP: &str = "\nInspect existing SpliceAI inputs; convert only confirmed native SpliceAI2\ntables. Header classification reads bounded headers, not variant rows or\nwhole-file CRCs. Explicit prepared-cache adoption streams original/copy SHA256\nonce and preserves the original receipt; it does not rerun the converter or\nqualify scientific use. Fresh installation guidance is assembly-specific;\nexisting compatible v1 VCFs also remain usable on GRCh38. Neither\nclassification nor registration qualifies SpliceAI2 for the unchanged\nACMG/CADD scientific contracts.\n\npositional arguments:\n  {install,adopt-prepared}\n    adopt-prepared      Copy a completed SpliceAI2 cache into managed storage\n                        without conversion\n\noptions:\n  -h, --help            show this help message and exit\n";
const INSTALL_HELP: &str = "\noptions:\n  -h, --help            show this help message and exit\n  --config CONFIG\n  --cache-dir CACHE_DIR\n  --snv SNV\n  --indel INDEL\n";
const ADOPT_HELP: &str = "\noptions:\n  -h, --help            show this help message and exit\n  --source-dir SOURCE_DIR\n  --output-dir OUTPUT_DIR\n  --cache-root CACHE_ROOT\n";
const COMMANDS: &str = "'install', 'adopt-prepared'";

/// One of the script's sub-commands.
struct SubCommand {
    name: &'static str,
    usage: &'static str,
    help: &'static str,
    /// The required options, without their dashes, in the order they are declared.
    options: &'static [&'static str],
}

const INSTALL: SubCommand = SubCommand {
    name: "install",
    usage: INSTALL_USAGE,
    help: INSTALL_HELP,
    options: &["config", "cache-dir", "snv", "indel"],
};

const ADOPT_PREPARED: SubCommand = SubCommand {
    name: "adopt-prepared",
    usage: ADOPT_USAGE,
    help: ADOPT_HELP,
    options: &["source-dir", "output-dir", "cache-root"],
};

/// `parser.error(message)`: the usage, then the message, exit code 2.
fn argparse_error(usage: &str, program: &str, message: &str) -> ! {
    eprintln!("{usage}");
    eprintln!("{program}: error: {message}");
    std::process::exit(2)
}

/// How argparse reads one argument string.
enum Token {
    /// A value or a command word.
    Plain,
    /// An option: the position of the option it names among the parser's option
    /// strings (`None` when the parser does not know it) and an explicit value.
    Option(Option<usize>, Option<String>),
}

/// argparse's `_negative_number_matcher`: `-5`, `-5.5` and `-.5`.
fn is_negative_number(text: &str) -> bool {
    let Some(body) = text.strip_prefix('-') else {
        return false;
    };
    let digits = |part: &str| !part.is_empty() && part.chars().all(|c| c.is_ascii_digit());
    match body.split_once('.') {
        Some((whole, fraction)) => (whole.is_empty() || digits(whole)) && digits(fraction),
        None => digits(body),
    }
}

/// `ArgumentParser._parse_optional` for a parser with these option strings
/// (`-h`, `--help` and the long options). An ambiguous abbreviation is an error message.
fn classify_token(text: &str, strings: &[String]) -> std::result::Result<Token, String> {
    if !text.starts_with('-') {
        return Ok(Token::Plain);
    }
    if let Some(found) = strings.iter().position(|name| name == text) {
        return Ok(Token::Option(Some(found), None));
    }
    if text.len() == 1 {
        return Ok(Token::Plain);
    }
    if let Some((name, value)) = text.split_once('=') {
        if let Some(found) = strings.iter().position(|known| known == name) {
            return Ok(Token::Option(Some(found), Some(value.to_owned())));
        }
    }
    // Abbreviations: a long option may be shortened, a short one may carry its value.
    let mut matches: Vec<(usize, Option<String>)> = Vec::new();
    if text.starts_with("--") {
        let (prefix, value) = match text.split_once('=') {
            Some((prefix, value)) => (prefix, Some(value.to_owned())),
            None => (text, None),
        };
        for (position, name) in strings.iter().enumerate() {
            if name.starts_with(prefix) {
                matches.push((position, value.clone()));
            }
        }
    } else {
        let short: String = text.chars().take(2).collect();
        let attached: String = text.chars().skip(2).collect();
        for (position, name) in strings.iter().enumerate() {
            if *name == short {
                matches.push((position, Some(attached.clone())));
            } else if name.starts_with(text) {
                matches.push((position, None));
            }
        }
    }
    if matches.len() > 1 {
        let names: Vec<&str> = matches
            .iter()
            .map(|(at, _)| strings[*at].as_str())
            .collect();
        return Err(format!(
            "ambiguous option: {text} could match {}",
            names.join(", ")
        ));
    }
    if let Some((position, value)) = matches.pop() {
        return Ok(Token::Option(Some(position), value));
    }
    if is_negative_number(text) || text.contains(' ') {
        return Ok(Token::Plain);
    }
    Ok(Token::Option(None, None))
}

/// `-h` was asked for (in a short or abbreviated spelling, with a value attached or not).
fn help_requested(text: &str, explicit: &Option<String>) -> bool {
    // `--help=x` is an error; `-hx` shows the help and keeps `-x` as an extra.
    explicit.is_none() || !text.starts_with("--")
}

fn print_help(usage: &str, body: &str) -> ! {
    print!("{usage}\n{body}");
    std::process::exit(0)
}

/// What the sub-command's parser found: its option values and the arguments it did not use.
struct SubParsed {
    values: HashMap<String, String>,
    extras: Vec<String>,
}

fn parse_sub_command(command: &SubCommand, args: &[String]) -> SubParsed {
    let program = format!("{PROGRAM} {}", command.name);
    let fail = |message: &str| -> ! { argparse_error(command.usage, &program, message) };
    let mut strings: Vec<String> = vec!["-h".to_owned(), "--help".to_owned()];
    strings.extend(command.options.iter().map(|name| format!("--{name}")));
    // First every argument is read as an option or a plain word; after `--`
    // everything is plain. An ambiguous abbreviation stops here.
    let mut kinds: Vec<Option<Token>> = Vec::with_capacity(args.len());
    let mut after_separator = false;
    for text in args {
        if after_separator {
            kinds.push(Some(Token::Plain));
        } else if text == "--" {
            after_separator = true;
            kinds.push(None);
        } else {
            match classify_token(text, &strings) {
                Ok(kind) => kinds.push(Some(kind)),
                Err(message) => fail(&message),
            }
        }
    }
    let mut values: HashMap<String, String> = HashMap::new();
    let mut extras: Vec<String> = Vec::new();
    let mut index = 0;
    while index < args.len() {
        let (position, explicit) = match &kinds[index] {
            Some(Token::Option(Some(position), explicit)) => (*position, explicit.clone()),
            Some(Token::Option(None, _)) | Some(Token::Plain) | None => {
                extras.push(args[index].clone());
                index += 1;
                continue;
            }
        };
        index += 1;
        if position < 2 {
            if help_requested(&args[index - 1], &explicit) {
                print_help(command.usage, command.help);
            }
            fail(&format!(
                "argument -h/--help: ignored explicit argument {}",
                python_str_repr(explicit.as_deref().unwrap_or(""))
            ));
        }
        let name = strings[position].clone();
        let value = match explicit {
            Some(value) => value,
            None => {
                // One plain word follows. An option does not take `--` as its separator
                // (argparse leaves it out of an option's pattern), so `--config -- x` fails.
                if index < args.len() && matches!(kinds[index], Some(Token::Plain)) {
                    let value = args[index].clone();
                    index += 1;
                    value
                } else {
                    fail(&format!("argument {name}: expected one argument"))
                }
            }
        };
        values.insert(name[2..].to_owned(), value);
    }
    let missing: Vec<String> = command
        .options
        .iter()
        .filter(|name| !values.contains_key(**name))
        .map(|name| format!("--{name}"))
        .collect();
    if !missing.is_empty() {
        fail(&format!(
            "the following arguments are required: {}",
            missing.join(", ")
        ));
    }
    SubParsed { values, extras }
}

fn finish(outcome: Done<Json>) -> Result<()> {
    match outcome {
        Ok(report) => {
            match json_text::dumps(&report, true) {
                Ok(text) => println!("{text}"),
                Err(error) => {
                    eprintln!("priva spliceai-cache-input: {error}");
                    std::process::exit(1);
                }
            }
            Ok(())
        }
        Err(Failure::Rejected(message)) => {
            eprintln!("SpliceAI cache input rejected: {message}");
            std::process::exit(2)
        }
        Err(Failure::Crash(message)) => {
            eprintln!("priva spliceai-cache-input: unhandled error: {message}");
            std::process::exit(1)
        }
        Err(Failure::Helper(code)) => {
            eprintln!("SpliceAI cache preparation failed (exit {code})");
            std::process::exit(code)
        }
    }
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let main_strings: Vec<String> = vec!["-h".to_owned(), "--help".to_owned()];
    let main_fail = |message: &str| -> ! { argparse_error(MAIN_USAGE, PROGRAM, message) };
    let invalid_choice = |word: &str| -> ! {
        main_fail(&format!(
            "argument command: invalid choice: {} (choose from {COMMANDS})",
            python_str_repr(word)
        ))
    };
    // Options before the command word are the main parser's: only the help is
    // known, anything else is set aside and reported once the command has run.
    let mut main_extras: Vec<String> = Vec::new();
    let mut found: Option<usize> = None;
    let mut index = 0;
    while index < args.len() {
        let text = &args[index];
        if text == "--" {
            if index + 1 < args.len() {
                invalid_choice("--");
            }
            break;
        }
        match classify_token(text, &main_strings) {
            Ok(Token::Plain) => {
                found = Some(index);
                break;
            }
            Ok(Token::Option(Some(_), explicit)) => {
                if help_requested(text, &explicit) {
                    print_help(MAIN_USAGE, MAIN_HELP);
                }
                main_fail(&format!(
                    "argument -h/--help: ignored explicit argument {}",
                    python_str_repr(explicit.as_deref().unwrap_or(""))
                ));
            }
            Ok(Token::Option(None, _)) => main_extras.push(text.clone()),
            Err(message) => main_fail(&message),
        }
        index += 1;
    }
    let Some(at) = found else {
        main_fail("the following arguments are required: command");
    };
    let command = match args[at].as_str() {
        "install" => INSTALL,
        "adopt-prepared" => ADOPT_PREPARED,
        other => invalid_choice(other),
    };
    let SubParsed { values, extras } = parse_sub_command(&command, &args[at + 1..]);
    main_extras.extend(extras);
    if !main_extras.is_empty() {
        main_fail(&format!(
            "unrecognized arguments: {}",
            main_extras.join(" ")
        ));
    }
    if command.name == "install" {
        let report = install(
            &values["config"],
            &values["cache-dir"],
            &values["snv"],
            &values["indel"],
        );
        finish(report.map(|report| report_json(&report)))
    } else {
        finish(adopt_completed_projection(
            Path::new(&values["source-dir"]),
            Path::new(&values["output-dir"]),
            Path::new(&values["cache-root"]),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;

    const VCF_HEADER: &str =
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    fn scratch(name: &str) -> PathBuf {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let unique = NEXT.fetch_add(1, Ordering::Relaxed);
        let folder = std::env::temp_dir().join(format!(
            "priva_spliceai_cache_input_{}_{unique}_{name}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&folder);
        fs::create_dir_all(&folder).unwrap();
        folder
    }

    fn write_text(folder: &Path, name: &str, text: &str, compressed: bool) -> PathBuf {
        let path = folder.join(name);
        if compressed {
            let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
            encoder.write_all(text.as_bytes()).unwrap();
            fs::write(&path, encoder.finish().unwrap()).unwrap();
        } else {
            fs::write(&path, text).unwrap();
        }
        path
    }

    fn header_of(text: &str) -> Done<Header> {
        let folder = scratch("header");
        let path = write_text(&folder, "scores", text, false);
        let result = inspect_header(&path);
        let _ = fs::remove_dir_all(&folder);
        result
    }

    fn with_metadata(metadata: &str) -> String {
        format!("##fileformat=VCFv4.2\n{metadata}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    }

    #[test]
    fn the_native_table_header_is_the_published_one() {
        let columns = native_columns();
        assert_eq!(columns.len(), 84);
        let line = format!("{}\n", columns.join("\t"));
        assert_eq!(
            sha256_hex(line.as_bytes()),
            "1feb799bed12dae1030e64b03930eab470340848aa4a7e0166072a140e7b0d66"
        );
    }

    #[test]
    fn headers_are_classified() {
        assert_eq!(
            header_of(VCF_HEADER).unwrap(),
            Header {
                format: "vcf_v1",
                assembly: None
            }
        );
        let native = format!("{}\n", native_columns().join("\t"));
        assert_eq!(
            header_of(&native).unwrap(),
            Header {
                format: "spliceai2_tsv",
                assembly: Some("GRCh38")
            }
        );
        for (metadata, format, assembly) in [
            (
                "##reference=file:///refs/hg19.fa\n",
                "vcf_v1",
                Some("GRCh37"),
            ),
            (
                "##contig=<ID=1,length=10,assembly=\"GRCh38.p14\">\n",
                "vcf_v1",
                Some("GRCh38"),
            ),
            (
                "##source=PriVA_SpliceAI2_projection\n",
                "vcf_spliceai2",
                Some("GRCh38"),
            ),
            ("##SpliceAI2Version=2\n", "vcf_spliceai2", Some("GRCh38")),
            (
                "##SpliceAIVersion=2.0\n##reference=GRCh38\n",
                "vcf_spliceai2",
                Some("GRCh38"),
            ),
        ] {
            assert_eq!(
                header_of(&with_metadata(metadata)).unwrap(),
                Header { format, assembly },
                "{metadata:?}"
            );
        }
    }

    #[test]
    fn bad_headers_are_rejected() {
        let columns = VCF_COLUMNS.join("\t");
        for text in [
            VCF_HEADER.replace("POS", "POSITION"),
            VCF_HEADER.split_once('\n').unwrap().1.to_owned(),
            with_metadata("##fileformat=VCFv4.2\n"),
            with_metadata("##INFO=<ID=x\n"),
            with_metadata("##INFO=<ID=x,Description=\"unclosed>\n"),
            with_metadata("##broken\n"),
            with_metadata("##SpliceAI2Version=1\n"),
            with_metadata("##SpliceAI2Version=2\n##SpliceAIVersion=1\n"),
            with_metadata("##SpliceAI2Version=2\n##assembly=hg19\n"),
            with_metadata("##reference=hg19\n##assembly=GRCh38\n"),
            format!("##fileformat=VCFv4.2\n{columns}\tFORMAT\n"),
            format!("##fileformat=VCFv4.2\n{columns}\tFORMAT\tS\tS\n"),
            "##fileformat=VCFv4.2\n1\t10\t.\tA\tC\t.\t.\t.\n".to_owned(),
            VCF_HEADER.trim_end().to_owned(),
            VCF_HEADER.replace("POS", "PO\u{0}S"),
        ] {
            assert!(
                matches!(header_of(&text), Err(Failure::Rejected(_))),
                "{text:?}"
            );
        }
    }

    #[test]
    fn gzip_is_found_by_its_magic_and_the_body_is_not_read() {
        let folder = scratch("gzip");
        let packed = write_text(&folder, "scores.anything", VCF_HEADER, true);
        assert_eq!(inspect_header(&packed).unwrap().format, "vcf_v1");
        let body = format!("{VCF_HEADER}NOT_A_VARIANT\n{}", "x".repeat(100_000));
        let long = write_text(&folder, "long", &body, false);
        assert_eq!(inspect_header(&long).unwrap().format, "vcf_v1");
        let truncated = folder.join("truncated");
        let mut bytes = fs::read(&packed).unwrap();
        bytes.truncate(bytes.len() - 6);
        fs::write(&truncated, bytes).unwrap();
        assert!(matches!(
            inspect_header(&truncated),
            Err(Failure::Rejected(_))
        ));
        let staged = folder.join("scores.partial.gz");
        fs::write(&staged, "x").unwrap();
        assert!(matches!(
            inspect_header(&staged),
            Err(Failure::Rejected(message)) if message.contains("Incomplete")
        ));
        let _ = fs::remove_dir_all(&folder);
    }

    #[test]
    fn header_limits_are_enforced() {
        let folder = scratch("limits");
        let path = write_text(&folder, "scores", VCF_HEADER, false);
        for (limits, expected) in [
            (
                Limits {
                    line: 10,
                    header: MAX_HEADER,
                    lines: MAX_LINES,
                },
                "bounded",
            ),
            (
                Limits {
                    line: MAX_LINE,
                    header: 50,
                    lines: MAX_LINES,
                },
                "bounded",
            ),
            (
                Limits {
                    line: MAX_LINE,
                    header: MAX_HEADER,
                    lines: 1,
                },
                "bounded",
            ),
        ] {
            match inspect_header_limited(&path, &limits) {
                Err(Failure::Rejected(message)) => assert!(message.contains(expected), "{message}"),
                other => panic!("{other:?}"),
            }
        }
        let _ = fs::remove_dir_all(&folder);
    }

    #[test]
    fn assembly_and_version_text_rules() {
        assert_eq!(assemblies_in("file:///refs/hg19.fa"), vec!["GRCh37"]);
        assert_eq!(assemblies_in("xhg19 hg38y GRCH38"), vec!["GRCh38"]);
        assert!(assemblies_in("hg1938").is_empty());
        assert_eq!(
            contig_assembly_values("ID=1,length=10,assembly=\"GRCh38.p14\",species=x"),
            vec!["\"GRCh38.p14\"".to_owned()]
        );
        assert_eq!(
            contig_assembly_values("assembly=hg19,ID=2"),
            vec!["hg19".to_owned()]
        );
        assert!(contig_assembly_values("ID=1,xassembly=hg19").is_empty());
        for (text, expected) in [
            ("PriVA_SpliceAI2_projection", true),
            ("SpliceAI 2", true),
            ("spliceai-v2.0", true),
            ("SpliceAI2x", false),
            ("SpliceAI v1", false),
            ("spliceai22", false),
            ("my spliceai_2_", true),
        ] {
            assert_eq!(mentions_spliceai2(text), expected, "{text}");
        }
        assert!(is_version_text("2.0.1", "2"));
        assert!(!is_version_text("1.0", "2"));
        assert!(is_version_text("1", "12"));
        assert!(!is_version_text("2.", "12"));
    }

    #[test]
    fn assemblies_have_aliases_and_a_python_style_message() {
        assert_eq!(canonical_assembly("hg19").unwrap(), "GRCh37");
        assert_eq!(canonical_assembly("GRCH38").unwrap(), "GRCh38");
        match canonical_assembly("hg37") {
            Err(Failure::Rejected(message)) => assert_eq!(message, "Unsupported assembly: 'hg37'"),
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn paths_are_read_the_way_pathlib_reads_them() {
        assert_eq!(path_name("/a/b/c.gz"), "c.gz");
        assert_eq!(path_name("/a/b/"), "b");
        assert_eq!(path_name("/"), "");
        assert_eq!(path_name("a/.."), "..");
        assert_eq!(
            lexical_normalize(Path::new("/a/./b/../c//d")),
            PathBuf::from("/a/c/d")
        );
        assert!(
            is_unfinished("X.PARTIAL.gz")
                && is_unfinished(".building-x")
                && !is_unfinished("build")
        );
    }

    #[test]
    fn the_receipt_of_a_copied_cache_moves_every_location() {
        let original = json_text::parse(
            "{\"outputs\": {\"snv\": {\"path\": \"/o/a.vcf.gz\", \"index_path\": \"/o/a.vcf.gz.tbi\"}}, \"gene_context\": {\"path\": \"/o/g.tsv\", \"sha256\": \"x\"}, \"other\": 1}",
            true,
        )
        .unwrap();
        let relocation = json_text::parse("{\"kind\": \"k\"}", true).unwrap();
        let moved = relocated_receipt(&original, Path::new("/new"), &relocation).unwrap();
        assert_eq!(
            json_text::dumps(&moved, false).unwrap(),
            "{\"outputs\": {\"snv\": {\"path\": \"/new/a.vcf.gz\", \"index_path\": \"/new/a.vcf.gz.tbi\"}}, \"gene_context\": {\"path\": \"/new/g.tsv\", \"sha256\": \"x\"}, \"other\": 1, \"relocation\": {\"kind\": \"k\"}}"
        );
        let without_context = json_text::parse("{\"outputs\": {}}", true).unwrap();
        assert!(matches!(
            relocated_receipt(&without_context, Path::new("/new"), &relocation),
            Err(Failure::Crash(_))
        ));
    }

    #[test]
    fn registration_replaces_only_the_two_values() {
        let folder = scratch("register");
        let snv = write_text(&folder, "snv.vcf", VCF_HEADER, false);
        let indel = write_text(&folder, "indel.vcf", VCF_HEADER, false);
        let config = folder.join("config.yaml");
        let original = "# before\r\nassembly: hg38\r\nspliceai_snv_prescore: \"old\" # snv comment\r\nother: [one, two] # unchanged\r\nspliceai_indel_prescore: old2 # indel comment\r\n";
        fs::write(&config, original).unwrap();
        fs::set_permissions(&config, fs::Permissions::from_mode(0o640)).unwrap();
        let before = identity(&fs::canonicalize(&config).unwrap()).unwrap();
        register_vcfs(&config, &snv, &indel, before, &[]).unwrap();
        let expected = original
            .replace("\"old\"", &json_text::quote(&path_text(&snv), false))
            .replace("old2", &json_text::quote(&path_text(&indel), false));
        assert_eq!(fs::read_to_string(&config).unwrap(), expected);
        assert_eq!(fs::metadata(&config).unwrap().mode() & 0o777, 0o640);
        // Registering the same values again changes nothing at all.
        let settled = identity(&fs::canonicalize(&config).unwrap()).unwrap();
        register_vcfs(&config, &snv, &indel, settled, &[]).unwrap();
        assert_eq!(identity(&config).unwrap(), settled);
        // Missing keys are added, empty ones filled, comments kept.
        fs::write(
            &config,
            "assembly: hg38\nspliceai_snv_prescore: # keep\nspliceai_indel_prescore:\n",
        )
        .unwrap();
        let again = identity(&fs::canonicalize(&config).unwrap()).unwrap();
        register_vcfs(&config, &snv, &indel, again, &[]).unwrap();
        let text = fs::read_to_string(&config).unwrap();
        assert!(
            text.contains("# keep")
                && text.contains(&path_text(&snv))
                && text.contains(&path_text(&indel))
        );
        fs::write(&config, "assembly: hg38").unwrap();
        let last = identity(&fs::canonicalize(&config).unwrap()).unwrap();
        register_vcfs(&config, &snv, &indel, last, &[]).unwrap();
        assert!(
            fs::read_to_string(&config)
                .unwrap()
                .starts_with("assembly: hg38\nspliceai_snv_prescore: ")
        );
        let _ = fs::remove_dir_all(&folder);
    }

    #[test]
    fn a_partial_copy_is_resumed_and_checked() {
        let folder = scratch("copy");
        let source = write_text(&folder, "member.bin", "0123456789abcdef", false);
        let sha = sha256_hex(b"0123456789abcdef");
        let target = folder.join("stage_member.bin");
        fs::write(folder.join("stage_member.bin.copying"), "01234").unwrap();
        let record = copy_verified_member(&source, &target, &sha).unwrap();
        assert_eq!(fs::read_to_string(&target).unwrap(), "0123456789abcdef");
        assert!(!folder.join("stage_member.bin.copying").exists());
        assert_eq!(
            record.get("bytes").map(|v| matches!(v, Json::Int(16))),
            Some(true)
        );
        // A second call finds the finished copy and verifies it.
        copy_verified_member(&source, &target, &sha).unwrap();
        // A partial prefix that differs from the original is refused.
        let other = folder.join("other_member.bin");
        fs::write(folder.join("other_member.bin.copying"), "xx").unwrap();
        assert!(matches!(
            copy_verified_member(&source, &other, &sha),
            Err(Failure::Rejected(message)) if message.contains("prefix differs")
        ));
        // A wrong checksum is refused.
        let wrong = folder.join("wrong_member.bin");
        assert!(matches!(
            copy_verified_member(&source, &wrong, &"0".repeat(64)),
            Err(Failure::Rejected(message)) if message.contains("SHA256 differs")
        ));
        let _ = fs::remove_dir_all(&folder);
    }
    #[test]
    fn command_lines_are_read_like_argparse() {
        let strings: Vec<String> = [
            "-h",
            "--help",
            "--config",
            "--cache-dir",
            "--snv",
            "--indel",
        ]
        .iter()
        .map(|name| name.to_string())
        .collect();
        let named = |text: &str| match classify_token(text, &strings) {
            Ok(Token::Option(Some(at), explicit)) => Some((at, explicit)),
            _ => None,
        };
        assert_eq!(named("--config"), Some((2, None)));
        assert_eq!(named("--conf"), Some((2, None)));
        assert_eq!(named("--snv=x=y"), Some((4, Some("x=y".to_owned()))));
        assert_eq!(named("--he"), Some((1, None)));
        assert_eq!(named("-hx"), Some((0, Some("x".to_owned()))));
        assert!(matches!(
            classify_token("--bogus", &strings),
            Ok(Token::Option(None, None))
        ));
        assert!(matches!(classify_token("-5", &strings), Ok(Token::Plain)));
        assert!(matches!(classify_token("-.5", &strings), Ok(Token::Plain)));
        assert!(matches!(classify_token("-x y", &strings), Ok(Token::Plain)));
        assert!(matches!(classify_token("-", &strings), Ok(Token::Plain)));
        assert!(matches!(classify_token("", &strings), Ok(Token::Plain)));
        assert_eq!(
            classify_token("--c", &strings).err().unwrap(),
            "ambiguous option: --c could match --config, --cache-dir"
        );
        let words =
            |items: &[&str]| -> Vec<String> { items.iter().map(|item| item.to_string()).collect() };
        let parsed = parse_sub_command(
            &INSTALL,
            &words(&[
                "--conf",
                "a",
                "--cache-dir=b",
                "--indel",
                "",
                "--snv",
                "c",
                "--",
                "stray",
                "--bogus",
            ]),
        );
        assert_eq!(parsed.values["config"], "a");
        assert_eq!(parsed.values["cache-dir"], "b");
        assert_eq!(parsed.values["snv"], "c");
        assert_eq!(parsed.values["indel"], "");
        assert_eq!(parsed.extras, words(&["--", "stray", "--bogus"]));
    }

    #[test]
    fn receipts_with_decimal_numbers_are_written_back() {
        let original = json_text::parse(
            "{\"outputs\": {}, \"gene_context\": {\"path\": \"/o/g.tsv\"}, \"thresholds\": {\"balanced\": 0.250, \"low\": 1e-5}}",
            true,
        )
        .unwrap();
        let relocation = json_text::parse("{\"kind\": \"k\"}", true).unwrap();
        let moved = relocated_receipt(&original, Path::new("/n"), &relocation).unwrap();
        let text = json_text::dumps(&moved, false).unwrap();
        assert!(
            text.contains("\"thresholds\": {\"balanced\": 0.25, \"low\": 1e-05}"),
            "{text}"
        );
    }

    #[test]
    fn damaged_headers_give_the_scripts_texts() {
        let folder = scratch("damaged");
        let packed = write_text(&folder, "ok.gz", VCF_HEADER, true);
        let good = fs::read(&packed).unwrap();
        let message = |name: &str, bytes: &[u8]| -> String {
            let path = folder.join(name);
            fs::write(&path, bytes).unwrap();
            match inspect_header(&path) {
                Err(Failure::Rejected(message)) => message,
                other => panic!("{name}: {other:?}"),
            }
        };
        let cut = message("cut.gz", &good[..good.len() - 6]);
        assert!(cut.starts_with("Corrupt/truncated cache header: "), "{cut}");
        assert!(
            cut.ends_with(&format!(": {COMPRESSED_FILE_ENDED}")),
            "{cut}"
        );
        let mut garbage = good.clone();
        garbage.extend_from_slice(b"garbage");
        assert_eq!(
            message("garbage.gz", &garbage),
            "Not a gzipped file (b'ga')"
        );
        let mut method = good.clone();
        method[2] = 9;
        assert_eq!(message("method.gz", &method), "Unknown compression method");
        let mut bytes = b"##fileformat=VCFv4.2\n##source=".to_vec();
        bytes.extend_from_slice(b"\xff\n");
        let coded = message("coded", &bytes);
        assert!(
            coded.ends_with(
                "'utf-8' codec can't decode byte 0xff in position 9: invalid start byte"
            ),
            "{coded}"
        );
        let _ = fs::remove_dir_all(&folder);
    }
}
