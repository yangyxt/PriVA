//! Native `scripts/spliceai2_to_vep_vcf.py`: `prepare` and `gene-context` for
//! the verified, annotation-only SpliceAI2 GRCh38 view. No ACMG policy is
//! activated; the original 84-column tables remain authoritative.
//!
//! Manifest, input and tool verification, the GENCODE48 gene context with the
//! VEP release's gene names, the gzip | projection | bgzip pipeline per variant
//! kind, the tabix index and first-record probe, the receipt and the
//! publication order follow the script. Printed lines and exit statuses follow
//! it too; a failure the script reported as `SpliceAI2 preparation failed:
//! message` exits 1 here as `priva: CATEGORY: message`. Differences, by design
//! of the Rust migration:
//! - The projection runs as this executable's `spliceai2-vcf`, the same
//!   conversion the script compiled from `spliceai2_converter.rs` with rustc at
//!   installation; `build-converter` is not provided. The receipt's
//!   `converter` names this executable (`role` `priva_builtin`) and
//!   `driver_sha256` is its SHA-256. `--converter-bin` still runs an explicit
//!   standalone converter's `project`.
//! - The subcommand comes first (the script also accepted bare prepare flags).
//! - `--jobs N` (default 1, which keeps the script's single pipeline) checks
//!   the source shards' SHA-256 and converts the shards of each kind N at a
//!   time, each into its own BGZF part, then joins the parts in chromosome
//!   order. The decompressed VCF and the receipt, apart from the output bytes
//!   and checksums, are the same as with one job.
use crate::builder_error;
use crate::cli::{Result, invalid};
use crate::installer_config_paths::files::{
    is_file, is_link, metadata, mkdir, native, os_error, rename,
};
use crate::json_text::{self, DumpsOptions, Json};
use crate::python_args::{parse_long_option_occurrences, require_options};
use crate::python_io::PythonGzipReader;
use crate::python_path::{absolute, repr};
use crate::python_subprocess::signal_repr;
use crate::python_text::is_python_space;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::ffi::CString;
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufRead, BufReader, Read, Write};
use std::os::fd::AsFd;
use std::os::unix::fs::MetadataExt;
use std::os::unix::process::ExitStatusExt;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

const HELP: &str = "usage: priva spliceai2-prepare {prepare,gene-context}\n\
       prepare --manifest MANIFEST --assembly ASSEMBLY --output-dir OUTPUT_DIR --receipt RECEIPT\n\
               [--converter-bin CONVERTER_BIN] [--threads THREADS] [--jobs JOBS]\n\
               [--validation-subset]\n\
       gene-context --regions REGIONS --gtf GTF --output OUTPUT [--vep-annotation VEP_ANNOTATION\n\
               --ensembl-release ENSEMBL_RELEASE]\n\
Prepare a verified, annotation-only SpliceAI2 view; never activate ACMG policy. The output and\n\
exit statuses of scripts/spliceai2_to_vep_vcf.py; the projection runs as priva spliceai2-vcf.";
const REVISION: &str = "cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9";
const REGION_COLUMNS: [&str; 5] = ["gene_id", "chrom", "gene_start", "gene_end", "strand"];
const CONTEXT_NAME: &str = "SpliceAI2_GRCh38.gene_context.tsv";

fn qualified_plugin(sha256: &str) -> Option<(&'static str, &'static str)> {
    match sha256 {
        "500a534c403bcbcc15753718768f0a1f352cc6b93f9119e5eba25745f7205b9a" => {
            Some(("113", "Ensembl get_matched_variant_alleles normalization"))
        }
        "d08ce7d8c2e6b3229cef638beff145e465e224f0e3fcd43fd49d501ebc2740b1" => Some((
            "116",
            "Exact position/REF/ALT after VF indel conversion and VF-strand orientation",
        )),
        _ => None,
    }
}
fn chroms() -> Vec<String> {
    (1..=22)
        .map(|n| format!("chr{n}"))
        .chain(["chrX".into(), "chrY".into()])
        .collect()
}
fn chrom_index(chrom: &str) -> usize {
    chroms()
        .iter()
        .position(|c| c == chrom)
        .unwrap_or(usize::MAX)
}
fn published(kind: &str, chrom: &str) -> bool {
    matches!(kind, "snv" | "indel")
        && chrom_index(chrom) != usize::MAX
        && (kind, chrom) != ("indel", "chrY")
}

fn error(category: &'static str, message: impl Into<String>) -> Box<dyn std::error::Error> {
    builder_error::category(category, invalid(message))
}
fn require(condition: bool, message: impl Into<String>) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(error("ValueError", message))
    }
}
fn key_error(key: &str) -> Box<dyn std::error::Error> {
    error("KeyError", repr(key))
}
fn say(line: &str) {
    println!("{line}");
    io::stdout().flush().ok();
}
fn text(value: &str) -> Json {
    Json::Str(value.to_owned())
}
fn object(items: Vec<(&str, Json)>) -> Json {
    Json::Object(items.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}

// ---------------------------------------------------------------- files

fn digest(path: &str) -> Result<String> {
    crate::installer_config_paths::files::sha256_file(Path::new(path))
}
type Identity = (u64, u64, u64, i128, i128);
fn identity(path: &str) -> Result<Identity> {
    let value = metadata(Path::new(path))?;
    Ok((
        value.dev(),
        value.ino(),
        value.len(),
        value.mtime() as i128 * 1_000_000_000 + value.mtime_nsec() as i128,
        value.ctime() as i128 * 1_000_000_000 + value.ctime_nsec() as i128,
    ))
}
fn partial(path: &str) -> bool {
    path.split('/')
        .any(|part| part.to_lowercase().contains(".partial"))
}

/// `complete_file(path)`: the resolved, complete, nonempty regular input.
fn complete_file(path: &str) -> Result<String> {
    require(
        path.starts_with('/'),
        format!("An absolute input path is required: {path}"),
    )?;
    require(!partial(path), format!("Incomplete input: {path}"))?;
    let resolved = crate::python_io::resolve_strict(&native(Path::new(path))?)
        .map_err(|(cause, at)| os_error(cause, &at, None))?
        .display()
        .to_string();
    require(
        !partial(&resolved),
        format!("Incomplete input target: {path}"),
    )?;
    require(
        is_file(Path::new(&resolved))? && metadata(Path::new(&resolved))?.len() > 0,
        format!("Missing/empty input: {path}"),
    )?;
    Ok(resolved)
}

/// `verified_file(item, identities)`: the pinned size and SHA-256.
fn verified_file(item: &Json, identities: &mut Vec<(String, Identity)>) -> Result<String> {
    let Json::Object(_) = item else {
        return Err(error("TypeError", "manifest entry is not a mapping"));
    };
    let raw = match item.get("path") {
        Some(Json::Str(path)) => path.clone(),
        Some(_) => {
            return Err(error(
                "TypeError",
                "expected str, bytes or os.PathLike object",
            ));
        }
        None => return Err(key_error("path")),
    };
    let path = complete_file(&crate::installer_prerequisites::paths::pathlib(&raw))?;
    let bytes = match item.get("bytes") {
        Some(Json::Int(value)) if *value > 0 => *value,
        _ => return Err(error("ValueError", "Expected positive input bytes")),
    };
    let sha256 = item
        .get("sha256")
        .and_then(Json::as_str)
        .unwrap_or_default();
    require(
        sha256.len() == 64
            && sha256
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b)),
        "Expected pinned SHA256",
    )?;
    let before = identity(&path)?;
    require(
        before.2 as i128 == bytes,
        format!("Input size mismatch: {path}"),
    )?;
    say(&format!(
        "Verifying {}",
        path.rsplit('/').next().unwrap_or_default()
    ));
    require(
        digest(&path)? == sha256,
        format!("Input SHA256 mismatch: {path}"),
    )?;
    require(
        identity(&path)? == before,
        format!("Input changed during verification: {path}"),
    )?;
    record(identities, &path, before);
    Ok(path)
}
fn record(identities: &mut Vec<(String, Identity)>, path: &str, value: Identity) {
    match identities.iter_mut().find(|(seen, _)| seen == path) {
        Some(slot) => slot.1 = value,
        None => identities.push((path.to_owned(), value)),
    }
}

/// Lines of `open(path, encoding="utf-8", newline="")` or its gzip text
/// equivalent: split at \n, \r and \r\n, terminators kept.
struct Lines {
    reader: BufReader<Box<dyn Read>>,
    path: String,
}
impl Lines {
    fn open(path: &str) -> Result<Self> {
        let file = File::open(native(Path::new(path))?)
            .map_err(|cause| os_error(cause, Path::new(path), None))?;
        let reader: Box<dyn Read> = if path.ends_with(".gz") {
            Box::new(PythonGzipReader::new(file))
        } else {
            Box::new(file)
        };
        Ok(Self {
            reader: BufReader::with_capacity(1 << 20, reader),
            path: path.to_owned(),
        })
    }
    fn next(&mut self) -> Result<Option<String>> {
        let mut line = Vec::new();
        loop {
            let available = self
                .reader
                .fill_buf()
                .map_err(|cause| os_error(cause, Path::new(&self.path), None))?;
            if available.is_empty() {
                break;
            }
            match available.iter().position(|b| matches!(b, b'\n' | b'\r')) {
                Some(at) => {
                    let carriage = available[at] == b'\r';
                    line.extend_from_slice(&available[..=at]);
                    self.reader.consume(at + 1);
                    if carriage {
                        let next = self
                            .reader
                            .fill_buf()
                            .map_err(|cause| os_error(cause, Path::new(&self.path), None))?;
                        if next.first() == Some(&b'\n') {
                            line.push(b'\n');
                            self.reader.consume(1);
                        }
                    }
                    break;
                }
                None => {
                    let length = available.len();
                    line.extend_from_slice(available);
                    self.reader.consume(length);
                }
            }
        }
        if line.is_empty() {
            return Ok(None);
        }
        String::from_utf8(line)
            .map(Some)
            .map_err(|failure| error("UnicodeDecodeError", failure.to_string()))
    }
}

/// `re.findall(r'(\w+) "([^"\n]*)";', text)`
fn gtf_attributes(text: &str) -> Vec<(String, String)> {
    let chars: Vec<char> = text.chars().collect();
    let word = |c: char| c.is_alphanumeric() || c == '_';
    let mut found = Vec::new();
    let mut at = 0;
    while at < chars.len() {
        if !word(chars[at]) {
            at += 1;
            continue;
        }
        let start = at;
        let mut end = at;
        while end < chars.len() && word(chars[end]) {
            end += 1;
        }
        // A match starting inside a word run ends where the run ends, so only
        // the run's first character can start one.
        let mut matched = None;
        if chars.get(end) == Some(&' ') && chars.get(end + 1) == Some(&'"') {
            let mut close = end + 2;
            while close < chars.len() && chars[close] != '"' && chars[close] != '\n' {
                close += 1;
            }
            if chars.get(close) == Some(&'"') && chars.get(close + 1) == Some(&';') {
                matched = Some(close);
            }
        }
        match matched {
            Some(close) => {
                found.push((
                    chars[start..end].iter().collect(),
                    chars[end + 2..close].iter().collect(),
                ));
                at = close + 2;
            }
            None => at = end,
        }
    }
    found
}

/// `urllib.parse.unquote(value)`: percent-decoded UTF-8, invalid bytes replaced.
fn unquote(value: &str) -> String {
    if !value.contains('%') {
        return value.to_owned();
    }
    let bytes = value.as_bytes();
    let mut out = Vec::new();
    let mut at = 0;
    while at < bytes.len() {
        if bytes[at] == b'%' && at + 2 < bytes.len() {
            let hex = std::str::from_utf8(&bytes[at + 1..at + 3]).ok();
            if let Some(value) = hex.and_then(|h| u8::from_str_radix(h, 16).ok()) {
                out.push(value);
                at += 3;
                continue;
            }
        }
        out.push(bytes[at]);
        at += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

/// `re.sub(r"\.[0-9]+(_PAR_Y)?$", "", gene_id)`
fn unversioned(gene_id: &str) -> String {
    let body = gene_id.strip_suffix("_PAR_Y").unwrap_or(gene_id);
    if let Some(dot) = body.rfind('.') {
        let digits = &body[dot + 1..];
        if !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()) {
            return body[..dot].to_owned();
        }
    }
    if body.len() != gene_id.len() {
        // `_PAR_Y` without a version stays.
        return gene_id.to_owned();
    }
    gene_id.to_owned()
}

fn fields(line: &str) -> Vec<&str> {
    line.trim_end_matches(['\r', '\n']).split('\t').collect()
}

type ReleaseNames = HashMap<(String, String), Vec<(String, String)>>;

/// `vep_release_names(annotation, ensembl_release)`
fn vep_release_names(annotation: &str, release: &str) -> Result<ReleaseNames> {
    let gff3 = annotation.ends_with(".gff3") || annotation.ends_with(".gff3.gz");
    require(
        gff3 || annotation.ends_with(".gtf") || annotation.ends_with(".gtf.gz"),
        "The VEP-release gene annotation must be a GTF or GFF3",
    )?;
    let mut names: ReleaseNames = HashMap::new();
    let mut description = String::new();
    let mut lines = Lines::open(annotation)?;
    let mut number = 0;
    while let Some(line) = lines.next()? {
        number += 1;
        if line.starts_with('#') {
            description.push_str(&line);
            continue;
        }
        let row = fields(&line);
        require(
            row.len() == 9,
            format!("Invalid VEP-release annotation row {number}"),
        )?;
        if row[2] != "gene" {
            continue;
        }
        let mut attributes: HashMap<String, String> = HashMap::new();
        if gff3 {
            for item in row[8].split(';').filter(|item| !item.is_empty()) {
                let Some((key, value)) = item.split_once('=') else {
                    return Err(error(
                        "ValueError",
                        "not enough values to unpack (expected 2, got 1)",
                    ));
                };
                attributes.insert(key.to_owned(), unquote(value));
            }
        } else {
            attributes.extend(gtf_attributes(row[8]));
        }
        let gene_id = attributes.get("gene_id").filter(|v| !v.is_empty());
        let symbol = attributes.get("gene_name").filter(|v| !v.is_empty());
        let (Some(gene_id), Some(symbol)) = (gene_id, symbol) else {
            return Err(error(
                "ValueError",
                format!("VEP-release annotation gene lacks ID/name at {number}"),
            ));
        };
        let set = names
            .entry((unversioned(gene_id), row[0].to_owned()))
            .or_default();
        let value = (row[6].to_owned(), symbol.clone());
        if !set.contains(&value) {
            set.push(value);
        }
    }
    require(
        description.contains("GRCh38") && description.contains(&format!("(Ensembl {release})")),
        format!(
            "VEP gene names require the GRCh38 GENCODE annotation of Ensembl release {release}"
        ),
    )?;
    Ok(names)
}

/// `int(text)` for a GTF or region coordinate.
fn integer(value: &str) -> Result<i64> {
    crate::python_text::python_int(value)
        .and_then(|v| i64::try_from(v).ok())
        .ok_or_else(|| {
            error(
                "ValueError",
                format!("invalid literal for int() with base 10: {}", repr(value)),
            )
        })
}
fn unsafe_symbol(name: &str) -> bool {
    name.chars()
        .any(|c| is_python_space(c) || matches!(c, '|' | ';' | ',' | '='))
}
/// csv.writer(delimiter="\t", lineterminator="\n") quoting.
fn csv_field(value: &str) -> String {
    if value.contains(['\t', '"', '\n', '\r']) {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_owned()
    }
}

type ContextRow = (String, String, i64, i64, String, String);

/// `gene_context(regions, gtf, destination, release_names)`
fn gene_context(
    regions: &str,
    gtf: &str,
    destination: &str,
    release_names: Option<&ReleaseNames>,
) -> Result<Json> {
    let mut genes: HashMap<String, (String, i64, i64, String, String, String)> = HashMap::new();
    let mut stable_ids: HashMap<String, Vec<String>> = HashMap::new();
    let mut description = String::new();
    let mut lines = Lines::open(gtf)?;
    let mut number = 0;
    while let Some(line) = lines.next()? {
        number += 1;
        if line.starts_with('#') {
            description.push_str(&line);
            continue;
        }
        let row = fields(&line);
        require(row.len() == 9, format!("Invalid GTF row {number}"))?;
        if row[2] != "gene" {
            continue;
        }
        let attributes = gtf_attributes(row[8]);
        for key in ["gene_id", "gene_name", "gene_type"] {
            require(
                attributes.iter().filter(|(k, _)| k == key).count() == 1,
                format!("Missing/duplicate GTF {key} at {number}"),
            )?;
        }
        let value = |key: &str| {
            attributes
                .iter()
                .rev()
                .find(|(k, _)| k == key)
                .map(|(_, v)| v.clone())
        };
        let (gene_id, symbol) = (
            value("gene_id").unwrap_or_default(),
            value("gene_name").unwrap_or_default(),
        );
        require(
            !gene_id.is_empty() && !symbol.is_empty(),
            format!("GTF gene lacks ID/name at {number}"),
        )?;
        require(
            !genes.contains_key(&gene_id),
            format!("Duplicate GTF gene ID: {gene_id}"),
        )?;
        genes.insert(
            gene_id.clone(),
            (
                row[0].to_owned(),
                integer(row[3])?,
                integer(row[4])?,
                row[6].to_owned(),
                symbol,
                value("gene_type").unwrap_or_default(),
            ),
        );
        if let Some(stable) = stable_id(&gene_id) {
            stable_ids.entry(stable).or_default().push(gene_id);
        }
    }
    require(
        description.contains("GRCh38") && description.contains("version 48"),
        "Gene mapping requires the original GRCh38 GENCODE version 48 GTF",
    )?;
    let mut rows: Vec<ContextRow> = Vec::new();
    let mut conventions: Vec<&str> = Vec::new();
    let (mut seen, mut resolved_seen) = (HashSet::new(), HashSet::new());
    let mut symbols: Vec<((String, String), String, usize)> = Vec::new();
    let mut stable_matches = 0;
    let mut reader = Lines::open(regions)?;
    let header = reader.next()?.unwrap_or_default();
    require(
        fields(&header) == REGION_COLUMNS,
        "Unexpected publisher scored-region schema",
    )?;
    while let Some(line) = reader.next()? {
        if line.trim_end_matches(['\r', '\n']).is_empty() {
            continue;
        }
        let row = fields(&line);
        require(
            row.len() == REGION_COLUMNS.len(),
            "Malformed scored-region row",
        )?;
        let gene_id = row[0].to_owned();
        require(
            !seen.contains(&gene_id),
            format!("Duplicate scored-region gene: {gene_id}"),
        )?;
        seen.insert(gene_id.clone());
        let mut resolved_id = gene_id.clone();
        if !genes.contains_key(&gene_id) {
            let candidates = stable_ids.get(&gene_id).cloned().unwrap_or_default();
            require(
                candidates.len() == 1,
                format!(
                    "No exact or unique GENCODE48 stable gene ID: {gene_id}; aliases/ambiguous versions are not allowed"
                ),
            )?;
            resolved_id = candidates[0].clone();
            stable_matches += 1;
        }
        require(
            !resolved_seen.contains(&resolved_id),
            format!("Multiple source IDs resolve to one GENCODE gene: {resolved_id}"),
        )?;
        resolved_seen.insert(resolved_id.clone());
        let (chrom, start, end, strand, name, gene_type) = genes[&resolved_id].clone();
        require(
            chrom == row[1] && strand == row[4],
            format!("Gene contig/strand mismatch: {gene_id}"),
        )?;
        require(
            gene_type == "protein_coding",
            format!("Scored region is not a GENCODE48 protein-coding gene: {gene_id}"),
        )?;
        let (raw_start, raw_end) = (integer(row[2])?, integer(row[3])?);
        require(
            raw_end == end && (raw_start == start || raw_start == start - 1),
            format!("Scored-region/GTF bounds differ: {gene_id}"),
        )?;
        let convention = if raw_start == start {
            "1-based-closed"
        } else {
            "0-based-half-open"
        };
        if !conventions.contains(&convention) {
            conventions.push(convention);
        }
        require(
            chrom_index(&chrom) != usize::MAX
                && (strand == "+" || strand == "-")
                && 0 < start
                && start <= end,
            format!("Unsupported gene context: {gene_id}"),
        )?;
        require(
            !unsafe_symbol(&name),
            format!("Gene symbol cannot be represented safely in the plugin: {gene_id}"),
        )?;
        let key = (chrom.clone(), name.clone());
        match symbols.iter_mut().find(|(k, _, _)| *k == key) {
            Some((_, prior_strand, count)) => {
                require(
                    *prior_strand == strand,
                    format!(
                        "Opposite-strand genes share a symbol on {chrom}: {name}; plugin matching would be ambiguous"
                    ),
                )?;
                *count += 1;
            }
            None => symbols.push((key, strand.clone(), 1)),
        }
        rows.push((resolved_id, chrom, start, end, strand, name));
    }
    require(
        !rows.is_empty() && conventions.len() == 1,
        "Empty or mixed-coordinate scored-region table",
    )?;
    let genes_scored = rows.len();
    let mut release = Json::Null;
    if let Some(release_names) = release_names {
        // A second row with the same gene ID, interval and strand carries the VEP
        // release's name. A name another scored gene uses on that chromosome is
        // skipped, so no score moves to a different gene's transcripts.
        let (mut same, mut not_found, mut strand_differs) = (0, 0, 0);
        let (mut renamed, mut skipped) = (Vec::new(), Vec::new());
        let held: HashSet<(String, String)> =
            rows.iter().map(|r| (r.1.clone(), r.5.clone())).collect();
        let mut proposed: BTreeMap<(String, String), Vec<ContextRow>> = BTreeMap::new();
        for row in &rows {
            let found = release_names.get(&(unversioned(&row.0), row.1.clone()));
            let Some(found) = found.filter(|found| found.len() == 1) else {
                not_found += 1;
                continue;
            };
            let (strand, name) = &found[0];
            if *strand != row.4 {
                strand_differs += 1;
            } else if *name == row.5 {
                same += 1;
            } else {
                proposed
                    .entry((row.1.clone(), name.clone()))
                    .or_default()
                    .push(row.clone());
            }
        }
        let triple = |owner: &ContextRow, name: &str| {
            Json::Array(vec![text(&owner.0), text(&owner.5), text(name)])
        };
        let mut added = Vec::new();
        for ((chrom, name), owners) in proposed {
            if owners.len() > 1
                || held.contains(&(chrom.clone(), name.clone()))
                || unsafe_symbol(&name)
            {
                skipped.extend(owners.iter().map(|owner| triple(owner, &name)));
                continue;
            }
            let owner = &owners[0];
            added.push((
                owner.0.clone(),
                owner.1.clone(),
                owner.2,
                owner.3,
                owner.4.clone(),
                name.clone(),
            ));
            renamed.push(triple(owner, &name));
        }
        rows.extend(added);
        release = object(vec![
            ("same_name", Json::Int(same)),
            ("not_found", Json::Int(not_found)),
            ("strand_differs", Json::Int(strand_differs)),
            ("renamed", Json::Array(renamed)),
            ("skipped_name_conflicts", Json::Array(skipped)),
        ]);
    }
    rows.sort_by(|a, b| {
        (chrom_index(&a.1), a.2, a.3, &a.0, &a.5).cmp(&(chrom_index(&b.1), b.2, b.3, &b.0, &b.5))
    });
    let mut body = String::from("gene_id\tchrom\tgene_start\tgene_end\tstrand\tgene_name\n");
    for row in &rows {
        body.push_str(
            &[
                csv_field(&row.0),
                csv_field(&row.1),
                row.2.to_string(),
                row.3.to_string(),
                csv_field(&row.4),
                csv_field(&row.5),
            ]
            .join("\t"),
        );
        body.push('\n');
    }
    let mut output = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(native(Path::new(destination))?)
        .map_err(|cause| os_error(cause, Path::new(destination), None))?;
    output
        .write_all(body.as_bytes())
        .map_err(|cause| os_error(cause, Path::new(destination), None))?;
    Ok(object(vec![
        ("genes", Json::Int(genes_scored as i128)),
        ("context_rows", Json::Int(rows.len() as i128)),
        ("vep_release_names", release),
        ("source_coordinates", text(conventions[0])),
        ("context_coordinates", text("1-based-closed")),
        (
            "gene_id_matching",
            text(if stable_matches > 0 {
                "unique-GENCODE48-stable-ID-with-exact-context"
            } else {
                "exact"
            }),
        ),
        ("stable_id_matches", Json::Int(stable_matches)),
        ("versioned_gene_ids_retained", Json::Bool(true)),
        (
            "same_strand_symbol_groups",
            Json::Int(symbols.iter().filter(|(_, _, n)| *n > 1).count() as i128),
        ),
        (
            "symbol_projection",
            text(
                "One unchanged prediction per distinct symbol per native row, under the GENCODE48 name and, when it differs, the VEP release's name of the same gene ID; all same-strand gene IDs/intervals retained in context. Opposite-strand same-symbol contexts rejected.",
            ),
        ),
        (
            "join",
            text("source variant position and model strand within verified scored-gene intervals"),
        ),
        (
            "scope",
            text(
                "Adapter-derived gene association, not a publisher gene-conditioned prediction. Predicted sites may lie outside the associated gene; no effects are filtered.",
            ),
        ),
    ]))
}

/// `re.fullmatch(r"(ENSG[0-9]{11})\.[0-9]+(_PAR_Y)?", gene_id)` groups joined.
fn stable_id(gene_id: &str) -> Option<String> {
    let (body, par) = match gene_id.strip_suffix("_PAR_Y") {
        Some(body) => (body, "_PAR_Y"),
        None => (gene_id, ""),
    };
    let (stem, version) = body.split_once('.')?;
    let valid = stem.len() == 15
        && stem.starts_with("ENSG")
        && stem[4..].bytes().all(|b| b.is_ascii_digit())
        && !version.is_empty()
        && version.bytes().all(|b| b.is_ascii_digit());
    valid.then(|| format!("{stem}{par}"))
}

// ---------------------------------------------------------------- tools

/// `shutil.which(name)`
fn which(name: &str) -> Option<String> {
    let executable = |candidate: &Path| {
        CString::new(candidate.as_os_str().as_encoded_bytes())
            .is_ok_and(|t| unsafe { libc::access(t.as_ptr(), libc::X_OK) } == 0)
            && candidate.exists()
            && !candidate.is_dir()
    };
    if name.contains('/') {
        return executable(Path::new(name)).then(|| name.to_owned());
    }
    let path = std::env::var_os("PATH").unwrap_or_else(|| "/bin:/usr/bin".into());
    std::env::split_paths(&path)
        .map(|folder| folder.join(name))
        .find(|candidate| executable(candidate))
        .map(|candidate| candidate.display().to_string())
}

fn called_process_error(
    arguments: &[String],
    status: std::process::ExitStatus,
) -> Box<dyn std::error::Error> {
    let shown = arguments
        .iter()
        .map(|a| repr(a))
        .collect::<Vec<_>>()
        .join(", ");
    let outcome = match (status.code(), status.signal()) {
        (Some(code), _) => format!("returned non-zero exit status {code}"),
        (None, Some(signal)) => format!("died with {}", signal_repr(signal)),
        (None, None) => format!("ended with {status}"),
    };
    error(
        "CalledProcessError",
        format!("Command '[{shown}]' {outcome}."),
    )
}

/// `subprocess.check_output(arguments, text=True)` (stderr as given).
fn check_output(arguments: &[String], merge_stderr: bool) -> Result<String> {
    let mut command = Command::new(&arguments[0]);
    command.args(&arguments[1..]).stdout(Stdio::piped());
    if merge_stderr {
        command.stderr(Stdio::piped());
    }
    let finished = command
        .output()
        .map_err(|cause| os_error(cause, Path::new(&arguments[0]), None))?;
    if !finished.status.success() {
        return Err(called_process_error(arguments, finished.status));
    }
    let mut bytes = finished.stdout;
    if merge_stderr {
        bytes.extend(finished.stderr);
    }
    let text = String::from_utf8(bytes).map_err(|f| error("UnicodeDecodeError", f.to_string()))?;
    Ok(text.replace("\r\n", "\n").replace('\r', "\n"))
}
fn run(arguments: &[String]) -> Result<()> {
    let status = Command::new(&arguments[0])
        .args(&arguments[1..])
        .status()
        .map_err(|cause| os_error(cause, Path::new(&arguments[0]), None))?;
    if status.success() {
        Ok(())
    } else {
        Err(called_process_error(arguments, status))
    }
}

/// Children of one conversion, ended as the script's `terminate` does.
struct Processes(Vec<Child>);
impl Drop for Processes {
    fn drop(&mut self) {
        for child in &mut self.0 {
            if matches!(child.try_wait(), Ok(None)) {
                unsafe { libc::kill(child.id() as i32, libc::SIGTERM) };
            }
        }
        for child in &mut self.0 {
            let _ = child.wait();
        }
    }
}

struct Entry {
    chrom: String,
    path: String,
    resolved: String,
}

/// One source shard, `gzip -dc | projection` into `sink`, and its checked
/// conversion report. Shards after the first omit the VCF header.
fn project(
    kind: &str,
    item: &Entry,
    header: bool,
    context: &str,
    converter: &[String],
    sink: Stdio,
    tools: &HashMap<&str, String>,
) -> Result<Json> {
    let mut processes = Processes(Vec::new());
    let mut source = Command::new(&tools["gzip"])
        .args(["-dc", &item.resolved])
        .stdout(Stdio::piped())
        .spawn()
        .map_err(|cause| os_error(cause, Path::new(&tools["gzip"]), None))?;
    let source_output = source.stdout.take().expect("piped stdout");
    let mut arguments = converter[1..].to_vec();
    arguments.extend(
        [
            "--genes",
            context,
            "--chrom",
            &item.chrom,
            "--kind",
            kind,
            "--input",
            "-",
            "--output",
            "-",
        ]
        .map(str::to_owned),
    );
    if !header {
        arguments.push("--no-header".into());
    }
    let projection = Command::new(&converter[0])
        .args(&arguments)
        .stdin(Stdio::from(source_output))
        .stdout(sink)
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|cause| os_error(cause, Path::new(&converter[0]), None))?;
    processes.0.push(source);
    let finished = projection
        .wait_with_output()
        .map_err(|cause| os_error(cause, Path::new(&converter[0]), None))?;
    let diagnostic = String::from_utf8_lossy(&finished.stderr).into_owned();
    if !finished.status.success() {
        let tail: String = diagnostic
            .chars()
            .rev()
            .take(4096)
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
            .collect();
        return Err(error(
            "ValueError",
            format!("Rust projection failed: {tail}"),
        ));
    }
    let source = processes.0.last_mut().expect("source process");
    let decompressed = source
        .wait()
        .map_err(|cause| os_error(cause, Path::new(&tools["gzip"]), None))?;
    require(
        decompressed.success(),
        format!("Source decompression failed: {}", item.path),
    )?;
    let report = json_text::parse(&diagnostic, false).map_err(|m| error("JSONDecodeError", m))?;
    let field = |key: &str| report.get(key).cloned().ok_or_else(|| key_error(key));
    require(
        matches!(field("module")?, Json::Str(ref v) if v == "spliceai2-vcf")
            && matches!(field("chrom")?, Json::Str(ref v) if *v == item.chrom)
            && matches!(field("kind")?, Json::Str(ref v) if v == kind),
        "Invalid Rust conversion report",
    )?;
    let count = |key: &str| -> Result<i128> {
        match field(key)? {
            Json::Int(value) => Ok(value),
            _ => Err(error("TypeError", "'>' not supported")),
        }
    };
    require(
        count("input_rows")? > 0 && count("output_records")? >= count("input_rows")?,
        "Empty/incomplete projection",
    )?;
    Ok(report)
}

/// A failure carried between threads: its category and message.
type Carried<T> = std::result::Result<T, (&'static str, String)>;
fn carry<T>(result: Result<T>) -> Carried<T> {
    builder_error::finish(result).map_err(|cause| {
        match cause.downcast_ref::<builder_error::RuntimeError>() {
            Some(failure) => (failure.category_name(), failure.message()),
            None => ("RuntimeError", cause.to_string()),
        }
    })
}
fn uncarry<T>(result: Carried<T>) -> Result<T> {
    result.map_err(|(category, message)| error(category, message))
}

/// Run `task` on each index below `count` with at most `jobs` threads; the
/// results come back in index order. After a failure no new index starts.
fn in_parallel<T: Send>(
    count: usize,
    jobs: usize,
    task: impl Fn(usize) -> Result<T> + Sync,
) -> Result<Vec<T>> {
    let next = AtomicUsize::new(0);
    let failed = AtomicBool::new(false);
    let slots: Mutex<Vec<Option<Carried<T>>>> = Mutex::new((0..count).map(|_| None).collect());
    std::thread::scope(|scope| {
        for _ in 0..jobs.min(count) {
            scope.spawn(|| {
                while !failed.load(Ordering::SeqCst) {
                    let index = next.fetch_add(1, Ordering::SeqCst);
                    if index >= count {
                        break;
                    }
                    let result = carry(task(index));
                    if result.is_err() {
                        failed.store(true, Ordering::SeqCst);
                    }
                    slots.lock().expect("result slots")[index] = Some(result);
                }
            });
        }
    });
    // Indexes start in order, so the first failure in index order comes
    // before any index that never started.
    let mut results = Vec::with_capacity(count);
    for slot in slots.into_inner().expect("result slots") {
        match slot {
            Some(result) => results.push(uncarry(result)?),
            None => break,
        }
    }
    require(results.len() == count, "Parallel work stopped early")?;
    Ok(results)
}

/// The 28-byte empty BGZF block that ends every bgzip output.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// `--jobs` above 1 (Rust only): each shard is projected and compressed into
/// its own BGZF part at the same time, and the parts are joined in
/// chromosome order without the end blocks of all but the last part. The
/// decompressed VCF is the same as the sequential pipeline's; the compressed
/// bytes are not, because the block boundaries differ.
fn convert_parts(
    kind: &str,
    entries: &[Entry],
    context: &str,
    target: &str,
    converter: &[String],
    threads: i64,
    jobs: usize,
    tools: &HashMap<&str, String>,
) -> Result<Vec<Json>> {
    let output = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(native(Path::new(target))?)
        .map_err(|cause| os_error(cause, Path::new(target), None))?;
    let bgzip = tools["bgzip"].clone();
    let parts: Vec<String> = (0..entries.len())
        .map(|index| format!("{target}.part{index}"))
        .collect();
    let statistics = in_parallel(entries.len(), jobs, |index| {
        let item = &entries[index];
        say(&format!("Projecting {kind} {}", item.chrom));
        let part = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(native(Path::new(&parts[index]))?)
            .map_err(|cause| os_error(cause, Path::new(&parts[index]), None))?;
        let mut processes = Processes(Vec::new());
        let mut compressor = Command::new(&bgzip)
            .args(["-c", "-@", &threads.to_string()])
            .stdin(Stdio::piped())
            .stdout(part)
            .spawn()
            .map_err(|cause| os_error(cause, Path::new(&bgzip), None))?;
        let sink = Stdio::from(compressor.stdin.take().expect("piped stdin"));
        processes.0.push(compressor);
        let report = project(kind, item, index == 0, context, converter, sink, tools)?;
        let status = processes.0[0]
            .wait()
            .map_err(|cause| os_error(cause, Path::new(&bgzip), None))?;
        require(status.success(), "BGZF compression failed")?;
        Ok(report)
    });
    let statistics = match statistics {
        Ok(statistics) => statistics,
        Err(failure) => {
            for part in &parts {
                let _ = fs::remove_file(part);
            }
            return Err(failure);
        }
    };
    let mut output = io::BufWriter::new(output);
    for (index, part) in parts.iter().enumerate() {
        let mut bytes = fs::read(part).map_err(|cause| os_error(cause, Path::new(part), None))?;
        if index + 1 < parts.len() {
            require(
                bytes.ends_with(&BGZF_EOF),
                format!("BGZF part lacks its end block: {part}"),
            )?;
            bytes.truncate(bytes.len() - BGZF_EOF.len());
        }
        output
            .write_all(&bytes)
            .map_err(|cause| os_error(cause, Path::new(target), None))?;
        fs::remove_file(part).map_err(|cause| os_error(cause, Path::new(part), None))?;
    }
    output
        .into_inner()
        .map_err(|cause| os_error(cause.into_error(), Path::new(target), None))?
        .sync_all()
        .map_err(|cause| os_error(cause, Path::new(target), None))?;
    Ok(statistics)
}

/// `convert_kind(kind, entries, context, target, converter, threads, tools)`
fn convert_kind(
    kind: &str,
    entries: &[Entry],
    context: &str,
    target: &str,
    converter: &[String],
    threads: i64,
    jobs: usize,
    tools: &HashMap<&str, String>,
) -> Result<Vec<Json>> {
    let statistics = if jobs > 1 {
        convert_parts(
            kind, entries, context, target, converter, threads, jobs, tools,
        )?
    } else {
        convert_sequentially(kind, entries, context, target, converter, threads, tools)?
    };
    index_and_probe(entries, target, tools)?;
    Ok(statistics)
}

/// The script's single `gzip | projection | bgzip` pipeline per kind.
fn convert_sequentially(
    kind: &str,
    entries: &[Entry],
    context: &str,
    target: &str,
    converter: &[String],
    threads: i64,
    tools: &HashMap<&str, String>,
) -> Result<Vec<Json>> {
    let mut statistics = Vec::new();
    let mut processes = Processes(Vec::new());
    let output = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(native(Path::new(target))?)
        .map_err(|cause| os_error(cause, Path::new(target), None))?;
    let bgzip = tools["bgzip"].clone();
    let mut compressor = Command::new(&bgzip)
        .args(["-c", "-@", &threads.to_string()])
        .stdin(Stdio::piped())
        .stdout(output)
        .spawn()
        .map_err(|cause| os_error(cause, Path::new(&bgzip), None))?;
    let compressor_input = compressor.stdin.take().expect("piped stdin");
    processes.0.push(compressor);
    for (index, item) in entries.iter().enumerate() {
        say(&format!("Projecting {kind} {}", item.chrom));
        let sink = compressor_input
            .as_fd()
            .try_clone_to_owned()
            .map_err(|cause| os_error(cause, Path::new(&bgzip), None))?;
        statistics.push(project(
            kind,
            item,
            index == 0,
            context,
            converter,
            Stdio::from(sink),
            tools,
        )?);
    }
    drop(compressor_input);
    let status = processes.0[0]
        .wait()
        .map_err(|cause| os_error(cause, Path::new(&bgzip), None))?;
    require(status.success(), "BGZF compression failed")?;
    Ok(statistics)
}

/// The tabix index, its contig list and the first-record probe.
fn index_and_probe(entries: &[Entry], target: &str, tools: &HashMap<&str, String>) -> Result<()> {
    let tabix = tools["tabix"].clone();
    run(&[tabix.clone(), "-p".into(), "vcf".into(), target.into()])?;
    let contigs: Vec<String> = check_output(&[tabix.clone(), "-l".into(), target.into()], false)?
        .lines()
        .map(str::to_owned)
        .collect();
    require(
        contigs == entries.iter().map(|e| e.chrom.clone()).collect::<Vec<_>>(),
        "Indexed contigs differ from converted inputs",
    )?;
    // Probe the first observed record, not an arbitrary possibly empty interval.
    let mut lines = Lines::open(target)?;
    let mut first = None;
    while let Some(line) = lines.next()? {
        if !line.starts_with('#') {
            first = Some(line.trim_end_matches('\n').to_owned());
            break;
        }
    }
    let first = first.ok_or_else(|| error("StopIteration", ""))?;
    let row = first.split('\t').collect::<Vec<_>>();
    let region = format!(
        "{}:{}-{}",
        row[0],
        row.get(1).unwrap_or(&""),
        row.get(1).unwrap_or(&"")
    );
    let query = check_output(&[tabix, target.into(), region], false)?;
    require(
        query.lines().any(|line| line == first),
        "Nonempty VCF/index query did not recover the first record",
    )?;
    Ok(())
}

// ---------------------------------------------------------------- prepare

struct Options {
    manifest: String,
    assembly: String,
    output_dir: String,
    receipt: String,
    converter_bin: Option<String>,
    threads: i64,
    jobs: usize,
    validation_subset: bool,
}

fn no_symlink_ancestors(path: &str) -> Result<bool> {
    let mut current = PathBuf::from(path);
    loop {
        if is_link(&current)? {
            return Ok(false);
        }
        if !current.pop() {
            return Ok(true);
        }
    }
}

fn prepare(options: &Options) -> Result<Json> {
    require(
        matches!(options.assembly.as_str(), "GRCh38" | "hg38"),
        "SpliceAI2 route supports GRCh38/hg38 only; hg19 remains v1",
    )?;
    require(options.threads > 0, "threads must be positive")?;
    let manifest_path = complete_file(&options.manifest)?;
    let manifest_before = identity(&manifest_path)?;
    let manifest_text = crate::installer_config_paths::files::read_text(Path::new(&manifest_path))?;
    let manifest =
        json_text::parse(&manifest_text, false).map_err(|m| error("JSONDecodeError", m))?;
    let manifest_sha = digest(&manifest_path)?;
    let get = |key: &str| manifest.get(key).cloned().unwrap_or(Json::Null);
    let is_text = |key: &str, value: &str| matches!(get(key), Json::Str(ref v) if v == value);
    require(
        matches!(get("schema_version"), Json::Int(1))
            && is_text("status", "downloads_complete")
            && matches!(get("confirmed_by"), Json::Str(ref v) if !v.trim_matches(is_python_space).is_empty()),
        "A confirmed, complete source manifest is required; downloading inventories are not inputs",
    )?;
    require(
        is_text("resource", "SpliceAI2")
            && is_text("version", "2.0")
            && is_text("assembly", "GRCh38")
            && is_text("source_revision", REVISION),
        "Wrong source version/assembly/revision",
    )?;
    require(
        matches!(get("gencode_gtf").get("release"), Some(Json::Str(v)) if v == "GENCODE48"),
        "Pinned GENCODE48 GTF is required",
    )?;
    require(
        matches!(get("vep_gene_annotation"), Json::Object(_)),
        "A pinned GENCODE GTF/GFF3 of VEP's Ensembl release is required so renamed genes keep their scores",
    )?;
    let files = match manifest.get("files") {
        Some(Json::Array(items)) => items.clone(),
        None => Vec::new(),
        Some(_) => return Err(error("TypeError", "manifest files must be a list")),
    };
    let mut pairs = Vec::new();
    for item in &files {
        let kind = item
            .get("kind")
            .and_then(Json::as_str)
            .ok_or_else(|| key_error("kind"))?;
        let chrom = item
            .get("chrom")
            .and_then(Json::as_str)
            .ok_or_else(|| key_error("chrom"))?;
        pairs.push((kind.to_owned(), chrom.to_owned()));
    }
    let unique: HashSet<_> = pairs.iter().cloned().collect();
    require(
        unique.len() == pairs.len()
            && unique.iter().all(|(k, c)| published(k, c))
            && !pairs.is_empty(),
        "Duplicate/unsupported source shards",
    )?;
    require(
        unique.len() == 47 || options.validation_subset,
        "Full cache requires all 47 published shards; chrY indels are not published",
    )?;
    let kinds: HashSet<&str> = pairs.iter().map(|(k, _)| k.as_str()).collect();
    require(
        !options.validation_subset || kinds.len() == 2,
        "Subset validation requires both variant types",
    )?;
    let driver = std::env::current_exe()
        .map_err(|cause| os_error(cause, Path::new("priva"), None))?
        .display()
        .to_string();
    let mut identities = vec![
        (manifest_path.clone(), manifest_before),
        (driver.clone(), identity(&driver)?),
    ];
    let driver_sha = digest(&driver)?;
    for (item, (kind, chrom)) in files.iter().zip(&pairs) {
        let expected = format!("SpliceAI2_GRCh38_{kind}_{chrom}.tsv.gz");
        let path = item
            .get("path")
            .and_then(Json::as_str)
            .ok_or_else(|| key_error("path"))?;
        require(
            path.rsplit('/').next() == Some(expected.as_str()),
            format!("Unexpected source name: {path}"),
        )?;
        complete_file(&crate::installer_prerequisites::paths::pathlib(path))?;
    }
    let regions = verified_file(&get("scored_regions"), &mut identities)?;
    let gtf = verified_file(&get("gencode_gtf"), &mut identities)?;
    let vep_annotation = verified_file(&get("vep_gene_annotation"), &mut identities)?;
    let plugin_item = get("vep_spliceai_plugin");
    let plugin_path = verified_file(&plugin_item, &mut identities)?;
    let plugin_sha = plugin_item
        .get("sha256")
        .and_then(Json::as_str)
        .unwrap_or_default()
        .to_owned();
    let Some((plugin_release, allele_matching)) = qualified_plugin(&plugin_sha) else {
        return Err(error("ValueError", "Unqualified official plugin revision"));
    };
    let plugin = object(vec![
        ("path", text(&plugin_path)),
        ("sha256", text(&plugin_sha)),
        ("release", text(plugin_release)),
        ("allele_matching", text(allele_matching)),
    ]);
    let mut tools: HashMap<&str, String> = HashMap::new();
    let mut tool_identities = Vec::new();
    for name in ["gzip", "bgzip", "tabix"] {
        let found = which(name)
            .ok_or_else(|| error("ValueError", format!("Missing mature command tool: {name}")))?;
        let path = crate::python_io::resolve_strict(Path::new(&found))
            .map_err(|(cause, at)| os_error(cause, &at, None))?
            .display()
            .to_string();
        tools.insert(name, path.clone());
        let before = identity(&path)?;
        record(&mut identities, &path, before);
        let version = check_output(&[path.clone(), "--version".into()], true)?
            .lines()
            .next()
            .unwrap_or_default()
            .to_owned();
        tool_identities.push((
            name.to_owned(),
            object(vec![
                ("path", text(&path)),
                ("sha256", text(&digest(&path)?)),
                ("version", text(&version)),
            ]),
        ));
    }
    let folder = absolute(Path::new(&options.output_dir))?;
    let receipt_path = absolute(Path::new(&options.receipt))?;
    for path in [&folder, &receipt_path] {
        require(
            !path
                .split('/')
                .any(|part| part.to_lowercase().contains(".partial") || part == ".."),
            "Incomplete/noncanonical output path",
        )?;
        require(
            no_symlink_ancestors(path)?,
            "Output paths must not traverse symlinks",
        )?;
    }
    let receipt_parent = crate::installer_prerequisites::paths::parent(&receipt_path);
    let resolved = |path: &str| -> Result<String> {
        Ok(
            crate::installer_config_paths::files::resolve(Path::new(path))?
                .display()
                .to_string(),
        )
    };
    require(
        resolved(&receipt_parent)? == resolved(&folder)?,
        "Receipt must be directly inside the output directory",
    )?;
    let (converter, mut converter_info) = match &options.converter_bin {
        Some(name) => {
            let found = which(name).ok_or_else(|| {
                error(
                    "ValueError",
                    "The explicit --converter-bin installation tool is not executable",
                )
            })?;
            let path = complete_file(&absolute(Path::new(&found))?)?;
            (
                vec![path.clone(), "project".to_owned()],
                vec![("role".to_owned(), text("explicit_converter_override"))],
            )
        }
        None => (
            vec![driver.clone(), "spliceai2-vcf".to_owned()],
            vec![("role".to_owned(), text("priva_builtin"))],
        ),
    };
    let converter_identity = identity(&converter[0])?;
    record(&mut identities, &converter[0], converter_identity);
    converter_info.push(("path".into(), text(&converter[0])));
    converter_info.push(("sha256".into(), text(&digest(&converter[0])?)));
    mkdir(Path::new(&folder))?;
    let mut final_paths = Vec::new();
    for kind in ["snv", "indel"] {
        for suffix in ["", ".tbi"] {
            final_paths.push(format!(
                "{folder}/SpliceAI2_GRCh38.{kind}.vep.vcf.gz{suffix}"
            ));
        }
    }
    final_paths.push(format!("{folder}/{CONTEXT_NAME}"));
    final_paths.push(receipt_path.clone());
    let lock_path = format!("{folder}/.install.lock");
    let lock = OpenOptions::new()
        .append(true)
        .create(true)
        .open(native(Path::new(&lock_path))?)
        .map_err(|cause| os_error(cause, Path::new(&lock_path), None))?;
    let _guard =
        crate::file_lock::ExclusiveLock::acquire(lock, crate::file_lock::Wait::Nonblocking)
            .map_err(|cause| {
                error(
                    "BlockingIOError",
                    format!(
                        "[Errno {}] Resource temporarily unavailable",
                        cause.raw_os_error().unwrap_or(11)
                    ),
                )
            })?;
    for path in &final_paths {
        require(
            !crate::installer_config_paths::files::exists(Path::new(path))?
                && !is_link(Path::new(path))?,
            "Refusing to replace existing prepared outputs",
        )?;
    }
    let stage = Building::new(&folder)?;
    let context = format!("{}/{CONTEXT_NAME}", stage.0);
    let names = vep_release_names(&vep_annotation, plugin_release)?;
    let mapping = gene_context(&regions, &gtf, &context, Some(&names))?;
    // With --jobs above 1 the shards' SHA-256 checks run at the same time;
    // their identities are recorded in manifest order either way.
    let checked = in_parallel(files.len(), options.jobs, |index| {
        let mut seen = Vec::new();
        let resolved = verified_file(&files[index], &mut seen)?;
        Ok((resolved, seen))
    })?;
    let mut entries: Vec<(String, Entry)> = Vec::new();
    for ((item, (kind, chrom)), (resolved, seen)) in files.iter().zip(&pairs).zip(checked) {
        for (path, value) in seen {
            record(&mut identities, &path, value);
        }
        let path = item
            .get("path")
            .and_then(Json::as_str)
            .unwrap_or_default()
            .to_owned();
        entries.push((
            kind.clone(),
            Entry {
                chrom: chrom.clone(),
                path,
                resolved,
            },
        ));
    }
    let (mut outputs, mut statistics) = (Vec::new(), Vec::new());
    for kind in ["snv", "indel"] {
        let mut selected: Vec<&Entry> = entries
            .iter()
            .filter(|(k, _)| k == kind)
            .map(|(_, e)| e)
            .collect();
        selected.sort_by_key(|entry| chrom_index(&entry.chrom));
        let selected: Vec<Entry> = selected
            .into_iter()
            .map(|e| Entry {
                chrom: e.chrom.clone(),
                path: e.path.clone(),
                resolved: e.resolved.clone(),
            })
            .collect();
        let name = format!("SpliceAI2_GRCh38.{kind}.vep.vcf.gz");
        let target = format!("{}/{name}", stage.0);
        let report = convert_kind(
            kind,
            &selected,
            &context,
            &target,
            &converter,
            options.threads,
            options.jobs,
            &tools,
        )?;
        statistics.push((kind.to_owned(), Json::Array(report)));
        outputs.push((
            kind.to_owned(),
            object(vec![
                ("path", text(&format!("{folder}/{name}"))),
                ("index_path", text(&format!("{folder}/{name}.tbi"))),
                (
                    "bytes",
                    Json::Int(metadata(Path::new(&target))?.len() as i128),
                ),
                ("sha256", text(&digest(&target)?)),
                ("index_sha256", text(&digest(&format!("{target}.tbi"))?)),
            ]),
        ));
    }
    for (path, before) in &identities {
        require(
            identity(path)? == *before,
            format!("Input changed while preparing cache: {path}"),
        )?;
    }
    let source_files = Json::Array(files.clone());
    let status = if options.validation_subset {
        "validated_subset_annotation_only"
    } else {
        "prepared_annotation_only"
    };
    let mut context_info = vec![
        ("path".to_owned(), text(&format!("{folder}/{CONTEXT_NAME}"))),
        ("sha256".to_owned(), text(&digest(&context)?)),
    ];
    if let Json::Object(items) = mapping {
        context_info.extend(items);
    }
    let receipt = object(vec![
        ("schema_version", Json::Int(1)),
        ("resource", text("SpliceAI2")),
        ("version", text("2.0")),
        ("assembly", text("GRCh38")),
        ("status", text(status)),
        (
            "coverage",
            text(if options.validation_subset {
                "validation_subset"
            } else {
                "complete_published"
            }),
        ),
        ("scientific_qualified", Json::Bool(false)),
        ("source_revision", text(REVISION)),
        ("source_manifest_sha256", text(&manifest_sha)),
        ("outputs", Json::Object(outputs)),
        ("statistics", Json::Object(statistics)),
        ("source_files", source_files),
        ("scored_regions", get("scored_regions")),
        ("gencode_gtf", get("gencode_gtf")),
        ("vep_gene_annotation", get("vep_gene_annotation")),
        ("gene_context", Json::Object(context_info)),
        ("vep_spliceai_plugin", plugin),
        ("converter", Json::Object(converter_info)),
        ("tools", Json::Object(tool_identities)),
        ("driver_sha256", text(&driver_sha)),
        (
            "slurm_job_id",
            std::env::var("SLURM_JOB_ID").map_or(Json::Null, |v| text(&v)),
        ),
        (
            "thresholds",
            object(vec![
                ("high_recall", Json::Float(0.1)),
                ("balanced", Json::Float(0.25)),
                ("high_precision", Json::Float(0.5)),
            ]),
        ),
        (
            "projection",
            text(
                "Unchanged top0 site delta/distance strings only, one allele/gene per VCF record. No cutoff filtering, score recomputation, masking, normalization or strand flipping. Original 84-column tables retain all top3/junction/reference/alternate/summary values.",
            ),
        ),
        (
            "coverage_limits",
            text(
                "GENCODE48 protein-coding scored regions; population-observed indels; no published chrY indel shard; the plugin cannot score transcripts whose gene has no name in the VEP cache. Missing coverage is not zero score.",
            ),
        ),
        (
            "license",
            object(vec![
                (
                    "source",
                    text(&format!(
                        "https://huggingface.co/datasets/illumina-ai/SpliceAI2-data/tree/{REVISION}"
                    )),
                ),
                (
                    "note",
                    text(
                        "Original and derived score resources retain Illumina SpliceAI2 terms; no redistribution or public service is enabled by this adapter.",
                    ),
                ),
            ]),
        ),
        (
            "activation",
            text(
                "Preparation only. Legacy PriVA ACMG/ClinVar/SpliceVault and CADD v1 feature contracts are not qualified for substitution.",
            ),
        ),
    ]);
    let receipt_name = receipt_path
        .rsplit('/')
        .next()
        .unwrap_or_default()
        .to_owned();
    let staged_receipt = format!("{}/{receipt_name}", stage.0);
    require(
        !final_paths[..final_paths.len() - 1]
            .iter()
            .any(|path| path.rsplit('/').next() == Some(receipt_name.as_str())),
        "Receipt name collides with an output",
    )?;
    let body = json_text::dumps_with_options(
        &receipt,
        DumpsOptions {
            indent: Some(2),
            ..DumpsOptions::default()
        },
    )
    .map_err(|m| error("ValueError", m))?;
    fs::write(native(Path::new(&staged_receipt))?, format!("{body}\n"))
        .map_err(|cause| os_error(cause, Path::new(&staged_receipt), None))?;
    // Receipt is published last; incomplete sets never have a readiness receipt.
    for path in &final_paths[..final_paths.len() - 1] {
        let name = path.rsplit('/').next().unwrap_or_default();
        rename(Path::new(&format!("{}/{name}", stage.0)), Path::new(path))?;
    }
    rename(Path::new(&staged_receipt), Path::new(&receipt_path))?;
    drop(stage);
    let printed = object(vec![
        ("status", text(status)),
        ("receipt", text(&receipt_path)),
    ]);
    say(&json_text::dumps(&printed, false).map_err(|m| error("ValueError", m))?);
    Ok(receipt)
}

/// `tempfile.TemporaryDirectory(prefix=".building-", dir=folder)`
struct Building(String);
impl Building {
    fn new(folder: &str) -> Result<Self> {
        let template = format!("{folder}/.building-XXXXXX");
        let mut raw = CString::new(template.clone())
            .map_err(|_| error("ValueError", "embedded null byte"))?
            .into_bytes_with_nul();
        if unsafe { libc::mkdtemp(raw.as_mut_ptr().cast()) }.is_null() {
            return Err(os_error(
                io::Error::last_os_error(),
                Path::new(&template),
                None,
            ));
        }
        Ok(Self(
            String::from_utf8_lossy(&raw[..raw.len() - 1]).into_owned(),
        ))
    }
}
impl Drop for Building {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

// ---------------------------------------------------------------- command line

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(command) = args.first() else {
        return Err(invalid(format!(
            "the following arguments are required: command\n{HELP}"
        )));
    };
    if matches!(command.as_str(), "-h" | "--help") {
        println!("{HELP}");
        return Ok(());
    }
    let (names, flags, required): (&[&str], &[&str], &[&str]) = match command.as_str() {
        "prepare" => (
            &[
                "manifest",
                "assembly",
                "output-dir",
                "receipt",
                "converter-bin",
                "threads",
                "jobs",
            ],
            &["validation-subset"],
            &["manifest", "assembly", "output-dir", "receipt"],
        ),
        "gene-context" => (
            &[
                "regions",
                "gtf",
                "output",
                "vep-annotation",
                "ensembl-release",
            ],
            &[],
            &["regions", "gtf", "output"],
        ),
        _ => {
            return Err(invalid(format!(
                "argument command: invalid choice: {} (choose from 'prepare', 'gene-context')",
                repr(command)
            )));
        }
    };
    let Some(occurrences) =
        parse_long_option_occurrences(&args[1..], names, flags).map_err(invalid)?
    else {
        println!("{HELP}");
        return Ok(());
    };
    let mut values = HashMap::new();
    let mut set_flags = Vec::new();
    for (name, value) in occurrences {
        if flags.contains(&name.as_str()) {
            set_flags.push(name);
        } else {
            values.insert(name, value);
        }
    }
    require_options(&values, required).map_err(invalid)?;
    let path = |name: &str| {
        values
            .get(name)
            .map(|v| crate::installer_prerequisites::paths::pathlib(v))
    };
    if command == "gene-context" {
        let threads_free = || -> Result<Json> {
            let annotation = path("vep-annotation");
            let release = values.get("ensembl-release");
            require(
                annotation.is_none() == release.is_none(),
                "--vep-annotation and --ensembl-release go together",
            )?;
            let names = match (annotation, release) {
                (Some(annotation), Some(release)) => Some(vep_release_names(&annotation, release)?),
                _ => None,
            };
            gene_context(
                &path("regions").unwrap_or_default(),
                &path("gtf").unwrap_or_default(),
                &path("output").unwrap_or_default(),
                names.as_ref(),
            )
        };
        let result = builder_error::finish(threads_free())?;
        say(&json_text::dumps(&result, false).map_err(|m| error("ValueError", m))?);
        return Ok(());
    }
    let threads = match values.get("threads") {
        None => 1,
        Some(value) => crate::python_text::python_int(value)
            .and_then(|v| i64::try_from(v).ok())
            .ok_or_else(|| {
                invalid(format!(
                    "argument --threads: invalid int value: {}",
                    repr(value)
                ))
            })?,
    };
    let jobs = match values.get("jobs") {
        None => 1,
        Some(value) => crate::python_text::python_int(value)
            .and_then(|v| usize::try_from(v).ok())
            .filter(|v| *v > 0)
            .ok_or_else(|| {
                invalid(format!(
                    "argument --jobs: invalid positive int value: {}",
                    repr(value)
                ))
            })?,
    };
    let options = Options {
        manifest: path("manifest").unwrap_or_default(),
        assembly: values["assembly"].clone(),
        output_dir: path("output-dir").unwrap_or_default(),
        receipt: path("receipt").unwrap_or_default(),
        converter_bin: values.get("converter-bin").cloned(),
        threads,
        jobs,
        validation_subset: set_flags.iter().any(|f| f == "validation-subset"),
    };
    builder_error::finish(prepare(&options)).map(|_| ())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identifiers_follow_the_script_patterns() {
        assert_eq!(unversioned("ENSG00000000003.15"), "ENSG00000000003");
        assert_eq!(unversioned("ENSG00000182378.15_PAR_Y"), "ENSG00000182378");
        assert_eq!(
            unversioned("ENSG00000182378_PAR_Y"),
            "ENSG00000182378_PAR_Y"
        );
        assert_eq!(
            stable_id("ENSG00000182378.15_PAR_Y").as_deref(),
            Some("ENSG00000182378_PAR_Y")
        );
        assert_eq!(stable_id("ENSG0000018237.15"), None);
        assert_eq!(unquote("A%2CB%3Dx%"), "A,B=x%");
    }

    #[test]
    fn gtf_attributes_follow_the_findall_pattern() {
        let found = gtf_attributes(
            r#"gene_id "ENSG1.1"; gene_type "protein_coding"; gene_name "A1BG"; level 2; tag "x""#,
        );
        assert_eq!(
            found,
            [
                ("gene_id".into(), "ENSG1.1".into()),
                ("gene_type".into(), "protein_coding".into()),
                ("gene_name".into(), "A1BG".into())
            ]
        );
    }
}
