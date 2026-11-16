//! The annotation hub: stored results of one annotation step, reused across
//! cohort runs. Layout: HUB/SOURCE/LINEAGE/CHROM/part-SESSION.parquet, where
//! SOURCE is the step (gnomad, clinvar, control, vep, cadd), LINEAGE names the
//! resources and code the results came from, and SESSION is the Nextflow
//! session that wrote the file. Files are only ever added, never edited.
//!
//! Each Parquet row is one variant: chrom, pos, ref, alt and value. For the
//! VCF steps the value is the INFO entries the step added (joined by `;`,
//! empty when it added none); for CADD it is the variant's output rows
//! (joined by newlines, empty when CADD wrote none). The step's own header
//! lines are kept once per lineage in the file metadata.
//!
//! - lookup: splits a preprocessed VCF into the variants the hub already
//!   holds (written as a Parquet file in VCF order) and the others (a VCF of
//!   whole records for the step to annotate). Files written by the session
//!   given with --exclude-session are not read, so a resumed run sees what
//!   its first launch saw.
//! - fill: the step's full output (a VCF with every preprocessed record, or a
//!   CADD table) from the stored values and the step's output for the new
//!   variants, in the preprocessed VCF's record order.
//! - add: stores the step's results for the new variants as new files.
//! - refresh-plan: for every step whose current folder (named by the lineages
//!   validation would record) is empty while earlier folders of the step hold
//!   files, writes a VCF of the variants those earlier folders hold, so the
//!   hub refresh run can annotate them into the current folder before the
//!   next cohort needs them. Earlier folders are never changed.
//!
//! Memory: Parquet files are written and read one batch of rows at a time,
//! and lookup holds the stored values of one chromosome at a time, so the
//! memory held follows the largest chromosome and not the genome. `add`
//! closes a chromosome's file when the next chromosome begins, so it needs
//! each chromosome's records together in the preprocessed VCF.
use crate::cli::{Result, invalid};
use crate::vcf_text::open_reader;
use polars_core::prelude::*;
use polars_io::parquet::write::{BatchedWriter, KeyValueMetadata};
use polars_io::prelude::*;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

const USAGE: &str = "usage: priva annotation-hub lookup --hub HUB --source SOURCE --lineage LINEAGE --vcf PREPROCESSED.vcf.gz --known KNOWN.parquet --new NEW.vcf.gz [--exclude-session SESSION]
       priva annotation-hub fill --source SOURCE --vcf PREPROCESSED.vcf.gz --known KNOWN.parquet --annotated STEP_OUTPUT --output OUTPUT
       priva annotation-hub add --hub HUB --source SOURCE --lineage LINEAGE --session SESSION --base NEW.vcf.gz --annotated STEP_OUTPUT
       priva annotation-hub refresh-plan --hub HUB --resources RESOURCES.json --legacy-root CHECKOUT [--resources-lock LOCK.json] --vcf VARIANTS.vcf.gz --plan PLAN.json
       priva annotation-hub lineages --resources RESOURCES.json --legacy-root CHECKOUT [--resources-lock LOCK.json]
SOURCE is gnomad, clinvar, control, vep or cadd; for cadd the step output is the CADD table.";
const SOURCES: [&str; 5] = ["gnomad", "clinvar", "control", "vep", "cadd"];
const HEADER_KEY: &str = "priva.hub.header";

type Key = (String, String, String, String);

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(command) = args.first() else {
        return Err(invalid(USAGE));
    };
    if command == "--help" || command == "-h" {
        println!("{USAGE}");
        return Ok(());
    }
    let mut options = HashMap::new();
    let mut rest = args[1..].iter();
    while let Some(option) = rest.next() {
        let value = rest
            .next()
            .ok_or_else(|| invalid(format!("{option} needs a value\n{USAGE}")))?;
        options.insert(option.trim_start_matches("--").to_owned(), value.clone());
    }
    let get = |name: &str| -> Result<&str> {
        options
            .get(name)
            .map(String::as_str)
            .ok_or_else(|| invalid(format!("--{name} is required\n{USAGE}")))
    };
    if command == "lineages" {
        // Each step's hub folder name with the code and resources it hashes.
        use crate::nextflow_task::contract::{json_document, read_json};
        let lock = match options.get("resources-lock") {
            Some(path) => Some(read_json(Path::new(path))?),
            None => None,
        };
        let explained = crate::nextflow_task::hub_lineage::explain(
            &read_json(Path::new(get("resources")?))?,
            Path::new(get("legacy-root")?),
            lock.as_ref(),
        )?;
        print!("{}", json_document(&explained)?);
        return Ok(());
    }
    if command == "refresh-plan" {
        return refresh_plan(
            Path::new(get("hub")?),
            Path::new(get("resources")?),
            Path::new(get("legacy-root")?),
            options.get("resources-lock").map(Path::new),
            Path::new(get("vcf")?),
            Path::new(get("plan")?),
        );
    }
    let source = get("source")?;
    if !SOURCES.contains(&source) {
        return Err(invalid(format!("unknown source {source:?}\n{USAGE}")));
    }
    match command.as_str() {
        "lookup" => lookup(
            Path::new(get("hub")?),
            source,
            get("lineage")?,
            options.get("exclude-session").map(String::as_str),
            Path::new(get("vcf")?),
            Path::new(get("known")?),
            Path::new(get("new")?),
        ),
        "fill" => fill(
            source,
            Path::new(get("vcf")?),
            Path::new(get("known")?),
            Path::new(get("annotated")?),
            Path::new(get("output")?),
        ),
        "add" => add(
            Path::new(get("hub")?),
            source,
            get("lineage")?,
            get("session")?,
            Path::new(get("base")?),
            Path::new(get("annotated")?),
        ),
        other => Err(invalid(format!("unknown command {other:?}\n{USAGE}"))),
    }
}

fn check_name(label: &str, value: &str) -> Result<()> {
    if value.is_empty()
        || !value
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'-' || b == b'_' || b == b'.')
        || value.starts_with('.')
    {
        return Err(invalid(format!(
            "{label} {value:?} must be letters, digits, '-', '_' or '.'"
        )));
    }
    Ok(())
}

/// A VCF's `##` lines, its `#CHROM` line and its record lines.
struct Vcf {
    meta: Vec<String>,
    chrom_line: String,
    records: Box<dyn Iterator<Item = std::io::Result<String>>>,
}

fn read_vcf(path: &Path) -> Result<Vcf> {
    let mut lines = open_reader(path)?.lines();
    let mut meta = Vec::new();
    for line in lines.by_ref() {
        let line = line?;
        if line.starts_with("##") {
            meta.push(line);
        } else if line.starts_with("#CHROM") {
            return Ok(Vcf {
                meta,
                chrom_line: line,
                records: Box::new(lines),
            });
        } else {
            break;
        }
    }
    Err(invalid(format!(
        "{}: no #CHROM header line",
        path.display()
    )))
}

fn key_of(line: &str) -> Result<Key> {
    let mut fields = line.splitn(6, '\t');
    let mut next = || fields.next().map(str::to_owned);
    match (next(), next(), next(), next(), next()) {
        (Some(chrom), Some(pos), Some(_id), Some(reference), Some(alt)) => {
            Ok((chrom, pos, reference, alt))
        }
        _ => Err(invalid(format!(
            "VCF record with fewer than 5 columns: {line}"
        ))),
    }
}

fn info_of(line: &str) -> Result<&str> {
    line.split('\t')
        .nth(7)
        .ok_or_else(|| invalid("VCF record with fewer than 8 columns"))
}

fn source_dir(hub: &Path, source: &str, lineage: &str) -> PathBuf {
    hub.join(source).join(lineage)
}

/// A part or known file is written and read this many rows (or this many
/// bytes of values) at a time.
const BATCH_ROWS: usize = 100_000;
const BATCH_BYTES: usize = 256 << 20;

/// The part files of one chromosome in a hub folder, by name, without the
/// file of the excluded session.
fn part_files(folder: &Path, chrom: &str, exclude: Option<&str>) -> Vec<PathBuf> {
    let skip = exclude.map(|s| format!("part-{s}.parquet"));
    let Ok(entries) = fs::read_dir(folder.join(chrom)) else {
        return Vec::new();
    };
    let mut parts: Vec<PathBuf> = entries
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| {
            let name = p.file_name().and_then(|n| n.to_str()).unwrap_or_default();
            name.starts_with("part-") && name.ends_with(".parquet") && skip.as_deref() != Some(name)
        })
        .collect();
    parts.sort();
    parts
}

/// A Parquet file of hub rows and its stored header, read one row group at a
/// time.
struct RowReader {
    path: PathBuf,
    header: Option<String>,
    groups: Vec<usize>,
    next: usize,
    offset: usize,
    batch: std::vec::IntoIter<(Key, String)>,
}

impl RowReader {
    fn open(path: &Path) -> Result<Self> {
        let mut reader = ParquetReader::new(File::open(path)?);
        let metadata = reader.get_metadata()?;
        let header = metadata
            .key_value_metadata
            .as_deref()
            .unwrap_or_default()
            .iter()
            .find(|kv| kv.key == HEADER_KEY)
            .and_then(|kv| kv.value.clone());
        let groups = metadata.row_groups.iter().map(|g| g.num_rows()).collect();
        Ok(Self {
            path: path.to_owned(),
            header,
            groups,
            next: 0,
            offset: 0,
            batch: Vec::new().into_iter(),
        })
    }

    /// The next row in the file's order, reading the next row group when the
    /// rows of the last one are used up.
    fn next_row(&mut self) -> Result<Option<(Key, String)>> {
        loop {
            if let Some(row) = self.batch.next() {
                return Ok(Some(row));
            }
            let Some(&rows) = self.groups.get(self.next) else {
                return Ok(None);
            };
            self.next += 1;
            if rows == 0 {
                continue;
            }
            let frame = ParquetReader::new(File::open(&self.path)?)
                .read_parallel(ParallelStrategy::RowGroups)
                .with_slice(Some((self.offset, rows)))
                .finish()?;
            self.offset += rows;
            let column = |name: &str| -> Result<Vec<String>> {
                Ok(frame
                    .column(name)?
                    .as_materialized_series()
                    .str()?
                    .into_iter()
                    .map(|v| v.unwrap_or_default().to_owned())
                    .collect())
            };
            let (c, p, r, a, v) = (
                column("chrom")?,
                column("pos")?,
                column("ref")?,
                column("alt")?,
                column("value")?,
            );
            let batch: Vec<(Key, String)> = c
                .into_iter()
                .zip(p)
                .zip(r)
                .zip(a)
                .zip(v)
                .map(|((((c, p), r), a), v)| ((c, p, r, a), v))
                .collect();
            self.batch = batch.into_iter();
        }
    }
}

/// Every stored value of one chromosome, the first file (by name) winning
/// when two sessions stored one variant.
fn load_chrom(folder: &Path, chrom: &str, exclude: Option<&str>) -> Result<HashMap<Key, String>> {
    let mut values = HashMap::new();
    for part in part_files(folder, chrom, exclude) {
        let mut reader = RowReader::open(&part)?;
        while let Some((key, value)) = reader.next_row()? {
            values.entry(key).or_insert(value);
        }
    }
    Ok(values)
}

/// The header a hub folder stores: that of the first part file that has one,
/// taking the chromosomes in the given order.
fn stored_header(
    folder: &Path,
    chroms: &[String],
    exclude: Option<&str>,
) -> Result<Option<String>> {
    for chrom in chroms {
        for part in part_files(folder, chrom, exclude) {
            if let Some(header) = RowReader::open(&part)?.header {
                return Ok(Some(header));
            }
        }
    }
    Ok(None)
}

fn frame(rows: &[(Key, String)]) -> Result<DataFrame> {
    let pick = |f: fn(&(Key, String)) -> &str| -> Vec<&str> { rows.iter().map(f).collect() };
    Ok(DataFrame::new(vec![
        Series::new("chrom".into(), pick(|r| &r.0.0)).into(),
        Series::new("pos".into(), pick(|r| &r.0.1)).into(),
        Series::new("ref".into(), pick(|r| &r.0.2)).into(),
        Series::new("alt".into(), pick(|r| &r.0.3)).into(),
        Series::new("value".into(), pick(|r| &r.1)).into(),
    ])?)
}

/// A Parquet file of hub rows, written one batch (one row group) at a time
/// under a staged name and renamed when it is complete.
struct RowWriter {
    path: PathBuf,
    header: Option<String>,
    writer: Option<BatchedWriter<File>>,
    rows: Vec<(Key, String)>,
    bytes: usize,
    batch_rows: usize,
    written: usize,
}

impl RowWriter {
    fn new(path: &Path, header: Option<&str>) -> Self {
        Self {
            path: path.to_owned(),
            header: header.map(str::to_owned),
            writer: None,
            rows: Vec::new(),
            bytes: 0,
            batch_rows: BATCH_ROWS,
            written: 0,
        }
    }

    fn staged(&self) -> PathBuf {
        self.path.with_extension("parquet.partial")
    }

    fn push(&mut self, key: Key, value: String) -> Result<()> {
        self.bytes += value.len();
        self.rows.push((key, value));
        if self.rows.len() >= self.batch_rows || self.bytes >= BATCH_BYTES {
            self.write_batch()?;
        }
        Ok(())
    }

    fn write_batch(&mut self) -> Result<()> {
        let data = frame(&self.rows)?;
        if self.writer.is_none() {
            let metadata = self.header.as_deref().map(|h| {
                KeyValueMetadata::from_static(vec![(HEADER_KEY.to_owned(), h.to_owned())])
            });
            self.writer = Some(
                ParquetWriter::new(File::create(self.staged())?)
                    .with_compression(ParquetCompression::Zstd(None))
                    .with_key_value_metadata(metadata)
                    .batched(data.schema())?,
            );
        }
        if let Some(writer) = self.writer.as_mut() {
            writer.write_batch(&data)?;
        }
        self.written += self.rows.len();
        self.rows.clear();
        self.bytes = 0;
        Ok(())
    }

    /// Writes the rows still held and the file's footer, and gives the
    /// number of rows in the file. A file without rows keeps its header.
    fn finish(mut self) -> Result<usize> {
        if !self.rows.is_empty() || self.writer.is_none() {
            self.write_batch()?;
        }
        if let Some(writer) = self.writer.take() {
            writer.finish()?;
        }
        File::open(self.staged())?.sync_all()?;
        fs::rename(self.staged(), &self.path)?;
        Ok(self.written)
    }
}

fn bgzip_and_index(plain: &Path, output: &Path) -> Result<()> {
    let status = Command::new("bgzip").arg("-f").arg(plain).status()?;
    if !status.success() {
        return Err(invalid(format!("bgzip failed on {}", plain.display())));
    }
    fs::rename(PathBuf::from(format!("{}.gz", plain.display())), output)?;
    let status = Command::new("tabix")
        .args(["-f", "-p", "vcf"])
        .arg(output)
        .status()?;
    if !status.success() {
        return Err(invalid(format!("tabix failed on {}", output.display())));
    }
    Ok(())
}

fn chroms_of(path: &Path) -> Result<Vec<String>> {
    let vcf = read_vcf(path)?;
    let mut seen = Vec::new();
    for line in vcf.records {
        let line = line?;
        let chrom = line.split('\t').next().unwrap_or_default();
        if seen.last().map(String::as_str) != Some(chrom) && !seen.iter().any(|c| c == chrom) {
            seen.push(chrom.to_owned());
        }
    }
    Ok(seen)
}

/// Whether a VCF's first record names its chromosome with "chr"; CADD
/// writes chromosome names without it.
fn uses_chr(path: &Path) -> Result<bool> {
    let mut vcf = read_vcf(path)?;
    match vcf.records.next() {
        Some(line) => Ok(line?.starts_with("chr")),
        None => Ok(true),
    }
}

fn lookup(
    hub: &Path,
    source: &str,
    lineage: &str,
    exclude: Option<&str>,
    vcf_path: &Path,
    known: &Path,
    new: &Path,
) -> Result<()> {
    check_name("lineage", lineage)?;
    if let Some(session) = exclude {
        check_name("session", session)?;
    }
    let folder = source_dir(hub, source, lineage);
    let header = stored_header(&folder, &chroms_of(vcf_path)?, exclude)?;
    let vcf = read_vcf(vcf_path)?;
    let plain = PathBuf::from(format!("{}.plain.vcf", new.display()));
    let mut out = BufWriter::new(File::create(&plain)?);
    for line in &vcf.meta {
        writeln!(out, "{line}")?;
    }
    writeln!(out, "{}", vcf.chrom_line)?;
    let mut hits = RowWriter::new(known, header.as_deref());
    // The stored values of the chromosome the records are in.
    let mut chrom: Option<String> = None;
    let mut stored = HashMap::new();
    let mut misses = 0usize;
    for line in vcf.records {
        let line = line?;
        let key = key_of(&line)?;
        if chrom.as_deref() != Some(key.0.as_str()) {
            // The last chromosome's values go before the next ones come.
            drop(std::mem::take(&mut stored));
            stored = load_chrom(&folder, &key.0, exclude)?;
            chrom = Some(key.0.clone());
        }
        match stored.get(&key) {
            Some(value) => hits.push(key, value.clone())?,
            None => {
                writeln!(out, "{line}")?;
                misses += 1;
            }
        }
    }
    out.flush()?;
    drop(out);
    bgzip_and_index(&plain, new)?;
    let hits = hits.finish()?;
    eprintln!("annotation-hub lookup {source}: {hits} known, {misses} new");
    Ok(())
}

/// The INFO entries `annotated` has beyond those of `base`, in its order.
fn added_entries(base: &str, annotated: &str) -> String {
    let keys: HashSet<&str> = base
        .split(';')
        .filter(|e| *e != ".")
        .map(|e| e.split_once('=').map_or(e, |(k, _)| k))
        .collect();
    annotated
        .split(';')
        .filter(|e| *e != "." && !e.is_empty())
        .filter(|e| !keys.contains(e.split_once('=').map_or(*e, |(k, _)| k)))
        .collect::<Vec<_>>()
        .join(";")
}

fn with_entries(record: &str, entries: &str) -> String {
    if entries.is_empty() {
        return record.to_owned();
    }
    let mut fields: Vec<&str> = record.split('\t').collect();
    let joined = if fields[7] == "." || fields[7].is_empty() {
        entries.to_owned()
    } else {
        format!("{};{entries}", fields[7])
    };
    fields[7] = &joined;
    fields.join("\t")
}

/// The CADD table's leading `#` lines and its rows grouped by variant, keyed
/// with the VCF's chromosome spelling.
fn cadd_rows(
    path: &Path,
    chr_prefix: bool,
) -> Result<(String, BTreeMap<Key, Vec<String>>, Vec<Key>)> {
    let mut header = Vec::new();
    let mut rows: BTreeMap<Key, Vec<String>> = BTreeMap::new();
    let mut order = Vec::new();
    for line in open_reader(path)?.lines() {
        let line = line?;
        if line.starts_with('#') && rows.is_empty() {
            header.push(line);
            continue;
        }
        let mut f = line.splitn(5, '\t');
        let (Some(c), Some(p), Some(r), Some(a)) = (f.next(), f.next(), f.next(), f.next()) else {
            return Err(invalid(format!(
                "CADD row with fewer than 4 columns in {}",
                path.display()
            )));
        };
        let chrom = if chr_prefix && !c.starts_with("chr") {
            format!("chr{c}")
        } else {
            c.to_owned()
        };
        let key = (chrom, p.to_owned(), r.to_owned(), a.to_owned());
        let entry = rows.entry(key.clone()).or_default();
        if entry.is_empty() {
            order.push(key);
        }
        entry.push(line);
    }
    Ok((header.join("\n"), rows, order))
}

fn fill(
    source: &str,
    vcf_path: &Path,
    known: &Path,
    annotated: &Path,
    output: &Path,
) -> Result<()> {
    let mut known_rows = RowReader::open(known)?;
    let stored_header = known_rows.header.clone();
    let mut next_known = known_rows.next_row()?;
    let vcf = read_vcf(vcf_path)?;
    let has_new = fs::metadata(annotated)
        .map(|m| m.len() > 0)
        .unwrap_or(false);
    if source == "cadd" {
        let (fresh_header, mut fresh, _) = if has_new {
            cadd_rows(annotated, uses_chr(vcf_path)?)?
        } else {
            (String::new(), BTreeMap::new(), Vec::new())
        };
        let header = if has_new && !fresh_header.is_empty() {
            fresh_header
        } else {
            stored_header.unwrap_or_default()
        };
        let mut out = BufWriter::new(File::create(output)?);
        if !header.is_empty() {
            writeln!(out, "{header}")?;
        }
        for line in vcf.records {
            let key = key_of(&line?)?;
            let value = match &next_known {
                Some((k, value)) if *k == key => {
                    let value = value.clone();
                    next_known = known_rows.next_row()?;
                    value
                }
                _ => fresh
                    .remove(&key)
                    .map(|rows| rows.join("\n"))
                    .unwrap_or_default(),
            };
            if !value.is_empty() {
                writeln!(out, "{value}")?;
            }
        }
        out.flush()?;
    } else {
        let mut fresh = if has_new {
            Some(read_vcf(annotated)?)
        } else {
            None
        };
        let base_lines: HashSet<&String> = vcf.meta.iter().collect();
        let added_meta: Vec<String> = match (&fresh, &stored_header) {
            (Some(f), _) => f
                .meta
                .iter()
                .filter(|l| !base_lines.contains(l))
                .cloned()
                .collect(),
            (None, Some(h)) => h
                .lines()
                .filter(|l| !l.is_empty())
                .map(str::to_owned)
                .collect(),
            (None, None) => Vec::new(),
        };
        let plain = PathBuf::from(format!("{}.plain.vcf", output.display()));
        let mut out = BufWriter::new(File::create(&plain)?);
        for line in vcf.meta.iter().chain(added_meta.iter()) {
            writeln!(out, "{line}")?;
        }
        writeln!(out, "{}", vcf.chrom_line)?;
        for line in vcf.records {
            let line = line?;
            let key = key_of(&line)?;
            match &next_known {
                Some((k, value)) if *k == key => {
                    writeln!(out, "{}", with_entries(&line, value))?;
                    next_known = known_rows.next_row()?;
                }
                _ => {
                    let Some(f) = fresh.as_mut() else {
                        return Err(invalid(format!(
                            "{source}: variant {key:?} is neither stored nor annotated"
                        )));
                    };
                    let annotated_line = f.records.next().ok_or_else(|| {
                        invalid(format!("{source}: the step's output ended before {key:?}"))
                    })??;
                    if key_of(&annotated_line)? != key {
                        return Err(invalid(format!(
                            "{source}: the step's output is out of order at {key:?}"
                        )));
                    }
                    writeln!(
                        out,
                        "{}",
                        with_entries(
                            &line,
                            &added_entries(info_of(&line)?, info_of(&annotated_line)?)
                        )
                    )?;
                }
            }
        }
        if let Some(f) = fresh.as_mut() {
            if f.records.next().is_some() {
                return Err(invalid(format!(
                    "{source}: the step's output has records the preprocessed VCF lacks"
                )));
            }
        }
        out.flush()?;
        drop(out);
        bgzip_and_index(&plain, output)?;
    }
    if next_known.is_some() {
        return Err(invalid(format!(
            "{source}: stored values out of the preprocessed VCF's order"
        )));
    }
    Ok(())
}

/// The part files one `add` writes: one per chromosome, written while the
/// chromosome's records pass and closed when the next chromosome begins.
struct PartStore<'a> {
    folder: &'a Path,
    source: &'a str,
    session: &'a str,
    header: &'a str,
    current: Option<(String, RowWriter)>,
    closed: HashSet<String>,
    stored: usize,
}

impl<'a> PartStore<'a> {
    fn new(folder: &'a Path, source: &'a str, session: &'a str, header: &'a str) -> Self {
        Self {
            folder,
            source,
            session,
            header,
            current: None,
            closed: HashSet::new(),
            stored: 0,
        }
    }

    fn push(&mut self, key: Key, value: String) -> Result<()> {
        if self.current.as_ref().map(|(chrom, _)| chrom) != Some(&key.0) {
            self.close()?;
            check_name("chromosome", &key.0)?;
            if self.closed.contains(&key.0) {
                return Err(invalid(format!(
                    "{}: the records of {} are not together in the preprocessed VCF",
                    self.source, key.0
                )));
            }
            let dir = self.folder.join(&key.0);
            fs::create_dir_all(&dir)?;
            let part = dir.join(format!("part-{}.parquet", self.session));
            self.current = Some((key.0.clone(), RowWriter::new(&part, Some(self.header))));
        }
        match self.current.as_mut() {
            Some((_, writer)) => writer.push(key, value),
            None => Ok(()),
        }
    }

    fn close(&mut self) -> Result<()> {
        if let Some((chrom, writer)) = self.current.take() {
            self.stored += writer.finish()?;
            self.closed.insert(chrom);
        }
        Ok(())
    }

    /// Closes the last chromosome's file and gives the number of variants
    /// stored.
    fn finish(mut self) -> Result<usize> {
        self.close()?;
        Ok(self.stored)
    }
}

fn add(
    hub: &Path,
    source: &str,
    lineage: &str,
    session: &str,
    base: &Path,
    annotated: &Path,
) -> Result<()> {
    check_name("lineage", lineage)?;
    check_name("session", session)?;
    let folder = source_dir(hub, source, lineage);
    let base_vcf = read_vcf(base)?;
    let stored;
    if source == "cadd" {
        let (header, mut rows, _) = cadd_rows(annotated, uses_chr(base)?)?;
        let mut store = PartStore::new(&folder, source, session, &header);
        for line in base_vcf.records {
            let key = key_of(&line?)?;
            let value = rows.remove(&key).map(|r| r.join("\n")).unwrap_or_default();
            store.push(key, value)?;
        }
        stored = store.finish()?;
    } else {
        let mut fresh = read_vcf(annotated)?;
        let base_lines: HashSet<&String> = base_vcf.meta.iter().collect();
        let header = fresh
            .meta
            .iter()
            .filter(|l| !base_lines.contains(l))
            .cloned()
            .collect::<Vec<_>>()
            .join("\n");
        let mut store = PartStore::new(&folder, source, session, &header);
        for line in base_vcf.records {
            let line = line?;
            let key = key_of(&line)?;
            let annotated_line = fresh.records.next().ok_or_else(|| {
                invalid(format!("{source}: the step's output ended before {key:?}"))
            })??;
            if key_of(&annotated_line)? != key {
                return Err(invalid(format!(
                    "{source}: the step's output is out of order at {key:?}"
                )));
            }
            let value = added_entries(info_of(&line)?, info_of(&annotated_line)?);
            store.push(key, value)?;
        }
        stored = store.finish()?;
    }
    eprintln!(
        "annotation-hub add {source}: {stored} variants stored under {}",
        folder.display()
    );
    Ok(())
}

fn subfolders(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .map(|entries| {
            entries
                .filter_map(|e| e.ok())
                .filter(|e| e.path().is_dir())
                .filter_map(|e| e.file_name().into_string().ok())
                .collect()
        })
        .unwrap_or_default();
    names.sort();
    names
}

fn has_parts(folder: &Path) -> bool {
    subfolders(folder).iter().any(|chrom| {
        fs::read_dir(folder.join(chrom))
            .map(|entries| {
                entries.filter_map(|e| e.ok()).any(|e| {
                    let name = e.file_name().to_string_lossy().into_owned();
                    name.starts_with("part-") && name.ends_with(".parquet")
                })
            })
            .unwrap_or(false)
    })
}

/// Chromosome order for the refresh VCF: 1-22, X, Y, M, then others by name.
fn chrom_rank(chrom: &str) -> (u32, String) {
    let bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    let rank = match bare {
        "X" => 23,
        "Y" => 24,
        "M" | "MT" => 25,
        other => other.parse::<u32>().unwrap_or(26),
    };
    (rank, chrom.to_owned())
}

fn refresh_plan(
    hub: &Path,
    resources: &Path,
    root: &Path,
    lock: Option<&Path>,
    vcf_out: &Path,
    plan_out: &Path,
) -> Result<()> {
    use crate::json_text::Json;
    use crate::nextflow_task::contract::{object, read_json, text, write_json};
    let manifest = read_json(resources)?;
    let lock = match lock {
        Some(path) => Some(read_json(path)?),
        None => None,
    };
    let lineages = crate::nextflow_task::hub_lineage::lineages(&manifest, root, lock.as_ref())?;
    let mut variants: BTreeSet<((u32, String), u64, String, String)> = BTreeSet::new();
    let mut steps = Vec::new();
    for source in SOURCES {
        let lineage = lineages
            .get(source)
            .and_then(Json::as_str)
            .unwrap_or("off")
            .to_owned();
        let current = source_dir(hub, source, &lineage);
        let earlier: Vec<String> = subfolders(&hub.join(source))
            .into_iter()
            .filter(|name| *name != lineage && has_parts(&hub.join(source).join(name)))
            .collect();
        let state = if lineage == "off" {
            "no hub folder"
        } else if has_parts(&current) {
            "current"
        } else if earlier.is_empty() {
            "nothing stored yet"
        } else {
            "refresh"
        };
        let mut count = 0usize;
        if state == "refresh" {
            for old in &earlier {
                let folder = hub.join(source).join(old);
                for stored_chrom in subfolders(&folder) {
                    let stored = load_chrom(&folder, &stored_chrom, None)?;
                    for (chrom, pos, reference, alt) in stored.into_keys() {
                        let position = pos.parse::<u64>().map_err(|_| {
                            invalid(format!("{source}: stored position {pos:?} is not a number"))
                        })?;
                        if variants.insert((chrom_rank(&chrom), position, reference, alt)) {
                            count += 1;
                        }
                    }
                }
            }
        }
        steps.push((
            source.to_owned(),
            object(vec![
                ("lineage", text(lineage.clone())),
                ("state", text(state)),
                (
                    "earlier",
                    Json::Array(earlier.iter().cloned().map(text).collect()),
                ),
                ("variants_added", Json::Int(count as i128)),
            ]),
        ));
    }
    let refresh: Vec<Json> = steps
        .iter()
        .filter(|(_, step)| step.get("state").and_then(Json::as_str) == Some("refresh"))
        .map(|(source, _)| text(source.clone()))
        .collect();
    let plain = PathBuf::from(format!("{}.plain.vcf", vcf_out.display()));
    let mut out = BufWriter::new(File::create(&plain)?);
    writeln!(out, "##fileformat=VCFv4.2")?;
    let mut contigs: Vec<&(u32, String)> = variants.iter().map(|v| &v.0).collect();
    contigs.dedup();
    for (_, chrom) in contigs {
        writeln!(out, "##contig=<ID={chrom}>")?;
    }
    writeln!(
        out,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    writeln!(
        out,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\thub_refresh"
    )?;
    for ((_, chrom), pos, reference, alt) in &variants {
        writeln!(
            out,
            "{chrom}\t{pos}\t.\t{reference}\t{alt}\t.\tPASS\t.\tGT\t0/1"
        )?;
    }
    out.flush()?;
    drop(out);
    bgzip_and_index(&plain, vcf_out)?;
    write_json(
        plan_out,
        &object(vec![
            ("refresh", Json::Array(refresh)),
            ("variants", Json::Int(variants.len() as i128)),
            ("steps", Json::Object(steps)),
        ]),
    )?;
    eprintln!(
        "annotation-hub refresh-plan: {} variants to refresh",
        variants.len()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn added_entries_are_the_new_keys_in_step_order() {
        assert_eq!(
            added_entries("DP=3;AC=1", "DP=3;AC=1;CSQ=x|y;AF_joint=0.1"),
            "CSQ=x|y;AF_joint=0.1"
        );
        assert_eq!(added_entries(".", "CSQ=a"), "CSQ=a");
        assert_eq!(added_entries("DP=3", "DP=3"), "");
        assert_eq!(
            with_entries("1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1", "CSQ=a"),
            "1\t10\t.\tA\tG\t.\tPASS\tCSQ=a\tGT\t0/1"
        );
        assert_eq!(
            with_entries("1\t10\t.\tA\tG\t.\tPASS\tDP=3", "CSQ=a"),
            "1\t10\t.\tA\tG\t.\tPASS\tDP=3;CSQ=a"
        );
    }

    #[test]
    fn names_with_path_parts_are_refused() {
        assert!(check_name("session", "abc-1_2.x").is_ok());
        for bad in ["", "../x", "a/b", ".hidden"] {
            assert!(check_name("session", bad).is_err(), "{bad}");
        }
    }

    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("priva_hub_{name}_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn key(chrom: &str, pos: &str) -> Key {
        (
            chrom.to_owned(),
            pos.to_owned(),
            "A".to_owned(),
            "G".to_owned(),
        )
    }

    #[test]
    fn rows_are_written_and_read_back_one_batch_at_a_time() {
        let dir = scratch("batches");
        let file = dir.join("rows.parquet");
        let rows: Vec<(Key, String)> = (1..=5)
            .map(|i| (key("chr1", &i.to_string()), format!("AF={i}")))
            .collect();
        let mut writer = RowWriter::new(&file, Some("##one\n##two"));
        writer.batch_rows = 2;
        for (key, value) in rows.clone() {
            writer.push(key, value).unwrap();
        }
        assert_eq!(writer.finish().unwrap(), 5);
        assert!(!file.with_extension("parquet.partial").exists());
        let mut reader = RowReader::open(&file).unwrap();
        assert_eq!(reader.groups, vec![2, 2, 1]);
        assert_eq!(reader.header.as_deref(), Some("##one\n##two"));
        let mut read = Vec::new();
        while let Some(row) = reader.next_row().unwrap() {
            read.push(row);
        }
        assert_eq!(read, rows);

        let empty = dir.join("empty.parquet");
        assert_eq!(RowWriter::new(&empty, Some("##h")).finish().unwrap(), 0);
        let mut reader = RowReader::open(&empty).unwrap();
        assert_eq!(reader.header.as_deref(), Some("##h"));
        assert!(reader.next_row().unwrap().is_none());
        fs::remove_dir_all(&dir).unwrap();
    }

    fn vcf(path: &Path, extra_meta: &[&str], records: &[(&str, &str, &str)]) {
        let mut text = String::from("##fileformat=VCFv4.2\n");
        for line in extra_meta {
            text.push_str(line);
            text.push('\n');
        }
        text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for (chrom, pos, info) in records {
            text.push_str(&format!("{chrom}\t{pos}\t.\tA\tG\t.\tPASS\t{info}\n"));
        }
        fs::write(path, text).unwrap();
    }

    #[test]
    fn add_writes_one_file_per_chromosome_as_the_records_pass() {
        let dir = scratch("add");
        let (hub, base, annotated) = (dir.join("hub"), dir.join("base.vcf"), dir.join("anno.vcf"));
        let meta = "##INFO=<ID=AF,Number=1,Type=Float,Description=\"x\">";
        vcf(
            &base,
            &[],
            &[
                ("chr1", "10", "DP=1"),
                ("chr1", "20", "."),
                ("chr2", "5", "DP=3"),
            ],
        );
        vcf(
            &annotated,
            &[meta],
            &[
                ("chr1", "10", "DP=1;AF=0.1"),
                ("chr1", "20", "."),
                ("chr2", "5", "DP=3;AF=0.3"),
            ],
        );
        add(&hub, "gnomad", "lineage", "s1", &base, &annotated).unwrap();
        let folder = source_dir(&hub, "gnomad", "lineage");
        assert_eq!(subfolders(&folder), vec!["chr1", "chr2"]);
        let chr1 = load_chrom(&folder, "chr1", None).unwrap();
        assert_eq!(chr1.len(), 2);
        assert_eq!(chr1[&key("chr1", "10")], "AF=0.1");
        assert_eq!(chr1[&key("chr1", "20")], "");
        let chr2 = load_chrom(&folder, "chr2", None).unwrap();
        assert_eq!(chr2[&key("chr2", "5")], "AF=0.3");
        assert!(load_chrom(&folder, "chr1", Some("s1")).unwrap().is_empty());
        let chroms = vec!["chr9".to_owned(), "chr2".to_owned()];
        assert_eq!(
            stored_header(&folder, &chroms, None).unwrap().as_deref(),
            Some(meta)
        );
        assert_eq!(stored_header(&folder, &chroms, Some("s1")).unwrap(), None);

        // A chromosome that returns after its file was closed is refused.
        let split = [("chr1", "10", "."), ("chr2", "5", "."), ("chr1", "20", ".")];
        vcf(&base, &[], &split);
        vcf(&annotated, &[meta], &split);
        let error = add(&hub, "gnomad", "lineage", "s2", &base, &annotated).unwrap_err();
        assert!(error.to_string().contains("not together"), "{error}");
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn fill_takes_stored_rows_across_batches_in_the_vcf_order() {
        let dir = scratch("fill");
        let (input, known, table, output) = (
            dir.join("in.vcf"),
            dir.join("known.parquet"),
            dir.join("new.tsv"),
            dir.join("out.tsv"),
        );
        vcf(
            &input,
            &[],
            &[("chr1", "10", "."), ("chr1", "20", "."), ("chr2", "5", ".")],
        );
        let header = "## CADD\n#Chrom\tPos\tRef\tAlt\tPHRED";
        let mut writer = RowWriter::new(&known, Some(header));
        writer.batch_rows = 1;
        writer
            .push(key("chr1", "10"), "1\t10\tA\tG\t1.5".to_owned())
            .unwrap();
        writer
            .push(key("chr2", "5"), "2\t5\tA\tG\t2.5".to_owned())
            .unwrap();
        writer.finish().unwrap();
        fs::write(&table, format!("{header}\n1\t20\tA\tG\t3.5\n")).unwrap();
        fill("cadd", &input, &known, &table, &output).unwrap();
        assert_eq!(
            fs::read_to_string(&output).unwrap(),
            format!("{header}\n1\t10\tA\tG\t1.5\n1\t20\tA\tG\t3.5\n2\t5\tA\tG\t2.5\n")
        );
        fs::remove_dir_all(&dir).unwrap();
    }
}
