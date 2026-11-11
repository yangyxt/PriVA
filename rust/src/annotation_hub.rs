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
use crate::cli::{Result, invalid};
use crate::vcf_text::open_reader;
use polars_core::prelude::*;
use polars_io::parquet::write::KeyValueMetadata;
use polars_io::prelude::*;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

const USAGE: &str = "usage: priva annotation-hub lookup --hub HUB --source SOURCE --lineage LINEAGE --vcf PREPROCESSED.vcf.gz --known KNOWN.parquet --new NEW.vcf.gz [--exclude-session SESSION]
       priva annotation-hub fill --source SOURCE --vcf PREPROCESSED.vcf.gz --known KNOWN.parquet --annotated STEP_OUTPUT --output OUTPUT
       priva annotation-hub add --hub HUB --source SOURCE --lineage LINEAGE --session SESSION --base NEW.vcf.gz --annotated STEP_OUTPUT
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

/// Every stored value of the chromosomes in `wanted`, the first file (by
/// name) winning when two sessions stored one variant, and the stored header.
fn load(
    folder: &Path,
    wanted: &HashSet<String>,
    exclude: Option<&str>,
) -> Result<(HashMap<Key, String>, Option<String>)> {
    let mut values = HashMap::new();
    let mut header = None;
    let skip = exclude.map(|s| format!("part-{s}.parquet"));
    for chrom in wanted {
        let dir = folder.join(chrom);
        let Ok(entries) = fs::read_dir(&dir) else {
            continue;
        };
        let mut parts: Vec<PathBuf> = entries
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| {
                let name = p.file_name().and_then(|n| n.to_str()).unwrap_or_default();
                name.starts_with("part-")
                    && name.ends_with(".parquet")
                    && skip.as_deref() != Some(name)
            })
            .collect();
        parts.sort();
        for part in parts {
            let mut reader = ParquetReader::new(File::open(&part)?);
            if header.is_none() {
                header = reader
                    .get_metadata()?
                    .key_value_metadata
                    .as_deref()
                    .unwrap_or_default()
                    .iter()
                    .find(|kv| kv.key == HEADER_KEY)
                    .and_then(|kv| kv.value.clone());
            }
            let frame = reader.finish()?;
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
            for i in 0..frame.height() {
                values
                    .entry((c[i].clone(), p[i].clone(), r[i].clone(), a[i].clone()))
                    .or_insert_with(|| v[i].clone());
            }
        }
    }
    Ok((values, header))
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

fn write_parquet(path: &Path, rows: &[(Key, String)], header: Option<&str>) -> Result<()> {
    let staged = path.with_extension("parquet.partial");
    let mut data = frame(rows)?;
    let metadata =
        header.map(|h| KeyValueMetadata::from_static(vec![(HEADER_KEY.to_owned(), h.to_owned())]));
    ParquetWriter::new(File::create(&staged)?)
        .with_compression(ParquetCompression::Zstd(None))
        .with_key_value_metadata(metadata)
        .finish(&mut data)?;
    File::open(&staged)?.sync_all()?;
    fs::rename(&staged, path)?;
    Ok(())
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
    let chroms: HashSet<String> = chroms_of(vcf_path)?.into_iter().collect();
    let (stored, header) = load(&source_dir(hub, source, lineage), &chroms, exclude)?;
    let vcf = read_vcf(vcf_path)?;
    let plain = PathBuf::from(format!("{}.plain.vcf", new.display()));
    let mut out = BufWriter::new(File::create(&plain)?);
    for line in &vcf.meta {
        writeln!(out, "{line}")?;
    }
    writeln!(out, "{}", vcf.chrom_line)?;
    let (mut hits, mut misses) = (Vec::new(), 0usize);
    for line in vcf.records {
        let line = line?;
        let key = key_of(&line)?;
        match stored.get(&key) {
            Some(value) => hits.push((key, value.clone())),
            None => {
                writeln!(out, "{line}")?;
                misses += 1;
            }
        }
    }
    out.flush()?;
    drop(out);
    bgzip_and_index(&plain, new)?;
    write_parquet(known, &hits, header.as_deref())?;
    eprintln!(
        "annotation-hub lookup {source}: {} known, {misses} new",
        hits.len()
    );
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
    let mut reader = ParquetReader::new(File::open(known)?);
    let stored_header = reader
        .get_metadata()?
        .key_value_metadata
        .as_deref()
        .unwrap_or_default()
        .iter()
        .find(|kv| kv.key == HEADER_KEY)
        .and_then(|kv| kv.value.clone());
    let frame = reader.finish()?;
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
    let mut known_rows = (0..frame.height()).map(|i| {
        (
            (c[i].clone(), p[i].clone(), r[i].clone(), a[i].clone()),
            v[i].clone(),
        )
    });
    let mut next_known = known_rows.next();
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
                    next_known = known_rows.next();
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
                    next_known = known_rows.next();
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
    let mut by_chrom: BTreeMap<String, Vec<(Key, String)>> = BTreeMap::new();
    let header;
    if source == "cadd" {
        let (h, mut rows, _) = cadd_rows(annotated, uses_chr(base)?)?;
        header = h;
        for line in base_vcf.records {
            let key = key_of(&line?)?;
            let value = rows.remove(&key).map(|r| r.join("\n")).unwrap_or_default();
            by_chrom
                .entry(key.0.clone())
                .or_default()
                .push((key, value));
        }
    } else {
        let mut fresh = read_vcf(annotated)?;
        let base_lines: HashSet<&String> = base_vcf.meta.iter().collect();
        header = fresh
            .meta
            .iter()
            .filter(|l| !base_lines.contains(l))
            .cloned()
            .collect::<Vec<_>>()
            .join("\n");
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
            by_chrom
                .entry(key.0.clone())
                .or_default()
                .push((key, value));
        }
    }
    let mut stored = 0;
    for (chrom, rows) in by_chrom {
        check_name("chromosome", &chrom)?;
        let dir = folder.join(&chrom);
        fs::create_dir_all(&dir)?;
        write_parquet(
            &dir.join(format!("part-{session}.parquet")),
            &rows,
            Some(&header),
        )?;
        stored += rows.len();
    }
    eprintln!(
        "annotation-hub add {source}: {stored} variants stored under {}",
        folder.display()
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
}
