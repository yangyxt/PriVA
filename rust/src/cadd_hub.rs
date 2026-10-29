//! Rust version of `scripts/cadd_hub_util.py`, the CADD score cache ("hub")
//! that annotation runs share.
//!
//! - match: split a VCF's variants (CHROM POS REF ALT, no header) into the hub
//!   rows that already score them (in hub order) and the positions still to
//!   score (in input order);
//! - update: append newly computed CADD rows whose variant the hub lacks;
//! - merge: join the cached and the newly computed rows into one table.
//!
//! Files are read and written the way the Python script's polars calls do
//! (checked against polars 1.27.1): a value spelled as one of `NULL_VALUES` is
//! missing and is written back empty; chromosome names in hub and CADD tables
//! (not in the position file) go through `data/liftover/GRC_to_ucsc.contig.map.txt`;
//! a value holding a quote, tab or line break is written quoted; a variant with
//! a missing CHROM, POS, REF or ALT matches nothing and is left out of both
//! match outputs; output is plain text even when the name ends in `.gz`. Before
//! replacing the hub or writing the merge, the Python script reads its own
//! output back with polars' type inference and gives up when that fails; the
//! same inference is applied here. Exit codes follow the Python script: 1 when
//! a step fails or merge finds no rows, else 0.
use crate::cli::{Result, invalid};
use crate::vcf_text::open_reader;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};

const USAGE: &str = "usage: priva cadd-hub match --hub HUB --pos POS --covered COVERED --uncovered UNCOVERED\n       priva cadd-hub update --new NEW --hub HUB\n       priva cadd-hub merge --covered COVERED --new NEW --output OUTPUT";
const OUTPUT_HEADER: [&str; 6] = ["#Chrom", "Pos", "Ref", "Alt", "FeatureID", "PHRED"];
const COLUMNS: [&str; 6] = ["Chrom", "Pos", "Ref", "Alt", "FeatureID", "PHRED"];
const NULL_VALUES: [&str; 10] = ["NA", "na", "NaN", "nan", "", ".", ",", ";", "NAN", "-"];
const CONTIG_MAP: &str = include_str!("../../data/liftover/GRC_to_ucsc.contig.map.txt");
/// Rows polars samples to choose column types: the position file, and the
/// read-back check of written files.
const POSITION_SAMPLE: usize = 100_000;
const CHECK_SAMPLE: usize = 10_000;

type Cell = Option<String>;

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
        options.insert(option.as_str(), value.as_str());
    }
    let allowed: &[&str] = match command.as_str() {
        "match" => &["--hub", "--pos", "--covered", "--uncovered", "--threads"],
        "update" => &["--new", "--hub", "--threads"],
        "merge" => &["--covered", "--new", "--output"],
        _ => return Err(invalid(format!("unknown command {command:?}\n{USAGE}"))),
    };
    if let Some(option) = options.keys().find(|option| !allowed.contains(*option)) {
        return Err(invalid(format!("unknown option {option}\n{USAGE}")));
    }
    let path = |name: &str| -> Result<PathBuf> {
        options
            .get(name)
            .map(|value| PathBuf::from(*value))
            .ok_or_else(|| invalid(format!("{name} is required\n{USAGE}")))
    };
    let outcome = match command.as_str() {
        "match" => match_variants(
            &path("--hub")?,
            &path("--pos")?,
            &path("--covered")?,
            &path("--uncovered")?,
        )
        .map(|_| true),
        "update" => update_hub(&path("--new")?, &path("--hub")?).map(|_| true),
        _ => merge_results(&path("--covered")?, &path("--new")?, &path("--output")?),
    };
    match outcome {
        Ok(true) => Ok(()),
        Ok(false) => std::process::exit(1),
        Err(error) => {
            eprintln!("priva cadd-hub {command}: {error}");
            std::process::exit(1)
        }
    }
}

/// Tab-separated fields; a field opening with a quote ends at the next lone
/// quote and `""` inside it is one quote, as polars reads them.
fn split_fields(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut rest = line;
    loop {
        if let Some(quoted) = rest.strip_prefix('"') {
            let (mut value, mut end) = (String::new(), quoted.len());
            let mut chars = quoted.char_indices().peekable();
            while let Some((index, character)) = chars.next() {
                if character != '"' {
                    value.push(character);
                } else if chars.peek().map(|(_, next)| *next) == Some('"') {
                    value.push('"');
                    chars.next();
                } else {
                    end = index + 1;
                    break;
                }
            }
            let after = &quoted[end..];
            let (tail, next) = after
                .split_once('\t')
                .map_or((after, None), |(a, b)| (a, Some(b)));
            value.push_str(tail);
            fields.push(value);
            match next {
                Some(next) => rest = next,
                None => return fields,
            }
        } else {
            match rest.split_once('\t') {
                Some((field, next)) => {
                    fields.push(field.to_owned());
                    rest = next;
                }
                None => {
                    fields.push(rest.to_owned());
                    return fields;
                }
            }
        }
    }
}

fn cell(value: String) -> Cell {
    (!NULL_VALUES.contains(&value.as_str())).then_some(value)
}

fn cell_text(value: &str) -> Cell {
    (!NULL_VALUES.contains(&value)).then(|| value.to_owned())
}

/// Calls `visit` with every non-empty line, without its line break; a file
/// without one fails like polars' "empty CSV".
fn each_line(path: &Path, mut visit: impl FnMut(&str) -> Result<()>) -> Result<()> {
    let mut reader = open_reader(path)?;
    let (mut buffer, mut seen) = (Vec::new(), false);
    loop {
        buffer.clear();
        if reader.read_until(b'\n', &mut buffer)? == 0 {
            break;
        }
        while matches!(buffer.last(), Some(b'\n' | b'\r')) {
            buffer.pop();
        }
        if buffer.is_empty() {
            continue;
        }
        let line = std::str::from_utf8(&buffer)
            .map_err(|_| invalid(format!("{}: not UTF-8 text", path.display())))?;
        seen = true;
        visit(line)?;
    }
    if seen {
        Ok(())
    } else {
        Err(invalid(format!("{}: empty CSV", path.display())))
    }
}

/// The wanted columns of a table with a header: `present[i]` says whether
/// `wanted[i]` exists, and absent columns hold `None` in every row.
struct Table {
    present: Vec<bool>,
    rows: Vec<Vec<Cell>>,
}

impl Table {
    fn complete(&self) -> bool {
        self.present.iter().all(|present| *present)
    }
}

/// Reads a hub or CADD table: `#Chrom` is renamed `Chrom` and mapped to UCSC
/// names; a table without a chromosome column fails, as in the Python script.
fn read_table(path: &Path, wanted: &[&str]) -> Result<Table> {
    let contigs = contig_map()?;
    let (mut header, mut width, mut chrom) = (false, 0, 0);
    let (mut index, mut rows, mut number) = (Vec::new(), Vec::new(), 1);
    each_line(path, |line| {
        if !header {
            let mut names = split_fields(line);
            if let Some(name) = names.iter_mut().find(|name| name.as_str() == "#Chrom") {
                *name = "Chrom".to_owned();
            }
            chrom = names
                .iter()
                .position(|name| name.as_str() == "Chrom")
                .ok_or_else(|| invalid(format!("{}: no Chrom column", path.display())))?;
            index = wanted
                .iter()
                .map(|name| names.iter().position(|column| column.as_str() == *name))
                .collect();
            (header, width) = (true, names.len());
            return Ok(());
        }
        number += 1;
        // Only a line holding a quote needs polars' quoted-field reading.
        let quoted;
        let fields: Vec<&str> = if line.contains('"') {
            quoted = split_fields(line);
            quoted.iter().map(String::as_str).collect()
        } else {
            line.split('\t').collect()
        };
        if fields.len() > width {
            return Err(invalid(format!(
                "{}: row {number} has more fields than the header",
                path.display()
            )));
        }
        let value = |column: usize| -> Cell {
            let raw = *fields.get(column)?;
            if column == chrom {
                cell_text(contigs.get(raw).copied().unwrap_or(raw))
            } else {
                cell_text(raw)
            }
        };
        let row: Vec<Cell> = index.iter().map(|column| column.and_then(value)).collect();
        rows.push(row);
        Ok(())
    })?;
    Ok(Table {
        present: index.iter().map(Option::is_some).collect(),
        rows,
    })
}

fn contig_map() -> Result<HashMap<&'static str, &'static str>> {
    let mut map = HashMap::new();
    for line in CONTIG_MAP
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
    {
        match line.split(' ').collect::<Vec<_>>()[..] {
            [key, value] => map.insert(key, value),
            _ => {
                return Err(invalid(format!(
                    "contig map line {line:?} is not two names"
                )));
            }
        };
    }
    Ok(map)
}

/// The position file (CHROM POS REF ALT, no header). POS becomes an integer,
/// as polars infers it, when every sampled value is one.
fn read_positions(path: &Path) -> Result<Vec<[Cell; 4]>> {
    let (mut rows, mut number) = (Vec::new(), 0);
    each_line(path, |line| {
        number += 1;
        let fields = split_fields(line);
        if fields.len() < 4 {
            return Err(invalid(format!(
                "{}: row {number} has fewer than 4 columns",
                path.display()
            )));
        }
        let mut fields = fields.into_iter().take(4).map(cell);
        rows.push(std::array::from_fn(|_| fields.next().flatten()));
        Ok(())
    })?;
    let integer = rows
        .iter()
        .take(POSITION_SAMPLE)
        .filter_map(|row| row[1].as_deref())
        .all(|pos| kind(pos) == Kind::Integer)
        && rows
            .iter()
            .take(POSITION_SAMPLE)
            .any(|row| row[1].is_some());
    if integer {
        for (number, row) in rows.iter_mut().enumerate() {
            if let Some(pos) = row[1].as_mut() {
                *pos = pos
                    .parse::<i64>()
                    .map_err(|_| {
                        invalid(format!(
                            "{}: row {}: POS is not an integer",
                            path.display(),
                            number + 1
                        ))
                    })?
                    .to_string();
            }
        }
    }
    Ok(rows)
}

fn key(chrom: &Cell, pos: &Cell, reference: &Cell, alt: &Cell) -> Option<String> {
    Some(format!(
        "{}_{}_{}_{}",
        chrom.as_deref()?,
        pos.as_deref()?,
        reference.as_deref()?,
        alt.as_deref()?
    ))
}

fn row_key(row: &[Cell]) -> Option<String> {
    key(&row[0], &row[1], &row[2], &row[3])
}

fn nonempty_file(path: &Path) -> bool {
    fs::metadata(path).is_ok_and(|meta| meta.len() > 0)
}

fn format_cell(value: &Cell) -> String {
    match value {
        None => String::new(),
        Some(text) if text.contains(['"', '\t', '\n', '\r']) => {
            format!("\"{}\"", text.replace('"', "\"\""))
        }
        Some(text) => text.clone(),
    }
}

fn write_rows<'a>(
    path: &Path,
    header: &[&str],
    rows: impl Iterator<Item = &'a Vec<Cell>>,
) -> Result<()> {
    let mut out = BufWriter::new(File::create(path)?);
    writeln!(out, "{}", header.join("\t"))?;
    for row in rows {
        let fields: Vec<String> = row.iter().map(format_cell).collect();
        writeln!(out, "{}", fields.join("\t"))?;
    }
    out.flush()?;
    Ok(())
}

/// Writes through a temporary file renamed into place, as the Python script does.
fn replace_with<'a>(
    path: &Path,
    header: &[&str],
    rows: impl Iterator<Item = &'a Vec<Cell>>,
) -> Result<()> {
    let temporary = PathBuf::from(format!("{}.temp.{}", path.display(), std::process::id()));
    write_rows(&temporary, header, rows)?;
    fs::rename(&temporary, path)?;
    Ok(())
}

fn write_uncovered<'a>(path: &Path, rows: impl Iterator<Item = &'a [Cell; 4]>) -> Result<()> {
    let mut out = BufWriter::new(File::create(path)?);
    for row in rows {
        // Python prints a missing value as None.
        let text = |value: &Cell| value.clone().unwrap_or_else(|| "None".to_owned());
        writeln!(out, "{}\t{}", text(&row[0]), text(&row[1]))?;
    }
    out.flush()?;
    Ok(())
}

fn match_variants(hub: &Path, pos: &Path, covered: &Path, uncovered: &Path) -> Result<()> {
    let variants = read_positions(pos)?;
    let none: [Vec<Cell>; 0] = [];
    let table = if nonempty_file(hub) {
        Some(read_table(hub, &COLUMNS)?).filter(Table::complete)
    } else {
        None
    };
    let Some(table) = table else {
        eprintln!(
            "CADD hub {} missing, empty or incomplete: all variants are uncovered",
            hub.display()
        );
        write_uncovered(uncovered, variants.iter())?;
        return write_rows(covered, &OUTPUT_HEADER, none.iter());
    };
    let wanted: HashSet<String> = variants.iter().filter_map(|row| row_key(row)).collect();
    let cached: HashSet<String> = table.rows.iter().filter_map(|row| row_key(row)).collect();
    let hits = table
        .rows
        .iter()
        .filter(|row| row_key(row).is_some_and(|key| wanted.contains(&key)));
    replace_with(covered, &OUTPUT_HEADER, hits)?;
    let misses: Vec<&[Cell; 4]> = variants
        .iter()
        .filter(|row| row_key(row.as_slice()).is_some_and(|key| !cached.contains(&key)))
        .collect();
    let temporary = PathBuf::from(format!(
        "{}.temp.{}",
        uncovered.display(),
        std::process::id()
    ));
    write_uncovered(&temporary, misses.iter().copied())?;
    fs::rename(&temporary, uncovered)?;
    eprintln!(
        "Found {} of {} distinct input variants in the CADD hub",
        wanted.iter().filter(|key| cached.contains(*key)).count(),
        wanted.len()
    );
    Ok(())
}

fn update_hub(new: &Path, hub: &Path) -> Result<()> {
    if !nonempty_file(new) {
        eprintln!(
            "New CADD file {} is empty or missing: no update",
            new.display()
        );
        return Ok(());
    }
    if let Some(folder) = hub.parent().filter(|folder| !folder.as_os_str().is_empty()) {
        fs::create_dir_all(folder)?;
    }
    let fresh = read_table(new, &COLUMNS)?;
    if !fresh.complete() {
        eprintln!(
            "New CADD file {} lacks a required column: hub not updated",
            new.display()
        );
        return Ok(());
    }
    if !nonempty_file(hub) {
        return replace_with(hub, &OUTPUT_HEADER, fresh.rows.iter());
    }
    let existing = read_table(hub, &COLUMNS)?;
    if !existing.complete() {
        eprintln!(
            "CADD hub {} lacks a required column: not updated",
            hub.display()
        );
        return Ok(());
    }
    let cached: HashSet<String> = existing
        .rows
        .iter()
        .filter_map(|row| row_key(row))
        .collect();
    let added: Vec<&Vec<Cell>> = fresh
        .rows
        .iter()
        .filter(|row| row_key(row).is_some_and(|key| !cached.contains(&key)))
        .collect();
    if added.is_empty() {
        eprintln!("No new variants to add to the CADD hub {}", hub.display());
        return Ok(());
    }
    let rows: Vec<&Vec<Cell>> = existing.rows.iter().chain(added.iter().copied()).collect();
    if let Err(error) = reads_back(&rows) {
        eprintln!("Error verifying the updated hub: {error}; hub not updated");
        return Ok(());
    }
    replace_with(hub, &OUTPUT_HEADER, rows.into_iter())?;
    eprintln!(
        "Added {} variants to the CADD hub {}",
        added.len(),
        hub.display()
    );
    Ok(())
}

fn merge_results(covered: &Path, new: &Path, output: &Path) -> Result<bool> {
    if let Some(folder) = output
        .parent()
        .filter(|folder| !folder.as_os_str().is_empty())
    {
        fs::create_dir_all(folder)?;
    }
    let none: [Vec<Cell>; 0] = [];
    let mut tables = Vec::new();
    for path in [covered, new] {
        if !nonempty_file(path) {
            continue;
        }
        match read_table(path, &COLUMNS) {
            Ok(table) if !table.rows.is_empty() => tables.push(table),
            Ok(_) => {}
            Err(error) => eprintln!("Could not read {}: {error}; skipped", path.display()),
        }
    }
    if tables.is_empty() {
        eprintln!("No CADD scores available to merge");
        write_rows(output, &OUTPUT_HEADER, none.iter())?;
        return Ok(false);
    }
    // Columns present in any input, in the standard order; rows of an input
    // lacking one hold it missing.
    let keep: Vec<usize> = (0..COLUMNS.len())
        .filter(|&column| tables.iter().any(|table| table.present[column]))
        .collect();
    let header: Vec<&str> = keep.iter().map(|&column| OUTPUT_HEADER[column]).collect();
    let rows: Vec<Vec<Cell>> = tables
        .iter()
        .flat_map(|table| table.rows.iter())
        .map(|row| keep.iter().map(|&column| row[column].clone()).collect())
        .collect();
    let all: Vec<&Vec<Cell>> = rows.iter().collect();
    if let Err(error) = reads_back(&all) {
        eprintln!("Error verifying the merged table: {error}");
        write_rows(output, &OUTPUT_HEADER, none.iter())?;
        return Ok(false);
    }
    replace_with(output, &header, rows.iter())?;
    eprintln!("Merged {} CADD rows into {}", rows.len(), output.display());
    Ok(true)
}

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
enum Kind {
    Boolean,
    Integer,
    Float,
    Text,
}

/// The type polars infers for one value.
fn kind(value: &str) -> Kind {
    let digits = |text: &str| !text.is_empty() && text.bytes().all(|b| b.is_ascii_digit());
    let unsigned = value.strip_prefix('-').unwrap_or(value);
    if value.eq_ignore_ascii_case("true") || value.eq_ignore_ascii_case("false") {
        return Kind::Boolean;
    }
    if digits(unsigned) {
        return Kind::Integer;
    }
    let signed = value.strip_prefix(['-', '+']).unwrap_or(value);
    let (mantissa, exponent) = match signed.split_once(['e', 'E']) {
        Some((mantissa, exponent)) => (mantissa, Some(exponent)),
        None => (signed, None),
    };
    let exponent_ok = exponent
        .is_none_or(|exponent| digits(exponent.strip_prefix(['-', '+']).unwrap_or(exponent)));
    let mantissa_ok = match mantissa.split_once('.') {
        Some((whole, fraction)) => {
            (whole.is_empty() || digits(whole))
                && (digits(fraction)
                    || (fraction.is_empty() && exponent.is_none() && digits(whole)))
        }
        None => digits(mantissa) && exponent.is_some(),
    };
    if (mantissa_ok && exponent_ok) || signed == "inf" {
        Kind::Float
    } else {
        Kind::Text
    }
}

/// Applies polars' read-back of a written table: column types come from the
/// first `CHECK_SAMPLE` rows, and a later value that does not parse as its
/// column's type makes the read fail.
fn reads_back(rows: &[&Vec<Cell>]) -> Result<()> {
    let columns = rows.first().map_or(0, |row| row.len());
    for column in 0..columns {
        let mut kinds = HashSet::new();
        for row in rows.iter().take(CHECK_SAMPLE) {
            if let Some(value) = &row[column] {
                kinds.insert(kind(value));
            }
        }
        let chosen = if kinds.is_empty() || kinds.contains(&Kind::Text) {
            Kind::Text
        } else if kinds.len() == 1 {
            *kinds.iter().next().unwrap()
        } else if !kinds.contains(&Kind::Boolean) {
            Kind::Float
        } else {
            Kind::Text
        };
        for (number, row) in rows.iter().enumerate().skip(CHECK_SAMPLE) {
            let Some(value) = &row[column] else { continue };
            let fits = match chosen {
                Kind::Text => true,
                Kind::Boolean => kind(value) == Kind::Boolean,
                Kind::Integer => value.parse::<i64>().is_ok(),
                Kind::Float => value.parse::<f64>().is_ok(),
            };
            if !fits {
                return Err(invalid(format!(
                    "row {}: {value:?} does not parse as the column's inferred {chosen:?}",
                    number + 1
                )));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fields_are_split_like_polars() {
        assert_eq!(split_fields("a\t\"b\"\"c\"\t"), ["a", "b\"c", ""]);
        assert_eq!(split_fields("x"), ["x"]);
        assert_eq!(format_cell(&Some("a\"b".into())), "\"a\"\"b\"");
        assert_eq!(format_cell(&cell("-".into())), "");
    }

    #[test]
    fn value_kinds_follow_polars_inference() {
        assert_eq!(kind("12"), Kind::Integer);
        assert_eq!(kind("-3"), Kind::Integer);
        assert_eq!(kind("12.5"), Kind::Float);
        assert_eq!(kind("1e-3"), Kind::Float);
        assert_eq!(kind(".5"), Kind::Float);
        assert_eq!(kind("True"), Kind::Boolean);
        assert_eq!(kind("ENST1"), Kind::Text);
        assert_eq!(kind("chr1"), Kind::Text);
    }

    #[test]
    fn missing_alleles_have_no_key() {
        let row = [
            Some("chr1".into()),
            Some("5".into()),
            Some("A".into()),
            cell(".".into()),
        ];
        assert_eq!(row_key(&row), None);
        assert!(contig_map().unwrap().get("MT").is_some());
    }
}
