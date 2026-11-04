//! Opt-in native SpliceVault_tsv_vcf_conversion.py. CrossMap remains upstream.
//!
//! Preserve the producer's stage order, chunk-local IDs, last-value dictionaries,
//! sequential publication and finally cleanup. A failed later chunk leaves the
//! completed prefix; unmapped CrossMap sidecars are deliberately not removed.
mod compression;
mod table;

use crate::cli::{Result, invalid};
use crate::pandas_table::{Cell, Column, ObjectValue, Table};
use crate::python_args::{parse_long_option_occurrences, require_options};
use crate::python_text::{decimal_digit_value, float_repr, python_int_is_valid, python_strip};
use polars_core::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

const HELP: &str = "usage: priva splicevault-convert --input_tsv TSV --output_tsv TSV --chain_file CHAIN --reference_fasta FASTA [--temp_dir ./temp] [--chunk_size 50000]\n\
Opt-in native SpliceVault table conversion for plain TSV input and pandas-inferred output compression (.gz, .bz2, .xz, .zst, .zip, .tar, .tar.gz, .tar.bz2, .tar.xz). Retains archive append behavior, including compressed-TAR failure on a second chunk. Runs the original upstream CrossMap vcf and bed commands, in order, once per original chunk. No download or installer default change.";

struct Options {
    input: PathBuf,
    output: PathBuf,
    chain: PathBuf,
    reference: PathBuf,
    temp: PathBuf,
    size: Option<usize>,
}

fn decimal_integer(value: &str) -> Option<String> {
    if !python_int_is_valid(value) {
        return None;
    }
    Some(
        python_strip(value)
            .chars()
            .filter(|c| *c != '_')
            .map(|c| decimal_digit_value(c).map_or(c, |d| char::from(b'0' + d as u8)))
            .collect(),
    )
}

// Python int() is not limited to the machine word used for VCF positions.
// BED receives its canonical decimal spelling and pos-1 even when upstream
// CrossMap will subsequently reject or leave that coordinate unmapped.
fn bed_integer(value: &str) -> Result<(String, String)> {
    let text =
        decimal_integer(value).ok_or_else(|| invalid("ValueError: invalid literal for int()"))?;
    let negative = text.starts_with('-');
    let digits = text.trim_start_matches(['-', '+']);
    if digits.len() > 4300 {
        return Err(invalid(
            "ValueError: Exceeds the limit for integer string conversion",
        ));
    }
    let digits = digits.trim_start_matches('0');
    let digits = if digits.is_empty() { "0" } else { digits };
    let negative = negative && digits != "0";
    let canonical = if negative {
        format!("-{digits}")
    } else {
        digits.to_owned()
    };
    if digits == "0" {
        return Ok((canonical, "-1".to_owned()));
    }
    let mut adjusted = digits.as_bytes().to_vec();
    if negative {
        let mut carry = true;
        for digit in adjusted.iter_mut().rev() {
            if *digit == b'9' {
                *digit = b'0';
            } else {
                *digit += 1;
                carry = false;
                break;
            }
        }
        if carry {
            adjusted.insert(0, b'1');
        }
        Ok((
            canonical,
            format!("-{}", String::from_utf8(adjusted).unwrap()),
        ))
    } else {
        for digit in adjusted.iter_mut().rev() {
            if *digit == b'0' {
                *digit = b'9';
            } else {
                *digit -= 1;
                break;
            }
        }
        let result = String::from_utf8(adjusted).unwrap();
        let result = result.trim_start_matches('0');
        Ok((
            canonical,
            if result.is_empty() {
                "0".to_owned()
            } else {
                result.to_owned()
            },
        ))
    }
}

fn arguments(args: &[String]) -> Result<Option<Options>> {
    let names = [
        "input_tsv",
        "output_tsv",
        "chain_file",
        "reference_fasta",
        "temp_dir",
        "chunk_size",
    ];
    // argparse applies an option's type immediately, before a later --help.
    for (index, token) in args.iter().enumerate() {
        if matches!(token.as_str(), "--help" | "-h") {
            break;
        }
        let (name, inline) = token
            .split_once('=')
            .map_or((token.as_str(), None), |(k, v)| (k, Some(v)));
        if name.starts_with("--") && "chunk_size".starts_with(&name[2..]) && !name[2..].is_empty() {
            if let Some(value) = inline.or_else(|| args.get(index + 1).map(String::as_str)) {
                if !value.starts_with("--") && decimal_integer(value).is_none() {
                    return Err(invalid(format!(
                        "argument --chunk_size: invalid int value: {value:?}"
                    )));
                }
            }
        }
    }
    let Some(occurrences) = parse_long_option_occurrences(args, &names, &[]).map_err(invalid)?
    else {
        return Ok(None);
    };
    for (name, value) in &occurrences {
        if name == "chunk_size" && decimal_integer(value).is_none() {
            return Err(invalid(format!(
                "argument --chunk_size: invalid int value: {value:?}"
            )));
        }
    }
    let values: HashMap<_, _> = occurrences.into_iter().collect();
    require_options(&values, &names[..4]).map_err(invalid)?;
    let size = decimal_integer(values.get("chunk_size").map_or("50000", String::as_str))
        .unwrap()
        .parse::<usize>()
        .ok()
        .filter(|n| *n > 0);
    Ok(Some(Options {
        input: PathBuf::from(&values["input_tsv"]),
        output: PathBuf::from(&values["output_tsv"]),
        chain: PathBuf::from(&values["chain_file"]),
        reference: PathBuf::from(&values["reference_fasta"]),
        temp: values
            .get("temp_dir")
            .map_or_else(|| PathBuf::from("./temp"), PathBuf::from),
        size,
    }))
}

fn uuid() -> Result<String> {
    let mut bytes = [0_u8; 16];
    File::open("/dev/urandom")?.read_exact(&mut bytes)?;
    bytes[6] = (bytes[6] & 15) | 0x40;
    bytes[8] = (bytes[8] & 63) | 0x80;
    let text: String = bytes.iter().map(|b| format!("{b:02x}")).collect();
    Ok(format!(
        "{}-{}-{}-{}-{}",
        &text[..8],
        &text[8..12],
        &text[12..16],
        &text[16..20],
        &text[20..]
    ))
}

fn scalar(chunk: &table::Chunk, column: usize, row: usize) -> String {
    if let Some(values) = chunk.unsigned.get(&column) {
        return values[row].to_string();
    }
    match chunk.table.cell(column, row) {
        Cell::Missing => "nan".into(),
        Cell::Text(s) => s.to_owned(),
        Cell::Int(n) => n.to_string(),
        Cell::Float(f) => float_repr(f),
        Cell::Bool(b) => if b { "True" } else { "False" }.into(),
    }
}

fn required(table: &Table, name: &str) -> Result<usize> {
    table
        .column_index(name)
        .ok_or_else(|| invalid(format!("KeyError: '{name}'")))
}

fn text_cell<'a>(chunk: &'a table::Chunk, col: usize, row: usize) -> Result<&'a str> {
    if chunk.unsigned.contains_key(&col) {
        return Err(invalid("TypeError: expected bytes or string, got integer"));
    }
    match chunk.table.cell(col, row) {
        Cell::Text(s) => Ok(s),
        _ => Err(invalid("TypeError: expected bytes or string")),
    }
}

fn record_pos(chunk: &table::Chunk, col: usize, row: usize) -> Result<i64> {
    let value = if let Some(values) = chunk.unsigned.get(&col) {
        i64::try_from(values[row])
            .map_err(|_| invalid("OverflowError: Python int too large to convert to C long"))?
    } else {
        match chunk.table.cell(col, row) {
            Cell::Int(n) => n,
            Cell::Bool(b) => i64::from(b),
            Cell::Float(n) if !n.is_finite() => {
                return Err(invalid(
                    "ValueError: cannot convert non-finite float to integer",
                ));
            }
            Cell::Float(n) if !(-((1_u64 << 63) as f64)..(1_u64 << 63) as f64).contains(&n) => {
                return Err(invalid(
                    "OverflowError: Python int too large to convert to C long",
                ));
            }
            Cell::Float(n) => n as i64,
            Cell::Missing => {
                return Err(invalid("ValueError: cannot convert float NaN to integer"));
            }
            Cell::Text(_) => return Err(invalid("TypeError: an integer is required")),
        }
    };
    if value < 1 {
        return Err(invalid("ValueError: Position must be positive"));
    }
    if value > i64::from(i32::MAX) {
        return Err(invalid("OverflowError: value too large to convert to int"));
    }
    Ok(value)
}

fn vcf_header(out: &mut impl Write, has_records: bool) -> Result<()> {
    writeln!(
        out,
        "##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"All filters passed\">"
    )?;
    for chrom in (1..=22)
        .map(|n| n.to_string())
        .chain(["X", "Y", "M"].map(str::to_owned))
    {
        writeln!(out, "##contig=<ID=chr{chrom}>")?;
    }
    // pysam new_record() registers END before it writes the first header.
    if has_records {
        writeln!(
            out,
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">"
        )?;
    }
    writeln!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
    Ok(())
}

fn write_vcf(path: &Path, chunk: &table::Chunk, ids: &[String], cols: [usize; 4]) -> Result<()> {
    let mut out = BufWriter::new(File::create(path)?);
    vcf_header(&mut out, !ids.is_empty())?;
    for (row, id) in ids.iter().enumerate() {
        let chrom = text_cell(chunk, cols[0], row)?;
        let known = chrom.strip_prefix("chr").is_some_and(|s| {
            matches!(s, "X" | "Y" | "M")
                || s.parse::<u8>()
                    .is_ok_and(|n| (1..=22).contains(&n) && n.to_string() == s)
        });
        if !known {
            return Err(invalid(format!(
                "ValueError: Invalid chromosome/contig {chrom}"
            )));
        }
        let pos = record_pos(chunk, cols[1], row)?;
        let reference = text_cell(chunk, cols[2], row)?.split(',').next().unwrap();
        if reference.is_empty() {
            return Err(invalid("ValueError: cannot set null allele"));
        }
        let alternate = text_cell(chunk, cols[3], row)?;
        if alternate.is_empty() {
            return Err(invalid("ValueError: cannot set null allele"));
        }
        let info = if alternate
            .split(',')
            .any(|allele| allele.starts_with('<') && allele.ends_with('>'))
        {
            format!("END={}", pos + reference.len() as i64 - 1)
        } else {
            ".".to_owned()
        };
        writeln!(
            out,
            "{chrom}\t{pos}\t{id}\t{reference}\t{alternate}\t.\t.\t{info}"
        )?;
    }
    out.flush()?;
    Ok(())
}

fn write_bed(path: &Path, chunk: &table::Chunk, ids: &[String]) -> Result<()> {
    let mut out = BufWriter::new(File::create(path)?);
    for (row, id) in ids.iter().enumerate() {
        let col = required(&chunk.table, "splicevault_site_loc")?;
        let site = text_cell(chunk, col, row)
            .map_err(|_| invalid("AttributeError: splice site has no attribute 'split'"))?;
        let fields: Vec<_> = site.split(':').collect();
        if fields.len() != 2 {
            return Err(invalid(
                "ValueError: splice site must unpack into chrom and position",
            ));
        }
        let (pos, start) = bed_integer(fields[1])?;
        // Each event row has its own site, so the lifted site is keyed by row.
        writeln!(out, "{}\t{start}\t{pos}\t{id}|{row}", fields[0])?;
    }
    out.flush()?;
    Ok(())
}

fn crossmap(
    options: &Options,
    vcf: &Path,
    lifted_vcf: &Path,
    bed: &Path,
    lifted_bed: &Path,
) -> Result<()> {
    let status = Command::new("CrossMap")
        .args(["vcf", "--chromid", "l"])
        .arg(&options.chain)
        .arg(vcf)
        .arg(&options.reference)
        .arg(lifted_vcf)
        .status()?;
    if !status.success() {
        return Err(invalid(format!(
            "CalledProcessError: CrossMap vcf returned {status}"
        )));
    }
    let unmapped = lifted_bed.to_string_lossy().replace(".bed", ".unmap.bed");
    let status = Command::new("CrossMap")
        .args(["bed", "--chromid", "l", "--unmap-file", &unmapped])
        .arg(&options.chain)
        .arg(bed)
        .arg(lifted_bed)
        .status()?;
    if !status.success() {
        return Err(invalid(format!(
            "CalledProcessError: CrossMap bed returned {status}"
        )));
    }
    Ok(())
}

struct Coordinate {
    chrom: String,
    pos: i64,
    reference: String,
    alternate: String,
}

fn mapped_vcf(path: &Path) -> Result<HashMap<String, Coordinate>> {
    let mut result = HashMap::new();
    let mut file = BufReader::new(File::open(path)?);
    let mut line = String::new();
    let mut header = false;
    while file.read_line(&mut line)? > 0 {
        if line.starts_with("#CHROM\t") {
            header = true;
        }
        if !line.starts_with('#') {
            if !header {
                return Err(invalid("ValueError: invalid VCF header"));
            }
            let fields: Vec<_> = line.trim_end_matches(['\r', '\n']).split('\t').collect();
            if fields.len() < 8 {
                return Err(invalid("OSError: unable to parse next VCF record"));
            }
            let pos = fields[1]
                .parse::<i64>()
                .map_err(|_| invalid("OSError: invalid VCF position"))?;
            if fields[4] == "." {
                return Err(invalid("TypeError: NoneType object is not subscriptable"));
            }
            result.insert(
                fields[2].to_owned(),
                Coordinate {
                    chrom: fields[0].to_owned(),
                    pos,
                    reference: fields[3].to_owned(),
                    alternate: fields[4].split(',').next().unwrap().to_owned(),
                },
            );
        }
        line.clear();
    }
    if !header {
        return Err(invalid("ValueError: invalid VCF header"));
    }
    Ok(result)
}

fn mapped_bed(path: &Path) -> Result<HashMap<String, String>> {
    let mut result = HashMap::new();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let fields: Vec<_> = python_strip(&line).split('\t').collect();
        if fields.len() != 4 {
            return Err(invalid(
                "ValueError: lifted BED row must unpack into four fields",
            ));
        }
        result.insert(fields[3].to_owned(), format!("{}:{}", fields[0], fields[2]));
    }
    Ok(result)
}

fn set_text(table: &mut Table, col: usize, row: usize, value: &str) -> Result<()> {
    match &mut table.columns[col] {
        Column::Text(values) => values[row] = Some(value.to_owned()),
        Column::Object(values) => values[row] = ObjectValue::Text(value.to_owned()),
        _ => {
            return Err(invalid(
                "Unexpected non-text dtype after a successful SpliceVault record write",
            ));
        }
    }
    Ok(())
}

fn set_pos(chunk: &mut table::Chunk, col: usize, row: usize, value: i64) -> Result<()> {
    match &mut chunk.table.columns[col] {
        Column::Int(values) => values[row] = value,
        Column::Float(values) => values[row] = value as f64,
        Column::Bool(values) => {
            let _ = values;
            let mut replacement = (0..chunk.table.rows)
                .map(|r| match chunk.table.cell(col, r) {
                    Cell::Bool(v) => ObjectValue::Bool(v),
                    _ => unreachable!(),
                })
                .collect::<Vec<_>>();
            replacement[row] = ObjectValue::Int(value);
            chunk.table.columns[col] = Column::Object(replacement);
        }
        Column::Object(values) => values[row] = ObjectValue::Int(value),
        _ => {
            return Err(invalid(
                "Unexpected position dtype after a successful SpliceVault record write",
            ));
        }
    }
    Ok(())
}

fn process(mut chunk: table::Chunk, options: &Options) -> Result<Table> {
    let id = uuid()?;
    let vcf = options.temp.join(format!("temp_grch38_{id}.vcf"));
    let lifted_vcf = options.temp.join(format!("temp_hg19_{id}.vcf"));
    let bed = options.temp.join(format!("temp_grch38_{id}.bed"));
    let lifted_bed = options.temp.join(format!("temp_hg19_{id}.bed"));
    let result = (|| {
        let mut cols = [0; 4];
        let mut ids = Vec::with_capacity(chunk.table.rows);
        if chunk.table.rows > 0 {
            for (target, name) in cols.iter_mut().zip(["#chrom", "pos", "ref", "alt"]) {
                *target = required(&chunk.table, name)?;
            }
            for row in 0..chunk.table.rows {
                ids.push(cols.map(|col| scalar(&chunk, col, row)).join("_"));
            }
        }
        chunk.table.set_column(
            "variant_id",
            Column::Text(ids.iter().cloned().map(Some).collect()),
        );
        write_vcf(&vcf, &chunk, &ids, cols)?;
        write_bed(&bed, &chunk, &ids)?;
        crossmap(options, &vcf, &lifted_vcf, &bed, &lifted_bed)?;
        let vcf_coords = mapped_vcf(&lifted_vcf)?;
        let bed_coords = mapped_bed(&lifted_bed)?;
        let mut unmapped = HashSet::new();
        for (row, id) in ids.iter().enumerate() {
            if let Some(value) = vcf_coords.get(id) {
                if chunk.duplicate_multi_index {
                    return Err(invalid(
                        "ValueError: Invalid call for scalar access (setting)!",
                    ));
                }
                for target in chunk.same_index[row].clone() {
                    set_text(&mut chunk.table, cols[0], target, &value.chrom)?;
                    set_pos(&mut chunk, cols[1], target, value.pos)?;
                    set_text(&mut chunk.table, cols[2], target, &value.reference)?;
                    set_text(&mut chunk.table, cols[3], target, &value.alternate)?;
                }
            } else {
                unmapped.insert(id.as_str());
            }
            if let Some(value) = bed_coords.get(&format!("{id}|{row}")) {
                if chunk.duplicate_multi_index {
                    return Err(invalid(
                        "ValueError: Invalid call for scalar access (setting)!",
                    ));
                }
                let col = required(&chunk.table, "splicevault_site_loc")?;
                for target in &chunk.same_index[row] {
                    set_text(&mut chunk.table, col, *target, value)?;
                }
            }
        }
        // Polars applies the original ID-membership filter, retaining every
        // duplicate occurrence and its source order after dictionary updates.
        let series = Series::new("variant_id".into(), &ids);
        let missing: Vec<_> = unmapped.into_iter().collect();
        let mask = crate::pandas_polars::isin_text(&series, &missing)?;
        let keep =
            BooleanChunked::from_iter_values("keep".into(), mask.into_no_null_iter().map(|v| !v));
        let indices = Series::new(
            "row".into(),
            (0..chunk.table.rows as u64).collect::<Vec<_>>(),
        )
        .filter(&keep)?;
        let rows: Vec<_> = indices
            .u64()?
            .into_no_null_iter()
            .map(|v| v as usize)
            .collect();
        let mut updated = chunk.table.take_rows(&rows);
        let id_column = required(&updated, "variant_id")?;
        updated.names.remove(id_column);
        updated.columns.remove(id_column);
        Ok(updated)
    })();
    // Python finally runs in exactly this order; its first cleanup failure
    // replaces a prior exception and prevents later cleanup operations.
    for path in [&vcf, &lifted_vcf, &bed, &lifted_bed] {
        if path.exists() {
            std::fs::remove_file(path)?;
        }
    }
    result
}

fn raw_header(path: &Path) -> Result<String> {
    let mut reader = BufReader::new(File::open(path)?);
    let mut raw = Vec::new();
    let mut pending_utf8 = Vec::new();
    loop {
        let buffer = reader.fill_buf()?;
        if buffer.is_empty() {
            break;
        }
        // TextIOWrapper.readline decodes its complete 8192-byte read-ahead
        // block before returning the first line. Invalid later bytes in that
        // same block therefore fail before the output is initially truncated.
        let mut decode = std::mem::take(&mut pending_utf8);
        decode.extend_from_slice(buffer);
        if let Err(error) = std::str::from_utf8(&decode) {
            if error.error_len().is_some() {
                return Err(invalid(
                    "UnicodeDecodeError: input header read-ahead is not UTF-8",
                ));
            }
            pending_utf8.extend_from_slice(&decode[error.valid_up_to()..]);
        }
        let end = buffer.iter().position(|b| matches!(*b, b'\n' | b'\r'));
        let length = end.map_or(buffer.len(), |i| i + 1);
        raw.extend_from_slice(&buffer[..length]);
        reader.consume(length);
        if end.is_some() {
            break;
        }
    }
    Ok(python_strip(
        std::str::from_utf8(&raw)
            .map_err(|_| invalid("UnicodeDecodeError: input header is not UTF-8"))?,
    )
    .to_owned())
}

fn publish(path: &Path, table: &Table, first: bool) -> Result<()> {
    let name = path.to_string_lossy().to_lowercase();
    if let Some(format) = compression::Format::infer(&name) {
        return compression::publish(path, table, first, format);
    }
    let file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(first)
        .append(!first)
        .open(path)?;
    let mut out = BufWriter::new(file);
    if name.ends_with(".gz") {
        let name = path.file_name().unwrap_or_default().to_string_lossy();
        let stem = name.strip_suffix(".gz").unwrap_or(&name);
        let latin: Option<Vec<_>> = stem
            .chars()
            .map(|c| u8::try_from(u32::from(c)).ok())
            .collect();
        let mut gzip = crate::gzip_writer::GzipWriter::new(
            out,
            crate::gzip_writer::GzipOptions {
                level: 9,
                mtime: SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs() as u32,
                filename: latin,
            },
        )?;
        table::write(&mut gzip, table, first)?;
        gzip.finish()?;
    } else {
        table::write(&mut out, table, first)?;
        out.flush()?;
    }
    Ok(())
}

fn run(options: Options) -> Result<()> {
    std::fs::create_dir_all(&options.temp)?;
    let header = raw_header(&options.input)?;
    {
        let mut out = File::create(&options.output)?;
        writeln!(out, "{header}")?;
    }
    let size = options
        .size
        .ok_or_else(|| invalid("ValueError: 'chunksize' must be an integer >=1"))?;
    let mut reader = table::Reader::new(
        BufReader::with_capacity(256 * 1024, File::open(&options.input)?),
        size,
    )?;
    let mut count = 0;
    while let Some(chunk) = reader.next()? {
        count += 1;
        eprintln!("Processing chunk {count}...");
        let updated = process(chunk, &options)?;
        publish(&options.output, &updated, count == 1)?;
        eprintln!("Completed chunk {count}");
    }
    Ok(())
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(options) = arguments(args)? else {
        println!("{HELP}");
        return Ok(());
    };
    match run(options) {
        Ok(()) => Ok(()),
        Err(error) => {
            eprintln!("priva splicevault-convert: {error}");
            std::process::exit(1)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bed_positions_use_python_decimal_not_machine_word_arithmetic() {
        for (input, end, start) in [
            (" ２_１ ", "21", "20"),
            ("-0", "0", "-1"),
            ("1000", "1000", "999"),
            (
                "-999999999999999999999999999999",
                "-999999999999999999999999999999",
                "-1000000000000000000000000000000",
            ),
        ] {
            assert_eq!(bed_integer(input).unwrap(), (end.into(), start.into()));
        }
        assert!(bed_integer("21.0").is_err());
    }
    #[test]
    fn invalid_integer_precedes_later_help() {
        let args = ["--chunk_size", "bad", "--help"].map(str::to_owned);
        assert!(arguments(&args).is_err());
        assert!(arguments(&["--help".into()]).unwrap().is_none());
    }
}
