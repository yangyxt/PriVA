//! Native `prepare_clingen_map.py`: pandas-compatible table inference, a
//! source-ordered ClinVar ID map, and direct schema-10 publication.
//!
//! Polars performs the dataframe filters and last-row reduction. Integer
//! conversion is explicit because pandas/NumPy and Polars differ for floats,
//! unsigned integers and object strings. No scientific evidence is inferred.
use crate::cache_convert::Provenance;
use crate::cache_convert::clingen_map::write_entries;
use crate::cli::{Result, invalid};
use crate::pandas_table::{self, Cell, Column as PandasColumn};
use crate::python_args::{Parsed, parse_long_options};
use crate::python_text::{decimal_digit_value, is_python_space, python_int_is_valid, python_strip};
use crate::rust_cache::ensure_distinct_paths;
use polars_core::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Read};
use std::path::{Path, PathBuf};
use std::time::Instant;

const ID: &str = "ClinVar Variation Id";
const CODE: &str = "Applied Evidence Codes (Met)";
const USAGE: &str = "usage: priva clingen-map build --clinvar-vcf VCF --clingen-tab TSV --output RKYV\n\
Build a schema-10 ClinGen map directly from the same inputs as prepare_clingen_map.py.\n\
The VCF may be plain text or BGZF. Evidence codes remain text or missing;\n\
duplicate ClinVar IDs and duplicate final variants keep their last source value.";

/// The shared reader deliberately refuses uint64 inference. This producer
/// needs that dtype for pandas' explicit `astype(int)` wrapping operation.
enum ProducerColumn {
    Pandas(PandasColumn),
    Unsigned(Vec<u64>),
}

impl ProducerColumn {
    fn len(&self) -> usize {
        match self {
            Self::Pandas(values) => values.len(),
            Self::Unsigned(values) => values.len(),
        }
    }

    fn series(&self, name: &str) -> Result<Series> {
        match self {
            Self::Pandas(values) => crate::pandas_polars::to_series(values, name),
            Self::Unsigned(values) => Ok(Series::new(name.into(), values)),
        }
    }

    fn integer_text(&self, row: usize) -> Result<String> {
        let number = match self {
            Self::Unsigned(values) => values[row] as i64,
            Self::Pandas(values) => match values.cell(row) {
                Cell::Int(number) => number,
                Cell::Bool(value) => i64::from(value),
                Cell::Float(number) => {
                    if !number.is_finite() {
                        return Err(invalid(
                            "IntCastingNaNError: cannot convert non-finite ClinVar IDs to integer",
                        ));
                    }
                    // NumPy's float64 -> platform int64 conversion truncates
                    // finite values, with INT64_MIN for out-of-range values.
                    if !(-((1_u64 << 63) as f64)..(1_u64 << 63) as f64).contains(&number) {
                        i64::MIN
                    } else {
                        number as i64
                    }
                }
                Cell::Text(text) => object_integer(text)?,
                Cell::Missing => {
                    return Err(invalid("a missing ClinVar ID reached the integer cast"));
                }
            },
        };
        Ok(number.to_string())
    }

    fn code(&self, row: usize) -> Result<Option<String>> {
        match self {
            Self::Pandas(values) => match values.cell(row) {
                Cell::Text(text) => Ok(Some(text.to_owned())),
                Cell::Missing => Ok(None),
                _ => Err(invalid(
                    "a surviving ClinGen code is not text or float NaN; schema 10 refuses it",
                )),
            },
            Self::Unsigned(_) => Err(invalid(
                "a surviving ClinGen code is uint64, not text or float NaN; schema 10 refuses it",
            )),
        }
    }
}

/// `object.astype(int)` uses Python's integer parser, including Unicode
/// decimal digits and underscores. It must not accept an object text `1.0`.
fn object_integer(text: &str) -> Result<i64> {
    if !python_int_is_valid(text) {
        return Err(invalid(
            "ValueError: invalid literal for int() in ClinVar Variation Id",
        ));
    }
    let ascii: String = text
        .chars()
        .map(|c| {
            if let Some(digit) = decimal_digit_value(c) {
                char::from(b'0' + digit as u8)
            } else if !c.is_ascii() && is_python_space(c) {
                ' '
            } else {
                c
            }
        })
        .filter(|c| *c != '_')
        .collect();
    ascii
        .trim_matches(|c| matches!(c, ' ' | '\t' | '\n' | '\r' | '\u{b}' | '\u{c}'))
        .parse::<i64>()
        .map_err(|_| invalid("OverflowError: Python int too large to convert to C long"))
}

fn infer(tokens: Vec<String>) -> Result<ProducerColumn> {
    match pandas_table::infer_column(tokens.clone(), true) {
        Ok(column) => Ok(ProducerColumn::Pandas(column)),
        Err(error)
            if error
                .to_string()
                .contains("outside the signed 64-bit integer range") =>
        {
            let mut numbers = Vec::new();
            let (mut unsigned_range, mut has_unsigned, mut any_negative) = (true, false, false);
            for token in tokens.iter().filter(|token| !pandas_table::is_na(token)) {
                let token = token
                    .trim_matches(|c| matches!(c, ' ' | '\t' | '\n' | '\r' | '\u{b}' | '\u{c}'));
                // pandas marks a leading minus before looking at any
                // following characters, including `-x` and huge values.
                if token.starts_with('-') {
                    any_negative = true;
                    continue;
                }
                let digits = token.strip_prefix('+').unwrap_or(token);
                let number = (!digits.is_empty()
                    && digits.bytes().all(|byte| byte.is_ascii_digit()))
                .then(|| digits.parse::<u64>().ok())
                .flatten();
                let Some(number) = number else {
                    unsigned_range = false;
                    break;
                };
                has_unsigned |= number > i64::MAX as u64;
                numbers.push(number);
            }
            if unsigned_range && has_unsigned {
                let any_missing = tokens.iter().any(|token| pandas_table::is_na(token));
                if !any_missing && !any_negative {
                    return Ok(ProducerColumn::Unsigned(numbers));
                }
                // pandas cancels uint64 conversion if it sees NA or a
                // negative integer: the original NA spellings stay text.
                return Ok(ProducerColumn::Pandas(PandasColumn::Text(
                    tokens.into_iter().map(Some).collect(),
                )));
            }
            // A true integer overflow, e.g. > UINT64_MAX, becomes object
            // text with normal NA handling, rather than rounded float64.
            Ok(ProducerColumn::Pandas(PandasColumn::Text(
                tokens
                    .into_iter()
                    .map(|token| {
                        if pandas_table::is_na(&token) {
                            None
                        } else {
                            Some(token)
                        }
                    })
                    .collect(),
            )))
        }
        Err(error) => Err(error),
    }
}

/// Read the two selected columns without losing their original NA spelling.
/// pandas may infer an implicit index from the first data row's extra leading
/// fields; that index is discarded by the producer's later `set_index`.
fn parse_curated(text: &str) -> Result<(ProducerColumn, ProducerColumn)> {
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    let mut names = None;
    let mut selected = None;
    let mut width = None;
    let mut ids = Vec::new();
    let mut codes = Vec::new();
    pandas_table::tokenize(text, None, |record| {
        if names.is_none() {
            let header = pandas_table::column_names(record);
            let id = header
                .iter()
                .position(|name| name == ID)
                .ok_or_else(|| invalid(format!("KeyError: '{ID}'")))?;
            let code = header
                .iter()
                .position(|name| name == CODE)
                .ok_or_else(|| invalid(format!("KeyError: '{CODE}'")))?;
            selected = Some((id, code));
            names = Some(header);
            return Ok(());
        }
        let named = names.as_ref().unwrap().len();
        let expected = *width.get_or_insert(record.len().max(named));
        if record.len() > expected {
            return Err(invalid(
                "ParserError: a ClinGen table row has too many fields",
            ));
        }
        let leading = expected - named;
        let (id, code) = selected.unwrap();
        ids.push(record.get(leading + id).cloned().unwrap_or_default());
        codes.push(record.get(leading + code).cloned().unwrap_or_default());
        Ok(())
    })?;
    if names.is_none() {
        return Err(invalid(
            "EmptyDataError: no columns to parse from the table",
        ));
    }
    Ok((infer(ids)?, infer(codes)?))
}

#[derive(Clone, Copy)]
enum FormatType {
    Integer,
    Float,
    Text,
}

fn format_definition(line: &[u8]) -> Option<(String, FormatType)> {
    let line = std::str::from_utf8(line).ok()?;
    let body = line.strip_prefix("##FORMAT=<")?.strip_suffix('>')?;
    let mut fields = Vec::new();
    let (mut start, mut quoted, mut escaped) = (0, false, false);
    for (at, character) in body.char_indices() {
        if escaped {
            escaped = false;
            continue;
        }
        if character == '\\' && quoted {
            escaped = true;
        } else if character == '"' {
            quoted = !quoted;
        } else if character == ',' && !quoted {
            fields.push(&body[start..at]);
            start = at + 1;
        }
    }
    fields.push(&body[start..]);
    let attribute = |name| {
        fields.iter().find_map(|field| {
            field
                .split_once('=')
                .filter(|(key, _)| *key == name)
                .map(|(_, value)| value.trim_matches('"'))
        })
    };
    let id = attribute("ID")?;
    let kind = match attribute("Type") {
        Some("Integer") => FormatType::Integer,
        Some("Float") => FormatType::Float,
        _ => FormatType::Text,
    };
    Some((id.to_owned(), kind))
}

fn integer_format(text: &str, negative: bool) -> bool {
    let digits = text
        .strip_prefix('+')
        .or_else(|| negative.then(|| text.strip_prefix('-')).flatten())
        .unwrap_or(text);
    !digits.is_empty() && digits.bytes().all(|byte| byte.is_ascii_digit())
}

/// pysam parses sample fields before exposing a record. Invalid GT/numeric
/// FORMAT data must therefore fail even though this producer only uses its ID.
fn check_samples<'a>(
    columns: &mut impl Iterator<Item = &'a [u8]>,
    samples: usize,
    types: &HashMap<String, FormatType>,
) -> Result<()> {
    // QUAL, FILTER and INFO are not used. A short record with no FORMAT is
    // accepted by pysam even when its header declared samples.
    let Some(format) = columns.nth(3) else {
        return Ok(());
    };
    let format = std::str::from_utf8(format).map_err(|_| invalid("VCF FORMAT is not UTF-8"))?;
    let names = format.split(':').collect::<Vec<_>>();
    let mut count = 0;
    for sample in columns {
        count += 1;
        let sample = std::str::from_utf8(sample).map_err(|_| invalid("VCF sample is not UTF-8"))?;
        let values = sample.split(':').collect::<Vec<_>>();
        if count > samples || sample.is_empty() || values.len() > names.len() {
            return Err(invalid(
                "OSError: invalid VCF sample field count (unable to parse next record)",
            ));
        }
        for (name, value) in names.iter().zip(values) {
            let valid = if *name == "GT" {
                value
                    .split(['/', '|'])
                    .all(|allele| allele == "." || integer_format(allele, false))
            } else {
                match types.get(*name).copied().unwrap_or(FormatType::Text) {
                    FormatType::Text => true,
                    FormatType::Integer => value
                        .split(',')
                        .all(|value| value == "." || integer_format(value, true)),
                    FormatType::Float => value.split(',').all(|value| {
                        value == "."
                            || value
                                .trim_start_matches(|c| {
                                    matches!(c, ' ' | '\t' | '\n' | '\r' | '\u{b}' | '\u{c}')
                                })
                                .parse::<f64>()
                                .is_ok()
                    }),
                }
            };
            if !valid {
                return Err(invalid(
                    "OSError: invalid VCF genotype or numeric FORMAT value (unable to parse next record)",
                ));
            }
        }
    }
    if count != samples {
        return Err(invalid(
            "OSError: missing VCF sample fields (unable to parse next record)",
        ));
    }
    Ok(())
}

/// Scan every VCF record, including when the table has no IDs. The original
/// constructs `record.alts[0]` before testing `record.id`; a missing ALT must
/// fail even for an unused ID. Only requested assignments need retained storage.
fn clinvar_map(
    reader: &mut dyn BufRead,
    wanted: &HashSet<String>,
) -> Result<(HashMap<String, String>, usize)> {
    let mut bytes = Vec::new();
    if reader.read_until(b'\n', &mut bytes)? == 0 || !bytes.starts_with(b"##fileformat=VCFv") {
        return Err(invalid("ValueError: input is not a VCF file"));
    }
    let mut header = false;
    let mut samples = 0;
    let mut format_types = HashMap::new();
    loop {
        bytes.clear();
        if reader.read_until(b'\n', &mut bytes)? == 0 {
            break;
        }
        while bytes
            .last()
            .is_some_and(|byte| matches!(byte, b'\r' | b'\n'))
        {
            bytes.pop();
        }
        if bytes.starts_with(b"##") {
            if let Some((name, kind)) = format_definition(&bytes) {
                format_types.entry(name).or_insert(kind);
            }
            continue;
        }
        let fields = bytes.split(|byte| *byte == b'\t').collect::<Vec<_>>();
        header = fields.len() >= 8
            && fields[..8]
                == [
                    b"#CHROM".as_slice(),
                    b"POS",
                    b"ID",
                    b"REF",
                    b"ALT",
                    b"QUAL",
                    b"FILTER",
                    b"INFO",
                ];
        if fields.len() > 8 {
            header &= fields[8] == b"FORMAT";
            samples = fields.len() - 9;
            let names = fields[9..].iter().copied().collect::<HashSet<_>>();
            header &= names.len() == samples && !names.contains(b"".as_slice());
        }
        break;
    }
    if !header {
        return Err(invalid(
            "ValueError: VCF does not have a valid column header",
        ));
    }
    let mut result = HashMap::new();
    let mut count = 0;
    loop {
        bytes.clear();
        if reader.read_until(b'\n', &mut bytes)? == 0 {
            break;
        }
        while bytes
            .last()
            .is_some_and(|byte| matches!(byte, b'\r' | b'\n'))
        {
            bytes.pop();
        }
        // Avoid decoding the large INFO annotation payload. When the header
        // declares samples, preserve pysam's FORMAT parse failures as well.
        let mut columns = bytes.split(|byte| *byte == b'\t');
        let chrom = columns.next().unwrap_or_default();
        let position = columns.next();
        let id = columns.next();
        let reference = columns.next();
        let alternate = columns.next();
        if samples > 0 {
            check_samples(&mut columns, samples, &format_types)?;
        }
        let alternate = alternate.filter(|allele| *allele != b".").ok_or_else(|| {
            invalid("TypeError: VCF record has no alternative allele (record.alts is None)")
        })?;
        let text = |value| std::str::from_utf8(value).map_err(|_| invalid("VCF key is not UTF-8"));
        let position = text(position.unwrap_or_default())?;
        let position: i64 = (if position.is_empty() { "0" } else { position })
            .parse()
            .ok()
            .filter(|value| *value >= 0)
            .ok_or_else(|| invalid("OSError: could not parse the VCF position (truncated file)"))?;
        let reference = reference.unwrap_or_default();
        let reference = if reference.is_empty() {
            "."
        } else {
            text(reference)?
        };
        let first_alt = alternate.split(|byte| *byte == b',').next().unwrap();
        let first_alt = if first_alt.is_empty() {
            "."
        } else {
            text(first_alt)?
        };
        let variant = format!("{}:{position}:{reference}-{first_alt}", text(chrom)?);
        if let Some(id) = id {
            for id in text(id)?.split(';').map(python_strip) {
                if !id.is_empty() && id != "." && wanted.contains(id) {
                    result.insert(id.to_owned(), variant.clone());
                }
            }
        }
        count += 1;
    }
    Ok((result, count))
}

fn curated_frame(ids: &ProducerColumn, codes: &ProducerColumn) -> Result<DataFrame> {
    let mut frame = DataFrame::new(vec![
        ids.series("id")?.into(),
        codes.series("code")?.into(),
        Column::new(
            "source_row".into(),
            (0..ids.len() as u64).collect::<Vec<_>>(),
        ),
    ])?;
    frame = frame.filter(&frame.column("id")?.is_not_null())?;
    let id_text = frame
        .column("source_row")?
        .u64()?
        .into_no_null_iter()
        .map(|row| ids.integer_text(row as usize))
        .collect::<Result<Vec<_>>>()?;
    frame.with_column(Column::new("id".into(), id_text.clone()))?;
    Ok(frame)
}

fn mapped_entries(
    mut frame: DataFrame,
    codes: &ProducerColumn,
    reader: &mut dyn BufRead,
) -> Result<(Vec<(String, Option<String>)>, usize)> {
    let id_text = frame
        .column("id")?
        .str()?
        .into_no_null_iter()
        .map(str::to_owned)
        .collect::<Vec<_>>();
    let wanted = id_text.iter().cloned().collect();
    let (mapping, records) = clinvar_map(reader, &wanted)?;
    let mapped = id_text
        .iter()
        .map(|id| mapping.get(id).map(String::as_str))
        .collect::<Vec<_>>();
    frame.with_column(Column::new("variant".into(), mapped))?;
    frame = frame.filter(&frame.column("variant")?.is_not_null())?;
    // Python's to_dict overwrites a duplicated variant with its last row.
    // Only the schema's explicitly sorted top-level key order is canonicalized.
    let last = frame.unique_stable(
        Some(&["variant".to_owned()]),
        UniqueKeepStrategy::Last,
        None,
    )?;
    let entries = last
        .column("variant")?
        .str()?
        .into_no_null_iter()
        .zip(last.column("source_row")?.u64()?.into_no_null_iter())
        .map(|(variant, row)| Ok((variant.to_owned(), codes.code(row as usize)?)))
        .collect::<Result<Vec<_>>>()?;
    Ok((entries, records))
}

#[cfg(test)]
fn reduce(
    ids: &ProducerColumn,
    codes: &ProducerColumn,
    reader: &mut dyn BufRead,
) -> Result<(Vec<(String, Option<String>)>, usize)> {
    mapped_entries(curated_frame(ids, codes)?, codes, reader)
}

fn input_provenance(table: &Path, vcf: &Path) -> Result<Vec<Provenance>> {
    let mut result = Vec::new();
    for (prefix, path) in [("source", table), ("clinvar_source", vcf)] {
        let basename = path
            .file_name()
            .and_then(|name| name.to_str())
            .ok_or_else(|| invalid("a ClinGen input needs a UTF-8 file name"))?;
        result.extend([
            Provenance {
                key: prefix.to_owned(),
                value: basename.to_owned(),
            },
            Provenance {
                key: format!("{prefix}_bytes"),
                value: std::fs::metadata(path)?.len().to_string(),
            },
            Provenance {
                key: format!("{prefix}_sha256"),
                value: crate::sha256::sha256_file(path)?,
            },
        ]);
    }
    result.push(Provenance {
        key: "converter".to_owned(),
        value: format!("priva clingen-map build {}", env!("CARGO_PKG_VERSION")),
    });
    Ok(result)
}

/// `pysam.VariantFile` accepts BGZF, but refuses ordinary gzip VCFs because
/// its record reader needs BGZF seeking. Table gzip is independently supported.
fn check_vcf_compression(path: &Path) -> Result<()> {
    let mut file = std::fs::File::open(path)?;
    let mut prefix = [0_u8; 12];
    let count = file.read(&mut prefix)?;
    if count >= 2 && prefix[..2] == [0x1f, 0x8b] {
        let fail = || {
            invalid(
                "NotImplementedError: seek not implemented for ordinary gzip VCF; BGZF is required",
            )
        };
        if count < prefix.len() || prefix[3] & 4 == 0 {
            return Err(fail());
        }
        let mut extra = vec![0; u16::from_le_bytes([prefix[10], prefix[11]]) as usize];
        file.read_exact(&mut extra)?;
        let mut at = 0;
        while at + 4 <= extra.len() {
            let size = u16::from_le_bytes([extra[at + 2], extra[at + 3]]) as usize;
            if at + 4 + size > extra.len() {
                break;
            }
            if extra[at..at + 2] == *b"BC" && size == 2 {
                return Ok(());
            }
            at += 4 + size;
        }
        return Err(fail());
    }
    Ok(())
}

pub(crate) fn build(clinvar: &Path, table: &Path, output: &Path) -> Result<String> {
    let start = Instant::now();
    let (ids, codes) = parse_curated(&pandas_table::read_table_text(table)?)?;
    let frame = curated_frame(&ids, &codes)?;
    ensure_distinct_paths(clinvar, output)?;
    ensure_distinct_paths(table, output)?;
    let provenance = input_provenance(table, clinvar)?;
    check_vcf_compression(clinvar)?;
    let mut reader = crate::vcf_text::open_reader(clinvar)?;
    let (entries, records) = mapped_entries(frame, &codes, reader.as_mut())?;
    if provenance != input_provenance(table, clinvar)? {
        return Err(invalid(
            "a ClinGen input changed during production; no archive was published",
        ));
    }
    let without_codes = entries.iter().filter(|(_, value)| value.is_none()).count();
    write_entries(output, &entries, provenance)?;
    Ok(format!(
        "{{\"cache\":\"clingen_map\",\"table_rows\":{},\"vcf_records\":{records},\"variants\":{},\"without_codes\":{without_codes},\"bytes\":{},\"elapsed_s\":{:.3}}}",
        ids.len(),
        entries.len(),
        std::fs::metadata(output)?.len(),
        start.elapsed().as_secs_f64()
    ))
}

fn arguments(args: &[String]) -> Result<Option<(PathBuf, PathBuf, PathBuf)>> {
    if args.is_empty() || args == ["--help"] || args == ["-h"] {
        return Ok(None);
    }
    if args[0] != "build" {
        return Err(invalid(format!("unknown ClinGen procedure {}", args[0])));
    }
    match parse_long_options(&args[1..], &["clinvar-vcf", "clingen-tab", "output"])
        .map_err(invalid)?
    {
        Parsed::Help => Ok(None),
        Parsed::Values(values) => {
            let required = |name: &str| {
                values
                    .get(name)
                    .map(PathBuf::from)
                    .ok_or_else(|| invalid(format!("the following argument is required: --{name}")))
            };
            Ok(Some((
                required("clinvar-vcf")?,
                required("clingen-tab")?,
                required("output")?,
            )))
        }
    }
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let parsed = arguments(args)?;
    let Some((vcf, table, output)) = parsed else {
        println!("{USAGE}");
        return Ok(());
    };
    match build(&vcf, &table, &output) {
        Ok(summary) => {
            println!("{summary}");
            Ok(())
        }
        Err(error) => {
            eprintln!("priva clingen-map: {error}");
            std::process::exit(1)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    fn run(table: &str, records: &str) -> Result<Vec<(String, Option<String>)>> {
        let (ids, codes) = parse_curated(&format!("{ID}\t{CODE}\n{table}"))?;
        let mut input = Cursor::new(format!("{HEADER}{records}"));
        let mut entries = reduce(&ids, &codes, &mut input)?.0;
        entries.sort_by(|a, b| a.0.cmp(&b.0));
        Ok(entries)
    }

    #[test]
    fn duplicate_ids_and_variants_keep_their_last_source_values() {
        let result = run(
            "1\tPS3\n2\tPM2\n1\tBA1\n3\tNA\n4\tPVS1\n",
            "1\t10\t1; 2 ;.;;3\tA\tT,C\t.\t.\t.\n1\t11\t1\tG\tC\t.\t.\t.\n",
        )
        .unwrap();
        assert_eq!(
            result,
            vec![
                ("1:10:A-T".into(), None),
                ("1:11:G-C".into(), Some("BA1".into()))
            ]
        );
    }

    #[test]
    fn pandas_ids_keep_missingness_integer_casts_and_object_rules() {
        let records =
            "1\t1\t1\tA\tT\n1\t2\t0\tA\tC\n1\t3\t-9223372036854775808\tA\tG\n1\t4\t1000\tG\tT\n";
        assert_eq!(
            run("1.9\tPS3\nNA\tPM2\n", records).unwrap(),
            vec![("1:1:A-T".into(), Some("PS3".into()))]
        );
        assert_eq!(
            run("True\tPS3\nFalse\tPM2\n\tPVS1\n", records)
                .unwrap()
                .len(),
            2
        );
        assert_eq!(
            run("1e30\tPS3\n", records).unwrap(),
            vec![("1:3:A-G".into(), Some("PS3".into()))]
        );
        assert_eq!(
            run("9223372036854775808\tPS3\n", records).unwrap(),
            vec![("1:3:A-G".into(), Some("PS3".into()))]
        );
        assert!(run("9223372036854775808\tPS3\nNA\tPM2\n", records).is_err());
        assert_eq!(
            run("1_000\tPS3\n", records).unwrap(),
            vec![("1:4:G-T".into(), Some("PS3".into()))]
        );
        assert!(run("1.0\tPS3\nx\tPM2\n", records).is_err());
        assert!(run("inf\tPS3\n", records).is_err());
        assert_eq!(object_integer("\u{a0}+１_２\u{2003}").unwrap(), 12);
        assert!(object_integer("1.2").is_err());
    }

    #[test]
    fn typed_codes_are_checked_only_after_mapping_and_last_row_reduction() {
        assert!(run("1\t2\n", "1\t1\t1\tA\tT\n").is_err());
        assert!(run("1\tTrue\n", "1\t1\t1\tA\tT\n").is_err());
        assert!(run("1\t2.0\n", "1\t1\t1\tA\tT\n").is_err());
        assert!(run("1\t2\n", "1\t1\t9\tA\tT\n").unwrap().is_empty());
        assert_eq!(
            run("1\tNA\n", "1\t1\t1\tA\tT\n").unwrap(),
            vec![("1:1:A-T".into(), None)]
        );
    }

    #[test]
    fn empty_tables_and_unused_ids_still_parse_every_vcf_record() {
        for rows in ["", "NA\tPM2\n", "1\tPS3\n"] {
            assert!(run(rows, "1\t1\t.\tA\t.\n").is_err());
            assert!(run(rows, "1\t1\t9\tA\tT\n\n").is_err());
            assert!(run(rows, "1\tx\t9\tA\tT\n").is_err());
        }
    }

    #[test]
    fn pysam_permissive_fields_and_first_alternate_are_preserved() {
        assert_eq!(
            run(
                "1\tPS3\n2\tPM2\n3\tPVS1\n",
                "1\t0\t1\t\t\n1\t1\t2\tA\t.,T\n1\t2\t3\tA\t<DEL>,C\tbad\n"
            )
            .unwrap(),
            vec![
                ("1:0:.-.".into(), Some("PS3".into())),
                ("1:1:A-.".into(), Some("PM2".into())),
                ("1:2:A-<DEL>".into(), Some("PVS1".into()))
            ]
        );
        assert!(run("1\tPS3\n", "1\t-1\t1\tA\tT\n").is_err());
        assert!(run("1\tPS3\n", "1\t1\t1\tA\n").is_err());
        assert_eq!(
            run("1\tPS3\n", "1\t\t1\tA\tT\n").unwrap(),
            vec![("1:0:A-T".into(), Some("PS3".into()))]
        );
    }

    #[test]
    fn malformed_sample_fields_fail_even_for_an_unused_clinvar_id() {
        let wanted = HashSet::from(["1".to_owned()]);
        for (kind, valid, invalid) in [
            ("GT", "0/1", "0//1"),
            ("DP", "+2", "2x"),
            ("AF", "1e999", "2x"),
        ] {
            let value_type = match kind {
                "GT" => "String",
                "DP" => "Integer",
                _ => "Float",
            };
            let header = HEADER.replace(
                "#CHROM",
                &format!(
                    "##FORMAT=<ID={kind},Number=1,Type={value_type},Description=\"Value\">\n#CHROM"
                ),
            );
            let header = format!("{}\tFORMAT\tS1\n", header.trim_end());
            let body = |value| format!("{header}1\t1\t9\tA\tT\t.\t.\t.\t{kind}\t{value}\n");
            assert!(
                clinvar_map(&mut Cursor::new(body(valid)), &wanted)
                    .unwrap()
                    .0
                    .is_empty()
            );
            assert!(clinvar_map(&mut Cursor::new(body(invalid)), &wanted).is_err());
            assert!(clinvar_map(&mut Cursor::new(body("")), &wanted).is_err());
        }
    }

    #[test]
    fn pandas_quoting_header_mangling_and_implicit_index_are_retained() {
        let table = format!("\u{feff}{ID}\t{CODE}\t{ID}\r\nindex\t1\t\"PS3\tPM2\"\tbad\r\n");
        let (ids, codes) = parse_curated(&table).unwrap();
        assert_eq!(ids.integer_text(0).unwrap(), "1");
        assert_eq!(codes.code(0).unwrap(), Some("PS3\tPM2".into()));
        assert!(parse_curated(&format!("{ID}\t{CODE}\n1\tPS3\n2\tPM2\textra\n")).is_err());
        assert!(parse_curated(&format!("{ID}\n1\n")).is_err());
    }

    #[test]
    fn input_failures_preserve_an_existing_archive_and_aliases_are_refused() {
        let directory =
            std::env::temp_dir().join(format!("priva_clingen_builder_{}", std::process::id()));
        std::fs::create_dir_all(&directory).unwrap();
        let table = directory.join("curated.tsv");
        let vcf = directory.join("clinvar.vcf");
        let output = directory.join("map.rkyv");
        std::fs::write(&table, format!("{ID}\t{CODE}\n1\tPS3\n")).unwrap();
        std::fs::write(&vcf, format!("{HEADER}1\t1\t1\tA\t.\n")).unwrap();
        std::fs::write(&output, b"previous archive").unwrap();
        assert!(build(&vcf, &table, &output).is_err());
        assert_eq!(std::fs::read(&output).unwrap(), b"previous archive");
        assert!(build(&vcf, &table, &table).is_err());
        assert!(build(&vcf, &table, &vcf).is_err());
        // The Python producer casts the table IDs before opening its VCF.
        std::fs::write(&table, format!("{ID}\t{CODE}\nx\tPS3\n")).unwrap();
        assert!(
            build(&directory.join("absent.vcf"), &table, &output)
                .unwrap_err()
                .to_string()
                .starts_with("ValueError")
        );
        std::fs::write(&table, format!("{ID}\t{CODE}\n1\tPS3\n")).unwrap();
        std::fs::write(&vcf, format!("{HEADER}1\t1\t1\tA\tT\n")).unwrap();
        build(&vcf, &table, &output).unwrap();
        assert_eq!(
            crate::cache_convert::clingen_map::ClingenCache::open(output.to_str().unwrap())
                .unwrap()
                .get("1:1:A-T"),
            Some(Some("PS3"))
        );
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn ordinary_gzip_vcf_is_refused_but_bgzf_and_gzip_tables_work() {
        use std::io::Write;
        let directory =
            std::env::temp_dir().join(format!("priva_clingen_compression_{}", std::process::id()));
        std::fs::create_dir_all(&directory).unwrap();
        let table = directory.join("curated.tsv.gz");
        let vcf = directory.join("clinvar.vcf.gz");
        let output = directory.join("map.rkyv");
        let mut gzip = flate2::write::GzEncoder::new(
            std::fs::File::create(&table).unwrap(),
            flate2::Compression::default(),
        );
        gzip.write_all(format!("{ID}\t{CODE}\n1\tPS3\n").as_bytes())
            .unwrap();
        gzip.finish().unwrap();
        let contents = format!("{HEADER}1\t1\t1\tA\tT\n");
        let mut gzip = flate2::write::GzEncoder::new(
            std::fs::File::create(&vcf).unwrap(),
            flate2::Compression::default(),
        );
        gzip.write_all(contents.as_bytes()).unwrap();
        gzip.finish().unwrap();
        assert!(
            build(&vcf, &table, &output)
                .unwrap_err()
                .to_string()
                .starts_with("NotImplementedError")
        );
        assert!(!output.exists());
        let mut bgzf = crate::vcf_text::create_writer(&vcf).unwrap();
        bgzf.write_line(contents.as_bytes()).unwrap();
        bgzf.finish().unwrap();
        build(&vcf, &table, &output).unwrap();
        assert_eq!(
            crate::cache_convert::clingen_map::ClingenCache::open(output.to_str().unwrap())
                .unwrap()
                .get("1:1:A-T"),
            Some(Some("PS3"))
        );
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn procedure_arguments_use_named_options_and_validate_required_inputs() {
        let args = [
            "build",
            "--clinvar-vcf=v.vcf",
            "--clingen-tab",
            "c.tsv",
            "--output",
            "c.rkyv",
        ]
        .map(str::to_owned);
        assert_eq!(
            arguments(&args).unwrap(),
            Some(("v.vcf".into(), "c.tsv".into(), "c.rkyv".into()))
        );
        assert!(arguments(&["build".into()]).is_err());
        assert!(arguments(&["unknown".into()]).is_err());
        assert_eq!(arguments(&["build".into(), "--help".into()]).unwrap(), None);
    }
}
