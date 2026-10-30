//! Rust version of `scripts/combine_annotations.py`, which turns an annotated
//! family VCF into PriVA's transcript-level table: one row per record and
//! Ensembl feature (CSQ Feature starting with `ENS`) with the record's gnomAD,
//! ClinVar, FILTER, VARIANT_SOURCE, control and PEXT values, the samples'
//! genotypes and allele depths, the feature's CADD PHRED score and the gene's
//! HPO annotation.
//!
//! The table is written as the Python script writes it with `-t 1`; with more
//! threads its worker processes return records in a varying order, while this
//! version always keeps the VCF's order. Behaviour reproduced from the script:
//! - values are read as pysam returns them: Float INFO values at single
//!   precision, a String with Number=1 split at commas only when it has one,
//!   `.` as None for numbers, an allele beyond the record's alleles missing;
//!   an INFO field the header does not define stops the run (pysam's "Invalid
//!   header"), except the sex-specific gnomAD fields, until a record uses it:
//!   htslib then defines it as one text value;
//! - CSQ fields of digits become integers and `""`, `.`, `-`, `NaN`, `nan`,
//!   `NA`, `N/A` missing values; a GT is joined with `|` only when every allele
//!   after the first is phased; a missing AD is rebuilt as `0,DP` for a 1/1
//!   genotype, else `.,.`;
//! - the script's column typing: rows are gathered in chunks of at least 20,000
//!   (whole records); in a chunk a column with any number becomes single
//!   precision unless it already holds text, and turns into a Python object
//!   column as soon as text arrives, its stored numbers keeping their single
//!   precision value; at the end a column without text becomes pandas Float32,
//!   or Int32 when every value is whole (out-of-range values wrap to
//!   -2147483648, as numpy's cast does);
//! - numbers are written as pandas writes them: Float32 in numpy's shortest
//!   single-precision form (`1e-05`, `0.1`), numbers in text columns in
//!   Python's `repr` (`123.0`, `0.10000000149011612` for a value stored at
//!   single precision); text holding a tab, quote or line break is quoted;
//! - the CADD and HPO tables are read as pandas reads them, plain or gzip;
//!   CADD PHRED is kept at single precision, chromosomes get a `chr` prefix,
//!   the first row of each variant and feature counts, and a regulatory
//!   feature's score goes to CADD_reg_phred;
//! - HPO annotations come from the disease-specific assertion table (or the
//!   older one-row-per-gene table): values stripped, pandas' missing-value
//!   spellings read as `-`, assertions written as Python's
//!   `json.dumps(..., ensure_ascii=True, separators=(",", ":"), sort_keys=True)`.
use crate::cli::{Result, invalid};
use crate::hpo_penetrance::{PENETRANCE_STATUS_BY_TERM, py_strip};
use crate::vcf_text::open_reader;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufWriter, Read, Write};
use std::path::Path;

const USAGE: &str = "usage: priva combine-annotations -i INPUT_VCF -c CADD_TABLE -p HPO_TABLE [-o OUTPUT_TABLE] [-t THREADS]";
const CHUNK_SIZE: usize = 20_000;
/// pandas' default missing-value spellings (`keep_default_na`).
const PANDAS_NA: [&str; 19] = [
    "", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN",
    "<NA>", "N/A", "NA", "NULL", "NaN", "None", "n/a", "nan", "null",
];
const HPO_COLUMNS: [&str; 8] = [
    "HPO_IDs",
    "HPO_terms",
    "HPO_sources",
    "HPO_gene_inheritance",
    "HPO_assertions",
    "HPO_scope_review_required",
    "HPO_scope_review_disease_ids",
    "HPO_scope_excluded_disease_ids",
];
const ASSERTION_REQUIRED: [&str; 6] = [
    "gene_symbol",
    "disease_id",
    "hpo_id",
    "frequency",
    "evidence",
    "reference",
];
const ASSERTION_SCOPE: [&str; 7] = [
    "mondo_id",
    "mondo_name",
    "disease_scope",
    "priva_scope",
    "scope_evidence",
    "scope_reference",
    "scope_review_status",
];
const INHERITANCE_TERMS: [(&str, &str); 13] = [
    ("HP:0000006", "Autosomal dominant inheritance"),
    ("HP:0000007", "Autosomal recessive inheritance"),
    ("HP:0001417", "X-linked inheritance"),
    ("HP:0001419", "X-linked recessive inheritance"),
    ("HP:0001423", "X-linked dominant inheritance"),
    ("HP:0001426", "Non-Mendelian inheritance"),
    ("HP:0001427", "Mitochondrial inheritance"),
    ("HP:0001450", "Y-linked inheritance"),
    ("HP:0010982", "Polygenic inheritance"),
    ("HP:0010983", "Oligogenic inheritance"),
    ("HP:0010984", "Digenic inheritance"),
    (
        "HP:0012275",
        "Autosomal dominant inheritance with maternal imprinting",
    ),
    ("HP:0034341", "Pseudoautosomal recessive inheritance"),
];
const ONSET_TERMS: [&str; 15] = [
    "HP:0030674",
    "HP:0011460",
    "HP:0011461",
    "HP:0003577",
    "HP:0003623",
    "HP:0003593",
    "HP:0011463",
    "HP:0003621",
    "HP:0410280",
    "HP:0003581",
    "HP:0011462",
    "HP:0003596",
    "HP:0003584",
    "HP:0034857",
    "HP:0003587",
];
/// Mild and Asymptomatic, kept in the flat annotation only.
const SEVERITY_TERMS: [&str; 2] = ["HP:0012829", "HP:0040007"];
const POPULATIONS: [&str; 9] = [
    "nfe",
    "eas",
    "afr",
    "amr",
    "asj",
    "fin",
    "sas",
    "mid",
    "remaining",
];

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args == ["--help"] || args == ["-h"] {
        println!("{USAGE}");
        return Ok(());
    }
    let mut options: HashMap<&str, String> = HashMap::new();
    let mut index = 0;
    while index < args.len() {
        let (option, inline) = match args[index].split_once('=') {
            Some((name, value)) if name.starts_with("--") => (name, Some(value.to_owned())),
            _ => (args[index].as_str(), None),
        };
        let name = match option {
            "-i" | "--input" => "input",
            "-o" | "--output" => "output",
            "-t" | "--threads" => "threads",
            "-c" | "--cadd" => "cadd",
            "-p" | "--hpo" => "hpo",
            _ => return Err(invalid(format!("unknown option {option:?}\n{USAGE}"))),
        };
        let value = match inline {
            Some(value) => value,
            None => {
                index += 1;
                args.get(index)
                    .cloned()
                    .ok_or_else(|| invalid(format!("{option} needs a value\n{USAGE}")))?
            }
        };
        options.insert(name, value);
        index += 1;
    }
    let mut take = |name: &str| {
        options
            .remove(name)
            .ok_or_else(|| invalid(format!("--{name} is required\n{USAGE}")))
    };
    let (input, cadd, hpo) = (take("input")?, take("cadd")?, take("hpo")?);
    // Records are always handled in VCF order; the thread count is accepted
    // for the caller's sake.
    if let Some(threads) = options.get("threads") {
        threads
            .parse::<i64>()
            .map_err(|_| invalid(format!("--threads must be an integer\n{USAGE}")))?;
    }
    let output = options.remove("output");
    // Like the Python script, a failure after the arguments exits with code 1.
    if let Err(error) = combine(
        Path::new(&input),
        Path::new(&cadd),
        Path::new(&hpo),
        output.as_deref().map(Path::new),
    ) {
        eprintln!("priva combine-annotations: {error}");
        std::process::exit(1);
    }
    Ok(())
}

fn combine(input: &Path, cadd: &Path, hpo: &Path, output: Option<&Path>) -> Result<()> {
    let cadd = read_cadd(cadd)?;
    eprintln!("INFO: CADD scores loaded for {} variants", cadd.len());
    let hpo = read_hpo(hpo)?;
    eprintln!("INFO: HPO annotations loaded for {} genes", hpo.len());
    let table = convert_vcf(input, &cadd, &hpo)?;
    eprintln!(
        "INFO: The annotation table has {} rows and {} columns",
        table.count,
        table.names.len()
    );
    if let Some(output) = output {
        table.write(output)?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Python values

/// A Python value as the script handles it.
#[derive(Clone, Debug, PartialEq)]
enum Py {
    None,
    Bool(bool),
    Int(i64),
    /// A CSQ field of ASCII digits, which the script turns into a Python int.
    Digits(String),
    Float(f64),
    Str(String),
    Tuple(Vec<Py>),
}

const NAN: Py = Py::Float(f64::NAN);

impl Py {
    /// The script's `na_value`.
    fn is_na(&self) -> bool {
        match self {
            Py::Float(x) => x.is_nan(),
            Py::None => true,
            Py::Str(s) => [
                "NaN", "nan", "na", "NA", "NAN", "None", "none", "", ".", "-",
            ]
            .contains(&s.as_str()),
            _ => false,
        }
    }

    /// Python `float(value)`, or `None` where it raises.
    fn to_float(&self) -> Option<f64> {
        match self {
            Py::Bool(b) => Some(if *b { 1.0 } else { 0.0 }),
            Py::Int(i) => Some(*i as f64),
            Py::Digits(d) => d.parse().ok(),
            Py::Float(x) => Some(*x),
            Py::Str(s) => py_float(s),
            Py::None | Py::Tuple(_) => None,
        }
    }
}

/// Python `value[0]`.
fn first(value: Py, key: &str) -> Result<Py> {
    match value {
        Py::Tuple(items) => items
            .into_iter()
            .next()
            .ok_or_else(|| invalid(format!("{key}: index out of range"))),
        Py::Str(s) => s
            .chars()
            .next()
            .map(|c| Py::Str(c.to_string()))
            .ok_or_else(|| invalid(format!("{key}: index out of range"))),
        _ => Err(invalid(format!("{key}: the value is not subscriptable"))),
    }
}

/// Python `",".join(value)`; a text is joined character by character.
fn comma_join(value: Py, key: &str) -> Result<String> {
    match value {
        Py::Tuple(items) => items
            .into_iter()
            .map(|item| match item {
                Py::Str(s) => Ok(s),
                _ => Err(invalid(format!("{key}: join expects text"))),
            })
            .collect::<Result<Vec<_>>>()
            .map(|items| items.join(",")),
        Py::Str(s) => Ok(s.chars().map(String::from).collect::<Vec<_>>().join(",")),
        _ => Err(invalid(format!("{key}: join expects text"))),
    }
}

/// Python `max` over values that are all numbers or all text.
fn py_max(items: Vec<Py>, key: &str) -> Result<Py> {
    let numeric = |item: &Py| {
        matches!(
            item,
            Py::Bool(_) | Py::Int(_) | Py::Digits(_) | Py::Float(_)
        )
    };
    let mut best: Option<Py> = None;
    for item in items {
        let replace = match &best {
            None => true,
            Some(current) if numeric(current) && numeric(&item) => {
                item.to_float().unwrap_or(f64::NAN) > current.to_float().unwrap_or(f64::NAN)
            }
            Some(Py::Str(current)) => match &item {
                Py::Str(text) => text > current,
                _ => return Err(invalid(format!("{key}: max over mixed values"))),
            },
            Some(_) => return Err(invalid(format!("{key}: max over mixed values"))),
        };
        if replace {
            best = Some(item);
        }
    }
    Ok(best.unwrap_or(NAN))
}

/// Python `float(text)`: surrounding whitespace, a sign, `_` between digits,
/// `inf`, `infinity` and `nan` in any case.
fn py_float(text: &str) -> Option<f64> {
    let text = py_strip(text);
    let (negative, body) = match text.as_bytes().first() {
        Some(b'-') => (true, &text[1..]),
        Some(b'+') => (false, &text[1..]),
        _ => (false, text),
    };
    let sign = if negative { -1.0 } else { 1.0 };
    match body.to_ascii_lowercase().as_str() {
        "inf" | "infinity" => return Some(sign * f64::INFINITY),
        "nan" => return Some(f64::NAN),
        _ => {}
    }
    let bytes = body.as_bytes();
    if bytes.is_empty() || !body.is_ascii() {
        return None;
    }
    for (index, byte) in bytes.iter().enumerate() {
        let allowed = match byte {
            b'0'..=b'9' | b'.' | b'e' | b'E' => true,
            b'+' | b'-' => index > 0 && matches!(bytes[index - 1], b'e' | b'E'),
            b'_' => {
                index > 0
                    && index + 1 < bytes.len()
                    && bytes[index - 1].is_ascii_digit()
                    && bytes[index + 1].is_ascii_digit()
            }
            _ => false,
        };
        if !allowed {
            return None;
        }
    }
    body.replace('_', "")
        .parse::<f64>()
        .ok()
        .map(|value| sign * value)
}

// ---------------------------------------------------------------------------
// Number formatting

/// Sign, shortest digits and decimal exponent from Rust's `{:e}` output.
fn shortest_digits(scientific: &str) -> (bool, String, i32) {
    let (negative, text) = match scientific.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, scientific),
    };
    let (mantissa, exponent) = text.split_once('e').unwrap_or((text, "0"));
    (
        negative,
        mantissa.replace('.', ""),
        exponent.parse().unwrap_or(0),
    )
}

fn place_digits(negative: bool, digits: &str, exponent: i32, scientific: bool) -> String {
    let mut out = String::new();
    if negative {
        out.push('-');
    }
    if scientific {
        out.push_str(&digits[..1]);
        if digits.len() > 1 {
            out.push('.');
            out.push_str(&digits[1..]);
        }
        out.push('e');
        out.push(if exponent < 0 { '-' } else { '+' });
        out.push_str(&format!("{:02}", exponent.abs()));
    } else if exponent >= 0 {
        let whole = exponent as usize + 1;
        if digits.len() <= whole {
            out.push_str(digits);
            out.push_str(&"0".repeat(whole - digits.len()));
            out.push_str(".0");
        } else {
            out.push_str(&digits[..whole]);
            out.push('.');
            out.push_str(&digits[whole..]);
        }
    } else {
        out.push_str("0.");
        out.push_str(&"0".repeat((-exponent - 1) as usize));
        out.push_str(digits);
    }
    out
}

fn special(value: f64) -> Option<String> {
    if value.is_nan() {
        Some("nan".to_owned())
    } else if value.is_infinite() {
        Some(if value > 0.0 { "inf" } else { "-inf" }.to_owned())
    } else {
        None
    }
}

/// Shortest digits with an exact tie broken towards the even digit, as
/// Python's `repr` and numpy do; Rust's shortest form can round such a tie
/// up. `exact` is the value (exactly representable as f64) and `round_trips`
/// checks a candidate spelling.
fn tie_to_even(
    exact: f64,
    digits: String,
    exponent: i32,
    round_trips: impl Fn(&str) -> bool,
) -> (String, i32) {
    let count = digits.len();
    let tied = |text: &str| {
        let (_, all, _) = shortest_digits(text);
        all.len() > count
            && all.as_bytes()[count] == b'5'
            && all[count + 1..].bytes().all(|b| b == b'0')
    };
    // A cheap probe first; only a possible tie pays for the full expansion.
    if !tied(&format!("{:.*e}", count + 9, exact.abs())) || !tied(&format!("{:.780e}", exact.abs()))
    {
        return (digits, exponent);
    }
    let (_, all, exact_exponent) = shortest_digits(&format!("{:.780e}", exact.abs()));
    let mut kept: Vec<u8> = all.as_bytes()[..count].to_vec();
    let mut new_exponent = exact_exponent;
    if (kept[count - 1] - b'0') % 2 == 1 {
        let mut index = count;
        loop {
            if index == 0 {
                kept.insert(0, b'1');
                kept.pop();
                new_exponent += 1;
                break;
            }
            index -= 1;
            if kept[index] == b'9' {
                kept[index] = b'0';
            } else {
                kept[index] += 1;
                break;
            }
        }
    }
    while kept.len() > 1 && kept.last() == Some(&b'0') {
        kept.pop();
    }
    let candidate = String::from_utf8(kept).unwrap_or_default();
    let sign = if exact < 0.0 { "-" } else { "" };
    let spelled = format!(
        "{sign}{}.{}e{new_exponent}",
        &candidate[..1],
        if candidate.len() > 1 {
            &candidate[1..]
        } else {
            "0"
        }
    );
    if round_trips(&spelled) {
        (candidate, new_exponent)
    } else {
        (digits, exponent)
    }
}

/// Python's `repr` of a float.
fn python_repr(value: f64) -> String {
    if let Some(text) = special(value) {
        return text;
    }
    let (negative, digits, exponent) = shortest_digits(&format!("{value:e}"));
    let (digits, exponent) = tie_to_even(value, digits, exponent, |text| {
        text.parse::<f64>() == Ok(value)
    });
    place_digits(negative, &digits, exponent, !(-4..16).contains(&exponent))
}

/// numpy's `repr` of a float32, as pandas writes Float32 values.
fn numpy_f32_repr(value: f32) -> String {
    if let Some(text) = special(value as f64) {
        return text;
    }
    let magnitude = (value as f64).abs();
    let positional = magnitude == 0.0 || (1e-4..1e16).contains(&magnitude);
    let (negative, digits, exponent) = shortest_digits(&format!("{value:e}"));
    let (digits, exponent) = tie_to_even(value as f64, digits, exponent, |text| {
        text.parse::<f32>() == Ok(value)
    });
    place_digits(negative, &digits, exponent, !positional)
}

/// numpy's float32 to int32 cast: out-of-range values become the minimum.
fn wrap_int32(value: f32) -> i64 {
    if (-2147483648.0..2147483648.0).contains(&value) {
        value as i64
    } else {
        i32::MIN as i64
    }
}

/// A field as pandas' `to_csv` writes it with a tab separator.
fn csv_field(text: &str) -> String {
    if text.contains(['\t', '"', '\n', '\r']) {
        format!("\"{}\"", text.replace('"', "\"\""))
    } else {
        text.to_owned()
    }
}

// ---------------------------------------------------------------------------
// pandas tables (CADD and HPO)

/// A tab-separated table read as pandas' C parser does: the first line holds
/// the column names, a field opening with `"` is quoted, blank lines are
/// skipped, and missing trailing fields are missing values.
struct Table {
    columns: Vec<String>,
    rows: Vec<Vec<Option<String>>>,
}

impl Table {
    fn read(path: &Path) -> Result<Table> {
        // Plain or gzip text, as pandas opens a `.gz` table.
        let mut text = String::new();
        open_reader(path)?.read_to_string(&mut text)?;
        let mut lines = text
            .split('\n')
            .map(|line| line.strip_suffix('\r').unwrap_or(line))
            .filter(|line| !line.is_empty());
        let columns: Vec<String> = match lines.next() {
            Some(line) => pandas_fields(line)?,
            None => return Err(invalid(format!("{}: no columns to parse", path.display()))),
        };
        let mut rows = Vec::new();
        for line in lines {
            let fields = pandas_fields(line)?;
            if fields.len() > columns.len() {
                return Err(invalid(format!(
                    "{}: expected {} fields, saw {}",
                    path.display(),
                    columns.len(),
                    fields.len()
                )));
            }
            let mut row: Vec<Option<String>> = fields.into_iter().map(Some).collect();
            row.resize(columns.len(), None);
            rows.push(row);
        }
        Ok(Table { columns, rows })
    }

    fn column(&self, name: &str) -> Option<usize> {
        self.columns.iter().position(|column| column == name)
    }

    fn require(&self, names: &[&str]) -> Result<Vec<usize>> {
        names
            .iter()
            .map(|name| {
                self.column(name).ok_or_else(|| {
                    invalid(format!(
                        "Usecols do not match columns, columns expected but not found: {name}"
                    ))
                })
            })
            .collect()
    }
}

fn pandas_fields(line: &str) -> Result<Vec<String>> {
    let mut fields = Vec::new();
    let mut chars = line.chars().peekable();
    loop {
        let mut field = String::new();
        if chars.peek() == Some(&'"') {
            chars.next();
            loop {
                match chars.next() {
                    Some('"') if chars.peek() == Some(&'"') => {
                        chars.next();
                        field.push('"');
                    }
                    Some('"') => break,
                    Some(c) => field.push(c),
                    None => {
                        return Err(invalid(
                            "a quoted field spanning lines is outside the validated contract",
                        ));
                    }
                }
            }
        }
        loop {
            match chars.next() {
                Some('\t') => {
                    fields.push(field);
                    break;
                }
                Some(c) => field.push(c),
                None => {
                    fields.push(field);
                    return Ok(fields);
                }
            }
        }
    }
}

/// pandas' missing-value spellings read as missing.
fn pandas_value(value: &Option<String>) -> Option<&str> {
    value.as_deref().filter(|value| !PANDAS_NA.contains(value))
}

/// variant ID (`chr1:100:A-G`) to the CADD PHRED score of each feature.
type CaddScores = HashMap<String, HashMap<String, f32>>;

fn read_cadd(path: &Path) -> Result<CaddScores> {
    let table = Table::read(path)?;
    let index = table.require(&["#Chrom", "Pos", "Ref", "Alt", "PHRED", "FeatureID"])?;
    let mut scores: CaddScores = HashMap::new();
    for row in &table.rows {
        let field = |column: usize| pandas_value(&row[index[column]]);
        // pandas writes a missing chromosome as "<NA>" once it is text.
        let chrom = field(0).unwrap_or("<NA>");
        let pos: i64 = field(1)
            .and_then(|pos| pos.trim().parse().ok())
            .filter(|pos| i32::try_from(*pos).is_ok())
            .ok_or_else(|| invalid("CADD Pos must be a 32-bit integer"))?;
        let phred = match field(4) {
            None => f32::NAN,
            Some(text) => py_float(text)
                .ok_or_else(|| invalid(format!("CADD PHRED {text:?} is not a number")))?
                as f32,
        };
        // A missing REF or ALT makes the variant ID missing: groupby drops it.
        let (Some(reference), Some(alternate)) = (field(2), field(3)) else {
            continue;
        };
        let chrom = if chrom.starts_with("chr") {
            chrom.to_owned()
        } else {
            format!("chr{chrom}")
        };
        let features = scores
            .entry(format!("{chrom}:{pos}:{reference}-{alternate}"))
            .or_default();
        // A missing FeatureID never matches a CSQ feature.
        if let Some(feature) = field(5) {
            features.entry(feature.to_owned()).or_insert(phred);
        }
    }
    Ok(scores)
}

/// Gene symbol to its eight HPO columns.
type HpoAnnotations = HashMap<String, Vec<Py>>;

fn read_hpo(path: &Path) -> Result<HpoAnnotations> {
    let table = Table::read(path)?;
    if ASSERTION_REQUIRED
        .iter()
        .all(|name| table.column(name).is_some())
    {
        assertion_annotations(&table)
    } else {
        legacy_annotations(&table)
    }
}

fn context_term(id: &str) -> bool {
    INHERITANCE_TERMS.iter().any(|(term, _)| *term == id)
        || ONSET_TERMS.contains(&id)
        || PENETRANCE_STATUS_BY_TERM
            .iter()
            .any(|(term, _)| *term == id)
        || SEVERITY_TERMS.contains(&id)
}

fn ordered_unique<'a>(values: impl Iterator<Item = &'a str>) -> Vec<&'a str> {
    let mut seen = Vec::new();
    for value in values {
        if !value.is_empty() && value != "-" && !seen.contains(&value) {
            seen.push(value);
        }
    }
    seen
}

fn assertion_annotations(table: &Table) -> Result<HpoAnnotations> {
    // Columns in the script's order: required ones, then the scope columns present.
    let names: Vec<&str> = ASSERTION_REQUIRED
        .iter()
        .chain(ASSERTION_SCOPE.iter())
        .copied()
        .filter(|name| table.column(name).is_some())
        .collect();
    let index = table.require(&names)?;
    let at = |name: &str| names.iter().position(|n| *n == name).unwrap();
    let mut seen = HashSet::new();
    let mut groups: Vec<(String, Vec<Vec<String>>)> = Vec::new();
    let mut group_of: HashMap<String, usize> = HashMap::new();
    for row in &table.rows {
        let values: Vec<String> = index
            .iter()
            .map(|&column| {
                let value = py_strip(pandas_value(&row[column]).unwrap_or("-"));
                if value.is_empty() { "-" } else { value }.to_owned()
            })
            .collect();
        let symbol = values[at("gene_symbol")].clone();
        if symbol == "-" || !context_term(&values[at("hpo_id")]) || !seen.insert(values.clone()) {
            continue;
        }
        let group = *group_of.entry(symbol.clone()).or_insert_with(|| {
            groups.push((symbol, Vec::new()));
            groups.len() - 1
        });
        groups[group].1.push(values);
    }
    let scoped = names.contains(&"priva_scope");
    let mut fields: Vec<&str> = vec!["disease_id", "hpo_id", "frequency", "evidence", "reference"];
    fields.extend(
        ASSERTION_SCOPE
            .iter()
            .copied()
            .filter(|name| names.contains(name)),
    );
    let mut sorted_fields = fields.clone();
    sorted_fields.sort();
    let mut annotations = HashMap::new();
    for (symbol, rows) in groups {
        let in_scope = |scope: &str| -> Vec<&Vec<String>> {
            rows.iter()
                .filter(|row| !scoped || row[at("priva_scope")] == scope)
                .collect()
        };
        let included = in_scope("include");
        let (review, excluded): (Vec<&Vec<String>>, Vec<&Vec<String>>) = if scoped {
            (in_scope("review"), in_scope("exclude"))
        } else {
            (Vec::new(), Vec::new())
        };
        let column = |rows: &[&Vec<String>], name: &str| -> Vec<String> {
            ordered_unique(rows.iter().map(|row| row[at(name)].as_str()))
                .into_iter()
                .map(str::to_owned)
                .collect()
        };
        let hpo_ids = column(&included, "hpo_id");
        let inheritance: Vec<&str> = ordered_unique(hpo_ids.iter().filter_map(|id| {
            INHERITANCE_TERMS
                .iter()
                .find(|(term, _)| *term == id.as_str())
                .map(|(_, mode)| *mode)
        }));
        let assertions: Vec<String> = rows
            .iter()
            .map(|row| {
                let items: Vec<String> = sorted_fields
                    .iter()
                    .map(|name| format!("{}:{}", json_string(name), json_string(&row[at(name)])))
                    .collect();
                format!("{{{}}}", items.join(","))
            })
            .collect();
        annotations.insert(
            symbol,
            vec![
                Py::Str(hpo_ids.join(";")),
                Py::Str(String::new()),
                Py::Str(column(&included, "disease_id").join(";")),
                Py::Str(inheritance.join(";")),
                Py::Str(format!("[{}]", assertions.join(","))),
                Py::Int(i64::from(!review.is_empty())),
                Py::Str(column(&review, "disease_id").join(";")),
                Py::Str(column(&excluded, "disease_id").join(";")),
            ],
        );
    }
    Ok(annotations)
}

fn legacy_annotations(table: &Table) -> Result<HpoAnnotations> {
    let index = table.require(&[
        "gene_symbol",
        "hpo_id",
        "hpo_name",
        "disease_id",
        "inheritance_modes",
    ])?;
    let mut annotations = HashMap::new();
    for row in &table.rows {
        let Some(symbol) = pandas_value(&row[index[0]]) else {
            continue;
        };
        if symbol == "-" {
            continue;
        }
        let values = [index[1], index[2], index[3], index[4]]
            .iter()
            .map(|&column| {
                pandas_value(&row[column])
                    .map(|value| Py::Str(value.to_owned()))
                    .ok_or_else(|| {
                        invalid("a missing value in the one-row-per-gene HPO table is outside the validated contract")
                    })
            })
            .collect::<Result<Vec<_>>>()?;
        let mut values = values;
        values.extend([
            Py::Str("[]".to_owned()),
            Py::Int(0),
            Py::Str(String::new()),
            Py::Str(String::new()),
        ]);
        if annotations.insert(symbol.to_owned(), values).is_some() {
            return Err(invalid(format!(
                "HPO table contains duplicate gene symbols: {symbol}"
            )));
        }
    }
    Ok(annotations)
}

/// A JSON string as Python's json.dumps writes it with `ensure_ascii=True`.
fn json_string(text: &str) -> String {
    let mut out = String::with_capacity(text.len() + 2);
    out.push('"');
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{8}' => out.push_str("\\b"),
            '\u{c}' => out.push_str("\\f"),
            c if (c as u32) < 0x20 || (c as u32) > 0x7e => {
                let mut units = [0u16; 2];
                for unit in c.encode_utf16(&mut units) {
                    out.push_str(&format!("\\u{unit:04x}"));
                }
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

// ---------------------------------------------------------------------------
// VCF records as pysam returns them

struct Definition {
    number: String,
    kind: String,
    description: String,
}

/// `key=value` pairs of a structured header line (`##INFO=<...>`).
fn structured(line: &str) -> HashMap<String, String> {
    let body = line
        .split_once('<')
        .map_or("", |(_, rest)| rest)
        .trim_end_matches('>');
    let mut pairs = HashMap::new();
    let mut chars = body.chars().peekable();
    while chars.peek().is_some() {
        let key: String = chars.by_ref().take_while(|c| *c != '=').collect();
        let mut value = String::new();
        if chars.peek() == Some(&'"') {
            chars.next();
            while let Some(c) = chars.next() {
                match c {
                    '\\' => {
                        if let Some(next) = chars.next() {
                            value.push(next);
                        }
                    }
                    '"' => break,
                    c => value.push(c),
                }
            }
            if chars.peek() == Some(&',') {
                chars.next();
            }
        } else {
            value = chars.by_ref().take_while(|c| *c != ',').collect();
        }
        pairs.insert(key.trim().to_owned(), value);
    }
    pairs
}

struct VcfHeader {
    infos: HashMap<String, Definition>,
    samples: Vec<String>,
}

impl VcfHeader {
    /// `vcf.header.info[key].description.split('Format: ')[1].split('|')`.
    fn csq_fields(&self, key: &str) -> Result<Vec<String>> {
        let definition = self
            .infos
            .get(key)
            .ok_or_else(|| invalid(format!("the VCF header does not define INFO {key}")))?;
        let parts: Vec<&str> = definition.description.split("Format: ").collect();
        let format = parts
            .get(1)
            .ok_or_else(|| invalid(format!("INFO {key}'s description has no \"Format: \"")))?;
        Ok(format.split('|').map(str::to_owned).collect())
    }
}

/// One record's INFO as pysam returns values.
struct Info<'a> {
    header: &'a VcfHeader,
    /// Key to raw text in record order (`None` for a key without `=`).
    raw: Vec<(&'a str, Option<&'a str>)>,
}

impl<'a> Info<'a> {
    /// `record.info.get(key)`: `Err` for a key the header does not define,
    /// `Ok(None)` for one the record lacks (a Flag is then False).
    fn get(&self, key: &str) -> Result<Option<Py>> {
        let definition = self.header.infos.get(key).ok_or_else(|| {
            invalid(format!(
                "INFO {key} is not defined in the VCF header (pysam: Invalid header)"
            ))
        })?;
        let raw = self.raw.iter().find(|(name, _)| *name == key);
        if definition.kind == "Flag" {
            return Ok(Some(Py::Bool(raw.is_some())));
        }
        match raw {
            None => Ok(None),
            Some((_, text)) => typed(definition, text.unwrap_or(".")).map(Some),
        }
    }

    fn get_or(&self, key: &str, default: Py) -> Result<Py> {
        Ok(self.get(key)?.unwrap_or(default))
    }

    /// `record.info.get(key, [np.nan])[0]`.
    fn first_or_nan(&self, key: &str) -> Result<Py> {
        first(self.get_or(key, Py::Tuple(vec![NAN]))?, key)
    }

    /// The script's `_safe_info`.
    fn safe(&self, key: &str, take_first: bool) -> Result<Py> {
        let value = match self.get(key) {
            Err(_) | Ok(None) | Ok(Some(Py::None)) => return Ok(NAN),
            Ok(Some(value)) => value,
        };
        if !take_first {
            return Ok(value);
        }
        let empty = match &value {
            Py::Tuple(items) => items.is_empty(),
            Py::Str(s) => s.is_empty(),
            _ => return Err(invalid(format!("{key}: the value has no length"))),
        };
        if empty { Ok(NAN) } else { first(value, key) }
    }
}

fn typed(definition: &Definition, text: &str) -> Result<Py> {
    let scalar = definition.number == "1";
    let numbers = |parse: &dyn Fn(&str) -> Option<Py>| -> Result<Py> {
        let items = text
            .split(',')
            .map(|item| {
                if item == "." {
                    Ok(Py::None)
                } else {
                    parse(item)
                        .ok_or_else(|| invalid(format!("INFO value {item:?} is not a number")))
                }
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(if scalar && items.len() == 1 {
            items.into_iter().next().unwrap()
        } else {
            Py::Tuple(items)
        })
    };
    match definition.kind.as_str() {
        "Integer" => numbers(&|item| item.parse::<i64>().ok().map(Py::Int)),
        // htslib keeps Float values at single precision.
        "Float" => numbers(&|item| {
            item.parse::<f64>()
                .ok()
                .map(|value| Py::Float(value as f32 as f64))
        }),
        _ if scalar && !text.contains(',') => Ok(Py::Str(text.to_owned())),
        _ => Ok(Py::Tuple(
            text.split(',')
                .map(|item| Py::Str(item.to_owned()))
                .collect(),
        )),
    }
}

/// GT text as pysam joins it: `|` only when every allele after the first is
/// phased, `.` for a missing allele or one beyond the record's alleles.
fn genotype_text(value: &str, alleles_in_record: i64) -> Result<(String, Vec<Option<i64>>)> {
    if value.is_empty() || value == "." {
        return Ok((".".to_owned(), vec![None]));
    }
    let alleles = value
        .split(['/', '|'])
        .map(|allele| match allele {
            "." => Ok(None),
            _ => allele
                .parse::<i64>()
                .map(|allele| (0..alleles_in_record).contains(&allele).then_some(allele))
                .map_err(|_| invalid(format!("GT {value:?} is not a genotype"))),
        })
        .collect::<Result<Vec<_>>>()?;
    let phased = !value.contains('/');
    let text = alleles
        .iter()
        .map(|allele| allele.map_or(".".to_owned(), |a| a.to_string()))
        .collect::<Vec<_>>()
        .join(if phased { "|" } else { "/" });
    Ok((text, alleles))
}

// ---------------------------------------------------------------------------
// Rows and the script's column typing

/// One value of a chunk after the script's collect step.
#[derive(Clone, Debug)]
enum Cell {
    None,
    Int(i64),
    Float(f64),
    Str(String),
}

/// A value of an object (text) column.
#[derive(Clone, Debug)]
enum Object {
    None,
    Float(f64),
    Str(String),
}

enum Column {
    Position(Vec<i64>),
    Single(Vec<f32>),
    Objects(Vec<Object>),
}

impl Column {
    fn len(&self) -> usize {
        match self {
            Column::Position(v) => v.len(),
            Column::Single(v) => v.len(),
            Column::Objects(v) => v.len(),
        }
    }

    /// Unassigned rows keep the array's initial value.
    fn fill_to(&mut self, length: usize) {
        match self {
            Column::Position(v) => v.resize(length, 0),
            Column::Single(v) => v.resize(length, f32::NAN),
            Column::Objects(v) => v.resize(length, Object::None),
        }
    }

    fn to_objects(&mut self) {
        if let Column::Single(values) = self {
            *self = Column::Objects(
                values
                    .iter()
                    .map(|value| Object::Float(*value as f64))
                    .collect(),
            );
        }
    }
}

/// The script's per-value collect step: missing values become NaN, `pos`
/// an integer, anything `float()` accepts a float (which marks the column as
/// numeric), a tuple its only item or its largest one.
fn collect(position: bool, raw: Option<&Py>, key: &str) -> Result<(Cell, bool)> {
    let value = match raw {
        Some(value) if !value.is_na() => value.clone(),
        _ => NAN,
    };
    if position && !value.is_na() {
        let integer = match &value {
            Py::Int(i) => Some(*i),
            Py::Digits(d) => d.parse().ok(),
            Py::Float(x) => Some(x.trunc() as i64),
            Py::Bool(b) => Some(i64::from(*b)),
            Py::Str(s) => py_strip(s).parse().ok(),
            _ => None,
        };
        if let Some(integer) = integer {
            return Ok((Cell::Int(integer), true));
        }
    }
    if let Some(x) = value.to_float() {
        return Ok((Cell::Float(x), !x.is_nan()));
    }
    let value = match value {
        Py::Tuple(items) => {
            let chosen = if items.len() == 1 {
                items.into_iter().next().unwrap()
            } else {
                let present: Vec<Py> = items.into_iter().filter(|item| *item != Py::None).collect();
                if present.is_empty() {
                    NAN
                } else {
                    py_max(present, key)?
                }
            };
            if let Some(x) = chosen.to_float() {
                return Ok((Cell::Float(x), false));
            }
            chosen
        }
        other => other,
    };
    Ok((
        match value {
            Py::Str(s) => Cell::Str(s),
            Py::None => Cell::None,
            other => return Err(invalid(format!("{key}: unexpected value {other:?}"))),
        },
        false,
    ))
}

/// The table being built: the script's `col_dicts`.
struct Annotations {
    names: Vec<String>,
    index: HashMap<String, usize>,
    columns: Vec<Column>,
    count: usize,
}

/// Rows gathered since the last chunk was stored.
struct Chunk {
    /// Columns of the chunk's first row, with each one's values.
    keys: Vec<String>,
    values: Vec<Vec<Cell>>,
    numeric: Vec<bool>,
    rows: usize,
}

impl Annotations {
    fn store(&mut self, chunk: &mut Chunk) -> Result<()> {
        let start = self.count;
        for (k, key) in chunk.keys.iter().enumerate() {
            let index = *self.index.get(key).ok_or_else(|| {
                invalid(format!(
                    "column {key} is not in the table (the script's KeyError)"
                ))
            })?;
            let cells = std::mem::take(&mut chunk.values[k]);
            let column = &mut self.columns[index];
            column.fill_to(start);
            // A numeric chunk turns an object column without text into float32.
            let single = match &*column {
                Column::Objects(values)
                    if chunk.numeric[k]
                        && key != "pos"
                        && values.iter().all(|value| !matches!(value, Object::Str(_))) =>
                {
                    Some(
                        values
                            .iter()
                            .map(|value| match value {
                                Object::Float(x) => *x as f32,
                                _ => f32::NAN,
                            })
                            .collect(),
                    )
                }
                _ => None,
            };
            if let Some(single) = single {
                *column = Column::Single(single);
            }
            if matches!(column, Column::Single(_))
                && cells.iter().any(|cell| matches!(cell, Cell::Str(_)))
            {
                column.to_objects();
            }
            match column {
                Column::Position(values) => values.extend(cells.iter().map(|cell| match cell {
                    Cell::Int(i) => *i,
                    Cell::Float(x) => *x as i64,
                    _ => 0,
                })),
                Column::Single(values) => values.extend(cells.iter().map(|cell| match cell {
                    Cell::Float(x) => *x as f32,
                    Cell::Int(i) => *i as f32,
                    _ => f32::NAN,
                })),
                Column::Objects(values) => {
                    values.extend(cells.into_iter().map(|cell| match cell {
                        Cell::Float(x) => Object::Float(x),
                        Cell::Int(i) => Object::Float(i as f64),
                        Cell::Str(s) => Object::Str(s),
                        Cell::None => Object::None,
                    }))
                }
            }
        }
        self.count += chunk.rows;
        chunk.keys.clear();
        chunk.values.clear();
        chunk.numeric.clear();
        chunk.rows = 0;
        Ok(())
    }

    fn write(&self, path: &Path) -> Result<()> {
        enum Final<'a> {
            Integers(&'a [i64]),
            Whole(Vec<f32>),
            Single(Vec<f32>),
            Objects(&'a [Object]),
        }
        let finals: Vec<Final> = self
            .columns
            .iter()
            .zip(&self.names)
            .map(|(column, name)| {
                let single = match column {
                    Column::Position(values) if name == "pos" => return Final::Integers(values),
                    Column::Position(values) => values.iter().map(|v| *v as f32).collect(),
                    Column::Single(values) => values.clone(),
                    Column::Objects(values) => {
                        if values.iter().any(|value| matches!(value, Object::Str(_))) {
                            return Final::Objects(values);
                        }
                        values
                            .iter()
                            .map(|value| match value {
                                Object::Float(x) => *x as f32,
                                _ => f32::NAN,
                            })
                            .collect()
                    }
                };
                let whole = single
                    .iter()
                    .filter(|value| !value.is_nan())
                    .all(|value| value.is_finite() && value.fract() == 0.0);
                if whole {
                    Final::Whole(single)
                } else {
                    Final::Single(single)
                }
            })
            .collect();
        let mut out = BufWriter::new(File::create(path)?);
        let header: Vec<String> = self.names.iter().map(|name| csv_field(name)).collect();
        writeln!(out, "{}", header.join("\t"))?;
        let mut line = String::new();
        for row in 0..self.count {
            line.clear();
            for (index, column) in finals.iter().enumerate() {
                if index > 0 {
                    line.push('\t');
                }
                match column {
                    Final::Integers(values) => {
                        line.push_str(&values.get(row).copied().unwrap_or(0).to_string())
                    }
                    Final::Whole(values) => {
                        let value = values.get(row).copied().unwrap_or(f32::NAN);
                        if !value.is_nan() {
                            line.push_str(&wrap_int32(value).to_string());
                        }
                    }
                    Final::Single(values) => {
                        let value = values.get(row).copied().unwrap_or(f32::NAN);
                        if !value.is_nan() {
                            line.push_str(&numpy_f32_repr(value));
                        }
                    }
                    Final::Objects(values) => match values.get(row).unwrap_or(&Object::None) {
                        Object::None => {}
                        Object::Float(x) if x.is_nan() => {}
                        Object::Float(x) => line.push_str(&python_repr(*x)),
                        Object::Str(s) => line.push_str(&csv_field(s)),
                    },
                }
            }
            line.push('\n');
            out.write_all(line.as_bytes())?;
        }
        out.flush()?;
        Ok(())
    }
}

/// An ordered row: a key keeps its first position when a later part of the
/// row (CSQ, genotypes, CADD, HPO) assigns it again, as a Python dict does.
#[derive(Clone, Default)]
struct Row {
    items: Vec<(String, Py)>,
    position: HashMap<String, usize>,
}

impl Row {
    fn set(&mut self, key: &str, value: Py) {
        match self.position.get(key) {
            Some(&index) => self.items[index].1 = value,
            None => {
                self.position.insert(key.to_owned(), self.items.len());
                self.items.push((key.to_owned(), value));
            }
        }
    }

    fn get(&self, key: &str) -> Option<&Py> {
        self.position.get(key).map(|&index| &self.items[index].1)
    }
}

fn convert_vcf(path: &Path, cadd: &CaddScores, hpo: &HpoAnnotations) -> Result<Annotations> {
    let mut reader = open_reader(path)?;
    let mut buffer = Vec::new();
    let mut infos = HashMap::new();
    let samples: Vec<String>;
    let mut number = 0usize;
    loop {
        buffer.clear();
        if reader.read_until(b'\n', &mut buffer)? == 0 {
            return Err(invalid("the VCF has no #CHROM line"));
        }
        number += 1;
        let line = std::str::from_utf8(&buffer)
            .map_err(|_| invalid(format!("line {number}: not UTF-8 text")))?
            .trim_end_matches(['\n', '\r']);
        if line.starts_with("#CHROM") {
            samples = line.split('\t').skip(9).map(str::to_owned).collect();
            break;
        }
        if let Some(rest) = line.strip_prefix("##INFO=") {
            let pairs = structured(rest);
            let field = |key: &str| pairs.get(key).cloned().unwrap_or_default();
            infos.insert(
                field("ID"),
                Definition {
                    number: field("Number"),
                    kind: field("Type"),
                    description: field("Description"),
                },
            );
        }
    }
    let mut header = VcfHeader { infos, samples };
    let csq_fields = header.csq_fields("CSQ")?;
    let feature_index = csq_fields.iter().position(|field| field == "Feature");
    if !csq_fields.iter().any(|field| field == "SYMBOL") {
        return Err(invalid("CSQ has no SYMBOL field"));
    }
    header.csq_fields("CLNCSQ")?;
    let variant_source = header.infos.contains_key("VARIANT_SOURCE");
    let control = header.infos.contains_key("control_AC");

    let mut table: Option<Annotations> = None;
    let mut chunk = Chunk {
        keys: Vec::new(),
        values: Vec::new(),
        numeric: Vec::new(),
        rows: 0,
    };
    loop {
        buffer.clear();
        if reader.read_until(b'\n', &mut buffer)? == 0 {
            break;
        }
        number += 1;
        let line = std::str::from_utf8(&buffer)
            .map_err(|_| invalid(format!("line {number}: not UTF-8 text")))?
            .trim_end_matches(['\n', '\r']);
        if line.is_empty() {
            continue;
        }
        // htslib defines an INFO key the header lacks when a record first
        // uses it ("assuming Type=String"); later reads then succeed.
        if let Some(info) = line.split('\t').nth(7).filter(|info| *info != ".") {
            for item in info.split(';') {
                let key = item.split_once('=').map_or(item, |(key, _)| key);
                if !header.infos.contains_key(key) {
                    header.infos.insert(
                        key.to_owned(),
                        Definition {
                            number: "1".to_owned(),
                            kind: "String".to_owned(),
                            description: "Dummy".to_owned(),
                        },
                    );
                }
            }
        }
        let rows = record_rows(
            line,
            &header,
            &csq_fields,
            feature_index,
            variant_source,
            control,
            cadd,
            hpo,
        )
        .map_err(|error| invalid(format!("line {number}: {error}")))?;
        if rows.is_empty() {
            continue;
        }
        let table = table.get_or_insert_with(|| {
            let names: Vec<String> = rows[0].items.iter().map(|(key, _)| key.clone()).collect();
            Annotations {
                index: names
                    .iter()
                    .enumerate()
                    .map(|(i, n)| (n.clone(), i))
                    .collect(),
                columns: names
                    .iter()
                    .map(|name| {
                        if name == "pos" {
                            Column::Position(Vec::new())
                        } else {
                            Column::Objects(Vec::new())
                        }
                    })
                    .collect(),
                names,
                count: 0,
            }
        });
        for row in &rows {
            if chunk.rows == 0 {
                chunk.keys = row.items.iter().map(|(key, _)| key.clone()).collect();
                chunk.values = vec![Vec::new(); chunk.keys.len()];
                chunk.numeric = vec![false; chunk.keys.len()];
            }
            for (k, key) in chunk.keys.iter().enumerate() {
                let (cell, numeric) = collect(key == "pos", row.get(key), key)?;
                chunk.values[k].push(cell);
                chunk.numeric[k] |= numeric;
            }
            chunk.rows += 1;
        }
        if chunk.rows >= CHUNK_SIZE {
            table.store(&mut chunk)?;
        }
    }
    // The script fails when no record has an Ensembl feature (col_dicts stays None).
    let mut table = table.ok_or_else(|| invalid("no record has an Ensembl feature in CSQ"))?;
    if chunk.rows > 0 {
        table.store(&mut chunk)?;
    }
    let count = table.count;
    for column in &mut table.columns {
        if column.len() < count {
            column.fill_to(count);
        }
    }
    Ok(table)
}

#[allow(clippy::too_many_arguments)]
fn record_rows(
    line: &str,
    header: &VcfHeader,
    csq_fields: &[String],
    feature_index: Option<usize>,
    variant_source: bool,
    control: bool,
    cadd: &CaddScores,
    hpo: &HpoAnnotations,
) -> Result<Vec<Row>> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 8 {
        return Err(invalid("fewer than 8 VCF columns"));
    }
    let (chrom, pos, reference, alts) = (fields[0], fields[1], fields[3], fields[4]);
    let pos: i64 = pos.parse().map_err(|_| invalid("POS is not an integer"))?;
    let alt = if alts == "." {
        ""
    } else {
        alts.split(',').next().unwrap_or("")
    };
    let info = Info {
        header,
        raw: if fields[7] == "." {
            Vec::new()
        } else {
            fields[7]
                .split(';')
                .map(|item| match item.split_once('=') {
                    Some((key, value)) => (key, Some(value)),
                    None => (item, None),
                })
                .collect()
        },
    };

    // extract_record_info, in the script's order (its errors stop the run).
    let mut base = Row::default();
    base.set("chrom", Py::Str(chrom.to_owned()));
    base.set("pos", Py::Int(pos));
    base.set("ref", Py::Str(reference.to_owned()));
    base.set("alt", Py::Str(alt.to_owned()));
    let _ac_joint = info.first_or_nan("AC_joint")?;
    let an_joint = info.get_or("AN_joint", NAN)?;
    let af_joint = info.first_or_nan("AF_joint")?;
    let _nhomalt_joint = info.first_or_nan("nhomalt_joint")?;
    let _ac_grpmax = info.first_or_nan("AC_grpmax_joint")?;
    let af_grpmax = info.first_or_nan("AF_grpmax_joint")?;
    let an_grpmax = info.first_or_nan("AN_grpmax_joint")?;
    let nhomalt_grpmax = info.first_or_nan("nhomalt_grpmax_joint")?;
    let sex = [
        info.safe("AC_joint_XX", true)?,
        info.safe("AN_joint_XX", false)?,
        info.safe("AF_joint_XX", true)?,
        info.safe("AC_joint_XY", true)?,
        info.safe("AN_joint_XY", false)?,
        info.safe("AF_joint_XY", true)?,
        info.safe("nhomalt_joint_XX", true)?,
        info.safe("nhomalt_joint_XY", true)?,
    ];
    let empty = || Py::Tuple(vec![Py::Str(String::new())]);
    comma_join(info.get_or("CLNDN", empty())?, "CLNDN")?;
    let clnsig = comma_join(info.get_or("CLNSIG", empty())?, "CLNSIG")?;
    let clnrevstat = comma_join(info.get_or("CLNREVSTAT", empty())?, "CLNREVSTAT")?;
    let csq = info.get_or("CSQ", empty())?;
    // pext: every PEXT_ key the record carries, first value of a tuple; the
    // script skips the rest after any failure.
    let mut pext = Vec::new();
    for (key, _) in &info.raw {
        if *key == "PEXT_MEAN" || key.starts_with("PEXT_") {
            let value = match info.get_or(key, NAN) {
                Ok(Py::Tuple(items)) => items.into_iter().next().unwrap_or(NAN),
                Ok(value) => value,
                Err(_) => break,
            };
            pext.push((key.to_string(), value));
        }
    }
    let source = if variant_source {
        Some(info.get_or("VARIANT_SOURCE", Py::Str(String::new()))?)
    } else {
        None
    };
    let controls = if control {
        Some([
            info.first_or_nan("control_AC")?,
            info.get_or("control_AN", NAN)?,
            info.first_or_nan("control_AF")?,
            info.get_or("control_nhomalt", NAN)?,
        ])
    } else {
        None
    };
    for population in POPULATIONS {
        info.first_or_nan(&format!("AF_joint_{population}_XX"))?;
        info.first_or_nan(&format!("AF_joint_{population}_XY"))?;
        info.first_or_nan(&format!("nhomalt_joint_{population}_XX"))?;
    }
    let filters = if fields[6] == "." {
        String::new()
    } else {
        fields[6].replace(';', ",")
    };
    let (genotypes, depths) = sample_columns(&fields, &header.samples)?;

    base.set("gnomAD_joint_AN", an_joint);
    base.set("gnomAD_joint_AF", af_joint);
    base.set("gnomAD_joint_AF_max", af_grpmax);
    base.set("gnomAD_joint_AN_max", an_grpmax);
    base.set("gnomAD_nhomalt_max", nhomalt_grpmax);
    let sex_columns = [
        "gnomAD_joint_AC_XX",
        "gnomAD_joint_AN_XX",
        "gnomAD_joint_AF_XX",
        "gnomAD_joint_AC_XY",
        "gnomAD_joint_AN_XY",
        "gnomAD_joint_AF_XY",
        "gnomAD_nhomalt_XX",
        "gnomAD_nhomalt_XY",
    ];
    for (column, value) in sex_columns.iter().zip(sex) {
        base.set(column, value);
    }
    base.set("CLNSIG", Py::Str(clnsig));
    base.set("CLNREVSTAT", Py::Str(clnrevstat));
    base.set("VCF_filters", Py::Str(filters));
    if let Some(source) = source {
        base.set("VARIANT_SOURCE", source);
    }
    if let Some(controls) = controls {
        for (column, value) in ["control_AC", "control_AN", "control_AF", "control_nhomalt"]
            .iter()
            .zip(controls)
        {
            base.set(column, value);
        }
    }
    for (key, value) in pext {
        base.set(&key, value);
    }

    // parse_csq_field: feature name to its fields, first position kept.
    let entries: Vec<String> = match csq {
        Py::Tuple(items) if !items.is_empty() => items
            .into_iter()
            .map(|item| match item {
                Py::Str(s) => Ok(s),
                _ => Err(invalid("CSQ is not text")),
            })
            .collect::<Result<_>>()?,
        Py::Str(s) if !s.is_empty() => s.chars().map(String::from).collect(),
        _ => {
            return Err(invalid(
                "CSQ is empty (the script's parse_csq_field returns None)",
            ));
        }
    };
    let feature_index = feature_index.ok_or_else(|| invalid("CSQ has no Feature field"))?;
    let mut features: Vec<(String, Vec<(String, Py)>)> = Vec::new();
    for entry in &entries {
        let mut values: Vec<&str> = entry.split('|').collect();
        if values.len() < csq_fields.len() {
            values.resize(csq_fields.len(), "");
        }
        let feature = values[feature_index];
        let mut dict = Row::default();
        for (name, value) in csq_fields.iter().zip(&values) {
            let value = if !value.is_empty() && value.bytes().all(|b| b.is_ascii_digit()) {
                Py::Digits((*value).to_owned())
            } else if ["", ".", "-", "NaN", "nan", "NA", "N/A"].contains(value) {
                NAN
            } else {
                Py::Str((*value).to_owned())
            };
            dict.set(name, value);
        }
        if !feature.is_empty() {
            match features.iter_mut().find(|(name, _)| *name == feature) {
                Some((_, existing)) => *existing = dict.items,
                None => features.push((feature.to_owned(), dict.items)),
            }
        }
    }

    let variant = format!("{chrom}:{pos}:{reference}-{alt}");
    let scores = cadd.get(&variant);
    let mut rows = Vec::new();
    for (feature, dict) in features {
        if !feature.starts_with("ENS") {
            continue;
        }
        let mut row = base.clone();
        for (name, value) in &dict {
            row.set(name, value.clone());
        }
        for (sample, genotype) in &genotypes {
            row.set(sample, genotype.clone());
        }
        for (column, depth) in &depths {
            row.set(column, depth.clone());
        }
        let as_py = |score: f32| Py::Float(score as f64);
        match scores.and_then(|scores| scores.get(&feature)) {
            None => {
                row.set("CADD_phred", NAN);
                row.set("CADD_reg_phred", NAN);
            }
            Some(&score) if feature.starts_with("ENST") => {
                row.set("CADD_phred", as_py(score));
                row.set("CADD_reg_phred", NAN);
            }
            Some(&score) => {
                row.set("CADD_reg_phred", as_py(score));
                row.set("CADD_phred", NAN);
            }
        }
        let symbol = dict
            .iter()
            .find(|(name, _)| name == "SYMBOL")
            .map(|(_, value)| value);
        let annotation = match symbol {
            Some(Py::Str(symbol)) => hpo.get(symbol),
            _ => None,
        };
        for (index, column) in HPO_COLUMNS.iter().enumerate() {
            row.set(
                column,
                annotation.map_or(NAN, |values| values[index].clone()),
            );
        }
        rows.push(row);
    }
    Ok(rows)
}

type SampleValues = Vec<(String, Py)>;

/// GT and AD of every sample, as extract_record_info builds them.
fn sample_columns(fields: &[&str], samples: &[String]) -> Result<(SampleValues, SampleValues)> {
    if samples.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }
    let format: Vec<&str> = fields.get(8).map_or(Vec::new(), |f| f.split(':').collect());
    let at = |key: &str| format.iter().position(|name| *name == key);
    let gt = at("GT").ok_or_else(|| invalid("the record has no GT (pysam: invalid FORMAT: GT)"))?;
    let alleles_in_record = match fields[4] {
        "." => 1,
        alts => 1 + alts.split(',').count() as i64,
    };
    let (ad, dp) = (at("AD"), at("DP"));
    let mut genotypes = Vec::new();
    let mut depths = Vec::new();
    for (index, sample) in samples.iter().enumerate() {
        let values: Vec<&str> = fields
            .get(9 + index)
            .map_or(Vec::new(), |value| value.split(':').collect());
        let value = |position: usize| values.get(position).copied().unwrap_or(".");
        let (text, alleles) = genotype_text(value(gt), alleles_in_record)?;
        genotypes.push((sample.clone(), Py::Str(text)));
        let depth = match ad {
            Some(ad) => value(ad)
                .split(',')
                .map(|count| match count {
                    "." => Ok(".".to_owned()),
                    _ => count
                        .parse::<i64>()
                        .map(|count| count.to_string())
                        .map_err(|_| invalid(format!("AD {count:?} is not an integer"))),
                })
                .collect::<Result<Vec<_>>>()?
                .join(","),
            None => {
                let depth = dp
                    .map(value)
                    .filter(|text| *text != ".")
                    .and_then(|text| text.parse::<i64>().ok())
                    .filter(|depth| *depth >= 0);
                match depth {
                    Some(depth) if alleles == [Some(1), Some(1)] => format!("0,{depth}"),
                    _ => ".,.".to_owned(),
                }
            }
        };
        depths.push((format!("{sample}_AD"), Py::Str(depth)));
    }
    Ok((genotypes, depths))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn floats_print_like_python_and_numpy() {
        let python = [
            (0.1, "0.1"),
            (1e-4, "0.0001"),
            (1e-5, "1e-05"),
            (1e16, "1e+16"),
            (9999999999999998.0, "9999999999999998.0"),
            (1e15, "1000000000000000.0"),
            (0.1 + 0.2, "0.30000000000000004"),
            (123.0, "123.0"),
            (-0.0, "-0.0"),
            (0.1f32 as f64, "0.10000000149011612"),
            (129.142857f32 as f64, "129.14285278320312"),
        ];
        for (value, text) in python {
            assert_eq!(python_repr(value), text);
        }
        let numpy = [
            (0.1f32, "0.1"),
            (5.0, "5.0"),
            (1e-5, "1e-05"),
            (123456789.0, "123456790.0"),
            (0.123456789, "0.12345679"),
            (1e16, "1e+16"),
            (2.5e-7, "2.5e-07"),
            (1e-4, "1e-04"),
            (0.00011, "0.00011"),
            (9.99e15, "9990000000000000.0"),
            (-0.0, "-0.0"),
            (3.4e38, "3.4e+38"),
        ];
        for (value, text) in numpy {
            assert_eq!(numpy_f32_repr(value), text);
        }
        assert_eq!(wrap_int32(3e9), -2147483648);
        assert_eq!(wrap_int32(2147483647.0), -2147483648);
        assert_eq!(wrap_int32(-7.0), -7);
    }

    #[test]
    fn python_float_accepts_what_python_accepts() {
        assert_eq!(py_float(" 1_000.5 "), Some(1000.5));
        assert_eq!(py_float("-Infinity"), Some(f64::NEG_INFINITY));
        assert_eq!(py_float("1e-3"), Some(0.001));
        assert_eq!(py_float(".5"), Some(0.5));
        assert!(py_float("1__0").is_none());
        assert!(py_float("_1").is_none());
        assert!(py_float("0/1").is_none());
        assert!(py_float("+-1").is_none());
        assert!(py_float("0x10").is_none());
        assert!(py_float("").is_none());
    }

    #[test]
    fn json_matches_python_ensure_ascii() {
        assert_eq!(
            json_string("a\"b\\c\u{e9}\u{1f600}\t"),
            "\"a\\\"b\\\\c\\u00e9\\ud83d\\ude00\\t\""
        );
    }

    #[test]
    fn csv_quotes_like_pandas() {
        assert_eq!(csv_field("a\tb"), "\"a\tb\"");
        assert_eq!(csv_field("[{\"a\":1}]"), "\"[{\"\"a\"\":1}]\"");
        assert_eq!(csv_field("plain"), "plain");
    }

    #[test]
    fn genotypes_join_like_pysam() {
        assert_eq!(genotype_text("0/1", 2).unwrap().0, "0/1");
        assert_eq!(genotype_text("1|0", 2).unwrap().0, "1|0");
        assert_eq!(genotype_text("0|1/2", 3).unwrap().0, "0/1/2");
        assert_eq!(genotype_text(".", 2).unwrap().0, ".");
        assert_eq!(genotype_text("./.", 2).unwrap().0, "./.");
        assert_eq!(genotype_text("1", 2).unwrap().0, "1");
        assert_eq!(genotype_text("1|1", 1).unwrap().0, ".|.");
        assert_eq!(genotype_text("0/2", 2).unwrap().1, [Some(0), None]);
    }

    #[test]
    fn collect_follows_the_script() {
        let cell = |value: Py| collect(false, Some(&value), "x").unwrap();
        assert!(matches!(cell(Py::Str("0.5".into())), (Cell::Float(x), true) if x == 0.5));
        assert!(matches!(cell(Py::Str("0/1".into())), (Cell::Str(_), false)));
        assert!(matches!(cell(Py::Str("-".into())), (Cell::Float(x), false) if x.is_nan()));
        assert!(matches!(cell(Py::Digits("12".into())), (Cell::Float(x), true) if x == 12.0));
        assert!(matches!(
            cell(Py::Tuple(vec![Py::None, Py::Int(4), Py::Int(2)])),
            (Cell::Float(x), false) if x == 4.0
        ));
        assert!(matches!(
            cell(Py::Tuple(vec![Py::None])),
            (Cell::None, false)
        ));
    }
}
