//! `scripts/acmg_pvs1_null_variant.py` (736 lines at main 6cacc18): PVS1 and its gene-level gate.
//!
//! `pvs1_criteria` is `PVS1_CRITERIA` with the helpers it calls (`span_functional_domains`,
//! `downstream_domain_impact`, `downstream_exon_patho_af`, `exon_patho_af_entry`,
//! `evaluate_start_lost_with_cds`, `identify_alternative_start_codon_genes`,
//! `load_non_lof_gene_set`); `summarize_clinvar_gene_pathogenicity` is the Python function of
//! the same name over the decoded ClinVar dictionaries of `acmg_ps1_pm5_known_variant`.
//!
//! The exon look-ups follow the fixed Python of main 6cacc18 (owner's decision of 2026-10-30):
//! an exon number is looked up under `"<exon>/<exon count>"`, the count taken from the
//! transcript's own first `N/N` key, for the exons downstream of a frameshift or stop-gained
//! variant, exon 1 of a 5' UTR variant and the exons a splicing variant skips; the look-ups by
//! VEP's `EXON` text are unchanged. The exon allele-frequency cache is the plain-float rebuild
//! (owner's decision), read with `pickle_dict_stream` for the table's transcripts only, as is
//! the transcript domain map.
//!
//! `df.apply(span_functional_domains, axis=1)` sees each row as Python objects (`df.values` of a
//! frame of mixed types: Python ints, floats, bools, text, NaN, `None`, `pd.NA`); `RowValue`
//! holds them with Python's truth and `str()` rules. The boolean column algebra is T6's
//! `pandas_masks`. The criterion uses T3's `variant_mechanism_masks` for `is_exact_gof` and the
//! hub's `gene_lof_mechanism_history` column.
//!
//! Python's `truncate_fraction` is `acmg_consequence::truncate_fraction` (thread T1).
//!
//! Refused (exit 1 with a message): a `span_functional_domains` value that is not a boolean
//! (numpy would make a float or text array of it and `~` would raise or act bitwise), text
//! with non-ASCII digits where Python's `int()`, `isdigit()` or `\d` would read them, a range
//! of more than a million exons, and the column types the shared mask rules refuse.
use super::acmg_consequence::truncate_fraction;
use super::acmg_pp3_bp4_bp7_insilico::harness_options;
use super::acmg_pp3_bp4_bp7_insilico::pandas_masks::{
    Mask, bool_input, bool_input_fillna_false, numeric_column, python_str_filled, text_contains,
    text_values,
};
use super::acmg_ps1_pm5_known_variant::{
    ClinvarKnownVariants, python_exception, table_transcripts,
};
use super::acmg_variant_mechanism::variant_mechanism_masks;
use super::alternative_start_codon::{
    AlternativeStart, AlternativeStartCodonDetector, adjust_pvs1_for_start_lost,
};
use super::pickle_dict_stream::{Item, visit_top_level_dict};
use super::types::{AcmgError, HarnessCall, StepOutput, not_converted};
use crate::cli::{Result, invalid};
use crate::interpro_mapping_io::{InterproMapping, load_interpro_mapping};
use crate::pandas_table::{Cell, Column, ObjectValue, ReadOptions, Table, read_table_with};
use crate::pickle_reader::{self, Value};
use crate::protein_domain_mapping::{DomainNormalizer, Functionality};
use crate::python_text::{
    decimal_digit_value, float_repr, python_int, python_int_is_valid, python_strip,
};
use std::collections::{HashMap, HashSet};
use std::io::Read;
use std::path::Path;
use std::rc::Rc;

/// `exon_patho_af_threshold` as PVS1_criteria passes it.
const EXON_PATHO_AF_THRESHOLD: f64 = 0.01;

fn python_error(kind: &'static str, message: impl Into<String>) -> Box<dyn std::error::Error> {
    Box::new(AcmgError::new(kind, message))
}

// ---------------------------------------------------------------- gene summary

/// The set `summarize_clinvar_gene_pathogenicity` returns: gene texts, and NaN when a
/// transcript whose `Gene` is missing qualifies (NaN is truthy, so Python keeps it).
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct PathogenicGenes {
    pub(crate) genes: HashSet<String>,
    pub(crate) holds_nan: bool,
}

/// `dict(zip(df['Feature'], df['Gene']))` for the text transcripts (only text can match the
/// ClinVar keys); the last row of a transcript wins. `None` is a missing (NaN) gene.
pub(crate) fn transcript_gene_map(table: &Table) -> Result<HashMap<String, Option<String>>> {
    let feature = table.column("Feature")?;
    let gene = table.column("Gene")?;
    let mut map = HashMap::new();
    for row in 0..table.rows {
        let Cell::Text(transcript) = feature.cell(row) else {
            continue;
        };
        let gene = match (gene, gene.cell(row)) {
            (Column::Object(values), _) if matches!(values[row], ObjectValue::None) => {
                return Err(invalid("Gene holds None; not reproduced"));
            }
            (_, Cell::Text(text)) => Some(text.to_owned()),
            (_, Cell::Missing) => None,
            (_, other) => {
                return Err(invalid(format!(
                    "Gene holds {other:?}; only text genes are reproduced"
                )));
            }
        };
        map.insert(transcript.to_owned(), gene);
    }
    Ok(map)
}

/// `summarize_clinvar_gene_pathogenicity(transcript_to_gene_map, clinvar_aa_dict,
/// clinvar_splice_dict)`: the genes with a pathogenic or likely pathogenic ClinVar variant at two
/// stars or more, in either dictionary. An empty gene text is falsy and skipped, as in Python.
pub(crate) fn summarize_clinvar_gene_pathogenicity(
    transcript_to_gene: &HashMap<String, Option<String>>,
    known: &ClinvarKnownVariants,
) -> PathogenicGenes {
    let mut found = PathogenicGenes::default();
    for (transcript, gene) in transcript_to_gene {
        if gene.as_deref() == Some("") {
            continue;
        }
        let (amino_acid, splice) = known.gene_level_pathogenic(transcript);
        if amino_acid == Some(true) || splice == Some(true) {
            match gene {
                Some(gene) => {
                    found.genes.insert(gene.clone());
                }
                None => found.holds_nan = true,
            }
        }
    }
    found
}

// ---------------------------------------------------------------- resources

/// The ClinVar exon allele-frequency cache for the table's transcripts:
/// `{transcript: {"<exon>/<count>": (af, ...)}}`, keeping the first item of each tuple.
#[derive(Debug, Default)]
pub(crate) struct ExonAlleleFrequencies {
    transcripts: HashMap<String, TranscriptExons>,
}

#[derive(Debug, Default)]
struct TranscriptExons {
    /// The exon count from the transcript's first key holding a `/`.
    count: Option<String>,
    first_items: HashMap<Rc<str>, f64>,
}

impl ExonAlleleFrequencies {
    /// Reads the pickle (gzip when the name ends in `.gz`) for the transcripts in `wanted`.
    pub(crate) fn load(path: &str, wanted: &HashSet<String>) -> Result<ExonAlleleFrequencies> {
        let mut cache = ExonAlleleFrequencies::default();
        let keep = |transcript: &str| wanted.contains(transcript);
        visit_top_level_dict(path, keep, |transcript, value| {
            let odd = |what: &str| {
                invalid(format!(
                    "{path}: transcript {transcript}: {what}; the exon allele-frequency cache is \
                     read for its installed shape only"
                ))
            };
            let Item::Dict(exons) = &value else {
                return Err(odd("the value is not a dict"));
            };
            let mut decoded = TranscriptExons::default();
            for (key, entry) in &exons.borrow().entries {
                if decoded.count.is_none() && key.contains('/') {
                    decoded.count = Some(key.split('/').nth(1).unwrap_or("").to_owned());
                }
                let first = match entry {
                    Item::Tuple(items) => items.first().cloned(),
                    Item::List(items) => items.borrow().first().cloned(),
                    _ => None,
                };
                let first = match first {
                    Some(Item::Float(number)) => number,
                    Some(Item::Int(number)) => number as f64,
                    _ => return Err(odd("an exon's value does not start with a number")),
                };
                decoded.first_items.insert(key.clone(), first);
            }
            cache.transcripts.insert(transcript.to_string(), decoded);
            Ok(())
        })?;
        Ok(cache)
    }

    /// `float(dict.get(transcript, {}).get(exon_text, (0,))[0])`.
    fn by_exon_text(&self, transcript: RowValue<'_>, exon: RowValue<'_>) -> f64 {
        let (RowValue::Text(transcript), RowValue::Text(exon)) = (transcript, exon) else {
            return 0.0;
        };
        self.transcripts
            .get(transcript)
            .and_then(|exons| exons.first_items.get(exon).copied())
            .unwrap_or(0.0)
    }

    /// `float(exon_patho_af_entry(dict, transcript, exon_number)[0])`.
    fn by_exon_number(&self, transcript: RowValue<'_>, exon_number: i128) -> f64 {
        let RowValue::Text(transcript) = transcript else {
            return 0.0;
        };
        let Some(exons) = self.transcripts.get(transcript) else {
            return 0.0;
        };
        let Some(count) = &exons.count else {
            return 0.0;
        };
        exons
            .first_items
            .get(format!("{exon_number}/{count}").as_str())
            .copied()
            .unwrap_or(0.0)
    }
}

/// `tranx_pp_domain_map[transcript].get('pp_domains', [])` for the table's transcripts (keys
/// without version): `(aa_start, domain path)` of each `(aa_start, aa_end, domain path)`.
#[derive(Debug, Default)]
pub(crate) struct ProteinDomains {
    transcripts: HashMap<String, Vec<(i64, Rc<str>)>>,
}

impl ProteinDomains {
    pub(crate) fn load(path: &str, wanted: &HashSet<String>) -> Result<ProteinDomains> {
        let mut map = ProteinDomains::default();
        let keep = |transcript: &str| wanted.contains(transcript);
        visit_top_level_dict(path, keep, |transcript, value| {
            let odd = |what: &str| {
                invalid(format!(
                    "{path}: transcript {transcript}: {what}; the transcript domain map is read \
                     for its installed shape only"
                ))
            };
            let Item::Dict(entry) = &value else {
                return Err(odd("the value is not a dict"));
            };
            let mut domains = Vec::new();
            if let Some(found) = entry.borrow().get("pp_domains") {
                let Item::List(items) = found else {
                    return Err(odd("pp_domains is not a list"));
                };
                for item in items.borrow().iter() {
                    match item {
                        Item::Tuple(parts) if parts.len() == 3 => match (&parts[0], &parts[2]) {
                            (Item::Int(start), Item::Str(path)) => {
                                domains.push((*start, path.clone()));
                            }
                            _ => return Err(odd("a pp_domains entry is not (int, _, str)")),
                        },
                        _ => return Err(odd("a pp_domains entry is not a 3-tuple")),
                    }
                }
            }
            map.transcripts.insert(transcript.to_string(), domains);
            Ok(())
        })?;
        Ok(map)
    }
}

/// The intolerant domain collection (`all_intolerant_domains`): its text items and whether
/// it is truthy (`if intol_domains:`).
#[derive(Clone, Debug, Default)]
pub(crate) struct IntolerantDomains {
    pub(crate) items: HashSet<String>,
    pub(crate) truthy: bool,
}

impl IntolerantDomains {
    /// The orchestrator's load of `--intolerant_domains_pkl` (through `mmap`; a `.gz` name maps
    /// the compressed bytes and fails in Python, so it is refused here).
    pub(crate) fn load(path: &str) -> Result<IntolerantDomains> {
        if path.ends_with(".gz") {
            return Err(invalid(format!(
                "{path}: Python maps the gzip bytes of a .gz intolerant domain pickle and fails"
            )));
        }
        let value = pickle_reader::load(&pickle_reader::read_file(Path::new(path), false)?)?;
        let items: Vec<Value> = match &value {
            Value::Set(set) => set.borrow().items.clone(),
            Value::FrozenSet(set) => set.items.clone(),
            Value::List(items) => items.borrow().clone(),
            other => {
                return Err(invalid(format!(
                    "{path} holds a {}; a set, frozenset or list of domains is expected",
                    other.type_name()
                )));
            }
        };
        Ok(IntolerantDomains {
            truthy: !items.is_empty(),
            items: items
                .iter()
                .filter_map(|item| item.as_str().map(str::to_owned))
                .collect(),
        })
    }

    fn contains(&self, text: &str) -> bool {
        self.items.contains(text)
    }
}

/// `load_non_lof_gene_set(path)`: the truthy `symbol` of every gene of the gzip JSON.
fn load_non_lof_gene_set(path: &str) -> Result<HashSet<String>> {
    let mut raw = Vec::new();
    std::fs::File::open(path)?.read_to_end(&mut raw)?;
    if raw.is_empty() {
        return Err(python_error(
            "JSONDecodeError",
            "Expecting value: line 1 column 1 (char 0)",
        ));
    }
    if !raw.starts_with(&[0x1f, 0x8b]) {
        return Err(python_error(
            "BadGzipFile",
            format!("Not a gzipped file ({path})"),
        ));
    }
    let mut plain = Vec::new();
    flate2::read::MultiGzDecoder::new(&raw[..]).read_to_end(&mut plain)?;
    let text = String::from_utf8(plain).map_err(|_| invalid(format!("{path} is not UTF-8")))?;
    // json.load keeps the last of repeated keys; the installed cache has none, so they are refused.
    let document = crate::json_text::parse(&text, true).map_err(invalid)?;
    let crate::json_text::Json::Object(genes) = document else {
        return Err(python_error(
            "AttributeError",
            "the non-LoF mechanism JSON is not an object ('.values')",
        ));
    };
    let mut symbols = HashSet::new();
    for (_, entry) in &genes {
        let crate::json_text::Json::Object(_) = entry else {
            return Err(python_error(
                "AttributeError",
                "a non-LoF mechanism entry is not an object ('.get')",
            ));
        };
        match entry.get("symbol") {
            None | Some(crate::json_text::Json::Null) => {}
            Some(crate::json_text::Json::Str(symbol)) => {
                if !symbol.is_empty() {
                    symbols.insert(symbol.clone());
                }
            }
            Some(other) => {
                return Err(invalid(format!(
                    "a non-LoF mechanism symbol is {other:?}, not text; not reproduced"
                )));
            }
        }
    }
    Ok(symbols)
}

/// `dict(zip(clingen["#Gene Symbol"], clingen["Haploinsufficiency Score"].astype(int)))` after
/// `dropna`; `None` where Python's `try` catches an exception (the scores then stay NaN).
fn clingen_scores(path: &str) -> Result<Option<HashMap<String, i64>>> {
    let Ok(table) = read_table_with(Path::new(path), &ReadOptions::default()) else {
        return Err(invalid(format!(
            "{path}: the ClinGen dosage file is not read alike by pandas and this reader; Python \
             would log a warning and go on without it"
        )));
    };
    let (Ok(symbols), Ok(scores)) = (
        table.column("#Gene Symbol"),
        table.column("Haploinsufficiency Score"),
    ) else {
        return Ok(None);
    };
    let mut map = HashMap::new();
    for row in 0..table.rows {
        let (symbol, score) = (symbols.cell(row), scores.cell(row));
        if symbol == Cell::Missing || score == Cell::Missing {
            continue;
        }
        let Cell::Text(symbol) = symbol else {
            return Err(invalid(
                "the ClinGen gene symbols are not text; not reproduced",
            ));
        };
        let score = match score {
            Cell::Int(number) => number,
            Cell::Bool(flag) => i64::from(flag),
            Cell::Float(number) if number.is_finite() && number.abs() < 9.0e18 => {
                number.trunc() as i64
            }
            Cell::Float(_) => return Ok(None),
            Cell::Text(text) => match python_int(text) {
                Some(number) if i64::try_from(number).is_ok() => number as i64,
                _ if python_int_is_valid(text) => {
                    return Err(invalid(
                        "a ClinGen score uses digits int() reads; not reproduced",
                    ));
                }
                _ => return Ok(None),
            },
            Cell::Missing => unreachable!("dropped above"),
        };
        map.insert(symbol.to_owned(), score);
    }
    Ok(Some(map))
}

// ---------------------------------------------------------------- row values

/// A value of one row as `df.apply(..., axis=1)` hands it to the function.
#[derive(Clone, Copy, Debug, PartialEq)]
enum RowValue<'a> {
    Text(&'a str),
    Int(i64),
    /// A float; NaN is how pandas holds a missing value in float and text columns.
    Float(f64),
    Bool(bool),
    PythonNone,
    PandasNa,
}

impl RowValue<'_> {
    fn of(column: &Column, row: usize) -> RowValue<'_> {
        match column {
            Column::Int(values) => RowValue::Int(values[row]),
            Column::Float(values) => RowValue::Float(values[row]),
            Column::Bool(values) => RowValue::Bool(values[row]),
            Column::BoolOrMissing(values) => {
                values[row].map_or(RowValue::Float(f64::NAN), RowValue::Bool)
            }
            Column::Text(values) => values[row]
                .as_deref()
                .map_or(RowValue::Float(f64::NAN), RowValue::Text),
            Column::Object(values) => match &values[row] {
                ObjectValue::Nan => RowValue::Float(f64::NAN),
                ObjectValue::None => RowValue::PythonNone,
                ObjectValue::Text(text) => RowValue::Text(text),
                ObjectValue::Int(number) => RowValue::Int(*number),
                ObjectValue::Float(number) => RowValue::Float(*number),
                ObjectValue::Bool(flag) => RowValue::Bool(*flag),
            },
            Column::NullableInt(values) => values[row].map_or(RowValue::PandasNa, RowValue::Int),
            Column::NullableBool(values) => values[row].map_or(RowValue::PandasNa, RowValue::Bool),
        }
    }

    /// `bool(value)`.
    fn truthy(self) -> Result<bool> {
        match self {
            RowValue::Text(text) => Ok(!text.is_empty()),
            RowValue::Int(number) => Ok(number != 0),
            RowValue::Float(number) => Ok(number != 0.0),
            RowValue::Bool(flag) => Ok(flag),
            RowValue::PythonNone => Ok(false),
            RowValue::PandasNa => Err(python_error(
                "TypeError",
                "boolean value of NA is ambiguous",
            )),
        }
    }

    /// `str(value)`.
    fn python_str(self) -> String {
        match self {
            RowValue::Text(text) => text.to_owned(),
            RowValue::Int(number) => number.to_string(),
            RowValue::Float(number) => float_repr(number),
            RowValue::Bool(flag) => (if flag { "True" } else { "False" }).to_owned(),
            RowValue::PythonNone => "None".to_owned(),
            RowValue::PandasNa => "<NA>".to_owned(),
        }
    }

    fn type_name(self) -> &'static str {
        match self {
            RowValue::Text(_) => "str",
            RowValue::Int(_) => "int",
            RowValue::Float(_) => "float",
            RowValue::Bool(_) => "bool",
            RowValue::PythonNone => "NoneType",
            RowValue::PandasNa => "NAType",
        }
    }
}

/// `int(text)`; `None` for Python's ValueError. Text that Python reads through digits of
/// other scripts is refused.
fn python_int_of(text: &str) -> Result<Option<i128>> {
    match python_int(text) {
        Some(number) => Ok(Some(number)),
        None if python_int_is_valid(text) => Err(invalid(format!(
            "int({text:?}) reads non-ASCII digits or a huge number in Python; not reproduced"
        ))),
        None => Ok(None),
    }
}

fn int_or_value_error(text: &str) -> Result<i128> {
    python_int_of(text)?.ok_or_else(|| {
        python_error(
            "ValueError",
            format!("invalid literal for int() with base 10: {text:?}"),
        )
    })
}

/// What one row needs from the loaded resources.
struct RowResources<'a> {
    normalizer: &'a DomainNormalizer,
    interpro: &'a InterproMapping,
    domains: &'a ProteinDomains,
    exon_af: &'a ExonAlleleFrequencies,
    intolerant: &'a IntolerantDomains,
}

impl RowResources<'_> {
    fn functional(&self, domain: &str) -> Result<bool> {
        match self
            .normalizer
            .interpret_functionality(domain, self.interpro)
        {
            Ok(found) => Ok(found == Functionality::Functional),
            Err(message) => Err(python_exception(&message)),
        }
    }
}

struct RowReader<'a> {
    table: &'a Table,
    row: usize,
}

impl<'a> RowReader<'a> {
    /// `row[name]`: KeyError when the table has no such column.
    fn get(&self, name: &str) -> Result<RowValue<'a>> {
        let column = self.table.column(name)?;
        Ok(RowValue::of(column, self.row))
    }
}

/// `downstream_domain_impact(row['Protein_position'], row['Feature'], ...,
/// domains=str(row["DOMAINS"]), intol_domains, ensg_id)`.
fn downstream_domain_impact(
    resources: &RowResources,
    protein_position: RowValue<'_>,
    transcript: RowValue<'_>,
    domains: &str,
    gene: &str,
) -> Result<bool> {
    for domain in domains.split('&') {
        if resources.intolerant.contains(&format!("{gene}:{domain}"))
            || resources.functional(domain)?
        {
            return Ok(true);
        }
    }
    let truncation = match protein_position {
        RowValue::Text(text) if !text.is_empty() => {
            let first = text.split('/').next().unwrap_or("");
            python_int_of(first.split('-').next().unwrap_or(""))?
        }
        _ => None,
    };
    let Some(truncation) = truncation else {
        return Ok(false);
    };
    let RowValue::Text(transcript) = transcript else {
        return Err(python_error(
            "AttributeError",
            format!(
                "'{}' object has no attribute 'split'",
                transcript.type_name()
            ),
        ));
    };
    let base = transcript.split('.').next().unwrap_or("");
    let Some(protein_domains) = resources.domains.transcripts.get(base) else {
        return Ok(false);
    };
    for (start, path) in protein_domains {
        if i128::from(*start) > truncation {
            let Some((_, domain)) = path.split_once(':') else {
                return Err(python_error("IndexError", "list index out of range"));
            };
            if resources.intolerant.contains(&format!("{gene}:{domain}"))
                || resources.functional(domain)?
            {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// `downstream_exon_patho_af(row, dict, logic="any", threshold)`.
fn downstream_exon_patho_af(
    resources: &RowResources,
    row: &RowReader,
    threshold: f64,
) -> Result<bool> {
    let RowValue::Text(exon) = row.get("EXON")? else {
        return Ok(false);
    };
    let (first, last) = if exon.contains('-') && exon.contains('/') {
        (
            exon.split('-').next().unwrap_or(""),
            exon.split('/').nth(1).unwrap_or(""),
        )
    } else if exon.contains('/') {
        (
            exon.split('/').next().unwrap_or(""),
            exon.split('/').nth(1).unwrap_or(""),
        )
    } else {
        return Err(python_error(
            "ValueError",
            format!("Invalid exon string: {exon}"),
        ));
    };
    let (first, last) = (int_or_value_error(first)?, int_or_value_error(last)?);
    if last - first > 1_000_000 {
        return Err(invalid(format!(
            "EXON {exon:?} spans more than a million exons; not reproduced"
        )));
    }
    let transcript = row.get("Feature")?;
    Ok((first..=last).any(|exon| resources.exon_af.by_exon_number(transcript, exon) < threshold))
}

/// The digit-only exon numbers of `str(row["splicing_affected_exons"])`.
fn skipped_exons(text: &str) -> Result<Vec<i128>> {
    let mut exons = Vec::new();
    for part in text.split(',') {
        let part = python_strip(part);
        if part.is_empty() {
            continue;
        }
        if !part.is_ascii() {
            if part
                .chars()
                .all(|c| c.is_numeric() || decimal_digit_value(c).is_some())
            {
                return Err(invalid(format!(
                    "splicing_affected_exons part {part:?} may be digits to str.isdigit(); not \
                     reproduced"
                )));
            }
            continue;
        }
        if part.bytes().all(|byte| byte.is_ascii_digit()) {
            exons.push(int_or_value_error(part)?);
        }
    }
    Ok(exons)
}

/// `span_functional_domains(row, ...)`: `(func_domain, exon_frequent_patho)`. `func_domain` is
/// the raw cell for 5' UTR and splicing variants, as Python takes it.
fn span_functional_domains<'a>(
    resources: &RowResources,
    row: &RowReader<'a>,
    threshold: f64,
) -> Result<(RowValue<'a>, bool)> {
    let gene = row.get("Gene")?.python_str();
    let consequence = match row.get("Consequence")? {
        RowValue::Text(text) => text,
        other => {
            return Err(python_error(
                "TypeError",
                format!("argument of type '{}' is not iterable", other.type_name()),
            ));
        }
    };
    if consequence.contains("frameshift") || consequence.contains("stop_gained") {
        let protein_position = row.get("Protein_position")?;
        let transcript = row.get("Feature")?;
        let domains = row.get("DOMAINS")?.python_str();
        let functional =
            downstream_domain_impact(resources, protein_position, transcript, &domains, &gene)?;
        let frequent = downstream_exon_patho_af(resources, row, threshold)?;
        for name in ["chrom", "pos", "Feature", "EXON"] {
            row.get(name)?;
        }
        return Ok((RowValue::Bool(functional), frequent));
    }
    let any_truthy = |names: [&str; 3]| -> Result<bool> {
        for name in names {
            if row.get(name)?.truthy()? {
                return Ok(true);
            }
        }
        Ok(false)
    };
    if any_truthy(["5UTR_lof", "5UTR_frameshift", "5UTR_len_changing"])? {
        let functional = row.get("5UTR_span_intol_domain")?;
        let af = resources.exon_af.by_exon_number(row.get("Feature")?, 1);
        for name in ["chrom", "pos", "Feature"] {
            row.get(name)?;
        }
        return Ok((functional, af < threshold));
    }
    if any_truthy([
        "splicing_lof",
        "splicing_frameshift",
        "splicing_len_changing",
    ])? {
        let functional = row.get("splicing_span_intol_domain")?;
        let exons = skipped_exons(&row.get("splicing_affected_exons")?.python_str())?;
        let transcript = row.get("Feature")?;
        if !exons.is_empty() {
            let afs: Vec<f64> = exons
                .iter()
                .map(|exon| resources.exon_af.by_exon_number(transcript, *exon))
                .collect();
            for name in ["chrom", "pos", "Feature"] {
                row.get(name)?;
            }
            return Ok((functional, afs.iter().any(|af| *af < threshold)));
        }
        for name in ["chrom", "pos"] {
            row.get(name)?;
        }
        let af = resources
            .exon_af
            .by_exon_text(row.get("Feature")?, row.get("EXON")?);
        return Ok((functional, af < threshold));
    }
    let mut functional = false;
    if let RowValue::Text(domains) = row.get("DOMAINS")? {
        let mut found = Vec::new();
        for domain in domains.split('&') {
            found.push(if resources.intolerant.truthy {
                resources.intolerant.contains(&format!("{gene}:{domain}"))
            } else {
                resources.functional(domain)?
            });
        }
        functional = found.into_iter().any(|flag| flag);
    }
    let af = resources
        .exon_af
        .by_exon_text(row.get("Feature")?, row.get("EXON")?);
    for name in ["chrom", "pos", "Feature", "EXON"] {
        row.get(name)?;
    }
    Ok((RowValue::Bool(functional), af < threshold))
}

/// `zip(*df.apply(span_functional_domains, axis=1, ...))` with `np.array` of each half; the
/// first half must hold booleans only.
fn functional_domains_of_rows(
    table: &Table,
    resources: &RowResources,
) -> Result<(Vec<bool>, Vec<bool>)> {
    if table.rows == 0 {
        return Err(python_error(
            "ValueError",
            "too many values to unpack (expected 2): df.apply on an empty table gives its column \
             names",
        ));
    }
    let mut functional = Vec::with_capacity(table.rows);
    let mut frequent = Vec::with_capacity(table.rows);
    let mut odd: Option<RowValue> = None;
    for row in 0..table.rows {
        let (domain, rare) = span_functional_domains(
            resources,
            &RowReader { table, row },
            EXON_PATHO_AF_THRESHOLD,
        )?;
        match domain {
            RowValue::Bool(flag) => functional.push(flag),
            other => {
                odd.get_or_insert(other);
                functional.push(false);
            }
        }
        frequent.push(rare);
    }
    if let Some(value) = odd {
        return Err(invalid(format!(
            "span_functional_domains gives a {} ({}) taken from 5UTR_span_intol_domain or \
             splicing_span_intol_domain; numpy's array of mixed values and its `~` are not \
             reproduced",
            value.type_name(),
            value.python_str()
        )));
    }
    Ok((functional, frequent))
}

// ---------------------------------------------------------------- PVS1_criteria

/// The arguments of `PVS1_criteria` besides the table. `row_labels` are the pandas index labels
/// of the rows when they are not 0, 1, 2, ... (the step's `read_table().drop_duplicates()`
/// keeps the labels of the rows it keeps): `evaluate_start_lost_with_cds` reads the LoF gate
/// by label (`lof_intolerant[iloc_idx]`).
pub(crate) struct Pvs1Arguments<'a> {
    pub(crate) clinvar_patho_exon_af_stat: &'a str,
    pub(crate) interpro_entry_map_pkl: &'a str,
    pub(crate) clinvar_pathogenic_genes: &'a PathogenicGenes,
    pub(crate) gene_dosage_sensitivity: Option<&'a str>,
    pub(crate) nonlof_mechanism_json: Option<&'a str>,
    pub(crate) intolerant_domains: &'a IntolerantDomains,
    pub(crate) tranx_exon_domain_map_pkl: Option<&'a str>,
    pub(crate) cds_fasta_path: Option<&'a str>,
    pub(crate) row_labels: Option<&'a [i64]>,
}

fn path_exists(path: Option<&str>) -> bool {
    path.is_some_and(|path| !path.is_empty() && Path::new(path).exists())
}

/// `df[name].isin(values)` for a gene or symbol column of text.
fn text_isin(
    table: &Table,
    name: &str,
    values: &HashSet<String>,
    nan_listed: bool,
) -> Result<Vec<bool>> {
    let column = table.column(name)?;
    (0..table.rows)
        .map(|row| match RowValue::of(column, row) {
            RowValue::Text(text) => Ok(values.contains(text)),
            RowValue::Float(number) if number.is_nan() => Ok(nan_listed),
            other => Err(invalid(format!(
                "{name} holds a {}; isin is reproduced for text and NaN only",
                other.type_name()
            ))),
        })
        .collect()
}

fn set_where(pvs1: &mut [i64], mask: &Mask, value: i64) -> Result<()> {
    for (slot, flag) in pvs1.iter_mut().zip(mask.plain()?) {
        if *flag {
            *slot = value;
        }
    }
    Ok(())
}

/// `pvs1_criteria[mask] = value`, the mask built from the strengths as they are before the assignment.
fn assign(pvs1: &mut [i64], value: i64, build: impl FnOnce(&[i64]) -> Mask) -> Result<()> {
    let mask = build(pvs1);
    set_where(pvs1, &mask, value)
}

fn below(pvs1: &[i64], limit: i64) -> Mask {
    Mask::Plain(pvs1.iter().map(|value| *value < limit).collect())
}

/// The three columns `evaluate_start_lost_with_cds` writes into the table, each from
/// `df['Feature'].map(...)` (pandas infers the type from the values).
fn alternative_start_columns(
    table: &Table,
    found: &HashMap<String, AlternativeStart>,
) -> Result<Vec<(String, Column)>> {
    let feature = table.column("Feature")?;
    let lookup = |row: usize| -> Option<&AlternativeStart> {
        match RowValue::of(feature, row) {
            RowValue::Text(text) => found.get(text),
            _ => None,
        }
    };
    let exists: Vec<bool> = (0..table.rows)
        .map(|row| lookup(row).is_some_and(|info| info.has_alternative))
        .collect();
    let truncation: Vec<Option<f64>> = (0..table.rows)
        .map(|row| lookup(row).and_then(|info| info.n_terminal_truncation_fraction))
        .collect();
    let position: Vec<Option<i64>> = (0..table.rows)
        .map(|row| lookup(row).and_then(|info| info.first_alt_position))
        .collect();
    let floats_or_none = |values: Vec<Option<f64>>| -> Column {
        if values.iter().all(Option::is_none) {
            Column::Object(vec![ObjectValue::None; values.len()])
        } else {
            Column::Float(
                values
                    .into_iter()
                    .map(|value| value.unwrap_or(f64::NAN))
                    .collect(),
            )
        }
    };
    let position = if position.iter().all(Option::is_some) {
        Column::Int(
            position
                .into_iter()
                .map(|value| value.unwrap_or(0))
                .collect(),
        )
    } else {
        floats_or_none(
            position
                .into_iter()
                .map(|value| value.map(|number| number as f64))
                .collect(),
        )
    };
    Ok(vec![
        ("alt_start_exists".to_owned(), Column::Bool(exists)),
        (
            "alt_start_truncation".to_owned(),
            floats_or_none(truncation),
        ),
        ("alt_start_position".to_owned(), position),
    ])
}

/// `evaluate_start_lost_with_cds(df, cds_fasta_path, intolerant_domains, lof_mechanism)`: the
/// start_lost strengths and the columns written into the table (none when no variant is
/// start_lost or the detector fails).
fn evaluate_start_lost_with_cds(
    table: &Table,
    cds_fasta_path: &str,
    intolerant_domains: &[bool],
    lof_mechanism: &[bool],
    row_labels: Option<&[i64]>,
) -> Result<(Vec<i64>, Vec<(String, Column)>)> {
    let rows = table.rows;
    let mut strengths = vec![0_i64; rows];
    let start_lost = text_contains(table, "Consequence", None, &["start_lost"])?[0].fillna(false);
    let start_lost = start_lost.plain()?.to_vec();
    if !start_lost.iter().any(|flag| *flag) {
        return Ok((strengths, Vec::new()));
    }
    let feature = table.column("Feature")?;
    // `start_lost_df['Feature'].dropna().unique()`; `find_alternative_start` splits each, so
    // a value that is not text fails there, once the detector has read the FASTA.
    let mut transcripts: Vec<&str> = Vec::new();
    let mut not_text: Option<&'static str> = None;
    let mut seen = HashSet::new();
    for row in (0..rows).filter(|row| start_lost[*row]) {
        match RowValue::of(feature, row) {
            RowValue::Text(text) => {
                if seen.insert(text) {
                    transcripts.push(text);
                }
            }
            RowValue::Float(number) if number.is_nan() => {}
            RowValue::PythonNone | RowValue::PandasNa => {}
            other => {
                not_text.get_or_insert(other.type_name());
            }
        }
    }
    let detector = match AlternativeStartCodonDetector::load(cds_fasta_path, &transcripts)? {
        Ok(detector) => detector,
        Err(_) => {
            for (slot, flag) in strengths.iter_mut().zip(&start_lost) {
                if *flag {
                    *slot = 2;
                }
            }
            return Ok((strengths, Vec::new()));
        }
    };
    table.column("Gene")?;
    if let Some(kind) = not_text {
        return Err(python_error(
            "AttributeError",
            format!("'{kind}' object has no attribute 'split'"),
        ));
    }
    let found: HashMap<String, AlternativeStart> = transcripts
        .iter()
        .map(|id| ((*id).to_owned(), detector.find_alternative_start(id)))
        .collect();
    let columns = alternative_start_columns(table, &found)?;
    for row in (0..rows).filter(|row| start_lost[*row]) {
        let info = match RowValue::of(feature, row) {
            RowValue::Text(text) => found.get(text),
            _ => None,
        };
        // `lof_intolerant[iloc_idx]` reads the Series by label.
        let lof = match row_labels {
            None => lof_mechanism[row],
            Some(labels) => {
                let label = row as i64;
                match labels.iter().position(|candidate| *candidate == label) {
                    Some(position) => lof_mechanism[position],
                    None => return Err(python_error("KeyError", format!("{label}"))),
                }
            }
        };
        strengths[row] = adjust_pvs1_for_start_lost(info, intolerant_domains[row], lof);
    }
    Ok((strengths, columns))
}

/// A hashable key of one cell of a group-by column, or `None` for a missing value (which
/// `groupby` drops).
fn group_key_part(value: RowValue<'_>, name: &str) -> Result<Option<String>> {
    Ok(match value {
        RowValue::Text(text) => Some(format!("s{text}")),
        RowValue::Int(number) => Some(format!("n{number}")),
        RowValue::Bool(flag) => Some(format!("n{}", i64::from(flag))),
        RowValue::Float(number) if number.is_nan() => None,
        RowValue::Float(number) if number.fract() == 0.0 && number.abs() < 9.0e15 => {
            Some(format!("n{}", number as i64))
        }
        RowValue::Float(number) => Some(format!("f{:016x}", number.to_bits())),
        RowValue::PythonNone | RowValue::PandasNa => {
            return Err(invalid(format!(
                "{name} holds None or pd.NA; groupby on it is not reproduced"
            )));
        }
    })
}

/// The deprecated start_lost path (no CDS FASTA): `alt_start_genes` from
/// `df.groupby(['chrom', 'pos', 'ref', 'alt', 'Gene']).apply(identify_alternative_start_codon_genes)`.
fn alternative_start_genes(table: &Table, start_lost: &[bool]) -> Result<PathogenicGenes> {
    let keys = ["chrom", "pos", "ref", "alt", "Gene"];
    let columns = keys
        .iter()
        .map(|name| table.column(name))
        .collect::<Result<Vec<&Column>>>()?;
    let mut groups: HashMap<Vec<String>, Vec<usize>> = HashMap::new();
    'rows: for row in 0..table.rows {
        let mut key = Vec::with_capacity(keys.len());
        for (name, column) in keys.iter().zip(&columns) {
            match group_key_part(RowValue::of(column, row), name)? {
                Some(part) => key.push(part),
                None => continue 'rows,
            }
        }
        groups.entry(key).or_default().push(row);
    }
    if groups.is_empty() {
        return Err(python_error(
            "AttributeError",
            "'DataFrame' object has no attribute 'unique'",
        ));
    }
    let mut biotype_coding: Option<Vec<bool>> = None;
    let mut genes = PathogenicGenes::default();
    for rows in groups.values() {
        let lost = rows.iter().filter(|row| start_lost[**row]).count();
        if lost == 0 || lost == rows.len() {
            genes.holds_nan = true;
            continue;
        }
        if biotype_coding.is_none() {
            let found = text_contains(table, "BIOTYPE", Some("."), &["protein_coding"])?;
            biotype_coding = Some(found[0].plain()?.to_vec());
        }
        let coding = biotype_coding.as_ref().expect("filled above");
        if rows.iter().any(|row| !start_lost[*row] && coding[*row]) {
            match RowValue::of(columns[4], rows[0]) {
                RowValue::Text(gene) => {
                    genes.genes.insert(gene.to_owned());
                }
                other => {
                    return Err(invalid(format!(
                        "Gene holds a {}; only text genes are reproduced",
                        other.type_name()
                    )));
                }
            }
        } else {
            genes.holds_nan = true;
        }
    }
    Ok(genes)
}

/// `PVS1_criteria(df, ...)`: the PVS1 strengths (int64, 0 to 4), the alt_start columns
/// `evaluate_start_lost_with_cds` adds to the table, and the `span_functional_domains` array
/// handed on as `locate_intol_domains`.
pub(crate) fn pvs1_criteria(table: &Table, arguments: &Pvs1Arguments) -> Result<StepOutput> {
    let rows = table.rows;
    let mechanism = variant_mechanism_masks(table)?;

    // 1. The gene's curated LoF history (`.fillna("").astype(str).str.lower().eq("true")`).
    let hpo_lof_history: Vec<bool> = match table.column_index("gene_lof_mechanism_history") {
        None => vec![false; rows],
        Some(index) => python_str_filled(&table.columns[index], "gene_lof_mechanism_history", "")?
            .iter()
            .map(|text| text.is_ascii() && text.eq_ignore_ascii_case("true"))
            .collect(),
    };
    // 2. ClinVar pathogenic gene.
    let genes = arguments.clinvar_pathogenic_genes;
    let clinvar_pathogenic = text_isin(table, "Gene", &genes.genes, genes.holds_nan)?;
    // 3 and 4. LOEUF and mean AlphaMissense.
    let loeuf_intolerant: Vec<bool> = numeric_column(table, "LOEUF")?
        .iter()
        .map(|value| !value.is_nan() && *value < 0.35)
        .collect();
    let am_intolerant: Vec<bool> = if table.column_index("Gene_avg_AM_score").is_some() {
        numeric_column(table, "Gene_avg_AM_score")?
            .iter()
            .map(|value| !value.is_nan() && *value > 0.6)
            .collect()
    } else {
        vec![false; rows]
    };
    // 5. ClinGen dosage: inside Python's try, so a failure leaves the scores NaN.
    let mut hi_scores: Vec<Option<i64>> = vec![None; rows];
    if let Some(path) = arguments.gene_dosage_sensitivity
        && path_exists(Some(path))
        && let Some(scores) = clingen_scores(path)?
        && let Ok(symbol) = table.column("SYMBOL")
    {
        for (row, slot) in hi_scores.iter_mut().enumerate() {
            // Only a text symbol can equal one of the dict's text keys.
            *slot = match RowValue::of(symbol, row) {
                RowValue::Text(text) => scores.get(text).copied(),
                _ => None,
            };
        }
    }
    let clingen_hi: Vec<bool> = hi_scores.iter().map(|score| *score == Some(3)).collect();
    let clingen_ar: Vec<bool> = hi_scores
        .iter()
        .map(|score| matches!(score, Some(30 | 40)))
        .collect();
    let non_lof_gene = match arguments.nonlof_mechanism_json {
        Some(path) if path_exists(Some(path)) => {
            let symbols = load_non_lof_gene_set(path)?;
            text_isin(table, "SYMBOL", &symbols, false)?
        }
        _ => vec![false; rows],
    };
    let lof_mechanism: Vec<bool> = (0..rows)
        .map(|row| {
            let blocked = non_lof_gene[row] && !hpo_lof_history[row] && !clingen_hi[row];
            (hpo_lof_history[row]
                || clinvar_pathogenic[row]
                || loeuf_intolerant[row]
                || am_intolerant[row]
                || clingen_hi[row]
                || clingen_ar[row])
                && !blocked
        })
        .collect();

    // The row-wise part, after the four loads in Python's order.
    let transcripts = if table.has_column("Feature") {
        table_transcripts(table)?
    } else {
        HashSet::new()
    };
    let exon_af = ExonAlleleFrequencies::load(arguments.clinvar_patho_exon_af_stat, &transcripts)?;
    let interpro = load_interpro_mapping(Path::new(arguments.interpro_entry_map_pkl))?;
    let Some(domain_map_path) = arguments.tranx_exon_domain_map_pkl else {
        return Err(python_error(
            "AttributeError",
            "'NoneType' object has no attribute 'endswith'",
        ));
    };
    let bases: HashSet<String> = transcripts
        .iter()
        .map(|id| id.split('.').next().unwrap_or("").to_owned())
        .collect();
    let domains = ProteinDomains::load(domain_map_path, &bases)?;
    let normalizer = DomainNormalizer::new()?;
    let resources = RowResources {
        normalizer: &normalizer,
        interpro: &interpro,
        domains: &domains,
        exon_af: &exon_af,
        intolerant: arguments.intolerant_domains,
    };
    let (intolerant, exon_rare) = functional_domains_of_rows(table, &resources)?;

    let lof = Mask::Plain(lof_mechanism.clone());
    let intol = Mask::Plain(intolerant.clone());
    let not_intol = Mask::Plain(intolerant.iter().map(|flag| !flag).collect());
    let rare = Mask::Plain(exon_rare);
    let mut pvs1 = vec![0_i64; rows];

    // Frameshift and nonsense variants.
    let consequence = text_contains(
        table,
        "Consequence",
        None,
        &["stop_gained", "frameshift", "inframe_deletion"],
    )?;
    let truncating = consequence[0]
        .fillna(false)
        .or(&consequence[1].fillna(false));
    let nmd_escaping = text_contains(table, "NMD", Some("."), &["escaping"])?.remove(0);
    let nmd_variants = truncating.and(&nmd_escaping.logical_not());
    assign(&mut pvs1, 4, |_| nmd_variants.and(&lof))?;
    let end_truncation = text_contains(table, "LoF_filter", Some("."), &["END_TRUNC"])?.remove(0);
    let escaping = truncating.and(&nmd_escaping.or(&end_truncation));
    let escaping_intol = escaping.and(&intol);
    assign(&mut pvs1, 3, |current| {
        escaping_intol.and(&lof).and(&below(current, 3))
    })?;
    let fraction = truncate_fraction(table)?;
    let long_truncation = Mask::Plain(fraction.iter().map(|value| *value >= 0.1).collect());
    let short_truncation = Mask::Plain(fraction.iter().map(|value| *value < 0.1).collect());
    let escaping_base = escaping.and(&lof).and(&not_intol).and(&rare);
    assign(&mut pvs1, 3, |current| {
        escaping_base.and(&long_truncation).and(&below(current, 3))
    })?;
    assign(&mut pvs1, 2, |current| {
        escaping_base.and(&short_truncation).and(&below(current, 2))
    })?;

    // In-frame deletions.
    let inframe = &consequence[2];
    let long_text = |name: &str| -> Result<Mask> {
        Ok(Mask::Plain(
            text_values(table, name, None)?
                .iter()
                .map(|text| text.as_ref().is_some_and(|text| text.chars().count() > 50))
                .collect(),
        ))
    };
    let exon_span = Mask::Plain(
        text_values(table, "EXON", Some(""))?
            .iter()
            .map(|text| exon_range_text(text.as_deref().unwrap_or("")))
            .collect::<Result<Vec<bool>>>()?,
    );
    let large = long_text("ref")?.or(&long_text("alt")?).or(&exon_span);
    assign(&mut pvs1, 3, |_| inframe.and(&intol).and(&large).and(&lof))?;
    let inframe_base = inframe.and(&lof).and(&not_intol).and(&rare).and(&large);
    assign(&mut pvs1, 3, |current| {
        inframe_base.and(&long_truncation).and(&below(current, 3))
    })?;
    assign(&mut pvs1, 2, |current| {
        inframe_base.and(&short_truncation).and(&below(current, 2))
    })?;

    // Splicing variants.
    let splicing_frameshift = bool_input_fillna_false(table, "splicing_frameshift")?;
    let last_exon_only = if table.column_index("splicing_last_exon_only").is_some() {
        bool_input_fillna_false(table, "splicing_last_exon_only")?
    } else {
        Mask::Plain(vec![false; rows])
    };
    let not_last_exon_only = last_exon_only.invert()?;
    assign(&mut pvs1, 4, |current| {
        splicing_frameshift
            .and(&not_last_exon_only)
            .and(&lof)
            .and(&below(current, 4))
    })?;
    let escape = splicing_frameshift.and(&last_exon_only);
    let span = bool_input_fillna_false(table, "splicing_span_intol_domain")?;
    let ten_percent = bool_input_fillna_false(table, "splicing_ten_percent_protein")?;
    assign(&mut pvs1, 3, |current| {
        escape.and(&span).and(&lof).and(&below(current, 3))
    })?;
    let not_span = span.invert()?;
    let not_ten_percent = ten_percent.invert()?;
    assign(&mut pvs1, 3, |current| {
        escape
            .and(&not_span)
            .and(&ten_percent)
            .and(&lof)
            .and(&rare)
            .and(&below(current, 3))
    })?;
    assign(&mut pvs1, 2, |current| {
        escape
            .and(&not_span)
            .and(&not_ten_percent)
            .and(&lof)
            .and(&rare)
            .and(&below(current, 2))
    })?;
    let lof_or_length =
        bool_input(table, "splicing_lof")?.or(&bool_input(table, "splicing_len_changing")?);
    let not_frameshift = splicing_frameshift.invert()?;
    let inframe_intol = lof_or_length
        .and(&not_frameshift)
        .and(&not_last_exon_only)
        .and(&bool_input(table, "splicing_span_intol_domain")?);
    assign(&mut pvs1, 3, |current| {
        inframe_intol.and(&lof).and(&below(current, 3))
    })?;
    let ten_percent_raw = bool_input(table, "splicing_ten_percent_protein")?;
    let splicing_base = lof_or_length.and(&not_frameshift).and(&not_span).and(&rare);
    assign(&mut pvs1, 3, |current| {
        splicing_base
            .and(&ten_percent_raw)
            .and(&lof)
            .and(&below(current, 3))
    })?;
    let not_ten_percent_raw = ten_percent_raw.invert()?;
    assign(&mut pvs1, 2, |current| {
        splicing_base
            .and(&not_ten_percent_raw)
            .and(&lof)
            .and(&below(current, 2))
    })?;

    // 5' UTR frameshifts (UTRAnnotator).
    let utr_frameshift = bool_input(table, "5UTR_frameshift")?;
    let utr_inframe_intol = bool_input(table, "5UTR_len_changing")?
        .and(&bool_input(table, "5UTR_span_intol_domain")?)
        .and(&utr_frameshift.invert()?);
    assign(&mut pvs1, 3, |current| {
        utr_inframe_intol.and(&lof).and(&below(current, 3))
    })?;
    assign(&mut pvs1, 4, |current| {
        utr_frameshift.and(&lof).and(&below(current, 4))
    })?;

    // start_lost variants.
    let mut columns = Vec::new();
    if path_exists(arguments.cds_fasta_path) {
        let (start_lost_strengths, added) = evaluate_start_lost_with_cds(
            table,
            arguments.cds_fasta_path.unwrap_or(""),
            &intolerant,
            &lof_mechanism,
            arguments.row_labels,
        )?;
        columns = added;
        let start_lost =
            text_contains(table, "Consequence", None, &["start_lost"])?[0].fillna(false);
        let start_lost = start_lost.plain()?;
        for row in 0..rows {
            if start_lost[row] && start_lost_strengths[row] > pvs1[row] {
                pvs1[row] = start_lost_strengths[row];
            }
        }
        for row in 0..rows {
            if start_lost[row] && pvs1[row] == 0 {
                pvs1[row] = start_lost_strengths[row];
            }
        }
    } else {
        let start_lost = text_contains(table, "Consequence", None, &["start_lost"])?.remove(0);
        let alt_genes = alternative_start_genes(table, start_lost.plain()?)?;
        let in_alt_genes = Mask::Plain(text_isin(
            table,
            "Gene",
            &alt_genes.genes,
            alt_genes.holds_nan,
        )?);
        let alt_start_losts = start_lost.and(&in_alt_genes);
        let not_alt_start_losts = alt_start_losts.invert()?;
        assign(&mut pvs1, 2, |current| {
            start_lost.and(&not_alt_start_losts).and(&below(current, 2))
        })?;
    }

    let exact_gof = mechanism.get("is_exact_gof")?;
    for (slot, gof) in pvs1.iter_mut().zip(exact_gof) {
        if *slot > 0 && *gof {
            *slot = 0;
        }
    }
    Ok(StepOutput {
        criteria: vec![("PVS1", Column::Int(pvs1))],
        columns,
        passed_on: vec![("locate_intol_domains".to_owned(), Column::Bool(intolerant))],
    })
}

/// `str.match(r'^\d+-\d+/\d+$')` on one text (`$` also before a final line break).
fn exon_range_text(text: &str) -> Result<bool> {
    if !text.is_ascii() {
        if text
            .chars()
            .any(|c| !c.is_ascii() && decimal_digit_value(c).is_some())
        {
            return Err(invalid(format!(
                "EXON {text:?} holds digits of another script, which Python's \\d reads; not \
                 reproduced"
            )));
        }
        return Ok(false);
    }
    let body = text.strip_suffix('\n').unwrap_or(text);
    let Some((range, count)) = body.split_once('/') else {
        return Ok(false);
    };
    let Some((first, last)) = range.split_once('-') else {
        return Ok(false);
    };
    let digits = |part: &str| !part.is_empty() && part.bytes().all(|byte| byte.is_ascii_digit());
    Ok(digits(first) && digits(last) && digits(count))
}

// ---------------------------------------------------------------- harness

/// `priva acmg-step acmg_pvs1_null_variant <function>`:
///
/// - `summarize_clinvar_gene_pathogenicity`, options `clinvar_aa_dict_pkl`,
///   `clinvar_splice_dict_pkl`, on a frame with `Feature` and `Gene`: a frame of `gene` (text,
///   sorted, a missing value last for NaN);
/// - `span_functional_domains`, options `clinvar_patho_exon_af_stat`, `interpro_entry_map_pkl`,
///   `tranx_exon_domain_map_pkl`, `intolerant_domains_pkl` (null for `[]`): a frame of
///   `func_domain` and `exon_frequent_patho` (bool), one row per row;
/// - `PVS1_criteria`, the same options and `clinvar_pathogenic_genes` (list of text, null for
///   NaN), `gene_dosage_sensitivity`, `nonlof_mechanism_json`, `cds_fasta_path` (each text or
///   null): a frame of `PVS1` (int64), `span_functional_domains` (bool) and, when Python adds
///   them to the table, `df_alt_start_exists`, `df_alt_start_truncation`, `df_alt_start_position`.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let table = &call.table;
    let options = &call.options;
    let text = |key: &str| harness_options::text(options, key);
    let intolerant = || -> Result<IntolerantDomains> {
        match text("intolerant_domains_pkl")? {
            Some(path) => IntolerantDomains::load(path),
            None => Ok(IntolerantDomains::default()),
        }
    };
    match function {
        "summarize_clinvar_gene_pathogenicity" => {
            let known = ClinvarKnownVariants::load(
                harness_options::required_text(options, "clinvar_aa_dict_pkl")?,
                harness_options::required_text(options, "clinvar_splice_dict_pkl")?,
                Some(&table_transcripts(table)?),
            )?;
            let found = summarize_clinvar_gene_pathogenicity(&transcript_gene_map(table)?, &known);
            let mut genes: Vec<Option<String>> = found.genes.into_iter().map(Some).collect();
            genes.sort();
            if found.holds_nan {
                genes.push(None);
            }
            let count = genes.len();
            Ok(harness_options::frame(
                vec![("gene".to_owned(), Column::Text(genes))],
                count,
            ))
        }
        "span_functional_domains" => {
            let transcripts = table_transcripts(table)?;
            let exon_af = ExonAlleleFrequencies::load(
                harness_options::required_text(options, "clinvar_patho_exon_af_stat")?,
                &transcripts,
            )?;
            let interpro = load_interpro_mapping(Path::new(harness_options::required_text(
                options,
                "interpro_entry_map_pkl",
            )?))?;
            let bases: HashSet<String> = transcripts
                .iter()
                .map(|id| id.split('.').next().unwrap_or("").to_owned())
                .collect();
            let domains = ProteinDomains::load(
                harness_options::required_text(options, "tranx_exon_domain_map_pkl")?,
                &bases,
            )?;
            let normalizer = DomainNormalizer::new()?;
            let intolerant = intolerant()?;
            let resources = RowResources {
                normalizer: &normalizer,
                interpro: &interpro,
                domains: &domains,
                exon_af: &exon_af,
                intolerant: &intolerant,
            };
            let (functional, frequent) = functional_domains_of_rows(table, &resources)?;
            Ok(harness_options::frame(
                vec![
                    ("func_domain".to_owned(), Column::Bool(functional)),
                    ("exon_frequent_patho".to_owned(), Column::Bool(frequent)),
                ],
                table.rows,
            ))
        }
        "PVS1_criteria" => {
            let genes = match options.get("clinvar_pathogenic_genes") {
                None | Some(crate::json_text::Json::Null) => PathogenicGenes::default(),
                Some(crate::json_text::Json::Array(items)) => {
                    let mut genes = PathogenicGenes::default();
                    for item in items {
                        match item {
                            crate::json_text::Json::Str(gene) => {
                                genes.genes.insert(gene.clone());
                            }
                            crate::json_text::Json::Null => genes.holds_nan = true,
                            _ => return Err(invalid("clinvar_pathogenic_genes holds a non-text")),
                        }
                    }
                    genes
                }
                Some(_) => return Err(invalid("clinvar_pathogenic_genes must be a list")),
            };
            let intolerant = intolerant()?;
            let arguments = Pvs1Arguments {
                clinvar_patho_exon_af_stat: harness_options::required_text(
                    options,
                    "clinvar_patho_exon_af_stat",
                )?,
                interpro_entry_map_pkl: harness_options::required_text(
                    options,
                    "interpro_entry_map_pkl",
                )?,
                clinvar_pathogenic_genes: &genes,
                gene_dosage_sensitivity: text("gene_dosage_sensitivity")?,
                nonlof_mechanism_json: text("nonlof_mechanism_json")?,
                intolerant_domains: &intolerant,
                tranx_exon_domain_map_pkl: text("tranx_exon_domain_map_pkl")?,
                cds_fasta_path: text("cds_fasta_path")?,
                row_labels: None,
            };
            let output = pvs1_criteria(table, &arguments)?;
            let mut columns = vec![
                (
                    "PVS1".to_owned(),
                    output.criterion("PVS1").cloned().expect("PVS1"),
                ),
                (
                    "span_functional_domains".to_owned(),
                    output
                        .passed_on("locate_intol_domains")
                        .cloned()
                        .expect("domains"),
                ),
            ];
            for (name, column) in output.columns {
                columns.push((format!("df_{name}"), column));
            }
            Ok(harness_options::frame(columns, table.rows))
        }
        _ => Err(not_converted("acmg_pvs1_null_variant", function)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exon_numbers_find_the_cache_key_of_their_transcript() {
        let mut exons = TranscriptExons {
            count: Some("36".into()),
            first_items: HashMap::new(),
        };
        for key in ["1/36", "34/36", "35-36/36"] {
            exons.first_items.insert(Rc::from(key), 0.05);
        }
        let mut cache = ExonAlleleFrequencies::default();
        cache.transcripts.insert("ENST1".into(), exons);
        let transcript = RowValue::Text("ENST1");
        assert_eq!(cache.by_exon_number(transcript, 34), 0.05);
        assert_eq!(cache.by_exon_number(transcript, 1), 0.05);
        assert_eq!(cache.by_exon_number(transcript, 2), 0.0);
        assert_eq!(cache.by_exon_number(RowValue::Text("ENST_absent"), 1), 0.0);
        assert_eq!(cache.by_exon_number(RowValue::Float(f64::NAN), 1), 0.0);
        assert_eq!(
            cache.by_exon_text(transcript, RowValue::Text("35-36/36")),
            0.05
        );
        assert_eq!(cache.by_exon_text(transcript, RowValue::Text("2/36")), 0.0);
        assert_eq!(
            cache.by_exon_text(transcript, RowValue::Float(f64::NAN)),
            0.0
        );
    }

    #[test]
    fn exon_range_pattern() {
        assert!(exon_range_text("3-5/10").unwrap());
        assert!(exon_range_text("3-5/10\n").unwrap());
        assert!(!exon_range_text("3/10").unwrap());
        assert!(!exon_range_text("").unwrap());
        assert!(!exon_range_text("3-5/10x").unwrap());
        assert!(exon_range_text("3-٥/10").is_err());
    }

    #[test]
    fn skipped_exons_keep_digit_parts() {
        assert_eq!(skipped_exons("35, 36,,x,007").unwrap(), vec![35, 36, 7]);
        assert_eq!(skipped_exons("nan").unwrap(), Vec::<i128>::new());
        assert_eq!(skipped_exons("5.0").unwrap(), Vec::<i128>::new());
        assert!(skipped_exons("²").is_err());
    }
}
