//! `scripts/gene_mechanism_common.py` (431 lines): the constants and text helpers the hub's
//! modules share.
//!
//! Thread T2 (gene-condition hub). Converted: the constants the ACMG step reads, the default data
//! files (found under `protein_domain_mapping::data_directory()`), and the helpers `_clean`,
//! `_bool_value`, `_split_annotation_tokens`, `_bounded_mechanism_score`,
//! `_exact_mechanisms_from_effect`, the `_norm*` key helpers, `_genomic_match_key(_from_text)`,
//! `_normalize_assembly`, `_requested_genomic_key_types`, `_normalize_sequence_mechanism`,
//! `_split_multi`, `_split_pmids`, `_safe_int`, `_safe_float` and `_deduplicate_by`. The DDG2P
//! constants serve only the hub's diagnostic summaries (group 4) and are not here.
//!
//! Python reads a value of any type through `str(value)`; here a table value is a
//! `pandas_table::Cell` and a cache value a `json_text::Json`. Lists and dicts where Python expects
//! a scalar are refused: Python would print their `repr()`, which the caches never hold.
use super::types::{AcmgError, HarnessCall, not_converted};
use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use crate::pandas_table::{Cell, Column, Table};
use crate::python_text::{decimal_digit_value, float_repr, parse_float, python_strip};
use std::collections::HashSet;
use std::hash::Hash;

/// The fields of the HGNC table the local resolver looks a query up in, best first.
pub(crate) const LOOKUP_FIELD_PRIORITY: [&str; 9] = [
    "symbol",
    "hgnc_id",
    "ensembl_gene_id",
    "entrez_id",
    "prev_symbol",
    "alias_symbol",
    "refseq_accession",
    "uniprot_ids",
    "mane_select",
];

pub(crate) const CANONICAL_MECHANISMS: [&str; 4] =
    ["LOF", "GOF", "DOMINANT_NEGATIVE", "TRIPLOSENSITIVITY"];
pub(crate) const UNRESOLVED_MECHANISM: &str = "UNRESOLVED";
/// Mechanisms an exact allele call can establish exclusively.
pub(crate) const EXACT_SEQUENCE_MECHANISMS: [&str; 3] = ["GOF", "DOMINANT_NEGATIVE", "LOF"];
/// The subset a curated allele match can establish, sorted as `sorted(...)` gives them.
pub(crate) const CURATED_EXACT_MECHANISMS: [&str; 2] = ["DOMINANT_NEGATIVE", "GOF"];
pub(crate) const VARIANT_MECHANISM_SCORE_KEYS: [&str; 3] = ["LOF", "GOF", "DOMINANT_NEGATIVE"];
/// Every source that may state a condition's mechanism.
pub(crate) const CONDITION_MECHANISM_SOURCES: [&str; 5] = [
    "G2P_DDG2P",
    "Orphadata",
    "ClinGen_haploinsufficiency",
    "GoFCards_exact+ClinVar_VCV",
    "deduced_from_inheritance",
];
pub(crate) const HPO_CONDITION_MECHANISM_SCHEMA_VERSION: &str = "1.2";

/// HPO inheritance terms and their labels.
pub(crate) const HPO_INHERITANCE_TERMS: [(&str, &str); 13] = [
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

/// The condition cache's inheritance keys and the HPO labels they stand for.
pub(crate) const HPO_CACHE_INHERITANCE_LABELS: [(&str, &str); 13] = [
    ("autosomal_dominant", "Autosomal dominant inheritance"),
    ("autosomal_recessive", "Autosomal recessive inheritance"),
    ("x_linked", "X-linked inheritance"),
    ("x_linked_recessive", "X-linked recessive inheritance"),
    ("x_linked_dominant", "X-linked dominant inheritance"),
    ("non_mendelian", "Non-Mendelian inheritance"),
    ("mitochondrial", "Mitochondrial inheritance"),
    ("y_linked", "Y-linked inheritance"),
    ("polygenic", "Polygenic inheritance"),
    ("oligogenic", "Oligogenic inheritance"),
    ("digenic", "Digenic inheritance"),
    (
        "autosomal_dominant_maternal_imprinting",
        "Autosomal dominant inheritance with maternal imprinting",
    ),
    (
        "pseudoautosomal_recessive",
        "Pseudoautosomal recessive inheritance",
    ),
];

/// The 15 per-row columns the hub writes after `var_plausible_patho_mechs`, in order.
pub(crate) const VARIANT_MECHANISM_OUTPUT_COLUMNS: [&str; 15] = [
    "variant_effect",
    "variant_lof_score",
    "variant_gof_score",
    "variant_dn_score",
    "variant_mechanism_exclusive",
    "variant_exact_mechanisms",
    "variant_mechanism_applicable",
    "variant_mechanism_uncertain",
    "variant_condition_ids",
    "variant_condition_histories",
    "variant_inheritance",
    "variant_inheritance_basis",
    "variant_x_linked",
    "variant_penetrance",
    "gene_lof_mechanism_history",
];

/// `DEFAULT_HGNC_TABLE`, below the data folder.
pub(crate) const DEFAULT_HGNC_TABLE: &str = "hgnc/hgnc_non_alt_loci_set_sha256_890b18743cc8.tsv";
/// `DEFAULT_CLINGEN_DOSAGE`: always the hg19 file, as in Python.
pub(crate) const DEFAULT_CLINGEN_DOSAGE: &str = "clingen/gene_dosage_sensitivity.hg19.tsv";
pub(crate) const DEFAULT_LOEUF_TABLE: &str = "loeuf/loeuf_dataset.tsv.gz";
pub(crate) const DEFAULT_HPO_CONDITION_MECHANISM_CACHE: &str =
    "patho_mechanism/hpo_condition_mechanism_cache.json.gz";
pub(crate) const DEFAULT_MECHANISM_JSON: &str =
    "patho_mechanism/gene_nonlof_mechanism_curated_assertions.json.gz";

/// A default file of the hub: `DATA_DIR / relative`, with the data folder of
/// `protein_domain_mapping::data_directory()` (`PRIVA_DATA_DIR`, else the source tree's `data/`).
pub(crate) fn default_data_path(relative: &str) -> String {
    crate::protein_domain_mapping::data_directory()
        .join(relative)
        .to_string_lossy()
        .into_owned()
}

/// A value Python handles in a way this conversion does not reproduce.
pub(crate) fn refused(what: impl std::fmt::Display) -> Box<dyn std::error::Error> {
    invalid(format!("refused (not reproduced): {what}"))
}

/// An exception Python raises; the step stops with it.
pub(crate) fn python_error(
    python_type: &'static str,
    message: impl Into<String>,
) -> Box<dyn std::error::Error> {
    Box::new(AcmgError::new(python_type, message))
}

/// Python's `str(value)` of a table value (`nan` for a missing one).
pub(crate) fn python_str(cell: Cell<'_>) -> String {
    match cell {
        Cell::Missing => "nan".to_owned(),
        Cell::Text(text) => text.to_owned(),
        Cell::Int(number) => number.to_string(),
        Cell::Float(number) => float_repr(number),
        Cell::Bool(flag) => (if flag { "True" } else { "False" }).to_owned(),
    }
}

/// `repr()` of a text, for the messages Python builds with `!r`. Only printable
/// ASCII is reproduced; other characters are refused.
pub(crate) fn python_text_repr(text: &str) -> Result<String> {
    if !text.chars().all(|c| (' '..='~').contains(&c)) {
        return Err(refused(format!("repr() of the text {text:?}")));
    }
    let quote = if text.contains('\'') && !text.contains('"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::with_capacity(text.len() + 2);
    out.push(quote);
    for c in text.chars() {
        if c == '\\' || c == quote {
            out.push('\\');
        }
        out.push(c);
    }
    out.push(quote);
    Ok(out)
}

/// `str(Path(text))`: repeated and trailing slashes and `.` parts dropped.
pub(crate) fn python_path_text(text: &str) -> String {
    let root = if text.starts_with("//") && !text.starts_with("///") {
        "//"
    } else if text.starts_with('/') {
        "/"
    } else {
        ""
    };
    let parts: Vec<&str> = text
        .split('/')
        .filter(|part| !part.is_empty() && *part != ".")
        .collect();
    if root.is_empty() && parts.is_empty() {
        ".".to_owned()
    } else {
        format!("{root}{}", parts.join("/"))
    }
}

/// `_clean`: the stripped text, or `""` for a missing value and the spellings
/// `nan`, `none` and `<na>` in any case. `None` stands for an absent value.
pub(crate) fn clean(cell: Option<Cell<'_>>) -> String {
    let Some(cell) = cell else {
        return String::new();
    };
    clean_text(&python_str(cell))
}

/// `_clean` of a text.
pub(crate) fn clean_text(text: &str) -> String {
    let text = python_strip(text);
    // Only ASCII letters lower-case to these spellings, so an ASCII comparison is Python's.
    if ["", "nan", "none", "<na>"]
        .iter()
        .any(|word| text.eq_ignore_ascii_case(word))
    {
        String::new()
    } else {
        text.to_owned()
    }
}

/// `type(value).__name__` of a JSON value, for the messages of Python's errors.
pub(crate) fn json_type_name(value: &Json) -> &'static str {
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

/// Python truthiness of a JSON value.
pub(crate) fn truthy(value: &Json) -> bool {
    match value {
        Json::Null => false,
        Json::Bool(flag) => *flag,
        Json::Int(number) => *number != 0,
        Json::Float(number) => *number != 0.0,
        Json::Str(text) => !text.is_empty(),
        Json::Array(items) => !items.is_empty(),
        Json::Object(items) => !items.is_empty(),
    }
}

/// A JSON scalar as the table value `_clean` and the resolver see; `None` for
/// JSON null or an absent key. Lists and dicts are refused (Python would use
/// their `repr()`), as are integers past 64 bits.
pub(crate) fn json_cell(value: Option<&Json>) -> Result<Option<Cell<'_>>> {
    Ok(match value {
        None | Some(Json::Null) => None,
        Some(Json::Bool(flag)) => Some(Cell::Bool(*flag)),
        Some(Json::Int(number)) => Some(Cell::Int(
            i64::try_from(*number).map_err(|_| refused(format!("the integer {number}")))?,
        )),
        Some(Json::Float(number)) if number.is_nan() => Some(Cell::Missing),
        Some(Json::Float(number)) => Some(Cell::Float(*number)),
        Some(Json::Str(text)) => Some(Cell::Text(text)),
        Some(other) => {
            return Err(refused(format!(
                "a {} where the cache holds text",
                json_type_name(other)
            )));
        }
    })
}

/// `_clean(value)` of a JSON value (`None` for an absent key).
pub(crate) fn clean_json(value: Option<&Json>) -> Result<String> {
    Ok(clean(json_cell(value)?))
}

/// `_bool_value`: True, a non-zero integer, or the texts `1`, `true`, `t`,
/// `yes`, `y` in any case. A float (other than NaN, which is false) goes
/// through its text, which is never one of those spellings.
pub(crate) fn bool_value(cell: Option<Cell<'_>>) -> bool {
    match cell {
        None | Some(Cell::Missing) => false,
        Some(Cell::Bool(flag)) => flag,
        Some(Cell::Int(number)) => number != 0,
        Some(other) => {
            let text = clean(Some(other));
            ["1", "true", "t", "yes", "y"]
                .iter()
                .any(|word| text.eq_ignore_ascii_case(word))
        }
    }
}

/// `_split_annotation_tokens`: the upper-cased tokens between `&`, `,`, `;` and
/// `|`, as a set (here a list without repeats, in first-seen order).
pub(crate) fn split_annotation_tokens(cell: Option<Cell<'_>>) -> Vec<String> {
    let mut tokens: Vec<String> = Vec::new();
    for token in clean(cell).split(['&', ',', ';', '|']) {
        let token = python_strip(token);
        if token.is_empty() {
            continue;
        }
        let token = token.to_uppercase();
        if !tokens.contains(&token) {
            tokens.push(token);
        }
    }
    tokens
}

/// `float(text)` of a cleaned text that is not empty: `None` for Python's
/// ValueError. Python also reads decimal digits of other scripts (`１２` is
/// 12.0): a text holding one is refused.
fn python_float(text: &str) -> Result<Option<f64>> {
    match parse_float(text) {
        Some(number) => Ok(Some(number)),
        None if text
            .chars()
            .any(|c| !c.is_ascii() && decimal_digit_value(c).is_some()) =>
        {
            Err(refused(format!(
                "float() of {text:?}, which holds a non-ASCII digit"
            )))
        }
        None => Ok(None),
    }
}

fn infinity_to_integer() -> Box<dyn std::error::Error> {
    python_error("OverflowError", "cannot convert float infinity to integer")
}

/// `_bounded_mechanism_score`: `int(float(_clean(value)))` when it is 0, 1 or 2,
/// else 0; text Python cannot read as a float gives 0. An infinite value raises
/// OverflowError in Python, and here.
pub(crate) fn bounded_mechanism_score(cell: Option<Cell<'_>>) -> Result<i64> {
    let number = match cell {
        None | Some(Cell::Missing) | Some(Cell::Bool(_)) => return Ok(0),
        Some(Cell::Int(number)) => return Ok(if (0..=2).contains(&number) { number } else { 0 }),
        Some(Cell::Float(number)) => number,
        Some(Cell::Text(text)) => {
            let text = clean_text(text);
            match python_float(&text)? {
                Some(number) => number,
                None => return Ok(0),
            }
        }
    };
    if number.is_nan() {
        return Ok(0);
    }
    if number.is_infinite() {
        return Err(infinity_to_integer());
    }
    let truncated = number.trunc();
    Ok(
        if truncated == 0.0 || truncated == 1.0 || truncated == 2.0 {
            truncated as i64
        } else {
            0
        },
    )
}

/// `_exact_mechanisms_from_effect`: the exact mechanisms an `exact_known_...`
/// effect names, as a set (here a list without repeats).
pub(crate) fn exact_mechanisms_from_effect(cleaned_effect: &str) -> Vec<String> {
    let Some(rest) = cleaned_effect.strip_prefix("exact_known_") else {
        return Vec::new();
    };
    let mut mechanisms: Vec<String> = Vec::new();
    for token in rest.split('+') {
        let mechanism = normalize_sequence_mechanism(&clean_text(token));
        if EXACT_SEQUENCE_MECHANISMS.contains(&mechanism.as_str())
            && !mechanisms.contains(&mechanism)
        {
            mechanisms.push(mechanism);
        }
    }
    mechanisms
}

/// `_norm`: the cleaned text, upper-cased.
pub(crate) fn norm(cell: Option<Cell<'_>>) -> String {
    clean(cell).to_uppercase()
}

/// `_norm_chrom` of a cleaned text: `chr` added when missing.
pub(crate) fn norm_chrom(cleaned: &str) -> String {
    if cleaned.is_empty() || cleaned.starts_with("chr") {
        cleaned.to_owned()
    } else {
        format!("chr{cleaned}")
    }
}

/// `_norm_chrom_key` of a cleaned text.
pub(crate) fn norm_chrom_key(cleaned: &str) -> String {
    if cleaned.is_empty() {
        return String::new();
    }
    // Only the ASCII letters lower-case to "chr", so the first three characters decide.
    let mut text = cleaned;
    let head: String = cleaned.chars().take(3).collect();
    if head.eq_ignore_ascii_case("chr") {
        text = &cleaned[head.len()..];
    }
    let upper = text.to_uppercase();
    if upper == "M" || upper == "MT" {
        "MT".to_owned()
    } else {
        upper
    }
}

/// `str(int(value))` for an integral float.
fn integral_float_text(value: f64) -> String {
    if value == 0.0 {
        "0".to_owned()
    } else {
        // Fixed-point formatting prints the exact integer value.
        format!("{value:.0}")
    }
}

/// `_norm_pos_key` of a cleaned text.
pub(crate) fn norm_pos_key(cleaned: &str) -> Result<String> {
    if cleaned.is_empty() {
        return Ok(String::new());
    }
    match python_float(cleaned)? {
        Some(number) if number.is_finite() && number.fract() == 0.0 => {
            Ok(integral_float_text(number))
        }
        _ => Ok(cleaned.to_owned()),
    }
}

/// `_norm_allele_key` of a cleaned text.
pub(crate) fn norm_allele_key(cleaned: &str) -> String {
    cleaned.to_uppercase()
}

/// `_genomic_match_key` of cleaned texts.
pub(crate) fn genomic_match_key(
    chrom: &str,
    pos: &str,
    ref_allele: &str,
    alt_allele: &str,
) -> Result<String> {
    let chrom_key = norm_chrom_key(chrom);
    let pos_key = norm_pos_key(pos)?;
    if chrom_key.is_empty() || pos_key.is_empty() {
        return Ok(String::new());
    }
    Ok([
        chrom_key,
        pos_key,
        norm_allele_key(ref_allele),
        norm_allele_key(alt_allele),
    ]
    .join("|"))
}

/// `_genomic_match_key_from_text` of a cleaned text (`chrom|pos|ref|alt`).
pub(crate) fn genomic_match_key_from_text(cleaned: &str) -> Result<String> {
    if cleaned.is_empty() {
        return Ok(String::new());
    }
    let parts: Vec<&str> = cleaned.split('|').collect();
    if parts.len() != 4 {
        return Ok(String::new());
    }
    genomic_match_key(
        &clean_text(parts[0]),
        &clean_text(parts[1]),
        &clean_text(parts[2]),
        &clean_text(parts[3]),
    )
}

/// `_normalize_assembly` of a cleaned text.
pub(crate) fn normalize_assembly(cleaned: &str) -> String {
    let lower = cleaned.to_lowercase();
    match lower.as_str() {
        "grch37" | "hg19" | "b37" => "hg19".to_owned(),
        "grch38" | "hg38" | "b38" => "hg38".to_owned(),
        _ => lower,
    }
}

/// `_requested_genomic_key_types` of a cleaned text.
pub(crate) fn requested_genomic_key_types(cleaned: &str) -> &'static [&'static str] {
    match cleaned.to_lowercase().as_str() {
        "" | "auto" => &["vcf", "genomic"],
        "vcf" => &["vcf"],
        "genomic" => &["genomic"],
        _ => &[],
    }
}

/// `_norm_hgvsp` of a cleaned text: the protein change after the last `:`,
/// without `p.`, upper-cased.
pub(crate) fn norm_hgvsp(cleaned: &str) -> String {
    if cleaned.is_empty() {
        return String::new();
    }
    let protein = python_strip(cleaned.rsplit(':').next().unwrap_or(""));
    protein.strip_prefix("p.").unwrap_or(protein).to_uppercase()
}

/// `_normalize_sequence_mechanism` of a cleaned text: runs of characters other
/// than `A-Z` and `0-9` (after upper-casing) become one `_`, outer `_` dropped,
/// then the aliases.
pub(crate) fn normalize_sequence_mechanism(cleaned: &str) -> String {
    let upper = cleaned.to_uppercase();
    let mut token = String::with_capacity(upper.len());
    let mut in_gap = false;
    for c in upper.chars() {
        if c.is_ascii_uppercase() || c.is_ascii_digit() {
            token.push(c);
            in_gap = false;
        } else if !in_gap {
            token.push('_');
            in_gap = true;
        }
    }
    let token = token.trim_matches('_');
    match token {
        "GAIN_OF_FUNCTION" | "ACTIVATING" => "GOF".to_owned(),
        "DN" | "DOMINANTNEGATIVE" => "DOMINANT_NEGATIVE".to_owned(),
        _ => token.to_owned(),
    }
}

/// `_split_multi`: the stripped, non-empty parts between `|`.
pub(crate) fn split_multi(cleaned: &str) -> Vec<String> {
    cleaned
        .split('|')
        .map(python_strip)
        .filter(|part| !part.is_empty())
        .map(str::to_owned)
        .collect()
}

/// `_split_pmids`: the stripped, non-empty parts between `;`, `|` and `,`.
pub(crate) fn split_pmids(cleaned: &str) -> Vec<String> {
    cleaned
        .split([';', '|', ','])
        .map(python_strip)
        .filter(|part| !part.is_empty())
        .map(str::to_owned)
        .collect()
}

/// `_safe_int`: `int(float(text))` of the cleaned value, `None` when it is
/// empty or not a number (or NaN). Infinity raises OverflowError as in Python;
/// a value past 127 bits is refused.
pub(crate) fn safe_int(cell: Option<Cell<'_>>) -> Result<Option<i128>> {
    let text = clean(cell);
    if text.is_empty() {
        return Ok(None);
    }
    let Some(number) = python_float(&text)? else {
        return Ok(None);
    };
    if number.is_nan() {
        return Ok(None);
    }
    if number.is_infinite() {
        return Err(infinity_to_integer());
    }
    let truncated = number.trunc();
    if truncated.abs() >= 1.7e38 {
        return Err(refused(format!("int() of {text:?} past 127 bits")));
    }
    // A float holds its integer value exactly, so the cast loses nothing.
    Ok(Some(truncated as i128))
}

/// `_safe_float`: `float(text)` of the cleaned value, NaN when it is empty or not a number.
pub(crate) fn safe_float(cell: Option<Cell<'_>>) -> Result<f64> {
    let text = clean(cell);
    if text.is_empty() {
        return Ok(f64::NAN);
    }
    Ok(python_float(&text)?.unwrap_or(f64::NAN))
}

/// `_deduplicate_by(records, fields)`: the first record of each key, in order.
/// `key` gives the compared fields as Python builds them (a bool kept as a bool,
/// anything else through `_clean`).
pub(crate) fn deduplicate_by<T, K: Hash + Eq>(records: Vec<T>, key: impl Fn(&T) -> K) -> Vec<T> {
    let mut seen: HashSet<K> = HashSet::new();
    records
        .into_iter()
        .filter(|record| seen.insert(key(record)))
        .collect()
}

/// One field of a `_deduplicate_by` key.
#[derive(Clone, Debug, Hash, PartialEq, Eq)]
pub(crate) enum KeyField {
    Flag(bool),
    Text(String),
}

/// `_deduplicate_by` on JSON records (dicts), for the harness.
fn deduplicate_json(records: &[Json], fields: &[String]) -> Result<Vec<Json>> {
    let mut keyed = Vec::with_capacity(records.len());
    for record in records {
        let Json::Object(_) = record else {
            return Err(python_error(
                "AttributeError",
                format!("'{}' object has no attribute 'get'", json_type_name(record)),
            ));
        };
        let mut key = Vec::with_capacity(fields.len());
        for field in fields {
            key.push(match record.get(field) {
                Some(Json::Bool(flag)) => KeyField::Flag(*flag),
                value => KeyField::Text(clean_json(value)?),
            });
        }
        keyed.push((key, record.clone()));
    }
    Ok(deduplicate_by(keyed, |(key, _)| key.clone())
        .into_iter()
        .map(|(_, record)| record)
        .collect())
}

// ------------------------------------------------------------------------------- harness

/// How the harness prints one Python result (the capture script's `show`):
/// `s:` text, `b:` bool, `i:` int, `f:` float repr, `None`, `j:` JSON of a list
/// (a set sorted first), `!Type: message` for an exception Python raised.
pub(crate) fn show_text(text: &str) -> String {
    format!("s:{text}")
}

pub(crate) fn show_list(items: &[String], sort: bool) -> String {
    let mut items = items.to_vec();
    if sort {
        items.sort();
    }
    let json = Json::Array(items.into_iter().map(Json::Str).collect());
    format!(
        "j:{}",
        json_text::dumps(&json, false).expect("texts always dump")
    )
}

pub(crate) fn show_json(value: &Json) -> Result<String> {
    Ok(format!(
        "j:{}",
        json_text::dumps(value, true).map_err(invalid)?
    ))
}

/// A row result, or the exception Python raises for that row; refusals stop the case.
pub(crate) fn show_outcome(outcome: Result<String>) -> Result<String> {
    match outcome {
        Ok(text) => Ok(text),
        Err(error) => match error.downcast_ref::<AcmgError>() {
            Some(python) => Ok(format!("!{}: {}", python.python_type, python.message)),
            None => Err(error),
        },
    }
}

/// The harness output: one text column `result`.
pub(crate) fn result_table(results: Vec<String>) -> Table {
    Table {
        names: vec!["result".to_owned()],
        rows: results.len(),
        columns: vec![Column::Text(results.into_iter().map(Some).collect())],
    }
}

/// A text cell of a harness frame holding JSON (`json.dumps` of a Python value).
pub(crate) fn json_column_value(table: &Table, name: &str, row: usize) -> Result<Json> {
    match table.column(name)?.cell(row) {
        Cell::Text(text) => json_text::parse(text, false).map_err(invalid),
        _ => Err(invalid(format!("the column {name} must hold JSON text"))),
    }
}

/// `priva acmg-step gene_mechanism_common <function>`: each row's `value` (and
/// the other inputs a function takes) given to the helper, each result shown as
/// the capture script shows Python's.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    const FUNCTIONS: [&str; 21] = [
        "_clean",
        "_bool_value",
        "_split_annotation_tokens",
        "_bounded_mechanism_score",
        "_exact_mechanisms_from_effect",
        "_norm",
        "_norm_chrom",
        "_norm_chrom_key",
        "_norm_pos_key",
        "_norm_allele_key",
        "_genomic_match_key_from_text",
        "_genomic_match_key",
        "_normalize_assembly",
        "_requested_genomic_key_types",
        "_norm_hgvsp",
        "_normalize_sequence_mechanism",
        "_split_multi",
        "_split_pmids",
        "_safe_int",
        "_safe_float",
        "_deduplicate_by",
    ];
    if !FUNCTIONS.contains(&function) {
        return Err(not_converted("gene_mechanism_common", function));
    }
    let table = &call.table;
    let value =
        |row: usize| -> Result<Option<Cell<'_>>> { Ok(Some(table.column("value")?.cell(row))) };
    let mut results = Vec::with_capacity(table.rows);
    for row in 0..table.rows {
        let outcome: Result<String> = match function {
            "_clean" => value(row).map(|cell| show_text(&clean(cell))),
            "_bool_value" => value(row).map(|cell| format!("b:{}", python_bool(bool_value(cell)))),
            "_split_annotation_tokens" => {
                value(row).map(|cell| show_list(&split_annotation_tokens(cell), true))
            }
            "_bounded_mechanism_score" => {
                value(row).and_then(|cell| Ok(format!("i:{}", bounded_mechanism_score(cell)?)))
            }
            "_exact_mechanisms_from_effect" => {
                value(row).map(|cell| show_list(&exact_mechanisms_from_effect(&clean(cell)), true))
            }
            "_norm" => value(row).map(|cell| show_text(&norm(cell))),
            "_norm_chrom" => value(row).map(|cell| show_text(&norm_chrom(&clean(cell)))),
            "_norm_chrom_key" => value(row).map(|cell| show_text(&norm_chrom_key(&clean(cell)))),
            "_norm_pos_key" => {
                value(row).and_then(|cell| Ok(show_text(&norm_pos_key(&clean(cell))?)))
            }
            "_norm_allele_key" => value(row).map(|cell| show_text(&norm_allele_key(&clean(cell)))),
            "_genomic_match_key_from_text" => value(row)
                .and_then(|cell| Ok(show_text(&genomic_match_key_from_text(&clean(cell))?))),
            "_genomic_match_key" => (|| {
                let text = |name: &str| -> Result<String> {
                    Ok(clean(Some(table.column(name)?.cell(row))))
                };
                Ok(show_text(&genomic_match_key(
                    &text("chrom")?,
                    &text("pos")?,
                    &text("ref")?,
                    &text("alt")?,
                )?))
            })(),
            "_normalize_assembly" => {
                value(row).map(|cell| show_text(&normalize_assembly(&clean(cell))))
            }
            "_requested_genomic_key_types" => value(row).map(|cell| {
                let types: Vec<String> = requested_genomic_key_types(&clean(cell))
                    .iter()
                    .map(|name| (*name).to_owned())
                    .collect();
                show_list(&types, false)
            }),
            "_norm_hgvsp" => value(row).map(|cell| show_text(&norm_hgvsp(&clean(cell)))),
            "_normalize_sequence_mechanism" => {
                value(row).map(|cell| show_text(&normalize_sequence_mechanism(&clean(cell))))
            }
            "_split_multi" => value(row).map(|cell| show_list(&split_multi(&clean(cell)), false)),
            "_split_pmids" => value(row).map(|cell| show_list(&split_pmids(&clean(cell)), false)),
            "_safe_int" => value(row).and_then(|cell| {
                Ok(match safe_int(cell)? {
                    Some(number) => format!("i:{number}"),
                    None => "None".to_owned(),
                })
            }),
            "_safe_float" => {
                value(row).and_then(|cell| Ok(format!("f:{}", float_repr(safe_float(cell)?))))
            }
            "_deduplicate_by" => (|| {
                let records = json_column_value(table, "records", row)?;
                let fields = json_column_value(table, "fields", row)?;
                let (Json::Array(records), Json::Array(fields)) = (records, fields) else {
                    return Err(invalid("records and fields must be JSON lists"));
                };
                let fields: Vec<String> = fields
                    .iter()
                    .map(|field| field.as_str().map(str::to_owned))
                    .collect::<Option<_>>()
                    .ok_or_else(|| invalid("fields must be texts"))?;
                show_json(&Json::Array(deduplicate_json(&records, &fields)?))
            })(),
            _ => return Err(not_converted("gene_mechanism_common", function)),
        };
        results.push(show_outcome(outcome)?);
    }
    Ok(result_table(results))
}

/// Python's spelling of a bool.
pub(crate) fn python_bool(flag: bool) -> &'static str {
    if flag { "True" } else { "False" }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clean_and_bool_value() {
        assert_eq!(clean(Some(Cell::Text("  <NA> "))), "");
        assert_eq!(clean(Some(Cell::Text("\u{1f}x\u{1c}"))), "x");
        assert_eq!(clean(Some(Cell::Float(1e16))), "1e+16");
        assert_eq!(clean(Some(Cell::Bool(true))), "True");
        assert_eq!(clean(Some(Cell::Missing)), "");
        assert_eq!(clean(None), "");
        assert!(bool_value(Some(Cell::Text(" Yes "))));
        assert!(bool_value(Some(Cell::Int(-3))));
        assert!(!bool_value(Some(Cell::Float(1.0))));
        assert!(!bool_value(Some(Cell::Missing)));
    }

    #[test]
    fn keys_and_splits() {
        assert_eq!(norm_chrom("7"), "chr7");
        assert_eq!(norm_chrom("chrX"), "chrX");
        assert_eq!(norm_chrom("CHR1"), "chrCHR1");
        assert_eq!(norm_chrom_key("Chr9"), "9");
        assert_eq!(norm_chrom_key("chrM"), "MT");
        assert_eq!(norm_pos_key("5073770.0").unwrap(), "5073770");
        assert_eq!(norm_pos_key("1e20").unwrap(), "100000000000000000000");
        assert_eq!(norm_pos_key("12a").unwrap(), "12a");
        assert_eq!(norm_pos_key("inf").unwrap(), "inf");
        assert_eq!(
            genomic_match_key_from_text(" chr1 | 5.0 |a|t").unwrap(),
            "1|5|A|T"
        );
        assert_eq!(genomic_match_key_from_text("1|2|3").unwrap(), "");
        assert_eq!(split_multi(" a| |b "), ["a", "b"]);
        assert_eq!(split_pmids("1;2|3, ,4"), ["1", "2", "3", "4"]);
        assert_eq!(norm_hgvsp("ENSP0001.1:p.Val617Phe"), "VAL617PHE");
        assert_eq!(normalize_assembly("GRCh38"), "hg38");
        assert_eq!(requested_genomic_key_types(""), &["vcf", "genomic"]);
        assert!(requested_genomic_key_types("x").is_empty());
        assert_eq!(normalize_sequence_mechanism("gain of function"), "GOF");
        assert_eq!(
            normalize_sequence_mechanism("Dominant-Negative"),
            "DOMINANT_NEGATIVE"
        );
        assert_eq!(normalize_sequence_mechanism("dn"), "DOMINANT_NEGATIVE");
        assert_eq!(
            exact_mechanisms_from_effect("exact_known_GOF+DN+x+GOF"),
            ["GOF", "DOMINANT_NEGATIVE"]
        );
        assert!(exact_mechanisms_from_effect("uncertain").is_empty());
    }

    #[test]
    fn numbers_follow_int_of_float() {
        assert_eq!(bounded_mechanism_score(Some(Cell::Text("1_0"))).unwrap(), 0);
        assert_eq!(bounded_mechanism_score(Some(Cell::Text("2.9"))).unwrap(), 2);
        assert_eq!(bounded_mechanism_score(Some(Cell::Int(3))).unwrap(), 0);
        assert_eq!(bounded_mechanism_score(Some(Cell::Bool(true))).unwrap(), 0);
        assert_eq!(
            bounded_mechanism_score(Some(Cell::Text("Infinity")))
                .unwrap_err()
                .to_string(),
            "OverflowError: cannot convert float infinity to integer"
        );
        assert!(bounded_mechanism_score(Some(Cell::Text("\u{ff12}"))).is_err());
        assert_eq!(safe_int(Some(Cell::Text(" 30 "))).unwrap(), Some(30));
        assert_eq!(safe_int(Some(Cell::Text("-2.7"))).unwrap(), Some(-2));
        assert_eq!(
            safe_int(Some(Cell::Text("Not yet evaluated"))).unwrap(),
            None
        );
        assert_eq!(safe_int(Some(Cell::Text("-nan"))).unwrap(), None);
        assert!(safe_int(Some(Cell::Text("inf"))).is_err());
        assert!(safe_float(Some(Cell::Text("x"))).unwrap().is_nan());
        assert_eq!(safe_float(Some(Cell::Text("0.3"))).unwrap(), 0.3);
    }

    #[test]
    fn deduplicate_keeps_the_first_of_each_key() {
        let records = json_text::parse(
            r#"[{"a": " x ", "b": true}, {"a": "x", "b": true}, {"a": "x", "b": "True"}, {"b": null}]"#,
            false,
        )
        .unwrap();
        let Json::Array(records) = records else {
            unreachable!()
        };
        let kept = deduplicate_json(&records, &["a".to_owned(), "b".to_owned()]).unwrap();
        assert_eq!(kept.len(), 3);
        assert_eq!(python_text_repr("it's").unwrap(), "\"it's\"".to_owned());
        assert_eq!(python_text_repr("a'\"b").unwrap(), r#"'a\'"b'"#);
    }
}
