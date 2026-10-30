//! `scripts/gene_mechanism_variants.py` (642 lines): `infer_query_variant_effect` and the exact
//! curated-allele matcher.
//!
//! Thread T3 (exact-allele matching and the scoring tail). Converted: `infer_query_variant_effect`,
//! which the hub (thread T2) calls for every row, and `ExactVariantMatcher` with its loaders and
//! `match_curated_nonlof_variant`, which `annotate_exact_nonlof_variants` uses. The two
//! compatibility views `match_gofcards_variant_gof*` are not on the step's path.
//!
//! The matcher also holds the hub's `_load_mechanisms` (`mechanisms_by_symbol`), which only it
//! uses in the step. Resolution of symbols goes through `types::SymbolResolver` (thread T2's hub).
//!
//! The text helpers of `gene_mechanism_common.py` this module needs (`_clean`, `_bool_value`,
//! `_split_annotation_tokens`, `_bounded_mechanism_score`) are kept here as private functions,
//! because `gene_mechanism_common` belongs to thread T2; they can move there at integration.
use super::types::{
    AcmgError, EffectInputs, HarnessCall, QueryEffect, SymbolResolver, not_converted,
};
use crate::cli::{Result, invalid};
use crate::json_text::{self, Json, quote};
use crate::pandas_table::{Cell, Column, Table};
use crate::python_text::{float_repr, parse_float, python_strip};
use std::collections::{HashMap, HashSet};
use std::io::Read;

/// `priva acmg-step gene_mechanism_variants <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    match function {
        "infer_query_variant_effect" => infer_frame(&call.table),
        "match_curated_nonlof_variant" => match_frame(call),
        _ => Err(not_converted("gene_mechanism_variants", function)),
    }
}

/// `pd.DataFrame([infer_query_variant_effect(row) for _, row in df.iterrows()])`.
fn infer_frame(table: &Table) -> Result<Table> {
    let column = |name: &str| table.column_index(name);
    let positions = [
        column("LoF"),
        column("Consequence"),
        column("NMD"),
        column("LoF_filter"),
        column("splicing_frameshift"),
        column("vep_consq_lof"),
        column("splicing_lof"),
        column("5UTR_lof"),
        column("variant_lof_score"),
        column("variant_gof_score"),
        column("variant_dn_score"),
        column("variant_gof_tag"),
    ];
    let mut effects = Vec::with_capacity(table.rows);
    for row in 0..table.rows {
        let cell = |index: usize| positions[index].map(|position| table.cell(position, row));
        let inputs = EffectInputs {
            lof: cell(0),
            consequence: cell(1),
            nmd: cell(2),
            lof_filter: cell(3),
            splicing_frameshift: cell(4),
            vep_consq_lof: cell(5),
            splicing_lof: cell(6),
            utr5_lof: cell(7),
            variant_lof_score: cell(8),
            variant_gof_score: cell(9),
            variant_dn_score: cell(10),
            variant_gof_tag: cell(11),
        };
        effects.push(infer_query_variant_effect(&inputs)?);
    }
    let texts = |pick: fn(&QueryEffect) -> &String| {
        Column::Text(effects.iter().map(|e| Some(pick(e).clone())).collect())
    };
    let numbers = |pick: fn(&QueryEffect) -> i64| Column::Int(effects.iter().map(pick).collect());
    let names = [
        "variant_effect",
        "variant_effect_evidence",
        "variant_effect_suppressed_evidence",
        "variant_lof_score",
        "variant_gof_score",
        "variant_dn_score",
        "variant_mechanism_exclusive",
        "variant_exact_mechanisms",
    ];
    // An empty list gives pandas a frame without columns.
    if effects.is_empty() {
        return Ok(Table {
            names: Vec::new(),
            columns: Vec::new(),
            rows: 0,
        });
    }
    let columns = vec![
        texts(|e| &e.variant_effect),
        texts(|e| &e.variant_effect_evidence),
        texts(|e| &e.variant_effect_suppressed_evidence),
        numbers(|e| e.variant_lof_score),
        numbers(|e| e.variant_gof_score),
        numbers(|e| e.variant_dn_score),
        Column::Bool(
            effects
                .iter()
                .map(|e| e.variant_mechanism_exclusive)
                .collect(),
        ),
        texts(|e| &e.variant_exact_mechanisms),
    ];
    Ok(Table {
        names: names.iter().map(|name| (*name).to_owned()).collect(),
        columns,
        rows: effects.len(),
    })
}

/// The symbol resolutions a Python capture recorded (`options["symbols"]`: cleaned
/// symbol to `[hub._try_resolve_symbol(s), hub.resolve_symbol(s)]`), for comparing
/// the matcher before thread T2's hub resolver is in the crate. Harness only.
pub(crate) struct CapturedResolver {
    symbols: HashMap<String, (String, String)>,
}

impl CapturedResolver {
    pub(crate) fn from_options(options: &Json) -> Result<CapturedResolver> {
        let mut symbols = HashMap::new();
        let Some(Json::Object(items)) = options.get("symbols") else {
            return Err(invalid(
                "options need \"symbols\": {symbol: [tried, resolved]}",
            ));
        };
        for (symbol, pair) in items {
            match pair {
                Json::Array(pair) if pair.len() == 2 => match (&pair[0], &pair[1]) {
                    (Json::Str(tried), Json::Str(resolved)) => {
                        symbols.insert(symbol.clone(), (tried.clone(), resolved.clone()));
                    }
                    _ => return Err(invalid("symbol resolutions must be texts")),
                },
                _ => return Err(invalid("each symbol needs [tried, resolved]")),
            }
        }
        Ok(CapturedResolver { symbols })
    }
}

impl SymbolResolver for CapturedResolver {
    fn resolve_symbol(&self, gene_symbol: Cell<'_>) -> Result<String> {
        let cleaned = clean(Some(gene_symbol));
        if cleaned.is_empty() {
            return Ok(String::new());
        }
        self.symbols
            .get(&cleaned)
            .map(|(_, resolved)| resolved.clone())
            .ok_or_else(|| invalid(format!("the symbol {cleaned:?} was not captured")))
    }

    fn try_resolve_symbol(&self, gene_symbol: Cell<'_>) -> String {
        let cleaned = clean(Some(gene_symbol));
        if cleaned.is_empty() {
            return cleaned;
        }
        match self.symbols.get(&cleaned) {
            Some((tried, _)) => tried.clone(),
            // Shows up as a difference rather than passing silently.
            None => format!("<uncaptured {cleaned}>"),
        }
    }
}

/// The text option `name` of a harness call.
pub(crate) fn text_option<'a>(options: &'a Json, name: &str) -> Result<&'a str> {
    options
        .get(name)
        .and_then(Json::as_str)
        .ok_or_else(|| invalid(format!("options need the text {name:?}")))
}

/// One `match_curated_nonlof_variant` call per row of the frame (columns `gene_symbol`, `hgvsp`,
/// `hgvsc`, `chrom`, `pos`, `ref`, `alt`, `assembly`, `key_type`; absent ones take the defaults).
fn match_frame(call: &HarnessCall) -> Result<Table> {
    let resolver = CapturedResolver::from_options(&call.options)?;
    let mut matcher =
        ExactVariantMatcher::new(text_option(&call.options, "mechanism_json")?, &resolver);
    let table = &call.table;
    let column = |name: &str| table.column_index(name);
    let positions = [
        "gene_symbol",
        "hgvsp",
        "hgvsc",
        "chrom",
        "pos",
        "ref",
        "alt",
        "assembly",
        "key_type",
    ]
    .map(column);
    let mut results = Vec::with_capacity(table.rows);
    for row in 0..table.rows {
        let cell = |index: usize| positions[index].map(|position| table.cell(position, row));
        let text = |index: usize, default: &str| cell(index).map_or(default.to_owned(), python_str);
        results.push(matcher.match_curated_nonlof_variant(
            cell(0).unwrap_or(Cell::Text("")),
            cell(1),
            cell(2),
            [cell(3), cell(4), cell(5), cell(6)],
            &text(7, "hg38"),
            &text(8, "auto"),
        )?);
    }
    if results.is_empty() {
        return Ok(Table {
            names: Vec::new(),
            columns: Vec::new(),
            rows: 0,
        });
    }
    let texts = |pick: &dyn Fn(&CuratedMatch) -> String| {
        Column::Text(results.iter().map(|result| Some(pick(result))).collect())
    };
    let score =
        |index: usize| Column::Int(results.iter().map(|result| result.scores[index]).collect());
    let names = [
        "match_route",
        "lof_score",
        "gof_score",
        "dn_score",
        "mechanism_exclusive",
        "exclusive_mechanisms",
        "gofcards_accession_id",
        "gofcards_variant_id",
        "matches",
    ];
    let columns = vec![
        texts(&|result| result.match_route.clone()),
        score(0),
        score(1),
        score(2),
        Column::Bool(
            results
                .iter()
                .map(CuratedMatch::mechanism_exclusive)
                .collect(),
        ),
        texts(&|result| result.mechanisms.join(";")),
        texts(&|result| result.gofcards_accession_id.clone()),
        texts(&|result| result.gofcards_variant_id.clone()),
        texts(&|result| result.evidence_json()),
    ];
    Ok(Table {
        names: names.iter().map(|name| (*name).to_owned()).collect(),
        columns,
        rows: results.len(),
    })
}

/// Python's `str(value)` of a table value (`nan` for a missing one).
fn python_str(cell: Cell<'_>) -> String {
    match cell {
        Cell::Missing => "nan".to_owned(),
        Cell::Text(text) => text.to_owned(),
        Cell::Int(number) => number.to_string(),
        Cell::Float(number) => float_repr(number),
        Cell::Bool(flag) => (if flag { "True" } else { "False" }).to_owned(),
    }
}

/// `gene_mechanism_common._clean`: the stripped text, or `""` for a missing
/// value and the spellings `nan`, `none` and `<na>` in any case.
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

/// `gene_mechanism_common._bool_value`: True, a non-zero integer, or the texts
/// `1`, `true`, `t`, `yes`, `y` in any case. A float (other than NaN, which is
/// false) goes through its text, which is never one of those spellings.
fn bool_value(cell: Option<Cell<'_>>) -> bool {
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

/// `gene_mechanism_common._split_annotation_tokens`: the upper-cased tokens
/// between `&`, `,`, `;` and `|`, as a set (here a list without repeats).
fn split_annotation_tokens(cell: Option<Cell<'_>>) -> Vec<String> {
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

/// `gene_mechanism_common._bounded_mechanism_score`: `int(float(_clean(value)))`
/// when it is 0, 1 or 2, else 0; text Python cannot read as a float gives 0.
/// An infinite value raises OverflowError in Python, and here.
fn bounded_mechanism_score(cell: Option<Cell<'_>>) -> Result<i64> {
    let overflow = || -> Box<dyn std::error::Error> {
        Box::new(AcmgError::new(
            "OverflowError",
            "cannot convert float infinity to integer",
        ))
    };
    let number = match cell {
        None | Some(Cell::Missing) | Some(Cell::Bool(_)) => return Ok(0),
        Some(Cell::Int(number)) => return Ok(if (0..=2).contains(&number) { number } else { 0 }),
        Some(Cell::Float(number)) => number,
        Some(Cell::Text(text)) => {
            let text = clean_text(text);
            match parse_float(&text) {
                Some(number) => number,
                // Python's float() also reads digits of other scripts: not reproduced.
                None if !text.is_ascii() => {
                    return Err(invalid(format!(
                        "refused (not reproduced): float() of the non-ASCII text {text:?}"
                    )));
                }
                None => return Ok(0),
            }
        }
    };
    if number.is_nan() {
        return Ok(0);
    }
    if number.is_infinite() {
        return Err(overflow());
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

/// The curated mechanisms a score of 2 can establish (`CURATED_EXACT_MECHANISMS`), sorted.
const CURATED_EXACT_MECHANISMS: [&str; 2] = ["DOMINANT_NEGATIVE", "GOF"];

/// `infer_query_variant_effect(row)` for one table row. The hub (T2) calls it
/// for every row.
pub(crate) fn infer_query_variant_effect(row: &EffectInputs<'_>) -> Result<QueryEffect> {
    let mut predicted: Vec<&'static str> = Vec::new();
    let loftee = split_annotation_tokens(row.lof);
    if loftee.iter().any(|token| token == "HC") {
        predicted.push("LOFTEE_HC");
    }
    if loftee.iter().any(|token| token == "OS") {
        predicted.push("LOFTEE_OS");
    }
    let consequence = clean(row.consequence).to_lowercase();
    let nmd = clean(row.nmd).to_lowercase();
    let lof_filter = clean(row.lof_filter).to_uppercase();
    let truncating = consequence.contains("stop_gained") || consequence.contains("frameshift");
    let escapes_nmd = nmd.contains("escaping") || lof_filter.contains("END_TRUNC");
    let splice_frameshift = bool_value(row.splicing_frameshift);
    if (truncating || splice_frameshift) && !escapes_nmd {
        predicted.push("NMD_PREDICTED_LOF");
    } else if truncating || splice_frameshift {
        predicted.push("NMD_ESCAPING_TRUNCATION");
    }
    if bool_value(row.vep_consq_lof) {
        predicted.push("VEP_LOF");
    }
    if bool_value(row.splicing_lof) {
        predicted.push("PRIVA_SPLICE_LOF");
    }
    if bool_value(row.utr5_lof) {
        predicted.push("PRIVA_5UTR_LOF");
    }

    // Scores in the order LOF, GOF, DOMINANT_NEGATIVE.
    let mut scores = [
        bounded_mechanism_score(row.variant_lof_score)?,
        bounded_mechanism_score(row.variant_gof_score)?,
        bounded_mechanism_score(row.variant_dn_score)?,
    ];
    if split_annotation_tokens(row.variant_gof_tag)
        .iter()
        .any(|token| token == "GOF")
    {
        scores[1] = 2;
    }
    // Sorted, as `sorted(curated_mechanisms)` gives them.
    let curated: Vec<&str> = CURATED_EXACT_MECHANISMS
        .into_iter()
        .filter(|mechanism| scores[if *mechanism == "GOF" { 1 } else { 2 }] == 2)
        .collect();
    let canonical = |mechanisms: &[&str]| -> Vec<String> {
        mechanisms
            .iter()
            .map(|mechanism| format!("CANONICAL_EXACT_{mechanism}"))
            .collect()
    };
    let predicted: Vec<String> = predicted.into_iter().map(str::to_owned).collect();
    let (effect, evidence, suppressed, exact): (String, Vec<String>, Vec<String>, Vec<&str>) =
        if predicted.iter().any(|token| token == "NMD_PREDICTED_LOF") {
            scores = [2, 0, 0];
            (
                "exact_known_LOF".to_owned(),
                predicted,
                canonical(&curated),
                vec!["LOF"],
            )
        } else if !curated.is_empty() {
            scores = [
                0,
                if curated.contains(&"GOF") { 2 } else { 0 },
                if curated.contains(&"DOMINANT_NEGATIVE") {
                    2
                } else {
                    0
                },
            ];
            (
                format!("exact_known_{}", curated.join("+")),
                canonical(&curated),
                predicted,
                curated.clone(),
            )
        } else if !predicted.is_empty() {
            let floor = if predicted.iter().any(|token| token == "LOFTEE_HC") {
                2
            } else {
                1
            };
            scores[0] = scores[0].max(floor);
            (
                "predicted_LOF_high_confidence".to_owned(),
                predicted,
                Vec::new(),
                Vec::new(),
            )
        } else {
            scores = scores.map(|score| score.max(1));
            ("uncertain".to_owned(), Vec::new(), Vec::new(), Vec::new())
        };
    Ok(QueryEffect {
        variant_effect: effect,
        variant_effect_evidence: evidence.join(";"),
        variant_effect_suppressed_evidence: suppressed.join(";"),
        variant_lof_score: scores[0],
        variant_gof_score: scores[1],
        variant_dn_score: scores[2],
        variant_mechanism_exclusive: !exact.is_empty(),
        // `exact` is already sorted (LOF alone, or DOMINANT_NEGATIVE before GOF).
        variant_exact_mechanisms: exact.join(";"),
    })
}

// ------------------------------------------------------------------------------------------
// The exact curated-allele matcher (`ExactVariantMatcher`, lines 178-642).

/// One curated allele as `_gofcards_match_entry` reports it: 18 text fields.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct MatchEntry {
    pub(crate) source: String,
    pub(crate) mechanism: String,
    pub(crate) canonical_mechanism_json: String,
    pub(crate) symbol: String,
    pub(crate) gofcards_variant_id: String,
    pub(crate) gofcards_allele_key: String,
    pub(crate) gofcards_accession_id: String,
    pub(crate) assembly: String,
    pub(crate) vep_transcript: String,
    pub(crate) hgvsp: String,
    pub(crate) hgvsc: String,
    pub(crate) chrom: String,
    pub(crate) pos: String,
    pub(crate) ref_allele: String,
    pub(crate) alt_allele: String,
    pub(crate) match_eligibility: String,
    pub(crate) pmids: String,
}

impl MatchEntry {
    /// `json.dumps(entry, sort_keys=True, separators=(",", ":"))` (keys in code-point order).
    fn json(&self) -> String {
        let fields: [(&str, &str); 18] = [
            ("HGVSc", &self.hgvsc),
            ("HGVSp", &self.hgvsp),
            ("alt", &self.alt_allele),
            ("assembly", &self.assembly),
            ("canonical_assertion_source", &self.source),
            ("canonical_mechanism_json", &self.canonical_mechanism_json),
            ("chrom", &self.chrom),
            ("gofcards_accession_id", &self.gofcards_accession_id),
            ("gofcards_allele_key", &self.gofcards_allele_key),
            ("gofcards_variant_id", &self.gofcards_variant_id),
            ("match_eligibility", &self.match_eligibility),
            ("mechanism", &self.mechanism),
            ("pmids", &self.pmids),
            ("pos", &self.pos),
            ("ref", &self.ref_allele),
            ("source", &self.source),
            ("symbol", &self.symbol),
            ("vep_transcript", &self.vep_transcript),
        ];
        let body: Vec<String> = fields
            .iter()
            .map(|(key, value)| format!("{}:{}", quote(key, true), quote(value, true)))
            .collect();
        format!("{{{}}}", body.join(","))
    }
}

/// What `match_curated_nonlof_variant` returns that the step reads.
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct CuratedMatch {
    pub(crate) match_route: String,
    /// Scores of LOF, GOF and DOMINANT_NEGATIVE (2 for a matched mechanism, else 0).
    pub(crate) scores: [i64; 3],
    /// `exclusive_mechanisms`: the matched mechanisms, sorted.
    pub(crate) mechanisms: Vec<String>,
    pub(crate) gofcards_accession_id: String,
    pub(crate) gofcards_variant_id: String,
    pub(crate) matches: Vec<MatchEntry>,
}

impl CuratedMatch {
    pub(crate) fn mechanism_exclusive(&self) -> bool {
        !self.mechanisms.is_empty()
    }

    /// `json.dumps(match["matches"], sort_keys=True, separators=(",", ":"))`.
    pub(crate) fn evidence_json(&self) -> String {
        let items: Vec<String> = self.matches.iter().map(MatchEntry::json).collect();
        format!("[{}]", items.join(","))
    }
}

/// A refusal: a value of a kind the cache never holds, where Python would raise or behave oddly.
fn odd(what: impl std::fmt::Display) -> Box<dyn std::error::Error> {
    invalid(format!("refused (not reproduced): {what}"))
}

/// Python truthiness of a JSON value.
fn truthy(value: &Json) -> bool {
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

/// `value.get(key)` on a dict; any other kind is the AttributeError Python raises.
fn get<'a>(value: &'a Json, key: &str, what: &str) -> Result<Option<&'a Json>> {
    match value {
        Json::Object(items) => Ok(items
            .iter()
            .rev()
            .find(|(name, _)| name == key)
            .map(|(_, v)| v)),
        _ => Err(odd(format!("{what} is not a JSON object"))),
    }
}

/// `(value or {})`: `None` stands for the empty dict.
fn or_empty(value: Option<&Json>) -> Option<&Json> {
    value.filter(|value| truthy(value))
}

/// `x.get(key)` where `x` is `(something or {})`.
fn get_in<'a>(value: Option<&'a Json>, key: &str, what: &str) -> Result<Option<&'a Json>> {
    match value {
        None => Ok(None),
        Some(value) => get(value, key, what),
    }
}

/// `(value or {}).items()` / `(value or []) ` iteration over a dict.
fn dict_items<'a>(value: Option<&'a Json>, what: &str) -> Result<&'a [(String, Json)]> {
    match or_empty(value) {
        None => Ok(&[]),
        Some(Json::Object(items)) => Ok(items),
        Some(_) => Err(odd(format!("{what} is not a JSON object"))),
    }
}

/// Iteration over `(value or [])` where only dict items are used: a list gives
/// its items, a dict or a text gives nothing usable (keys and characters are
/// never dicts), anything else is the TypeError Python raises.
fn list_items<'a>(value: Option<&'a Json>, what: &str) -> Result<&'a [Json]> {
    match or_empty(value) {
        None | Some(Json::Object(_)) | Some(Json::Str(_)) => Ok(&[]),
        Some(Json::Array(items)) => Ok(items),
        Some(_) => Err(odd(format!("{what} is not iterable"))),
    }
}

/// A JSON value as the table value Python's `_clean` and the resolver see.
fn json_cell(value: Option<&Json>) -> Result<Option<Cell<'_>>> {
    Ok(match value {
        None | Some(Json::Null) => None,
        Some(Json::Bool(flag)) => Some(Cell::Bool(*flag)),
        Some(Json::Int(number)) => Some(Cell::Int(
            i64::try_from(*number).map_err(|_| odd(format!("the integer {number}")))?,
        )),
        Some(Json::Float(number)) if number.is_nan() => Some(Cell::Missing),
        Some(Json::Float(number)) => Some(Cell::Float(*number)),
        Some(Json::Str(text)) => Some(Cell::Text(text)),
        Some(_) => return Err(odd("a list or dict where the cache holds text")),
    })
}

/// `_clean(value)` of a JSON value.
fn clean_json(value: Option<&Json>) -> Result<String> {
    Ok(clean(json_cell(value)?))
}

/// `str(Path(text))`: repeated and trailing slashes and `.` parts dropped.
fn python_path_text(text: &str) -> String {
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

/// `_norm_hgvsp`: the protein change after the last `:`, without `p.`, upper-cased.
pub(crate) fn norm_hgvsp(cleaned: &str) -> String {
    if cleaned.is_empty() {
        return String::new();
    }
    let protein = python_strip(cleaned.rsplit(':').next().unwrap_or(""));
    protein.strip_prefix("p.").unwrap_or(protein).to_uppercase()
}

/// `_norm_chrom_key` of a cleaned text.
fn norm_chrom_key(cleaned: &str) -> String {
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
fn norm_pos_key(cleaned: &str) -> Result<String> {
    if cleaned.is_empty() {
        return Ok(String::new());
    }
    match parse_float(cleaned) {
        Some(number) if number.is_finite() && number.fract() == 0.0 => {
            Ok(integral_float_text(number))
        }
        Some(_) => Ok(cleaned.to_owned()),
        None if !cleaned.is_ascii() => Err(odd(format!(
            "float() of the non-ASCII position {cleaned:?}"
        ))),
        None => Ok(cleaned.to_owned()),
    }
}

/// `_genomic_match_key` of cleaned texts.
fn genomic_match_key(chrom: &str, pos: &str, ref_allele: &str, alt_allele: &str) -> Result<String> {
    let chrom_key = norm_chrom_key(chrom);
    let pos_key = norm_pos_key(pos)?;
    if chrom_key.is_empty() || pos_key.is_empty() {
        return Ok(String::new());
    }
    Ok([
        chrom_key,
        pos_key,
        ref_allele.to_uppercase(),
        alt_allele.to_uppercase(),
    ]
    .join("|"))
}

/// `_normalize_assembly` of a cleaned text.
fn normalize_assembly(cleaned: &str) -> String {
    let lower = cleaned.to_lowercase();
    match lower.as_str() {
        "grch37" | "hg19" | "b37" => "hg19".to_owned(),
        "grch38" | "hg38" | "b38" => "hg38".to_owned(),
        _ => lower,
    }
}

/// `_requested_genomic_key_types` of a cleaned text.
fn requested_genomic_key_types(cleaned: &str) -> &'static [&'static str] {
    match cleaned.to_lowercase().as_str() {
        "" | "auto" => &["vcf", "genomic"],
        "vcf" => &["vcf"],
        "genomic" => &["genomic"],
        _ => &[],
    }
}

/// `_normalize_sequence_mechanism` of a cleaned text.
fn normalize_sequence_mechanism(cleaned: &str) -> String {
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

/// The genomic coordinates of one assembly block (`block["genomic"]`).
#[derive(Clone, Debug, Default)]
struct Genomic {
    /// `coords.get("chrom")` is truthy and `coords.get("pos")` is not None.
    indexable: bool,
    chrom: String,
    pos: String,
    ref_allele: String,
    alt_allele: String,
}

#[derive(Clone, Debug, Default)]
struct TranscriptView {
    name: String,
    /// `by_hgvsp`: (protein change, its first coding change or "").
    by_hgvsp: Vec<(String, String)>,
    /// `by_hgvsc`: (coding change, `detail.get("hgvsp") or ""`).
    by_hgvsc: Vec<(String, String)>,
}

#[derive(Clone, Debug, Default)]
struct AssemblyBlock {
    name: String,
    genomic: Genomic,
    transcripts: Vec<TranscriptView>,
}

/// One runtime-eligible curated allele (`_load_canonical_exact_nonlof_rows`).
#[derive(Clone, Debug, Default)]
struct CanonicalVariant {
    symbol: String,
    source: String,
    mechanism: String,
    variant_id: String,
    allele_key: String,
    accession: String,
    eligibility: String,
    pmids: String,
    assemblies: Vec<AssemblyBlock>,
}

impl CanonicalVariant {
    /// `_gofcards_match_entry(variant, assembly=..., transcript=..., hgvsp=..., hgvsc=...)`.
    fn entry(
        &self,
        block: &AssemblyBlock,
        transcript: &str,
        hgvsp: &str,
        hgvsc: &str,
    ) -> MatchEntry {
        MatchEntry {
            source: if self.source.is_empty() {
                "GoFCards".to_owned()
            } else {
                self.source.clone()
            },
            mechanism: if self.mechanism.is_empty() {
                "GOF".to_owned()
            } else {
                self.mechanism.clone()
            },
            canonical_mechanism_json: String::new(),
            symbol: self.symbol.clone(),
            gofcards_variant_id: self.variant_id.clone(),
            gofcards_allele_key: self.allele_key.clone(),
            gofcards_accession_id: self.accession.clone(),
            assembly: block.name.clone(),
            vep_transcript: transcript.to_owned(),
            hgvsp: hgvsp.to_owned(),
            hgvsc: hgvsc.to_owned(),
            chrom: block.genomic.chrom.clone(),
            pos: block.genomic.pos.clone(),
            ref_allele: block.genomic.ref_allele.clone(),
            alt_allele: block.genomic.alt_allele.clone(),
            match_eligibility: self.eligibility.clone(),
            pmids: self.pmids.clone(),
        }
    }
}

/// `_deduplicate_gofcards_matches`: sorted by (variant, assembly, transcript,
/// protein change), then the first of each (mechanism, source, variant).
fn deduplicate(mut rows: Vec<MatchEntry>) -> Vec<MatchEntry> {
    rows.sort_by(|a, b| {
        (
            &a.gofcards_variant_id,
            &a.assembly,
            &a.vep_transcript,
            &a.hgvsp,
        )
            .cmp(&(
                &b.gofcards_variant_id,
                &b.assembly,
                &b.vep_transcript,
                &b.hgvsp,
            ))
    });
    let mut seen = HashSet::new();
    rows.into_iter()
        .filter(|row| {
            seen.insert((
                row.mechanism.clone(),
                row.source.clone(),
                row.gofcards_variant_id.clone(),
            ))
        })
        .collect()
}

/// Reads the non-LOF mechanism JSON as `open_text` and `json.load` do.
pub(crate) fn read_mechanism_json(path: &str) -> Result<Json> {
    let mut raw = Vec::new();
    std::fs::File::open(path)
        .map_err(|error| odd(format!("reading the mechanism JSON {path}: {error}")))?
        .read_to_end(&mut raw)?;
    if path.ends_with(".gz") {
        let mut plain = Vec::new();
        flate2::read::MultiGzDecoder::new(&raw[..])
            .read_to_end(&mut plain)
            .map_err(|error| odd(format!("decompressing {path}: {error}")))?;
        raw = plain;
    }
    let text = String::from_utf8(raw).map_err(|_| odd(format!("{path} is not UTF-8")))?;
    if text.starts_with('\u{feff}') {
        return Err(odd(format!("{path} starts with a byte-order mark")));
    }
    json_text::parse(&text, false).map_err(|error| odd(format!("{path}: {error}")))
}

/// `GeneMechanismHub._load_mechanisms`: the cache's genes by symbol, each gene
/// also under its resolved symbol when that is new, in Python's dict order
/// (a repeated symbol keeps its first place and takes the later gene).
fn mechanisms_by_symbol<'a>(raw: &'a Json, resolver: &dyn SymbolResolver) -> Result<Vec<&'a Json>> {
    let Json::Object(genes) = raw else {
        return Err(odd("the mechanism JSON is not an object"));
    };
    let mut order: Vec<(String, &'a Json)> = Vec::new();
    let mut place: HashMap<String, usize> = HashMap::new();
    for (_, info) in genes {
        let symbol = clean_json(get(info, "symbol", "a gene record")?)?;
        if symbol.is_empty() {
            continue;
        }
        match place.get(&symbol) {
            Some(index) => order[*index].1 = info,
            None => {
                place.insert(symbol.clone(), order.len());
                order.push((symbol.clone(), info));
            }
        }
        let resolved = resolver.try_resolve_symbol(Cell::Text(&symbol));
        if !resolved.is_empty() && !place.contains_key(&resolved) {
            place.insert(resolved.clone(), order.len());
            order.push((resolved, info));
        }
    }
    Ok(order.into_iter().map(|(_, info)| info).collect())
}

/// One assembly block of a curated variant.
fn assembly_block(name: &str, block: &Json) -> Result<AssemblyBlock> {
    let coords = or_empty(get(block, "genomic", "an assembly block")?);
    let field = |key: &str| get_in(coords, key, "genomic coordinates");
    let chrom = field("chrom")?;
    let pos = field("pos")?;
    let genomic = Genomic {
        indexable: chrom.is_some_and(truthy) && !matches!(pos, None | Some(Json::Null)),
        chrom: clean_json(chrom)?,
        pos: clean_json(pos)?,
        ref_allele: clean_json(field("ref")?)?,
        alt_allele: clean_json(field("alt")?)?,
    };
    let mut transcripts = Vec::new();
    for (transcript, view) in dict_items(
        get(block, "transcripts", "an assembly block")?,
        "transcripts",
    )? {
        let mut by_hgvsp = Vec::new();
        for (hgvsp, coding) in dict_items(get(view, "by_hgvsp", "a transcript view")?, "by_hgvsp")?
        {
            let first = match coding {
                coding if !truthy(coding) => String::new(),
                Json::Array(items) => match &items[0] {
                    Json::Str(text) => text.clone(),
                    _ => return Err(odd("a coding change that is not text")),
                },
                Json::Str(text) => text.chars().next().map(String::from).unwrap_or_default(),
                _ => return Err(odd("by_hgvsp holds neither a list nor a text")),
            };
            by_hgvsp.push((hgvsp.clone(), first));
        }
        let mut by_hgvsc = Vec::new();
        for (hgvsc, detail) in dict_items(get(view, "by_hgvsc", "a transcript view")?, "by_hgvsc")?
        {
            let protein = match or_empty(get(detail, "hgvsp", "a by_hgvsc entry")?) {
                None => String::new(),
                Some(Json::Str(text)) => text.clone(),
                Some(_) => return Err(odd("a by_hgvsc protein change that is not text")),
            };
            by_hgvsc.push((hgvsc.clone(), protein));
        }
        transcripts.push(TranscriptView {
            name: transcript.clone(),
            by_hgvsp,
            by_hgvsc,
        });
    }
    Ok(AssemblyBlock {
        name: name.to_owned(),
        genomic,
        transcripts,
    })
}

/// `ExactVariantMatcher._load_canonical_exact_nonlof_rows`.
fn canonical_variants(
    genes: &[&Json],
    resolver: &dyn SymbolResolver,
) -> Result<Vec<CanonicalVariant>> {
    let mut entries = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();
    for gene in genes {
        let symbol = resolver.try_resolve_symbol(
            json_cell(get(gene, "symbol", "a gene record")?)?.unwrap_or(Cell::Missing),
        );
        if symbol.is_empty() {
            continue;
        }
        for wrapped in list_items(
            get(gene, "variant_level", "a gene record")?,
            "variant_level",
        )? {
            let Json::Object(sources) = wrapped else {
                continue;
            };
            for (source, assertion) in sources {
                if source == "ClinVar_VCV" || !matches!(assertion, Json::Object(_)) {
                    continue;
                }
                if let Some(status) = get(assertion, "exact_normalization_status", "an assertion")?
                {
                    if truthy(status) && status.as_str() != Some("matched_gene_concordant") {
                        continue;
                    }
                }
                let exact = get(assertion, "exact_normalized_variants", "an assertion")?;
                for variant in list_items(exact, "exact_normalized_variants")? {
                    if !matches!(variant, Json::Object(_)) {
                        continue;
                    }
                    let record = or_empty(get(variant, "record", "a variant")?);
                    let eligibility = or_empty(get_in(record, "eligibility", "a record")?);
                    let status = get_in(eligibility, "status", "an eligibility block")?;
                    if status.and_then(Json::as_str) != Some("eligible") {
                        continue;
                    }
                    let variant_id = clean_json(get(variant, "variant_id", "a variant")?)?;
                    let marker = format!("{source}|{symbol}|{variant_id}");
                    if variant_id.is_empty() || seen.contains(&marker) {
                        continue;
                    }
                    seen.insert(marker);
                    let record_source = or_empty(get_in(record, "source", "a record")?);
                    let annotations = or_empty(get_in(record, "annotations", "a record")?);
                    let mut pmids = Vec::new();
                    for item in list_items(get_in(record, "evidence", "a record")?, "evidence")? {
                        // `e.get("pmid", "")` on a non-dict item is Python's AttributeError.
                        match get(item, "pmid", "an evidence item")? {
                            Some(pmid) if truthy(pmid) => match pmid {
                                Json::Str(text) => pmids.push(text.clone()),
                                _ => return Err(odd("a PMID that is not text")),
                            },
                            _ => {}
                        }
                    }
                    let mut assemblies = Vec::new();
                    for (name, block) in
                        dict_items(get(variant, "assemblies", "a variant")?, "assemblies")?
                    {
                        assemblies.push(assembly_block(name, block)?);
                    }
                    entries.push(CanonicalVariant {
                        symbol: symbol.clone(),
                        source: clean_text(source),
                        mechanism: clean_json(get(assertion, "mechanism", "an assertion")?)?,
                        variant_id,
                        allele_key: clean_json(get_in(
                            record_source,
                            "gofcards_allele_key",
                            "a record source",
                        )?)?,
                        accession: clean_json(get_in(
                            annotations,
                            "clinvar_variation_id",
                            "annotations",
                        )?)?,
                        eligibility: clean_json(status)?,
                        pmids: pmids.join(";"),
                        assemblies,
                    });
                }
            }
        }
    }
    Ok(entries)
}

type Index = HashMap<Vec<String>, Vec<MatchEntry>>;

/// `ExactVariantMatcher`: the curated alleles of the non-LOF mechanism JSON,
/// indexed by protein change, coding change and genomic allele. The JSON is
/// read and each index built on first use, as in Python.
pub(crate) struct ExactVariantMatcher<'r> {
    mechanism_json: String,
    resolver: &'r dyn SymbolResolver,
    variants: Option<Vec<CanonicalVariant>>,
    by_hgvsp: Option<Index>,
    by_hgvsc: Option<Index>,
    by_genomic: Option<Index>,
}

impl<'r> ExactVariantMatcher<'r> {
    pub(crate) fn new(mechanism_json: &str, resolver: &'r dyn SymbolResolver) -> Self {
        ExactVariantMatcher {
            mechanism_json: mechanism_json.to_owned(),
            resolver,
            variants: None,
            by_hgvsp: None,
            by_hgvsc: None,
            by_genomic: None,
        }
    }

    fn variants(&mut self) -> Result<&[CanonicalVariant]> {
        if self.variants.is_none() {
            let raw = read_mechanism_json(&self.mechanism_json)?;
            let genes = mechanisms_by_symbol(&raw, self.resolver)?;
            self.variants = Some(canonical_variants(&genes, self.resolver)?);
        }
        Ok(self.variants.as_deref().unwrap_or_default())
    }

    /// The text `canonical_mechanism_json` carries in every entry.
    fn path_text(&self) -> String {
        python_path_text(&self.mechanism_json)
    }

    fn finish(&self, index: HashMap<Vec<String>, Vec<MatchEntry>>) -> Index {
        let path = self.path_text();
        index
            .into_iter()
            .map(|(key, rows)| {
                let rows = rows
                    .into_iter()
                    .map(|mut row| {
                        row.canonical_mechanism_json = path.clone();
                        row
                    })
                    .collect();
                (key, deduplicate(rows))
            })
            .collect()
    }

    fn hgvsp_index(&mut self) -> Result<&Index> {
        if self.by_hgvsp.is_none() {
            let mut index: Index = HashMap::new();
            for variant in self.variants()? {
                for block in &variant.assemblies {
                    for view in &block.transcripts {
                        for (hgvsp, coding) in &view.by_hgvsp {
                            let key = norm_hgvsp(&clean_text(hgvsp));
                            if key.is_empty() {
                                continue;
                            }
                            index
                                .entry(vec![variant.symbol.clone(), key])
                                .or_default()
                                .push(variant.entry(block, &view.name, hgvsp, coding));
                        }
                    }
                }
            }
            self.by_hgvsp = Some(self.finish(index));
        }
        Ok(self.by_hgvsp.as_ref().expect("built above"))
    }

    fn hgvsc_index(&mut self) -> Result<&Index> {
        if self.by_hgvsc.is_none() {
            let mut index: Index = HashMap::new();
            for variant in self.variants()? {
                for block in &variant.assemblies {
                    for view in &block.transcripts {
                        for (hgvsc, protein) in &view.by_hgvsc {
                            let key = clean_text(hgvsc).to_uppercase();
                            if key.is_empty() {
                                continue;
                            }
                            index
                                .entry(vec![variant.symbol.clone(), key])
                                .or_default()
                                .push(variant.entry(block, &view.name, protein, hgvsc));
                        }
                    }
                }
            }
            self.by_hgvsc = Some(self.finish(index));
        }
        Ok(self.by_hgvsc.as_ref().expect("built above"))
    }

    /// Keyed by (symbol, assembly, genomic key); Python files every entry under
    /// both key types `vcf` and `genomic`, so one list serves both.
    fn genomic_index(&mut self) -> Result<&Index> {
        if self.by_genomic.is_none() {
            let mut index: Index = HashMap::new();
            for variant in self.variants()? {
                for block in &variant.assemblies {
                    let coords = &block.genomic;
                    if !coords.indexable {
                        continue;
                    }
                    let key = genomic_match_key(
                        &coords.chrom,
                        &coords.pos,
                        &coords.ref_allele,
                        &coords.alt_allele,
                    )?;
                    if key.is_empty() {
                        continue;
                    }
                    index
                        .entry(vec![variant.symbol.clone(), block.name.clone(), key])
                        .or_default()
                        .push(variant.entry(block, "", "", ""));
                }
            }
            self.by_genomic = Some(self.finish(index));
        }
        Ok(self.by_genomic.as_ref().expect("built above"))
    }

    /// `match_curated_nonlof_variant(gene_symbol, hgvsp=, hgvsc=, chrom=, pos=, ref=, alt=,
    /// assembly=, key_type=)`. `None` stands for the default `""`.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn match_curated_nonlof_variant(
        &mut self,
        gene_symbol: Cell<'_>,
        hgvsp: Option<Cell<'_>>,
        hgvsc: Option<Cell<'_>>,
        genomic: [Option<Cell<'_>>; 4],
        assembly: &str,
        key_type: &str,
    ) -> Result<CuratedMatch> {
        let symbol = self.resolver.resolve_symbol(gene_symbol)?;
        let hgvsp_key = norm_hgvsp(&clean(hgvsp));
        let hgvsc_key = clean(hgvsc).rsplit(':').next().unwrap_or("").to_uppercase();
        let assembly_key = normalize_assembly(&clean_text(assembly));
        let [chrom, pos, ref_allele, alt_allele] = genomic;
        let genomic_key = genomic_match_key(
            &clean(chrom),
            &clean(pos),
            &clean(ref_allele),
            &clean(alt_allele),
        )?;
        let mut matches: Vec<MatchEntry> = Vec::new();
        let mut route = "";
        if !symbol.is_empty() && !hgvsp_key.is_empty() {
            let key = vec![symbol.clone(), hgvsp_key.clone()];
            matches = self.hgvsp_index()?.get(&key).cloned().unwrap_or_default();
            if !matches.is_empty() {
                route = "hgvsp";
            }
        }
        if !symbol.is_empty() && matches.is_empty() && !hgvsc_key.is_empty() {
            let key = vec![symbol.clone(), hgvsc_key.clone()];
            matches = self.hgvsc_index()?.get(&key).cloned().unwrap_or_default();
            if !matches.is_empty() {
                route = "hgvsc";
            }
        }
        if !symbol.is_empty()
            && matches.is_empty()
            && !genomic_key.is_empty()
            && !assembly_key.is_empty()
        {
            let key = vec![symbol.clone(), assembly_key.clone(), genomic_key.clone()];
            let found = self.genomic_index()?.get(&key).cloned().unwrap_or_default();
            for _ in requested_genomic_key_types(&clean_text(key_type)) {
                if !found.is_empty() {
                    matches.extend(found.iter().cloned());
                    route = "genomic";
                }
            }
            matches = deduplicate(matches);
        }
        let normalized: Vec<String> = matches
            .iter()
            .map(|entry| normalize_sequence_mechanism(&clean_text(&entry.mechanism)))
            .collect();
        let keep: Vec<bool> = normalized
            .iter()
            .map(|mechanism| ["GOF", "DOMINANT_NEGATIVE", "LOF"].contains(&mechanism.as_str()))
            .collect();
        let mut kept = Vec::new();
        let mut mechanisms: Vec<String> = Vec::new();
        for ((entry, mechanism), keep) in matches.into_iter().zip(normalized).zip(keep) {
            if keep {
                if !mechanisms.contains(&mechanism) {
                    mechanisms.push(mechanism);
                }
                kept.push(entry);
            }
        }
        mechanisms.sort();
        let mut scores = [0; 3];
        for mechanism in &mechanisms {
            match mechanism.as_str() {
                "LOF" => scores[0] = 2,
                "GOF" => scores[1] = 2,
                _ => scores[2] = 2,
            }
        }
        let sorted_ids = |pick: fn(&MatchEntry) -> &String| {
            let mut ids: Vec<String> = kept
                .iter()
                .map(|entry| clean_text(pick(entry)))
                .filter(|id| !id.is_empty())
                .collect();
            ids.sort();
            ids.dedup();
            ids.join(";")
        };
        Ok(CuratedMatch {
            match_route: route.to_owned(),
            scores,
            gofcards_accession_id: sorted_ids(|entry| &entry.gofcards_accession_id),
            gofcards_variant_id: sorted_ids(|entry| &entry.gofcards_variant_id),
            mechanisms,
            matches: kept,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn effect(inputs: EffectInputs<'_>) -> QueryEffect {
        infer_query_variant_effect(&inputs).unwrap()
    }

    #[test]
    fn decay_lof_outranks_a_curated_allele() {
        let call = effect(EffectInputs {
            consequence: Some(Cell::Text("Stop_Gained&splice_region_variant")),
            variant_gof_score: Some(Cell::Int(2)),
            variant_dn_score: Some(Cell::Text(" 2.7 ")),
            lof: Some(Cell::Text("HC")),
            ..EffectInputs::default()
        });
        assert_eq!(call.variant_effect, "exact_known_LOF");
        assert_eq!(call.variant_effect_evidence, "LOFTEE_HC;NMD_PREDICTED_LOF");
        assert_eq!(
            call.variant_effect_suppressed_evidence,
            "CANONICAL_EXACT_DOMINANT_NEGATIVE;CANONICAL_EXACT_GOF"
        );
        assert_eq!(
            (
                call.variant_lof_score,
                call.variant_gof_score,
                call.variant_dn_score
            ),
            (2, 0, 0)
        );
        assert!(call.variant_mechanism_exclusive);
        assert_eq!(call.variant_exact_mechanisms, "LOF");
    }

    #[test]
    fn curated_alleles_and_the_gof_tag() {
        let call = effect(EffectInputs {
            variant_gof_tag: Some(Cell::Text("x;gof")),
            vep_consq_lof: Some(Cell::Bool(true)),
            consequence: Some(Cell::Text("frameshift_variant")),
            nmd: Some(Cell::Text("NMD_escaping_variant")),
            ..EffectInputs::default()
        });
        assert_eq!(call.variant_effect, "exact_known_GOF");
        assert_eq!(call.variant_effect_evidence, "CANONICAL_EXACT_GOF");
        assert_eq!(
            call.variant_effect_suppressed_evidence,
            "NMD_ESCAPING_TRUNCATION;VEP_LOF"
        );
        assert_eq!(
            (
                call.variant_lof_score,
                call.variant_gof_score,
                call.variant_dn_score
            ),
            (0, 2, 0)
        );
        assert_eq!(call.variant_exact_mechanisms, "GOF");
    }

    #[test]
    fn predicted_and_uncertain_calls() {
        let call = effect(EffectInputs {
            splicing_lof: Some(Cell::Text(" Yes ")),
            variant_lof_score: Some(Cell::Float(1.0)),
            ..EffectInputs::default()
        });
        assert_eq!(call.variant_effect, "predicted_LOF_high_confidence");
        assert_eq!(call.variant_lof_score, 1);
        assert!(!call.variant_mechanism_exclusive);
        // A float 1.0 is not a true flag in _bool_value; NaN gives 0.
        let call = effect(EffectInputs {
            utr5_lof: Some(Cell::Float(1.0)),
            variant_gof_score: Some(Cell::Missing),
            variant_dn_score: Some(Cell::Text("nan")),
            ..EffectInputs::default()
        });
        assert_eq!(call.variant_effect, "uncertain");
        assert_eq!(
            (
                call.variant_lof_score,
                call.variant_gof_score,
                call.variant_dn_score
            ),
            (1, 1, 1)
        );
        assert_eq!(call.variant_effect_evidence, "");
    }

    #[test]
    fn scores_follow_int_of_float() {
        assert_eq!(bounded_mechanism_score(Some(Cell::Text("1_0"))).unwrap(), 0);
        assert_eq!(bounded_mechanism_score(Some(Cell::Text("2.9"))).unwrap(), 2);
        assert_eq!(
            bounded_mechanism_score(Some(Cell::Text("-0.5"))).unwrap(),
            0
        );
        assert_eq!(bounded_mechanism_score(Some(Cell::Text("abc"))).unwrap(), 0);
        assert_eq!(
            bounded_mechanism_score(Some(Cell::Text("-nan"))).unwrap(),
            0
        );
        assert_eq!(bounded_mechanism_score(Some(Cell::Int(3))).unwrap(), 0);
        assert_eq!(bounded_mechanism_score(Some(Cell::Bool(true))).unwrap(), 0);
        assert_eq!(
            bounded_mechanism_score(Some(Cell::Text("Infinity")))
                .unwrap_err()
                .to_string(),
            "OverflowError: cannot convert float infinity to integer"
        );
        assert!(bounded_mechanism_score(Some(Cell::Text("\u{ff12}"))).is_err());
    }

    #[test]
    fn matcher_text_rules() {
        assert_eq!(norm_hgvsp("ENSP0001.1:p.Val617Phe"), "VAL617PHE");
        assert_eq!(norm_hgvsp(""), "");
        assert_eq!(norm_chrom_key("Chr9"), "9");
        assert_eq!(norm_chrom_key("chrM"), "MT");
        assert_eq!(norm_pos_key("5073770.0").unwrap(), "5073770");
        assert_eq!(norm_pos_key("1e20").unwrap(), "100000000000000000000");
        assert_eq!(norm_pos_key("12a").unwrap(), "12a");
        assert_eq!(norm_pos_key("inf").unwrap(), "inf");
        assert_eq!(normalize_sequence_mechanism("gain of function"), "GOF");
        assert_eq!(
            normalize_sequence_mechanism("Dominant-Negative"),
            "DOMINANT_NEGATIVE"
        );
        assert_eq!(normalize_sequence_mechanism("dn"), "DOMINANT_NEGATIVE");
        assert_eq!(normalize_assembly("GRCh38"), "hg38");
        assert_eq!(requested_genomic_key_types(""), &["vcf", "genomic"]);
        assert!(requested_genomic_key_types("x").is_empty());
        assert_eq!(python_path_text("a//b/./c/"), "a/b/c");
        assert_eq!(python_path_text("//a"), "//a");
        assert_eq!(python_path_text("///a"), "/a");
        assert_eq!(python_path_text(""), ".");
    }

    struct Same;

    impl SymbolResolver for Same {
        fn resolve_symbol(&self, gene_symbol: Cell<'_>) -> Result<String> {
            Ok(clean(Some(gene_symbol)).to_uppercase())
        }

        fn try_resolve_symbol(&self, gene_symbol: Cell<'_>) -> String {
            clean(Some(gene_symbol)).to_uppercase()
        }
    }

    #[test]
    fn matcher_routes_and_evidence() {
        let path = std::env::temp_dir().join(format!("t3_matcher_{}.json", std::process::id()));
        std::fs::write(
            &path,
            r#"{"HGNC:1": {"symbol": "jak2", "variant_level": [{"GoFCards": {
                "mechanism": "gain of function",
                "exact_normalized_variants": [{"variant_id": "GOF1",
                  "record": {"eligibility": {"status": "eligible"}, "evidence": [{"pmid": "7"}, {}]},
                  "assemblies": {"hg19": {"genomic": {"chrom": "chr9", "pos": 5073770, "ref": "g", "alt": "T"},
                    "transcripts": {"ENST1": {"by_hgvsp": {"p.Val617Phe": ["c.1849G>T"]},
                                              "by_hgvsc": {"c.1849G>T": {"hgvsp": "p.Val617Phe"}}}}}}}]},
                "ClinVar_VCV": {"mechanism": "GOF"}}]}}"#,
        )
        .unwrap();
        let resolver = Same;
        let mut matcher = ExactVariantMatcher::new(path.to_str().unwrap(), &resolver);
        let by_protein = matcher
            .match_curated_nonlof_variant(
                Cell::Text("JAK2"),
                Some(Cell::Text("ENSP:p.val617phe")),
                None,
                [None; 4],
                "hg38",
                "auto",
            )
            .unwrap();
        assert_eq!(by_protein.match_route, "hgvsp");
        assert_eq!(by_protein.scores, [0, 2, 0]);
        assert_eq!(by_protein.mechanisms, vec!["GOF".to_owned()]);
        assert!(
            by_protein
                .evidence_json()
                .contains(r#""HGVSc":"c.1849G>T","HGVSp":"p.Val617Phe""#)
        );
        assert!(by_protein.evidence_json().contains(r#""pmids":"7""#));
        let by_allele = matcher
            .match_curated_nonlof_variant(
                Cell::Text("jak2"),
                None,
                None,
                [
                    Some(Cell::Text("9")),
                    Some(Cell::Int(5073770)),
                    Some(Cell::Text("G")),
                    Some(Cell::Text("t")),
                ],
                "GRCh37",
                "auto",
            )
            .unwrap();
        assert_eq!(by_allele.match_route, "genomic");
        assert_eq!(by_allele.gofcards_variant_id, "GOF1");
        let none = matcher
            .match_curated_nonlof_variant(
                Cell::Text("JAK2"),
                Some(Cell::Text("p.Val617Asp")),
                None,
                [None; 4],
                "hg38",
                "auto",
            )
            .unwrap();
        assert!(!none.mechanism_exclusive());
        std::fs::remove_file(&path).unwrap();
    }
}
