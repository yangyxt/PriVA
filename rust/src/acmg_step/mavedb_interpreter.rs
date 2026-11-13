//! `scripts/mavedb_interpreter.py` (842 lines): the assay-family rules that turn a MaveDB score
//! into loss, gain or neutral.
//!
//! `pick_interpreter_from_metadata` looks the URN up in the metadata table (the first row with
//! that `urn`), reads the quality flags and the neutral band (`parse_neutral_band_from_metadata`),
//! builds the assay description (`AssaySpec`) and runs the family's interpreter
//! (`interpret_score`). The PS3/BS3 code reads only the call (`MaveDB_PS3` for LoF or GoF,
//! `MaveDB_BS3` for Neutral); the standardized score and the reason are kept for the comparison
//! with Python. The free-text `notes` and the metadata echoes of the Python dict are not built.
//!
//! Python's `math.log(x, base)` is `log(x) / log(base)`, reproduced as `ln(x) / ln(base)`;
//! `math.log` of zero raises `ValueError: math domain error` (a Kd of infinity), which stops
//! the step in Python and here.
use super::acmg_pp3_bp4_bp7_insilico::harness_options;
use super::acmg_pp3_bp4_bp7_insilico::pandas_masks::python_str;
use super::types::{AcmgError, HarnessCall, not_converted};
use crate::cli::{Result, invalid};
use crate::json_text::Json;
use crate::pandas_table::{Cell, Column, Table};
use crate::python_text::{parse_float, python_strip};
use std::collections::HashMap;

/// `FAMILY_TO_SUBTYPE_ENUM`: the 13 assay families with an interpreter and the member names
/// of their subtype enums.
const FAMILY_SUBTYPES: [(&str, &[&str]); 13] = [
    (
        "Stability_K50",
        &["LOG10_K50_STANDARDIZED", "LOG10_K50_RELATIVE", "K50_RAW"],
    ),
    ("Stability_DG", &["DDG_KCAL_MOL", "DG_KCAL_MOL"]),
    (
        "Abundance",
        &["VAMP_SEQ", "CELL_SURFACE_TRAFF", "TOTAL_EXPRESSION"],
    ),
    (
        "Activity",
        &["ENZYME_ACTIVITY_REL", "REPORTER_ACTIVITY_REL"],
    ),
    ("Binding", &["KD", "PKD", "KA", "KD_APPARENT"]),
    ("Fitness", &["LOG2_ENRICHMENT"]),
    ("Splicing", &["DELTA_PSI_PERCENT", "DELTA_PSI_FRACTION"]),
    ("Electrophysiology", &["V1_2", "CURRENT_DENSITY"]),
    ("Complementation", &["YEAST", "MAMMALIAN_CELL"]),
    ("Localization", &["CELLULAR_LOCALIZATION"]),
    ("Aggregation", &["AGGREGATION_METRIC"]),
    ("Interaction", &["TWO_HYBRID"]),
    ("MPRA", &["REPORTER_ACTIVITY"]),
];

#[derive(Clone, Copy, Debug, PartialEq)]
enum Directionality {
    HigherIsFunction,
    LowerIsFunction,
    NonMonotonic,
}

impl Directionality {
    fn from_name(name: &str) -> Option<Directionality> {
        match name {
            "HIGHER_IS_FUNCTION" => Some(Directionality::HigherIsFunction),
            "LOWER_IS_FUNCTION" => Some(Directionality::LowerIsFunction),
            "NON_MONOTONIC" => Some(Directionality::NonMonotonic),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Scale {
    Centered0,
    WtEq1,
    Raw,
    Percent,
    Fraction,
}

impl Scale {
    fn from_name(name: &str) -> Option<Scale> {
        match name {
            "CENTERED_0" => Some(Scale::Centered0),
            "WT_EQ_1" => Some(Scale::WtEq1),
            "RAW" => Some(Scale::Raw),
            "PERCENT" => Some(Scale::Percent),
            "FRACTION" => Some(Scale::Fraction),
            _ => None,
        }
    }
}

/// `AssaySpec` without the descriptive fields (baseline, units), which only reach `notes`.
#[derive(Clone, Debug)]
struct AssaySpec {
    family: &'static str,
    subtype: &'static str,
    directionality: Directionality,
    scale: Scale,
    neutral_band: (f64, f64),
    flexibility_requirement: bool,
}

/// What `pick_interpreter_from_metadata` returns, as far as the step uses it.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Interpretation {
    pub(crate) call: &'static str,
    pub(crate) standardized_score: Option<f64>,
    pub(crate) reason: String,
    pub(crate) mavedb_ps3: bool,
    pub(crate) mavedb_bs3: bool,
}

fn not_available(reason: impl Into<String>) -> Interpretation {
    Interpretation {
        call: "NA",
        standardized_score: None,
        reason: reason.into(),
        mavedb_ps3: false,
        mavedb_bs3: false,
    }
}

/// `math.log(x, base)` for a float `x`: `log(x) / log(base)`; zero raises, a negative
/// number gives NaN in C and raises in Python (callers never pass one).
fn python_log(x: f64, base: f64) -> Result<f64> {
    if x == 0.0 || x < 0.0 {
        return Err(Box::new(AcmgError::new("ValueError", "math domain error")));
    }
    Ok(x.ln() / base.ln())
}

/// `_to_log2(x)`: `None` for `x <= 0`, else `math.log(x, 2)` (NaN stays NaN).
fn to_log2(x: f64) -> Option<f64> {
    if x <= 0.0 {
        return None;
    }
    python_log(x, 2.0).ok()
}

/// `BaseInterpreter.standardize`.
fn standardize(score: f64, spec: &AssaySpec) -> Option<f64> {
    match spec.scale {
        Scale::Centered0 | Scale::Percent | Scale::Fraction => Some(score),
        Scale::WtEq1 => to_log2(score),
        Scale::Raw => None,
    }
}

/// `BaseInterpreter.call_from_band` (comparisons with NaN are false, as in Python).
fn call_from_band(s: f64, spec: &AssaySpec) -> &'static str {
    let (lo, hi) = spec.neutral_band;
    if lo <= s && s <= hi {
        return "Neutral";
    }
    match spec.directionality {
        Directionality::HigherIsFunction => {
            if s > hi {
                "GoF"
            } else {
                "LoF"
            }
        }
        Directionality::LowerIsFunction => {
            if spec.flexibility_requirement {
                if s < lo { "GoF" } else { "LoF" }
            } else if s < lo {
                "LoF"
            } else {
                "GoF"
            }
        }
        Directionality::NonMonotonic => "Unknown",
    }
}

fn banded(s: f64, spec: &AssaySpec, reason: &str) -> Interpretation {
    Interpretation {
        call: call_from_band(s, spec),
        standardized_score: Some(s),
        reason: reason.to_owned(),
        mavedb_ps3: false,
        mavedb_bs3: false,
    }
}

/// `interpret_score(score, spec)`: the family's interpreter (Electrophysiology gets no
/// gene-specific rules, so it never calls).
fn interpret_score(score: f64, spec: &AssaySpec) -> Result<Interpretation> {
    // validate_spec: the subtype always belongs to the family here; the band must be ordered.
    let (lo, hi) = spec.neutral_band;
    if !(lo < hi) {
        return Ok(not_available("Invalid spec"));
    }
    let simple = |name: &str| -> Interpretation {
        match standardize(score, spec) {
            None => not_available(format!("{name} not standardized")),
            Some(s) => banded(s, spec, &format!("{name} band")),
        }
    };
    Ok(match spec.family {
        "Stability_K50" => {
            if matches!(
                spec.subtype,
                "LOG10_K50_STANDARDIZED" | "LOG10_K50_RELATIVE"
            ) {
                match standardize(score, spec) {
                    None => not_available("K50 not standardized"),
                    Some(s) => banded(s, spec, "K50 heuristic band"),
                }
            } else {
                not_available("K50 RAW requires transform")
            }
        }
        "Stability_DG" => match standardize(score, spec) {
            Some(s) => banded(s, spec, "ΔG/ΔΔG standardized"),
            None => {
                let s = if spec.subtype == "DDG_KCAL_MOL" {
                    -score
                } else {
                    score
                };
                banded(s, spec, "ΔG/ΔΔG transform")
            }
        },
        "Binding" => {
            let mut spec = spec.clone();
            let s = match standardize(score, &spec) {
                Some(s) => s,
                None => {
                    let s = if matches!(spec.subtype, "KD" | "KD_APPARENT") {
                        if score <= 0.0 {
                            return Ok(not_available("Kd <= 0 invalid"));
                        }
                        python_log(1.0 / score, 10.0)?
                    } else {
                        score
                    };
                    if spec.directionality == Directionality::NonMonotonic {
                        spec.directionality = Directionality::HigherIsFunction;
                    }
                    s
                }
            };
            if spec.directionality == Directionality::NonMonotonic {
                Interpretation {
                    standardized_score: Some(s),
                    ..not_available("Binding directionality unknown")
                }
            } else {
                banded(s, &spec, "Binding band")
            }
        }
        "Electrophysiology" => not_available("Gene-specific rule required"),
        "Abundance" => simple("Abundance"),
        "Activity" => simple("Activity"),
        "Fitness" => simple("Fitness"),
        "Splicing" => simple("Splicing"),
        "Complementation" => simple("Complementation"),
        "Localization" => simple("Localization"),
        "Aggregation" => simple("Aggregation"),
        "Interaction" => simple("Interaction"),
        "MPRA" => simple("MPRA"),
        other => unreachable!("family {other} is checked against FAMILY_SUBTYPES"),
    })
}

/// The metadata table with its first row per URN (`metadata_df[metadata_df['urn'] == urn].iloc[0]`).
pub(crate) struct Metadata<'a> {
    table: &'a Table,
    first_row: HashMap<&'a str, usize>,
}

impl<'a> Metadata<'a> {
    pub(crate) fn new(table: &'a Table) -> Metadata<'a> {
        let mut first_row = HashMap::new();
        if let Some(index) = table.column_index("urn") {
            for row in 0..table.rows {
                if let Cell::Text(urn) = table.columns[index].cell(row) {
                    first_row.entry(urn).or_insert(row);
                }
            }
        }
        Metadata { table, first_row }
    }

    /// The row of `urn`, if the table has it.
    pub(crate) fn row(&self, urn: &str) -> Option<usize> {
        self.first_row.get(urn).copied()
    }

    /// `meta.get(name)` of the row (`None` when the table has no such column).
    pub(crate) fn get(&self, row: usize, name: &str) -> Option<Cell<'a>> {
        self.table
            .column_index(name)
            .map(|index| self.table.columns[index].cell(row))
    }
}

/// Python's truthiness of `meta.get(name, default)`: NaN is true, empty text false.
fn truthy(value: Option<Cell<'_>>, default: bool) -> bool {
    match value {
        None => default,
        Some(Cell::Missing) => true,
        Some(Cell::Bool(flag)) => flag,
        Some(Cell::Int(number)) => number != 0,
        Some(Cell::Float(number)) => number != 0.0,
        Some(Cell::Text(text)) => !text.is_empty(),
    }
}

/// `f"{value}"` of a `meta.get(...)` result (`None` when the column is absent).
fn formatted(value: Option<Cell<'_>>) -> String {
    match value {
        None => "None".to_owned(),
        Some(cell) => python_str(cell),
    }
}

/// `float(text)`; ValueError is `None`. Text outside ASCII (where Python also reads other
/// scripts' digits) is refused.
pub(crate) fn python_float_of_text(text: &str) -> Result<Option<f64>> {
    if !text.is_ascii() {
        return Err(invalid(format!(
            "float({text:?}) of a non-ASCII text is not reproduced"
        )));
    }
    Ok(parse_float(text))
}

/// `float(value)` of a metadata value; ValueError is `None`.
fn python_float_of_cell(cell: Cell<'_>) -> Result<Option<f64>> {
    Ok(match cell {
        Cell::Missing => Some(f64::NAN),
        Cell::Float(number) => Some(number),
        Cell::Int(number) => Some(number as f64),
        Cell::Bool(flag) => Some(if flag { 1.0 } else { 0.0 }),
        Cell::Text(text) => python_float_of_text(text)?,
    })
}

/// `parse_neutral_band_from_metadata(meta)`: `(lo, hi)`, either `None` when unknown.
fn parse_neutral_band(metadata: &Metadata<'_>, row: usize) -> Result<(Option<f64>, Option<f64>)> {
    let raw = match metadata.get(row, "neutral_bands") {
        None => "[NA, NA]".to_owned(),
        Some(cell) => python_str(cell),
    };
    let band = raw.trim_matches([' ', '[', ']']);
    if !["[NA, NA]", "NA", "", "None", "nan"].contains(&band) {
        let mut parsed = Vec::new();
        let mut failed = false;
        for piece in band.split(',') {
            let piece = python_strip(piece);
            if ["NA", "nan", "None", ""].contains(&piece) {
                parsed.push(f64::NAN);
            } else {
                match python_float_of_text(piece)? {
                    Some(value) => parsed.push(value),
                    None => {
                        failed = true;
                        break;
                    }
                }
            }
        }
        if !failed && parsed.len() == 2 {
            let lo = (!parsed[0].is_nan()).then_some(parsed[0]);
            let hi = (!parsed[1].is_nan()).then_some(parsed[1]);
            let flexibility = truthy(metadata.get(row, "flexibility_requirement"), false);
            match (lo, hi) {
                (None, Some(hi)) => {
                    return Ok(if flexibility {
                        (Some(f64::NEG_INFINITY), Some(hi))
                    } else {
                        (Some(hi), Some(f64::INFINITY))
                    });
                }
                (Some(lo), None) => {
                    return Ok(if flexibility {
                        (Some(lo), Some(f64::INFINITY))
                    } else {
                        (Some(f64::NEG_INFINITY), Some(lo))
                    });
                }
                (Some(lo), Some(hi)) if lo < hi => return Ok((Some(lo), Some(hi))),
                _ => {}
            }
        }
    }
    if let (Some(q05), Some(q95)) = (metadata.get(row, "q05"), metadata.get(row, "q95")) {
        // float() of both; a ValueError skips the percentiles.
        if let (Some(q05), Some(q95)) = (python_float_of_cell(q05)?, python_float_of_cell(q95)?)
            && !(q05.is_nan() || q95.is_nan())
            && q05 < q95
        {
            return Ok((Some(q05), Some(q95)));
        }
    }
    Ok((None, None))
}

/// Python's `repr` of a text (`str(KeyError(name))`), for the reason of an unknown enum name.
fn python_text_repr(text: &str) -> String {
    let quote = if text.contains('\'') && !text.contains('"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::from(quote);
    for c in text.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c == quote => {
                out.push('\\');
                out.push(c);
            }
            c if (c as u32) < 0x20 || c as u32 == 0x7f => {
                out.push_str(&format!("\\x{:02x}", c as u32))
            }
            c => out.push(c),
        }
    }
    out.push(quote);
    out
}

/// `str(KeyError(value))` for an enum lookup of a metadata value.
fn key_error_text(value: Option<Cell<'_>>) -> String {
    match value {
        Some(Cell::Text(text)) => python_text_repr(text),
        other => formatted(other),
    }
}

/// `pick_interpreter_from_metadata(urn, score, metadata_df)`.
pub(crate) fn pick_interpreter_from_metadata(
    urn: &str,
    score: f64,
    metadata: &Metadata<'_>,
) -> Result<Interpretation> {
    let Some(row) = metadata.row(urn) else {
        return Ok(not_available("URN not found in metadata"));
    };
    if truthy(metadata.get(row, "ambiguous"), false) {
        return Ok(not_available("Dataset marked ambiguous"));
    }
    if truthy(metadata.get(row, "nonmonotonic"), false)
        || truthy(metadata.get(row, "needs_gene_specific_rules"), false)
    {
        return Ok(not_available(
            "Gene-specific rules required (not supported)",
        ));
    }
    let (Some(lo), Some(hi)) = parse_neutral_band(metadata, row)? else {
        return Ok(not_available("Unable to determine neutral band"));
    };
    let family_value = metadata.get(row, "assay_family");
    let subtype_value = metadata.get(row, "subtype_value");
    let Some((family, subtypes)) = FAMILY_SUBTYPES
        .iter()
        .find(|(name, _)| matches!(family_value, Some(Cell::Text(text)) if text == *name))
    else {
        return Ok(not_available(format!(
            "Unknown assay family: {}",
            formatted(family_value)
        )));
    };
    let Some(subtype) = subtypes
        .iter()
        .find(|name| matches!(subtype_value, Some(Cell::Text(text)) if text == **name))
    else {
        return Ok(not_available(format!(
            "Unknown subtype {} for {}",
            formatted(subtype_value),
            formatted(family_value)
        )));
    };
    let directionality_value = metadata
        .get(row, "directionality")
        .or(Some(Cell::Text("HIGHER_IS_FUNCTION")));
    let scale_value = metadata
        .get(row, "scale")
        .or(Some(Cell::Text("CENTERED_0")));
    let as_text = |value: Option<Cell<'_>>| match value {
        Some(Cell::Text(text)) => Some(text.to_owned()),
        _ => None,
    };
    let Some(directionality) = as_text(directionality_value)
        .as_deref()
        .and_then(Directionality::from_name)
    else {
        return Ok(not_available(format!(
            "Invalid enum value: {}",
            key_error_text(directionality_value)
        )));
    };
    let Some(scale) = as_text(scale_value).as_deref().and_then(Scale::from_name) else {
        return Ok(not_available(format!(
            "Invalid enum value: {}",
            key_error_text(scale_value)
        )));
    };
    let spec = AssaySpec {
        family,
        subtype,
        directionality,
        scale,
        neutral_band: (lo, hi),
        flexibility_requirement: truthy(metadata.get(row, "flexibility_requirement"), false),
    };
    let mut result = interpret_score(score, &spec)?;
    result.mavedb_ps3 = matches!(result.call, "LoF" | "GoF");
    result.mavedb_bs3 = result.call == "Neutral";
    Ok(result)
}

/// `priva acmg-step mavedb_interpreter pick_interpreter_from_metadata`: the input frame is the
/// metadata table; option `queries` is a list of `[urn, score]` with the score as the hex text
/// of its 64 bits. The output has one row per query: `urn`, `call`, `standardized_score` (NaN
/// for None), `standardized_is_none`, `reason`, `MaveDB_PS3`, `MaveDB_BS3`, `error` (the Python
/// exception type, or empty).
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    if function != "pick_interpreter_from_metadata" {
        return Err(not_converted("mavedb_interpreter", function));
    }
    let metadata = Metadata::new(&call.table);
    let Some(Json::Array(queries)) = call.options.get("queries") else {
        return Err(invalid(
            "option queries must be a list of [urn, score bits]",
        ));
    };
    let mut urns = Vec::new();
    let mut calls = Vec::new();
    let mut scores = Vec::new();
    let mut none = Vec::new();
    let mut reasons = Vec::new();
    let mut ps3 = Vec::new();
    let mut bs3 = Vec::new();
    let mut errors = Vec::new();
    for query in queries {
        let (Some(Json::Str(urn)), Some(Json::Str(bits))) =
            (query.get_index(0), query.get_index(1))
        else {
            return Err(invalid("each query is [urn, score bits]"));
        };
        let score = f64::from_bits(
            u64::from_str_radix(bits, 16).map_err(|_| invalid("score bits must be hex"))?,
        );
        urns.push(Some(urn.clone()));
        match pick_interpreter_from_metadata(urn, score, &metadata) {
            Ok(result) => {
                calls.push(Some(result.call.to_owned()));
                scores.push(result.standardized_score.unwrap_or(f64::NAN));
                none.push(result.standardized_score.is_none());
                reasons.push(Some(result.reason));
                ps3.push(result.mavedb_ps3);
                bs3.push(result.mavedb_bs3);
                errors.push(Some(String::new()));
            }
            Err(error) => {
                let kind = error
                    .downcast_ref::<AcmgError>()
                    .map(|python| python.python_type)
                    .ok_or(error.to_string())?;
                calls.push(Some(String::new()));
                scores.push(f64::NAN);
                none.push(true);
                reasons.push(Some(String::new()));
                ps3.push(false);
                bs3.push(false);
                errors.push(Some(kind.to_owned()));
            }
        }
    }
    let rows = urns.len();
    Ok(harness_options::frame(
        vec![
            ("urn".into(), Column::Text(urns)),
            ("call".into(), Column::Text(calls)),
            ("standardized_score".into(), Column::Float(scores)),
            ("standardized_is_none".into(), Column::Bool(none)),
            ("reason".into(), Column::Text(reasons)),
            ("MaveDB_PS3".into(), Column::Bool(ps3)),
            ("MaveDB_BS3".into(), Column::Bool(bs3)),
            ("error".into(), Column::Text(errors)),
        ],
        rows,
    ))
}

trait JsonIndex {
    fn get_index(&self, index: usize) -> Option<&Json>;
}

impl JsonIndex for Json {
    fn get_index(&self, index: usize) -> Option<&Json> {
        match self {
            Json::Array(items) => items.get(index),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn metadata_table(rows: &[[&str; 6]]) -> Table {
        let names = [
            "urn",
            "assay_family",
            "subtype_value",
            "directionality",
            "scale",
            "neutral_bands",
        ];
        Table {
            names: names.iter().map(|name| (*name).to_owned()).collect(),
            columns: (0..6)
                .map(|column| {
                    Column::Text(
                        rows.iter()
                            .map(|row| Some(row[column].to_owned()))
                            .collect(),
                    )
                })
                .collect(),
            rows: rows.len(),
        }
    }

    #[test]
    fn bands_and_calls_follow_the_python_rules() {
        let table = metadata_table(&[
            [
                "u1",
                "Fitness",
                "LOG2_ENRICHMENT",
                "HIGHER_IS_FUNCTION",
                "CENTERED_0",
                "[-0.5, 0.5]",
            ],
            [
                "u2",
                "Abundance",
                "VAMP_SEQ",
                "HIGHER_IS_FUNCTION",
                "WT_EQ_1",
                "[NA, 0.2]",
            ],
            ["u3", "Binding", "KD", "NON_MONOTONIC", "RAW", "[-1, 1]"],
            ["u4", "Other", "X", "HIGHER_IS_FUNCTION", "RAW", "[-1, 1]"],
            [
                "u5",
                "Fitness",
                "LOG2_ENRICHMENT",
                "SIDEWAYS",
                "RAW",
                "[-1, 1]",
            ],
        ]);
        let metadata = Metadata::new(&table);
        let pick = |urn: &str, score: f64| pick_interpreter_from_metadata(urn, score, &metadata);
        assert_eq!(pick("u1", 0.5).unwrap().call, "Neutral");
        assert!(pick("u1", 0.6).unwrap().mavedb_ps3);
        assert_eq!(pick("u1", f64::NAN).unwrap().call, "LoF");
        // [NA, 0.2] without flexibility: (0.2, inf); log2(1) = 0 is below the band.
        assert_eq!(pick("u2", 1.0).unwrap().call, "LoF");
        assert_eq!(
            pick("u2", 0.0).unwrap().reason,
            "Abundance not standardized"
        );
        let kd = pick("u3", 0.01).unwrap();
        assert_eq!(kd.call, "GoF");
        assert_eq!(
            kd.standardized_score,
            Some((1.0f64 / 0.01).ln() / 10f64.ln())
        );
        assert!(pick("u3", f64::INFINITY).is_err());
        assert_eq!(
            pick("u4", 0.0).unwrap().reason,
            "Unknown assay family: Other"
        );
        assert_eq!(
            pick("u5", 0.0).unwrap().reason,
            "Invalid enum value: 'SIDEWAYS'"
        );
        assert_eq!(pick("u9", 0.0).unwrap().reason, "URN not found in metadata");
    }
}
