//! `scripts/acmg_inheritance.py` (190 lines): the legacy gene-level inheritance flags.
//!
//! Thread T2 (gene-condition hub). Converted: `hpo_onset_modes`, `parse_hpo_inheritance` and
//! `identify_inheritance_mode_per_row`. In Python only the hub's diagnostic summary
//! (`GeneMechanismHub.known_inheritance_mode`, group 4) calls them; the ACMG step reads the hub's
//! per-variant columns instead.
//!
//! The row is a Python dict whose values can be of any type; `PyValue` keeps the kinds the
//! functions tell apart (`None` and NaN differ for `float()`).
use super::gene_mechanism_common::{
    json_column_value, python_error, python_str, python_text_repr, refused, result_table,
    show_json, show_outcome,
};
use super::types::{HarnessCall, not_converted};
use crate::cli::{Result, invalid};
use crate::hpo_penetrance::{Values, normalize_penetrance_hpo_ids};
use crate::json_text::Json;
use crate::pandas_table::{Cell, Column, ObjectValue, Table};
use crate::python_text::{parse_float, python_strip};

/// HPO terms of explicit incomplete, moderate or low penetrance and of late onset.
const INCOMPLETE_PENETRANCE_EQUIVALENT_TERMS: [&str; 4] =
    ["HP:0003829", "HP:4000159", "HP:4000160", "HP:0003584"];

/// A Python value of the row dict.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum PyValue<'a> {
    None,
    Nan,
    Text(&'a str),
    Int(i64),
    Float(f64),
    Bool(bool),
}

impl PyValue<'_> {
    fn truthy(&self) -> bool {
        match self {
            PyValue::None => false,
            PyValue::Nan => true,
            PyValue::Text(text) => !text.is_empty(),
            PyValue::Int(number) => *number != 0,
            PyValue::Float(number) => *number != 0.0,
            PyValue::Bool(flag) => *flag,
        }
    }

    /// `value == number` for a whole number.
    fn equals(&self, number: i64) -> bool {
        match self {
            PyValue::Int(value) => *value == number,
            PyValue::Float(value) => *value == number as f64,
            PyValue::Bool(flag) => i64::from(*flag) == number,
            _ => false,
        }
    }

    /// `str(value)`.
    fn text(&self) -> String {
        match self {
            PyValue::None => "None".to_owned(),
            PyValue::Nan => "nan".to_owned(),
            PyValue::Text(text) => (*text).to_owned(),
            PyValue::Int(number) => python_str(Cell::Int(*number)),
            PyValue::Float(number) => python_str(Cell::Float(*number)),
            PyValue::Bool(flag) => python_str(Cell::Bool(*flag)),
        }
    }

    /// `float(value)`, with Python's errors.
    fn to_float(self) -> Result<f64> {
        match self {
            PyValue::Nan => Ok(f64::NAN),
            PyValue::Int(number) => Ok(number as f64),
            PyValue::Float(number) => Ok(number),
            PyValue::Bool(flag) => Ok(if flag { 1.0 } else { 0.0 }),
            PyValue::None => Err(python_error(
                "TypeError",
                "float() argument must be a string or a real number, not 'NoneType'",
            )),
            PyValue::Text(text) => match parse_float(text) {
                Some(number) => Ok(number),
                None if !text.is_ascii() => {
                    Err(refused(format!("float() of the non-ASCII text {text:?}")))
                }
                None => Err(python_error(
                    "ValueError",
                    format!(
                        "could not convert string to float: {}",
                        python_text_repr(text)?
                    ),
                )),
            },
        }
    }
}

/// The row dict the functions read: `None` stands for an absent key.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct InheritanceRow<'a> {
    pub(crate) gene: Option<PyValue<'a>>,
    pub(crate) loeuf: Option<PyValue<'a>>,
    pub(crate) hpo_ids: Option<PyValue<'a>>,
    pub(crate) hpo_gene_inheritance: Option<PyValue<'a>>,
}

/// The flags `parse_hpo_inheritance` returns when the row has an HPO inheritance text.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct HpoInheritance {
    pub(crate) recessive: bool,
    pub(crate) dominant: bool,
    pub(crate) non_monogenic: bool,
    pub(crate) non_mendelian: bool,
    pub(crate) x_linked_recessive: bool,
    pub(crate) x_linked_dominant: bool,
    pub(crate) x_linked_unspecified: bool,
    pub(crate) incomplete_penetrance: bool,
}

/// What `parse_hpo_inheritance` returns: the incomplete-penetrance flag alone
/// when the row has no HPO inheritance text, else the dict of flags.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum ParsedInheritance {
    IncompletePenetranceOnly(bool),
    Flags(HpoInheritance),
}

/// `hpo_onset_modes(hpo_string)`: false when an incomplete-penetrance-equivalent
/// HPO term is among the `;`-separated terms of a text, else true.
pub(crate) fn hpo_onset_modes(value: PyValue<'_>) -> bool {
    let PyValue::Text(text) = value else {
        return true;
    };
    // A term matching `^HP:\d+$` (any case) is kept upper-cased; only the four
    // ASCII terms can then meet the set, so the digit rule needs no further check.
    !text.split(';').any(|term| {
        let upper = python_strip(term).to_uppercase();
        INCOMPLETE_PENETRANCE_EQUIVALENT_TERMS.contains(&upper.as_str())
    })
}

/// `parse_hpo_inheritance(row_dict)`.
pub(crate) fn parse_hpo_inheritance(row: &InheritanceRow<'_>) -> ParsedInheritance {
    let incomplete_penetrance = match row.hpo_ids {
        Some(PyValue::Text(text)) => {
            let terms: Vec<&str> = text.split(';').collect();
            normalize_penetrance_hpo_ids(Values::List(&terms)) == "incomplete"
        }
        _ => false,
    };
    let Some(PyValue::Text(text)) = row.hpo_gene_inheritance else {
        return ParsedInheritance::IncompletePenetranceOnly(incomplete_penetrance);
    };
    let modes: Vec<&str> = text.split(';').collect();
    let has = |label: &str| modes.contains(&label);
    ParsedInheritance::Flags(HpoInheritance {
        recessive: has("Autosomal recessive inheritance"),
        dominant: has("Autosomal dominant inheritance")
            || has("Autosomal dominant inheritance with maternal imprinting"),
        non_monogenic: has("Digenic inheritance")
            || has("Oligogenic inheritance")
            || has("Polygenic inheritance"),
        non_mendelian: has("Non-Mendelian inheritance"),
        x_linked_recessive: has("X-linked recessive inheritance"),
        x_linked_dominant: has("X-linked dominant inheritance"),
        x_linked_unspecified: has("X-linked inheritance"),
        incomplete_penetrance,
    })
}

/// `identify_inheritance_mode_per_row(row_dict, gene_mean_am_score,
/// clingen_curate_score, chromosome)`: (recessive, dominant, non-monogenic,
/// non-Mendelian, haploinsufficient, incomplete penetrance). The debug
/// messages Python builds read `row_dict['Gene']` (and `row_dict['LOEUF']`),
/// so their KeyErrors are raised too.
pub(crate) fn identify_inheritance_mode_per_row(
    row: &InheritanceRow<'_>,
    gene_mean_am_score: f64,
    clingen_curate_score: PyValue<'_>,
    chromosome: PyValue<'_>,
) -> Result<[bool; 6]> {
    let key_error = |name: &str| python_error("KeyError", format!("'{name}'"));
    let mut loeuf = row.loeuf.unwrap_or(PyValue::Float(0.6)).to_float()?;
    if loeuf.is_nan() {
        loeuf = 1.0;
    }
    let am = gene_mean_am_score;
    let mut haplo_insufficient =
        (loeuf <= 0.35 || am >= 0.564) && loeuf <= 0.7 && (am >= 0.5 || am.is_nan());
    let mut haplo_sufficient = !haplo_insufficient;

    let chromosome_text = python_strip(&chromosome.text()).to_lowercase();
    let first = chromosome_text.split(':').next().unwrap_or("");
    let normalized = first.strip_prefix("chr").unwrap_or(first);
    let sex_chromosome = normalized == "x" || normalized == "y";
    let mut clingen_recessive = None;
    let mut clingen_dominant = None;
    if clingen_curate_score.truthy() {
        if row.gene.is_none() {
            return Err(key_error("Gene"));
        }
        if clingen_curate_score.equals(3) {
            if !sex_chromosome {
                clingen_recessive = Some(false);
                clingen_dominant = Some(true);
            }
            haplo_insufficient = true;
        } else if clingen_curate_score.equals(30) || clingen_curate_score.equals(40) {
            if !sex_chromosome {
                clingen_recessive = Some(true);
            }
            haplo_insufficient = false;
            haplo_sufficient = true;
        }
    }

    let mut flags = match parse_hpo_inheritance(row) {
        ParsedInheritance::IncompletePenetranceOnly(incomplete) => {
            if row.gene.is_none() {
                return Err(key_error("Gene"));
            }
            if row.loeuf.is_none() {
                return Err(key_error("LOEUF"));
            }
            if sex_chromosome {
                return Ok([false, false, false, false, haplo_insufficient, incomplete]);
            }
            return Ok([
                clingen_recessive.unwrap_or(haplo_sufficient),
                clingen_dominant.unwrap_or(haplo_insufficient),
                false,
                false,
                haplo_insufficient,
                incomplete,
            ]);
        }
        ParsedInheritance::Flags(flags) => flags,
    };
    let x_specific =
        flags.x_linked_recessive || flags.x_linked_dominant || flags.x_linked_unspecified;
    if let Some(recessive) = clingen_recessive.filter(|_| !x_specific) {
        flags.recessive = recessive;
    }
    if let Some(dominant) = clingen_dominant.filter(|_| !x_specific) {
        flags.dominant = dominant;
    }
    if flags.recessive {
        haplo_insufficient = false;
    }
    Ok([
        flags.recessive,
        flags.dominant,
        flags.non_monogenic,
        flags.non_mendelian,
        haplo_insufficient,
        flags.incomplete_penetrance,
    ])
}

// ------------------------------------------------------------------------------- harness

/// The Python value of a frame cell (`None` and NaN kept apart).
fn py_value(column: &Column, row: usize) -> Result<PyValue<'_>> {
    Ok(match column {
        Column::Object(values) => match &values[row] {
            ObjectValue::None => PyValue::None,
            ObjectValue::Nan => PyValue::Nan,
            ObjectValue::Text(text) => PyValue::Text(text),
            ObjectValue::Int(number) => PyValue::Int(*number),
            ObjectValue::Float(number) if number.is_nan() => PyValue::Nan,
            ObjectValue::Float(number) => PyValue::Float(*number),
            ObjectValue::Bool(flag) => PyValue::Bool(*flag),
        },
        Column::NullableInt(_) | Column::NullableBool(_) => {
            return Err(refused("pd.NA in an inheritance row"));
        }
        other => match other.cell(row) {
            Cell::Missing => PyValue::Nan,
            Cell::Text(text) => PyValue::Text(text),
            Cell::Int(number) => PyValue::Int(number),
            Cell::Float(number) => PyValue::Float(number),
            Cell::Bool(flag) => PyValue::Bool(flag),
        },
    })
}

fn flags_json(flags: &HpoInheritance) -> Json {
    let field = |name: &str, flag: bool| (name.to_owned(), Json::Bool(flag));
    Json::Object(vec![
        field("hpo_recessive", flags.recessive),
        field("hpo_dominant", flags.dominant),
        field("hpo_non_monogenic", flags.non_monogenic),
        field("hpo_non_mendelian", flags.non_mendelian),
        field("hpo_x_linked_recessive", flags.x_linked_recessive),
        field("hpo_x_linked_dominant", flags.x_linked_dominant),
        field("hpo_x_linked_unspecified", flags.x_linked_unspecified),
        field("incomplete_penetrance", flags.incomplete_penetrance),
    ])
}

/// `priva acmg-step acmg_inheritance <function>`: `hpo_onset_modes` (column
/// `value`); `parse_hpo_inheritance` and `identify_inheritance_mode_per_row`
/// (the row dict's keys `Gene`, `LOEUF`, `HPO_IDs`, `HPO_gene_inheritance` as
/// columns, absent columns absent keys; the latter also `gene_mean_am_score`,
/// `clingen_curate_score`, `chromosome`). A JSON column `absent` may list keys
/// a row leaves out.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    const FUNCTIONS: [&str; 3] = [
        "hpo_onset_modes",
        "parse_hpo_inheritance",
        "identify_inheritance_mode_per_row",
    ];
    if !FUNCTIONS.contains(&function) {
        return Err(not_converted("acmg_inheritance", function));
    }
    let table = &call.table;
    let mut results = Vec::with_capacity(table.rows);
    for row in 0..table.rows {
        let value = |name: &str| -> Result<Option<PyValue<'_>>> {
            let absent = match table.column("absent") {
                Ok(_) => match json_column_value(table, "absent", row)? {
                    Json::Array(names) => names.iter().any(|item| item.as_str() == Some(name)),
                    _ => return Err(invalid("absent must be a JSON list")),
                },
                Err(_) => false,
            };
            match table.column(name) {
                Ok(column) if !absent => Ok(Some(py_value(column, row)?)),
                _ => Ok(None),
            }
        };
        let row_dict = || -> Result<InheritanceRow<'_>> {
            Ok(InheritanceRow {
                gene: value("Gene")?,
                loeuf: value("LOEUF")?,
                hpo_ids: value("HPO_IDs")?,
                hpo_gene_inheritance: value("HPO_gene_inheritance")?,
            })
        };
        let outcome: Result<String> = match function {
            "hpo_onset_modes" => (|| {
                let value = value("value")?.ok_or_else(|| invalid("the column value is needed"))?;
                Ok(format!(
                    "b:{}",
                    if hpo_onset_modes(value) {
                        "True"
                    } else {
                        "False"
                    }
                ))
            })(),
            "parse_hpo_inheritance" => (|| match parse_hpo_inheritance(&row_dict()?) {
                ParsedInheritance::IncompletePenetranceOnly(flag) => {
                    Ok(format!("b:{}", if flag { "True" } else { "False" }))
                }
                ParsedInheritance::Flags(flags) => show_json(&flags_json(&flags)),
            })(),
            "identify_inheritance_mode_per_row" => (|| {
                let am = match value("gene_mean_am_score")? {
                    Some(PyValue::Float(number)) => number,
                    Some(PyValue::Nan) => f64::NAN,
                    _ => return Err(invalid("gene_mean_am_score must be a float")),
                };
                let clingen = value("clingen_curate_score")?.unwrap_or(PyValue::None);
                let chromosome = value("chromosome")?.unwrap_or(PyValue::Text(""));
                let result =
                    identify_inheritance_mode_per_row(&row_dict()?, am, clingen, chromosome)?;
                show_json(&Json::Array(
                    result.iter().map(|flag| Json::Bool(*flag)).collect(),
                ))
            })(),
            _ => return Err(not_converted("acmg_inheritance", function)),
        };
        results.push(show_outcome(outcome)?);
    }
    Ok(result_table(results))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn onset_terms_block_unaffected_carriers() {
        assert!(!hpo_onset_modes(PyValue::Text("HP:0000001; hp:0003584 ")));
        assert!(hpo_onset_modes(PyValue::Text("HP:0003831")));
        assert!(hpo_onset_modes(PyValue::Nan));
        assert!(hpo_onset_modes(PyValue::Text("")));
    }

    #[test]
    fn inheritance_flags_and_modes() {
        let row = InheritanceRow {
            gene: Some(PyValue::Text("G")),
            loeuf: Some(PyValue::Float(0.2)),
            hpo_ids: Some(PyValue::Text("HP:0003829")),
            hpo_gene_inheritance: Some(PyValue::Text(
                "Autosomal recessive inheritance;X-linked inheritance",
            )),
        };
        let ParsedInheritance::Flags(flags) = parse_hpo_inheritance(&row) else {
            panic!("flags expected");
        };
        assert!(flags.recessive && flags.x_linked_unspecified && flags.incomplete_penetrance);
        assert_eq!(
            identify_inheritance_mode_per_row(
                &row,
                f64::NAN,
                PyValue::Int(3),
                PyValue::Text("chr1")
            )
            .unwrap(),
            [true, false, false, false, false, true]
        );
        let row = InheritanceRow {
            gene: Some(PyValue::Text("G")),
            loeuf: Some(PyValue::Nan),
            ..InheritanceRow::default()
        };
        assert_eq!(
            identify_inheritance_mode_per_row(&row, 0.6, PyValue::None, PyValue::Text("X"))
                .unwrap(),
            [false, false, false, false, false, false]
        );
        let row = InheritanceRow {
            gene: Some(PyValue::Text("G")),
            ..InheritanceRow::default()
        };
        assert_eq!(
            identify_inheritance_mode_per_row(&row, 0.6, PyValue::None, PyValue::Text("1"))
                .unwrap_err()
                .to_string(),
            "KeyError: 'LOEUF'"
        );
        let row = InheritanceRow {
            loeuf: Some(PyValue::None),
            ..InheritanceRow::default()
        };
        assert!(
            identify_inheritance_mode_per_row(&row, 0.6, PyValue::None, PyValue::Text("1"))
                .unwrap_err()
                .to_string()
                .starts_with("TypeError: float() argument")
        );
    }
}
