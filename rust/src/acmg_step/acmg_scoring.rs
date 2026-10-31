//! `scripts/acmg_scoring.py` (836 lines): the scoring tail: criteria matrix and ClinGen override,
//! the posterior per row, `compute_control_common`, and the glue that runs the existing posterior
//! and ranking modules in the Python order.
//!
//! Thread T3 (exact-allele matching and the scoring tail). Converted: `summarize_acmg_criteria`
//! (matrix construction and validation, the ClinGen map, the summary column and the override of
//! the matrix rows), through `crate::acmg_scoring`'s `create_criteria_summary` and
//! `translate_strength_level`; the orchestrator's posterior step (`criteria_matrix.apply(
//! calculate_posterior_probability, axis=1)` and the two columns it writes), through
//! `crate::acmg_scoring::calculate_posterior_probability`; and `compute_control_common`, which
//! calls `stats::binom_cdf` and `stats::beta_ppf` (thread T5; until they land, a row whose control
//! counts are tested stops at their placeholder); and `sort_and_rank_variants` (stages R1 to R8:
//! protein-coding filter, the proband from the pedigree, the control-common penalty, the
//! zygosity, inheritance and mechanism compatibility penalty from the mechanism masks, the
//! dispensable and relevant gene lists, the pext modulation through the earlier
//! `coerce_pext_value`, and the final order and rank).
use super::acmg_ps2_pm6_denovo::{family_rows, harness_ped, python_str};
use super::acmg_variant_mechanism::variant_mechanism_masks;
use super::stats;
use super::types::{AcmgError, CRITERIA_ORDER, HarnessCall, not_converted};
use crate::acmg_consequence::coerce_pext_value::{Value as PextValue, coerce_pext_value};
use crate::acmg_scoring::SourceError;
use crate::acmg_scoring::calculate_posterior_probability::{
    Probability, calculate_posterior_probability,
};
use crate::acmg_scoring::create_criteria_summary::create_criteria_summary;
use crate::acmg_scoring::translate_strength_level::translate_strength_level;
use crate::acmg_variant_mechanism::variant_mechanism_masks::Cell as MaskCell;
use crate::cli::{Result, invalid};
use crate::json_text::Json;
use crate::pandas_table::{Cell, Column, ObjectValue, Table};
use crate::pandas_table::{Header, ReadOptions, read_table_with};
use crate::pickle_reader::{self, Value};
use crate::python_text::python_strip;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::path::Path;

/// `priva acmg-step acmg_scoring <function>`: runs one converted function on a frame.
///
/// `summarize_acmg_criteria`: the frame holds the table's columns and the criteria as columns
/// named `criterion:NAME`; options `clingen_map_pkl` (text or null) and
/// `apply_clingen_override`. The result is the criteria matrix with the table's
/// `ACMG_criteria` appended as `__ACMG_criteria`.
/// `calculate_posterior_probability`: the frame is the criteria matrix; the result holds
/// `ACMG_quant_score` and `ACMG_class`.
/// `compute_control_common`: option `proband` (text or null); the result is the table.
/// `sort_and_rank_variants`: options `ped_frame` (a frame file, or absent for `ped_df=None`),
/// `fam_name`, `pext_tissues`, `pext_low_expression_cutoff`, `pext_penalty_floor`,
/// `pext_penalty_shape`, `relevant_gene_list`, `dispensable_gene_list` (null for None); the
/// result is the sorted table.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    match function {
        "summarize_acmg_criteria" => {
            let mut table = call.table.clone();
            let mut criteria = Vec::new();
            let mut keep = Vec::new();
            for (index, name) in table.names.iter().enumerate() {
                match name.strip_prefix("criterion:") {
                    Some(criterion) => {
                        let Some(known) = CRITERIA_ORDER.iter().find(|c| **c == criterion) else {
                            return Err(invalid(format!("unknown criterion {criterion}")));
                        };
                        criteria.push((*known, table.columns[index].clone()));
                    }
                    None => keep.push(index),
                }
            }
            table.names = keep
                .iter()
                .map(|index| table.names[*index].clone())
                .collect();
            table.columns = keep
                .iter()
                .map(|index| table.columns[*index].clone())
                .collect();
            let clingen = match call.options.get("clingen_map_pkl") {
                Some(Json::Str(path)) => Some(path.as_str()),
                _ => None,
            };
            let apply = !matches!(
                call.options.get("apply_clingen_override"),
                Some(Json::Bool(false))
            );
            let mut matrix = summarize_acmg_criteria(&mut table, &criteria, clingen, apply)?;
            let summary = table.column("ACMG_criteria")?.clone();
            matrix.set_column("__ACMG_criteria", summary);
            Ok(matrix)
        }
        "calculate_posterior_probability" => {
            let (scores, classes) = posterior_columns(&call.table)?;
            Ok(Table {
                names: vec!["ACMG_quant_score".into(), "ACMG_class".into()],
                columns: vec![scores, classes],
                rows: call.table.rows,
            })
        }
        "compute_control_common" => {
            let mut table = call.table.clone();
            let proband = match call.options.get("proband") {
                Some(Json::Str(name)) => Some(name.as_str()),
                _ => None,
            };
            compute_control_common(&mut table, proband)?;
            Ok(table)
        }
        "sort_and_rank_variants" => {
            let ped = match call.options.get("ped_frame") {
                Some(Json::Str(_)) => Some(harness_ped(call)?),
                _ => None,
            };
            let text = |name: &str| match call.options.get(name) {
                Some(Json::Str(value)) => Some(value.as_str()),
                _ => None,
            };
            let number = |name: &str, default: f64| match call.options.get(name) {
                Some(Json::Float(value)) => *value,
                Some(Json::Int(value)) => *value as f64,
                _ => default,
            };
            let options = RankingOptions {
                pext_tissues: text("pext_tissues").unwrap_or(""),
                pext_low_expression_cutoff: number("pext_low_expression_cutoff", 0.1),
                pext_penalty_floor: number("pext_penalty_floor", 0.8),
                pext_penalty_shape: number("pext_penalty_shape", 0.5),
                relevant_gene_list: text("relevant_gene_list"),
                dispensable_gene_list: text("dispensable_gene_list"),
            };
            sort_and_rank_variants(&call.table, ped.as_ref(), text("fam_name"), &options)
        }
        _ => Err(not_converted("acmg_scoring", function)),
    }
}

fn source_error(error: SourceError) -> Box<dyn std::error::Error> {
    Box::new(AcmgError::new(error.0, error.1))
}

/// A criterion strength as Python's `int(value)` reads a matrix cell: an integer,
/// or a boolean (numpy's True is 1). Other kinds are refused.
fn strength(column: &Column, row: usize, name: &str) -> Result<i64> {
    match column.cell(row) {
        Cell::Int(number) => Ok(number),
        Cell::Bool(flag) => Ok(i64::from(flag)),
        _ => Err(invalid(format!(
            "refused (not reproduced): criterion {name} holds a {} value",
            column.dtype_name()
        ))),
    }
}

/// `repr(path)` of a text path, as the FileNotFoundError message shows it.
fn python_repr_path(path: Option<&str>) -> Result<String> {
    match path {
        None => Ok("None".to_owned()),
        Some(text)
            if text
                .chars()
                .all(|c| (' '..='~').contains(&c) && c != '\'' && c != '\\') =>
        {
            Ok(format!("'{text}'"))
        }
        Some(text) => Err(invalid(format!(
            "refused (not reproduced): repr() of the ClinGen map path {text:?}"
        ))),
    }
}

/// The ClinGen map pickle: variant id to its evidence codes; an entry that is not
/// text (the map's float NaN for a curated variant without codes) is `None`.
fn load_clingen_map(path: &str) -> Result<HashMap<String, Option<String>>> {
    let bytes = pickle_reader::read_file(Path::new(path), path.ends_with(".gz"))?;
    let Value::Dict(dict) = pickle_reader::load(&bytes)? else {
        return Err(invalid(
            "refused (not reproduced): the ClinGen map is not a dict",
        ));
    };
    let mut map = HashMap::new();
    for (key, value) in &dict.borrow().entries {
        // Only text keys can equal a variant id.
        if let Value::Str(key) = key {
            let value = match value {
                Value::Str(text) => Some(text.to_string()),
                _ => None,
            };
            map.insert(key.to_string(), value);
        }
    }
    Ok(map)
}

/// `summarize_acmg_criteria(df, criteria_dict, clingen_map_pkl, apply_clingen_override)`:
/// writes `ACMG_criteria` into the table and returns the criteria matrix (the 28
/// criteria in `CRITERIA_ORDER`, absent ones as int64 zeros, then `variant_id`).
///
/// With a non-empty ClinGen map, each criterion column is rebuilt for the overridden
/// rows through a `.loc` assignment, which changes some columns to `object` the way
/// pandas does (see the comment there); the step's real run has all 28 as `object`.
pub(crate) fn summarize_acmg_criteria(
    table: &mut Table,
    criteria: &[(&str, Column)],
    clingen_map_pkl: Option<&str>,
    apply_clingen_override: bool,
) -> Result<Table> {
    let rows = table.rows;
    let mut columns: Vec<Column> = Vec::with_capacity(CRITERIA_ORDER.len());
    for name in CRITERIA_ORDER {
        let column = match criteria.iter().find(|(criterion, _)| *criterion == name) {
            Some((_, values)) => values.clone(),
            None => Column::Int(vec![0; rows]),
        };
        if column.len() != rows {
            return Err(invalid(format!(
                "refused (not reproduced): criterion {name} has {} values for {rows} rows",
                column.len()
            )));
        }
        columns.push(column);
    }
    let variant_id = table.column("variant_id")?.clone();
    if rows == 0 {
        return Err(invalid(
            "refused (not reproduced): summarize_acmg_criteria of an empty table",
        ));
    }
    for (name, column) in CRITERIA_ORDER.iter().zip(&columns) {
        let valid = match column {
            Column::Int(values) => values.iter().all(|value| (0..=5).contains(value)),
            Column::Bool(_) => true,
            _ => {
                return Err(invalid(format!(
                    "refused (not reproduced): criterion {name} is a {} column",
                    column.dtype_name()
                )));
            }
        };
        if !valid {
            return Err(Box::new(AcmgError::new(
                "ValueError",
                format!("Column {name} has non-integer values"),
            )));
        }
    }

    let clingen_map = if apply_clingen_override {
        let path = clingen_map_pkl.filter(|path| !path.is_empty() && Path::new(path).exists());
        let Some(path) = path else {
            return Err(Box::new(AcmgError::new(
                "FileNotFoundError",
                format!(
                    "ClinGen override is enabled but the evidence map is missing: {}. Build it \
                     with scripts/prepare_clingen_map.py, or pass apply_clingen_override=False to \
                     benchmark PriVA against ClinGen-curated variants.",
                    python_repr_path(clingen_map_pkl)?
                ),
            )));
        };
        load_clingen_map(path)?
    } else {
        HashMap::new()
    };

    let evidence: Vec<Option<&str>> = (0..rows)
        .map(|row| match variant_id.cell(row) {
            Cell::Text(id) => clingen_map.get(id).and_then(|value| value.as_deref()),
            _ => None,
        })
        .collect();
    let mut summaries: Vec<Option<String>> = Vec::with_capacity(rows);
    for (row, evidence) in evidence.iter().enumerate() {
        let strengths: Vec<(&str, i64)> = CRITERIA_ORDER
            .iter()
            .zip(&columns)
            .map(|(name, column)| Ok((*name, strength(column, row, name)?)))
            .collect::<Result<_>>()?;
        summaries.push(create_criteria_summary(&strengths, *evidence).map_err(source_error)?);
    }
    let summary_all_missing = summaries.iter().all(Option::is_none);
    let summary_column = if summary_all_missing {
        Column::Float(vec![f64::NAN; rows])
    } else {
        Column::Text(summaries.clone())
    };
    table.set_column("ACMG_criteria", summary_column);

    if !clingen_map.is_empty() {
        let overridden: Vec<bool> = evidence
            .iter()
            .map(|value| value.is_some_and(|text| !python_strip(text).is_empty()))
            .collect();
        let any = overridden.iter().any(|flag| *flag);
        for (name, column) in CRITERIA_ORDER.iter().zip(columns.iter_mut()) {
            let mut updates = Vec::new();
            for row in (0..rows).filter(|row| overridden[*row]) {
                updates.push((
                    row,
                    translate_strength_level(name, summaries[row].as_deref())
                        .map_err(source_error)?,
                ));
            }
            // pandas assigns `.map(...)` of the overridden rows: an int64 Series when a row is
            // overridden, else an empty Series of the summary column's type (object when any
            // summary is text, float64 when all are NaN). int64 holds int64 and an empty float64,
            // and becomes object for an empty object Series; bool becomes object in every case.
            let to_object = match &*column {
                Column::Int(_) => !any && !summary_all_missing,
                _ => true,
            };
            let mut values: Vec<ObjectValue> = match &*column {
                Column::Int(values) => values
                    .iter()
                    .map(|value| ObjectValue::Int(*value))
                    .collect(),
                Column::Bool(values) => values
                    .iter()
                    .map(|value| ObjectValue::Bool(*value))
                    .collect(),
                _ => unreachable!("criteria are int64 or bool here"),
            };
            for (row, value) in updates {
                values[row] = ObjectValue::Int(value);
            }
            *column = if to_object {
                Column::Object(values)
            } else {
                Column::Int(
                    values
                        .iter()
                        .map(|value| match value {
                            ObjectValue::Int(number) => *number,
                            _ => unreachable!("an int64 column"),
                        })
                        .collect(),
                )
            };
        }
    }

    let mut names: Vec<String> = CRITERIA_ORDER
        .iter()
        .map(|name| (*name).to_owned())
        .collect();
    names.push("variant_id".to_owned());
    columns.push(variant_id);
    Ok(Table {
        names,
        columns,
        rows,
    })
}

/// The orchestrator's posterior step: `zip(*criteria_matrix.apply(
/// calculate_posterior_probability, axis=1))` as the columns `ACMG_quant_score`
/// and `ACMG_class`. Only columns named with a leading `P` or `B` count, as in Python.
pub(crate) fn posterior_columns(matrix: &Table) -> Result<(Column, Column)> {
    if matrix.rows == 0 {
        return Err(Box::new(AcmgError::new(
            "ValueError",
            "not enough values to unpack (expected 2, got 0)",
        )));
    }
    let counted: Vec<(&str, &Column)> = matrix
        .names
        .iter()
        .zip(&matrix.columns)
        .filter(|(name, _)| name.starts_with('P') || name.starts_with('B'))
        .map(|(name, column)| (name.as_str(), column))
        .collect();
    let mut scores = Vec::with_capacity(matrix.rows);
    let mut classes = Vec::with_capacity(matrix.rows);
    for row in 0..matrix.rows {
        let mut strengths = Vec::with_capacity(counted.len());
        for (name, column) in &counted {
            strengths.push((*name, strength(column, row, name)?));
        }
        if strengths
            .iter()
            .any(|(_, value)| !(0..=5).contains(value) && *value > 0)
        {
            return Err(invalid(
                "refused (not reproduced): a strength outside 0..5 in the posterior",
            ));
        }
        let (score, class) = calculate_posterior_probability(&strengths);
        scores.push(score);
        classes.push(Some(class.to_owned()));
    }
    // A tuple of floats and the integer 0 (BA1) becomes float64; integers alone stay int64.
    let scores = if scores
        .iter()
        .all(|score| matches!(score, Probability::Integer(_)))
    {
        Column::Int(
            scores
                .iter()
                .map(|score| match score {
                    Probability::Integer(value) => *value,
                    Probability::Float(_) => unreachable!("all integers"),
                })
                .collect(),
        )
    } else {
        Column::Float(scores.iter().map(Probability::as_f64).collect())
    };
    Ok((scores, Column::Text(classes)))
}

/// `pd.to_numeric(df.get(name), errors="coerce")` for the column types the control
/// columns have (`None` when the column is absent: Python's `to_numeric(None)` is
/// reached only through `df.get` and would raise, so an absent column is refused
/// unless the four are all absent).
fn control_numbers(table: &Table, name: &str) -> Result<Option<Vec<f64>>> {
    let Ok(column) = table.column(name) else {
        return Ok(None);
    };
    match column {
        Column::Int(values) => Ok(Some(values.iter().map(|value| *value as f64).collect())),
        Column::Float(values) => Ok(Some(values.clone())),
        Column::NullableInt(values) => Ok(Some(
            values
                .iter()
                .map(|value| value.map_or(f64::NAN, |number| number as f64))
                .collect(),
        )),
        _ => Err(invalid(format!(
            "refused (not reproduced): control column {name} is {}",
            column.dtype_name()
        ))),
    }
}

/// `compute_control_common(df, proband, cutoff=0.01, conf_level=0.999)`: writes
/// `control_common` and `control_common_index`.
pub(crate) fn compute_control_common(table: &mut Table, proband: Option<&str>) -> Result<()> {
    const CUTOFF: f64 = 0.01;
    const CONF_LEVEL: f64 = 0.999;
    let rows = table.rows;
    let names = ["control_AC", "control_AN", "control_AF", "control_nhomalt"];
    if !names.iter().any(|name| table.has_column(name)) {
        table.set_column("control_common", Column::Bool(vec![false; rows]));
        table.set_column("control_common_index", Column::Float(vec![1.0; rows]));
        return Ok(());
    }
    let mut values = Vec::new();
    for name in names {
        match control_numbers(table, name)? {
            Some(numbers) => values.push(numbers),
            None => {
                return Err(invalid(format!(
                    "refused (not reproduced): {name} is absent while another control column exists"
                )));
            }
        }
    }
    let [ac, an, af, nhomo]: [Vec<f64>; 4] = values.try_into().expect("four columns");
    let valid_an: Vec<bool> = an
        .iter()
        .map(|value| !value.is_nan() && *value > 0.0)
        .collect();
    let is_hom: Vec<bool> = match proband.and_then(|name| table.column(name).ok()) {
        Some(column) => (0..rows)
            .map(|row| {
                let text = match column.cell(row) {
                    Cell::Text(text) => text.to_owned(),
                    Cell::Missing => "nan".to_owned(),
                    _ => {
                        return Err(invalid(
                            "refused (not reproduced): a genotype that is not text",
                        ));
                    }
                };
                let gt = text.split(':').next().unwrap_or("");
                Ok(["1/1", "1|1", "1"].contains(&gt))
            })
            .collect::<Result<_>>()?,
        None => vec![false; rows],
    };
    let use_nhomo: Vec<bool> = (0..rows)
        .map(|row| valid_an[row] && is_hom[row] && !nhomo[row].is_nan())
        .collect();
    let use_ac: Vec<bool> = (0..rows)
        .map(|row| valid_an[row] && !use_nhomo[row] && !ac[row].is_nan())
        .collect();
    let use_af: Vec<bool> = (0..rows)
        .map(|row| valid_an[row] && !use_nhomo[row] && !use_ac[row] && !af[row].is_nan())
        .collect();
    // `.clip(lower=...).astype(int)`: truncation of finite values; `round()` is half to even.
    let as_int = |value: f64| -> Result<i64> {
        if !value.is_finite() || value.abs() >= 9.2e18 {
            return Err(invalid(
                "refused (not reproduced): a control count that does not fit an integer",
            ));
        }
        Ok(value.trunc() as i64)
    };
    let mut common = vec![false; rows];
    for row in 0..rows {
        if use_nhomo[row] {
            let k = as_int(nhomo[row].max(0.0))?;
            let n = as_int((an[row] / 2.0).round_ties_even().max(1.0))?;
            common[row] = stats::binom_cdf(k as f64, n as f64, CUTOFF) > CONF_LEVEL;
        }
    }
    for row in 0..rows {
        if use_ac[row] {
            let k = as_int(ac[row].max(0.0))?;
            let n = as_int(an[row].round_ties_even().max(1.0))?;
            common[row] = stats::binom_cdf(k as f64, n as f64, CUTOFF) > CONF_LEVEL;
        }
    }
    for row in 0..rows {
        if use_af[row] {
            common[row] = af[row] > CUTOFF;
        }
    }
    // `common[m] = af[m] > cutoff` assigns a Series; when it covers only some rows pandas
    // aligns it with NaN elsewhere and `common` becomes an object column. numpy then
    // refuses that object array as an index as soon as one row is common.
    let af_partial = use_af.iter().any(|flag| *flag) && !use_af.iter().all(|flag| *flag);
    if af_partial {
        let indexed = (0..rows).any(|row| common[row] && !an[row].is_nan() && an[row] > 0.0);
        if indexed {
            return Err(Box::new(AcmgError::new(
                "IndexError",
                "arrays used as indices must be of integer (or boolean) type",
            )));
        }
        table.set_column(
            "control_common",
            Column::BoolOrMissing(common.iter().map(|flag| Some(*flag)).collect()),
        );
        table.set_column("control_common_index", Column::Float(vec![1.0; rows]));
        return Ok(());
    }
    table.set_column("control_common", Column::Bool(common.clone()));
    let mut index = Vec::with_capacity(rows);
    for row in 0..rows {
        let ac_safe = if ac[row].is_nan() { 0.0 } else { ac[row] };
        let an_safe = if an[row].is_nan() { 0.0 } else { an[row] };
        index.push(if common[row] {
            let confident = if an_safe > 0.0 {
                stats::beta_ppf(0.01, ac_safe + 1.0, an_safe - ac_safe + 1.0)
            } else {
                0.0
            };
            (1.0 - confident).clamp(0.01, 0.5)
        } else {
            1.0
        });
    }
    table.set_column("control_common_index", Column::Float(index));
    Ok(())
}

/// The keyword arguments of `sort_and_rank_variants` after `fam_name`, as the step
/// passes them (`None` for a gene list stands for Python's `None`: no list).
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct RankingOptions<'a> {
    pub(crate) pext_tissues: &'a str,
    pub(crate) pext_low_expression_cutoff: f64,
    pub(crate) pext_penalty_floor: f64,
    pub(crate) pext_penalty_shape: f64,
    pub(crate) relevant_gene_list: Option<&'a str>,
    pub(crate) dispensable_gene_list: Option<&'a str>,
}

fn python_error(
    python_type: &'static str,
    message: impl Into<String>,
) -> Box<dyn std::error::Error> {
    Box::new(AcmgError::new(python_type, message))
}

fn refuse(message: impl std::fmt::Display) -> Box<dyn std::error::Error> {
    invalid(format!("refused (not reproduced): {message}"))
}

/// A float64 column's values (NaN for missing); other column types are refused.
fn float_values(table: &Table, name: &str) -> Result<Vec<f64>> {
    match table.column(name)? {
        Column::Float(values) => Ok(values.clone()),
        other => Err(refuse(format!(
            "{name} is {} rather than float64",
            other.dtype_name()
        ))),
    }
}

/// `str(value).strip().lower()` of a pedigree cell.
fn sex_text(cell: Option<Cell<'_>>) -> String {
    let text = match cell {
        None => String::new(),
        Some(cell) => python_str(cell).unwrap_or_else(|| "nan".to_owned()),
    };
    python_strip(&text).to_lowercase()
}

/// R2: the proband, its sex flags and `proband_het` (`df[proband].str.count("1") == 1`).
fn pedigree_context(
    table: &Table,
    ped: Option<&Table>,
    fam_name: Option<&str>,
) -> Result<(Option<String>, bool, bool, Vec<bool>)> {
    let (Some(ped), Some(fam_name)) = (ped, fam_name) else {
        return Ok((None, false, false, vec![false; table.rows]));
    };
    let family = family_rows(ped, fam_name)?;
    let phenotypes = ped.column("Phenotype")?;
    let first = family.into_iter().find(|row| match phenotypes.cell(*row) {
        Cell::Text(text) => text == "2",
        Cell::Int(number) => number == 2,
        Cell::Float(number) => number == 2.0,
        _ => false,
    });
    let Some(row) = first else {
        return Err(python_error(
            "IndexError",
            "single positional indexer is out-of-bounds",
        ));
    };
    let proband = match ped.column("IndividualID")?.cell(row) {
        Cell::Text(text) => text.to_owned(),
        other => return Err(refuse(format!("a proband ID that is not text ({other:?})"))),
    };
    let sex = sex_text(ped.column("Sex").ok().map(|column| column.cell(row)));
    let male = ["1", "1.0", "m", "male"].contains(&sex.as_str());
    let female = ["2", "2.0", "f", "female"].contains(&sex.as_str());
    let not_text = || {
        python_error(
            "AttributeError",
            "Can only use .str accessor with string values!",
        )
    };
    let het = match table.column(&proband)? {
        Column::Text(values) => values
            .iter()
            .map(|value| {
                value
                    .as_deref()
                    .is_some_and(|text| text.matches('1').count() == 1)
            })
            .collect(),
        Column::Object(_) => {
            return Err(refuse(format!(
                "the proband column {proband} mixes kinds of values"
            )));
        }
        _ => return Err(not_text()),
    };
    Ok((Some(proband), male, female, het))
}

/// `pd.read_table(path, header=None, names=["SYMBOL"], comment="#")` as a set of
/// texts. A missing value in the list is refused (pandas would match missing symbols).
fn gene_list(path: &str) -> Result<std::collections::HashSet<String>> {
    let options = ReadOptions {
        header: Header::Names(vec!["SYMBOL".to_owned()]),
        comment: Some(b'#'),
        ..ReadOptions::default()
    };
    let list = read_table_with(Path::new(path), &options)?;
    let column = list.column("SYMBOL")?;
    let mut genes = std::collections::HashSet::new();
    for row in 0..list.rows {
        match column.cell(row) {
            Cell::Text(text) => {
                genes.insert(text.to_owned());
            }
            Cell::Missing => {
                return Err(refuse(format!("a missing value in the gene list {path}")));
            }
            // A number or boolean never equals a symbol text.
            _ => {}
        }
    }
    Ok(genes)
}

/// `df["SYMBOL"].isin(genes)`.
fn symbol_mask(table: &Table, genes: &std::collections::HashSet<String>) -> Result<Vec<bool>> {
    let symbols = table.column("SYMBOL")?;
    Ok((0..table.rows)
        .map(|row| matches!(symbols.cell(row), Cell::Text(text) if genes.contains(text)))
        .collect())
}

/// `Series.clip(lower, upper)`: NaN stays; a value that is not `>= lower` becomes `lower`,
/// then one that is not `<= upper` becomes `upper`.
fn clip(value: f64, lower: f64, upper: f64) -> f64 {
    if value.is_nan() {
        return value;
    }
    let value = if value >= lower { value } else { lower };
    if value <= upper { value } else { upper }
}

/// `np.power(series, shape)`: pandas hands a Series power to `operator.pow`, where numpy
/// 1.26.4 uses `sqrt`, the value itself or a square for the exponents 0.5, 1 and 2 (the
/// finding of `crate::acmg_scoring::ranking_pext`, checked there on 807,012 values).
fn series_power(base: f64, shape: f64) -> f64 {
    if shape == 1.0 {
        base
    } else if shape == 0.5 {
        base.sqrt()
    } else if shape == 2.0 {
        base * base
    } else {
        base.powf(std::hint::black_box(shape))
    }
}

/// `df[column].map(coerce_pext_value)` through the earlier scalar module.
fn pext_values(table: &Table, name: &str) -> Result<Vec<f64>> {
    let column = table.column(name)?;
    let value = |row: usize| -> Result<PextValue> {
        Ok(match column {
            Column::Bool(values) => PextValue::NumpyBool(values[row]),
            Column::NullableInt(_) | Column::NullableBool(_) => {
                return Err(refuse(format!(
                    "the pext column {name} is {}",
                    column.dtype_name()
                )));
            }
            Column::Float(values) => PextValue::Primitive(MaskCell::Float(values[row])),
            Column::Object(values) => match &values[row] {
                ObjectValue::None => PextValue::Primitive(MaskCell::Null),
                ObjectValue::Nan => PextValue::Primitive(MaskCell::Float(f64::NAN)),
                ObjectValue::Text(text) => PextValue::Primitive(MaskCell::Text(text.clone())),
                ObjectValue::Int(number) => PextValue::Primitive(MaskCell::Integer(*number)),
                ObjectValue::Float(number) => PextValue::Primitive(MaskCell::Float(*number)),
                ObjectValue::Bool(flag) => PextValue::Primitive(MaskCell::Bool(*flag)),
            },
            _ => match column.cell(row) {
                Cell::Missing => PextValue::Primitive(MaskCell::Float(f64::NAN)),
                Cell::Text(text) => PextValue::Primitive(MaskCell::Text(text.to_owned())),
                Cell::Int(number) => PextValue::Primitive(MaskCell::Integer(number)),
                Cell::Float(number) => PextValue::Primitive(MaskCell::Float(number)),
                Cell::Bool(flag) => PextValue::Primitive(MaskCell::Bool(flag)),
            },
        })
    };
    (0..table.rows)
        .map(|row| coerce_pext_value(&value(row)?).map_err(source_error))
        .collect()
}

/// One key column of the final sort and of the variant groups.
enum KeyValues {
    /// NaN is a missing value.
    Numbers(Vec<f64>),
    Texts(Vec<Option<String>>),
}

impl KeyValues {
    fn of(table: &Table, name: &str) -> Result<KeyValues> {
        Ok(match table.column(name)? {
            Column::Float(values) => KeyValues::Numbers(values.clone()),
            Column::Int(values) => {
                if values.iter().any(|value| value.unsigned_abs() > 1 << 53) {
                    return Err(refuse(format!("{name} holds an integer past 2**53")));
                }
                KeyValues::Numbers(values.iter().map(|value| *value as f64).collect())
            }
            Column::Bool(values) => KeyValues::Numbers(
                values
                    .iter()
                    .map(|value| f64::from(u8::from(*value)))
                    .collect(),
            ),
            Column::Text(values) => KeyValues::Texts(values.clone()),
            other => {
                return Err(refuse(format!(
                    "the sort key {name} is a {} column of mixed or nullable values",
                    other.dtype_name()
                )));
            }
        })
    }

    fn missing(&self, row: usize) -> bool {
        match self {
            KeyValues::Numbers(values) => values[row].is_nan(),
            KeyValues::Texts(values) => values[row].is_none(),
        }
    }

    /// Ascending order of two present values (-0.0 equals 0.0, as in pandas' categories).
    fn compare(&self, left: usize, right: usize) -> Ordering {
        match self {
            KeyValues::Numbers(values) => values[left]
                .partial_cmp(&values[right])
                .unwrap_or(Ordering::Equal),
            KeyValues::Texts(values) => values[left].cmp(&values[right]),
        }
    }

    /// The value as a group key (`-0.0` and `0.0` together).
    fn group_key(&self, row: usize) -> String {
        match self {
            KeyValues::Numbers(values) => {
                let value = if values[row] == 0.0 { 0.0 } else { values[row] };
                format!("n{:016x}", value.to_bits())
            }
            KeyValues::Texts(values) => format!("s{}", values[row].as_deref().unwrap_or("")),
        }
    }
}

/// `sort_and_rank_variants(df, ped_df, fam_name, ...)`: the protein-coding rows of the
/// table with the ranking columns, sorted by the best `sort_index` of each variant, then
/// `ACMG_quant_score`, `chrom`, `pos`, `ref`, `alt`, and ranked by `variant_id`.
pub(crate) fn sort_and_rank_variants(
    table: &Table,
    ped: Option<&Table>,
    fam_name: Option<&str>,
    options: &RankingOptions<'_>,
) -> Result<Table> {
    // R1: drop the stale columns, copy the score, keep the protein-coding rows.
    let mut df = table.clone();
    for stale in ["haplo_insuf_index", "zygosity_inheritance_mechanism_index"] {
        if let Some(index) = df.column_index(stale) {
            df.names.remove(index);
            df.columns.remove(index);
        }
    }
    let quant = df.column("ACMG_quant_score")?.clone();
    df.set_column("sort_index", quant);
    let biotype = df.column("BIOTYPE")?;
    let keep: Vec<bool> = (0..df.rows)
        .map(|row| biotype.cell(row) == Cell::Text("protein_coding"))
        .collect();
    let mut df = df.filter_rows(&keep);
    let rows = df.rows;

    // R2: the proband from the pedigree.
    let (proband, male, female, proband_het) = pedigree_context(&df, ped, fam_name)?;

    // Control-common penalty (R3 multiplies it in).
    compute_control_common(&mut df, proband.as_deref())?;
    let mut sort_index = float_values(&df, "sort_index")?;
    let control_index = float_values(&df, "control_common_index")?;
    for row in 0..rows {
        sort_index[row] *= control_index[row];
    }
    df.set_column("sort_index", Column::Float(sort_index.clone()));

    // R4: the zygosity, inheritance and mechanism compatibility penalty.
    let masks = variant_mechanism_masks(&df)?;
    let x_locus: Vec<bool> = match df.column("chrom") {
        Ok(chrom) => (0..rows)
            .map(|row| {
                let text = python_str(chrom.cell(row))
                    .unwrap_or_default()
                    .to_uppercase();
                text.strip_prefix("CHR").unwrap_or(&text) == "X"
            })
            .collect(),
        Err(_) => vec![false; rows],
    };
    let mask = |name: &str| masks.get(name);
    let (dom_lof, dom_unresolved) = (
        mask("has_dom_lof_history")?,
        mask("has_dom_unresolved_history")?,
    );
    let (xr_lof, xr_unresolved) = (
        mask("has_x_recessive_lof_history")?,
        mask("has_x_recessive_unresolved_history")?,
    );
    let (xd_lof, xd_unresolved) = (
        mask("has_x_dominant_lof_history")?,
        mask("has_x_dominant_unresolved_history")?,
    );
    let exact_gof = mask("is_exact_gof")?;
    let predicted_lof = mask("is_predicted_lof")?;
    let quant_values = float_values(&df, "ACMG_quant_score")?;
    let mut penalty_rows = vec![false; rows];
    let mut high_acmg_lof = vec![false; rows];
    for row in 0..rows {
        let compatible = dom_lof[row]
            || dom_unresolved[row]
            || (x_locus[row]
                && male
                && (xr_lof[row] || xr_unresolved[row] || xd_lof[row] || xd_unresolved[row]))
            || (x_locus[row] && female && (xd_lof[row] || xd_unresolved[row]));
        let incompatible = predicted_lof[row] && !exact_gof[row] && !compatible;
        high_acmg_lof[row] = quant_values[row] > 0.89;
        penalty_rows[row] = high_acmg_lof[row] && proband_het[row] && incompatible;
    }
    let compatibility: Vec<f64> = penalty_rows
        .iter()
        .map(|penalty| if *penalty { 0.9 } else { 1.0 })
        .collect();
    df.set_column(
        "zygosity_inheritance_mechanism_compatibility",
        Column::Float(compatibility.clone()),
    );

    // R5 and R6: the dispensable and relevant gene lists.
    for (list, column, weight) in [
        (options.dispensable_gene_list, "dispensable_gene_index", 0.9),
        (options.relevant_gene_list, "relevant_gene_index", 1.2),
    ] {
        let Some(path) = list else {
            continue;
        };
        let genes = gene_list(path)?;
        let listed = symbol_mask(&df, &genes)?;
        let index: Vec<f64> = listed
            .iter()
            .map(|listed| if *listed { weight } else { 1.0 })
            .collect();
        df.set_column(column, Column::Float(index.clone()));
        for row in 0..rows {
            sort_index[row] *= index[row];
        }
        df.set_column("sort_index", Column::Float(sort_index.clone()));
    }
    for row in 0..rows {
        if penalty_rows[row] {
            sort_index[row] *= compatibility[row];
        }
    }
    df.set_column("sort_index", Column::Float(sort_index.clone()));

    // R7: the pext modulation of high-scoring rows.
    let mut pext_column = None;
    let tissues: Vec<&str> = options
        .pext_tissues
        .split(',')
        .map(python_strip)
        .filter(|tissue| !tissue.is_empty())
        .collect();
    if let Some(first) = tissues.first() {
        let name = format!("PEXT_{}", first.replace('-', "_"));
        if df.has_column(&name) {
            pext_column = Some(name);
        }
    }
    if pext_column.is_none() && df.has_column("PEXT_MEAN") {
        pext_column = Some("PEXT_MEAN".to_owned());
    }
    if let Some(name) = pext_column {
        let cutoff = options.pext_low_expression_cutoff;
        let floor = options.pext_penalty_floor;
        let shape = options.pext_penalty_shape;
        if cutoff <= 0.0 {
            return Err(python_error(
                "ValueError",
                "pext_low_expression_cutoff must be > 0",
            ));
        }
        if !(0.0 <= floor && floor <= 1.0) {
            return Err(python_error(
                "ValueError",
                "pext_penalty_floor must be between 0 and 1",
            ));
        }
        if shape <= 0.0 {
            return Err(python_error("ValueError", "pext_penalty_shape must be > 0"));
        }
        let pext: Vec<f64> = pext_values(&df, &name)?
            .into_iter()
            .map(|value| clip(if value.is_nan() { 1.0 } else { value }, 0.0, 1.0))
            .collect();
        let pext_index: Vec<f64> = pext
            .iter()
            .map(|value| {
                let scaled = clip(value / cutoff, 0.0, 1.0);
                let modulated = floor + (1.0 - floor) * series_power(scaled, shape);
                if *value >= cutoff { 1.0 } else { modulated }
            })
            .collect();
        // On no rows, pandas' `map` of an object column gives an empty object Series, and
        // the arithmetic and `np.where` after it keep that type.
        let numeric = matches!(
            df.column(&name)?,
            Column::Float(_) | Column::Int(_) | Column::Bool(_)
        );
        let pext_column = if rows == 0 && !numeric {
            Column::Object(Vec::new())
        } else {
            Column::Float(pext_index.clone())
        };
        df.set_column("pext_sort_index", pext_column);
        for row in 0..rows {
            if high_acmg_lof[row] {
                sort_index[row] *= pext_index[row];
            }
        }
        df.set_column("sort_index", Column::Float(sort_index.clone()));
    }

    // R8: the best sort_index of each variant (groups by chrom, pos, ref, alt; a missing
    // key leaves NaN), the sort, and the rank by first appearance of variant_id.
    let keys: Vec<KeyValues> = ["chrom", "pos", "ref", "alt"]
        .iter()
        .map(|name| KeyValues::of(&df, name))
        .collect::<Result<_>>()?;
    let mut best: HashMap<Vec<String>, (usize, f64)> = HashMap::new();
    let mut group_of: Vec<Option<Vec<String>>> = Vec::with_capacity(rows);
    for row in 0..rows {
        if keys.iter().any(|key| key.missing(row)) {
            group_of.push(None);
            continue;
        }
        let key: Vec<String> = keys.iter().map(|key| key.group_key(row)).collect();
        // pandas' group max: the first value that is greater than every earlier one.
        let entry = best.entry(key.clone()).or_insert((0, f64::NEG_INFINITY));
        let value = sort_index[row];
        if !value.is_nan() {
            entry.0 += 1;
            if value > entry.1 {
                entry.1 = value;
            }
        }
        group_of.push(Some(key));
    }
    let max_score: Vec<f64> = group_of
        .iter()
        .map(|key| match key {
            Some(key) => match best[key] {
                (0, _) => f64::NAN,
                (_, value) => value,
            },
            None => f64::NAN,
        })
        .collect();
    let sort_keys = [
        (KeyValues::Numbers(max_score), false),
        (KeyValues::Numbers(quant_values), false),
    ];
    let mut order: Vec<usize> = (0..rows).collect();
    order.sort_by(|left, right| {
        let all = sort_keys
            .iter()
            .map(|(key, ascending)| (key, *ascending))
            .chain(keys.iter().map(|key| (key, true)));
        for (key, ascending) in all {
            let ordering = match (key.missing(*left), key.missing(*right)) {
                (true, true) => Ordering::Equal,
                (true, false) => Ordering::Greater,
                (false, true) => Ordering::Less,
                (false, false) => {
                    let ordering = key.compare(*left, *right);
                    if ascending {
                        ordering
                    } else {
                        ordering.reverse()
                    }
                }
            };
            if ordering != Ordering::Equal {
                return ordering;
            }
        }
        Ordering::Equal
    });
    if let Some(index) = df.column_index("max_variant_score") {
        df.names.remove(index);
        df.columns.remove(index);
    }
    let mut sorted = df.take_rows(&order);
    let identities: Vec<Option<String>> = match sorted.column("variant_id") {
        Ok(Column::Text(values)) => values.clone(),
        Ok(other) => {
            return Err(refuse(format!(
                "variant_id is a {} column",
                other.dtype_name()
            )));
        }
        Err(_) => {
            let parts: Vec<&Column> = ["chrom", "pos", "ref", "alt"]
                .iter()
                .map(|name| sorted.column(name))
                .collect::<Result<_>>()?;
            (0..sorted.rows)
                .map(|row| {
                    Some(
                        parts
                            .iter()
                            .map(|column| {
                                python_str(column.cell(row)).unwrap_or_else(|| "nan".to_owned())
                            })
                            .collect::<Vec<_>>()
                            .join(":"),
                    )
                })
                .collect()
        }
    };
    let mut codes: HashMap<String, i64> = HashMap::new();
    let ranks: Vec<i64> = identities
        .iter()
        .map(|identity| match identity {
            None => 0,
            Some(identity) => {
                let next = codes.len() as i64 + 1;
                *codes.entry(identity.clone()).or_insert(next)
            }
        })
        .collect();
    sorted.set_column("variant_rank", Column::Int(ranks));
    Ok(sorted)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table(columns: Vec<(&str, Column)>) -> Table {
        let rows = columns.first().map_or(0, |(_, column)| column.len());
        let (names, columns) = columns
            .into_iter()
            .map(|(name, column)| (name.to_owned(), column))
            .unzip();
        Table {
            names,
            columns,
            rows,
        }
    }

    #[test]
    fn summary_without_clingen_keeps_int_columns() {
        let mut input = table(vec![(
            "variant_id",
            Column::Text(vec![Some("1:1:A-T".into()), Some("1:2:C-G".into())]),
        )]);
        let criteria = vec![
            ("PVS1", Column::Int(vec![4, 0])),
            ("PM2", Column::Int(vec![1, 0])),
            ("BP7", Column::Bool(vec![false, false])),
        ];
        let matrix = summarize_acmg_criteria(&mut input, &criteria, None, false).unwrap();
        assert_eq!(matrix.names.len(), 29);
        assert_eq!(matrix.column("PS4").unwrap(), &Column::Int(vec![0, 0]));
        assert_eq!(
            input.column("ACMG_criteria").unwrap(),
            &Column::Text(vec![Some("PVS1;PM2_Supporting".into()), None])
        );
        let bad = vec![("PS1", Column::Int(vec![6, 0]))];
        let error = summarize_acmg_criteria(&mut input, &bad, None, false)
            .unwrap_err()
            .to_string();
        assert_eq!(error, "ValueError: Column PS1 has non-integer values");
        let error = summarize_acmg_criteria(&mut input, &criteria, Some("/no/such.pkl"), true)
            .unwrap_err()
            .to_string();
        assert!(error.starts_with(
            "FileNotFoundError: ClinGen override is enabled but the evidence map is missing: '/no/such.pkl'."
        ));
    }

    #[test]
    fn posterior_columns_follow_the_matrix() {
        let matrix = table(vec![
            ("PVS1", Column::Int(vec![4, 0])),
            ("BA1", Column::Int(vec![0, 5])),
            (
                "variant_id",
                Column::Text(vec![Some("a".into()), Some("b".into())]),
            ),
        ]);
        let (scores, classes) = posterior_columns(&matrix).unwrap();
        let Column::Float(scores) = scores else {
            panic!("float scores")
        };
        assert_eq!(scores[1], 0.0);
        assert_eq!(
            classes,
            Column::Text(vec![
                Some("Likely Pathogenic".into()),
                Some("Benign".into())
            ])
        );
    }

    #[test]
    fn control_common_without_control_columns() {
        let mut input = table(vec![("x", Column::Int(vec![1]))]);
        compute_control_common(&mut input, None).unwrap();
        assert_eq!(
            input.column("control_common").unwrap(),
            &Column::Bool(vec![false])
        );
        assert_eq!(
            input.column("control_common_index").unwrap(),
            &Column::Float(vec![1.0])
        );
    }

    fn texts(values: &[&str]) -> Column {
        Column::Text(
            values
                .iter()
                .map(|value| Some((*value).to_owned()))
                .collect(),
        )
    }

    #[test]
    fn ranking_sorts_by_the_best_score_of_each_variant() {
        let input = table(vec![
            ("chrom", texts(&["chr1", "chr1", "chr2", "chrX"])),
            ("pos", Column::Int(vec![5, 5, 3, 9])),
            ("ref", texts(&["A", "A", "C", "G"])),
            ("alt", texts(&["T", "T", "G", "A"])),
            (
                "BIOTYPE",
                texts(&[
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                    "lncRNA",
                ]),
            ),
            ("SYMBOL", texts(&["G1", "G1", "G2", "G3"])),
            (
                "ACMG_quant_score",
                Column::Float(vec![0.2, 0.95, 0.5, 0.99]),
            ),
            (
                "var_plausible_patho_mechs",
                texts(&["", "", "dominant_LOF", ""]),
            ),
            (
                "variant_effect",
                texts(&["uncertain", "uncertain", "uncertain", "uncertain"]),
            ),
            (
                "variant_id",
                texts(&["1:5:A-T", "1:5:A-T", "2:3:C-G", "X:9:G-A"]),
            ),
        ]);
        let options = RankingOptions {
            pext_tissues: "",
            pext_low_expression_cutoff: 0.1,
            pext_penalty_floor: 0.8,
            pext_penalty_shape: 0.5,
            relevant_gene_list: None,
            dispensable_gene_list: None,
        };
        let ranked = sort_and_rank_variants(&input, None, None, &options).unwrap();
        assert_eq!(ranked.rows, 3);
        assert_eq!(
            ranked.column("ACMG_quant_score").unwrap(),
            &Column::Float(vec![0.95, 0.2, 0.5])
        );
        assert_eq!(
            ranked.column("variant_rank").unwrap(),
            &Column::Int(vec![1, 1, 2])
        );
        assert_eq!(
            ranked.column("control_common_index").unwrap(),
            &Column::Float(vec![1.0; 3])
        );
        assert!(!ranked.has_column("max_variant_score"));
        let error = sort_and_rank_variants(&input, Some(&input), Some("F"), &options)
            .unwrap_err()
            .to_string();
        assert_eq!(error, "KeyError: '#FamilyID'");
    }
}
