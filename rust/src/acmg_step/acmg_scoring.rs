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
//! counts are tested stops at their placeholder).
//!
//! Not converted yet: `sort_and_rank_variants` (stages R1 to R8 around the existing
//! `crate::acmg_scoring` ranking modules).
use super::stats;
use super::types::{AcmgError, CRITERIA_ORDER, HarnessCall, not_converted};
use crate::acmg_scoring::SourceError;
use crate::acmg_scoring::calculate_posterior_probability::{
    Probability, calculate_posterior_probability,
};
use crate::acmg_scoring::create_criteria_summary::create_criteria_summary;
use crate::acmg_scoring::translate_strength_level::translate_strength_level;
use crate::cli::{Result, invalid};
use crate::json_text::Json;
use crate::pandas_table::{Cell, Column, ObjectValue, Table};
use crate::pickle_reader::{self, Value};
use crate::python_text::python_strip;
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
}
