//! `scripts/acmg_ps3_bs3_pp5_bp6_reported.py` (258 lines): PS3, BS3, PP5 and BP6 from ClinVar and
//! MaveDB.
//!
//! `ps3_bs3_criteria` is `PS3_BS3_criteria` with `mavedb_score_interpretation` and
//! `mavedb_interpretation_per_row`; `pp5_bp6_criteria` is `PP5_BP6_criteria`.
//!
//! With a MaveDB metadata table, Python runs `df.apply(mavedb_interpretation_per_row, axis=1)`
//! and the step goes on with the frame `apply` builds. That frame has three more columns
//! (`MaveDB_PS3`, `MaveDB_BS3`, `MaveDB_score_interpretation`), and every column has been
//! through `infer_objects` (the rows are handed out as `object` Series): an `object` column of
//! booleans without a missing value becomes `bool`, an all-missing interpretation column
//! becomes `float64`, as does any `object` column of only NaN; numbers, booleans and text
//! columns keep their type. `StepOutput.columns`
//! holds the columns whose type changes and the three new ones, to be set in that order.
//! Tables for which `apply` would turn numbers into other numbers (no `object` column, so the
//! rows are float Series), mixed `object` columns (which may hold `None`, which pandas keeps
//! apart from NaN), the nullable `Int64`/`boolean` columns and empty tables are refused.
use super::acmg_pp3_bp4_bp7_insilico::harness_options;
use super::acmg_pp3_bp4_bp7_insilico::pandas_masks::{
    Mask, high_confidence_review, python_str, strengths, text_contains,
};
use super::mavedb_interpreter::{Metadata, pick_interpreter_from_metadata, python_float_of_text};
use super::types::{AcmgError, HarnessCall, StepOutput, not_converted};
use crate::cli::{Result, invalid};
use crate::pandas_polars;
use crate::pandas_table::{Cell, Column, Table, read_table};
use std::path::Path;

fn index_error() -> Box<dyn std::error::Error> {
    Box::new(AcmgError::new("IndexError", "list index out of range"))
}

fn key_error(name: &str) -> Box<dyn std::error::Error> {
    Box::new(AcmgError::new("KeyError", format!("'{name}'")))
}

/// `str(row.get(name, default)).split('&')`.
fn split_field(table: &Table, name: &str, row: usize, default: &str) -> Vec<String> {
    let text = match table.column_index(name) {
        Some(index) => python_str(table.columns[index].cell(row)),
        None => default.to_owned(),
    };
    text.split('&').map(str::to_owned).collect()
}

/// One row of `mavedb_interpretation_per_row`: `(MaveDB_PS3, MaveDB_BS3, interpretation)`,
/// the interpretation `None` for NaN.
fn interpret_row(
    table: &Table,
    row: usize,
    metadata: &Metadata<'_>,
) -> Result<(bool, bool, Option<String>)> {
    let urn_cell = match table.column_index("MaveDB_urn") {
        Some(index) => table.columns[index].cell(row),
        None => Cell::Missing,
    };
    if urn_cell == Cell::Missing {
        return Ok((false, false, None));
    }
    let high_confs = split_field(table, "MaveDB_high_conf", row, "");
    let pvalues = split_field(table, "MaveDB_pvalue", row, "");
    let scores = split_field(table, "MaveDB_score", row, "");
    let urns: Vec<String> = python_str(urn_cell).split('&').map(str::to_owned).collect();
    let mut interpretations = Vec::new();
    // The assay family of each piece of evidence, in the order found.
    let mut ps3_evidence: Vec<String> = Vec::new();
    let mut bs3_evidence: Vec<String> = Vec::new();
    for (position, urn) in urns.iter().enumerate() {
        let Some(meta_row) = metadata.row(urn) else {
            continue;
        };
        let interpretation = metadata
            .get(meta_row, "rationale")
            .ok_or_else(|| key_error("rationale"))?;
        let family = metadata.get(meta_row, "assay_family");
        interpretations.push(python_str(interpretation));
        let score = scores.get(position).ok_or_else(index_error)?;
        high_confs.get(position).ok_or_else(index_error)?;
        let pvalue = pvalues.get(position).ok_or_else(index_error)?;
        let Some(score) = python_float_of_text(score)? else {
            continue;
        };
        let pvalue = python_float_of_text(pvalue)?.unwrap_or(f64::NAN);
        let family_key = || -> Result<String> {
            match family {
                None => Ok("Unknown".to_owned()),
                Some(Cell::Text(text)) => Ok(text.to_owned()),
                Some(other) => Err(invalid(format!(
                    "MaveDB metadata {urn}: assay_family {other:?} is not a text; \
                     grouping evidence by it is not reproduced"
                ))),
            }
        };
        if pvalue.is_nan() {
            let result = pick_interpreter_from_metadata(urn, score, metadata)?;
            if result.mavedb_ps3 {
                ps3_evidence.push(family_key()?);
            }
            if result.mavedb_bs3 {
                bs3_evidence.push(family_key()?);
            }
        } else if pvalue < 0.05 {
            ps3_evidence.push(family_key()?);
        } else if pvalue >= 0.05 {
            bs3_evidence.push(family_key()?);
        }
    }
    // Per-family resolution: a family with both kinds of evidence counts for the side with
    // more pieces of evidence (a tie counts for neither).
    let count = |evidence: &[String], family: &str| {
        evidence
            .iter()
            .filter(|item| item.as_str() == family)
            .count()
    };
    let decided = |own: &[String], other: &[String]| {
        own.iter()
            .any(|family| count(other, family) == 0 || count(own, family) > count(other, family))
    };
    let mavedb_ps3 = decided(&ps3_evidence, &bs3_evidence);
    let mavedb_bs3 = decided(&bs3_evidence, &ps3_evidence);
    Ok((mavedb_ps3, mavedb_bs3, Some(interpretations.join("&"))))
}

/// What `infer_objects` makes of a column after `df.apply(..., axis=1)`: `Ok(None)` when
/// its type stays.
fn after_apply(name: &str, column: &Column) -> Result<Option<Column>> {
    match column {
        Column::Int(_) | Column::Float(_) | Column::Bool(_) => Ok(None),
        // An object column of only NaN (the table's missing text values are NaN, not None)
        // becomes float64.
        Column::Text(values) if !values.is_empty() && values.iter().all(Option::is_none) => {
            Ok(Some(Column::Float(vec![f64::NAN; values.len()])))
        }
        Column::Text(_) => Ok(None),
        Column::BoolOrMissing(values) => {
            if values.iter().all(Option::is_some) {
                Ok(Some(Column::Bool(
                    values.iter().map(|value| value.unwrap_or(false)).collect(),
                )))
            } else if values.iter().all(Option::is_none) {
                Ok(Some(Column::Float(vec![f64::NAN; values.len()])))
            } else {
                Ok(None)
            }
        }
        other => Err(invalid(format!(
            "{name} is a pandas {} column ({}); df.apply(axis=1) is not reproduced for it",
            other.dtype_name(),
            match other {
                Column::Object(_) => "mixed values",
                _ => "nullable",
            }
        ))),
    }
}

/// `mavedb_score_interpretation(df, mavedb_metadata)`: the columns the returned frame has
/// that differ from `table` (re-typed existing columns, then the three MaveDB columns).
fn mavedb_score_interpretation(
    table: &Table,
    metadata_table: &Table,
) -> Result<Vec<(String, Column)>> {
    if table.rows == 0 {
        return Err(invalid(
            "df.apply(axis=1) on an empty table returns the table unchanged in pandas; not reproduced",
        ));
    }
    if !table
        .columns
        .iter()
        .any(|column| column.dtype_name() == "object")
    {
        return Err(invalid(
            "a table without an object column is handed to apply as float rows; not reproduced",
        ));
    }
    let mut changed = Vec::new();
    for (name, column) in table.names.iter().zip(&table.columns) {
        if let Some(column) = after_apply(name, column)? {
            changed.push((name.clone(), column));
        }
    }
    let metadata = Metadata::new(metadata_table);
    let mut ps3 = Vec::with_capacity(table.rows);
    let mut bs3 = Vec::with_capacity(table.rows);
    let mut interpretation = Vec::with_capacity(table.rows);
    for row in 0..table.rows {
        let (row_ps3, row_bs3, row_interpretation) = interpret_row(table, row, &metadata)?;
        ps3.push(row_ps3);
        bs3.push(row_bs3);
        interpretation.push(row_interpretation);
    }
    let interpretation = if interpretation.iter().all(Option::is_none) {
        Column::Float(vec![f64::NAN; table.rows])
    } else {
        Column::Text(interpretation)
    };
    let mut new_columns = vec![
        ("MaveDB_PS3".to_owned(), Column::Bool(ps3)),
        ("MaveDB_BS3".to_owned(), Column::Bool(bs3)),
        ("MaveDB_score_interpretation".to_owned(), interpretation),
    ];
    // An existing MaveDB column is overwritten in place; it is not re-typed separately.
    changed.retain(|(name, _)| !new_columns.iter().any(|(new, _)| new == name));
    changed.append(&mut new_columns);
    // The log line selects these columns of the new frame.
    for name in ["MaveDB_pvalue", "MaveDB_score", "MaveDB_high_conf"] {
        if !table.has_column(name) {
            return Err(Box::new(AcmgError::new(
                "KeyError",
                format!("\"['{name}'] not in index\""),
            )));
        }
    }
    Ok(changed)
}

/// Reads the MaveDB metadata as `pd.read_table(path, low_memory=False)` and keeps the first row
/// per `urn` (`drop_duplicates(subset=["urn"])`).
pub(crate) fn read_mavedb_metadata(path: &str) -> Result<Table> {
    let table = read_table(Path::new(path))?;
    if !table.has_column("urn") {
        return Err(key_error("urn"));
    }
    let rows = pandas_polars::drop_duplicates_first(&table, Some(&["urn"]))?;
    Ok(table.take_rows(&rows))
}

/// `PS3_BS3_criteria(df, mavedb_metadata_tsv)`: PS3 and BS3 (strength 3), `clinvar_patho` and
/// `clinvar_benign` (passed on to PS1/PM5, PM1 and PP5/BP6) and, with a metadata table, the
/// columns of the frame Python goes on with.
pub(crate) fn ps3_bs3_criteria(
    table: &Table,
    mavedb_metadata_tsv: Option<&str>,
) -> Result<StepOutput> {
    let pathogenic = text_contains(table, "CLNSIG", Some(""), &["Pathogenic"])?.remove(0);
    let high_confidence = Mask::Plain(high_confidence_review(table, "CLNREVSTAT")?);
    let clinvar_patho = pathogenic.and(&high_confidence);
    let benign = text_contains(table, "CLNSIG", Some(""), &["Benign"])?.remove(0);
    let clinvar_benign = benign.and(&high_confidence);
    let (columns, ps3, bs3) = match mavedb_metadata_tsv.filter(|path| !path.is_empty()) {
        Some(path) => {
            let metadata = read_mavedb_metadata(path)?;
            let columns = mavedb_score_interpretation(table, &metadata)?;
            let mavedb = |name: &str| match columns.iter().find(|(column, _)| column == name) {
                Some((_, Column::Bool(values))) => Mask::Plain(values.clone()),
                _ => unreachable!("the MaveDB columns are boolean"),
            };
            let ps3 = clinvar_patho.or(&mavedb("MaveDB_PS3"));
            let bs3 = clinvar_benign.or(&mavedb("MaveDB_BS3"));
            (columns, ps3, bs3)
        }
        None => (Vec::new(), clinvar_patho.clone(), clinvar_benign.clone()),
    };
    Ok(StepOutput {
        criteria: vec![
            ("PS3", Column::Int(strengths(&ps3, 3)?)),
            ("BS3", Column::Int(strengths(&bs3, 3)?)),
        ],
        columns,
        passed_on: vec![
            (
                "clinvar_patho".to_owned(),
                Column::Bool(clinvar_patho.plain()?.to_vec()),
            ),
            (
                "clinvar_benign".to_owned(),
                Column::Bool(clinvar_benign.plain()?.to_vec()),
            ),
        ],
    })
}

/// `PP5_BP6_criteria(df, clinvar_patho, clinvar_benign)`: ClinVar's (likely) pathogenic or
/// benign classifications not already counted by PS3/BS3 (strength 1).
pub(crate) fn pp5_bp6_criteria(
    table: &Table,
    clinvar_patho: &[bool],
    clinvar_benign: &[bool],
) -> Result<StepOutput> {
    if clinvar_patho.len() != table.rows || clinvar_benign.len() != table.rows {
        return Err(invalid(
            "clinvar_patho and clinvar_benign must cover every row",
        ));
    }
    let mut found = text_contains(table, "CLNSIG", Some(""), &["athogenic", "enign"])?;
    let benign = found.pop().expect("two patterns");
    let pathogenic = found.pop().expect("two patterns");
    let not = |values: &[bool]| Mask::Plain(values.iter().map(|value| !value).collect());
    let pp5 = pathogenic.and(&not(clinvar_patho));
    let bp6 = benign.and(&not(clinvar_benign));
    Ok(StepOutput {
        criteria: vec![
            ("PP5", Column::Int(strengths(&pp5, 1)?)),
            ("BP6", Column::Int(strengths(&bp6, 1)?)),
        ],
        ..StepOutput::default()
    })
}

/// `priva acmg-step acmg_ps3_bs3_pp5_bp6_reported <function>`:
///
/// - `PS3_BS3_criteria`, option `mavedb_metadata_tsv` (text or null): the frame Python returns
///   with `__PS3`, `__BS3`, `__clinvar_patho`, `__clinvar_benign` appended;
/// - `PP5_BP6_criteria`, options `clinvar_patho` and `clinvar_benign` (lists of booleans): a
///   frame of `PP5`, `BP6`.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let table = &call.table;
    match function {
        "PS3_BS3_criteria" => {
            let path = harness_options::text(&call.options, "mavedb_metadata_tsv")?;
            let output = ps3_bs3_criteria(table, path)?;
            let mut frame = table.clone();
            for (name, column) in &output.columns {
                frame.set_column(name, column.clone());
            }
            for (name, column) in output.criteria.iter().map(|(name, values)| (*name, values)) {
                frame.set_column(&format!("__{name}"), column.clone());
            }
            for (name, column) in &output.passed_on {
                frame.set_column(&format!("__{name}"), column.clone());
            }
            Ok(frame)
        }
        "PP5_BP6_criteria" => {
            let patho = harness_options::bool_list(&call.options, "clinvar_patho")?
                .ok_or_else(|| invalid("option clinvar_patho is required"))?;
            let benign = harness_options::bool_list(&call.options, "clinvar_benign")?
                .ok_or_else(|| invalid("option clinvar_benign is required"))?;
            let output = pp5_bp6_criteria(table, &patho, &benign)?;
            Ok(harness_options::criteria_frame(&output, table.rows))
        }
        _ => Err(not_converted("acmg_ps3_bs3_pp5_bp6_reported", function)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn text(values: &[Option<&str>]) -> Column {
        Column::Text(
            values
                .iter()
                .map(|value| value.map(str::to_owned))
                .collect(),
        )
    }

    #[test]
    fn clinvar_masks_need_two_stars() {
        let table = Table {
            names: vec!["CLNSIG".into(), "CLNREVSTAT".into()],
            columns: vec![
                text(&[
                    Some("Pathogenic"),
                    Some("Likely_pathogenic"),
                    Some("Benign"),
                    None,
                ]),
                text(&[
                    Some("reviewed_by_expert_panel"),
                    Some("practice_guideline"),
                    Some("criteria_provided,_single_submitter"),
                    None,
                ]),
            ],
            rows: 4,
        };
        let output = ps3_bs3_criteria(&table, None).unwrap();
        assert_eq!(
            output.criterion("PS3"),
            Some(&Column::Int(vec![3, 0, 0, 0]))
        );
        assert_eq!(
            output.criterion("BS3"),
            Some(&Column::Int(vec![0, 0, 0, 0]))
        );
        let patho = match output.passed_on("clinvar_patho") {
            Some(Column::Bool(values)) => values.clone(),
            other => panic!("{other:?}"),
        };
        let pp5 = pp5_bp6_criteria(&table, &patho, &[false; 4]).unwrap();
        assert_eq!(pp5.criterion("PP5"), Some(&Column::Int(vec![0, 1, 0, 0])));
        assert_eq!(pp5.criterion("BP6"), Some(&Column::Int(vec![0, 0, 1, 0])));
    }
}
