//! Explicit table boundary for ALREADY assigned criterion strengths.
//!
//! Keeps summarize_acmg_criteria's order: validate matrix → compute all summary
//! rows → rebuild overridden rows criterion-by-criterion → calculate posteriors.
//! Pickle loading and variant-to-ClinGen mapping remain upstream Python work.
use super::{
    CRITERIA_ORDER,
    calculate_posterior_probability::{Probability, calculate_posterior_probability},
    create_criteria_summary::create_criteria_summary,
    python_strip,
    translate_strength_level::translate_strength_level,
};
use crate::cli::{Result, invalid, required};
use polars_core::prelude::*;
use std::{collections::HashSet, fs::OpenOptions};

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let evidence_name = match args {
        [_, _, mode] if mode == "--no-clingen-override" => None,
        [_, _, mode, name] if mode == "--clingen-evidence-column" => Some(name.as_str()),
        _ => {
            return Err(invalid(
                "score-criteria requires INPUT.tsv NEW_OUTPUT.tsv and either --no-clingen-override or --clingen-evidence-column COLUMN",
            ));
        }
    };
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&args[0])?;
    let headers = reader.headers()?.clone();
    let unique: HashSet<_> = headers.iter().collect();
    if headers.is_empty() || unique.len() != headers.len() || headers.iter().any(str::is_empty) {
        return Err(invalid("empty or duplicate headers"));
    }
    required(&headers, "variant_id")?;
    let indices: Vec<usize> = CRITERIA_ORDER
        .iter()
        .map(|name| required(&headers, name))
        .collect::<Result<_>>()?;
    let evidence_index = evidence_name
        .map(|name| required(&headers, name))
        .transpose()?;
    if evidence_index.is_some_and(|index| indices.contains(&index)) {
        return Err(invalid(
            "ClinGen evidence column must not be a criterion-strength column",
        ));
    }
    let rows = reader
        .records()
        .collect::<std::result::Result<Vec<_>, _>>()?;
    if rows.is_empty() {
        // The reference's empty DataFrame.apply does not produce the summary
        // Series expected by assignment. Do not silently repair that failure.
        return Err(Box::new(super::SourceError(
            "ValueError",
            "empty criteria table is outside this qualified source boundary".into(),
        )));
    }
    let mut columns = Vec::new();
    for (&name, &index) in CRITERIA_ORDER.iter().zip(&indices) {
        let strengths: Vec<i64> = rows.iter().enumerate().map(|(row_number, row)| {
            let strength = row[index].parse::<i64>().map_err(|_| invalid(format!(
                "row {}: {name} must be an integer 0..5; missing/fractional values are outside this boundary", row_number + 2)))?;
            if !(0..=5).contains(&strength) {
                return Err(invalid(format!("row {}: {name} is outside 0..5", row_number + 2)));
            }
            Ok(strength)
        }).collect::<Result<_>>()?;
        columns.push(Column::new(name.into(), strengths));
    }
    let mut matrix = DataFrame::new(columns)?;
    let mut summaries: Vec<Option<String>> = Vec::with_capacity(rows.len());
    let mut overridden = Vec::with_capacity(rows.len());
    for (index, row) in rows.iter().enumerate() {
        let criteria = matrix_row(&matrix, index)?;
        let evidence = evidence_index.map(|column| &row[column]);
        summaries.push(create_criteria_summary(&criteria, evidence)?);
        overridden.push(evidence.is_some_and(|value| !python_strip(value).is_empty()));
    }
    // The reference loops columns outside rows. Preserve this failure/mutation
    // order even though a row-first loop could appear equivalent on good input.
    for &name in &CRITERIA_ORDER {
        let current = matrix.column(name)?.i64()?;
        let mut strengths = Vec::with_capacity(rows.len());
        for index in 0..rows.len() {
            strengths.push(if overridden[index] {
                translate_strength_level(name, summaries[index].as_deref())?
            } else {
                current.get(index).unwrap()
            });
        }
        matrix.with_column(Column::new(name.into(), strengths))?;
    }
    let mut scores = Vec::with_capacity(rows.len());
    for index in 0..rows.len() {
        scores.push(calculate_posterior_probability(&matrix_row(
            &matrix, index,
        )?));
    }

    let integer_score_column = scores
        .iter()
        .all(|(score, _)| matches!(score, Probability::Integer(_)));

    // All scientific validation finishes before any output is opened.
    let output = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&args[1])?;
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(output);
    let mut output_header = headers.clone();
    let derived_indices: Vec<usize> = ["ACMG_criteria", "ACMG_quant_score", "ACMG_class"]
        .iter()
        .map(|name| {
            if let Some(index) = output_header.iter().position(|value| value == *name) {
                index
            } else {
                let index = output_header.len();
                output_header.push_field(name);
                index
            }
        })
        .collect();
    writer.write_record(&output_header)?;
    for (row_index, source_row) in rows.iter().enumerate() {
        let mut row: Vec<String> = source_row.iter().map(str::to_owned).collect();
        row.resize(output_header.len(), String::new());
        for (&name, &column_index) in CRITERIA_ORDER.iter().zip(&indices) {
            // Keep original strings on unchanged rows; rewritten strengths use
            // canonical integer spelling, just as the explicit typed contract.
            if overridden[row_index] {
                row[column_index] = matrix
                    .column(name)?
                    .i64()?
                    .get(row_index)
                    .unwrap()
                    .to_string();
            }
        }
        row[derived_indices[0]] = summaries[row_index].clone().unwrap_or_default();
        let (score, class) = &scores[row_index];
        let probability = score.as_f64();
        row[derived_indices[1]] = if integer_score_column {
            "0".to_owned()
        } else if probability.fract() == 0.0 {
            format!("{probability:.1}")
        } else {
            probability.to_string()
        };
        row[derived_indices[2]] = (*class).to_owned();
        writer.write_record(row)?;
    }
    writer.flush()?;
    Ok(())
}

fn matrix_row(matrix: &DataFrame, index: usize) -> Result<Vec<(&'static str, i64)>> {
    CRITERIA_ORDER
        .iter()
        .map(|&name| Ok((name, matrix.column(name)?.i64()?.get(index).unwrap())))
        .collect()
}
