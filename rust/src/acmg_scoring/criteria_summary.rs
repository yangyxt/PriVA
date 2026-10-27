//! Non-statistical summary/replacement block of summarize_acmg_criteria.
//!
//! Matrix construction/validation, ClinGen loading and caller-owned pandas
//! buffers remain in Python. This procedure computes all summaries before
//! rebuilding strengths in the original criterion-then-row order. It does not
//! assign evidence, calculate posteriors, or invoke any ranking stage.
use super::{
    CRITERIA_ORDER, SourceError, create_criteria_summary::create_criteria_summary, python_strip,
    translate_strength_level::translate_strength_level,
};
use crate::cli::{Result, invalid};
use polars_core::prelude::*;
use std::fs::OpenOptions;
use std::path::Path;

pub struct SummaryPlan {
    pub summaries: Vec<Option<String>>,
    pub overridden: Vec<bool>,
    pub strengths: DataFrame,
}

/// Already validated primitive strengths; rows and criterion order are fixed.
pub fn criteria_summary(columns: &[Vec<i64>], evidence: &[Option<String>]) -> Result<SummaryPlan> {
    if columns.len() != CRITERIA_ORDER.len()
        || columns.iter().any(|column| column.len() != evidence.len())
        || evidence.is_empty()
    {
        return Err(invalid(
            "nonempty aligned 28-column criteria matrix required",
        ));
    }
    let mut frame_columns = Vec::with_capacity(columns.len());
    for (&name, values) in CRITERIA_ORDER.iter().zip(columns) {
        if values.iter().any(|value| !(0..=5).contains(value)) {
            return Err(invalid(format!(
                "{name} must contain integer strengths 0..5"
            )));
        }
        frame_columns.push(Column::new(name.into(), values.clone()));
    }
    let mut strengths = DataFrame::new(frame_columns)?;
    let mut summaries = Vec::with_capacity(evidence.len());
    let mut overridden = Vec::with_capacity(evidence.len());
    for (row, entry) in evidence.iter().enumerate() {
        let criteria = CRITERIA_ORDER
            .iter()
            .enumerate()
            .map(|(column, &name)| (name, columns[column][row]))
            .collect::<Vec<_>>();
        summaries.push(create_criteria_summary(&criteria, entry.as_deref())?);
        overridden.push(
            entry
                .as_deref()
                .is_some_and(|value| !python_strip(value).is_empty()),
        );
    }
    // Original Python completes every summary, then every row in PVS1, then
    // every row in PS1, etc. Malformed evidence must retain that failure order.
    for (column, &name) in CRITERIA_ORDER.iter().enumerate() {
        let mut values = columns[column].clone();
        for (row, replace) in overridden.iter().enumerate() {
            if *replace {
                values[row] = translate_strength_level(name, summaries[row].as_deref())?;
            }
        }
        strengths.with_column(Column::new(name.into(), values))?;
    }
    Ok(SummaryPlan {
        summaries,
        overridden,
        strengths,
    })
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args == ["--help"] {
        println!(
            "priva criteria-summary --input INPUT.csv --output NEW.csv\nAlready assigned integer strengths 0..5 and pre-resolved ClinGen strings only. All summaries precede ordered per-criterion replacements; no posterior, criteria assignment, pickle loading, or production activation. Unsupported Python objects/dtypes/alignment require the original Python function."
        );
        return Ok(());
    }
    let (mut input, mut output) = (None, None);
    for pair in args.chunks(2) {
        if pair.len() != 2 {
            return Err(invalid("named criteria-summary options require values"));
        }
        match pair[0].as_str() {
            "--input" if input.is_none() => input = Some(&pair[1]),
            "--output" if output.is_none() => output = Some(&pair[1]),
            _ => return Err(invalid("unknown or repeated criteria-summary option")),
        }
    }
    let input = input.ok_or_else(|| invalid("missing --input"))?;
    let output = output.ok_or_else(|| invalid("missing --output"))?;
    if input == output || Path::new(output).exists() {
        return Err(invalid("output must be a new separate file"));
    }
    let mut reader = csv::Reader::from_path(input)?;
    let expected: Vec<&str> = ["source_position"]
        .into_iter()
        .chain(CRITERIA_ORDER)
        .chain(["evidence_kind", "evidence"])
        .collect();
    if reader.headers()?.iter().ne(expected.iter().copied()) {
        return Err(invalid("invalid criteria-summary input header"));
    }
    let mut columns = vec![Vec::new(); CRITERIA_ORDER.len()];
    let mut evidence = Vec::new();
    for (position, row) in reader.records().enumerate() {
        let row = row?;
        if row[0] != position.to_string() {
            return Err(invalid(
                "source positions must be the complete original row order",
            ));
        }
        for (column, values) in columns.iter_mut().enumerate() {
            let value = row[column + 1].parse::<i64>()?;
            if !(0..=5).contains(&value) {
                return Err(invalid("criterion strengths must be integers 0..5"));
            }
            values.push(value);
        }
        evidence.push(match (&row[29], &row[30]) {
            ("none", "") => None,
            ("str", value) => Some(value.to_owned()),
            _ => return Err(invalid("invalid pre-resolved evidence encoding")),
        });
    }
    if evidence.is_empty() {
        return Err(Box::new(SourceError(
            "ValueError",
            "empty summary requires Python".into(),
        )));
    }
    let plan = criteria_summary(&columns, &evidence)?;
    // Complete the whole plan before opening output. create_new also closes
    // the exists/open race and refuses literal or hard-link input aliases.
    let file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)?;
    let mut writer = csv::Writer::from_writer(file);
    let header: Vec<&str> = ["source_position", "summary_kind", "summary", "overridden"]
        .into_iter()
        .chain(CRITERIA_ORDER)
        .collect();
    writer.write_record(header)?;
    for row in 0..evidence.len() {
        let mut record = vec![
            row.to_string(),
            if plan.summaries[row].is_some() {
                "str"
            } else {
                "nan"
            }
            .into(),
            plan.summaries[row].clone().unwrap_or_default(),
            u8::from(plan.overridden[row]).to_string(),
        ];
        for name in CRITERIA_ORDER {
            record.push(
                plan.strengths
                    .column(name)?
                    .i64()?
                    .get(row)
                    .unwrap()
                    .to_string(),
            );
        }
        writer.write_record(record)?;
    }
    writer.flush()?;
    Ok(())
}
