//! Posterior-only transport for already assigned, canonical integer criteria.
//!
//! The scalar implementation and its default model are unchanged. No summary,
//! replacement, matrix publication or ranking procedure executes here.
use super::{
    CRITERIA_ORDER,
    calculate_posterior_probability::{Probability, calculate_posterior_probability},
};
use crate::cli::{Result, invalid};
use std::{fs::OpenOptions, path::Path};

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args == ["--help"] {
        println!(
            "priva posterior-probability --input INPUT.csv --output NEW.csv\nAlready assigned canonical28 integer strengths0..5; default prior0.1, exp_base2, odds_pvst350 only. Output preserves ordinal, scalar kind, exact float64 bits and class. No summaries, criterion assignment, matrix writes, ranking or default activation."
        );
        return Ok(());
    }
    let (mut input, mut output) = (None, None);
    for pair in args.chunks(2) {
        if pair.len() != 2 {
            return Err(invalid(
                "named posterior-probability options require values",
            ));
        }
        match pair[0].as_str() {
            "--input" if input.is_none() => input = Some(&pair[1]),
            "--output" if output.is_none() => output = Some(&pair[1]),
            _ => return Err(invalid("unknown or repeated posterior-probability option")),
        }
    }
    let input = input.ok_or_else(|| invalid("missing --input"))?;
    let output = output.ok_or_else(|| invalid("missing --output"))?;
    if input == output || Path::new(output).exists() {
        return Err(invalid("output must be a new separate file"));
    }
    let mut reader = csv::Reader::from_path(input)?;
    let expected: Vec<_> = ["source_position"]
        .into_iter()
        .chain(CRITERIA_ORDER)
        .collect();
    if reader.headers()?.iter().ne(expected.iter().copied()) {
        return Err(invalid("invalid posterior-probability input header"));
    }
    let mut results = Vec::new();
    for (position, row) in reader.records().enumerate() {
        let row = row?;
        if row[0] != position.to_string() {
            return Err(invalid(
                "source positions must be the complete original row order",
            ));
        }
        let mut criteria = Vec::with_capacity(CRITERIA_ORDER.len());
        for (column, name) in CRITERIA_ORDER.iter().enumerate() {
            let value = match &row[column + 1] {
                "0" => 0,
                "1" => 1,
                "2" => 2,
                "3" => 3,
                "4" => 4,
                "5" => 5,
                _ => {
                    return Err(invalid(
                        "criterion strengths must have canonical integer spelling0..5",
                    ));
                }
            };
            criteria.push((*name, value));
        }
        results.push(calculate_posterior_probability(&criteria));
    }
    if results.is_empty() {
        return Err(invalid("empty posterior request requires original Python"));
    }
    // Complete every calculation before opening the new result. create_new
    // refuses existing outputs, aliases and the exists/open publication race.
    let file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)?;
    let mut writer = csv::Writer::from_writer(file);
    writer.write_record([
        "source_position",
        "probability_kind",
        "probability_value",
        "ACMG_class",
    ])?;
    for (position, (probability, class)) in results.into_iter().enumerate() {
        let (kind, value) = match probability {
            Probability::Float(value) => ("float64", format!("{:016x}", value.to_bits())),
            Probability::Integer(value) => ("int", value.to_string()),
        };
        writer.write_record([
            position.to_string(),
            kind.to_owned(),
            value,
            class.to_owned(),
        ])?;
    }
    writer.flush()?;
    Ok(())
}
