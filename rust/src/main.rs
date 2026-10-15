//! Experimental terminal ordering only: no criterion assignment or penalties.
//! Original TSV strings are preserved; only row order and variant_rank change.
use polars_core::prelude::*;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::OpenOptions;
use std::io;
use std::path::Path;

type Result<T> = std::result::Result<T, Box<dyn Error>>;

fn invalid(message: impl Into<String>) -> Box<dyn Error> {
    Box::new(io::Error::new(io::ErrorKind::InvalidInput, message.into()))
}

fn required(headers: &csv::StringRecord, name: &str) -> Result<usize> {
    headers
        .iter()
        .position(|item| item == name)
        .ok_or_else(|| invalid(format!("missing {name}")))
}

fn strings(rows: &[csv::StringRecord], index: usize, name: &str) -> Result<Vec<String>> {
    rows.iter()
        .enumerate()
        .map(|(i, row)| {
            let value = &row[index];
            // Missing-value spellings recognized by pandas must not become real IDs.
            if value.is_empty()
                || [
                    "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND",
                    "1.#QNAN", "<NA>", "N/A", "NA", "NULL", "NaN", "None", "n/a", "nan", "null",
                ]
                .contains(&value)
            {
                Err(invalid(format!(
                    "row {}: missing {name} is outside the validated contract",
                    i + 2
                )))
            } else {
                Ok(value.to_owned())
            }
        })
        .collect()
}

fn integer_values(values: &[String], name: &str) -> Result<Vec<i64>> {
    values
        .iter()
        .enumerate()
        .map(|(i, value)| {
            value
                .parse::<i64>()
                .map_err(|_| invalid(format!("row {}: {name} must be a signed integer", i + 2)))
        })
        .collect()
}

fn scores(rows: &[csv::StringRecord], index: usize, name: &str) -> Result<Vec<f64>> {
    rows.iter()
        .enumerate()
        .map(|(i, row)| {
            let value: f64 = row[index]
                .parse()
                .map_err(|_| invalid(format!("row {}: {name} must be finite f64", i + 2)))?;
            if !value.is_finite() {
                return Err(invalid(format!(
                    "row {}: nonfinite {name} is outside the validated contract",
                    i + 2
                )));
            }
            Ok(value)
        })
        .collect()
}

fn final_order(input: &Path, output: &Path, chrom_type: &str) -> Result<()> {
    if !["string", "integer"].contains(&chrom_type) {
        return Err(invalid(
            "--chrom-type must be string or integer, matching the legacy in-memory dtype",
        ));
    }
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(input)?;
    let headers = reader.headers()?.clone();
    let unique: HashSet<_> = headers.iter().collect();
    if headers.is_empty() || unique.len() != headers.len() || headers.iter().any(str::is_empty) {
        return Err(invalid("empty or duplicate headers"));
    }
    let ci = required(&headers, "chrom")?;
    let pi = required(&headers, "pos")?;
    let ri = required(&headers, "ref")?;
    let ai = required(&headers, "alt")?;
    let vi = required(&headers, "variant_id")?;
    let si = required(&headers, "sort_index")?;
    let qi = required(&headers, "ACMG_quant_score")?;
    let rank_index = headers.iter().position(|item| item == "variant_rank");
    let rows = reader
        .records()
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let chrom = strings(&rows, ci, "chrom")?;
    let pos = integer_values(&strings(&rows, pi, "pos")?, "pos")?;
    if pos.iter().any(|value| *value <= 0) {
        return Err(invalid("VCF position must be positive"));
    }
    let reference = strings(&rows, ri, "ref")?;
    let alternate = strings(&rows, ai, "alt")?;
    let variant = strings(&rows, vi, "variant_id")?;
    let sort_index = scores(&rows, si, "sort_index")?;
    let quant = scores(&rows, qi, "ACMG_quant_score")?;

    // Explicit dtypes avoid Polars inference changing chromosome/position order.
    let chrom_column = if chrom_type == "integer" {
        Column::new("chrom".into(), integer_values(&chrom, "chrom")?)
    } else {
        Column::new("chrom".into(), chrom.clone())
    };
    let mut frame = DataFrame::new(vec![
        chrom_column,
        Column::new("pos".into(), pos.clone()),
        Column::new("ref".into(), reference.clone()),
        Column::new("alt".into(), alternate.clone()),
        Column::new("sort_index".into(), sort_index),
        Column::new("ACMG_quant_score".into(), quant),
        Column::new(
            "source_row".into(),
            (0..rows.len() as u64).collect::<Vec<_>>(),
        ),
    ])?;

    // Polars performs the dataframe grouping and reduction. Broadcast the
    // result via an exact coordinate key; this never deduplicates source rows.
    #[allow(deprecated)]
    let maxima = frame
        .group_by(["chrom", "pos", "ref", "alt"])?
        .select(["sort_index"])
        .max()?;
    let max_chrom = maxima.column("chrom")?;
    let max_pos = maxima.column("pos")?.i64()?;
    let max_ref = maxima.column("ref")?.str()?;
    let max_alt = maxima.column("alt")?.str()?;
    let max_score = maxima.column("sort_index_max")?.f64()?;
    let mut by_coordinate = HashMap::new();
    for i in 0..maxima.height() {
        let chromosome = if chrom_type == "integer" {
            max_chrom.i64()?.get(i).unwrap().to_string()
        } else {
            max_chrom.str()?.get(i).unwrap().to_owned()
        };
        by_coordinate.insert(
            (
                chromosome,
                max_pos.get(i).unwrap(),
                max_ref.get(i).unwrap().to_owned(),
                max_alt.get(i).unwrap().to_owned(),
            ),
            max_score.get(i).unwrap(),
        );
    }
    let canonical_chrom = if chrom_type == "integer" {
        integer_values(&chrom, "chrom")?
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
    } else {
        chrom
    };
    let group_scores: Vec<f64> = (0..rows.len())
        .map(|i| {
            by_coordinate[&(
                canonical_chrom[i].clone(),
                pos[i],
                reference[i].clone(),
                alternate[i].clone(),
            )]
        })
        .collect();
    frame.with_column(Column::new("max_variant_score".into(), group_scores))?;

    // pandas' multikey lexsort preserves fully tied source rows. The explicit
    // final key makes that choice independent of Polars thread scheduling.
    let ordered = frame.sort(
        [
            "max_variant_score",
            "ACMG_quant_score",
            "chrom",
            "pos",
            "ref",
            "alt",
            "source_row",
        ],
        SortMultipleOptions::default()
            .with_order_descending_multi([true, true, false, false, false, false, false])
            .with_maintain_order(true),
    )?;

    // Never truncate an existing result or accidentally replace the input.
    let file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)?;
    let mut writer = csv::WriterBuilder::new().delimiter(b'\t').from_writer(file);
    let mut output_header = headers.clone();
    if rank_index.is_none() {
        output_header.push_field("variant_rank");
    }
    writer.write_record(&output_header)?;
    let mut rank_by_variant = HashMap::new();
    for row_number in ordered.column("source_row")?.u64()?.into_no_null_iter() {
        let index = row_number as usize;
        let next_rank = rank_by_variant.len() + 1;
        let rank = *rank_by_variant.entry(&variant[index]).or_insert(next_rank);
        let mut row: Vec<String> = rows[index].iter().map(str::to_owned).collect();
        if let Some(column) = rank_index {
            row[column] = rank.to_string();
        } else {
            row.push(rank.to_string());
        }
        writer.write_record(row)?;
    }
    writer.flush()?;
    Ok(())
}

fn run() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() == 2 && args[1] == "--version" {
        println!(
            "priva-tabular {} (polars-core 0.51.0; final-order experimental)",
            env!("CARGO_PKG_VERSION")
        );
        return Ok(());
    }
    if args.len() != 6 || args[1] != "final-order" || args[4] != "--chrom-type" {
        return Err(invalid(
            "usage: priva-tabular final-order INPUT.tsv NEW_OUTPUT.tsv --chrom-type string|integer\nComputes ordering and rank only; sort_index and ACMG judgement must already exist.",
        ));
    }
    final_order(Path::new(&args[2]), Path::new(&args[3]), &args[5])
}

fn main() {
    if let Err(error) = run() {
        eprintln!("priva-tabular: {error}");
        std::process::exit(2);
    }
}
