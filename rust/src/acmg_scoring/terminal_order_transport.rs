//! Explicit test/adapter transport. Original payload/index buffers stay caller owned.
use super::terminal_order::{Input, InputColumn, terminal_order};
use crate::acmg_variant_mechanism::transport::decode;
use crate::acmg_variant_mechanism::variant_mechanism_masks::Cell;
use crate::cli::{Result, invalid};
use std::fs::OpenOptions;
use std::path::Path;

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args == ["--help"] {
        println!(
            "priva terminal-order --input INPUT.csv --output NEW.csv --metadata-output NEW.csv --column-dtype POSITION=DTYPE [--axis-metadata OPAQUE]\nR8 only; scores must already be computed by unchanged Python. Typed cell pairs follow mechanism/condition transport; index kind/value are opaque. Output is a source-row permutation, original index tokens, temporary maxima and ranks; metadata gives exact output column positions/dtypes. The caller gathers original payload buffers/axes. No production selection or C1–R7 numerics. Unsupported objects/dtypes require Python. Existing final-order is unchanged."
        );
        return Ok(());
    }
    let (mut input, mut output, mut metadata) = (None, None, None);
    let mut axis = String::new();
    let mut dtypes = std::collections::HashMap::new();
    for pair in args.chunks(2) {
        if pair.len() != 2 {
            return Err(invalid("named options require values"));
        }
        match pair[0].as_str() {
            "--input" => input = Some(&pair[1]),
            "--output" => output = Some(&pair[1]),
            "--metadata-output" => metadata = Some(&pair[1]),
            "--axis-metadata" => axis = pair[1].clone(),
            "--column-dtype" => {
                let (position, dtype) = pair[1]
                    .split_once('=')
                    .ok_or_else(|| invalid("expected POSITION=DTYPE"))?;
                let position: usize = position.parse()?;
                if dtypes.insert(position, dtype.to_owned()).is_some() {
                    return Err(invalid("repeated dtype position"));
                }
            }
            _ => return Err(invalid("unknown terminal-order option")),
        }
    }
    let input = input.ok_or_else(|| invalid("missing --input"))?;
    let output = output.ok_or_else(|| invalid("missing --output"))?;
    let metadata = metadata.ok_or_else(|| invalid("missing --metadata-output"))?;
    if output == metadata
        || output == input
        || metadata == input
        || Path::new(output).exists()
        || Path::new(metadata).exists()
    {
        return Err(invalid("outputs must be new separate files"));
    }
    let mut reader = csv::Reader::from_path(input)?;
    let header = reader.headers()?.clone();
    if header.len() < 2
        || header.len() % 2 != 0
        || header.get(0) != Some("index_kind")
        || header.get(1) != Some("index_value")
    {
        return Err(invalid("invalid typed transport header"));
    }
    let mut columns = Vec::new();
    for (position, i) in (2..header.len()).step_by(2).enumerate() {
        let name = header[i + 1].to_owned();
        if header[i] != format!("{name}_kind") {
            return Err(invalid("mismatched kind header"));
        }
        columns.push(InputColumn {
            name,
            dtype: dtypes.remove(&position).unwrap_or("object".into()),
            values: Vec::new(),
        });
    }
    if !dtypes.is_empty() {
        return Err(invalid("dtype position outside columns"));
    }
    let mut index = Vec::new();
    for row in reader.records() {
        let row = row?;
        index.push([row[0].to_owned(), row[1].to_owned()]);
        for (position, column) in columns.iter_mut().enumerate() {
            let kind = &row[2 + 2 * position];
            if kind == "pd_na"
                && column.dtype == "object"
                && ["chrom", "pos", "ref", "alt"].contains(&column.name.as_str())
            {
                return Err(invalid(
                    "UnsupportedInput: object coordinate pd.NA requires Python backend",
                ));
            }
            column.values.push(decode(
                if kind == "pd_na" { "null" } else { kind },
                &row[3 + 2 * position],
            )?);
        }
    }
    let result = terminal_order(&Input {
        columns,
        index,
        axis_metadata: axis,
    })?;
    let mut writer = csv::Writer::from_writer(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output)?,
    );
    writer.write_record([
        "source_position",
        "index_kind",
        "index_value",
        "max_kind",
        "max_value",
        "variant_rank",
    ])?;
    for (row, &position) in result.source_positions.iter().enumerate() {
        let (kind, value) = match &result.maxima[position] {
            Cell::Float(value) => ("float", format!("{:016x}", value.to_bits())),
            Cell::Integer(value) => ("int", value.to_string()),
            Cell::Unsigned(value) if *value <= i64::MAX as u64 => ("int", value.to_string()),
            Cell::Unsigned(value) => ("uint", value.to_string()),
            Cell::Null => ("null", String::new()),
            _ => unreachable!(),
        };
        writer.write_record([
            position.to_string(),
            result.index[row][0].clone(),
            result.index[row][1].clone(),
            kind.into(),
            value,
            result.ranks[row].to_string(),
        ])?;
    }
    writer.flush()?;
    let mut writer = csv::Writer::from_writer(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(metadata)?,
    );
    writer.write_record(["record", "position", "name", "dtype"])?;
    for (i, (name, dtype)) in result.output_columns.iter().enumerate() {
        writer.write_record(["column", &i.to_string(), name, dtype])?;
    }
    writer.write_record([
        "temporary",
        &result.temporary_column_position.to_string(),
        "max_variant_score",
        &result.maximum_dtype,
    ])?;
    writer.write_record(["axis", "", &result.axis_metadata, "opaque"])?;
    writer.flush()?;
    Ok(())
}
