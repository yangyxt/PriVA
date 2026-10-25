//! Opt-in R1 primitive transport; output is a plan over original caller buffers.
use super::ranking_setup::ranking_setup;
use super::terminal_order::{Input, InputColumn};
use crate::acmg_variant_mechanism::transport::decode;
use crate::cli::{Result, invalid};
use std::collections::HashMap;
use std::fs::OpenOptions;
use std::path::Path;

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args == ["--help"] {
        println!(
            "priva rank-setup --input INPUT.csv --output NEW.csv --metadata-output NEW.csv [--column-dtype POSITION=DTYPE] [--axis-metadata OPAQUE]\nR1 copy/drop/score-column-copy/protein-coding filter only. No score arithmetic, PED, SciPy or production dispatch. Typed primitive cells, opaque index tokens, unique required/assigned labels; unrelated duplicate labels retained. Outputs row positions and column-copy/drop plan over original buffers. Unsupported objects/dtypes require Python."
        );
        return Ok(());
    }
    let (mut input, mut output, mut metadata) = (None, None, None);
    let mut axis = String::new();
    let mut dtypes = HashMap::new();
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
                if dtypes
                    .insert(position.parse::<usize>()?, dtype.to_owned())
                    .is_some()
                {
                    return Err(invalid("repeated dtype position"));
                }
            }
            _ => return Err(invalid("unknown rank-setup option")),
        }
    }
    let input = input.ok_or_else(|| invalid("missing --input"))?;
    let output = output.ok_or_else(|| invalid("missing --output"))?;
    let metadata = metadata.ok_or_else(|| invalid("missing --metadata-output"))?;
    if input == output
        || input == metadata
        || output == metadata
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
        for (i, column) in columns.iter_mut().enumerate() {
            let kind = &row[2 + 2 * i];
            column.values.push(decode(
                if kind == "pd_na" { "null" } else { kind },
                &row[3 + 2 * i],
            )?);
        }
    }
    let plan = ranking_setup(&Input {
        columns: columns.clone(),
        index,
        axis_metadata: axis,
    })?;
    let mut writer = csv::Writer::from_writer(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output)?,
    );
    writer.write_record(["source_position", "index_kind", "index_value"])?;
    for (row, &position) in plan.source_positions.iter().enumerate() {
        writer.write_record([
            position.to_string(),
            plan.index[row][0].clone(),
            plan.index[row][1].clone(),
        ])?;
    }
    writer.flush()?;
    let mut writer = csv::Writer::from_writer(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(metadata)?,
    );
    writer.write_record([
        "record",
        "position",
        "name",
        "dtype",
        "source_column",
        "operation",
    ])?;
    for (position, column) in plan.columns.iter().enumerate() {
        writer.write_record([
            "column",
            &position.to_string(),
            &column.name,
            &column.dtype,
            &column.source_column.to_string(),
            if column.copy_score {
                "copy_score"
            } else {
                "retain"
            },
        ])?;
    }
    for &position in &plan.dropped_columns {
        writer.write_record([
            "drop",
            &position.to_string(),
            &columns[position].name,
            &columns[position].dtype,
            &position.to_string(),
            "remove_legacy",
        ])?;
    }
    writer.write_record(["axis", "", &plan.axis_metadata, "opaque", "", ""])?;
    writer.flush()?;
    Ok(())
}
