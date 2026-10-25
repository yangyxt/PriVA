//! Optional procedure transport, not a production per-stage process pipeline.
use super::control_multiply::{RowIndex, control_multiply};
use super::pedigree_context_transport::read_frame;
use crate::acmg_variant_mechanism::variant_mechanism_masks::Cell;
use crate::cli::{Result, invalid};
use std::collections::HashMap;
use std::fs::OpenOptions;
use std::path::Path;

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args == ["--help"] {
        println!(
            "priva control-multiply --input INPUT.csv --output NEW.csv --metadata-output NEW.csv --row-index-kind flat|multi [--column-dtype POSITION=DTYPE] [--axis-metadata OPAQUE]\nOriginal L528 only after unchanged Python C1-C5. In-memory stage mutates native float64 sort_index, preserving columns/opaque axes. Typed CLI emits only modified score bits and target position. Other dtypes, two-NaN payloads, MultiIndex warnings and >1000000 rows defer to original L528 on current buffers; do not replay C1-C5/earlier stages. No production activation."
        );
        return Ok(());
    }
    let (mut input, mut output, mut metadata, mut index_kind) = (None, None, None, None);
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
            "--row-index-kind" => index_kind = Some(pair[1].as_str()),
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
            _ => return Err(invalid("unknown control-multiply option")),
        }
    }
    let input = input.ok_or_else(|| invalid("missing --input"))?;
    let output = output.ok_or_else(|| invalid("missing --output"))?;
    let metadata = metadata.ok_or_else(|| invalid("missing --metadata-output"))?;
    let row_index = match index_kind {
        Some("flat") => RowIndex::Flat,
        Some("multi") => RowIndex::MultiIndex,
        _ => return Err(invalid("--row-index-kind must be flat or multi")),
    };
    if output == input
        || metadata == input
        || output == metadata
        || Path::new(output).exists()
        || Path::new(metadata).exists()
    {
        return Err(invalid("outputs must be new separate files"));
    }
    let mut frame = read_frame(input, dtypes, axis)?;
    let target = match control_multiply(&mut frame, row_index) {
        Ok(target) => target,
        Err(error) => {
            if error.0 == "KeyError" {
                eprintln!("source_failure_line=528");
            }
            return Err(Box::new(error));
        }
    };
    let mut writer = csv::Writer::from_writer(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output)?,
    );
    writer.write_record(["source_position", "index_kind", "index_value", "float_bits"])?;
    for (i, index) in frame.index.iter().enumerate() {
        let Cell::Float(value) = frame.columns[target].values[i] else {
            unreachable!("validated output dtype")
        };
        writer.write_record([
            i.to_string(),
            index[0].clone(),
            index[1].clone(),
            format!("{:016x}", value.to_bits()),
        ])?;
    }
    writer.flush()?;
    let mut writer = csv::Writer::from_writer(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(metadata)?,
    );
    writer.write_record(["name", "value"])?;
    writer.write_record(["target_source_column_position", &target.to_string()])?;
    writer.write_record(["dtype", "float64"])?;
    writer.write_record(["axis", &frame.axis_metadata])?;
    writer.flush()?;
    Ok(())
}
