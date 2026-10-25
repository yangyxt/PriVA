//! Opt-in R2 context and raw Boolean/null-mask transport over caller buffers.
use super::pedigree_context::{Family, pedigree_context};
use super::terminal_order::{Input, InputColumn};
use crate::acmg_variant_mechanism::transport::decode;
use crate::cli::{Result, invalid};
use std::collections::HashMap;
use std::fs::OpenOptions;
use std::path::Path;

pub(crate) fn read_frame(
    path: &str,
    mut dtypes: HashMap<usize, String>,
    axis_metadata: String,
) -> Result<Input<[String; 2]>> {
    let mut reader = csv::Reader::from_path(path)?;
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
            let kind = &row[2 + position * 2];
            column.values.push(decode(
                if kind == "pd_na" { "null" } else { kind },
                &row[3 + position * 2],
            )?);
        }
    }
    Ok(Input {
        columns,
        index,
        axis_metadata,
    })
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args == ["--help"] {
        println!(
            "priva pedigree-context --input INPUT.csv --family-input SCALAR.csv [--pedigree-input PED.csv] --output NEW.csv --metadata-output NEW.csv [--column-dtype POSITION=DTYPE] [--pedigree-column-dtype POSITION=DTYPE] [--axis-metadata OPAQUE]\nOriginal R2 only: first affected pedigree row, sex flags, whole-sample-string count('1') mask. Output raw Boolean data plus null bits, opaque axes, selected source row and structured warning. Numeric nullable sample accessors, nontext selected IDs, ambiguous PED labels and unsupported objects/dtypes defer to Python. No control/scoring/production dispatch."
        );
        return Ok(());
    }
    let (mut input, mut family, mut pedigree, mut output, mut metadata) =
        (None, None, None, None, None);
    let mut frame_dtypes = HashMap::new();
    let mut ped_dtypes = HashMap::new();
    let mut axis = String::new();
    for pair in args.chunks(2) {
        if pair.len() != 2 {
            return Err(invalid("named options require values"));
        }
        match pair[0].as_str() {
            "--input" => input = Some(&pair[1]),
            "--family-input" => family = Some(&pair[1]),
            "--pedigree-input" => pedigree = Some(&pair[1]),
            "--output" => output = Some(&pair[1]),
            "--metadata-output" => metadata = Some(&pair[1]),
            "--axis-metadata" => axis = pair[1].clone(),
            "--column-dtype" | "--pedigree-column-dtype" => {
                let (position, dtype) = pair[1]
                    .split_once('=')
                    .ok_or_else(|| invalid("expected POSITION=DTYPE"))?;
                let dtypes = if pair[0] == "--column-dtype" {
                    &mut frame_dtypes
                } else {
                    &mut ped_dtypes
                };
                if dtypes
                    .insert(position.parse::<usize>()?, dtype.to_owned())
                    .is_some()
                {
                    return Err(invalid("repeated dtype position"));
                }
            }
            _ => return Err(invalid("unknown pedigree-context option")),
        }
    }
    let input = input.ok_or_else(|| invalid("missing --input"))?;
    let family = family.ok_or_else(|| invalid("missing --family-input"))?;
    let output = output.ok_or_else(|| invalid("missing --output"))?;
    let metadata = metadata.ok_or_else(|| invalid("missing --metadata-output"))?;
    if output == metadata
        || [Some(input), Some(family), pedigree]
            .into_iter()
            .flatten()
            .any(|path| path == output || path == metadata)
        || Path::new(output).exists()
        || Path::new(metadata).exists()
    {
        return Err(invalid("outputs must be new separate files"));
    }
    let frame = read_frame(input, frame_dtypes, axis)?;
    let pedigree = pedigree
        .map(|path| read_frame(path, ped_dtypes, String::new()))
        .transpose()?;
    let mut scalar = csv::Reader::from_path(family)?;
    let header = scalar.headers()?;
    if header.len() != 2 || header.get(0) != Some("kind") || header.get(1) != Some("value") {
        return Err(invalid("invalid family scalar header"));
    }
    let rows = scalar
        .records()
        .collect::<std::result::Result<Vec<_>, _>>()?;
    if rows.len() != 1 {
        return Err(invalid("one family scalar required"));
    }
    if rows[0][0].starts_with("unsupported_numpy_family_") && pedigree.is_some() {
        return Err(invalid(
            "UnsupportedInput: NumPy family scalar identity; retain Python backend",
        ));
    }
    let family = if &rows[0][0] == "null"
        || (rows[0][0].starts_with("unsupported_numpy_family_") && pedigree.is_none())
    {
        Family::None
    } else {
        Family::Value(decode(
            if &rows[0][0] == "pd_na" {
                "null"
            } else {
                &rows[0][0]
            },
            &rows[0][1],
        )?)
    };
    let context = match pedigree_context(&frame, pedigree.as_ref(), &family) {
        Ok(value) => value,
        Err(error) => {
            if let Some(warning) = error.warning {
                eprintln!("source_warning={}", warning.name());
            }
            if let Some(row) = error.pedigree_position {
                eprintln!("selected_pedigree_position={row}");
            }
            if let Some(line) = error.source_line {
                eprintln!("source_failure_line={line}");
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
    writer.write_record([
        "source_position",
        "index_kind",
        "index_value",
        "value",
        "is_null",
    ])?;
    let values = context.mask.column("value")?.bool()?;
    let missing = context.mask.column("is_null")?.bool()?;
    for (row, index) in context.index.iter().enumerate() {
        writer.write_record([
            row.to_string(),
            index[0].clone(),
            index[1].clone(),
            u8::from(values.get(row).unwrap()).to_string(),
            u8::from(missing.get(row).unwrap()).to_string(),
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
    for (name, value) in [
        (
            "proband_present",
            u8::from(context.proband.is_some()).to_string(),
        ),
        ("proband", context.proband.clone().unwrap_or_default()),
        (
            "selected_pedigree_position",
            context
                .pedigree_position
                .map(|row| row.to_string())
                .unwrap_or_default(),
        ),
        ("is_male", u8::from(context.is_male).to_string()),
        ("is_female", u8::from(context.is_female).to_string()),
        ("mask_dtype", context.mask_dtype.to_owned()),
        (
            "inherit_frame_metadata",
            u8::from(context.inherit_frame_metadata).to_string(),
        ),
        (
            "source_warning",
            context
                .warning
                .map(|warning| warning.name())
                .unwrap_or("none")
                .to_owned(),
        ),
        ("axis", context.axis_metadata),
    ] {
        writer.write_record([name, &value])?;
    }
    writer.flush()?;
    Ok(())
}
