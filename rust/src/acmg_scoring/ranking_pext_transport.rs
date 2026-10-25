//! Fixture-only scalar-origin transport. Ordered workflows use the in-memory API.
use super::control_multiply::RowIndex;
use super::pedigree_context_transport::read_frame;
use super::ranking_pext::{CoercionColumn, Modulation, Options, ranking_pext, select_column};
use crate::acmg_consequence::coerce_pext_value::Value;
use crate::acmg_variant_mechanism::transport::decode;
use crate::cli::{Result, invalid};
use std::collections::HashMap;
use std::fs::OpenOptions;
use std::path::Path;

fn unsupported(message: &str) -> Box<dyn std::error::Error> {
    invalid(format!(
        "UnsupportedInput: {message}; retain Python backend"
    ))
}
fn bit(value: &str) -> Result<bool> {
    match value {
        "0" => Ok(false),
        "1" => Ok(true),
        _ => Err(invalid("invalid Boolean storage bit")),
    }
}
pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args == ["--help"] {
        println!(
            "priva ranking-pext --input FRAME.csv --options-input OPTIONS.csv --coercion-input VALUES.csv --coercion-column-position N --mask-input MASK.csv --mask-dtype bool|boolean --row-index-kind flat|multi --output NEW.csv --metadata-output NEW.csv [--column-dtype POSITION=DTYPE] [--axis-metadata OPAQUE]\nOpt-in original R7 only, L621-665. Exact original map scalar origins, native float64 scores, finite valid float parameters, current R4 mask, separate elementary operations and source Series scalar power dispatch. Unsupported dtype/object/warning/numerical/axis routes retain original Python on current buffers. No production dispatch or prior-stage replay."
        );
        return Ok(());
    }
    let mut values = HashMap::new();
    let mut dtypes = HashMap::new();
    for pair in args.chunks(2) {
        if pair.len() != 2 {
            return Err(invalid("named options require values"));
        }
        if pair[0] == "--column-dtype" {
            let (position, dtype) = pair[1]
                .split_once('=')
                .ok_or_else(|| invalid("expected POSITION=DTYPE"))?;
            if dtypes
                .insert(position.parse::<usize>()?, dtype.to_owned())
                .is_some()
            {
                return Err(invalid("repeated dtype position"));
            }
        } else if ![
            "--input",
            "--options-input",
            "--coercion-input",
            "--coercion-column-position",
            "--mask-input",
            "--mask-dtype",
            "--row-index-kind",
            "--output",
            "--metadata-output",
            "--axis-metadata",
        ]
        .contains(&pair[0].as_str())
        {
            return Err(invalid("unknown ranking-pext option"));
        } else if values.insert(pair[0].as_str(), pair[1].as_str()).is_some() {
            return Err(invalid("repeated ranking-pext option"));
        }
    }
    let required = |name| {
        values
            .get(name)
            .copied()
            .ok_or_else(|| invalid(format!("missing {name}")))
    };
    let input = required("--input")?;
    let config = required("--options-input")?;
    let scalar_input = required("--coercion-input")?;
    let mask_input = required("--mask-input")?;
    let output = required("--output")?;
    let metadata = required("--metadata-output")?;
    if output == metadata
        || [input, config, scalar_input, mask_input]
            .iter()
            .any(|p| *p == output || *p == metadata)
        || Path::new(output).exists()
        || Path::new(metadata).exists()
    {
        return Err(invalid("outputs must be new separate files"));
    }
    let row_index = match required("--row-index-kind")? {
        "flat" => RowIndex::Flat,
        "multi" => RowIndex::MultiIndex,
        _ => return Err(invalid("invalid row-index kind")),
    };
    let mut frame = read_frame(
        input,
        dtypes,
        values
            .get("--axis-metadata")
            .copied()
            .unwrap_or("")
            .to_owned(),
    )?;
    let mut reader = csv::Reader::from_path(config)?;
    if reader.headers()?.iter().collect::<Vec<_>>() != ["name", "kind", "value"] {
        return Err(invalid("invalid options header"));
    }
    let mut options = HashMap::new();
    for row in reader.records() {
        let row = row?;
        if &row[1] == "pd_na" && &row[0] == "tissues" {
            return Err(unsupported("pd.NA tissue truth"));
        }
        let value = decode(if &row[1] == "pd_na" { "null" } else { &row[1] }, &row[2])?;
        if !["tissues", "cutoff", "floor", "shape"].contains(&&row[0])
            || options.insert(row[0].to_owned(), value).is_some()
        {
            return Err(invalid("invalid/repeated option scalar"));
        }
    }
    let mut take = |key| {
        options
            .remove(key)
            .ok_or_else(|| invalid(format!("missing option {key}")))
    };
    let options = Options {
        tissues: take("tissues")?,
        cutoff: take("cutoff")?,
        floor: take("floor")?,
        shape: take("shape")?,
    };
    let selection = select_column(&frame, &options.tissues)?;
    let mut mask = Modulation {
        index: Vec::new(),
        dtype: "bool",
        data: Vec::new(),
        null: Vec::new(),
    };
    let mut coercion = None;
    // Source skips every unused parameter/mask/score route without a column.
    if selection.column.is_some() {
        mask.dtype = match required("--mask-dtype")? {
            "bool" => "bool",
            "boolean" => "boolean",
            _ => return Err(unsupported("mask dtype")),
        };
        let mut reader = csv::Reader::from_path(mask_input)?;
        if reader.headers()?.iter().collect::<Vec<_>>()
            != [
                "index_kind",
                "index_value",
                "value_kind",
                "value",
                "raw_data",
                "is_null",
            ]
        {
            return Err(invalid("invalid mask header"));
        }
        for row in reader.records() {
            let row = row?;
            let na = bit(&row[5])?;
            let raw = bit(if mask.dtype == "boolean" {
                &row[4]
            } else {
                &row[3]
            })?;
            if mask.dtype == "bool" && (na || &row[2] != "bool")
                || mask.dtype == "boolean"
                    && (na && &row[2] != "pd_na"
                        || !na && (&row[2] != "bool" || bit(&row[3])? != raw))
            {
                return Err(unsupported("raw mask storage"));
            }
            mask.index.push([row[0].to_owned(), row[1].to_owned()]);
            mask.data.push(raw);
            mask.null.push(na);
        }
        let mut reader = csv::Reader::from_path(scalar_input)?;
        if reader.headers()?.iter().collect::<Vec<_>>()
            != ["source_position", "value_kind", "value"]
        {
            return Err(invalid("invalid coercion header"));
        }
        let mut boxed = Vec::new();
        for row in reader.records() {
            let row = row?;
            if row[0].parse::<usize>()? != boxed.len() {
                return Err(invalid("coercion order mismatch"));
            }
            boxed.push(match &row[1] {
                "numpy_bool" => Value::NumpyBool(bit(&row[2])?),
                "float32" => Value::Float32(f32::from_bits(u32::from_str_radix(&row[2], 16)?)),
                "bigint" => Value::BigInteger(row[2].to_owned()),
                kind => Value::Primitive(decode(
                    if kind == "pd_na" { "null" } else { kind },
                    &row[2],
                )?),
            });
        }
        coercion = Some(CoercionColumn {
            position: required("--coercion-column-position")?.parse()?,
            index: frame.index.clone(),
            values: boxed,
        });
    }
    let result = ranking_pext(&mut frame, &options, coercion.as_ref(), &mask, row_index)?;
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
        "selected",
        "pext_numeric_bits",
        "scaled_pext_bits",
        "power_bits",
        "product_bits",
        "sum_bits",
        "pext_sort_index_bits",
        "score_bits",
    ])?;
    if let Some(data) = &result.values {
        for (i, token) in result.index.iter().enumerate() {
            let mut row = vec![
                i.to_string(),
                token[0].clone(),
                token[1].clone(),
                u8::from(result.selected[i]).to_string(),
            ];
            for name in [
                "pext_numeric",
                "scaled_pext",
                "power",
                "product",
                "sum",
                "pext_sort_index",
                "sort_index",
            ] {
                row.push(format!(
                    "{:016x}",
                    data.column(name)?.f64()?.get(i).unwrap().to_bits()
                ));
            }
            writer.write_record(row)?;
        }
    }
    writer.flush()?;
    let mut writer = csv::Writer::from_writer(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(metadata)?,
    );
    writer.write_record(["name", "value"])?;
    let selection = result.selection;
    for (name, value) in [
        ("axis", result.axis_metadata),
        ("applied", result.values.is_some().to_string()),
        ("column", selection.column.unwrap_or_default()),
        ("used_tissue", selection.used_tissue.to_string()),
        (
            "tissue_list_present",
            selection.tissue_list.is_some().to_string(),
        ),
        (
            "tissue_count",
            selection
                .tissue_list
                .as_ref()
                .map_or(0, Vec::len)
                .to_string(),
        ),
        (
            "first_tissue_column_present",
            selection.first_tissue_column.is_some().to_string(),
        ),
        (
            "first_tissue_column",
            selection.first_tissue_column.unwrap_or_default(),
        ),
        (
            "source_column_position",
            result
                .source_column_position
                .map(|i| i.to_string())
                .unwrap_or_default(),
        ),
        (
            "target_column_position",
            result
                .target_column_position
                .map(|i| i.to_string())
                .unwrap_or_default(),
        ),
        (
            "score_column_position",
            result
                .score_column_position
                .map(|i| i.to_string())
                .unwrap_or_default(),
        ),
    ] {
        writer.write_record([name, &value])?;
    }
    for (i, tissue) in selection.tissue_list.unwrap_or_default().iter().enumerate() {
        writer.write_record([format!("tissue.{i}"), tissue.clone()])?;
    }
    writer.flush()?;
    Ok(())
}
