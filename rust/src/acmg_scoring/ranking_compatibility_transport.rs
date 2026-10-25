//! Optional R4 typed procedure boundary; ordered workflows call the Rust API.
use super::pedigree_context::Context;
use super::pedigree_context_transport::read_frame;
use super::ranking_compatibility::{MetadataOrigin, ROLES, ranking_compatibility};
use crate::acmg_variant_mechanism::variant_mechanism_masks::{Cell, NAMES};
use crate::cli::{Result, invalid};
use polars_core::prelude::*;
use std::collections::HashMap;
use std::fs::OpenOptions;
use std::path::Path;

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args == ["--help"] {
        println!(
            "priva ranking-compatibility --input INPUT.csv --heterozygosity-input MASK.csv --heterozygosity-dtype bool|boolean|object --proband-is-male true|false --proband-is-female true|false --output NEW.csv --metadata-output NEW.csv [--column-dtype POSITION=DTYPE] [--axis-metadata OPAQUE]\nOriginal R4 L537-597 only, after unchanged C1-C5/R3. Reuses mechanism helper and existing R2 context; returns ordered masks/raw null bits, metadata origins, source aliases, INFO counts and compatibility assignment. Unsupported representations defer to original R4 on current buffers. R6 multiplication and production dispatch remain unchanged."
        );
        return Ok(());
    }
    let (mut input, mut het_input, mut het_dtype, mut output, mut metadata) =
        (None, None, None, None, None);
    let (mut male, mut female) = (None, None);
    let mut dtypes = HashMap::new();
    let mut axis = String::new();
    for pair in args.chunks(2) {
        if pair.len() != 2 {
            return Err(invalid("named options require values"));
        }
        match pair[0].as_str() {
            "--input" => input = Some(&pair[1]),
            "--heterozygosity-input" => het_input = Some(&pair[1]),
            "--heterozygosity-dtype" => het_dtype = Some(pair[1].as_str()),
            "--proband-is-male" => male = Some(pair[1].parse::<bool>()?),
            "--proband-is-female" => female = Some(pair[1].parse::<bool>()?),
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
            _ => return Err(invalid("unknown ranking-compatibility option")),
        }
    }
    let input = input.ok_or_else(|| invalid("missing --input"))?;
    let het_input = het_input.ok_or_else(|| invalid("missing --heterozygosity-input"))?;
    let output = output.ok_or_else(|| invalid("missing --output"))?;
    let metadata = metadata.ok_or_else(|| invalid("missing --metadata-output"))?;
    if output == metadata
        || [input, het_input]
            .iter()
            .any(|p| *p == output || *p == metadata)
        || Path::new(output).exists()
        || Path::new(metadata).exists()
    {
        return Err(invalid("outputs must be new separate files"));
    }
    let mask_dtype = match het_dtype {
        Some("bool") => "bool",
        Some("boolean") => "boolean",
        Some("object") => "object",
        Some(_) => {
            return Err(invalid(
                "UnsupportedInput: R2 mask dtype; retain Python backend",
            ));
        }
        None => return Err(invalid("missing --heterozygosity-dtype")),
    };
    let mut frame = read_frame(input, dtypes, axis)?;
    let mut reader = csv::Reader::from_path(het_input)?;
    if reader.headers()?
        != [
            "index_kind",
            "index_value",
            "value_kind",
            "value",
            "raw_data",
            "is_null",
        ]
        .as_slice()
    {
        return Err(invalid("invalid heterozygosity transport header"));
    }
    let (mut index, mut data, mut nulls) = (Vec::new(), Vec::new(), Vec::new());
    for row in reader.records() {
        let row = row?;
        index.push([row[0].to_owned(), row[1].to_owned()]);
        if mask_dtype == "object" {
            return Err(invalid(
                "UnsupportedInput: nonempty object R2 mask; retain Python backend",
            ));
        }
        let boolean = |value: &str| -> Result<bool> {
            match value {
                "0" => Ok(false),
                "1" => Ok(true),
                _ => Err(invalid("invalid Boolean storage bit")),
            }
        };
        let missing = boolean(&row[5])?;
        let raw = if mask_dtype == "boolean" {
            boolean(&row[4])?
        } else {
            boolean(&row[3])?
        };
        if mask_dtype == "bool" && (missing || &row[2] != "bool")
            || mask_dtype == "boolean"
                && (!missing && (&row[2] != "bool" || boolean(&row[3])? != raw)
                    || missing && &row[2] != "pd_na")
        {
            return Err(invalid(
                "UnsupportedInput: R2 mask storage; retain Python backend",
            ));
        }
        data.push(raw);
        nulls.push(missing);
    }
    let context = Context {
        pedigree_position: None,
        proband: None,
        is_male: male.ok_or_else(|| invalid("missing --proband-is-male"))?,
        is_female: female.ok_or_else(|| invalid("missing --proband-is-female"))?,
        mask_dtype,
        mask: DataFrame::new(vec![
            Column::new("value".into(), data),
            Column::new("is_null".into(), nulls),
        ])?,
        index,
        axis_metadata: String::new(),
        inherit_frame_metadata: false,
        warning: None,
    };
    let result = match ranking_compatibility(&mut frame, &context) {
        Ok(value) => value,
        Err(error) => {
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
    let mut header = vec![
        "source_position",
        "index_kind",
        "index_value",
        "chrom_hex",
        "compatibility_bits",
    ];
    header.extend(NAMES);
    header.extend(ROLES.into_iter().filter(|name| *name != "chrom"));
    writer.write_record(header)?;
    for (i, index) in result.index.iter().enumerate() {
        let chrom = result.role("chrom")?.str()?.get(i).unwrap();
        let Cell::Float(compatibility) =
            frame.columns[result.compatibility_column_position].values[i]
        else {
            unreachable!("assigned float64")
        };
        let mut row = vec![
            i.to_string(),
            index[0].clone(),
            index[1].clone(),
            chrom
                .as_bytes()
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect(),
            format!("{:016x}", compatibility.to_bits()),
        ];
        for name in NAMES {
            row.push(
                u8::from(result.mechanism.frame.column(name)?.bool()?.get(i).unwrap()).to_string(),
            );
        }
        for name in ROLES.into_iter().filter(|name| *name != "chrom") {
            row.push(u8::from(result.role(name)?.bool()?.get(i).unwrap()).to_string());
        }
        writer.write_record(row)?;
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
            "target_source_column_position",
            result.compatibility_column_position.to_string(),
        ),
        ("penalty_dtype", result.penalty_dtype.to_owned()),
        ("axis", result.axis_metadata.clone()),
        ("log_message", result.log_message.clone()),
    ] {
        writer.write_record([name, &value])?;
    }
    for name in NAMES.into_iter().chain(ROLES) {
        writer.write_record([
            format!("{name}.series_name_kind"),
            if result.source_name(name).is_some() {
                "str"
            } else {
                "none"
            }
            .to_owned(),
        ])?;
        writer.write_record([
            format!("{name}.series_name"),
            result.source_name(name).unwrap_or("").to_owned(),
        ])?;
        writer.write_record([
            format!("{name}.metadata_origin"),
            match result.metadata_origin(name) {
                MetadataOrigin::SourceColumn(position) => format!("column:{position}"),
                MetadataOrigin::Fresh => "fresh".to_owned(),
            },
        ])?;
    }
    writer.flush()?;
    Ok(())
}
