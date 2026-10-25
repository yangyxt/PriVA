//! Opt-in fixture transport. Ordered workflows use the in-memory R4/R5/R6 APIs.
use super::control_multiply::RowIndex;
use super::pedigree_context_transport::read_frame;
use super::ranking_gene_weights::{GeneList, Penalty, ranking_gene_weights};
use crate::cli::{Result, invalid};
use std::collections::HashMap;
use std::fs::OpenOptions;
use std::path::Path;

fn selected_file(kind: &str, path: Option<&String>) -> Result<Option<String>> {
    match (kind, path) {
        ("none", None) => Ok(None),
        ("file", Some(path)) => {
            let bytes = std::fs::read(path).map_err(|_| {
                invalid("UnsupportedInput: selected gene-file I/O; retain Python backend")
            })?;
            String::from_utf8(bytes).map(Some).map_err(|_| {
                invalid("UnsupportedInput: gene-file decoding/compression; retain Python backend")
            })
        }
        ("unsupported", _) => Err(invalid(
            "UnsupportedInput: selected gene-file argument; retain Python backend",
        )),
        _ => Err(invalid(
            "gene selection requires --*-kind none|file|unsupported and --*-file only with file",
        )),
    }
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args == ["--help"] {
        println!(
            "priva ranking-gene-weights --input INPUT.csv --penalty-input MASK.csv --penalty-dtype bool|boolean|object --row-index-kind flat|multi --dispensable-kind none|file|unsupported --relevant-kind none|file|unsupported [--dispensable-file FILE] [--relevant-file FILE] --output NEW.csv --metadata-output NEW.csv [--column-dtype POSITION=DTYPE] [--axis-metadata OPAQUE]\nOriginal R5 then R6 L599-616, preserving three separate multiplications and the existing R4 mask. Only literal plain-text parser subset, native float64 and aligned flat Boolean masks are qualified. All other parser/error/dtype/object routes defer atomically to original L599-616 on current buffers. No earlier-stage replay, production dispatch or family activation."
        );
        return Ok(());
    }
    let (mut input, mut mask_input, mut dtype, mut index_kind, mut output, mut metadata) =
        (None, None, None, None, None, None);
    let (mut disp_kind, mut relevant_kind, mut disp_file, mut relevant_file) =
        (None, None, None, None);
    let (mut dtypes, mut axis) = (HashMap::new(), String::new());
    for pair in args.chunks(2) {
        if pair.len() != 2 {
            return Err(invalid("named options require values"));
        }
        match pair[0].as_str() {
            "--input" => input = Some(&pair[1]),
            "--penalty-input" => mask_input = Some(&pair[1]),
            "--penalty-dtype" => dtype = Some(pair[1].as_str()),
            "--row-index-kind" => index_kind = Some(pair[1].as_str()),
            "--dispensable-kind" => disp_kind = Some(pair[1].as_str()),
            "--relevant-kind" => relevant_kind = Some(pair[1].as_str()),
            "--dispensable-file" => disp_file = Some(&pair[1]),
            "--relevant-file" => relevant_file = Some(&pair[1]),
            "--output" => output = Some(&pair[1]),
            "--metadata-output" => metadata = Some(&pair[1]),
            "--axis-metadata" => axis = pair[1].clone(),
            "--column-dtype" => {
                let (position, value) = pair[1]
                    .split_once('=')
                    .ok_or_else(|| invalid("expected POSITION=DTYPE"))?;
                if dtypes
                    .insert(position.parse::<usize>()?, value.to_owned())
                    .is_some()
                {
                    return Err(invalid("repeated dtype position"));
                }
            }
            _ => return Err(invalid("unknown ranking-gene-weights option")),
        }
    }
    let input = input.ok_or_else(|| invalid("missing --input"))?;
    let mask_input = mask_input.ok_or_else(|| invalid("missing --penalty-input"))?;
    let output = output.ok_or_else(|| invalid("missing --output"))?;
    let metadata = metadata.ok_or_else(|| invalid("missing --metadata-output"))?;
    if output == metadata
        || [Some(input), Some(mask_input), disp_file, relevant_file]
            .iter()
            .flatten()
            .any(|path| *path == output || *path == metadata)
        || Path::new(output).exists()
        || Path::new(metadata).exists()
    {
        return Err(invalid("outputs must be new separate files"));
    }
    let row_index = match index_kind {
        Some("flat") => RowIndex::Flat,
        Some("multi") => RowIndex::MultiIndex,
        _ => return Err(invalid("--row-index-kind must be flat or multi")),
    };
    let dtype = match dtype {
        Some("bool") => "bool",
        Some("boolean") => "boolean",
        Some("object") => "object",
        Some(_) => {
            return Err(invalid(
                "UnsupportedInput: R4 mask dtype; retain Python backend",
            ));
        }
        None => return Err(invalid("missing --penalty-dtype")),
    };
    let mut frame = read_frame(input, dtypes, axis)?;
    let mut reader = csv::Reader::from_path(mask_input)?;
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
        return Err(invalid("invalid penalty transport header"));
    }
    let (mut index, mut data, mut null) = (Vec::new(), Vec::new(), Vec::new());
    for row in reader.records() {
        let row = row?;
        if dtype == "object" {
            return Err(invalid(
                "UnsupportedInput: nonempty object R4 mask; retain Python backend",
            ));
        }
        let bit = |value: &str| -> Result<bool> {
            match value {
                "0" => Ok(false),
                "1" => Ok(true),
                _ => Err(invalid("invalid Boolean storage bit")),
            }
        };
        let missing = bit(&row[5])?;
        let raw = bit(if dtype == "boolean" { &row[4] } else { &row[3] })?;
        if dtype == "bool" && (missing || &row[2] != "bool")
            || dtype == "boolean"
                && (!missing && (&row[2] != "bool" || bit(&row[3])? != raw)
                    || missing && &row[2] != "pd_na")
        {
            return Err(invalid(
                "UnsupportedInput: R4 mask storage; retain Python backend",
            ));
        }
        index.push([row[0].to_owned(), row[1].to_owned()]);
        data.push(raw);
        null.push(missing);
    }
    // Preserve selected-file order. A rejection has not changed current buffers.
    let dispensable = selected_file(
        disp_kind.ok_or_else(|| invalid("missing --dispensable-kind"))?,
        disp_file,
    )?;
    let relevant = selected_file(
        relevant_kind.ok_or_else(|| invalid("missing --relevant-kind"))?,
        relevant_file,
    )?;
    let result = ranking_gene_weights(
        &mut frame,
        &Penalty {
            index,
            dtype,
            data,
            null,
        },
        dispensable
            .as_deref()
            .map_or(GeneList::Disabled, GeneList::PlainText),
        relevant
            .as_deref()
            .map_or(GeneList::Disabled, GeneList::PlainText),
        row_index,
    )?;
    let mut writer = csv::Writer::from_writer(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output)?,
    );
    writer.write_record([
        "step",
        "source_position",
        "index_kind",
        "index_value",
        "selected",
        "multiplier_bits",
        "score_bits",
    ])?;
    for (step_position, step) in result.steps.iter().enumerate() {
        for (i, token) in result.index.iter().enumerate() {
            writer.write_record([
                step_position.to_string(),
                i.to_string(),
                token[0].clone(),
                token[1].clone(),
                u8::from(step.frame.column("selected")?.bool()?.get(i).unwrap()).to_string(),
                format!(
                    "{:016x}",
                    step.frame
                        .column("multiplier")?
                        .f64()?
                        .get(i)
                        .unwrap()
                        .to_bits()
                ),
                format!(
                    "{:016x}",
                    step.frame
                        .column("sort_index")?
                        .f64()?
                        .get(i)
                        .unwrap()
                        .to_bits()
                ),
            ])?;
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
    for (name, value) in [
        ("axis", result.axis_metadata),
        (
            "score_column_position",
            result.score_column_position.to_string(),
        ),
        ("step_count", result.steps.len().to_string()),
    ] {
        writer.write_record([name, &value])?;
    }
    for (i, step) in result.steps.iter().enumerate() {
        for (key, value) in [
            ("name", step.name.to_owned()),
            ("line", step.source_line.to_string()),
            (
                "multiplier_column_position",
                step.multiplier_column_position.to_string(),
            ),
        ] {
            writer.write_record([format!("step.{i}.{key}"), value])?;
        }
    }
    for (role, genes) in [
        ("dispensable", result.dispensable_genes),
        ("relevant", result.relevant_genes),
    ] {
        writer.write_record([format!("{role}.enabled"), genes.is_some().to_string()])?;
        if let Some(genes) = genes {
            writer.write_record([format!("{role}.count"), genes.len().to_string()])?;
            for (i, gene) in genes.into_iter().enumerate() {
                writer.write_record([format!("{role}.gene.{i}"), gene])?;
            }
        }
    }
    writer.flush()?;
    Ok(())
}
