//! Explicit typed CSV transport, not a replacement for pandas read_csv inference.
//! String cells are UTF-8 hex; floats are big-endian IEEE bits; null is separate.
use super::variant_mechanism_masks::*;
use crate::cli::{Result, invalid};
use std::fs::OpenOptions;
pub(crate) fn decode(kind: &str, value: &str) -> Result<Cell> {
    Ok(match kind {
        "null" => Cell::Null,
        "bool" => match value {
            "0" => Cell::Bool(false),
            "1" => Cell::Bool(true),
            _ => return Err(invalid("invalid bool")),
        },
        "int" => Cell::Integer(value.parse()?),
        "uint" => Cell::Unsigned(value.parse()?),
        "float" => Cell::Float(f64::from_bits(u64::from_str_radix(value, 16)?)),
        "str" => {
            if !value.is_ascii() || value.len() % 2 != 0 {
                return Err(invalid("invalid UTF8 hex"));
            }
            let bytes = (0..value.len())
                .step_by(2)
                .map(|i| u8::from_str_radix(&value[i..i + 2], 16))
                .collect::<std::result::Result<Vec<_>, _>>()?;
            Cell::Text(String::from_utf8(bytes)?)
        }
        _ => return Err(invalid("unsupported typed cell; retain Python backend")),
    })
}
pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args == ["--help"] {
        println!(
            "priva mechanism-masks --input INPUT.csv --output NEW.csv [--profile-dtype object|string|Float64|Int64|boolean] [--effect-dtype ...]\nOpt-in typed primitive-cell transport; preserves index tokens. No TSV dtype inference or production activation."
        );
        return Ok(());
    }
    let mut input = None;
    let mut output = None;
    let mut profile = Dtype::Object;
    let mut effect = Dtype::Object;
    let mut score_dtypes = [Dtype::Object, Dtype::Object, Dtype::Object];
    for pair in args.chunks(2) {
        if pair.len() != 2 {
            return Err(invalid("named options require values"));
        }
        let dtype = || match pair[1].as_str() {
            "object" => Ok(Dtype::Object),
            "string" => Ok(Dtype::String),
            "Float64" => Ok(Dtype::Nullable("Float64")),
            "Int64" => Ok(Dtype::Nullable("Int64")),
            "boolean" => Ok(Dtype::Nullable("boolean")),
            _ => Err(invalid("unsupported dtype metadata; retain Python backend")),
        };
        match pair[0].as_str() {
            "--input" => input = Some(&pair[1]),
            "--output" => output = Some(&pair[1]),
            "--profile-dtype" => profile = dtype()?,
            "--effect-dtype" => effect = dtype()?,
            "--lof-dtype" => score_dtypes[0] = dtype()?,
            "--gof-dtype" => score_dtypes[1] = dtype()?,
            "--dn-dtype" => score_dtypes[2] = dtype()?,
            _ => return Err(invalid("unknown mechanism-masks option")),
        }
    }
    let mut reader = csv::Reader::from_path(input.ok_or_else(|| invalid("missing --input"))?)?;
    let headers = reader.headers()?.clone();
    if headers.len() < 2
        || headers.get(0) != Some("index_kind")
        || headers.get(1) != Some("index_value")
        || headers.len() % 2 != 0
    {
        return Err(invalid("invalid typed transport header"));
    }
    let mut columns = Vec::new();
    for i in (2..headers.len()).step_by(2) {
        let name = headers[i + 1].to_owned();
        if headers[i] != format!("{name}_kind") {
            return Err(invalid("mismatched kind header"));
        }
        let dtype = match name.as_str() {
            "var_plausible_patho_mechs" => profile.clone(),
            "variant_effect" => effect.clone(),
            "variant_lof_score" => score_dtypes[0].clone(),
            "variant_gof_score" => score_dtypes[1].clone(),
            "variant_dn_score" => score_dtypes[2].clone(),
            _ => Dtype::Object,
        };
        columns.push(InputColumn {
            name,
            dtype,
            values: Vec::new(),
        });
    }
    let mut index = Vec::new();
    let mut labels = Vec::new();
    for row in reader.records() {
        let row = row?;
        index.push(decode(&row[0], &row[1])?);
        labels.push([row[0].to_string(), row[1].to_string()]);
        for (i, col) in columns.iter_mut().enumerate() {
            col.values.push(decode(&row[2 + 2 * i], &row[3 + 2 * i])?);
        }
    }
    let result = variant_mechanism_masks(&Input {
        columns,
        index,
        index_names: vec![None],
    })?;
    let file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output.ok_or_else(|| invalid("missing --output"))?)?;
    let mut writer = csv::Writer::from_writer(file);
    let mut header = vec!["index_kind", "index_value"];
    header.extend(NAMES);
    writer.write_record(header)?;
    for (i, label) in labels.iter().enumerate() {
        let mut row = label.to_vec();
        for c in result.frame.get_columns() {
            row.push(if c.bool()?.get(i).unwrap() { "1" } else { "0" }.into());
        }
        writer.write_record(row)?;
    }
    writer.flush()?;
    Ok(())
}

/// Separate opt-in condition boundary: the established mechanism route above is
/// unchanged. Column-position dtype overrides retain heterogeneous duplicate labels.
pub(crate) fn dispatch_condition(args: &[String]) -> Result<()> {
    use super::variant_condition_masks::{NAMES as CONDITION_NAMES, variant_condition_masks};
    if args == ["--help"] {
        println!(
            "priva condition-masks --input INPUT.csv --output NEW.csv [--column-dtype POSITION=object|string|Float64|Int64|boolean|float64|int64|uint64|bool] [--index-name NAME] [--index-dtype DTYPE] [--metadata-output NEW.csv]\nPositions are zero-based data columns after the typed index pair. UTF-8 string cells are hex, floats IEEE bits. Index labels/dtype are opaque; a single string index name is optional. No pandas inference, arbitrary objects, categorical/Arrow/time dtypes or MultiIndex reconstruction. Opt-in only; production Python is unchanged."
        );
        return Ok(());
    }
    let (mut input, mut output, mut metadata) = (None, None, None);
    let mut index_name = None;
    let mut index_dtype = None;
    let mut dtypes = Vec::new();
    for pair in args.chunks(2) {
        if pair.len() != 2 {
            return Err(invalid("named options require values"));
        }
        match pair[0].as_str() {
            "--input" => input = Some(&pair[1]),
            "--output" => output = Some(&pair[1]),
            "--metadata-output" => metadata = Some(&pair[1]),
            "--index-name" => index_name = Some(pair[1].clone()),
            "--index-dtype" => index_dtype = Some(pair[1].clone()),
            "--column-dtype" => {
                let (position, dtype) = pair[1]
                    .split_once('=')
                    .ok_or_else(|| invalid("expected POSITION=DTYPE"))?;
                let position: usize = position.parse()?;
                let dtype = match dtype {
                    "object" => Dtype::Object,
                    "string" => Dtype::String,
                    "Float64" => Dtype::Nullable("Float64"),
                    "Int64" => Dtype::Nullable("Int64"),
                    "boolean" => Dtype::Nullable("boolean"),
                    "float64" | "int64" | "uint64" | "bool" => Dtype::Numpy,
                    _ => return Err(invalid("unsupported dtype metadata; retain Python backend")),
                };
                if dtypes.iter().any(|(i, _)| *i == position) {
                    return Err(invalid("repeated column dtype position"));
                }
                dtypes.push((position, dtype));
            }
            _ => return Err(invalid("unknown condition-masks option")),
        }
    }
    let input = input.ok_or_else(|| invalid("missing --input"))?;
    let output = output.ok_or_else(|| invalid("missing --output"))?;
    if metadata.is_some_and(|m| m == input || m == output || std::path::Path::new(m).exists()) {
        return Err(invalid("metadata output must be a new separate file"));
    }
    let mut reader = csv::Reader::from_path(input)?;
    let headers = reader.headers()?.clone();
    if headers.len() < 2
        || headers.get(0) != Some("index_kind")
        || headers.get(1) != Some("index_value")
        || headers.len() % 2 != 0
    {
        return Err(invalid("invalid typed transport header"));
    }
    let mut columns = Vec::new();
    for (position, i) in (2..headers.len()).step_by(2).enumerate() {
        let name = headers[i + 1].to_owned();
        if headers[i] != format!("{name}_kind") {
            return Err(invalid("mismatched kind header"));
        }
        columns.push(InputColumn {
            name,
            dtype: dtypes
                .iter()
                .find(|(i, _)| *i == position)
                .map(|(_, d)| d.clone())
                .unwrap_or(Dtype::Object),
            values: Vec::new(),
        });
    }
    if dtypes.iter().any(|(i, _)| *i >= columns.len()) {
        return Err(invalid("column dtype position out of range"));
    }
    let (mut index, mut labels) = (Vec::new(), Vec::new());
    for row in reader.records() {
        let row = row?;
        index.push(decode(&row[0], &row[1])?);
        labels.push([row[0].to_owned(), row[1].to_owned()]);
        for (i, c) in columns.iter_mut().enumerate() {
            c.values.push(decode(&row[2 + 2 * i], &row[3 + 2 * i])?);
        }
    }
    let result = variant_condition_masks(&Input {
        columns,
        index,
        index_names: vec![index_name],
    })?;
    let file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)?;
    let mut writer = csv::Writer::from_writer(file);
    let mut header = vec!["index_kind", "index_value"];
    header.extend(CONDITION_NAMES);
    writer.write_record(header)?;
    for (i, label) in labels.iter().enumerate() {
        let mut row = label.to_vec();
        for c in result.frame.get_columns() {
            row.push(if c.bool()?.get(i).unwrap() { "1" } else { "0" }.into());
        }
        writer.write_record(row)?;
    }
    writer.flush()?;
    if let Some(path) = metadata {
        let file = OpenOptions::new().write(true).create_new(true).open(path)?;
        let mut writer = csv::Writer::from_writer(file);
        writer.write_record(["property", "position", "kind", "value"])?;
        let hex = |s: &str| {
            s.as_bytes()
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect::<String>()
        };
        for (i, name) in result.index_names.iter().enumerate() {
            writer.write_record([
                "index_name".into(),
                i.to_string(),
                if name.is_some() { "str" } else { "null" }.into(),
                name.as_deref().map(hex).unwrap_or_default(),
            ])?;
        }
        if let Some(dtype) = index_dtype {
            writer.write_record(["index_dtype".into(), "0".into(), "str".into(), hex(&dtype)])?;
        }
        for (i, name) in result.source_series_names.iter().enumerate() {
            writer.write_record([
                "series_name".into(),
                i.to_string(),
                if name.is_some() { "str" } else { "null" }.into(),
                name.map(hex).unwrap_or_default(),
            ])?;
            writer.write_record([
                "series_dtype".into(),
                i.to_string(),
                "str".into(),
                hex("bool"),
            ])?;
        }
        writer.flush()?;
    }
    Ok(())
}
