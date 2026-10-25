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
