//! Separate opt-in scalar transport; no pandas dtype inference or production hook.
use super::coerce_pext_value::{Value, coerce_values};
use crate::acmg_variant_mechanism::transport::decode;
use crate::cli::{Result, invalid};
use std::fs::OpenOptions;

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args == ["--help"] {
        println!(
            "priva coerce-pext --input INPUT.csv --output NEW.csv\nHeader: index_kind,index_value,value_kind,value. Strings are UTF-8 hex; f64/f32 are IEEE hex. Existing primitive kinds plus numpy_bool, float32 and decimal bigint. Output values are f64 IEEE bits (NaN is not null); index tokens/order are retained. Unsupported objects/arrays/complex/f16/longdouble/time/category/surrogate strings require unchanged Python. Opt-in only."
        );
        return Ok(());
    }
    let (mut input, mut output) = (None, None);
    for pair in args.chunks(2) {
        if pair.len() != 2 {
            return Err(invalid("named options require values"));
        }
        let target = match pair[0].as_str() {
            "--input" => &mut input,
            "--output" => &mut output,
            _ => return Err(invalid("unknown coerce-pext option")),
        };
        if target.is_some() {
            return Err(invalid("repeated coerce-pext option"));
        }
        *target = Some(&pair[1]);
    }
    let mut reader = csv::Reader::from_path(input.ok_or_else(|| invalid("missing --input"))?)?;
    if reader.headers()?.iter().collect::<Vec<_>>()
        != ["index_kind", "index_value", "value_kind", "value"]
    {
        return Err(invalid("invalid typed coerce-pext header"));
    }
    let mut labels = Vec::new();
    let mut values = Vec::new();
    for row in reader.records() {
        let row = row?;
        decode(&row[0], &row[1])?; // Validate opaque primitive labels; do not join/sort/cast them.
        labels.push([row[0].to_owned(), row[1].to_owned()]);
        let value = match &row[2] {
            "numpy_bool" => match &row[3] {
                "0" => Value::NumpyBool(false),
                "1" => Value::NumpyBool(true),
                _ => return Err(invalid("invalid numpy_bool")),
            },
            "float32" => Value::Float32(f32::from_bits(u32::from_str_radix(&row[3], 16)?)),
            "bigint" => Value::BigInteger(row[3].to_owned()),
            kind => Value::Primitive(decode(kind, &row[3])?),
        };
        values.push(value);
    }
    let result = coerce_values(&values, "value")?;
    let file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output.ok_or_else(|| invalid("missing --output"))?)?;
    let mut writer = csv::Writer::from_writer(file);
    writer.write_record(["index_kind", "index_value", "value_kind", "value"])?;
    for (label, value) in labels.iter().zip(result.f64()?.into_no_null_iter()) {
        writer.write_record([
            label[0].as_str(),
            label[1].as_str(),
            "float",
            &format!("{:016x}", value.to_bits()),
        ])?;
    }
    writer.flush()?;
    Ok(())
}
