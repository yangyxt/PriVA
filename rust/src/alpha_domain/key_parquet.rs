//! Explicit Parquet text/tuple tags and ordered, nullable typed tuple elements.
//! No repr/JSON stringification is used for scientific dictionary keys.
use super::key::Key;
use super::{Result, exception, parquet};
use crate::vcf_variant::Scalar;
use polars_core::prelude::*;

fn atom_type() -> DataType {
    DataType::Struct(vec![
        Field::new("kind".into(), DataType::String),
        Field::new("text".into(), DataType::String),
        Field::new("integer".into(), DataType::Int64),
        Field::new("float".into(), DataType::Float64),
        Field::new("boolean".into(), DataType::Boolean),
    ])
}
fn key_type() -> DataType {
    DataType::Struct(vec![
        Field::new("kind".into(), DataType::String),
        Field::new("text".into(), DataType::String),
        Field::new("elements".into(), DataType::List(Box::new(atom_type()))),
    ])
}
fn atoms(values: &[Scalar]) -> Result<Series> {
    let fields = [
        Series::new(
            "kind".into(),
            values
                .iter()
                .map(|v| match v {
                    Scalar::Missing => "none",
                    Scalar::Text(_) => "text",
                    Scalar::Int(_) => "int",
                    Scalar::Float(_) => "float",
                    Scalar::Bool(_) => "bool",
                })
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "text".into(),
            values
                .iter()
                .map(|v| {
                    if let Scalar::Text(v) = v {
                        Some(v.as_str())
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "integer".into(),
            values
                .iter()
                .map(|v| {
                    if let Scalar::Int(v) = v {
                        Some(*v)
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "float".into(),
            values
                .iter()
                .map(|v| {
                    if let Scalar::Float(v) = v {
                        Some(*v)
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "boolean".into(),
            values
                .iter()
                .map(|v| {
                    if let Scalar::Bool(v) = v {
                        Some(*v)
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>(),
        ),
    ];
    Ok(StructChunked::from_series("".into(), values.len(), fields.iter())?.into_series())
}
fn keys(values: &[Key]) -> Result<Series> {
    let parts = values
        .iter()
        .map(|v| {
            atoms(match v {
                Key::Tuple(items) => items,
                _ => &[],
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let parts = if parts.is_empty() {
        Series::new_empty("elements".into(), &DataType::List(Box::new(atom_type())))
    } else {
        Series::new("elements".into(), parts)
    };
    let fields = [
        Series::new(
            "kind".into(),
            values
                .iter()
                .map(|v| match v {
                    Key::Text(_) => "text",
                    Key::Tuple(_) => "tuple",
                })
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "text".into(),
            values
                .iter()
                .map(|v| {
                    if let Key::Text(v) = v {
                        Some(v.as_ref())
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>(),
        ),
        parts,
    ];
    Ok(StructChunked::from_series("".into(), values.len(), fields.iter())?.into_series())
}
pub(super) fn column<'a>(rows: impl Iterator<Item = &'a [Key]>) -> Result<Series> {
    let values = rows.map(keys).collect::<Result<Vec<_>>>()?;
    Ok(if values.is_empty() {
        Series::new_empty("typed_keys".into(), &DataType::List(Box::new(key_type())))
    } else {
        Series::new("typed_keys".into(), values)
    })
}
pub(super) fn decode(series: Series) -> Result<Vec<Key>> {
    if series.dtype() != &key_type() {
        return Err(exception(
            "ValueError",
            "invalid typed dictionary key schema",
        ));
    }
    // Fields are resolved once; tuple elements are offsets into one flat child.
    let series = series.rechunk();
    let structure = series.struct_()?;
    let kinds = structure.field_by_name("kind")?;
    let kinds = kinds.str()?;
    let texts = structure.field_by_name("text")?;
    let texts = texts.str()?;
    let (offsets, elements) = parquet::list_parts(&structure.field_by_name("elements")?, || {
        exception("ValueError", "null tuple elements")
    })?;
    let fields = elements.struct_()?;
    let tags = fields.field_by_name("kind")?;
    let tags = tags.str()?;
    let element_texts = fields.field_by_name("text")?;
    let element_texts = element_texts.str()?;
    let ints = fields.field_by_name("integer")?;
    let ints = ints.i64()?;
    let floats = fields.field_by_name("float")?;
    let floats = floats.f64()?;
    let booleans = fields.field_by_name("boolean")?;
    let booleans = booleans.bool()?;
    let mut result = Vec::with_capacity(series.len());
    for at in 0..series.len() {
        let tag = kinds
            .get(at)
            .ok_or_else(|| exception("ValueError", "null dictionary key tag"))?;
        let text = texts.get(at);
        let parts = offsets[at]..offsets[at + 1];
        let key = match (tag, text) {
            ("text", Some(text)) if parts.is_empty() => Key::Text(text.into()),
            ("tuple", None) => {
                let mut values = Vec::with_capacity(parts.len());
                for i in parts {
                    let value = match (
                        tags.get(i),
                        element_texts.get(i),
                        ints.get(i),
                        floats.get(i),
                        booleans.get(i),
                    ) {
                        (Some("none"), None, None, None, None) => Scalar::Missing,
                        (Some("text"), Some(v), None, None, None) => Scalar::Text(v.into()),
                        (Some("int"), None, Some(v), None, None) => Scalar::Int(v),
                        (Some("float"), None, None, Some(v), None) => Scalar::Float(v),
                        (Some("bool"), None, None, None, Some(v)) => Scalar::Bool(v),
                        _ => {
                            return Err(exception(
                                "ValueError",
                                "inconsistent tuple element tag/null fields",
                            ));
                        }
                    };
                    values.push(value);
                }
                Key::Tuple(values.into())
            }
            _ => {
                return Err(exception(
                    "ValueError",
                    "inconsistent dictionary key tag/null fields",
                ));
            }
        };
        result.push(key);
    }
    Ok(result)
}
