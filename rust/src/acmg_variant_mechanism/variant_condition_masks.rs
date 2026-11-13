//! `_variant_condition_masks`, unchanged Python source lines 267–340 at 09223ef.
//! Source SHA256 c7e8da80ca0ec47fb38a5fd1d965641c1791189cb0d25a651253f44508d0bf26.
//! Typed primitive cells only; this opt-in module does not select a production backend.
use super::variant_mechanism_masks::{Cell, Dtype, Input, InputColumn, Masks};
use crate::acmg_scoring::SourceError;
use polars_core::prelude::*;

pub const NAMES: [&str; 13] = [
    "has_recessive",
    "has_dominant",
    "has_x_linked_recessive",
    "has_x_linked_dominant",
    "has_x_linked_unspecified",
    "has_y_linked",
    "has_mitochondrial",
    "has_dominant_like",
    "has_mendelian",
    "has_non_monogenic",
    "has_non_mendelian",
    "has_incomplete_penetrance",
    "has_complete_penetrance",
];
const REQUIRED: [&str; 2] = ["variant_inheritance", "variant_penetrance"];

fn missing_columns(names: impl Fn(&str) -> bool) -> Result<(), SourceError> {
    // Python sorts the missing set before evaluating either fillna/astype phase.
    let missing: Vec<_> = REQUIRED.into_iter().filter(|n| !names(n)).collect();
    if !missing.is_empty() {
        return Err(SourceError(
            "KeyError",
            format!(
                "'variant-level condition annotations are required; missing columns: {}'",
                missing.join(", ")
            ),
        ));
    }
    Ok(())
}

fn lowercase(column: &InputColumn) -> Result<Vec<String>, SourceError> {
    if let Dtype::Nullable(dtype) = column.dtype {
        if !["Float64", "Int64", "boolean"].contains(&dtype) {
            return Err(SourceError(
                "UnfinishedParity",
                format!("unsupported dtype {dtype}; retain Python backend"),
            ));
        }
        if column
            .values
            .iter()
            .any(|v| matches!(v, Cell::Null) || matches!(v, Cell::Float(n) if n.is_nan()))
        {
            return Err(SourceError(
                "TypeError",
                format!("Invalid value '' for dtype {dtype}"),
            ));
        }
    }
    Ok(column
        .values
        .iter()
        .map(|v| match v {
            Cell::Text(s) => s.to_lowercase(),
            // fillna("").astype(str).str.lower(): every supported nontext scalar
            // spells something that cannot equal any inheritance/penetrance token.
            // Keeping it empty is predicate-equivalent; no string values are exported.
            _ => String::new(),
        })
        .collect())
}

fn convert(input: &Input<impl Clone>, name: &str) -> Result<Vec<String>, SourceError> {
    let mut columns = input.columns.iter().filter(|c| c.name == name);
    let first = columns.next().expect("required columns checked");
    let value = lowercase(first)?;
    if let Some(second) = columns.next() {
        // Duplicate labels select a DataFrame. All its fillna/astype phases run
        // before `.str` raises, including a later nullable-column fillna failure.
        lowercase(second)?;
        for column in columns {
            lowercase(column)?;
        }
        return Err(SourceError(
            "AttributeError",
            "'DataFrame' object has no attribute 'str'".into(),
        ));
    }
    Ok(value)
}

fn has_token(value: &str, token: &str) -> bool {
    // Exactly (?:^|;)TOKEN(?:;|$), with Python's non-MULTILINE `$` also
    // accepting the position before ONE final LF. No whitespace stripping.
    let mut parts = value.split(';').peekable();
    while let Some(part) = parts.next() {
        if part == token || (parts.peek().is_none() && part.strip_suffix('\n') == Some(token)) {
            return true;
        }
    }
    false
}

pub fn variant_condition_masks<I: Clone>(input: &Input<I>) -> Result<Masks<I>, SourceError> {
    missing_columns(|n| input.columns.iter().any(|c| c.name == n))?;
    let len = input.index.len();
    if input.columns.iter().any(|c| c.values.len() != len) {
        return Err(SourceError(
            "ValueError",
            "column/index length mismatch".into(),
        ));
    }
    // Preserve the source's column-phase order even when inheritance has no tags.
    let inheritance = convert(input, REQUIRED[0])?;
    let penetrance = convert(input, REQUIRED[1])?;
    let tokens = [
        "recessive",
        "dominant",
        "x_linked_recessive",
        "x_linked_dominant",
        "x_linked_unspecified",
        "y_linked",
        "mitochondrial",
        "polygenic",
        "digenic",
        "oligogenic",
        "non_mendelian",
    ];
    let token_masks: Vec<Vec<bool>> = tokens
        .iter()
        .map(|token| inheritance.iter().map(|s| has_token(s, token)).collect())
        .collect();
    let mut values = token_masks[..7].to_vec();
    for indices in [
        &[1, 3, 5, 6][..],
        &[0, 1, 2, 3, 4, 5, 6][..],
        &[7, 8, 9][..],
    ] {
        values.push(
            (0..len)
                .map(|row| indices.iter().any(|i| token_masks[*i][row]))
                .collect(),
        );
    }
    values.push(token_masks[10].clone());
    for token in ["incomplete", "complete"] {
        values.push(penetrance.iter().map(|s| has_token(s, token)).collect());
    }
    let frame = DataFrame::new(
        values
            .into_iter()
            .zip(NAMES)
            .map(|(v, n)| Column::new(n.into(), v))
            .collect(),
    )
    .map_err(|e| SourceError("PolarsError", e.to_string()))?;
    Ok(Masks {
        index: input.index.clone(),
        index_names: input.index_names.clone(),
        source_series_names: (0..13)
            .map(|i| Some(REQUIRED[usize::from(i >= 11)]))
            .collect(),
        frame,
    })
}

/// Polars primitive adapter with explicit pandas dtype metadata. Index labels and
/// names remain separate, retaining duplicates and caller-supplied tuple labels.
/// Unsupported dtypes (including empty/all-null ones) require the Python backend.
pub fn from_polars<I: Clone>(
    df: &DataFrame,
    index: Vec<I>,
    index_names: Vec<Option<String>>,
    column_dtypes: [Dtype; 2],
) -> Result<Masks<I>, SourceError> {
    missing_columns(|n| df.column(n).is_ok())?;
    let mut columns = Vec::new();
    for (i, name) in REQUIRED.into_iter().enumerate() {
        let c = df.column(name).expect("required columns checked");
        if !matches!(
            c.dtype(),
            DataType::Null
                | DataType::String
                | DataType::Boolean
                | DataType::Float64
                | DataType::Float32
                | DataType::Int64
                | DataType::Int32
                | DataType::Int16
                | DataType::Int8
                | DataType::UInt64
                | DataType::UInt32
                | DataType::UInt16
                | DataType::UInt8
        ) {
            return Err(SourceError(
                "UnfinishedParity",
                format!(
                    "{name}: dtype {:?} requires original Python; no production dispatch",
                    c.dtype()
                ),
            ));
        }
        let values = (0..c.len())
            .map(|row| {
                Ok(
                    match c
                        .get(row)
                        .map_err(|e| SourceError("PolarsError", e.to_string()))?
                    {
                        AnyValue::Null => Cell::Null,
                        AnyValue::String(s) => Cell::Text(s.into()),
                        AnyValue::StringOwned(s) => Cell::Text(s.to_string()),
                        AnyValue::Boolean(b) => Cell::Bool(b),
                        AnyValue::Float64(n) => Cell::Float(n),
                        AnyValue::Float32(n) => Cell::Float(n as f64),
                        AnyValue::Int64(n) => Cell::Integer(n),
                        AnyValue::Int32(n) => Cell::Integer(n as i64),
                        AnyValue::Int16(n) => Cell::Integer(n as i64),
                        AnyValue::Int8(n) => Cell::Integer(n as i64),
                        AnyValue::UInt64(n) => Cell::Unsigned(n),
                        AnyValue::UInt32(n) => Cell::Unsigned(n as u64),
                        AnyValue::UInt16(n) => Cell::Unsigned(n as u64),
                        AnyValue::UInt8(n) => Cell::Unsigned(n as u64),
                        _ => unreachable!("dtype checked"),
                    },
                )
            })
            .collect::<Result<Vec<_>, SourceError>>()?;
        let column = InputColumn {
            name: name.into(),
            dtype: column_dtypes[i].clone(),
            values,
        };
        lowercase(&column)?; // Earlier nullable failures precede later unsupported dtypes.
        columns.push(column);
    }
    variant_condition_masks(&Input {
        columns,
        index,
        index_names,
    })
}
