//! Original ranker L528 only, after the unchanged Python C1-C5 return.
//! An in-memory mutation for ordered workflow composition; CLI is a test boundary.
use super::SourceError;
use super::terminal_order::Input;
use crate::acmg_variant_mechanism::variant_mechanism_masks::Cell;
use polars_core::prelude::*;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RowIndex {
    /// Opaque flat index labels, including duplicates/gaps; never used as keys.
    Flat,
    /// pandas .loc can emit index-dependent warnings; retain original stage.
    MultiIndex,
}

fn unsupported(message: &str) -> SourceError {
    SourceError(
        "UnsupportedInput",
        format!("{message}; retain Python backend"),
    )
}

fn required<I>(frame: &Input<I>, name: &str) -> Result<usize, SourceError> {
    let mut columns = frame
        .columns
        .iter()
        .enumerate()
        .filter(|(_, column)| column.name == name);
    let (position, _) = columns
        .next()
        .ok_or_else(|| SourceError("KeyError", format!("'{name}'")))?;
    if columns.next().is_some() {
        return Err(unsupported("duplicate required multiplication label"));
    }
    Ok(position)
}

/// Mutate only the existing sort_index column and return its source position.
///
/// Initial contract: native float64 operands, flat opaque axes, at most one
/// million rows (the pinned pandas standard NumPy expression route). Other
/// dtypes/native routes/axis warnings/two-NaN payloads require original L528 on
/// current buffers.
/// Every rejection occurs before mutation; C1-C5 and prior stages are never run.
pub fn control_multiply<I>(
    frame: &mut Input<I>,
    row_index: RowIndex,
) -> Result<usize, SourceError> {
    if row_index == RowIndex::MultiIndex {
        return Err(unsupported("MultiIndex assignment warning route"));
    }
    // Preserve the source RHS's left lookup before its right lookup.
    let target = required(frame, "sort_index")?;
    let multiplier = required(frame, "control_common_index")?;
    let n = frame.index.len();
    if n > 1_000_000 {
        return Err(unsupported("unqualified pandas numexpr-size route"));
    }
    if frame.columns.iter().any(|column| column.values.len() != n) {
        return Err(unsupported("inconsistent column lengths"));
    }
    for position in [target, multiplier] {
        let column = &frame.columns[position];
        if column.dtype != "float64"
            || !column
                .values
                .iter()
                .all(|cell| matches!(cell, Cell::Float(_)))
        {
            return Err(unsupported("native float64 multiplication required"));
        }
    }
    // The frozen NumPy and Rust paths differ in payload/sign/quieting when both
    // operands are NaN. IEEE value equality does not establish the required bit
    // parity. Leave current buffers untouched for the original source stage.
    if frame.columns[target]
        .values
        .iter()
        .zip(&frame.columns[multiplier].values)
        .any(|(left, right)| matches!((left,right),(Cell::Float(a),Cell::Float(b)) if a.is_nan() && b.is_nan()))
    {
        return Err(unsupported("two-NaN operand payload route"));
    }
    let products = frame.columns[target]
        .values
        .iter()
        .zip(&frame.columns[multiplier].values)
        .map(|(left, right)| {
            let (Cell::Float(left), Cell::Float(right)) = (left, right) else {
                unreachable!("validated float64 operands")
            };
            left * right
        })
        .collect::<Vec<_>>();
    let column = Column::new("sort_index".into(), products);
    let values = column
        .f64()
        .map_err(|error| SourceError("PolarsError", error.to_string()))?
        .into_no_null_iter()
        .map(Cell::Float)
        .collect();
    frame.columns[target].values = values;
    Ok(target)
}
