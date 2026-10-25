//! Opt-in R8 boundary, acmg_scoring.py 09223ef lines668–697.
//!
//! Input is the actual post-R7 checkpoint: no scores/penalties are calculated.
//! Return a row permutation, original-axis tokens, temporary maxima, ranks and
//! output column plan. The caller retains/gathers its original payload buffers;
//! this module never reconstructs pandas objects or changes Python dispatch.
use super::SourceError;
use crate::acmg_variant_mechanism::variant_mechanism_masks::Cell;
use polars_core::prelude::*;
use std::collections::{HashMap, HashSet};

#[derive(Clone, Debug)]
pub struct InputColumn {
    pub name: String,
    pub dtype: String,
    pub values: Vec<Cell>,
}
#[derive(Clone, Debug)]
pub struct Input<I> {
    pub columns: Vec<InputColumn>,
    pub index: Vec<I>,
    /// Opaque caller-owned axis/attrs descriptor, echoed without interpretation.
    pub axis_metadata: String,
}
pub struct Ordering<I> {
    pub source_positions: Vec<usize>,
    pub index: Vec<I>,
    /// In INPUT order, retaining source integer/float and missing-value kinds.
    pub maxima: Vec<Cell>,
    pub maximum_dtype: String,
    /// In OUTPUT order, missing variant IDs have rank zero.
    pub ranks: Vec<i64>,
    pub output_columns: Vec<(String, String)>,
    pub temporary_column_position: usize,
    pub axis_metadata: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
enum Atom {
    Text(String),
    Int(i64),
    UInt(u64),
    Float(u64),
    Bool(bool),
}
fn unsupported(message: impl Into<String>) -> SourceError {
    SourceError(
        "UnsupportedInput",
        format!("{}; retain Python backend", message.into()),
    )
}
fn missing(value: &Cell) -> bool {
    matches!(value, Cell::Null) || matches!(value, Cell::Float(x) if x.is_nan())
}
fn atom(column: &InputColumn, value: &Cell) -> Result<Option<Atom>, SourceError> {
    if missing(value) {
        let valid = match (column.dtype.as_str(), value) {
            ("object", Cell::Null | Cell::Float(_)) => true,
            ("float64", Cell::Float(_)) => true,
            ("Float64" | "string" | "Int64" | "UInt64", Cell::Null) => true,
            _ => false,
        };
        if valid {
            return Ok(None);
        }
        return Err(unsupported(format!(
            "missing cell in {} dtype {}",
            column.name, column.dtype
        )));
    }
    Ok(Some(match (column.dtype.as_str(), value) {
        ("object" | "string", Cell::Text(x)) => Atom::Text(x.clone()),
        ("int64" | "Int64", Cell::Integer(x)) => Atom::Int(*x),
        ("uint64" | "UInt64", Cell::Unsigned(x)) => Atom::UInt(*x),
        ("uint64" | "UInt64", Cell::Integer(x)) if *x >= 0 => Atom::UInt(*x as u64),
        ("float64" | "Float64", Cell::Float(x)) => {
            Atom::Float(if *x == 0. { 0 } else { x.to_bits() })
        }
        ("bool", Cell::Bool(x)) => Atom::Bool(*x),
        _ => {
            return Err(unsupported(format!(
                "{} dtype/cell combination {}",
                column.name, column.dtype
            )));
        }
    }))
}
fn polars_column(column: &InputColumn) -> Result<Column, SourceError> {
    let values = column
        .values
        .iter()
        .map(|value| atom(column, value))
        .collect::<Result<Vec<_>, _>>()?;
    let name = column.name.clone().into();
    Ok(match column.dtype.as_str() {
        "object" | "string" => Column::new(
            name,
            values
                .into_iter()
                .map(|v| match v {
                    Some(Atom::Text(x)) => Some(x),
                    None => None,
                    _ => unreachable!(),
                })
                .collect::<Vec<_>>(),
        ),
        "int64" | "Int64" => Column::new(
            name,
            values
                .into_iter()
                .map(|v| match v {
                    Some(Atom::Int(x)) => Some(x),
                    None => None,
                    _ => unreachable!(),
                })
                .collect::<Vec<_>>(),
        ),
        "uint64" | "UInt64" => Column::new(
            name,
            values
                .into_iter()
                .map(|v| match v {
                    Some(Atom::UInt(x)) => Some(x),
                    None => None,
                    _ => unreachable!(),
                })
                .collect::<Vec<_>>(),
        ),
        "float64" | "Float64" => Column::new(
            name,
            values
                .into_iter()
                .map(|v| match v {
                    Some(Atom::Float(x)) => Some(f64::from_bits(x)),
                    None => None,
                    _ => unreachable!(),
                })
                .collect::<Vec<_>>(),
        ),
        "bool" => Column::new(
            name,
            values
                .into_iter()
                .map(|v| match v {
                    Some(Atom::Bool(x)) => Some(x),
                    None => None,
                    _ => unreachable!(),
                })
                .collect::<Vec<_>>(),
        ),
        _ => return Err(unsupported(format!("dtype {}", column.dtype))),
    })
}
fn require<'a, I>(input: &'a Input<I>, name: &str) -> Result<&'a InputColumn, SourceError> {
    input
        .columns
        .iter()
        .find(|column| column.name == name)
        .ok_or_else(|| SourceError("KeyError", format!("'{name}'")))
}

/// The supported absent-ID branch uses the actual colon string join when its
/// components have exact text/integer/bool spellings. Floating components use
/// injective typed tokens only when no string component contains ':'. Their
/// first-appearance equality is identical to pandas float astype(str), including
/// signed zero and NaN spelling collapse. No surrogate string is emitted.
/// Float components plus colon-bearing strings are explicitly unsupported.
fn absent_key(columns: &[&InputColumn], row: usize, has_float: bool) -> String {
    columns
        .iter()
        .map(|column| match &column.values[row] {
            Cell::Text(x) => x.clone(),
            Cell::Integer(x) => x.to_string(),
            Cell::Unsigned(x) => x.to_string(),
            Cell::Bool(x) => if *x { "True" } else { "False" }.into(),
            Cell::Null => match column.dtype.as_str() {
                "object" => "None",
                "float64" => "nan",
                _ => "<NA>",
            }
            .into(),
            Cell::Float(x) if x.is_nan() => "nan".into(),
            Cell::Float(x) => {
                debug_assert!(has_float);
                format!("F{:016x}", x.to_bits())
            }
        })
        .collect::<Vec<_>>()
        .join(":")
}

pub fn terminal_order<I: Clone>(input: &Input<I>) -> Result<Ordering<I>, SourceError> {
    // Preserve source lookup order before dtype validation. Upstream failures
    // belong to Python; this boundary begins at original line668 only.
    let coordinates = [
        require(input, "chrom")?,
        require(input, "pos")?,
        require(input, "ref")?,
        require(input, "alt")?,
    ];
    let scores = input
        .columns
        .iter()
        .find(|c| c.name == "sort_index")
        .ok_or_else(|| SourceError("KeyError", "'Column not found: sort_index'".into()))?;
    let n = input.index.len();
    let mut names = HashSet::new();
    for column in &input.columns {
        if !names.insert(&column.name) || column.values.len() != n {
            return Err(unsupported("duplicate columns or inconsistent lengths"));
        }
    }
    if !matches!(
        scores.dtype.as_str(),
        "float64" | "Float64" | "int64" | "Int64" | "uint64" | "UInt64"
    ) {
        return Err(unsupported(
            "sort_index must have a supported primitive 64-bit numeric dtype",
        ));
    }
    // Typed coordinate identity is independent of index labels and variant IDs.
    // pandas drops every group with any missing coordinate; Polars receives only
    // valid coordinates, represented by exact source-order group IDs.
    let mut group_ids = HashMap::<Vec<Atom>, u64>::new();
    let mut rows = Vec::with_capacity(n);
    let mut groups = Vec::new();
    let mut observations = Vec::new();
    let mut first_zero = HashMap::new();
    for row in 0..n {
        let key = coordinates
            .iter()
            .map(|c| atom(c, &c.values[row]))
            .collect::<Result<Vec<_>, _>>()?;
        atom(scores, &scores.values[row])?;
        if key.iter().any(Option::is_none) {
            rows.push(None);
            continue;
        }
        let next = group_ids.len() as u64;
        let id = *group_ids
            .entry(key.into_iter().map(Option::unwrap).collect())
            .or_insert(next);
        rows.push(Some(id));
        groups.push(id);
        // Keep the original zero sign for first tied maximum restoration below.
        let original_score = match &scores.values[row] {
            Cell::Float(x) if !x.is_nan() => Some(*x),
            _ => None,
        };
        if original_score == Some(0.) {
            first_zero.entry(id).or_insert(original_score.unwrap());
        }
        observations.push(scores.values[row].clone());
    }
    let grouped = DataFrame::new(vec![
        Column::new("group".into(), groups),
        polars_column(&InputColumn {
            name: "score".into(),
            dtype: scores.dtype.clone(),
            values: observations,
        })?,
    ])
    .map_err(|e| SourceError("PolarsError", e.to_string()))?;
    #[allow(deprecated)]
    let maxima = grouped
        .group_by(["group"])
        .and_then(|g| g.select(["score"]).max())
        .map_err(|e| SourceError("PolarsError", e.to_string()))?;
    let mut maximum_by_group = HashMap::new();
    let ids = maxima
        .column("group")
        .and_then(|c| c.u64())
        .map_err(|e| SourceError("PolarsError", e.to_string()))?;
    let values = maxima
        .column("score_max")
        .map_err(|e| SourceError("PolarsError", e.to_string()))?;
    for (position, id) in ids.into_no_null_iter().enumerate() {
        let value = match values
            .get(position)
            .map_err(|e| SourceError("PolarsError", e.to_string()))?
        {
            AnyValue::Float64(x) => Cell::Float(if x == 0. { first_zero[&id] } else { x }),
            AnyValue::Int64(x) => Cell::Integer(x),
            AnyValue::UInt64(x) => Cell::Unsigned(x),
            AnyValue::Null => Cell::Null,
            _ => return Err(unsupported("unexpected maximum dtype")),
        };
        maximum_by_group.insert(id, value);
    }
    // Actual pandas2.2.3 capture: nullable integer transform preserves its mask
    // and native dtype. NumPy integer transform promotes to float64 ONLY when
    // excluded coordinate groups require a missing result. Preserve that source
    // promotion after exact integer reduction; never cast before taking max.
    let promoted =
        matches!(scores.dtype.as_str(), "int64" | "uint64") && rows.iter().any(Option::is_none);
    let maximum_dtype = if promoted {
        "float64".into()
    } else {
        scores.dtype.clone()
    };
    let max_scores: Vec<Cell> = rows
        .iter()
        .map(|id| {
            let value = id
                .map(|id| maximum_by_group[&id].clone())
                .unwrap_or(Cell::Null);
            match value {
                Cell::Integer(x) if promoted => Cell::Float(x as f64),
                Cell::Unsigned(x) if promoted => Cell::Float(x as f64),
                Cell::Null if maximum_dtype == "float64" => Cell::Float(f64::NAN),
                value => value,
            }
        })
        .collect();
    // This lookup occurs only after grouping, as in original sort_values.
    let raw = require(input, "ACMG_quant_score")?;
    let mut sort_columns = vec![
        polars_column(&InputColumn {
            name: "max_variant_score".into(),
            dtype: maximum_dtype.clone(),
            values: max_scores.clone(),
        })?,
        polars_column(raw)?,
    ];
    for column in coordinates {
        sort_columns.push(polars_column(column)?);
    }
    sort_columns.push(Column::new(
        "source_position".into(),
        (0..n as u64).collect::<Vec<_>>(),
    ));
    let frame =
        DataFrame::new(sort_columns).map_err(|e| SourceError("PolarsError", e.to_string()))?;
    let sorted = frame
        .sort(
            [
                "max_variant_score",
                "ACMG_quant_score",
                "chrom",
                "pos",
                "ref",
                "alt",
                "source_position",
            ],
            SortMultipleOptions::default()
                .with_order_descending_multi([true, true, false, false, false, false, false])
                .with_nulls_last(true)
                .with_maintain_order(true),
        )
        .map_err(|e| SourceError("PolarsError", e.to_string()))?;
    let positions = sorted
        .column("source_position")
        .and_then(|c| c.u64())
        .map_err(|e| SourceError("PolarsError", e.to_string()))?
        .into_no_null_iter()
        .map(|x| x as usize)
        .collect::<Vec<_>>();
    let ids = input.columns.iter().find(|c| c.name == "variant_id");
    if ids.is_some_and(|c| {
        !matches!(
            c.dtype.as_str(),
            "object"
                | "string"
                | "int64"
                | "Int64"
                | "uint64"
                | "UInt64"
                | "float64"
                | "Float64"
                | "bool"
        )
    }) {
        return Err(unsupported("variant_id dtype"));
    }
    let has_float = coordinates
        .iter()
        .any(|c| matches!(c.dtype.as_str(), "float64" | "Float64"));
    if ids.is_none()
        && has_float
        && coordinates.iter().any(|c| {
            c.values
                .iter()
                .any(|v| matches!(v, Cell::Text(x) if x.contains(':')))
        })
    {
        return Err(unsupported(
            "absent variant_id with floating coordinates and colon-bearing text",
        ));
    }
    let mut rank_by_id = HashMap::new();
    let mut ranks = Vec::with_capacity(n);
    for &position in &positions {
        let key = if let Some(column) = ids {
            atom(column, &column.values[position])?
        } else {
            Some(Atom::Text(absent_key(&coordinates, position, has_float)))
        };
        let rank = if let Some(key) = key {
            let next = rank_by_id.len() as i64 + 1;
            *rank_by_id.entry(key).or_insert(next)
        } else {
            0
        };
        ranks.push(rank);
    }
    let mut output_columns = input
        .columns
        .iter()
        .filter(|c| c.name != "max_variant_score")
        .map(|c| {
            (
                c.name.clone(),
                if c.name == "variant_rank" {
                    "int64".into()
                } else {
                    c.dtype.clone()
                },
            )
        })
        .collect::<Vec<_>>();
    if !input.columns.iter().any(|c| c.name == "variant_rank") {
        output_columns.push(("variant_rank".into(), "int64".into()));
    }
    Ok(Ordering {
        index: positions.iter().map(|&i| input.index[i].clone()).collect(),
        source_positions: positions,
        maxima: max_scores,
        maximum_dtype,
        ranks,
        output_columns,
        temporary_column_position: input
            .columns
            .iter()
            .position(|c| c.name == "max_variant_score")
            .unwrap_or(input.columns.len()),
        axis_metadata: input.axis_metadata.clone(),
    })
}
