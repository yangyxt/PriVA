//! Original sort_and_rank_variants R1, 09223ef lines492–498 only.
//! Return a copy/drop/column-copy/filter plan over caller-owned primitive buffers.
//! No PED, score arithmetic, numerical package or later ranking stage is called.
use super::SourceError;
use super::terminal_order::Input;
use crate::acmg_variant_mechanism::variant_mechanism_masks::Cell;
use polars_core::prelude::*;

pub struct ColumnCopy {
    pub name: String,
    pub dtype: String,
    pub source_column: usize,
    pub copy_score: bool,
}
pub struct Setup<I> {
    pub source_positions: Vec<usize>,
    pub index: Vec<I>,
    pub columns: Vec<ColumnCopy>,
    /// Original column positions, grouped in the original legacy-drop order.
    pub dropped_columns: Vec<usize>,
    pub axis_metadata: String,
}
fn unsupported(message: impl Into<String>) -> SourceError {
    SourceError(
        "UnsupportedInput",
        format!("{}; retain Python backend", message.into()),
    )
}
fn valid_cell(dtype: &str, value: &Cell) -> bool {
    match (dtype, value) {
        ("object", _) => true,
        ("string", Cell::Text(_) | Cell::Null) => true,
        ("float64", Cell::Float(_)) => true,
        ("Float64", Cell::Float(x)) => !x.is_nan(),
        ("Float64", Cell::Null) => true,
        ("int64" | "Int64", Cell::Integer(_)) => true,
        ("Int64" | "UInt64" | "boolean", Cell::Null) => true,
        ("uint64" | "UInt64", Cell::Unsigned(_)) => true,
        ("uint64" | "UInt64", Cell::Integer(x)) => *x >= 0,
        ("bool" | "boolean", Cell::Bool(_)) => true,
        _ => false,
    }
}

pub fn ranking_setup<I: Clone>(input: &Input<I>) -> Result<Setup<I>, SourceError> {
    // The source copies first and removes legacy fields in this fixed order.
    // They cannot overlap a required label. Input buffers are only borrowed.
    let mut dropped = Vec::new();
    for stale in ["haplo_insuf_index", "zygosity_inheritance_mechanism_index"] {
        dropped.extend(
            input
                .columns
                .iter()
                .enumerate()
                .filter(|(_, c)| c.name == stale)
                .map(|(i, _)| i),
        );
    }
    let find = |name: &str| {
        input
            .columns
            .iter()
            .enumerate()
            .filter(|(_, c)| c.name == name)
            .map(|(i, _)| i)
            .collect::<Vec<_>>()
    };
    let acmg = find("ACMG_quant_score");
    if acmg.is_empty() {
        return Err(SourceError("KeyError", "'ACMG_quant_score'".into()));
    }
    if acmg.len() != 1 {
        return Err(unsupported("duplicate ACMG_quant_score label"));
    }
    let previous_sort = find("sort_index");
    if previous_sort.len() > 1 {
        return Err(unsupported("duplicate sort_index label"));
    }
    // Original score-column assignment precedes the BIOTYPE lookup, including
    // empty/all-filtered frames. This plan retains the source failure priority.
    let biotype = find("BIOTYPE");
    if biotype.is_empty() {
        return Err(SourceError("KeyError", "'BIOTYPE'".into()));
    }
    if biotype.len() != 1 {
        return Err(unsupported("duplicate BIOTYPE label"));
    }
    let n = input.index.len();
    for (position, column) in input.columns.iter().enumerate() {
        if column.values.len() != n {
            return Err(unsupported("inconsistent column/index lengths"));
        }
        if dropped.contains(&position) {
            continue;
        }
        if !matches!(
            column.dtype.as_str(),
            "object"
                | "string"
                | "float64"
                | "Float64"
                | "int64"
                | "Int64"
                | "uint64"
                | "UInt64"
                | "bool"
                | "boolean"
        ) || column
            .values
            .iter()
            .any(|value| !valid_cell(&column.dtype, value))
        {
            return Err(unsupported(format!(
                "{} dtype/cell combination {}",
                column.name, column.dtype
            )));
        }
    }
    let source_score = acmg[0];
    let mut columns = Vec::new();
    for (position, column) in input.columns.iter().enumerate() {
        if dropped.contains(&position) {
            continue;
        }
        let copy_score = column.name == "sort_index";
        columns.push(ColumnCopy {
            name: column.name.clone(),
            dtype: if copy_score {
                input.columns[source_score].dtype.clone()
            } else {
                column.dtype.clone()
            },
            source_column: if copy_score { source_score } else { position },
            copy_score,
        });
    }
    if previous_sort.is_empty() {
        columns.push(ColumnCopy {
            name: "sort_index".into(),
            dtype: input.columns[source_score].dtype.clone(),
            source_column: source_score,
            copy_score: true,
        });
    }
    // Equality is exact and case-sensitive. Every non-string primitive and
    // every supported missing value fails this text comparison; nullable mask
    // missing entries are excluded by pandas' original boolean indexing.
    let keep = input.columns[biotype[0]]
        .values
        .iter()
        .map(|value| matches!(value, Cell::Text(value) if value == "protein_coding"))
        .collect::<Vec<_>>();
    let frame = DataFrame::new(vec![
        Column::new("source_position".into(), (0..n as u64).collect::<Vec<_>>()),
        Column::new("keep".into(), keep),
    ])
    .map_err(|e| SourceError("PolarsError", e.to_string()))?;
    let retained = frame
        .filter(
            frame
                .column("keep")
                .and_then(|c| c.bool())
                .map_err(|e| SourceError("PolarsError", e.to_string()))?,
        )
        .map_err(|e| SourceError("PolarsError", e.to_string()))?;
    let positions = retained
        .column("source_position")
        .and_then(|c| c.u64())
        .map_err(|e| SourceError("PolarsError", e.to_string()))?
        .into_no_null_iter()
        .map(|value| value as usize)
        .collect::<Vec<_>>();
    Ok(Setup {
        index: positions.iter().map(|&i| input.index[i].clone()).collect(),
        source_positions: positions,
        columns,
        dropped_columns: dropped,
        axis_metadata: input.axis_metadata.clone(),
    })
}
