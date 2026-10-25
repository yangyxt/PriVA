//! Original R5 then R6, L599-616, after the unchanged C1-C5/R2/R3/R4 stages.
//! Independent typed procedure; caller retains source buffers/metadata.
use super::SourceError;
use super::control_multiply::RowIndex;
use super::ranking_compatibility::Compatibility;
use super::terminal_order::{Input, InputColumn};
use crate::acmg_variant_mechanism::variant_mechanism_masks::Cell;
use polars_core::prelude::*;
use std::collections::HashSet;

pub const DISPENSABLE: &str = "dispensable_gene_index";
pub const RELEVANT: &str = "relevant_gene_index";
pub const COMPATIBILITY: &str = "zygosity_inheritance_mechanism_compatibility";

/// None disables the source branch; an empty selected file enables it.
pub enum GeneList<'a> {
    Disabled,
    PlainText(&'a str),
}

pub struct Penalty<I> {
    pub index: Vec<I>,
    pub dtype: &'static str,
    /// Actual raw Boolean data and separate NA bits, both retained unchanged.
    pub data: Vec<bool>,
    pub null: Vec<bool>,
}
impl<I: Clone> Penalty<I> {
    /// Compose directly with the previously computed R4 result; no R2 replay.
    pub fn from_compatibility(value: &Compatibility<I>) -> Result<Self, SourceError> {
        let bits = |role| {
            value
                .role(role)
                .and_then(Column::bool)
                .map(|column| column.into_no_null_iter().collect())
                .map_err(|_| unsupported("R4 raw mask representation"))
        };
        Ok(Self {
            index: value.index.clone(),
            dtype: value.penalty_dtype,
            data: bits("compatibility_penalty_rows")?,
            null: bits("compatibility_penalty_rows_null")?,
        })
    }
}

pub struct Step {
    /// The literal source branch/statement, in execution order.
    pub name: &'static str,
    pub source_line: u32,
    pub multiplier_column_position: usize,
    /// selected, multiplier and sort_index in original row order.
    pub frame: DataFrame,
}
pub struct Weights<I> {
    pub steps: Vec<Step>,
    /// File-order strings, before Python's set construction; duplicates remain.
    pub dispensable_genes: Option<Vec<String>>,
    pub relevant_genes: Option<Vec<String>>,
    pub score_column_position: usize,
    pub index: Vec<I>,
    pub axis_metadata: String,
}

fn unsupported(message: &str) -> SourceError {
    SourceError(
        "UnsupportedInput",
        format!("{message}; retain Python backend"),
    )
}

/// Only a literal plain-text subset of pandas read_table is qualified here.
/// All inference-sensitive tokens, quoting, tabs, whitespace, non-ASCII and
/// compression/encoding routes retain the actual original pandas parser.
pub fn parse_gene_text(text: &str) -> Result<Vec<String>, SourceError> {
    let mut genes = Vec::new();
    for line in text.split('\n') {
        if !line.is_ascii() || line.contains(['\r', '\t', '"', '\0']) {
            return Err(unsupported("gene-file parser/encoding route"));
        }
        let field = line.split_once('#').map_or(line, |(left, _)| left);
        if field.is_empty() {
            continue;
        }
        if !field.as_bytes()[0].is_ascii_alphabetic() && !field.starts_with('_')
            || !field
                .bytes()
                .all(|byte| byte.is_ascii_alphanumeric() || b"_.-".contains(&byte))
            || [
                "NA", "N/A", "NaN", "nan", "NULL", "null", "None", "n/a", "True", "False", "TRUE",
                "FALSE", "true", "false", "inf", "Inf", "INF", "infinity", "Infinity", "INFINITY",
            ]
            .contains(&field)
            || ["inf", "infinity", "nan", "true", "false", "null", "none"]
                .iter()
                .any(|token| field.eq_ignore_ascii_case(token))
        {
            return Err(unsupported("gene-file dtype/NA inference or token route"));
        }
        genes.push(field.to_owned());
    }
    Ok(genes)
}

fn position<I>(frame: &Input<I>, name: &str) -> Result<usize, SourceError> {
    let positions: Vec<_> = frame
        .columns
        .iter()
        .enumerate()
        .filter(|(_, column)| column.name == name)
        .map(|(i, _)| i)
        .collect();
    match positions.as_slice() {
        [i] => Ok(*i),
        _ => Err(unsupported("missing/ambiguous R5-R6 consumed label")),
    }
}
fn floats(column: &InputColumn) -> Result<Vec<f64>, SourceError> {
    if column.dtype != "float64" {
        return Err(unsupported(
            "native float64 multiplication/assignment required",
        ));
    }
    column
        .values
        .iter()
        .map(|cell| match cell {
            Cell::Float(value) => Ok(*value),
            _ => Err(unsupported("native float64 storage required")),
        })
        .collect()
}
fn set_column<I>(frame: &mut Input<I>, name: &str, values: &[f64]) -> Result<usize, SourceError> {
    let mut found = frame
        .columns
        .iter()
        .enumerate()
        .filter(|(_, c)| c.name == name);
    let first = found.next().map(|(i, _)| i);
    if found.next().is_some() {
        return Err(unsupported("duplicate assigned gene-weight label"));
    }
    let column = InputColumn {
        name: name.to_owned(),
        dtype: "float64".into(),
        values: values.iter().copied().map(Cell::Float).collect(),
    };
    if let Some(i) = first {
        frame.columns[i] = column;
        Ok(i)
    } else {
        let i = frame.columns.len();
        frame.columns.push(column);
        Ok(i)
    }
}
fn multiply(
    scores: &[f64],
    multipliers: &[f64],
    selected: &[bool],
) -> Result<Vec<f64>, SourceError> {
    scores
        .iter()
        .zip(multipliers)
        .zip(selected)
        .map(|((&left, &right), &use_row)| {
            if !use_row {
                Ok(left)
            } else if left.is_nan() && right.is_nan() {
                Err(unsupported("two-NaN operand payload route"))
            } else {
                Ok(left * right)
            }
        })
        .collect()
}

/// Preflight and plan the complete block before changing current R5-entry
/// buffers. A rejection must run only original L599-616 on those unchanged
/// buffers; it must never replay prior stages or a multiplier already applied.
/// Qualified: <=1M rows, flat opaque axes, aligned R4 Boolean mask, native f64
/// scores/compatibility and literal text-file subset. No arbitrary objects,
/// dtype promotion, numexpr-size routes or parser/error approximation.
pub fn ranking_gene_weights<I: Clone + PartialEq>(
    frame: &mut Input<I>,
    penalty: &Penalty<I>,
    dispensable: GeneList<'_>,
    relevant: GeneList<'_>,
    row_index: RowIndex,
) -> Result<Weights<I>, SourceError> {
    let n = frame.index.len();
    if row_index != RowIndex::Flat
        || n > 1_000_000
        || penalty.index != frame.index
        || penalty.data.len() != n
        || penalty.null.len() != n
        || !matches!(penalty.dtype, "bool" | "boolean") && !(penalty.dtype == "object" && n == 0)
        || penalty.dtype != "boolean" && penalty.null.iter().any(|value| *value)
        || frame.columns.iter().any(|column| column.values.len() != n)
    {
        return Err(unsupported(
            "unaligned R4 mask/axis/storage or expression-size route",
        ));
    }
    let mut planned = frame.clone();
    let mut steps = Vec::new();
    let mut gene_lists = Vec::new();
    for (selection, name, factor, line) in [
        (dispensable, DISPENSABLE, 0.9, 604),
        (relevant, RELEVANT, 1.2, 611),
    ] {
        let GeneList::PlainText(text) = selection else {
            gene_lists.push(None);
            continue;
        };
        let genes = parse_gene_text(text)?;
        let set: HashSet<_> = genes.iter().collect();
        let symbols = &planned.columns[position(&planned, "SYMBOL")?];
        if !matches!(
            symbols.dtype.as_str(),
            "object" | "string" | "int64" | "uint64" | "float64" | "bool"
        ) {
            return Err(unsupported("SYMBOL dtype/membership representation"));
        }
        let selected = symbols
            .values
            .iter()
            .map(|cell| match (symbols.dtype.as_str(), cell) {
                ("object" | "string", Cell::Text(value)) => Ok(set.contains(value)),
                (
                    "object",
                    Cell::Null
                    | Cell::Integer(_)
                    | Cell::Unsigned(_)
                    | Cell::Float(_)
                    | Cell::Bool(_),
                )
                | ("string", Cell::Null)
                | ("int64", Cell::Integer(_))
                | ("uint64", Cell::Unsigned(_))
                | ("float64", Cell::Float(_))
                | ("bool", Cell::Bool(_)) => Ok(false),
                _ => Err(unsupported("SYMBOL object equality route")),
            })
            .collect::<Result<Vec<_>, _>>()?;
        let multipliers: Vec<_> = selected
            .iter()
            .map(|yes| if *yes { factor } else { 1.0 })
            .collect();
        let multiplier_position = set_column(&mut planned, name, &multipliers)?;
        let score_position = position(&planned, "sort_index")?;
        let scores = multiply(
            &floats(&planned.columns[score_position])?,
            &multipliers,
            &vec![true; n],
        )?;
        planned.columns[score_position].values = scores.iter().copied().map(Cell::Float).collect();
        steps.push(Step {
            name,
            source_line: line,
            multiplier_column_position: multiplier_position,
            frame: DataFrame::new(vec![
                Column::new("selected".into(), selected),
                Column::new("multiplier".into(), multipliers),
                Column::new("sort_index".into(), scores),
            ])
            .map_err(|_| unsupported("gene-step Polars storage"))?,
        });
        gene_lists.push(Some(genes));
    }
    // Literal R6 follows both gene branches. NA mask entries select no rows.
    let target = position(&planned, "sort_index")?;
    let compatibility = position(&planned, COMPATIBILITY)?;
    let multipliers = floats(&planned.columns[compatibility])?;
    let selected: Vec<_> = penalty
        .data
        .iter()
        .zip(&penalty.null)
        .map(|(&data, &null)| data && !null)
        .collect();
    let scores = multiply(&floats(&planned.columns[target])?, &multipliers, &selected)?;
    planned.columns[target].values = scores.iter().copied().map(Cell::Float).collect();
    steps.push(Step {
        name: COMPATIBILITY,
        source_line: 613,
        multiplier_column_position: compatibility,
        frame: DataFrame::new(vec![
            Column::new("selected".into(), selected),
            Column::new("multiplier".into(), multipliers),
            Column::new("sort_index".into(), scores),
        ])
        .map_err(|_| unsupported("compatibility-step Polars storage"))?,
    });
    frame.columns = planned.columns;
    let mut gene_lists = gene_lists.into_iter();
    Ok(Weights {
        steps,
        dispensable_genes: gene_lists.next().unwrap(),
        relevant_genes: gene_lists.next().unwrap(),
        score_column_position: target,
        index: frame.index.clone(),
        axis_metadata: frame.axis_metadata.clone(),
    })
}
