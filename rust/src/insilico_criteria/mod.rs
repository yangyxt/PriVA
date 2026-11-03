//! Column-wide PP3/BP4/BP7 calculation primitives for the exact b8a093a9 source.
//!
//! The host retains the original statement/getter/logger frontiers. Every
//! admitted table coercion, string predicate, mask, reduction and criterion
//! assignment runs here; a session has a bounded number of requests per source
//! procedure, never a process/request per variant. Unsupported mature types
//! remain explicit Python operations, not falsely reported native execution.
mod procedure_transport;
pub mod procedures;
pub(crate) mod transport;

use polars_core::chunked_array::ops::arity;
use polars_core::prelude::*;
use std::collections::{HashMap, HashSet};

pub type Result<T> = std::result::Result<T, String>;

#[derive(Clone, Debug, PartialEq)]
pub enum Cell {
    None,
    Na,
    Bool(bool),
    Int(i64),
    Float(f64),
    Text(String),
}

impl Cell {
    fn missing(&self) -> bool {
        matches!(self, Self::None | Self::Na)
            || matches!(self, Self::Float(value) if value.is_nan())
    }

    fn number(&self) -> Option<f64> {
        match self {
            Self::Bool(value) => Some(f64::from(u8::from(*value))),
            Self::Int(value) => Some(*value as f64),
            Self::Float(value) => Some(*value),
            _ => None,
        }
    }

    fn truth(&self) -> Result<Option<bool>> {
        Ok(match self {
            Self::Na => None,
            Self::None => Some(false),
            Self::Bool(value) => Some(*value),
            Self::Int(value) => Some(*value != 0),
            Self::Float(value) => Some(*value != 0.0),
            Self::Text(value) => Some(!value.is_empty()),
        })
    }
}

/// A typed whole column, ndarray, list, or scalar. Axis keys are lossless
/// host-encoded ordinary labels, not numeric positions masquerading as labels.
#[derive(Clone, Debug, PartialEq)]
pub struct Value {
    pub kind: String,
    pub dtype: String,
    pub index: Vec<String>,
    pub cells: Vec<Cell>,
}

impl Value {
    pub fn scalar(cell: Cell) -> Self {
        Self {
            kind: "scalar".into(),
            dtype: "scalar".into(),
            index: Vec::new(),
            cells: vec![cell],
        }
    }

    fn column(&self, cells: Vec<Cell>, dtype: &str) -> Result<Self> {
        if self.kind != "scalar" && self.cells.len() != cells.len() {
            return Err("column result length differs".into());
        }
        Ok(Self {
            kind: self.kind.clone(),
            dtype: dtype.into(),
            index: self.index.clone(),
            cells,
        })
    }

    fn text(&self) -> Result<&str> {
        match self.cells.as_slice() {
            [Cell::Text(value)] if self.kind == "scalar" => Ok(value),
            _ => Err("expected scalar text".into()),
        }
    }
}

fn numeric_column(cells: &[Cell]) -> Result<Float64Chunked> {
    let values = cells
        .iter()
        .map(|cell| {
            if cell.missing() {
                Ok(None)
            } else {
                cell.number()
                    .map(Some)
                    .ok_or_else(|| "unsupported numeric operand".into())
            }
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(Float64Chunked::from_iter_options(
        "numeric".into(),
        values.into_iter(),
    ))
}

fn boolean_column(cells: &[Cell]) -> Result<BooleanChunked> {
    let values = cells
        .iter()
        .map(|cell| match cell {
            Cell::Bool(value) => Ok(Some(*value)),
            Cell::None | Cell::Na => Ok(None),
            _ => Err("unsupported nonboolean operand".into()),
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(BooleanChunked::from_iter_options(
        "boolean".into(),
        values.into_iter(),
    ))
}

fn string_column(cells: &[Cell]) -> Result<StringChunked> {
    let values = cells
        .iter()
        .map(|cell| match cell {
            Cell::Text(value) => Ok(Some(value.as_str())),
            Cell::None | Cell::Na => Ok(None),
            _ => Err("unsupported nonstring operand".into()),
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(StringChunked::from_iter_options(
        "string".into(),
        values.into_iter(),
    ))
}

// These are explicit pandas null/dtype transport adapters AFTER a real typed
// Polars operation. Pandas ordinary float NaN comparisons are False (or True
// for !=), while nullable extension comparisons retain NA; Polars uses nulls.
fn boolean_cells(column: &BooleanChunked, nullable: bool, missing_default: bool) -> Vec<Cell> {
    column
        .into_iter()
        .map(|value| match value {
            Some(value) => Cell::Bool(value),
            None if nullable => Cell::Na,
            None => Cell::Bool(missing_default),
        })
        .collect()
}

fn numeric_cells(column: &Float64Chunked, nullable: bool) -> Vec<Cell> {
    column
        .into_iter()
        .map(|value| match value {
            Some(value) => Cell::Float(value),
            None if nullable => Cell::Na,
            None => Cell::Float(f64::NAN),
        })
        .collect()
}

fn bool_dtype(value: &Value) -> &'static str {
    if matches!(
        value.dtype.as_str(),
        "boolean" | "Int64" | "Float64" | "string"
    ) {
        "boolean"
    } else {
        "bool"
    }
}

fn string(cell: &Cell) -> Result<String> {
    Ok(match cell {
        Cell::None => "None".into(),
        Cell::Na => "<NA>".into(),
        Cell::Bool(true) => "True".into(),
        Cell::Bool(false) => "False".into(),
        Cell::Int(value) => value.to_string(),
        Cell::Float(value) if value.is_nan() => "nan".into(),
        Cell::Float(value) if *value == f64::INFINITY => "inf".into(),
        Cell::Float(value) if *value == f64::NEG_INFINITY => "-inf".into(),
        // Exact stringification of arbitrary floating object scalars is a
        // mature fallback. For ordinary source score ranges the shortest
        // decimal plus explicit .0 is identical and never rounded for scores.
        Cell::Float(value) if value.abs() < 1e16 && (*value == 0.0 || value.abs() >= 1e-4) => {
            let mut text = value.to_string();
            if !text.contains('.') {
                text.push_str(".0");
            }
            text
        }
        Cell::Float(_) => return Err("unsupported float stringification".into()),
        Cell::Text(value) => value.clone(),
    })
}

fn python_whitespace(value: char) -> bool {
    matches!(value as u32, 0x09..=0x0d | 0x1c..=0x20 | 0x85 | 0xa0 | 0x1680
        | 0x2000..=0x200a | 0x2028..=0x2029 | 0x202f | 0x205f | 0x3000)
}

/// pandas.to_numeric(errors='coerce') fixed lexical domain. Long significant
/// decimals are deliberately unsupported instead of substituting Rust's parser
/// for the pinned pandas parser's sometimes different rounding.
fn numeric(cell: &Cell) -> Result<Cell> {
    if cell.missing() {
        return Ok(Cell::Float(f64::NAN));
    }
    if let Cell::Text(value) = cell {
        let value = value.trim_matches(|c: char| c.is_ascii_whitespace());
        let digits = value.chars().filter(char::is_ascii_digit).count();
        if digits > 15
            && value
                .chars()
                .all(|c| c.is_ascii_digit() || "+-.eE".contains(c))
        {
            return Err("unsupported long numeric string".into());
        }
        let lower = value.to_ascii_lowercase();
        if ["inf", "+inf", "infinity", "+infinity"].contains(&lower.as_str()) {
            return Ok(Cell::Float(f64::INFINITY));
        }
        if ["-inf", "-infinity"].contains(&lower.as_str()) {
            return Ok(Cell::Float(f64::NEG_INFINITY));
        }
        if value.is_empty() || value.chars().any(|c| !c.is_ascii()) {
            return Ok(Cell::Float(f64::NAN));
        }
        return Ok(Cell::Float(match value.parse::<f64>() {
            Ok(number) if number.is_finite() => number,
            _ => f64::NAN,
        }));
    }
    cell.number()
        .map(Cell::Float)
        .ok_or_else(|| "unsupported numeric scalar".into())
}

fn aligned(left: &Value, right: &Value) -> Result<(Vec<String>, Vec<Cell>, Vec<Cell>)> {
    if left.kind == "scalar" {
        return Ok((
            right.index.clone(),
            vec![left.cells[0].clone(); right.cells.len()],
            right.cells.clone(),
        ));
    }
    if right.kind == "scalar" {
        return Ok((
            left.index.clone(),
            left.cells.clone(),
            vec![right.cells[0].clone(); left.cells.len()],
        ));
    }
    if left.cells.len() != right.cells.len() {
        return Err("unsupported unequal vector lengths".into());
    }
    if left.kind != "series" || right.kind != "series" {
        return Ok((left.index.clone(), left.cells.clone(), right.cells.clone()));
    }
    let a: HashMap<&str, usize> = left
        .index
        .iter()
        .enumerate()
        .map(|(i, key)| (key.as_str(), i))
        .collect();
    let b: HashMap<&str, usize> = right
        .index
        .iter()
        .enumerate()
        .map(|(i, key)| (key.as_str(), i))
        .collect();
    if a.len() != left.index.len()
        || b.len() != right.index.len()
        || a.len() != b.len()
        || a.keys().any(|key| !b.contains_key(key))
    {
        return Err("unsupported nonmatching/duplicate aligned index".into());
    }
    if left.index == right.index {
        return Ok((left.index.clone(), left.cells.clone(), right.cells.clone()));
    }
    // Wire labels carry a fixed sortable index-rank supplied by the raw host
    // ownership boundary. The native kernel, not pandas, performs alignment.
    let mut keys = left.index.clone();
    keys.sort();
    let left_cells = keys
        .iter()
        .map(|key| left.cells[a[key.as_str()]].clone())
        .collect();
    let right_cells = keys
        .iter()
        .map(|key| right.cells[b[key.as_str()]].clone())
        .collect();
    Ok((keys, left_cells, right_cells))
}

fn binary(op: &str, left: &Value, right: &Value) -> Result<Value> {
    let (index, a, b) = aligned(left, right)?;
    let prototype = if left.kind == "scalar" { right } else { left };
    let nullable = bool_dtype(left) == "boolean" || bool_dtype(right) == "boolean";
    let result = if op == "and" || op == "or" {
        let a = boolean_column(&a)?;
        let b = boolean_column(&b)?;
        if op == "and" { &a & &b } else { &a | &b }
    } else if a
        .iter()
        .chain(&b)
        .all(|cell| cell.missing() || cell.number().is_some())
    {
        let a = numeric_column(&a)?;
        let b = numeric_column(&b)?;
        match op {
            "eq" => a.equal(&b),
            "ne" => a.not_equal(&b),
            "lt" => a.lt(&b),
            "le" => a.lt_eq(&b),
            "gt" => a.gt(&b),
            "ge" => a.gt_eq(&b),
            _ => return Err("unknown numeric comparison".into()),
        }
    } else if a
        .iter()
        .chain(&b)
        .all(|cell| cell.missing() || matches!(cell, Cell::Text(_)))
    {
        let a = string_column(&a)?;
        let b = string_column(&b)?;
        match op {
            "eq" => a.equal(&b),
            "ne" => a.not_equal(&b),
            _ => return Err("unsupported ordered string comparison".into()),
        }
    } else {
        return Err("unsupported mixed object comparison".into());
    };
    let mut output = prototype.column(
        boolean_cells(&result, nullable, op == "ne"),
        if nullable { "boolean" } else { "bool" },
    )?;
    output.index = index;
    Ok(output)
}

/// Complete column-wide primitive dispatch. The phase ordering remains the
/// exact original function ordering, including both independent PVS1 stages.
pub fn calculate(op: &str, args: &[Value]) -> Result<Value> {
    let first = args.first().ok_or("missing operation operand")?;
    if ["and", "or", "eq", "ne", "lt", "le", "gt", "ge"].contains(&op) {
        if args.len() != 2 {
            return Err("binary operation arity".into());
        }
        return binary(op, first, &args[1]);
    }
    match op {
        "numeric" => {
            if matches!(
                first.dtype.as_str(),
                "bool" | "int64" | "float64" | "boolean" | "Int64" | "Float64"
            ) {
                return Ok(first.clone());
            }
            let nullable = matches!(
                first.dtype.as_str(),
                "string" | "Int64" | "Float64" | "boolean"
            );
            let cells = first
                .cells
                .iter()
                .map(numeric)
                .map(|cell| cell.map(|c| if nullable && c.missing() { Cell::Na } else { c }))
                .collect::<Result<Vec<_>>>()?;
            let parsed = numeric_column(&cells)?;
            first.column(
                numeric_cells(&parsed, nullable),
                if nullable { "Float64" } else { "float64" },
            )
        }
        "isna" | "notna" => {
            // Values are an ownership adapter; the actual validity predicate
            // is Polars' column operation, including all-null/empty columns.
            let validity = BooleanChunked::from_iter_options(
                "present".into(),
                first
                    .cells
                    .iter()
                    .map(|cell| if cell.missing() { None } else { Some(true) }),
            );
            let mask = if op == "isna" {
                validity.is_null()
            } else {
                validity.is_not_null()
            };
            first.column(boolean_cells(&mask, false, false), "bool")
        }
        "fillna" => {
            if args.len() != 2 || args[1].kind != "scalar" {
                return Err("fillna scalar required".into());
            }
            let replacement = &args[1].cells[0];
            let cells = match first.dtype.as_str() {
                "bool" | "boolean" => {
                    let Cell::Bool(fill) = replacement else {
                        return Err("unsupported Boolean fill dtype".into());
                    };
                    let filled = boolean_column(&first.cells)?
                        .fill_null_with_values(*fill)
                        .map_err(|error| error.to_string())?;
                    boolean_cells(&filled, first.dtype == "boolean", false)
                }
                "float64" | "Float64" => {
                    let fill = replacement
                        .number()
                        .ok_or("unsupported floating fill dtype")?;
                    let filled = numeric_column(&first.cells)?
                        .fill_null_with_values(fill)
                        .map_err(|error| error.to_string())?;
                    numeric_cells(&filled, first.dtype == "Float64")
                }
                "int64" | "Int64" => {
                    if first.cells.iter().any(Cell::missing) {
                        return Err("unsupported nullable integer fill".into());
                    }
                    first.cells.clone()
                }
                "object" | "string" => {
                    // Object downcasting/warnings remain mature and explicit.
                    // Canonical string fill is a real typed Polars operation.
                    let Cell::Text(fill) = replacement else {
                        return Err("unsupported object fill".into());
                    };
                    let original = string_column(&first.cells)?;
                    let filled: StringChunked =
                        arity::unary_elementwise(&original, |value: Option<&str>| {
                            Some(value.unwrap_or(fill.as_str()))
                        });
                    filled
                        .into_iter()
                        .map(|value| Cell::Text(value.unwrap().into()))
                        .collect()
                }
                _ => return Err("unsupported fill dtype".into()),
            };
            first.column(cells, &first.dtype)
        }
        "as_str" => {
            // Python scalar spelling is an explicit compatibility adapter.
            let values = first.cells.iter().map(string).collect::<Result<Vec<_>>>()?;
            let column = StringChunked::from_iter_values(
                "strings".into(),
                values.iter().map(String::as_str),
            );
            first.column(
                column
                    .into_iter()
                    .map(|value| Cell::Text(value.unwrap().into()))
                    .collect(),
                "object",
            )
        }
        "strip" => {
            let original = string_column(&first.cells)?;
            if original.null_count() != 0 {
                return Err("unsupported string strip null".into());
            }
            let stripped: StringChunked =
                arity::unary_elementwise(&original, |value: Option<&str>| {
                    value.map(|text| text.trim_matches(python_whitespace))
                });
            first.column(
                stripped
                    .into_iter()
                    .map(|value| Cell::Text(value.unwrap().into()))
                    .collect(),
                &first.dtype,
            )
        }
        "contains" => {
            if args.len() != 2 {
                return Err("contains arity".into());
            }
            let pattern = args[1].text()?;
            let original = string_column(&first.cells)?;
            if first.dtype != "string" && original.null_count() != 0 {
                return Err("unsupported object contains null".into());
            }
            let known = pattern == r"(?:^|[&,;|])(?:HC|OS)(?:[&,;|]|$)"
                || [
                    "missense_variant",
                    "stop",
                    "frameshift",
                    "splice",
                    "5_prime_UTR_variant",
                    "pathogenic",
                    "enign",
                    "athogenic",
                ]
                .contains(&pattern);
            if !known {
                return Err("unsupported regex pattern".into());
            }
            let result: BooleanChunked =
                arity::unary_elementwise(&original, |value: Option<&str>| {
                    value.map(|text| {
                        if pattern == r"(?:^|[&,;|])(?:HC|OS)(?:[&,;|]|$)" {
                            // Literal source delimiters; Python $ also matches before
                            // ONE final newline, not before two trailing newlines.
                            let text = text.strip_suffix('\n').unwrap_or(text);
                            text.split(['&', ',', ';', '|'])
                                .any(|part| part == "HC" || part == "OS")
                        } else {
                            text.contains(pattern)
                        }
                    })
                });
            first.column(
                boolean_cells(&result, first.dtype == "string", false),
                bool_dtype(first),
            )
        }
        "not" | "invert" => {
            if first.dtype == "object" && first.cells.iter().any(|cell| matches!(cell, Cell::Na)) {
                return Err("object pandas.NA logical evaluation remains mature".into());
            }
            let cells: Vec<Cell> = first
                .cells
                .iter()
                .map(|cell| {
                    if op == "invert" && !matches!(cell, Cell::Bool(_) | Cell::Na) {
                        return Err("unsupported invert scalar".into());
                    }
                    Ok(match cell.truth()? {
                        Some(value) => Cell::Bool(value),
                        None => Cell::Na,
                    })
                })
                .collect::<Result<Vec<_>>>()?;
            let result = !boolean_column(&cells)?;
            first.column(
                boolean_cells(&result, bool_dtype(first) == "boolean", false),
                bool_dtype(first),
            )
        }
        "sum" | "any" => {
            if op == "sum"
                && !matches!(first.dtype.as_str(), "bool" | "boolean" | "int64" | "Int64")
            {
                return Err("unsupported noninteger source reduction".into());
            }
            if matches!(first.dtype.as_str(), "bool" | "boolean") {
                let column = boolean_column(&first.cells)?;
                if op == "any" {
                    Ok(Value::scalar(Cell::Bool(column.any())))
                } else {
                    Ok(Value::scalar(Cell::Int(i64::from(
                        column.sum().unwrap_or(0),
                    ))))
                }
            } else if matches!(first.dtype.as_str(), "int64" | "Int64") {
                let values = first
                    .cells
                    .iter()
                    .map(|cell| match cell {
                        Cell::Int(value) => Ok(Some(*value)),
                        Cell::Na | Cell::None => Ok(None),
                        _ => Err("unsupported integer reduction cell".into()),
                    })
                    .collect::<Result<Vec<_>>>()?;
                let column = Int64Chunked::from_iter_options("integer".into(), values.into_iter());
                if op == "sum" {
                    Ok(Value::scalar(Cell::Int(column.sum().unwrap_or(0))))
                } else {
                    Ok(Value::scalar(Cell::Bool(column.not_equal(0i64).any())))
                }
            } else {
                Err("unsupported non-source reduction".into())
            }
        }
        "unique" | "sorted" => {
            let mut seen = Vec::new();
            let mut unique = HashSet::new();
            for cell in &first.cells {
                let Cell::Text(text) = cell else {
                    return Err("unsupported string uniqueness operand".into());
                };
                if op == "sorted" || unique.insert(text.as_str()) {
                    seen.push(cell.clone());
                }
            }
            if op == "sorted" {
                if seen.iter().any(|cell| !matches!(cell, Cell::Text(_))) {
                    return Err("unsupported sorted operand".into());
                }
                seen.sort_by(|a, b| {
                    if let (Cell::Text(a), Cell::Text(b)) = (a, b) {
                        a.cmp(b)
                    } else {
                        unreachable!()
                    }
                });
            }
            Ok(Value {
                kind: if op == "sorted" { "list" } else { "array" }.into(),
                dtype: first.dtype.clone(),
                index: Vec::new(),
                cells: seen,
            })
        }
        "map" => {
            if args.len() != 3 || args[1].cells.len() != args[2].cells.len() {
                return Err("map key/value arity".into());
            }
            let mut mapping = HashMap::new();
            for (key, value) in args[1].cells.iter().zip(&args[2].cells) {
                let Cell::Text(key) = key else {
                    return Err("unsupported nonstring mapping key".into());
                };
                if mapping.insert(key.as_str(), value).is_some() {
                    return Err("duplicate mapping key".into());
                }
            }
            let cells: Vec<Cell> = first
                .cells
                .iter()
                .map(|cell| {
                    if cell.missing() {
                        return Cell::Float(f64::NAN);
                    }
                    if let Cell::Text(key) = cell {
                        mapping
                            .get(key.as_str())
                            .map(|value| (*value).clone())
                            .unwrap_or(Cell::Float(f64::NAN))
                    } else {
                        Cell::Float(f64::NAN)
                    }
                })
                .collect();
            let dtype = if cells
                .iter()
                .all(|cell| cell.missing() || cell.number().is_some())
            {
                "float64"
            } else {
                "object"
            };
            first.column(cells, dtype)
        }
        "max_columns" => {
            if args
                .iter()
                .any(|arg| arg.cells.len() != first.cells.len() || arg.index != first.index)
            {
                return Err("unaligned maximum columns".into());
            }
            let mut result = numeric_column(&first.cells)?;
            for value in &args[1..] {
                let next = numeric_column(&value.cells)?;
                // Polars owns/chunk-aligns every complete input/output column.
                // This small explicit kernel is pandas skipna horizontal max,
                // avoiding a new Polars feature/dependency or row allocations.
                result =
                    arity::binary_elementwise(&result, &next, |a: Option<f64>, b: Option<f64>| {
                        match (a, b) {
                            (None, value) | (value, None) => value,
                            (Some(a), Some(b)) => Some(if b >= a { b } else { a }),
                        }
                    });
            }
            let nullable = args
                .iter()
                .any(|v| matches!(v.dtype.as_str(), "Float64" | "Int64"));
            first.column(
                numeric_cells(&result, nullable),
                if nullable { "Float64" } else { "float64" },
            )
        }
        "zeros" => {
            let [Cell::Int(size)] = first.cells.as_slice() else {
                return Err("integer zeros size required".into());
            };
            let size = usize::try_from(*size).map_err(|_| "negative zeros size")?;
            Ok(Value {
                kind: "array".into(),
                dtype: "int64".into(),
                index: Vec::new(),
                cells: vec![Cell::Int(0); size],
            })
        }
        "assign" => {
            if args.len() != 3
                || first.kind != "array"
                || first.cells.len() != args[1].cells.len()
                || args[2].kind != "scalar"
            {
                return Err("unsupported assignment shape".into());
            }
            if !matches!(args[1].dtype.as_str(), "bool" | "boolean")
                || args[1]
                    .cells
                    .iter()
                    .any(|cell| !matches!(cell, Cell::Bool(_)))
            {
                return Err("unsupported numpy assignment mask dtype/missingness".into());
            }
            let values = first
                .cells
                .iter()
                .map(|cell| match cell {
                    Cell::Int(value) => Ok(*value),
                    _ => Err("unsupported criterion array dtype".into()),
                })
                .collect::<Result<Vec<_>>>()?;
            let Cell::Int(assigned) = &args[2].cells[0] else {
                return Err("unsupported assignment scalar".into());
            };
            let assigned = *assigned;
            let column = Int64Chunked::from_vec("criterion".into(), values);
            let mask = boolean_column(&args[1].cells)?;
            let result: Int64Chunked = arity::binary_elementwise(
                &column,
                &mask,
                |value: Option<i64>, mask: Option<bool>| {
                    if mask == Some(true) {
                        Some(assigned)
                    } else {
                        value
                    }
                },
            );
            first.column(
                result.into_no_null_iter().map(Cell::Int).collect(),
                &first.dtype,
            )
        }
        _ => Err(format!("unsupported operation {op}")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn floats(values: &[f64]) -> Value {
        Value {
            kind: "series".into(),
            dtype: "float64".into(),
            index: (0..values.len()).map(|i| format!("i{i:020}")).collect(),
            cells: values.iter().copied().map(Cell::Float).collect(),
        }
    }
    fn texts(values: &[&str]) -> Value {
        Value {
            kind: "series".into(),
            dtype: "object".into(),
            index: (0..values.len()).map(|i| format!("i{i:020}")).collect(),
            cells: values.iter().map(|s| Cell::Text((*s).into())).collect(),
        }
    }

    #[test]
    fn strict_and_inclusive_boundaries_are_distinct() {
        let x = floats(&[0.79, 0.8, 0.81]);
        assert_eq!(
            calculate("gt", &[x.clone(), Value::scalar(Cell::Float(0.8))])
                .unwrap()
                .cells,
            vec![Cell::Bool(false), Cell::Bool(false), Cell::Bool(true)]
        );
        assert_eq!(
            calculate("lt", &[x, Value::scalar(Cell::Float(0.8))])
                .unwrap()
                .cells,
            vec![Cell::Bool(true), Cell::Bool(false), Cell::Bool(false)]
        );
        assert_eq!(
            calculate("ge", &[floats(&[19.9, 20.0]), Value::scalar(Cell::Int(20))])
                .unwrap()
                .cells,
            vec![Cell::Bool(false), Cell::Bool(true)]
        );
    }

    #[test]
    fn lof_delimiters_and_python_final_newline_remain_literal() {
        let args = [
            texts(&["HC", "OS;LC", " HC ", "XHC", "HC\n", "HC\n\n"]),
            Value::scalar(Cell::Text(r"(?:^|[&,;|])(?:HC|OS)(?:[&,;|]|$)".into())),
        ];
        assert_eq!(
            calculate("contains", &args).unwrap().cells,
            [true, true, false, false, true, false]
                .into_iter()
                .map(Cell::Bool)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn numeric_overflow_and_missing_are_not_infinity() {
        let got = calculate("numeric", &[texts(&["0.81", "1e309", "bad", "inf", "  "])]).unwrap();
        assert_eq!(got.cells[0], Cell::Float(0.81));
        assert!(got.cells[1].missing() && got.cells[2].missing() && got.cells[4].missing());
        assert_eq!(got.cells[3], Cell::Float(f64::INFINITY));
    }

    #[test]
    fn pvs1_alignment_differs_from_position() {
        let left = Value {
            kind: "series".into(),
            dtype: "bool".into(),
            index: vec!["a".into(), "b".into()],
            cells: vec![Cell::Bool(true), Cell::Bool(true)],
        };
        let right = Value {
            kind: "series".into(),
            dtype: "bool".into(),
            index: vec!["b".into(), "a".into()],
            cells: vec![Cell::Bool(false), Cell::Bool(true)],
        };
        assert_eq!(
            calculate("and", &[left.clone(), right.clone()])
                .unwrap()
                .cells,
            vec![Cell::Bool(true), Cell::Bool(false)]
        );
        let mut positional = right;
        positional.kind = "array".into();
        positional.index.clear();
        assert_eq!(
            calculate("and", &[left, positional]).unwrap().cells,
            vec![Cell::Bool(false), Cell::Bool(true)]
        );
    }

    #[test]
    fn all_missing_maximum_supplies_no_benign_evidence() {
        let got = calculate(
            "max_columns",
            &[floats(&[f64::NAN, 0.09]), floats(&[f64::NAN, 0.1])],
        )
        .unwrap();
        assert!(got.cells[0].missing());
        assert_eq!(got.cells[1], Cell::Float(0.1));
    }

    #[test]
    fn bp7_or_preserves_missing_disjunct() {
        let a = Value {
            kind: "array".into(),
            dtype: "bool".into(),
            index: vec![],
            cells: vec![Cell::Bool(false)],
        };
        let b = Value {
            cells: vec![Cell::Bool(true)],
            ..a.clone()
        };
        assert_eq!(
            calculate("or", &[a, b]).unwrap().cells,
            vec![Cell::Bool(true)]
        );
    }

    #[test]
    fn polars_kleene_boolean_masks_keep_source_null_policy() {
        let left = Value {
            kind: "array".into(),
            dtype: "boolean".into(),
            index: vec![],
            cells: vec![Cell::Na, Cell::Na, Cell::Bool(true)],
        };
        let right = Value {
            cells: vec![Cell::Bool(false), Cell::Bool(true), Cell::Na],
            ..left.clone()
        };
        assert_eq!(
            calculate("and", &[left.clone(), right.clone()])
                .unwrap()
                .cells,
            vec![Cell::Bool(false), Cell::Na, Cell::Na]
        );
        assert_eq!(
            calculate("or", &[left, right]).unwrap().cells,
            vec![Cell::Na, Cell::Bool(true), Cell::Bool(true)]
        );
    }

    #[test]
    fn duplicate_axis_is_rejected_before_alignment() {
        let a = Value {
            kind: "series".into(),
            dtype: "bool".into(),
            index: vec!["a".into(), "a".into()],
            cells: vec![Cell::Bool(true); 2],
        };
        let b = Value {
            index: vec!["a".into(), "b".into()],
            ..a.clone()
        };
        assert!(calculate("and", &[a.clone(), a.clone()]).is_err());
        assert!(calculate("and", &[a, b]).is_err());
    }

    #[test]
    fn warning_sorting_keeps_unicode_and_sorted_multiplicity() {
        let mut input = texts(&["λ", "a", "a"]);
        input.kind = "array".into();
        input.index.clear();
        assert_eq!(
            calculate("unique", &[input.clone()]).unwrap().cells,
            vec![Cell::Text("λ".into()), Cell::Text("a".into())]
        );
        assert_eq!(
            calculate("sorted", &[input]).unwrap().cells,
            vec![
                Cell::Text("a".into()),
                Cell::Text("a".into()),
                Cell::Text("λ".into())
            ]
        );
    }
}
