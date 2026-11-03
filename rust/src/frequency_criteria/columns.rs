//! Typed Polars kernels derived from the existing native column interface.
//! Reuse is source engineering, not inherited frequency qualification.
use super::values::{Cell, Value, Result};
use polars_core::chunked_array::ops::arity;
use polars_core::prelude::*;
use std::collections::{HashMap, HashSet};

pub(crate) fn numeric_column(cells: &[Cell]) -> Result<Float64Chunked> {
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

pub(crate) fn boolean_column(cells: &[Cell]) -> Result<BooleanChunked> {
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
pub(crate) fn boolean_cells(column: &BooleanChunked, nullable: bool, missing_default: bool) -> Vec<Cell> {
    column
        .into_iter()
        .map(|value| match value {
            Some(value) => Cell::Bool(value),
            None if nullable => Cell::Na,
            None => Cell::Bool(missing_default),
        })
        .collect()
}

pub(crate) fn numeric_cells(column: &Float64Chunked, nullable: bool) -> Vec<Cell> {
    column
        .into_iter()
        .map(|value| match value {
            Some(value) => Cell::Float(value),
            None if nullable => Cell::Na,
            None => Cell::Float(f64::NAN),
        })
        .collect()
}

pub(crate) fn bool_dtype(value: &Value) -> &'static str {
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
        let digits = value.strip_prefix('+').or_else(||value.strip_prefix('-')).unwrap_or(value);
        if !digits.is_empty() && digits.bytes().all(|c|c.is_ascii_digit()) {
            let integer=value.parse::<i64>().map_err(|_|"integer string outside admitted width")?;
            if !(-9_007_199_254_740_992..=9_007_199_254_740_992).contains(&integer) {
                return Err("integer string outside exact native-float domain".into());
            }
            return Ok(Cell::Int(integer));
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
        // pandas uses its own decimal parser. Even a short lexical decimal
        // is not admitted merely because Rust accepts it; retain that exact
        // mature coercion until its byte-level domain is separately qualified.
        return match value.parse::<f64>() {
            Ok(number) if number.is_finite() => Err("decimal-string coercion retains the mature parser".into()),
            _ => Ok(Cell::Float(f64::NAN)),
        };
    }
    match cell {
        Cell::Int(_) | Cell::Bool(_) | Cell::Float(_) => Ok(cell.clone()),
        _ => Err("unsupported numeric scalar".into()),
    }
}

pub(crate) fn aligned(left: &Value, right: &Value) -> Result<(Vec<String>, Vec<Cell>, Vec<Cell>)> {
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
        let index = if right.kind == "series" { &right.index } else { &left.index };
        return Ok((index.clone(), left.cells.clone(), right.cells.clone()));
    }
    if left.index != right.index {
        return Err("nonidentical pandas axes require the mature alignment operation".into());
    }
    Ok((left.index.clone(), left.cells.clone(), right.cells.clone()))
}

fn binary(op: &str, left: &Value, right: &Value) -> Result<Value> {
    let (index, a, b) = aligned(left, right)?;
    let prototype = if left.kind == "scalar" || right.kind == "series" && left.kind != "series" { right } else { left };
    let nullable = bool_dtype(left) == "boolean" || bool_dtype(right) == "boolean";
    let result = if op == "and" || op == "or" {
        let a = boolean_column(&a)?;
        let b = boolean_column(&b)?;
        if op == "and" {
            &a & &b
        } else {
            &a | &b
        }
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

/// Native column primitives. The frequency source interpreter retains statement order.
/// An unsupported operation is reported before its mature counterpart runs once.
pub fn calculate(op: &str, args: &[Value]) -> Result<Value> {
    if let Some(value) = super::binomial::calculate(op, args) { return value; }
    if let Some(value) = super::rules::calculate(op, args) { return value; }
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
            if cells.iter().all(|c|matches!(c,Cell::Int(_)) || nullable && matches!(c,Cell::Na)) {
                return first.column(cells,if nullable {"Int64"}else{"int64"});
            }
            if !cells.is_empty() && cells.iter().all(|c|matches!(c,Cell::Bool(_))) {
                return first.column(cells,if nullable {"boolean"}else{"bool"});
            }
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
            if size > 16_777_216 { return Err("bounded array size exceeded".into()); }
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

