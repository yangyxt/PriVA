//! Exact scalar boundary for `acmg_consequence.py::coerce_pext_value`, lines28–49.
//! File a9efd0964f25e84dc01ad276d4f0905a10492a5d1f94c2b277573ac94244d354,
//! function67d26c68d7149f47719c3b8db2289895c0c4f77fc5227cc4ef288fe972a55dbe.
//! The source docstring says finite, but its code does NOT filter nonfinite values.
use crate::acmg_scoring::{SourceError, python_strip};
use crate::acmg_variant_mechanism::variant_mechanism_masks::Cell;
use polars_core::prelude::*;

const PYTHON_NAN: f64 = f64::from_bits(0x7ff8_0000_0000_0000);

/// Primitive scalar identities are explicit: NumPy bool is not Python bool.
/// BigInteger is a decimal representation of a real integer, not string input.
/// Arbitrary objects/arrays, complex, f16/longdouble and surrogate strings require
/// the unchanged Python backend; they must not be coerced into these variants.
#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    Primitive(Cell),
    NumpyBool(bool),
    Float32(f32),
    BigInteger(String),
}

#[derive(Clone, Copy, Debug)]
pub enum BooleanOrigin {
    Python,
    Numpy,
}

// Derived from actual CPython3.12.3 unicodedata15.0.0 decimal(c)==0, frozen
// before implementation in the actual-Python oracle (68 blocks /680 digits).
// Deliberately do not accept decimal digits introduced by a later Unicode version.
const DECIMAL_ZEROES: [u32; 68] = [
    0x30, 0x660, 0x6f0, 0x7c0, 0x966, 0x9e6, 0xa66, 0xae6, 0xb66, 0xbe6, 0xc66, 0xce6, 0xd66,
    0xde6, 0xe50, 0xed0, 0xf20, 0x1040, 0x1090, 0x17e0, 0x1810, 0x1946, 0x19d0, 0x1a80, 0x1a90,
    0x1b50, 0x1bb0, 0x1c40, 0x1c50, 0xa620, 0xa8d0, 0xa900, 0xa9d0, 0xa9f0, 0xaa50, 0xabf0, 0xff10,
    0x104a0, 0x10d30, 0x11066, 0x110f0, 0x11136, 0x111d0, 0x112f0, 0x11450, 0x114d0, 0x11650,
    0x116c0, 0x11730, 0x118e0, 0x11950, 0x11c50, 0x11d50, 0x11da0, 0x11f50, 0x16a60, 0x16ac0,
    0x16b50, 0x1d7ce, 0x1d7d8, 0x1d7e2, 0x1d7ec, 0x1d7f6, 0x1e140, 0x1e2f0, 0x1e4f0, 0x1e950,
    0x1fbf0,
];

fn ascii_float(value: &str) -> Option<String> {
    let mut ascii = String::with_capacity(value.len());
    for c in value.chars() {
        if c.is_ascii() {
            ascii.push(c);
        } else if c.is_whitespace() {
            ascii.push(' ');
        } else {
            let code = c as u32;
            let zero = DECIMAL_ZEROES
                .iter()
                .find(|zero| code >= **zero && code - **zero < 10)?;
            ascii.push(char::from(b'0' + (code - zero) as u8));
        }
    }
    Some(ascii)
}

fn float_text(value: &str) -> Option<f64> {
    let ascii = ascii_float(value)?;
    let lower = ascii.to_ascii_lowercase();
    let (negative, body) = match lower.as_bytes().first() {
        Some(b'-') => (true, &lower[1..]),
        Some(b'+') => (false, &lower[1..]),
        _ => (false, lower.as_str()),
    };
    if matches!(body, "inf" | "infinity") {
        return Some(if negative {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        });
    }
    if body == "nan" {
        return Some(f64::from_bits(
            PYTHON_NAN.to_bits() | (u64::from(negative) << 63),
        ));
    }
    // CPython accepts underscores only between decimal digits, after Unicode
    // decimal transformation. Removing all underscores without this check is wrong.
    let bytes = ascii.as_bytes();
    let mut normalized = String::with_capacity(bytes.len());
    for (i, &c) in bytes.iter().enumerate() {
        if c == b'_' {
            if i == 0
                || i + 1 == bytes.len()
                || !bytes[i - 1].is_ascii_digit()
                || !bytes[i + 1].is_ascii_digit()
            {
                return None;
            }
        } else {
            normalized.push(char::from(c));
        }
    }
    let bytes = normalized.as_bytes();
    let mut i = usize::from(matches!(bytes.first(), Some(b'+' | b'-')));
    let mut digits = 0;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
        digits += 1;
    }
    if bytes.get(i) == Some(&b'.') {
        i += 1;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
            digits += 1;
        }
    }
    if digits == 0 {
        return None;
    }
    if matches!(bytes.get(i), Some(b'e' | b'E')) {
        i += 1;
        if matches!(bytes.get(i), Some(b'+' | b'-')) {
            i += 1;
        }
        let first_digit = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        if i == first_digit {
            return None;
        }
    }
    if i != bytes.len() {
        return None;
    }
    normalized.parse().ok()
}

fn integer_float(value: &str) -> Result<f64, SourceError> {
    let digits = value.strip_prefix(['+', '-']).unwrap_or(value);
    if digits.is_empty() || !digits.bytes().all(|c| c.is_ascii_digit()) {
        return Err(SourceError(
            "UnfinishedParity",
            "BigInteger needs an actual ASCII decimal integer; retain Python backend".into(),
        ));
    }
    if digits.bytes().all(|c| c == b'0') {
        return Ok(0.0); // Integer -0 is 0; it is not string '-0'.
    }
    let value: f64 = value.parse().map_err(|_| {
        SourceError(
            "UnfinishedParity",
            "integer parser rejected a validated decimal; retain Python backend".into(),
        )
    })?;
    if !value.is_finite() {
        return Err(SourceError(
            "OverflowError",
            "int too large to convert to float".into(),
        ));
    }
    Ok(value)
}

pub fn coerce_pext_value(value: &Value) -> Result<f64, SourceError> {
    Ok(match value {
        // Preserve pd.isna FIRST: even signed/payload numeric NaNs become np.nan.
        Value::Primitive(Cell::Null) => PYTHON_NAN,
        Value::Primitive(Cell::Float(v)) if v.is_nan() => PYTHON_NAN,
        Value::Float32(v) if v.is_nan() => PYTHON_NAN,
        Value::Primitive(Cell::Float(v)) => *v,
        Value::Float32(v) => f64::from(*v),
        Value::Primitive(Cell::Integer(v)) => *v as f64,
        Value::Primitive(Cell::Unsigned(v)) => *v as f64,
        Value::Primitive(Cell::Bool(v)) => f64::from(u8::from(*v)),
        Value::BigInteger(v) => integer_float(v)?,
        // np.bool_ fails the source isinstance test and str(True/False) is not numeric.
        Value::NumpyBool(_) => PYTHON_NAN,
        Value::Primitive(Cell::Text(value)) => {
            let mut best = None;
            for item in value.split(['&', ',']) {
                let item = python_strip(item);
                if matches!(item, "" | "." | "NA" | "NaN" | "nan") {
                    continue;
                }
                if let Some(candidate) = float_text(item) {
                    // Python max retains the FIRST value when > is false,
                    // including ties, signed zeros, and unordered NaN comparisons.
                    if best.is_none_or(|previous| candidate > previous) {
                        best = Some(candidate);
                    }
                }
            }
            best.unwrap_or(PYTHON_NAN)
        }
    })
}

pub fn coerce_values(values: &[Value], name: &str) -> Result<Series, SourceError> {
    let values = values
        .iter()
        .map(coerce_pext_value)
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Series::new(name.into(), values)) // Missing results are IEEE NaN, never Polars null.
}

/// Explicit primitive Polars adapter. Bool boxing must be selected by the caller;
/// no inference about pandas scalar origin, arbitrary objects, categorical/time/
/// decimal/list/struct dtypes or unsupported empty columns is claimed.
pub fn from_polars(column: &Column, bool_origin: BooleanOrigin) -> Result<Series, SourceError> {
    if !matches!(
        column.dtype(),
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
                "dtype {:?} needs unchanged Python; no production dispatch",
                column.dtype()
            ),
        ));
    }
    let mut values = Vec::with_capacity(column.len());
    for row in 0..column.len() {
        let cell = match column
            .get(row)
            .map_err(|e| SourceError("PolarsError", e.to_string()))?
        {
            AnyValue::Null => Value::Primitive(Cell::Null),
            AnyValue::String(v) => Value::Primitive(Cell::Text(v.into())),
            AnyValue::StringOwned(v) => Value::Primitive(Cell::Text(v.to_string())),
            AnyValue::Boolean(v) => match bool_origin {
                BooleanOrigin::Python => Value::Primitive(Cell::Bool(v)),
                BooleanOrigin::Numpy => Value::NumpyBool(v),
            },
            AnyValue::Float64(v) => Value::Primitive(Cell::Float(v)),
            AnyValue::Float32(v) => Value::Float32(v),
            AnyValue::Int64(v) => Value::Primitive(Cell::Integer(v)),
            AnyValue::Int32(v) => Value::Primitive(Cell::Integer(i64::from(v))),
            AnyValue::Int16(v) => Value::Primitive(Cell::Integer(i64::from(v))),
            AnyValue::Int8(v) => Value::Primitive(Cell::Integer(i64::from(v))),
            AnyValue::UInt64(v) => Value::Primitive(Cell::Unsigned(v)),
            AnyValue::UInt32(v) => Value::Primitive(Cell::Unsigned(u64::from(v))),
            AnyValue::UInt16(v) => Value::Primitive(Cell::Unsigned(u64::from(v))),
            AnyValue::UInt8(v) => Value::Primitive(Cell::Unsigned(u64::from(v))),
            _ => unreachable!("dtype checked"),
        };
        values.push(cell);
    }
    coerce_values(&values, column.name().as_str())
}
