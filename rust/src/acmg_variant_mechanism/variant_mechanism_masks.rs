//! Source lines468–591. Preserve column-phase evaluation, not row-wise shortcuts.
use crate::acmg_scoring::{SourceError, python_strip};
use polars_core::prelude::*;

#[derive(Clone, Debug, PartialEq)]
pub enum Cell {
    Null,
    Text(String),
    Integer(i64),
    Unsigned(u64),
    Float(f64),
    Bool(bool),
}
#[derive(Clone, Debug, PartialEq)]
pub enum Dtype {
    Object,
    String,
    Numpy,
    Nullable(&'static str),
}
#[derive(Clone, Debug)]
pub struct InputColumn {
    pub name: String,
    pub dtype: Dtype,
    pub values: Vec<Cell>,
}
#[derive(Clone, Debug)]
pub struct Input<I = Cell> {
    pub columns: Vec<InputColumn>,
    pub index: Vec<I>,
    pub index_names: Vec<Option<String>>,
}
/// Index labels are retained independently of row positions; duplicates are not joins.
pub struct Masks<I = Cell> {
    pub index: Vec<I>,
    pub index_names: Vec<Option<String>>,
    pub source_series_names: Vec<Option<&'static str>>,
    pub frame: DataFrame,
}
pub const NAMES: [&str; 27] = [
    "modern_profile",
    "is_exact_gof",
    "is_exact_dn",
    "is_exact_nonlof",
    "is_predicted_lof",
    "is_uncertain",
    "has_recessive_compatible",
    "has_dominant_compatible",
    "has_x_linked_recessive_compatible",
    "has_x_linked_dominant_compatible",
    "has_x_linked_unspecified_compatible",
    "has_rec_lof_history",
    "has_rec_gof_history",
    "has_rec_dn_history",
    "has_rec_unresolved_history",
    "has_dom_lof_history",
    "has_dom_gof_history",
    "has_dom_dn_history",
    "has_dom_unresolved_history",
    "has_x_recessive_lof_history",
    "has_x_recessive_gof_history",
    "has_x_recessive_dn_history",
    "has_x_recessive_unresolved_history",
    "has_x_dominant_lof_history",
    "has_x_dominant_gof_history",
    "has_x_dominant_dn_history",
    "has_x_dominant_unresolved_history",
];

fn missing(v: &Cell) -> bool {
    matches!(v, Cell::Null) || matches!(v,Cell::Float(x) if x.is_nan())
}
fn text(column: &InputColumn) -> Result<Vec<String>, SourceError> {
    if let Dtype::Nullable(dtype) = column.dtype {
        if column.values.iter().any(missing) {
            return Err(SourceError(
                "TypeError",
                format!("Invalid value '' for dtype {dtype}"),
            ));
        }
    }
    // Primitive non-string scalars cannot match any source tag/effect. Their
    // spelling still matters to modern_profile only through non-emptiness.
    Ok(column
        .values
        .iter()
        .map(|v| match v {
            v if missing(v) => String::new(),
            Cell::Text(s) => s.clone(),
            Cell::Integer(n) => n.to_string(),
            Cell::Unsigned(n) => n.to_string(),
            Cell::Float(n) => n.to_string(),
            Cell::Bool(v) => if *v { "True" } else { "False" }.into(),
            Cell::Null => unreachable!(),
        })
        .collect())
}
fn ascii_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
}
/// pandas2.2.3 floatify/precise_xstrtod semantics, not Rust's nearest-roundtrip parser.
/// Mirrors ordered17-digit accumulation, decimal scaling and ASCII-only whitespace.
/// Algorithm provenance: pandas2.2.3 (BSD3-Clause), PyData Development Team;
/// tokenizer.c (Lambda Foundry2012 / WarrenWeckesser textreader), precise_xstrtod;
/// pd_parser.c floatify; lib.pyx maybe_convert_numeric. No fast-math reassociation.
/// https://github.com/pandas-dev/pandas/blob/v2.2.3/pandas/_libs/src/parser/tokenizer.c
/// https://github.com/pandas-dev/pandas/blob/v2.2.3/pandas/_libs/lib.pyx
/// Original license: https://github.com/pandas-dev/pandas/blob/v2.2.3/LICENSE
/// Only the downstream >=1/==2 predicates consume this f64; dtype inference is
/// retained at the input-column boundary rather than inferred independently per row.
fn parse_numeric_text(value: &str) -> Option<(f64, Option<i128>)> {
    let mut maybe_int = true;
    // The original C parser sees the UTF-8 buffer up to its first NUL.
    let raw = value.split('\0').next().unwrap_or("");
    let inf = match raw.to_ascii_lowercase().as_str() {
        "inf" | "+inf" | "infinity" | "+infinity" => Some(f64::INFINITY),
        "-inf" | "-infinity" => Some(f64::NEG_INFINITY),
        _ => None,
    };
    if let Some(v) = inf {
        return Some((v, None));
    }
    let b = raw.as_bytes();
    let mut i = 0;
    while i < b.len() && ascii_space(b[i]) {
        i += 1;
    }
    let negative = i < b.len() && b[i] == b'-';
    if i < b.len() && matches!(b[i], b'+' | b'-') {
        i += 1;
    }
    let mut number = 0.0;
    let mut digits = 0;
    let mut exponent = 0i32;
    while i < b.len() && b[i].is_ascii_digit() {
        if digits < 17 {
            number = number * 10.0 + (b[i] - b'0') as f64;
            digits += 1;
        } else {
            exponent = exponent.saturating_add(1);
        }
        i += 1;
    }
    if i < b.len() && b[i] == b'.' {
        maybe_int = false;
        i += 1;
        while i < b.len() && b[i].is_ascii_digit() {
            if digits < 17 {
                number = number * 10.0 + (b[i] - b'0') as f64;
                digits += 1;
                exponent -= 1;
            }
            i += 1;
        }
    }
    if digits == 0 {
        return None;
    }
    if negative {
        number = -number;
    }
    if i < b.len() && matches!(b[i], b'e' | b'E') {
        maybe_int = false;
        i += 1;
        let minus = i < b.len() && b[i] == b'-';
        if i < b.len() && matches!(b[i], b'+' | b'-') {
            i += 1;
        }
        let mut n = 0i32;
        let mut count = 0;
        while i < b.len() && b[i].is_ascii_digit() && count < 17 {
            n = n.wrapping_mul(10).wrapping_add((b[i] - b'0') as i32);
            i += 1;
            count += 1;
        }
        if count == 0 {
            return None;
        }
        exponent = if minus {
            exponent.wrapping_sub(n)
        } else {
            exponent.wrapping_add(n)
        };
    }
    while i < b.len() && ascii_space(b[i]) {
        i += 1;
    }
    if i != b.len() || exponent > 308 {
        return None;
    }
    // Parsing a power literal gives the same binary constants as C's lookup table;
    // powi/powf instead may round differently before the final division.
    let power = |n: i32| format!("1e{n}").parse::<f64>().unwrap();
    if exponent > 0 {
        number *= power(exponent);
    } else if exponent < -616 {
        number = 0.0;
    } else if exponent < -308 {
        number /= power(-308 - exponent);
        number /= power(308);
    } else {
        number /= power(-exponent);
    }
    if !number.is_finite() {
        return None;
    }
    let integer = if maybe_int {
        // lib.maybe_convert_numeric calls int(original), not the NUL prefix.
        if value.contains('\0') {
            return None;
        }
        let original = value.trim_matches(|c: char| c.is_ascii() && ascii_space(c as u8));
        match original.parse::<i128>() {
            Ok(n) => Some(n),
            Err(_) => None, // Valid out-of-i128 integer strings remain the parsed float.
        }
    } else {
        None
    };
    Some((number, integer))
}

/// One-element object-column convenience; whole-frame evaluation uses numeric().
pub fn numeric_text(value: &str) -> f64 {
    numeric(
        Some(&InputColumn {
            name: String::new(),
            dtype: Dtype::Object,
            values: vec![Cell::Text(value.into())],
        }),
        1,
    )
    .expect("one-element object column cannot raise nullable fillna")[0]
}
fn numeric(column: Option<&InputColumn>, len: usize) -> Result<Vec<f64>, SourceError> {
    let Some(column) = column else {
        return Ok(vec![0.0; len]);
    };
    if column.dtype == Dtype::Nullable("boolean") && column.values.iter().any(missing) {
        return Err(SourceError(
            "TypeError",
            "Invalid value '0' for dtype boolean".into(),
        ));
    }
    let nullable = matches!(column.dtype, Dtype::String | Dtype::Nullable(_));
    let mut floats = Vec::with_capacity(len);
    let mut ints = Vec::with_capacity(len);
    let mut seen_float = false;
    let mut seen_signed = false;
    let mut seen_unsigned = false;
    for value in &column.values {
        let parsed = match value {
            Cell::Text(s) => parse_numeric_text(s),
            Cell::Null => None,
            Cell::Integer(n) => Some((*n as f64, Some(*n as i128))),
            Cell::Unsigned(n) => Some((*n as f64, Some(*n as i128))),
            Cell::Bool(v) => Some((if *v { 1.0 } else { 0.0 }, Some(if *v { 1 } else { 0 }))),
            Cell::Float(n) => {
                if n.is_nan() {
                    None
                } else {
                    Some((*n, None))
                }
            }
        };
        match parsed {
            None => {
                // pandas empty-string branch always calls seen.saw_null(),
                // unlike masked StringDtype NA or failed numeric coercion.
                if !nullable || matches!(value, Cell::Text(s) if s.is_empty()) {
                    seen_float = true;
                }
                floats.push(0.0);
                ints.push(0.0);
            }
            Some((float, integer)) => {
                floats.push(float);
                if let Some(n) = integer {
                    seen_signed |= n < 0;
                    seen_unsigned |= n > i64::MAX as i128;
                    seen_float |= n < i64::MIN as i128 || n > u64::MAX as i128;
                    ints.push(n as f64);
                } else {
                    seen_float = true;
                    ints.push(0.0);
                }
            }
        }
    }
    // pandas chooses the column array after processing every value. In a float
    // column it keeps floatify's values even for integer-looking strings.
    Ok(if seen_float || (seen_signed && seen_unsigned) {
        floats
    } else {
        ints
    })
}
/// Semicolon anchors plus Python's `$` just before ONE final newline.
fn tag(s: &str, base: &str, collapsed: bool) -> bool {
    let chunks: Vec<&str> = s.split(';').collect();
    chunks.iter().enumerate().any(|(i, t)| {
        let accepts = |t: &str| {
            t == base
                || (collapsed
                    && t.strip_prefix(base)
                        .is_some_and(|r| r.starts_with('_') && r.len() > 1))
        };
        accepts(t) || (i + 1 == chunks.len() && t.strip_suffix('\n').is_some_and(accepts))
    })
}
fn effect_match(s: &str, token: &str) -> bool {
    let right = |r: &str| r.is_empty() || r.starts_with('+') || r == "\n";
    if s.strip_prefix("exact_known_")
        .and_then(|r| r.strip_prefix(token))
        .is_some_and(right)
    {
        return true;
    }
    s.match_indices('+')
        .any(|(i, _)| s[i + 1..].strip_prefix(token).is_some_and(right))
}

pub fn variant_mechanism_masks<I: Clone>(input: &Input<I>) -> Result<Masks<I>, SourceError> {
    let lookup = |name: &str| input.columns.iter().find(|c| c.name == name);
    let absent: Vec<&str> = ["var_plausible_patho_mechs", "variant_effect"]
        .into_iter()
        .filter(|n| lookup(n).is_none())
        .collect();
    if !absent.is_empty() {
        return Err(SourceError(
            "KeyError",
            format!(
                "'variant-level mechanism annotations are required; missing columns: {}'",
                absent.join(", ")
            ),
        ));
    }
    let len = input.index.len();
    if input.columns.iter().any(|c| c.values.len() != len) {
        return Err(SourceError(
            "ValueError",
            "column/index length mismatch".into(),
        ));
    }
    let duplicates = |name: &str| input.columns.iter().filter(|c| c.name == name).count() > 1;
    let dataframe_str = || {
        SourceError(
            "AttributeError",
            "'DataFrame' object has no attribute 'str'".into(),
        )
    };
    if duplicates("var_plausible_patho_mechs") {
        // pandas selects a DataFrame here: fillna/astype run before .str fails.
        for column in input
            .columns
            .iter()
            .filter(|c| c.name == "var_plausible_patho_mechs")
        {
            text(column)?;
        }
        return Err(dataframe_str());
    }
    let profile = text(lookup("var_plausible_patho_mechs").unwrap())?;
    let modern: Vec<bool> = profile
        .iter()
        .map(|s| !python_strip(s).is_empty())
        .collect();
    let inheritance = [
        "recessive",
        "dominant",
        "x_linked_recessive",
        "x_linked_dominant",
        "x_linked_unspecified",
    ];
    let collapsed: Vec<Vec<bool>> = inheritance
        .iter()
        .map(|base| profile.iter().map(|s| tag(s, base, true)).collect())
        .collect();
    let exact: Vec<Vec<bool>> = inheritance[..4]
        .iter()
        .flat_map(|base| {
            ["_LOF", "_GOF", "_DN", ""].map(|suffix| {
                profile
                    .iter()
                    .map(|s| tag(s, &format!("{base}{suffix}"), false))
                    .collect()
            })
        })
        .collect();
    let effect = if duplicates("variant_effect") {
        // Effect fillna errors still precede whole-column score conversion.
        for column in input.columns.iter().filter(|c| c.name == "variant_effect") {
            text(column)?;
        }
        vec![String::new(); len]
    } else {
        text(lookup("variant_effect").unwrap())?
    };
    if duplicates("variant_lof_score") {
        return Err(SourceError(
            "TypeError",
            "arg must be a list, tuple, 1-d array, or Series".into(),
        ));
    }
    let lof = numeric(lookup("variant_lof_score"), len)?;
    if duplicates("variant_gof_score") {
        return Err(SourceError(
            "TypeError",
            "arg must be a list, tuple, 1-d array, or Series".into(),
        ));
    }
    let gof = numeric(lookup("variant_gof_score"), len)?;
    if duplicates("variant_dn_score") {
        return Err(SourceError(
            "TypeError",
            "arg must be a list, tuple, 1-d array, or Series".into(),
        ));
    }
    let dn = numeric(lookup("variant_dn_score"), len)?;
    if duplicates("variant_effect") {
        return Err(dataframe_str());
    }
    let exact_gof: Vec<bool> = (0..len)
        .map(|i| gof[i] == 2.0 || effect_match(&effect[i], "GOF"))
        .collect();
    let exact_dn: Vec<bool> = (0..len)
        .map(|i| dn[i] == 2.0 || effect_match(&effect[i], "DOMINANT_NEGATIVE"))
        .collect();
    let uncertain: Vec<bool> = effect.iter().map(|s| s == "uncertain").collect();
    let predicted: Vec<bool> = (0..len)
        .map(|i| !uncertain[i] && (lof[i] >= 1.0 || effect[i] == "predicted_LOF_high_confidence"))
        .collect();
    let nonlof: Vec<bool> = (0..len).map(|i| exact_gof[i] || exact_dn[i]).collect();
    let values = [
        vec![modern, exact_gof, exact_dn, nonlof, predicted, uncertain],
        collapsed,
        exact,
    ]
    .concat();
    let columns = values
        .into_iter()
        .zip(NAMES)
        .map(|(v, n)| Column::new(n.into(), v))
        .collect();
    let frame = DataFrame::new(columns).map_err(|e| SourceError("PolarsError", e.to_string()))?;
    let source_series_names = (0..27)
        .map(|i| match i {
            1..=4 => None,
            5 => Some("variant_effect"),
            _ => Some("var_plausible_patho_mechs"),
        })
        .collect();
    Ok(Masks {
        index: input.index.clone(),
        index_names: input.index_names.clone(),
        source_series_names,
        frame,
    })
}

/// Polars primitive-column adapter. `column_dtypes` explicitly distinguishes pandas
/// nullable extension columns from NumPy/object columns (Polars alone loses this).
/// Index labels remain separate: no sorting, deduplication, or positional reset.
pub fn from_polars<I: Clone>(
    df: &DataFrame,
    index: Vec<I>,
    index_names: Vec<Option<String>>,
    column_dtypes: [Dtype; 5],
) -> Result<Masks<I>, SourceError> {
    let mut columns = Vec::new();
    // Missing-column failures must precede unsupported dtype errors.
    let missing: Vec<&str> = ["var_plausible_patho_mechs", "variant_effect"]
        .into_iter()
        .filter(|n| df.column(n).is_err())
        .collect();
    if !missing.is_empty() {
        return Err(SourceError(
            "KeyError",
            format!(
                "'variant-level mechanism annotations are required; missing columns: {}'",
                missing.join(", ")
            ),
        ));
    }
    for (i, name) in [
        "var_plausible_patho_mechs",
        "variant_effect",
        "variant_lof_score",
        "variant_gof_score",
        "variant_dn_score",
    ]
    .iter()
    .enumerate()
    {
        let Ok(c) = df.column(name) else { continue };
        // Check the declared dtype even for empty/all-null columns. Otherwise
        // AnyValue iteration would silently accept an unsupported typed domain.
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
                match c
                    .get(row)
                    .map_err(|e| SourceError("PolarsError", e.to_string()))?
                {
                    AnyValue::Null => Ok(Cell::Null),
                    AnyValue::String(s) => Ok(Cell::Text(s.into())),
                    AnyValue::StringOwned(s) => Ok(Cell::Text(s.to_string())),
                    AnyValue::Boolean(b) => Ok(Cell::Bool(b)),
                    AnyValue::Float64(n) => Ok(Cell::Float(n)),
                    AnyValue::Float32(n) => Ok(Cell::Float(n as f64)),
                    AnyValue::Int64(n) => Ok(Cell::Integer(n)),
                    AnyValue::Int32(n) => Ok(Cell::Integer(n as i64)),
                    AnyValue::Int16(n) => Ok(Cell::Integer(n as i64)),
                    AnyValue::Int8(n) => Ok(Cell::Integer(n as i64)),
                    AnyValue::UInt64(n) => Ok(Cell::Unsigned(n)),
                    AnyValue::UInt32(n) => Ok(Cell::Unsigned(n as u64)),
                    AnyValue::UInt16(n) => Ok(Cell::Unsigned(n as u64)),
                    AnyValue::UInt8(n) => Ok(Cell::Unsigned(n as u64)),
                    _ => Err(SourceError(
                        "UnfinishedParity",
                        format!(
                            "{name}: dtype {:?} requires original Python; no production dispatch",
                            c.dtype()
                        ),
                    )),
                }
            })
            .collect::<Result<Vec<_>, SourceError>>()?;
        let column = InputColumn {
            name: name.to_string(),
            dtype: column_dtypes[i].clone(),
            values,
        };
        // Source fillna/astype phases precede extracting later score columns.
        // A later unsupported Polars dtype must not hide an earlier nullable error.
        if i < 2 {
            text(&column)?;
        } else {
            numeric(Some(&column), index.len())?;
        }
        columns.push(column);
    }
    variant_mechanism_masks(&Input {
        columns,
        index,
        index_names,
    })
}
