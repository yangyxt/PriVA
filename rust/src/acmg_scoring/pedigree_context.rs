//! Original sort_and_rank_variants R2, 09223ef lines501–524 only.
//! Primitive context/mask plan. No control/scoring numerics or dispatch changes.
use super::terminal_order::{Input, InputColumn};
use super::{SourceError, python_strip};
use crate::acmg_variant_mechanism::variant_mechanism_masks::Cell;
use polars_core::prelude::*;

#[derive(Clone, Debug)]
pub enum Family {
    /// Python None disables the PED branch; pd.NA is Value(Cell::Null) instead.
    None,
    Value(Cell),
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Warning {
    NoPedigree,
    UnknownSex,
}
impl Warning {
    pub fn name(self) -> &'static str {
        match self {
            Self::NoPedigree => "no_pedigree",
            Self::UnknownSex => "unknown_sex",
        }
    }
}
pub struct Context<I> {
    pub pedigree_position: Option<usize>,
    pub proband: Option<String>,
    pub is_male: bool,
    pub is_female: bool,
    pub mask_dtype: &'static str,
    /// Raw Boolean data and separate null flags; masked source bits are retained.
    pub mask: DataFrame,
    pub index: Vec<I>,
    pub axis_metadata: String,
    pub inherit_frame_metadata: bool,
    pub warning: Option<Warning>,
}
#[derive(Debug)]
pub struct Failure {
    pub error: SourceError,
    pub source_line: Option<u32>,
    pub pedigree_position: Option<usize>,
    pub warning: Option<Warning>,
}
impl std::fmt::Display for Failure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.error.fmt(f)
    }
}
impl std::error::Error for Failure {}

fn failure(
    kind: &'static str,
    message: impl Into<String>,
    line: Option<u32>,
    warning: Option<Warning>,
    row: Option<usize>,
) -> Failure {
    Failure {
        error: SourceError(kind, message.into()),
        source_line: line,
        pedigree_position: row,
        warning,
    }
}
fn unsupported(
    message: impl Into<String>,
    warning: Option<Warning>,
    row: Option<usize>,
) -> Failure {
    failure(
        "UnsupportedInput",
        format!("{}; retain Python backend", message.into()),
        None,
        warning,
        row,
    )
}
fn valid(column: &InputColumn, n: usize) -> bool {
    column.values.len() == n
        && column
            .values
            .iter()
            .all(|value| match (column.dtype.as_str(), value) {
                ("object", _) => true,
                ("string", Cell::Text(_) | Cell::Null) => true,
                ("float64", Cell::Float(_)) => true,
                ("Float64", Cell::Float(x)) => !x.is_nan(),
                ("Float64" | "Int64" | "UInt64" | "boolean", Cell::Null) => true,
                ("int64" | "Int64", Cell::Integer(_)) => true,
                ("uint64" | "UInt64", Cell::Unsigned(_)) => true,
                ("uint64" | "UInt64", Cell::Integer(x)) => *x >= 0,
                ("bool" | "boolean", Cell::Bool(_)) => true,
                _ => false,
            })
        && matches!(
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
        )
}
fn required<'a, I>(input: &'a Input<I>, name: &str, line: u32) -> Result<&'a InputColumn, Failure> {
    let mut columns = input.columns.iter().filter(|column| column.name == name);
    let first = columns
        .next()
        .ok_or_else(|| failure("KeyError", format!("'{name}'"), Some(line), None, None))?;
    if columns.next().is_some() {
        return Err(unsupported(
            format!("duplicate pedigree {name} label"),
            None,
            None,
        ));
    }
    if !valid(first, input.index.len()) {
        return Err(unsupported(
            format!("pedigree {name} dtype/cells"),
            None,
            None,
        ));
    }
    Ok(first)
}
fn integer(cell: &Cell) -> Option<i128> {
    match cell {
        Cell::Integer(x) => Some(*x as i128),
        Cell::Unsigned(x) => Some(*x as i128),
        Cell::Bool(x) => Some(*x as i128),
        _ => None,
    }
}
fn float_integer_equal(value: f64, other: i128) -> bool {
    // Python compares int/float numerically without rounding the int to f64.
    value.is_finite()
        && value.fract() == 0.
        && value >= i64::MIN as f64
        && value < 18446744073709551616.0
        && value as i128 == other
}
fn equal(left: &Cell, right: &Cell) -> bool {
    match (left, right) {
        (Cell::Text(a), Cell::Text(b)) => a == b,
        (Cell::Float(a), Cell::Float(b)) => a == b,
        (Cell::Float(value), other) | (other, Cell::Float(value)) => {
            integer(other).is_some_and(|other| float_integer_equal(*value, other))
        }
        _ => integer(left)
            .zip(integer(right))
            .is_some_and(|(a, b)| a == b),
    }
}
fn sex(cell: &Cell) -> (bool, bool) {
    match cell {
        Cell::Text(value) => match python_strip(value).to_ascii_lowercase().as_str() {
            "1" | "1.0" | "m" | "male" => (true, false),
            "2" | "2.0" | "f" | "female" => (false, true),
            _ => (false, false),
        },
        Cell::Integer(value) => (*value == 1, *value == 2),
        Cell::Unsigned(value) => (*value == 1, *value == 2),
        Cell::Float(value) => (*value == 1., *value == 2.),
        // Python str(True), str(None) and str(pd.NA) are unknown sex spellings.
        _ => (false, false),
    }
}
fn mask(values: Vec<bool>, missing: Vec<bool>) -> Result<DataFrame, Failure> {
    DataFrame::new(vec![
        Column::new("value".into(), values),
        Column::new("is_null".into(), missing),
    ])
    .map_err(|error| failure("PolarsError", error.to_string(), None, None, None))
}

pub fn pedigree_context<I: Clone>(
    frame: &Input<I>,
    pedigree: Option<&Input<I>>,
    family: &Family,
) -> Result<Context<I>, Failure> {
    let n = frame.index.len();
    let (Some(pedigree), Family::Value(family)) = (pedigree, family) else {
        return Ok(Context {
            pedigree_position: None,
            proband: None,
            is_male: false,
            is_female: false,
            mask_dtype: if n == 0 { "object" } else { "bool" },
            mask: mask(vec![false; n], vec![false; n])?,
            index: frame.index.clone(),
            axis_metadata: frame.axis_metadata.clone(),
            inherit_frame_metadata: false,
            warning: Some(Warning::NoPedigree),
        });
    };
    let families = required(pedigree, "#FamilyID", 503)?;
    if !["object", "string"].contains(&families.dtype.as_str()) {
        return Err(unsupported(
            "native numeric family-column NumPy comparison",
            None,
            None,
        ));
    }
    let phenotypes = required(pedigree, "Phenotype", 504)?;
    let matches = families
        .values
        .iter()
        .zip(&phenotypes.values)
        .map(|(value, phenotype)| {
            equal(value, family)
                && (matches!(phenotype,Cell::Text(value) if value=="2")
                    || equal(phenotype, &Cell::Integer(2)))
        })
        .collect::<Vec<_>>();
    let selected = DataFrame::new(vec![
        Column::new(
            "position".into(),
            (0..pedigree.index.len() as u64).collect::<Vec<_>>(),
        ),
        Column::new("affected".into(), matches),
    ])
    .map_err(|error| failure("PolarsError", error.to_string(), None, None, None))?;
    let selected = selected
        .filter(
            selected
                .column("affected")
                .and_then(|column| column.bool())
                .map_err(|error| failure("PolarsError", error.to_string(), None, None, None))?,
        )
        .map_err(|error| failure("PolarsError", error.to_string(), None, None, None))?;
    let row = selected
        .column("position")
        .and_then(|column| column.u64())
        .map_err(|error| failure("PolarsError", error.to_string(), None, None, None))?
        .into_no_null_iter()
        .next()
        .ok_or_else(|| {
            failure(
                "IndexError",
                "single positional indexer is out-of-bounds",
                Some(502),
                None,
                None,
            )
        })? as usize;
    let identifiers = required(pedigree, "IndividualID", 506)?;
    let Cell::Text(proband) = &identifiers.values[row] else {
        return Err(unsupported("nontext selected proband", None, Some(row)));
    };
    let mut sexes = pedigree
        .columns
        .iter()
        .filter(|column| column.name == "Sex");
    let first = sexes.next();
    if sexes.next().is_some() {
        return Err(unsupported("duplicate pedigree Sex label", None, Some(row)));
    }
    if first.is_some_and(|column| !valid(column, pedigree.index.len())) {
        return Err(unsupported("pedigree Sex dtype/cells", None, Some(row)));
    }
    let default = Cell::Text(String::new());
    let (is_male, is_female) = sex(first.map(|column| &column.values[row]).unwrap_or(&default));
    let warning = if is_male || is_female {
        None
    } else {
        Some(Warning::UnknownSex)
    };
    let mut samples = frame
        .columns
        .iter()
        .filter(|column| column.name == *proband);
    let sample = samples.next().ok_or_else(|| {
        if proband
            .chars()
            .all(|value| value.is_ascii_alphanumeric() || matches!(value, '_' | '-' | '.'))
        {
            failure(
                "KeyError",
                format!("'{proband}'"),
                Some(518),
                warning,
                Some(row),
            )
        } else {
            unsupported(
                "missing sample label requires Python repr",
                warning,
                Some(row),
            )
        }
    })?;
    if samples.next().is_some() {
        return Err(failure(
            "AttributeError",
            "'DataFrame' object has no attribute 'str'",
            Some(518),
            warning,
            Some(row),
        ));
    }
    if !valid(sample, n) {
        return Err(unsupported("sample dtype/cells", warning, Some(row)));
    }
    if ["Int64", "UInt64", "Float64", "boolean"].contains(&sample.dtype.as_str()) {
        // Actual failing StringMethods construction may warm the caller's cold
        // pandas dtype caches. Primitive CSV cannot attest those object effects.
        return Err(unsupported(
            "nullable numeric sample accessor/cache effects",
            warning,
            Some(row),
        ));
    }
    let mut text = false;
    let mut boolean = false;
    let mut numeric = false;
    let mut nonmissing = false;
    for value in &sample.values {
        match value {
            Cell::Null => {}
            Cell::Float(value) if value.is_nan() => {}
            Cell::Text(_) => {
                text = true;
                nonmissing = true;
            }
            Cell::Bool(_) => {
                boolean = true;
                nonmissing = true;
            }
            _ => {
                numeric = true;
                nonmissing = true;
            }
        }
    }
    if sample.dtype != "string"
        && (sample.dtype != "object" || !(text || !nonmissing || (boolean && numeric)))
    {
        return Err(failure(
            "AttributeError",
            "Can only use .str accessor with string values!",
            Some(518),
            warning,
            Some(row),
        ));
    }
    let nullable = sample.dtype == "string";
    let values = sample
        .values
        .iter()
        .map(|value| match value {
            Cell::Text(value) => value.bytes().filter(|&b| b == b'1').count() == 1,
            Cell::Null if nullable => true,
            _ => false,
        })
        .collect();
    let missing = sample
        .values
        .iter()
        .map(|value| nullable && matches!(value, Cell::Null))
        .collect();
    Ok(Context {
        pedigree_position: Some(row),
        proband: Some(proband.clone()),
        is_male,
        is_female,
        mask_dtype: if nullable { "boolean" } else { "bool" },
        mask: mask(values, missing)?,
        index: frame.index.clone(),
        axis_metadata: frame.axis_metadata.clone(),
        inherit_frame_metadata: true,
        warning,
    })
}
