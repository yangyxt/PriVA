//! Opt-in original R7 L621-665, after unchanged Python C1-C5 and ordered R1-R6.
//! Source function 8aa29d9a4f83ff5629bbdf5e6279a43a374d71ca87fe04fe9db88253cc57d2ca.
//! Caller retains the original pandas buffers, Series/axis metadata and logs.
use super::control_multiply::RowIndex;
use super::ranking_compatibility::Compatibility;
use super::terminal_order::{Input, InputColumn};
use super::{SourceError, python_strip};
use crate::acmg_consequence::coerce_pext_value::{Value, coerce_pext_value};
use crate::acmg_variant_mechanism::variant_mechanism_masks::Cell;
use polars_core::prelude::*;

pub const TARGET: &str = "pext_sort_index";

pub struct Options {
    pub tissues: Cell,
    pub cutoff: Cell,
    pub floor: Cell,
    pub shape: Cell,
}

/// Exact boxed arguments supplied by the original pandas Series.map boundary.
/// Scalar origins (notably Python/NumPy bool and f32) must never be inferred
/// from a flattened object column. Arbitrary objects need original Python.
pub struct CoercionColumn<I> {
    pub position: usize,
    pub index: Vec<I>,
    pub values: Vec<Value>,
}
pub struct Modulation<I> {
    pub index: Vec<I>,
    pub dtype: &'static str,
    pub data: Vec<bool>,
    pub null: Vec<bool>,
}
impl<I: Clone> Modulation<I> {
    /// R4's actual high-ACMG mask, not its separate compatibility-penalty mask.
    pub fn from_compatibility(value: &Compatibility<I>) -> Result<Self, SourceError> {
        let data: Vec<bool> = value
            .role("pext_modulation_rows")
            .and_then(Column::bool)
            .map_err(|_| unsupported("R4 PEXT mask representation"))?
            .into_no_null_iter()
            .collect();
        Ok(Self {
            index: value.index.clone(),
            dtype: "bool",
            null: vec![false; data.len()],
            data,
        })
    }
}

#[derive(Debug, PartialEq)]
pub struct Selection {
    pub column: Option<String>,
    /// None means the source never created this local; Some([]) is different.
    pub tissue_list: Option<Vec<String>>,
    pub first_tissue_column: Option<String>,
    pub used_tissue: bool,
}
pub struct Pext<I> {
    pub selection: Selection,
    /// Original order: coerce/fill/clip, division/clip, power, product, sum,
    /// where and finally masked score multiplication. NaN is never Polars null.
    pub values: Option<DataFrame>,
    pub selected: Vec<bool>,
    pub source_column_position: Option<usize>,
    pub target_column_position: Option<usize>,
    pub score_column_position: Option<usize>,
    pub index: Vec<I>,
    pub axis_metadata: String,
}
fn unsupported(message: &str) -> SourceError {
    SourceError(
        "UnsupportedInput",
        format!("{message}; retain Python backend"),
    )
}

pub fn select_column<I>(frame: &Input<I>, tissues: &Cell) -> Result<Selection, SourceError> {
    let mut result = Selection {
        column: None,
        tissue_list: None,
        first_tissue_column: None,
        used_tissue: false,
    };
    match tissues {
        Cell::Null => (),
        Cell::Text(text) => {
            if !text.is_empty() {
                let list: Vec<_> = text
                    .split(',')
                    .map(python_strip)
                    .filter(|s| !s.is_empty())
                    .map(str::to_owned)
                    .collect();
                if let Some(first) = list.first() {
                    let name = format!("PEXT_{}", first.replace('-', "_"));
                    if frame.columns.iter().any(|c| c.name == name) {
                        result.column = Some(name.clone());
                        result.used_tissue = true;
                    }
                    result.first_tissue_column = Some(name);
                }
                result.tissue_list = Some(list);
            }
        }
        _ => return Err(unsupported("tissue truth/split representation")),
    }
    if result.column.is_none() && frame.columns.iter().any(|c| c.name == "PEXT_MEAN") {
        result.column = Some("PEXT_MEAN".into());
    }
    Ok(result)
}
fn position<I>(frame: &Input<I>, name: &str) -> Result<usize, SourceError> {
    let found: Vec<_> = frame
        .columns
        .iter()
        .enumerate()
        .filter(|(_, c)| c.name == name)
        .map(|(i, _)| i)
        .collect();
    match found.as_slice() {
        [i] => Ok(*i),
        _ => Err(unsupported("missing/ambiguous consumed R7 label")),
    }
}
fn parameter(value: &Cell) -> Result<f64, SourceError> {
    match value {
        Cell::Float(value) if value.is_finite() => Ok(*value),
        _ => Err(unsupported("R7 parameter type/nonfinite numerical route")),
    }
}
fn clip(value: f64) -> f64 {
    // pandas where-based clip retains an in-range negative zero.
    if value < 0. {
        0.
    } else if value > 1. {
        1.
    } else {
        value
    }
}
#[inline(never)]
fn source_power(base: f64, shape: f64) -> f64 {
    // 25 October 2026 CT: R7 passes a Series to np.power. pandas2.2.3
    // dispatches to operator.pow; NumPy1.26.4 number.c fast_scalar_power uses
    // positive/sqrt/square at these scalar exponents. Direct ndarray np.power
    // differs. Both old mismatches and the actual Series/installed dispatch
    // proof (807,012 values, opt0/opt2) remain retained; no blanket sqrt
    // substitution or numerical tolerance is allowed.
    if shape == 1. {
        base
    } else if shape == 0.5 {
        base.sqrt()
    } else if shape == 2. {
        base * base
    } else {
        base.powf(std::hint::black_box(shape))
    }
}

/// Qualified initial boundary: 1..=1M rows, flat opaque axes, actual boxed
/// primitive map arguments, finite valid f64 parameters, native f64 scores and
/// aligned bool/nullable-Boolean R4 masks. No-column skips unused validation.
/// Plan every value/assignment before changing buffers. Rejection requires
/// original L621-665 once on the current entry buffers, without earlier replay.
pub fn ranking_pext<I: Clone + PartialEq>(
    frame: &mut Input<I>,
    options: &Options,
    coercion: Option<&CoercionColumn<I>>,
    mask: &Modulation<I>,
    row_index: RowIndex,
) -> Result<Pext<I>, SourceError> {
    let selection = select_column(frame, &options.tissues)?;
    let mut result = Pext {
        selection,
        values: None,
        selected: Vec::new(),
        source_column_position: None,
        target_column_position: None,
        score_column_position: None,
        index: frame.index.clone(),
        axis_metadata: frame.axis_metadata.clone(),
    };
    let Some(name) = result.selection.column.as_deref() else {
        return Ok(result);
    };
    // Source priority is cutoff, floor, shape; all source errors are deferred.
    let cutoff = parameter(&options.cutoff)?;
    if cutoff <= 0. {
        return Err(unsupported("original cutoff failure"));
    }
    let floor = parameter(&options.floor)?;
    if !(0. ..=1.).contains(&floor) {
        return Err(unsupported("original floor failure"));
    }
    let shape = parameter(&options.shape)?;
    if shape <= 0. {
        return Err(unsupported("original shape failure"));
    }
    let n = frame.index.len();
    if n == 0
        || n > 1_000_000
        || row_index != RowIndex::Flat
        || frame.columns.iter().any(|c| c.values.len() != n)
        || mask.index != frame.index
        || mask.data.len() != n
        || mask.null.len() != n
        || !matches!(mask.dtype, "bool" | "boolean")
        || mask.dtype == "bool" && mask.null.iter().any(|x| *x)
    {
        return Err(unsupported(
            "empty/dtype/axis/mask or expression-size route",
        ));
    }
    let source_position = position(frame, name)?;
    let score_position = position(frame, "sort_index")?;
    let target_positions: Vec<_> = frame
        .columns
        .iter()
        .enumerate()
        .filter(|(_, c)| c.name == TARGET)
        .map(|(i, _)| i)
        .collect();
    if target_positions.len() > 1 {
        return Err(unsupported("duplicate R7 assigned label"));
    }
    let target_position = target_positions
        .first()
        .copied()
        .unwrap_or(frame.columns.len());
    let coercion = coercion.ok_or_else(|| unsupported("missing original map scalar identities"))?;
    if coercion.position != source_position
        || coercion.index != frame.index
        || coercion.values.len() != n
    {
        return Err(unsupported("unaligned original map arguments"));
    }
    if frame.columns[score_position].dtype != "float64" {
        return Err(unsupported("score dtype/promotion/warning route"));
    }
    let scores = frame.columns[score_position]
        .values
        .iter()
        .map(|v| match v {
            Cell::Float(value) => Ok(*value),
            _ => Err(unsupported("native f64 score storage")),
        })
        .collect::<Result<Vec<_>, _>>()?;
    let numeric = coercion
        .values
        .iter()
        .map(|value| {
            coerce_pext_value(value)
                .map(|value| clip(if value.is_nan() { 1. } else { value }))
                .map_err(|_| unsupported("original scalar parser failure"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let scaled: Vec<_> = numeric.iter().map(|v| clip(*v / cutoff)).collect();
    // Preserve separate operations and the qualified Series power dispatch;
    // no FMA, algebraic combination or pruning of np.where's unused branch.
    let powers: Vec<_> = scaled.iter().map(|v| source_power(*v, shape)).collect();
    let products: Vec<_> = powers.iter().map(|v| (1. - floor) * *v).collect();
    let sums: Vec<_> = products.iter().map(|v| floor + *v).collect();
    let multipliers: Vec<_> = numeric
        .iter()
        .zip(&sums)
        .map(|(v, low)| if *v >= cutoff { 1. } else { *low })
        .collect();
    let selected: Vec<_> = mask
        .data
        .iter()
        .zip(&mask.null)
        .map(|(raw, na)| *raw && !*na)
        .collect();
    let scores: Vec<_> = scores
        .iter()
        .zip(&multipliers)
        .zip(&selected)
        .map(|((score, factor), yes)| if *yes { *score * *factor } else { *score })
        .collect();
    let values = DataFrame::new(vec![
        Column::new("pext_numeric".into(), numeric),
        Column::new("scaled_pext".into(), scaled),
        Column::new("power".into(), powers),
        Column::new("product".into(), products),
        Column::new("sum".into(), sums),
        Column::new(TARGET.into(), multipliers.clone()),
        Column::new("sort_index".into(), scores.clone()),
    ])
    .map_err(|_| unsupported("Polars value construction"))?;
    let target = InputColumn {
        name: TARGET.into(),
        dtype: "float64".into(),
        values: multipliers.into_iter().map(Cell::Float).collect(),
    };
    if target_position == frame.columns.len() {
        frame.columns.push(target)
    } else {
        frame.columns[target_position] = target
    }
    frame.columns[score_position].values = scores.into_iter().map(Cell::Float).collect();
    result.values = Some(values);
    result.selected = selected;
    result.source_column_position = Some(source_position);
    result.target_column_position = Some(target_position);
    result.score_column_position = Some(score_position);
    Ok(result)
}
