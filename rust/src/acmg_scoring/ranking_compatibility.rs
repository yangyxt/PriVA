//! Original R4 L537–597 only. R6 owns compatibility multiplication.
//! Consume current R2 context after unchanged Python C1–C5/R3; never replay it.
use super::SourceError;
use super::pedigree_context::Context;
use super::terminal_order::{Input, InputColumn};
use crate::acmg_variant_mechanism::variant_mechanism_masks::{
    self as mechanism, Cell, Dtype, Masks,
};
use polars_core::prelude::*;

pub const TARGET: &str = "zygosity_inheritance_mechanism_compatibility";
pub const ROLES: [&str; 11] = [
    "chrom",
    "x_locus",
    "has_single_alt_lof_compatible",
    "exact_gof",
    "asserted_lof_effect",
    "dn_without_hi_or_ambiguous",
    "heterozygous_lof_incompatible",
    "high_acmg_lof",
    "compatibility_penalty_rows",
    "compatibility_penalty_rows_null",
    "pext_modulation_rows",
];

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MetadataOrigin {
    /// Restore opaque attrs/flags from this original source Series.
    SourceColumn(usize),
    /// The source constructed a fresh default Series.
    Fresh,
}

#[derive(Debug)]
pub struct Failure {
    pub error: SourceError,
    pub source_line: Option<u32>,
}
impl std::fmt::Display for Failure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.error.fmt(f)
    }
}
impl std::error::Error for Failure {}
fn fail(error: SourceError, source_line: Option<u32>) -> Failure {
    Failure { error, source_line }
}
fn unsupported(message: &str) -> Failure {
    fail(
        SourceError(
            "UnsupportedInput",
            format!("{message}; retain Python backend"),
        ),
        None,
    )
}

pub struct Compatibility<I> {
    pub mechanism: Masks<I>,
    pub derived: DataFrame,
    pub index: Vec<I>,
    pub axis_metadata: String,
    pub compatibility_column_position: usize,
    pub penalty_dtype: &'static str,
    pub log_counts: [usize; 5],
    pub log_message: String,
    profile: usize,
    effect: usize,
    score: usize,
    chrom: Option<usize>,
}
impl<I> Compatibility<I> {
    /// Alias roles resolve to the SAME column, not independent copied masks.
    pub fn role(&self, name: &str) -> PolarsResult<&Column> {
        match name {
            "exact_gof" => self.mechanism.frame.column("is_exact_gof"),
            "asserted_lof_effect" => self.mechanism.frame.column("is_predicted_lof"),
            "pext_modulation_rows" => self.derived.column("high_acmg_lof"),
            _ => self.derived.column(name),
        }
    }
    pub fn source_name(&self, name: &str) -> Option<&'static str> {
        if let Some(i) = mechanism::NAMES.iter().position(|value| *value == name) {
            return self.mechanism.source_series_names[i];
        }
        match name {
            "chrom" | "x_locus" if self.chrom.is_some() => Some("chrom"),
            "dn_without_hi_or_ambiguous" => Some("var_plausible_patho_mechs"),
            "high_acmg_lof" | "pext_modulation_rows" => Some("ACMG_quant_score"),
            _ => None,
        }
    }
    pub fn metadata_origin(&self, name: &str) -> MetadataOrigin {
        let column = match name {
            "chrom" | "x_locus" => self.chrom,
            // pd.to_numeric constructs a fresh Series even for an existing
            // numeric column; it does not inherit its attrs/duplicate flag.
            "is_exact_gof" | "is_exact_nonlof" | "exact_gof" | "is_exact_dn" => None,
            "is_predicted_lof"
            | "is_uncertain"
            | "asserted_lof_effect"
            | "heterozygous_lof_incompatible" => Some(self.effect),
            "high_acmg_lof"
            | "pext_modulation_rows"
            | "compatibility_penalty_rows"
            | "compatibility_penalty_rows_null" => Some(self.score),
            _ => Some(self.profile),
        };
        column.map_or(MetadataOrigin::Fresh, MetadataOrigin::SourceColumn)
    }
}

fn dtype(value: &str) -> Option<Dtype> {
    match value {
        "object" => Some(Dtype::Object),
        "string" => Some(Dtype::String),
        "float64" | "int64" | "uint64" | "bool" => Some(Dtype::Numpy),
        // Nullable helper/cache routes remain on the original current buffers.
        _ => None,
    }
}

/// Native numeric raw scores; no lossy large-integer conversion is required
/// for comparison with the source's fixed 0.89 threshold.
fn high(column: &InputColumn) -> Result<Vec<bool>, Failure> {
    column
        .values
        .iter()
        .map(|value| match (column.dtype.as_str(), value) {
            ("float64", Cell::Float(value)) => Ok(*value > 0.89),
            ("int64", Cell::Integer(value)) => Ok(*value >= 1),
            ("uint64", Cell::Unsigned(value)) => Ok(*value >= 1),
            ("uint64", Cell::Integer(value)) if *value >= 0 => Ok(*value >= 1),
            ("bool", Cell::Bool(value)) => Ok(*value),
            _ => Err(unsupported("raw ACMG score representation/promotion")),
        })
        .collect()
}

/// Preserve the original mechanism-helper position, required-column errors,
/// Boolean operation order, INFO counts and two compatibility assignments.
/// Only aligned bool/nullable-Boolean R2 masks (or empty object masks), native
/// numeric scores, declared primitive mechanism columns and textual/native
/// integer/bool chromosomes are qualified. No arbitrary pandas objects/caches.
/// Every rejection precedes mutation; fallback executes original R4 once on
/// current buffers. Prior warnings/C1–C5/R3 and R6 multiplication never run here.
pub fn ranking_compatibility<I: Clone + PartialEq>(
    frame: &mut Input<I>,
    context: &Context<I>,
) -> Result<Compatibility<I>, Failure> {
    let n = frame.index.len();
    if frame.columns.iter().any(|column| column.values.len() != n)
        || context.index != frame.index
        || context.mask.height() != n
        || !matches!(context.mask_dtype, "bool" | "boolean")
            && !(context.mask_dtype == "object" && n == 0)
    {
        return Err(unsupported("unaligned/unsupported R2 context"));
    }
    let used = [
        "var_plausible_patho_mechs",
        "variant_effect",
        "variant_lof_score",
        "variant_gof_score",
        "variant_dn_score",
    ];
    let mut columns = Vec::new();
    for column in &frame.columns {
        if used.contains(&column.name.as_str()) {
            let dtype = dtype(&column.dtype)
                .ok_or_else(|| unsupported("mechanism dtype/cache representation"))?;
            columns.push(mechanism::InputColumn {
                name: column.name.clone(),
                dtype,
                values: column.values.clone(),
            });
        }
    }
    // The qualified original helper executes here, after C1–C5 and R3.
    let masks = mechanism::variant_mechanism_masks(&mechanism::Input {
        columns,
        index: frame.index.clone(),
        index_names: Vec::new(),
    })
    .map_err(|error| fail(error, Some(537)))?;
    let position = |name: &str| frame.columns.iter().position(|column| column.name == name);
    if ["chrom", "ACMG_quant_score", TARGET].iter().any(|name| {
        frame
            .columns
            .iter()
            .filter(|column| column.name == *name)
            .count()
            > 1
    }) {
        return Err(unsupported("ambiguous R4 column label"));
    }
    let chrom_position = position("chrom");
    let mut chrom = Vec::with_capacity(n);
    if let Some(i) = chrom_position {
        let column = &frame.columns[i];
        for value in &column.values {
            let text = match (column.dtype.as_str(), value) {
                ("object" | "string", Cell::Text(value)) => value.clone(),
                ("object" | "string", Cell::Null) => String::new(),
                ("object", Cell::Float(value)) if value.is_nan() => String::new(),
                ("int64", Cell::Integer(value)) => value.to_string(),
                ("uint64", Cell::Unsigned(value)) => value.to_string(),
                ("uint64", Cell::Integer(value)) if *value >= 0 => value.to_string(),
                ("bool", Cell::Bool(value)) => if *value { "True" } else { "False" }.into(),
                _ => return Err(unsupported("chromosome fillna/string representation")),
            };
            chrom.push(text);
        }
        if !matches!(
            column.dtype.as_str(),
            "object" | "string" | "int64" | "uint64" | "bool"
        ) {
            return Err(unsupported("chromosome dtype, including empty columns"));
        }
    } else {
        chrom = vec![String::new(); n];
    }
    let mask_data = context
        .mask
        .column("value")
        .and_then(Column::bool)
        .map_err(|_| unsupported("R2 mask data representation"))?;
    let mask_null = context
        .mask
        .column("is_null")
        .and_then(Column::bool)
        .map_err(|_| unsupported("R2 mask null representation"))?;
    if mask_data.null_count() != 0
        || mask_null.null_count() != 0
        || context.mask_dtype != "boolean" && mask_null.into_no_null_iter().any(|value| value)
    {
        return Err(unsupported("R2 raw mask storage"));
    }
    let get = |name: &str| -> Vec<bool> {
        masks
            .frame
            .column(name)
            .unwrap()
            .bool()
            .unwrap()
            .into_no_null_iter()
            .collect()
    };
    let dom_lof = get("has_dom_lof_history");
    let dom_bare = get("has_dom_unresolved_history");
    let xr_lof = get("has_x_recessive_lof_history");
    let xr_bare = get("has_x_recessive_unresolved_history");
    let xd_lof = get("has_x_dominant_lof_history");
    let xd_bare = get("has_x_dominant_unresolved_history");
    let dom_dn = get("has_dom_dn_history");
    let xd_dn = get("has_x_dominant_dn_history");
    let gof = get("is_exact_gof");
    let lof = get("is_predicted_lof");
    let modern = get("modern_profile");
    // The captured CPython witness has no non-ASCII character whose uppercase
    // is solely CHRX. ASCII upper therefore preserves this exact X predicate;
    // the original chromosome text itself is returned unchanged.
    let x: Vec<bool> = chrom
        .iter()
        .map(|value| {
            let upper = value.to_ascii_uppercase();
            upper.strip_prefix("CHR").unwrap_or(&upper) == "X"
        })
        .collect();
    let single: Vec<bool> = (0..n)
        .map(|i| {
            dom_lof[i]
                || dom_bare[i]
                || x[i] && context.is_male && (xr_lof[i] || xr_bare[i] || xd_lof[i] || xd_bare[i])
                || x[i] && context.is_female && (xd_lof[i] || xd_bare[i])
        })
        .collect();
    let dn: Vec<bool> = (0..n)
        .map(|i| {
            (dom_dn[i] || xd_dn[i]) && !(dom_lof[i] || xd_lof[i]) && !(dom_bare[i] || xd_bare[i])
        })
        .collect();
    let incompatible: Vec<bool> = (0..n).map(|i| lof[i] && !gof[i] && !single[i]).collect();
    // Original L584 follows the compatibility-history/effect masks above.
    let score = position("ACMG_quant_score").ok_or_else(|| {
        fail(
            SourceError("KeyError", "'ACMG_quant_score'".into()),
            Some(584),
        )
    })?;
    let high = high(&frame.columns[score])?;
    if !matches!(
        frame.columns[score].dtype.as_str(),
        "float64" | "int64" | "uint64" | "bool"
    ) {
        return Err(unsupported("raw score dtype, including empty columns"));
    }
    let het: Vec<bool> = mask_data.into_no_null_iter().collect();
    let null: Vec<bool> = mask_null.into_no_null_iter().collect();
    // pandas2.2.3 kleene_and preserves raw data beneath null flags. Native
    // high/incompatible operands unmask False; do not normalize hidden bits.
    let penalty: Vec<bool> = (0..n)
        .map(|i| high[i] && het[i] && incompatible[i])
        .collect();
    let penalty_null: Vec<bool> = (0..n)
        .map(|i| high[i] && null[i] && incompatible[i])
        .collect();
    let count = |values: &[bool]| values.iter().filter(|value| **value).count();
    let counts = [
        count(&modern),
        count(&lof),
        count(&dn),
        count(&incompatible),
        (0..n).filter(|i| penalty[*i] && !penalty_null[*i]).count(),
    ];
    let log_message = format!(
        "Zygosity/inheritance/mechanism compatibility signals for ranking: hub_profile_rows={}, asserted_lof_effect={}, dominant_DN_without_HI_or_ambiguous={}, heterozygous_lof_incompatible={}, compatibility_penalty_rows={}",
        counts[0], counts[1], counts[2], counts[3], counts[4]
    );
    let derived = DataFrame::new(vec![
        Column::new("chrom".into(), chrom),
        Column::new("x_locus".into(), x),
        Column::new("has_single_alt_lof_compatible".into(), single),
        Column::new("dn_without_hi_or_ambiguous".into(), dn),
        Column::new("heterozygous_lof_incompatible".into(), incompatible),
        Column::new("high_acmg_lof".into(), high),
        Column::new("compatibility_penalty_rows".into(), penalty.clone()),
        Column::new(
            "compatibility_penalty_rows_null".into(),
            penalty_null.clone(),
        ),
    ])
    .map_err(|error| fail(SourceError("PolarsError", error.to_string()), None))?;
    let profile = position("var_plausible_patho_mechs").unwrap();
    let effect = position("variant_effect").unwrap();
    let target = position(TARGET).unwrap_or(frame.columns.len());
    let initialized = InputColumn {
        name: TARGET.into(),
        dtype: "float64".into(),
        values: vec![Cell::Float(1.0); n],
    };
    if target == frame.columns.len() {
        frame.columns.push(initialized);
    } else {
        frame.columns[target] = initialized;
    }
    for i in 0..n {
        if penalty[i] && !penalty_null[i] {
            frame.columns[target].values[i] = Cell::Float(0.9);
        }
    }
    Ok(Compatibility {
        mechanism: masks,
        derived,
        index: frame.index.clone(),
        axis_metadata: frame.axis_metadata.clone(),
        compatibility_column_position: target,
        penalty_dtype: if context.mask_dtype == "boolean" {
            "boolean"
        } else {
            "bool"
        },
        log_counts: counts,
        log_message,
        profile,
        effect,
        score,
        chrom: chrom_position,
    })
}
