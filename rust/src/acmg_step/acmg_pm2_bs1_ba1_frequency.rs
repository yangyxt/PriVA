//! `scripts/acmg_pm2_bs1_ba1_frequency.py` (799 lines): PM2, BS1 and BA1 from gnomAD and the
//! ClinVar allele-frequency statistics.
//!
//! Thread T5 (allele frequency, carrier observation and statistics). Converted:
//! `_gene_max_patho_hom_freq`, `control_false_neg_rate`, `control_false_positive_rate`,
//! `gnomAD_rare_AF`, `_rare_af_binomial`, `_rare_hom_binomial`, `PM2_criteria`, `BS1_criteria`
//! and `BA1_criteria`, with SciPy's binomial tails from `stats` and the condition masks from
//! `acmg_variant_mechanism::variant_condition_masks` (T3).
//!
//! The Python code works on whole columns; every quantity here is computed per row with the
//! same missing-value rules (`fillna`, `isna`, NaN comparisons false, the nullable-boolean
//! "undecidable" outcome of the two binomial helpers). The numeric columns are read as pandas'
//! `to_numeric` leaves int64 and float64 columns; any other column type there is refused, as are
//! empty tables (pandas' `apply` changes its result type on them).
use super::acmg_variant_mechanism::variant_condition_masks;
use super::stats;
use super::types::{AcmgError, HarnessCall, NamedMasks, StepOutput, not_converted};
use crate::cli::{Result, invalid};
use crate::json_text::Json;
use crate::pandas_table::{Column, ObjectValue, Table};
use crate::pickle_reader::{self, Value};
use crate::python_text::float_repr;
use std::collections::HashMap;
use std::path::Path;

fn refused(message: impl Into<String>) -> Box<dyn std::error::Error> {
    invalid(format!("refused (not reproduced): {}", message.into()))
}

fn key_error(name: &str) -> Box<dyn std::error::Error> {
    Box::new(AcmgError::new("KeyError", format!("'{name}'")))
}

// ------------------------------------------------------------------------------------------------
// Column access with pandas' rules

/// `pd.to_numeric(column, errors="coerce")` for int64 and float64 columns (NaN is missing).
fn numeric(column: &Column, name: &str) -> Result<Vec<f64>> {
    match column {
        Column::Float(values) => Ok(values.clone()),
        Column::Int(values) => values
            .iter()
            .map(|value| {
                if value.unsigned_abs() > 1 << 53 {
                    Err(refused(format!("integer {value} in {name} beyond 2^53")))
                } else {
                    Ok(*value as f64)
                }
            })
            .collect(),
        other => Err(refused(format!(
            "{name} is a {} column of {}, not numbers",
            other.dtype_name(),
            column_kind(other)
        ))),
    }
}

fn column_kind(column: &Column) -> &'static str {
    match column {
        Column::Text(_) => "text",
        Column::Object(_) => "mixed values",
        Column::Bool(_) | Column::BoolOrMissing(_) | Column::NullableBool(_) => "booleans",
        Column::NullableInt(_) => "nullable integers",
        Column::Int(_) | Column::Float(_) => "numbers",
    }
}

/// `df[name]` as numbers (a missing column is pandas' KeyError).
fn required_numbers(table: &Table, name: &str) -> Result<Vec<f64>> {
    match table.column_index(name) {
        Some(index) => numeric(&table.columns[index], name),
        None => Err(key_error(name)),
    }
}

/// `pd.to_numeric(df.get(name, pd.Series(default, index=df.index)), errors="coerce")`.
fn numbers_or(table: &Table, name: &str, default: f64) -> Result<Vec<f64>> {
    match table.column_index(name) {
        Some(index) => numeric(&table.columns[index], name),
        None => Ok(vec![default; table.rows]),
    }
}

fn fill_nan(values: &[f64], fill: f64) -> Vec<f64> {
    values
        .iter()
        .map(|value| if value.is_nan() { fill } else { *value })
        .collect()
}

/// `series.clip(low, high)` (NaN stays NaN).
fn clip(value: f64, low: f64, high: f64) -> f64 {
    if value.is_nan() {
        value
    } else {
        value.max(low).min(high)
    }
}

/// `series.astype(int)` of float values: truncation; NaN or infinity raise pandas' error.
fn astype_int(values: &[f64]) -> Result<Vec<i64>> {
    values
        .iter()
        .map(|value| {
            if !value.is_finite() {
                return Err(Box::new(AcmgError::new(
                    "IntCastingNaNError",
                    "Cannot convert non-finite values (NA or inf) to integer",
                )) as Box<dyn std::error::Error>);
            }
            let truncated = value.trunc();
            if !(-9.223_372_036_854_776e18..9.223_372_036_854_776e18).contains(&truncated) {
                return Err(refused(format!(
                    "{value} is outside the 64-bit integer range"
                )));
            }
            Ok(truncated as i64)
        })
        .collect()
}

/// `df.get(name, pd.Series("", index=df.index)).fillna("").astype(str)`.
fn text_or_empty(table: &Table, name: &str) -> Result<Vec<String>> {
    let Some(index) = table.column_index(name) else {
        return Ok(vec![String::new(); table.rows]);
    };
    let python_bool = |flag: bool| if flag { "True" } else { "False" }.to_owned();
    let float_text = |value: f64| {
        if value.is_nan() {
            String::new()
        } else {
            float_repr(value)
        }
    };
    Ok(match &table.columns[index] {
        Column::Text(values) => values
            .iter()
            .map(|value| value.clone().unwrap_or_default())
            .collect(),
        Column::Float(values) => values.iter().map(|value| float_text(*value)).collect(),
        Column::Int(values) => values.iter().map(|value| value.to_string()).collect(),
        Column::Bool(values) => values.iter().map(|flag| python_bool(*flag)).collect(),
        Column::BoolOrMissing(values) => values
            .iter()
            .map(|flag| flag.map(python_bool).unwrap_or_default())
            .collect(),
        Column::Object(values) => values
            .iter()
            .map(|value| match value {
                ObjectValue::Nan | ObjectValue::None => String::new(),
                ObjectValue::Text(text) => text.clone(),
                ObjectValue::Int(number) => number.to_string(),
                ObjectValue::Float(number) => float_text(*number),
                ObjectValue::Bool(flag) => python_bool(*flag),
            })
            .collect(),
        other => {
            return Err(refused(format!(
                "{name} is a {} column; fillna(\"\") is not reproduced there",
                other.dtype_name()
            )));
        }
    })
}

/// `str.upper()` / `str.lower()` of ASCII text (other text is refused).
fn ascii_case(values: &[String], name: &str, upper: bool) -> Result<Vec<String>> {
    values
        .iter()
        .map(|value| {
            if !value.is_ascii() {
                return Err(refused(format!("non-ASCII text in {name}")));
            }
            Ok(if upper {
                value.to_ascii_uppercase()
            } else {
                value.to_ascii_lowercase()
            })
        })
        .collect()
}

/// The loci of BS1 and BS2: `chrom` without its `chr` prefix, upper case.
pub(crate) struct Loci {
    pub(crate) x: Vec<bool>,
    pub(crate) y: Vec<bool>,
    pub(crate) mitochondrial: Vec<bool>,
    pub(crate) autosomal: Vec<bool>,
}

pub(crate) fn loci(table: &Table) -> Result<Loci> {
    let chrom = ascii_case(&text_or_empty(table, "chrom")?, "chrom", true)?;
    let chrom: Vec<&str> = chrom
        .iter()
        .map(|value| value.strip_prefix("CHR").unwrap_or(value))
        .collect();
    let x: Vec<bool> = chrom.iter().map(|value| *value == "X").collect();
    let y: Vec<bool> = chrom.iter().map(|value| *value == "Y").collect();
    let mitochondrial: Vec<bool> = chrom
        .iter()
        .map(|value| *value == "M" || *value == "MT")
        .collect();
    let autosomal = (0..chrom.len())
        .map(|row| !(x[row] || y[row] || mitochondrial[row]))
        .collect();
    Ok(Loci {
        x,
        y,
        mitochondrial,
        autosomal,
    })
}

/// `_series_or_empty("HPO_scope_review_required").str.lower().isin({"1", "true", "yes"})`.
pub(crate) fn scope_review_required(table: &Table) -> Result<Vec<bool>> {
    let text = ascii_case(
        &text_or_empty(table, "HPO_scope_review_required")?,
        "HPO_scope_review_required",
        false,
    )?;
    Ok(text
        .iter()
        .map(|value| matches!(value.as_str(), "1" | "true" | "yes"))
        .collect())
}

fn mask<'a>(masks: &'a NamedMasks, name: &str) -> Result<&'a [bool]> {
    masks.get(name)
}

// ------------------------------------------------------------------------------------------------
// The binomial helpers

/// What a threshold argument of the two binomial helpers is.
enum Threshold<'a> {
    Scalar(f64),
    PerRow(&'a [f64]),
}

impl Threshold<'_> {
    fn at(&self, row: usize) -> f64 {
        match self {
            Threshold::Scalar(value) => *value,
            Threshold::PerRow(values) => values[row],
        }
    }
}

/// The shared front of `control_false_neg_rate` and `control_false_positive_rate`: the cleaned
/// thresholds, allele numbers and observed allele counts, and which rows can be tested.
struct Prepared {
    thresholds: Vec<f64>,
    an: Vec<i64>,
    ac_observed: Vec<i64>,
    valid: Vec<bool>,
}

fn prepare(af: &[f64], allele_numbers: &[f64], threshold: &Threshold) -> Result<Prepared> {
    let rows = af.len();
    let thresholds: Vec<f64> = (0..rows)
        .map(|row| clip(threshold.at(row), 0.0, 1.0))
        .collect();
    if thresholds.iter().any(|value| value.is_nan()) {
        return Err(Box::new(AcmgError::new(
            "ValueError",
            "Non-numeric values found in af_threshold Series/array after coercion.",
        )));
    }
    let mut an = astype_int(&fill_nan(allele_numbers, 0.0))?;
    for value in an.iter_mut() {
        if *value < 0 {
            *value = 0;
        }
    }
    let mut ac_observed = Vec::with_capacity(rows);
    let mut valid = Vec::with_capacity(rows);
    for row in 0..rows {
        let frequency = clip(af[row], 0.0, 1.0);
        let ac_float = if frequency.is_nan() {
            0.0
        } else {
            frequency * an[row] as f64
        };
        let mut observed = if ac_float.is_nan() {
            0
        } else {
            ac_float.round_ties_even() as i64
        };
        observed = observed.min(an[row]);
        if an[row] == 0 {
            observed = 0;
        }
        ac_observed.push(observed);
        valid.push(an[row] > 0 && !frequency.is_nan() && !allele_numbers[row].is_nan());
    }
    Ok(Prepared {
        thresholds,
        an,
        ac_observed,
        valid,
    })
}

/// `control_false_neg_rate(af, an, af_threshold, alpha)[1]`: true where the allele count is
/// significantly above the threshold, `None` where the test cannot be run.
fn control_false_neg_rate(
    af: &[f64],
    allele_numbers: &[f64],
    threshold: &Threshold,
    alpha: f64,
) -> Result<Vec<Option<bool>>> {
    let prepared = prepare(af, allele_numbers, threshold)?;
    Ok((0..af.len())
        .map(|row| {
            if !prepared.valid[row] {
                return None;
            }
            let k = prepared.ac_observed[row];
            let pvalue = if k == 0 {
                1.0
            } else {
                stats::binom_sf(
                    (k - 1) as f64,
                    prepared.an[row] as f64,
                    prepared.thresholds[row],
                )
            };
            Some(pvalue <= alpha)
        })
        .collect())
}

/// `control_false_positive_rate(af, an, af_threshold, alpha)[1]`: true where the allele count is
/// significantly below the threshold, `None` where the test cannot be run.
fn control_false_positive_rate(
    af: &[f64],
    allele_numbers: &[f64],
    threshold: &Threshold,
    alpha: f64,
) -> Result<Vec<Option<bool>>> {
    let prepared = prepare(af, allele_numbers, threshold)?;
    Ok((0..af.len())
        .map(|row| {
            if !prepared.valid[row] {
                return None;
            }
            let threshold = prepared.thresholds[row];
            let pvalue = if threshold == 0.0 {
                1.0
            } else if threshold > 0.0 {
                stats::binom_cdf(
                    prepared.ac_observed[row] as f64,
                    prepared.an[row] as f64,
                    threshold,
                )
            } else {
                f64::NAN
            };
            if pvalue.is_nan() {
                None
            } else {
                Some(pvalue <= alpha)
            }
        })
        .collect())
}

/// The gnomAD columns the PM2 helpers read.
struct Gnomad {
    af: Vec<f64>,
    an: Vec<f64>,
    af_max: Vec<f64>,
    an_max: Vec<f64>,
}

fn gnomad(table: &Table) -> Result<Gnomad> {
    Ok(Gnomad {
        af: required_numbers(table, "gnomAD_joint_AF")?,
        an: required_numbers(table, "gnomAD_joint_AN")?,
        af_max: required_numbers(table, "gnomAD_joint_AF_max")?,
        an_max: required_numbers(table, "gnomAD_joint_AN_max")?,
    })
}

/// `_rare_af_binomial(df, cutoff)` with a per-row cutoff (already `to_numeric`).
fn rare_af_binomial(table: &Table, cutoff: &[f64]) -> Result<Vec<bool>> {
    let threshold: Vec<f64> = cutoff
        .iter()
        .map(|value| {
            let clipped = clip(*value, 0.0, 1.0);
            if clipped.is_nan() { 0.0 } else { clipped }
        })
        .collect();
    let gnomad = gnomad(table)?;
    let rare_joint =
        control_false_positive_rate(&gnomad.af, &gnomad.an, &Threshold::PerRow(&threshold), 0.01)?;
    Ok((0..table.rows)
        .map(|row| {
            let rare = rare_joint[row].unwrap_or(false);
            let reliable_max_pop = gnomad.an_max[row] >= 1.0 / threshold[row];
            let too_common = reliable_max_pop && gnomad.af_max[row] > threshold[row];
            let own_af_missing = gnomad.af_max[row].is_nan() && gnomad.af[row].is_nan();
            (rare && !too_common) || own_af_missing
        })
        .collect())
}

/// `_rare_hom_binomial(df.loc[rows], cutoff.loc[rows])` for the chosen rows.
fn rare_hom_binomial(table: &Table, rows: &[usize], cutoff: &[f64]) -> Result<Vec<bool>> {
    let threshold: Vec<f64> = rows
        .iter()
        .map(|row| {
            let clipped = clip(cutoff[*row], 0.0, 1.0);
            if clipped.is_nan() { 0.0 } else { clipped }
        })
        .collect();
    let an = required_numbers(table, "gnomAD_joint_AN")?;
    let nhomalt_xx = required_numbers(table, "gnomAD_nhomalt_XX")?;
    let nhomalt_xy = required_numbers(table, "gnomAD_nhomalt_XY")?;
    let an_max = required_numbers(table, "gnomAD_joint_AN_max")?;
    let nhomalt_max = required_numbers(table, "gnomAD_nhomalt_max")?;
    let af = required_numbers(table, "gnomAD_joint_AF")?;
    let af_max = required_numbers(table, "gnomAD_joint_AF_max")?;
    let pick = |values: &[f64]| -> Vec<f64> { rows.iter().map(|row| values[*row]).collect() };
    let n_joint: Vec<f64> = pick(&an).iter().map(|value| value / 2.0).collect();
    let hom_joint: Vec<f64> = rows
        .iter()
        .map(|row| fill_nan(&[nhomalt_xx[*row]], 0.0)[0] + fill_nan(&[nhomalt_xy[*row]], 0.0)[0])
        .collect();
    let hom_freq_joint: Vec<f64> = (0..rows.len())
        .map(|at| hom_joint[at] / n_joint[at])
        .collect();
    let rare_joint = control_false_positive_rate(
        &hom_freq_joint,
        &n_joint,
        &Threshold::PerRow(&threshold),
        0.01,
    )?;
    Ok(rows
        .iter()
        .enumerate()
        .map(|(at, row)| {
            let rare = rare_joint[at].unwrap_or(false);
            let reliable_max_pop = an_max[*row] / 2.0 >= 1.0 / threshold[at];
            let hom_freq_max = nhomalt_max[*row] / (an_max[*row] / 2.0);
            let too_common = reliable_max_pop && hom_freq_max > threshold[at];
            let own_af_missing = af_max[*row].is_nan() && af[*row].is_nan();
            (rare && !too_common) || own_af_missing
        })
        .collect())
}

/// `gnomAD_rare_AF(df, cutoff)` with a float cutoff.
fn gnomad_rare_af(table: &Table, cutoff: f64) -> Result<Vec<bool>> {
    let gnomad = gnomad(table)?;
    Ok((0..table.rows)
        .map(|row| {
            let an_max = if gnomad.an_max[row].is_nan() {
                1_000_000.0
            } else {
                gnomad.an_max[row]
            };
            if an_max >= 1.0 / cutoff {
                fill_nan(&[gnomad.af_max[row]], 0.0)[0] <= cutoff
            } else {
                fill_nan(&[gnomad.af[row]], 0.0)[0] <= cutoff
            }
        })
        .collect())
}

// ------------------------------------------------------------------------------------------------
// The ClinVar comparators

const HIGH_CONFIDENCE_STATUS: [(&str, i64); 3] = [
    ("practice_guideline", 4),
    ("reviewed_by_expert_panel", 3),
    ("criteria_provided,_multiple_submitters,_no_conflicts", 2),
];

/// `_gene_max_patho_hom_freq(df)`: gene -> the largest homozygote frequency among the table's
/// 2-star pathogenic ClinVar variants of the gene (NaN when every such frequency is NaN).
fn gene_max_patho_hom_freq(table: &Table) -> Result<HashMap<String, f64>> {
    let clnsig = text_or_empty(table, "CLNSIG")?;
    let clnrev = text_or_empty(table, "CLNREVSTAT")?;
    let rows: Vec<usize> = (0..table.rows)
        .filter(|row| {
            clnsig[*row].contains("athogenic")
                && HIGH_CONFIDENCE_STATUS
                    .iter()
                    .any(|(status, stars)| *status == clnrev[*row] && *stars >= 2)
        })
        .collect();
    let mut result: HashMap<String, f64> = HashMap::new();
    if rows.is_empty() {
        return Ok(result);
    }
    let an_max = numbers_or(table, "gnomAD_joint_AN_max", f64::NAN)?;
    let hom_max = numbers_or(table, "gnomAD_nhomalt_max", f64::NAN)?;
    let an = numbers_or(table, "gnomAD_joint_AN", f64::NAN)?;
    let nhomalt_xx = numbers_or(table, "gnomAD_nhomalt_XX", 0.0)?;
    let nhomalt_xy = numbers_or(table, "gnomAD_nhomalt_XY", 0.0)?;
    let genes = text_or_empty(table, "Gene")?;
    // groupby(gene).max(): the largest non-missing value of each group, NaN when none.
    let mut seen: HashMap<String, bool> = HashMap::new();
    for row in rows {
        let n_max = an_max[row] / 2.0;
        let n_joint = an[row] / 2.0;
        let hom_joint = fill_nan(&[nhomalt_xx[row]], 0.0)[0] + fill_nan(&[nhomalt_xy[row]], 0.0)[0];
        let frequency = if n_max > 10_000.0 {
            hom_max[row] / n_max
        } else {
            hom_joint / n_joint
        };
        let gene = genes[row].clone();
        let entry = result.entry(gene.clone()).or_insert(f64::NAN);
        let has_value = seen.entry(gene).or_insert(false);
        if !frequency.is_nan() && (!*has_value || frequency > *entry) {
            *entry = frequency;
            *has_value = true;
        }
    }
    Ok(result)
}

/// A gene's value from `dict.get(gene, 0)` mapped over a column, as `Series.map` types it:
/// int64 when every value is an int, float64 when they mix ints, floats and None.
enum Mapped {
    Int(i64),
    Float(f64),
}

fn mapped_column(values: Vec<Mapped>) -> Column {
    if values.iter().all(|value| matches!(value, Mapped::Int(_))) {
        Column::Int(
            values
                .into_iter()
                .map(|value| match value {
                    Mapped::Int(number) => number,
                    Mapped::Float(_) => unreachable!(),
                })
                .collect(),
        )
    } else {
        Column::Float(
            values
                .into_iter()
                .map(|value| match value {
                    Mapped::Int(number) => number as f64,
                    Mapped::Float(number) => number,
                })
                .collect(),
        )
    }
}

/// The text of the `Gene` cell as a dictionary key (`None` for a missing value, which matches no
/// key of the pickles).
fn gene_keys(table: &Table, name: &str) -> Result<Vec<Option<String>>> {
    let index = table.column_index(name).ok_or_else(|| key_error(name))?;
    Ok(match &table.columns[index] {
        Column::Text(values) => values.clone(),
        Column::Float(values) if values.iter().all(|value| value.is_nan()) => {
            vec![None; table.rows]
        }
        Column::Object(values) => values
            .iter()
            .map(|value| match value {
                ObjectValue::Text(text) => Ok(Some(text.clone())),
                ObjectValue::Nan | ObjectValue::None => Ok(None),
                _ => Err(refused(format!("{name} holds a value that is not text"))),
            })
            .collect::<Result<_>>()?,
        other => {
            return Err(refused(format!(
                "{name} is a {} column of {}",
                other.dtype_name(),
                column_kind(other)
            )));
        }
    })
}

fn load_pickle(path: &str) -> Result<Value> {
    let bytes = pickle_reader::read_file(Path::new(path), path.ends_with(".gz"))?;
    pickle_reader::load(&bytes)
}

fn dict_entries(value: &Value, what: &str) -> Result<Vec<(Value, Value)>> {
    match value {
        Value::Dict(dict) => Ok(dict.borrow().entries.clone()),
        other => Err(refused(format!(
            "{what} is a {}, not a dict",
            other.type_name()
        ))),
    }
}

/// `clinvar_patho_af_dict.get(gene, {"af": 0}).get("af", 0)` for every gene of the pickle whose
/// value is not None.
fn gene_af_lookup(path: &str) -> Result<HashMap<String, Mapped>> {
    let stat = load_pickle(path)?;
    let mut lookup = HashMap::new();
    for (key, value) in dict_entries(&stat, "the ClinVar pathogenic AF statistics")? {
        if matches!(value, Value::None) {
            continue;
        }
        let Value::Str(gene) = key else {
            // A key that is not text matches no cell of the Gene column.
            continue;
        };
        let entry = match &value {
            Value::Dict(dict) => match dict.borrow().get(&Value::Str("af".into()))? {
                None => Mapped::Int(0),
                Some(Value::Int(number)) => Mapped::Int(*number),
                Some(Value::Float(number)) => Mapped::Float(*number),
                Some(Value::None) => Mapped::Float(f64::NAN),
                Some(other) => {
                    return Err(refused(format!(
                        "the af of gene {gene} is a {}",
                        other.type_name()
                    )));
                }
            },
            other => {
                return Err(refused(format!(
                    "the statistics of gene {gene} are a {}, not a dict",
                    other.type_name()
                )));
            }
        };
        lookup.insert(gene.to_string(), entry);
    }
    Ok(lookup)
}

/// The exon allele-frequency statistics: transcript -> exon -> (median, ?, max, ...).
struct ExonStats {
    by_transcript: HashMap<String, HashMap<String, (f64, f64)>>,
}

fn exon_stats(path: &str) -> Result<ExonStats> {
    let stat = load_pickle(path)?;
    let mut by_transcript = HashMap::new();
    for (key, value) in dict_entries(&stat, "the ClinVar exon AF statistics")? {
        let Value::Str(transcript) = key else {
            continue;
        };
        let mut exons = HashMap::new();
        for (exon, entry) in dict_entries(&value, "an exon AF entry")? {
            let Value::Str(exon) = exon else {
                continue;
            };
            let Value::Tuple(items) = &entry else {
                return Err(refused(format!(
                    "exon {exon} of {transcript} holds a {}, not a tuple",
                    entry.type_name()
                )));
            };
            let number = |at: usize| -> Result<f64> {
                match items.get(at) {
                    Some(Value::Float(number)) => Ok(*number),
                    Some(other) => Err(refused(format!(
                        "item {at} of exon {exon} of {transcript} is a {}",
                        other.type_name()
                    ))),
                    None => Err(refused(format!(
                        "exon {exon} of {transcript} has {} items",
                        items.len()
                    ))),
                }
            };
            exons.insert(exon.to_string(), (number(0)?, number(2)?));
        }
        by_transcript.insert(transcript.to_string(), exons);
    }
    Ok(ExonStats { by_transcript })
}

// ------------------------------------------------------------------------------------------------
// The criteria

/// `PM2_criteria(df, clinvar_patho_af_stat, clinvar_patho_exon_af_stat,
/// gnomAD_extreme_rare_threshold)`: PM2 strengths (0, 1, 2) and the four columns it writes, in
/// its order: `clinvar_patho_gene_max_af`, `clinvar_patho_gene_max_hom_freq`,
/// `exon_patho_median_af`, `exon_patho_max_af`.
pub(crate) fn pm2_criteria(
    table: &Table,
    clinvar_patho_af_stat: &str,
    clinvar_patho_exon_af_stat: &str,
    gnomad_extreme_rare_threshold: f64,
) -> Result<StepOutput> {
    if table.rows == 0 {
        return Err(refused(
            "an empty table (pandas' apply returns a frame there)",
        ));
    }
    let gene_af_by_gene = gene_af_lookup(clinvar_patho_af_stat)?;
    let genes = gene_keys(table, "Gene")?;
    let gene_af_column = mapped_column(
        genes
            .iter()
            .map(
                |gene| match gene.as_ref().and_then(|g| gene_af_by_gene.get(g)) {
                    Some(Mapped::Int(number)) => Mapped::Int(*number),
                    Some(Mapped::Float(number)) => Mapped::Float(*number),
                    None => Mapped::Int(0),
                },
            )
            .collect(),
    );
    let mut work = table.clone();
    work.set_column("clinvar_patho_gene_max_af", gene_af_column.clone());
    let hom_by_gene = gene_max_patho_hom_freq(&work)?;
    let hom_column = mapped_column(
        genes
            .iter()
            .map(
                |gene| match gene.as_ref().and_then(|g| hom_by_gene.get(g)) {
                    Some(value) => Mapped::Float(*value),
                    None => Mapped::Int(0),
                },
            )
            .collect(),
    );
    work.set_column("clinvar_patho_gene_max_hom_freq", hom_column.clone());

    let exons = exon_stats(clinvar_patho_exon_af_stat)?;
    let features = lookup_keys(&work, "Feature")?;
    let exon_keys = lookup_keys(&work, "EXON")?;
    let mut median = Vec::with_capacity(table.rows);
    let mut max = Vec::with_capacity(table.rows);
    for row in 0..table.rows {
        let found = match (&features[row], &exon_keys[row]) {
            (Some(feature), Some(exon)) => exons
                .by_transcript
                .get(feature)
                .and_then(|per_exon| per_exon.get(exon)),
            _ => None,
        };
        median.push(found.map_or(f64::NAN, |entry| entry.0));
        max.push(found.map_or(f64::NAN, |entry| entry.1));
    }
    work.set_column("exon_patho_median_af", Column::Float(median.clone()));
    work.set_column("exon_patho_max_af", Column::Float(max.clone()));

    let masks = variant_condition_masks(&work)?;
    let has_dominant = mask(&masks, "has_dominant")?;
    let has_recessive = mask(&masks, "has_recessive")?;

    let gene_af = fill_nan(&numeric(&gene_af_column, "clinvar_patho_gene_max_af")?, 0.0);
    let gene_hom = fill_nan(
        &numeric(&hom_column, "clinvar_patho_gene_max_hom_freq")?,
        0.0,
    );
    let pm2_af_threshold: Vec<f64> = gene_af
        .iter()
        .map(|value| {
            if *value > 0.0 {
                value / 10.0
            } else {
                gnomad_extreme_rare_threshold
            }
        })
        .collect();
    let pm2_hom_threshold: Vec<f64> = gene_hom
        .iter()
        .map(|value| if *value > 0.0 { value / 10.0 } else { f64::NAN })
        .collect();

    let pm2_af_rare = rare_af_binomial(&work, &pm2_af_threshold)?;
    let mut pm2_hom_rare = vec![false; table.rows];
    let recessive_rows: Vec<usize> = (0..table.rows)
        .filter(|row| has_recessive[*row] && gene_hom[*row] > 0.0)
        .collect();
    if !recessive_rows.is_empty() {
        let rare = rare_hom_binomial(&work, &recessive_rows, &pm2_hom_threshold)?;
        for (at, row) in recessive_rows.iter().enumerate() {
            pm2_hom_rare[*row] = rare[at];
        }
    }
    let gnomad_af_rare: Vec<bool> = (0..table.rows)
        .map(|row| {
            let use_af_route = has_dominant[row] || !has_recessive[row];
            (pm2_af_rare[row] && use_af_route) || (pm2_hom_rare[row] && has_recessive[row])
        })
        .collect();

    let pm2_supporting_exon = rare_af_binomial(&work, &max)?;
    let gene_af_raw = numeric(&gene_af_column, "clinvar_patho_gene_max_af")?;
    let pm2_supporting_gene = rare_af_binomial(&work, &gene_af_raw)?;
    let pm2_moderate_exon = rare_af_binomial(&work, &median)?;
    let gnomad_absent = gnomad_rare_af(&work, 1e-7)?;

    let mut pm2 = vec![0i64; table.rows];
    for row in 0..table.rows {
        let mut supporting = if max[row] == 0.0 || max[row].is_nan() {
            pm2_supporting_gene[row]
        } else {
            pm2_supporting_exon[row]
        };
        if gene_af_raw[row] == 0.0 || gene_af_raw[row].is_nan() {
            supporting = gnomad_af_rare[row];
        }
        let moderate = if median[row] == 0.0 || median[row].is_nan() {
            gnomad_absent[row]
        } else {
            pm2_moderate_exon[row] && gnomad_af_rare[row]
        };
        if supporting {
            pm2[row] = 1;
        }
        if moderate {
            pm2[row] = 2;
        }
    }
    Ok(StepOutput {
        criteria: vec![("PM2", Column::Int(pm2))],
        columns: vec![
            ("clinvar_patho_gene_max_af".to_owned(), gene_af_column),
            ("clinvar_patho_gene_max_hom_freq".to_owned(), hom_column),
            ("exon_patho_median_af".to_owned(), Column::Float(median)),
            ("exon_patho_max_af".to_owned(), Column::Float(max)),
        ],
        ..StepOutput::default()
    })
}

/// The cells of `Feature` or `EXON` as `dict.get` keys inside `df.apply(axis=1)`: text, or `None`
/// for a value that matches no key of the exon pickle (missing, a number).
fn lookup_keys(table: &Table, name: &str) -> Result<Vec<Option<String>>> {
    let index = table.column_index(name).ok_or_else(|| key_error(name))?;
    Ok(match &table.columns[index] {
        Column::Text(values) => values.clone(),
        Column::Float(_) | Column::Int(_) | Column::Bool(_) | Column::BoolOrMissing(_) => {
            vec![None; table.rows]
        }
        Column::Object(values) => values
            .iter()
            .map(|value| match value {
                ObjectValue::Text(text) => Some(text.clone()),
                _ => None,
            })
            .collect(),
        other => {
            return Err(refused(format!(
                "{name} is a {} column",
                other.dtype_name()
            )));
        }
    })
}

/// `BS1_criteria(df, expected_incidence, pm2_criteria, return_frequency_components=True)`:
/// BS1 strengths (0 or 3), with the three route masks passed on (`disease_incidence_bs1`,
/// `gene_pathogenic_af_bs1`, `gene_pathogenic_af_blocked_by_pm2`). BS1 writes no column.
pub(crate) fn bs1_criteria(
    table: &Table,
    expected_incidence: f64,
    pm2_criteria: Option<&[i64]>,
) -> Result<StepOutput> {
    let rows = table.rows;
    let pm2_present: Vec<bool> = match pm2_criteria {
        None => vec![false; rows],
        Some(values) => {
            if values.len() != rows {
                return Err(Box::new(AcmgError::new(
                    "ValueError",
                    "pm2_criteria must have the same number of rows as df",
                )));
            }
            values.iter().map(|value| *value > 0).collect()
        }
    };
    let loci = loci(table)?;
    let masks = variant_condition_masks(table)?;
    let scope_review = scope_review_required(table)?;

    let gene_af = match table.column_index("clinvar_patho_gene_max_af") {
        Some(index) => fill_nan(
            &numeric(&table.columns[index], "clinvar_patho_gene_max_af")?,
            0.0,
        ),
        None => vec![0.0; rows],
    };
    let gene_hom = match table.column_index("clinvar_patho_gene_max_hom_freq") {
        Some(index) => fill_nan(
            &numeric(&table.columns[index], "clinvar_patho_gene_max_hom_freq")?,
            0.0,
        ),
        None => {
            let hom_by_gene = gene_max_patho_hom_freq(table)?;
            let genes = gene_keys(table, "Gene")?;
            genes
                .iter()
                .map(|gene| {
                    let value = gene
                        .as_ref()
                        .and_then(|g| hom_by_gene.get(g))
                        .copied()
                        .unwrap_or(0.0);
                    if value.is_nan() { 0.0 } else { value }
                })
                .collect()
        }
    };
    let bs1_af_threshold: Vec<f64> = gene_af
        .iter()
        .map(|value| {
            if *value > 0.0 {
                *value
            } else {
                expected_incidence
            }
        })
        .collect();
    let bs1_hom_threshold: Vec<f64> = gene_hom
        .iter()
        .map(|value| {
            if *value > 0.0 {
                *value
            } else {
                expected_incidence
            }
        })
        .collect();

    let bonferroni_alpha = 0.01 / 2.0;
    let af = required_numbers(table, "gnomAD_joint_AF")?;
    let af_max = required_numbers(table, "gnomAD_joint_AF_max")?;
    let an = required_numbers(table, "gnomAD_joint_AN")?;
    let an_max = required_numbers(table, "gnomAD_joint_AN_max")?;
    let common_max = control_false_neg_rate(
        &af_max,
        &an_max,
        &Threshold::PerRow(&bs1_af_threshold),
        bonferroni_alpha,
    )?;
    let common_joint = control_false_neg_rate(
        &af,
        &an,
        &Threshold::PerRow(&bs1_af_threshold),
        bonferroni_alpha,
    )?;
    let max_af_larger_incidence: Vec<bool> = (0..rows)
        .map(|row| {
            let either = common_max[row].unwrap_or(false) || common_joint[row].unwrap_or(false);
            let both_na = common_max[row].is_none() && common_joint[row].is_none();
            either || (both_na && af[row] > bs1_af_threshold[row])
        })
        .collect();

    // The homozygote tests: counts and halved allele numbers as integers.
    let nhomalt_max = required_numbers(table, "gnomAD_nhomalt_max")?;
    let nhomalt_xx = required_numbers(table, "gnomAD_nhomalt_XX")?;
    let nhomalt_xy = required_numbers(table, "gnomAD_nhomalt_XY")?;
    let as_count = |values: Vec<f64>| -> Result<Vec<i64>> {
        astype_int(
            &values
                .iter()
                .map(|value| if value.is_nan() { 0.0 } else { value.max(0.0) })
                .collect::<Vec<f64>>(),
        )
    };
    let n_max = as_count(an_max.iter().map(|value| value / 2.0).collect())?;
    let x_max = as_count(nhomalt_max.clone())?;
    let n_joint = as_count(an.iter().map(|value| value / 2.0).collect())?;
    let x_joint = as_count(
        (0..rows)
            .map(|row| nhomalt_xx[row] + nhomalt_xy[row])
            .collect(),
    )?;
    let hom_bs1 = |n: &[i64], x: &[i64]| -> Vec<bool> {
        (0..rows)
            .map(|row| {
                let pvalue = if n[row] > 0 {
                    if x[row] == 0 {
                        1.0
                    } else {
                        stats::binom_sf(
                            (x[row] - 1).max(0) as f64,
                            n[row] as f64,
                            bs1_hom_threshold[row],
                        )
                    }
                } else {
                    f64::NAN
                };
                let hom_freq = if n[row] == 0 {
                    f64::NAN
                } else {
                    x[row] as f64 / n[row] as f64
                };
                hom_freq > 0.0001 && pvalue <= bonferroni_alpha
            })
            .collect()
    };
    let hom_max = hom_bs1(&n_max, &x_max);
    let hom_joint = hom_bs1(&n_joint, &x_joint);

    let has_recessive = mask(&masks, "has_recessive")?;
    let has_dominant = mask(&masks, "has_dominant")?;
    let has_x_recessive = mask(&masks, "has_x_linked_recessive")?;
    let has_x_dominant = mask(&masks, "has_x_linked_dominant")?;
    let has_mendelian = mask(&masks, "has_mendelian")?;
    let has_non_monogenic = mask(&masks, "has_non_monogenic")?;
    let has_non_mendelian = mask(&masks, "has_non_mendelian")?;
    let has_incomplete = mask(&masks, "has_incomplete_penetrance")?;
    let has_y_linked = mask(&masks, "has_y_linked")?;
    let has_mitochondrial = mask(&masks, "has_mitochondrial")?;

    let frequency_above = |af_column: &str, an_column: &str| -> Result<Vec<bool>> {
        let frequency = fill_nan(&numbers_or(table, af_column, 0.0)?, 0.0);
        let number = fill_nan(&numbers_or(table, an_column, 0.0)?, 0.0);
        let significant = control_false_neg_rate(
            &frequency,
            &number,
            &Threshold::Scalar(expected_incidence),
            0.01,
        )?;
        Ok((0..rows)
            .map(|row| significant[row].unwrap_or(frequency[row] > expected_incidence))
            .collect())
    };
    let x_xy_af_larger = frequency_above("gnomAD_joint_AF_XY", "gnomAD_joint_AN_XY")?;
    let x_xx_af_larger = frequency_above("gnomAD_joint_AF_XX", "gnomAD_joint_AN_XX")?;

    let greater_than_clinvar =
        control_false_neg_rate(&af_max, &an_max, &Threshold::PerRow(&gene_af), 0.01)?;
    let af_xy = fill_nan(&numbers_or(table, "gnomAD_joint_AF_XY", 0.0)?, 0.0);
    let af_xx = fill_nan(&numbers_or(table, "gnomAD_joint_AF_XX", 0.0)?, 0.0);

    let mut bs1 = vec![0i64; rows];
    let mut disease_incidence = vec![false; rows];
    let mut gene_pathogenic = vec![false; rows];
    let mut blocked = vec![false; rows];
    for row in 0..rows {
        let valid_model = has_mendelian[row]
            && !has_non_monogenic[row]
            && !has_non_mendelian[row]
            && !has_incomplete[row]
            && !scope_review[row];
        let max_af = max_af_larger_incidence[row];
        let autosomal_dominant = loci.autosomal[row] && has_dominant[row] && max_af && valid_model;
        let x_dominant_requirement = has_x_dominant[row] && !has_recessive[row];
        let x_linked_recessive =
            loci.x[row] && has_x_recessive[row] && x_xy_af_larger[row] && valid_model;
        let x_linked_dominant =
            loci.x[row] && x_dominant_requirement && x_xx_af_larger[row] && valid_model;
        let autosomal_recessive = loci.autosomal[row]
            && has_recessive[row]
            && (hom_max[row] || hom_joint[row])
            && valid_model;
        let y_linked = loci.y[row] && has_y_linked[row] && max_af && valid_model;
        let mitochondrial =
            loci.mitochondrial[row] && has_mitochondrial[row] && max_af && valid_model;
        let clinvar_route = match greater_than_clinvar[row] {
            Some(flag) => flag,
            None if gene_af[row] > 0.0 => af[row] > gene_af[row],
            None => false,
        } && gene_af[row] > 0.0;
        let basic = af_max[row] > 0.0001 || af[row] > 0.0001;
        disease_incidence[row] =
            ((autosomal_dominant || autosomal_recessive || y_linked || mitochondrial) && basic)
                || (x_linked_recessive && af_xy[row] > 0.0001)
                || (x_linked_dominant && af_xx[row] > 0.0001);
        let candidate = clinvar_route && valid_model && basic;
        blocked[row] = candidate && pm2_present[row];
        gene_pathogenic[row] = candidate && !pm2_present[row];
        if disease_incidence[row] || gene_pathogenic[row] {
            bs1[row] = 3;
        }
    }
    Ok(StepOutput {
        criteria: vec![("BS1", Column::Int(bs1))],
        passed_on: vec![
            (
                "disease_incidence_bs1".to_owned(),
                Column::Bool(disease_incidence),
            ),
            (
                "gene_pathogenic_af_bs1".to_owned(),
                Column::Bool(gene_pathogenic),
            ),
            (
                "gene_pathogenic_af_blocked_by_pm2".to_owned(),
                Column::Bool(blocked),
            ),
        ],
        ..StepOutput::default()
    })
}

/// `BA1_criteria(anno_df)`: 5 where the max-population frequency is significantly above 5%
/// (or, where that cannot be tested, the joint frequency is above 5%).
pub(crate) fn ba1_criteria(table: &Table) -> Result<StepOutput> {
    let af_max = required_numbers(table, "gnomAD_joint_AF_max")?;
    let an_max = required_numbers(table, "gnomAD_joint_AN_max")?;
    let common = control_false_neg_rate(&af_max, &an_max, &Threshold::Scalar(0.05), 0.01)?;
    let af = required_numbers(table, "gnomAD_joint_AF")?;
    let ba1 = (0..table.rows)
        .map(|row| {
            let flag = common[row].unwrap_or(af[row] > 0.05);
            if flag { 5 } else { 0 }
        })
        .collect();
    Ok(StepOutput {
        criteria: vec![("BA1", Column::Int(ba1))],
        ..StepOutput::default()
    })
}

// ------------------------------------------------------------------------------------------------
// Comparison harness

fn option_text<'a>(options: &'a Json, name: &str) -> Result<&'a str> {
    match options.get(name) {
        Some(Json::Str(text)) => Ok(text),
        _ => Err(invalid(format!("option {name} (text) is required"))),
    }
}

fn option_number(options: &Json, name: &str, default: f64) -> Result<f64> {
    match options.get(name) {
        None | Some(Json::Null) => Ok(default),
        Some(Json::Float(number)) => Ok(*number),
        Some(Json::Int(number)) => Ok(*number as f64),
        Some(Json::Str(text)) if text == "nan" => Ok(f64::NAN),
        _ => Err(invalid(format!("option {name} must be a number"))),
    }
}

/// `pm2_criteria` as a list of integers, or None.
pub(crate) fn option_ints(options: &Json, name: &str) -> Result<Option<Vec<i64>>> {
    match options.get(name) {
        None | Some(Json::Null) => Ok(None),
        Some(Json::Array(items)) => items
            .iter()
            .map(|item| match item {
                Json::Int(number) => i64::try_from(*number).map_err(|_| invalid("integer")),
                _ => Err(invalid(format!("{name} holds integers"))),
            })
            .collect::<Result<Vec<i64>>>()
            .map(Some),
        _ => Err(invalid(format!("option {name} must be a list"))),
    }
}

/// The input table with the step's written columns applied and its criteria and passed-on
/// values appended (`__NAME`), as the Python captures write them.
pub(crate) fn output_frame(table: &Table, output: StepOutput) -> Table {
    let mut frame = table.clone();
    for (name, column) in output.columns {
        frame.set_column(&name, column);
    }
    for (name, column) in output.criteria {
        frame.set_column(&format!("__{name}"), column);
    }
    for (name, column) in output.passed_on {
        frame.set_column(&format!("__{name}"), column);
    }
    frame
}

/// `priva acmg-step acmg_pm2_bs1_ba1_frequency <function>`:
/// - `PM2_criteria`, options `clinvar_patho_af_stat`, `clinvar_patho_exon_af_stat`,
///   `gnomAD_extreme_rare_threshold` (default 0.0001);
/// - `BS1_criteria`, options `expected_incidence` (default 0.0001), `pm2_criteria` (list or null);
/// - `BA1_criteria`.
///
/// The frame returned is the input with the columns written and `__PM2`, `__BS1` (and its three
/// route masks) or `__BA1` appended.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let options = &call.options;
    let output = match function {
        "PM2_criteria" => pm2_criteria(
            &call.table,
            option_text(options, "clinvar_patho_af_stat")?,
            option_text(options, "clinvar_patho_exon_af_stat")?,
            option_number(options, "gnomAD_extreme_rare_threshold", 0.0001)?,
        )?,
        "BS1_criteria" => {
            let pm2 = option_ints(options, "pm2_criteria")?;
            bs1_criteria(
                &call.table,
                option_number(options, "expected_incidence", 0.0001)?,
                pm2.as_deref(),
            )?
        }
        "BA1_criteria" => ba1_criteria(&call.table)?,
        _ => return Err(not_converted("acmg_pm2_bs1_ba1_frequency", function)),
    };
    Ok(output_frame(&call.table, output))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frame(columns: Vec<(&str, Column)>) -> Table {
        let rows = columns.first().map_or(0, |(_, column)| column.len());
        let (names, columns) = columns
            .into_iter()
            .map(|(name, column)| (name.to_owned(), column))
            .unzip();
        Table {
            names,
            columns,
            rows,
        }
    }

    #[test]
    fn prepared_counts_round_half_to_even_and_stay_within_an() {
        let prepared = prepare(
            &[0.25, 0.5, f64::NAN, 2.0, 0.1],
            &[10.0, 3.0, 100.0, 4.0, f64::NAN],
            &Threshold::Scalar(0.01),
        )
        .unwrap();
        // 2.5 -> 2, 1.5 -> 2, NaN -> 0, clipped AF 1 x 4 -> 4, AN missing -> 0.
        assert_eq!(prepared.ac_observed, vec![2, 2, 0, 4, 0]);
        assert_eq!(prepared.valid, vec![true, true, false, true, false]);
    }

    #[test]
    fn undecidable_rows_are_none() {
        let result = control_false_neg_rate(
            &[0.2, f64::NAN, 0.0],
            &[1000.0, 1000.0, 0.0],
            &Threshold::Scalar(0.05),
            0.01,
        )
        .unwrap();
        assert_eq!(result, vec![Some(true), None, None]);
        let error = control_false_neg_rate(&[0.1], &[10.0], &Threshold::Scalar(f64::NAN), 0.01)
            .unwrap_err()
            .to_string();
        assert!(error.contains("Non-numeric values"), "{error}");
    }

    #[test]
    fn ba1_falls_back_to_the_joint_frequency() {
        let table = frame(vec![
            (
                "gnomAD_joint_AF_max",
                Column::Float(vec![0.2, f64::NAN, f64::NAN, 0.06]),
            ),
            (
                "gnomAD_joint_AN_max",
                Column::Float(vec![5000.0, 100.0, f64::NAN, 50.0]),
            ),
            (
                "gnomAD_joint_AF",
                Column::Float(vec![0.0, 0.06, 0.01, 0.06]),
            ),
        ]);
        let output = ba1_criteria(&table).unwrap();
        assert_eq!(output.criteria[0].1, Column::Int(vec![5, 5, 0, 0]));
    }

    #[test]
    fn text_columns_follow_astype_str() {
        let table = frame(vec![(
            "HPO_scope_review_required",
            Column::Object(vec![
                ObjectValue::Int(1),
                ObjectValue::Float(1.0),
                ObjectValue::Text("Yes".into()),
                ObjectValue::Bool(true),
                ObjectValue::Nan,
            ]),
        )]);
        assert_eq!(
            scope_review_required(&table).unwrap(),
            vec![true, false, true, true, false]
        );
    }
}
