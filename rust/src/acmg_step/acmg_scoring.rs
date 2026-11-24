//! `scripts/acmg_scoring.py` (836 lines): the scoring tail: criteria matrix and ClinGen override,
//! the posterior per row, `compute_control_common`, and the glue that runs the existing posterior
//! and ranking modules in the Python order.
//!
//! Thread T3 (exact-allele matching and the scoring tail). Converted: `summarize_acmg_criteria`
//! (matrix construction and validation, the ClinGen map, the summary column and the override of
//! the matrix rows), through `crate::acmg_scoring`'s `create_criteria_summary` and
//! `translate_strength_level`; the orchestrator's posterior step (`criteria_matrix.apply(
//! calculate_posterior_probability, axis=1)` and the two columns it writes), through
//! `crate::acmg_scoring::calculate_posterior_probability`; and `compute_control_common`, which
//! calls `stats::binom_cdf` and `stats::beta_ppf` (thread T5; until they land, a row whose control
//! counts are tested stops at their placeholder); and `sort_and_rank_variants` (stages R1 to R8:
//! the protein-coding filter (`kept_rows`, with the Rust-only rows of listed genes and of high
//! SpliceAI2 scores), the proband from the pedigree, the control-common penalty, the
//! zygosity, inheritance and mechanism compatibility penalty from the mechanism masks, the
//! dispensable and relevant gene lists, the pext modulation through the earlier
//! `coerce_pext_value`, and the final order and rank).
use super::acmg_pp3_bp4_bp7_insilico::pandas_masks::numeric_column;
use super::acmg_ps2_pm6_denovo::{family_rows, harness_ped, python_str};
use super::acmg_variant_mechanism::variant_mechanism_masks;
use super::stats;
use super::types::{AcmgError, CRITERIA_ORDER, HarnessCall, not_converted};
use crate::acmg_consequence::coerce_pext_value::{Value as PextValue, coerce_pext_value};
use crate::acmg_scoring::SourceError;
use crate::acmg_scoring::calculate_posterior_probability::{
    Probability, calculate_posterior_probability,
};
use crate::acmg_scoring::create_criteria_summary::create_criteria_summary;
use crate::acmg_scoring::translate_strength_level::translate_strength_level;
use crate::acmg_variant_mechanism::variant_mechanism_masks::Cell as MaskCell;
use crate::cache_convert::clingen_map::ClingenCache;
use crate::cli::{Result, invalid};
use crate::json_text::Json;
use crate::pandas_table::{Cell, Column, ObjectValue, Table, TextColumn};
use crate::pandas_table::{Header, ReadOptions, UseColumns, read_table_with};
use crate::python_text::python_strip;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::path::Path;

/// `priva acmg-step acmg_scoring <function>`: runs one converted function on a frame.
///
/// `summarize_acmg_criteria`: the frame holds the table's columns and the criteria as columns
/// named `criterion:NAME`; options `clingen_map_pkl` (text or null) and
/// `apply_clingen_override`. The result is the criteria matrix with the table's
/// `ACMG_criteria` appended as `__ACMG_criteria`.
/// `calculate_posterior_probability`: the frame is the criteria matrix; the result holds
/// `ACMG_quant_score` and `ACMG_class`.
/// `compute_control_common`: option `proband` (text or null); the result is the table.
/// `sort_and_rank_variants`: options `ped_frame` (a frame file, or absent for `ped_df=None`),
/// `fam_name`, `pext_tissues`, `pext_low_expression_cutoff`, `pext_penalty_floor`,
/// `pext_penalty_shape`, `relevant_gene_list`, `dispensable_gene_list` (null for None), and
/// the Rust-only `keep_listed_genes` (default false) and `spliceai_keep_score` (default
/// none), without which the rows are Python's; the result is the sorted table.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    match function {
        "summarize_acmg_criteria" => {
            let mut table = call.table.clone();
            let mut criteria = Vec::new();
            let mut keep = Vec::new();
            for (index, name) in table.names.iter().enumerate() {
                match name.strip_prefix("criterion:") {
                    Some(criterion) => {
                        let Some(known) = CRITERIA_ORDER.iter().find(|c| **c == criterion) else {
                            return Err(invalid(format!("unknown criterion {criterion}")));
                        };
                        criteria.push((*known, table.columns[index].clone()));
                    }
                    None => keep.push(index),
                }
            }
            table.names = keep
                .iter()
                .map(|index| table.names[*index].clone())
                .collect();
            table.columns = keep
                .iter()
                .map(|index| table.columns[*index].clone())
                .collect();
            let clingen = match call.options.get("clingen_map_pkl") {
                Some(Json::Str(path)) => Some(path.as_str()),
                _ => None,
            };
            let apply = !matches!(
                call.options.get("apply_clingen_override"),
                Some(Json::Bool(false))
            );
            let mut matrix = summarize_acmg_criteria(&mut table, &criteria, clingen, apply)?;
            let summary = table.column("ACMG_criteria")?.clone();
            matrix.set_column("__ACMG_criteria", summary);
            Ok(matrix)
        }
        "calculate_posterior_probability" => {
            let (scores, classes) = posterior_columns(&call.table)?;
            Ok(Table {
                names: vec!["ACMG_quant_score".into(), "ACMG_class".into()],
                columns: vec![scores, classes],
                rows: call.table.rows,
            })
        }
        "compute_control_common" => {
            let mut table = call.table.clone();
            let proband = match call.options.get("proband") {
                Some(Json::Str(name)) => Some(name.as_str()),
                _ => None,
            };
            compute_control_common(&mut table, proband)?;
            Ok(table)
        }
        "sort_and_rank_variants" => {
            let ped = match call.options.get("ped_frame") {
                Some(Json::Str(_)) => Some(harness_ped(call)?),
                _ => None,
            };
            let text = |name: &str| match call.options.get(name) {
                Some(Json::Str(value)) => Some(value.as_str()),
                _ => None,
            };
            let number = |name: &str, default: f64| match call.options.get(name) {
                Some(Json::Float(value)) => *value,
                Some(Json::Int(value)) => *value as f64,
                _ => default,
            };
            let options = RankingOptions {
                pext_tissues: text("pext_tissues").unwrap_or(""),
                pext_low_expression_cutoff: number("pext_low_expression_cutoff", 0.1),
                pext_penalty_floor: number("pext_penalty_floor", 0.8),
                pext_penalty_shape: number("pext_penalty_shape", 0.5),
                relevant_gene_list: text("relevant_gene_list"),
                dispensable_gene_list: text("dispensable_gene_list"),
                phenotype_gene_ranking: text("phenotype_gene_ranking"),
                keep_listed_genes: matches!(
                    call.options.get("keep_listed_genes"),
                    Some(Json::Bool(true))
                ),
                spliceai_keep_score: match call.options.get("spliceai_keep_score") {
                    Some(Json::Float(value)) => Some(*value),
                    _ => None,
                },
            };
            sort_and_rank_variants(call.table.clone(), ped.as_ref(), text("fam_name"), &options)
        }
        _ => Err(not_converted("acmg_scoring", function)),
    }
}

fn source_error(error: SourceError) -> Box<dyn std::error::Error> {
    Box::new(AcmgError::new(error.0, error.1))
}

/// A criterion strength as Python's `int(value)` reads a matrix cell: an integer,
/// or a boolean (numpy's True is 1). Other kinds are refused.
fn strength(column: &Column, row: usize, name: &str) -> Result<i64> {
    match column.cell(row) {
        Cell::Int(number) => Ok(number),
        Cell::Bool(flag) => Ok(i64::from(flag)),
        _ => Err(invalid(format!(
            "refused (not reproduced): criterion {name} holds a {} value",
            column.dtype_name()
        ))),
    }
}

/// `repr(path)` of a text path, as the FileNotFoundError message shows it.
fn python_repr_path(path: Option<&str>) -> Result<String> {
    match path {
        None => Ok("None".to_owned()),
        Some(text)
            if text
                .chars()
                .all(|c| (' '..='~').contains(&c) && c != '\'' && c != '\\') =>
        {
            Ok(format!("'{text}'"))
        }
        Some(text) => Err(invalid(format!(
            "refused (not reproduced): repr() of the ClinGen map path {text:?}"
        ))),
    }
}

/// `summarize_acmg_criteria(df, criteria_dict, clingen_map_pkl, apply_clingen_override)`:
/// writes `ACMG_criteria` into the table and returns the criteria matrix (the 28
/// criteria in `CRITERIA_ORDER`, absent ones as int64 zeros, then `variant_id`).
///
/// With a non-empty ClinGen map, each criterion column is rebuilt for the overridden
/// rows through a `.loc` assignment, which changes some columns to `object` the way
/// pandas does (see the comment there); the step's real run has all 28 as `object`.
pub(crate) fn summarize_acmg_criteria(
    table: &mut Table,
    criteria: &[(&str, Column)],
    clingen_map_pkl: Option<&str>,
    apply_clingen_override: bool,
) -> Result<Table> {
    let rows = table.rows;
    let mut columns: Vec<Column> = Vec::with_capacity(CRITERIA_ORDER.len());
    for name in CRITERIA_ORDER {
        let column = match criteria.iter().find(|(criterion, _)| *criterion == name) {
            Some((_, values)) => values.clone(),
            None => Column::Int(vec![0; rows]),
        };
        if column.len() != rows {
            return Err(invalid(format!(
                "refused (not reproduced): criterion {name} has {} values for {rows} rows",
                column.len()
            )));
        }
        columns.push(column);
    }
    let variant_id = table.column("variant_id")?.clone();
    if rows == 0 {
        return Err(invalid(
            "refused (not reproduced): summarize_acmg_criteria of an empty table",
        ));
    }
    for (name, column) in CRITERIA_ORDER.iter().zip(&columns) {
        let valid = match column {
            Column::Int(values) => values.iter().all(|value| (0..=5).contains(value)),
            Column::Bool(_) => true,
            _ => {
                return Err(invalid(format!(
                    "refused (not reproduced): criterion {name} is a {} column",
                    column.dtype_name()
                )));
            }
        };
        if !valid {
            return Err(Box::new(AcmgError::new(
                "ValueError",
                format!("Column {name} has non-integer values"),
            )));
        }
    }

    let clingen_map = if apply_clingen_override {
        let path = clingen_map_pkl.filter(|path| !path.is_empty() && Path::new(path).exists());
        let Some(path) = path else {
            return Err(Box::new(AcmgError::new(
                "FileNotFoundError",
                format!(
                    "ClinGen override is enabled but the evidence map is missing: {}. Build it \
                     with scripts/prepare_clingen_map.py, or pass apply_clingen_override=False to \
                     benchmark PriVA against ClinGen-curated variants.",
                    python_repr_path(clingen_map_pkl)?
                ),
            )));
        };
        Some(ClingenCache::open(path)?)
    } else {
        None
    };

    let evidence: Vec<Option<&str>> = (0..rows)
        .map(|row| match variant_id.cell(row) {
            Cell::Text(id) => clingen_map
                .as_ref()
                .and_then(|cache| cache.get(id).flatten()),
            _ => None,
        })
        .collect();
    let mut summaries: Vec<Option<String>> = Vec::with_capacity(rows);
    for (row, evidence) in evidence.iter().enumerate() {
        let strengths: Vec<(&str, i64)> = CRITERIA_ORDER
            .iter()
            .zip(&columns)
            .map(|(name, column)| Ok((*name, strength(column, row, name)?)))
            .collect::<Result<_>>()?;
        summaries.push(create_criteria_summary(&strengths, *evidence).map_err(source_error)?);
    }
    let summary_all_missing = summaries.iter().all(Option::is_none);
    let summary_column = if summary_all_missing {
        Column::Float(vec![f64::NAN; rows])
    } else {
        Column::Text(summaries.iter().map(Option::as_deref).collect())
    };
    table.set_column("ACMG_criteria", summary_column);

    if clingen_map.as_ref().is_some_and(|cache| cache.len() != 0) {
        let overridden: Vec<bool> = evidence
            .iter()
            .map(|value| value.is_some_and(|text| !python_strip(text).is_empty()))
            .collect();
        let any = overridden.iter().any(|flag| *flag);
        for (name, column) in CRITERIA_ORDER.iter().zip(columns.iter_mut()) {
            let mut updates = Vec::new();
            for row in (0..rows).filter(|row| overridden[*row]) {
                updates.push((
                    row,
                    translate_strength_level(name, summaries[row].as_deref())
                        .map_err(source_error)?,
                ));
            }
            // pandas assigns `.map(...)` of the overridden rows: an int64 Series when a row is
            // overridden, else an empty Series of the summary column's type (object when any
            // summary is text, float64 when all are NaN). int64 holds int64 and an empty float64,
            // and becomes object for an empty object Series; bool becomes object in every case.
            let to_object = match &*column {
                Column::Int(_) => !any && !summary_all_missing,
                _ => true,
            };
            let mut values: Vec<ObjectValue> = match &*column {
                Column::Int(values) => values
                    .iter()
                    .map(|value| ObjectValue::Int(*value))
                    .collect(),
                Column::Bool(values) => values
                    .iter()
                    .map(|value| ObjectValue::Bool(*value))
                    .collect(),
                _ => unreachable!("criteria are int64 or bool here"),
            };
            for (row, value) in updates {
                values[row] = ObjectValue::Int(value);
            }
            *column = if to_object {
                Column::Object(values)
            } else {
                Column::Int(
                    values
                        .iter()
                        .map(|value| match value {
                            ObjectValue::Int(number) => *number,
                            _ => unreachable!("an int64 column"),
                        })
                        .collect(),
                )
            };
        }
    }

    let mut names: Vec<String> = CRITERIA_ORDER
        .iter()
        .map(|name| (*name).to_owned())
        .collect();
    names.push("variant_id".to_owned());
    columns.push(variant_id);
    Ok(Table {
        names,
        columns,
        rows,
    })
}

/// The orchestrator's posterior step: `zip(*criteria_matrix.apply(
/// calculate_posterior_probability, axis=1))` as the columns `ACMG_quant_score`
/// and `ACMG_class`. Only columns named with a leading `P` or `B` count, as in Python.
pub(crate) fn posterior_columns(matrix: &Table) -> Result<(Column, Column)> {
    if matrix.rows == 0 {
        return Err(Box::new(AcmgError::new(
            "ValueError",
            "not enough values to unpack (expected 2, got 0)",
        )));
    }
    let counted: Vec<(&str, &Column)> = matrix
        .names
        .iter()
        .zip(&matrix.columns)
        .filter(|(name, _)| name.starts_with('P') || name.starts_with('B'))
        .map(|(name, column)| (name.as_str(), column))
        .collect();
    let mut scores = Vec::with_capacity(matrix.rows);
    let mut classes = Vec::with_capacity(matrix.rows);
    for row in 0..matrix.rows {
        let mut strengths = Vec::with_capacity(counted.len());
        for (name, column) in &counted {
            strengths.push((*name, strength(column, row, name)?));
        }
        if strengths
            .iter()
            .any(|(_, value)| !(0..=5).contains(value) && *value > 0)
        {
            return Err(invalid(
                "refused (not reproduced): a strength outside 0..5 in the posterior",
            ));
        }
        let (score, class) = calculate_posterior_probability(&strengths);
        scores.push(score);
        classes.push(Some(class.to_owned()));
    }
    // A tuple of floats and the integer 0 (BA1) becomes float64; integers alone stay int64.
    let scores = if scores
        .iter()
        .all(|score| matches!(score, Probability::Integer(_)))
    {
        Column::Int(
            scores
                .iter()
                .map(|score| match score {
                    Probability::Integer(value) => *value,
                    Probability::Float(_) => unreachable!("all integers"),
                })
                .collect(),
        )
    } else {
        Column::Float(scores.iter().map(Probability::as_f64).collect())
    };
    Ok((scores, Column::Text(classes.into())))
}

/// `pd.to_numeric(df.get(name), errors="coerce")` for the column types the control
/// columns have (`None` when the column is absent: Python's `to_numeric(None)` is
/// reached only through `df.get` and would raise, so an absent column is refused
/// unless the four are all absent).
fn control_numbers(table: &Table, name: &str) -> Result<Option<Vec<f64>>> {
    let Ok(column) = table.column(name) else {
        return Ok(None);
    };
    match column {
        Column::Int(values) => Ok(Some(values.iter().map(|value| *value as f64).collect())),
        Column::Float(values) => Ok(Some(values.clone())),
        Column::NullableInt(values) => Ok(Some(
            values
                .iter()
                .map(|value| value.map_or(f64::NAN, |number| number as f64))
                .collect(),
        )),
        _ => Err(invalid(format!(
            "refused (not reproduced): control column {name} is {}",
            column.dtype_name()
        ))),
    }
}

/// `compute_control_common(df, proband, cutoff=0.01, conf_level=0.999)`: writes
/// `control_common` and `control_common_index`.
pub(crate) fn compute_control_common(table: &mut Table, proband: Option<&str>) -> Result<()> {
    const CUTOFF: f64 = 0.01;
    const CONF_LEVEL: f64 = 0.999;
    let rows = table.rows;
    let names = ["control_AC", "control_AN", "control_AF", "control_nhomalt"];
    if !names.iter().any(|name| table.has_column(name)) {
        table.set_column("control_common", Column::Bool(vec![false; rows]));
        table.set_column("control_common_index", Column::Float(vec![1.0; rows]));
        return Ok(());
    }
    let mut values = Vec::new();
    for name in names {
        match control_numbers(table, name)? {
            Some(numbers) => values.push(numbers),
            None => {
                return Err(invalid(format!(
                    "refused (not reproduced): {name} is absent while another control column exists"
                )));
            }
        }
    }
    let [ac, an, af, nhomo]: [Vec<f64>; 4] = values.try_into().expect("four columns");
    let valid_an: Vec<bool> = an
        .iter()
        .map(|value| !value.is_nan() && *value > 0.0)
        .collect();
    let is_hom: Vec<bool> = match proband.and_then(|name| table.column(name).ok()) {
        Some(column) => (0..rows)
            .map(|row| {
                let text = match column.cell(row) {
                    Cell::Text(text) => text.to_owned(),
                    Cell::Missing => "nan".to_owned(),
                    _ => {
                        return Err(invalid(
                            "refused (not reproduced): a genotype that is not text",
                        ));
                    }
                };
                let gt = text.split(':').next().unwrap_or("");
                Ok(["1/1", "1|1", "1"].contains(&gt))
            })
            .collect::<Result<_>>()?,
        None => vec![false; rows],
    };
    let use_nhomo: Vec<bool> = (0..rows)
        .map(|row| valid_an[row] && is_hom[row] && !nhomo[row].is_nan())
        .collect();
    let use_ac: Vec<bool> = (0..rows)
        .map(|row| valid_an[row] && !use_nhomo[row] && !ac[row].is_nan())
        .collect();
    let use_af: Vec<bool> = (0..rows)
        .map(|row| valid_an[row] && !use_nhomo[row] && !use_ac[row] && !af[row].is_nan())
        .collect();
    // `.clip(lower=...).astype(int)`: truncation of finite values; `round()` is half to even.
    let as_int = |value: f64| -> Result<i64> {
        if !value.is_finite() || value.abs() >= 9.2e18 {
            return Err(invalid(
                "refused (not reproduced): a control count that does not fit an integer",
            ));
        }
        Ok(value.trunc() as i64)
    };
    let mut common = vec![false; rows];
    for row in 0..rows {
        if use_nhomo[row] {
            let k = as_int(nhomo[row].max(0.0))?;
            let n = as_int((an[row] / 2.0).round_ties_even().max(1.0))?;
            common[row] = stats::binom_cdf(k as f64, n as f64, CUTOFF) > CONF_LEVEL;
        }
    }
    for row in 0..rows {
        if use_ac[row] {
            let k = as_int(ac[row].max(0.0))?;
            let n = as_int(an[row].round_ties_even().max(1.0))?;
            common[row] = stats::binom_cdf(k as f64, n as f64, CUTOFF) > CONF_LEVEL;
        }
    }
    for row in 0..rows {
        if use_af[row] {
            common[row] = af[row] > CUTOFF;
        }
    }
    // `common[m] = af[m] > cutoff` assigns a Series; when it covers only some rows pandas
    // aligns it with NaN elsewhere and `common` becomes an object column. numpy then
    // refuses that object array as an index as soon as one row is common.
    let af_partial = use_af.iter().any(|flag| *flag) && !use_af.iter().all(|flag| *flag);
    if af_partial {
        let indexed = (0..rows).any(|row| common[row] && !an[row].is_nan() && an[row] > 0.0);
        if indexed {
            return Err(Box::new(AcmgError::new(
                "IndexError",
                "arrays used as indices must be of integer (or boolean) type",
            )));
        }
        table.set_column(
            "control_common",
            Column::BoolOrMissing(common.iter().map(|flag| Some(*flag)).collect()),
        );
        table.set_column("control_common_index", Column::Float(vec![1.0; rows]));
        return Ok(());
    }
    table.set_column("control_common", Column::Bool(common.clone()));
    let mut index = Vec::with_capacity(rows);
    for row in 0..rows {
        let ac_safe = if ac[row].is_nan() { 0.0 } else { ac[row] };
        let an_safe = if an[row].is_nan() { 0.0 } else { an[row] };
        index.push(if common[row] {
            let confident = if an_safe > 0.0 {
                stats::beta_ppf(0.01, ac_safe + 1.0, an_safe - ac_safe + 1.0)
            } else {
                0.0
            };
            (1.0 - confident).clamp(0.01, 0.5)
        } else {
            1.0
        });
    }
    table.set_column("control_common_index", Column::Float(index));
    Ok(())
}

/// The keyword arguments of `sort_and_rank_variants` after `fam_name`, as the step
/// passes them (`None` for a gene list stands for Python's `None`: no list).
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct RankingOptions<'a> {
    pub(crate) pext_tissues: &'a str,
    pub(crate) pext_low_expression_cutoff: f64,
    pub(crate) pext_penalty_floor: f64,
    pub(crate) pext_penalty_shape: f64,
    pub(crate) relevant_gene_list: Option<&'a str>,
    pub(crate) dispensable_gene_list: Option<&'a str>,
    /// The family's gene ranking for the phenotype gene bonus (see
    /// `phenotype_gene_indices`). Python has no such argument: `None` ranks as Python does.
    pub(crate) phenotype_gene_ranking: Option<&'a str>,
    /// Keep the rows of every gene in the gene ranking and the relevant gene list, whatever
    /// the transcript's biotype (see `kept_rows`). Python has no such rule.
    pub(crate) keep_listed_genes: bool,
    /// Keep a variant without a protein-coding row when a row's SpliceAI score is at least
    /// this (see `kept_rows`; `None`: never). Python has no such rule.
    pub(crate) spliceai_keep_score: Option<f64>,
}

/// The SpliceAI2 score from which a variant without a protein-coding row is kept: between
/// the publisher's balanced (0.25) and high-precision (0.5) thresholds.
pub(crate) const SPLICEAI2_KEEP_SCORE: f64 = 0.3;

/// The weight of the phenotype gene bonus at the family's top gene, the same as the
/// relevant gene list's weight.
const PHENOTYPE_GENE_TOP_WEIGHT: f64 = 1.2;

/// `phenotype_gene_index` of a gene from its share of the family's top `total_score`:
/// `1 + 0.2 * clip(1 + log10(share) / 2, 0, 1)`. The top gene gets 1.2, a gene at a tenth
/// of the top score 1.1, and a gene at a hundredth of it or below 1.0 (no gain).
fn phenotype_gene_index(share: f64) -> f64 {
    if share.is_nan() || share <= 0.0 {
        return 1.0;
    }
    let level = (1.0 + share.log10() / 2.0).clamp(0.0, 1.0);
    1.0 + (PHENOTYPE_GENE_TOP_WEIGHT - 1.0) * level
}

/// The phenotype gene bonus of each gene in a gene ranking table written by
/// `hpo-gene-ranking rank-genes` (`<sample>.gene_ranked.tsv`) for the proband's HPO terms:
/// `phenotype_gene_index` of `total_score / top`, where `top` is the largest `total_score`
/// of the table. A gene listed twice keeps its larger bonus; a table whose top score is not
/// above zero gives no gene a gain. This is a Rust-only revision of the ranking.
fn phenotype_gene_indices(path: &str) -> Result<HashMap<String, f64>> {
    let options = ReadOptions {
        all_text: true,
        keep_default_na: false,
        use_columns: Some(UseColumns::Names(vec![
            "gene_symbol".to_owned(),
            "total_score".to_owned(),
        ])),
        ..ReadOptions::default()
    };
    let ranking = read_table_with(Path::new(path), &options)
        .map_err(|error| refuse(format!("the gene ranking {path}: {error}")))?;
    let symbols = ranking.column("gene_symbol")?;
    let scores = ranking.column("total_score")?;
    let mut genes = Vec::with_capacity(ranking.rows);
    let mut top = 0.0_f64;
    for row in 0..ranking.rows {
        let symbol = match symbols.cell(row) {
            Cell::Text(text) => text.trim(),
            _ => "",
        };
        if symbol.is_empty() {
            continue;
        }
        let score = match scores.cell(row) {
            Cell::Text(text) => text.trim().parse::<f64>().ok(),
            _ => None,
        }
        .filter(|score| !score.is_nan())
        .ok_or_else(|| {
            refuse(format!(
                "the gene ranking {path}: total_score of {symbol} is not a number"
            ))
        })?;
        top = top.max(score);
        genes.push((symbol.to_owned(), score));
    }
    let mut indices = HashMap::with_capacity(genes.len());
    for (symbol, score) in genes {
        let index = if top > 0.0 {
            phenotype_gene_index(score / top)
        } else {
            1.0
        };
        let kept = indices.entry(symbol).or_insert(index);
        if index > *kept {
            *kept = index;
        }
    }
    Ok(indices)
}

/// Each row's value in `values` for its `SYMBOL`, 1.0 for a symbol not listed.
fn symbol_values(table: &Table, values: &HashMap<String, f64>) -> Result<Vec<f64>> {
    let symbols = table.column("SYMBOL")?;
    Ok((0..table.rows)
        .map(|row| match symbols.cell(row) {
            Cell::Text(text) => values.get(text).copied().unwrap_or(1.0),
            _ => 1.0,
        })
        .collect())
}

fn python_error(
    python_type: &'static str,
    message: impl Into<String>,
) -> Box<dyn std::error::Error> {
    Box::new(AcmgError::new(python_type, message))
}

fn refuse(message: impl std::fmt::Display) -> Box<dyn std::error::Error> {
    invalid(format!("refused (not reproduced): {message}"))
}

/// A float64 column's values (NaN for missing); other column types are refused.
fn float_values(table: &Table, name: &str) -> Result<Vec<f64>> {
    match table.column(name)? {
        Column::Float(values) => Ok(values.clone()),
        other => Err(refuse(format!(
            "{name} is {} rather than float64",
            other.dtype_name()
        ))),
    }
}

/// `str(value).strip().lower()` of a pedigree cell.
fn sex_text(cell: Option<Cell<'_>>) -> String {
    let text = match cell {
        None => String::new(),
        Some(cell) => python_str(cell).unwrap_or_else(|| "nan".to_owned()),
    };
    python_strip(&text).to_lowercase()
}

/// R2: the proband, its sex flags and `proband_het` (`df[proband].str.count("1") == 1`).
fn pedigree_context(
    table: &Table,
    ped: Option<&Table>,
    fam_name: Option<&str>,
) -> Result<(Option<String>, bool, bool, Vec<bool>)> {
    let (Some(ped), Some(fam_name)) = (ped, fam_name) else {
        return Ok((None, false, false, vec![false; table.rows]));
    };
    let family = family_rows(ped, fam_name)?;
    let phenotypes = ped.column("Phenotype")?;
    let first = family.into_iter().find(|row| match phenotypes.cell(*row) {
        Cell::Text(text) => text == "2",
        Cell::Int(number) => number == 2,
        Cell::Float(number) => number == 2.0,
        _ => false,
    });
    let Some(row) = first else {
        return Err(python_error(
            "IndexError",
            "single positional indexer is out-of-bounds",
        ));
    };
    let proband = match ped.column("IndividualID")?.cell(row) {
        Cell::Text(text) => text.to_owned(),
        other => return Err(refuse(format!("a proband ID that is not text ({other:?})"))),
    };
    let sex = sex_text(ped.column("Sex").ok().map(|column| column.cell(row)));
    let male = ["1", "1.0", "m", "male"].contains(&sex.as_str());
    let female = ["2", "2.0", "f", "female"].contains(&sex.as_str());
    let not_text = || {
        python_error(
            "AttributeError",
            "Can only use .str accessor with string values!",
        )
    };
    let het = match table.column(&proband)? {
        Column::Text(values) => values
            .iter()
            .map(|value| value.is_some_and(|text| text.matches('1').count() == 1))
            .collect(),
        Column::Object(_) => {
            return Err(refuse(format!(
                "the proband column {proband} mixes kinds of values"
            )));
        }
        _ => return Err(not_text()),
    };
    Ok((Some(proband), male, female, het))
}

/// `pd.read_table(path, header=None, names=["SYMBOL"], comment="#")` as a set of
/// texts. A missing value in the list is refused (pandas would match missing symbols).
fn gene_list(path: &str) -> Result<HashSet<String>> {
    let options = ReadOptions {
        header: Header::Names(vec!["SYMBOL".to_owned()]),
        comment: Some(b'#'),
        ..ReadOptions::default()
    };
    let list = read_table_with(Path::new(path), &options)?;
    let column = list.column("SYMBOL")?;
    let mut genes = HashSet::new();
    for row in 0..list.rows {
        match column.cell(row) {
            Cell::Text(text) => {
                genes.insert(text.to_owned());
            }
            Cell::Missing => {
                return Err(refuse(format!("a missing value in the gene list {path}")));
            }
            // A number or boolean never equals a symbol text.
            _ => {}
        }
    }
    Ok(genes)
}

/// `df["SYMBOL"].isin(genes)`.
fn symbol_mask(table: &Table, genes: &HashSet<String>) -> Result<Vec<bool>> {
    let symbols = table.column("SYMBOL")?;
    Ok((0..table.rows)
        .map(|row| matches!(symbols.cell(row), Cell::Text(text) if genes.contains(text)))
        .collect())
}

/// `Series.clip(lower, upper)`: NaN stays; a value that is not `>= lower` becomes `lower`,
/// then one that is not `<= upper` becomes `upper`.
fn clip(value: f64, lower: f64, upper: f64) -> f64 {
    if value.is_nan() {
        return value;
    }
    let value = if value >= lower { value } else { lower };
    if value <= upper { value } else { upper }
}

/// `np.power(series, shape)`: pandas hands a Series power to `operator.pow`, where numpy
/// 1.26.4 uses `sqrt`, the value itself or a square for the exponents 0.5, 1 and 2 (the
/// finding of `crate::acmg_scoring::ranking_pext`, checked there on 807,012 values).
fn series_power(base: f64, shape: f64) -> f64 {
    if shape == 1.0 {
        base
    } else if shape == 0.5 {
        base.sqrt()
    } else if shape == 2.0 {
        base * base
    } else {
        base.powf(std::hint::black_box(shape))
    }
}

/// `df[column].map(coerce_pext_value)` through the earlier scalar module.
fn pext_values(table: &Table, name: &str) -> Result<Vec<f64>> {
    let column = table.column(name)?;
    let value = |row: usize| -> Result<PextValue> {
        Ok(match column {
            Column::Bool(values) => PextValue::NumpyBool(values[row]),
            Column::NullableInt(_) | Column::NullableBool(_) => {
                return Err(refuse(format!(
                    "the pext column {name} is {}",
                    column.dtype_name()
                )));
            }
            Column::Float(values) => PextValue::Primitive(MaskCell::Float(values[row])),
            Column::Object(values) => match &values[row] {
                ObjectValue::None => PextValue::Primitive(MaskCell::Null),
                ObjectValue::Nan => PextValue::Primitive(MaskCell::Float(f64::NAN)),
                ObjectValue::Text(text) => PextValue::Primitive(MaskCell::Text(text.clone())),
                ObjectValue::Int(number) => PextValue::Primitive(MaskCell::Integer(*number)),
                ObjectValue::Float(number) => PextValue::Primitive(MaskCell::Float(*number)),
                ObjectValue::Bool(flag) => PextValue::Primitive(MaskCell::Bool(*flag)),
            },
            _ => match column.cell(row) {
                Cell::Missing => PextValue::Primitive(MaskCell::Float(f64::NAN)),
                Cell::Text(text) => PextValue::Primitive(MaskCell::Text(text.to_owned())),
                Cell::Int(number) => PextValue::Primitive(MaskCell::Integer(number)),
                Cell::Float(number) => PextValue::Primitive(MaskCell::Float(number)),
                Cell::Bool(flag) => PextValue::Primitive(MaskCell::Bool(flag)),
            },
        })
    };
    (0..table.rows)
        .map(|row| coerce_pext_value(&value(row)?).map_err(source_error))
        .collect()
}

/// One key column of the final sort and of the variant groups.
enum KeyValues {
    /// NaN is a missing value.
    Numbers(Vec<f64>),
    Texts(TextColumn),
}

impl KeyValues {
    fn of(table: &Table, name: &str) -> Result<KeyValues> {
        Ok(match table.column(name)? {
            Column::Float(values) => KeyValues::Numbers(values.clone()),
            Column::Int(values) => {
                if values.iter().any(|value| value.unsigned_abs() > 1 << 53) {
                    return Err(refuse(format!("{name} holds an integer past 2**53")));
                }
                KeyValues::Numbers(values.iter().map(|value| *value as f64).collect())
            }
            Column::Bool(values) => KeyValues::Numbers(
                values
                    .iter()
                    .map(|value| f64::from(u8::from(*value)))
                    .collect(),
            ),
            Column::Text(values) => KeyValues::Texts(values.clone()),
            other => {
                return Err(refuse(format!(
                    "the sort key {name} is a {} column of mixed or nullable values",
                    other.dtype_name()
                )));
            }
        })
    }

    fn missing(&self, row: usize) -> bool {
        match self {
            KeyValues::Numbers(values) => values[row].is_nan(),
            KeyValues::Texts(values) => values.get(row).is_none(),
        }
    }

    /// Ascending order of two present values (-0.0 equals 0.0, as in pandas' categories).
    fn compare(&self, left: usize, right: usize) -> Ordering {
        match self {
            KeyValues::Numbers(values) => values[left]
                .partial_cmp(&values[right])
                .unwrap_or(Ordering::Equal),
            KeyValues::Texts(values) => values.get(left).cmp(&values.get(right)),
        }
    }

    /// The value as a group key (`-0.0` and `0.0` together).
    fn group_key(&self, row: usize) -> String {
        match self {
            KeyValues::Numbers(values) => {
                let value = if values[row] == 0.0 { 0.0 } else { values[row] };
                format!("n{:016x}", value.to_bits())
            }
            KeyValues::Texts(values) => format!("s{}", values.get(row).unwrap_or("")),
        }
    }
}

/// Each row's SpliceAI score: the largest of its four delta scores, NaN without one. A table
/// without these columns has no scores.
fn spliceai_scores(table: &Table) -> Result<Vec<f64>> {
    let mut best = vec![f64::NAN; table.rows];
    for column in ["AG", "AL", "DG", "DL"] {
        let name = format!("SpliceAI_pred_DS_{column}");
        if !table.has_column(&name) {
            return Ok(vec![f64::NAN; table.rows]);
        }
        for (best, value) in best.iter_mut().zip(numeric_column(table, &name)?) {
            if !value.is_nan() && !(*best >= value) {
                *best = value;
            }
        }
    }
    Ok(best)
}

/// A variant of each row, as a key: its `chrom`, `pos`, `ref` and `alt`.
struct VariantKeys(Vec<KeyValues>);

impl VariantKeys {
    fn of(table: &Table) -> Result<Self> {
        ["chrom", "pos", "ref", "alt"]
            .iter()
            .map(|name| KeyValues::of(table, name))
            .collect::<Result<_>>()
            .map(Self)
    }

    fn key(&self, row: usize) -> Vec<String> {
        self.0.iter().map(|key| key.group_key(row)).collect()
    }
}

/// The rows the ranking keeps: every protein-coding row, as Python does, and two Rust-only
/// additions (24 November 2026). Both add rows only where the protein-coding rows say
/// nothing, so a variant gains no second set of rows for a gene it already has.
///
/// - `keep_listed_genes`: the rows of a gene of the family's gene ranking (its HPO terms) or
///   of the user's relevant gene list, whatever the transcript's biotype, when the variant
///   overlaps the transcript and has no protein-coding row of that gene. A row whose only
///   consequence is upstream or downstream of the transcript (VEP gives those up to 5 kb
///   away) does not overlap it. Genes match by symbol, as the bonus and the list weights
///   do. Without either list nothing is added.
/// - `spliceai_keep_score`: the rows with at least this SpliceAI score, of a variant that
///   has no protein-coding row.
fn kept_rows(table: &Table, options: &RankingOptions<'_>) -> Result<Vec<bool>> {
    let biotype = table.column("BIOTYPE")?;
    let coding: Vec<bool> = (0..table.rows)
        .map(|row| biotype.cell(row) == Cell::Text("protein_coding"))
        .collect();
    let mut keep = coding.clone();
    if options.keep_listed_genes {
        let mut genes: HashSet<String> = HashSet::new();
        if let Some(path) = options.phenotype_gene_ranking {
            genes.extend(phenotype_gene_indices(path)?.into_keys());
        }
        if let Some(path) = options.relevant_gene_list {
            genes.extend(gene_list(path)?);
        }
        if !genes.is_empty() {
            let consequence = table.column("Consequence")?;
            let outside = |row: usize| {
                matches!(
                    consequence.cell(row),
                    Cell::Text("upstream_gene_variant" | "downstream_gene_variant")
                )
            };
            let listed = symbol_mask(table, &genes)?;
            let overlapping: Vec<usize> = (0..table.rows)
                .filter(|row| listed[*row] && !coding[*row] && !outside(*row))
                .collect();
            let mut added = 0;
            if !overlapping.is_empty() {
                let (variants, symbols) = (VariantKeys::of(table)?, table.column("SYMBOL")?);
                let variant_and_gene = |row: usize| -> Vec<String> {
                    let mut key = variants.key(row);
                    if let Cell::Text(symbol) = symbols.cell(row) {
                        key.push(symbol.to_string());
                    }
                    key
                };
                // Only the listed genes' rows are looked up among the coding rows.
                let mut has_coding: HashMap<Vec<String>, bool> = overlapping
                    .iter()
                    .map(|row| (variant_and_gene(*row), false))
                    .collect();
                for row in (0..table.rows).filter(|row| coding[*row] && listed[*row]) {
                    if let Some(found) = has_coding.get_mut(&variant_and_gene(row)) {
                        *found = true;
                    }
                }
                for row in overlapping {
                    if !has_coding[&variant_and_gene(row)] {
                        keep[row] = true;
                        added += 1;
                    }
                }
            }
            eprintln!(
                "INFO: {added} rows of {} listed genes kept for variants without a protein-coding row of the gene",
                genes.len()
            );
        }
    }
    if let Some(cutoff) = options.spliceai_keep_score {
        let scores = spliceai_scores(table)?;
        let high: Vec<usize> = (0..table.rows)
            .filter(|row| !coding[*row] && scores[*row] >= cutoff)
            .collect();
        if !high.is_empty() {
            let variants = VariantKeys::of(table)?;
            // Only the few variants with such a row are looked up among the coding rows.
            let mut has_coding: HashMap<Vec<String>, bool> =
                high.iter().map(|row| (variants.key(*row), false)).collect();
            for row in (0..table.rows).filter(|row| coding[*row]) {
                if let Some(found) = has_coding.get_mut(&variants.key(row)) {
                    *found = true;
                }
            }
            let mut added = 0;
            for row in high {
                if !keep[row] && !has_coding[&variants.key(row)] {
                    keep[row] = true;
                    added += 1;
                }
            }
            eprintln!(
                "INFO: {added} rows with a SpliceAI score of at least {cutoff} kept for variants without a protein-coding row"
            );
        }
    }
    Ok(keep)
}

/// `sort_and_rank_variants(df, ped_df, fam_name, ...)`: the rows `kept_rows` keeps (Python's
/// are the protein-coding rows) with the ranking columns, sorted by the best `sort_index` of
/// each variant, then `ACMG_quant_score`, `chrom`, `pos`, `ref`, `alt`, and ranked by
/// `variant_id`.
///
/// The table is taken, and its rows are kept and sorted column by column, so it
/// is held about once (the step's peak memory was here: three copies).
pub(crate) fn sort_and_rank_variants(
    table: Table,
    ped: Option<&Table>,
    fam_name: Option<&str>,
    options: &RankingOptions<'_>,
) -> Result<Table> {
    // R1: drop the stale columns, copy the score, keep the ranking's rows.
    let mut df = table;
    for stale in ["haplo_insuf_index", "zygosity_inheritance_mechanism_index"] {
        if let Some(index) = df.column_index(stale) {
            df.names.remove(index);
            df.columns.remove(index);
        }
    }
    let quant = df.column("ACMG_quant_score")?.clone();
    df.set_column("sort_index", quant);
    let keep = kept_rows(&df, options)?;
    let mut df = df.into_filtered(&keep);
    let rows = df.rows;

    // R2: the proband from the pedigree.
    let (proband, male, female, proband_het) = pedigree_context(&df, ped, fam_name)?;

    // Control-common penalty (R3 multiplies it in).
    compute_control_common(&mut df, proband.as_deref())?;
    let mut sort_index = float_values(&df, "sort_index")?;
    let control_index = float_values(&df, "control_common_index")?;
    for row in 0..rows {
        sort_index[row] *= control_index[row];
    }
    df.set_column("sort_index", Column::Float(sort_index.clone()));

    // R4: the zygosity, inheritance and mechanism compatibility penalty.
    let masks = variant_mechanism_masks(&df)?;
    let x_locus: Vec<bool> = match df.column("chrom") {
        Ok(chrom) => (0..rows)
            .map(|row| {
                let text = python_str(chrom.cell(row))
                    .unwrap_or_default()
                    .to_uppercase();
                text.strip_prefix("CHR").unwrap_or(&text) == "X"
            })
            .collect(),
        Err(_) => vec![false; rows],
    };
    let mask = |name: &str| masks.get(name);
    let (dom_lof, dom_unresolved) = (
        mask("has_dom_lof_history")?,
        mask("has_dom_unresolved_history")?,
    );
    let (xr_lof, xr_unresolved) = (
        mask("has_x_recessive_lof_history")?,
        mask("has_x_recessive_unresolved_history")?,
    );
    let (xd_lof, xd_unresolved) = (
        mask("has_x_dominant_lof_history")?,
        mask("has_x_dominant_unresolved_history")?,
    );
    let exact_gof = mask("is_exact_gof")?;
    let predicted_lof = mask("is_predicted_lof")?;
    let quant_values = float_values(&df, "ACMG_quant_score")?;
    let mut penalty_rows = vec![false; rows];
    let mut high_acmg_lof = vec![false; rows];
    for row in 0..rows {
        let compatible = dom_lof[row]
            || dom_unresolved[row]
            || (x_locus[row]
                && male
                && (xr_lof[row] || xr_unresolved[row] || xd_lof[row] || xd_unresolved[row]))
            || (x_locus[row] && female && (xd_lof[row] || xd_unresolved[row]));
        let incompatible = predicted_lof[row] && !exact_gof[row] && !compatible;
        high_acmg_lof[row] = quant_values[row] > 0.89;
        penalty_rows[row] = high_acmg_lof[row] && proband_het[row] && incompatible;
    }
    let compatibility: Vec<f64> = penalty_rows
        .iter()
        .map(|penalty| if *penalty { 0.9 } else { 1.0 })
        .collect();
    df.set_column(
        "zygosity_inheritance_mechanism_compatibility",
        Column::Float(compatibility.clone()),
    );

    // R5 and R6: the dispensable and relevant gene lists. The phenotype gene bonus (a
    // Rust-only revision) takes R6's place: a row gains the larger of its relevant gene
    // list weight and its bonus, never both.
    let phenotype_index = match options.phenotype_gene_ranking {
        Some(path) => {
            let indices = phenotype_gene_indices(path)?;
            eprintln!(
                "INFO: phenotype gene bonus from {path}: {} genes of the gene ranking",
                indices.len()
            );
            Some(symbol_values(&df, &indices)?)
        }
        None => None,
    };
    for (list, column, weight) in [
        (options.dispensable_gene_list, "dispensable_gene_index", 0.9),
        (options.relevant_gene_list, "relevant_gene_index", 1.2),
    ] {
        let Some(path) = list else {
            continue;
        };
        let genes = gene_list(path)?;
        let listed = symbol_mask(&df, &genes)?;
        let index: Vec<f64> = listed
            .iter()
            .map(|listed| if *listed { weight } else { 1.0 })
            .collect();
        df.set_column(column, Column::Float(index.clone()));
        let gain: Vec<f64> = match &phenotype_index {
            Some(bonus) if column == "relevant_gene_index" => index
                .iter()
                .zip(bonus)
                .map(|(listed, bonus)| listed.max(*bonus))
                .collect(),
            _ => index,
        };
        for row in 0..rows {
            sort_index[row] *= gain[row];
        }
        df.set_column("sort_index", Column::Float(sort_index.clone()));
    }
    if let Some(bonus) = phenotype_index {
        if options.relevant_gene_list.is_none() {
            for row in 0..rows {
                sort_index[row] *= bonus[row];
            }
            df.set_column("sort_index", Column::Float(sort_index.clone()));
        }
        df.set_column("phenotype_gene_index", Column::Float(bonus));
    }
    for row in 0..rows {
        if penalty_rows[row] {
            sort_index[row] *= compatibility[row];
        }
    }
    df.set_column("sort_index", Column::Float(sort_index.clone()));

    // R7: the pext modulation of high-scoring rows.
    let mut pext_column = None;
    let tissues: Vec<&str> = options
        .pext_tissues
        .split(',')
        .map(python_strip)
        .filter(|tissue| !tissue.is_empty())
        .collect();
    if let Some(first) = tissues.first() {
        let name = format!("PEXT_{}", first.replace('-', "_"));
        if df.has_column(&name) {
            pext_column = Some(name);
        }
    }
    if pext_column.is_none() && df.has_column("PEXT_MEAN") {
        pext_column = Some("PEXT_MEAN".to_owned());
    }
    if let Some(name) = pext_column {
        let cutoff = options.pext_low_expression_cutoff;
        let floor = options.pext_penalty_floor;
        let shape = options.pext_penalty_shape;
        if cutoff <= 0.0 {
            return Err(python_error(
                "ValueError",
                "pext_low_expression_cutoff must be > 0",
            ));
        }
        if !(0.0 <= floor && floor <= 1.0) {
            return Err(python_error(
                "ValueError",
                "pext_penalty_floor must be between 0 and 1",
            ));
        }
        if shape <= 0.0 {
            return Err(python_error("ValueError", "pext_penalty_shape must be > 0"));
        }
        let pext: Vec<f64> = pext_values(&df, &name)?
            .into_iter()
            .map(|value| clip(if value.is_nan() { 1.0 } else { value }, 0.0, 1.0))
            .collect();
        let pext_index: Vec<f64> = pext
            .iter()
            .map(|value| {
                let scaled = clip(value / cutoff, 0.0, 1.0);
                let modulated = floor + (1.0 - floor) * series_power(scaled, shape);
                if *value >= cutoff { 1.0 } else { modulated }
            })
            .collect();
        // On no rows, pandas' `map` of an object column gives an empty object Series, and
        // the arithmetic and `np.where` after it keep that type.
        let numeric = matches!(
            df.column(&name)?,
            Column::Float(_) | Column::Int(_) | Column::Bool(_)
        );
        let pext_column = if rows == 0 && !numeric {
            Column::Object(Vec::new())
        } else {
            Column::Float(pext_index.clone())
        };
        df.set_column("pext_sort_index", pext_column);
        for row in 0..rows {
            if high_acmg_lof[row] {
                sort_index[row] *= pext_index[row];
            }
        }
        df.set_column("sort_index", Column::Float(sort_index.clone()));
    }

    // R8: the best sort_index of each variant (groups by chrom, pos, ref, alt; a missing
    // key leaves NaN), the sort, and the rank by first appearance of variant_id.
    let keys: Vec<KeyValues> = ["chrom", "pos", "ref", "alt"]
        .iter()
        .map(|name| KeyValues::of(&df, name))
        .collect::<Result<_>>()?;
    let mut best: HashMap<Vec<String>, (usize, f64)> = HashMap::new();
    let mut group_of: Vec<Option<Vec<String>>> = Vec::with_capacity(rows);
    for row in 0..rows {
        if keys.iter().any(|key| key.missing(row)) {
            group_of.push(None);
            continue;
        }
        let key: Vec<String> = keys.iter().map(|key| key.group_key(row)).collect();
        // pandas' group max: the first value that is greater than every earlier one.
        let entry = best.entry(key.clone()).or_insert((0, f64::NEG_INFINITY));
        let value = sort_index[row];
        if !value.is_nan() {
            entry.0 += 1;
            if value > entry.1 {
                entry.1 = value;
            }
        }
        group_of.push(Some(key));
    }
    let max_score: Vec<f64> = group_of
        .iter()
        .map(|key| match key {
            Some(key) => match best[key] {
                (0, _) => f64::NAN,
                (_, value) => value,
            },
            None => f64::NAN,
        })
        .collect();
    let sort_keys = [
        (KeyValues::Numbers(max_score), false),
        (KeyValues::Numbers(quant_values), false),
    ];
    let mut order: Vec<usize> = (0..rows).collect();
    order.sort_by(|left, right| {
        let all = sort_keys
            .iter()
            .map(|(key, ascending)| (key, *ascending))
            .chain(keys.iter().map(|key| (key, true)));
        for (key, ascending) in all {
            let ordering = match (key.missing(*left), key.missing(*right)) {
                (true, true) => Ordering::Equal,
                (true, false) => Ordering::Greater,
                (false, true) => Ordering::Less,
                (false, false) => {
                    let ordering = key.compare(*left, *right);
                    if ascending {
                        ordering
                    } else {
                        ordering.reverse()
                    }
                }
            };
            if ordering != Ordering::Equal {
                return ordering;
            }
        }
        Ordering::Equal
    });
    if let Some(index) = df.column_index("max_variant_score") {
        df.names.remove(index);
        df.columns.remove(index);
    }
    let mut sorted = df.into_rows(&order);
    let identities: TextColumn = match sorted.column("variant_id") {
        Ok(Column::Text(values)) => values.clone(),
        Ok(other) => {
            return Err(refuse(format!(
                "variant_id is a {} column",
                other.dtype_name()
            )));
        }
        Err(_) => {
            let parts: Vec<&Column> = ["chrom", "pos", "ref", "alt"]
                .iter()
                .map(|name| sorted.column(name))
                .collect::<Result<_>>()?;
            (0..sorted.rows)
                .map(|row| {
                    Some(
                        parts
                            .iter()
                            .map(|column| {
                                python_str(column.cell(row)).unwrap_or_else(|| "nan".to_owned())
                            })
                            .collect::<Vec<_>>()
                            .join(":"),
                    )
                })
                .collect()
        }
    };
    let mut codes: HashMap<&str, i64> = HashMap::new();
    let ranks: Vec<i64> = identities
        .iter()
        .map(|identity| match identity {
            None => 0,
            Some(identity) => {
                let next = codes.len() as i64 + 1;
                *codes.entry(identity).or_insert(next)
            }
        })
        .collect();
    sorted.set_column("variant_rank", Column::Int(ranks));
    Ok(sorted)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table(columns: Vec<(&str, Column)>) -> Table {
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
    fn summary_without_clingen_keeps_int_columns() {
        let mut input = table(vec![(
            "variant_id",
            Column::Text(vec![Some("1:1:A-T".into()), Some("1:2:C-G".into())].into()),
        )]);
        let criteria = vec![
            ("PVS1", Column::Int(vec![4, 0])),
            ("PM2", Column::Int(vec![1, 0])),
            ("BP7", Column::Bool(vec![false, false])),
        ];
        let matrix = summarize_acmg_criteria(&mut input, &criteria, None, false).unwrap();
        assert_eq!(matrix.names.len(), 29);
        assert_eq!(matrix.column("PS4").unwrap(), &Column::Int(vec![0, 0]));
        assert_eq!(
            input.column("ACMG_criteria").unwrap(),
            &Column::Text(vec![Some("PVS1;PM2_Supporting".into()), None].into())
        );
        let bad = vec![("PS1", Column::Int(vec![6, 0]))];
        let error = summarize_acmg_criteria(&mut input, &bad, None, false)
            .unwrap_err()
            .to_string();
        assert_eq!(error, "ValueError: Column PS1 has non-integer values");
        let error = summarize_acmg_criteria(&mut input, &criteria, Some("/no/such.pkl"), true)
            .unwrap_err()
            .to_string();
        assert!(error.starts_with(
            "FileNotFoundError: ClinGen override is enabled but the evidence map is missing: '/no/such.pkl'."
        ));
    }

    #[test]
    fn native_clingen_keeps_overrides_missing_codes_and_nonempty_map_dtypes() {
        use crate::cache_convert::clingen_map::{ClingenMap, SCHEMA};
        let path =
            std::env::temp_dir().join(format!("priva_scoring_clingen_{}.rkyv", std::process::id()));
        let root = ClingenMap {
            provenance: vec![],
            variants: vec!["a".into(), "b".into()],
            codes: vec![Some("PS3".into()), None],
        };
        crate::rust_cache::write_cache(
            &path,
            SCHEMA,
            &rkyv::to_bytes::<rkyv::rancor::Error>(&root).unwrap(),
        )
        .unwrap();
        let mut input = table(vec![(
            "variant_id",
            Column::Text(vec![Some("a".into()), Some("b".into()), Some("missing".into())].into()),
        )]);
        let criteria = vec![("PVS1", Column::Int(vec![4, 4, 0]))];
        let matrix = summarize_acmg_criteria(&mut input, &criteria, path.to_str(), true).unwrap();
        assert_eq!(
            input.column("ACMG_criteria").unwrap(),
            &Column::Text(vec![Some("PS3".into()), Some("PVS1".into()), None].into())
        );
        assert_eq!(matrix.column("PVS1").unwrap(), &Column::Int(vec![0, 4, 0]));
        assert_eq!(matrix.column("PS3").unwrap(), &Column::Int(vec![3, 0, 0]));

        // A present NaN/no-code entry still makes the map nonempty. The original empty
        // override assignment changes integer criteria to object when a summary is text.
        let root = ClingenMap {
            provenance: vec![],
            variants: vec!["b".into()],
            codes: vec![None],
        };
        crate::rust_cache::write_cache(
            &path,
            SCHEMA,
            &rkyv::to_bytes::<rkyv::rancor::Error>(&root).unwrap(),
        )
        .unwrap();
        let mut input = table(vec![(
            "variant_id",
            Column::Text(vec![Some("b".into())].into()),
        )]);
        let matrix = summarize_acmg_criteria(
            &mut input,
            &[("PVS1", Column::Int(vec![4]))],
            path.to_str(),
            true,
        )
        .unwrap();
        assert!(
            matches!(matrix.column("PVS1").unwrap(), Column::Object(values) if values == &[ObjectValue::Int(4)])
        );
    }

    #[test]
    fn posterior_columns_follow_the_matrix() {
        let matrix = table(vec![
            ("PVS1", Column::Int(vec![4, 0])),
            ("BA1", Column::Int(vec![0, 5])),
            (
                "variant_id",
                Column::Text(vec![Some("a".into()), Some("b".into())].into()),
            ),
        ]);
        let (scores, classes) = posterior_columns(&matrix).unwrap();
        let Column::Float(scores) = scores else {
            panic!("float scores")
        };
        assert_eq!(scores[1], 0.0);
        assert_eq!(
            classes,
            Column::Text(vec![Some("Likely Pathogenic".into()), Some("Benign".into())].into())
        );
    }

    #[test]
    fn control_common_without_control_columns() {
        let mut input = table(vec![("x", Column::Int(vec![1]))]);
        compute_control_common(&mut input, None).unwrap();
        assert_eq!(
            input.column("control_common").unwrap(),
            &Column::Bool(vec![false])
        );
        assert_eq!(
            input.column("control_common_index").unwrap(),
            &Column::Float(vec![1.0])
        );
    }

    fn texts(values: &[&str]) -> Column {
        Column::Text(
            values
                .iter()
                .map(|value| Some((*value).to_owned()))
                .collect(),
        )
    }

    #[test]
    fn ranking_sorts_by_the_best_score_of_each_variant() {
        let input = table(vec![
            ("chrom", texts(&["chr1", "chr1", "chr2", "chrX"])),
            ("pos", Column::Int(vec![5, 5, 3, 9])),
            ("ref", texts(&["A", "A", "C", "G"])),
            ("alt", texts(&["T", "T", "G", "A"])),
            (
                "BIOTYPE",
                texts(&[
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                    "lncRNA",
                ]),
            ),
            ("SYMBOL", texts(&["G1", "G1", "G2", "G3"])),
            (
                "ACMG_quant_score",
                Column::Float(vec![0.2, 0.95, 0.5, 0.99]),
            ),
            (
                "var_plausible_patho_mechs",
                texts(&["", "", "dominant_LOF", ""]),
            ),
            (
                "variant_effect",
                texts(&["uncertain", "uncertain", "uncertain", "uncertain"]),
            ),
            (
                "variant_id",
                texts(&["1:5:A-T", "1:5:A-T", "2:3:C-G", "X:9:G-A"]),
            ),
        ]);
        let options = RankingOptions {
            pext_tissues: "",
            pext_low_expression_cutoff: 0.1,
            pext_penalty_floor: 0.8,
            pext_penalty_shape: 0.5,
            relevant_gene_list: None,
            dispensable_gene_list: None,
            phenotype_gene_ranking: None,
            keep_listed_genes: false,
            spliceai_keep_score: None,
        };
        let ranked = sort_and_rank_variants(input.clone(), None, None, &options).unwrap();
        assert_eq!(ranked.rows, 3);
        assert_eq!(
            ranked.column("ACMG_quant_score").unwrap(),
            &Column::Float(vec![0.95, 0.2, 0.5])
        );
        assert_eq!(
            ranked.column("variant_rank").unwrap(),
            &Column::Int(vec![1, 1, 2])
        );
        assert_eq!(
            ranked.column("control_common_index").unwrap(),
            &Column::Float(vec![1.0; 3])
        );
        assert!(!ranked.has_column("max_variant_score"));
        let error = sort_and_rank_variants(input.clone(), Some(&input), Some("F"), &options)
            .unwrap_err()
            .to_string();
        assert_eq!(error, "KeyError: '#FamilyID'");
    }

    #[test]
    fn phenotype_gene_index_follows_the_share_of_the_top_score() {
        let close = |value: f64, expected: f64| (value - expected).abs() < 1e-12;
        assert_eq!(phenotype_gene_index(1.0), 1.2);
        assert!(close(phenotype_gene_index(0.1), 1.1));
        assert!(close(phenotype_gene_index(10f64.powf(-0.5)), 1.15));
        assert!(close(phenotype_gene_index(10f64.powf(-1.5)), 1.05));
        assert_eq!(phenotype_gene_index(0.01), 1.0);
        assert_eq!(phenotype_gene_index(0.0001), 1.0);
        assert_eq!(phenotype_gene_index(0.0), 1.0);
        assert_eq!(phenotype_gene_index(f64::NAN), 1.0);
    }

    fn floats(table: &Table, name: &str) -> Vec<f64> {
        match table.column(name).unwrap() {
            Column::Float(values) => values.clone(),
            _ => panic!("{name} is not a float column"),
        }
    }

    #[test]
    fn phenotype_gene_bonus_raises_the_sort_index_of_ranked_genes() {
        // Two variants whose ACMG scores are close: G2 is the family's top gene, G1 has a
        // hundredth of its score.
        let input = table(vec![
            ("chrom", texts(&["chr1", "chr2"])),
            ("pos", Column::Int(vec![5, 3])),
            ("ref", texts(&["A", "C"])),
            ("alt", texts(&["T", "G"])),
            ("BIOTYPE", texts(&["protein_coding", "protein_coding"])),
            ("SYMBOL", texts(&["G1", "G2"])),
            ("ACMG_quant_score", Column::Float(vec![0.5, 0.45])),
            ("var_plausible_patho_mechs", texts(&["", ""])),
            ("variant_effect", texts(&["uncertain", "uncertain"])),
            ("variant_id", texts(&["1:5:A-T", "2:3:C-G"])),
        ]);
        let directory =
            std::env::temp_dir().join(format!("priva_phenotype_bonus_{}", std::process::id()));
        std::fs::create_dir_all(&directory).unwrap();
        let ranking = directory.join("F1.gene_ranked.tsv");
        std::fs::write(
            &ranking,
            "sample_id\tgene_symbol\tncbi_gene_ids\ttotal_score\tbase_term_score\n\
             F1\tG2\t2\t5.0\t1.0\n\
             F1\tG7\t7\t2.5\t1.0\n\
             F1\tG2\t2\t0.5\t1.0\n\
             F1\tG1\t1\t0.05\t1.0\n",
        )
        .unwrap();
        let relevant = directory.join("relevant_genes.txt");
        std::fs::write(&relevant, "G1\n").unwrap();
        let broken = directory.join("broken.gene_ranked.tsv");
        std::fs::write(
            &broken,
            "sample_id\tgene_symbol\ttotal_score\nF1\tG2\tabc\n",
        )
        .unwrap();
        let (ranking, relevant, broken) = (
            ranking.to_str().unwrap(),
            relevant.to_str().unwrap(),
            broken.to_str().unwrap(),
        );
        let options = |relevant_gene_list, phenotype_gene_ranking| RankingOptions {
            pext_tissues: "",
            pext_low_expression_cutoff: 0.1,
            pext_penalty_floor: 0.8,
            pext_penalty_shape: 0.5,
            relevant_gene_list,
            dispensable_gene_list: None,
            phenotype_gene_ranking,
            keep_listed_genes: false,
            spliceai_keep_score: None,
        };
        let rank = |relevant_gene_list, phenotype_gene_ranking| {
            sort_and_rank_variants(
                input.clone(),
                None,
                None,
                &options(relevant_gene_list, phenotype_gene_ranking),
            )
        };

        // Without the ranking: G1 first, and no new column.
        let plain = rank(None, None).unwrap();
        assert_eq!(
            plain.column("variant_id").unwrap(),
            &texts(&["1:5:A-T", "2:3:C-G"])
        );
        assert!(!plain.has_column("phenotype_gene_index"));
        let plain_index = floats(&plain, "sort_index");

        // With it: G2 gains 1.2 (the larger of its two rows), G1 nothing, so G2 goes first.
        let bonus = rank(None, Some(ranking)).unwrap();
        assert_eq!(
            bonus.column("variant_id").unwrap(),
            &texts(&["2:3:C-G", "1:5:A-T"])
        );
        assert_eq!(floats(&bonus, "phenotype_gene_index"), vec![1.2, 1.0]);
        let bonus_index = floats(&bonus, "sort_index");
        assert!((bonus_index[0] - plain_index[1] * 1.2).abs() < 1e-12);
        assert!((bonus_index[1] - plain_index[0]).abs() < 1e-12);

        // With a relevant gene list too: each row gains the larger weight, never both.
        let both = rank(Some(relevant), Some(ranking)).unwrap();
        assert_eq!(
            both.column("variant_id").unwrap(),
            &texts(&["1:5:A-T", "2:3:C-G"])
        );
        assert_eq!(floats(&both, "relevant_gene_index"), vec![1.2, 1.0]);
        assert_eq!(floats(&both, "phenotype_gene_index"), vec![1.0, 1.2]);
        let both_index = floats(&both, "sort_index");
        assert!((both_index[0] - plain_index[0] * 1.2).abs() < 1e-12);
        assert!((both_index[1] - plain_index[1] * 1.2).abs() < 1e-12);

        // The relevant gene list alone keeps Python's ranking.
        let listed = rank(Some(relevant), None).unwrap();
        assert!(!listed.has_column("phenotype_gene_index"));
        assert_eq!(floats(&listed, "relevant_gene_index"), vec![1.2, 1.0]);
        let listed_index = floats(&listed, "sort_index");
        assert_eq!(listed_index[0], both_index[0]);
        assert_eq!(listed_index[1], plain_index[1]);

        let error = rank(None, Some(broken)).unwrap_err().to_string();
        assert!(
            error.contains("total_score of G2 is not a number"),
            "{error}"
        );
        std::fs::remove_dir_all(&directory).unwrap();
    }

    #[test]
    fn listed_genes_and_high_spliceai_scores_keep_rows_beside_the_coding_ones() {
        // One row per transcript: variant, gene, biotype, consequence, SpliceAI score.
        let rows = [
            (
                "1:5:A:T",
                "SIRT4",
                "protein_coding",
                "upstream_gene_variant",
                0.9,
            ),
            (
                "1:5:A:T",
                "RNU4-2",
                "snRNA",
                "non_coding_transcript_exon_variant",
                0.0,
            ),
            ("1:5:A:T", "RNU4-1", "snRNA", "downstream_gene_variant", 0.0),
            (
                "2:3:C:G",
                "LNC1",
                "lncRNA",
                "intron_variant&non_coding_transcript_variant",
                0.3,
            ),
            (
                "2:3:C:G",
                "LNC2",
                "lncRNA",
                "intron_variant&non_coding_transcript_variant",
                f64::NAN,
            ),
            (
                "3:9:G:A",
                "G9",
                "retained_intron",
                "non_coding_transcript_exon_variant",
                0.29,
            ),
            ("4:1:T:C", "G4", "protein_coding", "intron_variant", 0.9),
            ("4:1:T:C", "G4", "retained_intron", "intron_variant", 0.9),
        ];
        let part = |index: usize| -> Vec<&str> {
            rows.iter()
                .map(|row| row.0.split(':').nth(index).unwrap())
                .collect()
        };
        let input = table(vec![
            ("chrom", texts(&part(0))),
            (
                "pos",
                Column::Int(part(1).iter().map(|pos| pos.parse().unwrap()).collect()),
            ),
            ("ref", texts(&part(2))),
            ("alt", texts(&part(3))),
            ("SYMBOL", texts(&rows.map(|row| row.1))),
            ("BIOTYPE", texts(&rows.map(|row| row.2))),
            ("Consequence", texts(&rows.map(|row| row.3))),
            // The score is the largest of the four columns, whichever holds it.
            ("SpliceAI_pred_DS_AG", Column::Float(vec![0.0; rows.len()])),
            (
                "SpliceAI_pred_DS_AL",
                Column::Float(rows.map(|row| row.4).to_vec()),
            ),
            (
                "SpliceAI_pred_DS_DG",
                Column::Float(vec![f64::NAN; rows.len()]),
            ),
            ("SpliceAI_pred_DS_DL", Column::Float(vec![0.1; rows.len()])),
            ("ACMG_quant_score", Column::Float(vec![0.5; rows.len()])),
            ("var_plausible_patho_mechs", texts(&[""; 8])),
            ("variant_effect", texts(&["uncertain"; 8])),
            ("variant_id", texts(&rows.map(|row| row.0))),
        ]);
        let directory =
            std::env::temp_dir().join(format!("priva_kept_rows_{}", std::process::id()));
        std::fs::create_dir_all(&directory).unwrap();
        let ranking = directory.join("F1.gene_ranked.tsv");
        std::fs::write(
            &ranking,
            "sample_id\tgene_symbol\ttotal_score\nF1\tRNU4-2\t5.0\nF1\tRNU4-1\t4.0\n",
        )
        .unwrap();
        let relevant = directory.join("relevant_genes.txt");
        std::fs::write(&relevant, "LNC2\nG4\n").unwrap();
        let (ranking, relevant) = (ranking.to_str().unwrap(), relevant.to_str().unwrap());
        // The kept rows' gene and biotype, in the input's order.
        let kept = |keep_listed_genes, relevant_gene_list, phenotype_gene_ranking, score| {
            let options = RankingOptions {
                pext_tissues: "",
                pext_low_expression_cutoff: 0.1,
                pext_penalty_floor: 0.8,
                pext_penalty_shape: 0.5,
                relevant_gene_list,
                dispensable_gene_list: None,
                phenotype_gene_ranking,
                keep_listed_genes,
                spliceai_keep_score: score,
            };
            let keep = kept_rows(&input, &options).unwrap();
            let ranked = sort_and_rank_variants(input.clone(), None, None, &options).unwrap();
            assert_eq!(ranked.rows, keep.iter().filter(|kept| **kept).count());
            rows.iter()
                .zip(keep)
                .filter(|(_, kept)| *kept)
                .map(|(row, _)| format!("{} {}", row.1, row.2))
                .collect::<Vec<_>>()
        };
        let coding = vec!["SIRT4 protein_coding", "G4 protein_coding"];

        // Python's rows: no lists and no SpliceAI score, or the rules off.
        assert_eq!(kept(true, None, None, None), coding);
        assert_eq!(kept(false, Some(relevant), Some(ranking), None), coding);

        // The gene ranking keeps the listed genes' overlapping rows: RNU4-2's exon row stays
        // beside SIRT4's, RNU4-1's downstream row does not.
        assert_eq!(
            kept(true, None, Some(ranking), None),
            vec!["SIRT4 protein_coding", "RNU4-2 snRNA", "G4 protein_coding"]
        );
        // Without a ranking the user's list does; with both, a gene in either counts. G4 is
        // listed too, but its variant has G4's protein-coding row, so its other row stays out.
        assert_eq!(
            kept(true, Some(relevant), None, None),
            vec!["SIRT4 protein_coding", "LNC2 lncRNA", "G4 protein_coding"]
        );
        assert_eq!(
            kept(true, Some(relevant), Some(ranking), None),
            vec![
                "SIRT4 protein_coding",
                "RNU4-2 snRNA",
                "LNC2 lncRNA",
                "G4 protein_coding"
            ]
        );

        // A score of 0.3 keeps the scored row of a variant without a coding row: 0.29 does
        // not, and the retained intron row of a variant with a coding row is not added.
        assert_eq!(
            kept(false, None, None, Some(SPLICEAI2_KEEP_SCORE)),
            vec!["SIRT4 protein_coding", "LNC1 lncRNA", "G4 protein_coding"]
        );

        // The kept row of a listed gene takes the gene's bonus under its own name (the rows
        // of one variant keep the table's order).
        let options = RankingOptions {
            pext_tissues: "",
            pext_low_expression_cutoff: 0.1,
            pext_penalty_floor: 0.8,
            pext_penalty_shape: 0.5,
            relevant_gene_list: None,
            dispensable_gene_list: None,
            phenotype_gene_ranking: Some(ranking),
            keep_listed_genes: true,
            spliceai_keep_score: Some(SPLICEAI2_KEEP_SCORE),
        };
        let ranked = sort_and_rank_variants(input.clone(), None, None, &options).unwrap();
        assert_eq!(
            ranked.column("SYMBOL").unwrap(),
            &texts(&["SIRT4", "RNU4-2", "LNC1", "G4"])
        );
        assert_eq!(
            floats(&ranked, "phenotype_gene_index"),
            vec![1.0, 1.2, 1.0, 1.0]
        );
        assert_eq!(
            ranked.column("variant_rank").unwrap(),
            &Column::Int(vec![1, 1, 2, 3])
        );
        std::fs::remove_dir_all(&directory).unwrap();
    }
}
