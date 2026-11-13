//! `scripts/acmg_pp1_bs4_bp2_pm3_family.py` (602 lines): PP1, BS4, BP2 and PM3.
//!
//! - `PP1_criteria`: when both `--pp1_vcf` and `--pp1_ped` are given (accepted
//!   and forwarded since main c529805), the co-segregation counts of
//!   `find_cosegregation_vars` keyed by `chrom:pos:ref-alt`, weighted by the
//!   variant's selected inheritance, give points: 1 to 1.9 supporting (1), 2
//!   to 3.9 moderate (2), 4 and more strong (3). The sums are made in the
//!   script's order (`count * 0.4` in floating point), so the cut-offs fall
//!   exactly where Python's do.
//! - `BS4_criteria`: genotype doses of affected and unaffected relatives of
//!   the family (count of `1` in the GT before the first `:`; a male's chrX or
//!   chrY carrier counts 2) against the inheritance gates; 3 wins over 1.
//! - `BP2_PM3_criteria`: `determine_phase` with the pathogenic partners
//!   `vep_consq_lof | splicing_lof | 5UTR_lof | (PS1 & PS3) | PVS1` (PS1 & PS3
//!   is numpy's bitwise and of the strengths, as in Python), then BP2 (1) and
//!   PM3 (2) from the cis and trans calls. The table gains `is_pathogenic` and
//!   the patients' `<sample>_phase_patho` columns, as in Python.
//!
//! The condition masks come from `acmg_variant_mechanism::variant_condition_masks`
//! (thread T3's adapter). Every function also has a `_with` form taking the
//! masks from elsewhere: the harness runs that form when the call names masks
//! the Python capture computed, and the plain form (T3's adapter) otherwise.
//!
//! Refused (exit 1): a `chrom` column that is not text without missing values
//! for PP1; LOF columns other than booleans for BP2/PM3 (pandas' `|` treats a
//! missing value on the left differently from one on the right); nullable
//! sample columns for BS4.
use super::acmg_ps2_pm6_denovo::{
    PedNumber, criteria_frame, family_rows, harness_ped, harness_text, ped_to_nullable_int,
    ped_to_numeric, python_str,
};
use super::acmg_variant_mechanism::variant_condition_masks;
use super::determine_phase::{determine_cis_trans_relationships, json_mask, written_columns};
use super::find_cosegregation_vars::find_cosegregating_variants;
use super::frame_file::read_frame;
use super::types::{AcmgError, HarnessCall, NamedMasks, StepOutput};
use crate::cli::{Result, invalid};
use crate::json_text::Json;
use crate::pandas_table::{Cell, Column, Table};
use std::collections::HashMap;
use std::path::Path;

/// Where the condition masks come from (`_variant_condition_masks(df)`).
pub(crate) type MaskSource<'a> = &'a dyn Fn(&Table) -> Result<NamedMasks>;

fn refused(message: impl Into<String>) -> Box<dyn std::error::Error> {
    invalid(format!("refused: {}", message.into()))
}

/// The masks the family criteria read, by their Python names.
struct Masks<'a> {
    all: &'a NamedMasks,
}

impl Masks<'_> {
    fn get(&self, name: &str) -> Result<&[bool]> {
        self.all.get(name)
    }

    /// `has_mendelian & ~has_non_monogenic & ~has_non_mendelian & ~has_incomplete_penetrance`.
    fn valid_model(&self) -> Result<Vec<bool>> {
        let (mendelian, non_monogenic, non_mendelian, incomplete) = (
            self.get("has_mendelian")?,
            self.get("has_non_monogenic")?,
            self.get("has_non_mendelian")?,
            self.get("has_incomplete_penetrance")?,
        );
        Ok((0..mendelian.len())
            .map(|row| {
                mendelian[row] && !non_monogenic[row] && !non_mendelian[row] && !incomplete[row]
            })
            .collect())
    }
}

/// The `chrom` text of every row, refusing what `.str.contains` would turn into NaN.
fn chrom_texts(table: &Table) -> Result<Vec<&str>> {
    let column = table.column("chrom")?;
    let Column::Text(values) = column else {
        return Err(refused(format!(
            "chrom of type {} (pandas' .str fails or gives NaN)",
            column.dtype_name()
        )));
    };
    values
        .iter()
        .map(|value| {
            value
                .as_deref()
                .ok_or_else(|| refused("a missing chrom value (str.contains gives NaN)"))
        })
        .collect()
}

/// `PP1_criteria(df, multi_fam_vcf, multi_fam_ped)` with the step's masks.
pub(crate) fn pp1_criteria(
    table: &Table,
    multi_fam_vcf: &str,
    multi_fam_ped: &str,
) -> Result<StepOutput> {
    pp1_criteria_with(
        table,
        multi_fam_vcf,
        multi_fam_ped,
        &variant_condition_masks,
    )
}

/// `PP1_criteria` with the condition masks taken from `masks`.
pub(crate) fn pp1_criteria_with(
    table: &Table,
    multi_fam_vcf: &str,
    multi_fam_ped: &str,
    masks: MaskSource<'_>,
) -> Result<StepOutput> {
    let mut pp1 = vec![0i64; table.rows];
    if !multi_fam_vcf.is_empty() && !multi_fam_ped.is_empty() {
        let found = find_cosegregating_variants(
            Path::new(multi_fam_vcf),
            Path::new(multi_fam_ped),
            "both",
        )?;
        let mut counts: HashMap<String, [i64; 8]> = HashMap::new();
        for (key, values) in &found {
            counts.insert(
                format!("{}:{}:{}-{}", key.chrom, key.pos, key.reference, key.alt),
                values.values(),
            );
        }
        let ids = table.column("variant_id")?;
        let per_row: Vec<[i64; 8]> = (0..table.rows)
            .map(|row| match ids.cell(row) {
                Cell::Text(id) => counts.get(id).copied().unwrap_or_default(),
                _ => [0; 8],
            })
            .collect();
        let all = masks(table)?;
        let masks = Masks { all: &all };
        let valid = masks.valid_model()?;
        let (recessive, x_recessive, x_dominant, dominant) = (
            masks.get("has_recessive")?,
            masks.get("has_x_linked_recessive")?,
            masks.get("has_x_linked_dominant")?,
            masks.get("has_dominant")?,
        );
        let chrom = chrom_texts(table)?;
        for row in 0..table.rows {
            let [
                affected,
                control,
                affected_hom,
                control_nonhom,
                male_affected,
                male_control,
                female_affected,
                female_control,
            ] = per_row[row];
            let encode = |flag: bool| i64::from(flag);
            let recessive_encode = encode(valid[row] && recessive[row]);
            let x_recessive_encode = encode(valid[row] && x_recessive[row]);
            let x_dominant_encode = encode(valid[row] && x_dominant[row] && !recessive[row]);
            let dominant_encode = encode(valid[row] && dominant[row] && !recessive[row]);
            let text = chrom[row];
            let autosomal =
                encode(!(text.contains('X') || text.contains('Y') || text.contains('M')));
            let x = encode(text.contains('X'));
            let recessive_points = (autosomal * recessive_encode * affected * 2) as f64
                + (autosomal * recessive_encode * control) as f64 * 0.4;
            let x_recessive_points = x * x_recessive_encode * (male_affected + male_control);
            let x_dominant_points = x
                * x_dominant_encode
                * (male_affected + male_control + female_affected + female_control);
            let dominant_points =
                autosomal * dominant_encode * (affected + control + affected_hom + control_nonhom);
            let points = recessive_points
                + x_recessive_points as f64
                + x_dominant_points as f64
                + dominant_points as f64;
            if (1.0..=1.9).contains(&points) {
                pp1[row] = 1;
            }
            if (2.0..=3.9).contains(&points) {
                pp1[row] = 2;
            }
            if points >= 4.0 {
                pp1[row] = 3;
            }
        }
    }
    Ok(StepOutput {
        criteria: vec![("PP1", Column::Int(pp1))],
        ..StepOutput::default()
    })
}

/// `BS4_criteria(df, ped_df, fam_name)` with the step's masks.
pub(crate) fn bs4_criteria(table: &Table, ped: &Table, fam_name: &str) -> Result<StepOutput> {
    bs4_criteria_with(table, ped, fam_name, &variant_condition_masks)
}

/// The text `df[col].fillna(".").astype(str)` gives a cell.
fn dotted_text(column: &Column, row: usize) -> Result<String> {
    match column {
        Column::NullableInt(_) | Column::NullableBool(_) => Err(refused(format!(
            "a nullable {} sample column (fillna(\".\") does not fit it)",
            column.dtype_name()
        ))),
        _ => Ok(python_str(column.cell(row)).unwrap_or_else(|| ".".to_owned())),
    }
}

fn bs4_output(rows: usize, values: Vec<i64>) -> StepOutput {
    debug_assert_eq!(values.len(), rows);
    StepOutput {
        criteria: vec![("BS4", Column::Int(values))],
        ..StepOutput::default()
    }
}

/// `BS4_criteria` with the condition masks taken from `masks`.
pub(crate) fn bs4_criteria_with(
    table: &Table,
    ped: &Table,
    fam_name: &str,
    masks: MaskSource<'_>,
) -> Result<StepOutput> {
    let rows = table.rows;
    let family = family_rows(ped, fam_name)?;
    if family.is_empty() {
        return Ok(bs4_output(rows, vec![0; rows]));
    }
    let phenotype_column = ped.column("Phenotype")?;
    let phenotypes: Vec<Option<i64>> = family
        .iter()
        .map(|row| ped_to_nullable_int(phenotype_column.cell(*row)))
        .collect::<Result<_>>()?;
    let ids_column = ped.column("IndividualID")?;
    let ids: Vec<String> = family
        .iter()
        .map(|row| python_str(ids_column.cell(*row)).unwrap_or_else(|| "nan".to_owned()))
        .collect();
    let sexes: Vec<PedNumber> = match ped.column_index("Sex") {
        Some(column) => family
            .iter()
            .map(|row| ped_to_numeric(ped.cell(column, *row)))
            .collect::<Result<_>>()?,
        None => vec![PedNumber::Missing; family.len()],
    };
    // dict(zip(...)): the last row of an ID gives its sex.
    let mut sex_by_sample: HashMap<&str, PedNumber> = HashMap::new();
    for (id, sex) in ids.iter().zip(&sexes) {
        sex_by_sample.insert(id.as_str(), *sex);
    }
    let pick = |phenotype: i64| -> Vec<&str> {
        ids.iter()
            .zip(&phenotypes)
            .filter(|(id, value)| **value == Some(phenotype) && table.has_column(id))
            .map(|(id, _)| id.as_str())
            .collect()
    };
    let patients = pick(2);
    let controls = pick(1);
    if patients.is_empty() {
        return Ok(bs4_output(rows, vec![0; rows]));
    }
    let sex_is = |sample: &str, number: i64| {
        sex_by_sample
            .get(sample)
            .is_some_and(|sex| sex.equals(number))
    };
    let all = masks(table)?;
    let masks = Masks { all: &all };
    let valid = masks.valid_model()?;
    let chrom: Vec<String> = match table.column_index("chrom") {
        Some(column) => (0..rows)
            .map(|row| python_str(table.cell(column, row)).unwrap_or_default())
            .collect(),
        None => vec![String::new(); rows],
    };
    let autosomal: Vec<bool> = chrom
        .iter()
        .map(|text| !(text.contains('X') || text.contains('Y') || text.contains('M')))
        .collect();
    let x_linked: Vec<bool> = chrom.iter().map(|text| text.contains('X')).collect();
    let y_linked: Vec<bool> = chrom.iter().map(|text| text.contains('Y')).collect();
    let (recessive, x_recessive, x_dominant, dominant, y_model, dominant_like) = (
        masks.get("has_recessive")?,
        masks.get("has_x_linked_recessive")?,
        masks.get("has_x_linked_dominant")?,
        masks.get("has_dominant")?,
        masks.get("has_y_linked")?,
        masks.get("has_dominant_like")?,
    );
    let hom_hemi_requirement: Vec<bool> = (0..rows)
        .map(|row| recessive[row] || x_recessive[row])
        .collect();
    let hom_hemi_only: Vec<bool> = (0..rows)
        .map(|row| hom_hemi_requirement[row] && !(dominant[row] || x_dominant[row]))
        .collect();
    let dominant_affected_vs_affected: Vec<bool> = (0..rows)
        .map(|row| dominant_like[row] && !hom_hemi_requirement[row])
        .collect();
    let dominant_no_biallelic: Vec<bool> = (0..rows)
        .map(|row| dominant_like[row] && !recessive[row])
        .collect();
    let applicable = |sample: &str| -> Vec<bool> {
        (0..rows)
            .map(|row| {
                let (x_model, y) = if sex_is(sample, 1) {
                    (x_recessive[row] || x_dominant[row], y_model[row])
                } else if sex_is(sample, 2) {
                    (x_dominant[row], false)
                } else {
                    (false, false)
                };
                autosomal[row] || (x_linked[row] && x_model) || (y_linked[row] && y)
            })
            .collect()
    };
    // Dose and callability of each sample (`_gt_alt_count`, `_effective_alt_state`, `_gt_called`).
    let mut doses: HashMap<&str, (Vec<i64>, Vec<bool>)> = HashMap::new();
    for sample in patients.iter().chain(&controls) {
        if doses.contains_key(sample) {
            continue;
        }
        let column = table.column(sample)?;
        let mut dose = Vec::with_capacity(rows);
        let mut called = Vec::with_capacity(rows);
        for row in 0..rows {
            let text = dotted_text(column, row)?;
            let gt = text.split(':').next().unwrap_or_default();
            let mut count = gt.matches('1').count() as i64;
            if sex_is(sample, 1) && (x_linked[row] || y_linked[row]) && count >= 1 {
                count = 2;
            }
            let mut is_called = !gt.contains('.') && gt.contains(['0', '1']);
            if sex_is(sample, 2) {
                is_called = is_called && !y_linked[row];
            }
            dose.push(count);
            called.push(is_called);
        }
        doses.insert(sample, (dose, called));
    }
    let mendelian_locus: Vec<bool> = (0..rows)
        .map(|row| autosomal[row] || x_linked[row] || y_linked[row])
        .collect();
    let mut final_criteria = vec![false; rows];
    let mut supporting = vec![false; rows];
    let applicable_of: HashMap<&str, Vec<bool>> = patients
        .iter()
        .chain(&controls)
        .map(|sample| (*sample, applicable(sample)))
        .collect();
    for (index, sample_a) in patients.iter().enumerate() {
        for sample_b in &patients[index + 1..] {
            let (a_dose, a_called) = &doses[sample_a];
            let (b_dose, b_called) = &doses[sample_b];
            let (a_ok, b_ok) = (&applicable_of[sample_a], &applicable_of[sample_b]);
            for row in 0..rows {
                let called = a_called[row] && b_called[row] && valid[row] && a_ok[row] && b_ok[row];
                if !(called && mendelian_locus[row]) {
                    continue;
                }
                let (a, b) = (a_dose[row], b_dose[row]);
                let hom_wt = (a >= 2 && b == 0) || (b >= 2 && a == 0);
                let hom_het = (a >= 2 && b == 1) || (b >= 2 && a == 1);
                let het_wt = (a == 1 && b == 0) || (b == 1 && a == 0);
                final_criteria[row] |= hom_wt || (het_wt && dominant_affected_vs_affected[row]);
                supporting[row] |= hom_het && hom_hemi_only[row];
            }
        }
    }
    for patient in &patients {
        let (p_dose, p_called) = &doses[patient];
        for control in &controls {
            let (c_dose, c_called) = &doses[control];
            let (p_ok, c_ok) = (&applicable_of[patient], &applicable_of[control]);
            for row in 0..rows {
                let called = p_called[row] && c_called[row] && valid[row] && p_ok[row] && c_ok[row];
                if !(called && mendelian_locus[row]) {
                    continue;
                }
                let (p, c) = (p_dose[row], c_dose[row]);
                final_criteria[row] |= p >= 1 && c >= 2;
                supporting[row] |=
                    ((p >= 2 && c == 1) || (p == 1 && c == 1)) && dominant_no_biallelic[row];
            }
        }
    }
    let values = (0..rows)
        .map(|row| {
            if final_criteria[row] {
                3
            } else if supporting[row] {
                1
            } else {
                0
            }
        })
        .collect();
    Ok(bs4_output(rows, values))
}

/// A strength or flag array as numpy integers (`True` is 1).
fn integers(column: &Column, name: &str) -> Result<Vec<i64>> {
    match column {
        Column::Int(values) => Ok(values.clone()),
        Column::Bool(values) => Ok(values.iter().map(|flag| i64::from(*flag)).collect()),
        other => Err(refused(format!(
            "{name} of type {} is not a criterion array",
            other.dtype_name()
        ))),
    }
}

/// A LOF flag column for `|`: booleans, or an object column of booleans without missing values.
fn lof_flags(table: &Table, name: &str) -> Result<Vec<bool>> {
    match table.column(name)? {
        _ if table.rows == 0 => Ok(Vec::new()),
        Column::Bool(values) => Ok(values.clone()),
        Column::BoolOrMissing(values) if values.iter().all(Option::is_some) => {
            Ok(values.iter().map(|flag| flag.unwrap_or(false)).collect())
        }
        other => Err(refused(format!(
            "{name} of type {} with missing or non-boolean values (pandas' | treats them by side)",
            other.dtype_name()
        ))),
    }
}

/// `BP2_PM3_criteria(df, ped_df, pm2, pvs1, ps1, ps3)` with the step's masks.
pub(crate) fn bp2_pm3_criteria(
    table: &Table,
    ped: &Table,
    pm2: &Column,
    pvs1: &Column,
    ps1: &Column,
    ps3: &Column,
) -> Result<StepOutput> {
    bp2_pm3_criteria_with(table, ped, pm2, pvs1, ps1, ps3, &variant_condition_masks)
}

/// `BP2_PM3_criteria` with the condition masks taken from `masks`.
pub(crate) fn bp2_pm3_criteria_with(
    table: &Table,
    ped: &Table,
    pm2: &Column,
    pvs1: &Column,
    ps1: &Column,
    ps3: &Column,
    masks: MaskSource<'_>,
) -> Result<StepOutput> {
    let rows = table.rows;
    let vep = lof_flags(table, "vep_consq_lof")?;
    let splicing = lof_flags(table, "splicing_lof")?;
    let utr5 = lof_flags(table, "5UTR_lof")?;
    let (ps1, ps3, pvs1, pm2) = (
        integers(ps1, "PS1")?,
        integers(ps3, "PS3")?,
        integers(pvs1, "PVS1")?,
        integers(pm2, "PM2")?,
    );
    for array in [&ps1, &ps3, &pvs1, &pm2] {
        if array.len() != rows {
            return Err(invalid("a criterion array does not cover every row"));
        }
    }
    let pathogenic: Vec<bool> = (0..rows)
        .map(|row| {
            vep[row] || splicing[row] || utr5[row] || ps1[row] & ps3[row] != 0 || pvs1[row] != 0
        })
        .collect();
    let phase = determine_cis_trans_relationships(table, &pathogenic, ped)?;
    let mut written = table.clone();
    for (name, column) in &phase.columns {
        written.set_column(name, column.clone());
    }
    let all = masks(&written)?;
    let masks = Masks { all: &all };
    let valid = masks.valid_model()?;
    let (recessive, x_recessive, dominant_like) = (
        masks.get("has_recessive")?,
        masks.get("has_x_linked_recessive")?,
        masks.get("has_dominant_like")?,
    );
    let mut bp2 = vec![0i64; rows];
    let mut pm3 = vec![0i64; rows];
    for row in 0..rows {
        let cis = phase.in_cis_pathogenic[row] > 0;
        let trans = phase.in_trans_pathogenic[row] > 0;
        let benign = ((trans && dominant_like[row] && !(recessive[row] || x_recessive[row]))
            || (cis && recessive[row]))
            && valid[row];
        if benign {
            bp2[row] = 1;
        }
        if trans && recessive[row] && pm2[row] > 0 && valid[row] {
            pm3[row] = 2;
        }
    }
    Ok(StepOutput {
        criteria: vec![("BP2", Column::Int(bp2)), ("PM3", Column::Int(pm3))],
        columns: phase.columns,
        passed_on: Vec::new(),
    })
}

/// Whether a harness call names the masks Python computed; without them the
/// harness runs the plain functions, which call T3's adapter.
fn has_python_masks(call: &HarnessCall) -> bool {
    call.options.get("condition_masks_frame").is_some()
        || call.options.get("condition_masks_error").is_some()
}

/// The condition masks a harness call provides in place of T3's adapter: the
/// frame Python computed (`condition_masks_frame`), or the KeyError message
/// Python raised (`condition_masks_error`).
fn harness_masks(call: &HarnessCall) -> Result<impl Fn(&Table) -> Result<NamedMasks>> {
    let frame = match call.options.get("condition_masks_frame") {
        Some(Json::Str(path)) => Some(read_frame(Path::new(path))?),
        _ => None,
    };
    let error = call
        .options
        .get("condition_masks_error")
        .and_then(Json::as_str)
        .map(str::to_owned);
    Ok(move |_: &Table| -> Result<NamedMasks> {
        if let Some(message) = &error {
            return Err(AcmgError::new("KeyError", message.clone()).into());
        }
        let Some(frame) = &frame else {
            return Err(invalid("option condition_masks_frame is required"));
        };
        let masks = frame
            .names
            .iter()
            .zip(&frame.columns)
            .map(|(name, column)| match column {
                Column::Bool(values) => Ok((name.clone(), values.clone())),
                _ => Err(invalid("condition mask frames hold bool columns")),
            })
            .collect::<Result<_>>()?;
        Ok(NamedMasks { masks })
    })
}

/// A criterion array given as a JSON list of integers or booleans.
fn json_array(call: &HarnessCall, name: &str) -> Result<Column> {
    match call.options.get(name) {
        Some(Json::Array(items)) if items.iter().all(|item| matches!(item, Json::Bool(_))) => {
            Ok(Column::Bool(json_mask(call.options.get(name), name)?))
        }
        Some(Json::Array(items)) => items
            .iter()
            .map(|item| match item {
                Json::Int(number) => i64::try_from(*number)
                    .map_err(|_| invalid(format!("option {name}: integer out of range"))),
                _ => Err(invalid(format!("option {name} holds integers"))),
            })
            .collect::<Result<Vec<i64>>>()
            .map(Column::Int),
        _ => Err(invalid(format!("option {name} (a list) is required"))),
    }
}

/// `priva acmg-step acmg_pp1_bs4_bp2_pm3_family <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let table = &call.table;
    match function {
        "PP1_criteria" => {
            let (vcf, ped) = (
                harness_text(call, "multi_fam_vcf")?,
                harness_text(call, "multi_fam_ped")?,
            );
            let output = if has_python_masks(call) {
                pp1_criteria_with(table, vcf, ped, &harness_masks(call)?)?
            } else {
                pp1_criteria(table, vcf, ped)?
            };
            criteria_frame(&output, table.rows)
        }
        "BS4_criteria" => {
            let ped = harness_ped(call)?;
            let fam_name = harness_text(call, "fam_name")?;
            let output = if has_python_masks(call) {
                bs4_criteria_with(table, &ped, fam_name, &harness_masks(call)?)?
            } else {
                bs4_criteria(table, &ped, fam_name)?
            };
            criteria_frame(&output, table.rows)
        }
        "BP2_PM3_criteria" => {
            let ped = harness_ped(call)?;
            let arrays = [
                json_array(call, "pm2")?,
                json_array(call, "pvs1")?,
                json_array(call, "ps1")?,
                json_array(call, "ps3")?,
            ];
            let [pm2, pvs1, ps1, ps3] = &arrays;
            let output = if has_python_masks(call) {
                bp2_pm3_criteria_with(table, &ped, pm2, pvs1, ps1, ps3, &harness_masks(call)?)?
            } else {
                bp2_pm3_criteria(table, &ped, pm2, pvs1, ps1, ps3)?
            };
            let mut frame = criteria_frame(&output, table.rows)?;
            for (name, column) in written_columns(table, &output.columns) {
                frame.names.push(name);
                frame.columns.push(column);
            }
            Ok(frame)
        }
        _ => Err(super::types::not_converted(
            "acmg_pp1_bs4_bp2_pm3_family",
            function,
        )),
    }
}

/// A stand-in for T3's `_variant_condition_masks`, for the unit tests only:
/// exact semicolon tokens of the lower-cased text columns.
#[cfg(test)]
pub(crate) fn token_masks(table: &Table) -> Result<NamedMasks> {
    let text = |name: &str| -> Result<Vec<String>> {
        let column = table.column(name)?;
        Ok((0..table.rows)
            .map(|row| {
                python_str(column.cell(row))
                    .unwrap_or_default()
                    .to_lowercase()
            })
            .collect())
    };
    let inheritance = text("variant_inheritance")?;
    let penetrance = text("variant_penetrance")?;
    let has = |values: &[String], token: &str| -> Vec<bool> {
        values
            .iter()
            .map(|value| value.split(';').any(|part| part == token))
            .collect()
    };
    let names = [
        "recessive",
        "dominant",
        "x_linked_recessive",
        "x_linked_dominant",
        "x_linked_unspecified",
        "y_linked",
        "mitochondrial",
        "polygenic",
        "digenic",
        "oligogenic",
        "non_mendelian",
    ];
    let tokens: HashMap<&str, Vec<bool>> = names
        .iter()
        .map(|name| (*name, has(&inheritance, name)))
        .collect();
    let any = |parts: &[&str]| -> Vec<bool> {
        (0..table.rows)
            .map(|row| parts.iter().any(|part| tokens[part][row]))
            .collect()
    };
    Ok(NamedMasks {
        masks: vec![
            ("has_recessive".into(), tokens["recessive"].clone()),
            ("has_dominant".into(), tokens["dominant"].clone()),
            (
                "has_x_linked_recessive".into(),
                tokens["x_linked_recessive"].clone(),
            ),
            (
                "has_x_linked_dominant".into(),
                tokens["x_linked_dominant"].clone(),
            ),
            (
                "has_x_linked_unspecified".into(),
                tokens["x_linked_unspecified"].clone(),
            ),
            ("has_y_linked".into(), tokens["y_linked"].clone()),
            ("has_mitochondrial".into(), tokens["mitochondrial"].clone()),
            (
                "has_dominant_like".into(),
                any(&["dominant", "x_linked_dominant", "y_linked", "mitochondrial"]),
            ),
            (
                "has_mendelian".into(),
                any(&[
                    "recessive",
                    "dominant",
                    "x_linked_recessive",
                    "x_linked_dominant",
                    "x_linked_unspecified",
                    "y_linked",
                    "mitochondrial",
                ]),
            ),
            (
                "has_non_monogenic".into(),
                any(&["polygenic", "digenic", "oligogenic"]),
            ),
            ("has_non_mendelian".into(), tokens["non_mendelian"].clone()),
            (
                "has_incomplete_penetrance".into(),
                has(&penetrance, "incomplete"),
            ),
            (
                "has_complete_penetrance".into(),
                has(&penetrance, "complete"),
            ),
        ],
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn texts(values: &[&str]) -> Column {
        Column::Text(
            values
                .iter()
                .map(|value| Some((*value).to_owned()))
                .collect(),
        )
    }

    fn table(columns: Vec<(&str, Column)>) -> Table {
        Table {
            rows: columns[0].1.len(),
            names: columns.iter().map(|(name, _)| (*name).to_owned()).collect(),
            columns: columns.into_iter().map(|(_, column)| column).collect(),
        }
    }

    fn pedigree() -> Table {
        table(vec![
            ("#FamilyID", texts(&["F", "F", "F", "F"])),
            ("IndividualID", texts(&["kid", "dad", "mum", "sib"])),
            ("PaternalID", texts(&["dad", "0", "0", "dad"])),
            ("MaternalID", texts(&["mum", "0", "0", "mum"])),
            ("Sex", texts(&["1", "1", "2", "2"])),
            ("Phenotype", texts(&["2", "1", "1", "2"])),
        ])
    }

    #[test]
    fn bs4_branches_follow_the_docstring() {
        let variants = table(vec![
            ("chrom", texts(&["chr1", "chr1", "chr1", "chrX", "chr2"])),
            ("kid", texts(&["1/1", "0/1", "0/1", "0/1", "0/1"])),
            ("sib", texts(&["0/0", "0/0", "0/1", "0/0", "0/1"])),
            ("dad", texts(&["0/0", "0/0", "0/1", "0/1", "1/1"])),
            ("mum", texts(&["0/0", "0/0", "0/0", "0/0", "./."])),
            (
                "variant_inheritance",
                texts(&[
                    "recessive",
                    "dominant",
                    "dominant",
                    "x_linked_recessive",
                    "recessive;polygenic",
                ]),
            ),
            ("variant_penetrance", texts(&["complete", "", "", "", ""])),
        ]);
        let output = bs4_criteria_with(&variants, &pedigree(), "F", &token_masks).unwrap();
        // Row 0: affected hom against affected WT. Row 1: affected het against affected WT
        // (dominant). Row 2: affected het against unaffected het (dominant): supporting.
        // Row 3: the healthy father's chrX carrier counts 2 against the affected son.
        // Row 4: polygenic history fails the gate.
        assert_eq!(
            output.criterion("BS4"),
            Some(&Column::Int(vec![3, 3, 1, 3, 0]))
        );
        assert_eq!(
            bs4_criteria_with(&variants, &pedigree(), "G", &token_masks)
                .unwrap()
                .criterion("BS4"),
            Some(&Column::Int(vec![0; 5]))
        );
    }

    #[test]
    fn bp2_and_pm3_follow_the_phase_calls() {
        let variants = table(vec![
            ("variant_id", texts(&["v1", "v2"])),
            ("Gene", texts(&["A", "A"])),
            ("chrom", texts(&["chr1", "chr1"])),
            ("pos", Column::Int(vec![1, 2])),
            ("ref", texts(&["A", "C"])),
            ("alt", texts(&["G", "T"])),
            ("kid", texts(&["0/1", "0/1"])),
            ("dad", texts(&["0/1", "0/0"])),
            ("mum", texts(&["0/0", "0/1"])),
            ("sib", texts(&["0/0", "0/0"])),
            ("vep_consq_lof", Column::Bool(vec![false, true])),
            ("splicing_lof", Column::Bool(vec![false, false])),
            ("5UTR_lof", Column::Bool(vec![false, false])),
            ("variant_inheritance", texts(&["recessive", "recessive"])),
            ("variant_penetrance", texts(&["", ""])),
        ]);
        let zeros = Column::Int(vec![0, 0]);
        let output = bp2_pm3_criteria_with(
            &variants,
            &pedigree(),
            &Column::Int(vec![2, 0]),
            &zeros,
            &zeros,
            &zeros,
            &token_masks,
        )
        .unwrap();
        assert_eq!(output.criterion("PM3"), Some(&Column::Int(vec![2, 0])));
        assert_eq!(output.criterion("BP2"), Some(&Column::Int(vec![0, 0])));
        assert_eq!(output.columns[0].0, "is_pathogenic");
        assert_eq!(output.columns[1].0, "kid_phase_patho");
    }

    #[test]
    fn pp1_without_family_data_is_zero() {
        let variants = table(vec![("variant_id", texts(&["v1"]))]);
        let output = pp1_criteria_with(&variants, "", "x.ped", &token_masks).unwrap();
        assert_eq!(output.criterion("PP1"), Some(&Column::Int(vec![0])));
    }
}
