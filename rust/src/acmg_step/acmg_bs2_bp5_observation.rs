//! `scripts/acmg_bs2_bp5_observation.py` (460 lines): BS2 and BP5.
//!
//! Thread T5 (allele frequency, carrier observation and statistics). Converted: `BS2_criteria`,
//! `extract_and_summarize_alt_disease_variants` and `BP5_criteria`. BP5 runs `bcftools` through
//! `/bin/sh` with exactly the Python command texts (the owner's choice for BP5: spawn the same
//! pipeline), so the alternative-disease VCF is read by the same program; the pipeline's output
//! is then typed as `pd.read_csv` types it.
use super::acmg_pm2_bs1_ba1_frequency::{loci, option_ints, output_frame, scope_review_required};
use super::acmg_variant_mechanism::variant_condition_masks;
use super::types::{AcmgError, HarnessCall, StepOutput, not_converted};
use crate::cli::{Result, invalid};
use crate::json_text::Json;
use crate::pandas_table::{Column, Header, ObjectValue, ReadOptions, Table, parse_table_with};
use crate::python_text::float_repr;
use std::collections::HashMap;
use std::process::Command;

fn refused(message: impl Into<String>) -> Box<dyn std::error::Error> {
    invalid(format!("refused (not reproduced): {}", message.into()))
}

/// `pd.to_numeric(df.get(name, pd.Series(0, index=df.index)), errors="coerce").fillna(0)`.
fn numbers_filled(table: &Table, name: &str) -> Result<Vec<f64>> {
    let Some(index) = table.column_index(name) else {
        return Ok(vec![0.0; table.rows]);
    };
    match &table.columns[index] {
        Column::Float(values) => Ok(values
            .iter()
            .map(|value| if value.is_nan() { 0.0 } else { *value })
            .collect()),
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
            "{name} is a {} column, not numbers",
            other.dtype_name()
        ))),
    }
}

/// `BS2_criteria(df, pm2_criteria)`: 3 where a healthy carrier at a dose the inherited model
/// forbids is observed in gnomAD or the control cohort and PM2 is not set.
pub(crate) fn bs2_criteria(table: &Table, pm2_criteria: &[i64]) -> Result<StepOutput> {
    let rows = table.rows;
    if pm2_criteria.len() != rows {
        return Err(Box::new(AcmgError::new(
            "ValueError",
            format!(
                "operands could not be broadcast together with shapes ({rows},) ({},)",
                pm2_criteria.len()
            ),
        )));
    }
    let loci = loci(table)?;
    let scope_review = scope_review_required(table)?;
    let masks = variant_condition_masks(table)?;
    let complete_penetrance = masks.get("has_complete_penetrance")?;
    let has_mendelian = masks.get("has_mendelian")?;
    let has_non_monogenic = masks.get("has_non_monogenic")?;
    let has_non_mendelian = masks.get("has_non_mendelian")?;
    let has_incomplete = masks.get("has_incomplete_penetrance")?;
    let has_recessive = masks.get("has_recessive")?;
    let has_dominant = masks.get("has_dominant")?;
    let has_y_linked = masks.get("has_y_linked")?;
    let has_mitochondrial = masks.get("has_mitochondrial")?;
    let has_x_recessive = masks.get("has_x_linked_recessive")?;
    let has_x_dominant = masks.get("has_x_linked_dominant")?;

    let af = numbers_filled(table, "gnomAD_joint_AF")?;
    let an = numbers_filled(table, "gnomAD_joint_AN")?;
    let nhomalt_xx = numbers_filled(table, "gnomAD_nhomalt_XX")?;
    let nhomalt_xy = numbers_filled(table, "gnomAD_nhomalt_XY")?;
    let ac_xx = numbers_filled(table, "gnomAD_joint_AC_XX")?;
    let af_xy = numbers_filled(table, "gnomAD_joint_AF_XY")?;
    let control_ac = numbers_filled(table, "control_AC")?;
    let control_nhomalt = numbers_filled(table, "control_nhomalt")?;
    let hom_carrier_af: f64 = 1e-3;
    let hom_carrier_af_complete: f64 = 3e-4;
    let squared = hom_carrier_af.powf(2.0);
    let squared_complete = hom_carrier_af_complete.powf(2.0);

    let mut bs2 = vec![0i64; rows];
    for row in 0..rows {
        let valid_model = has_mendelian[row] && !has_non_monogenic[row] && !has_non_mendelian[row];
        let baseline = valid_model && !has_incomplete[row] && !scope_review[row];
        let complete = baseline && complete_penetrance[row];
        let non_x_single_copy = (loci.autosomal[row] && has_dominant[row] && !has_recessive[row])
            || (loci.y[row] && has_y_linked[row])
            || (loci.mitochondrial[row] && has_mitochondrial[row]);
        let x_recessive = loci.x[row] && has_x_recessive[row];
        let x_dominant = loci.x[row] && has_x_dominant[row] && !has_recessive[row];

        let gnomad_ac = af[row] * an[row];
        let hom_hemi = nhomalt_xx[row] + nhomalt_xy[row];
        let y_allele = loci.y[row] && af_xy[row] > 0.0;
        let individuals = an[row] / 2.0;
        let hom_min = individuals * squared;
        let hom_min_complete = individuals * squared_complete;
        let control_carrier = control_ac[row] > 0.0;
        let control_hom = control_nhomalt[row] > 0.0;
        let hom_evidence = (hom_hemi > hom_min && hom_hemi > 0.0) || control_hom;
        let carrier_evidence = gnomad_ac > 10.0 || control_carrier || y_allele;
        let complete_hom = (hom_hemi > hom_min_complete && hom_hemi > 0.0) || control_hom;
        let complete_carrier = gnomad_ac > 3.0 || control_carrier || y_allele;
        let standard = baseline
            && ((has_recessive[row] && hom_evidence)
                || (non_x_single_copy && (carrier_evidence || hom_evidence))
                || (x_recessive && nhomalt_xy[row] > 10.0)
                || (x_dominant && ac_xx[row] > 10.0));
        let relaxed = complete
            && ((has_recessive[row] && complete_hom)
                || (x_recessive && nhomalt_xy[row] > 3.0)
                || (x_dominant && ac_xx[row] > 3.0)
                || (non_x_single_copy && (complete_carrier || complete_hom)));
        if (standard || relaxed) && pm2_criteria[row] == 0 {
            bs2[row] = 3;
        }
    }
    Ok(StepOutput {
        criteria: vec![("BS2", Column::Int(bs2))],
        ..StepOutput::default()
    })
}

// ------------------------------------------------------------------------------------------------
// The alternative-disease VCF through bcftools

/// The commands of `extract_and_summarize_alt_disease_variants`, as Python formats them.
fn extract_command(vcf: &str, chrom: Option<&str>) -> String {
    const TAIL: &str = "bcftools +fill-tags - -Ou -- -t AC_Het,AC_Hom,AC_Hemi | bcftools filter -e '(INFO/AC_Het == 0) && (INFO/AC_Hom == 0) && (INFO/AC_Hemi == 0)' -Ou | bcftools query -f \"%CHROM:%POS:%REF-%ALT\\t%INFO/AC_Het\\t%INFO/AC_Hom\\t%INFO/AC_Hemi\\n\"";
    match chrom {
        None => format!(
            "bcftools +fill-tags {vcf} -Ou -- -t AC_Het,AC_Hom,AC_Hemi | bcftools filter -e '(INFO/AC_Het == 0) && (INFO/AC_Hom == 0) && (INFO/AC_Hemi == 0)' -Ou | bcftools query -f \"%CHROM:%POS:%REF-%ALT\\t%INFO/AC_Het\\t%INFO/AC_Hom\\t%INFO/AC_Hemi\\n\""
        ),
        Some(chrom) => format!("bcftools view -Ou {vcf} {chrom} | {TAIL}"),
    }
}

/// `subprocess.run(command, shell=True, capture_output=True, text=True, check=True).stdout`.
fn run_shell(command: &str) -> Result<String> {
    let output = Command::new("/bin/sh").arg("-c").arg(command).output()?;
    if !output.status.success() {
        let code = output
            .status
            .code()
            .map_or_else(|| "a signal".to_owned(), |code| code.to_string());
        return Err(Box::new(AcmgError::new(
            "CalledProcessError",
            format!("Command '{command}' returned non-zero exit status {code}."),
        )));
    }
    String::from_utf8(output.stdout).map_err(|_| {
        Box::new(AcmgError::new(
            "UnicodeDecodeError",
            "the bcftools output is not UTF-8",
        )) as Box<dyn std::error::Error>
    })
}

/// The `> 0` test of one column of the bcftools table, as pandas types the column.
fn positive(column: &Column) -> Result<Vec<bool>> {
    match column {
        Column::Int(values) => Ok(values.iter().map(|value| *value > 0).collect()),
        Column::Float(values) => Ok(values.iter().map(|value| *value > 0.0).collect()),
        Column::Bool(values) => Ok(values.clone()),
        _ => Err(Box::new(AcmgError::new(
            "TypeError",
            "'>' not supported between instances of 'str' and 'int'",
        ))),
    }
}

/// `extract_and_summarize_alt_disease_variants(alt_disease_vcf, chrom)`: variant id ->
/// (`AC_Het > 0`, `AC_Hom > 0 or AC_Hemi > 0`), in the order of the bcftools output with the
/// last of repeated ids kept (as `to_dict` keeps it). Empty when bcftools prints nothing.
pub(crate) fn extract_and_summarize_alt_disease_variants(
    vcf: &str,
    chrom: Option<&str>,
) -> Result<Vec<(String, bool, bool)>> {
    let stdout = run_shell(&extract_command(vcf, chrom))?;
    if stdout.trim().is_empty() {
        return Ok(Vec::new());
    }
    let names = ["variant_id", "AC_Het", "AC_Hom", "AC_Hemi"].map(str::to_owned);
    let options = ReadOptions {
        header: Header::Names(names.to_vec()),
        ..ReadOptions::default()
    };
    let table = parse_table_with(&stdout, &options)
        .map_err(|error| refused(format!("the bcftools output as a table: {error}")))?;
    if table.names.len() != 4 {
        return Err(refused("the bcftools output does not have four columns"));
    }
    let het = positive(&table.columns[1])?;
    let hom = positive(&table.columns[2])?;
    let hemi = positive(&table.columns[3])?;
    let ids: Vec<String> = match &table.columns[0] {
        Column::Text(values) => values
            .iter()
            .map(|value| {
                value
                    .clone()
                    .ok_or_else(|| refused("a missing variant id in the bcftools output"))
            })
            .collect::<Result<_>>()?,
        other => {
            return Err(refused(format!(
                "the variant ids of the bcftools output read as {}",
                other.dtype_name()
            )));
        }
    };
    let mut order: Vec<String> = Vec::new();
    let mut values: HashMap<String, (bool, bool)> = HashMap::new();
    for row in 0..table.rows {
        let id = ids[row].clone();
        if !values.contains_key(&id) {
            order.push(id.clone());
        }
        values.insert(id, (het[row], hom[row] || hemi[row]));
    }
    Ok(order
        .into_iter()
        .map(|id| {
            let (het, hom) = values[&id];
            (id, het, hom)
        })
        .collect())
}

/// `str(value)` of one cell for the `variant_id` text (`astype(str)`).
fn cell_text(column: &Column, row: usize) -> Option<String> {
    match column {
        Column::Text(values) => Some(values[row].clone().unwrap_or_else(|| "nan".to_owned())),
        Column::Int(values) => Some(values[row].to_string()),
        Column::Float(values) => Some(if values[row].is_nan() {
            "nan".to_owned()
        } else {
            float_repr(values[row])
        }),
        Column::Bool(values) => Some(if values[row] { "True" } else { "False" }.to_owned()),
        _ => None,
    }
}

/// An object column of text (`ref`, `alt`): its values, `None` where missing.
fn text_values<'a>(table: &'a Table, name: &str) -> Result<&'a [Option<String>]> {
    match table.column(name)? {
        Column::Text(values) => Ok(values),
        other => Err(refused(format!(
            "{name} is a {} column, not text",
            other.dtype_name()
        ))),
    }
}

/// The `alt_disease_hets` / `alt_disease_homs` columns as BP5 reads them: kept when bool,
/// turned into bool by `astype(str).str.lower().isin({"true", "1", "yes"})` when object.
fn observed_column(table: &Table, name: &str) -> Result<(Vec<bool>, Option<Column>)> {
    let truthy = |text: &str| matches!(text.to_lowercase().as_str(), "true" | "1" | "yes");
    match table.column(name)? {
        Column::Bool(values) => Ok((values.clone(), None)),
        Column::BoolOrMissing(values) => {
            let flags: Vec<bool> = values.iter().map(|value| *value == Some(true)).collect();
            Ok((flags.clone(), Some(Column::Bool(flags))))
        }
        Column::Text(values) => {
            let flags: Vec<bool> = values
                .iter()
                .map(|value| value.as_deref().is_some_and(truthy))
                .collect();
            Ok((flags.clone(), Some(Column::Bool(flags))))
        }
        Column::Object(values) => {
            let flags: Vec<bool> = values
                .iter()
                .map(|value| match value {
                    ObjectValue::Text(text) => truthy(text),
                    ObjectValue::Bool(flag) => *flag,
                    ObjectValue::Int(number) => *number == 1,
                    ObjectValue::Float(number) => float_repr(*number) == "1",
                    ObjectValue::Nan | ObjectValue::None => false,
                })
                .collect();
            Ok((flags.clone(), Some(Column::Bool(flags))))
        }
        other => Err(refused(format!(
            "{name} is a {} column; BP5's `|` on it is not reproduced",
            other.dtype_name()
        ))),
    }
}

/// `BP5_criteria(df, alt_disease_vcf, n_processes)`: 1 where the variant was seen in a patient
/// whose disease has another molecular explanation, at a dose the inherited model makes
/// informative. The columns written: `variant_id`, `alt_disease_hets`, `alt_disease_homs` when
/// the VCF is scanned; the two observation columns turned into bool when they are reused as
/// object columns.
pub(crate) fn bp5_criteria(
    table: &Table,
    alt_disease_vcf: &str,
    n_processes: usize,
) -> Result<StepOutput> {
    let rows = table.rows;
    let mut written: Vec<(String, Column)> = Vec::new();
    let (hets, homs);
    if table.has_column("alt_disease_hets") && table.has_column("alt_disease_homs") {
        let (het_values, het_column) = observed_column(table, "alt_disease_hets")?;
        let (hom_values, hom_column) = observed_column(table, "alt_disease_homs")?;
        if let Some(column) = het_column {
            written.push(("alt_disease_hets".to_owned(), column));
        }
        if let Some(column) = hom_column {
            written.push(("alt_disease_homs".to_owned(), column));
        }
        hets = het_values;
        homs = hom_values;
    } else {
        let listing = run_shell(&format!(
            "bcftools query -f \"%CHROM\\n\" {alt_disease_vcf} | sort -u"
        ))?;
        if listing.trim().is_empty() {
            return Ok(StepOutput {
                criteria: vec![("BP5", Column::Int(vec![0; rows]))],
                ..StepOutput::default()
            });
        }
        let chroms: Vec<String> = listing
            .trim()
            .split('\n')
            .map(|chrom| chrom.trim().to_owned())
            .filter(|chrom| !chrom.is_empty())
            .collect();
        // pool.map over the chromosomes: results in chromosome order, the first failure raised.
        let workers = n_processes.max(1);
        let mut results: Vec<Result<Vec<(String, bool, bool)>>> = Vec::new();
        for batch in chroms.chunks(workers) {
            let batch_results: Vec<Result<Vec<(String, bool, bool)>>> =
                std::thread::scope(|scope| {
                    let handles: Vec<_> = batch
                        .iter()
                        .map(|chrom| {
                            scope.spawn(move || {
                                extract_and_summarize_alt_disease_variants(
                                    alt_disease_vcf,
                                    Some(chrom),
                                )
                                .map_err(|error| error.to_string())
                            })
                        })
                        .collect();
                    handles
                        .into_iter()
                        .map(|handle| {
                            handle
                                .join()
                                .expect("a bcftools worker panicked")
                                .map_err(|message| invalid(message))
                        })
                        .collect()
                });
            results.extend(batch_results);
        }
        let mut het_by_id: HashMap<String, bool> = HashMap::new();
        let mut hom_by_id: HashMap<String, bool> = HashMap::new();
        for result in results {
            for (id, het, hom) in result? {
                het_by_id.insert(id.clone(), het);
                hom_by_id.insert(id, hom);
            }
        }
        let chrom = table.column("chrom")?;
        let pos = table.column("pos")?;
        let reference = text_values(table, "ref")?;
        let alternative = text_values(table, "alt")?;
        let mut ids: Vec<Option<String>> = Vec::with_capacity(rows);
        for row in 0..rows {
            let chrom_text = cell_text(chrom, row)
                .ok_or_else(|| refused(format!("chrom is a {} column", chrom.dtype_name())))?;
            let pos_text = cell_text(pos, row)
                .ok_or_else(|| refused(format!("pos is a {} column", pos.dtype_name())))?;
            ids.push(match (&reference[row], &alternative[row]) {
                (Some(r), Some(a)) => Some(format!("{chrom_text}:{pos_text}:{r}-{a}")),
                _ => None,
            });
        }
        let lookup = |by_id: &HashMap<String, bool>| -> Vec<bool> {
            ids.iter()
                .map(|id| {
                    id.as_ref()
                        .and_then(|id| by_id.get(id))
                        .copied()
                        .unwrap_or(false)
                })
                .collect()
        };
        hets = lookup(&het_by_id);
        homs = lookup(&hom_by_id);
        written.push(("variant_id".to_owned(), Column::Text(ids)));
        written.push(("alt_disease_hets".to_owned(), Column::Bool(hets.clone())));
        written.push(("alt_disease_homs".to_owned(), Column::Bool(homs.clone())));
    }

    let masks = variant_condition_masks(table)?;
    let has_mendelian = masks.get("has_mendelian")?;
    let has_non_monogenic = masks.get("has_non_monogenic")?;
    let has_non_mendelian = masks.get("has_non_mendelian")?;
    let has_incomplete = masks.get("has_incomplete_penetrance")?;
    let has_dominant_like = masks.get("has_dominant_like")?;
    let has_recessive = masks.get("has_recessive")?;
    let has_x_recessive = masks.get("has_x_linked_recessive")?;
    let bp5 = (0..rows)
        .map(|row| {
            let eligible = has_mendelian[row]
                && !has_non_monogenic[row]
                && !has_non_mendelian[row]
                && !has_incomplete[row];
            let heterozygous_disease = has_dominant_like[row] && !has_recessive[row];
            let dominant = heterozygous_disease && eligible && (hets[row] || homs[row]);
            let hom_hemi_disease = has_recessive[row] || has_x_recessive[row];
            let recessive = hom_hemi_disease && eligible && homs[row];
            i64::from(dominant || recessive)
        })
        .collect();
    Ok(StepOutput {
        criteria: vec![("BP5", Column::Int(bp5))],
        columns: written,
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

/// `priva acmg-step acmg_bs2_bp5_observation <function>`:
/// - `BS2_criteria`, option `pm2_criteria` (list of integers);
/// - `BP5_criteria`, options `alt_disease_vcf`, `n_processes` (default 1);
/// - `extract_and_summarize_alt_disease_variants`, options `alt_disease_vcf`, `chrom` (text or
///   null): the frame of `variant_id`, `het`, `hom`.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let options = &call.options;
    match function {
        "BS2_criteria" => {
            let pm2 = option_ints(options, "pm2_criteria")?
                .ok_or_else(|| invalid("option pm2_criteria is required"))?;
            Ok(output_frame(&call.table, bs2_criteria(&call.table, &pm2)?))
        }
        "BP5_criteria" => {
            let processes = match options.get("n_processes") {
                Some(Json::Int(number)) => usize::try_from(*number).unwrap_or(1),
                _ => 1,
            };
            let output = bp5_criteria(
                &call.table,
                option_text(options, "alt_disease_vcf")?,
                processes,
            )?;
            Ok(output_frame(&call.table, output))
        }
        "extract_and_summarize_alt_disease_variants" => {
            let chrom = match options.get("chrom") {
                Some(Json::Str(text)) => Some(text.as_str()),
                _ => None,
            };
            let entries = extract_and_summarize_alt_disease_variants(
                option_text(options, "alt_disease_vcf")?,
                chrom,
            )?;
            Ok(Table {
                names: vec!["variant_id".into(), "het".into(), "hom".into()],
                rows: entries.len(),
                columns: vec![
                    Column::Text(entries.iter().map(|entry| Some(entry.0.clone())).collect()),
                    Column::Bool(entries.iter().map(|entry| entry.1).collect()),
                    Column::Bool(entries.iter().map(|entry| entry.2).collect()),
                ],
            })
        }
        _ => Err(not_converted("acmg_bs2_bp5_observation", function)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_commands_are_pythons_texts() {
        assert_eq!(
            extract_command("a.vcf.gz", Some("chr2")),
            "bcftools view -Ou a.vcf.gz chr2 | bcftools +fill-tags - -Ou -- -t AC_Het,AC_Hom,AC_Hemi | bcftools filter -e '(INFO/AC_Het == 0) && (INFO/AC_Hom == 0) && (INFO/AC_Hemi == 0)' -Ou | bcftools query -f \"%CHROM:%POS:%REF-%ALT\\t%INFO/AC_Het\\t%INFO/AC_Hom\\t%INFO/AC_Hemi\\n\""
        );
        assert!(extract_command("a.vcf.gz", None).starts_with("bcftools +fill-tags a.vcf.gz -Ou"));
    }

    #[test]
    fn text_counts_refuse_like_pandas() {
        let error = positive(&Column::Text(vec![Some("1,0".into())]))
            .unwrap_err()
            .to_string();
        assert!(error.contains("TypeError"), "{error}");
        assert_eq!(
            positive(&Column::Float(vec![f64::NAN, 2.0])).unwrap(),
            vec![false, true]
        );
    }
}
