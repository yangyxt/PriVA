//! `scripts/acmg_pm4_bp3_protein_length.py` (224 lines): PM4 and BP3 and the repeat-region overlap.
//!
//! `pm4_bp3_criteria` is `PM4_BP3_criteria`; `variants_in_repeats` is
//! `find_overlaps_bedtools_efficient(df, regions_file, method="all")`, which Python runs through
//! pybedtools as `bedtools intersect -a variants.bed -b regions -wa -f 1 -u`: the variant's
//! interval `[pos - 1, pos - 1 + len(ref))` must lie inside one region of the BED file. Here
//! that containment test is done natively: per chromosome, the regions sorted by start with the
//! running maximum of their ends; a variant is inside a region when the largest end among the
//! regions starting at or before the variant's start reaches the variant's end.
//!
//! Refused instead of guessed: a region of length zero (bedtools widens those), a region or a
//! BED line bedtools would reject, a variant position below 1, an empty or missing `ref`, and
//! names that `to_csv` would quote (tab, line break, double quote).
use super::acmg_pp3_bp4_bp7_insilico::harness_options;
use super::acmg_pp3_bp4_bp7_insilico::pandas_masks::{
    Mask, bool_input_fillna_false, strengths, text_contains,
};
use super::types::{HarnessCall, StepOutput, not_converted};
use crate::cli::{Result, invalid};
use crate::pandas_table::{Cell, Column, Table};
use flate2::read::MultiGzDecoder;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

/// One variant as a BED interval (`BedTool.from_dataframe` of chrom, start, end, variant_id).
struct VariantInterval {
    chrom: String,
    start: i64,
    end: i64,
    name: String,
}

/// The text `to_csv(na_rep='.')` writes for a text cell, refusing what it would quote.
fn bed_field(cell: Cell<'_>, column: &str) -> Result<String> {
    let text = match cell {
        Cell::Missing => ".".to_owned(),
        Cell::Text(text) => text.to_owned(),
        other => {
            return Err(invalid(format!(
                "{column} holds {other:?}; the BED file is written for text values only"
            )));
        }
    };
    if text.contains(['\t', '\n', '\r', '"']) {
        return Err(invalid(format!(
            "{column} value {text:?} would be quoted in the BED file pybedtools writes"
        )));
    }
    Ok(text)
}

fn variant_intervals(table: &Table) -> Result<Vec<VariantInterval>> {
    let chrom = table.column("chrom")?;
    let pos = table.column("pos")?;
    let variant_id = table.column("variant_id")?;
    let Column::Int(positions) = pos else {
        return Err(invalid(format!(
            "pos is a pandas {} column; the BED intervals are reproduced for int64 positions only",
            pos.dtype_name()
        )));
    };
    let reference = table.column("ref")?;
    (0..table.rows)
        .map(|row| {
            let start = positions[row] - 1;
            if start < 0 {
                return Err(invalid(format!(
                    "pos {} gives a negative BED start, which bedtools rejects",
                    positions[row]
                )));
            }
            let length = match reference.cell(row) {
                Cell::Text(text) if !text.is_empty() => text.chars().count() as i64,
                other => {
                    return Err(invalid(format!(
                        "ref {other:?} gives no interval length; not reproduced"
                    )));
                }
            };
            Ok(VariantInterval {
                chrom: bed_field(chrom.cell(row), "chrom")?,
                start,
                end: start + length,
                name: bed_field(variant_id.cell(row), "variant_id")?,
            })
        })
        .collect()
}

/// The regions of the chromosomes in `wanted`, sorted by start, with the running maximum of
/// the ends.
struct Regions {
    starts: Vec<i64>,
    max_end: Vec<i64>,
}

fn parse_coordinate(text: &str, line: usize) -> Result<i64> {
    if text.is_empty() || !text.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(invalid(format!(
            "repeat BED line {line}: coordinate {text:?} is not a non-negative integer"
        )));
    }
    text.parse::<i64>().map_err(|_| {
        invalid(format!(
            "repeat BED line {line}: coordinate {text:?} too large"
        ))
    })
}

fn read_regions(path: &str, wanted: &HashSet<&str>) -> Result<HashMap<String, Regions>> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 2];
    let gzip = file.read(&mut magic)? == 2 && magic == [0x1f, 0x8b];
    let file = File::open(path)?;
    let reader: Box<dyn BufRead> = if gzip {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let mut intervals: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
    for (number, line) in reader.split(b'\n').enumerate() {
        let line = line?;
        let line = line.strip_suffix(b"\r").unwrap_or(&line);
        if line.is_empty()
            || line.starts_with(b"#")
            || line.starts_with(b"track")
            || line.starts_with(b"browser")
        {
            continue;
        }
        let text = std::str::from_utf8(line)
            .map_err(|_| invalid(format!("repeat BED line {}: not UTF-8", number + 1)))?;
        let mut fields = text.split('\t');
        let (Some(chrom), Some(start), Some(end)) = (fields.next(), fields.next(), fields.next())
        else {
            return Err(invalid(format!(
                "repeat BED line {}: fewer than three fields",
                number + 1
            )));
        };
        if !wanted.contains(chrom) {
            continue;
        }
        let start = parse_coordinate(start, number + 1)?;
        let end = parse_coordinate(end, number + 1)?;
        if end <= start {
            return Err(invalid(format!(
                "repeat BED line {}: region of length {} (bedtools widens empty regions); \
                 not reproduced",
                number + 1,
                end - start
            )));
        }
        intervals
            .entry(chrom.to_owned())
            .or_default()
            .push((start, end));
    }
    Ok(intervals
        .into_iter()
        .map(|(chrom, mut regions)| {
            regions.sort_unstable();
            let mut running = i64::MIN;
            let max_end = regions
                .iter()
                .map(|(_, end)| {
                    running = running.max(*end);
                    running
                })
                .collect();
            let starts = regions.iter().map(|(start, _)| *start).collect();
            (chrom, Regions { starts, max_end })
        })
        .collect())
}

/// `find_overlaps_bedtools_efficient(df, regions_file, method="all")`: the variant IDs whose
/// interval lies inside one region.
pub(crate) fn variants_in_repeats(table: &Table, regions_file: &str) -> Result<HashSet<String>> {
    let variants = variant_intervals(table)?;
    let wanted: HashSet<&str> = variants
        .iter()
        .map(|variant| variant.chrom.as_str())
        .collect();
    let regions = if variants.is_empty() {
        // pybedtools still opens the regions file.
        File::open(regions_file)?;
        HashMap::new()
    } else {
        read_regions(regions_file, &wanted)?
    };
    let mut inside = HashSet::new();
    for variant in &variants {
        let Some(chromosome) = regions.get(&variant.chrom) else {
            continue;
        };
        // Regions starting at or before the variant's start.
        let count = chromosome
            .starts
            .partition_point(|start| *start <= variant.start);
        if count > 0 && chromosome.max_end[count - 1] >= variant.end {
            inside.insert(variant.name.clone());
        }
    }
    Ok(inside)
}

/// `PM4_BP3_criteria(df, pvs1_criteria, repeat_regions_file, loc_intol_domain)`: PM4 (strength
/// 2) and BP3 (strength 1); the column `variant_in_repeat_region` the step then sets.
pub(crate) fn pm4_bp3_criteria(
    table: &Table,
    pvs1_criteria: &[i64],
    repeat_regions_file: &str,
    loc_intol_domain: &[bool],
) -> Result<StepOutput> {
    let rows = table.rows;
    if pvs1_criteria.len() != rows || loc_intol_domain.len() != rows {
        return Err(invalid(
            "pvs1_criteria and loc_intol_domain must cover every row",
        ));
    }
    let in_repeat_ids = variants_in_repeats(table, repeat_regions_file)?;
    let variant_id = table.column("variant_id")?;
    let in_repeat: Vec<bool> = (0..rows)
        .map(|row| matches!(variant_id.cell(row), Cell::Text(id) if in_repeat_ids.contains(id)))
        .collect();
    let low_complexity = text_contains(table, "DOMAINS", Some("."), &["Low_complexity"])?.remove(0);
    let in_repetitive_region = Mask::Plain(in_repeat.clone()).or(&low_complexity);

    let mut filled = text_contains(
        table,
        "Consequence",
        Some(""),
        &[
            "inframe_deletion",
            "inframe_insertion",
            "frameshift",
            "stop_gained",
        ],
    )?
    .into_iter();
    let mut next = || filled.next().expect("four patterns");
    let (deletion, insertion, frameshift, stop_gained) = (next(), next(), next(), next());
    let inframe_indels = deletion.or(&insertion);
    let splicing_inframe = bool_input_fillna_false(table, "splicing_len_changing")?
        .and(&bool_input_fillna_false(table, "splicing_frameshift")?.invert()?);
    let utr_inframe = bool_input_fillna_false(table, "5UTR_len_changing")?
        .and(&bool_input_fillna_false(table, "5UTR_frameshift")?.invert()?);
    let vep_inframe = bool_input_fillna_false(table, "vep_consq_len_changing")?
        .and(&frameshift.invert()?)
        .and(&stop_gained.invert()?);
    let all_inframe = inframe_indels
        .or(&splicing_inframe)
        .or(&utr_inframe)
        .or(&vep_inframe);

    let stop_loss = text_contains(table, "Consequence", None, &["stop_lost"])?
        .remove(0)
        .fillna(false);
    let escaping = text_contains(table, "NMD", Some("."), &["escaping"])?.remove(0);
    let end_trunc = text_contains(table, "LoF_filter", Some("."), &["END_TRUNC"])?.remove(0);
    let nmd_escaping_trunc = frameshift.or(&stop_gained).and(&escaping.or(&end_trunc));

    let has_functional_domain = Mask::Plain(loc_intol_domain.to_vec())
        .or(&bool_input_fillna_false(
            table,
            "splicing_span_intol_domain",
        )?)
        .or(&bool_input_fillna_false(table, "5UTR_span_intol_domain")?)
        .or(&bool_input_fillna_false(table, "span_functional_domains")?);

    let pvs1 = |test: fn(i64) -> bool| {
        Mask::Plain(pvs1_criteria.iter().map(|value| test(*value)).collect())
    };
    let not_stop_loss = stop_loss.invert()?;
    let frameshift_with_pvs1 = frameshift.and(&pvs1(|value| value > 0)).and(&not_stop_loss);
    let strong_nmd = pvs1(|value| value >= 3).and(&not_stop_loss);
    let exclude_from_pm4 = frameshift_with_pvs1.or(&strong_nmd);

    let pm4_inframe_eligible = all_inframe
        .and(&exclude_from_pm4.invert()?)
        .and(&in_repetitive_region.invert()?.or(&has_functional_domain));
    let pm4_eligible = pm4_inframe_eligible.or(&stop_loss).or(&nmd_escaping_trunc);
    let bp3_eligible = all_inframe
        .and(&in_repetitive_region)
        .and(&has_functional_domain.invert()?)
        .and(&pm4_eligible.invert()?);
    let pm4 = strengths(&pm4_eligible.and(&pvs1(|value| value < 2)), 2)?;
    let bp3_mask = bp3_eligible.and(&Mask::Plain(pm4.iter().map(|value| *value == 0).collect()));
    let bp3 = strengths(&bp3_mask, 1)?;
    Ok(StepOutput {
        criteria: vec![("PM4", Column::Int(pm4)), ("BP3", Column::Int(bp3))],
        columns: vec![(
            "variant_in_repeat_region".to_owned(),
            Column::Bool(in_repeat),
        )],
        ..StepOutput::default()
    })
}

/// `priva acmg-step acmg_pm4_bp3_protein_length <function>`:
///
/// - `PM4_BP3_criteria`, options `pvs1_criteria`, `loc_intol_domain`, `repeat_regions_file`:
///   a frame of `PM4`, `BP3`, `variant_in_repeat_region`;
/// - `find_overlaps_bedtools_efficient`, option `regions_file`: a frame of the overlapping
///   `variant_id`s, sorted.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let table = &call.table;
    let options = &call.options;
    match function {
        "PM4_BP3_criteria" => {
            let pvs1 = harness_options::int_list(options, "pvs1_criteria")?
                .ok_or_else(|| invalid("option pvs1_criteria is required"))?;
            let loc = harness_options::bool_list(options, "loc_intol_domain")?
                .ok_or_else(|| invalid("option loc_intol_domain is required"))?;
            let path = harness_options::required_text(options, "repeat_regions_file")?;
            let output = pm4_bp3_criteria(table, &pvs1, path, &loc)?;
            let mut columns: Vec<(String, Column)> = output
                .criteria
                .iter()
                .map(|(name, values)| ((*name).to_owned(), values.clone()))
                .collect();
            columns.extend(output.columns.iter().cloned());
            Ok(harness_options::frame(columns, table.rows))
        }
        "find_overlaps_bedtools_efficient" => {
            let path = harness_options::required_text(options, "regions_file")?;
            let mut ids: Vec<String> = variants_in_repeats(table, path)?.into_iter().collect();
            ids.sort();
            let rows = ids.len();
            Ok(harness_options::frame(
                vec![(
                    "variant_id".into(),
                    Column::Text(ids.into_iter().map(Some).collect()),
                )],
                rows,
            ))
        }
        _ => Err(not_converted("acmg_pm4_bp3_protein_length", function)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn containment_needs_one_region_holding_the_whole_variant() {
        let directory = std::env::temp_dir().join(format!("priva_t6_bed_{}", std::process::id()));
        std::fs::create_dir_all(&directory).unwrap();
        let bed = directory.join("repeats.bed");
        std::fs::write(
            &bed,
            "chr1\t10\t20\nchr1\t15\t30\nchr1\t40\t41\nchr2\t0\t5\n",
        )
        .unwrap();
        let text =
            |values: &[&str]| Column::Text(values.iter().map(|v| Some((*v).to_owned())).collect());
        let table = Table {
            names: vec![
                "chrom".into(),
                "pos".into(),
                "ref".into(),
                "variant_id".into(),
            ],
            columns: vec![
                text(&["chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr3"]),
                // BED starts 10, 19, 25, 39, 40, 4, 0
                Column::Int(vec![11, 20, 26, 40, 41, 5, 1]),
                text(&["AAAAAAAAAAAAAAAAAAAA", "A", "AAAAAA", "A", "A", "A", "A"]),
                text(&["v1", "v2", "v3", "v4", "v5", "v6", "v7"]),
            ],
            rows: 7,
        };
        let mut inside: Vec<String> = variants_in_repeats(&table, bed.to_str().unwrap())
            .unwrap()
            .into_iter()
            .collect();
        inside.sort();
        // v1 [10, 30) spans two regions but lies in neither; v3 [25, 31) ends past 30;
        // v4 [39, 40) ends where [40, 41) starts.
        assert_eq!(inside, vec!["v2", "v5", "v6"]);
        std::fs::remove_dir_all(&directory).unwrap();
    }
}
