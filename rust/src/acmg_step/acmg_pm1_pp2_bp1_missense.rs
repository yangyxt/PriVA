//! `scripts/acmg_pm1_pp2_bp1_missense.py` (307 lines): PM1, PP2 and BP1.
//!
//! `pm1_criteria` is `PM1_criteria` with `locate_pm1_region`; `pp2_bp1_criteria` is
//! `PP2_BP1_criteria` with `analyze_bp1_pp2`. The per-gene PP2 test calls SciPy's
//! `chi2_contingency` (Yates-corrected) or `fisher_exact(alternative="greater")` on a 2x2 table
//! of counts; those are `stats::chi2_contingency_2x2` and `stats::fisher_exact_greater` (thread
//! T5). The p-value source is a parameter (`PValues`), so the comparison harness can hand in the
//! p-values SciPy gives until `stats` lands; the step uses `PValues::Scipy`.
//!
//! Python runs the genes in a process pool (`imap_unordered`) and keys the results by each
//! entry's `ensembl_id`; two genes whose entries name one `ensembl_id` with different results
//! would make the outcome depend on completion order, and are refused.
use super::acmg_pp3_bp4_bp7_insilico::harness_options;
use super::acmg_pp3_bp4_bp7_insilico::pandas_masks::{
    Mask, numeric_column, python_str_filled, strengths, text_contains, text_values,
};
use super::stats;
use super::types::{AcmgError, HarnessCall, StepOutput, not_converted};
use crate::cli::{Result, invalid};
use crate::pandas_table::{Cell, Column, Table, read_table};
use crate::pickle_reader::{self, Set, Value};
use std::collections::HashMap;
use std::ops::Deref;
use std::path::Path;

/// Loads a pickle as `pickle.load(gzip.open(path) if path.endswith(".gz") else open(path, "rb"))`.
pub(crate) fn load_pickle(path: &str) -> Result<Value> {
    pickle_reader::load(&pickle_reader::read_file(
        Path::new(path),
        path.ends_with(".gz"),
    )?)
}

fn type_error(message: impl Into<String>) -> Box<dyn std::error::Error> {
    Box::new(AcmgError::new("TypeError", message))
}

fn key_error(key: &str) -> Box<dyn std::error::Error> {
    Box::new(AcmgError::new("KeyError", format!("'{key}'")))
}

/// `dict.get(key)` of a pickled dict with a text key.
fn dict_get(value: &Value, key: &str) -> Result<Option<Value>> {
    match value {
        Value::Dict(dict) => Ok(dict.borrow().get(&Value::Str(key.into()))?.cloned()),
        other => Err(invalid(format!(
            "expected a pickled dict, found {}",
            other.type_name()
        ))),
    }
}

/// `item in container` for a pickled set, frozenset, list, tuple or dict.
fn python_contains(container: &Value, item: &Value) -> Result<bool> {
    match container {
        Value::Set(set) => set.borrow().contains(item),
        Value::FrozenSet(set) => set.contains(item),
        Value::Dict(dict) => Ok(dict.borrow().get(item)?.is_some()),
        Value::List(items) => {
            for candidate in items.borrow().iter() {
                if candidate.python_equals(item)? {
                    return Ok(true);
                }
            }
            Ok(false)
        }
        Value::Tuple(items) => {
            for candidate in items.iter() {
                if candidate.python_equals(item)? {
                    return Ok(true);
                }
            }
            Ok(false)
        }
        other => Err(invalid(format!(
            "`in` on a pickled {} is not reproduced",
            other.type_name()
        ))),
    }
}

/// A table value as the dict key Python looks up (`None` for a missing value, which as a
/// fresh NaN never finds a key).
fn key_of(cell: Cell<'_>) -> Option<Value> {
    match cell {
        Cell::Missing => None,
        Cell::Text(text) => Some(Value::Str(text.into())),
        Cell::Int(number) => Some(Value::Int(number)),
        Cell::Float(number) => Some(Value::Float(number)),
        Cell::Bool(flag) => Some(Value::Bool(flag)),
    }
}

/// The protein position `locate_pm1_region` reads from `Protein_position`: the digits before
/// the first `/`, or `None` (no position).
fn protein_position(cell: Option<Cell<'_>>) -> Result<Option<i64>> {
    let text = match cell {
        // `row.get('Protein_position', '')` of an absent column, or a falsy value.
        None | Some(Cell::Missing) | Some(Cell::Bool(_)) | Some(Cell::Float(_)) => return Ok(None),
        Some(Cell::Int(0)) => return Ok(None),
        Some(Cell::Int(number)) => number.to_string(),
        Some(Cell::Text(text)) => text.to_owned(),
    };
    let position = text.split('/').next().unwrap_or("");
    if position.is_empty()
        || position
            .bytes()
            .any(|byte| byte.is_ascii() && !byte.is_ascii_digit())
    {
        return Ok(None);
    }
    if !position.is_ascii() {
        return Err(invalid(format!(
            "Protein_position {text:?}: str.isdigit() and int() of non-ASCII digits are not reproduced"
        )));
    }
    position
        .parse::<i64>()
        .map(Some)
        .map_err(|_| invalid(format!("Protein_position {text:?} is beyond 64 bits")))
}

/// `locate_pm1_region(row, hcseeker_hotspots, rmc_constrained)`: `(in_hcseeker, in_rmc)`.
fn locate_pm1_region(
    feature: Option<Cell<'_>>,
    protein: Option<Cell<'_>>,
    hcseeker_hotspots: &Value,
    rmc_constrained: &Value,
) -> Result<(bool, bool)> {
    let Some(position) = protein_position(protein)? else {
        return Ok((false, false));
    };
    // `row.get('Feature', '')`.
    let enst = match feature {
        None => Some(Value::Str("".into())),
        Some(cell) => key_of(cell),
    };
    let position = Value::Int(position);
    let hit = |regions: &Value| -> Result<bool> {
        let Value::Dict(dict) = regions else {
            return Err(type_error(format!(
                "PM1 regions: a {} where a dict of transcripts is expected",
                regions.type_name()
            )));
        };
        let Some(enst) = &enst else {
            return Ok(false);
        };
        match dict.borrow().get(enst)? {
            Some(residues) => python_contains(residues, &position),
            None => Ok(false),
        }
    };
    Ok((hit(hcseeker_hotspots)?, hit(rmc_constrained)?))
}

/// `PM1_criteria(df, pvs1_criteria, loc_intol_domain, clinvar_patho, pm1_regions_pkl)`: PM1
/// (strength 2), and `loc_intol_domain` passed on unchanged (PM4/BP3 reads it).
pub(crate) fn pm1_criteria(
    table: &Table,
    pvs1_criteria: &[i64],
    loc_intol_domain: &[bool],
    pm1_regions_pkl: &str,
) -> Result<StepOutput> {
    let rows = table.rows;
    if pvs1_criteria.len() != rows || loc_intol_domain.len() != rows {
        return Err(invalid(
            "pvs1_criteria and loc_intol_domain must cover every row",
        ));
    }
    let missense = text_contains(table, "Consequence", None, &["missense_variant"])?.remove(0);
    // `missense_damaging` only reaches the log, but its columns are read (and can raise) here.
    text_values(table, "am_class", Some(""))?;
    numeric_column(table, "PrimateAI")?;
    let am_pathogenicity = numeric_column(table, "am_pathogenicity")?;
    let primateai = numeric_column(table, "PrimateAI")?;
    let am_says_benign = text_contains(table, "am_class", Some(""), &["benign"])?
        .remove(0)
        .or(&Mask::Plain(
            am_pathogenicity
                .iter()
                .map(|value| *value < 0.564)
                .collect(),
        ));
    let primateai_not_damaging =
        Mask::Plain(primateai.iter().map(|value| !(*value >= 0.8)).collect());
    let missense_benign = am_says_benign.and(&primateai_not_damaging);

    let regions = load_pickle(pm1_regions_pkl)?;
    if !matches!(regions, Value::Dict(_)) {
        return Err(invalid("the PM1 regions pickle is not a dict"));
    }
    let empty = || Value::Dict(Default::default());
    let hcseeker_hotspots = dict_get(&regions, "hcseeker_hotspots")?.unwrap_or_else(empty);
    let rmc_constrained = dict_get(&regions, "rmc_constrained")?.unwrap_or_else(empty);
    let feature = table.column_index("Feature");
    let protein = table.column_index("Protein_position");
    let mut hcseeker_hit = Vec::with_capacity(rows);
    let mut rmc_hit = Vec::with_capacity(rows);
    for row in 0..rows {
        let (in_hcseeker, in_rmc) = locate_pm1_region(
            feature.map(|index| table.columns[index].cell(row)),
            protein.map(|index| table.columns[index].cell(row)),
            &hcseeker_hotspots,
            &rmc_constrained,
        )?;
        hcseeker_hit.push(in_hcseeker);
        rmc_hit.push(in_rmc);
    }

    // `df["Consequence"].fillna("").astype(str)` then two literal searches.
    let consequence_text = python_str_filled(table.column("Consequence")?, "Consequence", "")?;
    let indels = Mask::Plain(
        consequence_text
            .iter()
            .map(|text| text.contains("inframe_deletion") || text.contains("inframe_insertion"))
            .collect(),
    );
    let eligible = missense.and(&missense_benign.invert()?).or(&indels);
    let pm1 = Mask::Plain(hcseeker_hit)
        .and(&eligible)
        .or(&Mask::Plain(rmc_hit).and(&eligible))
        .or(&Mask::Plain(loc_intol_domain.to_vec()).and(&eligible))
        .and(&Mask::Plain(
            pvs1_criteria.iter().map(|value| *value < 2).collect(),
        ));
    Ok(StepOutput {
        criteria: vec![("PM1", Column::Int(strengths(&pm1, 2)?))],
        passed_on: vec![(
            "loc_intol_domain".to_owned(),
            Column::Bool(loc_intol_domain.to_vec()),
        )],
        ..StepOutput::default()
    })
}

/// Which SciPy test `analyze_bp1_pp2` runs on its 2x2 table.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Pp2Test {
    ChiSquare,
    Fisher,
}

/// Where the PP2 p-values come from: SciPy's functions in `stats` (the step), or p-values
/// supplied by the comparison harness, keyed by the test and the table.
pub(crate) enum PValues {
    Scipy,
    Supplied(HashMap<(bool, [i64; 4]), f64>),
}

impl PValues {
    fn pvalue(&self, test: Pp2Test, table: [[i64; 2]; 2]) -> Result<f64> {
        match self {
            PValues::Scipy => Ok(match test {
                Pp2Test::ChiSquare => stats::chi2_contingency_2x2(table).1,
                Pp2Test::Fisher => stats::fisher_exact_greater(table).1,
            }),
            PValues::Supplied(values) => {
                let key = (
                    test == Pp2Test::ChiSquare,
                    [table[0][0], table[0][1], table[1][0], table[1][1]],
                );
                values.get(&key).copied().ok_or_else(|| {
                    invalid(format!("no p-value supplied for {test:?} on {table:?}"))
                })
            }
        }
    }
}

/// A pickled set or frozenset.
enum SetView<'a> {
    Shared(std::cell::Ref<'a, Set>),
    Frozen(&'a Set),
}

impl Deref for SetView<'_> {
    type Target = Set;
    fn deref(&self) -> &Set {
        match self {
            SetView::Shared(set) => set,
            SetView::Frozen(set) => set,
        }
    }
}

fn set_view<'a>(value: &'a Value, what: &str) -> Result<SetView<'a>> {
    match value {
        Value::Set(set) => Ok(SetView::Shared(set.borrow())),
        Value::FrozenSet(set) => Ok(SetView::Frozen(set)),
        other => Err(type_error(format!(
            "{what} is a pickled {}; the set arithmetic of analyze_bp1_pp2 needs sets",
            other.type_name()
        ))),
    }
}

/// `len(left & right)` (`any_of` holds the sets of a union on the right).
fn intersection_count(left: &Set, any_of: &[&Set]) -> Result<usize> {
    let mut count = 0;
    for item in &left.items {
        let mut found = false;
        for set in any_of {
            if set.contains(item)? {
                found = true;
                break;
            }
        }
        count += usize::from(found);
    }
    Ok(count)
}

/// `len(left - right)`.
fn difference_count(left: &Set, right: &Set) -> Result<usize> {
    Ok(left.items.len() - intersection_count(left, &[right])?)
}

/// `analyze_bp1_pp2(gene_stat_dict)`: `(ensembl_id, (pp2_granted, bp1_granted))`; the
/// `ensembl_id` is `None` when it is not a text (no table gene can then find the result).
pub(crate) fn analyze_bp1_pp2(
    stat: &Value,
    pvalues: &PValues,
) -> Result<(Option<String>, (bool, bool))> {
    let Value::Dict(dict) = stat else {
        return Err(type_error(format!(
            "a gene's ClinVar statistics are a pickled {}, not a dict",
            stat.type_name()
        )));
    };
    let dict = dict.borrow();
    let field = |key: &str| -> Result<Value> {
        dict.get(&Value::Str(key.into()))?
            .cloned()
            .ok_or_else(|| key_error(key))
    };
    let gene = match field("ensembl_id")? {
        Value::Str(text) => Some(text.to_string()),
        _ => None,
    };
    if dict.entries.len() <= 1 {
        return Ok((gene, (false, false)));
    }
    let patho_value = field("clinvar_pathogenic")?;
    let patho = set_view(&patho_value, "clinvar_pathogenic")?;
    let total_patho = patho.items.len();
    if total_patho == 0 {
        return Ok((gene, (false, false)));
    }
    let large_value = field("large_aachange")?;
    let nmd_value = field("putative_nmd_variants")?;
    let large = set_view(&large_value, "large_aachange")?;
    let nmd = set_view(&nmd_value, "putative_nmd_variants")?;
    let patho_large_aa_nmd = intersection_count(&patho, &[&*large, &*nmd])?;
    let bp1_granted = patho_large_aa_nmd as f64 / total_patho as f64 >= 0.7;

    let small_value = field("small_aachange")?;
    let small = set_view(&small_value, "small_aachange")?;
    let patho_small = intersection_count(&patho, &[&*small])? as i64;
    let patho_large = intersection_count(&patho, &[&*large])? as i64;
    let non_patho_small = difference_count(&small, &patho)? as i64;
    let non_patho_large = difference_count(&large, &patho)? as i64;
    let table = [
        [patho_small, patho_large],
        [non_patho_small, non_patho_large],
    ];
    let total_count = patho_small + patho_large + non_patho_small + non_patho_large;
    let min_cell = patho_small
        .min(patho_large)
        .min(non_patho_small)
        .min(non_patho_large);
    let pp2_fraction = patho_small as f64 / total_patho as f64;
    let pp2_granted = if total_count >= 5
        && table.iter().all(|row| row[0] + row[1] > 0)
        && (0..2).all(|column| table[0][column] + table[1][column] > 0)
    {
        let test = if total_count >= 20 && min_cell >= 5 {
            Pp2Test::ChiSquare
        } else {
            Pp2Test::Fisher
        };
        let pvalue = pvalues.pvalue(test, table)?;
        // `min(pvalue, 0.99999)` keeps a NaN p-value.
        let pvalue = if 0.99999 < pvalue { 0.99999 } else { pvalue };
        pvalue > 0.05 && pp2_fraction >= 0.4
    } else {
        pp2_fraction >= 0.5
    };
    Ok((gene, (pp2_granted, bp1_granted)))
}

/// The table's genes in `df['Gene'].dropna().unique()` order, as dict keys (text only: no
/// other kind of value can match the text keys of the statistics and of the results).
fn gene_results(
    table: &Table,
    clinvar_stats: &Value,
    pvalues: &PValues,
) -> Result<HashMap<String, (bool, bool)>> {
    let genes = table.column("Gene")?;
    let Value::Dict(stats_dict) = clinvar_stats else {
        return Err(invalid("the ClinVar gene statistics pickle is not a dict"));
    };
    let mut results: HashMap<String, (bool, bool)> = HashMap::new();
    let mut seen = std::collections::HashSet::new();
    for row in 0..table.rows {
        let Cell::Text(gene) = genes.cell(row) else {
            continue;
        };
        if !seen.insert(gene) {
            continue;
        }
        let stat = stats_dict.borrow().get(&Value::Str(gene.into()))?.cloned();
        let (key, result) = match stat {
            Some(stat) => analyze_bp1_pp2(&stat, pvalues)?,
            None => (Some(gene.to_owned()), (false, false)),
        };
        let Some(key) = key else {
            continue;
        };
        if let Some(previous) = results.insert(key.clone(), result)
            && previous != result
        {
            return Err(invalid(format!(
                "two genes' ClinVar statistics name ensembl_id {key} with different PP2/BP1 \
                 results; Python keeps whichever finishes last"
            )));
        }
    }
    Ok(results)
}

/// `PP2_BP1_criteria(df, clinvar_stats_dict, am_intol_domains_tsv)`: PP2 and BP1 (strength 1).
pub(crate) fn pp2_bp1_criteria(
    table: &Table,
    clinvar_stats: &Value,
    am_intol_domains_tsv: &str,
    pvalues: &PValues,
) -> Result<StepOutput> {
    let results = gene_results(table, clinvar_stats, pvalues)?;
    let genes = table.column("Gene")?;
    let gene_result = |row: usize| match genes.cell(row) {
        Cell::Text(gene) => results.get(gene).copied().unwrap_or((false, false)),
        _ => (false, false),
    };
    let gene_pp2 = Mask::Plain((0..table.rows).map(|row| gene_result(row).0).collect());
    let gene_bp1 = Mask::Plain((0..table.rows).map(|row| gene_result(row).1).collect());
    // `str.contains('missense_variant', na=False)`.
    let is_missense = text_contains(table, "Consequence", None, &["missense_variant"])?
        .remove(0)
        .fillna(false);
    let pp2 = gene_pp2.and(&is_missense);
    let bp1 = gene_bp1.and(&is_missense);

    let domains = read_table(Path::new(am_intol_domains_tsv))?;
    let domain_texts = text_values(&domains, "domain", None)?;
    let tolerant = domains.column("is_more_tolerant")?;
    if matches!(
        tolerant,
        Column::Object(_) | Column::NullableBool(_) | Column::NullableInt(_)
    ) {
        return Err(invalid(
            "is_more_tolerant is a mixed or nullable column; .any() is not reproduced there",
        ));
    }
    // `groupby(Ensembl_Gene_ID)` then `is_more_tolerant.any()` per gene (missing values skipped).
    let mut any_tolerant: HashMap<String, bool> = HashMap::new();
    for (row, domain) in domain_texts.iter().enumerate() {
        let Some(domain) = domain else {
            continue;
        };
        let gene = domain.split(':').next().unwrap_or("").to_owned();
        let truthy = match tolerant.cell(row) {
            Cell::Missing => false,
            Cell::Bool(flag) => flag,
            Cell::Int(number) => number != 0,
            Cell::Float(number) => number != 0.0,
            Cell::Text(text) => !text.is_empty(),
        };
        *any_tolerant.entry(gene).or_insert(false) |= truthy;
    }
    let pp2_am = Mask::Plain(
        (0..table.rows)
            .map(|row| match genes.cell(row) {
                Cell::Text(gene) => any_tolerant.get(gene).is_some_and(|tolerant| !tolerant),
                _ => false,
            })
            .collect(),
    );
    let gene_am = numeric_column(table, "Gene_avg_AM_score")?;
    let pp2_gene_am = Mask::Plain(
        gene_am
            .iter()
            .map(|value| if value.is_nan() { 0.0 } else { *value } > 0.564)
            .collect(),
    );
    let pp2 = pp2.or(&pp2_am).or(&pp2_gene_am).and(&is_missense);
    Ok(StepOutput {
        criteria: vec![
            ("PP2", Column::Int(strengths(&pp2, 1)?)),
            ("BP1", Column::Int(strengths(&bp1, 1)?)),
        ],
        ..StepOutput::default()
    })
}

/// The harness's `pvalues` option: `[[test, a, b, c, d, p_bits], ...]` with `test` "chi_square"
/// or "fisher" and the p-value as the hex text of its 64 bits.
fn supplied_pvalues(options: &crate::json_text::Json) -> Result<PValues> {
    use crate::json_text::Json;
    let Some(Json::Array(entries)) = options.get("pvalues") else {
        return Ok(PValues::Scipy);
    };
    let mut values = HashMap::new();
    for entry in entries {
        let Json::Array(items) = entry else {
            return Err(invalid("each p-value entry is a list"));
        };
        let (Some(Json::Str(test)), Some(Json::Str(bits))) = (items.first(), items.get(5)) else {
            return Err(invalid("a p-value entry is [test, a, b, c, d, p bits]"));
        };
        let mut cells = [0i64; 4];
        for (position, cell) in cells.iter_mut().enumerate() {
            let Some(Json::Int(number)) = items.get(position + 1) else {
                return Err(invalid("table cells must be integers"));
            };
            *cell = i64::try_from(*number).map_err(|_| invalid("table cell beyond 64 bits"))?;
        }
        let pvalue = f64::from_bits(
            u64::from_str_radix(bits, 16).map_err(|_| invalid("p-value bits must be hex"))?,
        );
        values.insert((test == "chi_square", cells), pvalue);
    }
    Ok(PValues::Supplied(values))
}

/// `priva acmg-step acmg_pm1_pp2_bp1_missense <function>`:
///
/// - `PM1_criteria`, options `pvs1_criteria` (ints), `loc_intol_domain` (booleans),
///   `pm1_regions_pkl`: a frame of `PM1` and `loc_intol_domain`;
/// - `PP2_BP1_criteria`, options `clinvar_gene_stat_pkl`, `am_intol_domains_tsv`, `pvalues`
///   (optional): a frame of `PP2`, `BP1`;
/// - `analyze_bp1_pp2`, options `clinvar_gene_stat_pkl`, `pvalues`: for each text of the input
///   column `Gene`, `analyze_bp1_pp2(stats.get(gene, {"ensembl_id": gene}))` as columns
///   `ensembl_id`, `pp2`, `bp1`.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let table = &call.table;
    let options = &call.options;
    match function {
        "PM1_criteria" => {
            let pvs1 = harness_options::int_list(options, "pvs1_criteria")?
                .ok_or_else(|| invalid("option pvs1_criteria is required"))?;
            let loc = harness_options::bool_list(options, "loc_intol_domain")?
                .ok_or_else(|| invalid("option loc_intol_domain is required"))?;
            let path = harness_options::required_text(options, "pm1_regions_pkl")?;
            let output = pm1_criteria(table, &pvs1, &loc, path)?;
            let mut columns: Vec<(String, Column)> = output
                .criteria
                .iter()
                .map(|(name, values)| ((*name).to_owned(), values.clone()))
                .collect();
            columns.extend(output.passed_on.iter().cloned());
            Ok(harness_options::frame(columns, table.rows))
        }
        "PP2_BP1_criteria" => {
            let stats = load_pickle(harness_options::required_text(
                options,
                "clinvar_gene_stat_pkl",
            )?)?;
            let tsv = harness_options::required_text(options, "am_intol_domains_tsv")?;
            let output = pp2_bp1_criteria(table, &stats, tsv, &supplied_pvalues(options)?)?;
            Ok(harness_options::criteria_frame(&output, table.rows))
        }
        "analyze_bp1_pp2" => {
            let stats = load_pickle(harness_options::required_text(
                options,
                "clinvar_gene_stat_pkl",
            )?)?;
            let pvalues = supplied_pvalues(options)?;
            let genes = table.column("Gene")?;
            let mut ids = Vec::new();
            let mut pp2 = Vec::new();
            let mut bp1 = Vec::new();
            for row in 0..table.rows {
                let Cell::Text(gene) = genes.cell(row) else {
                    return Err(invalid("analyze_bp1_pp2 cases name text genes"));
                };
                let stat = dict_get(&stats, gene)?;
                let (id, (granted_pp2, granted_bp1)) = match stat {
                    Some(stat) => analyze_bp1_pp2(&stat, &pvalues)?,
                    None => (Some(gene.to_owned()), (false, false)),
                };
                ids.push(id);
                pp2.push(granted_pp2);
                bp1.push(granted_bp1);
            }
            Ok(harness_options::frame(
                vec![
                    ("ensembl_id".into(), Column::Text(ids)),
                    ("pp2".into(), Column::Bool(pp2)),
                    ("bp1".into(), Column::Bool(bp1)),
                ],
                table.rows,
            ))
        }
        _ => Err(not_converted("acmg_pm1_pp2_bp1_missense", function)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn protein_positions_follow_python_isdigit() {
        assert_eq!(
            protein_position(Some(Cell::Text("12/300"))).unwrap(),
            Some(12)
        );
        assert_eq!(protein_position(Some(Cell::Text("007"))).unwrap(), Some(7));
        assert_eq!(
            protein_position(Some(Cell::Text("12-14/300"))).unwrap(),
            None
        );
        assert_eq!(protein_position(Some(Cell::Text("?-5/300"))).unwrap(), None);
        assert_eq!(protein_position(Some(Cell::Text("-"))).unwrap(), None);
        assert_eq!(protein_position(Some(Cell::Int(5))).unwrap(), Some(5));
        assert_eq!(protein_position(Some(Cell::Int(0))).unwrap(), None);
        assert_eq!(protein_position(Some(Cell::Float(5.0))).unwrap(), None);
        assert_eq!(protein_position(None).unwrap(), None);
        assert!(protein_position(Some(Cell::Text("\u{0661}2"))).is_err());
    }

    fn text_op(text: &str) -> Vec<u8> {
        let mut bytes = vec![0x8c, text.len() as u8];
        bytes.extend_from_slice(text.as_bytes());
        bytes
    }

    fn set_op(items: &[&str]) -> Vec<u8> {
        let mut bytes = vec![0x8f];
        if !items.is_empty() {
            bytes.push(0x28);
            for item in items {
                bytes.extend(text_op(item));
            }
            bytes.push(0x90);
        }
        bytes
    }

    #[test]
    fn pp2_fraction_branch_needs_no_statistics() {
        let mut bytes = vec![0x80, 0x04, 0x7d, 0x28];
        for (key, value) in [
            ("ensembl_id", text_op("G")),
            ("clinvar_pathogenic", set_op(&["a", "b"])),
            ("large_aachange", set_op(&["a"])),
            ("putative_nmd_variants", set_op(&[])),
            ("small_aachange", set_op(&["b"])),
        ] {
            bytes.extend(text_op(key));
            bytes.extend(value);
        }
        bytes.extend([0x75, 0x2e]);
        let stat = pickle_reader::load(&bytes).unwrap();
        // 2 pathogenic: 1 small, 1 large; under five counted, so the fraction decides PP2.
        let result = analyze_bp1_pp2(&stat, &PValues::Supplied(HashMap::new())).unwrap();
        assert_eq!(result, (Some("G".to_owned()), (true, false)));
    }
}
