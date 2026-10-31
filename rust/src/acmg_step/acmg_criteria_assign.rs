//! `scripts/acmg_criteria_assign.py` (711 lines on main): the ACMG step's orchestrator, as
//! `priva acmg-assign` with the Python option names.
//!
//! - The command line is the script's argparse parser: the same options (with `--pp1_vcf` and
//!   `--pp1_ped`, forwarded to PP1, and the branch's `--ranking_binary`, accepted and not needed
//!   once the step is Rust), defaults, the `store_true` flag `--ignore_clingen_override`, the
//!   `type=float` and `type=int` conversions of every occurrence, and exit status 2 for a usage
//!   error. A number written with non-ASCII digits or spaces, which Python's `float()` and
//!   `int()` may read, is refused with exit 1.
//! - `run_step` is `ACMG_criteria_assign` in the script's order: the table read and
//!   de-duplicated, `vep_consq_interpret`, the AlphaMissense mean scores, the pedigree (with the
//!   headerless fallback) and its proband, `Gene_avg_AM_score` and `variant_id`, the exact-allele
//!   columns, the hub, then the criteria, the PM2/BS1 audit columns, the criteria matrix
//!   (`<table without its last suffix>.acmg.tsv`, written before the posterior), the posterior,
//!   the ranking, and the table written over the input. A failure after the matrix is written
//!   leaves the matrix and the untouched input table, as in Python.
//! - The hub of thread T2, the criteria of threads T3, T6 and T7 and T3's ranking are called
//!   directly. Those of threads still converting (T4 PVS1/PS1/PM5, T5 PM2/BS1/BA1/BS2/BP5) go
//!   through the
//!   `pending` module below, one function each with the inputs the Python call takes; until
//!   their thread lands they stop the step with "not converted". The comparison harness can
//!   instead hand them the outputs a recording Python run wrote (`Source::Recorded`), so the
//!   converted parts are compared on the real table now. The ranking's control-common step
//!   needs T5's `stats` (`binom_cdf`, `beta_ppf`); in a recorded comparison its two columns
//!   come from the recording too.
//! - `read_ped`, `proband`, `am_score_dict`, `transcript_to_gene_map`, `gene_avg_am_score`,
//!   `variant_id` and `frequency_audit` are the script's inline code; refused where pandas'
//!   behaviour depends on odd column types (noted at each).
use super::acmg_consequence::vep_consq_interpret;
use super::acmg_pm1_pp2_bp1_missense::{load_pickle, pm1_criteria, pp2_bp1_criteria};
use super::acmg_pm4_bp3_protein_length::pm4_bp3_criteria;
use super::acmg_pp1_bs4_bp2_pm3_family::{bp2_pm3_criteria, bs4_criteria, pp1_criteria};
use super::acmg_pp3_bp4_bp7_insilico::pandas_masks::python_str;
use super::acmg_pp3_bp4_bp7_insilico::{bp7_criteria, pp3_bp4_criteria};
use super::acmg_ps2_pm6_denovo::ps2_pm6_criteria;
use super::acmg_ps3_bs3_pp5_bp6_reported::{pp5_bp6_criteria, ps3_bs3_criteria};
use super::acmg_scoring::{
    RankingOptions, posterior_columns, sort_and_rank_variants_with, summarize_acmg_criteria,
};
use super::acmg_variant_mechanism::annotate_exact_nonlof_variants;
use super::gene_mechanism_annotation::annotate_gene_mechanism_categories;
use super::gene_mechanism_common::clean;
use super::gene_mechanism_resources::{GeneMechanismHub, HubPaths};
use super::types::{
    AcmgError, HarnessCall, MATRIX_AUDIT_COLUMNS, StepArguments, StepOutput, SymbolResolver,
    not_converted,
};
use crate::cli::{Result, invalid};
use crate::json_text::Json;
use crate::pandas_polars;
use crate::pandas_table::{
    Cell, Column, Header, ObjectValue, ReadOptions, Table, UseColumns, read_table, read_table_with,
    write_table,
};
use crate::pickle_reader;
use crate::protein_domain_mapping::data_directory;
use crate::python_args::{self, parse_long_option_occurrences};
use crate::python_text::{parse_float, python_int};
use crate::spliceai_cutoffs::{Cutoffs, cutoffs_for_config};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

// ------------------------------------------------------------------------------ command line

/// How an option is converted (`type=str`, `type=float`, `type=int`, `action="store_true"`).
#[derive(Clone, Copy, PartialEq)]
enum Kind {
    Text,
    Float,
    Int,
    Flag,
}

/// The parser's options in the order the script adds them (argparse lists ambiguous
/// prefixes and missing required options in this order); `true` marks a required one.
const OPTIONS: [(&str, Kind, bool); 37] = [
    ("anno_table", Kind::Text, true),
    ("am_score_table", Kind::Text, true),
    ("ped_table", Kind::Text, false),
    ("fam_name", Kind::Text, false),
    ("clinvar_patho_af_stat", Kind::Text, true),
    ("clinvar_patho_exon_af_stat", Kind::Text, true),
    ("clinvar_aa_dict_pkl", Kind::Text, true),
    ("clinvar_splice_dict_pkl", Kind::Text, true),
    ("interpro_entry_map_pkl", Kind::Text, true),
    ("intolerant_domains_pkl", Kind::Text, true),
    ("gene_dosage_sensitivity", Kind::Text, true),
    ("am_intol_domains_tsv", Kind::Text, true),
    ("pm1_regions_pkl", Kind::Text, true),
    ("clinvar_gene_stat_pkl", Kind::Text, true),
    ("tranx_exon_domain_map_pkl", Kind::Text, true),
    ("am_score_vcf", Kind::Text, false),
    ("alt_disease_vcf", Kind::Text, false),
    ("clingen_map_pkl", Kind::Text, false),
    ("ignore_clingen_override", Kind::Flag, false),
    ("mavedb_metadata_tsv", Kind::Text, false),
    ("pp1_vcf", Kind::Text, false),
    ("pp1_ped", Kind::Text, false),
    ("cds_fasta_path", Kind::Text, false),
    ("pext_tissues", Kind::Text, false),
    ("pext_low_expression_cutoff", Kind::Float, false),
    ("pext_penalty_floor", Kind::Float, false),
    ("pext_penalty_shape", Kind::Float, false),
    ("repeat_region_file", Kind::Text, true),
    ("gnomAD_extreme_rare_threshold", Kind::Float, false),
    ("expected_incidence", Kind::Float, false),
    ("relevant_gene_list", Kind::Text, false),
    ("dispensable_gene_list", Kind::Text, false),
    ("hpo_condition_mechanism_json", Kind::Text, false),
    ("gene_nonlof_mechanism_json", Kind::Text, false),
    ("threads", Kind::Int, false),
    ("config", Kind::Text, false),
    ("ranking_binary", Kind::Text, false),
];

/// The argparse usage line under the program name `priva acmg-assign`.
fn usage() -> String {
    let mut parts = vec!["usage: priva acmg-assign [-h]".to_owned()];
    for (name, kind, required) in OPTIONS {
        let option = match kind {
            Kind::Flag => format!("--{name}"),
            _ => format!("--{name} {}", name.to_uppercase()),
        };
        parts.push(if required {
            option
        } else {
            format!("[{option}]")
        });
    }
    parts.join(" ")
}

/// Why the command line is not run.
#[derive(Debug, PartialEq)]
pub(crate) enum ArgumentError {
    /// `-h`/`--help`: the usage is printed and the program ends with status 0.
    Help,
    /// An argparse error: usage and message on stderr, exit status 2.
    Usage(String),
    /// A number Python may read but this parser does not: exit status 1.
    Refused(String),
}

/// A value of an option after argparse's conversion.
#[derive(Clone, Debug, PartialEq)]
enum Converted {
    Text(String),
    Float(f64),
    Int(i64),
    Flag,
}

fn convert(name: &str, kind: Kind, value: &str) -> std::result::Result<Converted, ArgumentError> {
    let refuse = |what: &str| {
        ArgumentError::Refused(format!(
            "refused (not reproduced): argument --{name}: {value:?} holds a non-ASCII \
             character, which Python's {what}() may read as a number"
        ))
    };
    match kind {
        Kind::Text => Ok(Converted::Text(value.to_owned())),
        Kind::Flag => Ok(Converted::Flag),
        Kind::Float => match parse_float(value) {
            Some(number) => Ok(Converted::Float(number)),
            None if !value.is_ascii() => Err(refuse("float")),
            None => Err(ArgumentError::Usage(format!(
                "argument --{name}: invalid float value: '{value}'"
            ))),
        },
        Kind::Int => match python_int(value) {
            Some(number) => i64::try_from(number).map(Converted::Int).map_err(|_| {
                ArgumentError::Refused(format!(
                    "refused (not reproduced): argument --{name}: {value} is beyond 64 bits"
                ))
            }),
            None if !value.is_ascii() => Err(refuse("int")),
            None => Err(ArgumentError::Usage(format!(
                "argument --{name}: invalid int value: '{value}'"
            ))),
        },
    }
}

/// A default under the repository's `data/` folder (`PRIVA_DATA_DIR`, else the source tree's).
fn data_file(parts: &[&str]) -> String {
    let mut path: PathBuf = data_directory();
    for part in parts {
        path.push(part);
    }
    path.to_string_lossy().into_owned()
}

/// `parser.parse_args(args)` of the script.
pub(crate) fn parse_arguments(
    args: &[String],
) -> std::result::Result<StepArguments, ArgumentError> {
    let names: Vec<&str> = OPTIONS
        .iter()
        .filter(|(_, kind, _)| *kind != Kind::Flag)
        .map(|(name, _, _)| *name)
        .collect();
    let flags: Vec<&str> = OPTIONS
        .iter()
        .filter(|(_, kind, _)| *kind == Kind::Flag)
        .map(|(name, _, _)| *name)
        .collect();
    let occurrences = match parse_long_option_occurrences(args, &names, &flags) {
        Ok(Some(occurrences)) => occurrences,
        Ok(None) => return Err(ArgumentError::Help),
        Err(message) => return Err(ArgumentError::Usage(message)),
    };
    // argparse converts every occurrence in command-line order; the last one is kept.
    let mut values: HashMap<&str, Converted> = HashMap::new();
    let mut refusal = None;
    for (name, value) in &occurrences {
        let (option, kind, _) = OPTIONS
            .iter()
            .find(|(option, _, _)| option == name)
            .expect("the parser only returns known options");
        match convert(option, *kind, value) {
            Ok(converted) => {
                values.insert(option, converted);
            }
            Err(ArgumentError::Refused(message)) => {
                refusal.get_or_insert(message);
            }
            Err(error) => return Err(error),
        }
    }
    let missing: Vec<&str> = OPTIONS
        .iter()
        .filter(|(name, _, required)| *required && !values.contains_key(name))
        .map(|(name, _, _)| *name)
        .collect();
    if !missing.is_empty() {
        let given: HashMap<String, String> = HashMap::new();
        let message =
            python_args::require_options(&given, &missing).expect_err("the options are missing");
        return Err(ArgumentError::Usage(message));
    }
    if let Some(message) = refusal {
        return Err(ArgumentError::Refused(message));
    }
    let text = |name: &str| match values.get(name) {
        Some(Converted::Text(value)) => Some(value.clone()),
        _ => None,
    };
    let float = |name: &str, default: f64| match values.get(name) {
        Some(Converted::Float(value)) => *value,
        _ => default,
    };
    let required = |name: &str| text(name).expect("checked above");
    Ok(StepArguments {
        anno_table: required("anno_table"),
        am_score_table: required("am_score_table"),
        ped_table: text("ped_table"),
        fam_name: text("fam_name"),
        clinvar_patho_af_stat: required("clinvar_patho_af_stat"),
        clinvar_patho_exon_af_stat: required("clinvar_patho_exon_af_stat"),
        clinvar_aa_dict_pkl: required("clinvar_aa_dict_pkl"),
        clinvar_splice_dict_pkl: required("clinvar_splice_dict_pkl"),
        interpro_entry_map_pkl: required("interpro_entry_map_pkl"),
        intolerant_domains_pkl: required("intolerant_domains_pkl"),
        gene_dosage_sensitivity: required("gene_dosage_sensitivity"),
        am_intol_domains_tsv: required("am_intol_domains_tsv"),
        pm1_regions_pkl: required("pm1_regions_pkl"),
        clinvar_gene_stat_pkl: required("clinvar_gene_stat_pkl"),
        tranx_exon_domain_map_pkl: required("tranx_exon_domain_map_pkl"),
        am_score_vcf: text("am_score_vcf"),
        alt_disease_vcf: text("alt_disease_vcf"),
        clingen_map_pkl: text("clingen_map_pkl"),
        ignore_clingen_override: values.contains_key("ignore_clingen_override"),
        mavedb_metadata_tsv: text("mavedb_metadata_tsv"),
        pp1_vcf: text("pp1_vcf").unwrap_or_default(),
        pp1_ped: text("pp1_ped").unwrap_or_default(),
        cds_fasta_path: text("cds_fasta_path"),
        pext_tissues: text("pext_tissues").unwrap_or_default(),
        pext_low_expression_cutoff: float("pext_low_expression_cutoff", 0.1),
        pext_penalty_floor: float("pext_penalty_floor", 0.8),
        pext_penalty_shape: float("pext_penalty_shape", 0.5),
        repeat_region_file: required("repeat_region_file"),
        gnomad_extreme_rare_threshold: float("gnomAD_extreme_rare_threshold", 0.0001),
        expected_incidence: float("expected_incidence", 0.0001),
        relevant_gene_list: text("relevant_gene_list"),
        dispensable_gene_list: Some(
            text("dispensable_gene_list")
                .unwrap_or_else(|| data_file(&["dispensable_genes", "dispensable_gene_list.txt"])),
        ),
        hpo_condition_mechanism_json: text("hpo_condition_mechanism_json").unwrap_or_else(|| {
            data_file(&["patho_mechanism", "hpo_condition_mechanism_cache.json.gz"])
        }),
        gene_nonlof_mechanism_json: text("gene_nonlof_mechanism_json").unwrap_or_else(|| {
            data_file(&[
                "patho_mechanism",
                "gene_nonlof_mechanism_curated_assertions.json.gz",
            ])
        }),
        threads: match values.get("threads") {
            Some(Converted::Int(value)) => *value,
            _ => 10,
        },
        ranking_binary: text("ranking_binary"),
        config: text("config"),
    })
}

// ------------------------------------------------------------------------------ the readers

const PED_NAMES: [&str; 6] = [
    "#FamilyID",
    "IndividualID",
    "PaternalID",
    "MaternalID",
    "Sex",
    "Phenotype",
];

/// The pedigree as the step reads it: `read_table(dtype=str, keep_default_na=False)`, and
/// when that has no `#FamilyID` column, again without a header, the first six fields named
/// `#FamilyID IndividualID PaternalID MaternalID Sex Phenotype`.
pub(crate) fn read_ped(path: &str) -> Result<Table> {
    let text_options = ReadOptions {
        all_text: true,
        keep_default_na: false,
        ..ReadOptions::default()
    };
    let first = read_table_with(Path::new(path), &text_options)?;
    if first.has_column("#FamilyID") {
        return Ok(first);
    }
    read_table_with(
        Path::new(path),
        &ReadOptions {
            header: Header::Names(PED_NAMES.iter().map(|name| (*name).to_owned()).collect()),
            use_columns: Some(UseColumns::Positions((0..6).collect())),
            ..text_options
        },
    )
}

/// The family's first affected member (`Phenotype == "2"`), or Python's ValueError.
pub(crate) fn proband(ped: &Table, fam_name: Option<&str>) -> Result<String> {
    let families = ped.column("#FamilyID")?;
    let phenotypes = ped.column("Phenotype")?;
    let individuals = ped.column("IndividualID")?;
    let found = (0..ped.rows).find(|row| {
        fam_name.is_some_and(|name| families.cell(*row) == Cell::Text(name))
            && phenotypes.cell(*row) == Cell::Text("2")
    });
    match found.map(|row| individuals.cell(row)) {
        Some(Cell::Text(individual)) => Ok(individual.to_owned()),
        Some(_) => Err(invalid(
            "refused (not reproduced): a pedigree read as text has a value that is not text",
        )),
        None => Err(Box::new(AcmgError::new(
            "ValueError",
            format!(
                "No affected individual found for family {}",
                fam_name.unwrap_or("None")
            ),
        ))),
    }
}

/// The text values of a column that must be text (missing values allowed or not).
fn text_column<'a>(
    table: &'a Table,
    name: &str,
    missing_allowed: bool,
) -> Result<Vec<Option<&'a str>>> {
    let column = table.column(name)?;
    let Column::Text(values) = column else {
        return Err(invalid(format!(
            "refused (not reproduced): {name} is a pandas {} column of other values than text",
            column.dtype_name()
        )));
    };
    if !missing_allowed && values.iter().any(Option::is_none) {
        return Err(invalid(format!(
            "refused (not reproduced): {name} has a missing value (Python keys a dict by NaN)"
        )));
    }
    Ok(values.iter().map(Option::as_deref).collect())
}

/// `dict(zip(am_score_df['transcript'], am_score_df['mean_am_pathogenicity']))` (the last
/// duplicate wins). Transcripts must be text and the scores floats.
pub(crate) fn am_score_dict(am_scores: &Table) -> Result<HashMap<String, f64>> {
    if am_scores.rows == 0 {
        am_scores.column("transcript")?;
        am_scores.column("mean_am_pathogenicity")?;
        return Ok(HashMap::new());
    }
    let transcripts = text_column(am_scores, "transcript", false)?;
    let scores = am_scores.column("mean_am_pathogenicity")?;
    let Column::Float(scores) = scores else {
        return Err(invalid(format!(
            "refused (not reproduced): mean_am_pathogenicity is a pandas {} column",
            scores.dtype_name()
        )));
    };
    Ok(transcripts
        .into_iter()
        .zip(scores)
        .map(|(transcript, score)| (transcript.unwrap_or_default().to_owned(), *score))
        .collect())
}

/// `dict(zip(anno_df['Feature'], anno_df['Gene']))`: each transcript in the order of its
/// first row, with the gene of its last row; `None` is NaN (one key, as pandas' NaN object).
pub(crate) fn transcript_to_gene_map(
    table: &Table,
) -> Result<Vec<(Option<String>, Option<String>)>> {
    let features = text_column(table, "Feature", true)?;
    let genes = text_column(table, "Gene", true)?;
    let mut order: Vec<(Option<String>, Option<String>)> = Vec::new();
    let mut place: HashMap<Option<&str>, usize> = HashMap::new();
    for (feature, gene) in features.into_iter().zip(genes) {
        match place.get(&feature) {
            Some(index) => order[*index].1 = gene.map(str::to_owned),
            None => {
                place.insert(feature, order.len());
                order.push((feature.map(str::to_owned), gene.map(str::to_owned)));
            }
        }
    }
    Ok(order)
}

/// `{g: am_score_dict[t] for t, g in transcript_to_gene_map.items() if t in am_score_dict}`.
pub(crate) fn gene_to_am_score_map(
    transcript_to_gene: &[(Option<String>, Option<String>)],
    am_scores: &HashMap<String, f64>,
) -> HashMap<Option<String>, f64> {
    let mut map = HashMap::new();
    for (transcript, gene) in transcript_to_gene {
        if let Some(score) = transcript.as_ref().and_then(|t| am_scores.get(t)) {
            map.insert(gene.clone(), *score);
        }
    }
    map
}

/// `anno_df["Gene"].map(gene_to_am_score_map)`: a float column, NaN where the gene has no
/// score (a NaN gene finds a NaN key, as in pandas' index).
pub(crate) fn gene_avg_am_score(
    table: &Table,
    map: &HashMap<Option<String>, f64>,
) -> Result<Column> {
    let genes = text_column(table, "Gene", true)?;
    Ok(Column::Float(
        genes
            .into_iter()
            .map(|gene| {
                map.get(&gene.map(str::to_owned))
                    .copied()
                    .unwrap_or(f64::NAN)
            })
            .collect(),
    ))
}

/// `anno_df["chrom"] + ":" + anno_df["pos"].astype(str) + ":" + anno_df["ref"] + "-" +
/// anno_df["alt"]`: missing where `chrom`, `ref` or `alt` is missing (pandas skips missing
/// values in `object` arithmetic; a `ref` or `alt` column of only missing values, read as
/// float64, makes every id missing). A `chrom` column that is not text, or a `ref` or `alt`
/// column of numbers, fails with Python's TypeError.
pub(crate) fn variant_id(table: &Table) -> Result<Column> {
    let type_error = |name: &str, column: &Column| -> Box<dyn std::error::Error> {
        Box::new(AcmgError::new(
            "TypeError",
            format!(
                "{name} is a pandas {} column: + with text fails",
                column.dtype_name()
            ),
        ))
    };
    // `+` on an object column of text and NaN, or a column of only missing values.
    let parts = |name: &str| -> Result<Vec<Option<&str>>> {
        let column = table.column(name)?;
        match column {
            Column::Text(values) => Ok(values.iter().map(Option::as_deref).collect()),
            Column::Float(values)
                if name != "chrom" && values.iter().all(|value| value.is_nan()) =>
            {
                Ok(vec![None; values.len()])
            }
            Column::Int(_) | Column::Float(_) | Column::Bool(_) => Err(type_error(name, column)),
            other => Err(invalid(format!(
                "refused (not reproduced): {name} is a pandas {} column of other values than text",
                other.dtype_name()
            ))),
        }
    };
    let chrom = parts("chrom")?;
    let pos = table.column("pos")?;
    if matches!(
        pos,
        Column::Object(_) | Column::NullableInt(_) | Column::NullableBool(_)
    ) {
        return Err(invalid(format!(
            "refused (not reproduced): pos is a pandas {} column",
            pos.dtype_name()
        )));
    }
    let reference = parts("ref")?;
    let alternate = parts("alt")?;
    Ok(Column::Text(
        (0..table.rows)
            .map(|row| match (chrom[row], reference[row], alternate[row]) {
                (Some(chrom), Some(reference), Some(alternate)) => Some(format!(
                    "{chrom}:{}:{reference}-{alternate}",
                    python_str(pos.cell(row))
                )),
                _ => None,
            })
            .collect(),
    ))
}

/// The integer values of a criterion array (numpy int64).
fn int_values(column: &Column, what: &str) -> Result<Vec<i64>> {
    match column {
        Column::Int(values) => Ok(values.clone()),
        other => Err(invalid(format!(
            "refused (not reproduced): {what} is a {} array where integers are expected",
            other.dtype_name()
        ))),
    }
}

/// The values of a boolean array (numpy bool).
fn bool_values(column: &Column, what: &str) -> Result<Vec<bool>> {
    match column {
        Column::Bool(values) => Ok(values.clone()),
        other => Err(invalid(format!(
            "refused (not reproduced): {what} is a {} array where booleans are expected",
            other.dtype_name()
        ))),
    }
}

/// `array > 0` of a criterion (an int64 or bool array).
fn positive(column: &Column, what: &str) -> Result<Vec<bool>> {
    match column {
        Column::Int(values) => Ok(values.iter().map(|value| *value > 0).collect()),
        Column::Bool(values) => Ok(values.clone()),
        other => Err(invalid(format!(
            "refused (not reproduced): {what} is a {} array",
            other.dtype_name()
        ))),
    }
}

/// The PM2 back-down for the disease-incidence BS1 route and the five audit columns
/// (script lines 412-431): `(PM2 after the back-down, [(audit column, 0/1 values)])`.
pub(crate) fn frequency_audit(
    pm2: &[i64],
    bs1: &Column,
    disease_incidence_bs1: &[bool],
    gene_pathogenic_af_bs1: &[bool],
    gene_pathogenic_af_blocked_by_pm2: &[bool],
) -> Result<(Vec<i64>, Vec<(&'static str, Column)>)> {
    let rows = pm2.len();
    if [
        bs1.len(),
        disease_incidence_bs1.len(),
        gene_pathogenic_af_bs1.len(),
        gene_pathogenic_af_blocked_by_pm2.len(),
    ]
    .iter()
    .any(|length| *length != rows)
    {
        return Err(invalid(
            "refused (not reproduced): the PM2 and BS1 arrays differ in length",
        ));
    }
    let backdown: Vec<bool> = (0..rows)
        .map(|row| pm2[row] > 0 && disease_incidence_bs1[row])
        .collect();
    let pm2_after: Vec<i64> = (0..rows)
        .map(|row| if backdown[row] { 0 } else { pm2[row] })
        .collect();
    let bs1_positive = positive(bs1, "BS1")?;
    let conflicts: Vec<bool> = (0..rows)
        .map(|row| pm2_after[row] > 0 && bs1_positive[row])
        .collect();
    let as_int = |flags: &[bool]| Column::Int(flags.iter().map(|flag| i64::from(*flag)).collect());
    let audit = vec![
        (MATRIX_AUDIT_COLUMNS[0], as_int(disease_incidence_bs1)),
        (MATRIX_AUDIT_COLUMNS[1], as_int(gene_pathogenic_af_bs1)),
        (
            MATRIX_AUDIT_COLUMNS[2],
            as_int(gene_pathogenic_af_blocked_by_pm2),
        ),
        (MATRIX_AUDIT_COLUMNS[3], as_int(&backdown)),
        (MATRIX_AUDIT_COLUMNS[4], as_int(&conflicts)),
    ];
    Ok((pm2_after, audit))
}

/// `".".join(anno_table.split(".")[:-1]) + ".acmg.tsv"`: the whole path text is split, so a
/// dot in a folder name counts too, and a name without a dot gives `.acmg.tsv`.
pub(crate) fn matrix_path(anno_table: &str) -> String {
    let parts: Vec<&str> = anno_table.split('.').collect();
    format!("{}.acmg.tsv", parts[..parts.len() - 1].join("."))
}

/// Writes a criterion function's table columns in its order (`df[name] = values`).
fn apply_columns(table: &mut Table, output: &StepOutput) -> Result<()> {
    for (name, column) in &output.columns {
        if column.len() != table.rows {
            return Err(invalid(format!(
                "the column {name} has {} values for {} rows",
                column.len(),
                table.rows
            )));
        }
        table.set_column(name, column.clone());
    }
    Ok(())
}

fn criterion(output: &StepOutput, name: &str) -> Result<Column> {
    output
        .criterion(name)
        .cloned()
        .ok_or_else(|| invalid(format!("the criterion function returned no {name}")))
}

fn passed_on(output: &StepOutput, name: &str) -> Result<Column> {
    output
        .passed_on(name)
        .cloned()
        .ok_or_else(|| invalid(format!("the criterion function passed on no {name}")))
}

/// `np.array([0] * len(anno_df))`.
fn zeros(rows: usize) -> Column {
    Column::Int(vec![0; rows])
}

// ------------------------------------------------------------------------------ the step

/// The resolver handed to `annotate_exact_nonlof_variants` when no row has a symbol, where
/// Python builds no hub and the function returns before matching.
struct NoSymbols;

impl SymbolResolver for NoSymbols {
    fn resolve_symbol(&self, _gene_symbol: Cell<'_>) -> Result<String> {
        Err(invalid("no HGNC resolver: the table has no gene symbol"))
    }

    fn try_resolve_symbol(&self, gene_symbol: Cell<'_>) -> String {
        clean(Some(gene_symbol))
    }
}

/// `ACMG_criteria_assign(...)` with the parsed command line, in the script's order; the
/// final table, which is also written over the input.
pub(crate) fn run_step(
    args: &StepArguments,
    cutoffs: &Cutoffs,
    source: &pending::Source<'_>,
) -> Result<Table> {
    let read = read_table(Path::new(&args.anno_table))?;
    let kept = pandas_polars::drop_duplicates_first(&read, None)?;
    let mut table = if kept.len() == read.rows {
        read
    } else {
        read.take_rows(&kept)
    };
    eprintln!(
        "INFO: {} rows after drop_duplicates() of {}",
        table.rows, args.anno_table
    );
    vep_consq_interpret(&mut table, args.threads)?;

    let am_scores = read_table(Path::new(&args.am_score_table))?;
    let (ped, proband) = match args.ped_table.as_deref() {
        Some(path) if !path.is_empty() => {
            let ped = read_ped(path)?;
            let proband = proband(&ped, args.fam_name.as_deref())?;
            (Some(ped), Some(proband))
        }
        _ => (None, None),
    };
    let family = match (&ped, args.fam_name.as_deref()) {
        (Some(ped), Some(fam_name)) => Some((ped, fam_name)),
        _ => None,
    };

    let am_score_dict = am_score_dict(&am_scores)?;
    let transcript_to_gene = transcript_to_gene_map(&table)?;
    let gene_scores = gene_to_am_score_map(&transcript_to_gene, &am_score_dict);
    let clinvar = pending::load_clinvar_dicts(
        source,
        &args.clinvar_aa_dict_pkl,
        &args.clinvar_splice_dict_pkl,
    )?;
    let pathogenic_genes =
        pending::summarize_clinvar_gene_pathogenicity(source, &transcript_to_gene, &clinvar)?;
    let gene_avg = gene_avg_am_score(&table, &gene_scores)?;
    table.set_column("Gene_avg_AM_score", gene_avg);
    let variant_ids = variant_id(&table)?;
    table.set_column("variant_id", variant_ids);
    // `annotate_exact_nonlof_variants` builds `GeneMechanismHub(mechanism_json=...,
    // use_hgnc_package=False)` (which reads the HGNC table) only when a row has a symbol.
    let needs_hub = table.column_index("SYMBOL").is_some_and(|index| {
        (0..table.rows).any(|row| !clean(Some(table.cell(index, row))).is_empty())
    });
    let hub = if needs_hub {
        Some(GeneMechanismHub::new(HubPaths::default())?)
    } else {
        None
    };
    let resolver: &dyn SymbolResolver = match &hub {
        Some(hub) => hub,
        None => &NoSymbols,
    };
    annotate_exact_nonlof_variants(&mut table, None, &args.gene_nonlof_mechanism_json, resolver)?;
    drop(hub);

    // The script maps the file into memory and unpickles it; a `.gz` name maps the
    // compressed bytes, which fail to unpickle.
    let intolerant_domains = pickle_reader::load(&pickle_reader::read_file(
        Path::new(&args.intolerant_domains_pkl),
        false,
    )?)?;

    // `annotate_gene_mechanism_categories(anno_df, condition_cache=..., symbol_col="SYMBOL",
    // use_hgnc_package=False)`: the other hub files are the repository's.
    table = annotate_gene_mechanism_categories(
        &table,
        HubPaths {
            condition_cache: args.hpo_condition_mechanism_json.clone(),
            ..HubPaths::default()
        },
        "SYMBOL",
        "var_plausible_patho_mechs",
    )?;
    table.column("variant_effect")?;

    let pvs1_output = pending::pvs1_criteria(
        source,
        &table,
        &pending::Pvs1Inputs {
            clinvar_patho_exon_af_stat: &args.clinvar_patho_exon_af_stat,
            interpro_entry_map_pkl: &args.interpro_entry_map_pkl,
            clinvar_pathogenic_genes: &pathogenic_genes,
            gene_dosage_sensitivity: &args.gene_dosage_sensitivity,
            nonlof_mechanism_json: &args.gene_nonlof_mechanism_json,
            intolerant_domains: &intolerant_domains,
            tranx_exon_domain_map_pkl: &args.tranx_exon_domain_map_pkl,
            proband_gt_col: proband.as_deref(),
            cds_fasta_path: args.cds_fasta_path.as_deref(),
        },
    )?;
    apply_columns(&mut table, &pvs1_output)?;
    let pvs1 = criterion(&pvs1_output, "PVS1")?;
    let pvs1_values = int_values(&pvs1, "pvs1_criteria")?;
    let locate_intol_domains = passed_on(&pvs1_output, "locate_intol_domains")?;
    table.set_column("span_functional_domains", locate_intol_domains.clone());
    let locate_intol_domains = bool_values(&locate_intol_domains, "locate_intol_domains")?;

    let ps3_bs3 = ps3_bs3_criteria(&table, args.mavedb_metadata_tsv.as_deref())?;
    apply_columns(&mut table, &ps3_bs3)?;
    let clinvar_patho = bool_values(&passed_on(&ps3_bs3, "clinvar_patho")?, "clinvar_patho")?;
    let clinvar_benign = bool_values(&passed_on(&ps3_bs3, "clinvar_benign")?, "clinvar_benign")?;

    let ps1_pm5 = pending::ps1_pm5_criteria(
        source,
        &table,
        &clinvar,
        &clinvar_patho,
        &pvs1_values,
        args.threads,
        cutoffs.known_pathogenic_comparison,
    )?;
    apply_columns(&mut table, &ps1_pm5)?;

    let (ps2, pm6) = match family {
        Some((ped, fam_name)) => {
            let output = ps2_pm6_criteria(&table, ped, fam_name)?;
            (criterion(&output, "PS2")?, criterion(&output, "PM6")?)
        }
        None => {
            eprintln!("WARNING: No ped_table provided, skip the PS2 and PM6 criteria");
            (zeros(table.rows), zeros(table.rows))
        }
    };

    let pm1 = pm1_criteria(
        &table,
        &pvs1_values,
        &locate_intol_domains,
        &args.pm1_regions_pkl,
    )?;
    let loc_intol_domain = bool_values(&passed_on(&pm1, "loc_intol_domain")?, "loc_intol_domain")?;

    let pm4_bp3 = pm4_bp3_criteria(
        &table,
        &pvs1_values,
        &args.repeat_region_file,
        &loc_intol_domain,
    )?;
    apply_columns(&mut table, &pm4_bp3)?;

    let clinvar_gene_stat = load_pickle(&args.clinvar_gene_stat_pkl)?;
    let pp2_bp1 = pp2_bp1_criteria(
        &table,
        &clinvar_gene_stat,
        &args.am_intol_domains_tsv,
        &pending::pp2_pvalues(source)?,
    )?;
    drop(clinvar_gene_stat);

    let pp3_bp4 = pp3_bp4_criteria(&table, Some(&pvs1_values), cutoffs.no_splicing_impact)?;
    let pp5_bp6 = pp5_bp6_criteria(&table, &clinvar_patho, &clinvar_benign)?;

    let ba1 = pending::ba1_criteria(source, &table)?;
    apply_columns(&mut table, &ba1)?;
    let pm2_output = pending::pm2_criteria(
        source,
        &table,
        &args.clinvar_patho_af_stat,
        &args.clinvar_patho_exon_af_stat,
        args.gnomad_extreme_rare_threshold,
    )?;
    apply_columns(&mut table, &pm2_output)?;
    let pm2_before = int_values(&criterion(&pm2_output, "PM2")?, "pm2_criteria")?;
    let bs1_output = pending::bs1_criteria(source, &table, args.expected_incidence, &pm2_before)?;
    apply_columns(&mut table, &bs1_output)?;
    let bs1 = criterion(&bs1_output, "BS1")?;
    let component =
        |name: &str| -> Result<Vec<bool>> { bool_values(&passed_on(&bs1_output, name)?, name) };
    let (pm2_values, audit) = frequency_audit(
        &pm2_before,
        &bs1,
        &component("disease_incidence_bs1")?,
        &component("gene_pathogenic_af_bs1")?,
        &component("gene_pathogenic_af_blocked_by_pm2")?,
    )?;
    for (name, column) in &audit {
        table.set_column(name, column.clone());
    }
    let pm2 = Column::Int(pm2_values.clone());

    let bs2 = pending::bs2_criteria(source, &table, &pm2_values)?;
    apply_columns(&mut table, &bs2)?;
    let pp1 = pp1_criteria(&table, &args.pp1_vcf, &args.pp1_ped)?;
    let bs4 = match family {
        Some((ped, fam_name)) => criterion(&bs4_criteria(&table, ped, fam_name)?, "BS4")?,
        None => {
            eprintln!("WARNING: No ped_table provided, skip the BS4 criteria");
            zeros(table.rows)
        }
    };
    let bp5 = match args.alt_disease_vcf.as_deref() {
        Some(vcf) if !vcf.is_empty() => {
            let output = pending::bp5_criteria(source, &table, vcf, args.threads)?;
            apply_columns(&mut table, &output)?;
            criterion(&output, "BP5")?
        }
        _ => zeros(table.rows),
    };
    let bp7 = bp7_criteria(&table)?;
    let ps1 = criterion(&ps1_pm5, "PS1")?;
    let ps3 = criterion(&ps3_bs3, "PS3")?;
    let (bp2, pm3) = match family {
        Some((ped, _)) => {
            let output = bp2_pm3_criteria(&table, ped, &pm2, &pvs1, &ps1, &ps3)?;
            apply_columns(&mut table, &output)?;
            (criterion(&output, "BP2")?, criterion(&output, "PM3")?)
        }
        None => {
            eprintln!("WARNING: No ped_table provided, skip the BP2 and PM3 criteria");
            (zeros(table.rows), zeros(table.rows))
        }
    };

    let criteria: Vec<(&str, Column)> = vec![
        ("PVS1", pvs1),
        ("PS1", ps1),
        ("PS2", ps2),
        ("PS3", ps3),
        ("PM1", criterion(&pm1, "PM1")?),
        ("PM2", pm2),
        ("PM3", pm3),
        ("PM4", criterion(&pm4_bp3, "PM4")?),
        ("PM5", criterion(&ps1_pm5, "PM5")?),
        ("PM6", pm6),
        ("PP1", criterion(&pp1, "PP1")?),
        ("PP2", criterion(&pp2_bp1, "PP2")?),
        ("PP3", criterion(&pp3_bp4, "PP3")?),
        ("PP5", criterion(&pp5_bp6, "PP5")?),
        ("BA1", criterion(&ba1, "BA1")?),
        ("BS1", bs1),
        ("BS2", criterion(&bs2, "BS2")?),
        ("BS3", criterion(&ps3_bs3, "BS3")?),
        ("BS4", bs4),
        ("BP1", criterion(&pp2_bp1, "BP1")?),
        ("BP2", bp2),
        ("BP3", criterion(&pm4_bp3, "BP3")?),
        ("BP4", criterion(&pp3_bp4, "BP4")?),
        ("BP5", bp5),
        ("BP6", criterion(&pp5_bp6, "BP6")?),
        ("BP7", criterion(&bp7, "BP7")?),
    ];
    let mut matrix = summarize_acmg_criteria(
        &mut table,
        &criteria,
        args.clingen_map_pkl.as_deref(),
        !args.ignore_clingen_override,
    )?;
    drop(criteria);
    for name in MATRIX_AUDIT_COLUMNS {
        let column = table.column(name)?.clone();
        matrix.set_column(name, column);
    }
    let output_matrix = matrix_path(&args.anno_table);
    write_table(Path::new(&output_matrix), &matrix)?;
    eprintln!("INFO: criteria matrix written to {output_matrix}");

    let (scores, classes) = posterior_columns(&matrix)?;
    drop(matrix);
    table.set_column("ACMG_quant_score", scores);
    table.set_column("ACMG_class", classes);

    let options = RankingOptions {
        pext_tissues: &args.pext_tissues,
        pext_low_expression_cutoff: args.pext_low_expression_cutoff,
        pext_penalty_floor: args.pext_penalty_floor,
        pext_penalty_shape: args.pext_penalty_shape,
        relevant_gene_list: args.relevant_gene_list.as_deref(),
        dispensable_gene_list: args.dispensable_gene_list.as_deref(),
    };
    let control_common = |rows: &mut Table, proband: Option<&str>| {
        pending::compute_control_common(source, rows, proband)
    };
    let ranked = sort_and_rank_variants_with(
        &table,
        ped.as_ref(),
        args.fam_name.as_deref(),
        &options,
        &control_common,
    )?;
    drop(table);
    write_table(Path::new(&args.anno_table), &ranked)?;
    eprintln!("INFO: annotated table written to {}", args.anno_table);
    Ok(ranked)
}

/// The step's exit for an error: Python's uncaught exception ends with status 1.
fn fail(error: impl std::fmt::Display) -> ! {
    eprintln!("priva acmg-assign: {error}");
    std::process::exit(1)
}

/// `priva acmg-assign`: the converted `acmg_criteria_assign.py` command line.
pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let parsed = match parse_arguments(args) {
        Ok(parsed) => parsed,
        Err(ArgumentError::Help) => {
            println!("{}", usage());
            return Ok(());
        }
        Err(ArgumentError::Usage(message)) => {
            eprintln!("{}", usage());
            eprintln!("priva acmg-assign: error: {message}");
            std::process::exit(2);
        }
        Err(ArgumentError::Refused(message)) => fail(message),
    };
    // `cutoffs_for_config(args.config)` is an argument of the call, so it runs first.
    let cutoffs = match cutoffs_for_config(parsed.config.as_deref().map(Path::new)) {
        Ok(cutoffs) => cutoffs,
        Err(failure) => fail(format!("SpliceAI cutoffs: {failure:?}")),
    };
    match run_step(&parsed, &cutoffs, &pending::Source::Unconverted) {
        Ok(_) => Ok(()),
        Err(error) => fail(error),
    }
}

// ------------------------------------------------------------------------------ pending

/// The calls into threads that are still converting their scripts (thread T4: PVS1, PS1/PM5
/// and the ClinVar dictionaries; T5: BA1, PM2, BS1, BS2, BP5 and `stats`, which
/// T3's `compute_control_common` in the ranking calls). Each function takes what the Python call
/// takes. Until the thread lands, `Source::Unconverted` stops the step with "not converted";
/// the merge replaces each body by the real call and removes the function's recorded branch.
///
/// `Source::Recorded(folder)` is for the comparison harness only (`priva acmg-step
/// acmg_criteria_assign ACMG_criteria_assign`): a Python run of the step with these functions
/// wrapped wrote, per function, its returned arrays (`<function>.returns.frame`) and the table
/// columns it added or changed (`<function>.columns.frame`) and the SciPy p-values of PP2
/// (`pvalues.json`).
pub(crate) mod pending {
    use super::super::acmg_pm1_pp2_bp1_missense::PValues;
    use super::super::frame_file::read_frame;
    use super::super::types::{StepOutput, not_converted};
    use super::Column;
    use crate::cli::{Result, invalid};
    use crate::json_text::{self, Json};
    use crate::pandas_table::Table;
    use crate::pickle_reader::Value;
    use std::collections::HashMap;
    use std::path::Path;

    pub(crate) enum Source<'a> {
        Unconverted,
        Recorded(&'a Path),
    }

    /// The two ClinVar dictionaries (`aa_change.stats.pkl`, `splice_change.stats.pkl`); T4
    /// decides their Rust shape (GROUP3_BASE.md, section 8.2).
    pub(crate) struct ClinvarDicts;

    /// `summarize_clinvar_gene_pathogenicity`'s set of gene IDs (T4).
    pub(crate) struct PathogenicGenes;

    /// What `PVS1_criteria` takes besides the table.
    pub(crate) struct Pvs1Inputs<'a> {
        pub(crate) clinvar_patho_exon_af_stat: &'a str,
        pub(crate) interpro_entry_map_pkl: &'a str,
        pub(crate) clinvar_pathogenic_genes: &'a PathogenicGenes,
        pub(crate) gene_dosage_sensitivity: &'a str,
        pub(crate) nonlof_mechanism_json: &'a str,
        pub(crate) intolerant_domains: &'a Value,
        pub(crate) tranx_exon_domain_map_pkl: &'a str,
        pub(crate) proband_gt_col: Option<&'a str>,
        pub(crate) cds_fasta_path: Option<&'a str>,
    }

    /// The returned arrays and changed columns a recording Python run wrote for `function`.
    fn recorded(folder: &Path, function: &str) -> Result<(Table, Table)> {
        let returns = read_frame(&folder.join(format!("{function}.returns.frame")))?;
        let columns = read_frame(&folder.join(format!("{function}.columns.frame")))?;
        Ok((returns, columns))
    }

    /// A `StepOutput` from a recording: `criteria` and `passed` name return frame columns.
    fn recorded_output(
        folder: &Path,
        function: &str,
        criteria: &[&'static str],
        passed: &[&str],
    ) -> Result<StepOutput> {
        let (returns, columns) = recorded(folder, function)?;
        let take = |name: &str| -> Result<Column> {
            returns.column(name).cloned().map_err(|_| {
                invalid(format!(
                    "the recording of {function} has no returned {name}"
                ))
            })
        };
        Ok(StepOutput {
            criteria: criteria
                .iter()
                .map(|name| Ok((*name, take(name)?)))
                .collect::<Result<_>>()?,
            columns: columns.names.into_iter().zip(columns.columns).collect(),
            passed_on: passed
                .iter()
                .map(|name| Ok(((*name).to_owned(), take(name)?)))
                .collect::<Result<_>>()?,
        })
    }

    pub(crate) fn load_clinvar_dicts(
        source: &Source<'_>,
        aa: &str,
        splice: &str,
    ) -> Result<ClinvarDicts> {
        match source {
            Source::Recorded(_) => Ok(ClinvarDicts),
            Source::Unconverted => {
                let _ = (aa, splice);
                Err(not_converted(
                    "acmg_criteria_assign",
                    "load of the ClinVar dictionaries (thread T4)",
                ))
            }
        }
    }

    pub(crate) fn summarize_clinvar_gene_pathogenicity(
        source: &Source<'_>,
        transcript_to_gene: &[(Option<String>, Option<String>)],
        clinvar: &ClinvarDicts,
    ) -> Result<PathogenicGenes> {
        match source {
            Source::Recorded(_) => Ok(PathogenicGenes),
            Source::Unconverted => {
                let _ = (transcript_to_gene, clinvar);
                Err(not_converted(
                    "acmg_pvs1_null_variant",
                    "summarize_clinvar_gene_pathogenicity",
                ))
            }
        }
    }

    /// `PVS1_criteria(df, ...)` (T4): `PVS1` and the passed-on `locate_intol_domains`, and the
    /// `alt_start_*` columns it may add.
    pub(crate) fn pvs1_criteria(
        source: &Source<'_>,
        table: &Table,
        inputs: &Pvs1Inputs<'_>,
    ) -> Result<StepOutput> {
        match source {
            Source::Recorded(folder) => recorded_output(
                folder,
                "PVS1_criteria",
                &["PVS1"],
                &["locate_intol_domains"],
            ),
            Source::Unconverted => {
                let _ = (table, inputs);
                Err(not_converted("acmg_pvs1_null_variant", "PVS1_criteria"))
            }
        }
    }

    /// `PS1_PM5_criteria(df, clinvar_aa_dict, clinvar_splice_dict, clinvar_patho, pvs1_criteria,
    /// threads, spliceai_cutoff=...)` (T4).
    pub(crate) fn ps1_pm5_criteria(
        source: &Source<'_>,
        table: &Table,
        clinvar: &ClinvarDicts,
        clinvar_patho: &[bool],
        pvs1_criteria: &[i64],
        threads: i64,
        spliceai_cutoff: f64,
    ) -> Result<StepOutput> {
        match source {
            Source::Recorded(folder) => {
                recorded_output(folder, "PS1_PM5_criteria", &["PS1", "PM5"], &[])
            }
            Source::Unconverted => {
                let _ = (
                    table,
                    clinvar,
                    clinvar_patho,
                    pvs1_criteria,
                    threads,
                    spliceai_cutoff,
                );
                Err(not_converted(
                    "acmg_ps1_pm5_known_variant",
                    "PS1_PM5_criteria",
                ))
            }
        }
    }

    /// `BA1_criteria(anno_df)` (T5).
    pub(crate) fn ba1_criteria(source: &Source<'_>, table: &Table) -> Result<StepOutput> {
        match source {
            Source::Recorded(folder) => recorded_output(folder, "BA1_criteria", &["BA1"], &[]),
            Source::Unconverted => {
                let _ = table;
                Err(not_converted("acmg_pm2_bs1_ba1_frequency", "BA1_criteria"))
            }
        }
    }

    /// `PM2_criteria(anno_df, clinvar_patho_af_stat, clinvar_patho_exon_af_stat,
    /// gnomAD_extreme_rare_threshold)` (T5), with the four columns it writes.
    pub(crate) fn pm2_criteria(
        source: &Source<'_>,
        table: &Table,
        clinvar_patho_af_stat: &str,
        clinvar_patho_exon_af_stat: &str,
        gnomad_extreme_rare_threshold: f64,
    ) -> Result<StepOutput> {
        match source {
            Source::Recorded(folder) => recorded_output(folder, "PM2_criteria", &["PM2"], &[]),
            Source::Unconverted => {
                let _ = (
                    table,
                    clinvar_patho_af_stat,
                    clinvar_patho_exon_af_stat,
                    gnomad_extreme_rare_threshold,
                );
                Err(not_converted("acmg_pm2_bs1_ba1_frequency", "PM2_criteria"))
            }
        }
    }

    /// `BS1_criteria(anno_df, expected_incidence=..., pm2_criteria=...,
    /// return_frequency_components=True)` (T5): `BS1` and the three components passed on.
    pub(crate) fn bs1_criteria(
        source: &Source<'_>,
        table: &Table,
        expected_incidence: f64,
        pm2_criteria: &[i64],
    ) -> Result<StepOutput> {
        match source {
            Source::Recorded(folder) => recorded_output(
                folder,
                "BS1_criteria",
                &["BS1"],
                &[
                    "disease_incidence_bs1",
                    "gene_pathogenic_af_bs1",
                    "gene_pathogenic_af_blocked_by_pm2",
                ],
            ),
            Source::Unconverted => {
                let _ = (table, expected_incidence, pm2_criteria);
                Err(not_converted("acmg_pm2_bs1_ba1_frequency", "BS1_criteria"))
            }
        }
    }

    /// `BS2_criteria(anno_df, pm2_criteria)` (T5), with PM2 after the BS1 back-down.
    pub(crate) fn bs2_criteria(
        source: &Source<'_>,
        table: &Table,
        pm2_criteria: &[i64],
    ) -> Result<StepOutput> {
        match source {
            Source::Recorded(folder) => recorded_output(folder, "BS2_criteria", &["BS2"], &[]),
            Source::Unconverted => {
                let _ = (table, pm2_criteria);
                Err(not_converted("acmg_bs2_bp5_observation", "BS2_criteria"))
            }
        }
    }

    /// `BP5_criteria(anno_df, alt_disease_vcf, n_processes=threads)` (T5), with the
    /// `alt_disease_*` columns and the rebuilt `variant_id`.
    pub(crate) fn bp5_criteria(
        source: &Source<'_>,
        table: &Table,
        alt_disease_vcf: &str,
        threads: i64,
    ) -> Result<StepOutput> {
        match source {
            Source::Recorded(folder) => recorded_output(folder, "BP5_criteria", &["BP5"], &[]),
            Source::Unconverted => {
                let _ = (table, alt_disease_vcf, threads);
                Err(not_converted("acmg_bs2_bp5_observation", "BP5_criteria"))
            }
        }
    }

    /// Where PP2's chi-square and Fisher p-values come from: SciPy's functions in `stats`
    /// (T5), or in a recorded comparison the p-values SciPy gave in the Python run.
    pub(crate) fn pp2_pvalues(source: &Source<'_>) -> Result<PValues> {
        match source {
            Source::Unconverted => Ok(PValues::Scipy),
            Source::Recorded(folder) => {
                let path = folder.join("pvalues.json");
                if !path.exists() {
                    return Ok(PValues::Supplied(HashMap::new()));
                }
                let text = std::fs::read_to_string(path)?;
                let Json::Array(entries) = json_text::parse(&text, true).map_err(invalid)? else {
                    return Err(invalid("pvalues.json is a list"));
                };
                let mut values = HashMap::new();
                for entry in entries {
                    let Json::Array(items) = entry else {
                        return Err(invalid("a p-value entry is a list"));
                    };
                    let (Some(Json::Str(test)), Some(Json::Str(bits))) =
                        (items.first(), items.get(5))
                    else {
                        return Err(invalid("a p-value entry is [test, a, b, c, d, p bits]"));
                    };
                    let mut cells = [0i64; 4];
                    for (position, cell) in cells.iter_mut().enumerate() {
                        let Some(Json::Int(number)) = items.get(position + 1) else {
                            return Err(invalid("table cells are integers"));
                        };
                        *cell = i64::try_from(*number)
                            .map_err(|_| invalid("a table cell beyond 64 bits"))?;
                    }
                    let pvalue = f64::from_bits(
                        u64::from_str_radix(bits, 16)
                            .map_err(|_| invalid("p-value bits are hex"))?,
                    );
                    values.insert((test == "chi_square", cells), pvalue);
                }
                Ok(PValues::Supplied(values))
            }
        }
    }

    /// `compute_control_common(df, proband)` inside `sort_and_rank_variants` (T3's, calling
    /// T5's `stats::binom_cdf` and `stats::beta_ppf`, still placeholders): in a recorded
    /// comparison the `control_common` and `control_common_index` columns the Python run
    /// wrote for the ranked rows. Once `stats` lands, the recorded branch goes.
    pub(crate) fn compute_control_common(
        source: &Source<'_>,
        rows: &mut Table,
        proband: Option<&str>,
    ) -> Result<()> {
        match source {
            Source::Unconverted => {
                super::super::acmg_scoring::compute_control_common(rows, proband)
            }
            Source::Recorded(folder) => {
                let (_, columns) = recorded(folder, "compute_control_common")?;
                if columns.rows != rows.rows {
                    return Err(invalid(format!(
                        "the recorded control-common columns have {} rows for {} ranked rows",
                        columns.rows, rows.rows
                    )));
                }
                for (name, column) in columns.names.iter().zip(columns.columns) {
                    rows.set_column(name, column);
                }
                Ok(())
            }
        }
    }
}

// ------------------------------------------------------------------------------ harness

fn text_option<'a>(options: &'a Json, name: &str) -> Result<&'a str> {
    options
        .get(name)
        .and_then(Json::as_str)
        .ok_or_else(|| invalid(format!("options need the text {name:?}")))
}

fn argv_option(options: &Json) -> Result<Vec<String>> {
    match options.get("argv") {
        Some(Json::Array(items)) => items
            .iter()
            .map(|item| {
                item.as_str()
                    .map(str::to_owned)
                    .ok_or_else(|| invalid("argv holds texts"))
            })
            .collect(),
        _ => Err(invalid("options need \"argv\": [...]")),
    }
}

fn one_row(names: Vec<String>, values: Vec<ObjectValue>) -> Table {
    Table {
        names,
        columns: values
            .into_iter()
            .map(|value| Column::Object(vec![value]))
            .collect(),
        rows: 1,
    }
}

/// `vars(parser.parse_args(argv))` as a one-row frame (the option names, in the parser's
/// order), or `__error` = `usage` / `help` where argparse exits.
fn parsed_frame(argv: &[String]) -> Result<Table> {
    let parsed = match parse_arguments(argv) {
        Ok(parsed) => parsed,
        Err(ArgumentError::Help) => {
            return Ok(one_row(
                vec!["__error".into()],
                vec![ObjectValue::Text("help".into())],
            ));
        }
        Err(ArgumentError::Usage(_)) => {
            return Ok(one_row(
                vec!["__error".into()],
                vec![ObjectValue::Text("usage".into())],
            ));
        }
        Err(ArgumentError::Refused(message)) => return Err(invalid(message)),
    };
    let text = |value: &Option<String>| match value {
        Some(text) => ObjectValue::Text(text.clone()),
        None => ObjectValue::None,
    };
    let given = |value: &str| ObjectValue::Text(value.to_owned());
    let values: Vec<(&str, ObjectValue)> = vec![
        ("anno_table", given(&parsed.anno_table)),
        ("am_score_table", given(&parsed.am_score_table)),
        ("ped_table", text(&parsed.ped_table)),
        ("fam_name", text(&parsed.fam_name)),
        (
            "clinvar_patho_af_stat",
            given(&parsed.clinvar_patho_af_stat),
        ),
        (
            "clinvar_patho_exon_af_stat",
            given(&parsed.clinvar_patho_exon_af_stat),
        ),
        ("clinvar_aa_dict_pkl", given(&parsed.clinvar_aa_dict_pkl)),
        (
            "clinvar_splice_dict_pkl",
            given(&parsed.clinvar_splice_dict_pkl),
        ),
        (
            "interpro_entry_map_pkl",
            given(&parsed.interpro_entry_map_pkl),
        ),
        (
            "intolerant_domains_pkl",
            given(&parsed.intolerant_domains_pkl),
        ),
        (
            "gene_dosage_sensitivity",
            given(&parsed.gene_dosage_sensitivity),
        ),
        ("am_intol_domains_tsv", given(&parsed.am_intol_domains_tsv)),
        ("pm1_regions_pkl", given(&parsed.pm1_regions_pkl)),
        (
            "clinvar_gene_stat_pkl",
            given(&parsed.clinvar_gene_stat_pkl),
        ),
        (
            "tranx_exon_domain_map_pkl",
            given(&parsed.tranx_exon_domain_map_pkl),
        ),
        ("am_score_vcf", text(&parsed.am_score_vcf)),
        ("alt_disease_vcf", text(&parsed.alt_disease_vcf)),
        ("clingen_map_pkl", text(&parsed.clingen_map_pkl)),
        (
            "ignore_clingen_override",
            ObjectValue::Bool(parsed.ignore_clingen_override),
        ),
        ("mavedb_metadata_tsv", text(&parsed.mavedb_metadata_tsv)),
        ("pp1_vcf", given(&parsed.pp1_vcf)),
        ("pp1_ped", given(&parsed.pp1_ped)),
        ("cds_fasta_path", text(&parsed.cds_fasta_path)),
        ("pext_tissues", given(&parsed.pext_tissues)),
        (
            "pext_low_expression_cutoff",
            ObjectValue::Float(parsed.pext_low_expression_cutoff),
        ),
        (
            "pext_penalty_floor",
            ObjectValue::Float(parsed.pext_penalty_floor),
        ),
        (
            "pext_penalty_shape",
            ObjectValue::Float(parsed.pext_penalty_shape),
        ),
        ("repeat_region_file", given(&parsed.repeat_region_file)),
        (
            "gnomAD_extreme_rare_threshold",
            ObjectValue::Float(parsed.gnomad_extreme_rare_threshold),
        ),
        (
            "expected_incidence",
            ObjectValue::Float(parsed.expected_incidence),
        ),
        ("relevant_gene_list", text(&parsed.relevant_gene_list)),
        ("dispensable_gene_list", text(&parsed.dispensable_gene_list)),
        (
            "hpo_condition_mechanism_json",
            given(&parsed.hpo_condition_mechanism_json),
        ),
        (
            "gene_nonlof_mechanism_json",
            given(&parsed.gene_nonlof_mechanism_json),
        ),
        ("threads", ObjectValue::Int(parsed.threads)),
        ("config", text(&parsed.config)),
        ("ranking_binary", text(&parsed.ranking_binary)),
    ];
    let (names, values) = values
        .into_iter()
        .map(|(name, value)| (name.to_owned(), value))
        .unzip();
    Ok(one_row(names, values))
}

/// `priva acmg-step acmg_criteria_assign <function>`:
///
/// - `parse_args`, option `argv`: `vars(args)` as one row (defaults that point into `data/`
///   are the data folder's), or `__error`;
/// - `read_ped`, option `ped_table`: the pedigree frame; `proband`, options `ped_table` and
///   `fam_name` (text or null): a frame of `proband`;
/// - `gene_avg_am_score`, option `am_score_table`, on the input table: a frame of
///   `Gene_avg_AM_score`;
/// - `variant_id`, on the input table: a frame of `variant_id`;
/// - `frequency_audit`, on a frame of `PM2`, `BS1`, `disease_incidence_bs1`,
///   `gene_pathogenic_af_bs1`, `gene_pathogenic_af_blocked_by_pm2`: a frame of `PM2` after the
///   back-down and the five audit columns;
/// - `ACMG_criteria_assign`, options `argv` and `substitutes` (the folder of a recording
///   Python run): the step with the pending functions taken from the recording; the result is
///   the final table (also written over `--anno_table`, with the matrix beside it).
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let options = &call.options;
    let single = |name: &str, column: Column| Table {
        names: vec![name.to_owned()],
        rows: column.len(),
        columns: vec![column],
    };
    match function {
        "parse_args" => parsed_frame(&argv_option(options)?),
        "read_ped" => read_ped(text_option(options, "ped_table")?),
        "proband" => {
            let ped = read_ped(text_option(options, "ped_table")?)?;
            let fam_name = options.get("fam_name").and_then(Json::as_str);
            Ok(single(
                "proband",
                Column::Text(vec![Some(proband(&ped, fam_name)?)]),
            ))
        }
        "gene_avg_am_score" => {
            let am_scores = read_table(Path::new(text_option(options, "am_score_table")?))?;
            let dict = am_score_dict(&am_scores)?;
            let pairs = transcript_to_gene_map(&call.table)?;
            let map = gene_to_am_score_map(&pairs, &dict);
            Ok(single(
                "Gene_avg_AM_score",
                gene_avg_am_score(&call.table, &map)?,
            ))
        }
        "variant_id" => Ok(single("variant_id", variant_id(&call.table)?)),
        "frequency_audit" => {
            let table = &call.table;
            let flags = |name: &str| bool_values(table.column(name)?, name);
            let (pm2, audit) = frequency_audit(
                &int_values(table.column("PM2")?, "PM2")?,
                table.column("BS1")?,
                &flags("disease_incidence_bs1")?,
                &flags("gene_pathogenic_af_bs1")?,
                &flags("gene_pathogenic_af_blocked_by_pm2")?,
            )?;
            let mut out = single("PM2", Column::Int(pm2));
            for (name, column) in audit {
                out.set_column(name, column);
            }
            Ok(out)
        }
        "ACMG_criteria_assign" => {
            let argv = argv_option(options)?;
            let parsed = match parse_arguments(&argv) {
                Ok(parsed) => parsed,
                Err(error) => return Err(invalid(format!("the command line fails: {error:?}"))),
            };
            let folder = PathBuf::from(text_option(options, "substitutes")?);
            let cutoffs = cutoffs_for_config(parsed.config.as_deref().map(Path::new))
                .map_err(|failure| invalid(format!("SpliceAI cutoffs: {failure:?}")))?;
            run_step(&parsed, &cutoffs, &pending::Source::Recorded(&folder))
        }
        _ => Err(not_converted("acmg_criteria_assign", function)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(items: &[&str]) -> Vec<String> {
        items.iter().map(|item| (*item).to_owned()).collect()
    }

    fn required_args() -> Vec<String> {
        let mut args = Vec::new();
        for (name, _, required) in OPTIONS {
            if required {
                args.push(format!("--{name}"));
                args.push(format!("{name}.value"));
            }
        }
        args
    }

    #[test]
    fn the_command_line_follows_argparse() {
        let parsed = parse_arguments(&required_args()).unwrap();
        assert_eq!(parsed.anno_table, "anno_table.value");
        assert_eq!(parsed.threads, 10);
        assert_eq!(parsed.gnomad_extreme_rare_threshold, 0.0001);
        assert!(!parsed.ignore_clingen_override);
        assert_eq!(parsed.ped_table, None);
        assert_eq!(parsed.pp1_vcf, "");
        let mut args = required_args();
        args.extend(strings(&[
            "--ignore",
            "--threads",
            " +4 ",
            "--expected_incidence=1e-5",
            "--pp1_vcf",
            "fam.vcf.gz",
            "--pp1_ped",
            "fam.ped",
            "--ranking_binary",
            "priva",
        ]));
        let parsed = parse_arguments(&args).unwrap();
        assert!(parsed.ignore_clingen_override);
        assert_eq!(parsed.threads, 4);
        assert_eq!(parsed.expected_incidence, 1e-5);
        assert_eq!(parsed.pp1_vcf, "fam.vcf.gz");
        assert_eq!(parsed.ranking_binary.as_deref(), Some("priva"));
        // Every occurrence is converted, also one a later occurrence replaces.
        let mut args = required_args();
        args.extend(strings(&["--threads", "x", "--threads", "4"]));
        assert!(matches!(
            parse_arguments(&args),
            Err(ArgumentError::Usage(_))
        ));
        let mut args = required_args();
        args.extend(strings(&[
            "--gnomAD_extreme_rare_threshold",
            "--threads",
            "4",
        ]));
        assert!(matches!(
            parse_arguments(&args),
            Err(ArgumentError::Usage(_))
        ));
        let mut args = required_args();
        args.extend(strings(&["--threads", "４"]));
        assert!(matches!(
            parse_arguments(&args),
            Err(ArgumentError::Refused(_))
        ));
        assert_eq!(
            parse_arguments(&strings(&["--anno_table", "a"])),
            Err(ArgumentError::Usage(
                "the following arguments are required: --am_score_table, \
                 --clinvar_patho_af_stat, --clinvar_patho_exon_af_stat, --clinvar_aa_dict_pkl, \
                 --clinvar_splice_dict_pkl, --interpro_entry_map_pkl, --intolerant_domains_pkl, \
                 --gene_dosage_sensitivity, --am_intol_domains_tsv, --pm1_regions_pkl, \
                 --clinvar_gene_stat_pkl, --tranx_exon_domain_map_pkl, --repeat_region_file"
                    .to_owned()
            ))
        );
        assert_eq!(parse_arguments(&strings(&["-h"])), Err(ArgumentError::Help));
    }

    #[test]
    fn the_matrix_path_splits_the_whole_path_text() {
        assert_eq!(matrix_path("/a/b/fam.tsv"), "/a/b/fam.acmg.tsv");
        assert_eq!(matrix_path("/a.b/fam"), "/a.acmg.tsv");
        assert_eq!(matrix_path("fam"), ".acmg.tsv");
        assert_eq!(matrix_path("fam.x.tsv"), "fam.x.acmg.tsv");
    }

    #[test]
    fn the_frequency_audit_backs_pm2_down() {
        let (pm2, audit) = frequency_audit(
            &[2, 1, 0, 2],
            &Column::Int(vec![1, 0, 1, 1]),
            &[true, false, true, false],
            &[false, true, false, false],
            &[false, false, true, false],
        )
        .unwrap();
        assert_eq!(pm2, vec![0, 1, 0, 2]);
        assert_eq!(audit[3].1, Column::Int(vec![1, 0, 0, 0]));
        assert_eq!(audit[4].1, Column::Int(vec![0, 0, 0, 1]));
    }

    #[test]
    fn the_am_scores_map_through_transcripts_to_genes() {
        let text = |values: &[Option<&str>]| {
            Column::Text(
                values
                    .iter()
                    .map(|value| value.map(str::to_owned))
                    .collect(),
            )
        };
        let table = Table {
            names: vec!["Feature".into(), "Gene".into()],
            columns: vec![
                text(&[Some("T1"), Some("T2"), Some("T1"), None, Some("T3")]),
                text(&[Some("G1"), Some("G2"), Some("G9"), Some("G3"), None]),
            ],
            rows: 5,
        };
        let pairs = transcript_to_gene_map(&table).unwrap();
        assert_eq!(pairs[0], (Some("T1".into()), Some("G9".into())));
        assert_eq!(pairs.len(), 4);
        let am = HashMap::from([("T1".to_owned(), 0.5), ("T3".to_owned(), 0.25)]);
        let map = gene_to_am_score_map(&pairs, &am);
        let Column::Float(values) = gene_avg_am_score(&table, &map).unwrap() else {
            panic!("a float column");
        };
        // G1 lost T1 to G9; the NaN gene finds T3's score through the NaN key.
        assert!(values[0].is_nan());
        assert!(values[1].is_nan());
        assert_eq!(values[2], 0.5);
        assert!(values[3].is_nan());
        assert_eq!(values[4], 0.25);
    }
}
