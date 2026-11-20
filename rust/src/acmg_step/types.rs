//! What the group 3 modules share. These types and signatures are fixed for
//! the seven threads: a thread that needs a change agrees it with the threads
//! that use the item, and the change lands as its own commit.
//!
//! The table is `pandas_table::Table`: each column keeps its pandas type, and
//! `pandas_polars` moves typed columns to polars where a polars call gives the
//! pandas result (`rust/GROUP3_BASE.md` lists, call by call, which pandas
//! calls go to polars and which stay own code).
use crate::cli::Result;
use crate::json_text::Json;
use crate::pandas_table::{Cell, Column, Table};
use std::fmt;

/// The 28 criteria in the column order of the criteria matrix
/// (`summarize_acmg_criteria`, `acmg_scoring.py`). PS4 and PP4 are always 0.
pub(crate) const CRITERIA_ORDER: [&str; 28] = [
    "PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4", "PM5", "PM6", "PP1", "PP2",
    "PP3", "PP4", "PP5", "BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP5",
    "BP6", "BP7",
];

/// The five audit columns copied into the criteria matrix after the criteria.
pub(crate) const MATRIX_AUDIT_COLUMNS: [&str; 5] = [
    "audit_BS1_disease_incidence_gate",
    "audit_BS1_gene_pathogenic_AF_gate",
    "audit_BS1_gene_pathogenic_AF_blocked_by_PM2",
    "audit_PM2_backdown_for_BS1_disease_incidence",
    "frequency_conflict_PM2_BS1",
];

/// What a criterion function hands back to the orchestrator (`acmg_criteria_assign`).
///
/// A criterion's values are a `Column` because the Python arrays differ in
/// type: most are int64 strengths 0 to 5 (`Column::Int`), some are numpy
/// booleans (`Column::Bool`), and the criteria matrix prints them as they are
/// (`2` or `True`). The orchestrator writes `columns` into the table in the
/// order given, as the Python function's `df[...] = ...` lines do.
#[derive(Clone, Debug, Default)]
pub(crate) struct StepOutput {
    /// `(criterion, values)` in the order the Python function returns them.
    pub(crate) criteria: Vec<(&'static str, Column)>,
    /// Columns the Python function adds to the table or overwrites, in its order.
    pub(crate) columns: Vec<(String, Column)>,
    /// Values handed to later criteria but not written into the table
    /// (for example `clinvar_patho` from PS3_BS3_criteria, `locate_intol_domains` from PVS1).
    pub(crate) passed_on: Vec<(String, Column)>,
}

impl StepOutput {
    pub(crate) fn criterion(&self, name: &str) -> Option<&Column> {
        self.criteria
            .iter()
            .find(|(criterion, _)| *criterion == name)
            .map(|(_, values)| values)
    }

    pub(crate) fn passed_on(&self, name: &str) -> Option<&Column> {
        self.passed_on
            .iter()
            .find(|(key, _)| key == name)
            .map(|(_, values)| values)
    }
}

/// Named boolean masks, as the Python helpers return them in a dict
/// (`_variant_condition_masks`, `_variant_mechanism_masks`).
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct NamedMasks {
    pub(crate) masks: Vec<(String, Vec<bool>)>,
}

impl NamedMasks {
    /// `masks[name]`; a missing name is the KeyError Python would raise.
    pub(crate) fn get(&self, name: &str) -> Result<&[bool]> {
        self.masks
            .iter()
            .find(|(key, _)| key == name)
            .map(|(_, mask)| mask.as_slice())
            .ok_or_else(|| AcmgError::new("KeyError", format!("'{name}'")).into())
    }
}

/// A Python exception the step would raise: its type name and message. The
/// step stops with exit status 1 and prints both, as Python's traceback ends.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct AcmgError {
    pub(crate) python_type: &'static str,
    pub(crate) message: String,
}

impl AcmgError {
    pub(crate) fn new(python_type: &'static str, message: impl Into<String>) -> AcmgError {
        AcmgError {
            python_type,
            message: message.into(),
        }
    }
}

impl fmt::Display for AcmgError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}: {}", self.python_type, self.message)
    }
}

impl std::error::Error for AcmgError {}

/// The error of a module or function that has not been converted yet.
pub(crate) fn not_converted(module: &str, function: &str) -> Box<dyn std::error::Error> {
    Box::new(AcmgError::new(
        "NotImplementedError",
        format!("{module}.{function} is not converted to Rust yet"),
    ))
}

/// The command line of `acmg_criteria_assign.py` after argparse, with the
/// Python option names and defaults. Paths stay text because the Python code
/// tests them as text (`path.endswith(".gz")`). The defaults that point into
/// the repository's `data/` folder (`--dispensable_gene_list`,
/// `--hpo_condition_mechanism_json`, `--gene_nonlof_mechanism_json`, and the
/// hub's HGNC, ClinGen hg19 and LOEUF files) are resolved from
/// `protein_domain_mapping::data_directory()`: the environment variable
/// `PRIVA_DATA_DIR`, else the `data/` folder of the source tree.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct StepArguments {
    pub(crate) anno_table: String,
    pub(crate) am_score_table: String,
    pub(crate) ped_table: Option<String>,
    pub(crate) fam_name: Option<String>,
    pub(crate) clinvar_patho_af_stat: String,
    pub(crate) clinvar_patho_exon_af_stat: String,
    pub(crate) clinvar_aa_dict_pkl: String,
    pub(crate) clinvar_splice_dict_pkl: String,
    pub(crate) interpro_entry_map_pkl: String,
    pub(crate) intolerant_domains_pkl: String,
    pub(crate) gene_dosage_sensitivity: String,
    pub(crate) am_intol_domains_tsv: String,
    pub(crate) pm1_regions_pkl: String,
    pub(crate) clinvar_gene_stat_pkl: String,
    pub(crate) tranx_exon_domain_map_pkl: String,
    /// Accepted and never read, as in Python.
    pub(crate) am_score_vcf: Option<String>,
    pub(crate) alt_disease_vcf: Option<String>,
    pub(crate) clingen_map_pkl: Option<String>,
    pub(crate) ignore_clingen_override: bool,
    pub(crate) mavedb_metadata_tsv: Option<String>,
    /// PP1's multi-family VCF and pedigree (default empty), forwarded to PP1_criteria.
    pub(crate) pp1_vcf: String,
    pub(crate) pp1_ped: String,
    pub(crate) cds_fasta_path: Option<String>,
    pub(crate) pext_tissues: String,
    pub(crate) pext_low_expression_cutoff: f64,
    pub(crate) pext_penalty_floor: f64,
    pub(crate) pext_penalty_shape: f64,
    pub(crate) repeat_region_file: String,
    pub(crate) gnomad_extreme_rare_threshold: f64,
    pub(crate) expected_incidence: f64,
    pub(crate) relevant_gene_list: Option<String>,
    pub(crate) dispensable_gene_list: Option<String>,
    pub(crate) hpo_condition_mechanism_json: String,
    pub(crate) gene_nonlof_mechanism_json: String,
    pub(crate) threads: i64,
    /// The Python step's opt-in Rust ranking; meaningless once the step is Rust.
    pub(crate) ranking_binary: Option<String>,
    pub(crate) config: Option<String>,
    /// Rust only (`--clinvar_pm5_index`): the PM5 index of confident pathogenic ClinVar protein
    /// changes that the owner's PM5 rules read (`acmg_pm5_residue_rules`).
    pub(crate) clinvar_pm5_index: Option<String>,
}

/// The values `infer_query_variant_effect` (gene_mechanism_variants.py, thread
/// T3) reads from one table row, for the hub (thread T2). `None` means the
/// column is absent (`row.get(name)` gives `None`); `Some(Cell::Missing)` is NaN.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct EffectInputs<'a> {
    pub(crate) lof: Option<Cell<'a>>,
    pub(crate) consequence: Option<Cell<'a>>,
    pub(crate) nmd: Option<Cell<'a>>,
    pub(crate) lof_filter: Option<Cell<'a>>,
    pub(crate) splicing_frameshift: Option<Cell<'a>>,
    pub(crate) vep_consq_lof: Option<Cell<'a>>,
    pub(crate) splicing_lof: Option<Cell<'a>>,
    /// The `5UTR_lof` column.
    pub(crate) utr5_lof: Option<Cell<'a>>,
    pub(crate) variant_lof_score: Option<Cell<'a>>,
    pub(crate) variant_gof_score: Option<Cell<'a>>,
    pub(crate) variant_dn_score: Option<Cell<'a>>,
    pub(crate) variant_gof_tag: Option<Cell<'a>>,
}

/// The dict `infer_query_variant_effect` returns.
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct QueryEffect {
    pub(crate) variant_effect: String,
    pub(crate) variant_effect_evidence: String,
    pub(crate) variant_effect_suppressed_evidence: String,
    pub(crate) variant_lof_score: i64,
    pub(crate) variant_gof_score: i64,
    pub(crate) variant_dn_score: i64,
    pub(crate) variant_mechanism_exclusive: bool,
    pub(crate) variant_exact_mechanisms: String,
}

/// The HGNC symbol resolver of `GeneMechanismHub` (gene_mechanism_resources.py,
/// thread T2), as the exact-allele matcher (thread T3) uses it.
pub(crate) trait SymbolResolver {
    /// `hub.resolve_symbol(gene_symbol)`: one current HGNC symbol, or the error Python raises.
    fn resolve_symbol(&self, gene_symbol: Cell<'_>) -> Result<String>;
    /// `ExactVariantMatcher._try_resolve_symbol`: the resolved symbol, or the
    /// cleaned input when resolution fails.
    fn try_resolve_symbol(&self, gene_symbol: Cell<'_>) -> String;
}

/// One call of `priva acmg-step <module> <function>`: the input frame and
/// the options (a JSON object written by the Python capture script).
#[derive(Clone, Debug)]
pub(crate) struct HarnessCall {
    pub(crate) table: Table,
    pub(crate) options: Json,
}

/// The dev harness entry every module provides: run `function` of the module
/// on the call and return the frame the Python capture wrote for it.
pub(crate) type HarnessEntry = fn(&str, &HarnessCall) -> Result<Table>;
