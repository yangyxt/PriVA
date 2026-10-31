//! The ACMG decision step (`scripts/acmg_criteria_assign.py` and the scripts it
//! imports), converted in group 3 by seven parallel threads. One module per
//! Python script, named like the script; `crate::acmg_consequence`,
//! `crate::acmg_scoring` and `crate::acmg_variant_mechanism` are the earlier
//! opt-in modules, which the new code calls and does not change.
//!
//! Shared by every thread: `types` (fixed signatures), `stats` (SciPy's four
//! functions), `frame_file` and `dev_harness` (`priva acmg-step`, which runs one
//! converted function on a pandas frame written by the Python capture script),
//! `resource_summary` (`priva acmg-step inspect`). The table layer is
//! `crate::pandas_table` and `crate::pandas_polars`. `rust/GROUP3_BASE.md`
//! explains the design and what each thread starts with.
// The fixed signatures in `types`, `stats` and the module stubs have their
// callers in the threads' later commits.
#![allow(dead_code)]
pub(crate) mod dev_harness;
pub(crate) mod frame_file;
pub(crate) mod resource_summary;
pub(crate) mod stats;
pub(crate) mod types;

// T1: table, command line and orchestrator.
pub(crate) mod acmg_consequence;
pub(crate) mod acmg_criteria_assign;
// T2: gene-condition hub.
pub(crate) mod acmg_inheritance;
pub(crate) mod clinvar_vcv;
pub(crate) mod gene_mechanism_annotation;
pub(crate) mod gene_mechanism_common;
pub(crate) mod gene_mechanism_conditions;
pub(crate) mod gene_mechanism_hub;
pub(crate) mod gene_mechanism_resources;
// T3: exact-allele matching and the scoring tail.
pub(crate) mod acmg_scoring;
pub(crate) mod acmg_variant_mechanism;
pub(crate) mod gene_mechanism_variants;
// T4: null variants and known variants.
pub(crate) mod acmg_ps1_pm5_known_variant;
pub(crate) mod acmg_pvs1_null_variant;
pub(crate) mod alternative_start_codon;
pub(crate) mod pickle_dict_stream;
// T5: allele frequency, carrier observation and statistics (and `stats`).
pub(crate) mod acmg_bs2_bp5_observation;
pub(crate) mod acmg_pm2_bs1_ba1_frequency;
// T6: missense, predictors, protein length and reported evidence.
pub(crate) mod acmg_pm1_pp2_bp1_missense;
pub(crate) mod acmg_pm4_bp3_protein_length;
pub(crate) mod acmg_pp3_bp4_bp7_insilico;
pub(crate) mod acmg_ps3_bs3_pp5_bp6_reported;
pub(crate) mod mavedb_interpreter;
// T7: family evidence.
pub(crate) mod acmg_pp1_bs4_bp2_pm3_family;
pub(crate) mod acmg_ps2_pm6_denovo;
pub(crate) mod determine_phase;
pub(crate) mod find_cosegregation_vars;
