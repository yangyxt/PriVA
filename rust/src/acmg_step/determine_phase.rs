//! `scripts/determine_phase.py` (521 lines): cis and trans calls for pairs of variants in a gene.
//!
//! Thread T7 (family evidence) fills this module. Python to convert: Genotype, Variant,
//! create_variants_from_gene_df, analyze_gene_variants, process_gene_wrapper, build_patient_info,
//! _can_resolve_phase, determine_cis_trans_relationships.
//!
//! The per-gene results are merged in the genes' groupby order (owner's decision 2026-10-30, as
//! main's Python now does with `pool.imap`): use `pandas_polars::first_appearance_groups`.
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step determine_phase <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("determine_phase", function))
}
