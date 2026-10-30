//! `scripts/acmg_scoring.py` (836 lines): the scoring tail: criteria matrix and ClinGen override,
//! `compute_control_common`, and the glue that runs the existing posterior and ranking modules in
//! the Python order.
//!
//! Thread T3 (exact-allele matching and the scoring tail) fills this module. Python to convert:
//! summarize_acmg_criteria (matrix construction and validation, ClinGen evidence lookup),
//! compute_control_common, sort_and_rank_variants; it reuses `crate::acmg_scoring` (summary
//! strings, posterior, ranking stages R1 to R8).
//!
//! Calls `stats::binom_cdf` and `stats::beta_ppf` (T5).
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step acmg_scoring <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("acmg_scoring", function))
}
