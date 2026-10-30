//! `scripts/acmg_pm1_pp2_bp1_missense.py` (307 lines): PM1, PP2 and BP1.
//!
//! Thread T6 (missense, predictors, protein length and reported evidence) fills this module. Python
//! to convert: locate_pm1_region, PM1_criteria, analyze_bp1_pp2, PP2_BP1_criteria.
//!
//! Calls `stats::chi2_contingency_2x2` and `stats::fisher_exact_greater` (T5).
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step acmg_pm1_pp2_bp1_missense <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("acmg_pm1_pp2_bp1_missense", function))
}
