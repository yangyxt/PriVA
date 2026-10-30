//! `scripts/acmg_pp3_bp4_bp7_insilico.py` (139 lines): PP3, BP4 and BP7 from the in-silico
//! predictors.
//!
//! Thread T6 (missense, predictors, protein length and reported evidence) fills this module. Python
//! to convert: PP3_BP4_criteria, BP7_criteria.
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step acmg_pp3_bp4_bp7_insilico <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("acmg_pp3_bp4_bp7_insilico", function))
}
