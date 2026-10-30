//! `scripts/acmg_pm4_bp3_protein_length.py` (224 lines): PM4 and BP3 and the repeat-region overlap.
//!
//! Thread T6 (missense, predictors, protein length and reported evidence) fills this module. Python
//! to convert: find_overlaps_bedtools_efficient, PM4_BP3_criteria.
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step acmg_pm4_bp3_protein_length <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("acmg_pm4_bp3_protein_length", function))
}
