//! `scripts/acmg_ps3_bs3_pp5_bp6_reported.py` (258 lines): PS3, BS3, PP5 and BP6 from ClinVar and
//! MaveDB.
//!
//! Thread T6 (missense, predictors, protein length and reported evidence) fills this module. Python
//! to convert: mavedb_interpretation_per_row, mavedb_score_interpretation, PS3_BS3_criteria,
//! PP5_BP6_criteria.
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step acmg_ps3_bs3_pp5_bp6_reported <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("acmg_ps3_bs3_pp5_bp6_reported", function))
}
