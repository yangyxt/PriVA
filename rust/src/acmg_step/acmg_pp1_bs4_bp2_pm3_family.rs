//! `scripts/acmg_pp1_bs4_bp2_pm3_family.py` (602 lines): PP1, BS4, BP2 and PM3.
//!
//! Thread T7 (family evidence) fills this module. Python to convert: PP1_criteria, BS4_criteria,
//! BP2_PM3_criteria.
//!
//! PP1 is reachable from the command line now (`--pp1_vcf`, `--pp1_ped`).
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step acmg_pp1_bs4_bp2_pm3_family <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("acmg_pp1_bs4_bp2_pm3_family", function))
}
