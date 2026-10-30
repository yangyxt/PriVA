//! `scripts/acmg_ps2_pm6_denovo.py` (215 lines): PS2 and PM6 from the family genotypes.
//!
//! Thread T7 (family evidence) fills this module. Python to convert: identify_fam_members,
//! determine_denovo_per_row, PS2_PM6_criteria.
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step acmg_ps2_pm6_denovo <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("acmg_ps2_pm6_denovo", function))
}
