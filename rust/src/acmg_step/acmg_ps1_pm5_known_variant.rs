//! `scripts/acmg_ps1_pm5_known_variant.py` (288 lines): PS1 and PM5 from the ClinVar amino-acid and
//! splice dictionaries.
//!
//! Thread T4 (null variants and known variants) fills this module. Python to convert:
//! check_splice_pathogenic, PS1_PM5_criteria.
//!
//! Uses `splicing_var_analysis::parse_hgvsc_splice_position` (group 2). The two ClinVar pickles are
//! large (see `rust/GROUP3_BASE.md`, probe results).
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step acmg_ps1_pm5_known_variant <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("acmg_ps1_pm5_known_variant", function))
}
