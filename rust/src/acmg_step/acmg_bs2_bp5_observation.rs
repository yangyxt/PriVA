//! `scripts/acmg_bs2_bp5_observation.py` (460 lines): BS2 and BP5 (BP5 runs `bcftools` through a
//! shell pipeline).
//!
//! Thread T5 (allele frequency, carrier observation and statistics) fills this module. Python to
//! convert: BS2_criteria, extract_and_summarize_alt_disease_variants, BP5_criteria.
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step acmg_bs2_bp5_observation <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("acmg_bs2_bp5_observation", function))
}
