//! `scripts/alternative_start_codon.py` (405 lines): the CDS FASTA index and the downstream-ATG
//! search for start_lost variants.
//!
//! Thread T4 (null variants and known variants) fills this module. Python to convert:
//! AlternativeStartCodonDetector, adjust_pvs1_for_start_lost,
//! identify_transcripts_with_alternative_starts, get_detector.
//!
//! The `.cds_cache.pkl.gz` next to the FASTA: see the group 3 plan, section 8.8.
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step alternative_start_codon <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("alternative_start_codon", function))
}
