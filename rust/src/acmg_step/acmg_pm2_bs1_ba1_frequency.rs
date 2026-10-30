//! `scripts/acmg_pm2_bs1_ba1_frequency.py` (799 lines): PM2, BS1 and BA1 from gnomAD and the
//! ClinVar allele-frequency statistics.
//!
//! Thread T5 (allele frequency, carrier observation and statistics) fills this module. Python to
//! convert: _gene_max_patho_hom_freq, control_false_neg_rate, control_false_positive_rate,
//! gnomAD_rare_AF, _rare_af_binomial, _rare_hom_binomial, PM2_criteria, BS1_criteria, BA1_criteria.
//!
//! Calls `stats::binom_sf` and `stats::binom_cdf` (this thread) and
//! `acmg_variant_mechanism::variant_condition_masks` (T3).
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step acmg_pm2_bs1_ba1_frequency <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("acmg_pm2_bs1_ba1_frequency", function))
}
