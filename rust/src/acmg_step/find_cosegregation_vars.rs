//! `scripts/find_cosegregation_vars.py` (228 lines): PP1 co-segregation counts from a multi-family
//! VCF.
//!
//! Thread T7 (family evidence) fills this module. Python to convert: parse_gt, _normalized_sex,
//! _normalized_chromosome, check_variant_presence, find_cosegregating_variants.
//!
//! Reads the VCF the way pysam does; reuse `vcf_text` and the genotype handling of
//! `pedigree_filter`.
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step find_cosegregation_vars <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("find_cosegregation_vars", function))
}
