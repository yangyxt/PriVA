//! `scripts/gene_mechanism_common.py` (431 lines): constants and text helpers of the hub.
//!
//! Thread T2 (gene-condition hub) fills this module. Python to convert: _clean, _bool_value,
//! _split_annotation_tokens, _bounded_mechanism_score, _exact_mechanisms_from_effect, the _norm*
//! key helpers, _split_multi, _split_pmids, _safe_int, _safe_float, _deduplicate_by.
//!
//! Uses `crate::hpo_penetrance` (group 1).
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step gene_mechanism_common <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("gene_mechanism_common", function))
}
