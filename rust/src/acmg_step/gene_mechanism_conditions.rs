//! `scripts/gene_mechanism_conditions.py` (1,060 lines): the condition cache: assertions,
//! inheritance and penetrance normalisation and histories.
//!
//! Thread T2 (gene-condition hub) fills this module. Python to convert: the functions on the step's
//! path (about 820 of 1,060 lines): condition_cache_context, condition_cache_mechanism_entries,
//! condition_cache_mechanism_assertions, normalize_inheritance, normalize_inheritance_histories,
//! gene_has_lof_mechanism_history, gene_inheritance_from_constraint, normalize_penetrance and their
//! helpers.
//!
//! Not needed by the step (group 4): build_hpo_condition_index,
//! enrich_condition_mechanism_assertion, gene_inheritance_consensus.
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step gene_mechanism_conditions <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("gene_mechanism_conditions", function))
}
