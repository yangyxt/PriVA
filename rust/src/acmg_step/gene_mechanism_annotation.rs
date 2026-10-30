//! `scripts/gene_mechanism_annotation.py` (411 lines): the per-row chain that writes the hub's 16
//! columns.
//!
//! Thread T2 (gene-condition hub) fills this module. Python to convert:
//! annotate_gene_mechanism_categories, select_condition_histories_for_variant,
//! _classify_variant_applicability, _compact_inheritance, _mechanism_profile_tag,
//! _compact_profile_tags.
//!
//! Calls `gene_mechanism_variants::infer_query_variant_effect` (T3) through the fixed signature
//! with `types::EffectInputs`; test against a stub until T3 lands.
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step gene_mechanism_annotation <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("gene_mechanism_annotation", function))
}
