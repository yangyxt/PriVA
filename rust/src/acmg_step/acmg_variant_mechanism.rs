//! `scripts/acmg_variant_mechanism.py` (591 lines): the exact-allele columns and the adapters to
//! the existing mask modules.
//!
//! Thread T3 (exact-allele matching and the scoring tail) fills this module. Python to convert:
//! _clean_text, _row_assembly_for_gof_lookup, _merge_semicolon_values,
//! _ensure_exact_nonlof_columns, _record_exact_nonlof_match, annotate_exact_nonlof_variants; and
//! the two adapters below to `crate::acmg_variant_mechanism` (the masks already in Rust).
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, NamedMasks, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step acmg_variant_mechanism <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("acmg_variant_mechanism", function))
}

/// `_variant_condition_masks(df)`: the 13 condition masks from
/// `variant_inheritance` and `variant_penetrance`, by their Python names. Used
/// by PVS1 (T4), PM2/BS1 and BS2/BP5 (T5) and PP1/BP2/PM3 (T7); this thread (T3)
/// fills it with an adapter to `crate::acmg_variant_mechanism::variant_condition_masks`.
pub(crate) fn variant_condition_masks(table: &Table) -> Result<NamedMasks> {
    let _ = table;
    Err(not_converted(
        "acmg_variant_mechanism",
        "_variant_condition_masks",
    ))
}

/// `_variant_mechanism_masks(df)`: the 27 mechanism masks, by their Python names;
/// an adapter to `crate::acmg_variant_mechanism::variant_mechanism_masks`.
pub(crate) fn variant_mechanism_masks(table: &Table) -> Result<NamedMasks> {
    let _ = table;
    Err(not_converted(
        "acmg_variant_mechanism",
        "_variant_mechanism_masks",
    ))
}
