//! `scripts/acmg_consequence.py` (196 lines): the consequence helpers not yet in Rust: the LoF and
//! length-changing flags per row, the HGVSp variant type and the truncated fraction
//! (`coerce_pext_value` is already in `crate::acmg_consequence`).
//!
//! Thread T1 (table, command line and orchestrator) fills this module. Python to convert:
//! vep_consq_interpret_per_row, vep_consq_interpret, get_variant_type, truncate_fraction.
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step acmg_consequence <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("acmg_consequence", function))
}
