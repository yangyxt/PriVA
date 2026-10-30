//! `scripts/gene_mechanism_hub.py` (541 lines): the facade that re-exports the hub.
//!
//! Thread T2 (gene-condition hub) fills this module. Python to convert: the 71-line facade
//! (resolve_gene_symbol and the match_* wrappers); its `main` diagnostic belongs to group 4.
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step gene_mechanism_hub <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("gene_mechanism_hub", function))
}
