//! `scripts/mavedb_interpreter.py` (842 lines): the assay-family rules that turn a MaveDB score
//! into loss, gain or neutral.
//!
//! Thread T6 (missense, predictors, protein length and reported evidence) fills this module. Python
//! to convert: the enums, AssaySpec, _to_log2, validate_spec, BaseInterpreter and the 13
//! interpreters, interpret_score, parse_neutral_band_from_metadata, pick_interpreter_from_metadata.
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step mavedb_interpreter <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("mavedb_interpreter", function))
}
