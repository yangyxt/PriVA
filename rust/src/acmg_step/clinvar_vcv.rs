//! `scripts/clinvar_vcv.py` (1,061 lines): the ClinVar VCV helpers; the step uses only `open_text`.
//!
//! Thread T2 (gene-condition hub) fills this module. Python to convert: open_text (26 lines); the
//! XML parser stays for group 4.
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step clinvar_vcv <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("clinvar_vcv", function))
}
