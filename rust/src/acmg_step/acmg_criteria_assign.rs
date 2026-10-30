//! `scripts/acmg_criteria_assign.py` (717 lines): the orchestrator: argparse, reading the inputs,
//! the call sequence of `ACMG_criteria_assign`, the criteria matrix `<table>.acmg.tsv` and the
//! rewritten table.
//!
//! Thread T1 (table, command line and orchestrator) fills this module. Python to convert:
//! ACMG_criteria_assign and the `__main__` argparse block.
//!
//! It owns `priva acmg-assign` (registered in `cli.rs`; `dispatch` below) with the Python option
//! names, `--pp1_vcf` and `--pp1_ped` included and forwarded to PP1_criteria (fixed on main in
//! c529805). The parsed options go into `types::StepArguments`. The table is read with
//! `pandas_table::read_table` and de-duplicated with `pandas_polars::drop_duplicates_first`; the
//! PED with `ReadOptions { all_text: true, keep_default_na: false, .. }` and the headerless
//! fallback (`Header::Names`, `UseColumns::Positions(0..6)`). Each criterion function returns a
//! `types::StepOutput`; this module writes its columns in the Python order. On a failure after the
//! matrix is written, the matrix stays and the input table is untouched.
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step acmg_criteria_assign <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("acmg_criteria_assign", function))
}

/// `priva acmg-assign`: the converted `acmg_criteria_assign.py` command line.
pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let _ = args;
    Err(not_converted("acmg_criteria_assign", "__main__"))
}
