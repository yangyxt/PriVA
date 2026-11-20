//! `scripts/gene_mechanism_hub.py` (541 lines): the facade that re-exports the hub.
//!
//! Thread T2 (gene-condition hub). The ACMG step imports from this facade only names defined in
//! the other hub modules, which are converted there: `GeneMechanismHub` and the HGNC resolver in
//! `gene_mechanism_resources`, `annotate_gene_mechanism_categories` in
//! `gene_mechanism_annotation`, the condition functions in `gene_mechanism_conditions`, the
//! constants and default files in `gene_mechanism_common`. Its own convenience functions
//! (`resolve_gene_symbol`, which builds the hub with the optional external HGNC package, and the
//! `match_*` wrappers) and its `main` diagnostic belong to group 4.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step gene_mechanism_hub <function>`: the facade has no function of its own here.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("gene_mechanism_hub", function))
}
