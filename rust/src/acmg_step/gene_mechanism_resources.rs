//! `scripts/gene_mechanism_resources.py` (1,028 lines): `GeneMechanismHub` and the local HGNC
//! resolver.
//!
//! Thread T2 (gene-condition hub) fills this module. Python to convert: _LocalHgncResolver,
//! GeneMechanismHub (constructor with use_hgnc_package=False, _load_condition_cache,
//! _load_mechanisms, _load_clingen, _load_loeuf, resolve_symbol).
//!
//! The hub implements `types::SymbolResolver`, which the exact-allele matcher (T3) calls. Its
//! default data files (HGNC, ClinGen hg19 dosage, LOEUF, dispensable genes) come from
//! `protein_domain_mapping::data_directory()`. The TSVs are read with
//! `pandas_table::read_table_with` (`all_text`, `use_columns`); a missing value then becomes `None`
//! and the Python `.fillna("")` is `unwrap_or_default`.
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step gene_mechanism_resources <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("gene_mechanism_resources", function))
}
