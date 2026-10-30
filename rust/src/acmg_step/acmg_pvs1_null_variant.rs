//! `scripts/acmg_pvs1_null_variant.py` (721 lines): PVS1 and its gene-level gate.
//!
//! Thread T4 (null variants and known variants) fills this module. Python to convert:
//! summarize_clinvar_gene_pathogenicity, identify_alternative_start_codon_genes,
//! evaluate_start_lost_with_cds, downstream_domain_impact, downstream_exon_patho_af,
//! span_functional_domains, load_non_lof_gene_set, PVS1_criteria.
//!
//! Uses `protein_domain_mapping::DomainNormalizer` and `interpro_mapping_io::load_interpro_mapping`
//! (group 2), `acmg_variant_mechanism::variant_mechanism_masks` (T3) and the hub's
//! `gene_lof_mechanism_history` column.
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step acmg_pvs1_null_variant <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("acmg_pvs1_null_variant", function))
}
