//! `scripts/gene_mechanism_variants.py` (642 lines): `infer_query_variant_effect` and the exact
//! curated-allele matcher.
//!
//! Thread T3 (exact-allele matching and the scoring tail) fills this module. Python to convert:
//! infer_query_variant_effect, ExactVariantMatcher (match_curated_nonlof_variant,
//! match_gofcards_variant_gof, match_gofcards_variant_gof_by_genomic and loaders).
//!
//! Until then it holds only the dev harness entry, which reports the function as not converted.
use super::types::{EffectInputs, HarnessCall, QueryEffect, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;

/// `priva acmg-step gene_mechanism_variants <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("gene_mechanism_variants", function))
}

/// `infer_query_variant_effect(row)` for one table row. The hub (T2) calls it
/// for every row; this thread (T3) fills it.
pub(crate) fn infer_query_variant_effect(row: &EffectInputs<'_>) -> Result<QueryEffect> {
    let _ = row;
    Err(not_converted(
        "gene_mechanism_variants",
        "infer_query_variant_effect",
    ))
}
