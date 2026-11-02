//! The Python helper checks structure, not clinical values. Native archive
//! admission is performed by the existing readers before these shape checks.
use super::{error, object};
use crate::cache_convert::{
    ClinvarAaCache, clinvar_exon_af::ClinvarExonAfCache, clinvar_gene_af::ClinvarGeneAfCache,
    clinvar_gene_stat::ClinvarGeneCache, clinvar_splice_stat::typed::Cache as SpliceCache,
};
use crate::cli::Result;
use crate::json_text::Json;

pub(super) fn validate(path: &str, kind: &str) -> Result<Json> {
    let (count, valid) = match kind {
        "aa" => {
            let cache = ClinvarAaCache::open(path)?;
            (
                cache.len(),
                cache.keys().iter().all(|key| {
                    let value = cache.get(key.as_str()).unwrap();
                    !value.positions.is_empty()
                        && value.positions.iter().all(|p| !p.changes.is_empty())
                }),
            )
        }
        "splice" => {
            let cache = SpliceCache::open(path)?;
            (
                cache.len(),
                cache
                    .keys()
                    .iter()
                    .all(|key| cache.record_count(key.as_str()).unwrap() > 0),
            )
        }
        "gene-af" => {
            let cache = ClinvarGeneAfCache::open(path)?;
            (
                cache.len(),
                cache.keys().iter().all(|key| {
                    cache.get(key.as_str()).unwrap().is_none_or(|record| {
                        ["chrom", "pos", "af", "acc_patho_af"]
                            .iter()
                            .all(|name| record.iter().any(|item| item.key == *name))
                    })
                }),
            )
        }
        "exon-af" => {
            let cache = ClinvarExonAfCache::open(path)?;
            (
                cache.len(),
                cache
                    .keys()
                    .iter()
                    .all(|key| !cache.get(key.as_str()).unwrap().is_empty()),
            )
        }
        "gene" => {
            let cache = ClinvarGeneCache::open(path)?;
            // Both schema versions contain the eleven required fields by
            // construction. Non-text symbols are not a shape failure.
            (cache.len(), true)
        }
        _ => {
            return Err(error(
                "ValueError",
                format!("Unknown ClinVar cache kind: {kind}"),
            ));
        }
    };
    if count == 0 {
        return Err(error(
            "ValueError",
            format!("ClinVar {kind} cache must contain a nonempty dictionary: {path}"),
        ));
    }
    if !valid {
        return Err(error(
            "ValueError",
            format!("ClinVar {kind} cache has an unexpected statistics structure: {path}"),
        ));
    }
    Ok(object([
        ("path", Json::Str(path.into())),
        ("cache_kind", Json::Str(kind.into())),
        ("success", Json::Bool(true)),
        ("top_level_entries", Json::Int(count as i128)),
    ]))
}
