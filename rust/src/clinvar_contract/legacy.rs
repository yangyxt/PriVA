//! `validate_cache(path, kind)` for the Python pickles of an existing
//! installation, which stay reusable for one transition. The restricted
//! pickle reader replaces `pickle.load`; a pickle it cannot read (one holding
//! numpy objects, for example) fails the check, so the stage rebuilds.
use super::{error, object};
use crate::cli::Result;
use crate::json_text::Json;
use crate::pickle_reader::{self, Dict, Value};
use std::path::Path;

/// `isinstance(value, dict) and value`, with the dict.
fn nonempty_dict(value: &Value, test: impl Fn(&Dict) -> bool) -> bool {
    match value {
        Value::Dict(dict) => {
            let dict = dict.borrow();
            !dict.entries.is_empty() && test(&dict)
        }
        _ => false,
    }
}

/// `isinstance(value, dict) and set(names) <= value.keys()`.
fn has_keys(value: &Value, names: &[&str]) -> bool {
    match value {
        Value::Dict(dict) => {
            let dict = dict.borrow();
            names
                .iter()
                .all(|name| matches!(dict.get(&Value::Str((*name).into())), Ok(Some(_))))
        }
        _ => false,
    }
}

fn shape(kind: &str, data: &Dict) -> Result<bool> {
    let values = || data.entries.iter().map(|(_, value)| value);
    Ok(match kind {
        "aa" => values().all(|positions| {
            nonempty_dict(positions, |positions| {
                positions.entries.iter().all(|(_, changes)| {
                    nonempty_dict(changes, |changes| {
                        changes
                            .entries
                            .iter()
                            .all(|(_, entry)| has_keys(entry, &["CLNSIG", "CLNREVSTAT"]))
                    })
                })
            })
        }),
        "splice" => values().all(|entries| match entries {
            Value::List(entries) => {
                let entries = entries.borrow();
                !entries.is_empty()
                    && entries.iter().all(|entry| {
                        has_keys(entry, &["chrom", "pos", "clinvar_sig", "splice_ai"])
                    })
            }
            _ => false,
        }),
        // The original maximum-AF rule can retain None when every AF is zero.
        "gene-af" => values().all(|value| {
            matches!(value, Value::None) || has_keys(value, &["chrom", "pos", "af", "acc_patho_af"])
        }),
        "exon-af" => values().all(|exons| {
            nonempty_dict(exons, |exons| {
                exons
                    .entries
                    .iter()
                    .all(|(_, value)| matches!(value, Value::Tuple(items) if items.len() == 5))
            })
        }),
        "gene" => values().all(|value| {
            has_keys(
                value,
                &[
                    "gene_symbol",
                    "ensembl_id",
                    "clinvar_variants",
                    "clinvar_pathogenic",
                    "clinvar_not_pathogenic",
                    "putative_nmd_variants",
                    "small_aachange",
                    "large_aachange",
                    "no_aachange",
                    "transcripts",
                    "domains",
                ],
            )
        }),
        _ => {
            return Err(error(
                "ValueError",
                format!("Unknown ClinVar cache kind: {kind}"),
            ));
        }
    })
}

pub(super) fn validate(path: &str, kind: &str) -> Result<Json> {
    // The script opens with gzip by the `.gz` name, not by the file's bytes.
    let bytes = pickle_reader::read_file(Path::new(path), path.ends_with(".gz"))?;
    let data = if kind == "gene" {
        pickle_reader::load_clinvar_gene(&bytes)?
    } else {
        pickle_reader::load(&bytes)?
    };
    let Value::Dict(data) = &data else {
        return Err(error(
            "ValueError",
            format!("ClinVar {kind} cache must contain a nonempty dictionary: {path}"),
        ));
    };
    let data = data.borrow();
    if data.entries.is_empty() {
        return Err(error(
            "ValueError",
            format!("ClinVar {kind} cache must contain a nonempty dictionary: {path}"),
        ));
    }
    if !shape(kind, &data)? {
        return Err(error(
            "ValueError",
            format!("ClinVar {kind} cache has an unexpected statistics structure: {path}"),
        ));
    }
    Ok(object([
        ("path", Json::Str(path.into())),
        ("cache_kind", Json::Str(kind.into())),
        ("success", Json::Bool(true)),
        ("top_level_entries", Json::Int(data.entries.len() as i128)),
    ]))
}
