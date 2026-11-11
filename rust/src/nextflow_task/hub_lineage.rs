//! The annotation hub folder name (lineage) of each annotation step: a hash
//! of the step's name, the text of the annotation script it runs, and the
//! resource records of every resource that step reads. A change to any of
//! them gives a new, empty hub folder, so stored results are never reused
//! with other resources or code. The control cohort step has no folder
//! ("off") when its counts come from the cohort itself (no control VCF).
use super::Result;
use super::contract::{json_digest, object, text};
use crate::json_text::Json;
use std::path::Path;

pub(crate) const SOURCES: [&str; 5] = ["gnomad", "clinvar", "control", "vep", "cadd"];

/// The resources the VEP step reads: its cache, plugins and their data.
const VEP_KEYS: &[&str] = &[
    "ref_genome",
    "reference_fai",
    "vep_cache_dir",
    "vep_plugins_dir",
    "vep_plugins_cachedir",
    "utr_annotator_file",
    "alphamissense_prescore",
    "spliceai_snv_prescore",
    "spliceai_indel_prescore",
    "primateai_prescore",
    "conservation_file",
    "splicevault_prescore",
    "loftee_repo",
    "human_ancestor_fasta",
    "loftee_conservation_file",
    "gerp_bigwig",
    "loeuf_used_by_shell",
    "mavedb_file",
    "mavedb_assay_metadata_used_by_shell",
    "pext_mean_bw",
    "pext_tissue_bw_dir",
];

fn reads(source: &str, key: &str) -> bool {
    match source {
        "gnomad" => key == "gnomad_vcf_chrX" || key.starts_with("gnomad_chr"),
        "clinvar" => key == "clinvar_vcf",
        "control" => key == "control_vcf",
        "vep" => VEP_KEYS.contains(&key),
        "cadd" => key == "cadd_base_dir" || key == "ref_genome",
        _ => false,
    }
}

/// {source: lineage} for the five annotation steps.
pub(crate) fn lineages(manifest: &Json, root: &Path) -> Result<Json> {
    let mut code = Vec::new();
    for name in ["annotation_vcf.sh", "common_bash_utils.sh"] {
        let path = root.join("scripts").join(name);
        code.push(text(if path.is_file() {
            crate::sha256::sha256_file(&path)?
        } else {
            String::new()
        }));
    }
    let resources: Vec<(String, Json)> = match manifest.get("resources") {
        Some(Json::Object(items)) => items.clone(),
        _ => Vec::new(),
    };
    let mut out = Vec::new();
    for source in SOURCES {
        let used: Vec<(String, Json)> = resources
            .iter()
            .filter(|(key, _)| reads(source, key.split('/').next().unwrap_or(key)))
            .cloned()
            .collect();
        let lineage = if source == "control" && used.is_empty() {
            "off".to_owned()
        } else {
            let payload = object(vec![
                ("source", text(source)),
                ("code", Json::Array(code.clone())),
                ("resources", Json::Object(used)),
            ]);
            json_digest(&payload)?[..16].to_owned()
        };
        out.push((source.to_owned(), text(lineage)));
    }
    Ok(Json::Object(out))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn manifest(entries: &[(&str, i128)]) -> Json {
        Json::Object(vec![(
            "resources".into(),
            Json::Object(
                entries
                    .iter()
                    .map(|(k, size)| {
                        let record = object(vec![("size", Json::Int(*size))]);
                        ((*k).to_owned(), record)
                    })
                    .collect(),
            ),
        )])
    }

    fn lineage(json: &Json, source: &str) -> String {
        json.get(source).and_then(Json::as_str).unwrap().to_owned()
    }

    #[test]
    fn a_resource_change_moves_only_the_steps_that_read_it() {
        let root = std::env::temp_dir();
        let of = |clinvar: i128, vep: i128| {
            let entries = [
                ("clinvar_vcf", clinvar),
                ("vep_cache_dir/a", vep),
                ("control_vcf", 3),
            ];
            lineages(&manifest(&entries), &root).unwrap()
        };
        let (before, after) = (of(1, 2), of(9, 2));
        assert_ne!(lineage(&before, "clinvar"), lineage(&after, "clinvar"));
        for source in ["gnomad", "control", "vep", "cadd"] {
            assert_eq!(
                lineage(&before, source),
                lineage(&after, source),
                "{source}"
            );
        }
        let vep_changed = of(1, 5);
        assert_ne!(lineage(&before, "vep"), lineage(&vep_changed, "vep"));
        assert_eq!(
            lineage(&before, "clinvar"),
            lineage(&vep_changed, "clinvar")
        );
    }

    #[test]
    fn control_counts_from_the_cohort_have_no_hub_folder() {
        let json = lineages(&manifest(&[("clinvar_vcf", 1)]), &std::env::temp_dir()).unwrap();
        assert_eq!(lineage(&json, "control"), "off");
        assert_eq!(lineage(&json, "clinvar").len(), 16);
    }
}
