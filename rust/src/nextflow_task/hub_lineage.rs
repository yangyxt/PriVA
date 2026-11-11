//! The annotation hub folder name (lineage) of each annotation step: a hash
//! of the step's name, the text of the annotation script it runs, and the
//! identity of every resource that step reads. A resource's identity leaves
//! its path out, so moving the checkout or the caches keeps the folder: a file
//! inside the PriVA checkout (bundled data) counts by its content, any other
//! by its size and recorded SHA-256, or its size and modification time. A change to any of
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

/// A resource record without its path: content digest for files inside the
/// checkout, else size with the recorded digest or the modification time.
fn identity(record: &Json, root: &Path) -> Result<Json> {
    let field = |name: &str| record.get(name).cloned().unwrap_or(Json::Null);
    let resolved = record
        .get("resolved_path")
        .and_then(Json::as_str)
        .map(Path::new);
    let canonical_root = root.canonicalize().unwrap_or_else(|_| root.to_path_buf());
    if let Some(path) = resolved.filter(|p| p.starts_with(&canonical_root) && p.is_file()) {
        return Ok(object(vec![(
            "sha256",
            text(crate::sha256::sha256_file(path)?),
        )]));
    }
    let digest = field("sha256");
    let marker = if matches!(digest, Json::Null) {
        field("mtime_ns")
    } else {
        digest
    };
    Ok(object(vec![
        ("type", field("type")),
        ("size", field("size")),
        ("marker", marker),
    ]))
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
        let mut used = Vec::new();
        for (key, record) in &resources {
            if reads(source, key.split('/').next().unwrap_or(key)) {
                used.push((key.clone(), identity(record, root)?));
            }
        }
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
    fn moving_a_resource_keeps_its_lineage() {
        let record = |path: &str| {
            object(vec![
                ("path", text(path)),
                ("resolved_path", text(path)),
                ("size", Json::Int(5)),
                ("mtime_ns", Json::Int(7)),
            ])
        };
        let at = |path: &str| {
            let resources = Json::Object(vec![("clinvar_vcf".into(), record(path))]);
            let manifest = Json::Object(vec![("resources".into(), resources)]);
            lineages(&manifest, &std::env::temp_dir()).unwrap()
        };
        assert_eq!(
            lineage(&at("/old/clinvar.vcf.gz"), "clinvar"),
            lineage(&at("/new/c.vcf.gz"), "clinvar")
        );
    }

    #[test]
    fn control_counts_from_the_cohort_have_no_hub_folder() {
        let json = lineages(&manifest(&[("clinvar_vcf", 1)]), &std::env::temp_dir()).unwrap();
        assert_eq!(lineage(&json, "control"), "off");
        assert_eq!(lineage(&json, "clinvar").len(), 16);
    }
}
