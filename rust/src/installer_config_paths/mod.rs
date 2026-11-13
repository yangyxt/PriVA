//! Source-mapped installer path planning and byte adoption, from
//! `scripts/installer_config_paths.py`. Planning keeps explicit sources distinct
//! from managed destinations. Copying never asserts scientific qualification.
mod adopt;
mod cli;
mod defaults;
pub(crate) mod files;
pub(crate) mod gnomad;
mod plan;
mod receipt_json;
mod values;
use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use crate::python_path::{absolute, path_value, repr, tempfile_parent};
pub(crate) use cli::dispatch;
use files::{is_file, metadata, os_error, resolve};
use std::path::{Path, PathBuf};
use values::{Config, Value};

pub(crate) const FOLDERS: &[(&str, &str, &str)] = &[
    ("vep_cache", "VEP_cache", "vep_cache_dir"),
    ("vep_plugin", "VEP_plugin", "vep_plugins_dir"),
    (
        "vep_plugin_cache",
        "VEP_plugin_cache",
        "vep_plugins_cachedir",
    ),
    ("gnomad", "gnomAD", "gnomad_vcf_dir"),
    ("clinvar", "ClinVar", "clinvar_vcf_dir"),
    ("cadd", "CADD", "cadd_parent_dir"),
];
const BUNDLED_PATHS: &[(&str, &str)] = &[
    ("conda_env_yaml", "acmg_conda_rough.yml"),
    (
        "hgnc_table",
        "data/hgnc/hgnc_non_alt_loci_set_sha256_890b18743cc8.tsv",
    ),
    (
        "dispensable_gene_list",
        "data/dispensable_genes/dispensable_gene_list.txt",
    ),
    (
        "gene_mechanism_builder_script",
        "scripts/build_gene_pathogenic_mechanism_cache.py",
    ),
    (
        "gene_nonlof_mechanism_builder_script",
        "scripts/build_gene_nonlof_mechanism_cache.py",
    ),
    (
        "gene_nonlof_mechanism_schema",
        "data/patho_mechanism/gene_nonlof_mechanism_curated_assertions.schema.json",
    ),
    (
        "hpo_condition_mechanism_builder_script",
        "scripts/build_hpo_condition_mechanism_cache.py",
    ),
    (
        "gofcards_clinvar_injector_script",
        "scripts/inject_clinvar_into_gofcards.py",
    ),
    (
        "gofcards_builder_script",
        "scripts/build_gofcards_exact_gof_cache.py",
    ),
    (
        "gofcards_mechanism_review_tsv",
        "data/gofcards/gofcards_mechanism_reviews.tsv",
    ),
    (
        "hpo_assertion_builder_script",
        "scripts/collapse_HPO_anno.py",
    ),
    (
        "mondo_scope_builder_script",
        "scripts/build_mondo_disease_scope.py",
    ),
    (
        "mondo_disease_scope_overrides",
        "data/mondo/disease_scope_overrides.tsv",
    ),
];
const BUNDLED_CACHES: &[(&str, &str)] = &[
    (
        "gofcards_exact_gof_cache",
        "data/gofcards/gofcards_exact_gof.json.gz",
    ),
    (
        "gene_nonlof_mechanism_json",
        "data/patho_mechanism/gene_nonlof_mechanism_curated_assertions.json.gz",
    ),
    (
        "hpo_condition_mechanism_json",
        "data/patho_mechanism/hpo_condition_mechanism_cache.json.gz",
    ),
    (
        "ddg2p_mechanism_evidence",
        "data/patho_mechanism/gene_pathogenic_mechanism_evidence.tsv",
    ),
];
const REPLACEMENT: &[&str] = &[
    "hgnc_table",
    "dispensable_gene_list",
    "gofcards_mechanism_review_tsv",
    "mondo_disease_scope_overrides",
    "gene_nonlof_mechanism_schema",
];
const GENERATED: &[&str] = &[
    "clinvar_vcf",
    "clinvar_aa_stat",
    "clinvar_splice_stat",
    "clinvar_patho_af_stat",
    "clinvar_patho_exon_af_stat",
    "clinvar_pm5_index",
    "clinvar_gene_stat",
    "alphamissense_vcf",
    "alphamissense_vep_vcf",
    "alphamissense_vep_base_vcf",
    "alphamissense_pd_stat",
    "alphamissense_tranx_domain_map",
    "alphamissense_intolerant_domains",
    "alphamissense_intolerant_motifs",
    "all_intolerant_domains",
    "pm1_regions_pkl",
    "interpro_mapping_pickle",
    "clingen_map",
    "gene_dosage_sensitivity",
    "gene_nonlof_mechanism_json",
    "hpo_condition_mechanism_json",
    "ddg2p_mechanism_evidence",
    "gofcards_exact_gof_cache",
    "hpo_assertions",
    "mondo_disease_scope_registry",
    "clinvar_vcv_match_audit",
];
const DIR_EXTRA: &[&str] = &[
    "loftee_parent_dir",
    "loftee_repo",
    "cadd_base_dir",
    "gene_mechanism_cache_dir",
    "gene_mechanism_raw_dir",
    "gofcards_workdir",
    "hpo_raw_dir",
    "pext_tissue_bw_dir",
];
const GENERATED_DIR_EXTRA: &[&str] = &[
    "loftee_parent_dir",
    "gene_mechanism_cache_dir",
    "gene_mechanism_raw_dir",
    "gofcards_workdir",
    "hpo_raw_dir",
];
const AUTO: &[&str] = &[
    "ref_genome",
    "gnomad_vcf_chrX",
    "alphamissense_prescore",
    "utr_annotator_file",
    "human_ancestor_fasta",
    "loftee_conservation_file",
    "gerp_bigwig",
    "conservation_file",
    "splicevault_prescore",
    "cds_fasta_file",
    "mondo_obo",
    "clinvar_vcv_xml",
    "mavedb_file",
    "pext_mean_bw",
    "pext_tissue_bw_dir",
    "loftee_repo",
    "cadd_base_dir",
];
const COMPANIONS: &[&str] = &[
    ".tbi",
    ".csi",
    ".fai",
    ".gzi",
    ".provenance.json",
    ".nm_to_enst.json",
    ".ready.json",
    ".source.json",
    ".annotation-base.json",
];
pub(crate) type Paths = Vec<(String, String)>;
fn error(kind: &'static str, message: impl Into<String>) -> Box<dyn std::error::Error> {
    crate::builder_error::category(
        kind,
        invalid(crate::python_unicode::display(&message.into())),
    )
}
fn require(ok: bool, message: impl Into<String>) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(error("ValueError", message))
    }
}
fn set(values: &mut Paths, key: &str, value: impl Into<String>) {
    let value = value.into();
    if let Some((_, v)) = values.iter_mut().find(|(k, _)| k == key) {
        *v = value;
    } else {
        values.push((key.into(), value));
    }
}
fn get<'a>(values: &'a Paths, key: &str) -> Result<&'a str> {
    values
        .iter()
        .find(|(k, _)| k == key)
        .map(|(_, v)| v.as_str())
        .ok_or_else(|| error("KeyError", repr(key)))
}
fn pathmap(values: &Paths) -> Json {
    Json::Object(
        values
            .iter()
            .map(|(k, v)| (k.clone(), Json::Str(v.clone())))
            .collect(),
    )
}
fn join(path: &Path, name: impl AsRef<Path>) -> PathBuf {
    path_value(&path.join(name).display().to_string())
}
fn name(path: &Path) -> String {
    match path.components().next_back() {
        Some(std::path::Component::Normal(n)) => n.to_string_lossy().into_owned(),
        Some(std::path::Component::ParentDir) => "..".into(),
        _ => String::new(),
    }
}
fn directory_key(key: &str) -> bool {
    FOLDERS.iter().any(|(_, _, k)| *k == key) || DIR_EXTRA.contains(&key)
}
fn generated_directory(key: &str) -> bool {
    FOLDERS.iter().any(|(_, _, k)| *k == key) || GENERATED_DIR_EXTRA.contains(&key)
}
fn bundled_path(key: &str) -> &'static str {
    BUNDLED_PATHS
        .iter()
        .chain(BUNDLED_CACHES)
        .find(|(k, _)| *k == key)
        .unwrap()
        .1
}
fn folder_name(s: &str) -> bool {
    s.bytes()
        .next()
        .is_some_and(|c| c.is_ascii_alphanumeric() || c == b'_')
        && s.bytes()
            .all(|c| c.is_ascii_alphanumeric() || b"_.-".contains(&c))
}
fn forced(key: &str) -> bool {
    let flags: &[&str] = match key {
        "gofcards_exact_gof_cache" => &[
            "PRIVA_FORCE_GOFCARDS_CACHE",
            "PRIVA_FORCE_GOFCARDS_CLINVAR_INJECTION",
        ],
        "gene_nonlof_mechanism_json" => &["PRIVA_FORCE_NONLOF_MECHANISM_CACHE"],
        "hpo_condition_mechanism_json" => &["PRIVA_FORCE_HPO_MECHANISM_JSON"],
        "ddg2p_mechanism_evidence" => &["PRIVA_FORCE_GENE_MECHANISM_CACHE"],
        _ => &[],
    };
    std::iter::once("PRIVA_FORCE_ALL_CACHES")
        .chain(flags.iter().copied())
        .any(|f| std::env::var(f).ok().as_deref() == Some("1"))
}
pub(crate) struct Context {
    pub(crate) helper_root: PathBuf,
    pub(crate) seed: u32,
}
impl Context {
    pub(crate) fn new(helper_root: &Path) -> Result<Self> {
        let seed=std::env::var("PYTHONHASHSEED").ok().filter(|v|!v.is_empty()&&v.bytes().all(|c|c.is_ascii_digit())).and_then(|v|v.parse().ok()).ok_or_else(||error("RuntimeError","Unsupported profile: set an explicit PYTHONHASHSEED in 0..4294967295; another Python process's random secret is unknown"))?;
        Ok(Self {
            helper_root: helper_root.to_path_buf(),
            seed,
        })
    }
    fn source_root(&self) -> Result<PathBuf> {
        Ok(
            resolve(&self.helper_root.join("scripts/installer_config_paths.py"))?
                .parent()
                .and_then(Path::parent)
                .ok_or_else(|| error("IndexError", "1"))?
                .to_path_buf(),
        )
    }
    fn import(&self, name: &str, expected: &[u8]) -> Result<()> {
        let path = self.helper_root.join("scripts").join(name);
        let actual = std::fs::read(files::native(&path)?).map_err(|e| os_error(e, &path, None))?;
        require(
            actual == expected,
            format!(
                "Helper source asset differs from the native implementation: {}",
                path.display()
            ),
        )
    }
    fn splice(&self) -> Result<()> {
        self.import(
            "spliceai_cache_input.py",
            include_bytes!("../../../scripts/spliceai_cache_input.py"),
        )
    }
}

/// `render_config(text, updates)` for another command's ordinary Unicode keys
/// and values, rendered in the same order and with the same quoting.
pub(crate) fn render_text(original: &str, updates: &[(String, String)]) -> Result<String> {
    use crate::python_unicode::encode;
    let updates: Paths = updates
        .iter()
        .map(|(key, value)| (encode(key), encode(value)))
        .collect();
    adopt::render_config(original, &updates)
}

/// The configuration replacement of `installer_prerequisites.manual_inputs`: a
/// written, flushed and synced `.NAME.` temporary beside the file, the file's
/// permission bits, a refusal when another writer changed the file, then one
/// rename. The temporary never survives.
pub(crate) fn replace_text_synced(
    destination: &Path,
    original: &str,
    updated: &str,
    message: &str,
) -> Result<()> {
    use std::io::Write;
    use std::os::unix::fs::{MetadataExt, PermissionsExt};
    let parent = destination
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let (file, temporary) = files::temporary(parent, &format!(".{}.", name(destination)))?;
    let result = (|| {
        let mut file = file;
        file.write_all(updated.as_bytes())
            .map_err(|e| os_error(e, &temporary, None))?;
        file.flush().map_err(|e| os_error(e, &temporary, None))?;
        file.sync_all().map_err(|e| os_error(e, &temporary, None))?;
        drop(file);
        std::fs::set_permissions(
            files::native(&temporary)?,
            std::fs::Permissions::from_mode(metadata(destination)?.mode() & 0o777),
        )
        .map_err(|e| os_error(e, &temporary, None))?;
        require(files::read_text(destination)? == original, message)?;
        files::rename(&temporary, destination)
    })();
    if files::exists(&temporary)? {
        files::unlink_missing(&temporary)?;
    }
    result
}
