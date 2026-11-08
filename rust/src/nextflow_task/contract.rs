//! Configuration, resource inventory and identity digests of the task contract.
use super::{Result, refuse};
use crate::json_text::{self, DumpsOptions, Json};
use crate::python_path::path_value;
use crate::yaml_config::{self, ConfigDocument, ConfigValue, Parsed, Scalar};
use std::fs;
use std::io::Read;
use std::os::unix::fs::{MetadataExt, PermissionsExt};
use std::path::{Path, PathBuf};

pub(crate) const PRIORITY_KEYS: &[&str] = &[
    "clinvar_aa_stat",
    "clinvar_splice_stat",
    "clinvar_pm5_index",
    "all_intolerant_domains",
    "alphamissense_intolerant_domains",
    "pm1_regions_pkl",
    "alphamissense_vcf",
    "alphamissense_tranx_domain_map",
    "clinvar_patho_af_stat",
    "clinvar_patho_exon_af_stat",
    "interpro_mapping_pickle",
    "clingen_map",
    "gene_dosage_sensitivity",
    "hpo_condition_mechanism_json",
    "gene_nonlof_mechanism_json",
];
const ANNOTATION_KEYS: &[&str] = &[
    "ref_genome",
    "gnomad_vcf_chrX",
    "clinvar_vcf",
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
    "cadd_base_dir",
];
pub(crate) const OPTIONAL_INPUTS: &[&str] = &[
    "control_vcf",
    "alt_disease_vcf",
    "pp1_vcf",
    "pp1_ped",
    "relevant_gene_list",
    "dispensable_gene_list",
    "cds_fasta_file",
];
/// Executables whose bytes are pinned by SHA-256 even when large caches use
/// stat reuse: `priva_binary` runs the converted pipeline steps from the
/// shell scripts; `priva_ranking_binary` is the optional ranking program.
pub(crate) const EXECUTABLE_KEYS: &[&str] = &["priva_binary", "priva_ranking_binary"];

/// The legacy assembly names and the one each stands for.
pub(crate) fn assembly(name: &str) -> Option<&'static str> {
    match name {
        "hg19" | "GRCh37" => Some("hg19"),
        "hg38" | "GRCh38" => Some("hg38"),
        _ => None,
    }
}

/// `^[A-Za-z0-9_./:+-]+$`.
fn safe_path(value: &str) -> bool {
    !value.is_empty()
        && value
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"_./:+-".contains(&b))
}

/// `^[A-Za-z0-9_][A-Za-z0-9_.-]*$`.
pub(crate) fn safe_id(value: &str) -> bool {
    let mut bytes = value.bytes();
    match bytes.next() {
        Some(first) if first.is_ascii_alphanumeric() || first == b'_' => {}
        _ => return false,
    }
    bytes.all(|b| b.is_ascii_alphanumeric() || b"_.-".contains(&b))
}

/// A legacy configuration as `yaml.load` with duplicate-key refusal reads it.
#[derive(Debug)]
pub(crate) struct Config {
    pub(crate) text: String,
    pub(crate) doc: ConfigDocument,
    /// The SpliceAI version of the configured score files, `v1` or
    /// `SpliceAI2`, when at least one of them is present.
    pub(crate) spliceai: Option<&'static str>,
}

impl Config {
    pub(crate) fn get(&self, key: &str) -> Option<&ConfigValue> {
        self.doc.get(key)
    }

    /// The value when it is a string, as `isinstance(value, str)`.
    pub(crate) fn str(&self, key: &str) -> Option<&str> {
        match self.doc.scalar(key) {
            Some(Scalar::Str(text)) => Some(text),
            _ => None,
        }
    }

    /// Python truthiness of `config.get(key)`.
    pub(crate) fn truthy(&self, key: &str) -> bool {
        match self.get(key) {
            None => false,
            Some(ConfigValue::Scalar(scalar)) => match scalar {
                Scalar::Null => false,
                Scalar::Bool(flag) => *flag,
                Scalar::Int(text) => text.trim_start_matches('-') != "0",
                Scalar::Str(text) => !text.is_empty(),
                Scalar::Other(text) => text
                    .replace('_', "")
                    .parse::<f64>()
                    .map_or(true, |value| value != 0.0),
            },
            Some(ConfigValue::Opaque(text)) => !matches!(text.trim(), "[]" | "{}" | ""),
        }
    }

    /// `str(config.get(key, default))` for scalar values.
    pub(crate) fn text_of(&self, key: &str, default: &str) -> String {
        match self.get(key) {
            None => default.to_owned(),
            Some(ConfigValue::Scalar(Scalar::Str(text))) => text.clone(),
            Some(ConfigValue::Scalar(scalar)) => scalar.python_repr(),
            Some(ConfigValue::Opaque(text)) => text.clone(),
        }
    }

    pub(crate) fn assembly(&self) -> &'static str {
        self.str("assembly").and_then(assembly).unwrap_or("hg19")
    }
}

pub(crate) fn read_text(path: &Path) -> Result<String> {
    let bytes = fs::read(path).map_err(|e| crate::python_path::os_error(e, path, None))?;
    String::from_utf8(bytes)
        .map_err(|_| format!("'utf-8' codec can't decode {}", path.display()).into())
}

pub(crate) fn read_config(path: &Path) -> Result<Config> {
    let text = read_text(path)?;
    let doc = match yaml_config::parse_config(&text) {
        Ok(Parsed::Mapping(doc)) => doc,
        Ok(_) => return refuse("Legacy config must be a YAML mapping"),
        Err(error) => return Err(error.to_string().into()),
    };
    for (index, entry) in doc.entries.iter().enumerate() {
        if doc.entries[..index].iter().any(|e| e.key == entry.key) {
            return refuse(format!("Duplicate YAML key: {}", entry.key));
        }
    }
    let mut config = Config {
        text,
        doc,
        spliceai: None,
    };
    if config.str("assembly").and_then(assembly).is_none() {
        return refuse("Config assembly must be hg19/hg38/GRCh37/GRCh38");
    }
    // Annotated imports skip the annotation shell, so inspect the actual
    // score files here too, including cached launches and tasks. SpliceAI v1
    // and SpliceAI2 each have their own cutoffs (scripts/spliceai_cutoffs.py);
    // the SNV and indel files must be the same version.
    let mut found: Vec<(&str, &'static str)> = Vec::new();
    for key in ["spliceai_snv_prescore", "spliceai_indel_prescore"] {
        let Some(value) = config.str(key) else {
            continue;
        };
        if value.is_empty() || !Path::new(value).is_file() {
            continue; // Keep existing missing-resource diagnostics at their stage.
        }
        use crate::spliceai_cache_input::{Failure, inspect_header};
        let detected = match inspect_header(Path::new(value)) {
            Ok(header) => header,
            Err(Failure::Rejected(m) | Failure::Crash(m)) => {
                return refuse(format!("Cannot verify {key} header: {m}"));
            }
            Err(Failure::Helper(code)) => {
                return refuse(format!("Cannot verify {key} header: exit {code}"));
            }
        };
        found.push((
            key,
            if detected.format == "vcf_v1" {
                "v1"
            } else {
                "SpliceAI2"
            },
        ));
    }
    if let [(_, snv), (_, indel)] = found.as_slice() {
        if snv != indel {
            return refuse(format!(
                "SpliceAI SNV ({snv}) and indel ({indel}) score files must be the same version; each version uses its own cutoffs"
            ));
        }
    }
    config.spliceai = found.first().map(|(_, version)| *version);
    Ok(config)
}

/// `Path.exists()`: false for the ordinary missing-path errors.
pub(crate) fn exists(path: &Path) -> bool {
    fs::metadata(path).is_ok()
}

/// A configured path: safe characters, absolute, present, not an LFS pointer.
pub(crate) fn checked_path(value: Option<&str>, key: &str) -> Result<PathBuf> {
    let Some(value) = value.filter(|v| safe_path(v)) else {
        return refuse(format!(
            "{key}: required path is empty or contains whitespace/shell metacharacters"
        ));
    };
    let path = path_value(value);
    if !path.is_absolute() {
        return refuse(format!("{key}: expected an absolute path, got {value}"));
    }
    if !exists(&path) {
        return refuse(format!("{key}: missing path {value}"));
    }
    if path.is_file() {
        let mut head = Vec::with_capacity(200);
        fs::File::open(&path)?.take(200).read_to_end(&mut head)?;
        if head.starts_with(b"version https://git-lfs.github.com/spec/v1") {
            return refuse(format!(
                "{key}: Git LFS pointer must be materialized: {value}"
            ));
        }
    }
    Ok(path)
}

fn config_path(config: &Config, key: &str) -> Result<PathBuf> {
    checked_path(config.str(key), key)
}

fn path_text(path: &Path) -> String {
    path.display().to_string()
}

/// `.fa.gz` → `.cds_cache.pkl.gz`, else the suffix is appended.
pub(crate) fn cds_cache_path(fasta: &str) -> PathBuf {
    match fasta.strip_suffix(".fa.gz") {
        Some(stem) => path_value(&format!("{stem}.cds_cache.pkl.gz")),
        None => path_value(&format!("{fasta}.cds_cache.pkl.gz")),
    }
}

/// Every resource the shell consumers read, by key, sorted like a dict
/// written with `sort_keys` (insertion order does not matter downstream).
pub(crate) fn runtime_resources(
    config: &Config,
    root: &Path,
    mode: &str,
    vcf_hub_mode: &str,
    cadd_hub_mode: &str,
) -> Result<Vec<(String, PathBuf)>> {
    let mut paths: Vec<(String, PathBuf)> = Vec::new();
    let put = |paths: &mut Vec<(String, PathBuf)>, key: String, path: PathBuf| {
        if let Some(slot) = paths.iter_mut().find(|(k, _)| *k == key) {
            slot.1 = path;
        } else {
            paths.push((key, path));
        }
    };
    let annotation: &[&str] = if mode == "raw" { ANNOTATION_KEYS } else { &[] };
    for key in PRIORITY_KEYS.iter().chain(annotation) {
        put(&mut paths, key.to_string(), config_path(config, key)?);
    }
    for key in OPTIONAL_INPUTS {
        if config.truthy(key) {
            put(&mut paths, key.to_string(), config_path(config, key)?);
        }
    }
    for key in EXECUTABLE_KEYS {
        // priva_binary is required: every phase after annotation runs it.
        if *key == "priva_binary" || config.truthy(key) {
            let binary = config_path(config, key)?;
            if !executable_file(&binary) {
                return refuse(format!("{key} must be a nonempty executable file"));
            }
            put(&mut paths, key.to_string(), binary);
        }
    }
    if config.truthy("pp1_vcf") != config.truthy("pp1_ped") {
        return refuse("pp1_vcf and pp1_ped must be supplied together");
    }
    let root_text = path_text(root);
    let under = |relative: &str| path_value(&format!("{root_text}/{relative}"));
    let hpo = match config.get("hpo_assertions") {
        None | Some(ConfigValue::Scalar(Scalar::Null)) => {
            under("data/hpo/genes_to_phenotype.assertions.tsv.gz")
        }
        Some(ConfigValue::Scalar(Scalar::Str(text))) if text.is_empty() || text == "null" => {
            under("data/hpo/genes_to_phenotype.assertions.tsv.gz")
        }
        _ => PathBuf::from(config.text_of("hpo_assertions", "")),
    };
    let selected_gene_stat = if config.truthy("clinvar_gene_stat") {
        config.text_of("clinvar_gene_stat", "")
    } else {
        String::new()
    };
    let gene_stat = if selected_gene_stat.contains('/') {
        PathBuf::from(crate::python_path::absolute(&path_value(
            &selected_gene_stat,
        ))?)
    } else {
        under(&format!("data/ClinVar/{selected_gene_stat}"))
    };
    let repeat_name = file_name(&config.text_of("repeat_region_file_name", ""));
    let mut fixed: Vec<(String, PathBuf)> = vec![
        ("hpo_assertions_used_by_shell".into(), hpo),
        (
            "alphamissense_mean_used_by_shell".into(),
            under("data/alphamissense/alphamissense_mean_score.tsv"),
        ),
        ("clinvar_gene_stat_used_by_shell".into(), gene_stat),
        (
            "repeat_regions_used_by_shell".into(),
            under(&format!("data/repeats/{repeat_name}")),
        ),
        (
            "hgnc_used_by_scientific_code".into(),
            under("data/hgnc/non_alt_loci_set.tsv"),
        ),
    ];
    if !config.truthy("dispensable_gene_list") {
        fixed.push((
            "dispensable_genes_default".into(),
            under("data/dispensable_genes/dispensable_gene_list.txt"),
        ));
    }
    if !config.truthy("clinvar_gene_stat") || !config.truthy("repeat_region_file_name") {
        return refuse("clinvar_gene_stat and repeat_region_file_name are required");
    }
    if mode == "raw" {
        fixed.push((
            "primary_contigs".into(),
            under("data/liftover/ucsc_GRC.primary.contigs.tsv"),
        ));
        fixed.push((
            "contig_name_mapping".into(),
            under("data/liftover/ucsc_to_GRC.contig.map.tsv"),
        ));
        if !config.truthy("loeuf_prescore") {
            return refuse("loeuf_prescore is required");
        }
        let Some(loeuf) = config.str("loeuf_prescore") else {
            return Err("expected str, bytes or os.PathLike object for loeuf_prescore".into());
        };
        fixed.push((
            "loeuf_used_by_shell".into(),
            under(&format!("data/loeuf/{}", file_name(loeuf))),
        ));
        // gnomAD is read chromosome by chromosome from the chrX template.
        let template = paths
            .iter()
            .find(|(k, _)| k == "gnomad_vcf_chrX")
            .map(|(_, p)| path_text(p))
            .unwrap_or_default();
        if !template.contains("chrX") {
            return refuse("gnomad_vcf_chrX must contain chrX for the legacy template");
        }
        let chromosomes = (1..=22).map(|n| n.to_string()).chain(["Y".to_owned()]);
        for chromosome in chromosomes {
            let key = format!("gnomad_chr{chromosome}");
            let value = template.replace("chrX", &format!("chr{chromosome}"));
            let path = checked_path(Some(&value), &key)?;
            put(&mut paths, key, path);
        }
        let reference = paths
            .iter()
            .find(|(k, _)| k == "ref_genome")
            .map(|(_, p)| path_text(p))
            .unwrap_or_default();
        fixed.push((
            "reference_fai".into(),
            path_value(&format!("{reference}.fai")),
        ));
        if config.assembly() == "hg38" {
            for key in ["gerp_bigwig", "mavedb_file"] {
                put(&mut paths, key.into(), config_path(config, key)?);
            }
        }
        if matches!(
            config.get("pext_enabled"),
            Some(ConfigValue::Scalar(Scalar::Bool(true)))
        ) {
            put(
                &mut paths,
                "pext_mean_bw".into(),
                config_path(config, "pext_mean_bw")?,
            );
            if config.truthy("pext_tissues") {
                put(
                    &mut paths,
                    "pext_tissue_bw_dir".into(),
                    config_path(config, "pext_tissue_bw_dir")?,
                );
            }
        }
        for (key, mode_value) in [
            ("hub_vcf_file", vcf_hub_mode),
            ("hub_cadd_file", cadd_hub_mode),
        ] {
            if mode_value != "off" && config.truthy(key) {
                let path = hub_source(config, key)?;
                let checked = checked_path(Some(&path_text(&path)), key)?;
                put(&mut paths, key.into(), checked);
            }
        }
    }
    if config.assembly() == "hg38" {
        fixed.push((
            "mavedb_assay_metadata_used_by_shell".into(),
            under("data/MaveDB/assay_mapping_result.neutral_bands.tsv"),
        ));
    }
    for (key, path) in fixed {
        let checked = checked_path(Some(&path_text(&path)), &key)?;
        put(&mut paths, key, checked);
    }
    if let Some((_, fasta)) = paths.iter().find(|(k, _)| k == "cds_fasta_file") {
        // Copy this mutable lazy cache into each prioritization task; don't
        // write a newly parsed pickle beside the shared FASTA.
        let cache = cds_cache_path(&path_text(fasta));
        if exists(&cache) {
            let checked = checked_path(Some(&path_text(&cache)), "cds_parsed_cache")?;
            put(&mut paths, "cds_parsed_cache".into(), checked);
        }
    }
    let indexed: Vec<(String, PathBuf)> = paths
        .iter()
        .filter(|(_, p)| {
            let text = path_text(p);
            text.ends_with(".vcf.gz") || text.ends_with(".vcf.bgz")
        })
        .cloned()
        .collect();
    for (key, path) in indexed {
        // Legacy tools rebuild missing tabix indexes; require a ready index.
        let index_key = format!("{key}_index");
        let checked = checked_path(Some(&format!("{}.tbi", path_text(&path))), &index_key)?;
        put(&mut paths, index_key, checked);
    }
    Ok(paths)
}

/// A hub file: absolute, or relative to `hub_file_dir`.
pub(crate) fn hub_source(config: &Config, key: &str) -> Result<PathBuf> {
    let value = config
        .str(key)
        .ok_or_else(|| format!("{key}: expected a path string"))?;
    let path = path_value(value);
    if path.is_absolute() {
        return Ok(path);
    }
    let directory = match config.get("hub_file_dir") {
        None => String::new(),
        Some(ConfigValue::Scalar(Scalar::Str(text))) => text.clone(),
        Some(_) => return Err("hub_file_dir: expected a path string".into()),
    };
    Ok(path_value(&format!("{directory}/{value}")))
}

/// `PurePath(text).name`.
pub(crate) fn file_name(text: &str) -> String {
    let path = path_value(text);
    path.file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default()
}

pub(crate) fn executable_file(path: &Path) -> bool {
    match fs::metadata(path) {
        Ok(meta) => meta.is_file() && meta.len() > 0 && access_executable(path),
        Err(_) => false,
    }
}

pub(crate) fn access_executable(path: &Path) -> bool {
    use std::ffi::CString;
    use std::os::unix::ffi::OsStrExt;
    let Ok(text) = CString::new(path.as_os_str().as_bytes()) else {
        return false;
    };
    // SAFETY: a valid NUL-terminated path; access() only reads it.
    unsafe { libc::access(text.as_ptr(), libc::X_OK) == 0 }
}

fn mtime_ns(meta: &fs::Metadata) -> i128 {
    meta.mtime() as i128 * 1_000_000_000 + meta.mtime_nsec() as i128
}

fn ctime_ns(meta: &fs::Metadata) -> i128 {
    meta.ctime() as i128 * 1_000_000_000 + meta.ctime_nsec() as i128
}

pub(crate) fn object(items: Vec<(&str, Json)>) -> Json {
    Json::Object(items.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}

pub(crate) fn text(value: impl Into<String>) -> Json {
    Json::Str(value.into())
}

pub(crate) fn sha256(path: &Path) -> Result<String> {
    crate::sha256::sha256_file(path).map_err(|e| crate::python_path::os_error(e, path, None))
}

/// `path`, `resolved_path`, `type`, `size`, `mtime_ns`, optional `sha256`.
pub(crate) fn inventory_record(path: &Path, checksum: bool) -> Result<Json> {
    let meta = fs::metadata(path).map_err(|e| crate::python_path::os_error(e, path, None))?;
    let resolved = fs::canonicalize(path)?;
    let mut items = vec![
        ("path", text(path_text(path))),
        ("resolved_path", text(path_text(&resolved))),
        (
            "type",
            text(if meta.is_dir() { "directory" } else { "file" }),
        ),
        ("size", Json::Int(meta.len() as i128)),
        ("mtime_ns", Json::Int(mtime_ns(&meta))),
    ];
    if checksum && meta.is_file() {
        items.push(("sha256", text(sha256(path)?)));
    }
    Ok(object(items))
}

/// Pin executable bytes even when large annotation caches use stat reuse.
fn executable_record(path: &Path, key: &str) -> Result<Json> {
    let identity = |meta: &fs::Metadata| {
        (
            meta.dev(),
            meta.ino(),
            meta.mode(),
            meta.size(),
            mtime_ns(meta),
            ctime_ns(meta),
            meta.nlink(),
        )
    };
    let capture = || -> Result<_> {
        Ok((
            identity(&fs::symlink_metadata(path)?),
            identity(&fs::metadata(path)?),
            fs::canonicalize(path)?,
        ))
    };
    let before = capture()?;
    if !executable_file(path) {
        return refuse(format!("{key} must be a nonempty executable file"));
    }
    let mut record = inventory_record(path, true)?;
    let after = capture()?;
    if before != after {
        return refuse(format!("{key} changed during identity capture"));
    }
    record.set("executable_mode", Json::Int((before.1.2 & 0o777) as i128));
    Ok(record)
}

/// Every entry below `directory`, without following directory symlinks,
/// ordered by path components as `sorted(path.rglob("*"))`.
fn tree(directory: &Path) -> Result<Vec<PathBuf>> {
    let mut found = Vec::new();
    let mut pending = vec![directory.to_path_buf()];
    while let Some(current) = pending.pop() {
        for entry in fs::read_dir(&current)? {
            let entry = entry?;
            let path = entry.path();
            if entry.file_type()?.is_dir() {
                pending.push(path.clone());
            }
            found.push(path);
        }
    }
    found.sort_by(|a, b| a.components().cmp(b.components()));
    Ok(found)
}

pub(crate) fn inventory(paths: &[(String, PathBuf)], checksum: bool) -> Result<Json> {
    let mut sorted: Vec<&(String, PathBuf)> = paths.iter().collect();
    sorted.sort_by(|a, b| a.0.cmp(&b.0));
    let mut records = Json::Object(Vec::new());
    for (key, path) in sorted {
        let record = if EXECUTABLE_KEYS.contains(&key.as_str()) {
            executable_record(path, key)?
        } else {
            inventory_record(path, checksum)?
        };
        records.set(key, record);
        if path.is_dir() {
            // Include every dependency below directory resources without
            // hashing TB of bytes by default. Symlinks are recorded and their
            // target stat included; directory symlinks are not traversed.
            for child in tree(path)? {
                let link = fs::symlink_metadata(&child)?.file_type().is_symlink();
                if link && child.is_dir() {
                    return refuse(format!(
                        "Directory symlink in resource tree needs explicit resolved config: {}",
                        child.display()
                    ));
                }
                if child.is_file() {
                    let relative = child.strip_prefix(path).unwrap_or(&child);
                    records.set(
                        &format!("{key}/{}", relative.display()),
                        inventory_record(&child, checksum)?,
                    );
                }
            }
        }
    }
    Ok(records)
}

/// `hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")))`.
pub(crate) fn json_digest(value: &Json) -> Result<String> {
    let encoded = json_text::dumps_with_options(
        value,
        DumpsOptions {
            sort_keys: true,
            compact: true,
            ..DumpsOptions::default()
        },
    )?;
    Ok(crate::sha256::sha256_hex(encoded.as_bytes()))
}

/// `json.dumps(value, indent=2, sort_keys=True) + "\n"`.
pub(crate) fn json_document(value: &Json) -> Result<String> {
    let mut encoded = json_text::dumps_with_options(
        value,
        DumpsOptions {
            indent: Some(2),
            sort_keys: true,
            ..DumpsOptions::default()
        },
    )?;
    encoded.push('\n');
    Ok(encoded)
}

pub(crate) fn write_json(path: &Path, value: &Json) -> Result<()> {
    fs::write(path, json_document(value)?).map_err(|e| crate::python_path::os_error(e, path, None))
}

pub(crate) fn read_json(path: &Path) -> Result<Json> {
    Ok(json_text::parse(&read_text(path)?, false)?)
}

/// The content digest of the source tree, as `main.nf` computes it.
pub(crate) fn source_digest(root: &Path) -> Result<String> {
    let mut files: Vec<PathBuf> = Vec::new();
    for folder in ["scripts", "bin", "modules", "conf"] {
        let directory = root.join(folder);
        if !directory.is_dir() {
            continue;
        }
        for path in tree(&directory)? {
            let text = path_text(&path);
            if path.is_file() && !text.ends_with(".pyc") && !text.contains("__pycache__") {
                files.push(path);
            }
        }
    }
    files.push(root.join("main.nf"));
    files.push(root.join("nextflow.config"));
    // Match Nextflow's full path-string ordering, which puts
    // conf/cadd.config before conf/cadd/profile.json.
    files.sort_by_key(|path| path_text(path));
    let mut digest = crate::sha256::Sha256::new();
    for path in files {
        let relative = path.strip_prefix(root).unwrap_or(&path);
        digest.update(path_text(relative).as_bytes());
        digest.update(b"\0");
        digest.update(&fs::read(&path).map_err(|e| crate::python_path::os_error(e, &path, None))?);
    }
    Ok(digest.finish_hex())
}

/// The keys `verify_manifest` compares with the manifest's own record.
pub(crate) fn verify_manifest(
    manifest: &Json,
    config_path: &Path,
    root: &Path,
    mode: &str,
    vcf_hub_mode: &str,
    cadd_hub_mode: &str,
) -> Result<Vec<(String, PathBuf)>> {
    let config = read_config(config_path)?;
    let Json::Object(items) = manifest else {
        return Err("Resource manifest must be a JSON object".into());
    };
    let empty = |value: Option<&Json>| match value {
        None | Some(Json::Null) => true,
        Some(Json::Object(items)) => items.is_empty(),
        Some(Json::Array(items)) => items.is_empty(),
        Some(Json::Str(text)) => text.is_empty(),
        Some(Json::Bool(flag)) => !flag,
        Some(Json::Int(value)) => *value == 0,
        Some(Json::Float(value)) => *value == 0.0,
    };
    let schema_one = manifest
        .get("schema_version")
        .is_some_and(|v| json_text::py_equal(v, &Json::Int(1)));
    if !schema_one || empty(manifest.get("resources")) {
        return refuse("Resource manifest is incomplete; create it with priva task manifest");
    }
    let content = Json::Object(
        items
            .iter()
            .filter(|(k, _)| k != "bundle_id")
            .cloned()
            .collect(),
    );
    let recorded = manifest.get("bundle_id").cloned().unwrap_or(Json::Null);
    if !json_text::py_equal(&recorded, &text(json_digest(&content)?)) {
        return refuse("Resource manifest bundle_id does not match its content");
    }
    let expected = [
        ("config_sha256", sha256(config_path)?),
        ("legacy_root", path_text(root)),
        ("mode", mode.to_owned()),
        ("assembly", config.assembly().to_owned()),
        ("vcf_hub_mode", vcf_hub_mode.to_owned()),
        ("cadd_hub_mode", cadd_hub_mode.to_owned()),
    ];
    for (key, value) in expected {
        let matches = manifest
            .get(key)
            .is_some_and(|recorded| json_text::py_equal(recorded, &text(value)));
        if !matches {
            return refuse(format!("Resource manifest mismatch: {key}"));
        }
    }
    let paths = runtime_resources(&config, root, mode, vcf_hub_mode, cadd_hub_mode)?;
    let current = inventory(&paths, false)?;
    let Some(Json::Object(recorded)) = manifest.get("resources") else {
        return refuse("Resource inventory membership changed; freeze and rebuild manifest");
    };
    let Json::Object(current_items) = &current else {
        unreachable!()
    };
    let mut current_keys: Vec<&str> = current_items.iter().map(|(k, _)| k.as_str()).collect();
    let mut recorded_keys: Vec<&str> = recorded.iter().map(|(k, _)| k.as_str()).collect();
    current_keys.sort_unstable();
    current_keys.dedup();
    recorded_keys.sort_unstable();
    recorded_keys.dedup();
    if current_keys != recorded_keys {
        return refuse("Resource inventory membership changed; freeze and rebuild manifest");
    }
    for (key, record) in current_items {
        let mut expected = manifest
            .get("resources")
            .and_then(|r| r.get(key))
            .cloned()
            .unwrap_or(Json::Null);
        if !EXECUTABLE_KEYS.contains(&key.as_str()) {
            if let Json::Object(fields) = &mut expected {
                fields.retain(|(name, _)| name != "sha256");
            }
        }
        if !json_text::py_equal(record, &expected) {
            return refuse(format!("Resource changed since manifest creation: {key}"));
        }
    }
    Ok(paths)
}

/// `shutil.copy2` then `chmod u+w`; returns the absolute destination.
pub(crate) fn copy_writable(source: &Path, destination: &Path) -> Result<PathBuf> {
    if let Some(parent) = destination.parent() {
        fs::create_dir_all(parent)?;
    }
    copy2(source, destination)?;
    let mut permissions = fs::metadata(destination)?.permissions();
    permissions.set_mode(permissions.mode() | 0o200);
    fs::set_permissions(destination, permissions)?;
    Ok(PathBuf::from(crate::python_path::absolute(destination)?))
}

/// `shutil.copy2`: bytes, permission bits and access/modification times.
pub(crate) fn copy2(source: &Path, destination: &Path) -> Result<()> {
    fs::copy(source, destination)
        .map_err(|e| crate::python_path::os_error(e, source, Some(destination)))?;
    let meta = fs::metadata(source)?;
    let times = fs::FileTimes::new()
        .set_accessed(meta.accessed()?)
        .set_modified(meta.modified()?);
    fs::OpenOptions::new()
        .write(true)
        .open(destination)
        .and_then(|file| file.set_times(times))
        .or_else(|_| {
            // A read-only copy cannot be opened for writing; use utimensat.
            set_times_by_path(destination, &meta)
        })?;
    Ok(())
}

fn set_times_by_path(path: &Path, meta: &fs::Metadata) -> std::io::Result<()> {
    use std::ffi::CString;
    use std::os::unix::ffi::OsStrExt;
    let text = CString::new(path.as_os_str().as_bytes())?;
    let times = [
        libc::timespec {
            tv_sec: meta.atime() as libc::time_t,
            tv_nsec: meta.atime_nsec() as _,
        },
        libc::timespec {
            tv_sec: meta.mtime() as libc::time_t,
            tv_nsec: meta.mtime_nsec() as _,
        },
    ];
    // SAFETY: valid path and a two-element timespec array.
    let status = unsafe { libc::utimensat(libc::AT_FDCWD, text.as_ptr(), times.as_ptr(), 0) };
    if status == 0 {
        Ok(())
    } else {
        Err(std::io::Error::last_os_error())
    }
}

/// Re-spell every plain top-level scalar the way `yaml.safe_dump` would
/// (so a YAML 1.2 reader such as `yq` sees the values PyYAML loaded), apply
/// `updates` in insertion order, and keep comments and other values as text.
pub(crate) fn render_config(config: &Config, updates: &[(String, Json)]) -> Result<String> {
    let mut all = Json::Object(Vec::new());
    for entry in &config.doc.entries {
        let ConfigValue::Scalar(scalar) = &entry.value else {
            continue;
        };
        let value = match scalar {
            Scalar::Null => Json::Null,
            Scalar::Bool(flag) => Json::Bool(*flag),
            Scalar::Int(digits) => match digits.parse::<i128>() {
                Ok(number) => Json::Int(number),
                Err(_) => continue,
            },
            Scalar::Str(text) => Json::Str(text.clone()),
            Scalar::Other(_) => continue,
        };
        all.set(&entry.key, value);
    }
    for (key, value) in updates {
        if matches!(config.get(key), Some(ConfigValue::Opaque(_))) {
            return Err(format!("{key}: cannot replace a nested configuration value").into());
        }
        all.set(key, value.clone());
    }
    let Json::Object(all) = all else {
        unreachable!()
    };
    crate::config_text::render(&config.text, &all, crate::config_text::Writer::Reference)
        .map_err(|e| e.message.into())
}
