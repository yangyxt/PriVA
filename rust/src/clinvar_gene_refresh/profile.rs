//! The resource policy, the scientific profile and its signature, and the
//! unqualified installation-candidate receipts.
use super::*;
use crate::installer_prerequisites::paths::{exists, is_dir, is_file, is_relative_to};
use crate::python_path::absolute;
use crate::spliceai_cache_input::{self as spliceai, Failure};
use std::process::{Command, Stdio};
use std::sync::OnceLock;

const NCBI: &str = "https://ftp.ncbi.nlm.nih.gov";
const PROFILE_KEYS: [&str; 25] = [
    "assembly",
    "vep_installed_version",
    "ref_genome",
    "vep_cache_dir",
    "vep_plugins_dir",
    "vep_plugins_cachedir",
    "utr_annotator_file",
    "loeuf_prescore",
    "alphamissense_prescore",
    "spliceai_snv_prescore",
    "spliceai_indel_prescore",
    "primateai_prescore",
    "splicevault_prescore",
    "conservation_file",
    "loftee_repo",
    "human_ancestor_fasta",
    "loftee_conservation_file",
    "gerp_bigwig",
    "mavedb_file",
    "gnomad_vcf_chrX",
    "pext_enabled",
    "pext_mean_bw",
    "pext_tissue_bw_dir",
    "pext_tissues",
    "priva_python_hash_seed",
];
const SCIENTIFIC_SOURCES: [&str; 7] = [
    "clinvar_stat_variants.py",
    "combine_annotations.py",
    "annotation_vcf.sh",
    "clinvar_annotation_contract.py",
    "common_bash_utils.sh",
    "stat_aachange_clinvar.py",
    "hpo_penetrance.py",
];
// v2 binds the installer by the ClinVar stage bodies below instead of the whole
// install_utils.sh, so v1 and v2 profiles can never share a signature.
const PROFILE_SCHEMA: &str = "priva_clinvar_scientific_profile_v2";
const INSTALLER: &str = "install_utils.sh";
// The install_utils.sh functions the refresh runs to build the ClinVar
// annotation and its statistics; see the list in clinvar_gene_refresh.py.
const INSTALLER_PROCEDURES: [&str; 4] = [
    "ClinVar_VCF_deploy",
    "_ClinVar_VCF_deploy",
    "ClinVar_AA_stat",
    "ClinVar_Gene_stat_regenerate",
];
const PLUGINS: [&str; 9] = [
    "UTRAnnotator",
    "LOEUF",
    "AlphaMissense",
    "SpliceAI",
    "PrimateAI",
    "SpliceVault",
    "Conservation",
    "NMD",
    "MaveDB",
];
const PACKAGES: [&str; 10] = [
    "python-",
    "pandas-",
    "numpy-",
    "pysam-",
    "pyarrow-",
    "ensembl-vep-",
    "htslib-",
    "bcftools-",
    "samtools-",
    "perl-",
];

/// Whether a failure is one of the classes `except (OSError, ValueError,
/// KeyError, TypeError)` catches.
fn policy_failure(category: &str) -> bool {
    matches!(
        category,
        "OSError"
            | "FileNotFoundError"
            | "PermissionError"
            | "IsADirectoryError"
            | "NotADirectoryError"
            | "ValueError"
            | "UnicodeDecodeError"
            | "JSONDecodeError"
            | "KeyError"
            | "TypeError"
    )
}

/// `re.findall(r"\{([^}]+)\}", text)`
fn braced_names(text: &str) -> Vec<&str> {
    let mut names = Vec::new();
    let mut at = 0;
    while let Some(offset) = text[at..].find('{') {
        let start = at + offset;
        match text[start + 1..].find('}') {
            Some(length) if length > 0 => {
                names.push(&text[start + 1..start + 1 + length]);
                at = start + 2 + length;
            }
            _ => at = start + 1,
        }
    }
    names
}

/// `resource_policy(source_root)`: the active source and refresh choices.
pub(super) fn resource_policy(source_root: &str) -> Result<Json> {
    let path = join(source_root, "cache_resources.json");
    let loaded = (|| -> Result<Json> {
        let value = read_json(&path)?;
        if !matches!(value, Json::Object(_))
            || !lookup(&value, "schema_version")?
                .is_some_and(|version| json_text::py_equal(version, &Json::Int(1)))
        {
            return Err(value_error("expected schema_version 1"));
        }
        let resources = item(&value, "resources")?;
        let refresh = item(item(resources, "clinvar_gene_statistics")?, "refresh")?;
        let months = item(refresh, "threshold_months")?;
        let automatic = item(refresh, "automatic")?;
        if !matches!(automatic, Json::Bool(_))
            || lookup(refresh, "policy")?.and_then(Json::as_str)
                != Some("latest_filedate_vs_bundled_mtime")
        {
            return Err(value_error(
                "expected boolean automatic and policy latest_filedate_vs_bundled_mtime",
            ));
        }
        let latest = item(
            item(item(resources, "clinvar_vcf")?, "source")?,
            "latest_url_template",
        )?;
        let summary = item(
            item(item(resources, "clinvar_variant_summary")?, "source")?,
            "url",
        )?;
        let mane = item(item(item(resources, "mane_select")?, "source")?, "url")?;
        if !matches!(months, Json::Int(number) if *number > 0) {
            return Err(value_error("threshold_months must be a positive integer"));
        }
        let template_valid = latest.as_str().is_some_and(|text| {
            text.contains("{assembly}")
                && braced_names(text) == ["assembly"]
                && text.ends_with("/clinvar.vcf.gz")
        });
        if !template_valid {
            return Err(value_error(
                "latest_url_template must contain {assembly} and end in /clinvar.vcf.gz",
            ));
        }
        for url in [latest, summary, mane] {
            if !url
                .as_str()
                .is_some_and(|text| text.starts_with(&format!("{NCBI}/")))
            {
                return Err(value_error(
                    "ClinVar/MANE source URLs must use official NCBI HTTPS",
                ));
            }
        }
        Ok(object([
            ("threshold_months", months.clone()),
            ("automatic", automatic.clone()),
            ("latest_url_template", latest.clone()),
            ("summary_url", summary.clone()),
            ("mane_url", mane.clone()),
        ]))
    })();
    let fields = match loaded {
        Ok(fields) => fields,
        Err(failure) => {
            let (category, message) = exception(failure.as_ref());
            if !policy_failure(category) {
                return Err(failure);
            }
            return Err(value_error(format!(
                "Cannot load ClinVar refresh policy from {path}: {message}"
            )));
        }
    };
    let mut policy = object([
        ("registry", string(path.clone())),
        ("registry_sha256", string(sha256(&path)?)),
    ]);
    if let Json::Object(items) = fields {
        for (key, value) in items {
            policy.set(&key, value);
        }
    }
    Ok(policy)
}

fn line_start(text: &str, at: usize) -> bool {
    at == 0 || text.as_bytes()[at - 1] == b'\n'
}
/// The end of the first line that is exactly `}` at its start, from `from`.
fn closing_brace(text: &str, from: usize) -> Option<usize> {
    text[from..]
        .match_indices('}')
        .map(|(offset, _)| from + offset)
        .find(|at| line_start(text, *at))
        .map(|at| at + 1)
}
/// `re.findall(r"^function NAME ?\(\) \{\n.*?^\}", text, re.M | re.S)`
fn procedure_matches<'a>(text: &'a str, name: &str) -> Vec<&'a str> {
    let headers = [
        format!("function {name}() {{\n"),
        format!("function {name} () {{\n"),
    ];
    let mut found = Vec::new();
    let mut at = 0;
    while at <= text.len() {
        let next = text[at..]
            .match_indices("function ")
            .map(|(offset, _)| at + offset)
            .find(|start| line_start(text, *start));
        let Some(start) = next else { break };
        let matched = headers.iter().find_map(|header| {
            text[start..]
                .starts_with(header.as_str())
                .then(|| closing_brace(text, start + header.len()))
                .flatten()
        });
        match matched {
            Some(end) => {
                found.push(&text[start..end]);
                at = end;
            }
            None => at = start + 1,
        }
    }
    found
}
/// `procedure_text(path, name)`: one shell function exactly as written.
fn procedure_text(path: &str, name: &str) -> Result<String> {
    let text = read_text(path)?;
    let matches = procedure_matches(&text, name);
    require(
        matches.len() == 1,
        format!(
            "Cannot identify exactly one ClinVar shell procedure {name} in {path} (found {}); \
             update INSTALLER_PROCEDURES with the stage it was renamed to",
            matches.len()
        ),
    )?;
    Ok(matches[0].to_owned())
}
/// `installer_procedures(source_root)`
fn installer_procedures(source_root: &str) -> Result<Json> {
    let path = absolute(Path::new(&join(
        source_root,
        &format!("scripts/{INSTALLER}"),
    )))?;
    if !exists(&path)? {
        return Ok(object([
            ("path", string(path)),
            ("exists", Json::Bool(false)),
        ]));
    }
    let mut hashes = Json::Object(Vec::new());
    for name in INSTALLER_PROCEDURES {
        let text = procedure_text(&path, name)?;
        hashes.set(name, string(crate::sha256::sha256_hex(text.as_bytes())));
    }
    Ok(object([
        ("path", string(path)),
        ("exists", Json::Bool(true)),
        ("sha256", hashes),
    ]))
}

/// The Python runtime the script would run in: `sys.executable`,
/// `sys.version` and `sys.prefix` of `python` from PATH.
fn python_runtime() -> Result<(String, String, String)> {
    static OBSERVED: OnceLock<(String, String, String)> = OnceLock::new();
    if let Some(found) = OBSERVED.get() {
        return Ok(found.clone());
    }
    const OBSERVER: &str =
        "import json,sys\nprint(json.dumps([sys.executable, sys.version, sys.prefix]))";
    let output = Command::new("python")
        .args(["-c", OBSERVER])
        .stdin(Stdio::null())
        .stderr(Stdio::inherit())
        .output()
        .map_err(|e| os_error(e, Path::new("python"), None))?;
    let text = String::from_utf8_lossy(&output.stdout);
    let parsed = json_text::parse(text.trim_end(), false).ok();
    let values: Option<Vec<String>> = match (&output.status.success(), parsed) {
        (true, Some(Json::Array(items))) if items.len() == 3 => items
            .iter()
            .map(|v| v.as_str().map(str::to_owned))
            .collect(),
        _ => None,
    };
    let Some(values) = values else {
        return Err(error(
            "RuntimeError",
            format!(
                "Cannot observe the Python runtime of the ClinVar scientific profile ({})",
                output.status
            ),
        ));
    };
    let found = (values[0].clone(), values[1].clone(), values[2].clone());
    Ok(OBSERVED.get_or_init(|| found).clone())
}

/// `shutil.which(name)` as the string the script records.
fn which(name: &str) -> Option<String> {
    let path = std::env::var("PATH").unwrap_or_else(|_| "/bin:/usr/bin".into());
    let mut seen = Vec::new();
    for folder in path.split(':') {
        if seen.contains(&folder) {
            continue;
        }
        seen.push(folder);
        let candidate = if folder.is_empty() {
            name.to_owned()
        } else if folder.ends_with('/') {
            format!("{folder}{name}")
        } else {
            format!("{folder}/{name}")
        };
        let executable = CString::new(candidate.as_bytes())
            .is_ok_and(|text| unsafe { libc::access(text.as_ptr(), libc::X_OK) } == 0);
        if fs::metadata(&candidate).is_ok_and(|info| !info.is_dir()) && executable {
            return Some(candidate);
        }
    }
    None
}

/// `sorted(Path(folder).glob(pattern))` for a one-level wildcard.
pub(super) fn glob(folder: &str, matches: impl Fn(&str) -> bool) -> Vec<String> {
    let Ok(entries) = fs::read_dir(folder) else {
        return Vec::new();
    };
    let mut names: Vec<String> = entries
        .filter_map(|entry| entry.ok())
        .filter_map(|entry| entry.file_name().to_str().map(str::to_owned))
        .filter(|name| matches(name))
        .collect();
    names.sort();
    names.into_iter().map(|name| join(folder, &name)).collect()
}

struct Recorder {
    files: Json,
    source_data: String,
}

impl Recorder {
    /// `record(path, code=...)`
    fn record(&mut self, path: &str, code: bool) -> Result<()> {
        let path = absolute(Path::new(path))?;
        let resolved = crate::installer_prerequisites::paths::resolve(&path)?;
        let present = exists(&path)?;
        let mut value = object([
            ("path", string(path.clone())),
            ("resolved", string(resolved.clone())),
            ("exists", Json::Bool(present)),
        ]);
        if present {
            let info = stat(&path)?;
            value.set("directory", Json::Bool(is_dir(&path)?));
            let file = is_file(&path)?;
            if file {
                value.set("size", Json::Int(info.len() as i128));
                value.set(
                    "mtime_ns",
                    Json::Int(nanoseconds(info.mtime(), info.mtime_nsec())),
                );
            }
            let bundled = is_relative_to(&resolved, &self.source_data);
            if file && (code || (!bundled && info.len() <= 1024 * 1024)) {
                value.set("sha256", string(sha256(&path)?));
            }
        }
        self.files.set(&path, value);
        Ok(())
    }
}

/// `scientific_profile(config, source_root, policy)`
pub(super) fn scientific_profile(config: &Json, source_root: &str, policy: &Json) -> Result<Json> {
    let mut values = Json::Object(
        PROFILE_KEYS
            .iter()
            .map(|key| Ok(((*key).to_owned(), get_or_null(config, key)?)))
            .collect::<Result<_>>()?,
    );
    // ClinVar_Gene_stat_regenerate reads an unset seed as 0.
    let seed = values
        .get("priva_python_hash_seed")
        .cloned()
        .unwrap_or(Json::Null);
    values.set(
        "priva_python_hash_seed",
        string(if empty(&seed) {
            "0".into()
        } else {
            py_str(&seed)
        }),
    );
    let pext = config
        .get("pext_enabled")
        .map_or_else(|| "false".to_owned(), py_str)
        .to_lowercase()
        == "true";
    if !pext && let Json::Object(items) = &mut values {
        items.retain(|(key, _)| {
            !matches!(
                key.as_str(),
                "pext_mean_bw" | "pext_tissue_bw_dir" | "pext_tissues"
            )
        });
    }
    let mut recorder = Recorder {
        files: Json::Object(Vec::new()),
        source_data: join(
            &crate::installer_prerequisites::paths::resolve(source_root)?,
            "data",
        ),
    };
    let entries = values.as_object().cloned().unwrap_or_default();
    for (key, value) in &entries {
        if matches!(
            key.as_str(),
            "assembly" | "vep_installed_version" | "pext_enabled" | "pext_tissues"
        ) || empty(value)
        {
            continue;
        }
        if let Json::Str(text) = value {
            let mut path = path_of(&string(text.clone()))?;
            if key == "loeuf_prescore" {
                // Match the existing VEP shell's bundled basename resolution.
                path = join(&join(source_root, "data/loeuf"), name_of(&path));
            }
            recorder.record(&path, false)?;
            if is_file(&path)? {
                for suffix in [".tbi", ".csi", ".fai"] {
                    let index = format!("{path}{suffix}");
                    if exists(&index)? {
                        recorder.record(&index, false)?;
                    }
                }
            }
        }
    }
    for name in SCIENTIFIC_SOURCES {
        recorder.record(&join(source_root, &format!("scripts/{name}")), true)?;
    }
    let plugins = get_or_null(config, "vep_plugins_dir")?;
    if !empty(&plugins) {
        let folder = path_of(&plugins)?;
        for name in PLUGINS {
            recorder.record(&join(&folder, &format!("{name}.pm")), true)?;
        }
    }
    let loftee = get_or_null(config, "loftee_repo")?;
    if !empty(&loftee) {
        let found = glob(&path_of(&loftee)?, |name| {
            name.len() >= 3 && (name.ends_with(".pm") || name.ends_with(".pl"))
        });
        for path in found {
            recorder.record(&path, true)?;
        }
    }
    let gnomad = get_or_null(config, "gnomad_vcf_chrX")?;
    if !empty(&gnomad) {
        let template = py_str(item(config, "gnomad_vcf_chrX")?);
        let chromosomes = (1..=22)
            .map(|n| n.to_string())
            .chain(["X".into(), "Y".into()]);
        for chrom in chromosomes {
            let path = pathlib(&template.replace(".chrX.vcf", &format!(".chr{chrom}.vcf")));
            recorder.record(&path, false)?;
            recorder.record(&format!("{path}.tbi"), false)?;
        }
    }
    let cache = get_or_null(config, "vep_cache_dir")?;
    if !empty(&cache) {
        let assembly = match item(config, "assembly")?.as_str() {
            Some("hg19" | "GRCh37") => "GRCh37",
            _ => "GRCh38",
        };
        let release = py_str(item(config, "vep_installed_version")?);
        let folder = join(
            &path_of(&cache)?,
            &format!("homo_sapiens_merged/{release}_{assembly}"),
        );
        for name in ["info.txt", "chr_synonyms.txt"] {
            recorder.record(&join(&folder, name), false)?;
        }
    }
    if pext && !empty(&get_or_null(config, "pext_tissue_bw_dir")?) {
        let mut tissues = get_or_null(config, "pext_tissues")?;
        if !truth(&tissues) {
            tissues = Json::Array(Vec::new());
        }
        if let Json::Str(text) = &tissues {
            tissues = safe_load(text)?;
            if let Json::Str(text) = &tissues {
                tissues = Json::Array(
                    text.split(',')
                        .map(|name| name.trim())
                        .filter(|name| !name.is_empty())
                        .map(string)
                        .collect(),
                );
            }
        }
        if let Json::Array(names) = &tissues {
            let folder = path_of(item(config, "pext_tissue_bw_dir")?)?;
            for tissue in names {
                let file = format!("{}.bw", py_str(tissue));
                recorder.record(&join(&join(&folder, "tissues"), &file), false)?;
            }
        }
    }
    if let Some(executable) = which("vep") {
        recorder.record(&executable, true)?;
    }
    let (executable, version, prefix) = python_runtime()?;
    let records = glob(&join(&prefix, "conda-meta"), |name| {
        name.len() >= 5 && name.ends_with(".json")
    });
    for path in records {
        if PACKAGES
            .iter()
            .any(|prefix| name_of(&path).starts_with(prefix))
        {
            recorder.record(&path, true)?;
        }
    }
    let mut source_policy = Json::Object(Vec::new());
    for key in [
        "latest_url_template",
        "summary_url",
        "mane_url",
        "automatic",
        "threshold_months",
    ] {
        source_policy.set(key, item(policy, key)?.clone());
    }
    Ok(object([
        ("schema", string(PROFILE_SCHEMA)),
        ("config", values),
        ("files", recorder.files),
        ("installer_procedures", installer_procedures(source_root)?),
        (
            "python",
            object([
                ("executable", string(executable)),
                ("version", string(version)),
            ]),
        ),
        ("source_policy", source_policy),
    ]))
}

/// `profile_signature(profile)`: deterministic, without observation times
/// for files whose content hash is recorded.
pub(super) fn profile_signature(profile: Option<&Json>) -> Result<String> {
    let Some(profile @ Json::Object(_)) = profile else {
        return Ok(String::new());
    };
    let mut stable = profile.clone();
    if let Json::Object(items) = &mut stable
        && let Some((_, files)) = items.iter_mut().find(|(key, _)| key == "files")
    {
        let Json::Object(files) = files else {
            return Err(error(
                "AttributeError",
                format!("'{}' object has no attribute 'values'", type_name(files)),
            ));
        };
        for (_, value) in files.iter_mut() {
            if let Json::Object(fields) = value
                && fields.iter().any(|(key, _)| key == "sha256")
            {
                fields.retain(|(key, _)| key != "mtime_ns");
            }
        }
    }
    let text = json_text::dumps_with_options(
        &stable,
        DumpsOptions {
            sort_keys: true,
            compact: true,
            ..DumpsOptions::default()
        },
    )
    .map_err(value_error)?;
    Ok(crate::sha256::sha256_hex(text.as_bytes()))
}

/// `candidate_route(source_root)`: source wiring only; never qualifies a profile.
pub(super) fn candidate_route(source_root: &str) -> Result<Json> {
    let script = join(source_root, "scripts/annotation_vcf.sh");
    let text = read_text(&script)?;
    let start = |name: &str| {
        let header = format!("function {name}() {{");
        text.match_indices(&header)
            .map(|(at, _)| at)
            .find(|at| line_start(&text, *at))
            .map(|at| (at, at + header.len()))
    };
    for name in [
        "anno_VEP_data",
        "anno_VEP_cache_candidate",
        "_priva_anno_VEP_data",
        "priva_check_spliceai_versions",
    ] {
        require(
            start(name).is_some(),
            format!("ClinVar candidate preparation entry point is missing: {name}"),
        )?;
    }
    let guard = text
        .match_indices("function priva_check_spliceai_versions() {")
        .map(|(at, _)| at)
        .filter(|at| line_start(&text, *at))
        .find_map(|at| {
            closing_brace(
                &text,
                at + "function priva_check_spliceai_versions() {".len(),
            )
            .map(|end| &text[at..end])
        });
    let Some(guard) = guard else {
        return Err(error("TypeError", "'NoneType' object is not subscriptable"));
    };
    Ok(object([
        ("kind", string(CANDIDATE_KIND)),
        (
            "status",
            string("preparation_route_available_inputs_not_checked"),
        ),
        ("scientific_qualified", Json::Bool(false)),
        ("source", identity(&script)?),
        (
            "helper_source",
            identity(&join(source_root, "scripts/clinvar_gene_refresh.py"))?,
        ),
        (
            "runtime_guard_sha256",
            string(crate::sha256::sha256_hex(guard.as_bytes())),
        ),
    ]))
}

fn candidate_receipt_path(vcf: &str, pending: bool) -> String {
    format!(
        "{vcf}{}",
        if pending {
            ".installation-candidate.pending.json"
        } else {
            ".installation-candidate.ready.json"
        }
    )
}

/// A SpliceAI cache-input failure as the Python exception it stands for.
fn spliceai_failure(failure: Failure) -> Box<dyn std::error::Error> {
    match failure {
        Failure::Rejected(message) => {
            let category = match message.split_once(']').map(|(code, _)| code) {
                Some("[Errno 2") => "FileNotFoundError",
                Some("[Errno 13") => "PermissionError",
                Some("[Errno 20") => "NotADirectoryError",
                Some("[Errno 21") => "IsADirectoryError",
                Some(code) if code.starts_with("[Errno ") => "OSError",
                _ => "ValueError",
            };
            error(category, message)
        }
        Failure::Crash(message) => {
            const KNOWN: [&str; 5] = [
                "TypeError",
                "KeyError",
                "AttributeError",
                "UnicodeDecodeError",
                "JSONDecodeError",
            ];
            match message
                .split_once(": ")
                .and_then(|(name, rest)| KNOWN.iter().find(|k| **k == name).map(|k| (*k, rest)))
            {
                Some((category, rest)) => error(category, rest),
                None => error("RuntimeError", message),
            }
        }
        Failure::Helper(code) => error("RuntimeError", format!("helper exit status {code}")),
    }
}

/// `candidate_context(config_path, source_root, input_vcf)`
fn candidate_context(config_path: &str, source_root: &str, input_vcf: &str) -> Result<Json> {
    let config = load_config(config_path)?;
    let root_value = get_or_null(&config, "cache_root")?;
    let root_text = if empty(&root_value) {
        None
    } else {
        Some(path_of(&root_value)?).filter(|path| path.starts_with('/'))
    };
    let Some(root_text) = root_text else {
        return Err(value_error(
            "ClinVar candidate preparation requires an absolute configured cache_root",
        ));
    };
    let root = crate::installer_prerequisites::paths::resolve(&root_text)?;
    let Some(assembly) = assembly_of(&get_or_null(&config, "assembly")?)? else {
        return Err(value_error(
            "ClinVar candidate preparation requires a configured assembly",
        ));
    };
    let folder = path_of(item(&config, "clinvar_vcf_dir")?)?;
    let selected = selected_vcf(&config, &folder, ucsc(assembly))?;
    let vcf = absolute(Path::new(input_vcf))?;
    let resolve = crate::installer_prerequisites::paths::resolve;
    let outside = !folder.starts_with('/')
        || !is_relative_to(&resolve(&folder)?, &root)
        || !is_relative_to(&vcf, &absolute(Path::new(&root_text))?)
        || !is_relative_to(&resolve(&vcf)?, &resolve(&folder)?)
        || resolve(&vcf)? != resolve(&selected)?;
    require(
        !outside,
        "ClinVar candidate output must be the selected managed VCF inside cache_root/clinvar_vcf_dir",
    )?;
    // The original shell replaces the first '.vcf' suffix; require its output to
    // be exactly this owned file rather than allow an untracked alternate path.
    require(
        vcf.ends_with(".vcf.gz") && format!("{}.vcf.gz", first_vcf_stem(&vcf)) == vcf,
        "ClinVar candidate preparation requires a managed .vcf.gz path without earlier .vcf components",
    )?;
    nonempty(&vcf, "managed ClinVar candidate input")?;
    let snv = get_or_null(&config, "spliceai_snv_prescore")?;
    let indel = get_or_null(&config, "spliceai_indel_prescore")?;
    require(
        !empty(&snv) && !empty(&indel),
        "ClinVar candidate preparation requires both selected SpliceAI VCFs",
    )?;
    let (snv, indel) = (path_of(&snv)?, path_of(&indel)?);
    let format = spliceai::pair_format(assembly, Path::new(&snv), Path::new(&indel))
        .map_err(spliceai_failure)?;
    require(
        matches!(format, "vcf_v1" | "vcf_spliceai2"),
        "ClinVar candidate preparation requires installed SpliceAI VCFs, not native publisher TSVs",
    )?;
    let mut adapter = Json::Null;
    if format == "vcf_spliceai2" {
        let receipt = join(&parent_of(&resolve(&snv)?), "adapter.ready.json");
        adapter = spliceai::within_completed_projection(Path::new(&snv), Path::new(&indel), || {
            identity(&receipt)
        })
        .map_err(spliceai_failure)??;
    }
    let config_sha256 = sha256(config_path)?;
    let route = candidate_route(source_root)?;
    let policy = resource_policy(source_root)?;
    Ok(object([
        ("kind", string(CANDIDATE_KIND)),
        ("config_path", string(absolute(Path::new(config_path))?)),
        ("config_sha256", string(config_sha256)),
        ("source_root", string(absolute(Path::new(source_root))?)),
        ("input_vcf", string(vcf)),
        ("assembly", string(assembly)),
        ("spliceai_format", string(format)),
        ("adapter_receipt", adapter),
        ("route", route),
        (
            "scientific_profile",
            scientific_profile(&config, source_root, &policy)?,
        ),
        ("scientific_qualified", Json::Bool(false)),
    ]))
}

/// `str(path).split(".vcf", 1)[0]`
pub(super) fn first_vcf_stem(path: &str) -> &str {
    path.split_once(".vcf").map_or(path, |(stem, _)| stem)
}

/// `candidate_begin(config_path, source_root, input_vcf)`
pub(super) fn candidate_begin(
    config_path: &str,
    source_root: &str,
    input_vcf: &str,
) -> Result<Json> {
    let mut value = candidate_context(config_path, source_root, input_vcf)?;
    let vcf = text_of(item(&value, "input_vcf")?);
    let info = stat(&vcf)?;
    value.set("status", string("preparing_annotation_only"));
    value.set("started_utc", string(utc_now()));
    value.set(
        "input_before",
        object([
            ("path", string(vcf.clone())),
            ("size", Json::Int(info.len() as i128)),
            (
                "mtime_ns",
                Json::Int(nanoseconds(info.mtime(), info.mtime_nsec())),
            ),
        ]),
    );
    write_json(&candidate_receipt_path(&vcf, true), &value)?;
    Ok(value)
}

/// `candidate_complete(config_path, source_root, input_vcf)`
pub(super) fn candidate_complete(
    config_path: &str,
    source_root: &str,
    input_vcf: &str,
) -> Result<Json> {
    let vcf = absolute(Path::new(input_vcf))?;
    let pending = candidate_receipt_path(&vcf, true);
    let before = read_json(&pending)?;
    let current = candidate_context(config_path, source_root, &vcf)?;
    let status = lookup(&before, "status")?.cloned().unwrap_or(Json::Null);
    let changed = !json_text::py_equal(&status, &string("preparing_annotation_only"))
        || current
            .as_object()
            .into_iter()
            .flatten()
            .any(|(key, value)| {
                !json_text::py_equal(before.get(key).unwrap_or(&Json::Null), value)
            });
    require(
        !changed,
        "ClinVar candidate configuration/scientific inputs changed during preparation",
    )?;
    let contract = crate::clinvar_contract::header::validate(&vcf, "vep")?;
    let mut value = before;
    value.set("status", string("prepared_annotation_only"));
    value.set("completed_utc", string(utc_now()));
    value.set("stage", string("vep_before_gnomad_enrichment"));
    value.set("output", identity(&vcf)?);
    value.set("output_index", identity(&format!("{vcf}.tbi"))?);
    value.set("annotation_contract", contract);
    value.set(
        "qualification",
        string(
            "installation candidate only; ACMG/ClinVar/CADD and public/family qualification remain pending",
        ),
    );
    write_json(&candidate_receipt_path(&vcf, false), &value)?;
    unlink(&pending)?;
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn braces_follow_the_findall_pattern() {
        assert_eq!(braced_names("a/{assembly}/b"), ["assembly"]);
        assert_eq!(braced_names("{}{assembly}"), ["assembly"]);
        assert_eq!(braced_names("{{assembly}}"), ["{assembly"]);
        assert!(braced_names("{open").is_empty());
    }

    #[test]
    fn shell_procedures_end_at_the_first_closing_line() {
        let text = "x\nfunction a() {\n  b\n}\nfunction a () {\n}\nfunction ab() {\n}\n";
        assert_eq!(
            procedure_matches(text, "a"),
            ["function a() {\n  b\n}", "function a () {\n}"]
        );
        assert_eq!(procedure_matches("function a() {\nx", "a").len(), 0);
    }
}
