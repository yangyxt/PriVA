//! The four task-contract commands: manifest, validate, launch and run.
use super::contract::{
    self, Config, OPTIONAL_INPUTS, assembly, cds_cache_path, copy_writable, copy2, exists,
    inventory, inventory_record, json_digest, object, read_config, read_json, read_text,
    render_config, runtime_resources, safe_id, sha256, source_digest, text, verify_manifest,
    verify_manifest_for, write_json,
};
use super::{Options, Result, refuse};
use crate::json_text::Json;
use crate::python_path::path_value;
use crate::yaml_config::{ConfigValue, Scalar};
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

fn resolve(value: &str) -> Result<PathBuf> {
    Ok(fs::canonicalize(path_value(value))?)
}

fn display(path: &Path) -> String {
    path.display().to_string()
}

pub(crate) fn build_manifest(options: &Options) -> Result<()> {
    let config_path = PathBuf::from(options.value("config"));
    let config = read_config(&config_path)?;
    let root = contract::checked_path(Some(options.value("legacy-root")), "legacy_root")?;
    let root = fs::canonicalize(root)?;
    let mode = options.get("mode").unwrap_or("raw");
    let vcf_hub = options.get("vcf-hub-mode").unwrap_or("off");
    let cadd_hub = options.get("cadd-hub-mode").unwrap_or("copy");
    let checksum = options.flag("checksum");
    let paths = runtime_resources(&config, &root, mode, vcf_hub, cadd_hub)?;
    let mut payload = object(vec![
        ("schema_version", Json::Int(1)),
        ("assembly", text(config.assembly())),
        ("mode", text(mode)),
        ("config_sha256", text(sha256(&config_path)?)),
        ("legacy_root", text(display(&root))),
        (
            "identity_kind",
            text(if checksum {
                "sha256"
            } else {
                "size-mtime-inventory+executable-sha256"
            }),
        ),
        ("vcf_hub_mode", text(vcf_hub)),
        ("cadd_hub_mode", text(cadd_hub)),
        ("resources", inventory(&paths, checksum)?),
    ]);
    let id = json_digest(&payload)?;
    payload.set("bundle_id", text(id));
    write_json(Path::new(options.value("output")), &payload)
}

/// The opt-in read-only view is part of the run identity. This detects an
/// omitted profile; it is not an attestation of scientific qualification.
fn require_task_boundary(boundary: &str) -> Result<()> {
    if boundary != "off" && boundary != "bwrap" {
        return refuse("Unsupported task boundary");
    }
    if boundary == "bwrap"
        && std::env::var("PRIVA_TASK_BOUNDARY").ok().as_deref() != Some("bwrap-readonly-v1")
    {
        return refuse("The requested read-only task boundary is not active");
    }
    Ok(())
}

/// `shutil.which`: the first executable, non-directory PATH entry.
fn which(tool: &str) -> Option<PathBuf> {
    let search = std::env::var("PATH").unwrap_or_default();
    search.split(':').find_map(|directory| {
        let directory = if directory.is_empty() { "." } else { directory };
        let candidate = Path::new(directory).join(tool);
        let usable = fs::metadata(&candidate).is_ok_and(|m| !m.is_dir())
            && contract::access_executable(&candidate);
        usable.then_some(candidate)
    })
}

/// Run `tool --version` with a 30 s limit; return the first 4096 characters
/// of standard output followed by standard error, and the exit code.
fn version_output(tool: &str) -> Result<(String, i32)> {
    let mut child = Command::new(tool)
        .arg("--version")
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;
    let mut stdout = child.stdout.take().unwrap();
    let mut stderr = child.stderr.take().unwrap();
    let out = std::thread::spawn(move || {
        let mut bytes = Vec::new();
        stdout.read_to_end(&mut bytes).map(|_| bytes)
    });
    let err = std::thread::spawn(move || {
        let mut bytes = Vec::new();
        stderr.read_to_end(&mut bytes).map(|_| bytes)
    });
    let started = Instant::now();
    let status = loop {
        if let Some(status) = child.try_wait()? {
            break status;
        }
        if started.elapsed() > Duration::from_secs(30) {
            let _ = child.kill();
            let _ = child.wait();
            return Err(
                format!("Command '['{tool}', '--version']' timed out after 30 seconds").into(),
            );
        }
        std::thread::sleep(Duration::from_millis(20));
    };
    let out = out.join().map_err(|_| "reader failed")??;
    let err = err.join().map_err(|_| "reader failed")??;
    let combined = String::from_utf8(out)? + &String::from_utf8(err)?;
    let code = status.code().unwrap_or_else(|| {
        use std::os::unix::process::ExitStatusExt;
        -status.signal().unwrap_or(0)
    });
    Ok((combined.chars().take(4096).collect(), code))
}

fn software_readiness(mode: &str, priva_binary: Option<&Path>) -> Result<Json> {
    let mut required = vec!["bash", "bcftools", "tabix", "yq", "parallel", "tsv-pretty"];
    // PriVA's own steps run in priva. Only raw-mode annotation runs Python:
    // scripts/cadd_conda.py finds CADD's Conda, and CADD.sh runs Snakemake.
    if mode == "raw" {
        required.extend(["vep", "samtools", "conda", "python", "snakemake"]);
    }
    let missing: Vec<&str> = required
        .iter()
        .copied()
        .filter(|tool| which(tool).is_none())
        .collect();
    if !missing.is_empty() {
        return refuse(format!(
            "Runtime tools missing from PATH: {}",
            missing.join(", ")
        ));
    }
    let mut versions = Json::Object(Vec::new());
    for tool in &required {
        let found = which(tool).unwrap();
        let mut entry = object(vec![(
            "executable",
            text(display(&fs::canonicalize(&found).unwrap_or(found))),
        )]);
        if ["python", "bcftools", "tabix", "yq", "samtools", "conda"].contains(tool) {
            let (output, code) = version_output(tool)?;
            entry.set("version_output", text(output));
            entry.set("version_exit_code", Json::Int(code as i128));
        }
        versions.set(tool, entry);
    }
    // The modules cadd_conda.py and CADD's Snakemake import, in the python
    // the shell runs (not necessarily python3).
    if mode == "raw" {
        let status = Command::new("python")
            .args(["-c", "import packaging, snakemake"])
            .status()?;
        if !status.success() {
            return Err(called_process_error(
                &["python", "-c", "import packaging, snakemake"],
                status,
            ));
        }
    }
    if let Some(binary) = priva_binary {
        // The same check as priva_resolve_binary in scripts/common_bash_utils.sh.
        let output = Command::new(binary)
            .arg("--version")
            .stdin(Stdio::null())
            .output()?;
        let stdout = String::from_utf8_lossy(&output.stdout).into_owned();
        if !output.status.success() || !stdout.starts_with("priva ") {
            return refuse("priva_binary is not a working PriVA program (priva --version failed)");
        }
        versions.set(
            "priva",
            object(vec![
                ("executable", text(display(&fs::canonicalize(binary)?))),
                (
                    "version_output",
                    text(stdout.chars().take(4096).collect::<String>()),
                ),
                ("version_exit_code", Json::Int(0)),
            ]),
        );
    }
    Ok(versions)
}

fn called_process_error(
    command: &[&str],
    status: std::process::ExitStatus,
) -> Box<dyn std::error::Error> {
    use std::os::unix::process::ExitStatusExt;
    let listed = command
        .iter()
        .map(|part| crate::python_path::repr(part))
        .collect::<Vec<_>>()
        .join(", ");
    match (status.code(), status.signal()) {
        (Some(code), _) => {
            format!("Command '[{listed}]' returned non-zero exit status {code}.").into()
        }
        (None, Some(signal)) => format!(
            "Command '[{listed}]' died with {}.",
            crate::python_subprocess::signal_repr(signal)
        )
        .into(),
        _ => format!("Command '[{listed}]' failed.").into(),
    }
}

/// The sample columns of the input VCF.
fn vcf_samples(path: &Path) -> Result<Vec<String>> {
    let file = fs::File::open(path).map_err(|e| crate::python_path::os_error(e, path, None))?;
    let mut lines = crate::python_lines::Lines::new(flate2::read::MultiGzDecoder::new(
        std::io::BufReader::new(file),
    ));
    while let Some(line) = lines.next()? {
        if line.starts_with("#CHROM\t") {
            let fields: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
            if fields.len() < 10 || fields[8] != "FORMAT" {
                return refuse("Family input VCF must have FORMAT and sample columns");
            }
            return Ok(fields[9..].iter().map(|s| s.to_string()).collect());
        }
        if !line.starts_with('#') {
            break;
        }
    }
    refuse("Input VCF has no #CHROM header")
}

/// `str.splitlines()` boundaries.
fn split_lines(text: &str) -> Vec<&str> {
    let mut lines = Vec::new();
    let mut start = 0;
    let mut chars = text.char_indices().peekable();
    while let Some((at, c)) = chars.next() {
        let boundary = matches!(
            c,
            '\n' | '\r'
                | '\u{0b}'
                | '\u{0c}'
                | '\u{1c}'
                | '\u{1d}'
                | '\u{1e}'
                | '\u{85}'
                | '\u{2028}'
                | '\u{2029}'
        );
        if boundary {
            lines.push(&text[start..at]);
            let mut end = at + c.len_utf8();
            if c == '\r' {
                if let Some((next_at, '\n')) = chars.peek().copied() {
                    chars.next();
                    end = next_at + 1;
                }
            }
            start = end;
        }
    }
    if start < text.len() {
        lines.push(&text[start..]);
    }
    lines
}

fn python_strip(text: &str) -> &str {
    text.trim_matches(crate::python_text::is_python_space)
}

/// Each family's sample identifiers, in PED order, for every family in
/// order of first appearance, or only `family` when one is selected.
pub(crate) fn validate_ped(
    path: &Path,
    family: Option<&str>,
    samples: &[String],
) -> Result<Vec<(String, Vec<String>)>> {
    let content = read_text(path)?;
    let rows: Vec<Vec<&str>> = split_lines(&content)
        .into_iter()
        .filter(|line| {
            !python_strip(line).is_empty()
                && !line
                    .trim_start_matches(crate::python_text::is_python_space)
                    .starts_with('#')
        })
        .map(|line| line.split('\t').collect())
        .collect();
    if rows.is_empty() || rows.iter().any(|row| row.len() < 6) {
        return refuse(
            "Legacy Python requires a tab-delimited PED with at least six columns and one person",
        );
    }
    let mut order: Vec<&str> = Vec::new();
    for row in &rows {
        if !order.contains(&row[0]) {
            order.push(row[0]);
        }
    }
    if let Some(family) = family {
        if !order.contains(&family) {
            return refuse("The selected family is absent from the PED");
        }
        order = vec![family];
    }
    if order.iter().any(|id| !safe_id(id)) {
        return refuse("PED family identifier contains unsupported shell characters");
    }
    let ids: Vec<&str> = rows.iter().map(|row| row[1]).collect();
    let unique = |items: &[&str]| {
        let mut sorted = items.to_vec();
        sorted.sort_unstable();
        sorted.windows(2).all(|pair| pair[0] != pair[1])
    };
    let sample_refs: Vec<&str> = samples.iter().map(String::as_str).collect();
    if !unique(&ids) || !unique(&sample_refs) {
        return refuse("Duplicate sample identifier in PED or VCF");
    }
    if ids.iter().any(|id| !safe_id(id)) {
        return refuse("PED sample identifier contains unsupported shell characters");
    }
    if ids.iter().any(|id| !sample_refs.contains(id)) {
        return refuse("Every PED sample must be present in the VCF");
    }
    // Missing affection is neither affected nor unaffected in the legacy science.
    if rows
        .iter()
        .any(|row| !["0", "1", "2"].contains(&row[4]) || !["-9", "0", "1", "2"].contains(&row[5]))
    {
        return refuse("PED requires sex 0/1/2 and phenotype -9/0/1/2");
    }
    let mut families = Vec::new();
    for id in order {
        let selected: Vec<&Vec<&str>> = rows.iter().filter(|row| row[0] == id).collect();
        if selected[0][5] != "2" {
            return refuse(format!(
                "Legacy PriVA requires the affected proband first (family {id}); PED order is never rewritten"
            ));
        }
        families.push((
            id.to_owned(),
            selected.iter().map(|row| row[1].to_owned()).collect(),
        ));
    }
    Ok(families)
}

fn validate_cadd(path: &Path) -> Result<()> {
    let file = fs::File::open(path).map_err(|e| crate::python_path::os_error(e, path, None))?;
    let mut lines = crate::python_lines::Lines::new(file);
    let first = lines.next()?.unwrap_or_default();
    let fields: Vec<&str> = first.trim_end_matches('\n').split('\t').collect();
    let required = ["#Chrom", "Pos", "Ref", "Alt", "PHRED", "FeatureID"];
    let mut sorted = fields.clone();
    sorted.sort_unstable();
    let duplicate = sorted.windows(2).any(|pair| pair[0] == pair[1]);
    if !required.iter().all(|name| fields.contains(name)) || duplicate {
        return refuse(
            "CADD first row must have unique #Chrom, Pos, Ref, Alt, PHRED, FeatureID columns; remove no rows implicitly",
        );
    }
    Ok(())
}

pub(crate) fn validate(options: &Options) -> Result<()> {
    let mode = options.value("mode");
    let stub = options.flag("stub");
    let boundary = options.get("task-boundary").unwrap_or("off");
    require_task_boundary(boundary)?;
    let config_path = PathBuf::from(options.value("config"));
    let config = read_config(&config_path)?;
    let family = options.get("family").filter(|v| !v.is_empty());
    let base = options.value("base");
    if family.is_some_and(|f| !safe_id(f)) || !safe_id(base) {
        return refuse("Unsafe family/base identifier");
    }
    let ped = options
        .get("ped")
        .filter(|v| !v.is_empty())
        .map(PathBuf::from);
    if family.is_some() && ped.is_none() {
        return refuse("A selected family needs its PED");
    }
    let requested = assembly(options.value("assembly"))
        .ok_or_else(|| format!("KeyError: '{}'", options.value("assembly")))?;
    if config.assembly() != requested {
        return refuse("Requested assembly disagrees with legacy config");
    }
    let cadd = options.get("cadd").filter(|v| !v.is_empty());
    if mode == "annotated" && cadd.is_none() {
        return refuse("Annotated mode requires the paired CADD TSV");
    }
    let vcf = PathBuf::from(options.value("vcf"));
    let samples = vcf_samples(&vcf)?;
    // No PED: one run over every sample without family filtering, as the
    // shell workflow's no-pedigree mode.
    let families = match &ped {
        Some(ped) => validate_ped(ped, family, &samples)?,
        None => Vec::new(),
    };
    if !Path::new(options.value("vcf-index")).is_file() {
        return refuse("Input VCF index is missing");
    }
    if let Some(cadd) = cadd {
        validate_cadd(Path::new(cadd))?;
    }
    let resources = PathBuf::from(options.value("resources"));
    let manifest = read_json(&resources)?;
    let root = resolve(options.value("legacy-root"))?;
    if source_digest(&root)? != options.value("source-digest") {
        return refuse("Source changed after Nextflow launch; freeze the checkout before retrying");
    }
    let vcf_hub = options.value("vcf-hub-mode");
    let cadd_hub = options.value("cadd-hub-mode");
    let mut versions = Json::Object(Vec::new());
    if !stub {
        let paths = verify_manifest(&manifest, &config_path, &root, mode, vcf_hub, cadd_hub)?;
        let binary = paths
            .iter()
            .find(|(k, _)| k == "priva_binary")
            .map(|(_, p)| p.clone());
        versions = software_readiness(mode, binary.as_deref())?;
        let status = Command::new("bcftools")
            .args(["index", "--nrecords"])
            .arg(&vcf)
            .stdout(Stdio::null())
            .status()?;
        if !status.success() {
            let vcf_text = display(&vcf);
            return Err(called_process_error(
                &["bcftools", "index", "--nrecords", &vcf_text],
                status,
            ));
        }
    }
    let bundle = Path::new("bundle");
    fs::create_dir(bundle).map_err(|e| crate::python_path::os_error(e, bundle, None))?;
    for (source, name) in [
        (config_path.as_path(), "legacy.yaml"),
        (resources.as_path(), "resources.json"),
    ] {
        copy2(source, &bundle.join(name))?;
    }
    if let Some(ped) = &ped {
        copy2(ped, &bundle.join("family.ped"))?;
    }
    let samples_json =
        |items: &[String]| Json::Array(items.iter().cloned().map(Json::Str).collect());
    let bundle_id = match &manifest {
        Json::Object(_) => manifest.get("bundle_id").cloned().unwrap_or(text("stub")),
        _ => return Err("Resource manifest must be a JSON object".into()),
    };
    let provided = match cadd {
        Some(cadd) => text(sha256(Path::new(cadd))?),
        None => Json::Null,
    };
    let record = object(vec![
        ("schema_version", Json::Int(1)),
        (
            "families",
            Json::Array(families.iter().map(|(id, _)| text(id.clone())).collect()),
        ),
        ("base", text(base)),
        ("assembly", text(requested)),
        ("mode", text(mode)),
        ("legacy_root", text(display(&root))),
        ("source_digest", text(options.value("source-digest"))),
        ("software_id", text(options.value("software-id"))),
        ("resource_bundle_id", bundle_id),
        ("task_boundary", text(boundary)),
        ("vcf_hub_mode", text(vcf_hub)),
        ("cadd_hub_mode", text(cadd_hub)),
        ("samples", samples_json(&samples)),
        (
            "family_samples",
            Json::Object(
                families
                    .iter()
                    .map(|(id, members)| (id.clone(), samples_json(members)))
                    .collect(),
            ),
        ),
        ("stub", Json::Bool(stub)),
        ("input_sha256", text(sha256(&vcf)?)),
        (
            "ped_sha256",
            match &ped {
                Some(ped) => text(sha256(ped)?),
                None => Json::Null,
            },
        ),
        ("provided_cadd_sha256", provided),
    ]);
    write_json(&bundle.join("contract.json"), &record)?;
    write_json(&bundle.join("software-observed.json"), &versions)
}

/// Keep scientific worker settings distinct from the task allocation.
fn task_thread_settings(
    config: &Config,
    phase: &str,
    mode: &str,
    cpus: i64,
) -> Result<Vec<(String, Json)>> {
    if cpus < 1 {
        return refuse("Allocated cpus must be a positive integer");
    }
    // `config.get(key, cpus)`: the configured value, else the allocation.
    let setting = |key: &str| -> (Json, String, Option<u128>) {
        match config.get(key) {
            None => (
                Json::Int(cpus as i128),
                cpus.to_string(),
                Some(cpus as u128),
            ),
            Some(ConfigValue::Scalar(scalar)) => {
                let shown = match scalar {
                    Scalar::Str(text) => text.clone(),
                    other => other.python_repr(),
                };
                let digits = match scalar {
                    Scalar::Int(text) | Scalar::Str(text)
                        if !text.is_empty() && text.bytes().all(|b| b.is_ascii_digit()) =>
                    {
                        text.parse::<u128>().ok()
                    }
                    _ => None,
                };
                let json = match scalar {
                    Scalar::Null => Json::Null,
                    Scalar::Bool(flag) => Json::Bool(*flag),
                    Scalar::Int(text) => text.parse::<i128>().map_or(Json::Null, Json::Int),
                    Scalar::Str(text) => Json::Str(text.clone()),
                    Scalar::Other(text) => Json::Str(text.clone()),
                };
                (json, shown, digits)
            }
            Some(ConfigValue::Opaque(text)) => (Json::Null, text.clone(), None),
        }
    };
    let mut settings = Vec::new();
    for key in ["threads", "threads_per_fam"] {
        // Only keys absent from the configuration are written; present
        // values are already in the text the task configuration keeps.
        if config.get(key).is_none() {
            settings.push((key.to_owned(), setting(key).0));
        }
    }
    // VEP runs in the whole raw annotation and in its own step. The other
    // annotation steps, and annotated mode's import, need no CPU minimum.
    let vep = phase == "annotation-vep" || (phase == "annotation" && mode == "raw");
    if phase.starts_with("annotation") && !vep {
        return Ok(settings);
    }
    let key = if vep { "threads" } else { "threads_per_fam" };
    if vep && config.get(key).is_none() {
        return refuse("Raw annotation requires explicit legacy-config threads");
    }
    let (_, shown, digits) = setting(key);
    let Some(value) = digits.filter(|v| *v >= 1) else {
        return refuse(format!("Legacy-config {key} must be a positive integer"));
    };
    // Preserve the fixed scientific profile: allocating 24 CPUs must not
    // rewrite a qualified --fork 22 invocation. VEP 116 can spawn N+1
    // children besides its parent; reserve that headroom for raw annotation.
    let required = value + if vep { 2 } else { 0 };
    if (cpus as u128) < required {
        return refuse(format!(
            "{phase} needs at least {required} allocated CPUs to preserve {key}={shown}"
        ));
    }
    Ok(settings)
}

pub(crate) fn launch(options: &Options) -> Result<()> {
    let mode = options.value("mode");
    let stub = options.flag("stub");
    let root = resolve(options.value("legacy-root"))?;
    if source_digest(&root)? != options.value("source-digest") {
        return refuse("Source changed while launching");
    }
    let config_path = PathBuf::from(options.value("config"));
    let config = read_config(&config_path)?;
    let requested = assembly(options.value("assembly"))
        .ok_or_else(|| format!("KeyError: '{}'", options.value("assembly")))?;
    if config.assembly() != requested {
        return refuse("Requested assembly disagrees with legacy config");
    }
    let annotation_cpus = options.int("annotation-cpus")?;
    let family_cpus = options.int("family-cpus")?;
    if !stub {
        task_thread_settings(&config, "annotation", mode, annotation_cpus)?;
        task_thread_settings(&config, "filtration", mode, family_cpus)?;
    }
    let resources = PathBuf::from(options.value("resources"));
    let manifest = read_json(&resources)?;
    let vcf_hub = options.value("vcf-hub-mode");
    let cadd_hub = options.value("cadd-hub-mode");
    if !stub {
        verify_manifest(&manifest, &config_path, &root, mode, vcf_hub, cadd_hub)?;
    }
    let cadd = options.get("cadd").filter(|v| !v.is_empty());
    let provided = match cadd {
        Some(cadd) => inventory_record(&resolve(cadd)?, false)?,
        None => Json::Null,
    };
    let fingerprint = object(vec![
        ("schema_version", Json::Int(1)),
        ("source_digest", text(options.value("source-digest"))),
        ("software_id", text(options.value("software-id"))),
        (
            "family",
            options
                .get("family")
                .filter(|v| !v.is_empty())
                .map_or(Json::Null, text),
        ),
        ("assembly", text(requested)),
        ("mode", text(mode)),
        ("config_sha256", text(sha256(&config_path)?)),
        (
            "ped_sha256",
            match options.get("ped").filter(|v| !v.is_empty()) {
                Some(ped) => text(sha256(Path::new(ped))?),
                None => Json::Null,
            },
        ),
        ("manifest_sha256", text(sha256(&resources)?)),
        (
            "input",
            inventory_record(&resolve(options.value("vcf"))?, false)?,
        ),
        (
            "input_index",
            inventory_record(&resolve(options.value("vcf-index"))?, false)?,
        ),
        ("provided_cadd", provided),
        ("vcf_hub_mode", text(vcf_hub)),
        ("cadd_hub_mode", text(cadd_hub)),
        ("annotation_cpus", Json::Int(annotation_cpus as i128)),
        ("family_cpus", Json::Int(family_cpus as i128)),
        ("stub", Json::Bool(stub)),
        (
            "task_boundary",
            text(options.get("task-boundary").unwrap_or("off")),
        ),
    ]);
    let marker = object(vec![
        ("run_identity", text(json_digest(&fingerprint)?)),
        ("contract", fingerprint),
    ]);
    let outdir = PathBuf::from(crate::python_path::absolute(&path_value(
        options.value("outdir"),
    ))?);
    fs::create_dir_all(&outdir)?;
    let outdir = fs::canonicalize(&outdir)?;
    let path = outdir.join(".priva-run.json");
    if exists(&path) {
        if !crate::json_text::py_equal(&read_json(&path)?, &marker) {
            return refuse(
                "Output directory belongs to different inputs/source/resources/settings; use a new --outdir",
            );
        }
    } else {
        // Nextflow may create trace/report directories before the workflow.
        for child in fs::read_dir(&outdir)? {
            if child?.file_name() != "pipeline_info" {
                return refuse(
                    "Output directory is not empty and has no matching PriVA run identity",
                );
            }
        }
        use std::io::Write;
        let mut handle = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .map_err(|e| crate::python_path::os_error(e, &path, None))?;
        handle.write_all(contract::json_document(&marker)?.as_bytes())?;
    }
    Ok(())
}

fn private_config(
    bundle: &Path,
    record: &Json,
    phase: &str,
    cpus: i64,
) -> Result<(PathBuf, PathBuf)> {
    let config = read_config(&bundle.join("legacy.yaml"))?;
    let mode = record
        .get("mode")
        .and_then(Json::as_str)
        .unwrap_or_default();
    let thread_settings = task_thread_settings(&config, phase, mode, cpus)?;
    let result = PathBuf::from(crate::python_path::absolute(Path::new("result"))?);
    if !result.is_dir() {
        fs::create_dir(&result)?;
    }
    let temp = PathBuf::from(crate::python_path::absolute(Path::new("tmp"))?);
    if !temp.is_dir() {
        fs::create_dir(&temp)?;
    }
    let field = |key: &str| record.get(key).and_then(Json::as_str).unwrap_or_default();
    // An empty ped_file is the shell scripts' no-pedigree mode.
    let ped = if bundle.join("family.ped").is_file() {
        display(&copy_writable(
            &bundle.join("family.ped"),
            &result.join("family.ped"),
        )?)
    } else {
        String::new()
    };
    let mut updates: Vec<(String, Json)> = vec![
        (
            "input_vcf".into(),
            text(display(&result.join(format!("{}.vcf.gz", field("base"))))),
        ),
        ("output_dir".into(), text(display(&result))),
        ("tmp_dir".into(), text(display(&temp))),
        ("ped_file".into(), text(ped)),
    ];
    updates.extend(thread_settings);
    updates.extend([
        ("assembly".into(), text(field("assembly"))),
        ("hub_vcf_file".into(), text("")),
        ("hub_cadd_file".into(), text("")),
        (
            "hub_file_dir".into(),
            text(display(&result.join("private_hubs"))),
        ),
    ]);
    for key in OPTIONAL_INPUTS {
        if !config.truthy(key)
            || (!matches!(phase, "prioritization" | "prioritization-acmg")
                && *key == "cds_fasta_file")
        {
            continue;
        }
        let source = config
            .str(key)
            .ok_or_else(|| format!("{key}: expected a path string"))?;
        let name = contract::file_name(source);
        let destination = copy_writable(
            Path::new(source),
            &result.join("private_inputs").join(key).join(&name),
        )?;
        updates.push((key.to_string(), text(display(&destination))));
        if source.ends_with(".vcf.gz") || source.ends_with(".vcf.bgz") {
            copy_writable(
                Path::new(&format!("{source}.tbi")),
                Path::new(&format!("{}.tbi", display(&destination))),
            )?;
        }
        if *key == "cds_fasta_file" {
            let original_cache = cds_cache_path(source);
            if exists(&original_cache) {
                copy_writable(&original_cache, &cds_cache_path(&display(&destination)))?;
            }
        }
    }
    if phase.starts_with("annotation") && mode == "raw" {
        for (key, option) in [
            ("hub_vcf_file", "vcf_hub_mode"),
            ("hub_cadd_file", "cadd_hub_mode"),
        ] {
            if field(option) == "off" || !config.truthy(key) {
                continue;
            }
            let source = contract::hub_source(&config, key)?;
            let name = contract::file_name(&display(&source));
            let destination = copy_writable(&source, &result.join("private_hubs").join(name))?;
            updates.push((key.into(), text(display(&destination))));
            if key == "hub_vcf_file" {
                copy_writable(
                    Path::new(&format!("{}.tbi", display(&source))),
                    Path::new(&format!("{}.tbi", display(&destination))),
                )?;
            }
        }
    }
    let path = PathBuf::from(crate::python_path::absolute(Path::new("task.yaml"))?);
    fs::write(&path, render_config(&config, &updates)?)?;
    Ok((path, temp))
}

#[cfg(test)]
pub(crate) fn private_config_for_test(
    bundle: &Path,
    record: &Json,
    phase: &str,
    cpus: i64,
) -> Result<PathBuf> {
    private_config(bundle, record, phase, cpus).map(|(path, _)| path)
}

const STUB_CADD: &str =
    "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\tFeatureID\n1\t10\tA\tG\t1.0\t10\tSTUB\n";

/// One function of annotation_vcf.sh, with the arguments main_workflow
/// gives it (read with the script's own read_yaml, so a missing key is
/// `null` as there). Sourced from this file, the script defines its
/// functions and runs nothing itself.
/// Arguments: SOURCE_ROOT CONFIG VCF STEP [OUTPUT].
/// One step of prioritization_vcf_per_fam.sh, as main_prioritization runs
/// it (threads from threads_per_fam). Sourced from this file, the script
/// defines its functions and runs nothing itself. A family with no variant
/// left gets empty tables, as main_prioritization writes them.
/// Arguments: SOURCE_ROOT CONFIG STEP ARGUMENTS...
const PRIORITIZATION_STEP: &str = r#"#!/usr/bin/env bash
source "$1/scripts/prioritization_vcf_per_fam.sh" || exit 1
config=$2 step=$3
shift 3
threads=$(read_yaml "$config" threads_per_fam)
case $step in
    combine)
        vcf=$1 cadd=$2 table=${1/.vcf*/.tsv}
        if [[ $(count_vcf_records "$vcf") -lt 1 ]] && check_vcf_validity "$vcf" 0; then : > "$table"; exit 0; fi
        [[ $(wc -l < "$cadd") -le 1 ]] && { echo "The CADD table $cadd is empty" >&2; exit 0; }
        prepare_combined_tab "$vcf" "$cadd" "$threads" "$config" && [[ -f $table ]] ;;
    splicing)
        [[ -s $1 ]] || exit 0
        interpret_splicing_annotations "$1" "$config" "$threads" ;;
    utr)
        [[ -s $1 ]] || exit 0
        interpret_utr_annotations "$1" "$config" "$threads" ;;
    acmg)
        table=$1 family=${2:-}
        [[ -s $table ]] || { : > "${table%.tsv}.acmg.tsv"; exit 0; }
        acmg_assign_table -i "$table" -c "$config" -t "$threads" ${family:+-f "$family"} ;;
    *)
        echo "Unknown prioritization step: $step" >&2; exit 2 ;;
esac
"#;

const ANNOTATION_STEP: &str = r#"#!/usr/bin/env bash
source "$1/scripts/annotation_vcf.sh" || exit 1
config=$2 vcf=$3 step=$4 output=${5:-}
value() { read_yaml "$config" "$1"; }
case $step in
    preprocess)
        preprocess_vcf -i "$vcf" -o "$output" -r "$(value ref_genome)" -t "$(value threads)" ;;
    gnomad)
        anno_agg_gnomAD_data "$vcf" "$(value threads)" "$(value assembly)" "$(value gnomad_vcf_chrX)" ;;
    clinvar)
        anno_clinvar_data "$vcf" "$(value clinvar_vcf)" "$(value threads)" ;;
    control)
        optional=$(value control_vcf_optional)
        anno_control_vcf_allele "$vcf" "$(value control_vcf)" "$(value ped_file)" "$(value threads)" "${optional:-false}" ;;
    vep)
        anno_VEP_data --input_vcf "$vcf" --config "$config" --assembly "$(value assembly)"             --ref_genome "$(value ref_genome)" --vep_cache_dir "$(value vep_cache_dir)"             --vep_plugins_dir "$(value vep_plugins_dir)"             --vep_plugins_cachedir "$(value vep_plugins_cachedir)" --threads "$(value threads)" ;;
    cadd)
        Calculate_CADD "$vcf" "$config" "$output" ;;
    *)
        echo "Unknown annotation step: $step" >&2; exit 2 ;;
esac
"#;

pub(crate) fn run(options: &Options) -> Result<()> {
    let bundle = PathBuf::from(options.value("bundle"));
    let record = read_json(&bundle.join("contract.json"))?;
    let field = |key: &str| record.get(key).and_then(Json::as_str).unwrap_or_default();
    let stub = options.flag("stub");
    let phase = options.value("phase");
    let boundary = record
        .get("task_boundary")
        .and_then(Json::as_str)
        .unwrap_or("off");
    require_task_boundary(boundary)?;
    let recorded_stub = match record.get("stub") {
        Some(Json::Bool(flag)) => *flag,
        Some(Json::Null) | None => false,
        Some(other) => !crate::json_text::py_equal(other, &Json::Int(0)),
    };
    if recorded_stub != stub {
        return refuse("Stub and real task contracts cannot be mixed");
    }
    let cpus = options.int("cpus")?;
    if cpus < 1 {
        return refuse("cpus must be positive");
    }
    let root = PathBuf::from(field("legacy_root"));
    if source_digest(&root)? != field("source_digest") {
        return refuse("Scientific or adapter source changed after validation; task refused");
    }
    let mut verified = None;
    if !stub {
        let manifest = read_json(&bundle.join("resources.json"))?;
        let paths = verify_manifest_for(
            &manifest,
            &bundle.join("legacy.yaml"),
            &root,
            field("mode"),
            field("vcf_hub_mode"),
            field("cadd_hub_mode"),
            phase,
        )?;
        task_thread_settings(
            &read_config(&bundle.join("legacy.yaml"))?,
            phase,
            field("mode"),
            cpus,
        )?;
        verified = Some((manifest, paths));
    }
    let base = field("base").to_owned();
    // The contract lists every family the run covers; each family task
    // names one. A contract without families is the no-pedigree run.
    let families: Vec<String> = match record.get("families") {
        Some(Json::Array(items)) => items
            .iter()
            .filter_map(|item| item.as_str().map(str::to_owned))
            .collect(),
        _ => Vec::new(),
    };
    let family = options
        .get("family")
        .filter(|v| !v.is_empty())
        .map(str::to_owned);
    // One function of the annotation script per task; only raw mode runs them.
    let step = phase.strip_prefix("annotation-");
    if step.is_some() && field("mode") != "raw" {
        return refuse("Annotation steps run only in raw mode");
    }
    match (&family, phase) {
        (Some(_), annotation) if annotation.starts_with("annotation") => {
            return refuse("Annotation covers the whole cohort; it takes no --family");
        }
        (Some(id), _) if !families.contains(id) => {
            return refuse(format!(
                "Family {id} is not in this run's validated families"
            ));
        }
        (None, family_phase)
            if (family_phase == "filtration" || family_phase.starts_with("prioritization"))
                && !families.is_empty() =>
        {
            return refuse("This run has PED families; name one with --family");
        }
        _ => {}
    }
    let family_base = match &family {
        Some(id) => format!("{base}.anno.{id}.filtered"),
        None => format!("{base}.anno.filtered"),
    };
    let result = PathBuf::from(crate::python_path::absolute(Path::new("result"))?);
    fs::create_dir(&result).map_err(|e| crate::python_path::os_error(e, &result, None))?;
    let name = match phase {
        "annotation" if field("mode") == "annotated" => format!("{base}.anno"),
        "annotation" | "annotation-preprocess" => base.clone(),
        // The other steps work on (a copy of) the preprocessed VCF.
        _ if step.is_some() => format!("{base}.anno"),
        "filtration" => format!("{base}.anno"),
        _ => family_base.clone(),
    };
    // The prioritization steps after combine read and rewrite the family
    // table instead of a VCF; their input is the table itself.
    let table_phase = matches!(
        phase,
        "prioritization-splicing"
            | "prioritization-utr"
            | "prioritization-join"
            | "prioritization-acmg"
    );
    let input_vcf = if table_phase {
        copy_writable(
            Path::new(options.value("table")),
            &result.join(format!("{family_base}.tsv")),
        )?
    } else {
        copy_writable(
            Path::new(options.value("vcf")),
            &result.join(format!("{name}.vcf.gz")),
        )?
    };
    let input_index = PathBuf::from(format!("{}.tbi", display(&input_vcf)));
    if !table_phase {
        copy_writable(Path::new(options.value("vcf-index")), &input_index)?;
    }
    let cadd = options.get("cadd").filter(|v| !v.is_empty());
    let cadd_path = match cadd {
        Some(cadd) => {
            let cadd_base = if phase.starts_with("prioritization") {
                family_base.clone()
            } else {
                format!("{base}.anno")
            };
            Some(copy_writable(
                Path::new(cadd),
                &result.join(format!("{cadd_base}.cadd.tsv")),
            )?)
        }
        None if matches!(
            phase,
            "filtration" | "prioritization" | "prioritization-combine" | "annotation-merge"
        ) =>
        {
            return refuse("Filtration/prioritization require the paired CADD TSV");
        }
        None => None,
    };
    if stub {
        // Tiny syntactic fixtures; no scientific acceptance is claimed.
        match phase {
            "annotation" | "annotation-preprocess" => {
                let destination = result.join(format!("{base}.anno.vcf.gz"));
                if input_vcf != destination {
                    copy_writable(&input_vcf, &destination)?;
                    copy_writable(
                        &input_index,
                        Path::new(&format!("{}.tbi", display(&destination))),
                    )?;
                }
                if phase == "annotation" {
                    fs::write(result.join(format!("{base}.anno.cadd.tsv")), STUB_CADD)?;
                }
            }
            "annotation-cadd" => {
                fs::write(result.join(format!("{base}.anno.cadd.tsv")), STUB_CADD)?
            }
            // The copied input already is the step's output.
            _ if step.is_some() => {}
            "prioritization-combine" => fs::write(
                result.join(format!("{family_base}.tsv")),
                "Chr\tStart\tRef\tAlt\tFeature\n1\t10\tA\tG\tSTUB\n",
            )?,
            "prioritization-splicing" | "prioritization-utr" | "prioritization-join" => {}
            "filtration" => {
                let destination = result.join(format!("{family_base}.vcf.gz"));
                copy_writable(&input_vcf, &destination)?;
                copy_writable(
                    &input_index,
                    Path::new(&format!("{}.tbi", display(&destination))),
                )?;
                copy_writable(
                    cadd_path.as_deref().unwrap(),
                    &result.join(format!("{family_base}.cadd.tsv")),
                )?;
            }
            _ => {
                fs::write(
                    result.join(format!("{family_base}.tsv")),
                    "Chr\tStart\tRef\tAlt\tFeature\tsort_index\tvariant_rank\tACMG_class\tACMG_criteria\n1\t10\tA\tG\tSTUB\t1\t1\tSTUB_ONLY\tSTUB_ONLY\n",
                )?;
                fs::write(
                    result.join(format!("{family_base}.acmg.tsv")),
                    "variant_id\tSTUB_ONLY\n1:10:A:G\t1\n",
                )?;
            }
        }
    } else {
        let (config_path, temp) = private_config(&bundle, &record, phase, cpus)?;
        let paths = &verified.as_ref().unwrap().1;
        let binary = paths
            .iter()
            .find(|(k, _)| k == "priva_binary")
            .map(|(_, p)| display(p))
            .unwrap_or_default();
        let script = |name: &str| display(&root.join("scripts").join(name));
        let input = display(&input_vcf);
        let config_text = display(&config_path);
        let step_command = |step: &str, output: Option<String>| -> Result<Vec<String>> {
            fs::write("annotation_step.sh", ANNOTATION_STEP)?;
            let mut command = vec![
                "bash".into(),
                "annotation_step.sh".into(),
                display(&root),
                config_text.clone(),
                input.clone(),
                step.into(),
            ];
            command.extend(output);
            Ok(command)
        };
        let anno = |suffix: &str| display(&result.join(format!("{base}.anno{suffix}")));
        let prioritization_step = |step: &str, arguments: Vec<String>| -> Result<Vec<String>> {
            fs::write("prioritization_step.sh", PRIORITIZATION_STEP)?;
            let mut command = vec![
                "bash".into(),
                "prioritization_step.sh".into(),
                display(&root),
                config_text.clone(),
                step.into(),
            ];
            command.extend(arguments);
            Ok(command)
        };
        let command: Option<Vec<String>> = match phase {
            "annotation-preprocess" => Some(step_command("preprocess", Some(anno(".vcf.gz")))?),
            "annotation-cadd" => Some(step_command("cadd", Some(anno(".cadd.tsv")))?),
            "annotation-merge" => {
                let mut steps = Vec::new();
                for (name, option) in [
                    ("gnomAD", "gnomad-vcf"),
                    ("ClinVar", "clinvar-vcf"),
                    ("control cohort", "control-vcf"),
                    ("VEP", "vep-vcf"),
                ] {
                    let Some(path) = options.get(option).filter(|v| !v.is_empty()) else {
                        return refuse(format!("annotation-merge needs --{option}"));
                    };
                    steps.push((name, PathBuf::from(path)));
                }
                let steps: Vec<(&str, &Path)> =
                    steps.iter().map(|(n, p)| (*n, p.as_path())).collect();
                super::info_merge::merge(&input_vcf, &steps, &input_vcf)?;
                None
            }
            _ if step.is_some() => Some(step_command(step.unwrap(), None)?),
            "prioritization-combine" => Some(prioritization_step(
                "combine",
                vec![input.clone(), display(cadd_path.as_deref().unwrap())],
            )?),
            "prioritization-splicing" => {
                Some(prioritization_step("splicing", vec![input.clone()])?)
            }
            "prioritization-utr" => Some(prioritization_step("utr", vec![input.clone()])?),
            "prioritization-join" => {
                let (Some(splicing), Some(utr)) = (
                    options.get("splicing-table").filter(|v| !v.is_empty()),
                    options.get("utr-table").filter(|v| !v.is_empty()),
                ) else {
                    return refuse("prioritization-join needs --splicing-table and --utr-table");
                };
                super::table_join::join(
                    &input_vcf,
                    Path::new(splicing),
                    Path::new(utr),
                    &input_vcf,
                )?;
                None
            }
            "prioritization-acmg" => {
                let mut arguments = vec![input.clone()];
                arguments.extend(family.clone());
                Some(prioritization_step("acmg", arguments)?)
            }
            "annotation" if field("mode") == "raw" => Some(vec![
                "bash".into(),
                script("annotation_vcf.sh"),
                "main_workflow".into(),
                "-i".into(),
                input.clone(),
                "--config".into(),
                config_text.clone(),
            ]),
            // Without a family the scripts run their no-pedigree mode.
            "filtration" => {
                let mut command: Vec<String> = vec![
                    "bash".into(),
                    script("filtration_vcf_per_fam.sh"),
                    "-v".into(),
                    input.clone(),
                ];
                if let Some(id) = &family {
                    command.extend(["-f".into(), id.clone()]);
                }
                command.extend(["-c".into(), config_text.clone()]);
                Some(command)
            }
            "prioritization" => {
                let mut command: Vec<String> = vec![
                    "bash".into(),
                    script("prioritization_vcf_per_fam.sh"),
                    "-i".into(),
                    input.clone(),
                ];
                if let Some(id) = &family {
                    command.extend(["-f".into(), id.clone()]);
                }
                command.extend([
                    "-c".into(),
                    config_text.clone(),
                    "-t".into(),
                    display(cadd_path.as_deref().unwrap()),
                ]);
                Some(command)
            }
            _ => None, // An explicit preannotated import keeps the input bytes.
        };
        if let Some(command) = command {
            let threads = cpus.to_string();
            let temp_text = display(&temp);
            // The pinned program, whatever PRIVA_BINARY the caller exported.
            let status = Command::new(&command[0])
                .args(&command[1..])
                .env("PRIVA_BINARY", &binary)
                .env("TMPDIR", &temp_text)
                .env("TMP", &temp_text)
                .env("TEMP", &temp_text)
                .env("PYTHONDONTWRITEBYTECODE", "1")
                .env("OMP_NUM_THREADS", &threads)
                .env("OPENBLAS_NUM_THREADS", &threads)
                .env("MKL_NUM_THREADS", &threads)
                .env("NUMEXPR_NUM_THREADS", &threads)
                .status()?;
            if !status.success() {
                let parts: Vec<&str> = command.iter().map(String::as_str).collect();
                return Err(called_process_error(&parts, status));
            }
        }
    }
    let suffixes: &[&str] = match phase {
        "annotation" | "annotation-merge" => {
            &[".anno.vcf.gz", ".anno.vcf.gz.tbi", ".anno.cadd.tsv"]
        }
        "annotation-cadd" => &[".anno.cadd.tsv"],
        _ if step.is_some() => &[".anno.vcf.gz", ".anno.vcf.gz.tbi"],
        "filtration" => &[".vcf.gz", ".vcf.gz.tbi", ".cadd.tsv"],
        "prioritization-combine"
        | "prioritization-splicing"
        | "prioritization-utr"
        | "prioritization-join" => &[".tsv"],
        _ => &[".tsv", ".acmg.tsv"],
    };
    let stem = if phase.starts_with("annotation") {
        &base
    } else {
        &family_base
    };
    let outputs: Vec<PathBuf> = suffixes
        .iter()
        .map(|suffix| result.join(format!("{stem}{suffix}")))
        .collect();
    for path in &outputs {
        if !path.is_file() {
            return refuse(format!(
                "Legacy task returned without required output: {}",
                path.display()
            ));
        }
    }
    if let Some((manifest, _)) = &verified {
        verify_manifest_for(
            manifest,
            &bundle.join("legacy.yaml"),
            &root,
            field("mode"),
            field("vcf_hub_mode"),
            field("cadd_hub_mode"),
            phase,
        )?;
    }
    let mut provenance = record.clone();
    provenance.set("phase", text(phase));
    provenance.set("family", family.clone().map_or(Json::Null, text));
    provenance.set("cpus", Json::Int(cpus as i128));
    let mut listed = Vec::new();
    for path in &outputs {
        listed.push(object(vec![
            (
                "name",
                text(path.file_name().unwrap().to_string_lossy().into_owned()),
            ),
            ("sha256", text(sha256(path)?)),
            ("size", Json::Int(fs::metadata(path)?.len() as i128)),
        ]));
    }
    provenance.set("outputs", Json::Array(listed));
    provenance.set(
        "cache_isolation",
        text("shared resources read; mutable VCF/CADD hubs, optional cohort VCFs and CDS lazy cache copied into task"),
    );
    // One provenance file per family task, so the published names differ.
    let name = match &family {
        Some(id) => format!("{phase}.{id}.provenance.json"),
        None => format!("{phase}.provenance.json"),
    };
    write_json(Path::new(&name), &provenance)
}
