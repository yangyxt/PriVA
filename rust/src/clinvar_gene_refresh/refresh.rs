//! The refresh plan, completed-receipt reuse and the locked installation.
use super::deploy::{build_lock, run_original};
use super::profile::{glob, profile_signature, resource_policy, scientific_profile};
use super::remote::{acquire_inputs, latest_metadata, local_file_date};
use super::*;
use crate::config_text::{self, Writer};
use crate::installer_prerequisites::paths::{exists, is_dir, is_file};
use crate::python_path::absolute;

fn os_failure(category: &str) -> bool {
    matches!(
        category,
        "OSError"
            | "FileNotFoundError"
            | "PermissionError"
            | "IsADirectoryError"
            | "NotADirectoryError"
    )
}
fn same(left: &Json, right: &Json) -> bool {
    json_text::py_equal(left, right)
}
fn receipt_path(output: &str) -> String {
    format!("{output}.ready.json")
}
/// `render_config(text, updates)` of installer_config_paths.py.
fn render_config(text: &str, updates: &[(String, Json)]) -> Result<String> {
    config_text::render(text, updates, Writer::Paths)
        .map_err(|failure| error(failure.kind, failure.message))
}
fn updates_of(value: &Json) -> Result<Vec<(String, Json)>> {
    match value {
        Json::Object(items) => Ok(items.clone()),
        other => Err(error(
            "AttributeError",
            format!("'{}' object has no attribute 'items'", type_name(other)),
        )),
    }
}

/// `clinvar_stats_path(config, source_root)` of installer_prerequisites.py.
fn clinvar_stats_path(config: &Json, source_root: &str) -> Result<String> {
    let configured = get_or_null(config, "clinvar_gene_stat")?;
    let selected = if truth(&configured) {
        py_str(&configured)
    } else {
        "clinvar_2star_stats.pkl.gz".to_owned()
    };
    let path = pathlib(&selected);
    Ok(if path.starts_with('/') || selected.contains('/') {
        path
    } else {
        join(&join(source_root, "data/ClinVar"), &path)
    })
}

/// `baseline_metadata(path)`: existence, size and mtime only.
fn baseline_metadata(path: &str) -> Result<Json> {
    let present = exists(path)?;
    let mut value = object([("path", string(path)), ("exists", Json::Bool(present))]);
    if !present {
        return Ok(value);
    }
    require(
        is_file(path)?,
        format!("Bundled ClinVar gene statistics is not a file: {path}"),
    )?;
    let info = stat(path)?;
    let (seconds, micros) = float_timestamp(info.mtime(), info.mtime_nsec());
    value.set("size", Json::Int(info.len() as i128));
    value.set(
        "mtime_ns",
        Json::Int(nanoseconds(info.mtime(), info.mtime_nsec())),
    );
    value.set("mtime_utc", string(utc_text(seconds, micros)));
    Ok(value)
}

/// `owned_receipt(output, assembly=..., vep_major=...)`
fn owned_receipt(output: &str, assembly: &str, major: &str) -> Result<Option<Json>> {
    let path = receipt_path(output);
    if !exists(&path)? {
        return Ok(None);
    }
    let value = read_json(&path)?;
    if !same(&get_or_null(&value, "kind")?, &string(RECEIPT_KIND)) {
        return Ok(None);
    }
    let matches = same(&get_or_null(&value, "status")?, &string("complete"))
        && same(&get_or_null(&value, "assembly")?, &string(assembly))
        && same(&get_or_null(&value, "vep_major")?, &string(major))
        && absolute(Path::new(&path_of(item(item(&value, "output")?, "path")?)?))?
            == absolute(Path::new(output))?;
    require(
        matches,
        format!("Completed ClinVar refresh receipt does not match the selected profile: {path}"),
    )?;
    require(
        same(&identity(output)?, item(&value, "output")?),
        format!("Generated ClinVar statistics differ from their completed receipt: {output}"),
    )?;
    Ok(Some(value))
}

fn suffix(path: &str) -> &str {
    let name = name_of(path);
    match name.rfind('.') {
        Some(at) if at > 0 && at < name.len() - 1 => &name[at..],
        _ => "",
    }
}

/// `selected_receipt(selected, root, ...)`: the owned refresh receipt for the
/// selected statistics, following a cache-install copy to its source.
fn selected_receipt(
    selected: &str,
    root: &str,
    assembly: &str,
    major: &str,
) -> Result<Option<Json>> {
    let known = owned_receipt(selected, assembly, major)?;
    let record_path = format!("{selected}.source.json");
    if known.is_some() || suffix(selected) != ".rkyv" || !is_file(&record_path)? {
        return Ok(known);
    }
    let record = read_json(&record_path)?;
    let resolve = crate::installer_prerequisites::paths::resolve;
    let archive = lookup(&record, "archive")?.map_or_else(String::new, py_str);
    if !same(&get_or_null(&record, "key")?, &string("clinvar_gene_stat"))
        || resolve(&pathlib(&archive))? != resolve(selected)?
    {
        return Ok(None);
    }
    let origin = get_or_null(&record, "origin")?;
    if matches!(origin, Json::Null) {
        return owned_receipt(&path_of(item(&record, "source")?)?, assembly, major);
    }
    if !same(&origin, &string("direct")) {
        return Ok(None);
    }
    let wanted = get_or_null(&record, "archive_sha256")?;
    let mut matches: Vec<Json> = Vec::new();
    let base = join(root, &format!("gene_statistics/{assembly}"));
    for first in glob(&base, |_| true) {
        if !is_dir(&first)? {
            continue;
        }
        for second in glob(&first, |_| true) {
            let path = join(&second, "clinvar_2star_stats.rkyv.ready.json");
            if !is_dir(&second)? || !exists(&path)? {
                continue;
            }
            let output = path.strip_suffix(".ready.json").unwrap_or(&path);
            let value = match owned_receipt(output, assembly, major) {
                Ok(value) => value,
                // Another profile's or a changed completion is never adopted.
                Err(failure)
                    if matches!(exception(failure.as_ref()).0, "ValueError" | "KeyError")
                        || os_failure(exception(failure.as_ref()).0) =>
                {
                    continue;
                }
                Err(failure) => return Err(failure),
            };
            if let Some(value) = value
                && same(item(item(&value, "output")?, "sha256")?, &wanted)
            {
                matches.push(value);
            }
        }
    }
    let key = |value: &Json| -> (String, String) {
        (
            value.get("data_date").map_or_else(String::new, text_of),
            value.get("completed_utc").map_or_else(String::new, text_of),
        )
    };
    let mut best: Option<Json> = None;
    for value in matches {
        if best.as_ref().is_none_or(|known| key(&value) > key(known)) {
            best = Some(value);
        }
    }
    Ok(best)
}

/// `validate_reuse(prior, plan)`: a ready statistic alone cannot authorize
/// stale annotation/input reuse.
fn validate_reuse(prior: &Json, plan: &Json) -> Result<()> {
    require(
        profile_signature(lookup(prior, "scientific_profile")?)?
            == profile_signature(Some(item(plan, "scientific_profile")?))?,
        "ClinVar refresh scientific configuration/input/source identity changed; \
         the previous statistic was not selected for this different profile",
    )?;
    require(
        same(
            &get_or_null(item(prior, "latest")?, "publisher_md5")?,
            &get_or_null(item(plan, "latest")?, "publisher_md5")?,
        ),
        "Official ClinVar bytes changed within the recorded fileDate; no stale result was selected",
    )?;
    let artifacts = get_or_null(prior, "completed_artifacts")?;
    let artifacts = match artifacts {
        Json::Null => Json::Object(Vec::new()),
        other => other,
    };
    for (label, recorded) in updates_of(&artifacts)? {
        let path = item(&recorded, "path")?;
        let current = match identity(&path_of(path)?) {
            Ok(value) => value,
            Err(failure) if os_failure(exception(failure.as_ref()).0) => {
                return Err(value_error(format!(
                    "Completed ClinVar artifact is missing: {label}: {}",
                    py_str(path)
                )));
            }
            Err(failure) => return Err(failure),
        };
        require(
            same(&current, &recorded),
            format!(
                "Completed ClinVar artifact changed: {label}: {}",
                py_str(path)
            ),
        )?;
    }
    let mut labels: Vec<String> = updates_of(&artifacts)?
        .into_iter()
        .map(|(k, _)| k)
        .collect();
    labels.sort();
    require(
        labels == ["aa", "splice", "vcf", "vcf_index"],
        "ClinVar ready receipt does not bind all completed annotation/AA/splice artifacts",
    )?;
    for (_, recorded) in updates_of(item(item(prior, "inputs")?, "sources")?)? {
        let current = identity(&path_of(item(&recorded, "path")?)?)?;
        for key in ["path", "size", "sha256"] {
            require(
                same(item(&current, key)?, item(&recorded, key)?),
                "A pinned ClinVar/MANE regeneration input changed; no stale result was selected",
            )?;
        }
    }
    if truth(item(plan, "explicit_vcf")?) {
        let current = identity(&path_of(item(plan, "input_vcf")?)?)?;
        require(
            same(&current, item(item(prior, "completed_artifacts")?, "vcf")?),
            "Explicit ClinVar VCF identity differs from the completed refresh, even though fileDate \
             matches; the explicit selection was not replaced",
        )?;
    }
    Ok(())
}

/// `plan_refresh(config_path, source_root)`: bundled metadata only; NCBI is
/// probed for a default or owned selection.
pub(super) fn plan_refresh(config_path: &str, source_root: &str) -> Result<Json> {
    let config = load_config(config_path)?;
    let Some(assembly) = assembly_of(&get_or_null(&config, "assembly")?)? else {
        return Err(value_error(
            "ClinVar refresh requires hg19/GRCh37 or hg38/GRCh38",
        ));
    };
    let ucsc = ucsc(assembly);
    let major = vep_major(&config);
    require(
        positive(&major),
        "ClinVar refresh requires the configured numeric VEP major",
    )?;
    let baseline = join(source_root, "data/ClinVar/clinvar_2star_stats.pkl.gz");
    let selected = clinvar_stats_path(&config, source_root)?;
    let is_baseline = absolute(Path::new(&selected))? == absolute(Path::new(&baseline))?;
    if !is_baseline {
        nonempty(&selected, "selected ClinVar gene statistics")?;
    }
    let root_value = get_or_null(&config, "clinvar_vcf_dir")?;
    let known = if is_baseline {
        None
    } else if !empty(&root_value) {
        selected_receipt(&selected, &path_of(&root_value)?, assembly, &major)?
    } else {
        owned_receipt(&selected, assembly, &major)?
    };
    let mut plan = object([
        ("kind", string(RECEIPT_KIND)),
        ("config_path", string(config_path)),
        ("source_root", string(source_root)),
        ("config_sha256", string(sha256(config_path)?)),
        ("assembly", string(assembly)),
        ("ucsc_assembly", string(ucsc)),
        ("vep_major", string(major.clone())),
        ("selected_stat", string(selected.clone())),
        ("action", string("preserve")),
        (
            "reason",
            string("explicitly selected statistics are retained"),
        ),
    ]);
    if !is_baseline && known.is_none() {
        return Ok(plan);
    }
    let policy = resource_policy(source_root)?;
    let mut baseline_state = baseline_metadata(&baseline)?;
    let missing = !baseline_state.get("size").is_some_and(truth);
    let months = match item(&policy, "threshold_months")? {
        Json::Int(months) => *months,
        _ => 0,
    };
    if !truth(item(&policy, "automatic")?) && !missing {
        plan.set("resource_policy", policy);
        plan.set(
            "reason",
            string(
                "automatic ClinVar gene-statistics refresh is disabled in the resource registry",
            ),
        );
        return Ok(plan);
    }
    let boundary = if missing {
        None
    } else {
        let stamp = text_of(item(&baseline_state, "mtime_utc")?);
        Some(add_calendar_months(parse_date(&stamp[..10])?, months)?)
    };
    if let Some(boundary) = boundary {
        baseline_state.set("refresh_boundary_utc_date", string(date_text(boundary)));
    }
    plan.set("resource_policy", policy.clone());
    plan.set("baseline", baseline_state);
    plan.set("baseline_missing", Json::Bool(missing));
    // The install workflow's download step (D0) names the release it acquired;
    // the plan then judges that release instead of asking NCBI again.
    let release_inputs = get_or_null(&config, "clinvar_release_inputs")?;
    let latest = if empty(&release_inputs) {
        latest_metadata(assembly, &text_of(item(&policy, "latest_url_template")?))?
    } else {
        let release = read_json(&path_of(&release_inputs)?)?;
        require(
            same(&get_or_null(&release, "assembly")?, &string(assembly)),
            "The ClinVar download step fetched another assembly",
        )?;
        plan.set("release_inputs", release_inputs.clone());
        item(&release, "latest")?.clone()
    };
    let file_date = text_of(item(&latest, "file_date")?);
    plan.set("latest", latest.clone());
    if let Some(boundary) = boundary
        && parse_date(&file_date)? < boundary
    {
        plan.set(
            "reason",
            string(format!(
                "latest ClinVar data date is before the {months}-calendar-month boundary"
            )),
        );
        return Ok(plan);
    }
    let root = path_of(item(&config, "clinvar_vcf_dir")?)?;
    let profile = scientific_profile(&config, source_root, &policy)?;
    let signature = profile_signature(Some(&profile))?;
    let release = join(
        &root,
        &format!("gene_statistics/{assembly}/{}", file_date.replace('-', "")),
    );
    let mut build = join(&release, &format!("vep{major}-{}", &signature[..16]));
    // priva clinvar-gene-stats writes the archive; cache-install publish adopts it.
    let mut output = join(&build, "clinvar_2star_stats.rkyv");
    let mut prior = owned_receipt(&output, assembly, &major)?;
    // Reuse an existing unsuffixed completion only for its exact profile.
    // Unsuffixed workspaces predate the native builder and hold a Python pickle.
    let legacy_build = join(&release, &format!("vep{major}"));
    let legacy_output = join(&legacy_build, "clinvar_2star_stats.pkl.gz");
    if prior.is_none() && is_file(&receipt_path(&legacy_output))? {
        let legacy = owned_receipt(&legacy_output, assembly, &major)?;
        if let Some(legacy) = legacy
            && profile_signature(lookup(&legacy, "scientific_profile")?)? == signature
        {
            (build, output, prior) = (legacy_build, legacy_output, Some(legacy));
        }
    }
    let input_vcf = selected_vcf(&config, &root, ucsc)?;
    let default_vcf = join(&root, &format!("clinvar.{ucsc}.vep.vcf.gz"));
    let known = known.unwrap_or_else(|| Json::Object(Vec::new()));
    let prior_vcf = match lookup(&known, "config_updates")? {
        Some(updates) => lookup(updates, "clinvar_vcf")?.cloned(),
        None => None,
    }
    .filter(|value| !matches!(value, Json::Null));
    let input_absolute = absolute(Path::new(&input_vcf))?;
    let managed = input_absolute == absolute(Path::new(&default_vcf))?
        || match &prior_vcf {
            Some(value) => {
                input_absolute == absolute(Path::new(&path_of(value)?))?
                    && same(&get_or_null(&known, "managed_root")?, &string(root.clone()))
            }
            None => false,
        };
    if !managed && date_text(local_file_date(&input_vcf)?) != file_date {
        return Err(value_error(format!(
            "Explicit ClinVar VCF {input_vcf} has fileDate {}; latest official ClinVar is {file_date}. \
             No replacement was selected; do not combine a newer summary with the explicitly \
             selected older/different release.",
            date_text(local_file_date(&input_vcf)?)
        )));
    }
    let reuse = prior.as_ref().is_some_and(|prior| {
        same(
            &get_or_null(prior, "data_date").unwrap_or(Json::Null),
            &string(file_date.clone()),
        )
    });
    plan.set("action", string(if reuse { "reuse" } else { "rebuild" }));
    plan.set(
        "reason",
        string(if missing {
            "missing/empty bundled ClinVar gene statistics require managed reconstruction"
                .to_owned()
        } else {
            format!("latest ClinVar data date reaches the {months}-calendar-month boundary")
        }),
    );
    plan.set("managed_root", string(root));
    plan.set("build_dir", string(build));
    plan.set("output_path", string(output));
    plan.set("input_vcf", string(input_vcf.clone()));
    plan.set("explicit_vcf", Json::Bool(!managed));
    plan.set(
        "explicit_vcf_identity",
        if managed {
            Json::Null
        } else {
            identity(&input_vcf)?
        },
    );
    plan.set("scientific_profile", profile);
    plan.set("scientific_profile_sha256", string(signature));
    if reuse && let Some(prior) = &prior {
        validate_reuse(prior, &plan)?;
    }
    Ok(plan)
}

/// `require_workspace_owner(plan)`: an interrupted annotation is reusable
/// only within its recorded profile.
fn require_workspace_owner(plan: &Json) -> Result<()> {
    let folder = path_of(item(plan, "build_dir")?)?;
    let previous_path = join(&folder, "refresh.plan.json");
    if exists(&previous_path)? {
        let previous = read_json(&previous_path)?;
        let latest = get_or_null(&previous, "latest")?;
        let latest = if matches!(latest, Json::Null) && previous.get("latest").is_none() {
            Json::Object(Vec::new())
        } else {
            latest
        };
        let plan_latest = item(plan, "latest")?;
        let owned = same(&get_or_null(&previous, "kind")?, &string(RECEIPT_KIND))
            && profile_signature(lookup(&previous, "scientific_profile")?)?
                == profile_signature(Some(item(plan, "scientific_profile")?))?
            && same(
                &get_or_null(&previous, "explicit_vcf_identity")?,
                item(plan, "explicit_vcf_identity")?,
            )
            && same(
                &get_or_null(&latest, "file_date")?,
                item(plan_latest, "file_date")?,
            )
            && same(
                &get_or_null(&latest, "publisher_md5")?,
                item(plan_latest, "publisher_md5")?,
            );
        require(
            owned,
            format!(
                "Occupied incomplete ClinVar workspace belongs to a different or unknown scientific \
                 profile/input; its annotation and previous plan were retained without reuse: {folder}"
            ),
        )?;
    } else {
        let occupied = exists(&join(&folder, "install.config.yaml"))?
            || exists(&path_of(item(plan, "output_path")?)?)?
            || !glob(&folder, |name| {
                name.starts_with("clinvar.") && name[8..].contains(".vep")
            })
            .is_empty();
        require(
            !occupied,
            format!(
                "Occupied ClinVar annotation workspace has no owned scientific-profile plan; \
                 it was retained without reuse: {folder}"
            ),
        )?;
    }
    Ok(())
}

/// `complete_refresh(plan, working_config, inputs)`: the receipt, then the
/// original configuration's new selection.
fn complete_refresh(plan: &Json, working: &str, inputs: &Json) -> Result<()> {
    let original = path_of(item(plan, "config_path")?)?;
    require(
        same(&string(sha256(&original)?), item(plan, "config_sha256")?),
        "The original configuration changed during ClinVar regeneration; no new selection was written",
    )?;
    let recorded = item(plan, "baseline")?;
    let mut before = Json::Object(Vec::new());
    for (key, value) in updates_of(recorded)? {
        if key != "refresh_boundary_utc_date" {
            before.set(&key, value);
        }
    }
    require(
        same(
            &baseline_metadata(&path_of(item(recorded, "path")?)?)?,
            &before,
        ),
        "Bundled ClinVar metadata changed during regeneration; no new selection was written",
    )?;
    let current = load_config(working)?;
    let source_root = path_of(item(plan, "source_root")?)?;
    let policy = resource_policy(&source_root)?;
    let now = scientific_profile(&load_config(&original)?, &source_root, &policy)?;
    require(
        profile_signature(Some(&now))?
            == profile_signature(Some(item(plan, "scientific_profile")?))?,
        "ClinVar scientific inputs/source changed during regeneration; no new selection was written",
    )?;
    let ucsc = text_of(item(plan, "ucsc_assembly")?);
    let vcf = selected_vcf(
        &current,
        &path_of(item(&current, "clinvar_vcf_dir")?)?,
        &ucsc,
    )?;
    let latest = item(plan, "latest")?;
    require(
        same(
            &string(date_text(local_file_date(&vcf)?)),
            item(latest, "file_date")?,
        ),
        "Completed annotation fileDate differs from the detected latest ClinVar release",
    )?;
    let mut updates = Vec::new();
    for key in ["clinvar_vcf", "clinvar_aa_stat", "clinvar_splice_stat"] {
        updates.push((key.to_owned(), string(py_str(item(&current, key)?))));
    }
    let output = path_of(item(plan, "output_path")?)?;
    updates.push((
        "clinvar_gene_stat".to_owned(),
        item(plan, "output_path")?.clone(),
    ));
    let mut updates_json = Json::Object(Vec::new());
    for (key, value) in &updates {
        updates_json.set(key, value.clone());
    }
    let completed = object([
        ("vcf", identity(&vcf)?),
        ("vcf_index", identity(&format!("{vcf}.tbi"))?),
        (
            "aa",
            identity(&path_of(item(&current, "clinvar_aa_stat")?)?)?,
        ),
        (
            "splice",
            identity(&path_of(item(&current, "clinvar_splice_stat")?)?)?,
        ),
    ]);
    let mut sources = Json::Object(Vec::new());
    for name in [
        "clinvar_stat_variants.py",
        "combine_annotations.py",
        "annotation_vcf.sh",
    ] {
        sources.set(
            name,
            string(sha256(&join(&source_root, &format!("scripts/{name}")))?),
        );
    }
    let value = object([
        ("kind", string(RECEIPT_KIND)),
        ("status", string("complete")),
        ("data_date", item(latest, "file_date")?.clone()),
        ("scientific_qualified", Json::Bool(false)),
        ("assembly", item(plan, "assembly")?.clone()),
        ("vep_major", item(plan, "vep_major")?.clone()),
        ("baseline", recorded.clone()),
        ("managed_root", item(plan, "managed_root")?.clone()),
        ("latest", latest.clone()),
        ("inputs", inputs.clone()),
        ("resource_policy", item(plan, "resource_policy")?.clone()),
        (
            "scientific_profile",
            item(plan, "scientific_profile")?.clone(),
        ),
        (
            "scientific_profile_sha256",
            item(plan, "scientific_profile_sha256")?.clone(),
        ),
        ("completed_artifacts", completed),
        ("output", identity(&output)?),
        ("annotation", identity(&vcf)?),
        ("completed_utc", string(utc_now())),
        ("config_updates", updates_json),
        ("source_root", item(plan, "source_root")?.clone()),
        ("scientific_sources", sources),
        (
            "qualification",
            string("software refresh; public/family acceptance remains separate"),
        ),
    ]);
    write_json(&receipt_path(&output), &value)?;
    atomic_text(&original, &render_config(&read_text(&original)?, &updates)?)
}

/// `install_refresh(config_path, source_root)`
/// The raw VCF a ClinVar download step (install step D0) acquired, when the
/// configuration names its `clinvar_release_inputs`.
fn downloaded_vcf(config_path: &str) -> Result<Option<String>> {
    let value = get_or_null(&load_config(config_path)?, "clinvar_release_inputs")?;
    if empty(&value) {
        return Ok(None);
    }
    let release = read_json(&path_of(&value)?)?;
    Ok(Some(path_of(item(item(item(&release, "sources")?, "vcf")?, "path")?)?))
}

pub(super) fn install_refresh(config_path: &str, source_root: &str) -> Result<Json> {
    // Preflight selected inputs before metadata access or any acquisition.
    let prerequisites =
        crate::installer_prerequisites::Config::load(&read_text(config_path)?, true)?;
    crate::installer_prerequisites::checks::validate_main(
        &prerequisites,
        source_root,
        false,
        true,
    )?;
    let mut plan = plan_refresh(config_path, source_root)?;
    let action = text_of(item(&plan, "action")?);
    println!(
        "ClinVar gene-statistics refresh: {action}; {}",
        text_of(item(&plan, "reason")?)
    );
    io::stdout().flush().ok();
    if action == "preserve" {
        // A downloaded release (install step D0) is deployed instead of the
        // newest one.
        let downloaded = downloaded_vcf(config_path)?;
        let mut deploy = vec![config_path];
        if let Some(source) = &downloaded {
            deploy.extend(["--source-vcf", source.as_str()]);
        }
        run_original(source_root, "ClinVar_VCF_deploy", &deploy)?;
        for function in ["ClinVar_AA_stat", "ClinVar_Gene_stat"] {
            run_original(source_root, function, &[config_path])?;
        }
        return Ok(plan);
    }
    let folder = path_of(item(&plan, "build_dir")?)?;
    let _lock = build_lock(&folder)?;
    let latest_date = item(item(&plan, "latest")?, "file_date")?.clone();
    let prior = owned_receipt(
        &path_of(item(&plan, "output_path")?)?,
        &text_of(item(&plan, "assembly")?),
        &text_of(item(&plan, "vep_major")?),
    )?;
    if let Some(prior) = prior
        && same(&get_or_null(&prior, "data_date")?, &latest_date)
    {
        validate_reuse(&prior, &plan)?;
        require(
            same(&string(sha256(config_path)?), item(&plan, "config_sha256")?),
            "Configuration changed while waiting for the ClinVar refresh lock",
        )?;
        let updates = updates_of(item(&prior, "config_updates")?)?;
        atomic_text(
            config_path,
            &render_config(&read_text(config_path)?, &updates)?,
        )?;
        plan.set("action", string("reuse"));
        return Ok(plan);
    }
    require_workspace_owner(&plan)?;
    write_json(&join(&folder, "refresh.plan.json"), &plan)?;
    let inputs = acquire_inputs(&plan)?;
    let working = join(&folder, "install.config.yaml");
    let explicit = truth(item(&plan, "explicit_vcf")?);
    let ucsc = text_of(item(&plan, "ucsc_assembly")?);
    let mut updates = Vec::new();
    if !explicit {
        let vcf = join(&folder, &format!("clinvar.{ucsc}.vep.vcf.gz"));
        updates.push(("clinvar_vcf_dir".to_owned(), string(folder.clone())));
        updates.push(("clinvar_vcf".to_owned(), string(vcf)));
        for (key, suffix) in [
            ("clinvar_aa_stat", "aa_change.stats.rkyv"),
            ("clinvar_splice_stat", "splice_change.stats.rkyv"),
        ] {
            updates.push((
                key.to_owned(),
                string(join(&folder, &format!("clinvar.{ucsc}.vep.{suffix}"))),
            ));
        }
    }
    atomic_text(
        &working,
        &render_config(&read_text(config_path)?, &updates)?,
    )?;
    let source = text_of(item(item(item(&inputs, "sources")?, "vcf")?, "path")?);
    let mut arguments = vec![working.as_str()];
    if !explicit {
        arguments.extend(["--source-vcf", source.as_str()]);
    }
    run_original(source_root, "ClinVar_VCF_deploy", &arguments)?;
    let annotated_config = load_config(&working)?;
    let annotated = selected_vcf(
        &annotated_config,
        &path_of(item(&load_config(&working)?, "clinvar_vcf_dir")?)?,
        &ucsc,
    )?;
    require(
        same(
            &string(date_text(local_file_date(&annotated)?)),
            &latest_date,
        ),
        "Managed ClinVar annotation has an older/different fileDate; refusing reuse",
    )?;
    run_original(source_root, "ClinVar_AA_stat", &[&working])?;
    let sources = item(&inputs, "sources")?;
    let summary = text_of(item(item(sources, "summary")?, "path")?);
    let mane = text_of(item(item(sources, "mane")?, "path")?);
    let output = text_of(item(&plan, "output_path")?);
    run_original(
        source_root,
        "ClinVar_Gene_stat_regenerate",
        &[&working, &summary, &mane, &output],
    )?;
    complete_refresh(&plan, &working, &inputs)?;
    Ok(plan)
}
