//! The locked ClinVar annotation deployment: ownership, the bounded VEP header
//! check, and retention of every overwrite member before a cold rebuild.
//! Nothing here deletes annotation bytes; retained members are moved aside
//! under a journal written before the first rename.
use super::profile::first_vcf_stem;
use super::*;
use crate::file_lock::{ExclusiveLock, Wait};
use crate::installer_prerequisites::paths::{
    abspath, exists, is_dir, is_file, is_relative_to, is_symlink, lexists,
};
use crate::python_path::absolute;
use std::os::unix::process::ExitStatusExt;
use std::process::Command;

const RECOVERY_KIND: &str = "priva_clinvar_annotation_recovery_v1";
const IDENTITY_FIELDS: [&str; 7] = ["dev", "ino", "mode", "nlink", "size", "mtime_ns", "ctime_ns"];

fn resolve(path: &str) -> Result<String> {
    crate::installer_prerequisites::paths::resolve(path)
}

/// `run_original(source_root, function, *arguments)`
pub(super) fn run_original(source_root: &str, function: &str, arguments: &[&str]) -> Result<()> {
    let mut command = vec![
        "bash".to_owned(),
        join(source_root, "scripts/install_utils.sh"),
        function.to_owned(),
    ];
    command.extend(arguments.iter().map(|value| (*value).to_owned()));
    io::stdout().flush().ok();
    let status = Command::new(&command[0])
        .args(&command[1..])
        .status()
        .map_err(|e| os_error(e, Path::new("bash"), None))?;
    let listed = format!(
        "[{}]",
        command.iter().map(|value| repr(value)).collect::<Vec<_>>().join(", ")
    );
    match (status.code(), status.signal()) {
        (Some(0), _) => Ok(()),
        (Some(code), _) => Err(error(
            "CalledProcessError",
            format!("Command '{listed}' returned non-zero exit status {code}."),
        )),
        (None, Some(signal)) => Err(error(
            "CalledProcessError",
            format!(
                "Command '{listed}' died with {}.",
                crate::python_subprocess::signal_repr(signal)
            ),
        )),
        _ => Err(error("CalledProcessError", format!("Command '{listed}' failed."))),
    }
}

/// `build_lock(folder)`: an exclusive lock on `folder/.refresh.lock`.
pub(super) fn build_lock(folder: &str) -> Result<ExclusiveLock> {
    mkdirs(folder)?;
    let path = join(folder, ".refresh.lock");
    let file = fs::OpenOptions::new()
        .append(true)
        .create(true)
        .open(&path)
        .map_err(|e| os_error(e, Path::new(&path), None))?;
    ExclusiveLock::acquire(file, Wait::Blocking).map_err(|e| os_error(e, Path::new(&path), None))
}

/// `deployment_context(config_path)`: ownership without inferring a profile
/// from old file contents.
#[derive(Clone)]
pub(super) struct Context {
    config: Json,
    config_path: String,
    folder: String,
    target: String,
    root: Option<String>,
    assembly: &'static str,
    ucsc: &'static str,
    major: String,
    owned: bool,
}

impl Context {
    fn same(&self, other: &Context) -> bool {
        json_text::py_equal(&self.config, &other.config)
            && (
                &self.config_path,
                &self.folder,
                &self.target,
                &self.root,
                self.assembly,
                &self.major,
                self.owned,
            ) == (
                &other.config_path,
                &other.folder,
                &other.target,
                &other.root,
                other.assembly,
                &other.major,
                other.owned,
            )
    }
    fn root(&self) -> Result<String> {
        let Some(root) = &self.root else {
            return Err(error(
                "AttributeError",
                "'NoneType' object has no attribute 'resolve'",
            ));
        };
        resolve(root)
    }
    fn target_name(&self) -> &str {
        name_of(&self.target)
    }
}

pub(super) fn deployment_context(config_path: &str) -> Result<Context> {
    let config_path = absolute(Path::new(config_path))?;
    let config = load_config(&config_path)?;
    let assembly = assembly_of(&get_or_null(&config, "assembly")?)?;
    let major = vep_major(&config);
    let Some(assembly) = assembly.filter(|_| positive(&major)) else {
        return Err(value_error(
            "ClinVar deployment requires a configured assembly and numeric VEP major",
        ));
    };
    let folder = path_of(item(&config, "clinvar_vcf_dir")?)?;
    require(
        folder.starts_with('/'),
        "ClinVar deployment requires an absolute clinvar_vcf_dir",
    )?;
    let folder = abspath(&folder)?;
    let ucsc = ucsc(assembly);
    let target = abspath(&selected_vcf(&config, &folder, ucsc)?)?;
    let owned = target == join(&folder, &format!("clinvar.{ucsc}.vep.vcf.gz"));
    require(
        !is_symlink(&target)?,
        format!("Linked ClinVar annotation preserved without replacement: {target}"),
    )?;
    let root = match config.get("cache_root") {
        Some(Json::Str(text)) if text.starts_with('/') => Some(pathlib(text)),
        _ => None,
    };
    if owned {
        let inside = match &root {
            Some(root) => is_relative_to(&resolve(&folder)?, &resolve(root)?),
            None => false,
        };
        require(
            inside,
            "Owned ClinVar output must remain inside the configured cache_root; prior paths preserved",
        )?;
        require(
            format!("{}.vcf.gz", first_vcf_stem(&target)) == target,
            "Owned ClinVar output cannot have an earlier .vcf path component",
        )?;
    }
    Ok(Context {
        config,
        config_path,
        folder,
        target,
        root,
        assembly,
        ucsc,
        major,
        owned,
    })
}

/// `deployment_header(context, input_file)`
pub(super) fn deployment_header(context: &Context, input: Option<&str>) -> Result<Json> {
    crate::alpha_annotation_base::header::check_resource(
        Path::new(input.unwrap_or(&context.target)),
        &string(context.assembly),
        &string(context.major.clone()),
        "ClinVar",
    )
}

/// The failures a header check may report for a stable, owned, obsolete or
/// unreadable annotation: `except (ValueError, EOFError, gzip.BadGzipFile, zlib.error)`.
fn obsolete(failure: &(dyn std::error::Error + 'static)) -> bool {
    matches!(
        exception(failure).0,
        "ValueError"
            | "UnicodeDecodeError"
            | "AnnotationHeaderChangedError"
            | "EOFError"
            | "BadGzipFile"
            | "error"
    )
}

/// `retirement_members(context)`: fixed names replaced by acquisition/VEP.
fn retirement_members(context: &Context) -> Vec<String> {
    let (folder, target) = (&context.folder, &context.target);
    let prefix = first_vcf_stem(context.target_name());
    let payloads = [
        target.clone(),
        join(folder, &format!("clinvar.{}.vcf.gz", context.assembly)),
        join(folder, &format!("clinvar.{}.vcf.gz", context.ucsc)),
    ];
    let mut members = Vec::new();
    for path in &payloads {
        for suffix in ["", ".tbi", ".csi"] {
            members.push(format!("{path}{suffix}"));
        }
    }
    members.extend(payloads.iter().map(|path| without_suffix(path)));
    for suffix in [
        ".installation-candidate.pending.json",
        ".installation-candidate.ready.json",
        ".gnomAD.publish.pending",
    ] {
        members.push(format!("{target}{suffix}"));
    }
    for suffix in [".gnomAD.log", ".vep.stats.html"] {
        members.push(join(folder, &format!("{prefix}{suffix}")));
    }
    unique(members)
}
/// `Path(str(path)[:-3])`
fn without_suffix(path: &str) -> String {
    let end = path.char_indices().rev().nth(2).map_or(0, |(at, _)| at);
    pathlib(&path[..end])
}
fn unique(paths: Vec<String>) -> Vec<String> {
    let mut seen = Vec::new();
    for path in paths {
        if !seen.contains(&path) {
            seen.push(path);
        }
    }
    seen
}

/// `shard_retirement_members(context, chromosomes)`: the writer's literal names.
fn shard_retirement_members(context: &Context, chromosomes: &[String]) -> Result<Vec<String>> {
    let prefix = first_vcf_stem(context.target_name());
    let mut members = Vec::new();
    for chrom in chromosomes {
        require(
            !chrom.is_empty() && !chrom.contains(['/', '\0', '\r', '\n', '\t']),
            "Unsafe literal ClinVar shard name; all paths preserved",
        )?;
        for ending in [".vcf.gz", ".gnomAD.vcf.gz"] {
            let shard = join(&context.folder, &format!("{prefix}.{chrom}{ending}"));
            for suffix in ["", ".tbi", ".csi"] {
                members.push(format!("{shard}{suffix}"));
            }
            members.push(without_suffix(&shard));
            if ending == ".gnomAD.vcf.gz" {
                members.push(format!("{shard}.gnomAD.publish.pending"));
            }
        }
    }
    Ok(unique(members))
}

/// `retirement_identity(path)`: raw lstat fields of a regular file.
fn retirement_identity(path: &str) -> Result<Json> {
    let info = lstat(path)?;
    require(
        info.file_type().is_file(),
        format!("Nonregular/linked ClinVar retirement member preserved: {path}"),
    )?;
    let values = [
        info.dev() as i128,
        info.ino() as i128,
        info.mode() as i128,
        info.nlink() as i128,
        info.len() as i128,
        nanoseconds(info.mtime(), info.mtime_nsec()),
        nanoseconds(info.ctime(), info.ctime_nsec()),
    ];
    Ok(Json::Object(
        IDENTITY_FIELDS
            .iter()
            .zip(values)
            .map(|(name, value)| ((*name).to_owned(), string(value.to_string())))
            .collect(),
    ))
}

/// `persisted_identity_matches(recorded, observed)`: a retry observation
/// compared without the node-dependent device number.
fn persisted_identity_matches(recorded: &Json, observed: &Json) -> bool {
    let fields = |value: &Json| -> Option<Vec<(String, String)>> {
        let Json::Object(items) = value else {
            return None;
        };
        let mut found: Vec<(String, String)> = items
            .iter()
            .map(|(key, value)| value.as_str().map(|text| (key.clone(), text.to_owned())))
            .collect::<Option<_>>()?;
        found.sort();
        let mut names: Vec<&str> = found.iter().map(|(key, _)| key.as_str()).collect();
        names.dedup();
        let mut expected = IDENTITY_FIELDS.to_vec();
        expected.sort();
        (names == expected).then_some(found)
    };
    let (Some(recorded), Some(observed)) = (fields(recorded), fields(observed)) else {
        return false;
    };
    let digits = |values: &[(String, String)]| {
        values.iter().any(|(key, value)| {
            key == "dev" && !value.is_empty() && value.bytes().all(|b| b.is_ascii_digit())
        })
    };
    digits(&recorded)
        && digits(&observed)
        && recorded
            .iter()
            .zip(&observed)
            .all(|(left, right)| left.0 == "dev" || left == right)
}

/// `record_retry_observation(pending, protected_sources=..., rebuilt_target=...)`
fn record_retry_observation(
    pending: &mut Json,
    protected: &Json,
    rebuilt: Option<&Json>,
) -> Result<()> {
    if pending.get("retry_observations").is_none() {
        pending.set("retry_observations", Json::Array(Vec::new()));
    }
    require(
        matches!(pending.get("retry_observations"), Some(Json::Array(_))),
        "Malformed ClinVar retry observations; previous identities preserved",
    )?;
    let mut value = object([
        ("observed_utc", string(utc_now())),
        ("phase", item(pending, "phase")?.clone()),
        ("protected_sources", protected.clone()),
    ]);
    if let Some(rebuilt) = rebuilt {
        value.set("rebuilt_target", rebuilt.clone());
    }
    let Some(Json::Array(observations)) = value_mut(pending, "retry_observations") else {
        return Err(value_error(
            "Malformed ClinVar retry observations; previous identities preserved",
        ));
    };
    observations.push(value);
    Ok(())
}

/// `retirement_inventory(context, source_vcf, members=..., protect_target=...)`:
/// foreign aliases are rejected before any archive or move.
fn retirement_inventory(
    context: &Context,
    source_vcf: Option<&str>,
    members: Option<Vec<String>>,
    protect_target: bool,
) -> Result<Vec<(String, Json)>> {
    let root = context.root()?;
    let mut protected = vec![context.config_path.clone()];
    if protect_target {
        let target = &context.target;
        protected.extend([target.clone(), format!("{target}.tbi"), format!("{target}.csi")]);
    }
    if let Some(source) = source_vcf {
        let source = absolute(Path::new(source))?;
        protected.extend([source.clone(), format!("{source}.tbi"), format!("{source}.csi")]);
    }
    let mut moves = Vec::new();
    for member in members.unwrap_or_else(|| retirement_members(context)) {
        require(
            is_relative_to(&resolve(&parent_of(&member))?, &root),
            format!("ClinVar retirement parent escapes cache_root; preserved: {member}"),
        )?;
        if !lexists(&member) {
            continue;
        }
        let before = retirement_identity(&member)?;
        for item in &protected {
            let aliased = resolve(&member)? == resolve(item)?
                || (exists(item)? && same_file(&member, item)?);
            require(
                !aliased,
                format!(
                    "ClinVar retirement would overwrite a selected config/source; preserved: {member}"
                ),
            )?;
        }
        moves.push((member, before));
    }
    Ok(moves)
}

fn deployment_pending_path(context: &Context) -> Result<String> {
    Ok(join(
        &resolve(&context.folder)?,
        &format!(
            ".clinvar-deploy-locks/{}/retirement.pending.json",
            context.target_name()
        ),
    ))
}

/// `deployment_pending(context)`: this target's recovery journal, if any.
fn deployment_pending(context: &Context) -> Result<Option<Json>> {
    let path = deployment_pending_path(context)?;
    require(!is_symlink(&path)?, "Linked ClinVar recovery journal preserved")?;
    if !exists(&path)? {
        return Ok(None);
    }
    let value = read_json(&path)?;
    let expected = [
        ("kind", string(RECOVERY_KIND)),
        ("target", string(context.target.clone())),
        ("config", string(context.config_path.clone())),
        ("assembly", string(context.assembly)),
        ("vep_major", string(context.major.clone())),
    ];
    let mut matches = true;
    for (key, wanted) in &expected {
        let found = get_or_null(&value, key)?;
        matches &= json_text::py_equal(&found, wanted);
    }
    let phase = get_or_null(&value, "phase")?;
    let known = matches!(
        phase.as_str(),
        Some("awaiting_vep" | "retiring_shards" | "shards_retained")
    );
    require(
        matches && known,
        "ClinVar recovery journal belongs to a different target/configuration; preserved",
    )?;
    Ok(Some(value))
}

/// `pinned_source_identity(context, source_vcf)`: a supplied source cannot
/// also occupy the controlled overwrite namespace.
fn pinned_source_identity(context: &Context, source_vcf: Option<&str>) -> Result<Json> {
    let mut result = Json::Object(Vec::new());
    let Some(source) = source_vcf else {
        return Ok(result);
    };
    let prefix = format!("{}.", first_vcf_stem(context.target_name()));
    let source = absolute(Path::new(source))?;
    let folder = resolve(&context.folder)?;
    for path in [source.clone(), format!("{source}.tbi"), format!("{source}.csi")] {
        if !lexists(&path) {
            continue;
        }
        let info = retirement_identity(&path)?;
        let resolved = resolve(&path)?;
        let aliased = info.get("nlink").and_then(Json::as_str) != Some("1")
            || (parent_of(&resolved) == folder && name_of(&resolved).starts_with(&prefix));
        require(
            !aliased,
            format!(
                "Pinned ClinVar source/index may alias a managed overwrite member; preserved: {path}"
            ),
        )?;
        result.set(&path, info);
    }
    Ok(result)
}

/// `retire_owned_cache(config_path, reason=..., ...)`: retain the exact
/// overwrite members. The caller holds the distinct deployment lock.
pub(super) fn retire_owned_cache(
    config_path: &str,
    reason: &str,
    source_vcf: Option<&str>,
    expected_target: Option<&Json>,
    expected_members: &[(String, Json)],
) -> Result<Option<Json>> {
    let context = deployment_context(config_path)?;
    require(
        context.owned,
        "Explicit ClinVar annotations are preserved; automatic retirement is limited to the managed default",
    )?;
    let pending = deployment_pending(&context)?;
    if let Some(pending) = &pending {
        require(
            get_or_null(pending, "phase")?.as_str() == Some("awaiting_vep"),
            "ClinVar shard retirement is pending; the new VEP output was preserved",
        )?;
    }
    let changed = || {
        error(
            "AnnotationHeaderChangedError",
            "ClinVar annotation changed after its header check; preserved",
        )
    };
    if lexists(&context.target) {
        let before = retirement_identity(&context.target)?;
        if expected_target.is_some_and(|expected| !json_text::py_equal(&before, expected)) {
            return Err(changed());
        }
        if let Err(failure) = deployment_header(&context, None) {
            if exception(failure.as_ref()).0 == "AnnotationHeaderChangedError"
                || !obsolete(failure.as_ref())
            {
                return Err(failure);
            }
            // Only a stable, owned unreadable/mismatched header may retire.
        }
        if !json_text::py_equal(&retirement_identity(&context.target)?, &before) {
            return Err(changed());
        }
    } else {
        require(
            expected_target.is_none(),
            "ClinVar annotation disappeared after its header check; no retirement started",
        )?;
    }
    let moves = retirement_inventory(&context, source_vcf, None, false)?;
    let protected = pinned_source_identity(&context, source_vcf)?;
    for (path, expected) in expected_members {
        require(
            json_text::py_equal(&retirement_identity(path)?, expected),
            "ClinVar repair source changed after its header check; no retirement started",
        )?;
    }
    let journal = deployment_pending_path(&context)?;
    match pending {
        None => {
            let value = object([
                ("kind", string(RECOVERY_KIND)),
                ("phase", string("awaiting_vep")),
                ("target", string(context.target.clone())),
                ("config", string(context.config_path.clone())),
                ("assembly", string(context.assembly)),
                ("vep_major", string(context.major.clone())),
                ("started_utc", string(utc_now())),
                ("scientific_qualified", Json::Bool(false)),
                ("protected_sources", protected),
            ]);
            write_json(&journal, &value)?;
        }
        Some(mut pending) => {
            let previous = get_or_null(&pending, "protected_sources")?;
            let unchanged = match (&previous, &protected) {
                (Json::Object(before), Json::Object(now)) => {
                    let mut keys: Vec<&String> = before.iter().map(|(key, _)| key).collect();
                    let mut current: Vec<&String> = now.iter().map(|(key, _)| key).collect();
                    keys.sort();
                    current.sort();
                    keys == current
                        && now.iter().all(|(path, value)| {
                            previous
                                .get(path)
                                .is_some_and(|recorded| persisted_identity_matches(recorded, value))
                        })
                }
                _ => false,
            };
            require(
                unchanged,
                "Pinned ClinVar source/index identity changed during recovery; preserved",
            )?;
            record_retry_observation(&mut pending, &protected, None)?;
            write_json(&journal, &pending)?;
        }
    }
    retain_members(&context, moves, reason, "primary")
}

/// `retain_members(context, moves, reason=..., phase=...)`: journal the exact
/// overwrite set before moving it; bytes are never deleted.
fn retain_members(
    context: &Context,
    mut moves: Vec<(String, Json)>,
    reason: &str,
    phase: &str,
) -> Result<Option<Json>> {
    if moves.is_empty() {
        return Ok(None);
    }
    let root = context.root()?;
    let parent = join(&root, ".installer-retained/clinvar-annotations");
    for directory in [parent_of(&parent), parent.clone()] {
        let unsafe_directory =
            is_symlink(&directory)? || (exists(&directory)? && !is_dir(&directory)?);
        require(
            !unsafe_directory,
            format!("Unsafe ClinVar retention directory preserved: {directory}"),
        )?;
    }
    mkdirs(&parent)?;
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs() as i128;
    let stamp = utc_text(now, 0);
    let prefix = format!(
        "{}{}{}T{}{}{}Z-",
        &stamp[..4],
        &stamp[5..7],
        &stamp[8..10],
        &stamp[11..13],
        &stamp[14..16],
        &stamp[17..19]
    );
    let archive = mkdtemp(&parent, &prefix)?;
    // Move the main annotation last; every overwrite member was preflighted.
    moves.sort_by_key(|(path, _)| *path == context.target);
    let journal = join(&archive, "retirement.json");
    let listed = moves
        .iter()
        .map(|(path, before)| {
            object([
                ("source", string(path.clone())),
                ("retained", string(join(&archive, name_of(path)))),
                ("before", before.clone()),
                ("moved", Json::Bool(false)),
            ])
        })
        .collect();
    let mut value = object([
        ("kind", string("priva_clinvar_annotation_retirement_v1")),
        ("status", string("retiring")),
        ("requested_assembly", string(context.assembly)),
        ("requested_vep_major", string(context.major.clone())),
        ("selected_vcf", string(context.target.clone())),
        ("reason", string(reason)),
        ("created_utc", string(utc_now())),
        ("scientific_qualified", Json::Bool(false)),
        ("phase", string(phase)),
        ("moves", Json::Array(listed)),
    ]);
    write_json(&journal, &value)?;
    let moved = (|| -> Result<()> {
        for index in 0..moves.len() {
            let (source, before) = &moves[index];
            let retained = join(&archive, name_of(source));
            let current = retirement_identity(source)?;
            require(
                json_text::py_equal(&current, before) && !lexists(&retained),
                "ClinVar retirement member changed or destination occupied; no rebuild started",
            )?;
            rename(source, &retained)?;
            let after = retirement_identity(&retained)?;
            if let Some(Json::Array(entries)) = value_mut(&mut value, "moves") {
                entries[index].set("moved", Json::Bool(true));
                entries[index].set("after", after.clone());
            }
            let kept = IDENTITY_FIELDS[..6].iter().all(|field| {
                after.get(field).and_then(Json::as_str) == before.get(field).and_then(Json::as_str)
            });
            require(
                kept,
                "ClinVar member changed while being retained; no rebuild started",
            )?;
            write_json(&journal, &value)?;
        }
        value.set("status", string("retained_before_rebuild"));
        value.set("completed_utc", string(utc_now()));
        write_json(&journal, &value)
    })();
    if let Err(failure) = moved {
        let (category, message) = exception(failure.as_ref());
        value.set("status", string("retirement_incomplete"));
        value.set("error", string(format!("{category}: {message}")));
        write_json(&journal, &value)?;
        return Err(failure);
    }
    println!("Retained previous managed ClinVar generation: {journal}");
    io::stdout().flush().ok();
    Ok(Some(value))
}

fn value_mut<'a>(value: &'a mut Json, key: &str) -> Option<&'a mut Json> {
    match value {
        Json::Object(items) => items.iter_mut().find(|(name, _)| name == key).map(|(_, v)| v),
        _ => None,
    }
}

/// `retire_rebuilt_shards(config_path, input_vcf, chromosomes)`: after the
/// unchanged gnomAD query/selection, before mtime-only shard reuse.
pub(super) fn retire_rebuilt_shards(
    config_path: &str,
    input_vcf: &str,
    chromosomes: &[String],
) -> Result<Option<Json>> {
    let context = deployment_context(config_path)?;
    let pending = deployment_pending(&context)?;
    require(
        context.owned && absolute(Path::new(input_vcf))? == context.target,
        "ClinVar recovery context does not match this selected input; no shards moved",
    )?;
    let Some(mut pending) = pending else {
        return Ok(None);
    };
    let journal = deployment_pending_path(&context)?;
    let before = retirement_identity(&context.target)?;
    deployment_header(&context, None)?;
    require(
        json_text::py_equal(&retirement_identity(&context.target)?, &before),
        "New ClinVar VEP output changed before shard retirement; preserved",
    )?;
    let mut observed = Json::Object(Vec::new());
    let sources = item(&pending, "protected_sources")?.clone();
    let Json::Object(sources) = sources else {
        return Err(error(
            "AttributeError",
            format!("'{}' object has no attribute 'items'", type_name(&sources)),
        ));
    };
    for (path, recorded) in &sources {
        let current = retirement_identity(path)?;
        require(
            persisted_identity_matches(recorded, &current),
            "Pinned ClinVar source/index changed before shard retirement; preserved",
        )?;
        observed.set(path, current);
    }
    let phase = text_of(item(&pending, "phase")?);
    let listed = Json::Array(chromosomes.iter().cloned().map(string).collect());
    if matches!(phase.as_str(), "retiring_shards" | "shards_retained") {
        let unchanged = persisted_identity_matches(item(&pending, "rebuilt_target")?, &before)
            && json_text::py_equal(item(&pending, "chromosomes")?, &listed);
        require(
            unchanged,
            "New ClinVar VEP output/selected contigs changed during shard retirement; preserved",
        )?;
        if phase == "shards_retained" {
            record_retry_observation(&mut pending, &observed, Some(&before))?;
            write_json(&journal, &pending)?;
            return Ok(Some(pending));
        }
    } else {
        pending.set("phase", string("retiring_shards"));
        pending.set("rebuilt_target", before.clone());
        pending.set("chromosomes", listed);
    }
    let members = shard_retirement_members(&context, chromosomes)?;
    let moves = retirement_inventory(&context, None, Some(members), true)?;
    for (member, _) in &moves {
        for (path, _) in &sources {
            require(
                !same_file(member, path)?,
                "A ClinVar shard aliases a supplied source/index; all shards preserved",
            )?;
        }
    }
    record_retry_observation(&mut pending, &observed, Some(&before))?;
    write_json(&journal, &pending)?;
    retain_members(
        &context,
        moves,
        "cold VEP rebuild: retire mtime-only gnomAD reuse members",
        "shards",
    )?;
    require(
        json_text::py_equal(&retirement_identity(&context.target)?, &before),
        "New ClinVar VEP output changed during shard retirement; preserved",
    )?;
    pending.set("phase", string("shards_retained"));
    pending.set("shards_retained_utc", string(utc_now()));
    write_json(&journal, &pending)?;
    Ok(Some(pending))
}

/// `deploy_cache(config_path, source_root, source_vcf=...)`: one target is
/// serialized independently of an enclosing gene-refresh lock.
pub(super) fn deploy_cache(config_path: &str, source_root: &str, source_vcf: Option<&str>) -> Result<()> {
    let context = deployment_context(config_path)?;
    let mut arguments = vec![config_path];
    if let Some(source) = source_vcf {
        arguments.extend(["--source-vcf", source]);
    }
    if !context.owned {
        if let Err(failure) = deployment_header(&context, None) {
            if !obsolete(failure.as_ref()) {
                return Err(failure);
            }
            return Err(value_error(format!(
                "Explicit ClinVar annotation incompatible or unverifiable; preserved {}: {}",
                context.target,
                exception(failure.as_ref()).1
            )));
        }
        return run_original(source_root, "_ClinVar_VCF_deploy", &arguments);
    }
    // 2026-10-23 CT: never reuse the enclosing refresh workspace's .refresh.lock.
    let lock_folder = join(
        &resolve(&context.folder)?,
        &format!(".clinvar-deploy-locks/{}", context.target_name()),
    );
    let lock_path = join(&lock_folder, ".refresh.lock");
    let linked = is_symlink(&parent_of(&lock_folder))?
        || is_symlink(&lock_folder)?
        || is_symlink(&lock_path)?
        || (exists(&lock_path)? && (!is_file(&lock_path)? || stat(&lock_path)?.nlink() != 1));
    require(!linked, "Linked ClinVar deployment lock directory preserved")?;
    retirement_inventory(&context, source_vcf, None, false)?;
    pinned_source_identity(&context, source_vcf)?;
    let _lock = build_lock(&lock_folder)?;
    // Waiting must not permit a different configuration/selection to inherit
    // the first target's lock or its managed ownership.
    require(
        deployment_context(config_path)?.same(&context),
        "ClinVar deployment configuration changed while waiting for its lock",
    )?;
    retirement_inventory(&context, source_vcf, None, false)?;
    deployment_pending(&context)?;
    if exists(&context.target)? {
        let before = retirement_identity(&context.target)?;
        let mut checked = Vec::new();
        let header = (|| -> Result<()> {
            deployment_header(&context, None)?;
            // The shared validity helper can repair a damaged gzip from a
            // newer plain sibling. Check that possible source beforehand,
            // retaining both generations if its provenance is incompatible.
            let plain = without_suffix(&context.target);
            if exists(&plain)? {
                let newer = {
                    let (a, b) = (stat(&plain)?, stat(&context.target)?);
                    nanoseconds(a.mtime(), a.mtime_nsec()) > nanoseconds(b.mtime(), b.mtime_nsec())
                };
                if newer {
                    checked.push((plain.clone(), retirement_identity(&plain)?));
                    deployment_header(&context, Some(&plain))?;
                }
            }
            Ok(())
        })();
        if let Err(failure) = header {
            if exception(failure.as_ref()).0 == "AnnotationHeaderChangedError"
                || !obsolete(failure.as_ref())
            {
                return Err(failure);
            }
            retire_owned_cache(
                config_path,
                &exception(failure.as_ref()).1,
                source_vcf,
                Some(&before),
                &checked,
            )?;
        }
    }
    run_original(source_root, "_ClinVar_VCF_deploy", &arguments)?;
    deployment_header(&deployment_context(config_path)?, None)?;
    let pending = deployment_pending_path(&context)?;
    if exists(&pending)? {
        unlink(&pending)?;
    }
    Ok(())
}
