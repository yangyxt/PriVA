//! `run`, `adopt`, `publish` and the `call-args` helper of the step driver.
use super::plan::normalize_pin;
use super::{
    Options, Result, absolute, clear, compact, copy_preserving, delete_key, fail, hostname,
    identity, now_iso, object, pretty, read_json, realpath_m, set_key, set_writable, step, steps,
    string, strings, text, today_utc, truthy, write_file, yaml_value,
};
use crate::json_text::{self, Json};
use std::collections::BTreeSet;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// Runs the installer functions of a step in one shell, as the shell runner
/// did: the installer sourced once, each call's arguments read from the
/// private configuration as the calls before it left it, pipefail off while an
/// installer function runs (they were written without it: checks such as
/// `bcftools view FILE | head` end their writer early). Its own variables have
/// long names so that an installer function's global cannot change them.
const DRIVER: &str = r#"set -o pipefail
_install_task_scripts=$1 _install_task_priva=$2 _install_task_registry=$3 _install_task_step=$4
_install_task_config=$5 _install_task_root=$6 _install_task_call=0
source "$_install_task_scripts/install_utils.sh" || { echo "install_task: cannot source install_utils.sh" >&2; exit 1; }
while :; do
    "$_install_task_priva" install-task call-args --registry "$_install_task_registry" --step "$_install_task_step" \
        --index "$_install_task_call" --config "$_install_task_config" --source-root "$_install_task_root" \
        > "$TMPDIR/install_task_call.args"
    _install_task_rc=$?
    ((_install_task_rc == 3)) && break
    ((_install_task_rc == 0)) || exit 1
    mapfile -d '' -t _install_task_args < "$TMPDIR/install_task_call.args"
    rm -f -- "$TMPDIR/install_task_call.args"
    echo "== $(date -Is) ${_install_task_args[*]}"
    _install_task_rc=0
    set +o pipefail; "${_install_task_args[@]}" || _install_task_rc=$?; set -o pipefail
    ((_install_task_rc == 0)) || { echo "install_task: ${_install_task_args[0]} failed (exit $_install_task_rc)" >&2; exit 1; }
    _install_task_call=$((_install_task_call + 1))
done
"#;

/// Environment for the programs a step starts.
struct Environment(Vec<(String, String)>);

impl Environment {
    fn apply<'a>(&self, command: &'a mut Command) -> &'a mut Command {
        for (name, value) in &self.0 {
            command.env(name, value);
        }
        command
    }
}

/// A receipt's keys as (KEY, VALUE). A step that reuses a file shipped with
/// the repository records the file's path in the source tree it ran from, and
/// lists the key under "shipped" with the path inside the tree. From another
/// source tree (a new checkout) the recorded path is no longer the shipped
/// file to the installer's functions, so the reader names the file in its own
/// tree (F7 after F5 in job 1546329).
fn receipt_keys(receipt: &Json, root: &Path) -> Vec<(String, String)> {
    let root = realpath_m(root);
    let shipped = receipt.get("shipped");
    let Some(keys) = receipt.get("keys").and_then(Json::as_object) else {
        return Vec::new();
    };
    keys.iter()
        .map(|(key, value)| {
            let inside = shipped.and_then(|s| s.get(key)).filter(|v| truthy(Some(v)));
            let value = match inside {
                Some(relative) => format!("{}/{}", root.display(), json_scalar_text(relative)),
                None => json_scalar_text(value),
            };
            (key.clone(), value)
        })
        .collect()
}

/// jq's `"\(.)"` of a value: a string as itself, anything else as JSON.
fn json_scalar_text(value: &Json) -> String {
    match value {
        Json::Str(text) => text.clone(),
        other => compact(other).unwrap_or_default(),
    }
}

/// The folder of a shipped file the step may reuse in place
/// (`priva_bundled_cache_reuse --key KEY --path PATH` without an assembly).
fn bundled_reuse(key: &str, selected: &str, root: &Path) -> bool {
    let relative = match key {
        "gofcards_exact_gof_cache" => "gofcards/gofcards_exact_gof.json.gz",
        "gene_nonlof_mechanism_json" => {
            "patho_mechanism/gene_nonlof_mechanism_curated_assertions.json.gz"
        }
        "hpo_condition_mechanism_json" => "patho_mechanism/hpo_condition_mechanism_cache.json.gz",
        "ddg2p_mechanism_evidence" => "patho_mechanism/gene_pathogenic_mechanism_evidence.tsv",
        "clingen_erepo" => "clingen/clingen_erepo.tsv",
        // The ClinGen keys need an assembly, which a copy does not name.
        _ => return false,
    };
    if selected.is_empty() {
        return false;
    }
    let canonical = absolute(&realpath_m(root).join("data").join(relative));
    let selected = Path::new(selected);
    absolute(selected) == canonical
        && realpath_m(selected) == canonical
        && std::fs::metadata(selected).is_ok_and(|m| m.is_file() && m.len() > 0)
}

fn current_dir() -> Result<PathBuf> {
    std::env::current_dir().map_err(|e| fail(format!("cannot read the working folder: {e}")))
}

fn make_dir(path: &Path) -> Result<()> {
    std::fs::create_dir_all(path)
        .map_err(|e| fail(format!("cannot create {}: {e}", path.display())))
}

fn this_program() -> Result<PathBuf> {
    std::env::current_exe().map_err(|e| fail(format!("cannot find the priva program: {e}")))
}

fn nonempty_env(name: &str) -> Option<String> {
    std::env::var(name).ok().filter(|v| !v.is_empty())
}

pub(crate) fn run(options: &Options) -> Result<()> {
    let registry_path = absolute(&options.path("registry")?);
    let plan = read_json(&options.path("plan")?)?;
    let id = options.one("step")?.to_owned();
    let lineage = options.one("lineage")?.to_owned();
    let folder = options.path("folder")?;
    let receipt_path = options.path("receipt")?;
    let root = options.source_root()?;
    let registry = read_json(&registry_path)?;
    let entry = step(&registry, &id)?;
    let stub = options.stub;
    if folder.join("receipt.json").exists() {
        return Err(fail(format!(
            "{} already holds a receipt",
            folder.display()
        )));
    }
    let here = current_dir()?;
    let work = folder.join("data");
    clear(&work)?;
    make_dir(&work)?;
    make_dir(&here.join("tmp"))?;
    let private = here.join("private.yaml");
    let base = text(&plan, "base_config");
    std::fs::copy(&base, &private).map_err(|e| fail(format!("cannot copy {base}: {e}")))?;
    let work_text = work.to_string_lossy().into_owned();
    set_key(&private, "cache_root", &work_text)?;
    set_key(&private, "tmp_dir", &here.join("tmp").to_string_lossy())?;
    if let Some(program) = nonempty_env("PRIVA_BINARY") {
        set_key(&private, "priva_binary", &program)?;
    }
    if let Some(threads) = nonempty_env("INSTALL_THREADS") {
        set_key(&private, "threads", &threads)?;
    }
    let mut environment = Environment(vec![(
        "TMPDIR".into(),
        here.join("tmp").to_string_lossy().into_owned(),
    )]);

    // Files you supply for another step (its user_inputs) reach the steps that
    // need them through that step's receipt. Resolving folders would copy them
    // into this step's folder: 47 GB for the SpliceAI2 pair in every running
    // step (job 1542391). They are left out of that and put back as configured.
    let own: BTreeSet<String> = ["user_inputs", "outputs", "optional_outputs", "copies"]
        .iter()
        .flat_map(|field| strings(entry, field))
        .collect();
    let supplied_keys: BTreeSet<String> = steps(&registry)?
        .iter()
        .flat_map(|s| strings(s, "user_inputs"))
        .collect();
    let mut supplied = Vec::new();
    for key in supplied_keys.iter().filter(|k| !own.contains(*k)) {
        let value = yaml_value(&private, key)?;
        if value.is_empty() {
            continue;
        }
        delete_key(&private, key)?;
        supplied.push((key.clone(), value));
    }
    if !stub {
        resolve_folders(&private, &root, &here.join("resolve.log"), &environment)?;
    }
    for (key, value) in &supplied {
        set_key(&private, key, value)?;
    }
    for input in options.all("input") {
        let receipt = read_json(Path::new(input))?;
        for (key, value) in receipt_keys(&receipt, &root) {
            set_key(&private, &key, &value)?;
        }
    }
    for key in strings(entry, "copies") {
        let value = yaml_value(&private, &key)?;
        if !Path::new(&value).is_file() {
            return Err(fail(format!("{key} ({value}) is not a file to copy")));
        }
        // A file shipped with the repository is reused where it is and never
        // rewritten. A copy is no longer that file to the step's function,
        // which then looks for build files nobody made (F7 in job 1545925).
        if !stub && bundled_reuse(&key, &value, &root) {
            continue;
        }
        let copied = work.join("copied");
        make_dir(&copied)?;
        let source = Path::new(&value);
        let name = source
            .file_name()
            .ok_or_else(|| fail(format!("{value} has no file name")))?;
        copy_preserving(source, &copied.join(name))?;
        for index in [format!("{value}.tbi"), format!("{value}.csi")] {
            let index = Path::new(&index);
            if index.is_file() {
                copy_preserving(index, &copied.join(index.file_name().unwrap_or_default()))?;
            }
        }
        set_key(&private, &key, &copied.join(name).to_string_lossy())?;
    }
    std::fs::copy(&private, here.join("private.before.yaml"))
        .map_err(|e| fail(format!("cannot copy {}: {e}", private.display())))?;
    environment
        .0
        .push(("PRIVA_CACHE_INSTALL_INTERNAL".into(), "1".into()));
    environment.0.push((
        "PRIVA_CACHE_INSTALL_BUILD_CONFIG".into(),
        private.to_string_lossy().into_owned(),
    ));

    if stub {
        for key in strings(entry, "outputs") {
            let path = work.join("stub").join(&key);
            if key.ends_with("_dir") || key == "loftee_repo" {
                make_dir(&path)?;
            } else {
                make_dir(&work.join("stub"))?;
                write_file(&path, &format!("{id} {key}\n"))?;
            }
            set_key(&private, &key, &path.to_string_lossy())?;
        }
    } else {
        let scripts = realpath_m(&root).join("scripts");
        let status = environment
            .apply(
                Command::new("bash")
                    .arg("-c")
                    .arg(DRIVER)
                    .arg("install-task"),
            )
            .arg(&scripts)
            .arg(this_program()?)
            .arg(&registry_path)
            .arg(&id)
            .arg(&private)
            .arg(&root)
            .status()
            .map_err(|e| fail(format!("cannot run bash: {e}")))?;
        if !status.success() {
            return Err(fail(format!("step {id} failed")));
        }
    }

    let mut keys: Vec<(String, Json)> = Vec::new();
    let mut missing = String::new();
    for key in strings(entry, "outputs") {
        let value = yaml_value(&private, &key)?;
        if value.is_empty() || !Path::new(&value).exists() {
            missing.push(' ');
            missing.push_str(&key);
        }
        keys.push((key, string(value)));
    }
    if !missing.is_empty() {
        return Err(fail(format!("outputs missing after the step:{missing}")));
    }
    // An optional output is recorded only when the step made it.
    for key in strings(entry, "optional_outputs") {
        let value = yaml_value(&private, &key)?;
        if !value.is_empty() && Path::new(&value).exists() {
            match keys.iter_mut().find(|(name, _)| *name == key) {
                Some(slot) => slot.1 = string(value),
                None => keys.push((key, string(value))),
            }
        }
    }
    let files: Vec<Json> = keys
        .iter()
        .map(|(_, v)| string(identity(&json_scalar_text(v))))
        .collect();
    // Outputs that lie in the source tree and outside this step's folder are
    // files shipped with the repository (see receipt_keys).
    let tree = realpath_m(&root);
    let inside = realpath_m(&folder);
    let mut shipped = Vec::new();
    for (key, value) in &keys {
        let value = realpath_m(Path::new(&json_scalar_text(value)));
        if value.starts_with(&tree)
            && value != tree
            && !(value.starts_with(&inside) && value != inside)
        {
            let relative = value
                .strip_prefix(&tree)
                .unwrap_or(&value)
                .to_string_lossy()
                .into_owned();
            shipped.push((key.clone(), string(relative)));
        }
    }
    let mut upstream = Vec::new();
    for input in options.all("input") {
        let receipt = read_json(Path::new(input))?;
        upstream.push(object(vec![
            (
                "id".into(),
                receipt.get("id").cloned().unwrap_or(Json::Null),
            ),
            (
                "lineage".into(),
                receipt.get("lineage").cloned().unwrap_or(Json::Null),
            ),
        ]));
    }
    let diff = Command::new("diff")
        .arg(here.join("private.before.yaml"))
        .arg(&private)
        .stderr(Stdio::null())
        .output()
        .map(|o| o.stdout)
        .unwrap_or_default();
    std::fs::write(here.join("config_changes.diff"), diff)?;
    if !truthy(entry.get("writable")) {
        set_writable(&work, false);
    }
    let receipt = object(vec![
        ("schema".into(), Json::Int(1)),
        ("id".into(), string(&id)),
        ("name".into(), string(text(entry, "name"))),
        ("lineage".into(), string(lineage)),
        ("folder".into(), string(folder.to_string_lossy())),
        ("keys".into(), object(keys)),
        ("shipped".into(), object(shipped)),
        ("files".into(), Json::Array(files)),
        ("inputs".into(), Json::Array(upstream)),
        ("finished".into(), string(now_iso())),
        ("host".into(), string(hostname())),
        ("stub".into(), Json::Bool(stub)),
    ]);
    write_file(&receipt_path, &pretty(&receipt)?)
}

/// The step's folders under its own cache_root (`priva_configure_cache_paths`
/// and `priva_configure_portable_paths`); the output goes to resolve.log and
/// is shown when a command fails.
fn resolve_folders(
    private: &Path,
    root: &Path,
    log: &Path,
    environment: &Environment,
) -> Result<()> {
    let base = realpath_m(root);
    let mut file = std::fs::File::create(log)
        .map_err(|e| fail(format!("cannot write {}: {e}", log.display())))?;
    for mode in ["managed-roots", "resolve"] {
        let output = environment
            .apply(&mut Command::new(this_program()?))
            .args(["installer-config-paths", mode, "--config"])
            .arg(private)
            .arg("--source-root")
            .arg(&base)
            .env("PYTHONHASHSEED", "0")
            .env("PRIVA_DATA_DIR", base.join("data"))
            .env("PRIVA_SCRIPT_DIR", base.join("scripts"))
            .output()
            .map_err(|e| fail(format!("cannot run priva installer-config-paths: {e}")))?;
        file.write_all(&output.stdout)?;
        file.write_all(&output.stderr)?;
        if !output.status.success() {
            let _ = std::io::stderr().write_all(&std::fs::read(log).unwrap_or_default());
            let what = if mode == "resolve" {
                "paths"
            } else {
                "folders"
            };
            return Err(fail(format!("resolving {what} failed")));
        }
    }
    Ok(())
}

/// call-args: the arguments of a step's call INDEX, each ended by a NUL, with
/// its tokens filled in from the configuration: `{config}` is the
/// configuration, `{script_dir}` the scripts folder, `{KEY}` the value of KEY.
/// Exit status 3 when the step has no such call.
pub(crate) fn call_args(options: &Options) -> Result<()> {
    let registry = read_json(&options.path("registry")?)?;
    let entry = step(&registry, options.one("step")?)?;
    let config = options.path("config")?;
    let root = options.source_root()?;
    let index: usize = options
        .one("index")?
        .parse()
        .map_err(|_| fail("--index must be a number"))?;
    let calls = match entry.get("calls") {
        Some(Json::Array(calls)) => calls,
        _ => return Err(fail("the step has no calls")),
    };
    let Some(Json::Array(tokens)) = calls.get(index) else {
        std::process::exit(3);
    };
    let mut out = Vec::new();
    for token in tokens {
        let token = json_scalar_text(token);
        let value = match token.as_str() {
            "{config}" => config.to_string_lossy().into_owned(),
            "{script_dir}" => realpath_m(&root)
                .join("scripts")
                .to_string_lossy()
                .into_owned(),
            t if t.len() >= 2 && t.starts_with('{') && t.ends_with('}') => {
                yaml_value(&config, &t[1..t.len() - 1])?
            }
            t => t.to_owned(),
        };
        out.extend_from_slice(value.as_bytes());
        out.push(0);
    }
    std::io::stdout().write_all(&out)?;
    Ok(())
}

/// When FILE lies in the data folder of another PriVA source tree and this
/// source tree ships the same bytes at the same place, the place inside the
/// tree (data/...). An existing installation names such a file by its path in
/// the source tree it was installed from; adopted as it is, the path is not
/// the shipped file to the installer's functions of this tree (see
/// receipt_keys), so the receipt lists the key under "shipped".
fn shipped_twin(file: &str, root: &Path) -> Option<String> {
    let value = realpath_m(Path::new(file)).to_string_lossy().into_owned();
    let root = realpath_m(root).to_string_lossy().into_owned();
    if !Path::new(&value).is_file() || value.starts_with(&format!("{root}/")) {
        return None;
    }
    let mut prefix = value.clone();
    while let Some(at) = prefix.rfind("/data/") {
        prefix.truncate(at);
        let relative = value[prefix.len() + 1..].to_owned();
        let tree = Path::new(&prefix);
        let is_tree = tree.join("scripts").join("install_utils.sh").is_file()
            || tree.join("conf").join("install_steps.json").is_file();
        let ours = Path::new(&root).join(&relative);
        if is_tree && ours.is_file() && same_bytes(Path::new(&value), &ours) {
            return Some(relative);
        }
    }
    None
}

fn same_bytes(a: &Path, b: &Path) -> bool {
    match (std::fs::read(a), std::fs::read(b)) {
        (Ok(x), Ok(y)) => x == y,
        _ => false,
    }
}

pub(crate) fn adopt(options: &Options) -> Result<()> {
    let plan = read_json(&options.path("plan")?)?;
    let id = options.one("step")?;
    let lineage = options.one("lineage")?;
    let receipt_path = options.path("receipt")?;
    let root = options.source_root()?;
    let adopted = plan
        .get("steps")
        .and_then(|s| s.get(id))
        .and_then(|s| s.get("adopt"))
        .filter(|a| !matches!(a, Json::Null))
        .ok_or_else(|| fail(format!("{id} is not adopted in the plan")))?;
    let keys = adopted
        .get("keys")
        .cloned()
        .unwrap_or(Json::Object(Vec::new()));
    let mut shipped = Vec::new();
    for (key, value) in keys.as_object().into_iter().flatten() {
        if let Some(relative) = shipped_twin(&json_scalar_text(value), &root) {
            shipped.push((key.clone(), string(relative)));
        }
    }
    let receipt = object(vec![
        ("schema".into(), Json::Int(1)),
        ("id".into(), string(id)),
        ("lineage".into(), string(lineage)),
        ("adopted".into(), Json::Bool(true)),
        ("keys".into(), keys),
        ("shipped".into(), object(shipped)),
        (
            "files".into(),
            adopted.get("identity").cloned().unwrap_or(Json::Null),
        ),
        ("inputs".into(), Json::Array(Vec::new())),
        ("finished".into(), string(now_iso())),
    ]);
    write_file(&receipt_path, &pretty(&receipt)?)
}

/// publish: one configuration from your base configuration and every step's
/// keys (in registry order, so a later step's key wins), native statistics
/// archives converted under F/data, and resources.lock.json naming every
/// step's lineage and keys.
pub(crate) fn publish(options: &Options) -> Result<()> {
    let registry = read_json(&options.path("registry")?)?;
    let plan = read_json(&options.path("plan")?)?;
    let lineage = options.one("lineage")?;
    let folder = options.path("folder")?;
    let root = options.source_root()?;
    let here = current_dir()?;
    let work = folder.join("data");
    clear(&work)?;
    make_dir(&work)?;
    make_dir(&here.join("tmp"))?;
    let merged = work.join("config.yaml");
    let base = text(&plan, "base_config");
    std::fs::copy(&base, &merged).map_err(|e| fail(format!("cannot copy {base}: {e}")))?;
    set_key(&merged, "cache_root", &work.to_string_lossy())?;
    let mut by_id: Vec<(String, Json)> = Vec::new();
    for input in options.all("input") {
        let receipt = read_json(Path::new(input))?;
        let id = text(&receipt, "id");
        by_id.retain(|(name, _)| *name != id);
        by_id.push((id, receipt));
    }
    let mut published = Vec::new();
    for entry in steps(&registry)? {
        let id = text(entry, "id");
        let receipt = by_id
            .iter()
            .find(|(name, _)| *name == id)
            .map(|(_, r)| r)
            .ok_or_else(|| fail(format!("no receipt for step {id}")))?;
        for (key, value) in receipt_keys(receipt, &root) {
            set_key(&merged, &key, &value)?;
        }
        let adopted = match receipt.get("adopted") {
            None | Some(Json::Null) | Some(Json::Bool(false)) => Json::Bool(false),
            Some(value) => value.clone(),
        };
        published.push((
            id,
            object(vec![
                (
                    "name".into(),
                    receipt.get("name").cloned().unwrap_or(Json::Null),
                ),
                (
                    "lineage".into(),
                    receipt.get("lineage").cloned().unwrap_or(Json::Null),
                ),
                ("adopted".into(), adopted),
                (
                    "keys".into(),
                    receipt.get("keys").cloned().unwrap_or(Json::Null),
                ),
                (
                    "receipt".into(),
                    receipt.get("folder").cloned().unwrap_or(Json::Null),
                ),
            ]),
        ));
    }
    if let Some(program) = nonempty_env("PRIVA_BINARY") {
        set_key(&merged, "priva_binary", &program)?;
    }
    // The native files are written below cache_root in the folders these keys
    // name; the steps' own folders lie outside it, so they are removed here.
    // vep_plugins_cachedir is restored below: the VEP step checks it exists.
    for key in [
        "clinvar_vcf_dir",
        "vep_plugins_cachedir",
        "gnomad_vcf_dir",
        "cadd_parent_dir",
        "loftee_parent_dir",
        "hpo_raw_dir",
        "gene_mechanism_cache_dir",
        "gene_mechanism_raw_dir",
        "gofcards_workdir",
    ] {
        delete_key(&merged, key)?;
    }
    make_dir(&work.join("VEP_plugin_cache"))?;
    if !options.stub {
        let program = nonempty_env("PRIVA_BINARY")
            .ok_or_else(|| fail("PRIVA_BINARY must name the priva program"))?;
        for stage in ["prepare", "publish"] {
            let status = Command::new(&program)
                .args(["cache-install", stage, "--config"])
                .arg(&merged)
                .arg("--build-config")
                .arg(work.join("build.yaml"))
                .arg("--source-root")
                .arg(realpath_m(&root))
                .env("TMPDIR", here.join("tmp"))
                .status()
                .map_err(|e| fail(format!("cannot run {program}: {e}")))?;
            if !status.success() {
                return Err(fail(format!("cache-install {stage} failed")));
            }
        }
    }
    set_key(
        &merged,
        "vep_plugins_cachedir",
        &work.join("VEP_plugin_cache").to_string_lossy(),
    )?;
    std::fs::copy(&merged, here.join("config.yaml"))
        .map_err(|e| fail(format!("cannot copy {}: {e}", merged.display())))?;
    // Pinned steps keep the release they were built from until --refresh. A
    // release already pinned keeps the day it was first taken.
    if let Some(pins) = options.get("pins").filter(|p| !p.is_empty()) {
        let pins = Path::new(pins);
        let mut current = if pins.is_file() {
            read_json(pins)?
        } else {
            Json::Object(Vec::new())
        };
        let today = today_utc();
        if let Some(planned) = plan.get("steps").and_then(Json::as_object) {
            for (id, entry) in planned {
                let release = entry.get("release").cloned().unwrap_or(Json::Null);
                if !truthy(entry.get("pinned"))
                    || matches!(release, Json::Null)
                    || matches!(&release, Json::Str(t) if t == "initial")
                {
                    continue;
                }
                let old = normalize_pin(current.get(id));
                let keep = !matches!(old, Json::Null)
                    && old
                        .get("release")
                        .is_some_and(|r| json_text::py_equal(r, &release))
                    && old
                        .get("recorded")
                        .is_some_and(|r| !matches!(r, Json::Null));
                let value = if keep {
                    old
                } else {
                    object(vec![
                        ("release".into(), release),
                        ("recorded".into(), string(today.clone())),
                    ])
                };
                current.set(id, value);
            }
        }
        let partial = PathBuf::from(format!("{}.partial", pins.display()));
        write_file(&partial, &(compact(&current)? + "\n"))?;
        std::fs::rename(&partial, pins)
            .map_err(|e| fail(format!("cannot write {}: {e}", pins.display())))?;
    }
    let lock = object(vec![
        ("schema".into(), Json::Int(1)),
        ("lineage".into(), string(lineage)),
        ("published".into(), string(now_iso())),
        (
            "priva_identity".into(),
            string(text(&plan, "priva_identity")),
        ),
        ("steps".into(), object(published)),
    ]);
    write_file(&here.join("resources.lock.json"), &pretty(&lock)?)
}
