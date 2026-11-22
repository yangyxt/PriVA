//! `adopt` and `publish`, and the receipt helpers `record` shares.
use super::plan::normalize_pin;
use super::{
    Options, Result, clear, compact, delete_key, fail, now_iso, object, pretty, read_json,
    realpath_m, set_key, steps, string, text, today_utc, truthy, write_file,
};
use crate::json_text::{self, Json};
use std::path::{Path, PathBuf};
use std::process::Command;

/// A receipt's keys as (KEY, VALUE). A step that reuses a file shipped with
/// the repository records the file's path in the source tree it ran from, and
/// lists the key under "shipped" with the path inside the tree. From another
/// source tree (a new checkout) the recorded path is no longer the shipped
/// file to the steps that read it, so the reader names the file in its own
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

fn current_dir() -> Result<PathBuf> {
    std::env::current_dir().map_err(|e| fail(format!("cannot read the working folder: {e}")))
}

fn make_dir(path: &Path) -> Result<()> {
    std::fs::create_dir_all(path)
        .map_err(|e| fail(format!("cannot create {}: {e}", path.display())))
}

fn nonempty_env(name: &str) -> Option<String> {
    std::env::var(name).ok().filter(|v| !v.is_empty())
}

/// What a step downstream keys on: the receipt's content, or its lineage
/// for a receipt without one (adopted, or written before contents).
pub(crate) fn receipt_content(receipt: &Json) -> String {
    receipt
        .get("content")
        .and_then(Json::as_str)
        .or_else(|| receipt.get("lineage").and_then(Json::as_str))
        .unwrap_or_default()
        .to_owned()
}

/// When FILE lies in the data folder of another PriVA source tree and this
/// source tree ships the same bytes at the same place, the place inside the
/// tree (data/...). An existing installation names such a file by its path in
/// the source tree it was installed from; adopted as it is, the path is not
/// the shipped file to this tree's steps (see receipt_keys), so the receipt
/// lists the key under "shipped". A source tree is a checkout with the step
/// registry, or an older one with the shell installer (scripts/install_utils.sh).
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
            "publish_digest".into(),
            string(text(&plan, "publish_digest")),
        ),
        ("steps".into(), object(published)),
    ]);
    write_file(&here.join("resources.lock.json"), &pretty(&lock)?)
}
