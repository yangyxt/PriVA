//! `record`: the receipt of a step that runs as its own Nextflow process
//! (registry `module`, `modules/local/install_steps/`).
//!
//! Such a step writes its results into `data/` of its task's work folder and
//! nowhere else. Its script runs the tools' own commands (curl, git, tabix,
//! samtools) and `priva` for PriVA's checks; its last line is
//!
//! ```text
//! priva install-task record --registry R --plan plan.json --step ID --lineage L
//!     --folder RESOURCE_ROOT/ID/L --receipt receipt.json --source-root CHECKOUT
//!     [--input inputs/N/receipt.json]... [--key KEY=data/PATH]... [--stub]
//! ```
//!
//! It checks that every output key of the step is given and exists, and
//! writes `receipt.json` naming each key where storeDir moves it when the
//! task succeeds: `data/PATH` becomes `FOLDER/data/PATH`. A key given as a
//! file staged into the task (a file you supply, kept as it is) names that
//! file's own path. The receipt has the fields `install-task run` writes, so
//! the steps after it, publish and adopt read it the same way; its content is
//! the same listing of `data/` (`content.rs`).
//!
//! The files under `data/` are made read-only. The folders stay writable:
//! storeDir moves `data/` into the step's folder with `mv`, and a read-only
//! folder cannot be moved to another parent by a user who is not root.
//!
//! A tool that names its own files by absolute path in a JSON file it
//! writes (the ClinVar download's inputs.ready.json, the reference's
//! .source.json) would name the task's work folder, which is gone once
//! storeDir has moved the results. `--relocate data/FILE` rewrites every
//! string of such a file that starts with this task's `data/` path to the
//! same path under `FOLDER/data/`, before the content is hashed.
//!
//! `--stub` writes `data/stub/KEY` for every output key not given, as
//! `install-task run --stub` does.
use super::content;
use super::run::receipt_content;
use super::{
    Options, Result, absolute, fail, hostname, now_iso, object, pretty, read_json, realpath_m,
    step, steps, string, strings, text, truthy, write_file,
};
use crate::json_text::Json;
use std::os::unix::fs::MetadataExt;
use std::path::{Path, PathBuf};

pub(crate) fn record(options: &Options) -> Result<()> {
    let registry = read_json(&absolute(&options.path("registry")?))?;
    let plan = read_json(&options.path("plan")?)?;
    let id = options.one("step")?.to_owned();
    let lineage = options.one("lineage")?.to_owned();
    let folder = absolute(&options.path("folder")?);
    let receipt_path = options.path("receipt")?;
    let root = options.source_root()?;
    let entry = step(&registry, &id)?;
    if entry.get("module").and_then(Json::as_str).is_none() {
        return Err(fail(format!(
            "step {id} runs installer functions (install-task run), not its own process"
        )));
    }
    if plan.get("steps").and_then(|s| s.get(&id)).is_none() {
        return Err(fail(format!("step {id} is not in the plan")));
    }
    let here = absolute(&std::env::current_dir()?);
    let work = here.join("data");
    std::fs::create_dir_all(&work)
        .map_err(|e| fail(format!("cannot create {}: {e}", work.display())))?;

    for file in options.all("relocate") {
        relocate(&absolute(&here.join(file)), &work, &folder.join("data"))?;
    }

    let outputs = strings(entry, "outputs");
    let optional = strings(entry, "optional_outputs");
    let mut given: Vec<(String, PathBuf)> = Vec::new();
    for pair in options.all("key") {
        let Some((key, value)) = pair.split_once('=') else {
            return Err(fail(format!("--key {pair}: use KEY=PATH")));
        };
        if !outputs.iter().chain(&optional).any(|k| k == key) {
            return Err(fail(format!(
                "{key} is not an output of step {id} (registry outputs and optional_outputs)"
            )));
        }
        if given.iter().any(|(k, _)| k == key) {
            return Err(fail(format!("--key {key} is given twice")));
        }
        if !value.is_empty() {
            given.push((key.to_owned(), absolute(&here.join(value))));
        }
    }
    if options.stub {
        for key in &outputs {
            if given.iter().any(|(k, _)| k == key) {
                continue;
            }
            let path = work.join("stub").join(key);
            if key.ends_with("_dir") || key == "loftee_repo" {
                std::fs::create_dir_all(&path)?;
            } else {
                std::fs::create_dir_all(work.join("stub"))?;
                write_file(&path, &format!("{id} {key}\n"))?;
            }
            given.push((key.clone(), path));
        }
    }

    // Each key as (KEY, the path in this task, the path it will have).
    let mut keys: Vec<(String, PathBuf, PathBuf)> = Vec::new();
    let mut missing = String::new();
    for key in &outputs {
        match given.iter().find(|(k, _)| k == key) {
            Some((_, path)) if std::fs::symlink_metadata(path).is_ok() => {
                keys.push((key.clone(), path.clone(), destination(path, &here, &folder)));
            }
            _ => {
                missing.push(' ');
                missing.push_str(key);
            }
        }
    }
    if !missing.is_empty() {
        return Err(fail(format!("outputs missing after the step:{missing}")));
    }
    // An optional output is recorded only when the step made it.
    for key in &optional {
        if let Some((_, path)) = given.iter().find(|(k, _)| k == key) {
            if std::fs::metadata(path).is_ok() && !keys.iter().any(|(k, _, _)| k == key) {
                keys.push((key.clone(), path.clone(), destination(path, &here, &folder)));
            }
        }
    }
    for (key, path, _) in &keys {
        if std::fs::metadata(path).is_err() {
            return Err(fail(format!("{key} ({}) does not exist", path.display())));
        }
    }

    // Outputs that lie in the source tree and outside the step's folder are
    // files shipped with the repository (install-task run's receipt_keys
    // reads them back).
    let tree = realpath_m(&root);
    let mut shipped = Vec::new();
    for (key, _, value) in &keys {
        if value.starts_with(&folder) {
            continue;
        }
        let value = realpath_m(value);
        if value.starts_with(&tree) && value != tree {
            let relative = value.strip_prefix(&tree).unwrap_or(&value);
            shipped.push((key.clone(), string(relative.to_string_lossy())));
        }
    }

    let mut upstream = Vec::new();
    let mut inputs = Vec::new();
    for input in options.all("input") {
        let receipt = read_json(Path::new(input))?;
        let mut fields = vec![
            (
                "id".into(),
                receipt.get("id").cloned().unwrap_or(Json::Null),
            ),
            (
                "lineage".into(),
                receipt.get("lineage").cloned().unwrap_or(Json::Null),
            ),
        ];
        if let Some(content) = receipt.get("content") {
            fields.push(("content".into(), content.clone()));
        }
        upstream.push(object(fields));
        if let Some(folder) = receipt.get("folder").and_then(Json::as_str) {
            inputs.push(content::Input {
                id: text(&receipt, "id"),
                folder: realpath_m(Path::new(folder)),
                content: receipt_content(&receipt),
            });
        }
    }

    // As install-task run: only a step that others need keys them by content.
    // The listing names this task's data folder as "own", as run names the
    // step's folder, so the same results give the same content either way.
    let needed = steps(&registry)?
        .iter()
        .any(|s| strings(s, "needs").iter().any(|n| *n == id));
    let content = if needed {
        let places = content::Places {
            own: realpath_m(&here),
            tree: tree.clone(),
            inputs: &inputs,
        };
        let pairs: Vec<(String, String)> = keys
            .iter()
            .map(|(k, path, _)| (k.clone(), path.to_string_lossy().into_owned()))
            .collect();
        let listing = content::listing(&work, &pairs, &places, threads())?;
        std::fs::write(here.join("content.txt"), &listing)?;
        Some(content::content(&listing))
    } else {
        None
    };
    if !truthy(entry.get("writable")) {
        files_read_only(&work);
    }
    let files: Vec<Json> = keys
        .iter()
        .map(|(_, path, value)| string(moved_identity(path, value)))
        .collect();
    let mut fields = vec![
        ("schema".into(), Json::Int(1)),
        ("id".into(), string(&id)),
        ("name".into(), string(text(entry, "name"))),
        ("lineage".into(), string(lineage)),
        ("folder".into(), string(folder.to_string_lossy())),
        (
            "keys".into(),
            object(
                keys.iter()
                    .map(|(k, _, v)| (k.clone(), string(v.to_string_lossy())))
                    .collect(),
            ),
        ),
        ("shipped".into(), object(shipped)),
    ];
    if let Some(content) = content {
        fields.push(("content".into(), string(content)));
    }
    fields.extend([
        ("files".into(), Json::Array(files)),
        ("inputs".into(), Json::Array(upstream)),
        ("finished".into(), string(now_iso())),
        ("host".into(), string(hostname())),
        ("stub".into(), Json::Bool(options.stub)),
    ]);
    write_file(&receipt_path, &pretty(&object(fields))?)
}

/// Rewrites the strings of a JSON file that name a path below `from` to the
/// same path below `to`, keeping the order of keys; a file written over
/// several lines is written back indented, a one-line file on one line.
fn relocate(file: &Path, from: &Path, to: &Path) -> Result<()> {
    let text = std::fs::read_to_string(file)
        .map_err(|e| fail(format!("--relocate {}: {e}", file.display())))?;
    let mut value = crate::json_text::parse(&text, false)
        .map_err(|e| fail(format!("--relocate {}: {e}", file.display())))?;
    let from_text = format!("{}/", from.display());
    let to_text = format!("{}/", to.display());
    let mut changed = 0usize;
    fn walk(value: &mut Json, from: &str, to: &str, changed: &mut usize) {
        match value {
            Json::Str(s) if s == from.trim_end_matches('/') => {
                *s = to.trim_end_matches('/').to_owned();
                *changed += 1;
            }
            Json::Str(s) if s.starts_with(from) => {
                *s = format!("{to}{}", &s[from.len()..]);
                *changed += 1;
            }
            Json::Array(items) => items.iter_mut().for_each(|v| walk(v, from, to, changed)),
            Json::Object(fields) => fields
                .iter_mut()
                .for_each(|(_, v)| walk(v, from, to, changed)),
            _ => {}
        }
    }
    walk(&mut value, &from_text, &to_text, &mut changed);
    if changed > 0 {
        let mut out = if text.trim_end().contains('\n') {
            pretty(&value)?
        } else {
            super::compact(&value)?
        };
        // The file keeps its last line as it was: with or without a newline.
        if text.ends_with('\n') && !out.ends_with('\n') {
            out.push('\n');
        }
        if !text.ends_with('\n') {
            out.truncate(out.trim_end_matches('\n').len());
        }
        let partial = PathBuf::from(format!("{}.relocate", file.display()));
        std::fs::write(&partial, out)?;
        std::fs::rename(&partial, file)?;
    }
    Ok(())
}

/// Where a key's file will be: a path under this task's `data/` moves into
/// the step's folder; a file staged into the task (a link to a file you
/// supply) is named by its own path; any other path stays as it is.
fn destination(path: &Path, here: &Path, folder: &Path) -> PathBuf {
    if let Ok(relative) = path.strip_prefix(here.join("data")) {
        return folder.join("data").join(relative);
    }
    if path.starts_with(here) {
        return realpath_m(path);
    }
    path.to_path_buf()
}

/// The identity `install-task run` records (`PATH SIZE MTIME`), with the path
/// the file will have: storeDir's `mv` keeps its size and modification time.
fn moved_identity(path: &Path, value: &Path) -> String {
    match std::fs::metadata(path) {
        Ok(meta) if meta.is_dir() => format!("{} folder {}", value.display(), meta.mtime()),
        Ok(meta) => format!("{} {} {}", value.display(), meta.size(), meta.mtime()),
        Err(_) => format!("missing {}", value.display()),
    }
}

/// Every file below `path` loses its write bits; folders and links are left.
fn files_read_only(path: &Path) {
    use std::os::unix::fs::PermissionsExt;
    let Ok(entries) = std::fs::read_dir(path) else {
        return;
    };
    for entry in entries.flatten() {
        let child = entry.path();
        let Ok(meta) = std::fs::symlink_metadata(&child) else {
            continue;
        };
        if meta.is_dir() {
            files_read_only(&child);
        } else if meta.is_file() {
            let mode = meta.permissions().mode() & !0o222;
            let _ = std::fs::set_permissions(&child, std::fs::Permissions::from_mode(mode));
        }
    }
}

fn threads() -> usize {
    std::env::var("INSTALL_THREADS")
        .ok()
        .and_then(|t| t.parse().ok())
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch(name: &str) -> PathBuf {
        let folder = std::env::temp_dir().join(format!("record-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&folder);
        std::fs::create_dir_all(&folder).unwrap();
        folder
    }

    #[test]
    fn relocate_rewrites_only_paths_below_the_task_data_folder() {
        let folder = scratch("relocate");
        let file = folder.join("ready.json");
        std::fs::write(
            &file,
            "{\n  \"path\": \"/w/1/data/A/t.gz\",\n  \"other\": \"/w/1/datafile\",\n  \"n\": [\"/w/1/data\", 3]\n}\n",
        )
        .unwrap();
        relocate(&file, Path::new("/w/1/data"), Path::new("/r/B3/L/data")).unwrap();
        assert_eq!(
            std::fs::read_to_string(&file).unwrap(),
            "{\n  \"path\": \"/r/B3/L/data/A/t.gz\",\n  \"other\": \"/w/1/datafile\",\n  \"n\": [\n    \"/r/B3/L/data\",\n    3\n  ]\n}\n"
        );
        let line = folder.join("line.json");
        std::fs::write(&line, "{\"output\": \"/w/1/data/ref.fa\"}").unwrap();
        relocate(&line, Path::new("/w/1/data"), Path::new("/r/A1/L/data")).unwrap();
        assert_eq!(
            std::fs::read_to_string(&line).unwrap(),
            "{\"output\":\"/r/A1/L/data/ref.fa\"}"
        );
        // json.dump(indent=2) writes no newline after the closing brace.
        let open = folder.join("manifest.json");
        std::fs::write(&open, "{\n  \"raw_dir\": \"/w/1/data/raw\"\n}").unwrap();
        relocate(&open, Path::new("/w/1/data"), Path::new("/r/F10/L/data")).unwrap();
        assert_eq!(
            std::fs::read_to_string(&open).unwrap(),
            "{\n  \"raw_dir\": \"/r/F10/L/data/raw\"\n}"
        );
        let _ = std::fs::remove_dir_all(&folder);
    }

    #[test]
    fn destination_moves_data_and_resolves_staged_links() {
        let here = scratch("destination");
        let supplied = here.join("yours.fa");
        std::fs::write(&supplied, ">chr1\nA\n").unwrap();
        std::fs::create_dir_all(here.join("staged")).unwrap();
        std::os::unix::fs::symlink(&supplied, here.join("staged/yours.fa")).unwrap();
        let folder = Path::new("/r/A1/L");
        assert_eq!(
            destination(&here.join("data/reference/x.fa"), &here, folder),
            PathBuf::from("/r/A1/L/data/reference/x.fa")
        );
        assert_eq!(
            destination(&here.join("staged/yours.fa"), &here, folder),
            realpath_m(&supplied)
        );
        assert_eq!(
            destination(Path::new("/elsewhere/f"), &here, folder),
            PathBuf::from("/elsewhere/f")
        );
        let _ = std::fs::remove_dir_all(&here);
    }
}
