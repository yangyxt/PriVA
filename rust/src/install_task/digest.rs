//! A step's code digest. Every step runs as its own Nextflow process
//! (registry `module`), so the digest is its registry entry, the text of its
//! module, the content of the scripts and shipped data files that text names,
//! and the Rust files behind every `priva` command it runs. Lists are in byte
//! order, so every machine computes the same digest whatever its language
//! setting.
//!
//! A step counts only the Rust code it runs: the files the command table hands
//! each of its `priva` commands to, and the files those name (the reach the
//! annotation hub lineage uses, `nextflow_task/rust_reach.rs`). A command
//! chosen at run time (`"$PRIVA_BIN" "$command"`) cannot be read from the
//! text, so a step whose module has such a line must name its commands in the
//! registry (`priva_commands`).
use super::{Options, Result, compact, fail, read_json, step, strings};
use crate::json_text::Json;
use crate::nextflow_task::rust_reach;
use crate::sha256::{Sha256, sha256_file};
use std::collections::{BTreeMap, BTreeSet, HashSet, VecDeque};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// What the steps' digests read: the checkout's scripts, data and Rust.
pub(crate) struct StepCode {
    scripts: PathBuf,
    root: PathBuf,
}

impl StepCode {
    pub(crate) fn load(root: &Path) -> StepCode {
        StepCode {
            scripts: root.join("scripts"),
            root: root.to_path_buf(),
        }
    }

    /// The named script files that exist in the scripts folder, and the
    /// Python files of that folder they import (and those import).
    fn with_imports(&self, named: Vec<String>) -> Vec<String> {
        let named = sorted(named, true);
        let mut queue: VecDeque<String> = named
            .into_iter()
            .filter(|w| self.scripts.join(w).is_file())
            .collect();
        let mut seen: HashSet<String> = HashSet::new();
        let mut found = Vec::new();
        while let Some(name) = queue.pop_front() {
            if !seen.insert(name.clone()) {
                continue;
            }
            found.push(name.clone());
            if !name.ends_with(".py") {
                continue;
            }
            let content = std::fs::read(self.scripts.join(&name)).unwrap_or_default();
            let modules: BTreeSet<String> = python_imports(&content).into_iter().collect();
            let modules = sorted(modules.into_iter().collect(), true);
            for module in modules {
                let file = format!("{module}.py");
                if self.scripts.join(&file).is_file() && !seen.contains(&file) {
                    queue.push_back(file);
                }
            }
        }
        sorted(found, true)
    }

    /// The digest of a step: its registry entry, the text of its Nextflow
    /// module, the content of the scripts and shipped data files that text
    /// names (`${projectDir}/scripts/NAME`, `${projectDir}/data/PATH`), and
    /// the Rust files of the `priva` commands it runs. The shared helpers of
    /// `modules/local/install_steps/common.nf` only write the adopt and
    /// record lines of `install-task`, the install runner, which no step's
    /// digest counts.
    pub(crate) fn digest(&self, entry: &Json) -> Result<String> {
        let module = process_module(entry).ok_or_else(|| {
            fail(format!(
                "step {} names no module (modules/local/install_steps/)",
                super::text(entry, "id")
            ))
        })?;
        let text = self.module_text(entry, module)?;
        let mut hash = Sha256::new();
        hash.update(compact(entry)?.as_bytes());
        hash.update(b"\n");
        hash.update(format!("== process {module}\n").as_bytes());
        hash.update(text.as_bytes());
        for script in self.process_scripts(&text) {
            let sum = sha256_file(&self.scripts.join(&script))?;
            hash.update(format!("{script} {sum}  -\n").as_bytes());
        }
        for data in self.process_data(&text) {
            let sum = sha256_file(&self.root.join("data").join(&data))?;
            hash.update(format!("data {data} {sum}  -\n").as_bytes());
        }
        let commands = self.process_commands(entry, &text)?;
        for command in &commands {
            hash.update(format!("priva {command}\n").as_bytes());
        }
        for (file, sum) in self.rust_files(&commands)? {
            hash.update(format!("rust {file} {sum}\n").as_bytes());
        }
        Ok(hash.finish_hex())
    }

    fn module_text(&self, entry: &Json, module: &str) -> Result<String> {
        let path = self.root.join(module);
        let text = std::fs::read_to_string(&path).map_err(|e| {
            fail(format!(
                "step {}: cannot read its module {}: {e}",
                super::text(entry, "id"),
                path.display()
            ))
        })?;
        // Steps that share a process name it (registry "process").
        let name = match entry.get("process").and_then(Json::as_str) {
            Some(name) => name.to_owned(),
            None => super::text(entry, "name"),
        };
        let process = format!("process {name} ");
        if !text.lines().any(|line| line.starts_with(&process)) {
            return Err(fail(format!(
                "step {}: {module} defines no process {name}",
                super::text(entry, "id"),
            )));
        }
        Ok(text)
    }

    /// The scripts of the checkout a module runs (`${projectDir}/scripts/NAME`)
    /// and the Python files they import.
    fn process_scripts(&self, text: &str) -> Vec<String> {
        self.with_imports(project_paths(text, "scripts"))
    }

    /// The shipped data files a module names: `${projectDir}/data/PATH`, with
    /// `${assembly}` read as hg19 and as hg38, tracked by git in a checkout.
    fn process_data(&self, text: &str) -> Vec<String> {
        let git = self.root.join(".git").exists();
        let mut found = Vec::new();
        for name in sorted(project_paths(text, "data"), true) {
            for assembly in ["hg19", "hg38"] {
                let path = name.replace("${assembly}", assembly);
                if path.contains('$') || !self.root.join("data").join(&path).is_file() {
                    continue;
                }
                if git && !git_tracked(&self.root, &format!("data/{path}")) {
                    continue;
                }
                found.push(path);
            }
        }
        sorted(found, true)
    }

    /// The `priva` commands a module runs (`"$PRIVA_BIN" COMMAND`) and those
    /// its registry entry lists, without `install-task`, the runner.
    pub(crate) fn process_commands(&self, entry: &Json, text: &str) -> Result<Vec<String>> {
        let table = rust_reach::commands(&self.root);
        let listed = strings(entry, "priva_commands");
        for command in &listed {
            if !table.contains_key(command) {
                return Err(fail(format!(
                    "step {}: priva_commands names {command}, which is not a priva command",
                    super::text(entry, "id")
                )));
            }
        }
        let mut found: BTreeSet<String> = listed.iter().cloned().collect();
        for line in text.lines() {
            if dynamic_command(line) && listed.is_empty() {
                return Err(fail(format!(
                    "step {}: its module picks a priva command at run time; list the commands in the registry's priva_commands",
                    super::text(entry, "id")
                )));
            }
            for command in called_commands(line) {
                if table.contains_key(&command) {
                    found.insert(command);
                }
            }
        }
        found.remove("install-task");
        Ok(found.into_iter().collect())
    }

    /// {file: SHA-256 of its counted lines} for the Rust files behind the
    /// commands.
    pub(crate) fn rust_files(&self, commands: &[String]) -> Result<BTreeMap<String, String>> {
        let table = rust_reach::commands(&self.root);
        let mut files = BTreeMap::new();
        for command in commands {
            if let Some(entry) = table.get(command) {
                files.extend(rust_reach::command_files(&self.root, entry)?);
            }
        }
        Ok(files)
    }
}

/// The Nextflow module of a step that runs as its own process.
pub(crate) fn process_module(entry: &Json) -> Option<&str> {
    entry.get("module").and_then(Json::as_str)
}

fn git_tracked(root: &Path, path: &str) -> bool {
    Command::new("git")
        .arg("-C")
        .arg(root)
        .args(["ls-files", "--error-unmatch", "--"])
        .arg(path)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

/// The paths a module names below `${projectDir}/FOLDER/`.
pub(crate) fn project_paths(text: &str, folder: &str) -> Vec<String> {
    let marker = format!("projectDir}}/{folder}/");
    let mut found = Vec::new();
    let mut rest = text;
    while let Some(at) = rest.find(&marker) {
        let after = &rest[at + marker.len()..];
        let end = after
            .find(|c: char| {
                !(c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '/' | '-' | '$' | '{' | '}'))
            })
            .unwrap_or(after.len());
        let name = after[..end].trim_end_matches(['.', '/']);
        if !name.is_empty() {
            found.push(name.to_owned());
        }
        rest = &after[end..];
    }
    found
}

/// The commands a line runs: `"${PRIVA_BIN}" COMMAND` and the forms
/// `rust_reach` reads, and `"$PRIVA_CACHE_INSTALL_BINARY" COMMAND`.
pub(super) fn called_commands(line: &str) -> Vec<String> {
    let mut found = rust_reach::called_commands(line);
    let mut rest = line;
    while let Some(at) = rest.find("PRIVA_CACHE_INSTALL_BINARY") {
        let after = rest[at + "PRIVA_CACHE_INSTALL_BINARY".len()..]
            .trim_start_matches(['}', '"', '\'', ' ']);
        let end = after
            .find(|c: char| !(c.is_ascii_lowercase() || c.is_ascii_digit() || c == '-'))
            .unwrap_or(after.len());
        if after.starts_with(|c: char| c.is_ascii_lowercase()) {
            found.push(after[..end].to_owned());
        }
        rest = &rest[at + 1..];
    }
    found
}

/// A line that runs `priva` with a command held in a variable.
pub(super) fn dynamic_command(line: &str) -> bool {
    ["PRIVA_BIN}\" \"$", "PRIVA_BIN\" \"$"]
        .iter()
        .any(|form| line.contains(form))
}

fn is_word_start(b: u8) -> bool {
    b.is_ascii_alphabetic() || b == b'_'
}

fn is_word(b: u8) -> bool {
    b.is_ascii_alphanumeric() || b == b'_'
}

fn is_space(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// The modules a Python file imports at the start of a line:
/// `^[[:space:]]*(from|import)[[:space:]]+NAME`.
fn python_imports(content: &[u8]) -> Vec<String> {
    let mut found = Vec::new();
    for line in content.split(|b| *b == b'\n') {
        let start = line.iter().take_while(|b| is_space(**b)).count();
        let rest = &line[start..];
        let Some(rest) = rest
            .strip_prefix(b"from")
            .or_else(|| rest.strip_prefix(b"import"))
        else {
            continue;
        };
        let spaces = rest.iter().take_while(|b| is_space(**b)).count();
        if spaces == 0 {
            continue;
        }
        let rest = &rest[spaces..];
        if !rest.first().is_some_and(|b| is_word_start(*b)) {
            continue;
        }
        let end = rest.iter().position(|b| !is_word(*b)).unwrap_or(rest.len());
        found.push(String::from_utf8_lossy(&rest[..end]).into_owned());
    }
    found
}

/// reach: what a step's code digest covers, to see why a step ran again.
pub(crate) fn reach(options: &Options) -> Result<()> {
    let registry = read_json(&options.path("registry")?)?;
    let id = options.one("step")?;
    let code = StepCode::load(&options.source_root()?);
    let entry = step(&registry, id)?;
    let module = process_module(entry).ok_or_else(|| {
        fail(format!(
            "step {id} names no module (modules/local/install_steps/)"
        ))
    })?;
    let text = code.module_text(entry, module)?;
    println!("module {module}");
    for script in code.process_scripts(&text) {
        println!("script {script}");
    }
    for data in code.process_data(&text) {
        println!("data {data}");
    }
    let commands = code.process_commands(entry, &text)?;
    for command in &commands {
        println!("priva {command}");
    }
    for file in code.rust_files(&commands)?.keys() {
        println!("rust {file}");
    }
    println!("digest {}", code.digest(entry)?);
    Ok(())
}

/// Byte order, once each with `unique`.
fn sorted(mut items: Vec<String>, unique: bool) -> Vec<String> {
    items.sort();
    if unique {
        items.dedup();
    }
    items
}
