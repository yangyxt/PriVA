//! The annotation hub folder name (lineage) of each annotation step: a hash
//! of the step's name, the code that step runs, and the identity of every
//! resource that step reads.
//!
//! The code of a step is only what that step runs: its own lines of the
//! annotation step wrapper, the shell functions those lines call and every
//! function they call in turn (in annotation_vcf.sh and common_bash_utils.sh),
//! the script-level variables those functions name, and the files under
//! scripts/ and data/ those functions name. Comment lines and blank lines do
//! not count. A change inside one step's function therefore gives a new
//! folder for that step only; a change to a shared function gives a new
//! folder for each step that calls it. Where those functions run a `priva`
//! command, the Rust files behind that command count too (rust_reach.rs).
//!
//! A resource's identity leaves its path out, so moving the checkout or the
//! caches keeps the folder: a file inside the PriVA checkout (bundled data)
//! counts by its content, any other by its size and recorded SHA-256, or its
//! size and modification time. With the install workflow's
//! resources.lock.json, a resource an install step produced counts by that
//! step's lineage instead, so a hub folder traces back through the install
//! steps to the downloads they came from. A change to any of them gives a
//! new, empty hub folder, so stored results are never reused with other
//! resources or code. The control cohort step has no folder ("off") when its
//! counts come from the cohort itself (no control VCF).
use super::Result;
use super::contract::{json_digest, object, text};
use crate::json_text::Json;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

/// The shell files whose functions the annotation steps run.
const SHELL_FILES: [&str; 2] = ["annotation_vcf.sh", "common_bash_utils.sh"];

pub(crate) const SOURCES: [&str; 5] = ["gnomad", "clinvar", "control", "vep", "cadd"];

/// The resources the VEP step reads: its cache, plugins and their data.
const VEP_KEYS: &[&str] = &[
    "ref_genome",
    "reference_fai",
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
    "gerp_bigwig",
    "loeuf_used_by_shell",
    "mavedb_file",
    "mavedb_assay_metadata_used_by_shell",
    "pext_mean_bw",
    "pext_tissue_bw_dir",
];

fn reads(source: &str, key: &str) -> bool {
    match source {
        "gnomad" => key == "gnomad_vcf_chrX" || key.starts_with("gnomad_chr"),
        "clinvar" => key == "clinvar_vcf",
        "control" => key == "control_vcf",
        "vep" => VEP_KEYS.contains(&key),
        "cadd" => key == "cadd_base_dir" || key == "ref_genome",
        _ => false,
    }
}

/// A resource record without its path: content digest for files inside the
/// checkout, else size with the recorded digest or the modification time.
fn identity(record: &Json, root: &Path) -> Result<Json> {
    let field = |name: &str| record.get(name).cloned().unwrap_or(Json::Null);
    let resolved = record
        .get("resolved_path")
        .and_then(Json::as_str)
        .map(Path::new);
    let canonical_root = root.canonicalize().unwrap_or_else(|_| root.to_path_buf());
    if let Some(path) = resolved.filter(|p| p.starts_with(&canonical_root) && p.is_file()) {
        return Ok(object(vec![(
            "sha256",
            text(crate::sha256::sha256_file(path)?),
        )]));
    }
    let digest = field("sha256");
    let marker = if matches!(digest, Json::Null) {
        field("mtime_ns")
    } else {
        digest
    };
    Ok(object(vec![
        ("type", field("type")),
        ("size", field("size")),
        ("marker", marker),
    ]))
}

/// The name a shell function definition line gives, for `function name {`,
/// `function name() {` and `name () {` at the start of a line.
fn function_name(line: &str) -> Option<&str> {
    let is_name = |c: char| c.is_ascii_alphanumeric() || c == '_';
    let (rest, keyword) = match line.strip_prefix("function ") {
        Some(rest) => (rest.trim_start(), true),
        None => (line, false),
    };
    let end = rest.find(|c: char| !is_name(c)).unwrap_or(rest.len());
    let (name, tail) = rest.split_at(end);
    if name.is_empty() || name.starts_with(|c: char| c.is_ascii_digit()) {
        return None;
    }
    let tail = tail.trim_start();
    let (tail, parentheses) = match tail.strip_prefix("()") {
        Some(tail) => (tail.trim_start(), true),
        None => (tail, false),
    };
    (tail.starts_with('{') && (keyword || parentheses)).then_some(name)
}

/// The variable a script-level assignment line sets (`NAME=...`,
/// `declare -a NAME=(`), if any.
fn assigned_name(line: &str) -> Option<&str> {
    let mut rest = line;
    for prefix in ["declare ", "export ", "readonly "] {
        if let Some(after) = rest.strip_prefix(prefix) {
            rest = after.trim_start();
            while rest.starts_with('-') {
                rest = rest.split_once(' ')?.1.trim_start();
            }
        }
    }
    let (name, _) = rest.split_once('=')?;
    let named = !name.is_empty()
        && !name.starts_with(|c: char| c.is_ascii_digit())
        && name.chars().all(|c| c.is_ascii_alphanumeric() || c == '_');
    named.then_some(name)
}

/// The names in a line: runs of letters, digits and `_`.
fn names(line: &str) -> impl Iterator<Item = &str> {
    line.split(|c: char| !(c.is_ascii_alphanumeric() || c == '_'))
        .filter(|word| !word.is_empty())
}

/// The names in a line that may be a function the line runs. In a message (a
/// log, echo or printf line) only a name where a command can start counts:
/// after `$(`, a backquote, `|`, `&` or `;`. Read whole, a message that
/// mentions a function would tie the step to it. The install steps' code
/// digest (`called_words` in scripts/install_task.sh) reads lines the same way.
fn run_names(line: &str) -> Vec<&str> {
    let trimmed = line.trim_start();
    let message = ["log", "echo", "printf"].iter().any(|command| {
        trimmed
            .strip_prefix(command)
            .is_some_and(|rest| rest.is_empty() || rest.starts_with(char::is_whitespace))
    });
    if !message {
        return names(line).collect();
    }
    let mut found = Vec::new();
    let mut rest = line;
    while let Some(at) = rest.find(['`', '|', '&', ';', '(']) {
        let opener = rest.as_bytes()[at];
        let after = &rest[at + 1..];
        rest = after;
        if opener == b'(' && !line[..line.len() - after.len() - 1].ends_with('$') {
            continue;
        }
        if let Some(name) = names(after.trim_start()).next() {
            if after.trim_start().starts_with(name) {
                found.push(name);
            }
        }
    }
    found
}

/// The path-like words of a line: runs of letters, digits and `_ . / -`.
fn path_words(line: &str) -> impl Iterator<Item = &str> {
    line.split(|c: char| !(c.is_ascii_alphanumeric() || "_./-".contains(c)))
        .filter(|word| !word.is_empty())
}

/// The shell code the annotation steps can run, by unit: every function, every
/// script-level variable, and the remaining script-level lines.
#[derive(Default)]
struct Shell {
    functions: BTreeMap<String, Vec<String>>,
    /// Script-level variable -> the lines that assign it.
    variables: BTreeMap<String, Vec<String>>,
    /// Script-level lines outside every function and assignment.
    script_level: Vec<String>,
}

impl Shell {
    fn read(root: &Path) -> Shell {
        let mut shell = Shell::default();
        for name in SHELL_FILES {
            if let Ok(content) = std::fs::read_to_string(root.join("scripts").join(name)) {
                shell.add(&content);
            }
        }
        shell
    }

    /// Sorts the lines of one shell file into units. Comment lines and blank
    /// lines are left out; a function ends at the first line starting with
    /// `}`, an array assignment at the first line starting with `)`. The
    /// block that runs only when the file is executed itself (not sourced,
    /// as every task does) is left out: `if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]`
    /// to its `fi`.
    fn add(&mut self, content: &str) {
        #[derive(PartialEq)]
        enum Inside {
            Nothing,
            Function,
            Array,
            DirectRun,
        }
        let mut inside = Inside::Nothing;
        let mut unit = String::new();
        for line in content.lines() {
            let line = line.trim_end();
            let counted = !line.trim_start().is_empty() && !line.trim_start().starts_with('#');
            match inside {
                Inside::DirectRun => {
                    if line == "fi" {
                        inside = Inside::Nothing;
                    }
                    continue;
                }
                Inside::Function | Inside::Array => {
                    let (units, close) = match inside {
                        Inside::Function => (&mut self.functions, '}'),
                        _ => (&mut self.variables, ')'),
                    };
                    if counted {
                        units.entry(unit.clone()).or_default().push(line.into());
                    }
                    if line.starts_with(close) {
                        inside = Inside::Nothing;
                    }
                    continue;
                }
                Inside::Nothing => {}
            }
            if !counted {
                continue;
            }
            if line.starts_with("if ") && line.contains("${#BASH_SOURCE[@]}\" -eq 1") {
                inside = Inside::DirectRun;
            } else if let Some(name) = function_name(line) {
                unit = name.to_owned();
                // A later definition of the same name replaces the earlier.
                self.functions.insert(unit.clone(), vec![line.into()]);
                inside = Inside::Function;
            } else if let Some(name) = assigned_name(line) {
                unit = name.to_owned();
                self.variables
                    .entry(unit.clone())
                    .or_default()
                    .push(line.into());
                if line.ends_with('(') {
                    inside = Inside::Array;
                }
            } else {
                self.script_level.push(line.into());
            }
        }
    }

    /// The functions and script-level variables a step runs: those its
    /// wrapper lines name, and those they name in turn.
    fn reached(&self, lines: &[&str]) -> (BTreeSet<&str>, BTreeSet<&str>) {
        let (mut functions, mut variables) = (BTreeSet::new(), BTreeSet::new());
        let mut pending: Vec<&str> = lines.to_vec();
        while let Some(line) = pending.pop() {
            for name in run_names(line) {
                if let Some((name, unit)) = self.functions.get_key_value(name) {
                    if functions.insert(name.as_str()) {
                        pending.extend(unit.iter().map(String::as_str));
                    }
                }
            }
            for name in names(line) {
                if let Some((name, unit)) = self.variables.get_key_value(name) {
                    if variables.insert(name.as_str()) {
                        pending.extend(unit.iter().map(String::as_str));
                    }
                }
            }
        }
        (functions, variables)
    }

    /// What one step's lineage records of its code.
    fn step_code(&self, source: &str, root: &Path) -> Result<Json> {
        let (wrapper, own) = wrapper_lines(source);
        let mut lines: Vec<&str> = wrapper.clone();
        lines.extend(&own);
        let (functions, variables) = self.reached(&lines);
        let digest = |lines: &[&str]| text(crate::sha256::sha256_hex(lines.join("\n").as_bytes()));
        // The `priva` commands those lines run, each with its Rust files.
        let mut called = BTreeSet::new();
        let function_lines = functions
            .iter()
            .flat_map(|name| self.functions[*name].iter().map(String::as_str));
        for line in lines.iter().copied().chain(function_lines) {
            called.extend(super::rust_reach::called_commands(line));
        }
        let table = super::rust_reach::commands(root);
        if !called.is_empty() && !root.join("rust/src/cli.rs").is_file() {
            return Err(format!(
                "{}: rust/src/cli.rs is needed to count the Rust code the {source} step runs",
                root.display()
            )
            .into());
        }
        let mut rust = Vec::new();
        for command in called {
            let Some(entry) = table.get(&command) else {
                continue;
            };
            let files = super::rust_reach::command_files(root, entry)?
                .into_iter()
                .map(|(file, digest)| (file, text(digest)))
                .collect();
            rust.push((command, Json::Object(files)));
        }
        let mut files = BTreeSet::new();
        let mut record = |units: &BTreeMap<String, Vec<String>>, reached: &BTreeSet<&str>| {
            let mut out = Vec::new();
            for name in reached {
                let unit: Vec<&str> = units[*name].iter().map(String::as_str).collect();
                for line in &unit {
                    named_files(line, root, &mut files);
                }
                out.push(((*name).to_owned(), digest(&unit)));
            }
            Json::Object(out)
        };
        let functions = record(&self.functions, &functions);
        let variables = record(&self.variables, &variables);
        imported_python(root, &mut files);
        let mut contents = Vec::new();
        for path in files {
            let content = crate::sha256::sha256_file(&root.join(&path))?;
            contents.push((path, text(content)));
        }
        let script_level: Vec<&str> = self.script_level.iter().map(String::as_str).collect();
        Ok(object(vec![
            ("wrapper", digest(&wrapper)),
            ("step", digest(&own)),
            ("script_level", digest(&script_level)),
            ("functions", functions),
            ("variables", variables),
            ("files", Json::Object(contents)),
            ("rust", Json::Object(rust)),
        ]))
    }
}

/// The lines of the annotation step wrapper every step runs (before its
/// `case`), and the lines of this step's own branch.
fn wrapper_lines(source: &str) -> (Vec<&'static str>, Vec<&'static str>) {
    let mut lines = super::tasks::ANNOTATION_STEP.lines();
    let wrapper: Vec<&str> = lines
        .by_ref()
        .take_while(|line| !line.starts_with("case "))
        .collect();
    let label = format!("{source})");
    let mut own = Vec::new();
    let mut found = false;
    for line in lines {
        if !found && line.trim() != label {
            continue;
        }
        found = true;
        own.push(line);
        if line.trim_end().ends_with(";;") {
            break;
        }
    }
    (wrapper, own)
}

/// Adds the files of the checkout a line names, as paths relative to the
/// checkout: a script file by its name (`NAME.py`, `.sh`, `.R`, `.pl`,
/// `.bash`) or as a path below a folder variable (`${SCRIPT_DIR}/NAME`), when
/// scripts/ holds it; and a shipped data file as `.../data/PATH`, when the
/// checkout holds it. The shared shell files count by function instead. A
/// path put together while the step runs is not found.
fn named_files(line: &str, root: &Path, files: &mut BTreeSet<String>) {
    for word in path_words(line) {
        let below = word.strip_prefix('/');
        let base = word.rsplit('/').next().unwrap_or(word);
        let script = matches!(
            base.rsplit_once('.'),
            Some((_, "py" | "sh" | "R" | "pl" | "bash"))
        );
        for name in [below.unwrap_or(""), if script { base } else { "" }] {
            let shared = SHELL_FILES.contains(&name) || name == "install_utils.sh";
            if !name.is_empty() && !shared && root.join("scripts").join(name).is_file() {
                files.insert(format!("scripts/{name}"));
            }
        }
        if let Some(at) = word.find("data/") {
            let path = &word[at..];
            if (at == 0 || word[..at].ends_with('/')) && root.join(path).is_file() {
                files.insert(path.to_owned());
            }
        }
    }
}

/// Adds the Python files of scripts/ that the named Python scripts import,
/// and those import.
fn imported_python(root: &Path, files: &mut BTreeSet<String>) {
    let mut pending: Vec<String> = files.iter().cloned().collect();
    while let Some(path) = pending.pop() {
        if !path.ends_with(".py") {
            continue;
        }
        let Ok(content) = std::fs::read_to_string(root.join(&path)) else {
            continue;
        };
        for line in content.lines() {
            let mut parts = line.split_whitespace();
            if !matches!(parts.next(), Some("from" | "import")) {
                continue;
            }
            let Some(module) = parts.next().and_then(|name| names(name).next()) else {
                continue;
            };
            let imported = format!("scripts/{module}.py");
            if root.join(&imported).is_file() && files.insert(imported.clone()) {
                pending.push(imported);
            }
        }
    }
}

/// What each of the five annotation steps hashes into its lineage; None for
/// a step without a hub folder.
fn records(
    manifest: &Json,
    root: &Path,
    lock: Option<&Json>,
) -> Result<Vec<(&'static str, Option<Json>)>> {
    let shell = Shell::read(root);
    let resources: Vec<(String, Json)> = match manifest.get("resources") {
        Some(Json::Object(items)) => items.clone(),
        _ => Vec::new(),
    };
    // Config key -> (install step, its lineage), from the lock's step keys.
    let mut produced = std::collections::HashMap::new();
    if let Some(Json::Object(steps)) = lock.and_then(|l| l.get("steps")) {
        for (step, entry) in steps {
            let lineage = entry
                .get("lineage")
                .and_then(Json::as_str)
                .unwrap_or_default()
                .to_owned();
            if let Some(Json::Object(keys)) = entry.get("keys") {
                for (key, _) in keys {
                    produced.insert(key.clone(), (step.clone(), lineage.clone()));
                }
            }
        }
    }
    let mut out = Vec::new();
    for source in SOURCES {
        let mut used = Vec::new();
        let mut from_install = std::collections::BTreeSet::new();
        for (key, record) in &resources {
            let top = key.split('/').next().unwrap_or(key);
            if !reads(source, top) {
                continue;
            }
            match produced.get(top) {
                // One record per resource: every file of a folder shares it.
                Some((step, lineage)) => {
                    if from_install.insert(top.to_owned()) {
                        let record = object(vec![
                            ("install_step", text(step.clone())),
                            ("lineage", text(lineage.clone())),
                        ]);
                        used.push((top.to_owned(), record));
                    }
                }
                None => used.push((key.clone(), identity(record, root)?)),
            }
        }
        let record = if source == "control" && used.is_empty() {
            None
        } else {
            Some(object(vec![
                ("source", text(source)),
                ("code", shell.step_code(source, root)?),
                ("resources", Json::Object(used)),
            ]))
        };
        out.push((source, record));
    }
    Ok(out)
}

fn lineage_of(record: &Option<Json>) -> Result<String> {
    Ok(match record {
        Some(record) => json_digest(record)?[..16].to_owned(),
        None => "off".to_owned(),
    })
}

/// {source: lineage} for the five annotation steps.
pub(crate) fn lineages(manifest: &Json, root: &Path, lock: Option<&Json>) -> Result<Json> {
    let mut out = Vec::new();
    for (source, record) in records(manifest, root, lock)? {
        out.push((source.to_owned(), text(lineage_of(&record)?)));
    }
    Ok(Json::Object(out))
}

/// {source: {lineage, code, resources}}: each step's lineage with the code
/// record it hashes and the names of the resources it hashes, for review.
pub(crate) fn explain(manifest: &Json, root: &Path, lock: Option<&Json>) -> Result<Json> {
    let mut out = Vec::new();
    for (source, record) in records(manifest, root, lock)? {
        let mut entry = vec![("lineage", text(lineage_of(&record)?))];
        if let Some(record) = &record {
            entry.push(("code", record.get("code").cloned().unwrap_or(Json::Null)));
            let names = match record.get("resources") {
                Some(Json::Object(items)) => items.iter().map(|(k, _)| text(k.clone())).collect(),
                _ => Vec::new(),
            };
            entry.push(("resources", Json::Array(names)));
        }
        out.push((source.to_owned(), object(entry)));
    }
    Ok(Json::Object(out))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn manifest(entries: &[(&str, i128)]) -> Json {
        Json::Object(vec![(
            "resources".into(),
            Json::Object(
                entries
                    .iter()
                    .map(|(k, size)| {
                        let record = object(vec![("size", Json::Int(*size))]);
                        ((*k).to_owned(), record)
                    })
                    .collect(),
            ),
        )])
    }

    fn lineage(json: &Json, source: &str) -> String {
        json.get(source).and_then(Json::as_str).unwrap().to_owned()
    }

    /// A checkout with the two shell files, for the code part of a lineage.
    fn checkout(name: &str, annotation: &str, common: &str) -> std::path::PathBuf {
        let root =
            std::env::temp_dir().join(format!("priva-lineage-{}-{name}", std::process::id()));
        let scripts = root.join("scripts");
        std::fs::create_dir_all(&scripts).unwrap();
        std::fs::create_dir_all(root.join("data")).unwrap();
        std::fs::write(scripts.join("annotation_vcf.sh"), annotation).unwrap();
        std::fs::write(scripts.join("common_bash_utils.sh"), common).unwrap();
        root
    }

    const COMMON: &str = "\
HELPER=\"$(dirname \"$0\")/helper.py\"
declare -a TISSUES=(
    \"Liver\"
)
function count_records {
    wc -l \"$1\"
}
function read_yaml () {
    cat \"$1\"
}
function shared_check {
    count_records checked
    log \"see main_workflow, or $(named_in_message)\"
}
function named_in_message {
    true
}
if [[ \"${#BASH_SOURCE[@]}\" -eq 1 ]]; then
    \"$@\"
fi
";

    fn annotation(clinvar: &str, cadd: &str) -> String {
        format!(
            "\
source common_bash_utils.sh
function anno_clinvar_data {{
    shared_check
    {clinvar}
}}
function anno_agg_gnomAD_data() {{
    count_records gnomad
}}
function Calculate_CADD {{
    python \"$HELPER\" {cadd}
    cut -f1 ${{BASE_DIR}}/data/map.tsv
}}
function anno_VEP_data {{
    echo \"${{TISSUES[@]}}\"
}}
function main_workflow {{
    anno_clinvar_data; anno_agg_gnomAD_data; Calculate_CADD; anno_VEP_data
}}
if [[ \"${{#BASH_SOURCE[@]}}\" -eq 1 ]]; then
    main_workflow \"$@\"
fi
"
        )
    }

    fn code_lineages(root: &Path) -> Json {
        lineages(&manifest(&[("control_vcf", 3)]), root, None).unwrap()
    }

    fn moved(before: &Json, after: &Json) -> Vec<&'static str> {
        SOURCES
            .into_iter()
            .filter(|source| lineage(before, source) != lineage(after, source))
            .collect()
    }

    #[test]
    fn a_change_inside_one_function_moves_only_the_steps_that_run_it() {
        let root = checkout("one", &annotation("echo a", "--one"), COMMON);
        let before = code_lineages(&root);
        // Calculate_CADD alone: only the CADD step.
        std::fs::write(
            root.join("scripts/annotation_vcf.sh"),
            annotation("echo a", "--two"),
        )
        .unwrap();
        assert_eq!(moved(&before, &code_lineages(&root)), ["cadd"]);
        // anno_clinvar_data alone: only the ClinVar step.
        std::fs::write(
            root.join("scripts/annotation_vcf.sh"),
            annotation("echo b", "--one"),
        )
        .unwrap();
        assert_eq!(moved(&before, &code_lineages(&root)), ["clinvar"]);
        // main_workflow and the direct-run block belong to no step.
        let other = annotation("echo a", "--one").replace(
            "anno_clinvar_data; anno_agg",
            "echo start; anno_clinvar_data; anno_agg",
        );
        std::fs::write(root.join("scripts/annotation_vcf.sh"), other).unwrap();
        assert!(moved(&before, &code_lineages(&root)).is_empty());
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn a_shared_function_moves_each_step_that_calls_it() {
        let root = checkout("shared", &annotation("echo a", "--one"), COMMON);
        let before = code_lineages(&root);
        let write = |common: String| {
            std::fs::write(root.join("scripts/common_bash_utils.sh"), common).unwrap();
            code_lineages(&root)
        };
        // shared_check is called by the ClinVar function only.
        let after = write(COMMON.replace("count_records checked", "count_records twice"));
        assert_eq!(moved(&before, &after), ["clinvar"]);
        // count_records is called by gnomAD directly and by ClinVar through
        // shared_check.
        let after = write(COMMON.replace("wc -l", "wc -c"));
        assert_eq!(moved(&before, &after), ["gnomad", "clinvar"]);
        // A function a message runs counts; one it only mentions does not.
        let after = write(COMMON.replace("    true", "    false"));
        assert_eq!(moved(&before, &after), ["clinvar"]);
        // read_yaml is called by the wrapper of every step.
        let after = write(COMMON.replace("cat \"$1\"", "cat -- \"$1\""));
        assert_eq!(moved(&before, &after).len(), 5);
        // A script-level variable counts for the steps that name it.
        let after = write(COMMON.replace("\"Liver\"", "\"Lung\""));
        assert_eq!(moved(&before, &after), ["vep"]);
        // Comment lines and blank lines do not count.
        let after = write(COMMON.replace(
            "    count_records checked",
            "    # why\n\n    count_records checked",
        ));
        assert!(moved(&before, &after).is_empty());
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn a_named_file_counts_for_the_steps_whose_code_names_it() {
        let root = checkout("files", &annotation("echo a", "--one"), COMMON);
        std::fs::write(root.join("scripts/helper.py"), "import shared_part\n").unwrap();
        std::fs::write(root.join("scripts/shared_part.py"), "X = 1\n").unwrap();
        std::fs::write(root.join("scripts/unused.py"), "X = 1\n").unwrap();
        std::fs::write(root.join("data/map.tsv"), "a\tb\n").unwrap();
        let mut before = code_lineages(&root);
        for (path, content) in [
            ("scripts/helper.py", "import shared_part\nY = 2\n"),
            ("scripts/shared_part.py", "X = 2\n"),
            ("data/map.tsv", "a\tc\n"),
        ] {
            std::fs::write(root.join(path), content).unwrap();
            let after = code_lineages(&root);
            assert_eq!(moved(&before, &after), ["cadd"], "{path}");
            before = after;
        }
        std::fs::write(root.join("scripts/unused.py"), "X = 2\n").unwrap();
        assert!(moved(&before, &code_lineages(&root)).is_empty());
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn a_rust_file_moves_only_the_steps_whose_code_runs_its_command() {
        let script = annotation("echo a", "--one").replace(
            "count_records gnomad",
            "\"${PRIVA_BIN}\" count-records gnomad",
        );
        let root = checkout("rust", &script, COMMON);
        let src = root.join("rust/src");
        std::fs::create_dir_all(&src).unwrap();
        let table = "pub(crate) fn run() {
    if first() == Some(\"count-records\") {
        return crate::count::dispatch();
    }
    if first() == Some(\"other\") {
        return crate::other::dispatch();
    }
}
";
        std::fs::write(src.join("cli.rs"), table).unwrap();
        std::fs::write(src.join("count.rs"), "pub(crate) fn dispatch() {}\n").unwrap();
        std::fs::write(src.join("other.rs"), "pub(crate) fn dispatch() {}\n").unwrap();
        let before = code_lineages(&root);
        // A command no step runs, then a test of the command gnomAD runs.
        std::fs::write(
            src.join("other.rs"),
            "pub(crate) fn dispatch() -> u8 { 1 }\n",
        )
        .unwrap();
        let tested = "pub(crate) fn dispatch() {}\n#[cfg(test)]\nmod tests {\n    fn t() {}\n}\n";
        std::fs::write(src.join("count.rs"), tested).unwrap();
        assert!(moved(&before, &code_lineages(&root)).is_empty());
        std::fs::write(
            src.join("count.rs"),
            "pub(crate) fn dispatch() -> u8 { 1 }\n",
        )
        .unwrap();
        assert_eq!(moved(&before, &code_lineages(&root)), ["gnomad"]);
        // Without the Rust source the folder cannot be named.
        std::fs::remove_dir_all(&src).unwrap();
        assert!(lineages(&manifest(&[("control_vcf", 3)]), &root, None).is_err());
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn shell_definition_lines_are_recognised() {
        for (line, name) in [
            ("function anno_x {", Some("anno_x")),
            ("function anno_x() {", Some("anno_x")),
            ("function normalize_vcf () {", Some("normalize_vcf")),
            ("anno_x () {", Some("anno_x")),
            ("anno_x() {", Some("anno_x")),
            ("anno_x {", None),
            ("if [[ -z $x ]]; then", None),
            ("function {", None),
        ] {
            assert_eq!(function_name(line), name, "{line}");
        }
        for (line, name) in [
            ("BASE_DIR=\"$(dirname x)\"", Some("BASE_DIR")),
            ("declare -a PEXT=(", Some("PEXT")),
            ("export -x NAME=1", Some("NAME")),
            ("if [[ -z $TMPDIR ]]; then TMPDIR=/tmp; fi", None),
            ("log \"a=b\"", None),
        ] {
            assert_eq!(assigned_name(line), name, "{line}");
        }
        for (line, expected) in [
            (
                "    anno_x \"$vcf\" | sort_it",
                vec!["anno_x", "vcf", "sort_it"],
            ),
            ("    log \"run main_install first\"", vec![]),
            (
                "    echo \"$(count_it x) and `other y`; third\"",
                vec!["count_it", "other", "third"],
            ),
            ("    logger (x)", vec!["logger", "x"]),
        ] {
            assert_eq!(run_names(line), expected, "{line}");
        }
    }

    #[test]
    fn a_resource_change_moves_only_the_steps_that_read_it() {
        let root = std::env::temp_dir();
        let of = |clinvar: i128, vep: i128| {
            let entries = [
                ("clinvar_vcf", clinvar),
                ("vep_cache_dir/a", vep),
                ("control_vcf", 3),
            ];
            lineages(&manifest(&entries), &root, None).unwrap()
        };
        let (before, after) = (of(1, 2), of(9, 2));
        assert_ne!(lineage(&before, "clinvar"), lineage(&after, "clinvar"));
        for source in ["gnomad", "control", "vep", "cadd"] {
            assert_eq!(
                lineage(&before, source),
                lineage(&after, source),
                "{source}"
            );
        }
        let vep_changed = of(1, 5);
        assert_ne!(lineage(&before, "vep"), lineage(&vep_changed, "vep"));
        assert_eq!(
            lineage(&before, "clinvar"),
            lineage(&vep_changed, "clinvar")
        );
    }

    #[test]
    fn moving_a_resource_keeps_its_lineage() {
        let record = |path: &str| {
            object(vec![
                ("path", text(path)),
                ("resolved_path", text(path)),
                ("size", Json::Int(5)),
                ("mtime_ns", Json::Int(7)),
            ])
        };
        let at = |path: &str| {
            let resources = Json::Object(vec![("clinvar_vcf".into(), record(path))]);
            let manifest = Json::Object(vec![("resources".into(), resources)]);
            lineages(&manifest, &std::env::temp_dir(), None).unwrap()
        };
        assert_eq!(
            lineage(&at("/old/clinvar.vcf.gz"), "clinvar"),
            lineage(&at("/new/c.vcf.gz"), "clinvar")
        );
    }

    #[test]
    fn with_a_lock_install_lineages_name_the_resources() {
        let entries = [
            ("vep_cache_dir/a", 1),
            ("vep_cache_dir/b", 2),
            ("clinvar_vcf", 3),
        ];
        let lock = |vep: &str| {
            let step = |lineage: &str, key: &str| {
                let keys = Json::Object(vec![(key.to_owned(), text("/x"))]);
                object(vec![("lineage", text(lineage)), ("keys", keys)])
            };
            let steps = vec![
                ("A2".to_owned(), step(vep, "vep_cache_dir")),
                ("D1".to_owned(), step("d1", "clinvar_vcf")),
            ];
            Json::Object(vec![("steps".into(), Json::Object(steps))])
        };
        let root = std::env::temp_dir();
        let (one, two) = (lock("v1"), lock("v2"));
        let first = lineages(&manifest(&entries), &root, Some(&one)).unwrap();
        let second = lineages(&manifest(&entries), &root, Some(&two)).unwrap();
        assert_ne!(lineage(&first, "vep"), lineage(&second, "vep"));
        assert_eq!(lineage(&first, "clinvar"), lineage(&second, "clinvar"));
        // The folder's file records no longer count once the lock names it.
        let resized = [
            ("vep_cache_dir/a", 9),
            ("vep_cache_dir/b", 2),
            ("clinvar_vcf", 3),
        ];
        let third = lineages(&manifest(&resized), &root, Some(&one)).unwrap();
        assert_eq!(lineage(&first, "vep"), lineage(&third, "vep"));
    }

    #[test]
    fn control_counts_from_the_cohort_have_no_hub_folder() {
        let json = lineages(
            &manifest(&[("clinvar_vcf", 1)]),
            &std::env::temp_dir(),
            None,
        )
        .unwrap();
        assert_eq!(lineage(&json, "control"), "off");
        assert_eq!(lineage(&json, "clinvar").len(), 16);
    }
}
