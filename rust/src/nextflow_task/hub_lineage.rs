//! The annotation hub folder name (lineage) of each annotation step: a hash
//! of the step's name, the code that step runs, and the identity of every
//! resource that step reads.
//!
//! The code of a step is only what that step runs: its Rust file
//! (rust/src/pipeline_steps/STEP.rs) and every file of the crate that one
//! names in turn (rust_reach.rs), the `priva` commands those files run, each
//! with the Rust files behind it, and the files of the checkout's data/
//! folder they name. Comment lines, blank lines and test modules do not
//! count. A change to one step's file therefore gives a new folder for that
//! step only; a change to a file several steps share (pipeline_steps/tools.rs)
//! gives a new folder for each of them.
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

/// The string literals that follow `marker` in a text, with white space
/// between them allowed: `.priva("NAME"` gives NAME.
fn literals_after<'a>(text: &'a str, marker: &str) -> Vec<&'a str> {
    let mut found = Vec::new();
    let mut from = 0;
    while let Some(offset) = text[from..].find(marker) {
        from += offset + marker.len();
        let Some(rest) = text[from..].trim_start().strip_prefix('"') else {
            continue;
        };
        if let Some(end) = rest.find('"') {
            found.push(&rest[..end]);
        }
    }
    found
}

/// The `priva` commands a step's Rust text runs: `context.priva("NAME", ..)`,
/// `context.priva_succeeds("NAME", ..)`, and the command a tool is handed the
/// program for (`"PRIVA_CADD_PRESCORE_BIN"` is `cadd-prescore`).
fn run_commands(text: &str) -> BTreeSet<String> {
    let mut found: BTreeSet<String> = [".priva(", ".priva_succeeds("]
        .into_iter()
        .flat_map(|marker| literals_after(text, marker))
        .map(str::to_owned)
        .collect();
    let mut rest = text;
    while let Some(at) = rest.find("\"PRIVA_") {
        rest = &rest[at + "\"PRIVA_".len()..];
        let Some(end) = rest.find('"') else {
            break;
        };
        if let Some(tool) = rest[..end].strip_suffix("_BIN") {
            if tool.chars().all(|c| c.is_ascii_uppercase() || c == '_') {
                found.insert(tool.to_ascii_lowercase().replace('_', "-"));
            }
        }
    }
    found
}

/// The files of data/ a step's Rust text names: `context.data("PATH")`.
fn data_files(text: &str) -> BTreeSet<String> {
    literals_after(text, ".data(")
        .into_iter()
        .map(|path| format!("data/{path}"))
        .collect()
}

/// What one step's lineage records of its code.
fn step_code(source: &str, root: &Path, table: &BTreeMap<String, String>) -> Result<Json> {
    let entry = format!("pipeline_steps/{source}.rs");
    if !root.join("rust/src").join(&entry).is_file() {
        return Err(format!(
            "{}: rust/src/{entry} is needed to count the code the {source} step runs",
            root.display()
        )
        .into());
    }
    let texts = super::rust_reach::reached_texts(root, &entry)?;
    let (mut commands, mut data) = (BTreeSet::new(), BTreeSet::new());
    for text in texts.values() {
        commands.extend(run_commands(text));
        data.extend(data_files(text));
    }
    let digests = |files: BTreeMap<String, String>| {
        Json::Object(files.into_iter().map(|(f, d)| (f, text(d))).collect())
    };
    let step = digests(super::rust_reach::command_files(root, &entry)?);
    if !commands.is_empty() && table.is_empty() {
        return Err(format!(
            "{}: rust/src/cli.rs is needed to count the Rust code the {source} step runs",
            root.display()
        )
        .into());
    }
    let mut rust = Vec::new();
    for command in commands {
        let file = table.get(&command).ok_or_else(|| {
            format!("the {source} step runs `priva {command}`, which rust/src/cli.rs does not name")
        })?;
        rust.push((
            command,
            digests(super::rust_reach::command_files(root, file)?),
        ));
    }
    let mut files = Vec::new();
    for path in data {
        let file = root.join(&path);
        let content = if file.is_file() {
            text(crate::sha256::sha256_file(&file)?)
        } else {
            Json::Null
        };
        files.push((path, content));
    }
    Ok(object(vec![
        ("step", step),
        ("commands", Json::Object(rust)),
        ("data", Json::Object(files)),
    ]))
}

/// What each of the five annotation steps hashes into its lineage; None for
/// a step without a hub folder.
fn records(
    manifest: &Json,
    root: &Path,
    lock: Option<&Json>,
) -> Result<Vec<(&'static str, Option<Json>)>> {
    let table = super::rust_reach::commands(root);
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
                ("code", step_code(source, root, &table)?),
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

    const TABLE: &str = "pub(crate) fn run() {
    if first() == Some(\"spliceai-cutoffs\") {
        return crate::splice::dispatch();
    }
    if first() == Some(\"cadd-prescore\") {
        return crate::prescore::dispatch();
    }
    if first() == Some(\"other\") {
        return crate::other::dispatch();
    }
}
";

    const PLAIN_STEP: &str = "use super::tools::shared;\npub(crate) fn run() {\n    shared()\n}\n";

    const VEP_STEP: &str = "use super::annotation::assembly;
use super::tools::shared;
pub(crate) fn run(context: &Context) -> bool {
    assembly();
    shared();
    context.priva_succeeds(\"spliceai-cutoffs\", &[])
}
";

    const CADD_STEP: &str = "use super::annotation::assembly;
use super::tools::shared;
pub(crate) fn run(context: &Context) {
    assembly();
    shared();
    let map = context.data(\"liftover/map.tsv\");
    command.env(\"PRIVA_CADD_PRESCORE_BIN\", context.binary);
}
";

    /// A checkout with a command table, the five step files and the files
    /// they share, for the code part of a lineage.
    fn checkout(name: &str) -> std::path::PathBuf {
        let root = std::env::temp_dir().join(format!(
            "priva-lineage-{}-{name}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        let _ = std::fs::remove_dir_all(&root);
        let src = root.join("rust/src");
        let steps = src.join("pipeline_steps");
        std::fs::create_dir_all(&steps).unwrap();
        std::fs::create_dir_all(root.join("data/liftover")).unwrap();
        for (path, content) in [
            ("rust/src/cli.rs", TABLE),
            ("rust/src/splice.rs", "pub(crate) fn dispatch() {}\n"),
            ("rust/src/prescore.rs", "pub(crate) fn dispatch() {}\n"),
            ("rust/src/other.rs", "pub(crate) fn dispatch() {}\n"),
            (
                "rust/src/pipeline_steps/mod.rs",
                "pub(crate) mod annotation;\npub(crate) mod tools;\n",
            ),
            (
                "rust/src/pipeline_steps/tools.rs",
                "pub(crate) fn shared() {}\n",
            ),
            (
                "rust/src/pipeline_steps/annotation.rs",
                "pub(crate) fn assembly() {}\n",
            ),
            ("rust/src/pipeline_steps/gnomad.rs", PLAIN_STEP),
            ("rust/src/pipeline_steps/clinvar.rs", PLAIN_STEP),
            ("rust/src/pipeline_steps/control.rs", PLAIN_STEP),
            ("rust/src/pipeline_steps/vep.rs", VEP_STEP),
            ("rust/src/pipeline_steps/cadd.rs", CADD_STEP),
            ("data/liftover/map.tsv", "1\tchr1\n"),
            ("data/liftover/unused.tsv", "1\tchr1\n"),
        ] {
            std::fs::write(root.join(path), content).unwrap();
        }
        root
    }

    fn code_lineages(root: &Path) -> Json {
        lineages(&manifest(&[("control_vcf", 3)]), root, None).unwrap()
    }

    /// The steps whose folder moves when FILE of the checkout becomes
    /// CONTENT; the file is put back afterwards.
    fn moved_by(root: &Path, file: &str, content: &str) -> Vec<&'static str> {
        let before = code_lineages(root);
        let path = root.join(file);
        let old = std::fs::read_to_string(&path).unwrap();
        std::fs::write(&path, content).unwrap();
        let after = code_lineages(root);
        std::fs::write(&path, old).unwrap();
        SOURCES
            .into_iter()
            .filter(|source| lineage(&before, source) != lineage(&after, source))
            .collect()
    }

    #[test]
    fn a_change_to_one_step_file_moves_only_that_step() {
        let root = checkout("one");
        let step = "rust/src/pipeline_steps/clinvar.rs";
        let changed = PLAIN_STEP.replace("shared()", "shared();\n    shared()");
        assert_eq!(moved_by(&root, step, &changed), ["clinvar"]);
        let cadd = CADD_STEP.replace("assembly();", "assembly();\n    assembly();");
        assert_eq!(
            moved_by(&root, "rust/src/pipeline_steps/cadd.rs", &cadd),
            ["cadd"]
        );
        // Comment lines, blank lines and tests do not count.
        let quiet =
            format!("// why\n\n{PLAIN_STEP}#[cfg(test)]\nmod tests {{\n    fn t() {{}}\n}}\n");
        assert!(moved_by(&root, step, &quiet).is_empty());
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn a_shared_file_moves_each_step_that_names_it() {
        let root = checkout("shared");
        let changed = "pub(crate) fn shared() -> u8 { 1 }\n";
        let all = moved_by(&root, "rust/src/pipeline_steps/tools.rs", changed);
        assert_eq!(all, SOURCES);
        let changed = "pub(crate) fn assembly() -> u8 { 1 }\n";
        let named = moved_by(&root, "rust/src/pipeline_steps/annotation.rs", changed);
        assert_eq!(named, ["vep", "cadd"]);
        // The module list names no step's code.
        let listed = "pub(crate) mod annotation;\npub(crate) mod cadd;\npub(crate) mod tools;\n";
        assert!(moved_by(&root, "rust/src/pipeline_steps/mod.rs", listed).is_empty());
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn a_data_file_counts_for_the_steps_that_name_it() {
        let root = checkout("data");
        assert_eq!(
            moved_by(&root, "data/liftover/map.tsv", "1\tchr2\n"),
            ["cadd"]
        );
        assert!(moved_by(&root, "data/liftover/unused.tsv", "1\tchr2\n").is_empty());
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn a_command_file_moves_only_the_steps_that_run_it() {
        let root = checkout("commands");
        let changed = "pub(crate) fn dispatch() -> u8 { 1 }\n";
        assert_eq!(moved_by(&root, "rust/src/splice.rs", changed), ["vep"]);
        assert_eq!(moved_by(&root, "rust/src/prescore.rs", changed), ["cadd"]);
        assert!(moved_by(&root, "rust/src/other.rs", changed).is_empty());
        // A command the table does not name, or no table, names no folder.
        let unknown = VEP_STEP.replace("spliceai-cutoffs", "spliceai-guess");
        std::fs::write(root.join("rust/src/pipeline_steps/vep.rs"), unknown).unwrap();
        assert!(lineages(&manifest(&[]), &root, None).is_err());
        std::fs::write(root.join("rust/src/pipeline_steps/vep.rs"), VEP_STEP).unwrap();
        std::fs::remove_file(root.join("rust/src/cli.rs")).unwrap();
        assert!(lineages(&manifest(&[]), &root, None).is_err());
        // Without a step's file the folder cannot be named.
        std::fs::remove_dir_all(root.join("rust/src")).unwrap();
        assert!(lineages(&manifest(&[]), &root, None).is_err());
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn commands_and_data_files_are_read_from_the_step_text() {
        let text = "context.priva(\n        \"variant-ids\",\n        args(&[]),\n    )?;
let ok = context.priva_succeeds(\"spliceai-cutoffs\", &[]);
context.priva(command, args);
command.env(\"PRIVA_CADD_PRESCORE_BIN\", binary).env(\"PRIVA_DATA_DIR\", data);
let map = context.data(\"liftover/map.tsv\");
let named = context.data(&format!(\"loeuf/{name}\"));";
        let commands: Vec<String> = run_commands(text).into_iter().collect();
        assert_eq!(
            commands,
            ["cadd-prescore", "spliceai-cutoffs", "variant-ids"]
        );
        let data: Vec<String> = data_files(text).into_iter().collect();
        assert_eq!(data, ["data/liftover/map.tsv"]);
    }

    #[test]
    fn each_step_of_this_checkout_has_a_lineage() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let explained = explain(&manifest(&[("control_vcf", 3)]), &root, None).unwrap();
        let commands = |source: &str| -> Vec<String> {
            match explained
                .get(source)
                .and_then(|e| e.get("code"))
                .and_then(|c| c.get("commands"))
            {
                Some(Json::Object(items)) => items.iter().map(|(k, _)| k.clone()).collect(),
                _ => Vec::new(),
            }
        };
        assert_eq!(commands("vep"), ["spliceai-cutoffs"]);
        assert_eq!(commands("cadd"), ["cadd-prescore"]);
        for source in ["gnomad", "clinvar", "control"] {
            assert!(commands(source).is_empty(), "{source}");
        }
    }

    #[test]
    fn a_resource_change_moves_only_the_steps_that_read_it() {
        let root = checkout("resources");
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
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn moving_a_resource_keeps_its_lineage() {
        let root = checkout("moving");
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
            lineages(&manifest, &root, None).unwrap()
        };
        assert_eq!(
            lineage(&at("/old/clinvar.vcf.gz"), "clinvar"),
            lineage(&at("/new/c.vcf.gz"), "clinvar")
        );
        std::fs::remove_dir_all(root).unwrap();
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
        let root = checkout("lock");
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
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn control_counts_from_the_cohort_have_no_hub_folder() {
        let root = checkout("control");
        let json = lineages(&manifest(&[("clinvar_vcf", 1)]), &root, None).unwrap();
        assert_eq!(lineage(&json, "control"), "off");
        assert_eq!(lineage(&json, "clinvar").len(), 16);
        std::fs::remove_dir_all(root).unwrap();
    }
}
