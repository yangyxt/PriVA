//! `priva hpo-gene-ranking` (src/hpo_gene_ranking) against the Python scripts
//! of llm_gene_reranker: runs each command on the made-up inputs in
//! tests/hpo_gene_ranking/tests/fixtures and compares its files with what the
//! Python scripts wrote for the same inputs (tests/hpo_gene_ranking/tests/expected,
//! made by llm_gene_reranker's rust/tests/python_reference/make_expected.sh).
//! The folder mirrors llm_gene_reranker's rust/ folder, so the relative paths
//! the expected files record are unchanged.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Command;

const BINARY: &str = env!("CARGO_BIN_EXE_priva");

/// The folder the commands run in, as llm_gene_reranker's rust/ folder.
fn manifest_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/hpo_gene_ranking")
}

fn output_dir(name: &str) -> PathBuf {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR"))
        .join("hpo_gene_ranking_parity")
        .join(name);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create output folder");
    dir
}

fn run(arguments: &[&str]) -> std::process::Output {
    let output = Command::new(BINARY)
        .arg("hpo-gene-ranking")
        .args(arguments)
        .current_dir(manifest_dir())
        .output()
        .expect("run hpo-gene-ranking");
    assert!(
        output.status.success(),
        "hpo-gene-ranking {:?} failed:\n{}",
        arguments,
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn files_in(dir: &Path) -> BTreeSet<String> {
    std::fs::read_dir(dir)
        .unwrap_or_else(|error| panic!("{}: {error}", dir.display()))
        .map(|entry| {
            entry
                .expect("entry")
                .file_name()
                .to_string_lossy()
                .into_owned()
        })
        .collect()
}

/// Every file of `expected` must exist in `actual` with the same bytes, and no other file.
fn assert_same_files(expected: &Path, actual: &Path) {
    let expected_files = files_in(expected);
    let actual_files = files_in(actual);
    assert_eq!(
        expected_files,
        actual_files,
        "file names differ in {}",
        actual.display()
    );
    for name in &expected_files {
        let left = std::fs::read(expected.join(name)).expect("read expected");
        let right = std::fs::read(actual.join(name)).expect("read actual");
        if left != right {
            let left = String::from_utf8_lossy(&left);
            let right = String::from_utf8_lossy(&right);
            let line = left
                .lines()
                .zip(right.lines())
                .position(|(a, b)| a != b)
                .unwrap_or(left.lines().count().min(right.lines().count()));
            panic!(
                "{name} differs from the Python output at line {}:\npython: {:?}\nrust:   {:?}",
                line + 1,
                left.lines().nth(line),
                right.lines().nth(line)
            );
        }
    }
}

fn expected(path: &str) -> String {
    format!("tests/expected/{path}")
}

fn reference_arguments() -> Vec<String> {
    let pairs = [
        (
            "--nodes-tsv",
            expected("hpo_graph/hp_owl_hpo_only.nodes.tsv"),
        ),
        (
            "--edges-tsv",
            expected("hpo_graph/hp_owl_hpo_only.edges.tsv"),
        ),
        (
            "--phenotype-to-genes",
            "tests/fixtures/raw/phenotype_to_genes.txt".to_owned(),
        ),
        ("--go-to-genes", expected("prepared/go_to_genes.tsv")),
        (
            "--pr-to-gene-or-protein",
            expected("prepared/pr_to_gene_or_protein.tsv"),
        ),
        (
            "--mp-to-human-orthologs",
            expected("prepared/mp_to_human_orthologs.tsv"),
        ),
        (
            "--orpha-disease-to-genes",
            expected("prepared/orpha_disease_to_genes.tsv"),
        ),
        (
            "--uniprot-evidence",
            expected("prepared/uniprot_human_reviewed_evidence.tsv"),
        ),
        (
            "--human-gene-info",
            expected("prepared/human_gene_info.tsv"),
        ),
        (
            "--reactome-to-genes",
            expected("prepared/reactome_to_ncbi_genes.tsv"),
        ),
        (
            "--reactome-relations",
            "tests/fixtures/raw/ReactomePathwaysRelation.txt".to_owned(),
        ),
        ("--gene2go", "tests/fixtures/raw/gene2go".to_owned()),
        ("--go-obo", "tests/fixtures/raw/go-basic.obo".to_owned()),
    ];
    let mut arguments = vec![
        "rank-genes".to_owned(),
        "--graph-tool".to_owned(),
        String::new(),
    ];
    for (option, value) in pairs {
        arguments.push(option.to_owned());
        arguments.push(value);
    }
    arguments
}

/// Runs rank-genes as `make_expected.sh` ran the Python script for case `name`.
fn rank_case(name: &str, extra: &[&str]) {
    let out = output_dir(name);
    let mut arguments = reference_arguments();
    arguments.push("--output-dir".to_owned());
    arguments.push(out.to_string_lossy().into_owned());
    arguments.extend(extra.iter().map(|value| value.to_string()));
    let refs: Vec<&str> = arguments.iter().map(String::as_str).collect();
    let output = run(&refs);
    assert_same_files(&manifest_dir().join(expected(name)), &out);
    let stdout: String = String::from_utf8_lossy(&output.stdout)
        .lines()
        .filter(|line| {
            line.starts_with("external_immune_genes=") || line.starts_with("distance_reference=")
        })
        .map(|line| format!("{line}\n"))
        .collect();
    let python_stdout =
        std::fs::read_to_string(manifest_dir().join(expected(&format!("{name}.stdout.txt"))))
            .expect("read Python stdout");
    assert_eq!(
        python_stdout, stdout,
        "printed summary lines differ for {name}"
    );
}

#[test]
fn rank_all_samples() {
    rank_case(
        "rank_all",
        &[
            "--patient-hpo-tsv",
            "tests/fixtures/patient_hpos.tsv",
            "--all-samples",
        ],
    );
}

#[test]
fn rank_manual_hp_ids() {
    rank_case(
        "rank_hp_ids",
        &[
            "--hp-id",
            "HP:0002721",
            "--hp-id",
            "HP:0002090",
            "--hp-id",
            "HP:0000964",
            "--hp-id",
            "HP:0002721",
            "--hp-id",
            "HP:0000006",
            "--synthetic-sample-id",
            "MANUAL",
        ],
    );
}

#[test]
fn rank_with_other_options() {
    rank_case(
        "rank_options",
        &[
            "--patient-hpo-tsv",
            "tests/fixtures/patient_hpos.tsv",
            "--sample-id",
            "S2",
            "--sample-id",
            "S1",
            "--sample-id",
            "MISSING",
            "--distance-reference",
            "2.5",
            "--weight-model",
            "tier_sqrt",
            "--primary-max-gene-set",
            "2",
            "--weak-max-gene-set",
            "4",
            "--parent-propagation-depth-floor",
            "-1",
            "--include-secondary-always",
            "--no-external-immune",
        ],
    );
}

#[test]
fn rank_first_sample_by_default() {
    rank_case(
        "rank_first",
        &["--patient-hpo-tsv", "tests/fixtures/patient_hpos.tsv"],
    );
}

#[test]
fn rank_other_sources_in_parallel() {
    rank_case(
        "rank_sources",
        &[
            "--patient-hpo-tsv",
            "tests/fixtures/patient_hpos.tsv",
            "--direct-source",
            "rag_direct",
            "--direct-source",
            "propagated",
            "--all-samples",
            "--jobs",
            "3",
            "--primary-max-gene-set",
            "3",
            "--weak-max-gene-set",
            "5",
        ],
    );
}

#[test]
fn prepared_tables_match_python() {
    let out = output_dir("prepared");
    let out_text = out.to_string_lossy().into_owned();
    run(&[
        "prepare-gene-evidence",
        "--raw-dir",
        "tests/fixtures/raw",
        "--out-dir",
        &out_text,
    ]);
    // The manifest records the output folder; put the Python one in its place.
    let manifest = std::fs::read_to_string(out.join("manifest.json")).expect("manifest");
    std::fs::write(
        out.join("manifest.json"),
        manifest.replace(&out_text, "tests/expected/prepared"),
    )
    .expect("rewrite manifest");
    assert_same_files(&manifest_dir().join("tests/expected/prepared"), &out);
}

/// Reads a tab-separated table (with Python csv quoting) as rows of fields.
fn read_table(path: &Path) -> Vec<Vec<String>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .flexible(true)
        .from_path(path)
        .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    reader
        .records()
        .map(|record| record.expect("record").iter().map(str::to_owned).collect())
        .collect()
}

#[test]
// Fails in llm_gene_reranker too (master bb2df76): for HP:0001250, whose
// logical definition is a blank node (rdf:nodeID genid1) written after the
// class, the program follows the node and lists UBERON:0001016, while the
// stored Python table of the fixture lists no linked term. The stored table
// needs to be made again (or the program changed) in llm_gene_reranker first.
#[ignore = "stored Python table of the blank-node fixture disagrees; see the comment above"]
fn hpo_graph_tables_match_python() {
    let out = output_dir("hpo_graph");
    let prefix = out.join("hp_owl_hpo_only").to_string_lossy().into_owned();
    run(&[
        "build-hpo-graph",
        "--hpo-owl",
        "tests/fixtures/hp_fixture.owl",
        "--hpo-only",
        "--output-prefix",
        &prefix,
    ]);
    let python_dir = manifest_dir().join("tests/expected/hpo_graph");
    let edges_python =
        std::fs::read(python_dir.join("hp_owl_hpo_only.edges.tsv")).expect("python edges");
    let edges_rust = std::fs::read(out.join("hp_owl_hpo_only.edges.tsv")).expect("rust edges");
    assert!(
        edges_python == edges_rust,
        "edges table differs from the Python output"
    );
    // The Rust nodes table leaves out three owlready2-rendered columns; every
    // column it does write must equal the Python column of the same name.
    let python = read_table(&python_dir.join("hp_owl_hpo_only.nodes.tsv"));
    let rust = read_table(&out.join("hp_owl_hpo_only.nodes.tsv"));
    assert_eq!(python.len(), rust.len(), "number of node rows");
    let python_header = &python[0];
    for (rust_column, name) in rust[0].iter().enumerate() {
        let python_column = python_header
            .iter()
            .position(|column| column == name)
            .unwrap_or_else(|| panic!("Python nodes table has no column {name}"));
        for (row, (left, right)) in python.iter().zip(rust.iter()).enumerate().skip(1) {
            assert_eq!(
                left[python_column], right[rust_column],
                "nodes table row {row} column {name}"
            );
        }
    }
    let left_out: Vec<&String> = python_header
        .iter()
        .filter(|name| !rust[0].contains(name))
        .collect();
    assert_eq!(
        left_out,
        vec![
            "logical_definition_text",
            "annotation_property_count",
            "annotation_value_count"
        ],
        "columns the Rust table leaves out"
    );
    // Both tables end their lines with "\r\n", as Python's csv module does by default.
    let nodes = std::fs::read(out.join("hp_owl_hpo_only.nodes.tsv")).expect("rust nodes");
    assert!(nodes.ends_with(b"\r\n"));
}

#[test]
fn term_names_give_the_same_ranking_as_their_ids() {
    let by_name = output_dir("terms_by_name");
    let by_id = output_dir("terms_by_id");
    let mut arguments = reference_arguments();
    arguments.extend(
        [
            "--hpo-term-table",
            "tests/fixtures/hpo_terms.tsv",
            "--all-samples",
            "--output-dir",
            &by_name.to_string_lossy(),
        ]
        .iter()
        .map(|value| value.to_string()),
    );
    run(&arguments.iter().map(String::as_str).collect::<Vec<_>>());
    let mut arguments = reference_arguments();
    arguments.extend(
        [
            "--hp-id",
            "HP:0002721",
            "--hp-id",
            "HP:0002090",
            "--hp-id",
            "HP:0000964",
            "--synthetic-sample-id",
            "T1",
            "--distance-reference",
            "0",
            "--output-dir",
            &by_id.to_string_lossy(),
        ]
        .iter()
        .map(|value| value.to_string()),
    );
    run(&arguments.iter().map(String::as_str).collect::<Vec<_>>());
    // The term table holds T1 and T2, so its distance reference spans both;
    // T2 has a single frontier term, so it adds no pair and the reference is T1's.
    for table in [
        "gene_ranked",
        "gene_score_contributions",
        "p_pairs",
        "gene_summary",
    ] {
        let left =
            std::fs::read_to_string(by_name.join(format!("T1.{table}.tsv"))).expect("by name");
        let right = std::fs::read_to_string(by_id.join(format!("T1.{table}.tsv"))).expect("by id");
        assert_eq!(left, right, "T1.{table}.tsv differs between names and ids");
    }
    let frontier =
        std::fs::read_to_string(by_name.join("T1.frontier_p_terms.tsv")).expect("frontier");
    assert!(
        frontier.contains("eczematous DERMATITIS"),
        "the typed name is kept as the original entry"
    );
    assert!(frontier.contains("manual_hp_term"));
}

#[test]
fn unclear_term_names_are_refused() {
    let out = output_dir("unclear_terms");
    let mut arguments = reference_arguments();
    arguments.extend(
        [
            "--hp-term",
            "Dermatitis",
            "--hp-term",
            "No such phenotype",
            "--output-dir",
            &out.to_string_lossy(),
        ]
        .iter()
        .map(|value| value.to_string()),
    );
    let output = Command::new(BINARY)
        .arg("hpo-gene-ranking")
        .args(&arguments)
        .current_dir(manifest_dir())
        .output()
        .expect("run");
    assert!(!output.status.success());
    let message = String::from_utf8_lossy(&output.stderr);
    assert!(
        message.contains(
            "\"Dermatitis\" is the exact synonym of several terms: HP:0000964, HP:0011123"
        ),
        "{message}"
    );
    assert!(
        message.contains("\"No such phenotype\" is neither an HPO id"),
        "{message}"
    );
}
