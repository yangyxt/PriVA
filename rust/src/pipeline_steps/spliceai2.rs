//! The SpliceAI2 step of a family's prioritization (optional): the SpliceAI2
//! model's scores for the insertions, deletions and multi-base substitutions
//! that the published SpliceAI2 tables do not hold, written into the combined
//! family table before the splicing step and the ACMG step read them
//! (`crate::spliceai2_scores`).
//!
//! The model is Illumina's `spliceai2` package, run by
//! scripts/spliceai2_score_variants.py with the Python of the setting
//! `spliceai2_python` and the weights in `spliceai2_model_dir`. Nobody can
//! install the weights for the user: they are given out only to a person who
//! has accepted the publisher's terms. Without `spliceai2_model_dir` this step
//! does not run.
//!
//! With an annotation hub, each variant is scored once: the hub folder
//! spliceai2/LINEAGE holds the scores as tables that are only ever added, and
//! LINEAGE counts what decides a score: the scoring script, the package's
//! version, the weights and the reference genome. New scores are stored after
//! every `CHUNK` variants, so a task that is stopped keeps what it has done.
use super::tools::{Context, Result, is_file, log, not_empty, random_tag, refuse};
use crate::spliceai2_scores::{SCORE_HEADER, VARIANT_HEADER};
use std::fs;
use std::os::unix::fs::MetadataExt;
use std::path::{Path, PathBuf};

/// The CPUs each model process works with.
const THREADS_PER_PROCESS: usize = 8;
/// The variants scored between two additions to the hub.
const CHUNK: usize = 100;

/// What the step did, for the task's record.
#[derive(Debug, PartialEq)]
pub(crate) struct Summary {
    /// Variants (one per strand) whose scores the published tables lack.
    pub(crate) needed: usize,
    /// Of them, those the annotation hub held.
    pub(crate) stored: usize,
    /// Those the model scored now.
    pub(crate) scored: usize,
    /// The table's rows that got scores.
    pub(crate) rows: usize,
    /// The hub folder of the scores.
    pub(crate) lineage: String,
}

/// The gene context the SpliceAI2 score files were prepared with, beside them
/// (`priva spliceai2-prepare`); `None` for a file that is not one of them.
pub(crate) fn gene_context(snv_prescore: &str) -> Option<String> {
    snv_prescore
        .strip_suffix(".snv.vep.vcf.gz")
        .map(|stem| format!("{stem}.gene_context.tsv"))
}

/// The weights: every model.ckpt one folder below the model folder, by name.
fn weights(model_dir: &Path) -> Result<Vec<PathBuf>> {
    let mut found: Vec<PathBuf> = fs::read_dir(model_dir)
        .map_err(|e| format!("{}: {e}", model_dir.display()))?
        .filter_map(|entry| entry.ok().map(|entry| entry.path().join("model.ckpt")))
        .filter(|path| path.is_file())
        .collect();
    found.sort();
    if found.is_empty() {
        return refuse(format!(
            "spliceai2_model_dir {} holds no FOLDER/model.ckpt",
            model_dir.display()
        ));
    }
    Ok(found)
}

/// The hub folder's name: a hash of what decides a score.
fn lineage(script: &Path, package: &str, weights: &[PathBuf], reference: &str) -> Result<String> {
    let genome = fs::metadata(reference).map_err(|e| format!("{reference}: {e}"))?;
    let mut counted = format!(
        "{SCORE_HEADER}\n{}\n{package}\n{} {}.{}\n",
        crate::sha256::sha256_file(script)?,
        genome.len(),
        genome.mtime(),
        genome.mtime_nsec()
    );
    for path in weights {
        counted.push_str(&format!("{}\n", crate::sha256::sha256_file(path)?));
    }
    Ok(crate::sha256::sha256_hex(counted.as_bytes())[..16].to_owned())
}

/// Adds a table of scores to the hub under a name made from its content. A hub
/// that cannot be written keeps the task going: its table has the scores.
fn store(scores: &Path, folder: &Path) {
    let add = || -> Result<()> {
        fs::create_dir_all(folder)?;
        let name = &crate::sha256::sha256_file(scores)?[..16];
        let target = folder.join(format!("scores-{name}.tsv"));
        if !target.is_file() {
            let staged = folder.join(format!("scores-{name}.{}.partial", random_tag()));
            fs::copy(scores, &staged)?;
            fs::rename(&staged, &target)?;
        }
        Ok(())
    };
    if let Err(error) = add() {
        log(&format!(
            "WARNING: the SpliceAI2 scores were not added to {}: {error}",
            folder.display()
        ));
    }
}

/// Fills the empty SpliceAI columns of the table's variants; `cpus` is what the
/// model may use and `hub` the annotation hub folder of the run, if any.
pub(crate) fn run(
    context: &Context,
    table: &str,
    hub: Option<&Path>,
    cpus: usize,
) -> Result<Summary> {
    let python = context.value("spliceai2_python");
    let model_dir = PathBuf::from(context.value("spliceai2_model_dir"));
    let reference = context.value("ref_genome");
    let Some(genes) = gene_context(&context.value("spliceai_snv_prescore")) else {
        return refuse(
            "spliceai2_model_dir needs SpliceAI2 score files: spliceai_snv_prescore is not one prepared by priva spliceai2-prepare",
        );
    };
    for (key, path) in [
        ("spliceai2_python", &python),
        ("ref_genome", &reference),
        ("the gene context of spliceai_snv_prescore", &genes),
    ] {
        if !is_file(path) {
            return refuse(format!("{key}: {path} is not a file"));
        }
    }
    let script = context.root.join("scripts/spliceai2_score_variants.py");
    let script_text = script.display().to_string();
    let weights = weights(&model_dir)?;
    let package = context.first_line(&python, &[&script_text, "--version"]);
    if !package.starts_with("spliceai2 ") {
        return refuse(format!(
            "spliceai2_python {python} cannot run the spliceai2 package ({script_text} --version)"
        ));
    }
    let mut summary = Summary {
        needed: 0,
        stored: 0,
        scored: 0,
        rows: 0,
        lineage: lineage(&script, &package, &weights, &reference)?,
    };
    if !not_empty(table) {
        return Ok(summary);
    }
    let folder = context
        .temp_dir()
        .join(format!("spliceai2.{}", random_tag()));
    fs::create_dir_all(&folder)?;
    let hub_folder = hub.map(|hub| hub.join("spliceai2").join(&summary.lineage));
    let (known, new) = (folder.join("known.tsv"), folder.join("new.tsv"));
    let (table_path, genes_path) = (Path::new(table), Path::new(&genes));
    let (needed, stored) = crate::spliceai2_scores::missing(
        table_path,
        genes_path,
        &new,
        hub_folder.as_deref().map(|hub| (hub, known.as_path())),
    )?;
    (summary.needed, summary.stored) = (needed, stored);
    log(&format!(
        "SpliceAI2: {needed} variants without a published score, {stored} of them in the annotation hub"
    ));
    let mut scores = Vec::new();
    if hub_folder.is_some() {
        scores.push(known);
    }
    let text = fs::read_to_string(&new)?;
    let variants: Vec<&str> = text.lines().skip(1).collect();
    let processes = (cpus / THREADS_PER_PROCESS).max(1);
    for (index, chunk) in variants.chunks(CHUNK).enumerate() {
        let list = folder.join(format!("variants.{index}.tsv"));
        let output = folder.join(format!("scores.{index}.tsv"));
        fs::write(&list, format!("{VARIANT_HEADER}\n{}\n", chunk.join("\n")))?;
        context.run(
            &python,
            &[
                script_text.clone(),
                "--variants".into(),
                list.display().to_string(),
                "--fasta".into(),
                reference.clone(),
                "--model-dir".into(),
                model_dir.display().to_string(),
                "--output".into(),
                output.display().to_string(),
                "--processes".into(),
                processes.to_string(),
                "--threads".into(),
                (cpus / processes).max(1).to_string(),
            ],
        )?;
        if !output.is_file() {
            return refuse(format!("{script_text} wrote no scores"));
        }
        summary.scored += fs::read_to_string(&output)?
            .lines()
            .count()
            .saturating_sub(1);
        if let Some(folder) = &hub_folder {
            store(&output, folder);
        }
        scores.push(output);
    }
    let unscored = needed - stored - summary.scored;
    if unscored > 0 {
        log(&format!(
            "SpliceAI2: {unscored} variants lie too close to a chromosome end for the model and stay without a score"
        ));
    }
    if !scores.is_empty() {
        let filled = folder.join("filled.tsv");
        summary.rows = crate::spliceai2_scores::fill(table_path, genes_path, &scores, &filled)?;
        // The temporary folder may be on another file system than the table.
        fs::copy(&filled, table)?;
    }
    let _ = fs::remove_dir_all(&folder);
    log(&format!(
        "SpliceAI2: {} variants scored, {} rows filled",
        summary.scored, summary.rows
    ));
    Ok(summary)
}

#[cfg(test)]
mod tests {
    use super::super::tools::Settings;
    use super::*;
    use std::os::unix::fs::PermissionsExt;

    #[test]
    fn the_gene_context_lies_beside_the_score_files() {
        assert_eq!(
            gene_context("/c/SpliceAI2_GRCh38.snv.vep.vcf.gz").as_deref(),
            Some("/c/SpliceAI2_GRCh38.gene_context.tsv")
        );
        assert_eq!(gene_context("/c/spliceai_scores.raw.snv.hg38.vcf.gz"), None);
    }

    /// A stand-in for the spliceai2 package's Python: it scores every variant
    /// of the list with the same values and notes how many it was given.
    const FAKE_PYTHON: &str = r#"#!/bin/bash
if [ "$2" = --version ]; then echo "spliceai2 2.0"; exit 0; fi
while [ $# -gt 0 ]; do
    case "$1" in --variants) list=$2;; --output) output=$2;; esac
    shift
done
printf 'chrom\tpos\tref\talt\tstrand\tDS_AG\tDS_AL\tDS_DG\tDS_DL\tDP_AG\tDP_AL\tDP_DG\tDP_DL\n' > "$output"
tail -n +2 "$list" | sed 's/$/\t0.00\t0.31\t0.02\t0.00\t1\t-2\t3\t4/' >> "$output"
wc -l < "$list" >> CALLS
"#;

    #[test]
    fn variants_are_scored_once_and_taken_from_the_hub_afterwards() {
        let root =
            std::env::temp_dir().join(format!("priva_spliceai2_step_{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(root.join("scripts")).unwrap();
        fs::create_dir_all(root.join("models/one")).unwrap();
        let calls = root.join("calls.txt");
        let write = |name: &str, text: &str| {
            let path = root.join(name);
            fs::write(&path, text).unwrap();
            path.display().to_string()
        };
        write("scripts/spliceai2_score_variants.py", "# stand-in\n");
        write("models/one/model.ckpt", "weights");
        let reference = write("reference.fa", ">chr1\nACGT\n");
        let python = write(
            "python",
            &FAKE_PYTHON.replace("CALLS", &calls.display().to_string()),
        );
        fs::set_permissions(&python, fs::Permissions::from_mode(0o755)).unwrap();
        write(
            "S.gene_context.tsv",
            "gene_id\tchrom\tgene_start\tgene_end\tstrand\tgene_name\nENSG01.1\tchr1\t1\t9999\t+\tONE\n",
        );
        let settings = Settings::parse(&format!(
            "spliceai2_python: {python}\nspliceai2_model_dir: {}\nref_genome: {reference}\nspliceai_snv_prescore: {}/S.snv.vep.vcf.gz\n",
            root.join("models").display(),
            root.display()
        ));
        let env = vec![("TMPDIR".to_owned(), root.display().to_string())];
        let config = root.join("task.yaml");
        let context = Context {
            root: &root,
            binary: "priva",
            config_path: &config,
            settings: &settings,
            env: &env,
        };
        let header = "chrom\tpos\tref\talt\tSYMBOL\tGene\tSpliceAI_pred_DP_AG\tSpliceAI_pred_DP_AL\tSpliceAI_pred_DP_DG\tSpliceAI_pred_DP_DL\tSpliceAI_pred_DS_AG\tSpliceAI_pred_DS_AL\tSpliceAI_pred_DS_DG\tSpliceAI_pred_DS_DL\tSpliceAI_pred_SYMBOL";
        // 150 deletions, so two runs of the model, and one single-base change.
        let mut text = format!("{header}\nchr1\t5\tA\tG\tONE\tENSG01\t\t\t\t\t\t\t\t\t\n");
        for pos in 100..250 {
            text.push_str(&format!(
                "chr1\t{pos}\tAC\tA\tONE\tENSG01\t\t\t\t\t\t\t\t\t\n"
            ));
        }
        let table = write("table.tsv", &text);
        let again = write("again.tsv", &text);
        let hub = root.join("hub");
        let first = run(&context, &table, Some(&hub), 24).unwrap();
        assert_eq!(
            (first.needed, first.stored, first.scored, first.rows),
            (150, 0, 150, 150)
        );
        // 100 and 50 variants, each list with its header line.
        assert_eq!(fs::read_to_string(&calls).unwrap(), "101\n51\n");
        let filled = fs::read_to_string(&table).unwrap();
        assert_eq!(
            filled
                .matches("\t1\t-2\t3\t4\t0.00\t0.31\t0.02\t0.00\tONE\n")
                .count(),
            150
        );
        assert!(filled.contains("chr1\t5\tA\tG\tONE\tENSG01\t\t\t\t\t\t\t\t\t\n"));
        let second = run(&context, &again, Some(&hub), 24).unwrap();
        assert_eq!(
            (second.needed, second.stored, second.scored, second.rows),
            (150, 150, 0, 150)
        );
        assert_eq!(second.lineage, first.lineage);
        assert_eq!(fs::read_to_string(&again).unwrap(), filled);
        assert_eq!(fs::read_to_string(&calls).unwrap(), "101\n51\n");
        // Other weights give another, empty hub folder; a table without
        // variants needs no model.
        write("models/one/model.ckpt", "other weights");
        let empty = write("empty.tsv", "");
        let none = run(&context, &empty, Some(&hub), 24).unwrap();
        assert_eq!((none.needed, none.rows), (0, 0));
        assert_ne!(none.lineage, first.lineage);
        fs::remove_dir_all(&root).unwrap();
    }
}
