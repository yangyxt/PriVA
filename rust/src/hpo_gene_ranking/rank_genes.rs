//! `rank-genes`: the command-line part of `patient_frontier_gene_evidence.py`
//! (inputs, sample choice, distance reference, files), plus HPO term names as
//! input.

use crate::hpo_gene_ranking::evidence::{EvidencePaths, EvidenceStore};
use crate::hpo_gene_ranking::immune::external_immune_genes;
use crate::hpo_gene_ranking::ontology::Ontology;
use crate::hpo_gene_ranking::python_text::{general, path_join, path_text, strip, with_suffix};
use crate::hpo_gene_ranking::ranking::{
    self, DirectTerm, FRONTIER_FIELDNAMES, GENE_EVIDENCE_FIELDNAMES, GENE_RANK_FIELDNAMES,
    GENE_SCORE_CONTRIBUTION_FIELDNAMES, GENE_SUMMARY_FIELDNAMES, PAIR_FIELDNAMES,
    RUN_SUMMARY_FIELDNAMES, Settings,
};
use crate::hpo_gene_ranking::tsv::{Record, Table, write_records};
use crate::hpo_gene_ranking::{Options, Result};
use indexmap::IndexMap;
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::sync::atomic::{AtomicUsize, Ordering};

const HELP: &str = "\
hpo-gene-ranking rank-genes: rank genes for each sample's HPO terms

Input terms (one of):
  --hp-id HP:0000001            an HPO id of one sample; repeat for more
  --hp-term TEXT                an HPO id or an exact HPO term name or exact
                                synonym (any letter case); repeat for more
  --synthetic-sample-id NAME    sample name for --hp-id/--hp-term (SYNTHETIC_HP_TERMS)
  --patient-hpo-tsv PATH        table with sample_id, hp_id, hp_name, source
                                [, original_entry] (the HPO extraction output)
  --hpo-term-table PATH         table with sample_id and hpo_term (an id or name)
  --sample-id ID                sample of a table to rank; repeat for more
  --all-samples                 rank every sample of the table
                                (without --sample-id or --all-samples: the first sample)
  --direct-source NAME          source values of --patient-hpo-tsv to use; repeat
                                (rag_direct, embedded_id, text_match)

Reference tables (each defaults to its place under --reference-dir):
  --reference-dir DIR           folder with raw/, prepared/ and hpo_graph/
  --nodes-tsv --edges-tsv       hpo_graph/hp_owl_hpo_only.nodes.tsv, .edges.tsv
  --phenotype-to-genes          raw/phenotype_to_genes.txt
  --go-to-genes --pr-to-gene-or-protein --mp-to-human-orthologs
  --orpha-disease-to-genes --uniprot-evidence --human-gene-info
  --reactome-to-genes           prepared/<the Python default names>
  --reactome-relations --gene2go --go-obo
                                raw/ReactomePathwaysRelation.txt, raw/gene2go, raw/go-basic.obo
  --no-external-immune          do not use Reactome and GO immune gene sets

Scoring (as the Python script):
  --primary-max-gene-set N (100) --weak-max-gene-set N (300)
  --weight-model inverse_gene_set|tier_sqrt   --distance-reference X (0: computed)
  --parent-propagation-depth-floor N (3; -1 turns it off) --include-secondary-always

Output:
  --output-dir DIR              per-sample tables and run_summary.tsv
  --jobs N                      samples ranked at the same time (1)";

const VALUED: &[&str] = &[
    "patient-hpo-tsv",
    "hpo-term-table",
    "graph-tool",
    "nodes-tsv",
    "edges-tsv",
    "phenotype-to-genes",
    "go-to-genes",
    "pr-to-gene-or-protein",
    "mp-to-human-orthologs",
    "orpha-disease-to-genes",
    "uniprot-evidence",
    "human-gene-info",
    "reactome-to-genes",
    "reactome-relations",
    "gene2go",
    "go-obo",
    "reference-dir",
    "output-dir",
    "sample-id",
    "jobs",
    "hp-id",
    "hp-term",
    "synthetic-sample-id",
    "random-samples",
    "random-seed",
    "direct-source",
    "primary-max-gene-set",
    "weak-max-gene-set",
    "weight-model",
    "distance-bonus-scale",
    "cross-module-bonus-scale",
    "distance-reference",
    "parent-propagation-depth-floor",
];
const FLAGS: &[&str] = &[
    "no-external-immune",
    "all-samples",
    "include-secondary-always",
];

/// Where each reference table sits under `--reference-dir`.
const REFERENCE_LAYOUT: [(&str, &str); 13] = [
    ("nodes-tsv", "hpo_graph/hp_owl_hpo_only.nodes.tsv"),
    ("edges-tsv", "hpo_graph/hp_owl_hpo_only.edges.tsv"),
    ("phenotype-to-genes", "raw/phenotype_to_genes.txt"),
    ("go-to-genes", "prepared/go_to_genes.tsv"),
    (
        "pr-to-gene-or-protein",
        "prepared/pr_to_gene_or_protein.tsv",
    ),
    (
        "mp-to-human-orthologs",
        "prepared/mp_to_human_orthologs.tsv",
    ),
    (
        "orpha-disease-to-genes",
        "prepared/orpha_disease_to_genes.tsv",
    ),
    (
        "uniprot-evidence",
        "prepared/uniprot_human_reviewed_evidence.tsv",
    ),
    ("human-gene-info", "prepared/human_gene_info.tsv"),
    ("reactome-to-genes", "prepared/reactome_to_ncbi_genes.tsv"),
    ("reactome-relations", "raw/ReactomePathwaysRelation.txt"),
    ("gene2go", "raw/gene2go"),
    ("go-obo", "raw/go-basic.obo"),
];

fn reference_path(options: &Options, name: &str) -> Result<String> {
    if let Some(path) = options.value(name) {
        return Ok(path.to_owned());
    }
    let relative = REFERENCE_LAYOUT
        .iter()
        .find(|(option, _)| *option == name)
        .map(|(_, relative)| *relative)
        .expect("every reference option has a place");
    match options.value("reference-dir") {
        Some(directory) => {
            let path = path_join(directory, relative);
            // A reference folder may hold gene2go only compressed; the reader takes either.
            let compressed = format!("{path}.gz");
            if !std::path::Path::new(&path).exists() && std::path::Path::new(&compressed).exists() {
                Ok(compressed)
            } else {
                Ok(path)
            }
        }
        None => Err(format!("--{name} or --reference-dir is required")),
    }
}

/// Exact names and exact synonyms of current (not obsolete) HPO terms, by lower-case text.
struct TermNames {
    by_name: BTreeMap<String, BTreeSet<String>>,
    by_synonym: BTreeMap<String, BTreeSet<String>>,
}

impl TermNames {
    fn load(nodes_tsv: &str) -> Result<TermNames> {
        let (table, _) = Table::read(
            nodes_tsv,
            &["hp_id", "name", "exact_synonyms", "obsolete"],
            &["exact_synonyms", "obsolete"],
        )?;
        let mut by_name: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
        let mut by_synonym: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
        for row in &table.rows {
            if row[3] == "true" {
                continue;
            }
            by_name
                .entry(strip(&row[1]).to_lowercase())
                .or_default()
                .insert(row[0].clone());
            for synonym in row[2].split(" | ") {
                let synonym = strip(synonym).to_lowercase();
                if !synonym.is_empty() {
                    by_synonym
                        .entry(synonym)
                        .or_default()
                        .insert(row[0].clone());
                }
            }
        }
        Ok(TermNames {
            by_name,
            by_synonym,
        })
    }

    /// An HPO id as given (HP:0001234, also hp:0001234 or HP_0001234), or
    /// the one term whose name, else whose exact synonym, is the text.
    fn resolve(&self, text: &str) -> std::result::Result<(String, bool), String> {
        let trimmed = strip(text);
        let bytes = trimmed.as_bytes();
        if bytes.len() == 10
            && trimmed[..2].eq_ignore_ascii_case("HP")
            && (bytes[2] == b':' || bytes[2] == b'_')
            && bytes[3..].iter().all(u8::is_ascii_digit)
        {
            return Ok((format!("HP:{}", &trimmed[3..]), true));
        }
        let key = trimmed.to_lowercase();
        for (table, kind) in [(&self.by_name, "name"), (&self.by_synonym, "exact synonym")] {
            if let Some(ids) = table.get(&key) {
                if ids.len() == 1 {
                    return Ok((ids.iter().next().expect("one id").clone(), false));
                }
                let listed: Vec<&str> = ids.iter().map(String::as_str).collect();
                return Err(format!(
                    "\"{trimmed}\" is the {kind} of several terms: {}",
                    listed.join(", ")
                ));
            }
        }
        Err(format!(
            "\"{trimmed}\" is neither an HPO id nor the exact name or exact synonym of a current HPO term"
        ))
    }
}

fn term_rows(
    ontology: &Ontology,
    ids: &[String],
    texts: &[String],
    names: Option<&TermNames>,
) -> Result<Vec<DirectTerm>> {
    let mut rows: Vec<DirectTerm> = ids
        .iter()
        .map(|id| DirectTerm {
            hp_id: id.clone(),
            source: "manual_hp_id".to_owned(),
            original_entry: if ontology.has_name(id) {
                ontology.name_of(id).to_owned()
            } else {
                id.clone()
            },
        })
        .collect();
    let mut problems = Vec::new();
    for text in texts {
        let names = names.expect("names are loaded for --hp-term");
        match names.resolve(text) {
            Ok((id, true)) => rows.push(DirectTerm {
                original_entry: if ontology.has_name(&id) {
                    ontology.name_of(&id).to_owned()
                } else {
                    id.clone()
                },
                hp_id: id,
                source: "manual_hp_id".to_owned(),
            }),
            Ok((id, false)) => rows.push(DirectTerm {
                hp_id: id,
                source: "manual_hp_term".to_owned(),
                original_entry: strip(text).to_owned(),
            }),
            Err(problem) => problems.push(problem),
        }
    }
    if problems.is_empty() {
        Ok(rows)
    } else {
        Err(problems.join("\n"))
    }
}

/// `load_patient_direct_terms`: rows whose source is a direct source, per sample.
fn load_patient_terms(
    path: &str,
    direct_sources: &HashSet<String>,
) -> Result<IndexMap<String, Vec<DirectTerm>>> {
    let (table, present) = Table::read(
        path,
        &["sample_id", "hp_id", "hp_name", "source", "original_entry"],
        &["original_entry"],
    )?;
    let _ = present;
    let mut by_sample: IndexMap<String, Vec<DirectTerm>> = IndexMap::new();
    for row in table.rows {
        if !direct_sources.contains(&row[3]) {
            continue;
        }
        by_sample
            .entry(row[0].clone())
            .or_default()
            .push(DirectTerm {
                hp_id: row[1].clone(),
                source: row[3].clone(),
                original_entry: row[4].clone(),
            });
    }
    Ok(by_sample)
}

fn load_term_table(
    path: &str,
    ontology: &Ontology,
    names: &TermNames,
) -> Result<IndexMap<String, Vec<DirectTerm>>> {
    let (table, _) = Table::read(path, &["sample_id", "hpo_term"], &[])?;
    let mut by_sample: IndexMap<String, Vec<DirectTerm>> = IndexMap::new();
    let mut problems = Vec::new();
    for row in table.rows {
        match term_rows(ontology, &[], std::slice::from_ref(&row[1]), Some(names)) {
            Ok(rows) => by_sample.entry(row[0].clone()).or_default().extend(rows),
            Err(problem) => problems.push(format!("sample {}: {problem}", row[0])),
        }
    }
    if problems.is_empty() {
        Ok(by_sample)
    } else {
        Err(problems.join("\n"))
    }
}

fn choose_samples(
    by_sample: &IndexMap<String, Vec<DirectTerm>>,
    requested: Vec<String>,
    all: bool,
) -> Vec<String> {
    if !requested.is_empty() {
        return requested;
    }
    let mut ids: Vec<String> = by_sample.keys().cloned().collect();
    ids.sort();
    if all {
        ids
    } else {
        ids.into_iter().take(1).collect()
    }
}

fn distance_reference(
    by_sample: &IndexMap<String, Vec<DirectTerm>>,
    samples: &[String],
    ontology: &Ontology,
) -> f64 {
    let mut largest = 0.0f64;
    for sample in samples {
        let ids: Vec<String> = by_sample
            .get(sample)
            .map(|rows| rows.iter().map(|row| row.hp_id.clone()).collect())
            .unwrap_or_default();
        let p_ids: Vec<String> = ids
            .into_iter()
            .filter(|id| ontology.is_p_term(id))
            .collect();
        let frontier = ontology.frontier_terms(&p_ids);
        for (i, a) in frontier.iter().enumerate() {
            for b in &frontier[i + 1..] {
                let distance = ontology.pair_shortest_distance(a, b) as f64;
                if distance > largest {
                    largest = distance;
                }
            }
        }
    }
    if largest > 1.0 { largest } else { 1.0 }
}

fn write_sample(output_dir: &str, sample_id: &str, tables: &ranking::SampleTables) -> Result<()> {
    let prefix = path_join(output_dir, sample_id);
    let files: [(&str, &[&str], &Vec<Record>); 6] = [
        (
            ".frontier_p_terms.tsv",
            FRONTIER_FIELDNAMES,
            &tables.frontier,
        ),
        (".p_pairs.tsv", PAIR_FIELDNAMES, &tables.pairs),
        (
            ".frontier_gene_evidence.tsv",
            GENE_EVIDENCE_FIELDNAMES,
            &tables.evidence,
        ),
        (
            ".gene_summary.tsv",
            GENE_SUMMARY_FIELDNAMES,
            &tables.summary,
        ),
        (".gene_ranked.tsv", GENE_RANK_FIELDNAMES, &tables.ranked),
        (
            ".gene_score_contributions.tsv",
            GENE_SCORE_CONTRIBUTION_FIELDNAMES,
            &tables.contributions,
        ),
    ];
    for (suffix, fieldnames, records) in files {
        write_records(&with_suffix(&prefix, suffix), fieldnames, records, "\n")?;
    }
    Ok(())
}

pub fn run(arguments: &[String]) -> Result<()> {
    let Some(options) = Options::parse(arguments, VALUED, FLAGS)? else {
        println!("{HELP}");
        return Ok(());
    };
    if options
        .value("graph-tool")
        .is_some_and(|value| !value.is_empty())
    {
        return Err(
            "--graph-tool (a graph-tool .gt file) is not read; give --nodes-tsv and --edges-tsv \
                    from the same build_hpo_owl_graph_tool.py or build-hpo-graph run"
                .to_owned(),
        );
    }
    if options.integer("random-samples", 0)? > 0 {
        return Err("--random-samples is not converted; name the samples with --sample-id or use --all-samples".to_owned());
    }
    let weight_model = options
        .value("weight-model")
        .unwrap_or("inverse_gene_set")
        .to_owned();
    if weight_model != "inverse_gene_set" && weight_model != "tier_sqrt" {
        return Err(format!(
            "--weight-model must be inverse_gene_set or tier_sqrt, not {weight_model}"
        ));
    }
    let output_dir = path_text(options.required("output-dir")?);
    std::fs::create_dir_all(&output_dir).map_err(|error| format!("{output_dir}: {error}"))?;
    let jobs = options.integer("jobs", 1)?.max(1) as usize;

    let nodes_tsv = reference_path(&options, "nodes-tsv")?;
    let ontology = Ontology::from_exports(&nodes_tsv, &reference_path(&options, "edges-tsv")?)?;
    let store = EvidenceStore::load(&EvidencePaths {
        nodes_tsv: nodes_tsv.clone(),
        phenotype_to_genes: reference_path(&options, "phenotype-to-genes")?,
        go_to_genes: reference_path(&options, "go-to-genes")?,
        pr_to_gene_or_protein: reference_path(&options, "pr-to-gene-or-protein")?,
        mp_to_human_orthologs: reference_path(&options, "mp-to-human-orthologs")?,
        orpha_disease_to_genes: reference_path(&options, "orpha-disease-to-genes")?,
        uniprot_evidence: reference_path(&options, "uniprot-evidence")?,
        human_gene_info: reference_path(&options, "human-gene-info")?,
    })?;
    let immune_genes: HashSet<String> = if options.flag("no-external-immune") {
        HashSet::new()
    } else {
        let genes = external_immune_genes(
            &reference_path(&options, "reactome-to-genes")?,
            &reference_path(&options, "reactome-relations")?,
            &reference_path(&options, "gene2go")?,
            &reference_path(&options, "go-obo")?,
            store.ncbi_to_symbol(),
        )?;
        println!("external_immune_genes={}", genes.len());
        genes
    };

    let hp_ids = options.all("hp-id");
    let hp_terms = options.all("hp-term");
    let (by_sample, samples): (IndexMap<String, Vec<DirectTerm>>, Vec<String>) =
        if !hp_ids.is_empty() || !hp_terms.is_empty() {
            let names = if hp_terms.is_empty() {
                None
            } else {
                Some(TermNames::load(&nodes_tsv)?)
            };
            let sample = options
                .value("synthetic-sample-id")
                .unwrap_or("SYNTHETIC_HP_TERMS")
                .to_owned();
            let rows = term_rows(&ontology, &hp_ids, &hp_terms, names.as_ref())?;
            (IndexMap::from([(sample.clone(), rows)]), vec![sample])
        } else if let Some(path) = options.value("hpo-term-table") {
            let names = TermNames::load(&nodes_tsv)?;
            let by_sample = load_term_table(path, &ontology, &names)?;
            let samples = choose_samples(
                &by_sample,
                options.all("sample-id"),
                options.flag("all-samples"),
            );
            (by_sample, samples)
        } else {
            let path = options
                .value("patient-hpo-tsv")
                .ok_or("give --hp-id, --hp-term, --hpo-term-table or --patient-hpo-tsv")?;
            let mut sources: HashSet<String> = options.all("direct-source").into_iter().collect();
            if sources.is_empty() {
                sources = ["rag_direct", "embedded_id", "text_match"]
                    .iter()
                    .map(|s| s.to_string())
                    .collect();
            }
            let by_sample = load_patient_terms(path, &sources)?;
            let samples = choose_samples(
                &by_sample,
                options.all("sample-id"),
                options.flag("all-samples"),
            );
            (by_sample, samples)
        };

    let given_reference = options.float("distance-reference", 0.0)?;
    let reference = if given_reference > 0.0 {
        given_reference
    } else {
        distance_reference(&by_sample, &samples, &ontology)
    };
    println!("distance_reference={}", general(reference));
    let settings = Settings {
        include_secondary_always: options.flag("include-secondary-always"),
        primary_max_gene_set: options.integer("primary-max-gene-set", 100)?,
        weak_max_gene_set: options.integer("weak-max-gene-set", 300)?,
        weight_model,
        distance_reference: reference,
        parent_depth_floor: options.integer("parent-propagation-depth-floor", 3)?,
    };

    let empty: Vec<DirectTerm> = Vec::new();
    let run_one = |sample: &String| -> Result<Record> {
        let rows = by_sample.get(sample).unwrap_or(&empty);
        let tables =
            ranking::process_sample(sample, rows, &ontology, &store, &settings, &immune_genes);
        write_sample(&output_dir, sample, &tables)?;
        Ok(tables.run_summary)
    };
    let summaries: Vec<Record> = if jobs > 1 && samples.len() > 1 {
        let next = AtomicUsize::new(0);
        let results: std::sync::Mutex<Vec<Option<Result<Record>>>> =
            std::sync::Mutex::new((0..samples.len()).map(|_| None).collect());
        std::thread::scope(|scope| {
            for _ in 0..jobs.min(samples.len()) {
                scope.spawn(|| {
                    loop {
                        let index = next.fetch_add(1, Ordering::SeqCst);
                        if index >= samples.len() {
                            break;
                        }
                        let result = run_one(&samples[index]);
                        results.lock().expect("lock")[index] = Some(result);
                    }
                });
            }
        });
        results
            .into_inner()
            .expect("lock")
            .into_iter()
            .map(|result| result.expect("every sample ran"))
            .collect::<Result<Vec<_>>>()?
    } else {
        samples.iter().map(run_one).collect::<Result<Vec<_>>>()?
    };
    write_records(
        &path_join(&output_dir, "run_summary.tsv"),
        RUN_SUMMARY_FIELDNAMES,
        &summaries,
        "\n",
    )?;
    println!("output_dir={output_dir}");
    println!("samples={}", samples.join(","));
    for summary in &summaries {
        let value = |name: &str| {
            summary
                .iter()
                .find(|(key, _)| *key == name)
                .map(|(_, value)| value.render())
                .unwrap_or_default()
        };
        println!(
            "{}\tfrontier={}\tpairs={}\tgene_rows={}\tgenes={}",
            value("sample_id"),
            value("frontier_p_terms"),
            value("pair_rows"),
            value("gene_evidence_rows"),
            value("gene_summary_rows")
        );
    }
    Ok(())
}
