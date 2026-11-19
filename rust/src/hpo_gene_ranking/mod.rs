//! `priva hpo-gene-ranking`: ranks genes by how well their known phenotypes
//! fit a patient's HPO terms. Taken from the `hpo-gene-ranking` program of the
//! llm_gene_reranker repository (rust/, commit bb2df76) with its commands and
//! options unchanged; only this file's entry point differs (priva dispatches
//! `priva hpo-gene-ranking COMMAND ...` here). It is the Rust version of three
//! llm_gene_reranker scripts:
//!
//! | Command | Python script it replaces |
//! |---|---|
//! | `build-hpo-graph` | `scripts/build_hpo_owl_graph_tool.py --hpo-only` (the nodes and edges tables) |
//! | `prepare-gene-evidence` | `scripts/prepare_gene_evidence_resources.py` |
//! | `rank-genes` | `scripts/patient_frontier_gene_evidence.py` |

pub(crate) mod evidence;
pub(crate) mod gene_evidence;
pub(crate) mod hpo_graph;
pub(crate) mod immune;
pub(crate) mod ontology;
pub(crate) mod python_text;
pub(crate) mod rank_genes;
pub(crate) mod ranking;
pub(crate) mod tsv;
pub(crate) mod xml_tree;

use std::collections::HashMap;

pub type Result<T> = std::result::Result<T, String>;

/// The llm_gene_reranker version and commit this copy was taken from.
const UPSTREAM_VERSION: &str = "0.1.0, llm_gene_reranker bb2df76";

const USAGE: &str = "\
priva hpo-gene-ranking: ranks genes by how well their known phenotypes fit a patient's HPO terms

Commands, in the order a full run uses them:
  prepare-gene-evidence  raw public files -> term-to-gene tables
                         (replaces scripts/prepare_gene_evidence_resources.py)
  build-hpo-graph        hp.owl -> HPO nodes and edges tables
                         (replaces the tables of scripts/build_hpo_owl_graph_tool.py --hpo-only)
  rank-genes             HPO terms or ids of each sample -> ranked genes
                         (replaces scripts/patient_frontier_gene_evidence.py)

Run `priva hpo-gene-ranking COMMAND --help` for the options of a command.";

/// `priva hpo-gene-ranking ARGUMENTS`: the program's own commands and options.
pub(crate) fn dispatch(arguments: &[String]) -> std::result::Result<(), Box<dyn std::error::Error>> {
    let result = match arguments.first().map(String::as_str) {
        Some("rank-genes") => rank_genes::run(&arguments[1..]),
        Some("prepare-gene-evidence") => gene_evidence::run(&arguments[1..]),
        Some("build-hpo-graph") => hpo_graph::run(&arguments[1..]),
        Some("--version") => {
            println!("hpo-gene-ranking {} (in priva {})", UPSTREAM_VERSION, env!("CARGO_PKG_VERSION"));
            Ok(())
        }
        None | Some("-h") | Some("--help") | Some("help") => {
            println!("{USAGE}");
            Ok(())
        }
        Some(other) => Err(format!("unknown command {other}\n\n{USAGE}")),
    };
    if let Err(error) = result {
        eprintln!("hpo-gene-ranking: {error}");
        std::process::exit(1);
    }
    Ok(())
}

/// Parsed `--name value`, `--name=value` and `--flag` options.
pub struct Options {
    values: HashMap<String, Vec<String>>,
    flags: Vec<String>,
}

impl Options {
    /// `valued` options take a value (and may repeat); `flags` take none.
    /// Returns None when `--help` was asked for.
    pub fn parse(arguments: &[String], valued: &[&str], flags: &[&str]) -> Result<Option<Options>> {
        let mut values: HashMap<String, Vec<String>> = HashMap::new();
        let mut set_flags = Vec::new();
        let mut index = 0;
        while index < arguments.len() {
            let argument = &arguments[index];
            if argument == "-h" || argument == "--help" {
                return Ok(None);
            }
            let Some(body) = argument.strip_prefix("--") else {
                return Err(format!("unexpected argument {argument}"));
            };
            let (name, inline) = match body.split_once('=') {
                Some((name, value)) => (name, Some(value.to_owned())),
                None => (body, None),
            };
            if flags.contains(&name) {
                if inline.is_some() {
                    return Err(format!("--{name} takes no value"));
                }
                set_flags.push(name.to_owned());
            } else if valued.contains(&name) {
                let value = match inline {
                    Some(value) => value,
                    None => {
                        index += 1;
                        arguments
                            .get(index)
                            .cloned()
                            .ok_or_else(|| format!("--{name} needs a value"))?
                    }
                };
                values.entry(name.to_owned()).or_default().push(value);
            } else {
                return Err(format!("unknown option --{name}"));
            }
            index += 1;
        }
        Ok(Some(Options {
            values,
            flags: set_flags,
        }))
    }

    /// The last value given (argparse keeps the last of repeated options).
    pub fn value(&self, name: &str) -> Option<&str> {
        self.values
            .get(name)
            .and_then(|list| list.last())
            .map(String::as_str)
    }

    pub fn all(&self, name: &str) -> Vec<String> {
        self.values.get(name).cloned().unwrap_or_default()
    }

    pub fn flag(&self, name: &str) -> bool {
        self.flags.iter().any(|flag| flag == name)
    }

    pub fn required(&self, name: &str) -> Result<&str> {
        self.value(name)
            .ok_or_else(|| format!("--{name} is required"))
    }

    pub fn integer(&self, name: &str, default: i64) -> Result<i64> {
        match self.value(name) {
            None => Ok(default),
            Some(text) => text
                .trim()
                .parse()
                .map_err(|_| format!("--{name}: invalid integer {text}")),
        }
    }

    pub fn float(&self, name: &str, default: f64) -> Result<f64> {
        match self.value(name) {
            None => Ok(default),
            Some(text) => text
                .trim()
                .parse()
                .map_err(|_| format!("--{name}: invalid number {text}")),
        }
    }
}
