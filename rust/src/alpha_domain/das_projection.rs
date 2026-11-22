//! DAS inputs projected from the final typed statistics model, in the order of
//! `am_pick_intolerant_domains.py`. Domain score arrays follow
//! `visualize_domain_distribution`: depth-first member order, at least six
//! `max_distribution` values, and a repeated joined path keeps its first
//! position with its last value. Common-missense summaries follow
//! `filter_common_missense`: the last residue dictionary for a joined path wins,
//! and the flat `gnomAD_distribution` fallback applies only when no residue
//! dictionary exists anywhere in the model.
use super::model::{Distribution, Member, Node, Samples};
use super::{EXTRACTOR, PRODUCER, Result, exception, parquet, parse_texts, store};
use crate::das_common::{self, Summary};
use std::collections::HashMap;
use std::path::Path;

/// `visualize_domain_distribution` passes this literal minimum.
const MIN_VARIANTS: usize = 6;
/// `analyze_domain_data` defaults; only the common AF has a source option.
pub(crate) const AF_COMMON: f64 = 0.05;
pub(crate) const AF_RARE: f64 = 0.01;
pub(crate) const RESIDUE_FRACTION: f64 = 0.05;

fn invalid(message: &str) -> Box<dyn std::error::Error> {
    exception("ValueError", message)
}

pub(crate) struct Projection {
    pub(crate) scores: Vec<(String, Vec<f64>)>,
    score_index: HashMap<String, usize>,
    residues: HashMap<String, Summary>,
    flat: HashMap<String, Vec<f64>>,
    any_residue_data: bool,
    af_common: f64,
}

impl Projection {
    /// The summary `filter_common_missense` records for one requested domain.
    pub(crate) fn common(&self, domain: &str) -> Summary {
        if let Some(summary) = self.residues.get(domain) {
            return summary.clone();
        }
        das_common::filter(
            None::<Vec<Vec<f64>>>,
            self.flat.get(domain).map(|values| values.iter().copied()),
            self.any_residue_data,
            self.af_common,
            AF_RARE,
            RESIDUE_FRACTION,
        )
    }

    fn visit(&mut self, node: &Node, path: &mut Vec<String>) {
        for (key, member) in &node.members {
            match member {
                Member::Array(Distribution::Maximum, values) if key == "max_distribution" => {
                    if values.len() >= MIN_VARIANTS {
                        let joined = path.join(":");
                        match self.score_index.get(&joined) {
                            Some(&at) => self.scores[at].1 = values.clone(),
                            None => {
                                self.score_index.insert(joined.clone(), self.scores.len());
                                self.scores.push((joined, values.clone()));
                            }
                        }
                    }
                }
                Member::Array(Distribution::Af, values) if key == "gnomAD_distribution" => {
                    // Flat arrays matter only while no residue dictionary exists.
                    if !self.any_residue_data {
                        self.flat.insert(path.join(":"), values.clone());
                    }
                }
                Member::Residues(residues) => {
                    if !self.any_residue_data {
                        self.any_residue_data = true;
                        self.flat = HashMap::new();
                    }
                    // Repeats of one value cannot change a minimum or maximum,
                    // so each run contributes its value once.
                    let lists = residues
                        .items
                        .iter()
                        .map(|(_, samples)| run_values(samples));
                    let summary = das_common::filter(
                        Some(lists),
                        None::<Vec<f64>>,
                        true,
                        self.af_common,
                        AF_RARE,
                        RESIDUE_FRACTION,
                    );
                    self.residues.insert(path.join(":"), summary);
                }
                Member::Node(child) if !key.contains("distribution") => {
                    path.push(key.clone());
                    self.visit(child, path);
                    path.pop();
                }
                _ => {}
            }
        }
    }
}

fn run_values(samples: &Samples) -> Vec<f64> {
    let mut values = Vec::new();
    samples.visit(|value, _| values.push(value));
    values
}

pub(crate) fn project(
    statistics: &Path,
    source_builder: &str,
    af_common: f64,
) -> Result<Projection> {
    let metadata = parquet::metadata(statistics)?;
    for (key, expected) in [
        ("priva.kind", "statistics"),
        ("priva.stage", "final"),
        ("priva.producer", PRODUCER),
        ("priva.extractor", EXTRACTOR),
        ("priva.builder_source_sha256", source_builder),
    ] {
        if metadata.get(key).map(String::as_str) != Some(expected) {
            return Err(invalid("typed statistics source profile differs"));
        }
    }
    let genes = parse_texts(
        metadata
            .get("priva.keys")
            .ok_or_else(|| invalid("ordered gene inventory absent"))?,
    )?;
    let mut projection = Projection {
        scores: Vec::new(),
        score_index: HashMap::new(),
        residues: HashMap::new(),
        flat: HashMap::new(),
        any_residue_data: false,
        af_common,
    };
    let mut visited = 0;
    store::scan_genes(statistics, |gene, ordinal, node| {
        if ordinal != visited as u64 || genes.get(visited) != Some(&gene) {
            return Err(invalid("typed gene inventory or order differs"));
        }
        let mut path = vec![gene];
        for (database, member) in &node.members {
            // The source iterates each database value with `.items()`.
            let Member::Node(source) = member else {
                return Err(invalid("gene source is not a domain dictionary"));
            };
            path.push(database.clone());
            projection.visit(source, &mut path);
            path.pop();
        }
        visited += 1;
        Ok(())
    })?;
    if visited != genes.len() {
        return Err(invalid("typed model gene coverage differs"));
    }
    Ok(projection)
}
