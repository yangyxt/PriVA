//! Scoring of one sample (`process_sample`, `build_gene_summary` and
//! `build_gene_ranking` in `patient_frontier_gene_evidence.py`).

use crate::hpo_gene_ranking::evidence::{EvidenceRow, EvidenceStore};
use crate::hpo_gene_ranking::ontology::{Ontology, PairSummary};
use crate::hpo_gene_ranking::python_text::{fixed, general, round, sum_floats};
use crate::hpo_gene_ranking::tsv::{Cell, Record};
use indexmap::IndexMap;
use std::collections::{BTreeSet, HashMap, HashSet};

pub const FRONTIER_FIELDNAMES: &[&str] = &[
    "sample_id",
    "hp_id",
    "hp_name",
    "shortest_p_depth",
    "longest_p_depth",
    "original_entries",
    "sources",
    "gene_evidence_rows",
    "direct_hpo_gene_count",
    "goa_gene_count",
    "pr_gene_count",
    "mp_gene_count",
    "orpha_gene_count",
    "term_gene_set_size_for_scoring",
    "term_gene_breadth_tier",
    "term_candidate_role",
    "term_weight",
    "all_gene_symbols",
];

pub const PAIR_FIELDNAMES: &[&str] = &[
    "sample_id",
    "term1_id",
    "term1_name",
    "term1_shortest_p_depth",
    "term1_longest_p_depth",
    "term2_id",
    "term2_name",
    "term2_shortest_p_depth",
    "term2_longest_p_depth",
    "shortest_distance",
    "nearest_common_ancestor_id",
    "nearest_common_ancestor_name",
    "nearest_common_ancestor_gene_count",
    "nearest_common_ancestor_descendant_leaf_count",
    "nearest_common_ancestor_specificity",
    "nearest_common_ancestor_specificity_source",
    "distance_specificity_factor",
    "all_tied_common_ancestors",
    "nearest_informative_common_ancestors",
    "mica_id",
    "mica_name",
    "mica_tied_common_ancestors",
    "same_top_P_branch",
    "term1_top_P_branches",
    "term2_top_P_branches",
    "shared_top_P_branches",
    "cross_module",
    "cross_module_pair_key",
];

pub const GENE_EVIDENCE_FIELDNAMES: &[&str] = &[
    "sample_id",
    "hp_id",
    "hp_name",
    "gene_symbol",
    "ncbi_gene_id",
    "uniprot_accession",
    "candidate_source",
    "candidate_status",
    "external_id",
    "external_name",
    "gene_set_size",
    "source_quality",
    "term_gene_breadth_tier",
    "term_candidate_role",
    "term_weight",
    "evidence_detail",
    "evidence_row_count",
    "propagated_from_hp_id",
    "propagated_from_hp_name",
];

pub const GENE_SUMMARY_FIELDNAMES: &[&str] = &[
    "sample_id",
    "gene_symbol",
    "ncbi_gene_ids",
    "frontier_hp_count",
    "frontier_hp_ids",
    "frontier_hp_names",
    "candidate_sources",
    "candidate_statuses",
    "external_ids",
    "evidence_row_count",
];

pub const GENE_RANK_FIELDNAMES: &[&str] = &[
    "sample_id",
    "gene_symbol",
    "ncbi_gene_ids",
    "total_score",
    "base_term_score",
    "distance_bonus_raw_score",
    "distance_bonus_score",
    "distance_bonus_cap",
    "cross_module_bonus_raw_score",
    "cross_module_bonus_score",
    "cross_module_bonus_cap",
    "cross_module_pair_count",
    "coherence_coefficient",
    "coherence_bonus_score",
    "coherence_pair_count",
    "coherence_cross_module_pair_count",
    "distance_reference",
    "immune_context_penalty",
    "raw_total_score_before_immune_penalty",
    "candidate_gene_role",
    "immune_supported",
    "immune_support_source",
    "immune_candidate_generating_hp_count",
    "nonimmune_candidate_generating_hp_count",
    "candidate_generating_hp_count",
    "support_only_hp_count",
    "frontier_hp_count",
    "candidate_generating_hp_ids",
    "candidate_generating_hp_names",
    "support_only_hp_ids",
    "support_only_hp_names",
    "all_frontier_hp_ids",
    "all_frontier_hp_names",
    "max_pair_distance_bonus",
    "bonus_pair_details",
    "max_cross_module_bonus",
    "cross_module_bonus_details",
    "candidate_sources",
    "candidate_statuses",
    "evidence_row_count",
];

pub const GENE_SCORE_CONTRIBUTION_FIELDNAMES: &[&str] = &[
    "sample_id",
    "gene_symbol",
    "contribution_type",
    "selected_for_hp_gene",
    "contribution_score",
    "raw_contribution_score",
    "hp_id",
    "hp_name",
    "hp2_id",
    "hp2_name",
    "candidate_source",
    "candidate_status",
    "external_id",
    "external_name",
    "gene_set_size",
    "source_quality",
    "term_weight",
    "shortest_distance",
    "nearest_common_ancestor_id",
    "nearest_common_ancestor_name",
    "nearest_common_ancestor_specificity",
    "nearest_common_ancestor_specificity_source",
    "distance_specificity_factor",
    "cross_module",
    "module_pair_key",
    "cross_module_bonus_factor",
    "pair_bonus_scale_due_to_cap",
    "distance_reference",
    "distance_factor",
    "ancestor_bonus_factor",
    "pair_topology_factor",
    "pair_source_quality",
    "coherence_bucket",
    "immune_context_penalty",
    "reason",
];

pub const RUN_SUMMARY_FIELDNAMES: &[&str] = &[
    "sample_id",
    "direct_terms",
    "direct_p_terms",
    "frontier_p_terms",
    "pair_rows",
    "gene_evidence_rows",
    "gene_summary_rows",
];

const CANDIDATE_GENERATING_ROLES: [&str; 3] = [
    "primary_candidate_generator",
    "weak_candidate_generator",
    "broad_low_weight_candidate_generator",
];

const STRONG_COHERENCE_STATUS_SOURCE_PAIRS: [(&str, &str); 7] = [
    (
        "direct_responsible_gene_candidate",
        "HPO_PHENOTYPE_TO_GENES",
    ),
    (
        "functional_support_for_direct_candidate",
        "HPO_LINKED_GO_TO_GOA_EXPERIMENTAL",
    ),
    (
        "secondary_candidate_from_goa",
        "HPO_LINKED_GO_TO_GOA_EXPERIMENTAL",
    ),
    (
        "protein_bridge_support_for_direct_candidate",
        "HPO_LINKED_PR_TO_PRO_UNIPROT",
    ),
    (
        "secondary_candidate_from_pr_proximal_molecule",
        "HPO_LINKED_PR_TO_PRO_UNIPROT",
    ),
    (
        "mouse_ortholog_support_for_direct_candidate",
        "HPO_XREF_MP_TO_MGI_HUMAN_ORTHOLOG",
    ),
    (
        "secondary_candidate_from_mp_mouse_ortholog",
        "HPO_XREF_MP_TO_MGI_HUMAN_ORTHOLOG",
    ),
];

/// Every source in SOURCE_QUALITY has quality 1.0, as has any other source.
fn source_quality(_source: &str) -> f64 {
    1.0
}

/// The options of one ranking run.
pub struct Settings {
    pub include_secondary_always: bool,
    pub primary_max_gene_set: i64,
    pub weak_max_gene_set: i64,
    /// "inverse_gene_set" or "tier_sqrt".
    pub weight_model: String,
    pub distance_reference: f64,
    pub parent_depth_floor: i64,
}

/// One input row of a sample: an HPO id with where it came from.
#[derive(Clone)]
pub struct DirectTerm {
    pub hp_id: String,
    pub source: String,
    pub original_entry: String,
}

/// Every table of one sample, ready to be written.
pub struct SampleTables {
    pub frontier: Vec<Record>,
    pub pairs: Vec<Record>,
    pub evidence: Vec<Record>,
    pub summary: Vec<Record>,
    pub ranked: Vec<Record>,
    pub contributions: Vec<Record>,
    pub run_summary: Record,
}

pub fn classify_gene_set_size(
    gene_set_size: i64,
    primary_max: i64,
    weak_max: i64,
    include_broad_candidates: bool,
) -> (&'static str, &'static str) {
    if gene_set_size <= 0 {
        ("no_gene_annotation", "no_candidate")
    } else if gene_set_size <= primary_max {
        ("well_specified", "primary_candidate_generator")
    } else if gene_set_size <= weak_max {
        ("broad_weak", "weak_candidate_generator")
    } else if include_broad_candidates {
        ("very_broad", "broad_low_weight_candidate_generator")
    } else {
        ("too_broad", "support_only")
    }
}

pub fn term_weight(
    gene_set_size: i64,
    role: &str,
    primary_max: i64,
    weight_model: &str,
    quality: f64,
) -> f64 {
    if gene_set_size <= 0 || role == "no_candidate" {
        return 0.0;
    }
    if weight_model == "inverse_gene_set" {
        return round(quality / gene_set_size as f64, 9);
    }
    let base = (primary_max as f64 / gene_set_size as f64).sqrt();
    let multiplier = match role {
        "primary_candidate_generator" => 1.0,
        "weak_candidate_generator" => 0.25,
        _ => 0.05,
    };
    round(quality * multiplier * base, 6)
}

fn evidence_record(sample_id: &str, row: &EvidenceRow) -> Record {
    vec![
        ("sample_id", Cell::from(sample_id)),
        ("hp_id", Cell::from(&row.hp_id)),
        ("hp_name", Cell::from(&row.hp_name)),
        ("gene_symbol", Cell::from(&row.gene_symbol)),
        ("ncbi_gene_id", Cell::from(&row.ncbi_gene_id)),
        ("uniprot_accession", Cell::from(&row.uniprot_accession)),
        ("candidate_source", Cell::from(row.candidate_source)),
        ("candidate_status", Cell::from(row.candidate_status)),
        ("external_id", Cell::from(&row.external_id)),
        ("external_name", Cell::from(&row.external_name)),
        ("gene_set_size", Cell::Int(row.gene_set_size)),
        ("source_quality", Cell::Float(row.source_quality)),
        (
            "term_gene_breadth_tier",
            Cell::from(row.term_gene_breadth_tier),
        ),
        ("term_candidate_role", Cell::from(row.term_candidate_role)),
        ("term_weight", Cell::Float(row.term_weight)),
        ("evidence_detail", Cell::from(&row.evidence_detail)),
        ("evidence_row_count", Cell::Int(row.evidence_row_count)),
        (
            "propagated_from_hp_id",
            Cell::from(&row.propagated_from_hp_id),
        ),
        (
            "propagated_from_hp_name",
            Cell::from(&row.propagated_from_hp_name),
        ),
    ]
}

fn joined_set<'a>(values: impl Iterator<Item = &'a str>, skip_empty: bool) -> String {
    values
        .filter(|value| !skip_empty || !value.is_empty())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>()
        .join("|")
}

fn build_gene_summary(sample_id: &str, rows: &[EvidenceRow]) -> Vec<Record> {
    let mut grouped: std::collections::BTreeMap<&str, Vec<&EvidenceRow>> = Default::default();
    for row in rows {
        grouped
            .entry(row.gene_symbol.as_str())
            .or_default()
            .push(row);
    }
    grouped
        .into_iter()
        .map(|(gene_symbol, rows)| {
            vec![
                ("sample_id", Cell::from(sample_id)),
                ("gene_symbol", Cell::from(gene_symbol)),
                (
                    "ncbi_gene_ids",
                    Cell::from(joined_set(
                        rows.iter().map(|r| r.ncbi_gene_id.as_str()),
                        true,
                    )),
                ),
                (
                    "frontier_hp_count",
                    Cell::from(
                        rows.iter()
                            .map(|r| r.hp_id.as_str())
                            .collect::<HashSet<_>>()
                            .len(),
                    ),
                ),
                (
                    "frontier_hp_ids",
                    Cell::from(joined_set(rows.iter().map(|r| r.hp_id.as_str()), false)),
                ),
                (
                    "frontier_hp_names",
                    Cell::from(joined_set(rows.iter().map(|r| r.hp_name.as_str()), false)),
                ),
                (
                    "candidate_sources",
                    Cell::from(joined_set(rows.iter().map(|r| r.candidate_source), false)),
                ),
                (
                    "candidate_statuses",
                    Cell::from(joined_set(rows.iter().map(|r| r.candidate_status), false)),
                ),
                (
                    "external_ids",
                    Cell::from(joined_set(
                        rows.iter().map(|r| r.external_id.as_str()),
                        true,
                    )),
                ),
                ("evidence_row_count", Cell::from(rows.len())),
            ]
        })
        .collect()
}

fn coherence_pair_bucket(pair: &PairSummary) -> String {
    if pair.cross_module {
        return if pair.cross_module_pair_key.is_empty() {
            String::new()
        } else {
            format!("cross_module:{}", pair.cross_module_pair_key)
        };
    }
    if pair.nearest_common_ancestor_id.is_empty() {
        return String::new();
    }
    format!(
        "within_module:{}:{}",
        pair.shared_top_p_branches, pair.nearest_common_ancestor_id
    )
}

/// A pair of frontier terms with the specificity values `annotate_pair_specificity` adds.
pub struct ScoredPair {
    pub summary: PairSummary,
    pub specificity: f64,
    pub specificity_source: &'static str,
}

/// One retained pair coherence candidate of a gene.
struct CoherenceCandidate {
    record: Record,
    raw_contribution: f64,
    cross_module: bool,
    bucket: String,
    detail: String,
}

fn base_contribution(
    sample_id: &str,
    gene_symbol: &str,
    row: &EvidenceRow,
    selected: bool,
) -> Record {
    let weight = round(row.term_weight, 9);
    vec![
        ("sample_id", Cell::from(sample_id)),
        ("gene_symbol", Cell::from(gene_symbol)),
        ("contribution_type", Cell::from("term_base")),
        (
            "selected_for_hp_gene",
            Cell::from(if selected { "true" } else { "false" }),
        ),
        (
            "contribution_score",
            Cell::Float(if selected { weight } else { 0.0 }),
        ),
        ("raw_contribution_score", Cell::Float(weight)),
        ("hp_id", Cell::from(&row.hp_id)),
        ("hp_name", Cell::from(&row.hp_name)),
        ("hp2_id", Cell::from("")),
        ("hp2_name", Cell::from("")),
        ("candidate_source", Cell::from(row.candidate_source)),
        ("candidate_status", Cell::from(row.candidate_status)),
        ("external_id", Cell::from(&row.external_id)),
        ("external_name", Cell::from(&row.external_name)),
        ("gene_set_size", Cell::Int(row.gene_set_size)),
        ("source_quality", Cell::Float(row.source_quality)),
        ("term_weight", Cell::Float(row.term_weight)),
        (
            "reason",
            Cell::from(if selected {
                "selected max source contribution for this gene within this HP term; weight=1/source_gene_set_size"
            } else {
                "not added to base score because another source gave an equal or stronger contribution for this gene within this HP term"
            }),
        ),
    ]
}

fn is_candidate_generating(role: &str) -> bool {
    CANDIDATE_GENERATING_ROLES.contains(&role)
}

fn is_strong_coherence_link(row: &EvidenceRow) -> bool {
    STRONG_COHERENCE_STATUS_SOURCE_PAIRS.contains(&(row.candidate_status, row.candidate_source))
}

fn build_gene_ranking(
    sample_id: &str,
    rows: &[EvidenceRow],
    pairs: &[ScoredPair],
    immune_hp_ids: &HashSet<String>,
    distance_reference: f64,
    external_immune_genes: &HashSet<String>,
) -> (Vec<Record>, Vec<Record>) {
    let mut pair_by_terms: HashMap<(String, String), &ScoredPair> = HashMap::new();
    for pair in pairs {
        let (a, b) = (&pair.summary.term1_id, &pair.summary.term2_id);
        let key = if a <= b {
            (a.clone(), b.clone())
        } else {
            (b.clone(), a.clone())
        };
        pair_by_terms.insert(key, pair);
    }
    let mut grouped: IndexMap<&str, Vec<usize>> = IndexMap::new();
    for (index, row) in rows.iter().enumerate() {
        grouped
            .entry(row.gene_symbol.as_str())
            .or_default()
            .push(index);
    }
    let distance_reference = distance_reference.max(1.0);
    let mut ranked: Vec<(bool, f64, String, Record)> = Vec::new();
    let mut contributions: Vec<Record> = Vec::new();
    for (gene_symbol, members) in &grouped {
        // by_hp: the first row with the largest weight per HP term, in first-seen order.
        let mut by_hp: IndexMap<&str, usize> = IndexMap::new();
        for &index in members {
            let row = &rows[index];
            match by_hp.get(row.hp_id.as_str()) {
                Some(&current) if row.term_weight <= rows[current].term_weight => {}
                _ => {
                    by_hp.insert(row.hp_id.as_str(), index);
                }
            }
        }
        let selected: HashSet<usize> = by_hp.values().copied().collect();
        for &index in members {
            contributions.push(base_contribution(
                sample_id,
                gene_symbol,
                &rows[index],
                selected.contains(&index),
            ));
        }
        let base_score = sum_floats(
            &by_hp
                .values()
                .map(|&index| rows[index].term_weight)
                .collect::<Vec<_>>(),
        );
        let candidate_hps: Vec<(&str, usize)> = by_hp
            .iter()
            .filter(|(_, index)| is_candidate_generating(rows[**index].term_candidate_role))
            .map(|(hp, index)| (*hp, *index))
            .collect();
        let immune_candidates = candidate_hps
            .iter()
            .filter(|(hp, _)| immune_hp_ids.contains(*hp))
            .count();
        let nonimmune_candidates = candidate_hps.len() - immune_candidates;
        let support_hps: Vec<(&str, usize)> = by_hp
            .iter()
            .filter(|(_, index)| rows[**index].term_candidate_role == "support_only")
            .map(|(hp, index)| (*hp, *index))
            .collect();
        // The strongest strong-coherence row per HP term: largest (quality, weight), first kept on ties.
        let mut strong: HashMap<&str, usize> = HashMap::new();
        for &index in members {
            let row = &rows[index];
            if !is_candidate_generating(row.term_candidate_role) || !is_strong_coherence_link(row) {
                continue;
            }
            let better = match strong.get(row.hp_id.as_str()) {
                None => true,
                Some(&current) => {
                    let current_row = &rows[current];
                    (row.source_quality, row.term_weight)
                        > (current_row.source_quality, current_row.term_weight)
                }
            };
            if better {
                strong.insert(row.hp_id.as_str(), index);
            }
        }
        let mut strong_items: Vec<(&str, usize)> = strong.into_iter().collect();
        strong_items.sort_by(|a, b| a.0.cmp(b.0));
        let mut candidates: IndexMap<String, CoherenceCandidate> = IndexMap::new();
        for (i, &(hp_a, index_a)) in strong_items.iter().enumerate() {
            for &(hp_b, index_b) in &strong_items[i + 1..] {
                let key = if hp_a <= hp_b {
                    (hp_a.to_owned(), hp_b.to_owned())
                } else {
                    (hp_b.to_owned(), hp_a.to_owned())
                };
                let Some(pair) = pair_by_terms.get(&key) else {
                    continue;
                };
                let bucket = coherence_pair_bucket(&pair.summary);
                if bucket.is_empty() {
                    continue;
                }
                let (row_a, row_b) = (&rows[index_a], &rows[index_b]);
                let shortest_distance = pair.summary.shortest_distance as f64;
                let specificity = pair.specificity;
                let distance_factor = shortest_distance / distance_reference;
                let ancestor_bonus_factor = 1.0 + specificity;
                let pair_topology_factor = (1.0 + distance_factor) * ancestor_bonus_factor;
                let pair_source_quality = (row_a.source_quality * row_b.source_quality).sqrt();
                let pair_coherence = pair_source_quality * pair_topology_factor;
                if pair_coherence <= 0.0 {
                    continue;
                }
                let detail = format!(
                    "{hp_a} + {hp_b}: bucket={bucket}, distance={}, D_ref={}, distance_factor={}, \
                     ancestor_specificity={}, ancestor_bonus_factor={}, pair_topology_factor={}, \
                     pair_source_quality={}, coherence={}",
                    general(shortest_distance),
                    general(distance_reference),
                    fixed(distance_factor, 6),
                    fixed(specificity, 6),
                    fixed(ancestor_bonus_factor, 6),
                    fixed(pair_topology_factor, 6),
                    fixed(pair_source_quality, 6),
                    fixed(pair_coherence, 6),
                );
                let raw = round(pair_coherence, 9);
                let ancestor_id = pair.summary.nearest_common_ancestor_id.clone();
                let ancestor_name = pair
                    .summary
                    .record
                    .iter()
                    .find(|(key, _)| *key == "nearest_common_ancestor_name")
                    .map(|(_, value)| value.render())
                    .unwrap_or_default();
                let record: Record = vec![
                    ("sample_id", Cell::from(sample_id)),
                    ("gene_symbol", Cell::from(*gene_symbol)),
                    ("contribution_type", Cell::from("pair_coherence_bonus")),
                    ("selected_for_hp_gene", Cell::from("true")),
                    ("contribution_score", Cell::Float(0.0)),
                    ("raw_contribution_score", Cell::Float(raw)),
                    ("hp_id", Cell::from(hp_a)),
                    ("hp_name", Cell::from(&row_a.hp_name)),
                    ("hp2_id", Cell::from(hp_b)),
                    ("hp2_name", Cell::from(&row_b.hp_name)),
                    ("candidate_source", Cell::from("")),
                    ("candidate_status", Cell::from("")),
                    ("external_id", Cell::from("")),
                    ("external_name", Cell::from("")),
                    ("gene_set_size", Cell::from("")),
                    ("source_quality", Cell::Float(round(pair_source_quality, 9))),
                    ("term_weight", Cell::from("")),
                    ("shortest_distance", Cell::Float(shortest_distance)),
                    ("nearest_common_ancestor_id", Cell::from(ancestor_id)),
                    ("nearest_common_ancestor_name", Cell::from(ancestor_name)),
                    (
                        "nearest_common_ancestor_specificity",
                        Cell::Float(specificity),
                    ),
                    (
                        "nearest_common_ancestor_specificity_source",
                        Cell::from(pair.specificity_source),
                    ),
                    (
                        "distance_specificity_factor",
                        Cell::Float(round(pair_topology_factor, 9)),
                    ),
                    (
                        "cross_module",
                        Cell::from(if pair.summary.cross_module {
                            "true"
                        } else {
                            "false"
                        }),
                    ),
                    (
                        "module_pair_key",
                        Cell::from(&pair.summary.cross_module_pair_key),
                    ),
                    ("cross_module_bonus_factor", Cell::from("")),
                    ("pair_bonus_scale_due_to_cap", Cell::from("")),
                    (
                        "distance_reference",
                        Cell::Float(round(distance_reference, 9)),
                    ),
                    ("distance_factor", Cell::Float(round(distance_factor, 9))),
                    (
                        "ancestor_bonus_factor",
                        Cell::Float(round(ancestor_bonus_factor, 9)),
                    ),
                    (
                        "pair_topology_factor",
                        Cell::Float(round(pair_topology_factor, 9)),
                    ),
                    (
                        "pair_source_quality",
                        Cell::Float(round(pair_source_quality, 9)),
                    ),
                    ("coherence_bucket", Cell::from(bucket.clone())),
                    ("immune_context_penalty", Cell::from("")),
                    (
                        "reason",
                        Cell::from(
                            "same gene has strong evidence for both frontier P terms; unified pair coherence uses \
                             source_quality * (1 + distance/D_ref) * (1 + ancestor specificity), with one best pair \
                             retained per coherence bucket",
                        ),
                    ),
                ];
                // The new unrounded coherence is compared with the kept, rounded one.
                let replace = match candidates.get(&bucket) {
                    None => true,
                    Some(current) => pair_coherence > current.raw_contribution,
                };
                if replace {
                    candidates.insert(
                        bucket.clone(),
                        CoherenceCandidate {
                            record,
                            raw_contribution: raw,
                            cross_module: pair.summary.cross_module,
                            bucket,
                            detail,
                        },
                    );
                }
            }
        }
        let raw_coherence_coefficient = sum_floats(
            &candidates
                .values()
                .map(|candidate| candidate.raw_contribution)
                .collect::<Vec<_>>(),
        );
        let coherence_bonus = base_score * raw_coherence_coefficient;
        let mut cross_module_pairs = 0i64;
        let mut max_pair_bonus = 0.0f64;
        let mut max_cross_module_bonus = 0.0f64;
        let mut bonus_details = Vec::new();
        let mut cross_module_details = Vec::new();
        let mut ordered: Vec<CoherenceCandidate> = candidates.into_values().collect();
        ordered.sort_by(|a, b| a.bucket.cmp(&b.bucket));
        let coherence_pair_count = ordered.len();
        for mut candidate in ordered {
            let pair_coherence = candidate.raw_contribution;
            for (key, value) in candidate.record.iter_mut() {
                if *key == "contribution_score" {
                    *value = Cell::Float(round(base_score * pair_coherence, 9));
                }
            }
            if pair_coherence > max_pair_bonus {
                max_pair_bonus = pair_coherence;
            }
            bonus_details.push(candidate.detail.clone());
            if candidate.cross_module {
                cross_module_pairs += 1;
                if pair_coherence > max_cross_module_bonus {
                    max_cross_module_bonus = pair_coherence;
                }
                cross_module_details.push(candidate.detail.clone());
            }
            contributions.push(candidate.record);
        }
        let raw_total_score = base_score * (1.0 + raw_coherence_coefficient);
        let hpo_immune_supported = immune_candidates > 0;
        let external_immune_supported = external_immune_genes.contains(*gene_symbol);
        let immune_supported = hpo_immune_supported || external_immune_supported;
        let immune_context_penalty =
            if !immune_hp_ids.is_empty() && !candidate_hps.is_empty() && !immune_supported {
                0.5
            } else {
                1.0
            };
        let total_score = raw_total_score * immune_context_penalty;
        if immune_context_penalty < 1.0 {
            contributions.push(vec![
                ("sample_id", Cell::from(sample_id)),
                ("gene_symbol", Cell::from(*gene_symbol)),
                ("contribution_type", Cell::from("immune_context_penalty")),
                ("selected_for_hp_gene", Cell::from("true")),
                ("contribution_score", Cell::Float(round(total_score - raw_total_score, 9))),
                ("raw_contribution_score", Cell::Float(round(raw_total_score, 9))),
                ("immune_context_penalty", Cell::Float(immune_context_penalty)),
                (
                    "reason",
                    Cell::from(
                        "patient has immune frontier P terms, but this gene is supported only by non-immune frontier P terms",
                    ),
                ),
            ]);
        }
        let roles: Vec<&str> = by_hp
            .values()
            .map(|&index| rows[index].term_candidate_role)
            .collect();
        let gene_role = if roles.contains(&"primary_candidate_generator") {
            "primary_candidate"
        } else if roles.contains(&"weak_candidate_generator") {
            "weak_candidate"
        } else if roles.contains(&"broad_low_weight_candidate_generator") {
            "broad_low_weight_candidate"
        } else {
            "support_only"
        };
        let mut all_hp: Vec<(&str, usize)> =
            by_hp.iter().map(|(hp, index)| (*hp, *index)).collect();
        all_hp.sort_by(|a, b| a.0.cmp(b.0));
        let ids = |items: &[(&str, usize)]| {
            items
                .iter()
                .map(|(hp, _)| *hp)
                .collect::<Vec<_>>()
                .join("|")
        };
        let names = |items: &[(&str, usize)]| {
            items
                .iter()
                .map(|(_, index)| rows[*index].hp_name.as_str())
                .collect::<Vec<_>>()
                .join("|")
        };
        let member_rows = || members.iter().map(|&index| &rows[index]);
        let immune_support_source = {
            let mut text = String::new();
            if hpo_immune_supported {
                text.push_str("hpo");
            }
            if hpo_immune_supported && external_immune_supported {
                text.push('+');
            }
            if external_immune_supported {
                text.push_str("external");
            }
            if text.is_empty() {
                "none".to_owned()
            } else {
                text
            }
        };
        let total_rounded = round(total_score, 6);
        // Python's sum() of no pairs is the integer 0, which prints without ".0".
        let coefficient_cell = if coherence_pair_count == 0 {
            Cell::Int(0)
        } else {
            Cell::Float(round(raw_coherence_coefficient, 6))
        };
        let record: Record = vec![
            ("sample_id", Cell::from(sample_id)),
            ("gene_symbol", Cell::from(*gene_symbol)),
            (
                "ncbi_gene_ids",
                Cell::from(joined_set(
                    member_rows().map(|r| r.ncbi_gene_id.as_str()),
                    true,
                )),
            ),
            ("total_score", Cell::Float(total_rounded)),
            ("base_term_score", Cell::Float(round(base_score, 6))),
            ("distance_bonus_raw_score", coefficient_cell.clone()),
            (
                "distance_bonus_score",
                Cell::Float(round(coherence_bonus, 6)),
            ),
            ("distance_bonus_cap", Cell::Float(round(0.0, 6))),
            ("cross_module_bonus_raw_score", Cell::Float(round(0.0, 6))),
            ("cross_module_bonus_score", Cell::Float(round(0.0, 6))),
            ("cross_module_bonus_cap", Cell::Float(round(0.0, 6))),
            ("cross_module_pair_count", Cell::Int(cross_module_pairs)),
            ("coherence_coefficient", coefficient_cell),
            (
                "coherence_bonus_score",
                Cell::Float(round(coherence_bonus, 6)),
            ),
            ("coherence_pair_count", Cell::from(coherence_pair_count)),
            (
                "coherence_cross_module_pair_count",
                Cell::Int(cross_module_pairs),
            ),
            (
                "distance_reference",
                Cell::Float(round(distance_reference, 6)),
            ),
            (
                "immune_context_penalty",
                Cell::Float(round(immune_context_penalty, 6)),
            ),
            (
                "raw_total_score_before_immune_penalty",
                Cell::Float(round(raw_total_score, 6)),
            ),
            ("candidate_gene_role", Cell::from(gene_role)),
            (
                "immune_supported",
                Cell::from(if immune_supported { "true" } else { "false" }),
            ),
            ("immune_support_source", Cell::from(immune_support_source)),
            (
                "immune_candidate_generating_hp_count",
                Cell::from(immune_candidates),
            ),
            (
                "nonimmune_candidate_generating_hp_count",
                Cell::from(nonimmune_candidates),
            ),
            (
                "candidate_generating_hp_count",
                Cell::from(candidate_hps.len()),
            ),
            ("support_only_hp_count", Cell::from(support_hps.len())),
            ("frontier_hp_count", Cell::from(by_hp.len())),
            (
                "candidate_generating_hp_ids",
                Cell::from(ids(&candidate_hps)),
            ),
            (
                "candidate_generating_hp_names",
                Cell::from(names(&candidate_hps)),
            ),
            ("support_only_hp_ids", Cell::from(ids(&support_hps))),
            ("support_only_hp_names", Cell::from(names(&support_hps))),
            ("all_frontier_hp_ids", Cell::from(ids(&all_hp))),
            ("all_frontier_hp_names", Cell::from(names(&all_hp))),
            (
                "max_pair_distance_bonus",
                Cell::Float(round(max_pair_bonus, 6)),
            ),
            ("bonus_pair_details", Cell::from(bonus_details.join(" | "))),
            (
                "max_cross_module_bonus",
                Cell::Float(round(max_cross_module_bonus, 6)),
            ),
            (
                "cross_module_bonus_details",
                Cell::from(cross_module_details.join(" | ")),
            ),
            (
                "candidate_sources",
                Cell::from(joined_set(member_rows().map(|r| r.candidate_source), false)),
            ),
            (
                "candidate_statuses",
                Cell::from(joined_set(member_rows().map(|r| r.candidate_status), false)),
            ),
            ("evidence_row_count", Cell::from(members.len())),
        ];
        ranked.push((
            gene_role == "support_only",
            total_rounded,
            gene_symbol.to_string(),
            record,
        ));
    }
    ranked.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then_with(|| {
                (-a.1)
                    .partial_cmp(&-b.1)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| a.2.cmp(&b.2))
    });
    (
        ranked.into_iter().map(|(_, _, _, record)| record).collect(),
        contributions,
    )
}

/// `process_sample` without the file writing.
pub fn process_sample(
    sample_id: &str,
    direct_rows: &[DirectTerm],
    ontology: &Ontology,
    store: &EvidenceStore,
    settings: &Settings,
    external_immune_genes: &HashSet<String>,
) -> SampleTables {
    let direct_ids: Vec<String> = direct_rows.iter().map(|row| row.hp_id.clone()).collect();
    let direct_p_ids: Vec<String> = direct_ids
        .iter()
        .filter(|id| ontology.is_p_term(id))
        .cloned()
        .collect();
    let frontier = ontology.frontier_terms(&direct_p_ids);
    let frontier_set: HashSet<&String> = frontier.iter().collect();
    let mut originals: HashMap<&str, Vec<&str>> = HashMap::new();
    let mut sources: HashMap<&str, BTreeSet<&str>> = HashMap::new();
    for row in direct_rows {
        if frontier_set.contains(&row.hp_id) {
            originals
                .entry(row.hp_id.as_str())
                .or_default()
                .push(row.original_entry.as_str());
            sources
                .entry(row.hp_id.as_str())
                .or_default()
                .insert(row.source.as_str());
        }
    }
    let include_broad = settings.weight_model == "inverse_gene_set";
    let mut frontier_records = Vec::new();
    let mut all_rows: Vec<EvidenceRow> = Vec::new();
    for hp_id in &frontier {
        let hp_name = ontology.name_of(hp_id);
        let gene_rows = store.gene_evidence_for_frontier_term(
            hp_id,
            hp_name,
            settings.include_secondary_always,
            ontology,
            settings.parent_depth_floor,
        );
        let count_genes = |test: &dyn Fn(&EvidenceRow) -> bool| {
            gene_rows
                .iter()
                .filter(|row| test(row))
                .map(|row| row.gene_symbol.as_str())
                .collect::<HashSet<_>>()
                .len()
        };
        let direct_hpo_gene_count =
            count_genes(&|row| row.candidate_source == crate::hpo_gene_ranking::evidence::HPO_SOURCE);
        let mut positive_sizes = Vec::new();
        let mut best: Option<(&'static str, &'static str, f64)> = None;
        for row in &gene_rows {
            let size = row.gene_set_size;
            if size > 0 {
                positive_sizes.push(size);
            }
            let (tier, role) = classify_gene_set_size(
                size,
                settings.primary_max_gene_set,
                settings.weak_max_gene_set,
                include_broad,
            );
            let quality = source_quality(row.candidate_source);
            let weight = term_weight(
                size,
                role,
                settings.primary_max_gene_set,
                &settings.weight_model,
                quality,
            );
            if best.is_none_or(|(_, _, best_weight)| weight > best_weight) {
                best = Some((tier, role, weight));
            }
            let mut scored = row.clone();
            scored.source_quality = quality;
            scored.term_gene_breadth_tier = tier;
            scored.term_candidate_role = role;
            scored.term_weight = weight;
            all_rows.push(scored);
        }
        let scoring_size = positive_sizes.iter().copied().min().unwrap_or(0);
        let (tier, role, weight) = best.unwrap_or_else(|| {
            let (tier, role) = classify_gene_set_size(
                0,
                settings.primary_max_gene_set,
                settings.weak_max_gene_set,
                include_broad,
            );
            (tier, role, 0.0)
        });
        let original_text = originals
            .get(hp_id.as_str())
            .map(|entries| joined_set(entries.iter().copied(), true))
            .unwrap_or_default();
        let source_text = sources
            .get(hp_id.as_str())
            .map(|set| set.iter().copied().collect::<Vec<_>>().join("|"))
            .unwrap_or_default();
        frontier_records.push(vec![
            ("sample_id", Cell::from(sample_id)),
            ("hp_id", Cell::from(hp_id)),
            ("hp_name", Cell::from(hp_name)),
            (
                "shortest_p_depth",
                Cell::Int(ontology.shortest_p_depth(hp_id)),
            ),
            (
                "longest_p_depth",
                Cell::Int(ontology.longest_p_depth(hp_id)),
            ),
            ("original_entries", Cell::from(original_text)),
            ("sources", Cell::from(source_text)),
            ("gene_evidence_rows", Cell::from(gene_rows.len())),
            ("direct_hpo_gene_count", Cell::from(direct_hpo_gene_count)),
            (
                "goa_gene_count",
                Cell::from(count_genes(&|row| {
                    row.candidate_source.starts_with("HPO_LINKED_GO_TO_GOA")
                })),
            ),
            (
                "pr_gene_count",
                Cell::from(count_genes(&|row| {
                    row.candidate_source == crate::hpo_gene_ranking::evidence::PR_SOURCE
                })),
            ),
            (
                "mp_gene_count",
                Cell::from(count_genes(&|row| {
                    row.candidate_source == crate::hpo_gene_ranking::evidence::MP_SOURCE
                })),
            ),
            (
                "orpha_gene_count",
                Cell::from(count_genes(&|row| {
                    row.candidate_source == crate::hpo_gene_ranking::evidence::ORPHA_SOURCE
                })),
            ),
            ("term_gene_set_size_for_scoring", Cell::Int(scoring_size)),
            ("term_gene_breadth_tier", Cell::from(tier)),
            ("term_candidate_role", Cell::from(role)),
            ("term_weight", Cell::Float(weight)),
            (
                "all_gene_symbols",
                Cell::from(joined_set(
                    gene_rows.iter().map(|row| row.gene_symbol.as_str()),
                    false,
                )),
            ),
        ]);
    }

    let mut pairs: Vec<ScoredPair> = Vec::new();
    for (i, hp_a) in frontier.iter().enumerate() {
        for hp_b in &frontier[i + 1..] {
            let summary = ontology.pair_summary(hp_a, hp_b);
            let ancestor = summary.nearest_common_ancestor_id.clone();
            let specificity = ontology.graph_specificity(&ancestor);
            pairs.push(ScoredPair {
                summary,
                specificity,
                specificity_source: "hpo_graph_descendant_leaf",
            });
        }
    }
    let pair_records: Vec<Record> = pairs
        .iter()
        .map(|pair| {
            let ancestor = &pair.summary.nearest_common_ancestor_id;
            let mut record: Record = vec![("sample_id", Cell::from(sample_id))];
            record.extend(pair.summary.record.iter().cloned());
            record.push((
                "nearest_common_ancestor_gene_count",
                Cell::Int(store.hpo_gene_count(ancestor)),
            ));
            record.push((
                "nearest_common_ancestor_descendant_leaf_count",
                Cell::from(ontology.descendant_leaf_count(ancestor)),
            ));
            record.push((
                "nearest_common_ancestor_specificity",
                Cell::Float(pair.specificity),
            ));
            record.push((
                "nearest_common_ancestor_specificity_source",
                Cell::from(pair.specificity_source),
            ));
            record.push((
                "distance_specificity_factor",
                Cell::Float(round(
                    pair.summary.shortest_distance as f64 * (1.0 + pair.specificity),
                    6,
                )),
            ));
            record
        })
        .collect();

    let evidence_records: Vec<Record> = all_rows
        .iter()
        .map(|row| evidence_record(sample_id, row))
        .collect();
    let summary = build_gene_summary(sample_id, &all_rows);
    let immune_hp_ids: HashSet<String> = frontier
        .iter()
        .filter(|id| ontology.is_immune_p_term(id))
        .cloned()
        .collect();
    let (ranked, contributions) = build_gene_ranking(
        sample_id,
        &all_rows,
        &pairs,
        &immune_hp_ids,
        settings.distance_reference,
        external_immune_genes,
    );
    let distinct_p: HashSet<&String> = direct_p_ids.iter().collect();
    let run_summary: Record = vec![
        ("sample_id", Cell::from(sample_id)),
        ("direct_terms", Cell::from(direct_ids.len())),
        ("direct_p_terms", Cell::from(distinct_p.len())),
        ("frontier_p_terms", Cell::from(frontier.len())),
        ("pair_rows", Cell::from(pair_records.len())),
        ("gene_evidence_rows", Cell::from(evidence_records.len())),
        ("gene_summary_rows", Cell::from(summary.len())),
    ];
    SampleTables {
        frontier: frontier_records,
        pairs: pair_records,
        evidence: evidence_records,
        summary,
        ranked,
        contributions,
        run_summary,
    }
}
