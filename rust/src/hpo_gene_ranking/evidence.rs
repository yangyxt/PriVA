//! Gene evidence for one frontier term (`EvidenceStore` in
//! `patient_frontier_gene_evidence.py`): direct HPO phenotype-to-gene
//! annotations, genes of parent terms, and GO, PR, MP and ORPHA links of the
//! term in the HPO graph.

use crate::hpo_gene_ranking::Result;
use crate::hpo_gene_ranking::ontology::Ontology;
use crate::hpo_gene_ranking::python_text::{split_first_word, strip};
use crate::hpo_gene_ranking::tsv::Table;
use indexmap::IndexMap;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

pub const HPO_SOURCE: &str = "HPO_PHENOTYPE_TO_GENES";
pub const GOA_EXPERIMENTAL: &str = "HPO_LINKED_GO_TO_GOA_EXPERIMENTAL";
pub const GOA_NON_EXPERIMENTAL: &str = "HPO_LINKED_GO_TO_GOA_NON_EXPERIMENTAL";
pub const PR_SOURCE: &str = "HPO_LINKED_PR_TO_PRO_UNIPROT";
pub const MP_SOURCE: &str = "HPO_XREF_MP_TO_MGI_HUMAN_ORTHOLOG";
pub const ORPHA_SOURCE: &str = "HPO_XREF_ORPHA_TO_DISEASE_GENES";

const EXPERIMENTAL_GOA_CODES: [&str; 11] = [
    "EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP",
];

/// One gene evidence row of one frontier term, with the scoring fields that
/// `process_sample` adds.
#[derive(Clone, Debug)]
pub struct EvidenceRow {
    pub hp_id: String,
    pub hp_name: String,
    pub gene_symbol: String,
    pub ncbi_gene_id: String,
    pub uniprot_accession: String,
    pub candidate_source: &'static str,
    pub candidate_status: &'static str,
    pub external_id: String,
    pub external_name: String,
    pub gene_set_size: i64,
    pub evidence_detail: String,
    pub evidence_row_count: i64,
    pub propagated_from_hp_id: String,
    pub propagated_from_hp_name: String,
    pub source_quality: f64,
    pub term_gene_breadth_tier: &'static str,
    pub term_candidate_role: &'static str,
    pub term_weight: f64,
}

impl EvidenceRow {
    fn new(
        hp_id: &str,
        hp_name: &str,
        candidate_source: &'static str,
        candidate_status: &'static str,
    ) -> Self {
        EvidenceRow {
            hp_id: hp_id.to_owned(),
            hp_name: hp_name.to_owned(),
            gene_symbol: String::new(),
            ncbi_gene_id: String::new(),
            uniprot_accession: String::new(),
            candidate_source,
            candidate_status,
            external_id: String::new(),
            external_name: String::new(),
            gene_set_size: 0,
            evidence_detail: String::new(),
            evidence_row_count: 0,
            propagated_from_hp_id: String::new(),
            propagated_from_hp_name: String::new(),
            source_quality: 1.0,
            term_gene_breadth_tier: "",
            term_candidate_role: "",
            term_weight: 0.0,
        }
    }
}

/// The links of one HP term read from the nodes table.
#[derive(Default)]
struct NodeLinks {
    linked_go_terms: String,
    linked_other_ontology_terms: String,
    xrefs: String,
}

/// Paths of the tables the evidence store reads.
pub struct EvidencePaths {
    pub nodes_tsv: String,
    pub phenotype_to_genes: String,
    pub go_to_genes: String,
    pub pr_to_gene_or_protein: String,
    pub mp_to_human_orthologs: String,
    pub orpha_disease_to_genes: String,
    pub uniprot_evidence: String,
    pub human_gene_info: String,
}

pub struct EvidenceStore {
    nodes: HashMap<String, NodeLinks>,
    /// phenotype_to_genes columns: hpo_id, gene_symbol, ncbi_gene_id, disease_id.
    hpo: Table,
    hpo_has_disease: bool,
    hpo_rows: HashMap<String, Vec<usize>>,
    hpo_gene_count: HashMap<String, i64>,
    /// go_to_genes columns: go_id, gene_symbol, db, db_object_id, evidence_code.
    goa: Table,
    goa_rows: HashMap<String, Vec<usize>>,
    /// pr_to_gene_or_protein columns: pr_id, pr_name, xref_db, uniprot_accession, ncbi_gene_id.
    pr: Table,
    pr_rows: HashMap<String, Vec<usize>>,
    /// mp_to_human_orthologs columns: mp_id, human_gene_symbol, human_entrez_gene_id, mouse_gene_symbol, pubmed_ids.
    mp: Table,
    mp_rows: HashMap<String, Vec<usize>>,
    /// orpha columns: orpha_disease_id, gene_symbol, disease_name, association_status, association_type.
    orpha: Table,
    orpha_rows: HashMap<String, Vec<usize>>,
    symbol_to_ncbi: IndexMap<String, String>,
    ncbi_to_symbol: HashMap<String, String>,
    first_symbol_of_ncbi: HashMap<String, String>,
    uniprot_to_gene: HashMap<String, (String, String)>,
}

/// `_split_linked_terms`: "ID name | ID name" into (id, name) pairs.
fn split_linked_terms(value: &str) -> Vec<(String, String)> {
    value
        .split('|')
        .map(strip)
        .filter(|part| !part.is_empty())
        .map(|part| {
            let (id, name) = split_first_word(part);
            (id.to_owned(), name.to_owned())
        })
        .collect()
}

/// `_split_xrefs`.
fn split_xrefs(value: &str) -> Vec<String> {
    value
        .split('|')
        .map(strip)
        .filter(|part| !part.is_empty())
        .map(str::to_owned)
        .collect()
}

fn joined_sorted<'a>(values: impl Iterator<Item = &'a str>) -> String {
    values
        .filter(|value| !value.is_empty())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>()
        .join("|")
}

impl EvidenceStore {
    pub fn load(paths: &EvidencePaths) -> Result<EvidenceStore> {
        let (node_table, _) = Table::read(
            &paths.nodes_tsv,
            &[
                "hp_id",
                "linked_go_terms",
                "linked_other_ontology_terms",
                "xrefs",
            ],
            &["linked_go_terms", "linked_other_ontology_terms", "xrefs"],
        )?;
        let mut nodes: HashMap<String, NodeLinks> = HashMap::new();
        for row in node_table.rows {
            // `_node_row` takes the first row of an id.
            nodes.entry(row[0].clone()).or_insert(NodeLinks {
                linked_go_terms: row[1].clone(),
                linked_other_ontology_terms: row[2].clone(),
                xrefs: row[3].clone(),
            });
        }
        // Only the GO, PR, MP and ORPHA ids that some HP term links to are
        // ever looked up, so other rows are not kept.
        let mut go_ids: HashSet<String> = HashSet::new();
        let mut other_ids: HashSet<String> = HashSet::new();
        for links in nodes.values() {
            for (id, _) in split_linked_terms(&links.linked_go_terms) {
                go_ids.insert(id);
            }
            for (id, _) in split_linked_terms(&links.linked_other_ontology_terms) {
                other_ids.insert(id);
            }
            for id in split_xrefs(&links.xrefs) {
                other_ids.insert(id);
            }
        }

        let (hpo, hpo_present) = Table::read(
            &paths.phenotype_to_genes,
            &["hpo_id", "gene_symbol", "ncbi_gene_id", "disease_id"],
            &["disease_id"],
        )?;
        let hpo_rows = hpo.index_by(0);
        let mut genes_per_term: HashMap<&str, HashSet<&str>> = HashMap::new();
        for row in &hpo.rows {
            genes_per_term
                .entry(row[0].as_str())
                .or_default()
                .insert(row[1].as_str());
        }
        let hpo_gene_count = genes_per_term
            .into_iter()
            .map(|(term, genes)| (term.to_owned(), genes.len() as i64))
            .collect();

        let mut goa = Table::read(
            &paths.go_to_genes,
            &[
                "go_id",
                "gene_symbol",
                "db",
                "db_object_id",
                "evidence_code",
            ],
            &[],
        )?
        .0;
        goa.rows.retain(|row| go_ids.contains(&row[0]));
        let goa_rows = goa.index_by(0);

        let mut pr = Table::read(
            &paths.pr_to_gene_or_protein,
            &[
                "pr_id",
                "pr_name",
                "xref_db",
                "uniprot_accession",
                "ncbi_gene_id",
            ],
            &[],
        )?
        .0;
        pr.rows.retain(|row| other_ids.contains(&row[0]));
        let pr_rows = pr.index_by(0);

        let mut mp = Table::read(
            &paths.mp_to_human_orthologs,
            &[
                "mp_id",
                "human_gene_symbol",
                "human_entrez_gene_id",
                "mouse_gene_symbol",
                "pubmed_ids",
            ],
            &[],
        )?
        .0;
        mp.rows.retain(|row| other_ids.contains(&row[0]));
        let mp_rows = mp.index_by(0);

        let mut orpha = Table::read(
            &paths.orpha_disease_to_genes,
            &[
                "orpha_disease_id",
                "gene_symbol",
                "disease_name",
                "association_status",
                "association_type",
            ],
            &[],
        )?
        .0;
        orpha.rows.retain(|row| other_ids.contains(&row[0]));
        let orpha_rows = orpha.index_by(0);

        let (uniprot, _) = Table::read(
            &paths.uniprot_evidence,
            &["uniprot_accession", "primary_gene_symbol", "ncbi_gene_ids"],
            &[],
        )?;
        let mut uniprot_to_gene = HashMap::new();
        for row in &uniprot.rows {
            let ids = row[2].trim_matches(';');
            let first = if ids.is_empty() {
                ""
            } else {
                ids.split(';').next().unwrap_or("")
            };
            uniprot_to_gene.insert(row[0].clone(), (row[1].clone(), first.to_owned()));
        }

        let (gene_info, _) = Table::read(
            &paths.human_gene_info,
            &["gene_symbol", "ncbi_gene_id"],
            &[],
        )?;
        let mut symbol_to_ncbi: IndexMap<String, String> = IndexMap::new();
        let mut first_symbol_of_ncbi: HashMap<String, String> = HashMap::new();
        for row in &gene_info.rows {
            symbol_to_ncbi.insert(row[0].clone(), row[1].clone());
            first_symbol_of_ncbi
                .entry(row[1].clone())
                .or_insert_with(|| row[0].clone());
        }
        let mut ncbi_to_symbol = HashMap::new();
        for (symbol, ncbi) in &symbol_to_ncbi {
            ncbi_to_symbol.insert(ncbi.clone(), symbol.clone());
        }

        Ok(EvidenceStore {
            nodes,
            hpo,
            hpo_has_disease: hpo_present[3],
            hpo_rows,
            hpo_gene_count,
            goa,
            goa_rows,
            pr,
            pr_rows,
            mp,
            mp_rows,
            orpha,
            orpha_rows,
            symbol_to_ncbi,
            ncbi_to_symbol,
            first_symbol_of_ncbi,
            uniprot_to_gene,
        })
    }

    /// NCBI gene id to symbol, as `_ncbi_to_symbol` builds it.
    pub fn ncbi_to_symbol(&self) -> &HashMap<String, String> {
        &self.ncbi_to_symbol
    }

    /// Number of distinct genes annotated to a term in phenotype_to_genes.
    pub fn hpo_gene_count(&self, hp_id: &str) -> i64 {
        self.hpo_gene_count.get(hp_id).copied().unwrap_or(0)
    }

    fn hpo_gene_evidence(
        &self,
        hp_id: &str,
        hp_name: &str,
        propagated_from: Option<(&str, &str)>,
    ) -> Vec<EvidenceRow> {
        let Some(rows) = self.hpo_rows.get(hp_id) else {
            return Vec::new();
        };
        let gene_set_size = rows
            .iter()
            .map(|&row| self.hpo.rows[row][1].as_str())
            .collect::<HashSet<_>>()
            .len() as i64;
        // groupby(["gene_symbol", "ncbi_gene_id"]) with sorted keys.
        let mut groups: BTreeMap<(&str, &str), Vec<usize>> = BTreeMap::new();
        for &row in rows {
            let values = &self.hpo.rows[row];
            groups
                .entry((values[1].as_str(), values[2].as_str()))
                .or_default()
                .push(row);
        }
        let status = if propagated_from.is_some() {
            "parent_propagated_gene_candidate"
        } else {
            "direct_responsible_gene_candidate"
        };
        groups
            .into_iter()
            .map(|((gene_symbol, ncbi_gene_id), group)| {
                let mut row = EvidenceRow::new(hp_id, hp_name, HPO_SOURCE, status);
                row.gene_symbol = gene_symbol.to_owned();
                row.ncbi_gene_id = ncbi_gene_id.to_owned();
                row.external_id = hp_id.to_owned();
                row.external_name = hp_name.to_owned();
                row.gene_set_size = gene_set_size;
                if self.hpo_has_disease {
                    let diseases: BTreeSet<&str> = group
                        .iter()
                        .map(|&r| self.hpo.rows[r][3].as_str())
                        .collect();
                    row.evidence_detail = diseases.into_iter().collect::<Vec<_>>().join("|");
                }
                row.evidence_row_count = group.len() as i64;
                if let Some((from_id, from_name)) = propagated_from {
                    row.propagated_from_hp_id = from_id.to_owned();
                    row.propagated_from_hp_name = from_name.to_owned();
                }
                row
            })
            .collect()
    }

    /// The status of a secondary-source gene given the direct genes, or None
    /// when the row is left out.
    fn secondary_status(
        gene_symbol: &str,
        direct: &HashSet<String>,
        include_secondary_always: bool,
        support: &'static str,
        exploratory: &'static str,
        secondary: &'static str,
    ) -> Option<&'static str> {
        if direct.is_empty() {
            Some(secondary)
        } else if direct.contains(gene_symbol) {
            Some(support)
        } else if include_secondary_always {
            Some(exploratory)
        } else {
            None
        }
    }

    fn goa_gene_evidence(
        &self,
        hp_id: &str,
        hp_name: &str,
        direct: &HashSet<String>,
        include_secondary_always: bool,
    ) -> Vec<EvidenceRow> {
        let Some(node) = self.nodes.get(hp_id) else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for (go_id, go_name) in split_linked_terms(&node.linked_go_terms) {
            let Some(rows) = self.goa_rows.get(&go_id) else {
                continue;
            };
            let mut grouped: BTreeMap<(&'static str, String, String), Vec<usize>> = BTreeMap::new();
            for &row in rows {
                let values = &self.goa.rows[row];
                let (raw_symbol, db, object_id) = (&values[1], &values[2], &values[3]);
                let (gene_symbol, ncbi_gene_id) =
                    if let Some(ncbi) = self.symbol_to_ncbi.get(raw_symbol) {
                        (raw_symbol.clone(), ncbi.clone())
                    } else if db == "UniProtKB" && self.uniprot_to_gene.contains_key(object_id) {
                        self.uniprot_to_gene[object_id].clone()
                    } else {
                        continue;
                    };
                if gene_symbol.is_empty() {
                    continue;
                }
                let source = if EXPERIMENTAL_GOA_CODES.contains(&values[4].as_str()) {
                    GOA_EXPERIMENTAL
                } else {
                    GOA_NON_EXPERIMENTAL
                };
                grouped
                    .entry((source, gene_symbol, ncbi_gene_id))
                    .or_default()
                    .push(row);
            }
            let mut genes_per_source: HashMap<&'static str, HashSet<&str>> = HashMap::new();
            for (source, gene_symbol, _) in grouped.keys() {
                genes_per_source
                    .entry(*source)
                    .or_default()
                    .insert(gene_symbol.as_str());
            }
            for ((source, gene_symbol, ncbi_gene_id), group) in &grouped {
                let Some(status) = Self::secondary_status(
                    gene_symbol,
                    direct,
                    include_secondary_always,
                    "functional_support_for_direct_candidate",
                    "exploratory_goa_candidate",
                    "secondary_candidate_from_goa",
                ) else {
                    continue;
                };
                let mut row = EvidenceRow::new(hp_id, hp_name, *source, status);
                row.gene_symbol = gene_symbol.clone();
                row.ncbi_gene_id = ncbi_gene_id.clone();
                // Every UniProtKB object id of the group, an empty one included.
                let accessions: BTreeSet<&str> = group
                    .iter()
                    .filter(|&&r| self.goa.rows[r][2] == "UniProtKB")
                    .map(|&r| self.goa.rows[r][3].as_str())
                    .collect();
                row.uniprot_accession = accessions.into_iter().collect::<Vec<_>>().join("|");
                row.external_id = go_id.clone();
                row.external_name = go_name.clone();
                row.gene_set_size = genes_per_source[source].len() as i64;
                let codes: BTreeSet<&str> = group
                    .iter()
                    .map(|&r| self.goa.rows[r][4].as_str())
                    .collect();
                row.evidence_detail = format!(
                    "evidence_codes={}",
                    codes.into_iter().collect::<Vec<_>>().join("|")
                );
                row.evidence_row_count = group.len() as i64;
                out.push(row);
            }
        }
        out
    }

    fn pr_gene_evidence(
        &self,
        hp_id: &str,
        hp_name: &str,
        direct: &HashSet<String>,
        include_secondary_always: bool,
    ) -> Vec<EvidenceRow> {
        let Some(node) = self.nodes.get(hp_id) else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for (pr_id, pr_name) in split_linked_terms(&node.linked_other_ontology_terms) {
            if !pr_id.starts_with("PR:") {
                continue;
            }
            let Some(rows) = self.pr_rows.get(&pr_id) else {
                continue;
            };
            let distinct = |column: usize| {
                rows.iter()
                    .map(|&r| self.pr.rows[r][column].as_str())
                    .filter(|value| !value.is_empty())
                    .collect::<HashSet<_>>()
                    .len() as i64
            };
            let gene_set_size = distinct(3).max(distinct(4));
            for &r in rows {
                let values = &self.pr.rows[r];
                let accession = &values[3];
                let mut ncbi_gene_id = values[4].clone();
                let mut gene_symbol;
                if !accession.is_empty() {
                    let (symbol, mapped) = self
                        .uniprot_to_gene
                        .get(accession)
                        .cloned()
                        .unwrap_or_default();
                    gene_symbol = symbol;
                    if ncbi_gene_id.is_empty() {
                        ncbi_gene_id = mapped;
                    }
                } else {
                    gene_symbol = String::new();
                    if !ncbi_gene_id.is_empty()
                        && let Some(first) = self.first_symbol_of_ncbi.get(&ncbi_gene_id)
                    {
                        gene_symbol = first.clone();
                    }
                }
                if gene_symbol.is_empty() && !ncbi_gene_id.is_empty() {
                    gene_symbol = self
                        .ncbi_to_symbol
                        .get(&ncbi_gene_id)
                        .cloned()
                        .unwrap_or_default();
                }
                if gene_symbol.is_empty() {
                    continue;
                }
                let Some(status) = Self::secondary_status(
                    &gene_symbol,
                    direct,
                    include_secondary_always,
                    "protein_bridge_support_for_direct_candidate",
                    "exploratory_pr_candidate",
                    "secondary_candidate_from_pr_proximal_molecule",
                ) else {
                    continue;
                };
                let mut row = EvidenceRow::new(hp_id, hp_name, PR_SOURCE, status);
                row.gene_symbol = gene_symbol;
                row.ncbi_gene_id = ncbi_gene_id;
                row.uniprot_accession = accession.clone();
                row.external_id = pr_id.clone();
                row.external_name = pr_name.clone();
                row.gene_set_size = gene_set_size;
                row.evidence_detail = format!("pr_name={};xref_db={}", values[1], values[2]);
                row.evidence_row_count = 1;
                out.push(row);
            }
        }
        out
    }

    /// MP or ORPHA ids among a term's xrefs and linked terms, sorted.
    fn linked_ids_with_prefix(node: &NodeLinks, prefix: &str) -> Vec<String> {
        let mut ids: BTreeSet<String> = BTreeSet::new();
        for id in split_xrefs(&node.xrefs) {
            if id.starts_with(prefix) {
                ids.insert(id);
            }
        }
        for (id, _) in split_linked_terms(&node.linked_other_ontology_terms) {
            if id.starts_with(prefix) {
                ids.insert(id);
            }
        }
        ids.into_iter().collect()
    }

    fn mp_gene_evidence(
        &self,
        hp_id: &str,
        hp_name: &str,
        direct: &HashSet<String>,
        include_secondary_always: bool,
    ) -> Vec<EvidenceRow> {
        let Some(node) = self.nodes.get(hp_id) else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for mp_id in Self::linked_ids_with_prefix(node, "MP:") {
            let Some(rows) = self.mp_rows.get(&mp_id) else {
                continue;
            };
            let mapped: Vec<usize> = rows
                .iter()
                .copied()
                .filter(|&r| !self.mp.rows[r][1].is_empty())
                .collect();
            if mapped.is_empty() {
                continue;
            }
            let gene_set_size = mapped
                .iter()
                .map(|&r| self.mp.rows[r][1].as_str())
                .collect::<HashSet<_>>()
                .len() as i64;
            let mut groups: BTreeMap<(&str, &str), Vec<usize>> = BTreeMap::new();
            for &r in &mapped {
                groups
                    .entry((self.mp.rows[r][1].as_str(), self.mp.rows[r][2].as_str()))
                    .or_default()
                    .push(r);
            }
            for ((gene_symbol, ncbi_gene_id), group) in groups {
                let Some(status) = Self::secondary_status(
                    gene_symbol,
                    direct,
                    include_secondary_always,
                    "mouse_ortholog_support_for_direct_candidate",
                    "exploratory_mp_ortholog_candidate",
                    "secondary_candidate_from_mp_mouse_ortholog",
                ) else {
                    continue;
                };
                let mouse = joined_sorted(group.iter().map(|&r| self.mp.rows[r][3].as_str()));
                let pubmed: Vec<&str> = group
                    .iter()
                    .map(|&r| self.mp.rows[r][4].as_str())
                    .filter(|value| !value.is_empty())
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .take(20)
                    .collect();
                let mut row = EvidenceRow::new(hp_id, hp_name, MP_SOURCE, status);
                row.gene_symbol = gene_symbol.to_owned();
                row.ncbi_gene_id = ncbi_gene_id.to_owned();
                row.external_id = mp_id.clone();
                row.gene_set_size = gene_set_size;
                row.evidence_detail =
                    format!("mouse_genes={mouse};pubmed_ids={}", pubmed.join("|"));
                row.evidence_row_count = group.len() as i64;
                out.push(row);
            }
        }
        out
    }

    fn orpha_gene_evidence(
        &self,
        hp_id: &str,
        hp_name: &str,
        direct: &HashSet<String>,
        include_secondary_always: bool,
    ) -> Vec<EvidenceRow> {
        let Some(node) = self.nodes.get(hp_id) else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for orpha_id in Self::linked_ids_with_prefix(node, "ORPHA:") {
            let Some(rows) = self.orpha_rows.get(&orpha_id) else {
                continue;
            };
            let mapped: Vec<usize> = rows
                .iter()
                .copied()
                .filter(|&r| !self.orpha.rows[r][1].is_empty())
                .collect();
            if mapped.is_empty() {
                continue;
            }
            let gene_set_size = mapped
                .iter()
                .map(|&r| self.orpha.rows[r][1].as_str())
                .collect::<HashSet<_>>()
                .len() as i64;
            let mut groups: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
            for &r in &mapped {
                groups
                    .entry(self.orpha.rows[r][1].as_str())
                    .or_default()
                    .push(r);
            }
            for (gene_symbol, group) in groups {
                let Some(status) = Self::secondary_status(
                    gene_symbol,
                    direct,
                    include_secondary_always,
                    "orpha_disease_support_for_direct_candidate",
                    "exploratory_orpha_disease_candidate",
                    "secondary_candidate_from_orpha_disease",
                ) else {
                    continue;
                };
                let mut row = EvidenceRow::new(hp_id, hp_name, ORPHA_SOURCE, status);
                row.gene_symbol = gene_symbol.to_owned();
                row.ncbi_gene_id = self
                    .symbol_to_ncbi
                    .get(gene_symbol)
                    .cloned()
                    .unwrap_or_default();
                row.external_id = orpha_id.clone();
                row.external_name =
                    joined_sorted(group.iter().map(|&r| self.orpha.rows[r][2].as_str()));
                row.gene_set_size = gene_set_size;
                row.evidence_detail = format!(
                    "association_status={};association_type={}",
                    joined_sorted(group.iter().map(|&r| self.orpha.rows[r][3].as_str())),
                    joined_sorted(group.iter().map(|&r| self.orpha.rows[r][4].as_str())),
                );
                row.evidence_row_count = group.len() as i64;
                out.push(row);
            }
        }
        out
    }

    /// `gene_evidence_for_frontier_term`.
    pub fn gene_evidence_for_frontier_term(
        &self,
        hp_id: &str,
        hp_name: &str,
        include_secondary_always: bool,
        ontology: &Ontology,
        parent_depth_floor: i64,
    ) -> Vec<EvidenceRow> {
        let mut evidence = self.hpo_gene_evidence(hp_id, hp_name, None);
        let mut direct: HashSet<String> =
            evidence.iter().map(|row| row.gene_symbol.clone()).collect();
        if parent_depth_floor >= 0 {
            for (ancestor, _) in ontology.p_ancestors_above_depth_floor(hp_id, parent_depth_floor) {
                let ancestor_id = &ontology.ids[ancestor];
                let ancestor_name = ontology.name_at(ancestor);
                for row in
                    self.hpo_gene_evidence(ancestor_id, ancestor_name, Some((hp_id, hp_name)))
                {
                    if !direct.contains(&row.gene_symbol) {
                        direct.insert(row.gene_symbol.clone());
                        evidence.push(row);
                    }
                }
            }
        }
        evidence.extend(self.goa_gene_evidence(hp_id, hp_name, &direct, include_secondary_always));
        evidence.extend(self.pr_gene_evidence(hp_id, hp_name, &direct, include_secondary_always));
        evidence.extend(self.mp_gene_evidence(hp_id, hp_name, &direct, include_secondary_always));
        evidence.extend(self.orpha_gene_evidence(
            hp_id,
            hp_name,
            &direct,
            include_secondary_always,
        ));
        evidence
    }
}
