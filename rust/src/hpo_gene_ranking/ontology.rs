//! The HPO graph as `patient_frontier_gene_evidence.py` uses it
//! (`HpoOntology.from_exports`): term names from the nodes table, parent and
//! child links from the edges table, and the depth, frontier, common-ancestor
//! and specificity questions asked of it.

use crate::hpo_gene_ranking::Result;
use crate::hpo_gene_ranking::python_text::round;
use crate::hpo_gene_ranking::tsv::Table;
use std::collections::{HashMap, HashSet, VecDeque};
use std::sync::Mutex;

pub const ALL_ROOT: &str = "HP:0000001";
pub const P_ROOT: &str = "HP:0000118";
pub const IMMUNE_ROOT: &str = "HP:0002715";

pub struct Ontology {
    /// Every term id that appears in the nodes or edges table.
    pub ids: Vec<String>,
    index: HashMap<String, usize>,
    /// Name per term; `None` for ids only seen in the edges table.
    names: Vec<Option<String>>,
    /// Parents and children, each sorted by id text without repeats.
    parents: Vec<Vec<usize>>,
    children: Vec<Vec<usize>>,
    is_p: Vec<bool>,
    is_immune: Vec<bool>,
    shortest_p_depth: Vec<i64>,
    longest_p_depth: Vec<i64>,
    top_p_branches: HashSet<usize>,
    p_leaf_count: usize,
    leaf_counts: Mutex<HashMap<usize, usize>>,
}

/// Distances from a term up to each of its ancestors (the term itself at 0).
pub type AncestorDistances = HashMap<usize, u32>;

impl Ontology {
    pub fn from_exports(nodes_path: &str, edges_path: &str) -> Result<Ontology> {
        let (nodes, _) = Table::read(nodes_path, &["hp_id", "name"], &[])?;
        let (edges, _) = Table::read(edges_path, &["child_hp_id", "parent_hp_id"], &[])?;
        let mut ids: Vec<String> = Vec::new();
        let mut index: HashMap<String, usize> = HashMap::new();
        let mut intern = |id: &str, ids: &mut Vec<String>| -> usize {
            if let Some(&position) = index.get(id) {
                return position;
            }
            ids.push(id.to_owned());
            index.insert(id.to_owned(), ids.len() - 1);
            ids.len() - 1
        };
        let mut named: Vec<(usize, String)> = Vec::new();
        for row in &nodes.rows {
            let position = intern(&row[0], &mut ids);
            named.push((position, row[1].clone()));
        }
        let mut links: Vec<(usize, usize)> = Vec::new();
        for row in &edges.rows {
            let child = intern(&row[0], &mut ids);
            let parent = intern(&row[1], &mut ids);
            links.push((child, parent));
        }
        drop(intern);
        let count = ids.len();
        let mut names: Vec<Option<String>> = vec![None; count];
        // dict(zip(hp_id, name)): a repeated id keeps its last name.
        for (position, name) in named {
            names[position] = Some(name);
        }
        let mut parents: Vec<Vec<usize>> = vec![Vec::new(); count];
        let mut children: Vec<Vec<usize>> = vec![Vec::new(); count];
        for (child, parent) in links {
            parents[child].push(parent);
            children[parent].push(child);
        }
        for list in parents.iter_mut().chain(children.iter_mut()) {
            list.sort_by(|a, b| ids[*a].cmp(&ids[*b]));
            list.dedup();
        }
        let mut ontology = Ontology {
            ids,
            index,
            names,
            parents,
            children,
            is_p: vec![false; count],
            is_immune: vec![false; count],
            shortest_p_depth: vec![-1; count],
            longest_p_depth: vec![-1; count],
            top_p_branches: HashSet::new(),
            p_leaf_count: 0,
            leaf_counts: Mutex::new(HashMap::new()),
        };
        ontology.precompute();
        Ok(ontology)
    }

    fn precompute(&mut self) {
        let count = self.ids.len();
        // A term is a P term when it has a name and HP:0000118 is among its
        // ancestors (itself included); the same holds for the immune root.
        let below = |root: Option<usize>, children: &Vec<Vec<usize>>| -> (Vec<bool>, Vec<i64>) {
            let mut reached = vec![false; count];
            let mut distance = vec![-1i64; count];
            if let Some(root) = root {
                reached[root] = true;
                distance[root] = 0;
                let mut queue = VecDeque::from([root]);
                while let Some(current) = queue.pop_front() {
                    for &child in &children[current] {
                        if !reached[child] {
                            reached[child] = true;
                            distance[child] = distance[current] + 1;
                            queue.push_back(child);
                        }
                    }
                }
            }
            (reached, distance)
        };
        let p_root = self.index.get(P_ROOT).copied();
        let immune_root = self.index.get(IMMUNE_ROOT).copied();
        let (below_p, p_distance) = below(p_root, &self.children);
        let (below_immune, _) = below(immune_root, &self.children);
        for term in 0..count {
            let named = self.names[term].is_some();
            self.is_p[term] = named && below_p[term];
            self.is_immune[term] = named && below_immune[term];
            if self.is_p[term] {
                self.shortest_p_depth[term] = p_distance[term];
            }
        }
        if let Some(root) = p_root {
            self.shortest_p_depth[root] = 0;
        }
        let mut memo: Vec<Option<i64>> = vec![None; count];
        for term in 0..count {
            self.longest_depth(term, p_root, &mut memo);
        }
        for term in 0..count {
            self.longest_p_depth[term] = memo[term].unwrap_or(-1);
        }
        match p_root {
            Some(root) => {
                self.top_p_branches = self.children[root].iter().copied().collect();
                // p_leaf_count: descendants of the root (named or not) without P children.
                self.p_leaf_count = (0..count)
                    .filter(|&term| below_p[term])
                    .filter(|&term| !self.children[term].iter().any(|&child| self.is_p[child]))
                    .count();
            }
            // Python counts the absent root itself as one leaf.
            None => self.p_leaf_count = 1,
        }
    }

    fn longest_depth(
        &self,
        term: usize,
        p_root: Option<usize>,
        memo: &mut Vec<Option<i64>>,
    ) -> i64 {
        if let Some(value) = memo[term] {
            return value;
        }
        let value = if Some(term) == p_root {
            0
        } else if !self.is_p[term] {
            -1
        } else {
            let mut best = -1i64;
            for &parent in &self.parents[term] {
                if self.is_p[parent] {
                    let depth = self.longest_depth(parent, p_root, memo);
                    if depth >= 0 {
                        best = best.max(depth);
                    }
                }
            }
            if best >= 0 { best + 1 } else { -1 }
        };
        memo[term] = Some(value);
        value
    }

    pub fn position(&self, id: &str) -> Option<usize> {
        self.index.get(id).copied()
    }

    /// `ontology.names.get(id, "")`.
    pub fn name_of(&self, id: &str) -> &str {
        self.position(id)
            .and_then(|position| self.names[position].as_deref())
            .unwrap_or("")
    }

    pub fn name_at(&self, term: usize) -> &str {
        self.names[term].as_deref().unwrap_or("")
    }

    pub fn has_name(&self, id: &str) -> bool {
        self.position(id)
            .is_some_and(|position| self.names[position].is_some())
    }

    pub fn is_p_term(&self, id: &str) -> bool {
        self.position(id)
            .is_some_and(|position| self.is_p[position])
    }

    pub fn is_immune_p_term(&self, id: &str) -> bool {
        self.position(id)
            .is_some_and(|position| self.is_immune[position])
    }

    pub fn shortest_p_depth(&self, id: &str) -> i64 {
        if id == P_ROOT {
            return 0;
        }
        self.position(id)
            .map_or(-1, |position| self.shortest_p_depth[position])
    }

    pub fn longest_p_depth(&self, id: &str) -> i64 {
        if id == P_ROOT {
            return 0;
        }
        self.position(id)
            .map_or(-1, |position| self.longest_p_depth[position])
    }

    /// `ancestor_distances(hp_id)`: empty for an id without a name.
    pub fn ancestor_distances(&self, id: &str) -> AncestorDistances {
        let mut distances = AncestorDistances::new();
        let Some(start) = self.position(id) else {
            return distances;
        };
        if self.names[start].is_none() {
            return distances;
        }
        distances.insert(start, 0);
        let mut queue = VecDeque::from([start]);
        while let Some(current) = queue.pop_front() {
            let next = distances[&current] + 1;
            for &parent in &self.parents[current] {
                let better = distances.get(&parent).is_none_or(|&known| next < known);
                if better {
                    distances.insert(parent, next);
                    queue.push_back(parent);
                }
            }
        }
        distances
    }

    /// `p_ancestors_above_depth_floor`: P ancestors (breadth first from the
    /// sorted parents) whose shortest P depth is at least `floor`, deepest
    /// first, ties kept in the order found.
    pub fn p_ancestors_above_depth_floor(&self, id: &str, floor: i64) -> Vec<(usize, i64)> {
        let mut result = Vec::new();
        if floor < 0 || !self.has_name(id) {
            return result;
        }
        let start = self.position(id).expect("a named term has a position");
        let mut visited: HashSet<usize> = HashSet::from([start]);
        let mut queue: VecDeque<usize> = self.parents[start].iter().copied().collect();
        while let Some(current) = queue.pop_front() {
            if !visited.insert(current) {
                continue;
            }
            if !self.is_p[current] {
                continue;
            }
            let depth = self.shortest_p_depth(&self.ids[current]);
            if depth < floor {
                continue;
            }
            result.push((current, depth));
            for &parent in &self.parents[current] {
                if !visited.contains(&parent) {
                    queue.push_back(parent);
                }
            }
        }
        result.sort_by_key(|&(_, depth)| -depth);
        result
    }

    /// The first-level P branches among a term's ancestors, sorted by id.
    pub fn top_p_branches_of(&self, distances: &AncestorDistances) -> Vec<usize> {
        let mut branches: Vec<usize> = distances
            .keys()
            .copied()
            .filter(|term| self.top_p_branches.contains(term))
            .collect();
        branches.sort_by(|a, b| self.ids[*a].cmp(&self.ids[*b]));
        branches
    }

    /// `descendant_leaf_count`: P descendants (the term included) that have no P child.
    pub fn descendant_leaf_count(&self, id: &str) -> usize {
        let Some(start) = self.position(id) else {
            return 0;
        };
        if !self.is_p[start] {
            return 0;
        }
        if let Some(&known) = self.leaf_counts.lock().expect("lock").get(&start) {
            return known;
        }
        let mut leaves = 0usize;
        let mut seen: HashSet<usize> = HashSet::from([start]);
        let mut queue = VecDeque::from([start]);
        while let Some(current) = queue.pop_front() {
            let mut has_p_child = false;
            for &child in &self.children[current] {
                if self.is_p[child] {
                    has_p_child = true;
                    if seen.insert(child) {
                        queue.push_back(child);
                    }
                }
            }
            if !has_p_child {
                leaves += 1;
            }
        }
        self.leaf_counts.lock().expect("lock").insert(start, leaves);
        leaves
    }

    /// `graph_specificity`: 1 for a P leaf, falling with the number of P leaves below.
    pub fn graph_specificity(&self, id: &str) -> f64 {
        let total = self.p_leaf_count;
        let leaves = self.descendant_leaf_count(id);
        if total <= 1 || leaves == 0 {
            return 0.0;
        }
        let total = total as f64;
        let score = -(leaves as f64 / total).ln() / -(1.0 / total).ln();
        // max(0.0, min(1.0, score)): Python keeps the first argument unless the
        // second is strictly smaller (min) or larger (max), so -0.0 becomes 0.0.
        let capped = if score < 1.0 { score } else { 1.0 };
        let floored = if capped > 0.0 { capped } else { 0.0 };
        round(floored, 6)
    }

    /// `frontier_terms`: the distinct P terms, sorted, minus those that are an
    /// ancestor of another listed P term.
    pub fn frontier_terms(&self, ids: &[String]) -> Vec<String> {
        let mut direct: Vec<String> = ids
            .iter()
            .filter(|id| self.is_p_term(id))
            .cloned()
            .collect();
        direct.sort();
        direct.dedup();
        let ancestors: Vec<AncestorDistances> = direct
            .iter()
            .map(|id| self.ancestor_distances(id))
            .collect();
        let positions: Vec<usize> = direct
            .iter()
            .map(|id| self.position(id).expect("a P term has a position"))
            .collect();
        let mut frontier = Vec::new();
        for (index, id) in direct.iter().enumerate() {
            let covered = (0..direct.len())
                .any(|other| other != index && ancestors[other].contains_key(&positions[index]));
            if !covered {
                frontier.push(id.clone());
            }
        }
        frontier
    }
}

/// `pair_summary` of two frontier terms (plus the fields that
/// `annotate_pair_specificity` adds later).
pub struct PairSummary {
    pub term1_id: String,
    pub term2_id: String,
    pub shortest_distance: i64,
    pub nearest_common_ancestor_id: String,
    pub cross_module: bool,
    pub cross_module_pair_key: String,
    pub shared_top_p_branches: String,
    pub record: crate::hpo_gene_ranking::tsv::Record,
}

fn joined_terms(ontology: &Ontology, terms: &[usize]) -> String {
    terms
        .iter()
        .map(|&term| format!("{} {}", ontology.ids[term], ontology.name_at(term)))
        .collect::<Vec<_>>()
        .join(" | ")
}

impl Ontology {
    /// Only the shortest distance of `pair_summary` (for the distance reference).
    pub fn pair_shortest_distance(&self, a: &str, b: &str) -> i64 {
        let distances_a = self.ancestor_distances(a);
        let distances_b = self.ancestor_distances(b);
        self.common_ancestors(&distances_a, &distances_b)
            .iter()
            .map(|term| (distances_a[term] + distances_b[term]) as i64)
            .min()
            .expect("two P terms share HP:0000118")
    }

    fn common_ancestors(&self, a: &AncestorDistances, b: &AncestorDistances) -> Vec<usize> {
        let common: Vec<usize> = a
            .keys()
            .copied()
            .filter(|term| b.contains_key(term))
            .collect();
        let common_p: Vec<usize> = common
            .iter()
            .copied()
            .filter(|&term| self.is_p[term])
            .collect();
        if common_p.is_empty() {
            common
        } else {
            common_p
        }
    }

    pub fn pair_summary(&self, a: &str, b: &str) -> PairSummary {
        use crate::hpo_gene_ranking::tsv::Cell;
        let distances_a = self.ancestor_distances(a);
        let distances_b = self.ancestor_distances(b);
        let common = self.common_ancestors(&distances_a, &distances_b);
        let sum = |term: usize| distances_a[&term] + distances_b[&term];
        let shortest = common
            .iter()
            .map(|&term| sum(term))
            .min()
            .expect("a shared ancestor");
        let mut nearest: Vec<usize> = common
            .iter()
            .copied()
            .filter(|&term| sum(term) == shortest)
            .collect();
        nearest.sort_by(|x, y| {
            (-self.longest_p_depth[*x], &self.ids[*x])
                .cmp(&(-self.longest_p_depth[*y], &self.ids[*y]))
        });
        let uninformative = |term: usize| self.ids[term] == ALL_ROOT || self.ids[term] == P_ROOT;
        let nearest_informative: Vec<usize> = nearest
            .iter()
            .copied()
            .filter(|&term| !uninformative(term))
            .collect();
        let selected_nearest = nearest_informative.first().copied().unwrap_or(nearest[0]);
        let informative_common: Vec<usize> = common
            .iter()
            .copied()
            .filter(|&term| !uninformative(term))
            .collect();
        let mica_pool = if informative_common.is_empty() {
            common.clone()
        } else {
            informative_common
        };
        let mica_depth = mica_pool
            .iter()
            .map(|&term| self.longest_p_depth[term])
            .max()
            .expect("a pool");
        let mut mica: Vec<usize> = mica_pool
            .iter()
            .copied()
            .filter(|&term| self.longest_p_depth[term] == mica_depth)
            .collect();
        mica.sort_by(|x, y| (sum(*x), &self.ids[*x]).cmp(&(sum(*y), &self.ids[*y])));
        let top_a = self.top_p_branches_of(&distances_a);
        let top_b = self.top_p_branches_of(&distances_b);
        let shared: Vec<usize> = top_a
            .iter()
            .copied()
            .filter(|term| top_b.contains(term))
            .collect();
        let same_top = !shared.is_empty();
        let ids_of = |terms: &[usize]| {
            terms
                .iter()
                .map(|&t| self.ids[t].clone())
                .collect::<Vec<_>>()
                .join(" | ")
        };
        let mut signatures = [ids_of(&top_a), ids_of(&top_b)];
        signatures.sort();
        let tied = |terms: &[usize]| {
            terms
                .iter()
                .map(|&term| {
                    format!(
                        "{} {} ({}+{})",
                        self.ids[term],
                        self.name_at(term),
                        distances_a[&term],
                        distances_b[&term]
                    )
                })
                .collect::<Vec<_>>()
                .join(" | ")
        };
        let mica_text = mica
            .iter()
            .map(|&term| {
                format!(
                    "{} {} depth={}",
                    self.ids[term],
                    self.name_at(term),
                    self.longest_p_depth[term]
                )
            })
            .collect::<Vec<_>>()
            .join(" | ");
        let cross_module_pair_key = if same_top {
            String::new()
        } else {
            signatures.join(" <> ")
        };
        let shared_text = joined_terms(self, &shared);
        let record: crate::hpo_gene_ranking::tsv::Record = vec![
            ("term1_id", Cell::from(a)),
            ("term1_name", Cell::from(self.name_of(a))),
            (
                "term1_shortest_p_depth",
                Cell::Int(self.shortest_p_depth(a)),
            ),
            ("term1_longest_p_depth", Cell::Int(self.longest_p_depth(a))),
            ("term2_id", Cell::from(b)),
            ("term2_name", Cell::from(self.name_of(b))),
            (
                "term2_shortest_p_depth",
                Cell::Int(self.shortest_p_depth(b)),
            ),
            ("term2_longest_p_depth", Cell::Int(self.longest_p_depth(b))),
            ("shortest_distance", Cell::Int(shortest as i64)),
            (
                "nearest_common_ancestor_id",
                Cell::from(self.ids[selected_nearest].as_str()),
            ),
            (
                "nearest_common_ancestor_name",
                Cell::from(self.name_at(selected_nearest)),
            ),
            ("all_tied_common_ancestors", Cell::from(tied(&nearest))),
            (
                "nearest_informative_common_ancestors",
                Cell::from(tied(&nearest_informative)),
            ),
            ("mica_id", Cell::from(self.ids[mica[0]].as_str())),
            ("mica_name", Cell::from(self.name_at(mica[0]))),
            ("mica_tied_common_ancestors", Cell::from(mica_text)),
            (
                "same_top_P_branch",
                Cell::from(if same_top { "true" } else { "false" }),
            ),
            (
                "term1_top_P_branches",
                Cell::from(joined_terms(self, &top_a)),
            ),
            (
                "term2_top_P_branches",
                Cell::from(joined_terms(self, &top_b)),
            ),
            ("shared_top_P_branches", Cell::from(shared_text.clone())),
            (
                "cross_module",
                Cell::from(if same_top { "false" } else { "true" }),
            ),
            (
                "cross_module_pair_key",
                Cell::from(cross_module_pair_key.clone()),
            ),
        ];
        PairSummary {
            term1_id: a.to_owned(),
            term2_id: b.to_owned(),
            shortest_distance: shortest as i64,
            nearest_common_ancestor_id: self.ids[selected_nearest].clone(),
            cross_module: !same_top,
            cross_module_pair_key,
            shared_top_p_branches: shared_text,
            record,
        }
    }
}
