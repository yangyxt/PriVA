//! `build_external_immune_genes`: genes in Reactome's immune system pathway
//! (R-HSA-168256 and every pathway below it) or annotated in gene2go to the
//! GO immune system process or defense response terms or their descendants.

use crate::hpo_gene_ranking::Result;
use crate::hpo_gene_ranking::python_text::strip;
use crate::hpo_gene_ranking::tsv::{Table, for_each_line};
use std::collections::{HashMap, HashSet, VecDeque};

const REACTOME_IMMUNE_ROOT: &str = "R-HSA-168256";
const GO_IMMUNE_ROOTS: [&str; 2] = ["GO:0002376", "GO:0006952"];

/// `line.split(separator)[1]`: the text between the first and the second
/// occurrence of the separator.
fn second_piece<'a>(line: &'a str, separator: &str) -> &'a str {
    let after = line.split_once(separator).map_or("", |(_, rest)| rest);
    after
        .split_once(separator)
        .map_or(after, |(piece, _)| piece)
}

/// `text.split(" !")[0].strip()`.
fn before_comment(text: &str) -> &str {
    strip(text.split_once(" !").map_or(text, |(piece, _)| piece))
}

pub fn external_immune_genes(
    reactome_genes_tsv: &str,
    reactome_relations_txt: &str,
    gene2go_txt: &str,
    go_obo_txt: &str,
    ncbi_to_symbol: &HashMap<String, String>,
) -> Result<HashSet<String>> {
    let mut symbols = HashSet::new();

    // Reactome: breadth first from the immune system pathway.
    let mut children_of: HashMap<String, Vec<String>> = HashMap::new();
    for_each_line(reactome_relations_txt, true, |line| {
        let line = line.strip_suffix('\n').unwrap_or(line);
        if line.is_empty() {
            return Ok(());
        }
        let mut fields = line.split('\t');
        let parent = fields.next().unwrap_or("");
        let child = fields.next().unwrap_or("");
        children_of
            .entry(parent.to_owned())
            .or_default()
            .push(child.to_owned());
        Ok(())
    })?;
    let mut immune_pathways: HashSet<String> = HashSet::from([REACTOME_IMMUNE_ROOT.to_owned()]);
    let mut queue = VecDeque::from([REACTOME_IMMUNE_ROOT.to_owned()]);
    while let Some(pathway) = queue.pop_front() {
        for child in children_of.get(&pathway).into_iter().flatten() {
            if !immune_pathways.contains(child) && child.starts_with("R-HSA-") {
                immune_pathways.insert(child.clone());
                queue.push_back(child.clone());
            }
        }
    }
    let (reactome, _) = Table::read(reactome_genes_tsv, &["reactome_id", "ncbi_gene_id"], &[])?;
    for row in &reactome.rows {
        if immune_pathways.contains(&row[0])
            && let Some(symbol) = ncbi_to_symbol.get(&row[1])
            && !symbol.is_empty()
        {
            symbols.insert(symbol.clone());
        }
    }

    // GO: breadth first from the immune roots through is_a and part_of.
    let mut go_children: HashMap<String, HashSet<String>> = HashMap::new();
    let mut current: Option<String> = None;
    for_each_line(go_obo_txt, true, |line| {
        let line = strip(line);
        if line.starts_with("id: GO:") {
            current = Some(second_piece(line, "id: ").to_owned());
        } else if let Some(term) = current.as_ref().filter(|term| !term.is_empty()) {
            if line.starts_with("is_a: GO:") {
                let parent = before_comment(second_piece(line, "is_a: "));
                go_children
                    .entry(parent.to_owned())
                    .or_default()
                    .insert(term.clone());
            } else if line.starts_with("relationship: part_of GO:") {
                let parent = before_comment(second_piece(line, "part_of "));
                go_children
                    .entry(parent.to_owned())
                    .or_default()
                    .insert(term.clone());
            }
        }
        Ok(())
    })?;
    let mut go_immune: HashSet<String> = HashSet::new();
    for root in GO_IMMUNE_ROOTS {
        go_immune.insert(root.to_owned());
        let mut queue = VecDeque::from([root.to_owned()]);
        while let Some(term) = queue.pop_front() {
            for child in go_children.get(&term).into_iter().flatten() {
                if go_immune.insert(child.clone()) {
                    queue.push_back(child.clone());
                }
            }
        }
    }
    for_each_line(gene2go_txt, true, |line| {
        if line.starts_with('#') {
            return Ok(());
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 3
            && parts[0] == "9606"
            && go_immune.contains(parts[2])
            && let Some(symbol) = ncbi_to_symbol.get(parts[1])
            && !symbol.is_empty()
        {
            symbols.insert(symbol.clone());
        }
        Ok(())
    })?;
    Ok(symbols)
}
