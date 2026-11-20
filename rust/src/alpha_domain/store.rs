//! Model-specific flattening. Ordinals, not Parquet physical row sorting,
//! define the source's complete gene and nested member iteration order.
use super::model::{
    Distribution, Member, Minimum, Model, Node, Residues, Samples, Scores, Transcript,
};
use super::parquet::{self, Row, Writer};
use super::{Result, exception, unsupported};
use std::collections::BTreeSet;
use std::path::Path;

fn base(
    gene: &str,
    gene_order: u64,
    path: &[String],
    member_order: usize,
    key: &str,
    kind: &str,
) -> Row {
    Row {
        gene: gene.into(),
        gene_order,
        path: path.to_vec(),
        member_order: member_order as u64,
        key: key.into(),
        kind: kind.into(),
        ..Row::default()
    }
}
fn write_node(
    writer: &mut Writer,
    gene: &str,
    ordinal: u64,
    path: &mut Vec<String>,
    node: &Node,
    raw: bool,
) -> Result<()> {
    for (order, (key, value)) in node.members.iter().enumerate() {
        match value {
            Member::Node(child) => {
                path.push(key.clone());
                writer.push(base(
                    gene,
                    ordinal,
                    path,
                    order,
                    key,
                    if path.len() == 1 { "source" } else { "domain" },
                ))?;
                write_node(writer, gene, ordinal, path, child, raw)?;
                path.pop();
            }
            Member::Transcript(transcript) => {
                let mut row = base(gene, ordinal, path, order, key, "transcript");
                row.keys = transcript.exons.iter().cloned().collect();
                match transcript.minimum {
                    Minimum::Unknown => row.minimum_float = Some(f64::INFINITY),
                    Minimum::Position(value) => row.minimum_int = Some(value),
                }
                row.maximum = Some(transcript.maximum);
                writer.push(row)?;
            }
            Member::Scores(mode, scores) => {
                if mode.name() != key {
                    return Err(unsupported("distribution member name/type differ"));
                }
                let mut row = base(
                    gene,
                    ordinal,
                    path,
                    order,
                    key,
                    if raw { "scores" } else { "array" },
                );
                if raw {
                    row.typed_keys = scores.items.iter().map(|(k, _)| k.clone()).collect();
                }
                row.values = scores.items.iter().map(|(_, v)| *v).collect();
                writer.push(row)?;
            }
            Member::Array(mode, values) => {
                if raw || mode.name() != key {
                    return Err(exception(
                        "ValueError",
                        "finalized arrays cannot become raw merge checkpoints",
                    ));
                }
                let mut row = base(gene, ordinal, path, order, key, "array");
                row.values = values.clone();
                writer.push(row)?;
            }
            Member::Residues(residues) => {
                writer.push(base(gene, ordinal, path, order, key, "residues"))?;
                for (item_order, (key, samples)) in residues.items.iter().enumerate() {
                    let mut row = base(gene, ordinal, path, order, "", "residue");
                    row.typed_keys.push(key.clone());
                    row.item_order = item_order as u64;
                    samples.visit(|value, count| {
                        row.values.push(value);
                        row.counts.push(count);
                    });
                    writer.push(row)?;
                }
            }
        }
    }
    Ok(())
}
pub(super) fn write_gene(
    writer: &mut Writer,
    gene: &str,
    ordinal: u64,
    node: &Node,
    raw: bool,
) -> Result<()> {
    writer.push(base(gene, ordinal, &[], 0, gene, "gene"))?;
    write_node(writer, gene, ordinal, &mut vec![], node, raw)
}
pub(super) fn write_model(writer: &mut Writer, model: &Model, raw: bool) -> Result<()> {
    for (ordinal, (gene, value)) in model.genes.members.iter().enumerate() {
        let Member::Node(node) = value else {
            return Err(unsupported("gene value is not a dictionary"));
        };
        write_gene(writer, gene, ordinal as u64, node, raw)?;
    }
    Ok(())
}

fn parent<'a>(node: &'a mut Node, path: &[String]) -> Result<&'a mut Node> {
    let mut node = node;
    for key in path {
        node = node.existing_child(key)?;
    }
    Ok(node)
}
fn insert(node: &mut Node, row: Row, raw: bool) -> Result<()> {
    if row.kind == "source" || row.kind == "domain" {
        if row.path.last() != Some(&row.key) {
            return Err(exception("ValueError", "domain path/key differ"));
        }
        let node = parent(node, &row.path[..row.path.len() - 1])?;
        if node.members.len() as u64 != row.member_order {
            return Err(exception(
                "ValueError",
                "domain member order is discontinuous",
            ));
        }
        return node.append_unique(row.key, Member::Node(Node::default()));
    }
    let node = parent(node, &row.path)?;
    if row.kind == "residue" {
        let (_, Member::Residues(residues)) = node
            .members
            .get_mut(row.member_order as usize)
            .ok_or_else(|| exception("ValueError", "residue parent missing"))?
        else {
            return Err(exception("ValueError", "residue parent has another type"));
        };
        if row.item_order != residues.items.len() as u64
            || row.counts.len() != row.values.len()
            || row.values.is_empty()
            || row.typed_keys.len() != 1
            || !row.key.is_empty()
            || !row.keys.is_empty()
        {
            return Err(exception("ValueError", "invalid residue sequence/order"));
        }
        let mut samples = None;
        for (value, count) in row.values.into_iter().zip(row.counts) {
            if count == 0 {
                return Err(exception("ValueError", "zero residue repeat count"));
            }
            match &mut samples {
                Some(samples) => Samples::push_run(samples, value, count)?,
                None => samples = Some(Samples::One(value, count)),
            }
        }
        return residues
            .append_unique(row.typed_keys.into_iter().next().unwrap(), samples.unwrap());
    }
    if node.members.len() as u64 != row.member_order {
        return Err(exception("ValueError", "member order is discontinuous"));
    }
    let value = match row.kind.as_str() {
        "transcript" => {
            let exons: BTreeSet<_> = row.keys.iter().cloned().collect();
            if exons.len() != row.keys.len() || !row.keys.windows(2).all(|v| v[0] < v[1]) {
                return Err(exception(
                    "ValueError",
                    "exon set must be sorted and unique",
                ));
            }
            let minimum = match (row.minimum_int, row.minimum_float) {
                (Some(value), None) => Minimum::Position(value),
                (None, Some(value)) if value == f64::INFINITY => Minimum::Unknown,
                _ => return Err(exception("ValueError", "invalid typed transcript minimum")),
            };
            Member::Transcript(Transcript {
                exons,
                minimum,
                maximum: row
                    .maximum
                    .ok_or_else(|| exception("ValueError", "transcript maximum absent"))?,
            })
        }
        "scores" if raw => {
            let mode = Distribution::from_name(&row.key)
                .ok_or_else(|| exception("ValueError", "unknown distribution name"))?;
            if row.typed_keys.len() != row.values.len() || !row.keys.is_empty() {
                return Err(exception("ValueError", "raw score keys/values differ"));
            }
            let mut scores = Scores::default();
            for (key, value) in row.typed_keys.into_iter().zip(row.values) {
                scores.append_unique(key, value)?;
            }
            Member::Scores(mode, scores)
        }
        "array" if !raw => Member::Array(
            Distribution::from_name(&row.key)
                .ok_or_else(|| exception("ValueError", "unknown distribution name"))?,
            row.values,
        ),
        "residues" if row.key == "gnomAD_residue_distribution" => {
            Member::Residues(Residues::default())
        }
        _ => {
            return Err(exception(
                "ValueError",
                "unknown or stage-incompatible model row",
            ));
        }
    };
    node.append_unique(row.key, value)
}

/// A single gene at a time; a whole-genome model is never required for merge.
pub(super) fn scan_genes(
    path: &Path,
    mut visit: impl FnMut(String, u64, Node) -> Result<()>,
) -> Result<()> {
    let metadata = parquet::metadata(path)?;
    let raw = match metadata.get("priva.stage").map(String::as_str) {
        Some("raw") => true,
        Some("final") => false,
        _ => return Err(exception("ValueError", "invalid model stage")),
    };
    let mut pending: Option<(String, u64, Node)> = None;
    let mut seen = BTreeSet::new();
    parquet::scan(path, |row| {
        if row.kind == "gene" {
            if let Some((gene, ordinal, node)) = pending.take() {
                visit(gene, ordinal, node)?;
            }
            if row.key != row.gene || !row.path.is_empty() || !seen.insert(row.gene.clone()) {
                return Err(exception("ValueError", "duplicate/invalid gene marker"));
            }
            pending = Some((row.gene, row.gene_order, Node::default()));
            return Ok(());
        }
        let (gene, ordinal, node) = pending
            .as_mut()
            .ok_or_else(|| exception("ValueError", "model row precedes gene marker"))?;
        if *gene != row.gene || *ordinal != row.gene_order {
            return Err(exception("ValueError", "model row gene identity differs"));
        }
        insert(node, row, raw)
    })?;
    if let Some((gene, ordinal, node)) = pending {
        visit(gene, ordinal, node)?;
    }
    Ok(())
}
