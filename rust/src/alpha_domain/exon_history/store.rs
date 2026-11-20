//! Ordered history rows use the existing typed Parquet transport, with a
//! distinct semantic kind/version. They contain no score or residue arrays.
use super::super::model::Minimum;
use super::super::parquet::{self, Row, Writer};
use super::model::{Entry, Model, Node, Transcript};
use super::{Result, invalid};
use std::path::Path;

fn base(gene: &str, ordinal: u64, path: &[String], order: usize, key: &str, kind: &str) -> Row {
    Row {
        gene: gene.into(),
        gene_order: ordinal,
        path: path.to_vec(),
        member_order: order as u64,
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
) -> Result<()> {
    for (order, (key, entry)) in node.members.iter().enumerate() {
        match entry {
            Entry::Branch(child) => {
                path.push(key.clone());
                writer.push(base(gene, ordinal, path, order, key, "history-branch"))?;
                write_node(writer, gene, ordinal, path, child)?;
                path.pop();
            }
            Entry::Transcript(value) => {
                let mut row = base(gene, ordinal, path, order, key, "history-transcript");
                row.keys = value.sequence();
                match value.minimum {
                    Minimum::Unknown => row.minimum_float = Some(f64::INFINITY),
                    Minimum::Position(value) => row.minimum_int = Some(value),
                }
                row.maximum = Some(value.maximum);
                writer.push(row)?;
            }
            Entry::Marker => {
                writer.push(base(gene, ordinal, path, order, key, "history-marker"))?
            }
        }
    }
    Ok(())
}
pub(super) fn write_gene(writer: &mut Writer, gene: &str, ordinal: u64, node: &Node) -> Result<()> {
    writer.push(base(gene, ordinal, &[], 0, gene, "history-gene"))?;
    write_node(writer, gene, ordinal, &mut vec![], node)
}
pub(super) fn write_model(writer: &mut Writer, model: &Model) -> Result<()> {
    for (ordinal, (gene, entry)) in model.genes.members.iter().enumerate() {
        let Entry::Branch(node) = entry else {
            return Err(invalid("history gene is not a dictionary"));
        };
        write_gene(writer, gene, ordinal as u64, node)?;
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
fn unused(row: &Row) -> Result<()> {
    if row.item_order != 0
        || !row.values.is_empty()
        || !row.counts.is_empty()
        || !row.typed_keys.is_empty()
        || row.auxiliary_int.is_some()
    {
        return Err(invalid("history row contains a score/residue payload"));
    }
    if row.kind != "history-transcript"
        && (!row.keys.is_empty()
            || row.minimum_int.is_some()
            || row.minimum_float.is_some()
            || row.maximum.is_some())
    {
        return Err(invalid("non-transcript history row contains exon fields"));
    }
    Ok(())
}
fn insert(node: &mut Node, row: Row) -> Result<()> {
    unused(&row)?;
    let (path, entry) = match row.kind.as_str() {
        "history-branch" => {
            if row.path.last() != Some(&row.key) {
                return Err(invalid("history branch path/key differ"));
            }
            (
                &row.path[..row.path.len() - 1],
                Entry::Branch(Node::default()),
            )
        }
        "history-transcript" => {
            let minimum = match (row.minimum_int, row.minimum_float) {
                (Some(value), None) => Minimum::Position(value),
                (None, Some(value)) if value.to_bits() == f64::INFINITY.to_bits() => {
                    Minimum::Unknown
                }
                _ => return Err(invalid("invalid history minimum type/value")),
            };
            let maximum = row
                .maximum
                .ok_or_else(|| invalid("history maximum absent"))?;
            if row.path.len() < 2 || !row.key.starts_with("ENST") {
                return Err(invalid("invalid history transcript path"));
            }
            (
                &row.path[..],
                Entry::Transcript(Transcript::from_sequence(row.keys, minimum, maximum)?),
            )
        }
        "history-marker" => {
            if row.path.len() < 2
                || ![
                    "distribution",
                    "gnomAD_distribution",
                    "gnomAD_residue_distribution",
                    "min_distribution",
                    "max_distribution",
                ]
                .contains(&row.key.as_str())
            {
                return Err(invalid("invalid history distribution marker"));
            }
            (&row.path[..], Entry::Marker)
        }
        _ => return Err(invalid("unknown history row kind")),
    };
    let node = parent(node, path)?;
    if node.members.len() as u64 != row.member_order {
        return Err(invalid("history member order is discontinuous"));
    }
    node.append(row.key, entry)
}
pub(super) fn scan_genes(
    path: &Path,
    mut visit: impl FnMut(String, u64, Node) -> Result<()>,
) -> Result<()> {
    let metadata = parquet::metadata(path)?;
    let keys = super::super::parse_texts(
        metadata
            .get("priva.keys")
            .ok_or_else(|| invalid("history inventory absent"))?,
    )?;
    let mut current: Option<(String, u64, Node)> = None;
    let mut observed = 0;
    parquet::scan(path, |row| {
        if row.kind == "history-gene" {
            unused(&row)?;
            if row.gene_order != observed as u64
                || keys.get(observed) != Some(&row.gene)
                || !row.path.is_empty()
                || row.member_order != 0
                || row.key != row.gene
            {
                return Err(invalid("history gene inventory/order differs"));
            }
            if let Some((gene, ordinal, node)) = current.take() {
                visit(gene, ordinal, node)?;
            }
            observed += 1;
            current = Some((row.gene, row.gene_order, Node::default()));
        } else {
            let (gene, ordinal, node) = current
                .as_mut()
                .ok_or_else(|| invalid("history row precedes gene"))?;
            if row.gene != *gene || row.gene_order != *ordinal {
                return Err(invalid("history row gene differs"));
            }
            insert(node, row)?;
        }
        Ok(())
    })?;
    if let Some((gene, ordinal, node)) = current {
        visit(gene, ordinal, node)?;
    }
    if observed != keys.len() {
        return Err(invalid("history gene coverage differs"));
    }
    Ok(())
}
