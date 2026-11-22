//! Exon construction history only. Score/residue payloads remain in the normal
//! Alpha model; marker entries retain their insertion positions and collisions.
use super::super::model::{Minimum, Observation};
use super::super::{Result, unsupported};
use crate::python_set_order::StringSet;
use std::collections::{HashMap, HashSet};

pub(super) struct Transcript {
    pub exons: Vec<String>,
    seen: HashSet<String>,
    pub minimum: Minimum,
    pub maximum: i64,
    pub table: Option<StringSet>,
}
impl Transcript {
    pub fn from_sequence(exons: Vec<String>, minimum: Minimum, maximum: i64) -> Result<Self> {
        let seen: HashSet<_> = exons.iter().cloned().collect();
        if seen.len() != exons.len() {
            return Err(unsupported("duplicate exon history member"));
        }
        Ok(Self {
            exons,
            seen,
            minimum,
            maximum,
            table: None,
        })
    }
    fn new() -> Self {
        Self {
            exons: vec![],
            seen: HashSet::new(),
            minimum: Minimum::Unknown,
            maximum: 0,
            table: None,
        }
    }
    fn observe(&mut self, exon: &str, position: Option<i64>) {
        if self.seen.insert(exon.to_owned()) {
            self.exons.push(exon.to_owned());
        }
        if let Some(value) = position {
            if matches!(self.minimum, Minimum::Unknown)
                || matches!(self.minimum, Minimum::Position(old) if value < old)
            {
                self.minimum = Minimum::Position(value);
            }
            if value > self.maximum {
                self.maximum = value;
            }
        }
    }
    pub fn prepare(&mut self, seed: u32, worker_transport: bool) {
        let mut table = StringSet::new(seed);
        for exon in self.exons.drain(..) {
            table.insert(exon);
        }
        self.seen.clear();
        if worker_transport {
            // Pool.map returns each worker value through pickle's iterable
            // ADDITEMS reconstruction in the parent process, once.
            let mut parent = StringSet::new(seed);
            for exon in table.into_order() {
                parent.insert(exon);
            }
            table = parent;
        }
        self.table = Some(table);
    }
    fn merge(&mut self, other: Self) -> Result<()> {
        let left = self
            .table
            .as_mut()
            .ok_or_else(|| unsupported("unprepared target exon history"))?;
        let right = other
            .table
            .as_ref()
            .ok_or_else(|| unsupported("unprepared source exon history"))?;
        left.update_set(right).map_err(unsupported)?;
        if let Minimum::Position(value) = other.minimum {
            if matches!(self.minimum, Minimum::Unknown)
                || matches!(self.minimum, Minimum::Position(old) if value < old)
            {
                self.minimum = Minimum::Position(value);
            }
        }
        if other.maximum > self.maximum {
            self.maximum = other.maximum;
        }
        Ok(())
    }
    pub fn sequence(&self) -> Vec<String> {
        self.table.as_ref().map_or_else(
            || self.exons.clone(),
            |table| table.iter().map(str::to_owned).collect(),
        )
    }
}

pub(super) enum Entry {
    Branch(Node),
    Transcript(Transcript),
    Marker,
}
#[derive(Default)]
pub(super) struct Node {
    pub members: Vec<(String, Entry)>,
    index: HashMap<String, usize>,
}
impl Node {
    pub fn append(&mut self, key: String, value: Entry) -> Result<()> {
        if self.index.insert(key.clone(), self.members.len()).is_some() {
            return Err(unsupported("duplicate exon-history path member"));
        }
        self.members.push((key, value));
        Ok(())
    }
    fn member(&mut self, key: &str, create: impl FnOnce() -> Entry) -> &mut Entry {
        let at = *self.index.entry(key.to_owned()).or_insert_with(|| {
            self.members.push((key.to_owned(), create()));
            self.members.len() - 1
        });
        &mut self.members[at].1
    }
    pub fn child(&mut self, key: &str) -> Result<&mut Self> {
        match self.member(key, || Entry::Branch(Node::default())) {
            Entry::Branch(node) => Ok(node),
            _ => Err(unsupported(
                "domain hierarchy collides with a transcript or distribution",
            )),
        }
    }
    pub fn existing_child(&mut self, key: &str) -> Result<&mut Self> {
        let at = *self
            .index
            .get(key)
            .ok_or_else(|| unsupported("missing exon-history parent"))?;
        match &mut self.members[at].1 {
            Entry::Branch(node) => Ok(node),
            _ => Err(unsupported("exon-history parent is not a dictionary")),
        }
    }
    fn observe(&mut self, transcript: &str, exon: &str, position: Option<i64>) -> Result<()> {
        match self.member(transcript, || Entry::Transcript(Transcript::new())) {
            Entry::Transcript(item) => item.observe(exon, position),
            _ => {
                return Err(unsupported(
                    "transcript name collides with domain hierarchy",
                ));
            }
        }
        // These five assignments follow the same transcript mutation in the
        // original collector. Their numeric contents never affect exon sets.
        for name in [
            "distribution",
            "gnomAD_distribution",
            "gnomAD_residue_distribution",
            "min_distribution",
            "max_distribution",
        ] {
            if !matches!(self.member(name, || Entry::Marker), Entry::Marker) {
                return Err(unsupported(
                    "distribution name collides with domain hierarchy",
                ));
            }
        }
        Ok(())
    }
    pub fn prepare(&mut self, seed: u32, worker_transport: bool) {
        for (_, member) in &mut self.members {
            match member {
                Entry::Branch(node) => node.prepare(seed, worker_transport),
                Entry::Transcript(transcript) => transcript.prepare(seed, worker_transport),
                Entry::Marker => {}
            }
        }
    }
    fn merge(&mut self, other: Self, depth: usize) -> Result<()> {
        for (key, value) in other.members {
            let Some(&at) = self.index.get(&key) else {
                self.append(key, value)?;
                continue;
            };
            match (&mut self.members[at].1, value) {
                (Entry::Branch(left), Entry::Branch(right))
                    if depth < 2 || !key.starts_with("ENST") =>
                {
                    left.merge(right, depth + 1)?
                }
                (Entry::Transcript(left), Entry::Transcript(right))
                    if depth >= 2 && key.starts_with("ENST") =>
                {
                    left.merge(right)?
                }
                (Entry::Marker, Entry::Marker) if depth >= 2 && !key.starts_with("ENST") => {}
                _ => {
                    return Err(unsupported(
                        "incompatible cross-chromosome domain member layout",
                    ));
                }
            }
        }
        Ok(())
    }
    pub fn merge_gene(&mut self, other: Self) -> Result<()> {
        self.merge(other, 1)
    }
}

#[derive(Default)]
pub(super) struct Model {
    pub genes: Node,
}
impl Model {
    pub fn observe(&mut self, item: Observation<'_>) -> Result<()> {
        for entry in item.domains.split('&') {
            let parts: Vec<_> = entry.split(':').collect();
            if ["PDB-ENSP_mappings", "AFDB-ENSP_mappings", "MobiDB_lite"].contains(&parts[0]) {
                continue;
            }
            // Every prefix visits every ancestor, in the original source order.
            for length in 1..parts.len() {
                let mut node = self.genes.child(item.gene)?.child(parts[0])?;
                for level in &parts[1..=length] {
                    node = node.child(level)?;
                    node.observe(item.transcript, item.exon, item.position)?;
                }
            }
        }
        Ok(())
    }
    pub fn prepare(&mut self, seed: u32, worker_transport: bool) {
        self.genes.prepare(seed, worker_transport);
    }
}
