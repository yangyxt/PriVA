//! Ordered mutations from `stat_protein_domain_amscores.py`.
//! Repeated residue observations are losslessly run-length encoded by f64 bits;
//! they are never deduplicated and signed zeros/NaN payloads remain distinct.
use super::key::Key;
use super::{Result, exception, unsupported};
use std::collections::{BTreeSet, HashMap};

#[derive(Clone, Debug)]
pub(super) enum Minimum {
    Unknown,
    Position(i64),
}
#[derive(Clone, Debug)]
pub(super) struct Transcript {
    pub exons: BTreeSet<String>,
    pub minimum: Minimum,
    pub maximum: i64,
}
impl Transcript {
    fn new() -> Self {
        Self {
            exons: BTreeSet::new(),
            minimum: Minimum::Unknown,
            maximum: 0,
        }
    }
    fn observe(&mut self, exon: &str, position: Option<i64>) {
        self.exons.insert(exon.to_owned());
        if let Some(position) = position {
            if matches!(self.minimum, Minimum::Unknown)
                || matches!(self.minimum, Minimum::Position(previous) if position < previous)
            {
                self.minimum = Minimum::Position(position);
            }
            if position > self.maximum {
                self.maximum = position;
            }
        }
    }
    fn merge(&mut self, other: Self) {
        self.exons.extend(other.exons);
        if let Minimum::Position(value) = other.minimum {
            if matches!(self.minimum, Minimum::Unknown)
                || matches!(self.minimum, Minimum::Position(previous) if value < previous)
            {
                self.minimum = Minimum::Position(value);
            }
        }
        if other.maximum > self.maximum {
            self.maximum = other.maximum;
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Distribution {
    Score,
    Af,
    Minimum,
    Maximum,
}
impl Distribution {
    pub fn name(self) -> &'static str {
        match self {
            Self::Score => "distribution",
            Self::Af => "gnomAD_distribution",
            Self::Minimum => "min_distribution",
            Self::Maximum => "max_distribution",
        }
    }
    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "distribution" => Some(Self::Score),
            "gnomAD_distribution" => Some(Self::Af),
            "min_distribution" => Some(Self::Minimum),
            "max_distribution" => Some(Self::Maximum),
            _ => None,
        }
    }
    fn replaces(self, old: f64, new: f64) -> bool {
        match self {
            Self::Score => false,
            Self::Af | Self::Maximum => new > old,
            Self::Minimum => new < old,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub(super) struct Scores {
    pub items: Vec<(Key, f64)>,
    index: HashMap<Key, usize>,
}
impl Scores {
    pub fn insert(&mut self, key: Key, value: f64, mode: Distribution) {
        if let Some(&index) = self.index.get(&key) {
            if mode.replaces(self.items[index].1, value) {
                self.items[index].1 = value;
            }
        } else {
            self.index.insert(key.clone(), self.items.len());
            self.items.push((key, value));
        }
    }
    pub fn append_unique(&mut self, key: Key, value: f64) -> Result<()> {
        if self.index.contains_key(&key) {
            return Err(exception("ValueError", "duplicate raw distribution key"));
        }
        self.index.insert(key.clone(), self.items.len());
        self.items.push((key, value));
        Ok(())
    }
}

#[derive(Clone, Debug)]
pub(super) enum Samples {
    One(f64, u64),
    Many(Vec<(f64, u64)>),
}
impl Samples {
    pub fn new(value: f64) -> Self {
        Self::One(value, 1)
    }
    pub fn visit(&self, mut visit: impl FnMut(f64, u64)) {
        match self {
            Self::One(value, count) => visit(*value, *count),
            Self::Many(runs) => {
                for &(value, count) in runs {
                    visit(value, count);
                }
            }
        }
    }
    pub fn push_run(&mut self, value: f64, count: u64) -> Result<()> {
        if count == 0 {
            return Err(exception(
                "ValueError",
                "a residue repeat count must be positive",
            ));
        }
        match self {
            Self::One(old, old_count) if old.to_bits() == value.to_bits() => {
                *old_count = old_count
                    .checked_add(count)
                    .ok_or_else(|| unsupported("residue sample count exceeds u64"))?;
            }
            Self::One(old, old_count) => {
                *self = Self::Many(vec![(*old, *old_count), (value, count)]);
            }
            Self::Many(runs) => {
                let last = runs.last_mut().unwrap();
                if last.0.to_bits() == value.to_bits() {
                    last.1 = last
                        .1
                        .checked_add(count)
                        .ok_or_else(|| unsupported("residue sample count exceeds u64"))?;
                } else {
                    runs.push((value, count));
                }
            }
        }
        Ok(())
    }
    fn extend(&mut self, other: Self) -> Result<()> {
        match other {
            Self::One(value, count) => self.push_run(value, count),
            Self::Many(runs) => {
                for (value, count) in runs {
                    self.push_run(value, count)?;
                }
                Ok(())
            }
        }
    }
}

#[derive(Clone, Debug, Default)]
pub(super) struct Residues {
    pub items: Vec<(Key, Samples)>,
    index: HashMap<Key, usize>,
}
impl Residues {
    pub fn insert(&mut self, key: Key, samples: Samples) -> Result<()> {
        if let Some(&index) = self.index.get(&key) {
            self.items[index].1.extend(samples)?;
        } else {
            self.index.insert(key.clone(), self.items.len());
            self.items.push((key, samples));
        }
        Ok(())
    }
    pub fn append_unique(&mut self, key: Key, samples: Samples) -> Result<()> {
        if self.index.contains_key(&key) {
            return Err(exception("ValueError", "duplicate residue key"));
        }
        self.index.insert(key.clone(), self.items.len());
        self.items.push((key, samples));
        Ok(())
    }
}

#[derive(Clone, Debug)]
pub(super) enum Member {
    Node(Node),
    Transcript(Transcript),
    Scores(Distribution, Scores),
    Residues(Residues),
    Array(Distribution, Vec<f64>),
}
#[derive(Clone, Debug, Default)]
pub(super) struct Node {
    pub members: Vec<(String, Member)>,
    index: HashMap<String, usize>,
}
impl Node {
    pub fn member(&mut self, key: &str, new: impl FnOnce() -> Member) -> &mut Member {
        let index = *self.index.entry(key.to_owned()).or_insert_with(|| {
            let index = self.members.len();
            self.members.push((key.to_owned(), new()));
            index
        });
        &mut self.members[index].1
    }
    pub fn append_unique(&mut self, key: String, value: Member) -> Result<()> {
        if self.index.contains_key(&key) {
            return Err(exception("ValueError", "duplicate ordered model member"));
        }
        self.index.insert(key.clone(), self.members.len());
        self.members.push((key, value));
        Ok(())
    }
    pub fn child(&mut self, key: &str) -> Result<&mut Node> {
        match self.member(key, || Member::Node(Node::default())) {
            Member::Node(node) => Ok(node),
            _ => Err(unsupported(
                "domain hierarchy collides with a transcript or distribution",
            )),
        }
    }
    pub fn existing_child(&mut self, key: &str) -> Result<&mut Node> {
        let index = self
            .index
            .get(key)
            .copied()
            .ok_or_else(|| exception("ValueError", "missing parent domain row"))?;
        match &mut self.members[index].1 {
            Member::Node(node) => Ok(node),
            _ => Err(exception("ValueError", "parent is not a domain")),
        }
    }
    fn observe(
        &mut self,
        transcript: &str,
        exon: &str,
        position: Option<i64>,
        variant: &Key,
        residue: &Key,
        score: f64,
        af: f64,
    ) -> Result<()> {
        match self.member(transcript, || Member::Transcript(Transcript::new())) {
            Member::Transcript(item) => item.observe(exon, position),
            _ => {
                return Err(unsupported(
                    "transcript name collides with domain hierarchy",
                ));
            }
        }
        for (mode, key, value) in [
            (Distribution::Score, variant, score),
            (Distribution::Af, variant, af),
        ] {
            match self.member(mode.name(), || Member::Scores(mode, Scores::default())) {
                Member::Scores(_, items) => items.insert(key.clone(), value, mode),
                _ => {
                    return Err(unsupported(
                        "distribution name collides with domain hierarchy",
                    ));
                }
            }
        }
        match self.member("gnomAD_residue_distribution", || {
            Member::Residues(Residues::default())
        }) {
            Member::Residues(items) => items.insert(residue.clone(), Samples::new(af))?,
            _ => {
                return Err(unsupported(
                    "residue distribution name collides with domain hierarchy",
                ));
            }
        }
        for mode in [Distribution::Minimum, Distribution::Maximum] {
            match self.member(mode.name(), || Member::Scores(mode, Scores::default())) {
                Member::Scores(_, items) => items.insert(residue.clone(), score, mode),
                _ => {
                    return Err(unsupported(
                        "distribution name collides with domain hierarchy",
                    ));
                }
            }
        }
        Ok(())
    }
    pub fn merge(&mut self, other: Self) -> Result<()> {
        for (key, value) in other.members {
            if !self.index.contains_key(&key) {
                self.append_unique(key, value)?;
                continue;
            }
            let old = &mut self.members[self.index[&key]].1;
            // Literal `_merge_domain_level` dispatch is key-first: ENST takes
            // precedence even if a malformed hierarchy uses that prefix.
            match (old, value) {
                (Member::Node(left), Member::Node(right)) if !key.starts_with("ENST") => {
                    left.merge(right)?
                }
                (Member::Transcript(left), Member::Transcript(right))
                    if key.starts_with("ENST") =>
                {
                    left.merge(right)
                }
                (Member::Scores(mode, left), Member::Scores(other_mode, right))
                    if *mode == other_mode =>
                {
                    for (k, v) in right.items {
                        left.insert(k, v, *mode);
                    }
                }
                (Member::Residues(left), Member::Residues(right)) => {
                    for (k, v) in right.items {
                        left.insert(k, v)?;
                    }
                }
                _ => {
                    return Err(unsupported(
                        "incompatible cross-chromosome domain member layout",
                    ));
                }
            }
        }
        Ok(())
    }
    fn merge_branches(&mut self, other: Self, depth: usize) -> Result<()> {
        for (key, value) in other.members {
            if !self.index.contains_key(&key) {
                self.append_unique(key, value)?;
                continue;
            }
            match (&mut self.members[self.index[&key]].1, value) {
                (Member::Node(left), Member::Node(right)) if depth > 1 => {
                    left.merge_branches(right, depth - 1)?
                }
                (Member::Node(left), Member::Node(right)) => left.merge(right)?,
                _ => return Err(unsupported("gene/source level is not a dictionary")),
            }
        }
        Ok(())
    }
    pub fn merge_gene(&mut self, other: Self) -> Result<()> {
        self.merge_branches(other, 1)
    }
}

#[derive(Default)]
pub(super) struct Model {
    pub genes: Node,
}
pub(super) struct Observation<'a> {
    pub gene: &'a str,
    pub domains: &'a str,
    pub transcript: &'a str,
    pub exon: &'a str,
    pub position: Option<i64>,
    pub variant: &'a Key,
    pub residue: &'a Key,
    pub score: f64,
    pub af: f64,
}
impl Model {
    pub fn observe(&mut self, observation: Observation<'_>) -> Result<()> {
        let Observation {
            gene,
            domains,
            transcript,
            exon,
            position,
            variant,
            residue,
            score,
            af,
        } = observation;
        for entry in domains.split('&') {
            let parts: Vec<_> = entry.split(':').collect();
            if ["PDB-ENSP_mappings", "AFDB-ENSP_mappings", "MobiDB_lite"].contains(&parts[0]) {
                continue;
            }
            // The original emits every prefix, then visits every level of that
            // prefix. Repeated parent visits intentionally repeat AF samples.
            for length in 1..parts.len() {
                let mut node = self.genes.child(gene)?.child(parts[0])?;
                for level in &parts[1..=length] {
                    node = node.child(level)?;
                    node.observe(transcript, exon, position, variant, residue, score, af)?;
                }
            }
        }
        Ok(())
    }
    pub fn merge(&mut self, other: Self) -> Result<()> {
        self.genes.merge_branches(other.genes, 2)
    }
}
