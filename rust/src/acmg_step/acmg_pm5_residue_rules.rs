//! PM5 by the owner's rules of 30 October 2026, which apply to the Rust step only: Python's
//! PM5 (`PS1_PM5_criteria`) is unchanged, so the Rust PM5 differs from it on purpose. The design
//! is `PriVA_upgrade/PriVA_PM5_Rust_Design_20261030.md`:
//!
//! | Confident pathogenic ClinVar variant | Query variants that get PM5 |
//! | --- | --- |
//! | 1. missense, not gain of function | missense at the same residue (another amino acid); in-frame deletion covering the residue; frameshift or stop-gain at or upstream of the residue |
//! | 1b. missense, gain of function (GoFCards) | missense at the same residue; in-frame deletion covering it |
//! | 2. in-frame deletion | in-frame deletion spanning all of its residues; frameshift or stop-gain upstream of it |
//! | 3. frameshift escaping NMD | frameshift at or upstream of it |
//! | 3b. stop-gain escaping NMD | frameshift or stop-gain at or upstream of it |
//! | 4. in-frame insertion or duplication | the same type covering all of its residues |
//! | 5. truncating variant that triggers NMD, and any other change | not counted |
//!
//! "Confident" is PS1's rule (CLNSIG holding `Pathogenic` with two stars or more, or exactly
//! `Likely_pathogenic` with three or more). An identical protein change (the same HGVSp) is PS1's
//! business, not PM5's. A truncating query (frameshift, stop-gain) gets PM5 only when PVS1 is not
//! applied to it, and PS1 takes precedence (`ps1_pm5_criteria`).
//!
//! The Python amino-acid cache holds neither consequence, NMD status, residue span nor gain of
//! function, so the facts come from an index of the VEP-annotated ClinVar VCF, written once by
//! `build_pm5_index` (the records and transcript entries the Python cache builder,
//! `stat_aachange_clinvar.py`, takes: a CSQ entry with a transcript, a protein position and an
//! HGVSp) and read for the table's transcripts by `Pm5Index::read`. Gain of function is the
//! curated GOF of the non-LoF mechanism cache (`ExactVariantMatcher`, GoFCards alleles).
//!
//! A variant's protein change is read from VEP's `Consequence` (`&`-joined terms, taken in the
//! order frameshift, stop_gained, inframe_deletion, inframe_insertion, missense_variant), its
//! `HGVSp` (a `delins` is "another change"; an in-frame insertion is a duplication when its HGVSp
//! says `dup`, else an insertion) and `Protein_position` (`first[-last][/length]`; a `?` makes the
//! change unplaceable). The residues of insertions and duplications are VEP's `Protein_position`
//! for both the ClinVar variant and the query.
use super::gene_mechanism_variants::ExactVariantMatcher;
use super::types::SymbolResolver;
use crate::cli::{Result, invalid};
use crate::pandas_table::Cell;
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

/// The kinds of protein change the rules tell apart.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ChangeKind {
    Missense,
    InframeDeletion,
    Frameshift,
    StopGain,
    Insertion,
    Duplication,
}

impl ChangeKind {
    pub(crate) fn truncating(self) -> bool {
        matches!(self, ChangeKind::Frameshift | ChangeKind::StopGain)
    }
}

/// A protein change placed on its residues (first to last, 1-based).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ProteinChange {
    pub(crate) kind: ChangeKind,
    pub(crate) first: i64,
    pub(crate) last: i64,
}

/// The residues of a `Protein_position` text: `150`, `150-152`, either with `/length`.
fn residues(protein_position: &str) -> Option<(i64, i64)> {
    let span = protein_position.split('/').next()?;
    let number = |text: &str| -> Option<i64> {
        (!text.is_empty() && text.len() < 10 && text.bytes().all(|byte| byte.is_ascii_digit()))
            .then(|| text.parse().ok())
            .flatten()
    };
    match span.split_once('-') {
        Some((first, last)) => {
            let (first, last) = (number(first)?, number(last)?);
            (first <= last).then_some((first, last))
        }
        None => number(span).map(|first| (first, first)),
    }
}

/// The protein change of a variant from its VEP `Consequence`, `Protein_position` and `HGVSp`;
/// `None` for a change the rules do not count or cannot place.
pub(crate) fn protein_change(
    consequence: &str,
    protein_position: &str,
    hgvsp: &str,
) -> Option<ProteinChange> {
    let terms: Vec<&str> = consequence.split('&').collect();
    let has = |term: &str| terms.contains(&term);
    let protein = hgvsp.rsplit(':').next().unwrap_or("");
    let kind = if has("frameshift_variant") {
        ChangeKind::Frameshift
    } else if has("stop_gained") {
        ChangeKind::StopGain
    } else if protein.contains("delins") {
        return None;
    } else if has("inframe_deletion") {
        ChangeKind::InframeDeletion
    } else if has("inframe_insertion") {
        if protein.contains("dup") {
            ChangeKind::Duplication
        } else if protein.contains("ins") {
            ChangeKind::Insertion
        } else {
            return None;
        }
    } else if has("missense_variant") {
        ChangeKind::Missense
    } else {
        return None;
    };
    let (first, last) = residues(protein_position)?;
    if kind == ChangeKind::Missense && first != last {
        return None;
    }
    Some(ProteinChange { kind, first, last })
}

/// One confident pathogenic ClinVar entry of a transcript.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Pm5Entry {
    pub(crate) change: ProteinChange,
    pub(crate) hgvsp: String,
    /// VEP's NMD plugin marks the variant `NMD_escaping_variant`.
    pub(crate) escapes_nmd: bool,
    /// For the gain-of-function look-up.
    pub(crate) symbol: String,
    pub(crate) hgvsc: String,
    pub(crate) genomic: [String; 4],
}

/// The columns of the index file, in order.
const INDEX_COLUMNS: [&str; 13] = [
    "transcript",
    "protein_position",
    "hgvsp",
    "consequence",
    "nmd",
    "symbol",
    "hgvsc",
    "chrom",
    "pos",
    "ref",
    "alt",
    "clnsig",
    "clnrevstat",
];

/// The review statuses and their stars (`high_confidence_status` of the step).
fn review_stars(status: &str) -> i64 {
    match status {
        "practice_guideline" => 4,
        "reviewed_by_expert_panel" => 3,
        "criteria_provided,_multiple_submitters,_no_conflicts" => 2,
        _ => 0,
    }
}

/// PS1's confidence rule for one ClinVar record (CLNSIG and CLNREVSTAT as the VCF writes them).
pub(crate) fn confident_pathogenic(clnsig: &str, clnrevstat: &str) -> bool {
    let stars = review_stars(clnrevstat);
    (clnsig.contains("Pathogenic") && stars >= 2) || (clnsig == "Likely_pathogenic" && stars >= 3)
}

fn reader(path: &str) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).map_err(|error| invalid(format!("{path}: {error}")))?;
    let source: Box<dyn Read> = if path.ends_with(".gz") {
        Box::new(MultiGzDecoder::new(BufReader::with_capacity(1 << 20, file)))
    } else {
        Box::new(file)
    };
    Ok(Box::new(BufReader::with_capacity(1 << 20, source)))
}

/// Writes the index of the confident pathogenic protein changes of the VEP-annotated ClinVar
/// VCF at `vcf` to `output` (gzip when the name ends in `.gz`): one row per CSQ entry with a
/// transcript, a protein position and an HGVSp, as the Python cache builder takes them.
/// Returns the number of rows.
pub(crate) fn build_pm5_index(vcf: &str, output: &str) -> Result<usize> {
    let mut input = reader(vcf)?;
    let file = File::create(output).map_err(|error| invalid(format!("{output}: {error}")))?;
    let mut out: Box<dyn Write> = if output.ends_with(".gz") {
        Box::new(GzEncoder::new(BufWriter::new(file), Compression::default()))
    } else {
        Box::new(BufWriter::new(file))
    };
    writeln!(out, "{}", INDEX_COLUMNS.join("\t"))?;
    let mut fields: Option<HashMap<String, usize>> = None;
    let mut line = Vec::new();
    let mut rows = 0_usize;
    loop {
        line.clear();
        if input.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        let text = std::str::from_utf8(&line)
            .map_err(|_| invalid(format!("{vcf}: a line is not UTF-8 text")))?
            .trim_end_matches(['\n', '\r']);
        if let Some(header) = text.strip_prefix("##INFO=<ID=CSQ,") {
            let format = header
                .split("Format: ")
                .nth(1)
                .ok_or_else(|| invalid(format!("{vcf}: the CSQ header has no Format")))?
                .trim_end_matches(">")
                .trim_end_matches('"');
            fields = Some(
                format
                    .split('|')
                    .enumerate()
                    .map(|(index, name)| (name.to_owned(), index))
                    .collect(),
            );
            continue;
        }
        if text.starts_with('#') {
            continue;
        }
        let Some(fields) = &fields else {
            return Err(invalid(format!(
                "{vcf}: records come before the CSQ header"
            )));
        };
        let columns: Vec<&str> = text.splitn(9, '\t').collect();
        if columns.len() < 8 || !columns[7].contains("athogenic") {
            continue;
        }
        let mut clnsig = None;
        let mut clnrevstat = None;
        let mut csq = None;
        for item in columns[7].split(';') {
            match item.split_once('=') {
                Some(("CLNSIG", value)) => clnsig = Some(value),
                Some(("CLNREVSTAT", value)) => clnrevstat = Some(value),
                Some(("CSQ", value)) => csq = Some(value),
                _ => {}
            }
        }
        let (Some(clnsig), Some(clnrevstat), Some(csq)) = (clnsig, clnrevstat, csq) else {
            continue;
        };
        if !confident_pathogenic(clnsig, clnrevstat) {
            continue;
        }
        let index = |name: &str| -> Result<usize> {
            fields
                .get(name)
                .copied()
                .ok_or_else(|| invalid(format!("{vcf}: CSQ has no {name}")))
        };
        let wanted = [
            index("Feature")?,
            index("Protein_position")?,
            index("HGVSp")?,
            index("Consequence")?,
            index("NMD")?,
            index("SYMBOL")?,
            index("HGVSc")?,
        ];
        let widest = *wanted.iter().max().expect("seven fields");
        for entry in csq.split(',') {
            let parts: Vec<&str> = entry.split('|').collect();
            if parts.len() <= widest {
                continue;
            }
            let values = wanted.map(|position| parts[position]);
            if values[0].is_empty() || values[1].is_empty() || values[2].is_empty() {
                continue;
            }
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                values.join("\t"),
                columns[0],
                columns[1],
                columns[3],
                columns[4],
                clnsig,
                clnrevstat
            )?;
            rows += 1;
        }
    }
    out.flush()?;
    drop(out);
    Ok(rows)
}

/// The confident pathogenic ClinVar protein changes of the table's transcripts.
#[derive(Debug, Default)]
pub(crate) struct Pm5Index {
    transcripts: HashMap<String, Vec<Pm5Entry>>,
}

impl Pm5Index {
    /// Reads the index written by `build_pm5_index` for the transcripts in `wanted`, keeping the
    /// entries whose change the rules count and can place.
    pub(crate) fn read(path: &str, wanted: &HashSet<String>) -> Result<Pm5Index> {
        let mut input = reader(path)?;
        let mut line = String::new();
        input.read_line(&mut line)?;
        if line.trim_end_matches(['\n', '\r']) != INDEX_COLUMNS.join("\t") {
            return Err(invalid(format!(
                "{path} is not a PM5 index (its first line is not the index header)"
            )));
        }
        let mut index = Pm5Index::default();
        loop {
            line.clear();
            if input.read_line(&mut line)? == 0 {
                break;
            }
            let values: Vec<&str> = line.trim_end_matches(['\n', '\r']).split('\t').collect();
            if values.len() != INDEX_COLUMNS.len() {
                return Err(invalid(format!(
                    "{path}: a row has {} fields",
                    values.len()
                )));
            }
            if !wanted.contains(values[0]) {
                continue;
            }
            if !confident_pathogenic(values[11], values[12]) {
                return Err(invalid(format!(
                    "{path}: a row is not a confident pathogenic ClinVar variant"
                )));
            }
            let Some(change) = protein_change(values[3], values[1], values[2]) else {
                continue;
            };
            index
                .transcripts
                .entry(values[0].to_owned())
                .or_default()
                .push(Pm5Entry {
                    change,
                    hgvsp: values[2].to_owned(),
                    escapes_nmd: values[4]
                        .split('&')
                        .any(|flag| flag == "NMD_escaping_variant"),
                    symbol: values[5].to_owned(),
                    hgvsc: values[6].to_owned(),
                    genomic: [values[7], values[8], values[9], values[10]].map(str::to_owned),
                });
        }
        Ok(index)
    }

    pub(crate) fn entries(&self, transcript: &str) -> &[Pm5Entry] {
        self.transcripts
            .get(transcript)
            .map(Vec::as_slice)
            .unwrap_or_default()
    }

    #[cfg(test)]
    pub(crate) fn from_entries(entries: Vec<(&str, Pm5Entry)>) -> Pm5Index {
        let mut index = Pm5Index::default();
        for (transcript, entry) in entries {
            index
                .transcripts
                .entry(transcript.to_owned())
                .or_default()
                .push(entry);
        }
        index
    }
}

/// Whether a confident pathogenic missense is a curated gain-of-function allele.
pub(crate) trait GainOfFunction {
    fn is_gain_of_function(&mut self, entry: &Pm5Entry) -> bool;
}

/// Gain of function from the non-LoF mechanism cache: `match_curated_nonlof_variant` scores
/// GOF 2. An entry the matcher cannot look up (for example a symbol the HGNC resolver refuses)
/// counts as not gain of function; `unresolved` counts them.
pub(crate) struct CuratedGainOfFunction<'r> {
    matcher: ExactVariantMatcher<'r>,
    assembly: String,
    known: HashMap<String, bool>,
    pub(crate) unresolved: usize,
}

impl<'r> CuratedGainOfFunction<'r> {
    pub(crate) fn new(
        mechanism_json: &str,
        resolver: &'r dyn SymbolResolver,
        assembly: &str,
    ) -> CuratedGainOfFunction<'r> {
        CuratedGainOfFunction {
            matcher: ExactVariantMatcher::new(mechanism_json, resolver),
            assembly: assembly.to_owned(),
            known: HashMap::new(),
            unresolved: 0,
        }
    }
}

impl GainOfFunction for CuratedGainOfFunction<'_> {
    fn is_gain_of_function(&mut self, entry: &Pm5Entry) -> bool {
        let key = format!(
            "{}\t{}\t{}\t{}",
            entry.symbol,
            entry.hgvsp,
            entry.hgvsc,
            entry.genomic.join("\t")
        );
        if let Some(found) = self.known.get(&key) {
            return *found;
        }
        let [chrom, pos, ref_allele, alt_allele] = &entry.genomic;
        let found = match self.matcher.match_curated_nonlof_variant(
            Cell::Text(&entry.symbol),
            Some(Cell::Text(&entry.hgvsp)),
            Some(Cell::Text(&entry.hgvsc)),
            [
                Some(Cell::Text(chrom)),
                Some(Cell::Text(pos)),
                Some(Cell::Text(ref_allele)),
                Some(Cell::Text(alt_allele)),
            ],
            &self.assembly,
            "auto",
        ) {
            Ok(found) => found.scores[1] == 2,
            Err(_) => {
                self.unresolved += 1;
                false
            }
        };
        self.known.insert(key, found);
        found
    }
}

/// Gain of function given as a list of HGVSp (comparisons and tests).
pub(crate) struct ListedGainOfFunction(pub(crate) HashSet<String>);

impl GainOfFunction for ListedGainOfFunction {
    fn is_gain_of_function(&mut self, entry: &Pm5Entry) -> bool {
        self.0.contains(&entry.hgvsp)
    }
}

/// The rule under which a ClinVar variant gives PM5.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Pm5Rule {
    /// Rules 1 and 1b alike: a missense or an in-frame deletion query at the residue.
    MissenseAtResidue,
    /// Rule 1 for a truncating query: the missense is not gain of function.
    MissenseDownstreamOfTruncation,
    InframeDeletion,
    FrameshiftEscapingNmd,
    StopGainEscapingNmd,
    InsertionOrDuplication,
}

impl Pm5Rule {
    pub(crate) fn label(self) -> &'static str {
        match self {
            Pm5Rule::MissenseAtResidue => "1/1b missense at the residue",
            Pm5Rule::MissenseDownstreamOfTruncation => "1 missense (not GoF) at or downstream",
            Pm5Rule::InframeDeletion => "2 in-frame deletion",
            Pm5Rule::FrameshiftEscapingNmd => "3 frameshift escaping NMD",
            Pm5Rule::StopGainEscapingNmd => "3b stop-gain escaping NMD",
            Pm5Rule::InsertionOrDuplication => "4 in-frame insertion or duplication",
        }
    }
}

/// The rule under which `entry` gives PM5 to the query change, if any; gain of function is
/// looked up only when it decides (a missense with a truncating query).
pub(crate) fn pm5_rule(
    entry: &Pm5Entry,
    query: &ProteinChange,
    query_hgvsp: Option<&str>,
    gain_of_function: &mut dyn GainOfFunction,
) -> Option<Pm5Rule> {
    if query_hgvsp == Some(entry.hgvsp.as_str()) {
        return None;
    }
    let (first, last) = (entry.change.first, entry.change.last);
    let covers = |from: i64, to: i64| query.first <= from && query.last >= to;
    match (entry.change.kind, query.kind) {
        (ChangeKind::Missense, ChangeKind::Missense) => {
            (query.first == first).then_some(Pm5Rule::MissenseAtResidue)
        }
        (ChangeKind::Missense, ChangeKind::InframeDeletion) => {
            covers(first, first).then_some(Pm5Rule::MissenseAtResidue)
        }
        (ChangeKind::Missense, ChangeKind::Frameshift | ChangeKind::StopGain) => {
            (query.first <= first && !gain_of_function.is_gain_of_function(entry))
                .then_some(Pm5Rule::MissenseDownstreamOfTruncation)
        }
        (ChangeKind::InframeDeletion, ChangeKind::InframeDeletion) => {
            covers(first, last).then_some(Pm5Rule::InframeDeletion)
        }
        (ChangeKind::InframeDeletion, ChangeKind::Frameshift | ChangeKind::StopGain) => {
            (query.first < first).then_some(Pm5Rule::InframeDeletion)
        }
        (ChangeKind::Frameshift, ChangeKind::Frameshift) if entry.escapes_nmd => {
            (query.first <= first).then_some(Pm5Rule::FrameshiftEscapingNmd)
        }
        (ChangeKind::StopGain, ChangeKind::Frameshift | ChangeKind::StopGain)
            if entry.escapes_nmd =>
        {
            (query.first <= first).then_some(Pm5Rule::StopGainEscapingNmd)
        }
        (ChangeKind::Insertion, ChangeKind::Insertion)
        | (ChangeKind::Duplication, ChangeKind::Duplication) => {
            covers(first, last).then_some(Pm5Rule::InsertionOrDuplication)
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(consequence: &str, position: &str, hgvsp: &str, escapes_nmd: bool) -> Pm5Entry {
        Pm5Entry {
            change: protein_change(consequence, position, hgvsp).expect("a counted change"),
            hgvsp: hgvsp.to_owned(),
            escapes_nmd,
            symbol: "G".into(),
            hgvsc: String::new(),
            genomic: Default::default(),
        }
    }

    fn rule(
        entry: &Pm5Entry,
        consequence: &str,
        position: &str,
        hgvsp: &str,
        gof: bool,
    ) -> Option<Pm5Rule> {
        let query = protein_change(consequence, position, hgvsp).expect("a counted query");
        let listed: HashSet<String> = if gof {
            [entry.hgvsp.clone()].into()
        } else {
            HashSet::new()
        };
        pm5_rule(
            entry,
            &query,
            Some(hgvsp),
            &mut ListedGainOfFunction(listed),
        )
    }

    #[test]
    fn classifies_vep_changes() {
        let change =
            |c: &str, p: &str, h: &str| protein_change(c, p, h).map(|c| (c.kind, c.first, c.last));
        assert_eq!(
            change("missense_variant", "150/500", "p.Arg150Cys"),
            Some((ChangeKind::Missense, 150, 150))
        );
        assert_eq!(
            change(
                "missense_variant",
                "150-151/500",
                "p.Arg150_Gly151delinsCysSer"
            ),
            None
        );
        assert_eq!(
            change("inframe_deletion", "150-152/500", "p.Arg150_Gly152del"),
            Some((ChangeKind::InframeDeletion, 150, 152))
        );
        assert_eq!(
            change("inframe_insertion", "150-151/500", "p.Gly151dup"),
            Some((ChangeKind::Duplication, 150, 151))
        );
        assert_eq!(
            change("inframe_insertion", "150-151/500", "p.Arg150_Gly151insTrp"),
            Some((ChangeKind::Insertion, 150, 151))
        );
        assert_eq!(
            change("stop_gained&frameshift_variant", "150/500", "p.Arg150fs"),
            Some((ChangeKind::Frameshift, 150, 150))
        );
        assert_eq!(
            change(
                "stop_gained&splice_region_variant",
                "150/500",
                "p.Arg150Ter"
            ),
            Some((ChangeKind::StopGain, 150, 150))
        );
        assert_eq!(change("synonymous_variant", "150/500", "p.Arg150="), None);
        assert_eq!(change("frameshift_variant", "?-150/500", "p.X"), None);
        assert_eq!(change("missense_variant", "", "p.Arg150Cys"), None);
    }

    #[test]
    fn every_rule_of_the_design_table() {
        let missense = entry("missense_variant", "150/500", "p.Arg150Cys", false);
        // 1: missense at the residue, deletion covering it, truncation at or upstream.
        assert_eq!(
            rule(
                &missense,
                "missense_variant",
                "150/500",
                "p.Arg150His",
                false
            ),
            Some(Pm5Rule::MissenseAtResidue)
        );
        assert_eq!(
            rule(
                &missense,
                "missense_variant",
                "150/500",
                "p.Arg150Cys",
                false
            ),
            None
        );
        assert_eq!(
            rule(
                &missense,
                "missense_variant",
                "151/500",
                "p.Gly151Cys",
                false
            ),
            None
        );
        assert_eq!(
            rule(
                &missense,
                "inframe_deletion",
                "149-150/500",
                "p.Lys149_Arg150del",
                false
            ),
            Some(Pm5Rule::MissenseAtResidue)
        );
        assert_eq!(
            rule(
                &missense,
                "inframe_deletion",
                "151-152/500",
                "p.Gly151_Ala152del",
                false
            ),
            None
        );
        assert_eq!(
            rule(
                &missense,
                "frameshift_variant",
                "150/500",
                "p.Arg150fs",
                false
            ),
            Some(Pm5Rule::MissenseDownstreamOfTruncation)
        );
        assert_eq!(
            rule(&missense, "stop_gained", "20/500", "p.Trp20Ter", false),
            Some(Pm5Rule::MissenseDownstreamOfTruncation)
        );
        assert_eq!(
            rule(&missense, "stop_gained", "151/500", "p.Trp151Ter", false),
            None
        );
        // 1b: a gain-of-function missense counts for missense and covering deletions only.
        assert_eq!(
            rule(
                &missense,
                "missense_variant",
                "150/500",
                "p.Arg150His",
                true
            ),
            Some(Pm5Rule::MissenseAtResidue)
        );
        assert_eq!(
            rule(
                &missense,
                "frameshift_variant",
                "100/500",
                "p.Leu100fs",
                true
            ),
            None
        );
        assert_eq!(
            rule(&missense, "stop_gained", "100/500", "p.Leu100Ter", true),
            None
        );
        // 2: deletions spanning all residues; truncation strictly upstream.
        let deletion = entry(
            "inframe_deletion",
            "150-152/500",
            "p.Arg150_Gly152del",
            false,
        );
        assert_eq!(
            rule(
                &deletion,
                "inframe_deletion",
                "149-153/500",
                "p.Lys149_Ala153del",
                false
            ),
            Some(Pm5Rule::InframeDeletion)
        );
        assert_eq!(
            rule(
                &deletion,
                "inframe_deletion",
                "151-153/500",
                "p.Gly151_Ala153del",
                false
            ),
            None
        );
        assert_eq!(
            rule(
                &deletion,
                "frameshift_variant",
                "149/500",
                "p.Lys149fs",
                false
            ),
            Some(Pm5Rule::InframeDeletion)
        );
        assert_eq!(
            rule(&deletion, "stop_gained", "150/500", "p.Arg150Ter", false),
            None
        );
        assert_eq!(
            rule(
                &deletion,
                "missense_variant",
                "151/500",
                "p.Gly151Cys",
                false
            ),
            None
        );
        // 3 and 3b: only NMD-escaping truncations count.
        let frameshift = entry("frameshift_variant", "480/500", "p.Leu480fs", true);
        assert_eq!(
            rule(
                &frameshift,
                "frameshift_variant",
                "470/500",
                "p.Ser470fs",
                false
            ),
            Some(Pm5Rule::FrameshiftEscapingNmd)
        );
        assert_eq!(
            rule(&frameshift, "stop_gained", "470/500", "p.Ser470Ter", false),
            None
        );
        assert_eq!(
            rule(
                &frameshift,
                "frameshift_variant",
                "490/500",
                "p.Ser490fs",
                false
            ),
            None
        );
        let stop = entry("stop_gained", "480/500", "p.Trp480Ter", true);
        assert_eq!(
            rule(&stop, "stop_gained", "480/500", "p.Trp480Ter2", false),
            Some(Pm5Rule::StopGainEscapingNmd)
        );
        assert_eq!(
            rule(&stop, "frameshift_variant", "300/500", "p.Ser300fs", false),
            Some(Pm5Rule::StopGainEscapingNmd)
        );
        // 5: an NMD-triggering truncation gives nothing.
        let triggering = entry("stop_gained", "480/500", "p.Trp480Ter", false);
        assert_eq!(
            rule(&triggering, "stop_gained", "300/500", "p.Ser300Ter", false),
            None
        );
        // 4: the same type, covering all residues.
        let duplication = entry("inframe_insertion", "150-151/500", "p.Gly151dup", false);
        assert_eq!(
            rule(
                &duplication,
                "inframe_insertion",
                "150-151/500",
                "p.Arg150_Gly151dup",
                false
            ),
            Some(Pm5Rule::InsertionOrDuplication)
        );
        assert_eq!(
            rule(
                &duplication,
                "inframe_insertion",
                "150-151/500",
                "p.Arg150_Gly151insTrp",
                false
            ),
            None
        );
        assert_eq!(
            rule(
                &duplication,
                "inframe_insertion",
                "151-152/500",
                "p.Gly152dup",
                false
            ),
            None
        );
    }

    #[test]
    fn confidence_is_ps1s_rule() {
        let two = "criteria_provided,_multiple_submitters,_no_conflicts";
        assert!(confident_pathogenic("Pathogenic/Likely_pathogenic", two));
        assert!(!confident_pathogenic("Likely_pathogenic", two));
        assert!(confident_pathogenic(
            "Likely_pathogenic",
            "reviewed_by_expert_panel"
        ));
        assert!(!confident_pathogenic(
            "Likely_pathogenic|risk_factor",
            "reviewed_by_expert_panel"
        ));
        assert!(!confident_pathogenic(
            "Pathogenic",
            "criteria_provided,_single_submitter"
        ));
    }
}
