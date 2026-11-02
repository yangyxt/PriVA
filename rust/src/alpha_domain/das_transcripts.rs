//! Retained-input DAS transcript mapping. The schema-2 model has canonical
//! exon membership, so a separately bound serialized-set history is required.
//! This does not establish fresh native Alpha-to-DAS installation parity.
use super::model::{Member, Minimum, Node};
use super::{EXTRACTOR, PRODUCER, Result, exception, parquet, parse_texts, store};
use crate::cache_convert::Provenance;
use crate::cache_convert::transcript_domains::{PlainTranscript, write_entries};
use crate::python_set_order::StringSet;
use crate::rust_cache::ensure_distinct_paths;
use crate::sha256::sha256_file;
use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

fn invalid(message: &str) -> Box<dyn std::error::Error> {
    exception("ValueError", message)
}

struct History {
    input: BufReader<File>,
    seed: u32,
    sets: usize,
}
impl History {
    fn next(&mut self, path: &[String], expected: &BTreeSet<String>) -> Result<Vec<String>> {
        let mut line = String::new();
        if self.input.read_line(&mut line)? == 0 {
            return Err(invalid("exon history ended before the typed model"));
        }
        let words: Vec<_> = line.split_whitespace().collect();
        let number = |at: usize| -> Result<usize> {
            words
                .get(at)
                .and_then(|s| s.parse().ok())
                .ok_or_else(|| invalid("invalid exon-history length"))
        };
        let decode = |word: &str| -> Result<String> {
            if word == "-" {
                return Ok(String::new());
            }
            if !word.is_ascii() || !word.len().is_multiple_of(2) {
                return Err(invalid("invalid exon-history UTF-8 hex"));
            }
            let bytes = (0..word.len())
                .step_by(2)
                .map(|at| u8::from_str_radix(&word[at..at + 2], 16))
                .collect::<std::result::Result<Vec<_>, _>>()
                .map_err(|_| invalid("invalid exon-history hex"))?;
            String::from_utf8(bytes).map_err(|_| invalid("invalid exon-history UTF-8"))
        };
        let fields = number(0)?;
        if fields != path.len() || fields + 2 > words.len() {
            return Err(invalid("exon-history path length differs"));
        }
        for (wanted, word) in path.iter().zip(&words[1..fields + 1]) {
            if *wanted != decode(word)? {
                return Err(invalid("exon-history path or source order differs"));
            }
        }
        let count = number(fields + 1)?;
        if count != words.len() - fields - 2 {
            return Err(invalid("exon-history member count differs"));
        }
        let mut membership = BTreeSet::new();
        let mut loaded = StringSet::new(self.seed);
        for word in &words[fields + 2..] {
            let text = decode(word)?;
            if !membership.insert(text.clone()) {
                return Err(invalid("duplicate serialized exon set member"));
            }
            loaded.insert(text);
        }
        if &membership != expected {
            return Err(invalid("exon history and typed model memberships differ"));
        }
        self.sets += 1;
        Ok(loaded.into_order())
    }
}

fn load_orders(
    node: &Node,
    path: &mut Vec<String>,
    history: &mut History,
    orders: &mut HashMap<Vec<String>, Vec<String>>,
) -> Result<()> {
    for (key, member) in &node.members {
        path.push(key.clone());
        match member {
            Member::Node(child) => load_orders(child, path, history, orders)?,
            Member::Transcript(transcript) => {
                path.push("exons".into());
                let values = history.next(path, &transcript.exons)?;
                path.pop();
                if orders.insert(path.clone(), values).is_some() {
                    return Err(invalid("duplicate transcript metadata path"));
                }
            }
            _ => {}
        }
        path.pop();
    }
    Ok(())
}

#[derive(Default)]
struct Accumulated {
    positions: Vec<(i64, i64, String)>,
    position_index: HashMap<String, usize>,
    exons: Vec<(String, BTreeSet<String>)>,
    exon_index: HashMap<String, usize>,
}
impl Accumulated {
    fn add(&mut self, domain: &str, minimum: &Minimum, maximum: i64, exons: &[String]) {
        if let Minimum::Position(start) = *minimum {
            if start <= maximum {
                let position = *self
                    .position_index
                    .entry(domain.to_owned())
                    .or_insert_with(|| {
                        self.positions.push((start, maximum, domain.to_owned()));
                        self.positions.len() - 1
                    });
                self.positions[position] = (start, maximum, domain.to_owned());
            }
        }
        for exon in exons {
            let exon = exon.split('/').next().unwrap();
            let index = *self.exon_index.entry(exon.to_owned()).or_insert_with(|| {
                self.exons.push((exon.to_owned(), BTreeSet::new()));
                self.exons.len() - 1
            });
            self.exons[index].1.insert(domain.to_owned());
        }
    }
    fn finish(mut self) -> PlainTranscript {
        // Python's sort is stable and compares only aa_start, leaving ties in
        // first-domain-assignment order. Dictionary overwrites retain position.
        self.positions.sort_by_key(|value| value.0);
        PlainTranscript {
            pp_domains: self.positions,
            exon_domains: self
                .exons
                .into_iter()
                .map(|(key, values)| (key, values.into_iter().collect()))
                .collect(),
        }
    }
}

fn collect(
    node: &Node,
    path: &mut Vec<String>,
    orders: &HashMap<Vec<String>, Vec<String>>,
    transcripts: &mut HashMap<String, Accumulated>,
) -> Result<()> {
    for (key, member) in &node.members {
        if let Member::Node(child) = member {
            path.push(key.clone());
            if child.members.iter().any(|(name, _)| name == "distribution") {
                let domain = path.join(":");
                // Source processing stops at this domain, even when it has
                // nested domain dictionaries: it visits only its own ENSTs.
                for (transcript, member) in &child.members {
                    if !transcript.starts_with("ENST") {
                        continue;
                    }
                    if let Member::Transcript(value) = member {
                        path.push(transcript.clone());
                        let exons = orders
                            .get(path)
                            .ok_or_else(|| invalid("transcript history absent"))?;
                        if matches!(value.minimum, Minimum::Position(start) if start <= value.maximum)
                            || !exons.is_empty()
                        {
                            transcripts.entry(transcript.clone()).or_default().add(
                                &domain,
                                &value.minimum,
                                value.maximum,
                                exons,
                            );
                        }
                        path.pop();
                    }
                }
            } else {
                collect(child, path, orders, transcripts)?;
            }
            path.pop();
        }
    }
    Ok(())
}

pub(crate) fn build(
    statistics: &Path,
    exon_history: &Path,
    source_builder: &str,
    seed: u32,
    output: &Path,
) -> Result<(usize, usize)> {
    ensure_distinct_paths(statistics, output)?;
    ensure_distinct_paths(exon_history, output)?;
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
    let identities = [statistics, exon_history]
        .map(|path| std::fs::metadata(path))
        .into_iter()
        .collect::<std::io::Result<Vec<_>>>()?;
    let mut history = History {
        input: BufReader::new(File::open(exon_history)?),
        seed,
        sets: 0,
    };
    let mut transcripts = HashMap::new();
    let mut visited = 0;
    // Reuse the existing complete typed/ordinal validator. Only one compressed
    // gene model is held at a time; residue sample counts are never expanded.
    store::scan_genes(statistics, |gene, ordinal, node| {
        if ordinal != visited as u64 || genes.get(visited) != Some(&gene) {
            return Err(invalid("typed gene inventory or order differs"));
        }
        let mut path = vec![gene];
        let mut orders = HashMap::new();
        load_orders(&node, &mut path, &mut history, &mut orders)?;
        for (database, member) in &node.members {
            let Member::Node(source) = member else {
                return Err(invalid("gene source is not a domain dictionary"));
            };
            path.push(database.clone());
            collect(source, &mut path, &orders, &mut transcripts)?;
            path.pop();
        }
        visited += 1;
        Ok(())
    })?;
    if visited != genes.len() || history.input.read_line(&mut String::new())? != 0 {
        return Err(invalid("typed model/exon history coverage differs"));
    }
    let mut provenance = vec![
        Provenance {
            key: "procedure".into(),
            value: "retained-input-transcript-map-v1".into(),
        },
        Provenance {
            key: "python_hash_seed".into(),
            value: seed.to_string(),
        },
        Provenance {
            key: "statistics_builder_sha256".into(),
            value: source_builder.to_owned(),
        },
    ];
    for ((path, before), label) in [statistics, exon_history]
        .into_iter()
        .zip(identities)
        .zip(["statistics", "exon_history"])
    {
        let digest = sha256_file(path)?;
        let after = std::fs::metadata(path)?;
        if before.len() != after.len() || before.modified()? != after.modified()? {
            return Err(invalid("input changed during transcript projection"));
        }
        provenance.push(Provenance {
            key: format!("{label}_sha256"),
            value: digest,
        });
    }
    let count = write_entries(
        output,
        transcripts
            .into_iter()
            .map(|(key, value)| (key, value.finish()))
            .collect(),
        provenance,
    )?;
    Ok((count, history.sets))
}
