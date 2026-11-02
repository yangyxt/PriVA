//! Schema 3, version 1: AlphaMissense transcript-to-domain lookup.
//!
//! Python: `{str: {"pp_domains": [(int, int, str), ...],
//! "exon_domains": {str: [str, ...]}}}`. The two record keys are required;
//! protein tuples, domain lists and the free-key exon dictionary retain their order
//! and multiplicity. Transcript lookup keys alone are sorted. No absent field is
//! replaced with an empty value and no unrecognized field is discarded.
use super::Provenance;
use super::common::{ExactValue, StringPool, provenance, write_entry};
use crate::acmg_step::pickle_dict_stream::{Item, visit_top_level_dict};
use crate::cli::{Result, invalid};
use crate::rust_cache::{
    MappedCache, Schema, ensure_distinct_paths, find_sorted, write_cache_checked,
};
use rkyv::{Archive, Deserialize, Serialize, rancor};
use std::collections::HashSet;
use std::io::Write;
use std::path::Path;
use std::time::Instant;

pub(crate) const ALPHAMISSENSE_TRANX_DOMAIN_MAP: Schema = Schema {
    name: "alphamissense_tranx_domain_map",
    id: 3,
    version: 1,
};

/// One free-key exon-dictionary entry. `domains` is a Python list of text.
#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct ExonDomains {
    pub(crate) exon: u32,
    pub(crate) domains: Vec<u32>,
}

/// A Python 3-tuple `(int start, int end, str domain)`, in list order.
#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct ProteinDomain {
    pub(crate) start: i64,
    pub(crate) end: i64,
    pub(crate) domain: u32,
}

/// Fixed-key Python record; both fields are required.
#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct TranscriptDomains {
    pub(crate) exon_domains: Vec<ExonDomains>,
    pub(crate) pp_domains: Vec<ProteinDomain>,
}

#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct TranscriptDomainMap {
    pub(crate) provenance: Vec<Provenance>,
    pub(crate) texts: Vec<String>,
    /// UTF-8 byte sorted, unique, one key per corresponding entry.
    pub(crate) transcripts: Vec<String>,
    pub(crate) entries: Vec<TranscriptDomains>,
}

#[derive(Debug, PartialEq)]
struct PlainTranscript {
    exon_domains: Vec<(String, Vec<String>)>,
    pp_domains: Vec<(i64, i64, String)>,
}

fn plain_transcript(path: &str, transcript: &str, value: &Item) -> Result<PlainTranscript> {
    let odd = |what: &str| invalid(format!("{path}: transcript {transcript}: {what}; refused"));
    let Item::Dict(record) = value else {
        return Err(odd("the value is not a dict"));
    };
    let record = record.borrow();
    if record.len() != 2
        || record.get("pp_domains").is_none()
        || record.get("exon_domains").is_none()
    {
        return Err(odd(
            "the record must contain exactly pp_domains and exon_domains",
        ));
    }
    let Some(Item::List(protein)) = record.get("pp_domains") else {
        return Err(odd("pp_domains is not a list"));
    };
    let mut pp_domains = Vec::with_capacity(protein.borrow().len());
    for item in protein.borrow().iter() {
        let Item::Tuple(parts) = item else {
            return Err(odd("a pp_domains item is not a tuple"));
        };
        let [Item::Int(start), Item::Int(end), Item::Str(domain)] = parts.as_slice() else {
            return Err(odd("a pp_domains item is not a 3-tuple of (int, int, str)"));
        };
        pp_domains.push((*start, *end, domain.to_string()));
    }
    let Some(Item::Dict(exons)) = record.get("exon_domains") else {
        return Err(odd("exon_domains is not a dict"));
    };
    let mut exon_domains = Vec::with_capacity(exons.borrow().len());
    for (exon, domains) in &exons.borrow().entries {
        let Item::List(domains) = domains else {
            return Err(odd("an exon_domains value is not a list"));
        };
        let domains = domains
            .borrow()
            .iter()
            .map(|domain| match domain {
                Item::Str(domain) => Ok(domain.to_string()),
                _ => Err(odd("an exon_domains list item is not text")),
            })
            .collect::<Result<Vec<_>>>()?;
        exon_domains.push((exon.to_string(), domains));
    }
    Ok(PlainTranscript {
        exon_domains,
        pp_domains,
    })
}

fn read_pickle(path: &str) -> Result<Vec<(String, PlainTranscript)>> {
    let mut entries = Vec::new();
    let distinct = visit_top_level_dict(
        path,
        |_| true,
        |transcript, value| {
            entries.push((
                transcript.to_string(),
                plain_transcript(path, &transcript, &value)?,
            ));
            Ok(())
        },
    )?;
    if distinct != entries.len() {
        return Err(invalid(format!(
            "{path}: a transcript is set twice; refused"
        )));
    }
    Ok(entries)
}

fn build(
    mut entries: Vec<(String, PlainTranscript)>,
    provenance: Vec<Provenance>,
) -> Result<Vec<u8>> {
    entries.sort_by(|a, b| a.0.cmp(&b.0));
    let mut pool = StringPool::default();
    let mut transcripts = Vec::with_capacity(entries.len());
    let mut records = Vec::with_capacity(entries.len());
    for (transcript, entry) in entries {
        transcripts.push(transcript);
        let exon_domains = entry
            .exon_domains
            .into_iter()
            .map(|(exon, domains)| {
                Ok(ExonDomains {
                    exon: pool.intern(&exon)?,
                    domains: domains
                        .into_iter()
                        .map(|domain| pool.intern(&domain))
                        .collect::<Result<_>>()?,
                })
            })
            .collect::<Result<_>>()?;
        let pp_domains = entry
            .pp_domains
            .into_iter()
            .map(|(start, end, domain)| {
                Ok(ProteinDomain {
                    start,
                    end,
                    domain: pool.intern(&domain)?,
                })
            })
            .collect::<Result<_>>()?;
        records.push(TranscriptDomains {
            exon_domains,
            pp_domains,
        });
    }
    Ok(rkyv::to_bytes::<rancor::Error>(&TranscriptDomainMap {
        provenance,
        texts: pool.into_texts(),
        transcripts,
        entries: records,
    })?
    .to_vec())
}

/// Read-only map with bytechecked pointers, sorted lookup keys and checked pool indexes.
/// Lookup and text access borrow the map; opening does not copy its string pool.
pub(crate) struct TranscriptDomainCache {
    mapped: MappedCache<TranscriptDomainMap>,
}

impl TranscriptDomainCache {
    pub(crate) fn open(path: &str) -> Result<Self> {
        let mapped =
            MappedCache::<TranscriptDomainMap>::open(path, ALPHAMISSENSE_TRANX_DOMAIN_MAP)?;
        let root = mapped.root();
        if root.transcripts.len() != root.entries.len()
            || root
                .transcripts
                .windows(2)
                .any(|pair| pair[0].as_str() >= pair[1].as_str())
        {
            return Err(invalid(format!(
                "{path}: transcript keys are not sorted and unique, or not one per entry"
            )));
        }
        let count = root.texts.len();
        let in_range = |id: u32| (id as usize) < count;
        for entry in root.entries.iter() {
            let mut exons = HashSet::new();
            for exon in entry.exon_domains.iter() {
                if !in_range(exon.exon.to_native())
                    || exon.domains.iter().any(|id| !in_range(id.to_native()))
                {
                    return Err(invalid(format!(
                        "{path}: an exon or domain text index is out of range"
                    )));
                }
                if !exons.insert(root.texts[exon.exon.to_native() as usize].as_str()) {
                    return Err(invalid(format!(
                        "{path}: a transcript has a repeated exon key"
                    )));
                }
            }
            if entry
                .pp_domains
                .iter()
                .any(|domain| !in_range(domain.domain.to_native()))
            {
                return Err(invalid(format!(
                    "{path}: a protein-domain text index is out of range"
                )));
            }
        }
        Ok(Self { mapped })
    }

    pub(crate) fn len(&self) -> usize {
        self.mapped.root().transcripts.len()
    }

    pub(crate) fn get(&self, transcript: &str) -> Option<&ArchivedTranscriptDomains> {
        let root = self.mapped.root();
        find_sorted(&root.transcripts, transcript).map(|at| &root.entries[at])
    }

    pub(crate) fn text(&self, id: &rkyv::Archived<u32>) -> &str {
        self.mapped.root().texts[id.to_native() as usize].as_str()
    }

    fn plain(&self, transcript: &str) -> Option<PlainTranscript> {
        self.get(transcript).map(|entry| PlainTranscript {
            exon_domains: entry
                .exon_domains
                .iter()
                .map(|exon| {
                    (
                        self.text(&exon.exon).to_owned(),
                        exon.domains
                            .iter()
                            .map(|id| self.text(id).to_owned())
                            .collect(),
                    )
                })
                .collect(),
            pp_domains: entry
                .pp_domains
                .iter()
                .map(|domain| {
                    (
                        domain.start.to_native(),
                        domain.end.to_native(),
                        self.text(&domain.domain).to_owned(),
                    )
                })
                .collect(),
        })
    }
}

fn check(pickle: &str, output: &Path, count: usize) -> Result<()> {
    let path = output
        .to_str()
        .ok_or_else(|| invalid("the cache path is not UTF-8"))?;
    let cache = TranscriptDomainCache::open(path)?;
    if cache.len() != count {
        return Err(invalid("the converted transcript count differs"));
    }
    let mut compared = 0;
    let distinct = visit_top_level_dict(
        pickle,
        |_| true,
        |transcript, value| {
            let expected = plain_transcript(pickle, &transcript, &value)?;
            if cache.plain(&transcript).as_ref() != Some(&expected) {
                return Err(invalid(format!(
                    "{path}: transcript {transcript} differs from the pickle"
                )));
            }
            compared += 1;
            Ok(())
        },
    )?;
    if distinct != count || compared != count {
        return Err(invalid("the pickle changed while it was converted"));
    }
    Ok(())
}

pub(crate) fn convert(pickle: &str, output: &Path) -> Result<String> {
    ensure_distinct_paths(Path::new(pickle), output)?;
    let start = Instant::now();
    let entries = read_pickle(pickle)?;
    let count = entries.len();
    let protein_domains: usize = entries
        .iter()
        .map(|(_, entry)| entry.pp_domains.len())
        .sum();
    let exons: usize = entries
        .iter()
        .map(|(_, entry)| entry.exon_domains.len())
        .sum();
    let archive = build(entries, provenance(pickle)?)?;
    write_cache_checked(output, ALPHAMISSENSE_TRANX_DOMAIN_MAP, &archive, |staged| {
        check(pickle, staged, count)
    })?;
    Ok(format!(
        "{{\"cache\":\"alphamissense_tranx_domain_map\",\"transcripts\":{count},\"exons\":{exons},\"protein_domains\":{protein_domains},\"checked_transcripts\":{count},\"bytes\":{},\"elapsed_s\":{:.3}}}",
        std::fs::metadata(output)?.len(),
        start.elapsed().as_secs_f64()
    ))
}

/// Complete typed values for an independent Python comparison, one sorted transcript per line.
pub(crate) fn export(path: &str, out: &mut dyn Write) -> Result<usize> {
    let cache = TranscriptDomainCache::open(path)?;
    for (key, entry) in cache
        .mapped
        .root()
        .transcripts
        .iter()
        .zip(cache.mapped.root().entries.iter())
    {
        let pp_domains = ExactValue::List(
            entry
                .pp_domains
                .iter()
                .map(|domain| {
                    ExactValue::Tuple(vec![
                        ExactValue::Int(domain.start.to_native()),
                        ExactValue::Int(domain.end.to_native()),
                        ExactValue::Text(cache.text(&domain.domain).to_owned()),
                    ])
                })
                .collect(),
        );
        let exon_domains = ExactValue::Dict(
            entry
                .exon_domains
                .iter()
                .map(|exon| {
                    (
                        ExactValue::Text(cache.text(&exon.exon).to_owned()),
                        ExactValue::List(
                            exon.domains
                                .iter()
                                .map(|id| ExactValue::Text(cache.text(id).to_owned()))
                                .collect(),
                        ),
                    )
                })
                .collect(),
        );
        write_entry(
            out,
            &ExactValue::Text(key.as_str().to_owned()),
            &ExactValue::Dict(vec![
                (ExactValue::Text("pp_domains".to_owned()), pp_domains),
                (ExactValue::Text("exon_domains".to_owned()), exon_domains),
            ]),
        )?;
    }
    Ok(cache.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    // Python protocol 4: tx2 has out-of-coordinate-order tuples, repeated list
    // entries and exon insertion order 2,1; tx1 has two empty containers.
    const GOOD: &[u8] = b"\x80\x04\x95\x71\x00\x00\x00\x00\x00\x00\x00\x7d\x94\x28\x8c\x03\x74\x78\x32\x94\x7d\x94\x28\x8c\x0a\x70\x70\x5f\x64\x6f\x6d\x61\x69\x6e\x73\x94\x5d\x94\x28\x4b\x07\x4b\x0b\x8c\x02\x64\x32\x94\x87\x94\x4b\x02\x4b\x09\x8c\x02\x64\x31\x94\x87\x94\x68\x06\x65\x8c\x0c\x65\x78\x6f\x6e\x5f\x64\x6f\x6d\x61\x69\x6e\x73\x94\x7d\x94\x28\x8c\x01\x32\x94\x5d\x94\x28\x68\x05\x68\x05\x68\x07\x65\x8c\x01\x31\x94\x5d\x94\x75\x75\x8c\x03\x74\x78\x31\x94\x7d\x94\x28\x68\x03\x5d\x94\x68\x09\x7d\x94\x75\x75\x2e";

    fn temporary(name: &str, bytes: &[u8]) -> String {
        let path = std::env::temp_dir().join(format!(
            "priva_transcript_domains_{}_{name}",
            std::process::id()
        ));
        std::fs::write(&path, bytes).unwrap();
        path.to_str().unwrap().to_owned()
    }

    #[test]
    fn exact_roundtrip_preserves_tuples_list_multiplicity_and_exon_order() {
        let pickle = temporary("good.pkl", GOOD);
        let output = temporary("good.rkyv", b"");
        convert(&pickle, Path::new(&output)).unwrap();
        let cache = TranscriptDomainCache::open(&output).unwrap();
        assert_eq!(cache.len(), 2);
        assert!(cache.get("missing").is_none());
        let entry = cache.plain("tx2").unwrap();
        assert_eq!(
            entry.exon_domains,
            vec![
                ("2".into(), vec!["d2".into(), "d2".into(), "d1".into()]),
                ("1".into(), vec![])
            ]
        );
        assert_eq!(
            entry.pp_domains,
            vec![
                (7, 11, "d2".into()),
                (2, 9, "d1".into()),
                (7, 11, "d2".into())
            ]
        );
        let mut exported = Vec::new();
        assert_eq!(export(&output, &mut exported).unwrap(), 2);
        let text = String::from_utf8(exported).unwrap();
        assert!(text.find("tx1").unwrap() < text.find("tx2").unwrap());
        let again = temporary("again.rkyv", b"");
        convert(&pickle, Path::new(&again)).unwrap();
        assert_eq!(
            std::fs::read(&output).unwrap(),
            std::fs::read(&again).unwrap()
        );
    }

    #[test]
    fn refuses_coercions_unknown_fields_and_duplicate_transcripts() {
        let cases: &[&[u8]] = &[
            // pp_domains tuple rather than list.
            b"\x80\x04}\x8c\x01t}(\x8c\x0app_domains)\x8c\x0cexon_domains}us.",
            // pp_domains item list rather than 3-tuple.
            b"\x80\x04}\x8c\x01t}(\x8c\x0app_domains]]a\x8c\x0cexon_domains}us.",
            // bool is not an integer endpoint.
            b"\x80\x04}\x8c\x01t}(\x8c\x0app_domains](\x88K\x02\x8c\x01dta\x8c\x0cexon_domains}us.",
            // exon_domains values must be lists.
            b"\x80\x04}\x8c\x01t}(\x8c\x0app_domains]\x8c\x0cexon_domains}\x8c\x011)sus.",
            // An unrecognized record field.
            b"\x80\x04}\x8c\x01t}(\x8c\x0app_domains]\x8c\x0cexon_domains}\x8c\x01xNus.",
            // Required field missing.
            b"\x80\x04}\x8c\x01t}\x8c\x0app_domains]ss.",
            // Same top key assigned twice; both valid records would decode.
            b"\x80\x04}(\x8c\x01t}(\x8c\x0app_domains]\x8c\x0cexon_domains}u\x8c\x01t}(\x8c\x0app_domains]\x8c\x0cexon_domains}uu.",
        ];
        let output = temporary("refuse.rkyv", b"unchanged");
        for (at, bytes) in cases.iter().enumerate() {
            let pickle = temporary(&format!("bad{at}.pkl"), bytes);
            assert!(convert(&pickle, Path::new(&output)).is_err(), "case {at}");
            assert_eq!(std::fs::read(&output).unwrap(), b"unchanged");
        }
    }

    #[test]
    fn reader_rejects_semantic_archive_damage() {
        let output = temporary("bad_index.rkyv", b"");
        let root = TranscriptDomainMap {
            provenance: vec![],
            texts: vec!["exon".into()],
            transcripts: vec!["t".into()],
            entries: vec![TranscriptDomains {
                exon_domains: vec![ExonDomains {
                    exon: 0,
                    domains: vec![1],
                }],
                pp_domains: vec![],
            }],
        };
        let archive = rkyv::to_bytes::<rancor::Error>(&root).unwrap();
        crate::rust_cache::write_cache(
            Path::new(&output),
            ALPHAMISSENSE_TRANX_DOMAIN_MAP,
            &archive,
        )
        .unwrap();
        assert!(
            TranscriptDomainCache::open(&output)
                .err()
                .unwrap()
                .to_string()
                .contains("out of range")
        );
    }
}
