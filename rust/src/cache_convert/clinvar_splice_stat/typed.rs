//! Schema 2 v2 retains source-supported missing fields and per-record SpliceAI
//! dictionaries. The v1 structs and encoder remain unchanged. No score coercion
//! happens here; consumers decide when a field is actually evaluated.
use super::{
    ClinvarSpliceCache, ClinvarSpliceStat, PlainRecord, RECORD_FIELDS, ai_names_valid, encode, odd,
    text,
};
use crate::acmg_step::pickle_dict_stream::{Item, TextDict, visit_top_level_dict};
use crate::cache_convert::Provenance;
use crate::cache_convert::common::{ExactValue, StringPool, write_entry};
use crate::cli::{Result, invalid};
use crate::rust_cache::{MappedCache, Schema, find_sorted, write_cache_checked};
use rkyv::{Archive, Deserialize, Serialize, rancor};
use std::collections::HashSet;
use std::io::{Read, Write};
use std::path::Path;

pub(crate) const SCHEMA: Schema = Schema {
    name: "clinvar_splice_stat",
    id: 2,
    version: 2,
};

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Record {
    pub(crate) chrom: String,
    pub(crate) pos: i64,
    pub(crate) reference: String,
    pub(crate) alternative: Option<String>,
    pub(crate) exon: Option<String>,
    pub(crate) intron: Option<String>,
    pub(crate) hgvsc: Option<String>,
    pub(crate) consequence: Option<String>,
    pub(crate) clinvar_sig: String,
    pub(crate) clinvar_review: String,
    pub(crate) splice_ai: Vec<(String, String)>,
}

#[derive(Archive, Serialize, Deserialize, Debug)]
struct EncodedRecord {
    chrom: u32,
    pos: i64,
    reference: String,
    alternative: Option<String>,
    exon: Option<u32>,
    intron: Option<u32>,
    hgvsc: Option<String>,
    consequence: Option<u32>,
    clinvar_sig: u32,
    clinvar_review: u32,
    splice_ai: Vec<(u32, u32)>,
}
#[derive(Archive, Serialize, Deserialize, Debug)]
struct Root {
    provenance: Vec<Provenance>,
    texts: Vec<String>,
    transcripts: Vec<String>,
    entries: Vec<Vec<EncodedRecord>>,
}

fn nullable(dict: &TextDict, key: &str, file: &str, transcript: &str) -> Result<Option<String>> {
    match dict.get(key) {
        Some(Item::None) => Ok(None),
        Some(Item::Str(s)) => Ok(Some(s.to_string())),
        _ => Err(odd(
            file,
            transcript,
            format!("{key} is absent or not str/None"),
        )),
    }
}

pub(super) fn plain_transcript(file: &str, transcript: &str, input: &Item) -> Result<Vec<Record>> {
    let Item::List(records) = input else {
        return Err(odd(file, transcript, "value is not a list"));
    };
    records
        .borrow()
        .iter()
        .map(|input| {
            let Item::Dict(input) = input else {
                return Err(odd(file, transcript, "record is not a dict"));
            };
            let input = input.borrow();
            if input.len() != RECORD_FIELDS.len()
                || input
                    .entries
                    .iter()
                    .any(|(k, _)| !RECORD_FIELDS.contains(&k.as_ref()))
            {
                return Err(odd(
                    file,
                    transcript,
                    "record must have exactly the eleven known fields",
                ));
            }
            let pos = match input.get("pos") {
                Some(Item::Int(i)) => *i,
                _ => return Err(odd(file, transcript, "pos is not an int")),
            };
            let Some(Item::Dict(ai)) = input.get("splice_ai") else {
                return Err(odd(file, transcript, "splice_ai is not a dict"));
            };
            let splice_ai = ai
                .borrow()
                .entries
                .iter()
                .map(|(key, value)| {
                    if !key.starts_with("SpliceAI_pred_") {
                        return Err(odd(file, transcript, "unknown splice_ai field prefix"));
                    }
                    let Item::Str(value) = value else {
                        return Err(odd(file, transcript, "splice_ai value is not str"));
                    };
                    Ok((key.to_string(), value.to_string()))
                })
                .collect::<Result<Vec<_>>>()?;
            Ok(Record {
                chrom: text(&input, "chrom", file, transcript)?,
                pos,
                reference: text(&input, "ref", file, transcript)?,
                alternative: nullable(&input, "alt", file, transcript)?,
                exon: nullable(&input, "exon", file, transcript)?,
                intron: nullable(&input, "intron", file, transcript)?,
                hgvsc: nullable(&input, "hgvsc", file, transcript)?,
                consequence: nullable(&input, "consequence", file, transcript)?,
                clinvar_sig: text(&input, "clinvar_sig", file, transcript)?,
                clinvar_review: text(&input, "clinvar_review", file, transcript)?,
                splice_ai,
            })
        })
        .collect()
}

/// Check v1 representability without normalizing values or dictionary order.
pub(super) fn legacy_names(records: &[Record], names: &mut Option<Vec<String>>) -> bool {
    for record in records {
        if [
            &record.alternative,
            &record.exon,
            &record.intron,
            &record.hgvsc,
            &record.consequence,
        ]
        .iter()
        .any(|s| s.is_none())
        {
            return false;
        }
        if !record.splice_ai.is_empty() {
            let present: Vec<_> = record.splice_ai.iter().map(|(k, _)| k.clone()).collect();
            if !ai_names_valid(&present) {
                return false;
            }
            match names {
                Some(previous) if *previous != present => return false,
                None => *names = Some(present),
                _ => {}
            }
        }
    }
    true
}

fn from_v1(record: PlainRecord, names: &[String]) -> Record {
    Record {
        chrom: record.chrom,
        pos: record.pos,
        reference: record.reference,
        alternative: Some(record.alternative),
        exon: Some(record.exon),
        intron: Some(record.intron),
        hgvsc: Some(record.hgvsc),
        consequence: Some(record.consequence),
        clinvar_sig: record.clinvar_sig,
        clinvar_review: record.clinvar_review,
        splice_ai: names.iter().cloned().zip(record.splice_ai).collect(),
    }
}

pub(crate) enum Cache {
    V1(ClinvarSpliceCache),
    V2(TypedCache),
}
pub(crate) struct TypedCache {
    mapped: MappedCache<Root>,
}

impl Cache {
    pub(crate) fn open(path: &str) -> Result<Self> {
        let mut file = std::fs::File::open(path)?;
        let mut header = [0; 20];
        let v2 = file.read_exact(&mut header).is_ok()
            && &header[..8] == crate::rust_cache::MAGIC
            && u32::from_le_bytes(header[12..16].try_into().unwrap()) == SCHEMA.id
            && u32::from_le_bytes(header[16..20].try_into().unwrap()) == SCHEMA.version;
        if !v2 {
            return ClinvarSpliceCache::open(path).map(Self::V1);
        }
        let mapped = MappedCache::<Root>::open_file(&file, path, SCHEMA)?;
        let root = mapped.root();
        if root.transcripts.len() != root.entries.len()
            || root
                .transcripts
                .windows(2)
                .any(|p| p[0].as_str() >= p[1].as_str())
        {
            return Err(invalid(
                "splice v2 transcripts must be sorted, unique and one per entry",
            ));
        }
        if root
            .texts
            .iter()
            .map(|s| s.as_str())
            .collect::<HashSet<_>>()
            .len()
            != root.texts.len()
        {
            return Err(invalid("splice v2 text pool contains duplicates"));
        }
        let valid = |v: &rkyv::Archived<u32>| (v.to_native() as usize) < root.texts.len();
        for record in root.entries.iter().flat_map(|e| e.iter()) {
            if ![&record.chrom, &record.clinvar_sig, &record.clinvar_review]
                .into_iter()
                .all(valid)
                || ![&record.exon, &record.intron, &record.consequence]
                    .into_iter()
                    .filter_map(|o| o.as_ref())
                    .all(valid)
                || !record
                    .splice_ai
                    .iter()
                    .all(|pair| valid(&pair.0) && valid(&pair.1))
            {
                return Err(invalid("splice v2 text index is out of range"));
            }
            let names: Vec<_> = record
                .splice_ai
                .iter()
                .map(|pair| root.texts[pair.0.to_native() as usize].as_str())
                .collect();
            if names.iter().any(|k| !k.starts_with("SpliceAI_pred_"))
                || names.iter().collect::<HashSet<_>>().len() != names.len()
            {
                return Err(invalid("splice v2 has duplicate or unknown splice_ai keys"));
            }
        }
        Ok(Self::V2(TypedCache { mapped }))
    }
    pub(crate) fn version(&self) -> u32 {
        match self {
            Self::V1(_) => 1,
            Self::V2(_) => 2,
        }
    }
    pub(crate) fn keys(&self) -> &[rkyv::string::ArchivedString] {
        match self {
            Self::V1(c) => c.keys(),
            Self::V2(c) => &c.mapped.root().transcripts,
        }
    }
    pub(crate) fn len(&self) -> usize {
        self.keys().len()
    }
    pub(crate) fn plain(&self, key: &str) -> Option<Vec<Record>> {
        match self {
            Self::V1(cache) => {
                let names: Vec<_> = cache
                    .splice_ai_names()
                    .iter()
                    .map(|s| s.as_str().to_owned())
                    .collect();
                cache
                    .plain(key)
                    .map(|records| records.into_iter().map(|r| from_v1(r, &names)).collect())
            }
            Self::V2(cache) => {
                let root = cache.mapped.root();
                let text = |id: &rkyv::Archived<u32>| {
                    root.texts[id.to_native() as usize].as_str().to_owned()
                };
                find_sorted(&root.transcripts, key).map(|i| {
                    root.entries[i]
                        .iter()
                        .map(|r| Record {
                            chrom: text(&r.chrom),
                            pos: r.pos.to_native(),
                            reference: r.reference.as_str().to_owned(),
                            alternative: r.alternative.as_ref().map(|s| s.as_str().to_owned()),
                            exon: r.exon.as_ref().map(text),
                            intron: r.intron.as_ref().map(text),
                            hgvsc: r.hgvsc.as_ref().map(|s| s.as_str().to_owned()),
                            consequence: r.consequence.as_ref().map(text),
                            clinvar_sig: text(&r.clinvar_sig),
                            clinvar_review: text(&r.clinvar_review),
                            splice_ai: r
                                .splice_ai
                                .iter()
                                .map(|pair| (text(&pair.0), text(&pair.1)))
                                .collect(),
                        })
                        .collect()
                })
            }
        }
    }
}

pub(super) fn exact(records: Vec<Record>) -> ExactValue {
    let text = ExactValue::Text;
    let nullable = |s: Option<String>| s.map(ExactValue::Text).unwrap_or(ExactValue::None);
    ExactValue::List(
        records
            .into_iter()
            .map(|r| {
                ExactValue::Dict(vec![
                    (text("chrom".into()), text(r.chrom)),
                    (text("pos".into()), ExactValue::Int(r.pos)),
                    (text("ref".into()), text(r.reference)),
                    (text("alt".into()), nullable(r.alternative)),
                    (text("exon".into()), nullable(r.exon)),
                    (text("intron".into()), nullable(r.intron)),
                    (text("hgvsc".into()), nullable(r.hgvsc)),
                    (text("consequence".into()), nullable(r.consequence)),
                    (text("clinvar_sig".into()), text(r.clinvar_sig)),
                    (text("clinvar_review".into()), text(r.clinvar_review)),
                    (
                        text("splice_ai".into()),
                        ExactValue::Dict(
                            r.splice_ai
                                .into_iter()
                                .map(|(k, v)| (text(k), text(v)))
                                .collect(),
                        ),
                    ),
                ])
            })
            .collect(),
    )
}

pub(super) fn export(path: &str, out: &mut dyn Write) -> Result<usize> {
    let cache = Cache::open(path)?;
    for key in cache.keys() {
        write_entry(
            out,
            &ExactValue::Text(key.as_str().to_owned()),
            &exact(cache.plain(key.as_str()).unwrap()),
        )?;
    }
    Ok(cache.len())
}

pub(super) fn pickle_version(path: &str) -> Result<u32> {
    let mut names = None;
    let mut legacy = true;
    visit_top_level_dict(
        path,
        |_| true,
        |key, value| {
            let records = plain_transcript(path, &key, &value)?;
            legacy &= legacy_names(&records, &mut names);
            Ok(())
        },
    )?;
    Ok(if legacy { 1 } else { 2 })
}

pub(crate) fn write_records(
    output: &Path,
    entries: Vec<(String, Vec<Record>)>,
    provenance: Vec<Provenance>,
) -> Result<usize> {
    write_checked(output, entries, provenance, |_| Ok(()))
}
pub(super) fn write_checked(
    output: &Path,
    entries: Vec<(String, Vec<Record>)>,
    provenance: Vec<Provenance>,
    check: impl FnOnce(&Cache) -> Result<()>,
) -> Result<usize> {
    let mut names = None;
    let legacy = entries.iter().all(|(_, r)| legacy_names(r, &mut names));
    let count = entries.len();
    let mut pool = StringPool::default();
    let (schema, bytes) = if legacy {
        let mut encoded = Vec::with_capacity(count);
        for (key, records) in &entries {
            let records = records
                .iter()
                .map(|r| {
                    encode(
                        PlainRecord {
                            chrom: r.chrom.clone(),
                            pos: r.pos,
                            reference: r.reference.clone(),
                            alternative: r.alternative.clone().unwrap(),
                            exon: r.exon.clone().unwrap(),
                            intron: r.intron.clone().unwrap(),
                            hgvsc: r.hgvsc.clone().unwrap(),
                            consequence: r.consequence.clone().unwrap(),
                            clinvar_sig: r.clinvar_sig.clone(),
                            clinvar_review: r.clinvar_review.clone(),
                            splice_ai: r.splice_ai.iter().map(|(_, v)| v.clone()).collect(),
                        },
                        &mut pool,
                    )
                })
                .collect::<Result<Vec<_>>>()?;
            encoded.push((key.clone(), records));
        }
        encoded.sort_by(|a, b| a.0.cmp(&b.0));
        let (transcripts, entries) = encoded.into_iter().unzip();
        let value = ClinvarSpliceStat {
            provenance,
            texts: pool.into_texts(),
            splice_ai_names: names.unwrap_or_default(),
            transcripts,
            entries,
        };
        (
            super::SCHEMA,
            rkyv::to_bytes::<rancor::Error>(&value)?.to_vec(),
        )
    } else {
        let mut encoded = Vec::with_capacity(count);
        for (key, records) in &entries {
            let mut encoded_records = Vec::with_capacity(records.len());
            for r in records {
                encoded_records.push(EncodedRecord {
                    chrom: pool.intern(&r.chrom)?,
                    pos: r.pos,
                    reference: r.reference.clone(),
                    alternative: r.alternative.clone(),
                    exon: r.exon.as_deref().map(|s| pool.intern(s)).transpose()?,
                    intron: r.intron.as_deref().map(|s| pool.intern(s)).transpose()?,
                    hgvsc: r.hgvsc.clone(),
                    consequence: r
                        .consequence
                        .as_deref()
                        .map(|s| pool.intern(s))
                        .transpose()?,
                    clinvar_sig: pool.intern(&r.clinvar_sig)?,
                    clinvar_review: pool.intern(&r.clinvar_review)?,
                    splice_ai: r
                        .splice_ai
                        .iter()
                        .map(|(k, v)| Ok((pool.intern(k)?, pool.intern(v)?)))
                        .collect::<Result<_>>()?,
                });
            }
            encoded.push((key.clone(), encoded_records));
        }
        encoded.sort_by(|a, b| a.0.cmp(&b.0));
        let (transcripts, entries) = encoded.into_iter().unzip();
        let value = Root {
            provenance,
            texts: pool.into_texts(),
            transcripts,
            entries,
        };
        (SCHEMA, rkyv::to_bytes::<rancor::Error>(&value)?.to_vec())
    };
    write_cache_checked(output, schema, &bytes, |staged| {
        let cache = Cache::open(&staged.to_string_lossy())?;
        if cache.len() != count
            || entries
                .iter()
                .any(|(k, v)| cache.plain(k).as_ref() != Some(v))
        {
            return Err(invalid("direct splice values differ after writing"));
        }
        check(&cache)
    })?;
    Ok(count)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn path(label: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "splice-v2-{label}-{}-{}.rkyv",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    }
    fn record() -> Record {
        Record {
            chrom: "1".into(),
            pos: 42,
            reference: "A".into(),
            alternative: Some("T".into()),
            exon: Some("1/2".into()),
            intron: Some("".into()),
            hgvsc: Some("t:c.1+1A>T".into()),
            consequence: Some("splice_region_variant".into()),
            clinvar_sig: "Pathogenic".into(),
            clinvar_review: "reviewed_by_expert_panel".into(),
            splice_ai: vec![],
        }
    }
    #[test]
    fn version_selection_preserves_v1_and_nullable_pair_order() {
        let old = path("v1");
        let entries = vec![("T".into(), vec![record()])];
        write_records(&old, entries.clone(), vec![]).unwrap();
        let cache = Cache::open(old.to_str().unwrap()).unwrap();
        assert_eq!(cache.version(), 1);
        assert_eq!(cache.plain("T").unwrap(), entries[0].1);
        assert!(ClinvarSpliceCache::open(old.to_str().unwrap()).is_ok());
        let mut nullable = record();
        nullable.alternative = None;
        nullable.hgvsc = None;
        nullable.consequence = None;
        nullable.splice_ai = vec![
            ("SpliceAI_pred_DS_DL".into(), ".9".into()),
            ("SpliceAI_pred_NEW".into(), "literal".into()),
        ];
        let entries = vec![("T".into(), vec![nullable, record()])];
        let new = path("v2");
        write_records(&new, entries.clone(), vec![]).unwrap();
        let cache = Cache::open(new.to_str().unwrap()).unwrap();
        assert_eq!(cache.version(), 2);
        assert_eq!(cache.plain("T").unwrap(), entries[0].1);
        assert!(ClinvarSpliceCache::open(new.to_str().unwrap()).is_err());
        std::fs::remove_file(old).unwrap();
        std::fs::remove_file(new).unwrap();
    }
    #[test]
    fn invalid_v2_dictionary_keeps_existing_output() {
        let output = path("invalid");
        std::fs::write(&output, b"preserved").unwrap();
        for keys in [
            vec!["other"],
            vec!["SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AG"],
        ] {
            let mut value = record();
            value.splice_ai = keys.into_iter().map(|k| (k.into(), ".5".into())).collect();
            assert!(write_records(&output, vec![("T".into(), vec![value])], vec![]).is_err());
            assert_eq!(std::fs::read(&output).unwrap(), b"preserved");
        }
        std::fs::remove_file(output).unwrap();
    }
    #[test]
    fn invalid_v2_text_index_and_future_schema_are_refused() {
        let output = path("invalid-id");
        let root = Root {
            provenance: vec![],
            texts: vec!["1".into()],
            transcripts: vec!["T".into()],
            entries: vec![vec![EncodedRecord {
                chrom: 5,
                pos: 1,
                reference: "A".into(),
                alternative: None,
                exon: None,
                intron: None,
                hgvsc: None,
                consequence: None,
                clinvar_sig: 0,
                clinvar_review: 0,
                splice_ai: vec![],
            }]],
        };
        let bytes = rkyv::to_bytes::<rancor::Error>(&root).unwrap();
        crate::rust_cache::write_cache(&output, SCHEMA, &bytes).unwrap();
        assert!(Cache::open(output.to_str().unwrap()).is_err());
        crate::rust_cache::write_cache(
            &output,
            Schema {
                version: 3,
                ..SCHEMA
            },
            &bytes,
        )
        .unwrap();
        assert!(Cache::open(output.to_str().unwrap()).is_err());
        std::fs::remove_file(output).unwrap();
    }
}
