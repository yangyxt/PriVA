//! Schema 2, version 1: `{transcript: [splice_record, ...]}`. The list order and
//! multiplicity are retained. The eleven fixed record keys may have any Python
//! insertion order and export in `RECORD_FIELDS` order. `splice_ai` is either
//! `{}` or the complete nine-key dictionary; its insertion order is retained in
//! `splice_ai_names` and must be the same in every nonempty record.
use super::Provenance;
use super::common::{ExactValue, StringPool, provenance, write_entry};
use crate::acmg_step::pickle_dict_stream::{Item, TextDict, visit_top_level_dict};
use crate::cli::{Result, invalid};
use crate::rust_cache::{
    MappedCache, Schema, ensure_distinct_paths, find_sorted, write_cache_checked,
};
use rkyv::{Archive, Deserialize, Serialize, rancor};
use std::collections::HashSet;
use std::io::Write;
use std::path::Path;
use std::time::Instant;

pub(crate) const SCHEMA: Schema = Schema {
    name: "clinvar_splice_stat",
    id: 2,
    version: 1,
};
const RECORD_FIELDS: [&str; 11] = [
    "chrom",
    "pos",
    "ref",
    "alt",
    "exon",
    "intron",
    "hgvsc",
    "consequence",
    "clinvar_sig",
    "clinvar_review",
    "splice_ai",
];
const SPLICE_AI_FIELDS: [&str; 9] = [
    "SpliceAI_pred_SYMBOL",
    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL",
    "SpliceAI_pred_DP_AG",
    "SpliceAI_pred_DP_AL",
    "SpliceAI_pred_DP_DG",
    "SpliceAI_pred_DP_DL",
];

#[derive(Archive, Serialize, Deserialize, Debug)]
pub(crate) struct SpliceRecord {
    pub(crate) chrom: u32,
    pub(crate) pos: i64,
    pub(crate) reference: String,
    pub(crate) alternative: String,
    pub(crate) exon: u32,
    pub(crate) intron: u32,
    pub(crate) hgvsc: String,
    pub(crate) consequence: u32,
    pub(crate) clinvar_sig: u32,
    pub(crate) clinvar_review: u32,
    pub(crate) splice_ai: Vec<u32>,
}

#[derive(Archive, Serialize, Deserialize, Debug)]
pub(crate) struct ClinvarSpliceStat {
    pub(crate) provenance: Vec<Provenance>,
    pub(crate) texts: Vec<String>,
    pub(crate) splice_ai_names: Vec<String>,
    pub(crate) transcripts: Vec<String>,
    pub(crate) entries: Vec<Vec<SpliceRecord>>,
}

#[derive(Debug, PartialEq)]
struct PlainRecord {
    chrom: String,
    pos: i64,
    reference: String,
    alternative: String,
    exon: String,
    intron: String,
    hgvsc: String,
    consequence: String,
    clinvar_sig: String,
    clinvar_review: String,
    splice_ai: Vec<String>,
}

fn odd(file: &str, transcript: &str, what: impl std::fmt::Display) -> Box<dyn std::error::Error> {
    invalid(format!("{file}: transcript {transcript}: {what}; refused"))
}

fn text(dict: &TextDict, field: &str, file: &str, transcript: &str) -> Result<String> {
    match dict.get(field) {
        Some(Item::Str(value)) => Ok(value.to_string()),
        Some(value) => Err(odd(
            file,
            transcript,
            format!("{field} is a {}, not str", value.type_name()),
        )),
        None => Err(odd(file, transcript, format!("a record lacks {field}"))),
    }
}

fn ai_names_valid(names: &[String]) -> bool {
    names.len() == 9
        && names.iter().collect::<HashSet<_>>().len() == 9
        && names
            .iter()
            .all(|name| SPLICE_AI_FIELDS.contains(&name.as_str()))
}

fn plain_transcript(
    file: &str,
    transcript: &str,
    value: &Item,
    names: &mut Option<Vec<String>>,
) -> Result<Vec<PlainRecord>> {
    let Item::List(records) = value else {
        return Err(odd(
            file,
            transcript,
            format!("the value is a {}, not list", value.type_name()),
        ));
    };
    records.borrow().iter().map(|item| {
        let Item::Dict(record) = item else {
            return Err(odd(file, transcript, "a splice record is not a dict"));
        };
        let record = record.borrow();
        if record.len() != RECORD_FIELDS.len() || record.entries.iter().any(|(key, _)| !RECORD_FIELDS.contains(&key.as_ref())) {
            return Err(odd(file, transcript, "a splice record does not have exactly the eleven known fields"));
        }
        let pos = match record.get("pos") {
            Some(Item::Int(number)) => *number,
            Some(other) => return Err(odd(file, transcript, format!("pos is a {}, not int", other.type_name()))),
            None => return Err(odd(file, transcript, "a record lacks pos")),
        };
        let Some(Item::Dict(ai)) = record.get("splice_ai") else {
            return Err(odd(file, transcript, "splice_ai is not a dict"));
        };
        let ai = ai.borrow();
        let mut values = Vec::with_capacity(ai.len());
        if !ai.entries.is_empty() {
            let present: Vec<String> = ai.entries.iter().map(|(name, _)| name.to_string()).collect();
            if !ai_names_valid(&present) {
                return Err(odd(file, transcript, "splice_ai has a partial or unknown key set; expected {} or all nine known keys"));
            }
            match names {
                Some(expected) if *expected != present => return Err(odd(file, transcript, "splice_ai dictionaries have inconsistent insertion order")),
                None => *names = Some(present),
                _ => {},
            }
            for (_, value) in &ai.entries {
                let Item::Str(value) = value else {
                    return Err(odd(file, transcript, format!("a splice_ai value is a {}, not str", value.type_name())));
                };
                values.push(value.to_string());
            }
        }
        Ok(PlainRecord {
            chrom: text(&record, "chrom", file, transcript)?, pos,
            reference: text(&record, "ref", file, transcript)?,
            alternative: text(&record, "alt", file, transcript)?,
            exon: text(&record, "exon", file, transcript)?,
            intron: text(&record, "intron", file, transcript)?,
            hgvsc: text(&record, "hgvsc", file, transcript)?,
            consequence: text(&record, "consequence", file, transcript)?,
            clinvar_sig: text(&record, "clinvar_sig", file, transcript)?,
            clinvar_review: text(&record, "clinvar_review", file, transcript)?, splice_ai: values,
        })
    }).collect()
}

fn encode(record: PlainRecord, pool: &mut StringPool) -> Result<SpliceRecord> {
    Ok(SpliceRecord {
        chrom: pool.intern(&record.chrom)?,
        pos: record.pos,
        reference: record.reference,
        alternative: record.alternative,
        exon: pool.intern(&record.exon)?,
        intron: pool.intern(&record.intron)?,
        hgvsc: record.hgvsc,
        consequence: pool.intern(&record.consequence)?,
        clinvar_sig: pool.intern(&record.clinvar_sig)?,
        clinvar_review: pool.intern(&record.clinvar_review)?,
        splice_ai: record
            .splice_ai
            .iter()
            .map(|value| pool.intern(value))
            .collect::<Result<_>>()?,
    })
}

/// A checked map; string lookups borrow the pool directly without copying it.
pub(crate) struct ClinvarSpliceCache {
    mapped: MappedCache<ClinvarSpliceStat>,
}

impl ClinvarSpliceCache {
    pub(crate) fn open(path: &str) -> Result<Self> {
        let mapped = MappedCache::<ClinvarSpliceStat>::open(path, SCHEMA)?;
        let root = mapped.root();
        if root.transcripts.len() != root.entries.len()
            || root
                .transcripts
                .windows(2)
                .any(|p| p[0].as_str() >= p[1].as_str())
        {
            return Err(invalid(format!(
                "{path}: transcripts must be sorted, unique and one per entry"
            )));
        }
        let names: Vec<String> = root
            .splice_ai_names
            .iter()
            .map(|s| s.as_str().to_owned())
            .collect();
        if !names.is_empty() && !ai_names_valid(&names) {
            return Err(invalid(format!("{path}: invalid splice_ai_names")));
        }
        let count = root.texts.len();
        let unique = root
            .texts
            .iter()
            .map(|s| s.as_str())
            .collect::<HashSet<_>>()
            .len();
        if count != unique {
            return Err(invalid(format!("{path}: text pool contains duplicates")));
        }
        let valid_id = |id: &rkyv::Archived<u32>| (id.to_native() as usize) < count;
        if !root
            .entries
            .iter()
            .flat_map(|records| records.iter())
            .all(|record| {
                [
                    &record.chrom,
                    &record.exon,
                    &record.intron,
                    &record.consequence,
                    &record.clinvar_sig,
                    &record.clinvar_review,
                ]
                .into_iter()
                .all(valid_id)
                    && (record.splice_ai.is_empty()
                        || record.splice_ai.len() == names.len() && !names.is_empty())
                    && record.splice_ai.iter().all(valid_id)
            })
        {
            return Err(invalid(format!(
                "{path}: a text index or splice_ai count is invalid"
            )));
        }
        Ok(Self { mapped })
    }
    pub(crate) fn len(&self) -> usize {
        self.mapped.root().transcripts.len()
    }
    pub(crate) fn get(
        &self,
        transcript: &str,
    ) -> Option<&rkyv::vec::ArchivedVec<ArchivedSpliceRecord>> {
        let root = self.mapped.root();
        find_sorted(&root.transcripts, transcript).map(|at| &root.entries[at])
    }
    pub(crate) fn text(&self, id: &rkyv::Archived<u32>) -> &str {
        self.mapped.root().texts[id.to_native() as usize].as_str()
    }
    pub(crate) fn splice_ai_names(&self) -> &[rkyv::string::ArchivedString] {
        &self.mapped.root().splice_ai_names
    }
    fn plain(&self, transcript: &str) -> Option<Vec<PlainRecord>> {
        self.get(transcript).map(|records| {
            records
                .iter()
                .map(|r| PlainRecord {
                    chrom: self.text(&r.chrom).to_owned(),
                    pos: r.pos.to_native(),
                    reference: r.reference.as_str().to_owned(),
                    alternative: r.alternative.as_str().to_owned(),
                    exon: self.text(&r.exon).to_owned(),
                    intron: self.text(&r.intron).to_owned(),
                    hgvsc: r.hgvsc.as_str().to_owned(),
                    consequence: self.text(&r.consequence).to_owned(),
                    clinvar_sig: self.text(&r.clinvar_sig).to_owned(),
                    clinvar_review: self.text(&r.clinvar_review).to_owned(),
                    splice_ai: r
                        .splice_ai
                        .iter()
                        .map(|id| self.text(id).to_owned())
                        .collect(),
                })
                .collect()
        })
    }
}

fn exact(records: Vec<PlainRecord>, names: &[String]) -> ExactValue {
    let text = |value: String| ExactValue::Text(value);
    ExactValue::List(
        records
            .into_iter()
            .map(|r| {
                ExactValue::Dict(vec![
                    (text("chrom".into()), text(r.chrom)),
                    (text("pos".into()), ExactValue::Int(r.pos)),
                    (text("ref".into()), text(r.reference)),
                    (text("alt".into()), text(r.alternative)),
                    (text("exon".into()), text(r.exon)),
                    (text("intron".into()), text(r.intron)),
                    (text("hgvsc".into()), text(r.hgvsc)),
                    (text("consequence".into()), text(r.consequence)),
                    (text("clinvar_sig".into()), text(r.clinvar_sig)),
                    (text("clinvar_review".into()), text(r.clinvar_review)),
                    (
                        text("splice_ai".into()),
                        ExactValue::Dict(
                            names
                                .iter()
                                .cloned()
                                .zip(r.splice_ai)
                                .map(|(k, v)| (text(k), text(v)))
                                .collect(),
                        ),
                    ),
                ])
            })
            .collect(),
    )
}

#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn export(path: &str, out: &mut dyn Write) -> Result<usize> {
    let cache = ClinvarSpliceCache::open(path)?;
    let names: Vec<String> = cache
        .splice_ai_names()
        .iter()
        .map(|s| s.as_str().to_owned())
        .collect();
    for transcript in cache.mapped.root().transcripts.iter() {
        let key = transcript.as_str();
        write_entry(
            out,
            &ExactValue::Text(key.to_owned()),
            &exact(cache.plain(key).unwrap(), &names),
        )?;
    }
    Ok(cache.len())
}

pub(crate) fn convert(pickle: &str, output: &Path) -> Result<String> {
    ensure_distinct_paths(Path::new(pickle), output)?;
    let started = Instant::now();
    let mut pool = StringPool::default();
    let mut names = None;
    let mut entries = Vec::new();
    let mut records = 0_usize;
    let distinct = visit_top_level_dict(
        pickle,
        |_| true,
        |key, value| {
            let plain = plain_transcript(pickle, &key, &value, &mut names)?;
            records += plain.len();
            entries.push((
                key.to_string(),
                plain
                    .into_iter()
                    .map(|r| encode(r, &mut pool))
                    .collect::<Result<Vec<_>>>()?,
            ));
            Ok(())
        },
    )?;
    if distinct != entries.len() {
        return Err(invalid(format!(
            "{pickle}: a transcript is set twice; refused"
        )));
    }
    entries.sort_by(|a, b| a.0.cmp(&b.0));
    let (transcripts, entries) = entries.into_iter().unzip();
    let value = ClinvarSpliceStat {
        provenance: provenance(pickle)?,
        texts: pool.into_texts(),
        splice_ai_names: names.unwrap_or_default(),
        transcripts,
        entries,
    };
    let count = value.transcripts.len();
    let text_count = value.texts.len();
    let bytes = rkyv::to_bytes::<rancor::Error>(&value)?;
    write_cache_checked(output, SCHEMA, &bytes, |staged| {
        let cache = ClinvarSpliceCache::open(&staged.to_string_lossy())?;
        if cache.len() != count {
            return Err(invalid("converted transcript count differs"));
        }
        let mut names = if cache.splice_ai_names().is_empty() {
            None
        } else {
            Some(
                cache
                    .splice_ai_names()
                    .iter()
                    .map(|s| s.as_str().to_owned())
                    .collect(),
            )
        };
        let mut compared = 0_usize;
        let distinct = visit_top_level_dict(
            pickle,
            |_| true,
            |key, input| {
                let expected = plain_transcript(pickle, &key, &input, &mut names)?;
                if cache.plain(&key).as_ref() != Some(&expected) {
                    return Err(invalid(format!(
                        "{pickle}: transcript {key} differs after conversion"
                    )));
                }
                compared += 1;
                Ok(())
            },
        )?;
        if distinct != count || compared != count {
            return Err(invalid("the pickle changed while converted"));
        }
        Ok(())
    })?;
    Ok(format!(
        "{{\"cache\":\"clinvar_splice_stat\",\"transcripts\":{count},\"records\":{records},\"texts\":{text_count},\"checked_transcripts\":{count},\"bytes\":{},\"elapsed_s\":{:.3}}}",
        std::fs::metadata(output)?.len(),
        started.elapsed().as_secs_f64()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{cell::RefCell, rc::Rc};
    fn dict(entries: Vec<(&str, Item)>) -> Item {
        let mut dict = TextDict::default();
        dict.entries = entries
            .into_iter()
            .map(|(key, value)| (Rc::from(key), value))
            .collect();
        Item::Dict(Rc::new(RefCell::new(dict)))
    }
    fn input_record(ai: Item) -> Item {
        dict(
            RECORD_FIELDS
                .iter()
                .map(|field| {
                    (
                        *field,
                        match *field {
                            "pos" => Item::Int(9),
                            "splice_ai" => ai.clone(),
                            _ => Item::Str(Rc::from("")),
                        },
                    )
                })
                .collect(),
        )
    }
    fn ai_input(names: &[&str]) -> Item {
        dict(
            names
                .iter()
                .map(|name| (*name, Item::Str(Rc::from(""))))
                .collect(),
        )
    }
    fn input_list(record: Item) -> Item {
        Item::List(Rc::new(RefCell::new(vec![record])))
    }
    fn path(name: &str) -> String {
        std::env::temp_dir()
            .join(format!("priva_splice_{}_{name}", std::process::id()))
            .to_str()
            .unwrap()
            .to_owned()
    }
    fn record() -> PlainRecord {
        PlainRecord {
            chrom: "chr1".into(),
            pos: 9,
            reference: "A".into(),
            alternative: "T".into(),
            exon: "2/3".into(),
            intron: "".into(),
            hgvsc: "c.9A>T".into(),
            consequence: "missense_variant".into(),
            clinvar_sig: "Pathogenic".into(),
            clinvar_review: "reviewed_by_expert_panel".into(),
            splice_ai: vec!["".into(); 9],
        }
    }
    fn root() -> ClinvarSpliceStat {
        let mut pool = StringPool::default();
        let a = encode(record(), &mut pool).unwrap();
        let b = encode(record(), &mut pool).unwrap();
        ClinvarSpliceStat {
            provenance: vec![],
            texts: pool.into_texts(),
            splice_ai_names: SPLICE_AI_FIELDS.iter().map(|s| (*s).to_owned()).collect(),
            transcripts: vec!["a".into(), "z".into()],
            entries: vec![vec![a, b], vec![]],
        }
    }
    fn save(name: &str, value: &ClinvarSpliceStat) -> String {
        let path = path(name);
        crate::rust_cache::write_cache(
            Path::new(&path),
            SCHEMA,
            &rkyv::to_bytes::<rancor::Error>(value).unwrap(),
        )
        .unwrap();
        path
    }
    #[test]
    fn typed_lookup_preserves_list_multiplicity_and_pool() {
        let path = save("good", &root());
        let cache = ClinvarSpliceCache::open(&path).unwrap();
        assert_eq!(cache.len(), 2);
        assert!(cache.get("missing").is_none());
        assert_eq!(cache.plain("a").unwrap(), vec![record(), record()]);
        assert_eq!(cache.get("z").unwrap().len(), 0);
        let mut out = Vec::new();
        assert_eq!(export(&path, &mut out).unwrap(), 2);
        assert!(
            std::str::from_utf8(&out)
                .unwrap()
                .contains("SpliceAI_pred_SYMBOL")
        );
        std::fs::remove_file(path).unwrap();
    }
    #[test]
    fn refuses_wrong_shapes_missing_fields_and_inconsistent_ai_order() {
        let value = input_list(input_record(ai_input(&SPLICE_AI_FIELDS)));
        let mut names = None;
        assert_eq!(
            plain_transcript("x", "t", &value, &mut names)
                .unwrap()
                .len(),
            1
        );
        let mut reverse = SPLICE_AI_FIELDS;
        reverse.reverse();
        assert!(
            plain_transcript(
                "x",
                "t",
                &input_list(input_record(ai_input(&reverse))),
                &mut names
            )
            .is_err()
        );
        assert!(
            plain_transcript(
                "x",
                "t",
                &input_list(input_record(ai_input(&SPLICE_AI_FIELDS[..8]))),
                &mut None
            )
            .is_err()
        );
        assert!(
            plain_transcript(
                "x",
                "t",
                &input_list(input_record(dict(vec![("other", Item::Str(Rc::from("")))]))),
                &mut None
            )
            .is_err()
        );
        assert!(plain_transcript("x", "t", &Item::Tuple(Rc::new(vec![])), &mut None).is_err());
        for change in [0, 1, 2, 3] {
            let record = input_record(dict(vec![]));
            let Item::Dict(dict) = &record else {
                unreachable!()
            };
            match change {
                0 => dict.borrow_mut().entries[1].1 = Item::Bool(true),
                1 => dict.borrow_mut().entries[0].1 = Item::None,
                2 => {
                    dict.borrow_mut().entries.pop();
                }
                _ => dict
                    .borrow_mut()
                    .entries
                    .push((Rc::from("unknown"), Item::None)),
            }
            assert!(plain_transcript("x", "t", &input_list(record), &mut None).is_err());
        }
        let empty =
            plain_transcript("x", "t", &input_list(input_record(dict(vec![]))), &mut None).unwrap();
        assert!(empty[0].splice_ai.is_empty());
    }
    #[test]
    fn refuses_invalid_semantic_indices_names_and_counts() {
        for (name, alter) in [
            ("index", 0),
            ("count", 1),
            ("keys", 2),
            ("names", 3),
            ("pool", 4),
            ("order", 5),
        ] {
            let mut value = root();
            match alter {
                0 => value.entries[0][0].chrom = u32::MAX,
                1 => {
                    value.entries[0][0].splice_ai.pop();
                }
                2 => {
                    value.entries.pop();
                }
                3 => value.splice_ai_names[0] = "unknown".into(),
                4 => value.texts.push(value.texts[0].clone()),
                _ => value.transcripts.swap(0, 1),
            }
            let path = save(name, &value);
            assert!(ClinvarSpliceCache::open(&path).is_err(), "{name}");
            std::fs::remove_file(path).unwrap();
        }
    }
    #[test]
    fn refuses_wrong_schema_and_damaged_archive() {
        let path = save("corrupt", &root());
        let mut bytes = std::fs::read(&path).unwrap();
        bytes[12] = 3;
        std::fs::write(&path, &bytes).unwrap();
        assert!(ClinvarSpliceCache::open(&path).is_err());
        bytes[12] = 2;
        let end = bytes.len();
        bytes[end - 8..].fill(0xff);
        std::fs::write(&path, &bytes).unwrap();
        assert!(ClinvarSpliceCache::open(&path).is_err());
        std::fs::remove_file(path).unwrap();
    }
    #[test]
    fn conversion_sorts_top_keys_and_preserves_destination_on_failure() {
        let input = path("convert.pkl");
        let output = path("convert.rkyv");
        std::fs::write(&input, b"\x80\x04}\x8c\x01z]s\x8c\x01a]s.").unwrap();
        convert(&input, Path::new(&output)).unwrap();
        let cache = ClinvarSpliceCache::open(&output).unwrap();
        assert_eq!(cache.mapped.root().transcripts[0].as_str(), "a");
        drop(cache);
        let before = std::fs::read(&output).unwrap();
        // Tuple instead of a transcript's required list.
        std::fs::write(&input, b"\x80\x04}\x8c\x01z)s.").unwrap();
        assert!(convert(&input, Path::new(&output)).is_err());
        assert_eq!(std::fs::read(&output).unwrap(), before);
        // Duplicate top-level assignment is refused even when its values agree.
        std::fs::write(&input, b"\x80\x04}\x8c\x01z]s\x8c\x01z]s.").unwrap();
        assert!(convert(&input, Path::new(&output)).is_err());
        assert_eq!(std::fs::read(&output).unwrap(), before);
        std::fs::remove_file(input).unwrap();
        std::fs::remove_file(output).unwrap();
    }
}
