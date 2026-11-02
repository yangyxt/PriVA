//! Schema 4, version 1: `{transcript: {exon_key: (float, float, float,
//! float, (float, ...))}}`. Both Python sequences are tuples, not lists;
//! exon dictionary order, repeated samples, NaN payloads and signed zero survive.
use super::Provenance;
use super::common::{ExactValue, provenance, write_entry};
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

pub(crate) const SCHEMA: Schema = Schema {
    name: "clinvar_patho_exon_af_stat",
    id: 4,
    version: 1,
};

/// One ordered dict entry; the value is the Python five-tuple, with a tuple at item4.
#[derive(Archive, Serialize, Deserialize, Debug)]
pub(crate) struct ExonAf {
    pub(crate) exon_key: String,
    pub(crate) item0: f64,
    pub(crate) item1: f64,
    pub(crate) item2: f64,
    pub(crate) item3: f64,
    pub(crate) item4: Vec<f64>,
}

#[derive(Archive, Serialize, Deserialize, Debug)]
pub(crate) struct ClinvarExonAfStat {
    pub(crate) provenance: Vec<Provenance>,
    pub(crate) transcripts: Vec<String>,
    pub(crate) entries: Vec<Vec<ExonAf>>,
}

fn odd(file: &str, transcript: &str, what: impl std::fmt::Display) -> Box<dyn std::error::Error> {
    invalid(format!("{file}: transcript {transcript}: {what}; refused"))
}

fn plain_transcript(file: &str, transcript: &str, input: &Item) -> Result<Vec<ExonAf>> {
    let Item::Dict(exons) = input else {
        return Err(odd(
            file,
            transcript,
            format!("the value is a {}, not dict", input.type_name()),
        ));
    };
    exons
        .borrow()
        .entries
        .iter()
        .map(|(key, value)| {
            let Item::Tuple(items) = value else {
                return Err(odd(
                    file,
                    transcript,
                    format!("exon {key} holds a {}, not a five-tuple", value.type_name()),
                ));
            };
            if items.len() != 5 {
                return Err(odd(
                    file,
                    transcript,
                    format!("exon {key} has {} items, expected five", items.len()),
                ));
            }
            let float = |item: &Item, at: &str| -> Result<f64> {
                match item {
                    Item::Float(value) => Ok(*value),
                    other => Err(odd(
                        file,
                        transcript,
                        format!("exon {key} item {at} is a {}, not float", other.type_name()),
                    )),
                }
            };
            let Item::Tuple(samples) = &items[4] else {
                return Err(odd(
                    file,
                    transcript,
                    format!("exon {key} item 4 is a {}, not tuple", items[4].type_name()),
                ));
            };
            Ok(ExonAf {
                exon_key: key.to_string(),
                item0: float(&items[0], "0")?,
                item1: float(&items[1], "1")?,
                item2: float(&items[2], "2")?,
                item3: float(&items[3], "3")?,
                item4: samples
                    .iter()
                    .map(|item| float(item, "4 tuple member"))
                    .collect::<Result<_>>()?,
            })
        })
        .collect()
}

fn exact_value(
    key: String,
    first: [f64; 4],
    samples: impl Iterator<Item = f64>,
) -> (ExactValue, ExactValue) {
    let mut items: Vec<ExactValue> = first
        .into_iter()
        .map(|f| ExactValue::Float(f.to_bits()))
        .collect();
    items.push(ExactValue::Tuple(
        samples.map(|f| ExactValue::Float(f.to_bits())).collect(),
    ));
    (ExactValue::Text(key), ExactValue::Tuple(items))
}

fn exact_plain(entries: &[ExonAf]) -> ExactValue {
    ExactValue::Dict(
        entries
            .iter()
            .map(|e| {
                exact_value(
                    e.exon_key.clone(),
                    [e.item0, e.item1, e.item2, e.item3],
                    e.item4.iter().copied(),
                )
            })
            .collect(),
    )
}

/// Zero-copy ordered exon entries, validated once when mapped.
pub(crate) struct ClinvarExonAfCache {
    mapped: MappedCache<ClinvarExonAfStat>,
}
impl ClinvarExonAfCache {
    pub(crate) fn open(path: &str) -> Result<Self> {
        let mapped = MappedCache::<ClinvarExonAfStat>::open(path, SCHEMA)?;
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
        for exons in root.entries.iter() {
            if exons
                .iter()
                .map(|e| e.exon_key.as_str())
                .collect::<HashSet<_>>()
                .len()
                != exons.len()
            {
                return Err(invalid(format!("{path}: a transcript repeats an exon key")));
            }
        }
        Ok(Self { mapped })
    }
    pub(crate) fn len(&self) -> usize {
        self.mapped.root().transcripts.len()
    }
    pub(crate) fn keys(&self) -> &[rkyv::string::ArchivedString] {
        &self.mapped.root().transcripts
    }
    pub(crate) fn get(&self, transcript: &str) -> Option<&rkyv::vec::ArchivedVec<ArchivedExonAf>> {
        let root = self.mapped.root();
        find_sorted(&root.transcripts, transcript).map(|at| &root.entries[at])
    }
    fn exact(&self, transcript: &str) -> Option<ExactValue> {
        self.get(transcript).map(|entries| {
            ExactValue::Dict(
                entries
                    .iter()
                    .map(|e| {
                        exact_value(
                            e.exon_key.as_str().to_owned(),
                            [
                                e.item0.to_native(),
                                e.item1.to_native(),
                                e.item2.to_native(),
                                e.item3.to_native(),
                            ],
                            e.item4.iter().map(|v| v.to_native()),
                        )
                    })
                    .collect(),
            )
        })
    }
}

#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn export(path: &str, out: &mut dyn Write) -> Result<usize> {
    let cache = ClinvarExonAfCache::open(path)?;
    for key in cache.mapped.root().transcripts.iter() {
        write_entry(
            out,
            &ExactValue::Text(key.as_str().to_owned()),
            &cache.exact(key.as_str()).unwrap(),
        )?;
    }
    Ok(cache.len())
}

pub(crate) fn convert(pickle: &str, output: &Path) -> Result<String> {
    ensure_distinct_paths(Path::new(pickle), output)?;
    let started = Instant::now();
    let mut entries = Vec::new();
    let mut exons = 0_usize;
    let mut samples = 0_usize;
    let distinct = visit_top_level_dict(
        pickle,
        |_| true,
        |key, input| {
            let value = plain_transcript(pickle, &key, &input)?;
            exons += value.len();
            samples += value.iter().map(|e| e.item4.len()).sum::<usize>();
            entries.push((key.to_string(), value));
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
    let value = ClinvarExonAfStat {
        provenance: provenance(pickle)?,
        transcripts,
        entries,
    };
    let count = value.transcripts.len();
    let bytes = rkyv::to_bytes::<rancor::Error>(&value)?;
    write_cache_checked(output, SCHEMA, &bytes, |staged| {
        let cache = ClinvarExonAfCache::open(&staged.to_string_lossy())?;
        if cache.len() != count {
            return Err(invalid("converted transcript count differs"));
        }
        let mut compared = 0_usize;
        let distinct = visit_top_level_dict(
            pickle,
            |_| true,
            |key, input| {
                let expected = exact_plain(&plain_transcript(pickle, &key, &input)?);
                if cache.exact(&key).as_ref() != Some(&expected) {
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
        "{{\"cache\":\"clinvar_patho_exon_af_stat\",\"transcripts\":{count},\"exons\":{exons},\"samples\":{samples},\"checked_transcripts\":{count},\"bytes\":{},\"elapsed_s\":{:.3}}}",
        std::fs::metadata(output)?.len(),
        started.elapsed().as_secs_f64()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::acmg_step::pickle_dict_stream::TextDict;
    use std::{cell::RefCell, rc::Rc};
    fn dict(entries: Vec<(&str, Item)>) -> Item {
        let mut d = TextDict::default();
        d.entries = entries.into_iter().map(|(k, v)| (Rc::from(k), v)).collect();
        Item::Dict(Rc::new(RefCell::new(d)))
    }
    fn tuple(items: Vec<Item>) -> Item {
        Item::Tuple(Rc::new(items))
    }
    fn tuple_value() -> Item {
        tuple(vec![
            Item::Float(-0.0),
            Item::Float(f64::from_bits(0x7ff8_0000_0000_0042)),
            Item::Float(f64::INFINITY),
            Item::Float(f64::NEG_INFINITY),
            tuple(vec![Item::Float(1.0), Item::Float(1.0), Item::Float(-0.0)]),
        ])
    }
    fn path(name: &str) -> String {
        std::env::temp_dir()
            .join(format!("priva_exon_{}_{name}", std::process::id()))
            .to_str()
            .unwrap()
            .to_owned()
    }
    fn root() -> ClinvarExonAfStat {
        ClinvarExonAfStat {
            provenance: vec![],
            transcripts: vec!["a".into(), "z".into()],
            entries: vec![
                plain_transcript(
                    "x",
                    "a",
                    &dict(vec![("2/3", tuple_value()), ("1/3", tuple_value())]),
                )
                .unwrap(),
                vec![],
            ],
        }
    }
    fn save(name: &str, value: &ClinvarExonAfStat) -> String {
        let p = path(name);
        crate::rust_cache::write_cache(
            Path::new(&p),
            SCHEMA,
            &rkyv::to_bytes::<rancor::Error>(value).unwrap(),
        )
        .unwrap();
        p
    }
    #[test]
    fn preserves_all_float_bits_order_and_tuple_shapes() {
        let value = root();
        let p = save("float", &value);
        let c = ClinvarExonAfCache::open(&p).unwrap();
        assert_eq!(c.len(), 2);
        assert_eq!(c.exact("a"), Some(exact_plain(&value.entries[0])));
        assert!(c.get("missing").is_none());
        let exons = c.get("a").unwrap();
        assert_eq!(exons[0].exon_key.as_str(), "2/3");
        assert_eq!(exons[1].exon_key.as_str(), "1/3");
        assert_eq!(exons[0].item0.to_native().to_bits(), (-0.0f64).to_bits());
        assert_eq!(exons[0].item1.to_native().to_bits(), 0x7ff8_0000_0000_0042);
        let mut output = Vec::new();
        assert_eq!(export(&p, &mut output).unwrap(), 2);
        std::fs::remove_file(p).unwrap();
    }
    #[test]
    fn refuses_lists_ints_bools_none_wrong_arity_and_nonfloat_samples() {
        let list = Item::List(Rc::new(RefCell::new(vec![])));
        for value in [
            list.clone(),
            Item::None,
            Item::Bool(false),
            tuple(vec![]),
            tuple(vec![
                Item::Int(0),
                Item::Float(0.0),
                Item::Float(0.0),
                Item::Float(0.0),
                tuple(vec![]),
            ]),
            tuple(vec![
                Item::Float(0.0),
                Item::Float(0.0),
                Item::Float(0.0),
                Item::Float(0.0),
                list,
            ]),
            tuple(vec![
                Item::Float(0.0),
                Item::Float(0.0),
                Item::Float(0.0),
                Item::Float(0.0),
                tuple(vec![Item::Bool(true)]),
            ]),
        ] {
            assert!(plain_transcript("x", "t", &dict(vec![("1/1", value)])).is_err());
        }
    }
    #[test]
    fn refuses_corrupt_schema_duplicate_exons_and_top_level_counts() {
        for kind in 0..4 {
            let mut value = root();
            match kind {
                0 => value.transcripts.swap(0, 1),
                1 => {
                    value.entries.pop();
                }
                2 => value.entries[0][1].exon_key = "2/3".into(),
                _ => value.transcripts[1] = "a".into(),
            };
            let p = save(&format!("bad{kind}"), &value);
            assert!(ClinvarExonAfCache::open(&p).is_err());
            std::fs::remove_file(p).unwrap();
        }
        let p = save("schema", &root());
        let mut bytes = std::fs::read(&p).unwrap();
        bytes[12] = 6;
        std::fs::write(&p, &bytes).unwrap();
        assert!(ClinvarExonAfCache::open(&p).is_err());
        bytes[12] = 4;
        let end = bytes.len();
        bytes[end - 8..].fill(0xff);
        std::fs::write(&p, &bytes).unwrap();
        assert!(ClinvarExonAfCache::open(&p).is_err());
        std::fs::remove_file(p).unwrap();
    }
    #[test]
    fn converts_pickle_float_bits_and_keeps_failed_destination() {
        // Python protocol 4: two transcripts, one exon five-tuple with -0, a
        // NaN payload, infinities and duplicate sample values in the inner tuple.
        const PICKLE: &[u8] = b"\x80\x04\x95\\\x00\x00\x00\x00\x00\x00\x00}\x94(\x8c\x01z\x94}\x94\x8c\x032/3\x94(G\x80\x00\x00\x00\x00\x00\x00\x00G\x7f\xf8\x00\x00\x00\x00\x00BG\x7f\xf0\x00\x00\x00\x00\x00\x00G\xff\xf0\x00\x00\x00\x00\x00\x00G?\xf0\x00\x00\x00\x00\x00\x00G?\xf0\x00\x00\x00\x00\x00\x00G\x80\x00\x00\x00\x00\x00\x00\x00\x87\x94t\x94s\x8c\x01a\x94}\x94u.";
        let input = path("convert.pkl");
        let output = path("convert.rkyv");
        std::fs::write(&input, PICKLE).unwrap();
        convert(&input, Path::new(&output)).unwrap();
        let c = ClinvarExonAfCache::open(&output).unwrap();
        assert_eq!(c.mapped.root().transcripts[0].as_str(), "a");
        assert_eq!(
            c.exact("z"),
            Some(exact_plain(
                &plain_transcript("x", "z", &dict(vec![("2/3", tuple_value())])).unwrap()
            ))
        );
        drop(c);
        let before = std::fs::read(&output).unwrap();
        std::fs::write(&input, b"\x80\x04}\x8c\x01z]s.").unwrap();
        assert!(convert(&input, Path::new(&output)).is_err());
        assert_eq!(std::fs::read(&output).unwrap(), before);
        std::fs::remove_file(input).unwrap();
        std::fs::remove_file(output).unwrap();
    }
}
