//! Schema 6, version 1: `{gene: None | {known_field: str | int | float}}`.
//! The twelve known builder field names may be absent; present entries retain
//! their dictionary insertion order and exact scalar kind. `None`, an empty
//! dictionary and a missing top-level gene remain distinct.
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
    name: "clinvar_patho_af_stat",
    id: 6,
    version: 1,
};
const FIELDS: [&str; 12] = [
    "chrom",
    "pos",
    "ref",
    "alt",
    "af",
    "hgvsp",
    "impact",
    "consequence",
    "clnsig",
    "clnrevstat",
    "gene_symbol",
    "acc_patho_af",
];

#[derive(Archive, Serialize, Deserialize, Debug, Clone)]
pub(crate) enum GeneAfValue {
    Text(String),
    Int(i64),
    Float(f64),
}
#[derive(Archive, Serialize, Deserialize, Debug, Clone)]
pub(crate) struct GeneAfEntry {
    pub(crate) key: String,
    pub(crate) value: GeneAfValue,
}
#[derive(Archive, Serialize, Deserialize, Debug)]
pub(crate) struct ClinvarGeneAfStat {
    pub(crate) provenance: Vec<Provenance>,
    pub(crate) genes: Vec<String>,
    pub(crate) entries: Vec<Option<Vec<GeneAfEntry>>>,
}

fn plain_gene(file: &str, gene: &str, input: &Item) -> Result<Option<Vec<GeneAfEntry>>> {
    if matches!(input, Item::None) {
        return Ok(None);
    }
    let odd = |what: String| invalid(format!("{file}: gene {gene}: {what}; refused"));
    let Item::Dict(record) = input else {
        return Err(odd(format!(
            "value is a {}, not None or dict",
            input.type_name()
        )));
    };
    record
        .borrow()
        .entries
        .iter()
        .map(|(key, input)| {
            if !FIELDS.contains(&key.as_ref()) {
                return Err(odd(format!("unknown field {key}")));
            }
            let value = match input {
                Item::Str(s) => GeneAfValue::Text(s.to_string()),
                Item::Int(i) => GeneAfValue::Int(*i),
                Item::Float(f) => GeneAfValue::Float(*f),
                other => {
                    return Err(odd(format!(
                        "field {key} is a {}, not str/int/float",
                        other.type_name()
                    )));
                }
            };
            Ok(GeneAfEntry {
                key: key.to_string(),
                value,
            })
        })
        .collect::<Result<Vec<_>>>()
        .map(Some)
}

fn exact_plain(record: Option<&[GeneAfEntry]>) -> ExactValue {
    match record {
        None => ExactValue::None,
        Some(entries) => ExactValue::Dict(
            entries
                .iter()
                .map(|entry| {
                    (
                        ExactValue::Text(entry.key.clone()),
                        match &entry.value {
                            GeneAfValue::Text(s) => ExactValue::Text(s.clone()),
                            GeneAfValue::Int(i) => ExactValue::Int(*i),
                            GeneAfValue::Float(f) => ExactValue::Float(f.to_bits()),
                        },
                    )
                })
                .collect(),
        ),
    }
}

/// A checked, mapped gene dictionary. `get` returns None for a missing gene,
/// Some(None) for Python None and Some(Some(entries)) for the ordered record.
pub(crate) struct ClinvarGeneAfCache {
    mapped: MappedCache<ClinvarGeneAfStat>,
}
impl ClinvarGeneAfCache {
    pub(crate) fn open(path: &str) -> Result<Self> {
        let mapped = MappedCache::<ClinvarGeneAfStat>::open(path, SCHEMA)?;
        let root = mapped.root();
        if root.genes.len() != root.entries.len()
            || root
                .genes
                .windows(2)
                .any(|p| p[0].as_str() >= p[1].as_str())
        {
            return Err(invalid(format!(
                "{path}: genes must be sorted, unique and one per entry"
            )));
        }
        for entry in root.entries.iter() {
            if let Some(record) = entry.as_ref() {
                let names: HashSet<&str> = record.iter().map(|e| e.key.as_str()).collect();
                if names.len() != record.len() || names.iter().any(|name| !FIELDS.contains(name)) {
                    return Err(invalid(format!(
                        "{path}: a gene record has duplicate or unknown fields"
                    )));
                }
            }
        }
        Ok(Self { mapped })
    }
    pub(crate) fn len(&self) -> usize {
        self.mapped.root().genes.len()
    }
    pub(crate) fn keys(&self) -> &[rkyv::string::ArchivedString] {
        &self.mapped.root().genes
    }
    pub(crate) fn get(
        &self,
        gene: &str,
    ) -> Option<Option<&rkyv::vec::ArchivedVec<ArchivedGeneAfEntry>>> {
        let root = self.mapped.root();
        find_sorted(&root.genes, gene).map(|at| root.entries[at].as_ref())
    }
    pub(crate) fn owned_entries(&self) -> Vec<(String, Option<Vec<GeneAfEntry>>)> {
        self.keys()
            .iter()
            .map(|gene| {
                let entries = self.get(gene.as_str()).unwrap().map(|record| {
                    record
                        .iter()
                        .map(|entry| GeneAfEntry {
                            key: entry.key.as_str().to_owned(),
                            value: match &entry.value {
                                ArchivedGeneAfValue::Text(s) => {
                                    GeneAfValue::Text(s.as_str().to_owned())
                                }
                                ArchivedGeneAfValue::Int(i) => GeneAfValue::Int(i.to_native()),
                                ArchivedGeneAfValue::Float(f) => GeneAfValue::Float(f.to_native()),
                            },
                        })
                        .collect()
                });
                (gene.as_str().to_owned(), entries)
            })
            .collect()
    }
    fn exact(&self, gene: &str) -> Option<ExactValue> {
        self.get(gene).map(|record| match record {
            None => ExactValue::None,
            Some(entries) => ExactValue::Dict(
                entries
                    .iter()
                    .map(|entry| {
                        (
                            ExactValue::Text(entry.key.as_str().to_owned()),
                            match &entry.value {
                                ArchivedGeneAfValue::Text(s) => {
                                    ExactValue::Text(s.as_str().to_owned())
                                }
                                ArchivedGeneAfValue::Int(i) => ExactValue::Int(i.to_native()),
                                ArchivedGeneAfValue::Float(f) => {
                                    ExactValue::Float(f.to_native().to_bits())
                                }
                            },
                        )
                    })
                    .collect(),
            ),
        })
    }
}

/// Direct producer boundary for the existing schema. Only lookup keys are
/// sorted; None, record field order, scalar kinds and float bits are retained.
pub(crate) fn write_entries(
    output: &Path,
    mut entries: Vec<(String, Option<Vec<GeneAfEntry>>)>,
    provenance: Vec<Provenance>,
) -> Result<usize> {
    entries.sort_by(|a, b| a.0.cmp(&b.0));
    let (genes, entries) = entries.into_iter().unzip();
    let value = ClinvarGeneAfStat {
        provenance,
        genes,
        entries,
    };
    let count = value.genes.len();
    let bytes = rkyv::to_bytes::<rancor::Error>(&value)?;
    write_cache_checked(output, SCHEMA, &bytes, |staged| {
        let cache = ClinvarGeneAfCache::open(&staged.to_string_lossy())?;
        if cache.len() != count {
            return Err(invalid("native gene count differs"));
        }
        for (key, expected) in value.genes.iter().zip(&value.entries) {
            if cache.exact(key).as_ref() != Some(&exact_plain(expected.as_deref())) {
                return Err(invalid("native gene value differs after writing"));
            }
        }
        Ok(())
    })?;
    Ok(count)
}

#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn export(path: &str, out: &mut dyn Write) -> Result<usize> {
    let cache = ClinvarGeneAfCache::open(path)?;
    for key in cache.mapped.root().genes.iter() {
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
    let mut fields = 0_usize;
    let mut none = 0_usize;
    let distinct = visit_top_level_dict(
        pickle,
        |_| true,
        |key, input| {
            let record = plain_gene(pickle, &key, &input)?;
            match &record {
                None => none += 1,
                Some(record) => fields += record.len(),
            };
            entries.push((key.to_string(), record));
            Ok(())
        },
    )?;
    if distinct != entries.len() {
        return Err(invalid(format!("{pickle}: a gene is set twice; refused")));
    }
    entries.sort_by(|a, b| a.0.cmp(&b.0));
    let (genes, entries) = entries.into_iter().unzip();
    let value = ClinvarGeneAfStat {
        provenance: provenance(pickle)?,
        genes,
        entries,
    };
    let count = value.genes.len();
    let bytes = rkyv::to_bytes::<rancor::Error>(&value)?;
    write_cache_checked(output, SCHEMA, &bytes, |staged| {
        let cache = ClinvarGeneAfCache::open(&staged.to_string_lossy())?;
        if cache.len() != count {
            return Err(invalid("converted gene count differs"));
        }
        let mut compared = 0_usize;
        let distinct = visit_top_level_dict(
            pickle,
            |_| true,
            |key, input| {
                let expected = plain_gene(pickle, &key, &input)?;
                if cache.exact(&key).as_ref() != Some(&exact_plain(expected.as_deref())) {
                    return Err(invalid(format!(
                        "{pickle}: gene {key} differs after conversion"
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
        "{{\"cache\":\"clinvar_patho_af_stat\",\"genes\":{count},\"fields\":{fields},\"none\":{none},\"checked_genes\":{count},\"bytes\":{},\"elapsed_s\":{:.3}}}",
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
    fn record() -> Vec<GeneAfEntry> {
        vec![
            GeneAfEntry {
                key: "af".into(),
                value: GeneAfValue::Int(0),
            },
            GeneAfEntry {
                key: "acc_patho_af".into(),
                value: GeneAfValue::Float(-0.0),
            },
            GeneAfEntry {
                key: "pos".into(),
                value: GeneAfValue::Float(f64::from_bits(0x7ff8_0000_0000_1234)),
            },
            GeneAfEntry {
                key: "chrom".into(),
                value: GeneAfValue::Text("chr1".into()),
            },
        ]
    }
    fn root() -> ClinvarGeneAfStat {
        ClinvarGeneAfStat {
            provenance: vec![],
            genes: vec!["a".into(), "empty".into(), "none".into()],
            entries: vec![Some(record()), Some(vec![]), None],
        }
    }
    fn path(name: &str) -> String {
        std::env::temp_dir()
            .join(format!("priva_geneaf_{}_{name}", std::process::id()))
            .to_str()
            .unwrap()
            .to_owned()
    }
    fn save(name: &str, value: &ClinvarGeneAfStat) -> String {
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
    fn distinguishes_missing_none_empty_and_preserves_types_bits_order() {
        let p = save("types", &root());
        let c = ClinvarGeneAfCache::open(&p).unwrap();
        assert_eq!(c.len(), 3);
        assert!(c.get("absent").is_none());
        assert!(matches!(c.get("none"), Some(None)));
        assert_eq!(c.get("empty").unwrap().unwrap().len(), 0);
        assert_eq!(c.exact("a"), Some(exact_plain(Some(&record()))));
        let mut out = Vec::new();
        assert_eq!(export(&p, &mut out).unwrap(), 3);
        std::fs::remove_file(p).unwrap();
    }
    #[test]
    fn refuses_unknown_fields_containers_bool_and_inner_none() {
        for value in [
            Item::Bool(true),
            Item::None,
            Item::Tuple(Rc::new(vec![])),
            Item::List(Rc::new(RefCell::new(vec![]))),
            dict(vec![]),
        ] {
            assert!(plain_gene("x", "g", &dict(vec![("af", value)])).is_err());
        }
        assert!(plain_gene("x", "g", &dict(vec![("unknown", Item::Int(1))])).is_err());
        assert!(plain_gene("x", "g", &Item::Bool(false)).is_err());
        let expected = plain_gene(
            "x",
            "g",
            &dict(vec![
                ("af", Item::Int(0)),
                ("chrom", Item::Str(Rc::from("chr1"))),
                ("pos", Item::Float(-0.0)),
            ]),
        )
        .unwrap();
        assert_eq!(expected.as_ref().unwrap().len(), 3);
        assert_eq!(expected.as_ref().unwrap()[1].key, "chrom");
    }
    #[test]
    fn rejects_invalid_fields_key_counts_and_corruption() {
        for kind in 0..5 {
            let mut value = root();
            match kind {
                0 => value.genes.swap(0, 1),
                1 => value.genes[1] = "a".into(),
                2 => {
                    value.entries.pop();
                }
                3 => value.entries[0].as_mut().unwrap()[1].key = "af".into(),
                _ => value.entries[0].as_mut().unwrap()[0].key = "extra".into(),
            };
            let p = save(&format!("bad{kind}"), &value);
            assert!(ClinvarGeneAfCache::open(&p).is_err());
            std::fs::remove_file(p).unwrap();
        }
        let p = save("schema", &root());
        let mut bytes = std::fs::read(&p).unwrap();
        bytes[12] = 4;
        std::fs::write(&p, &bytes).unwrap();
        assert!(ClinvarGeneAfCache::open(&p).is_err());
        bytes[12] = 6;
        let end = bytes.len();
        bytes[end - 8..].fill(0xff);
        std::fs::write(&p, &bytes).unwrap();
        assert!(ClinvarGeneAfCache::open(&p).is_err());
        std::fs::remove_file(p).unwrap();
    }
    #[test]
    fn converts_pickle_without_scalar_coercion_and_preserves_failed_destination() {
        const PICKLE:&[u8]=b"\x80\x04\x95\\\x00\x00\x00\x00\x00\x00\x00}\x94(\x8c\x01z\x94}\x94(\x8c\x02af\x94K\x00\x8c\x0cacc_patho_af\x94G\x80\x00\x00\x00\x00\x00\x00\x00\x8c\x03pos\x94G\x7f\xf8\x00\x00\x00\x00\x00B\x8c\x05chrom\x94\x8c\x04chr1\x94u\x8c\x04none\x94N\x8c\x05empty\x94}\x94u.";
        let input = path("convert.pkl");
        let output = path("convert.rkyv");
        std::fs::write(&input, PICKLE).unwrap();
        convert(&input, Path::new(&output)).unwrap();
        let c = ClinvarGeneAfCache::open(&output).unwrap();
        assert_eq!(c.mapped.root().genes[0].as_str(), "empty");
        let expected = vec![
            GeneAfEntry {
                key: "af".into(),
                value: GeneAfValue::Int(0),
            },
            GeneAfEntry {
                key: "acc_patho_af".into(),
                value: GeneAfValue::Float(-0.0),
            },
            GeneAfEntry {
                key: "pos".into(),
                value: GeneAfValue::Float(f64::from_bits(0x7ff8_0000_0000_0042)),
            },
            GeneAfEntry {
                key: "chrom".into(),
                value: GeneAfValue::Text("chr1".into()),
            },
        ];
        assert_eq!(c.exact("z"), Some(exact_plain(Some(&expected))));
        assert!(matches!(c.get("none"), Some(None)));
        drop(c);
        let before = std::fs::read(&output).unwrap();
        std::fs::write(&input, b"\x80\x04}\x8c\x01z\x88s.").unwrap();
        assert!(convert(&input, Path::new(&output)).is_err());
        assert_eq!(std::fs::read(&output).unwrap(), before);
        std::fs::remove_file(input).unwrap();
        std::fs::remove_file(output).unwrap();
    }
}
