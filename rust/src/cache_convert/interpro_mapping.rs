//! Schema 9, version 1: the complete InterPro mapping written by
//! `protein_domain_mapping.parse_interpro_xml_for_mapping`.
//!
//! Python keys are exact `(database: str, accession: str)` tuples; values are lists of exact
//! six-tuples `(ipr_id, kind, short_name, name, go_terms, abstract_text)`. The five text fields
//! and all members of the GO three-tuples are pooled. `go_terms` is a list, each member is a
//! tuple of three texts; tuple/list substitutions are refused. Entries, repeated entries and
//! GO terms retain their order and multiplicity. Only the top-level keys are sorted.
use super::Provenance;
use super::common::{ExactValue, StringPool, provenance, write_entry};
use crate::cli::{Result, invalid};
use crate::pickle_reader::{self, Value};
use crate::rust_cache::{MappedCache, Schema, ensure_distinct_paths, write_cache_checked};
use rkyv::{Archive, Deserialize, Serialize, rancor};
use std::io::Write;
use std::path::Path;
use std::time::Instant;

pub(crate) const SCHEMA: Schema = Schema {
    name: "interpro_mapping_pickle",
    id: 9,
    version: 1,
};

/// A Python `(str, str)` tuple, sorted lexicographically by the two UTF-8 texts.
#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct InterproKey {
    pub(crate) database: String,
    pub(crate) accession: String,
}

/// A Python `(go_id: str, category: str, description: str)` tuple; text-pool indexes.
#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct GoTerm {
    pub(crate) id: u32,
    pub(crate) category: u32,
    pub(crate) description: u32,
}

/// The producer's six-tuple; the five scalar fields are texts and `go_terms` is a list.
#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct InterproEntry {
    pub(crate) id: u32,
    pub(crate) kind: u32,
    pub(crate) short_name: u32,
    pub(crate) name: u32,
    pub(crate) go_terms: Vec<GoTerm>,
    pub(crate) abstract_text: u32,
}

#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct InterproMapping {
    pub(crate) provenance: Vec<Provenance>,
    pub(crate) texts: Vec<String>,
    /// Strictly sorted and unique; `entries[i]` is the Python list at `keys[i]`.
    pub(crate) keys: Vec<InterproKey>,
    pub(crate) entries: Vec<Vec<InterproEntry>>,
}

#[derive(Debug, PartialEq)]
pub(crate) struct PlainEntry {
    pub(crate) id: String,
    pub(crate) kind: String,
    pub(crate) short_name: String,
    pub(crate) name: String,
    pub(crate) go_terms: Vec<(String, String, String)>,
    pub(crate) abstract_text: String,
}

pub(crate) type PlainMapping = Vec<((String, String), Vec<PlainEntry>)>;

fn text(path: &str, value: &Value) -> Result<String> {
    match value {
        Value::Str(text) => Ok(text.to_string()),
        other => Err(invalid(format!(
            "{path}: an InterPro text is a {}; refused",
            other.type_name()
        ))),
    }
}

fn plain_key(path: &str, value: &Value) -> Result<(String, String)> {
    let Value::Tuple(parts) = value else {
        return Err(invalid(format!(
            "{path}: an InterPro key is not a two-text tuple; refused"
        )));
    };
    if parts.len() != 2 {
        return Err(invalid(format!(
            "{path}: an InterPro key is not a two-text tuple; refused"
        )));
    }
    Ok((text(path, &parts[0])?, text(path, &parts[1])?))
}

fn plain_entry(path: &str, value: &Value) -> Result<PlainEntry> {
    let odd = || {
        invalid(format!(
            "{path}: an InterPro entry is not an exact six-tuple; refused"
        ))
    };
    let Value::Tuple(parts) = value else {
        return Err(odd());
    };
    if parts.len() != 6 {
        return Err(odd());
    }
    let Value::List(terms) = &parts[4] else {
        return Err(invalid(format!(
            "{path}: InterPro GO terms are not a list; refused"
        )));
    };
    let mut go_terms = Vec::with_capacity(terms.borrow().len());
    for term in terms.borrow().iter() {
        let Value::Tuple(term) = term else {
            return Err(invalid(format!(
                "{path}: an InterPro GO term is not a three-text tuple; refused"
            )));
        };
        if term.len() != 3 {
            return Err(invalid(format!(
                "{path}: an InterPro GO term is not a three-text tuple; refused"
            )));
        }
        go_terms.push((
            text(path, &term[0])?,
            text(path, &term[1])?,
            text(path, &term[2])?,
        ));
    }
    Ok(PlainEntry {
        id: text(path, &parts[0])?,
        kind: text(path, &parts[1])?,
        short_name: text(path, &parts[2])?,
        name: text(path, &parts[3])?,
        go_terms,
        abstract_text: text(path, &parts[5])?,
    })
}

fn plain_list(path: &str, value: &Value) -> Result<Vec<PlainEntry>> {
    let Value::List(list) = value else {
        return Err(invalid(format!(
            "{path}: the InterPro value is not a list; refused"
        )));
    };
    list.borrow()
        .iter()
        .map(|entry| plain_entry(path, entry))
        .collect()
}

fn plain_mapping(path: &str, value: &Value) -> Result<PlainMapping> {
    let Value::Dict(dict) = value else {
        return Err(invalid(format!(
            "{path}: the InterPro mapping is not a dict; refused"
        )));
    };
    let dict = dict.borrow();
    if dict.duplicate_assignments != 0 {
        return Err(invalid(format!(
            "{path}: an InterPro key is set twice; refused"
        )));
    }
    dict.entries
        .iter()
        .map(|(key, entries)| Ok((plain_key(path, key)?, plain_list(path, entries)?)))
        .collect()
}

fn read_pickle(path: &str) -> Result<PlainMapping> {
    let bytes = pickle_reader::read_file_by_magic(Path::new(path))?;
    plain_mapping(path, &pickle_reader::load(&bytes)?)
}

fn build(mut entries: PlainMapping, provenance: Vec<Provenance>) -> Result<Vec<u8>> {
    entries.sort_by(|left, right| left.0.cmp(&right.0));
    if entries.windows(2).any(|pair| pair[0].0 == pair[1].0) {
        return Err(invalid("the InterPro keys are not unique"));
    }
    let mut pool = StringPool::default();
    let mut keys = Vec::with_capacity(entries.len());
    let mut lists = Vec::with_capacity(entries.len());
    for ((database, accession), entries) in entries {
        keys.push(InterproKey {
            database,
            accession,
        });
        let mut list = Vec::with_capacity(entries.len());
        for entry in entries {
            list.push(InterproEntry {
                id: pool.intern(&entry.id)?,
                kind: pool.intern(&entry.kind)?,
                short_name: pool.intern(&entry.short_name)?,
                name: pool.intern(&entry.name)?,
                go_terms: entry
                    .go_terms
                    .iter()
                    .map(|(id, category, description)| {
                        Ok(GoTerm {
                            id: pool.intern(id)?,
                            category: pool.intern(category)?,
                            description: pool.intern(description)?,
                        })
                    })
                    .collect::<Result<_>>()?,
                abstract_text: pool.intern(&entry.abstract_text)?,
            });
        }
        lists.push(list);
    }
    Ok(rkyv::to_bytes::<rancor::Error>(&InterproMapping {
        provenance,
        texts: pool.into_texts(),
        keys,
        entries: lists,
    })?
    .to_vec())
}

/// Checked mmap access to the complete mapping, including fields not used by today's rules.
pub(crate) struct InterproCache {
    mapped: MappedCache<InterproMapping>,
}

impl InterproCache {
    pub(crate) fn open(path: &str) -> Result<Self> {
        let mapped = MappedCache::<InterproMapping>::open(path, SCHEMA)?;
        let root = mapped.root();
        if root.keys.len() != root.entries.len()
            || root.keys.windows(2).any(|pair| {
                (pair[0].database.as_str(), pair[0].accession.as_str())
                    >= (pair[1].database.as_str(), pair[1].accession.as_str())
            })
        {
            return Err(invalid(format!(
                "{path}: the InterPro keys are not sorted and unique, or not one per entry"
            )));
        }
        let mut distinct_texts = std::collections::HashSet::with_capacity(root.texts.len());
        if !root
            .texts
            .iter()
            .all(|text| distinct_texts.insert(text.as_str()))
        {
            return Err(invalid(format!(
                "{path}: the InterPro text pool has duplicate values"
            )));
        }
        let texts = root.texts.len();
        let valid = |id: &rkyv::Archived<u32>| (id.to_native() as usize) < texts;
        if !root.entries.iter().all(|list| {
            list.iter().all(|entry| {
                [
                    &entry.id,
                    &entry.kind,
                    &entry.short_name,
                    &entry.name,
                    &entry.abstract_text,
                ]
                .into_iter()
                .all(valid)
                    && entry.go_terms.iter().all(|term| {
                        [&term.id, &term.category, &term.description]
                            .into_iter()
                            .all(valid)
                    })
            })
        }) {
            return Err(invalid(format!(
                "{path}: an InterPro text index is out of range"
            )));
        }
        Ok(Self { mapped })
    }

    pub(crate) fn len(&self) -> usize {
        self.mapped.root().keys.len()
    }

    pub(crate) fn keys(&self) -> impl Iterator<Item = (&str, &str)> {
        self.mapped
            .root()
            .keys
            .iter()
            .map(|key| (key.database.as_str(), key.accession.as_str()))
    }

    pub(crate) fn get(
        &self,
        database: &str,
        accession: &str,
    ) -> Option<&rkyv::vec::ArchivedVec<ArchivedInterproEntry>> {
        let root = self.mapped.root();
        self.position(database, accession)
            .map(|at| &root.entries[at])
    }

    /// The stable slot for a key, for consumers that retain a small lazily decoded view.
    pub(crate) fn position(&self, database: &str, accession: &str) -> Option<usize> {
        self.mapped
            .root()
            .keys
            .binary_search_by(|key| {
                (key.database.as_str(), key.accession.as_str()).cmp(&(database, accession))
            })
            .ok()
    }

    /// Access a validated pooled string without allocating.
    pub(crate) fn text(&self, id: &rkyv::Archived<u32>) -> &str {
        self.mapped.root().texts[id.to_native() as usize].as_str()
    }

    pub(crate) fn plain(&self, database: &str, accession: &str) -> Option<Vec<PlainEntry>> {
        self.get(database, accession).map(|list| {
            list.iter()
                .map(|entry| PlainEntry {
                    id: self.text(&entry.id).to_owned(),
                    kind: self.text(&entry.kind).to_owned(),
                    short_name: self.text(&entry.short_name).to_owned(),
                    name: self.text(&entry.name).to_owned(),
                    go_terms: entry
                        .go_terms
                        .iter()
                        .map(|term| {
                            (
                                self.text(&term.id).to_owned(),
                                self.text(&term.category).to_owned(),
                                self.text(&term.description).to_owned(),
                            )
                        })
                        .collect(),
                    abstract_text: self.text(&entry.abstract_text).to_owned(),
                })
                .collect()
        })
    }
}

fn exact(entries: Vec<PlainEntry>) -> ExactValue {
    ExactValue::List(
        entries
            .into_iter()
            .map(|entry| {
                ExactValue::Tuple(vec![
                    ExactValue::Text(entry.id),
                    ExactValue::Text(entry.kind),
                    ExactValue::Text(entry.short_name),
                    ExactValue::Text(entry.name),
                    ExactValue::List(
                        entry
                            .go_terms
                            .into_iter()
                            .map(|(id, category, description)| {
                                ExactValue::Tuple(vec![
                                    ExactValue::Text(id),
                                    ExactValue::Text(category),
                                    ExactValue::Text(description),
                                ])
                            })
                            .collect(),
                    ),
                    ExactValue::Text(entry.abstract_text),
                ])
            })
            .collect(),
    )
}

pub(crate) fn export(path: &str, out: &mut dyn Write) -> Result<usize> {
    let cache = InterproCache::open(path)?;
    for (database, accession) in cache.keys() {
        write_entry(
            out,
            &ExactValue::Tuple(vec![
                ExactValue::Text(database.to_owned()),
                ExactValue::Text(accession.to_owned()),
            ]),
            &exact(cache.plain(database, accession).unwrap()),
        )?;
    }
    Ok(cache.len())
}

pub(crate) fn convert(pickle: &str, output: &Path) -> Result<String> {
    ensure_distinct_paths(Path::new(pickle), output)?;
    let start = Instant::now();
    let entries = read_pickle(pickle)?;
    let count = entries.len();
    let records: usize = entries.iter().map(|(_, list)| list.len()).sum();
    let go_terms: usize = entries
        .iter()
        .flat_map(|(_, list)| list.iter())
        .map(|entry| entry.go_terms.len())
        .sum();
    let archive = build(entries, provenance(pickle)?)?;
    let converted_s = start.elapsed().as_secs_f64();
    let check = Instant::now();
    write_cache_checked(output, SCHEMA, &archive, |staged| {
        let cache = InterproCache::open(&staged.to_string_lossy())?;
        let expected = read_pickle(pickle)?;
        if cache.len() != count || expected.len() != count {
            return Err(invalid("the InterPro mapping has a different key count"));
        }
        for ((database, accession), list) in expected {
            if cache.plain(&database, &accession) != Some(list) {
                return Err(invalid("an InterPro entry differs from the pickle"));
            }
        }
        Ok(())
    })?;
    Ok(format!(
        "{{\"cache\":\"interpro_mapping_pickle\",\"keys\":{count},\"records\":{records},\"go_terms\":{go_terms},\"bytes\":{},\"converted_s\":{converted_s:.3},\"checked_keys\":{count},\"check_s\":{:.3}}}",
        std::fs::metadata(output)?.len(),
        check.elapsed().as_secs_f64()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;
    use std::rc::Rc;

    fn text(value: &str) -> Value {
        Value::Str(value.into())
    }
    fn tuple(items: Vec<Value>) -> Value {
        Value::Tuple(Rc::new(items))
    }
    fn list(items: Vec<Value>) -> Value {
        Value::List(Rc::new(RefCell::new(items)))
    }
    fn entry() -> Value {
        tuple(vec![
            text("IPR:1"),
            text("Domain"),
            text("short"),
            text("name"),
            list(vec![tuple(vec![
                text("GO:1"),
                text("function"),
                text("binding"),
            ])]),
            text("complete abstract \u{00e9}\ntext"),
        ])
    }
    fn temporary(name: &str) -> String {
        std::env::temp_dir()
            .join(format!(
                "priva_interpro_cache_{}_{name}",
                std::process::id()
            ))
            .to_str()
            .unwrap()
            .to_owned()
    }

    #[test]
    fn retains_all_six_fields_and_list_multiplicity() {
        let plain = plain_list("fixture", &list(vec![entry(), entry()])).unwrap();
        assert_eq!(plain.len(), 2);
        let bytes = build(vec![(("PFAM".into(), "PF1".into()), plain)], vec![]).unwrap();
        let path = temporary("good.rkyv");
        crate::rust_cache::write_cache(Path::new(&path), SCHEMA, &bytes).unwrap();
        let cache = InterproCache::open(&path).unwrap();
        assert_eq!(cache.len(), 1);
        assert!(cache.get("PFAM", "missing").is_none());
        let list = cache.plain("PFAM", "PF1").unwrap();
        assert_eq!(list[0], list[1]);
        assert_eq!(list[0].abstract_text, "complete abstract \u{00e9}\ntext");
        assert_eq!(
            list[0].go_terms,
            [("GO:1".into(), "function".into(), "binding".into())]
        );
        let mut export_bytes = Vec::new();
        assert_eq!(export(&path, &mut export_bytes).unwrap(), 1);
        assert!(!export_bytes.is_empty());
    }

    #[test]
    fn refuses_wrong_shapes_at_every_level() {
        assert!(plain_key("f", &list(vec![text("PFAM"), text("PF1")])).is_err());
        assert!(plain_key("f", &tuple(vec![text("PFAM")])).is_err());
        assert!(plain_key("f", &tuple(vec![text("PFAM"), Value::Int(1)])).is_err());
        assert!(plain_list("f", &tuple(vec![entry()])).is_err());
        assert!(plain_entry("f", &list(vec![])).is_err());
        assert!(plain_entry("f", &tuple(vec![text("x"); 5])).is_err());
        assert!(plain_entry("f", &tuple(vec![text("x"); 7])).is_err());
        let Value::Tuple(good) = entry() else {
            unreachable!()
        };
        for replacement in [
            tuple(vec![]),
            list(vec![list(vec![text("x"); 3])]),
            list(vec![tuple(vec![text("x"); 2])]),
            list(vec![tuple(vec![text("x"), text("y"), Value::None])]),
        ] {
            let mut fields = good.as_ref().clone();
            fields[4] = replacement;
            assert!(plain_entry("f", &tuple(fields)).is_err());
        }
        for index in [0, 1, 2, 3, 5] {
            let mut fields = good.as_ref().clone();
            fields[index] = Value::Int(1);
            assert!(plain_entry("f", &tuple(fields)).is_err());
        }
        // The restricted reader retains duplicate-assignment metadata rather than silently
        // losing this distinction when it constructs the final Python dict value.
        let duplicate =
            pickle_reader::load(b"\x80\x04}(\x8c\x01D\x8c\x01A\x86]\x8c\x01D\x8c\x01A\x86]u.")
                .unwrap();
        assert!(
            plain_mapping("f", &duplicate)
                .unwrap_err()
                .to_string()
                .contains("twice")
        );
    }

    fn open_raw(name: &str, value: &InterproMapping) -> Result<InterproCache> {
        let path = temporary(name);
        crate::rust_cache::write_cache(
            Path::new(&path),
            SCHEMA,
            &rkyv::to_bytes::<rancor::Error>(value)?,
        )?;
        InterproCache::open(&path)
    }

    #[test]
    fn validates_key_order_lengths_and_every_pool_index() {
        let key = |db: &str, acc: &str| InterproKey {
            database: db.into(),
            accession: acc.into(),
        };
        let mut root = InterproMapping {
            provenance: vec![],
            texts: vec!["x".into()],
            keys: vec![],
            entries: vec![],
        };
        assert_eq!(open_raw("empty", &root).unwrap().len(), 0);
        root.keys = vec![key("Z", "a"), key("A", "a")];
        root.entries = vec![vec![], vec![]];
        assert!(open_raw("unsorted", &root).is_err());
        root.keys = vec![key("A", "a"), key("A", "a")];
        assert!(open_raw("duplicates", &root).is_err());
        root.keys.pop();
        assert!(open_raw("unequal", &root).is_err());
        root.entries.pop();
        for index in 0..8 {
            let mut ids = [0_u32; 8];
            ids[index] = 1;
            root.entries[0] = vec![InterproEntry {
                id: ids[0],
                kind: ids[1],
                short_name: ids[2],
                name: ids[3],
                abstract_text: ids[4],
                go_terms: vec![GoTerm {
                    id: ids[5],
                    category: ids[6],
                    description: ids[7],
                }],
            }];
            assert!(open_raw(&format!("bad_index_{index}"), &root).is_err());
        }
        root.entries[0].clear();
        root.texts.push("x".into());
        assert!(open_raw("duplicate_pool", &root).is_err());
    }
}
