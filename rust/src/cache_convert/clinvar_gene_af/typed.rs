//! Schema 6 v2 adds Boolean scalars without changing the v1 archive layout.
//! Values are retained as Python emits them; arithmetic belongs to consumers.
use super::{FIELDS, GeneAfEntry, GeneAfValue};
use crate::acmg_step::pickle_dict_stream::Item;
use crate::cache_convert::Provenance;
use crate::cache_convert::common::ExactValue;
use crate::cli::{Result, invalid};
use crate::rust_cache::Schema;
use rkyv::{Archive, Deserialize, Serialize};

pub(super) const SCHEMA: Schema = Schema {
    name: "clinvar_patho_af_stat",
    id: 6,
    version: 2,
};

#[derive(Archive, Serialize, Deserialize, Debug, Clone)]
pub(crate) enum Value {
    Text(String),
    Int(i64),
    Float(f64),
    Bool(bool),
}
#[derive(Archive, Serialize, Deserialize, Debug, Clone)]
pub(crate) struct Entry {
    pub(crate) key: String,
    pub(crate) value: Value,
}
#[derive(Archive, Serialize, Deserialize, Debug)]
pub(super) struct Root {
    pub(super) provenance: Vec<Provenance>,
    pub(super) genes: Vec<String>,
    pub(super) entries: Vec<Option<Vec<Entry>>>,
}

#[derive(Clone, Copy)]
pub(crate) enum ValueRef<'a> {
    Text(&'a str),
    Int(i64),
    Float(f64),
    Bool(bool),
}
impl ValueRef<'_> {
    pub(crate) fn owned(self) -> Value {
        match self {
            Self::Text(s) => Value::Text(s.to_owned()),
            Self::Int(i) => Value::Int(i),
            Self::Float(f) => Value::Float(f),
            Self::Bool(b) => Value::Bool(b),
        }
    }
    pub(super) fn exact(self) -> ExactValue {
        match self {
            Self::Text(s) => ExactValue::Text(s.to_owned()),
            Self::Int(i) => ExactValue::Int(i),
            Self::Float(f) => ExactValue::Float(f.to_bits()),
            Self::Bool(b) => ExactValue::Bool(b),
        }
    }
}
impl Value {
    pub(super) fn legacy(&self) -> Option<GeneAfValue> {
        Some(match self {
            Self::Text(s) => GeneAfValue::Text(s.clone()),
            Self::Int(i) => GeneAfValue::Int(*i),
            Self::Float(f) => GeneAfValue::Float(*f),
            Self::Bool(_) => return None,
        })
    }
    fn exact(&self) -> ExactValue {
        match self {
            Self::Text(s) => ExactValue::Text(s.clone()),
            Self::Int(i) => ExactValue::Int(*i),
            Self::Float(f) => ExactValue::Float(f.to_bits()),
            Self::Bool(b) => ExactValue::Bool(*b),
        }
    }
}

/// A borrowed record in either schema. Neither lookup copies its strings.
#[derive(Clone, Copy)]
pub(crate) enum RecordRef<'a> {
    V1(&'a rkyv::vec::ArchivedVec<super::ArchivedGeneAfEntry>),
    V2(&'a rkyv::vec::ArchivedVec<ArchivedEntry>),
}
pub(crate) struct EntryRef<'a> {
    pub(crate) key: &'a str,
    pub(crate) value: ValueRef<'a>,
}
impl<'a> RecordRef<'a> {
    pub(crate) fn len(self) -> usize {
        match self {
            Self::V1(v) => v.len(),
            Self::V2(v) => v.len(),
        }
    }
    pub(crate) fn iter(self) -> impl Iterator<Item = EntryRef<'a>> {
        (0..self.len()).map(move |i| match self {
            Self::V1(v) => EntryRef {
                key: v[i].key.as_str(),
                value: match &v[i].value {
                    super::ArchivedGeneAfValue::Text(s) => ValueRef::Text(s.as_str()),
                    super::ArchivedGeneAfValue::Int(n) => ValueRef::Int(n.to_native()),
                    super::ArchivedGeneAfValue::Float(n) => ValueRef::Float(n.to_native()),
                },
            },
            Self::V2(v) => EntryRef {
                key: v[i].key.as_str(),
                value: match &v[i].value {
                    ArchivedValue::Text(s) => ValueRef::Text(s.as_str()),
                    ArchivedValue::Int(n) => ValueRef::Int(n.to_native()),
                    ArchivedValue::Float(n) => ValueRef::Float(n.to_native()),
                    ArchivedValue::Bool(b) => ValueRef::Bool(*b),
                },
            },
        })
    }
}

pub(super) fn legacy(
    entries: &[(String, Option<Vec<Entry>>)],
) -> Option<Vec<(String, Option<Vec<GeneAfEntry>>)>> {
    entries
        .iter()
        .map(|(key, record)| {
            let record = match record {
                None => None,
                Some(record) => Some(
                    record
                        .iter()
                        .map(|entry| {
                            Some(GeneAfEntry {
                                key: entry.key.clone(),
                                value: entry.value.legacy()?,
                            })
                        })
                        .collect::<Option<Vec<_>>>()?,
                ),
            };
            Some((key.clone(), record))
        })
        .collect()
}

pub(super) fn plain_gene(file: &str, gene: &str, input: &Item) -> Result<Option<Vec<Entry>>> {
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
                Item::Str(s) => Value::Text(s.to_string()),
                Item::Int(i) => Value::Int(*i),
                Item::Float(f) => Value::Float(*f),
                Item::Bool(b) => Value::Bool(*b),
                other => {
                    return Err(odd(format!(
                        "field {key} is a {}, not str/int/float/bool",
                        other.type_name()
                    )));
                }
            };
            Ok(Entry {
                key: key.to_string(),
                value,
            })
        })
        .collect::<Result<Vec<_>>>()
        .map(Some)
}

pub(super) fn exact_plain(record: Option<&[Entry]>) -> ExactValue {
    match record {
        None => ExactValue::None,
        Some(entries) => ExactValue::Dict(
            entries
                .iter()
                .map(|e| (ExactValue::Text(e.key.clone()), e.value.exact()))
                .collect(),
        ),
    }
}
