//! Schema 8, version 1: the DAS intolerant-domain Python `set[str]`.
//!
//! Set order is not meaningful; values are stored once, in UTF-8 byte order. Lists,
//! tuples and frozen sets are different Python types and are refused rather than
//! silently coerced to a set. The archive contains no Python object tree.
use super::Provenance;
use super::common::{ExactValue, provenance, write_entry};
use crate::cli::{Result, invalid};
use crate::pickle_reader::{self, Value};
use crate::rust_cache::{
    MappedCache, Schema, ensure_distinct_paths, find_sorted, write_cache_checked,
};
use rkyv::{Archive, Deserialize, Serialize, rancor};
use std::io::Write;
use std::path::Path;
use std::time::Instant;

pub(crate) const ALL_INTOLERANT_DOMAINS: Schema = Schema {
    name: "all_intolerant_domains",
    id: 8,
    version: 1,
};

#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct IntolerantDomainSet {
    pub(crate) provenance: Vec<Provenance>,
    /// The Python set's texts, sorted by UTF-8 bytes and unique.
    pub(crate) domains: Vec<String>,
}

pub(super) fn plain_set(path: &str, value: &Value) -> Result<Vec<String>> {
    let Value::Set(set) = value else {
        return Err(invalid(format!(
            "{path}: all_intolerant_domains is a {}, not a set; refused",
            value.type_name()
        )));
    };
    let mut domains = set
        .borrow()
        .items
        .iter()
        .map(|value| match value {
            Value::Str(text) => Ok(text.to_string()),
            _ => Err(invalid(format!(
                "{path}: an intolerant-domain set member is a {}, not text; refused",
                value.type_name()
            ))),
        })
        .collect::<Result<Vec<_>>>()?;
    domains.sort();
    if domains.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(invalid(format!(
            "{path}: repeated domain set member; refused"
        )));
    }
    Ok(domains)
}

fn read_pickle(path: &str) -> Result<Vec<String>> {
    let bytes = pickle_reader::read_file_by_magic(Path::new(path))?;
    plain_set(path, &pickle_reader::load(&bytes)?)
}

/// The validated set stays memory mapped; membership is a binary search.
pub(crate) struct IntolerantDomainCache {
    mapped: MappedCache<IntolerantDomainSet>,
}

impl IntolerantDomainCache {
    pub(crate) fn open(path: &str) -> Result<Self> {
        let mapped = MappedCache::<IntolerantDomainSet>::open(path, ALL_INTOLERANT_DOMAINS)?;
        if mapped
            .root()
            .domains
            .windows(2)
            .any(|pair| pair[0].as_str() >= pair[1].as_str())
        {
            return Err(invalid(format!(
                "{path}: domain set values are not sorted and unique"
            )));
        }
        Ok(Self { mapped })
    }

    pub(crate) fn len(&self) -> usize {
        self.mapped.root().domains.len()
    }

    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn contains(&self, domain: &str) -> bool {
        find_sorted(&self.mapped.root().domains, domain).is_some()
    }

    pub(crate) fn domains(&self) -> impl ExactSizeIterator<Item = &str> {
        self.mapped
            .root()
            .domains
            .iter()
            .map(|domain| domain.as_str())
    }
}

pub(crate) fn convert(pickle: &str, output: &Path) -> Result<String> {
    ensure_distinct_paths(Path::new(pickle), output)?;
    let start = Instant::now();
    let domains = read_pickle(pickle)?;
    let count = domains.len();
    let root = IntolerantDomainSet {
        provenance: provenance(pickle)?,
        domains,
    };
    let archive = rkyv::to_bytes::<rancor::Error>(&root)?;
    write_cache_checked(output, ALL_INTOLERANT_DOMAINS, &archive, |staged| {
        let path = staged
            .to_str()
            .ok_or_else(|| invalid("the cache path is not UTF-8"))?;
        let cache = IntolerantDomainCache::open(path)?;
        let expected = read_pickle(pickle)?;
        if cache.len() != count
            || expected.len() != count
            || !cache.domains().eq(expected.iter().map(String::as_str))
        {
            return Err(invalid(format!(
                "{path}: the converted intolerant-domain set differs from the pickle"
            )));
        }
        Ok(())
    })?;
    Ok(format!(
        "{{\"cache\":\"all_intolerant_domains\",\"domains\":{count},\"checked_domains\":{count},\"bytes\":{},\"elapsed_s\":{:.3}}}",
        std::fs::metadata(output)?.len(),
        start.elapsed().as_secs_f64()
    ))
}

/// A set is one full exported value, keyed by None to distinguish it from dict entries.
pub(crate) fn export(path: &str, out: &mut dyn Write) -> Result<usize> {
    let cache = IntolerantDomainCache::open(path)?;
    write_entry(
        out,
        &ExactValue::None,
        &ExactValue::Set(
            cache
                .domains()
                .map(|text| ExactValue::Text(text.to_owned()))
                .collect(),
        ),
    )?;
    Ok(1)
}

#[cfg(test)]
mod tests {
    use super::*;

    const GOOD: &[u8] = b"\x80\x04\x95\x14\x00\x00\x00\x00\x00\x00\x00\x8f\x94\x28\x8c\x02\xc3\xa9\x94\x8c\x02\x64\x31\x94\x8c\x02\x64\x32\x94\x90\x2e";

    fn temporary(name: &str, bytes: &[u8]) -> String {
        let path = std::env::temp_dir().join(format!(
            "priva_intolerant_domains_{}_{name}",
            std::process::id()
        ));
        std::fs::write(&path, bytes).unwrap();
        path.to_str().unwrap().to_owned()
    }

    #[test]
    fn maps_exact_set_and_empty_set() {
        for (at, source) in [GOOD, b"\x80\x04\x8f."].iter().enumerate() {
            let pickle = temporary(&format!("good{at}.pkl"), source);
            let output = temporary(&format!("good{at}.rkyv"), b"");
            convert(&pickle, Path::new(&output)).unwrap();
            let cache = IntolerantDomainCache::open(&output).unwrap();
            assert_eq!(cache.len(), if at == 0 { 3 } else { 0 });
            assert_eq!(cache.contains("d2"), at == 0);
            assert!(!cache.contains("absent"));
            if at == 0 {
                assert_eq!(cache.domains().collect::<Vec<_>>(), ["d1", "d2", "é"]);
            }
            let mut exported = Vec::new();
            assert_eq!(export(&output, &mut exported).unwrap(), 1);
            let again = temporary(&format!("again{at}.rkyv"), b"");
            convert(&pickle, Path::new(&again)).unwrap();
            assert_eq!(
                std::fs::read(&output).unwrap(),
                std::fs::read(&again).unwrap()
            );
        }
    }

    #[test]
    fn refuses_other_container_types_and_nontext_members_without_replacing_output() {
        let output = temporary("refuse.rkyv", b"unchanged");
        for (at, source) in [
            b"\x80\x04].".as_slice(),
            b"\x80\x04(\x91.",
            b"\x80\x04\x8f(K\x01\x90.",
        ]
        .iter()
        .enumerate()
        {
            let pickle = temporary(&format!("bad{at}.pkl"), source);
            assert!(convert(&pickle, Path::new(&output)).is_err());
            assert_eq!(std::fs::read(&output).unwrap(), b"unchanged");
        }
    }

    #[test]
    fn reader_rejects_unsorted_and_repeated_set_values() {
        let output = temporary("bad_order.rkyv", b"");
        for domains in [vec!["b".into(), "a".into()], vec!["a".into(), "a".into()]] {
            let archive = rkyv::to_bytes::<rancor::Error>(&IntolerantDomainSet {
                provenance: vec![],
                domains,
            })
            .unwrap();
            crate::rust_cache::write_cache(Path::new(&output), ALL_INTOLERANT_DOMAINS, &archive)
                .unwrap();
            assert!(IntolerantDomainCache::open(&output).is_err());
        }
    }
}
