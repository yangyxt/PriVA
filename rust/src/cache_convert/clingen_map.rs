//! Schema 10, version 1: the Python ClinGen map `{variant_id: str | float('nan')}`.
//!
//! `prepare_clingen_map.py` writes text evidence codes, or a float NaN for a curated variant
//! without codes. Only that NaN is normalized to `None` in `codes`. Python `None`, finite
//! floats and all other scalar/container types are refused. Top-level keys are sorted by
//! bytes; a missing key remains distinguishable from a present key without evidence codes.
use super::Provenance;
use super::common::{ExactValue, provenance, write_entry};
use crate::acmg_step::pickle_dict_stream::{Item, visit_top_level_dict};
use crate::cli::{Result, invalid};
use crate::rust_cache::{
    MappedCache, Schema, ensure_distinct_paths, find_sorted, write_cache_checked,
};
use rkyv::{Archive, Deserialize, Serialize, rancor};
use std::io::Write;
use std::path::Path;
use std::time::Instant;

pub(crate) const SCHEMA: Schema = Schema {
    name: "clingen_map",
    id: 10,
    version: 1,
};

#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct ClingenMap {
    pub(crate) provenance: Vec<Provenance>,
    /// Python dict keys, sorted by bytes and unique; parallel to `codes`.
    pub(crate) variants: Vec<String>,
    /// Python text, or `None` for the producer's float NaN only.
    pub(crate) codes: Vec<Option<String>>,
}

pub(crate) fn plain_code(path: &str, value: &Item) -> Result<Option<String>> {
    match value {
        Item::Str(text) => Ok(Some(text.to_string())),
        Item::Float(number) if number.is_nan() => Ok(None),
        other => Err(invalid(format!(
            "{path}: a ClinGen code is a {}, not text or float NaN; refused",
            other.type_name()
        ))),
    }
}

fn read_pickle(path: &str) -> Result<Vec<(String, Option<String>)>> {
    let mut entries = Vec::new();
    let distinct = visit_top_level_dict(
        path,
        |_| true,
        |variant, value| {
            entries.push((variant.to_string(), plain_code(path, &value)?));
            Ok(())
        },
    )?;
    if distinct != entries.len() {
        return Err(invalid(format!(
            "{path}: a ClinGen variant is set twice in the pickle; refused"
        )));
    }
    Ok(entries)
}

fn build(
    mut entries: Vec<(String, Option<String>)>,
    provenance: Vec<Provenance>,
) -> Result<Vec<u8>> {
    entries.sort_by(|left, right| left.0.cmp(&right.0));
    if entries.windows(2).any(|pair| pair[0].0 == pair[1].0) {
        return Err(invalid("the ClinGen variant keys are not unique"));
    }
    let (variants, codes) = entries.into_iter().unzip();
    Ok(rkyv::to_bytes::<rancor::Error>(&ClingenMap {
        provenance,
        variants,
        codes,
    })?
    .to_vec())
}

/// A checked mmap view of the ClinGen map. Lookups do not copy the key/code vectors.
pub(crate) struct ClingenCache {
    mapped: MappedCache<ClingenMap>,
}

impl ClingenCache {
    pub(crate) fn open(path: &str) -> Result<Self> {
        let mapped = MappedCache::<ClingenMap>::open(path, SCHEMA)?;
        let root = mapped.root();
        if root.variants.len() != root.codes.len()
            || root
                .variants
                .windows(2)
                .any(|pair| pair[0].as_str() >= pair[1].as_str())
        {
            return Err(invalid(format!(
                "{path}: the ClinGen variants are not sorted and unique, or not one per code"
            )));
        }
        Ok(Self { mapped })
    }

    pub(crate) fn len(&self) -> usize {
        self.mapped.root().variants.len()
    }

    pub(crate) fn keys(&self) -> impl Iterator<Item = &str> {
        self.mapped.root().variants.iter().map(|key| key.as_str())
    }

    /// `None` is an absent variant; `Some(None)` is the producer's NaN/no-code value.
    pub(crate) fn get(&self, variant: &str) -> Option<Option<&str>> {
        let root = self.mapped.root();
        find_sorted(&root.variants, variant)
            .map(|at| root.codes[at].as_ref().map(|code| code.as_str()))
    }

    /// Plain values for independent exporters; the only normalization is NaN to `None`.
    pub(crate) fn plain(&self, variant: &str) -> Option<Option<String>> {
        self.get(variant).map(|code| code.map(str::to_owned))
    }
}

pub(crate) fn convert(pickle: &str, output: &Path) -> Result<String> {
    ensure_distinct_paths(Path::new(pickle), output)?;
    let start = Instant::now();
    let entries = read_pickle(pickle)?;
    let count = entries.len();
    let without_codes = entries.iter().filter(|(_, code)| code.is_none()).count();
    let archive = build(entries, provenance(pickle)?)?;
    let converted_s = start.elapsed().as_secs_f64();
    let check = Instant::now();
    let mut compared = 0_usize;
    write_cache_checked(output, SCHEMA, &archive, |staged| {
        let cache = ClingenCache::open(&staged.to_string_lossy())?;
        if cache.len() != count {
            return Err(invalid("the ClinGen map has a different key count"));
        }
        let distinct = visit_top_level_dict(
            pickle,
            |_| true,
            |variant, value| {
                if cache.plain(&variant) != Some(plain_code(pickle, &value)?) {
                    return Err(invalid("a ClinGen code differs from the pickle"));
                }
                compared += 1;
                Ok(())
            },
        )?;
        if compared != count || distinct != count {
            return Err(invalid("the ClinGen pickle changed while it was converted"));
        }
        Ok(())
    })?;
    Ok(format!(
        "{{\"cache\":\"clingen_map\",\"variants\":{count},\"without_codes\":{without_codes},\"bytes\":{},\"converted_s\":{converted_s:.3},\"checked_keys\":{compared},\"check_s\":{:.3}}}",
        std::fs::metadata(output)?.len(),
        check.elapsed().as_secs_f64()
    ))
}

pub(crate) fn export(path: &str, out: &mut dyn Write) -> Result<usize> {
    let cache = ClingenCache::open(path)?;
    for key in cache.keys() {
        // `None` is the explicit schema-10 normalization of a present float NaN.
        let value = match cache.get(key).unwrap() {
            Some(code) => ExactValue::Text(code.to_owned()),
            None => ExactValue::None,
        };
        write_entry(out, &ExactValue::Text(key.to_owned()), &value)?;
    }
    Ok(cache.len())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rust_cache::write_cache;
    use std::io::Write;

    fn temporary(name: &str, bytes: &[u8]) -> String {
        let path =
            std::env::temp_dir().join(format!("priva_clingen_cache_{}_{name}", std::process::id()));
        std::fs::File::create(&path)
            .unwrap()
            .write_all(bytes)
            .unwrap();
        path.to_str().unwrap().to_owned()
    }

    // Protocol 4: {"z": NaN (payload 0x123), "a": "", "b": "PS3, PM2"}.
    const GOOD: &[u8] = b"\x80\x04}(\x8c\x01zG\x7f\xf8\x00\x00\x00\x00\x01\x23\x8c\x01a\x8c\x00\x8c\x01b\x8c\x08PS3, PM2u.";

    #[test]
    fn converts_text_and_nan_but_distinguishes_a_missing_variant() {
        let pickle = temporary("good.pkl", GOOD);
        let output = temporary("good.rkyv", b"");
        let summary = convert(&pickle, Path::new(&output)).unwrap();
        assert!(summary.contains("\"variants\":3"), "{summary}");
        assert!(summary.contains("\"without_codes\":1"), "{summary}");
        let cache = ClingenCache::open(&output).unwrap();
        assert_eq!(cache.keys().collect::<Vec<_>>(), ["a", "b", "z"]);
        assert_eq!(cache.get("a"), Some(Some("")));
        assert_eq!(cache.get("b"), Some(Some("PS3, PM2")));
        assert_eq!(cache.get("z"), Some(None));
        assert_eq!(cache.get("missing"), None);
        let again = temporary("again.rkyv", b"");
        convert(&pickle, Path::new(&again)).unwrap();
        assert_eq!(
            std::fs::read(output).unwrap(),
            std::fs::read(again).unwrap()
        );
    }

    #[test]
    fn refuses_every_other_value_type_and_repeated_top_keys() {
        for value in [
            Item::None,
            Item::Bool(true),
            Item::Int(1),
            Item::Float(0.0),
            Item::Float(f64::INFINITY),
            Item::List(Default::default()),
            Item::Tuple(Default::default()),
            Item::Dict(Default::default()),
        ] {
            assert!(plain_code("test", &value).is_err());
        }
        let duplicate = temporary(
            "duplicate.pkl",
            b"\x80\x04}(\x8c\x01a\x8c\x01x\x8c\x01a\x8c\x01yu.",
        );
        assert!(
            read_pickle(&duplicate)
                .unwrap_err()
                .to_string()
                .contains("twice")
        );
    }

    fn open_raw(
        name: &str,
        variants: Vec<String>,
        codes: Vec<Option<String>>,
    ) -> Result<ClingenCache> {
        let archive = rkyv::to_bytes::<rancor::Error>(&ClingenMap {
            provenance: vec![],
            variants,
            codes,
        })?;
        let path = temporary(name, b"");
        write_cache(Path::new(&path), SCHEMA, &archive)?;
        ClingenCache::open(&path)
    }

    #[test]
    fn rejects_invalid_key_vectors_in_a_structurally_valid_archive() {
        assert!(open_raw("unsorted", vec!["z".into(), "a".into()], vec![None, None]).is_err());
        assert!(open_raw("repeated", vec!["a".into(), "a".into()], vec![None, None]).is_err());
        assert!(open_raw("unequal", vec!["a".into()], vec![]).is_err());
        assert_eq!(open_raw("empty", vec![], vec![]).unwrap().len(), 0);
    }
}
