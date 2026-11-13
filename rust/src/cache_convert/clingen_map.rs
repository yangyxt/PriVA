//! Schema 10: unchanged version 1 for text/NaN; version 2 for typed scalar codes.
//!
//! The v1 representation and its NaN-to-missing normalization are frozen. V2
//! additionally retains integer, Boolean and non-NaN float codes. Python None
//! and containers remain refused. Lookup keys are sorted by bytes; missing keys
//! remain distinct from present non-string codes that cannot override criteria.
use super::Provenance;
use super::common::{ExactValue, provenance, write_entry};
use crate::acmg_step::pickle_dict_stream::{
    Item, visit_top_level_dict_with_uint64 as visit_top_level_dict,
};
use crate::cli::{Result, invalid};
use crate::rust_cache::{
    MappedCache, Schema, ensure_distinct_paths, find_sorted, write_cache_checked,
};
use rkyv::{Archive, Deserialize, Serialize, rancor};
use std::io::{Read, Write};
use std::path::Path;
use std::time::Instant;

mod typed;
pub(crate) use typed::Code;

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

fn read_pickle(path: &str) -> Result<Vec<(String, Code)>> {
    let mut entries = Vec::new();
    let distinct = visit_top_level_dict(
        path,
        |_| true,
        |variant, value| {
            entries.push((variant.to_string(), Code::from_item(path, &value)?));
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

/// Managed installation names must identify the actual schema version before
/// publication. This reads only the small completed map, never its VCF source.
pub(crate) fn pickle_schema_version(path: &str) -> Result<u32> {
    Ok(
        if read_pickle(path)?
            .iter()
            .all(|(_, code)| matches!(code, Code::Text(_) | Code::Missing))
        {
            SCHEMA.version
        } else {
            typed::SCHEMA.version
        },
    )
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

/// Publish a producer's already reduced text/NaN values using the same schema
/// and complete mmap lookup check as the pickle converter. Duplicate keys are
/// an error here; source-specific last-wins reduction belongs to the producer.
pub(crate) fn write_entries(
    output: &Path,
    entries: &[(String, Option<String>)],
    provenance: Vec<Provenance>,
) -> Result<()> {
    let archive = build(entries.to_vec(), provenance)?;
    write_cache_checked(output, SCHEMA, &archive, |staged| {
        let cache = ClingenCache::open(&staged.to_string_lossy())?;
        if cache.len() != entries.len() {
            return Err(invalid("the ClinGen map has a different key count"));
        }
        for (variant, code) in entries {
            if cache.get(variant) != Some(code.as_deref()) {
                return Err(invalid("a ClinGen code differs from the producer"));
            }
        }
        Ok(())
    })
}

fn archive_for_codes(
    entries: &[(String, Code)],
    provenance: Vec<Provenance>,
) -> Result<(Schema, Vec<u8>)> {
    let legacy: Option<Vec<_>> = entries
        .iter()
        .map(|(key, code)| code.legacy().map(|value| (key.clone(), value)))
        .collect();
    if let Some(entries) = legacy {
        // Preserve the complete version-1 archive bytes for every old input.
        return Ok((SCHEMA, build(entries, provenance)?));
    }
    let mut entries = entries.to_vec();
    entries.sort_by(|a, b| a.0.cmp(&b.0));
    if entries.windows(2).any(|pair| pair[0].0 == pair[1].0) {
        return Err(invalid("the ClinGen variant keys are not unique"));
    }
    let (variants, codes) = entries.into_iter().unzip();
    Ok((
        typed::SCHEMA,
        rkyv::to_bytes::<rancor::Error>(&typed::ClingenMap {
            provenance,
            variants,
            codes,
        })?
        .to_vec(),
    ))
}

/// Select v2 only when a surviving code cannot be represented by v1.
pub(crate) fn write_codes(
    output: &Path,
    entries: &[(String, Code)],
    provenance: Vec<Provenance>,
) -> Result<()> {
    let legacy: Option<Vec<_>> = entries
        .iter()
        .map(|(key, code)| code.legacy().map(|value| (key.clone(), value)))
        .collect();
    if let Some(entries) = legacy {
        return write_entries(output, &entries, provenance);
    }
    let (schema, archive) = archive_for_codes(entries, provenance)?;
    write_cache_checked(output, schema, &archive, |staged| {
        let cache = ClingenCache::open(&staged.to_string_lossy())?;
        if cache.len() != entries.len()
            || entries
                .iter()
                .any(|(key, code)| cache.code(key).as_ref() != Some(code))
        {
            return Err(invalid("a typed ClinGen code differs from the producer"));
        }
        Ok(())
    })
}

enum Mapping {
    V1(MappedCache<ClingenMap>),
    V2(MappedCache<typed::ClingenMap>),
}

/// A checked mmap view of the ClinGen map. Lookups do not copy the key/code vectors.
pub(crate) struct ClingenCache {
    mapped: Mapping,
}

impl ClingenCache {
    pub(crate) fn open(path: &str) -> Result<Self> {
        let mut file = std::fs::File::open(path)
            .map_err(|error| invalid(format!("FileNotFoundError: {path}: {error}")))?;
        let mut prefix = [0_u8; 20];
        let v2 = file.read_exact(&mut prefix).is_ok()
            && &prefix[..8] == crate::rust_cache::MAGIC
            && u32::from_le_bytes(prefix[12..16].try_into().unwrap()) == SCHEMA.id
            && u32::from_le_bytes(prefix[16..20].try_into().unwrap()) == typed::SCHEMA.version;
        let mapped = if v2 {
            Mapping::V2(MappedCache::<typed::ClingenMap>::open_file(
                &file,
                path,
                typed::SCHEMA,
            )?)
        } else {
            Mapping::V1(MappedCache::<ClingenMap>::open_file(&file, path, SCHEMA)?)
        };
        let cache = Self { mapped };
        let code_count = match &cache.mapped {
            Mapping::V1(mapped) => mapped.root().codes.len(),
            Mapping::V2(mapped) => {
                if mapped.root().codes.iter().any(|code| !code.valid()) {
                    return Err(invalid(format!(
                        "{path}: invalid typed ClinGen scalar representation"
                    )));
                }
                mapped.root().codes.len()
            }
        };
        if cache.variants().len() != code_count
            || cache
                .variants()
                .windows(2)
                .any(|pair| pair[0].as_str() >= pair[1].as_str())
        {
            return Err(invalid(format!(
                "{path}: the ClinGen variants are not sorted and unique, or not one per code"
            )));
        }
        Ok(cache)
    }

    fn variants(&self) -> &[rkyv::string::ArchivedString] {
        match &self.mapped {
            Mapping::V1(mapped) => &mapped.root().variants,
            Mapping::V2(mapped) => &mapped.root().variants,
        }
    }

    pub(crate) fn len(&self) -> usize {
        self.variants().len()
    }

    pub(crate) fn schema_version(&self) -> u32 {
        match &self.mapped {
            Mapping::V1(_) => SCHEMA.version,
            Mapping::V2(_) => typed::SCHEMA.version,
        }
    }

    pub(crate) fn keys(&self) -> impl Iterator<Item = &str> {
        self.variants().iter().map(|key| key.as_str())
    }

    /// Runtime text view: absent key is None; every present non-string is
    /// Some(None). Python likewise overrides criteria only for nonempty strings.
    pub(crate) fn get(&self, variant: &str) -> Option<Option<&str>> {
        let at = find_sorted(self.variants(), variant)?;
        Some(match &self.mapped {
            Mapping::V1(mapped) => mapped.root().codes[at].as_ref().map(|code| code.as_str()),
            Mapping::V2(mapped) => mapped.root().codes[at].text(),
        })
    }

    /// Plain values for independent exporters; the only normalization is NaN to `None`.
    pub(crate) fn code(&self, variant: &str) -> Option<Code> {
        let at = find_sorted(self.variants(), variant)?;
        Some(match &self.mapped {
            Mapping::V1(mapped) => mapped.root().codes[at]
                .as_ref()
                .map(|code| code.as_str().to_owned())
                .into(),
            Mapping::V2(mapped) => mapped.root().codes[at].plain(),
        })
    }
}

pub(crate) fn convert(pickle: &str, output: &Path) -> Result<String> {
    ensure_distinct_paths(Path::new(pickle), output)?;
    let start = Instant::now();
    let entries = read_pickle(pickle)?;
    let count = entries.len();
    let without_codes = entries.iter().filter(|(_, code)| code.is_missing()).count();
    let (schema, archive) = archive_for_codes(&entries, provenance(pickle)?)?;
    let converted_s = start.elapsed().as_secs_f64();
    let check = Instant::now();
    let mut compared = 0_usize;
    write_cache_checked(output, schema, &archive, |staged| {
        let cache = ClingenCache::open(&staged.to_string_lossy())?;
        if cache.len() != count {
            return Err(invalid("the ClinGen map has a different key count"));
        }
        let distinct = visit_top_level_dict(
            pickle,
            |_| true,
            |variant, value| {
                if cache.code(&variant) != Some(Code::from_item(pickle, &value)?) {
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
        let value = cache.code(key).unwrap().exact();
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

    #[test]
    fn typed_codes_keep_types_bits_and_string_only_runtime_overrides() {
        let entries = vec![
            ("text".into(), Code::Text("PS3".into())),
            ("missing".into(), Code::Missing),
            ("signed".into(), Code::Integer(i128::from(i64::MIN))),
            ("unsigned".into(), Code::Integer(i128::from(u64::MAX))),
            ("bool".into(), Code::Boolean(true)),
            ("negative_zero".into(), Code::Float((-0.0_f64).to_bits())),
            ("infinity".into(), Code::Float(f64::INFINITY.to_bits())),
        ];
        let path = temporary("typed.rkyv", b"previous");
        write_codes(Path::new(&path), &entries, vec![]).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        assert_eq!(u32::from_le_bytes(bytes[16..20].try_into().unwrap()), 2);
        let cache = ClingenCache::open(&path).unwrap();
        for (key, code) in &entries {
            assert_eq!(cache.code(key).as_ref(), Some(code));
            assert_eq!(
                cache.get(key),
                Some(if key == "text" { Some("PS3") } else { None })
            );
        }
        assert_eq!(cache.get("absent"), None);
        let mut exported = Vec::new();
        export(&path, &mut exported).unwrap();
        let text = String::from_utf8(exported).unwrap();
        assert!(text.contains("[\"int\",\"18446744073709551615\"]"));
        assert!(text.contains("[\"float\",\"8000000000000000\"]"));
        assert!(text.contains("[\"bool\",true]"));
    }

    #[test]
    fn old_values_use_the_exact_v1_layout_and_invalid_v2_preserves_output() {
        let old = vec![("a".into(), Some("PS3".into())), ("b".into(), None)];
        let codes: Vec<_> = old
            .iter()
            .map(|(key, code): &(String, Option<String>)| (key.clone(), code.clone().into()))
            .collect();
        let (schema, bytes) = archive_for_codes(&codes, vec![]).unwrap();
        assert_eq!(schema, SCHEMA);
        assert_eq!(bytes, build(old, vec![]).unwrap());
        let path = temporary("invalid-typed.rkyv", b"previous");
        for code in [
            Code::Integer(i128::from(u64::MAX) + 1),
            Code::Float(f64::NAN.to_bits()),
        ] {
            assert!(write_codes(Path::new(&path), &[("a".into(), code)], vec![]).is_err());
            assert_eq!(std::fs::read(&path).unwrap(), b"previous");
        }
        for (variants, codes) in [
            (
                vec!["z".into(), "a".into()],
                vec![Code::Integer(1), Code::Integer(2)],
            ),
            (vec!["a".into()], vec![]),
            (
                vec!["a".into(), "a".into()],
                vec![Code::Integer(1), Code::Integer(2)],
            ),
        ] {
            let bytes = rkyv::to_bytes::<rancor::Error>(&typed::ClingenMap {
                provenance: vec![],
                variants,
                codes,
            })
            .unwrap();
            write_cache(Path::new(&path), typed::SCHEMA, &bytes).unwrap();
            assert!(ClingenCache::open(&path).is_err());
        }
    }
}
