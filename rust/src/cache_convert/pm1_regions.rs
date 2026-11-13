//! Schema 7, version 1: the three PM1 lookup arms and their producer metadata.
//!
//! Python root: four required fixed keys, `das_intolerant_domains: set[str]`,
//! `hcseeker_hotspots: dict[str, set[int]]`, `rmc_constrained: dict[str, set[int]]`
//! and `metadata: dict`. Transcript lookup keys and set members are sorted;
//! metadata keeps field presence, its exact scalar types and insertion order.
//! The closed metadata enum describes the producer's declared fields. Unknown
//! keys and wrong types are errors, never discarded metadata.
use super::Provenance;
use super::common::{ExactValue, provenance, write_entry};
use super::intolerant_domains::plain_set;
use crate::cli::{Result, invalid};
use crate::pickle_reader::{self, Dict, Value};
use crate::rust_cache::{MappedCache, Schema, ensure_distinct_paths, write_cache_checked};
use rkyv::{Archive, Deserialize, Serialize, rancor};
use std::collections::HashSet;
use std::io::Write;
use std::path::Path;
use std::time::Instant;

pub(crate) const PM1_REGIONS_PKL: Schema = Schema {
    name: "pm1_regions_pkl",
    id: 7,
    version: 1,
};

/// One transcript's Python `set[int]`, sorted numerically with no duplicates.
#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct TranscriptResidues {
    pub(crate) transcript: String,
    pub(crate) residues: Vec<i64>,
}

/// Allowed producer metadata keys and Python value types. Enum order in the vector
/// is original metadata insertion order; absent keys remain absent. Float fields
/// are `f64` with exact bits, including signed zero, infinity and NaN payloads.
#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) enum Pm1Metadata {
    Assembly(String),
    RmcOeThreshold(f64),
    RmcPBhThreshold(f64),
    NDasIntolerantDomains(i64),
    NHcseekerHotspotTranscripts(i64),
    NHcseekerHotspotResidues(i64),
    NRmcConstrainedTranscripts(i64),
    NRmcConstrainedResidues(i64),
}

#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct Pm1Regions {
    pub(crate) provenance: Vec<Provenance>,
    pub(crate) das_intolerant_domains: Vec<String>,
    /// UTF-8 byte-sorted, unique transcript lookup keys.
    pub(crate) hcseeker_hotspots: Vec<TranscriptResidues>,
    pub(crate) rmc_constrained: Vec<TranscriptResidues>,
    pub(crate) metadata: Vec<Pm1Metadata>,
}

fn required<'a>(path: &str, dict: &'a Dict, key: &str) -> Result<&'a Value> {
    dict.get(&Value::Str(key.into()))?
        .ok_or_else(|| invalid(format!("{path}: missing required PM1 field {key}; refused")))
}

fn residue_map(path: &str, name: &str, value: &Value) -> Result<Vec<TranscriptResidues>> {
    let Value::Dict(dict) = value else {
        return Err(invalid(format!("{path}: {name} is not a dict; refused")));
    };
    let mut result = Vec::with_capacity(dict.borrow().entries.len());
    for (key, value) in &dict.borrow().entries {
        let Value::Str(transcript) = key else {
            return Err(invalid(format!(
                "{path}: {name} has a nontext transcript key; refused"
            )));
        };
        let Value::Set(set) = value else {
            return Err(invalid(format!(
                "{path}: {name} transcript {transcript} has a {}, not a set; refused",
                value.type_name()
            )));
        };
        let mut residues = set
            .borrow()
            .items
            .iter()
            .map(|value| match value {
                Value::Int(residue) => Ok(*residue),
                _ => Err(invalid(format!(
                    "{path}: {name} transcript {transcript} has a {}, not an int residue; refused",
                    value.type_name()
                ))),
            })
            .collect::<Result<Vec<_>>>()?;
        residues.sort_unstable();
        if residues.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(invalid(format!(
                "{path}: {name} has duplicate residue set members; refused"
            )));
        }
        result.push(TranscriptResidues {
            transcript: transcript.to_string(),
            residues,
        });
    }
    result.sort_by(|a, b| a.transcript.cmp(&b.transcript));
    Ok(result)
}

fn metadata(path: &str, value: &Value) -> Result<Vec<Pm1Metadata>> {
    let Value::Dict(dict) = value else {
        return Err(invalid(format!(
            "{path}: PM1 metadata is not a dict; refused"
        )));
    };
    dict.borrow().entries.iter().map(|(key, value)| {
        Ok(match (key.as_str(), value) {
            (Some("assembly"), Value::Str(value)) => Pm1Metadata::Assembly(value.to_string()),
            (Some("rmc_oe_threshold"), Value::Float(value)) => Pm1Metadata::RmcOeThreshold(*value),
            (Some("rmc_p_bh_threshold"), Value::Float(value)) => Pm1Metadata::RmcPBhThreshold(*value),
            (Some("n_das_intolerant_domains"), Value::Int(value)) => Pm1Metadata::NDasIntolerantDomains(*value),
            (Some("n_hcseeker_hotspot_transcripts"), Value::Int(value)) => Pm1Metadata::NHcseekerHotspotTranscripts(*value),
            (Some("n_hcseeker_hotspot_residues"), Value::Int(value)) => Pm1Metadata::NHcseekerHotspotResidues(*value),
            (Some("n_rmc_constrained_transcripts"), Value::Int(value)) => Pm1Metadata::NRmcConstrainedTranscripts(*value),
            (Some("n_rmc_constrained_residues"), Value::Int(value)) => Pm1Metadata::NRmcConstrainedResidues(*value),
            (Some(key), _) => return Err(invalid(format!("{path}: unknown PM1 metadata field {key:?} or undeclared value type {}; refused", value.type_name()))),
            _ => return Err(invalid(format!("{path}: PM1 metadata has a nontext key; refused"))),
        })
    }).collect()
}

fn read_pickle(path: &str) -> Result<Pm1Regions> {
    let bytes = pickle_reader::read_file_by_magic(Path::new(path))?;
    let value = pickle_reader::load(&bytes)?;
    let Value::Dict(dict) = value else {
        return Err(invalid(format!("{path}: PM1 root is not a dict; refused")));
    };
    let dict = dict.borrow();
    if dict.duplicate_assignments != 0 {
        return Err(invalid(format!(
            "{path}: a PM1 top-level key is assigned twice; refused"
        )));
    }
    if dict.entries.len() != 4 {
        return Err(invalid(format!(
            "{path}: PM1 root must have exactly das_intolerant_domains, hcseeker_hotspots, rmc_constrained and metadata; refused"
        )));
    }
    Ok(Pm1Regions {
        provenance: vec![],
        das_intolerant_domains: plain_set(path, required(path, &dict, "das_intolerant_domains")?)?,
        hcseeker_hotspots: residue_map(
            path,
            "hcseeker_hotspots",
            required(path, &dict, "hcseeker_hotspots")?,
        )?,
        rmc_constrained: residue_map(
            path,
            "rmc_constrained",
            required(path, &dict, "rmc_constrained")?,
        )?,
        metadata: metadata(path, required(path, &dict, "metadata")?)?,
    })
}

fn metadata_plain(value: &ArchivedPm1Metadata) -> (String, ExactValue) {
    match value {
        ArchivedPm1Metadata::Assembly(value) => (
            "assembly".to_owned(),
            ExactValue::Text(value.as_str().to_owned()),
        ),
        ArchivedPm1Metadata::RmcOeThreshold(value) => (
            "rmc_oe_threshold".to_owned(),
            ExactValue::Float(value.to_native().to_bits()),
        ),
        ArchivedPm1Metadata::RmcPBhThreshold(value) => (
            "rmc_p_bh_threshold".to_owned(),
            ExactValue::Float(value.to_native().to_bits()),
        ),
        ArchivedPm1Metadata::NDasIntolerantDomains(value) => (
            "n_das_intolerant_domains".to_owned(),
            ExactValue::Int(value.to_native()),
        ),
        ArchivedPm1Metadata::NHcseekerHotspotTranscripts(value) => (
            "n_hcseeker_hotspot_transcripts".to_owned(),
            ExactValue::Int(value.to_native()),
        ),
        ArchivedPm1Metadata::NHcseekerHotspotResidues(value) => (
            "n_hcseeker_hotspot_residues".to_owned(),
            ExactValue::Int(value.to_native()),
        ),
        ArchivedPm1Metadata::NRmcConstrainedTranscripts(value) => (
            "n_rmc_constrained_transcripts".to_owned(),
            ExactValue::Int(value.to_native()),
        ),
        ArchivedPm1Metadata::NRmcConstrainedResidues(value) => (
            "n_rmc_constrained_residues".to_owned(),
            ExactValue::Int(value.to_native()),
        ),
    }
}

fn plain_metadata(value: &Pm1Metadata) -> (String, ExactValue) {
    match value {
        Pm1Metadata::Assembly(value) => ("assembly".to_owned(), ExactValue::Text(value.to_owned())),
        Pm1Metadata::RmcOeThreshold(value) => (
            "rmc_oe_threshold".to_owned(),
            ExactValue::Float(value.to_bits()),
        ),
        Pm1Metadata::RmcPBhThreshold(value) => (
            "rmc_p_bh_threshold".to_owned(),
            ExactValue::Float(value.to_bits()),
        ),
        Pm1Metadata::NDasIntolerantDomains(value) => (
            "n_das_intolerant_domains".to_owned(),
            ExactValue::Int(*value),
        ),
        Pm1Metadata::NHcseekerHotspotTranscripts(value) => (
            "n_hcseeker_hotspot_transcripts".to_owned(),
            ExactValue::Int(*value),
        ),
        Pm1Metadata::NHcseekerHotspotResidues(value) => (
            "n_hcseeker_hotspot_residues".to_owned(),
            ExactValue::Int(*value),
        ),
        Pm1Metadata::NRmcConstrainedTranscripts(value) => (
            "n_rmc_constrained_transcripts".to_owned(),
            ExactValue::Int(*value),
        ),
        Pm1Metadata::NRmcConstrainedResidues(value) => (
            "n_rmc_constrained_residues".to_owned(),
            ExactValue::Int(*value),
        ),
    }
}

fn map_valid(regions: &[ArchivedTranscriptResidues]) -> bool {
    !regions
        .windows(2)
        .any(|pair| pair[0].transcript.as_str() >= pair[1].transcript.as_str())
        && regions.iter().all(|entry| {
            !entry
                .residues
                .windows(2)
                .any(|pair| pair[0].to_native() >= pair[1].to_native())
        })
}

/// The PM1 regions stay mapped; per-transcript and per-residue membership use binary searches.
pub(crate) struct Pm1RegionCache {
    mapped: MappedCache<Pm1Regions>,
}

// Native runtime consumers can use these borrowed lookup APIs without deserializing.
#[cfg_attr(not(test), allow(dead_code))]
impl Pm1RegionCache {
    pub(crate) fn open(path: &str) -> Result<Self> {
        let mapped = MappedCache::<Pm1Regions>::open(path, PM1_REGIONS_PKL)?;
        let root = mapped.root();
        if root
            .das_intolerant_domains
            .windows(2)
            .any(|pair| pair[0].as_str() >= pair[1].as_str())
            || !map_valid(&root.hcseeker_hotspots)
            || !map_valid(&root.rmc_constrained)
        {
            return Err(invalid(format!(
                "{path}: PM1 lookup keys or residue sets are not sorted and unique"
            )));
        }
        let mut keys = HashSet::new();
        for item in root.metadata.iter() {
            if !keys.insert(metadata_plain(item).0) {
                return Err(invalid(format!("{path}: duplicate PM1 metadata key")));
            }
        }
        Ok(Self { mapped })
    }

    pub(crate) fn hcseeker_residues(&self, transcript: &str) -> Option<&[rkyv::Archived<i64>]> {
        let regions = &self.mapped.root().hcseeker_hotspots;
        regions
            .binary_search_by(|entry| entry.transcript.as_str().cmp(transcript))
            .ok()
            .map(|at| regions[at].residues.as_slice())
    }

    pub(crate) fn rmc_residues(&self, transcript: &str) -> Option<&[rkyv::Archived<i64>]> {
        let regions = &self.mapped.root().rmc_constrained;
        regions
            .binary_search_by(|entry| entry.transcript.as_str().cmp(transcript))
            .ok()
            .map(|at| regions[at].residues.as_slice())
    }

    pub(crate) fn contains_hcseeker(&self, transcript: &str, residue: i64) -> bool {
        self.hcseeker_residues(transcript).is_some_and(|set| {
            set.binary_search_by_key(&residue, |value| value.to_native())
                .is_ok()
        })
    }

    pub(crate) fn contains_rmc(&self, transcript: &str, residue: i64) -> bool {
        self.rmc_residues(transcript).is_some_and(|set| {
            set.binary_search_by_key(&residue, |value| value.to_native())
                .is_ok()
        })
    }

    pub(crate) fn contains_das(&self, domain: &str) -> bool {
        crate::rust_cache::find_sorted(&self.mapped.root().das_intolerant_domains, domain).is_some()
    }

    pub(crate) fn metadata(&self) -> &[ArchivedPm1Metadata] {
        &self.mapped.root().metadata
    }
}

fn check(pickle: &str, output: &Path) -> Result<()> {
    let path = output
        .to_str()
        .ok_or_else(|| invalid("the cache path is not UTF-8"))?;
    let cache = Pm1RegionCache::open(path)?;
    let expected = read_pickle(pickle)?;
    if !matches_expected(&cache, &expected) {
        return Err(invalid(format!(
            "{path}: the converted PM1 values differ from the pickle"
        )));
    }
    Ok(())
}

fn matches_expected(cache: &Pm1RegionCache, expected: &Pm1Regions) -> bool {
    let root = cache.mapped.root();
    let same_map = |actual: &[ArchivedTranscriptResidues], expected: &[TranscriptResidues]| {
        actual.len() == expected.len()
            && actual.iter().zip(expected).all(|(a, b)| {
                a.transcript.as_str() == b.transcript
                    && a.residues
                        .iter()
                        .map(|value| value.to_native())
                        .eq(b.residues.iter().copied())
            })
    };
    root.das_intolerant_domains
        .iter()
        .map(|value| value.as_str())
        .eq(expected.das_intolerant_domains.iter().map(String::as_str))
        && same_map(&root.hcseeker_hotspots, &expected.hcseeker_hotspots)
        && same_map(&root.rmc_constrained, &expected.rmc_constrained)
        && root
            .metadata
            .iter()
            .map(metadata_plain)
            .eq(expected.metadata.iter().map(plain_metadata))
}

/// A native producer writes the frozen schema directly and checks every value
/// through its runtime reader before atomic publication. No pickle is emitted.
pub(crate) fn write_produced(output: &Path, root: &Pm1Regions) -> Result<()> {
    let archive = rkyv::to_bytes::<rancor::Error>(root)?;
    write_cache_checked(output, PM1_REGIONS_PKL, &archive, |staged| {
        let path = staged
            .to_str()
            .ok_or_else(|| invalid("PM1 output path is not UTF-8"))?;
        let cache = Pm1RegionCache::open(path)?;
        if !matches_expected(&cache, root) {
            return Err(invalid(
                "the staged PM1 cache differs from its native producer",
            ));
        }
        Ok(())
    })
}

pub(crate) fn convert(pickle: &str, output: &Path) -> Result<String> {
    ensure_distinct_paths(Path::new(pickle), output)?;
    let start = Instant::now();
    let mut root = read_pickle(pickle)?;
    let das = root.das_intolerant_domains.len();
    let hcseeker = root.hcseeker_hotspots.len();
    let rmc = root.rmc_constrained.len();
    let metadata = root.metadata.len();
    root.provenance = provenance(pickle)?;
    let archive = rkyv::to_bytes::<rancor::Error>(&root)?;
    write_cache_checked(output, PM1_REGIONS_PKL, &archive, |staged| {
        check(pickle, staged)
    })?;
    Ok(format!(
        "{{\"cache\":\"pm1_regions_pkl\",\"das_domains\":{das},\"hcseeker_transcripts\":{hcseeker},\"rmc_transcripts\":{rmc},\"metadata_fields\":{metadata},\"checked_root_fields\":4,\"bytes\":{},\"elapsed_s\":{:.3}}}",
        std::fs::metadata(output)?.len(),
        start.elapsed().as_secs_f64()
    ))
}

pub(crate) fn export(path: &str, out: &mut dyn Write) -> Result<usize> {
    let cache = Pm1RegionCache::open(path)?;
    let root = cache.mapped.root();
    let das = ExactValue::Set(
        root.das_intolerant_domains
            .iter()
            .map(|domain| ExactValue::Text(domain.as_str().to_owned()))
            .collect(),
    );
    let regions = |map: &[ArchivedTranscriptResidues]| {
        ExactValue::Dict(
            map.iter()
                .map(|entry| {
                    (
                        ExactValue::Text(entry.transcript.as_str().to_owned()),
                        ExactValue::Set(
                            entry
                                .residues
                                .iter()
                                .map(|value| ExactValue::Int(value.to_native()))
                                .collect(),
                        ),
                    )
                })
                .collect(),
        )
    };
    let metadata = ExactValue::Dict(
        root.metadata
            .iter()
            .map(|item| {
                let (key, value) = metadata_plain(item);
                (ExactValue::Text(key), value)
            })
            .collect(),
    );
    for (key, value) in [
        ("das_intolerant_domains", das),
        ("hcseeker_hotspots", regions(&root.hcseeker_hotspots)),
        ("metadata", metadata),
        ("rmc_constrained", regions(&root.rmc_constrained)),
    ] {
        write_entry(out, &ExactValue::Text(key.to_owned()), &value)?;
    }
    Ok(4)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Protocol 4, all producer metadata fields, an explicit NaN payload and -0.0.
    // Hotspot transcript order is tx2,tx1; tx1 is an empty set.
    const GOOD: &[u8] = b"\x80\x04\x95\x7e\x01\x00\x00\x00\x00\x00\x00\x7d\x94\x28\x8c\x16\x64\x61\x73\x5f\x69\x6e\x74\x6f\x6c\x65\x72\x61\x6e\x74\x5f\x64\x6f\x6d\x61\x69\x6e\x73\x94\x8f\x94\x28\x8c\x02\x64\x31\x94\x8c\x02\x64\x32\x94\x90\x8c\x11\x68\x63\x73\x65\x65\x6b\x65\x72\x5f\x68\x6f\x74\x73\x70\x6f\x74\x73\x94\x7d\x94\x28\x8c\x03\x74\x78\x32\x94\x8f\x94\x28\x4b\x01\x4b\x03\x4a\xfe\xff\xff\xff\x90\x8c\x03\x74\x78\x31\x94\x8f\x94\x75\x8c\x0f\x72\x6d\x63\x5f\x63\x6f\x6e\x73\x74\x72\x61\x69\x6e\x65\x64\x94\x7d\x94\x8c\x03\x74\x78\x33\x94\x8f\x94\x28\x4b\x04\x4b\x07\x90\x73\x8c\x08\x6d\x65\x74\x61\x64\x61\x74\x61\x94\x7d\x94\x28\x8c\x08\x61\x73\x73\x65\x6d\x62\x6c\x79\x94\x8c\x04\x68\x67\x33\x38\x94\x8c\x10\x72\x6d\x63\x5f\x6f\x65\x5f\x74\x68\x72\x65\x73\x68\x6f\x6c\x64\x94\x47\x7f\xf8\x00\x00\x00\x00\x00\x42\x8c\x12\x72\x6d\x63\x5f\x70\x5f\x62\x68\x5f\x74\x68\x72\x65\x73\x68\x6f\x6c\x64\x94\x47\x80\x00\x00\x00\x00\x00\x00\x00\x8c\x18\x6e\x5f\x64\x61\x73\x5f\x69\x6e\x74\x6f\x6c\x65\x72\x61\x6e\x74\x5f\x64\x6f\x6d\x61\x69\x6e\x73\x94\x4b\x02\x8c\x1e\x6e\x5f\x68\x63\x73\x65\x65\x6b\x65\x72\x5f\x68\x6f\x74\x73\x70\x6f\x74\x5f\x74\x72\x61\x6e\x73\x63\x72\x69\x70\x74\x73\x94\x4b\x02\x8c\x1b\x6e\x5f\x68\x63\x73\x65\x65\x6b\x65\x72\x5f\x68\x6f\x74\x73\x70\x6f\x74\x5f\x72\x65\x73\x69\x64\x75\x65\x73\x94\x4b\x03\x8c\x1d\x6e\x5f\x72\x6d\x63\x5f\x63\x6f\x6e\x73\x74\x72\x61\x69\x6e\x65\x64\x5f\x74\x72\x61\x6e\x73\x63\x72\x69\x70\x74\x73\x94\x4b\x01\x8c\x1a\x6e\x5f\x72\x6d\x63\x5f\x63\x6f\x6e\x73\x74\x72\x61\x69\x6e\x65\x64\x5f\x72\x65\x73\x69\x64\x75\x65\x73\x94\x4b\x02\x75\x75\x2e";

    fn temporary(name: &str, bytes: &[u8]) -> String {
        let path =
            std::env::temp_dir().join(format!("priva_pm1_regions_{}_{name}", std::process::id()));
        std::fs::write(&path, bytes).unwrap();
        path.to_str().unwrap().to_owned()
    }

    #[test]
    fn maps_residue_membership_and_preserves_metadata_float_bits() {
        let pickle = temporary("good.pkl", GOOD);
        let output = temporary("good.rkyv", b"");
        convert(&pickle, Path::new(&output)).unwrap();
        let cache = Pm1RegionCache::open(&output).unwrap();
        assert!(cache.contains_das("d2"));
        assert!(!cache.contains_das("absent"));
        assert!(cache.contains_hcseeker("tx2", -2));
        assert!(cache.contains_hcseeker("tx2", 3));
        assert!(!cache.contains_hcseeker("tx2", 2));
        assert!(!cache.contains_hcseeker("tx1", 1));
        assert!(!cache.contains_hcseeker("absent", 1));
        assert_eq!(cache.hcseeker_residues("tx1").unwrap().len(), 0);
        assert!(cache.contains_rmc("tx3", 7));
        assert!(!cache.contains_rmc("tx2", 3));
        assert_eq!(cache.metadata().len(), 8);
        let metadata: Vec<_> = cache.metadata().iter().map(metadata_plain).collect();
        assert_eq!(metadata[1].1, ExactValue::Float(0x7ff8_0000_0000_0042));
        assert_eq!(metadata[2].1, ExactValue::Float(0x8000_0000_0000_0000));
        let mut exported = Vec::new();
        assert_eq!(export(&output, &mut exported).unwrap(), 4);
        let text = String::from_utf8(exported).unwrap();
        assert!(text.contains("7ff8000000000042"));
        assert!(text.contains("8000000000000000"));
        let again = temporary("again.rkyv", b"");
        convert(&pickle, Path::new(&again)).unwrap();
        assert_eq!(
            std::fs::read(&output).unwrap(),
            std::fs::read(&again).unwrap()
        );
    }

    #[test]
    fn refuses_unknown_metadata_types_and_residue_container_coercions() {
        for bytes in [
            b"\x80\x04}\x8c\x01xNs.".as_slice(),
            b"\x80\x04}\x8c\x10rmc_oe_thresholdK\x01s.",
            b"\x80\x04}\x8c\x18n_das_intolerant_domains\x88s.",
            b"\x80\x04}\x8c\x08assemblyNs.",
        ] {
            let value = pickle_reader::load(bytes).unwrap();
            assert!(metadata("fixture", &value).is_err());
        }
        for bytes in [
            b"\x80\x04}\x8c\x01t]s.".as_slice(),
            b"\x80\x04}\x8c\x01t\x8f(\x88\x90s.",
            b"\x80\x04}\x8c\x01t(K\x01\x91s.",
        ] {
            let value = pickle_reader::load(bytes).unwrap();
            assert!(residue_map("fixture", "hotspots", &value).is_err());
        }
    }

    #[test]
    fn refuses_unknown_fields_and_duplicate_top_assignments_before_publication() {
        let output = temporary("refuse.rkyv", b"unchanged");
        let mut unknown = GOOD.to_vec();
        let at = unknown
            .windows(8)
            .position(|part| part == b"assembly")
            .unwrap();
        unknown[at..at + 8].copy_from_slice(b"unknownx");
        let mut duplicate = GOOD[..GOOD.len() - 1].to_vec();
        duplicate.extend_from_slice(b"\x8c\x08metadata}s.");
        for (at, bytes) in [unknown, duplicate].iter().enumerate() {
            let pickle = temporary(&format!("bad{at}.pkl"), bytes);
            assert!(convert(&pickle, Path::new(&output)).is_err());
            assert_eq!(std::fs::read(&output).unwrap(), b"unchanged");
        }
    }

    #[test]
    fn reader_rejects_unsorted_residues_and_duplicate_metadata_keys() {
        let output = temporary("bad_schema.rkyv", b"");
        let mut root = Pm1Regions {
            provenance: vec![],
            das_intolerant_domains: vec![],
            hcseeker_hotspots: vec![TranscriptResidues {
                transcript: "t".into(),
                residues: vec![2, 1],
            }],
            rmc_constrained: vec![],
            metadata: vec![],
        };
        let archive = rkyv::to_bytes::<rancor::Error>(&root).unwrap();
        crate::rust_cache::write_cache(Path::new(&output), PM1_REGIONS_PKL, &archive).unwrap();
        assert!(Pm1RegionCache::open(&output).is_err());
        root.hcseeker_hotspots.clear();
        root.metadata = vec![
            Pm1Metadata::Assembly("a".into()),
            Pm1Metadata::Assembly("b".into()),
        ];
        let archive = rkyv::to_bytes::<rancor::Error>(&root).unwrap();
        crate::rust_cache::write_cache(Path::new(&output), PM1_REGIONS_PKL, &archive).unwrap();
        assert!(Pm1RegionCache::open(&output).is_err());
    }
}
