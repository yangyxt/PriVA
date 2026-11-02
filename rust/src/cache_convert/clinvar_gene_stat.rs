//! Schema 5, version 1: `clinvar_stat_variants.per_gene_stat` and its complete variant sets.
//!
//! Python: `{gene: {gene_symbol: str, ensembl_id: str, <seven named sets of str>,
//! domains: {domain: <seven named sets>}, transcripts: {transcript: <seven named sets>}}}`.
//! The eleven gene-record fields and seven set-record fields are exact; absent/extra fields,
//! non-text IDs and set/list/tuple substitutions are refused. The seven collections are Python
//! `set`, not `frozenset`. Domain/transcript dict order is retained. Only gene keys and set
//! members are sorted. A sorted shared variant-ID pool makes sorted `u32` IDs also text-sorted.
use super::Provenance;
use super::common::{ExactValue, provenance, write_entry};
use crate::cli::{Result, invalid};
use crate::pickle_reader::{self, Value};
use crate::rust_cache::{
    MappedCache, Schema, ensure_distinct_paths, find_sorted, write_cache_checked,
};
use rkyv::{Archive, Deserialize, Serialize, rancor};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::io::Write;
use std::path::Path;
use std::time::Instant;

pub(crate) const SCHEMA: Schema = Schema {
    name: "clinvar_gene_stat",
    id: 5,
    version: 1,
};

/// Stable order of the seven fixed record fields, used by plain values and exact export.
pub(crate) const SET_NAMES: [&str; 7] = [
    "clinvar_variants",
    "clinvar_pathogenic",
    "clinvar_not_pathogenic",
    "putative_nmd_variants",
    "small_aachange",
    "large_aachange",
    "no_aachange",
];

#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct VariantSets {
    pub(crate) clinvar_variants: Vec<u32>,
    pub(crate) clinvar_pathogenic: Vec<u32>,
    pub(crate) clinvar_not_pathogenic: Vec<u32>,
    pub(crate) putative_nmd_variants: Vec<u32>,
    pub(crate) small_aachange: Vec<u32>,
    pub(crate) large_aachange: Vec<u32>,
    pub(crate) no_aachange: Vec<u32>,
}

impl ArchivedVariantSets {
    /// All seven sets, in `SET_NAMES` order; each is sorted and unique.
    pub(crate) fn members(&self) -> [&[rkyv::Archived<u32>]; 7] {
        [
            &self.clinvar_variants,
            &self.clinvar_pathogenic,
            &self.clinvar_not_pathogenic,
            &self.putative_nmd_variants,
            &self.small_aachange,
            &self.large_aachange,
            &self.no_aachange,
        ]
    }
}

#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct Domain {
    pub(crate) domain: String,
    pub(crate) sets: VariantSets,
}

#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct Transcript {
    pub(crate) transcript: String,
    pub(crate) sets: VariantSets,
}

#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct Gene {
    pub(crate) ensembl_id: String,
    pub(crate) gene_symbol: String,
    pub(crate) sets: VariantSets,
    /// Python free-key dicts, in insertion order.
    pub(crate) domains: Vec<Domain>,
    pub(crate) transcripts: Vec<Transcript>,
}

#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct ClinvarGeneStat {
    pub(crate) provenance: Vec<Provenance>,
    /// Sorted and unique. Every set member is an index into this pool.
    pub(crate) variants: Vec<String>,
    /// Sorted and unique, parallel to `entries`.
    pub(crate) genes: Vec<String>,
    pub(crate) entries: Vec<Gene>,
}

/// Seven Python sets represented as sorted text vectors in `SET_NAMES` order.
pub(crate) type PlainSets = [Vec<String>; 7];

#[derive(Debug, PartialEq)]
pub(crate) struct PlainGene {
    pub(crate) ensembl_id: String,
    pub(crate) gene_symbol: String,
    pub(crate) sets: PlainSets,
    pub(crate) domains: Vec<(String, PlainSets)>,
    pub(crate) transcripts: Vec<(String, PlainSets)>,
}

pub(crate) type PlainStats = Vec<(String, PlainGene)>;

fn text(path: &str, value: &Value) -> Result<String> {
    match value {
        Value::Str(text) => Ok(text.to_string()),
        other => Err(invalid(format!(
            "{path}: a ClinVar gene text is a {}; refused",
            other.type_name()
        ))),
    }
}

fn check_fields(path: &str, entries: &[(Value, Value)], expected: &[&str]) -> Result<()> {
    if entries.len() != expected.len()
        || expected.iter().any(|name| {
            entries
                .iter()
                .filter(|(key, _)| matches!(key, Value::Str(key) if key.as_ref() == *name))
                .count()
                != 1
        })
    {
        return Err(invalid(format!(
            "{path}: a ClinVar gene record has missing or extra fields; refused"
        )));
    }
    Ok(())
}

fn field<'a>(entries: &'a [(Value, Value)], name: &str) -> &'a Value {
    // `check_fields` proved the complete exact key set first.
    &entries
        .iter()
        .find(|(key, _)| matches!(key, Value::Str(text) if text.as_ref() == name))
        .unwrap()
        .1
}

fn plain_set(path: &str, value: &Value) -> Result<Vec<String>> {
    let Value::Set(set) = value else {
        return Err(invalid(format!(
            "{path}: a ClinVar variant collection is a {}, not a set; refused",
            value.type_name()
        )));
    };
    let mut members = set
        .borrow()
        .items
        .iter()
        .map(|item| text(path, item))
        .collect::<Result<Vec<_>>>()?;
    members.sort();
    if members.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(invalid(format!(
            "{path}: a ClinVar variant set has duplicate members; refused"
        )));
    }
    Ok(members)
}

fn plain_sets(path: &str, entries: &[(Value, Value)]) -> Result<PlainSets> {
    let values = SET_NAMES
        .iter()
        .map(|name| plain_set(path, field(entries, name)))
        .collect::<Result<Vec<_>>>()?;
    Ok(values.try_into().unwrap())
}

fn plain_nested(path: &str, value: &Value) -> Result<Vec<(String, PlainSets)>> {
    let Value::Dict(dict) = value else {
        return Err(invalid(format!(
            "{path}: ClinVar domains/transcripts are not a dict; refused"
        )));
    };
    let dict = dict.borrow();
    if dict.duplicate_assignments != 0 {
        return Err(invalid(format!(
            "{path}: a ClinVar domain/transcript key is set twice; refused"
        )));
    }
    let mut result = Vec::with_capacity(dict.entries.len());
    for (key, value) in &dict.entries {
        let Value::Dict(sets) = value else {
            return Err(invalid(format!(
                "{path}: ClinVar domain/transcript sets are not a dict; refused"
            )));
        };
        let sets = sets.borrow();
        check_fields(path, &sets.entries, &SET_NAMES)?;
        if sets.duplicate_assignments != 0 {
            return Err(invalid(format!(
                "{path}: a ClinVar set field is set twice; refused"
            )));
        }
        result.push((text(path, key)?, plain_sets(path, &sets.entries)?));
    }
    Ok(result)
}

fn plain_gene(path: &str, value: &Value) -> Result<PlainGene> {
    let Value::Dict(dict) = value else {
        return Err(invalid(format!(
            "{path}: a ClinVar gene value is not a dict; refused"
        )));
    };
    let dict = dict.borrow();
    let mut names = SET_NAMES.to_vec();
    names.extend_from_slice(&["ensembl_id", "gene_symbol", "domains", "transcripts"]);
    check_fields(path, &dict.entries, &names)?;
    if dict.duplicate_assignments != 0 {
        return Err(invalid(format!(
            "{path}: a ClinVar gene record field is set twice; refused"
        )));
    }
    Ok(PlainGene {
        ensembl_id: text(path, field(&dict.entries, "ensembl_id"))?,
        gene_symbol: text(path, field(&dict.entries, "gene_symbol"))?,
        sets: plain_sets(path, &dict.entries)?,
        domains: plain_nested(path, field(&dict.entries, "domains"))?,
        transcripts: plain_nested(path, field(&dict.entries, "transcripts"))?,
    })
}

fn plain_stats(path: &str, value: &Value) -> Result<PlainStats> {
    let Value::Dict(dict) = value else {
        return Err(invalid(format!(
            "{path}: the ClinVar gene statistics are not a dict; refused"
        )));
    };
    let dict = dict.borrow();
    if dict.duplicate_assignments != 0 {
        return Err(invalid(format!(
            "{path}: a ClinVar gene key is set twice; refused"
        )));
    }
    dict.entries
        .iter()
        .map(|(key, value)| Ok((text(path, key)?, plain_gene(path, value)?)))
        .collect()
}

fn read_pickle(path: &str) -> Result<PlainStats> {
    let bytes = pickle_reader::read_file(Path::new(path), path.ends_with(".gz"))?;
    plain_stats(path, &pickle_reader::load(&bytes)?)
}

fn encode_sets(sets: PlainSets, ids: &HashMap<&str, u32>) -> VariantSets {
    let encode = |set: Vec<String>| set.iter().map(|variant| ids[variant.as_str()]).collect();
    let [
        variants,
        pathogenic,
        not_pathogenic,
        nmd,
        small,
        large,
        none,
    ] = sets;
    VariantSets {
        clinvar_variants: encode(variants),
        clinvar_pathogenic: encode(pathogenic),
        clinvar_not_pathogenic: encode(not_pathogenic),
        putative_nmd_variants: encode(nmd),
        small_aachange: encode(small),
        large_aachange: encode(large),
        no_aachange: encode(none),
    }
}

fn build(mut entries: PlainStats, provenance: Vec<Provenance>) -> Result<Vec<u8>> {
    entries.sort_by(|left, right| left.0.cmp(&right.0));
    if entries.windows(2).any(|pair| pair[0].0 == pair[1].0) {
        return Err(invalid("the ClinVar gene keys are not unique"));
    }
    let mut variants = BTreeSet::new();
    for (_, gene) in &entries {
        for sets in std::iter::once(&gene.sets)
            .chain(gene.domains.iter().map(|(_, sets)| sets))
            .chain(gene.transcripts.iter().map(|(_, sets)| sets))
        {
            for set in sets {
                if set.windows(2).any(|pair| pair[0] >= pair[1]) {
                    return Err(invalid("a plain ClinVar set is not sorted and unique"));
                }
                variants.extend(set.iter().cloned());
            }
        }
    }
    let variants: Vec<String> = variants.into_iter().collect();
    if variants.len() > u32::MAX as usize {
        return Err(invalid("too many ClinVar variants for the u32 pool"));
    }
    let ids: HashMap<&str, u32> = variants
        .iter()
        .enumerate()
        .map(|(at, variant)| (variant.as_str(), at as u32))
        .collect();
    let mut genes = Vec::with_capacity(entries.len());
    let mut archived = Vec::with_capacity(entries.len());
    for (key, gene) in entries {
        genes.push(key);
        archived.push(Gene {
            ensembl_id: gene.ensembl_id,
            gene_symbol: gene.gene_symbol,
            sets: encode_sets(gene.sets, &ids),
            domains: gene
                .domains
                .into_iter()
                .map(|(domain, sets)| Domain {
                    domain,
                    sets: encode_sets(sets, &ids),
                })
                .collect(),
            transcripts: gene
                .transcripts
                .into_iter()
                .map(|(transcript, sets)| Transcript {
                    transcript,
                    sets: encode_sets(sets, &ids),
                })
                .collect(),
        });
    }
    drop(ids);
    Ok(rkyv::to_bytes::<rancor::Error>(&ClinvarGeneStat {
        provenance,
        variants,
        genes,
        entries: archived,
    })?
    .to_vec())
}

/// Checked mmap view. The complete sets remain archived; no whole-map materialization is needed.
pub(crate) struct ClinvarGeneCache {
    mapped: MappedCache<ClinvarGeneStat>,
}

impl ClinvarGeneCache {
    pub(crate) fn open(path: &str) -> Result<Self> {
        let mapped = MappedCache::<ClinvarGeneStat>::open(path, SCHEMA)?;
        let root = mapped.root();
        let ordered = |keys: &[rkyv::string::ArchivedString]| {
            keys.windows(2)
                .all(|pair| pair[0].as_str() < pair[1].as_str())
        };
        if root.genes.len() != root.entries.len()
            || !ordered(&root.genes)
            || !ordered(&root.variants)
        {
            return Err(invalid(format!(
                "{path}: ClinVar genes/variants are not sorted and unique, or not one entry per gene"
            )));
        }
        let count = root.variants.len();
        let valid_sets = |sets: &ArchivedVariantSets| {
            sets.members().iter().all(|ids| {
                ids.iter().all(|id| (id.to_native() as usize) < count)
                    && ids
                        .windows(2)
                        .all(|pair| pair[0].to_native() < pair[1].to_native())
            })
        };
        for entry in root.entries.iter() {
            let mut domains = HashSet::with_capacity(entry.domains.len());
            let mut transcripts = HashSet::with_capacity(entry.transcripts.len());
            if !valid_sets(&entry.sets)
                || !entry.domains.iter().all(|domain| {
                    domains.insert(domain.domain.as_str()) && valid_sets(&domain.sets)
                })
                || !entry.transcripts.iter().all(|transcript| {
                    transcripts.insert(transcript.transcript.as_str())
                        && valid_sets(&transcript.sets)
                })
            {
                return Err(invalid(format!(
                    "{path}: ClinVar sets have an invalid index/order/duplicate, or an inner dict has duplicate keys"
                )));
            }
        }
        Ok(Self { mapped })
    }

    pub(crate) fn len(&self) -> usize {
        self.mapped.root().genes.len()
    }

    pub(crate) fn keys(&self) -> impl Iterator<Item = &str> {
        self.mapped.root().genes.iter().map(|key| key.as_str())
    }

    pub(crate) fn get(&self, gene: &str) -> Option<&ArchivedGene> {
        let root = self.mapped.root();
        find_sorted(&root.genes, gene).map(|at| &root.entries[at])
    }

    pub(crate) fn variant(&self, id: &rkyv::Archived<u32>) -> &str {
        self.mapped.root().variants[id.to_native() as usize].as_str()
    }

    fn plain_sets(&self, sets: &ArchivedVariantSets) -> PlainSets {
        sets.members()
            .map(|ids| ids.iter().map(|id| self.variant(id).to_owned()).collect())
    }

    pub(crate) fn plain(&self, gene: &str) -> Option<PlainGene> {
        self.get(gene).map(|entry| PlainGene {
            ensembl_id: entry.ensembl_id.as_str().to_owned(),
            gene_symbol: entry.gene_symbol.as_str().to_owned(),
            sets: self.plain_sets(&entry.sets),
            domains: entry
                .domains
                .iter()
                .map(|domain| {
                    (
                        domain.domain.as_str().to_owned(),
                        self.plain_sets(&domain.sets),
                    )
                })
                .collect(),
            transcripts: entry
                .transcripts
                .iter()
                .map(|transcript| {
                    (
                        transcript.transcript.as_str().to_owned(),
                        self.plain_sets(&transcript.sets),
                    )
                })
                .collect(),
        })
    }
}

fn exact_sets(sets: PlainSets) -> Vec<(ExactValue, ExactValue)> {
    SET_NAMES
        .into_iter()
        .zip(sets)
        .map(|(name, members)| {
            (
                ExactValue::Text(name.to_owned()),
                ExactValue::Set(members.into_iter().map(ExactValue::Text).collect()),
            )
        })
        .collect()
}

fn exact(gene: PlainGene) -> ExactValue {
    let mut fields = vec![
        (
            ExactValue::Text("ensembl_id".to_owned()),
            ExactValue::Text(gene.ensembl_id),
        ),
        (
            ExactValue::Text("gene_symbol".to_owned()),
            ExactValue::Text(gene.gene_symbol),
        ),
    ];
    fields.extend(exact_sets(gene.sets));
    for (name, entries) in [("domains", gene.domains), ("transcripts", gene.transcripts)] {
        fields.push((
            ExactValue::Text(name.to_owned()),
            ExactValue::Dict(
                entries
                    .into_iter()
                    .map(|(key, sets)| (ExactValue::Text(key), ExactValue::Dict(exact_sets(sets))))
                    .collect(),
            ),
        ));
    }
    ExactValue::Dict(fields)
}

pub(crate) fn export(path: &str, out: &mut dyn Write) -> Result<usize> {
    let cache = ClinvarGeneCache::open(path)?;
    for key in cache.keys() {
        write_entry(
            out,
            &ExactValue::Text(key.to_owned()),
            &exact(cache.plain(key).unwrap()),
        )?;
    }
    Ok(cache.len())
}

pub(crate) fn convert(pickle: &str, output: &Path) -> Result<String> {
    ensure_distinct_paths(Path::new(pickle), output)?;
    let start = Instant::now();
    let entries = read_pickle(pickle)?;
    let count = entries.len();
    let domains: usize = entries.iter().map(|(_, gene)| gene.domains.len()).sum();
    let transcripts: usize = entries.iter().map(|(_, gene)| gene.transcripts.len()).sum();
    let memberships: usize = entries
        .iter()
        .flat_map(|(_, gene)| {
            std::iter::once(&gene.sets)
                .chain(gene.domains.iter().map(|(_, sets)| sets))
                .chain(gene.transcripts.iter().map(|(_, sets)| sets))
        })
        .flat_map(|sets| sets.iter())
        .map(Vec::len)
        .sum();
    let archive = build(entries, provenance(pickle)?)?;
    let converted_s = start.elapsed().as_secs_f64();
    let check = Instant::now();
    write_cache_checked(output, SCHEMA, &archive, |staged| {
        let cache = ClinvarGeneCache::open(&staged.to_string_lossy())?;
        let expected = read_pickle(pickle)?;
        if cache.len() != count || expected.len() != count {
            return Err(invalid(
                "the ClinVar gene statistics have a different key count",
            ));
        }
        for (gene, value) in expected {
            if cache.plain(&gene) != Some(value) {
                return Err(invalid("a ClinVar gene record differs from the pickle"));
            }
        }
        Ok(())
    })?;
    Ok(format!(
        "{{\"cache\":\"clinvar_gene_stat\",\"genes\":{count},\"domains\":{domains},\"transcripts\":{transcripts},\"memberships\":{memberships},\"bytes\":{},\"converted_s\":{converted_s:.3},\"checked_keys\":{count},\"check_s\":{:.3}}}",
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
    fn dict(entries: Vec<(Value, Value)>) -> Value {
        let mut value = crate::pickle_reader::Dict::default();
        value.entries = entries;
        Value::Dict(Rc::new(RefCell::new(value)))
    }
    fn sets() -> Value {
        dict(
            SET_NAMES
                .iter()
                .map(|name| (text(name), Value::Set(Default::default())))
                .collect(),
        )
    }
    fn gene() -> Value {
        let Value::Dict(sets) = sets() else {
            unreachable!()
        };
        let mut entries = sets.borrow().entries.clone();
        entries.extend([
            (text("ensembl_id"), text("ENSG1")),
            (text("gene_symbol"), text("Gene\u{00e9}")),
            (
                text("domains"),
                dict(vec![(text("z"), self::sets()), (text("a"), self::sets())]),
            ),
            (
                text("transcripts"),
                dict(vec![
                    (text("ENST2"), self::sets()),
                    (text("ENST1"), self::sets()),
                ]),
            ),
        ]);
        dict(entries)
    }
    fn temporary(name: &str) -> String {
        std::env::temp_dir()
            .join(format!("priva_gene_cache_{}_{name}", std::process::id()))
            .to_str()
            .unwrap()
            .to_owned()
    }

    #[test]
    fn keeps_free_dict_order_and_all_seven_sets() {
        let mut gene = plain_gene("fixture", &gene()).unwrap();
        for index in 0..7 {
            gene.sets[index] = vec![format!("v{index}"), "z".into()];
        }
        gene.domains[0].1[0] = vec!["domain_variant".into()];
        gene.transcripts[1].1[6] = vec!["tranx_variant".into()];
        let archive = build(vec![("ENSG1".into(), gene)], vec![]).unwrap();
        let path = temporary("good.rkyv");
        crate::rust_cache::write_cache(Path::new(&path), SCHEMA, &archive).unwrap();
        let cache = ClinvarGeneCache::open(&path).unwrap();
        let actual = cache.plain("ENSG1").unwrap();
        assert_eq!(
            actual
                .domains
                .iter()
                .map(|(key, _)| key.as_str())
                .collect::<Vec<_>>(),
            ["z", "a"]
        );
        assert_eq!(
            actual
                .transcripts
                .iter()
                .map(|(key, _)| key.as_str())
                .collect::<Vec<_>>(),
            ["ENST2", "ENST1"]
        );
        for index in 0..7 {
            assert_eq!(actual.sets[index], [format!("v{index}"), "z".into()]);
        }
        assert_eq!(actual.domains[0].1[0], ["domain_variant"]);
        assert_eq!(actual.transcripts[1].1[6], ["tranx_variant"]);
        assert!(cache.get("missing").is_none());
        let mut export_bytes = Vec::new();
        assert_eq!(export(&path, &mut export_bytes).unwrap(), 1);
    }

    #[test]
    fn refuses_missing_extra_and_wrongly_typed_fields() {
        for value in [Value::None, Value::Int(1), Value::List(Default::default())] {
            assert!(plain_gene("f", &value).is_err());
        }
        let Value::Dict(good) = gene() else {
            unreachable!()
        };
        let good = good.borrow();
        let mut missing = good.entries.clone();
        missing.pop();
        assert!(plain_gene("f", &dict(missing)).is_err());
        let mut extra = good.entries.clone();
        extra.push((text("extra"), Value::None));
        assert!(plain_gene("f", &dict(extra)).is_err());
        for replacement in [
            Value::List(Default::default()),
            Value::Tuple(Default::default()),
            Value::FrozenSet(Default::default()),
            Value::None,
        ] {
            let mut entries = good.entries.clone();
            entries[0].1 = replacement;
            assert!(plain_gene("f", &dict(entries)).is_err());
        }
        let bad_set = Value::Set(Default::default());
        let Value::Set(set) = &bad_set else {
            unreachable!()
        };
        set.borrow_mut().items.push(Value::Int(1));
        assert!(plain_set("f", &bad_set).is_err());
        let duplicate = pickle_reader::load(b"\x80\x04}(\x8c\x01aN\x8c\x01aNu.").unwrap();
        assert!(
            plain_stats("f", &duplicate)
                .unwrap_err()
                .to_string()
                .contains("twice")
        );
    }

    fn blank_sets() -> VariantSets {
        VariantSets {
            clinvar_variants: vec![],
            clinvar_pathogenic: vec![],
            clinvar_not_pathogenic: vec![],
            putative_nmd_variants: vec![],
            small_aachange: vec![],
            large_aachange: vec![],
            no_aachange: vec![],
        }
    }
    fn raw_gene() -> Gene {
        Gene {
            ensembl_id: "id".into(),
            gene_symbol: "symbol".into(),
            sets: blank_sets(),
            domains: vec![],
            transcripts: vec![],
        }
    }
    fn open_raw(name: &str, value: &ClinvarGeneStat) -> Result<ClinvarGeneCache> {
        let path = temporary(name);
        crate::rust_cache::write_cache(
            Path::new(&path),
            SCHEMA,
            &rkyv::to_bytes::<rancor::Error>(value)?,
        )?;
        ClinvarGeneCache::open(&path)
    }

    fn mutable_sets(sets: &mut VariantSets) -> [&mut Vec<u32>; 7] {
        [
            &mut sets.clinvar_variants,
            &mut sets.clinvar_pathogenic,
            &mut sets.clinvar_not_pathogenic,
            &mut sets.putative_nmd_variants,
            &mut sets.small_aachange,
            &mut sets.large_aachange,
            &mut sets.no_aachange,
        ]
    }

    #[test]
    fn checks_every_set_in_gene_domain_and_transcript_scopes() {
        for scope in 0..3 {
            for set_at in 0..7 {
                let mut entry = raw_gene();
                let sets = match scope {
                    0 => &mut entry.sets,
                    1 => {
                        entry.domains.push(Domain {
                            domain: "d".into(),
                            sets: blank_sets(),
                        });
                        &mut entry.domains[0].sets
                    }
                    _ => {
                        entry.transcripts.push(Transcript {
                            transcript: "t".into(),
                            sets: blank_sets(),
                        });
                        &mut entry.transcripts[0].sets
                    }
                };
                mutable_sets(sets)[set_at].push(1);
                let root = ClinvarGeneStat {
                    provenance: vec![],
                    variants: vec!["v".into()],
                    genes: vec!["g".into()],
                    entries: vec![entry],
                };
                assert!(open_raw(&format!("bad_{scope}_{set_at}"), &root).is_err());
            }
        }
    }

    #[test]
    fn validates_all_key_vectors_nested_duplicates_and_set_indices() {
        let mut root = ClinvarGeneStat {
            provenance: vec![],
            variants: vec!["a".into(), "z".into()],
            genes: vec![],
            entries: vec![],
        };
        assert_eq!(open_raw("empty", &root).unwrap().len(), 0);
        root.genes = vec!["z".into(), "a".into()];
        root.entries = vec![raw_gene(), raw_gene()];
        assert!(open_raw("unsorted", &root).is_err());
        root.genes = vec!["a".into(), "a".into()];
        assert!(open_raw("repeated", &root).is_err());
        root.genes.pop();
        assert!(open_raw("unequal", &root).is_err());
        root.entries.pop();
        root.variants.reverse();
        assert!(open_raw("unsorted_pool", &root).is_err());
        root.variants = vec!["a".into(), "a".into()];
        assert!(open_raw("repeated_pool", &root).is_err());
        root.variants = vec!["a".into(), "z".into()];
        root.entries[0].sets.clinvar_variants = vec![2];
        assert!(open_raw("out_of_range", &root).is_err());
        root.entries[0].sets.clinvar_variants = vec![1, 0];
        assert!(open_raw("unsorted_set", &root).is_err());
        root.entries[0].sets.clinvar_variants = vec![0, 0];
        assert!(open_raw("repeated_set", &root).is_err());
        root.entries[0].sets.clinvar_variants.clear();
        root.entries[0].domains = vec![
            Domain {
                domain: "same".into(),
                sets: blank_sets(),
            },
            Domain {
                domain: "same".into(),
                sets: blank_sets(),
            },
        ];
        assert!(open_raw("repeated_domain", &root).is_err());
        root.entries[0].domains.clear();
        root.entries[0].transcripts = vec![
            Transcript {
                transcript: "same".into(),
                sets: blank_sets(),
            },
            Transcript {
                transcript: "same".into(),
                sets: blank_sets(),
            },
        ];
        assert!(open_raw("repeated_transcript", &root).is_err());
    }
}
