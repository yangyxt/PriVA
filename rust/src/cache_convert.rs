//! `priva cache-convert`: one-time conversion of an installed Python pickle cache to the branch's
//! Rust-native format (`rust/CACHE_FORMATS.md`), with a check that the converted file answers
//! every key of the pickle with the same value. Prototype: `clinvar_aa_stat` only.
//!
//! `clinvar_aa_stat` (schema 1, version 1) holds the ClinVar amino-acid dictionary
//! `{transcript: {protein_pos: {hgvsp: {"CLNSIG": [str], "CLNREVSTAT": [str]}}}}`:
//! transcripts sorted by their bytes (Python's top-level order is not kept; no reader uses it),
//! everything below in the pickle's dict and list order, and the CLNSIG and CLNREVSTAT texts
//! stored once each in `texts` and referenced by index. Any other shape is refused: a tuple
//! where a list stands, other keys than CLNSIG and CLNREVSTAT, text that is not text.
use crate::acmg_step::pickle_dict_stream::{Item, visit_top_level_dict};
use crate::cli::{Result, invalid};
use crate::rust_cache::{
    MappedCache, Schema, ensure_distinct_paths, find_sorted, write_cache_checked,
};
use rkyv::{Archive, Deserialize, Serialize, rancor};
use std::collections::HashMap;
use std::path::Path;
use std::time::Instant;

pub(crate) mod clinvar_splice_stat;
pub(crate) mod common;
pub(crate) mod intolerant_domains;
pub(crate) mod pm1_regions;
pub(crate) mod transcript_domains;

pub(crate) const CLINVAR_AA_STAT: Schema = Schema {
    name: "clinvar_aa_stat",
    id: 1,
    version: 1,
};

/// One `key: value` line of where the cache came from.
#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct Provenance {
    pub(crate) key: String,
    pub(crate) value: String,
}

/// `{"CLNSIG": [...], "CLNREVSTAT": [...]}` of one HGVSp, as indexes into `texts`.
#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct AaChange {
    pub(crate) hgvsp: String,
    pub(crate) clnsig: Vec<u32>,
    pub(crate) clnrevstat: Vec<u32>,
}

#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct AaPosition {
    pub(crate) protein_pos: String,
    pub(crate) changes: Vec<AaChange>,
}

#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct AaTranscript {
    pub(crate) positions: Vec<AaPosition>,
}

#[derive(Archive, Serialize, Deserialize, Debug, PartialEq)]
pub(crate) struct ClinvarAaStat {
    pub(crate) provenance: Vec<Provenance>,
    pub(crate) texts: Vec<String>,
    /// Sorted by bytes, unique; `entries[i]` belongs to `transcripts[i]`.
    pub(crate) transcripts: Vec<String>,
    pub(crate) entries: Vec<AaTranscript>,
}

/// One transcript as plain values, the unit the converter and its check compare.
type PlainChange = (String, Vec<String>, Vec<String>);
type PlainTranscript = Vec<(String, Vec<PlainChange>)>;

fn text_list(file: &str, transcript: &str, item: Option<&Item>, name: &str) -> Result<Vec<String>> {
    match item {
        Some(Item::List(items)) => items
            .borrow()
            .iter()
            .map(|item| match item {
                Item::Str(text) => Ok(text.to_string()),
                other => Err(invalid(format!(
                    "{file}: transcript {transcript}: an item of {name} is a {}; refused",
                    other.type_name()
                ))),
            })
            .collect(),
        Some(other) => Err(invalid(format!(
            "{file}: transcript {transcript}: {name} is a {}, not a list; refused",
            other.type_name()
        ))),
        None => Err(invalid(format!(
            "{file}: transcript {transcript}: an entry lacks {name}; refused"
        ))),
    }
}

/// Decodes one transcript of the pickle, refusing every shape but the installed one.
fn plain_transcript(file: &str, transcript: &str, value: &Item) -> Result<PlainTranscript> {
    let odd = |what: &str| invalid(format!("{file}: transcript {transcript}: {what}; refused"));
    let Item::Dict(positions) = value else {
        return Err(odd(&format!(
            "the value is a {}, not a dict",
            value.type_name()
        )));
    };
    let mut plain = Vec::new();
    for (position, changes) in &positions.borrow().entries {
        let Item::Dict(changes) = changes else {
            return Err(odd("a protein position maps to a non-dict"));
        };
        let mut list = Vec::new();
        for (hgvsp, entry) in &changes.borrow().entries {
            let Item::Dict(entry) = entry else {
                return Err(odd("an HGVSp maps to a non-dict"));
            };
            let entry = entry.borrow();
            if entry.len() != 2 {
                return Err(odd("an entry holds other keys than CLNSIG and CLNREVSTAT"));
            }
            list.push((
                hgvsp.to_string(),
                text_list(file, transcript, entry.get("CLNSIG"), "CLNSIG")?,
                text_list(file, transcript, entry.get("CLNREVSTAT"), "CLNREVSTAT")?,
            ));
        }
        plain.push((position.to_string(), list));
    }
    Ok(plain)
}

fn plain_from_archive(entry: &ArchivedAaTranscript, texts: &[String]) -> PlainTranscript {
    let text = |id: &rkyv::Archived<u32>| texts[id.to_native() as usize].clone();
    entry
        .positions
        .iter()
        .map(|position| {
            (
                position.protein_pos.as_str().to_owned(),
                position
                    .changes
                    .iter()
                    .map(|change| {
                        (
                            change.hgvsp.as_str().to_owned(),
                            change.clnsig.iter().map(text).collect(),
                            change.clnrevstat.iter().map(text).collect(),
                        )
                    })
                    .collect(),
            )
        })
        .collect()
}

/// Every transcript of the pickle, decoded, in pickle order; a key set twice is refused.
fn read_pickle(path: &str) -> Result<Vec<(String, PlainTranscript)>> {
    let mut entries = Vec::new();
    let distinct = visit_top_level_dict(
        path,
        |_| true,
        |transcript, value| {
            entries.push((
                transcript.to_string(),
                plain_transcript(path, &transcript, &value)?,
            ));
            Ok(())
        },
    )?;
    if distinct != entries.len() {
        return Err(invalid(format!(
            "{path}: a transcript is set twice in the pickle; refused"
        )));
    }
    Ok(entries)
}

/// The archive bytes for these transcripts (any order; sorted here).
fn build(
    mut entries: Vec<(String, PlainTranscript)>,
    provenance: Vec<Provenance>,
) -> Result<Vec<u8>> {
    entries.sort_by(|a, b| a.0.cmp(&b.0));
    let mut texts: Vec<String> = Vec::new();
    let mut ids: HashMap<String, u32> = HashMap::new();
    let mut id = |text: String| -> u32 {
        if let Some(&id) = ids.get(&text) {
            return id;
        }
        let id = texts.len() as u32;
        texts.push(text.clone());
        ids.insert(text, id);
        id
    };
    let mut transcripts = Vec::with_capacity(entries.len());
    let mut archived = Vec::with_capacity(entries.len());
    for (transcript, positions) in entries {
        transcripts.push(transcript);
        archived.push(AaTranscript {
            positions: positions
                .into_iter()
                .map(|(protein_pos, changes)| AaPosition {
                    protein_pos,
                    changes: changes
                        .into_iter()
                        .map(|(hgvsp, clnsig, clnrevstat)| AaChange {
                            hgvsp,
                            clnsig: clnsig.into_iter().map(&mut id).collect(),
                            clnrevstat: clnrevstat.into_iter().map(&mut id).collect(),
                        })
                        .collect(),
                })
                .collect(),
        });
    }
    let value = ClinvarAaStat {
        provenance,
        texts,
        transcripts,
        entries: archived,
    };
    Ok(rkyv::to_bytes::<rancor::Error>(&value)?.to_vec())
}

/// The converted ClinVar amino-acid cache, mapped and validated.
pub(crate) struct ClinvarAaCache {
    mapped: MappedCache<ClinvarAaStat>,
    texts: Vec<String>,
}

impl ClinvarAaCache {
    pub(crate) fn open(path: &str) -> Result<ClinvarAaCache> {
        let mapped = MappedCache::<ClinvarAaStat>::open(path, CLINVAR_AA_STAT)?;
        let root = mapped.root();
        if root.transcripts.len() != root.entries.len()
            || root
                .transcripts
                .windows(2)
                .any(|pair| pair[0].as_str() >= pair[1].as_str())
        {
            return Err(invalid(format!(
                "{path}: the transcripts are not sorted and unique, or not one per entry"
            )));
        }
        let texts: Vec<String> = root.texts.iter().map(|t| t.as_str().to_owned()).collect();
        let count = texts.len();
        let in_range = |ids: &rkyv::vec::ArchivedVec<rkyv::Archived<u32>>| {
            ids.iter().all(|id| (id.to_native() as usize) < count)
        };
        if !root.entries.iter().all(|entry| {
            entry.positions.iter().all(|position| {
                position
                    .changes
                    .iter()
                    .all(|change| in_range(&change.clnsig) && in_range(&change.clnrevstat))
            })
        }) {
            return Err(invalid(format!("{path}: a text index is out of range")));
        }
        Ok(ClinvarAaCache { mapped, texts })
    }

    pub(crate) fn len(&self) -> usize {
        self.mapped.root().transcripts.len()
    }

    pub(crate) fn get(&self, transcript: &str) -> Option<&ArchivedAaTranscript> {
        let root = self.mapped.root();
        find_sorted(&root.transcripts, transcript).map(|at| &root.entries[at])
    }

    /// The text an index of `clnsig` or `clnrevstat` stands for.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn text(&self, id: &rkyv::Archived<u32>) -> &str {
        &self.texts[id.to_native() as usize]
    }

    fn plain(&self, transcript: &str) -> Option<PlainTranscript> {
        self.get(transcript)
            .map(|entry| plain_from_archive(entry, &self.texts))
    }
}

/// Converts the pickle at `pickle` to `output`, then reads the pickle again and checks that the
/// converted file gives every transcript the same value and holds no other. Returns a summary.
pub(crate) fn convert_clinvar_aa_stat(pickle: &str, output: &Path) -> Result<String> {
    ensure_distinct_paths(Path::new(pickle), output)?;
    let start = Instant::now();
    let entries = read_pickle(pickle)?;
    let positions: usize = entries.iter().map(|(_, value)| value.len()).sum();
    let changes: usize = entries
        .iter()
        .flat_map(|(_, value)| value.iter().map(|(_, list)| list.len()))
        .sum();
    let metadata = std::fs::metadata(pickle)?;
    let provenance = vec![
        Provenance {
            key: "source".to_owned(),
            value: pickle.to_owned(),
        },
        Provenance {
            key: "source_bytes".to_owned(),
            value: metadata.len().to_string(),
        },
        Provenance {
            key: "converter".to_owned(),
            value: format!("priva cache-convert {}", env!("CARGO_PKG_VERSION")),
        },
    ];
    let count = entries.len();
    let archive = build(entries, provenance)?;
    let converted_s = start.elapsed().as_secs_f64();

    // The check reads the pickle a second time, independently of the conversion.
    let check = Instant::now();
    let mut compared = 0_usize;
    let mut text_count = 0;
    write_cache_checked(output, CLINVAR_AA_STAT, &archive, |staged| {
        let staged_text = staged.to_string_lossy();
        let cache = ClinvarAaCache::open(&staged_text)?;
        if cache.len() != count {
            return Err(invalid(format!(
                "{} transcripts converted, {count} in the pickle",
                cache.len()
            )));
        }
        text_count = cache.texts.len();
        let distinct = visit_top_level_dict(
            pickle,
            |_| true,
            |transcript, value| {
                let expected = plain_transcript(pickle, &transcript, &value)?;
                if cache.plain(&transcript).as_ref() != Some(&expected) {
                    return Err(invalid(format!(
                        "transcript {transcript} differs from the pickle"
                    )));
                }
                compared += 1;
                Ok(())
            },
        )?;
        if compared != count || distinct != count {
            return Err(invalid("the pickle changed while it was converted"));
        }
        Ok(())
    })?;
    Ok(format!(
        "{{\"cache\": \"clinvar_aa_stat\", \"transcripts\": {count}, \"positions\": {positions}, \
         \"changes\": {changes}, \"texts\": {}, \"bytes\": {}, \"converted_s\": {converted_s:.1}, \
         \"checked_transcripts\": {compared}, \"check_s\": {:.1}}}",
        text_count,
        std::fs::metadata(output)?.len(),
        check.elapsed().as_secs_f64()
    ))
}

/// Full typed inspection, including every position/change and list occurrence.
fn export_clinvar_aa_stat(path: &str, output: &mut dyn std::io::Write) -> Result<usize> {
    use common::ExactValue as V;
    let cache = ClinvarAaCache::open(path)?;
    for key in cache.mapped.root().transcripts.iter() {
        let positions = cache.plain(key.as_str()).unwrap();
        let value = V::Dict(
            positions
                .into_iter()
                .map(|(position, changes)| {
                    (
                        V::Text(position),
                        V::Dict(
                            changes
                                .into_iter()
                                .map(|(change, clnsig, review)| {
                                    (
                                        V::Text(change),
                                        V::Dict(vec![
                                            (
                                                V::Text("CLNSIG".into()),
                                                V::List(clnsig.into_iter().map(V::Text).collect()),
                                            ),
                                            (
                                                V::Text("CLNREVSTAT".into()),
                                                V::List(review.into_iter().map(V::Text).collect()),
                                            ),
                                        ]),
                                    )
                                })
                                .collect(),
                        ),
                    )
                })
                .collect(),
        );
        common::write_entry(output, &V::Text(key.as_str().to_owned()), &value)?;
    }
    Ok(cache.len())
}

/// Inspect one cache without Python. `--output -` streams tagged rows for an independent
/// comparator; a file destination must be new and is removed if exporting fails.
pub(crate) fn dispatch_export(args: &[String]) -> Result<()> {
    use std::io::Write;
    const HELP: &str = "usage: priva cache-export CACHE --input CACHE --output JSONL_OR_DASH";
    let [cache, first, first_value, second, second_value] = args else {
        return Err(invalid(HELP));
    };
    let mut input = None;
    let mut output = None;
    for (name, value) in [(first, first_value), (second, second_value)] {
        match name.as_str() {
            "--input" if input.is_none() => input = Some(value),
            "--output" if output.is_none() => output = Some(value),
            _ => return Err(invalid(HELP)),
        }
    }
    let (Some(input), Some(output)) = (input, output) else {
        return Err(invalid(HELP));
    };
    let exporter = match cache.as_str() {
        "clinvar_aa_stat" => export_clinvar_aa_stat,
        "alphamissense_tranx_domain_map" => transcript_domains::export,
        "all_intolerant_domains" => intolerant_domains::export,
        "pm1_regions_pkl" => pm1_regions::export,
        _ => return Err(invalid(format!("cache-export: unsupported cache {cache}"))),
    };
    if output == "-" {
        let stdout = std::io::stdout();
        let mut writer = std::io::BufWriter::new(stdout.lock());
        exporter(input, &mut writer)?;
        writer.flush()?;
    } else {
        ensure_distinct_paths(Path::new(input), Path::new(output))?;
        let file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output)?;
        let result = (|| -> Result<()> {
            let mut writer = std::io::BufWriter::new(file);
            exporter(input, &mut writer)?;
            writer.flush()?;
            writer.get_ref().sync_all()?;
            Ok(())
        })();
        if result.is_err() {
            let _ = std::fs::remove_file(output);
        }
        result?;
    }
    Ok(())
}

const USAGE: &str = "usage: priva cache-convert CACHE --pickle PICKLE --output CACHE\nCACHE: clinvar_aa_stat, clinvar_splice_stat";

/// `priva cache-convert CACHE --pickle PICKLE --output OUTPUT`.
pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let [cache, first, first_value, second, second_value] = args else {
        return Err(invalid(USAGE));
    };
    let mut pickle = None;
    let mut output = None;
    for (name, value) in [(first, first_value), (second, second_value)] {
        match name.as_str() {
            "--pickle" if pickle.is_none() => pickle = Some(value),
            "--output" if output.is_none() => output = Some(value),
            _ => return Err(invalid(USAGE)),
        }
    }
    let (Some(pickle), Some(output)) = (pickle, output) else {
        return Err(invalid(USAGE));
    };
    if output == pickle {
        return Err(invalid("the output must not be the pickle itself"));
    }
    match cache.as_str() {
        "clinvar_aa_stat" => {
            println!("{}", convert_clinvar_aa_stat(pickle, Path::new(output))?);
            Ok(())
        }
        "alphamissense_tranx_domain_map" => {
            println!("{}", transcript_domains::convert(pickle, Path::new(output))?);
            Ok(())
        }
        "all_intolerant_domains" => {
            println!("{}", intolerant_domains::convert(pickle, Path::new(output))?);
            Ok(())
        }
        "pm1_regions_pkl" => {
            println!("{}", pm1_regions::convert(pickle, Path::new(output))?);
            Ok(())
        }
        "clinvar_splice_stat" => {
            println!(
                "{}",
                clinvar_splice_stat::convert(pickle, Path::new(output))?
            );
            Ok(())
        }
        other => Err(invalid(format!(
            "cache-convert: {other} is not converted yet (prototype: clinvar_aa_stat)\n{USAGE}"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn temporary(name: &str, bytes: &[u8]) -> String {
        let path =
            std::env::temp_dir().join(format!("priva_cache_convert_{}_{name}", std::process::id()));
        std::fs::File::create(&path)
            .unwrap()
            .write_all(bytes)
            .unwrap();
        path.to_str().unwrap().to_owned()
    }

    // pickle.dumps({"ENST2": {"10": {"p.A10V": {"CLNSIG": ["Pathogenic"], "CLNREVSTAT":
    //   ["criteria_provided,_multiple_submitters,_no_conflicts"]}, "p.A10T": {"CLNSIG":
    //   ["Benign", "Pathogenic"], "CLNREVSTAT": ["no_assertion_criteria_provided",
    //   "reviewed_by_expert_panel"]}}, "7-8": {"p.K7_L8del": {"CLNSIG": ["Pathogenic"],
    //   "CLNREVSTAT": []}}}, "ENST1": {"3": {"p.M3I": {"CLNSIG": ["Likely_pathogenic"],
    //   "CLNREVSTAT": ["criteria_provided,_single_submitter"]}}}}, protocol=5)
    const GOOD: &[u8] = b"\x80\x05\x95n\x01\x00\x00\x00\x00\x00\x00}\x94(\x8c\x05ENST2\x94}\x94(\x8c\x0210\x94}\x94(\x8c\x06p.A10V\x94}\x94(\x8c\x06CLNSIG\x94]\x94\x8c\x0aPathogenic\x94a\x8c\x0aCLNREVSTAT\x94]\x94\x8c4criteria_provided,_multiple_submitters,_no_conflicts\x94au\x8c\x06p.A10T\x94}\x94(h\x07]\x94(\x8c\x06Benign\x94h\x09eh\x0a]\x94(\x8c\x1eno_assertion_criteria_provided\x94\x8c\x18reviewed_by_expert_panel\x94euu\x8c\x037-8\x94}\x94\x8c\x0ap.K7_L8del\x94}\x94(h\x07]\x94h\x09ah\x0a]\x94usu\x8c\x05ENST1\x94}\x94\x8c\x013\x94}\x94\x8c\x05p.M3I\x94}\x94(h\x07]\x94\x8c\x11Likely_pathogenic\x94ah\x0a]\x94\x8c#criteria_provided,_single_submitter\x94aussu.";
    // {"ENST1": {"3": {"p.M3I": {"CLNSIG": ("Pathogenic",), "CLNREVSTAT": ["x"]}}}}, protocol=4
    const TUPLE: &[u8] = b"\x80\x04\x95N\x00\x00\x00\x00\x00\x00\x00}\x94\x8c\x05ENST1\x94}\x94\x8c\x013\x94}\x94\x8c\x05p.M3I\x94}\x94(\x8c\x06CLNSIG\x94\x8c\x0aPathogenic\x94\x85\x94\x8c\x0aCLNREVSTAT\x94]\x94\x8c\x01x\x94ausss.";
    // {"ENST1": {"3": {"p.M3I": {"CLNSIG": ["P"], "CLNREVSTAT": ["x"], "X": []}}}}, protocol=4
    const EXTRA: &[u8] = b"\x80\x04\x95L\x00\x00\x00\x00\x00\x00\x00}\x94\x8c\x05ENST1\x94}\x94\x8c\x013\x94}\x94\x8c\x05p.M3I\x94}\x94(\x8c\x06CLNSIG\x94]\x94\x8c\x01P\x94a\x8c\x0aCLNREVSTAT\x94]\x94\x8c\x01x\x94a\x8c\x01X\x94]\x94usss.";

    #[test]
    fn converts_and_answers_every_key_as_the_pickle() {
        let pickle = temporary("good.pkl", GOOD);
        let output = temporary("good.rkyv", b"");
        let summary = convert_clinvar_aa_stat(&pickle, Path::new(&output)).unwrap();
        assert!(summary.contains("\"transcripts\": 2"), "{summary}");
        assert!(summary.contains("\"changes\": 4"), "{summary}");
        assert!(summary.contains("\"checked_transcripts\": 2"), "{summary}");
        let cache = ClinvarAaCache::open(&output).unwrap();
        assert_eq!(cache.len(), 2);
        assert!(cache.get("ENST3").is_none());
        assert!(cache.get("").is_none());
        let enst2 = cache.plain("ENST2").unwrap();
        let positions: Vec<&str> = enst2.iter().map(|(pos, _)| pos.as_str()).collect();
        assert_eq!(positions, ["10", "7-8"]);
        assert_eq!(
            enst2[0].1[1],
            (
                "p.A10T".to_owned(),
                vec!["Benign".to_owned(), "Pathogenic".to_owned()],
                vec![
                    "no_assertion_criteria_provided".to_owned(),
                    "reviewed_by_expert_panel".to_owned()
                ]
            )
        );
        assert_eq!(
            enst2[1].1[0],
            (
                "p.K7_L8del".to_owned(),
                vec!["Pathogenic".to_owned()],
                vec![]
            )
        );
        let change = &cache.get("ENST1").unwrap().positions[0].changes[0];
        assert_eq!(cache.text(&change.clnsig[0]), "Likely_pathogenic");
        // The same pickle gives the same bytes again.
        let again = temporary("good2.rkyv", b"");
        convert_clinvar_aa_stat(&pickle, Path::new(&again)).unwrap();
        assert_eq!(
            std::fs::read(&output).unwrap(),
            std::fs::read(&again).unwrap()
        );
    }

    #[test]
    fn refuses_other_shapes() {
        let output = temporary("odd.rkyv", b"");
        let error = convert_clinvar_aa_stat(&temporary("tuple.pkl", TUPLE), Path::new(&output))
            .unwrap_err();
        assert!(error.to_string().contains("CLNSIG is a tuple"), "{error}");
        let error = convert_clinvar_aa_stat(&temporary("extra.pkl", EXTRA), Path::new(&output))
            .unwrap_err();
        assert!(error.to_string().contains("other keys"), "{error}");
    }

    #[test]
    fn refuses_pickles_damaged_and_foreign_caches() {
        let pickle = temporary("refuse.pkl", GOOD);
        let error = ClinvarAaCache::open(&pickle).err().unwrap();
        assert!(error.to_string().contains("cache-convert"), "{error}");
        let output = temporary("refuse.rkyv", b"");
        convert_clinvar_aa_stat(&pickle, Path::new(&output)).unwrap();
        let mut bytes = std::fs::read(&output).unwrap();
        // A wrong schema id.
        let mut foreign = bytes.clone();
        foreign[12] = 2;
        let error = ClinvarAaCache::open(&temporary("foreign.rkyv", &foreign))
            .err()
            .unwrap();
        assert!(error.to_string().contains("schema 2"), "{error}");
        // Truncated.
        let error = ClinvarAaCache::open(&temporary("short.rkyv", &bytes[..bytes.len() - 8]))
            .err()
            .unwrap();
        assert!(error.to_string().contains("archive length"), "{error}");
        // The root pointer at the end of the archive overwritten.
        let end = bytes.len();
        for byte in &mut bytes[end - 8..] {
            *byte = 0xff;
        }
        let error = ClinvarAaCache::open(&temporary("damaged.rkyv", &bytes))
            .err()
            .unwrap();
        assert!(error.to_string().contains("damaged"), "{error}");
    }

    #[test]
    fn command_line() {
        assert!(dispatch(&[]).is_err());
        let args = |a: &[&str]| a.iter().map(|s| s.to_string()).collect::<Vec<_>>();
        let error =
            dispatch(&args(&["unknown_cache", "--pickle", "a", "--output", "b"])).unwrap_err();
        assert!(error.to_string().contains("not converted yet"), "{error}");
        assert!(
            dispatch(&args(&[
                "clinvar_aa_stat",
                "--pickle",
                "a",
                "--pickle",
                "b"
            ]))
            .is_err()
        );
        assert!(
            dispatch(&args(&[
                "clinvar_aa_stat",
                "--pickle",
                "a",
                "--output",
                "a"
            ]))
            .is_err()
        );
    }
}
