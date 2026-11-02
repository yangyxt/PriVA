//! Native installer path of prepare_pm1_regions.py with --skip_das.
//!
//! Keep mapping -> existing DAS -> HCSeeker -> RMC -> bundle order. HCSeeker
//! and RMC use Python text parsing, not pandas. Transcript multiplicity is
//! retained until the original set unions; the native schema sorts lookup
//! keys and set members without changing their values.
use crate::cache_convert::Provenance;
use crate::cache_convert::intolerant_domains::IntolerantDomainCache;
use crate::cache_convert::pm1_regions::{
    Pm1Metadata, Pm1Regions, TranscriptResidues, write_produced,
};
use crate::cli::{Result, invalid};
use crate::python_args::{Parsed, parse_long_options};
use crate::python_text::{decimal_digit_value, python_int_is_valid, python_strip};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

const HELP: &str = "usage: priva pm1-regions build --assembly hg19|hg38 --nm-to-enst TSV --intolerant-domains RKYV --hcseeker TSV --rmc TSV --output RKYV [--rmc-oe-threshold 0.4] [--rmc-p-bh-threshold 0.05]\n\
Build the original --skip_das PM1 bundle from explicit prepared inputs. No DAS calculation, source download or pickle output.";

type Regions = BTreeMap<String, BTreeSet<i64>>;

/// Python open(..., newline=None) translates CRLF and bare CR before a script
/// strips its newline. These simple TSV readers do not interpret CSV quoting.
fn text(path: &Path) -> Result<String> {
    Ok(std::fs::read_to_string(path)?
        .replace("\r\n", "\n")
        .replace('\r', "\n"))
}

fn decimal_text(value: &str) -> String {
    value
        .chars()
        .map(|c| decimal_digit_value(c).map_or(c, |d| char::from(b'0' + d as u8)))
        .collect()
}

fn integer(value: &str) -> Result<i64> {
    if !python_int_is_valid(value) {
        return Err(invalid(
            "ValueError: invalid literal for int() in a PM1 residue",
        ));
    }
    let normalized = decimal_text(python_strip(value)).replace('_', "");
    normalized
        .parse()
        .map_err(|_| invalid("PM1 residue is outside the native schema's signed 64-bit range"))
}

fn float(value: &str) -> Result<f64> {
    crate::python_text::parse_float(&decimal_text(value))
        .ok_or_else(|| invalid("ValueError: could not convert PM1 text to float"))
}

fn nm_mapping(input: &str) -> HashMap<String, Vec<String>> {
    let mut result: HashMap<String, Vec<String>> = HashMap::new();
    for line in input.split_terminator('\n').skip(1) {
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }
        let transcript = python_strip(fields[1]);
        let refseq = python_strip(fields[4]);
        if !refseq.starts_with("NM_") {
            continue;
        }
        let name = refseq.split('.').next().unwrap();
        let entries = result.entry(name.to_owned()).or_default();
        if !entries.iter().any(|known| known == transcript) {
            entries.push(transcript.to_owned());
        }
    }
    result
}

fn add_region(regions: &mut Regions, transcript: &str, start: i64, end: i64) {
    // Even an empty Python range creates this transcript's defaultdict entry.
    let positions = regions.entry(transcript.to_owned()).or_default();
    if start <= end {
        positions.extend(start..=end);
    }
}

fn hotspots(input: &str, mapping: &HashMap<String, Vec<String>>) -> Result<Regions> {
    let mut result = Regions::new();
    for line in input.split_terminator('\n').skip(1) {
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() < 7 || python_strip(fields[0]) != "hotspot" {
            continue;
        }
        let isoform = python_strip(fields[4]);
        // Parse positions before checking mapping, as the original does.
        let start = integer(fields[5])?;
        let end = integer(fields[6])?;
        let name = isoform.split('.').next().unwrap();
        if let Some(transcripts) = mapping.get(name) {
            for transcript in transcripts {
                add_region(&mut result, transcript, start, end);
            }
        }
    }
    Ok(result)
}

fn aa_position(value: &str) -> Result<i64> {
    let mut digits = String::new();
    for c in value.chars() {
        match decimal_digit_value(c) {
            Some(d) => digits.push(char::from(b'0' + d as u8)),
            None if !digits.is_empty() => break,
            None => (),
        }
    }
    if digits.is_empty() {
        return Err(invalid(format!(
            "ValueError: Cannot parse AA position from '{value}'"
        )));
    }
    integer(&digits)
}

fn constrained(input: &str, oe_threshold: f64, p_threshold: f64) -> Result<Regions> {
    let mut lines = input.split_terminator('\n');
    let header = lines.next().unwrap_or("");
    let columns: HashMap<_, _> = header
        .split('\t')
        .enumerate()
        .map(|(at, name)| (name, at))
        .collect();
    let mut records = Vec::new();
    for line in lines {
        let fields: Vec<_> = line.split('\t').collect();
        let field = |key: &str| -> Result<&str> {
            let at = columns
                .get(key)
                .ok_or_else(|| invalid(format!("KeyError: {key}")))?;
            fields
                .get(*at)
                .copied()
                .ok_or_else(|| invalid("IndexError: list index out of range"))
        };
        let transcript = python_strip(field("transcript")?).to_owned();
        let start_raw = python_strip(field("start_aa_raw")?);
        let end_raw = python_strip(field("end_aa_raw")?);
        let oe = float(field("oe")?)?;
        let p = float(field("p")?)?;
        let start = aa_position(start_raw)?;
        let end = aa_position(end_raw)?;
        records.push((transcript, start, end, oe, p));
    }
    if records.is_empty() {
        return Ok(Regions::new());
    }
    let p_values: Vec<_> = records.iter().map(|r| r.4).collect();
    let adjusted = crate::fdr::benjamini_hochberg(&p_values)
        .map_err(|e| invalid(format!("ValueError: {e}")))?;
    let mut result = Regions::new();
    for ((transcript, start, end, oe, _), p) in records.into_iter().zip(adjusted) {
        if oe < oe_threshold && p < p_threshold {
            add_region(&mut result, &transcript, start, end);
        }
    }
    Ok(result)
}

fn residues(regions: Regions) -> Vec<TranscriptResidues> {
    regions
        .into_iter()
        .map(|(transcript, positions)| TranscriptResidues {
            transcript,
            residues: positions.into_iter().collect(),
        })
        .collect()
}

struct Inputs {
    assembly: String,
    mapping: PathBuf,
    das: PathBuf,
    hcseeker: PathBuf,
    rmc: PathBuf,
    output: PathBuf,
    oe_threshold: f64,
    p_threshold: f64,
}

fn build(input: Inputs) -> Result<()> {
    let paths = [
        ("nm_to_enst", &input.mapping),
        ("all_intolerant_domains", &input.das),
        ("hcseeker", &input.hcseeker),
        ("rmc", &input.rmc),
    ];
    let mut provenance = Vec::new();
    let mut hashes = Vec::new();
    let mut capture = |name: &str, path: &Path| -> Result<()> {
        crate::rust_cache::ensure_distinct_paths(path, &input.output)?;
        let hash = crate::sha256::sha256_file(path)?;
        hashes.push(hash.clone());
        provenance.push(Provenance {
            key: format!("{name}_sha256"),
            value: hash,
        });
        Ok(())
    };
    // Do not pre-open later sources: the original reads the mapping, then DAS,
    // then HCSeeker, then RMC. An earlier parse failure must win over a missing
    // later input. Record identity at each stage and recheck before publication.
    capture(paths[0].0, paths[0].1)?;
    let mapping = nm_mapping(&text(&input.mapping)?);
    capture(paths[1].0, paths[1].1)?;
    let das = IntolerantDomainCache::open(
        input
            .das
            .to_str()
            .ok_or_else(|| invalid("DAS path is not UTF-8"))?,
    )?;
    let domains: Vec<_> = das.domains().map(str::to_owned).collect();
    capture(paths[2].0, paths[2].1)?;
    let hcseeker = residues(hotspots(&text(&input.hcseeker)?, &mapping)?);
    capture(paths[3].0, paths[3].1)?;
    let rmc = residues(constrained(
        &text(&input.rmc)?,
        input.oe_threshold,
        input.p_threshold,
    )?);
    let total =
        |regions: &[TranscriptResidues]| regions.iter().map(|r| r.residues.len() as i64).sum();
    let metadata = vec![
        Pm1Metadata::Assembly(input.assembly),
        Pm1Metadata::RmcOeThreshold(input.oe_threshold),
        Pm1Metadata::RmcPBhThreshold(input.p_threshold),
        Pm1Metadata::NDasIntolerantDomains(domains.len() as i64),
        Pm1Metadata::NHcseekerHotspotTranscripts(hcseeker.len() as i64),
        Pm1Metadata::NHcseekerHotspotResidues(total(&hcseeker)),
        Pm1Metadata::NRmcConstrainedTranscripts(rmc.len() as i64),
        Pm1Metadata::NRmcConstrainedResidues(total(&rmc)),
    ];
    let root = Pm1Regions {
        provenance,
        das_intolerant_domains: domains,
        hcseeker_hotspots: hcseeker,
        rmc_constrained: rmc,
        metadata,
    };
    for ((_, path), expected) in paths.into_iter().zip(hashes) {
        if crate::sha256::sha256_file(path)? != expected {
            return Err(invalid(
                "a PM1 source changed during production; no cache published",
            ));
        }
    }
    write_produced(&input.output, &root)?;
    println!(
        "{{\"cache\":\"pm1_regions_pkl\",\"das_domains\":{},\"hcseeker_transcripts\":{},\"rmc_transcripts\":{}}}",
        root.das_intolerant_domains.len(),
        root.hcseeker_hotspots.len(),
        root.rmc_constrained.len()
    );
    Ok(())
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args.is_empty() || matches!(args[0].as_str(), "--help" | "-h") {
        println!("{HELP}");
        return Ok(());
    }
    if args[0] != "build" {
        return Err(invalid("PM1 procedure must be build"));
    }
    let Parsed::Values(values) = parse_long_options(
        &args[1..],
        &[
            "assembly",
            "nm-to-enst",
            "intolerant-domains",
            "hcseeker",
            "rmc",
            "output",
            "rmc-oe-threshold",
            "rmc-p-bh-threshold",
        ],
    )
    .map_err(invalid)?
    else {
        println!("{HELP}");
        return Ok(());
    };
    let required = |key: &str| {
        values
            .get(key)
            .cloned()
            .ok_or_else(|| invalid(format!("missing --{key}")))
    };
    let assembly = required("assembly")?;
    if assembly != "hg19" && assembly != "hg38" {
        return Err(invalid("assembly must be hg19 or hg38"));
    }
    let inputs = Inputs {
        assembly,
        mapping: required("nm-to-enst")?.into(),
        das: required("intolerant-domains")?.into(),
        hcseeker: required("hcseeker")?.into(),
        rmc: required("rmc")?.into(),
        output: required("output")?.into(),
        oe_threshold: float(values.get("rmc-oe-threshold").map_or("0.4", String::as_str))?,
        p_threshold: float(
            values
                .get("rmc-p-bh-threshold")
                .map_or("0.05", String::as_str),
        )?,
    };
    match build(inputs) {
        Ok(()) => Ok(()),
        Err(error) => {
            eprintln!("priva pm1-regions: {error}");
            std::process::exit(1)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mapping_multiplicity_and_empty_hotspots_are_preserved() {
        let mapping =
            nm_mapping("h\na\tENST2\tb\tc\tNM_1.2\na\tENST1\tb\tc\tNM_1.3\na\tENST2\tb\tc\tNM_1\n");
        assert_eq!(mapping["NM_1"], vec!["ENST2", "ENST1"]);
        let regions = hotspots("h\nhotspot\ta\tb\tc\tNM_1\t4\t2\n", &mapping).unwrap();
        assert_eq!(regions.len(), 2);
        assert!(regions.values().all(BTreeSet::is_empty));
        assert!(hotspots("h\nhotspot\ta\tb\tc\tUNKNOWN\tx\t2\n", &mapping).is_err());
    }

    #[test]
    fn rmc_uses_strict_thresholds_and_first_unicode_decimal_run() {
        assert_eq!(aa_position("Arg９７and100").unwrap(), 97);
        assert_eq!(aa_position("-Arg1.5").unwrap(), 1);
        let rows = "transcript\tstart_aa_raw\tend_aa_raw\toe\tp\nENST1\tArg1\tArg2\t0.4\t0.01\nENST2\tArg3\tArg4\t0.3\t0.025\n";
        let regions = constrained(rows, 0.4, 0.05).unwrap();
        assert!(!regions.contains_key("ENST1"));
        assert_eq!(regions["ENST2"], BTreeSet::from([3, 4]));
        assert!(constrained(rows, f64::NAN, 0.05).unwrap().is_empty());
        assert!(constrained("nonsense\n", 0.4, 0.05).unwrap().is_empty());
    }
}
