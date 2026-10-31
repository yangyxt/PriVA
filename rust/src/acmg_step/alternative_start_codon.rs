//! `scripts/alternative_start_codon.py` (405 lines): the CDS FASTA index and the downstream-ATG
//! search for start_lost variants, and the PVS1 strength they lead to.
//!
//! `AlternativeStartCodonDetector` is the Python class as PVS1 uses it (`get_detector` then
//! `find_alternative_start`); `adjust_pvs1_for_start_lost` is the Python function.
//!
//! The Python class reads `<fasta>.cds_cache.pkl.gz` when that file is newer than the FASTA,
//! and writes it when it is missing (group 3 plan, section 8.8). The cache is the FASTA's own
//! parse, so Rust reads the FASTA directly and never reads or writes the cache: the results
//! differ only if a cache newer than its FASTA holds other sequences. The command line of the
//! script (which the installer runs to pre-build the cache) stays in Python.
//!
//! What Python's `get_detector` can fail on while reading the FASTA (a missing file, a broken
//! gzip stream, text that is not UTF-8) is a `DetectorFailure`, which PVS1 catches as Python
//! does. Only the sequences of the transcripts asked for are kept.
use super::acmg_pp3_bp4_bp7_insilico::harness_options;
use super::types::{HarnessCall, not_converted};
use crate::cli::{Result, invalid};
use crate::pandas_table::{Cell, Column, Table};
use crate::python_text::python_strip;
use flate2::read::MultiGzDecoder;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Read;

/// The dict `find_alternative_start` returns.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct AlternativeStart {
    pub(crate) has_alternative: bool,
    pub(crate) first_alt_position: Option<i64>,
    pub(crate) all_alt_positions: Vec<i64>,
    pub(crate) n_terminal_truncation_fraction: Option<f64>,
    pub(crate) total_protein_length: Option<i64>,
    pub(crate) error: Option<String>,
}

/// An exception Python's detector raises while reading the FASTA; PVS1 catches it
/// (`except Exception`) and falls back to PVS1_Moderate for every start_lost variant.
#[derive(Debug)]
pub(crate) struct DetectorFailure(pub(crate) String);

/// `AlternativeStartCodonDetector`, holding the CDS of the transcripts asked for.
pub(crate) struct AlternativeStartCodonDetector {
    /// `cds_sequences`: the last sequence of each header ID (with version).
    cds_sequences: HashMap<String, String>,
    /// `transcript_to_key`: the first header ID of each ID without version.
    transcript_to_key: HashMap<String, String>,
}

/// `transcript_id.split('.')[0]`.
fn without_version(transcript_id: &str) -> &str {
    transcript_id.split('.').next().unwrap_or("")
}

impl AlternativeStartCodonDetector {
    /// `get_detector(cds_fasta_path)`, keeping the sequences `find_alternative_start` can
    /// reach for the transcripts in `wanted`: those whose header ID, without its version,
    /// is the version-less form of one of them.
    pub(crate) fn load(
        cds_fasta_path: &str,
        wanted: &[&str],
    ) -> Result<std::result::Result<AlternativeStartCodonDetector, DetectorFailure>> {
        let bases: HashSet<&str> = wanted.iter().map(|id| without_version(id)).collect();
        let text = match read_fasta_text(cds_fasta_path) {
            Ok(text) => text,
            Err(failure) => return Ok(Err(failure)),
        };
        let mut detector = AlternativeStartCodonDetector {
            cds_sequences: HashMap::new(),
            transcript_to_key: HashMap::new(),
        };
        // `for line in f: line = line.strip()`: lines end at "\n" after the text mode turned
        // "\r\n" and "\r" into "\n"; a header without any sequence line is skipped, a header
        // followed only by blank lines keeps an empty sequence.
        let mut header: Option<&str> = None;
        let mut sequence: Option<Vec<&str>> = None;
        for raw_line in text.split_inclusive('\n') {
            let line = python_strip(raw_line);
            if let Some(rest) = line.strip_prefix('>') {
                if let (Some(name), Some(parts)) = (header, sequence.take()) {
                    detector.add_entry(name, &parts, &bases)?;
                }
                header = (!rest.is_empty()).then_some(rest);
                sequence = None;
            } else {
                sequence.get_or_insert_with(Vec::new).push(line);
            }
        }
        if let (Some(name), Some(parts)) = (header, sequence) {
            detector.add_entry(name, &parts, &bases)?;
        }
        Ok(Ok(detector))
    }

    /// `_process_fasta_entry`, for the transcripts that can be asked for.
    fn add_entry(&mut self, header: &str, parts: &[&str], bases: &HashSet<&str>) -> Result<()> {
        let transcript_with_version = header
            .split(crate::python_text::is_python_space)
            .find(|part| !part.is_empty())
            .unwrap_or("");
        let transcript_id = without_version(transcript_with_version);
        if !bases.contains(transcript_id) {
            return Ok(());
        }
        let sequence: String = parts.concat();
        if !sequence.is_ascii() {
            return Err(invalid(format!(
                "the CDS of {transcript_with_version} holds non-ASCII text; Python's upper() and \
                 codon slicing of it are not reproduced"
            )));
        }
        self.cds_sequences.insert(
            transcript_with_version.to_owned(),
            sequence.to_ascii_uppercase(),
        );
        self.transcript_to_key
            .entry(transcript_id.to_owned())
            .or_insert_with(|| transcript_with_version.to_owned());
        Ok(())
    }

    /// `find_alternative_start(transcript_id)`. Only transcripts given to `load` can be found.
    pub(crate) fn find_alternative_start(&self, transcript_id: &str) -> AlternativeStart {
        let sequence = match self.cds_sequences.get(transcript_id) {
            Some(sequence) => Some(sequence),
            None => self
                .transcript_to_key
                .get(without_version(transcript_id))
                .and_then(|full_key| self.cds_sequences.get(full_key)),
        };
        let Some(sequence) = sequence else {
            return AlternativeStart {
                has_alternative: false,
                first_alt_position: None,
                all_alt_positions: Vec::new(),
                n_terminal_truncation_fraction: None,
                total_protein_length: None,
                error: Some(format!("Transcript {transcript_id} not found in CDS FASTA")),
            };
        };
        let bytes = sequence.as_bytes();
        let total_protein_length = (bytes.len() / 3) as i64;
        // `range(3, len(cds_seq) - 2, 3)`: in-frame codons after the first.
        let all_alt_positions: Vec<i64> = (3..bytes.len().saturating_sub(2))
            .step_by(3)
            .filter(|start| &bytes[*start..*start + 3] == b"ATG")
            .map(|start| (start / 3) as i64 + 1)
            .collect();
        match all_alt_positions.first() {
            Some(&first) => AlternativeStart {
                has_alternative: true,
                first_alt_position: Some(first),
                n_terminal_truncation_fraction: Some(python_round_4(
                    (first - 1) as f64 / total_protein_length as f64,
                )),
                all_alt_positions,
                total_protein_length: Some(total_protein_length),
                error: None,
            },
            None => AlternativeStart {
                has_alternative: false,
                first_alt_position: None,
                all_alt_positions,
                n_terminal_truncation_fraction: None,
                total_protein_length: Some(total_protein_length),
                error: None,
            },
        }
    }
}

/// Python's `round(x, 4)` for `x >= 0`: the 4-digit decimal nearest to the exact value of `x`,
/// ties to even, read back as a float. `format!` rounds the exact value correctly; an exact
/// tie needs `x * 10000` to end in `.5`, which for a binary fraction happens only at odd
/// multiples of 1/32, and those are settled here to the even neighbour.
fn python_round_4(value: f64) -> f64 {
    let thirty_seconds = value * 32.0;
    if thirty_seconds.fract() == 0.0 && thirty_seconds < 1.0e12 && thirty_seconds as i64 % 2 == 1 {
        // value * 10000 = (j * 625) / 2 with j * 625 odd.
        let lower = (thirty_seconds as i64 * 625 - 1) / 2;
        let even = if lower % 2 == 0 { lower } else { lower + 1 };
        return format!("{even}e-4").parse().expect("a decimal reads back");
    }
    format!("{value:.4}")
        .parse()
        .expect("a formatted float reads back")
}

/// The FASTA as Python's text mode reads it (`gzip.open(path, 'rt')` when the name ends in
/// `.gz`, else `open(path, 'rt')`), with "\r\n" and "\r" turned into "\n".
fn read_fasta_text(path: &str) -> std::result::Result<String, DetectorFailure> {
    let failure = |what: String| DetectorFailure(format!("{path}: {what}"));
    let mut raw = Vec::new();
    File::open(path)
        .and_then(|mut file| file.read_to_end(&mut raw))
        .map_err(|error| failure(error.to_string()))?;
    if path.ends_with(".gz") {
        let mut plain = Vec::new();
        MultiGzDecoder::new(&raw[..])
            .read_to_end(&mut plain)
            .map_err(|error| failure(error.to_string()))?;
        raw = plain;
    }
    let text = String::from_utf8(raw).map_err(|_| failure("not UTF-8 text".to_owned()))?;
    Ok(if text.contains('\r') {
        text.replace("\r\n", "\n").replace('\r', "\n")
    } else {
        text
    })
}

/// `adjust_pvs1_for_start_lost(info, spans_critical_domain, gene_is_lof_intolerant)`;
/// `None` is the empty dict PVS1 passes for a transcript it did not look up.
pub(crate) fn adjust_pvs1_for_start_lost(
    info: Option<&AlternativeStart>,
    spans_critical_domain: bool,
    gene_is_lof_intolerant: bool,
) -> i64 {
    if !gene_is_lof_intolerant {
        return 0;
    }
    // `info.get('error')` is truthy for a non-empty message.
    if info
        .and_then(|info| info.error.as_deref())
        .is_some_and(|error| !error.is_empty())
    {
        return 2;
    }
    let has_alternative = info.is_some_and(|info| info.has_alternative);
    if !has_alternative {
        return 4;
    }
    let truncation = info
        .and_then(|info| info.n_terminal_truncation_fraction)
        .unwrap_or(0.0);
    match (truncation > 0.10, spans_critical_domain) {
        (true, true) => 3,
        (true, false) => 2,
        (false, true) => 2,
        (false, false) => 1,
    }
}

/// `priva acmg-step alternative_start_codon <function>`:
///
/// - `find_alternative_start`, option `cds_fasta_path`, on a frame with a text column
///   `transcript`: one row per transcript with `has_alternative` (bool), `first_alt_position`
///   and `total_protein_length` (int64, -1 for None), `all_alt_positions` (the positions
///   joined by commas), `n_terminal_truncation_fraction` (float64, NaN for None) and `error`
///   (text, empty for None); a FASTA the detector cannot read gives the error `DetectorFailure`;
/// - `adjust_pvs1_for_start_lost` on a frame of `has_info`, `has_alternative`, `spans`,
///   `lof_intolerant` (bool), `truncation` (float64, NaN for None) and `error` (text, missing
///   for None): a frame of `strength` (int64).
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let table = &call.table;
    match function {
        "find_alternative_start" => {
            let path = harness_options::required_text(&call.options, "cds_fasta_path")?;
            let transcripts = text_cells(table, "transcript")?;
            let wanted: Vec<&str> = transcripts.iter().map(String::as_str).collect();
            let detector = AlternativeStartCodonDetector::load(path, &wanted)?
                .map_err(|failure| invalid(format!("DetectorFailure: {}", failure.0)))?;
            let found: Vec<AlternativeStart> = wanted
                .iter()
                .map(|id| detector.find_alternative_start(id))
                .collect();
            Ok(harness_options::frame(
                vec![
                    (
                        "has_alternative".to_owned(),
                        Column::Bool(found.iter().map(|info| info.has_alternative).collect()),
                    ),
                    (
                        "first_alt_position".to_owned(),
                        Column::Int(
                            found
                                .iter()
                                .map(|info| info.first_alt_position.unwrap_or(-1))
                                .collect(),
                        ),
                    ),
                    (
                        "all_alt_positions".to_owned(),
                        Column::Text(
                            found
                                .iter()
                                .map(|info| {
                                    Some(
                                        info.all_alt_positions
                                            .iter()
                                            .map(i64::to_string)
                                            .collect::<Vec<_>>()
                                            .join(","),
                                    )
                                })
                                .collect(),
                        ),
                    ),
                    (
                        "n_terminal_truncation_fraction".to_owned(),
                        Column::Float(
                            found
                                .iter()
                                .map(|info| info.n_terminal_truncation_fraction.unwrap_or(f64::NAN))
                                .collect(),
                        ),
                    ),
                    (
                        "total_protein_length".to_owned(),
                        Column::Int(
                            found
                                .iter()
                                .map(|info| info.total_protein_length.unwrap_or(-1))
                                .collect(),
                        ),
                    ),
                    (
                        "error".to_owned(),
                        Column::Text(
                            found
                                .iter()
                                .map(|info| Some(info.error.clone().unwrap_or_default()))
                                .collect(),
                        ),
                    ),
                ],
                table.rows,
            ))
        }
        "adjust_pvs1_for_start_lost" => {
            let flag = |name: &str| -> Result<Vec<bool>> {
                match table.column(name)? {
                    Column::Bool(values) => Ok(values.clone()),
                    _ => Err(invalid(format!("{name} must be a bool column"))),
                }
            };
            let (has_info, has_alternative, spans, lof) = (
                flag("has_info")?,
                flag("has_alternative")?,
                flag("spans")?,
                flag("lof_intolerant")?,
            );
            let truncation = match table.column("truncation")? {
                Column::Float(values) => values.clone(),
                _ => return Err(invalid("truncation must be a float64 column")),
            };
            let error = table.column("error")?;
            let strengths = (0..table.rows)
                .map(|row| {
                    let info = AlternativeStart {
                        has_alternative: has_alternative[row],
                        first_alt_position: None,
                        all_alt_positions: Vec::new(),
                        n_terminal_truncation_fraction: (!truncation[row].is_nan())
                            .then_some(truncation[row]),
                        total_protein_length: None,
                        error: match error.cell(row) {
                            Cell::Text(text) => Some(text.to_owned()),
                            _ => None,
                        },
                    };
                    adjust_pvs1_for_start_lost(has_info[row].then_some(&info), spans[row], lof[row])
                })
                .collect();
            Ok(harness_options::frame(
                vec![("strength".to_owned(), Column::Int(strengths))],
                table.rows,
            ))
        }
        _ => Err(not_converted("alternative_start_codon", function)),
    }
}

fn text_cells(table: &Table, name: &str) -> Result<Vec<String>> {
    let column = table.column(name)?;
    (0..table.rows)
        .map(|row| match column.cell(row) {
            Cell::Text(text) => Ok(text.to_owned()),
            _ => Err(invalid(format!("{name} must hold text in every row"))),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn fasta(name: &str, text: &str) -> String {
        let path =
            std::env::temp_dir().join(format!("priva_alt_start_{}_{name}.fa", std::process::id()));
        File::create(&path)
            .unwrap()
            .write_all(text.as_bytes())
            .unwrap();
        path.to_str().unwrap().to_owned()
    }

    #[test]
    fn finds_in_frame_atg_after_the_first_codon() {
        let path = fasta(
            "basic",
            ">ENST1.2 cds chromosome:GRCh38:1:1:30:1\natgAAAATG\nCCCATGTAA\n>ENST2.1 cds\nATGCCC\n\
             >ENST1.3 other\nATGTTT\r\n",
        );
        let detector = AlternativeStartCodonDetector::load(&path, &["ENST1", "ENST2.1", "ENST9"])
            .unwrap()
            .unwrap();
        let info = detector.find_alternative_start("ENST1.2");
        assert_eq!(info.all_alt_positions, vec![3, 5]);
        assert_eq!(info.total_protein_length, Some(6));
        assert_eq!(info.n_terminal_truncation_fraction, Some(0.3333));
        // Without a version the first header of the transcript is used.
        assert_eq!(detector.find_alternative_start("ENST1"), info);
        // ENST1.3 is a later header with its own version.
        assert!(!detector.find_alternative_start("ENST1.3").has_alternative);
        let none = detector.find_alternative_start("ENST2.1");
        assert_eq!((none.has_alternative, none.error), (false, None));
        let missing = detector.find_alternative_start("ENST9");
        assert_eq!(
            missing.error.as_deref(),
            Some("Transcript ENST9 not found in CDS FASTA")
        );
    }

    #[test]
    fn strength_tree() {
        let info = |has_alternative, truncation: Option<f64>| AlternativeStart {
            has_alternative,
            first_alt_position: None,
            all_alt_positions: Vec::new(),
            n_terminal_truncation_fraction: truncation,
            total_protein_length: None,
            error: None,
        };
        assert_eq!(adjust_pvs1_for_start_lost(None, true, false), 0);
        assert_eq!(adjust_pvs1_for_start_lost(None, false, true), 4);
        assert_eq!(
            adjust_pvs1_for_start_lost(Some(&info(true, Some(0.2))), true, true),
            3
        );
        assert_eq!(
            adjust_pvs1_for_start_lost(Some(&info(true, Some(0.1))), true, true),
            2
        );
        assert_eq!(
            adjust_pvs1_for_start_lost(Some(&info(true, None)), false, true),
            1
        );
        let mut failed = info(true, Some(0.5));
        failed.error = Some("x".into());
        assert_eq!(adjust_pvs1_for_start_lost(Some(&failed), true, true), 2);
        failed.error = Some(String::new());
        assert_eq!(adjust_pvs1_for_start_lost(Some(&failed), true, true), 3);
    }

    #[test]
    fn rounds_like_python() {
        assert_eq!(python_round_4(0.12345), 0.1235);
        assert_eq!(python_round_4(0.00005), 0.0001);
        assert_eq!(python_round_4(2.0 / 3.0), 0.6667);
        // Exact ties go to the even digit.
        assert_eq!(python_round_4(1.0 / 32.0), 0.0312);
        assert_eq!(python_round_4(3.0 / 32.0), 0.0938);
        assert_eq!(python_round_4(5.0 / 32.0), 0.1562);
        assert_eq!(python_round_4(0.0), 0.0);
    }
}
