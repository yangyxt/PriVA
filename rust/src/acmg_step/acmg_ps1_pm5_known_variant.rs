//! `scripts/acmg_ps1_pm5_known_variant.py` (288 lines): PS1 and PM5 from the ClinVar amino-acid and
//! splice dictionaries.
//!
//! `ps1_pm5_criteria` is `PS1_PM5_criteria`, `check_splice_pathogenic` the per-row splice
//! fallback. The two ClinVar pickles are read once by `ClinvarKnownVariants::load` (the step
//! loads them once and shares them with `summarize_clinvar_gene_pathogenicity` of PVS1).
//!
//! The pickles hold 6.2 M (amino-acid) and 5.0 M (splice) ClinVar records; read whole they take
//! 2.5 GB and 7.6 GB in Python. They are decoded with `pickle_dict_stream` straight into what the
//! criteria ask of them, for the transcripts of the table only: PS1/PM5 look a transcript up by
//! the table's `Feature`, and the gene summary only counts transcripts the table maps to a gene.
//! Every amino-acid entry is still checked, because PS1_PM5_criteria reads all of them. Shapes
//! other than the installed ones (probe of 2026-10-30: text keys, lists of text, an int `pos`,
//! a dict of text for `splice_ai`) are refused with a message naming the transcript.
//!
//! Python's `get_variant_type` (acmg_consequence.py, thread T1) is copied here as
//! `variant_type_of_hgvsp` until that module lands; the integration replaces the copy with it.
use super::acmg_pp3_bp4_bp7_insilico::harness_options;
use super::acmg_pp3_bp4_bp7_insilico::pandas_masks::{high_confidence_review, python_str};
use super::pickle_dict_stream::{Item, TextDict, visit_top_level_dict};
use super::types::{AcmgError, HarnessCall, StepOutput, not_converted};
use crate::cli::{Result, invalid};
use crate::pandas_table::{Cell, Column, ObjectValue, Table};
use crate::python_text::{parse_float, python_strip};
use crate::splicing_var_analysis::parse_hgvsc_splice_position;
use std::collections::{HashMap, HashSet};
use std::rc::Rc;

/// The default `high_confidence_status` of PS1_PM5_criteria and
/// summarize_clinvar_gene_pathogenicity: ClinVar review statuses and their stars.
pub(crate) fn review_stars(status: &str) -> i64 {
    match status {
        "practice_guideline" => 4,
        "reviewed_by_expert_panel" => 3,
        "criteria_provided,_multiple_submitters,_no_conflicts" => 2,
        _ => 0,
    }
}

/// `qualifies(clnsig, revstat)` of summarize_clinvar_gene_pathogenicity: CLNSIG split on
/// `, / | ;` holds `pathogenic` or `likely_pathogenic` (case and white space aside), and the
/// review status has two stars or more. The texts are ASCII (checked when decoding).
pub(crate) fn gene_level_pathogenic(clnsig: &str, revstat: &str) -> bool {
    let tokens = clnsig
        .split([',', '/', '|', ';'])
        .map(|token| python_strip(token));
    let listed = tokens.filter(|token| !token.is_empty()).any(|token| {
        token.eq_ignore_ascii_case("pathogenic") || token.eq_ignore_ascii_case("likely_pathogenic")
    });
    listed && review_stars(revstat) >= 2
}

/// What PS1, PM5 and the gene summary need of one transcript of the amino-acid dict.
#[derive(Debug, Default)]
struct AminoAcidTranscript {
    /// Some entry passes the gene summary's test.
    gene_level_pathogenic: bool,
    /// The protein position keys (`p in clinvar_aa_dict[f]`).
    positions: HashSet<Rc<str>>,
    /// `(position, HGVSp)` of high-confidence pathogenic entries (`ps1_set`).
    same_change: HashSet<(Rc<str>, Rc<str>)>,
    /// `(position, variant type)` of the same entries (`pm5_residues`).
    same_residue: HashSet<(Rc<str>, &'static str)>,
}

/// One record of the splice dict, as check_splice_pathogenic reads it.
#[derive(Clone, Debug)]
struct SpliceRecord {
    chrom: Rc<str>,
    pos: i64,
    /// `abs(len(ref) - len(alt))`, at least 1 (lengths in code points).
    length_change: i64,
    hgvsc: Rc<str>,
    exon: Rc<str>,
    intron: Rc<str>,
    clinvar_sig: Rc<str>,
    /// `splice_ai.get(name)` for DS_AG, DS_AL, DS_DG, DS_DL (`None`: absent).
    delta_scores: [Option<Rc<str>>; 4],
}

#[derive(Debug, Default)]
struct SpliceTranscript {
    gene_level_pathogenic: bool,
    records: Vec<SpliceRecord>,
}

/// The two ClinVar dictionaries, decoded for the transcripts of one table.
#[derive(Debug, Default)]
pub(crate) struct ClinvarKnownVariants {
    /// `len(clinvar_aa_dict)` and `len(clinvar_splice_dict)` (for the logs).
    pub(crate) amino_acid_transcripts: usize,
    pub(crate) splice_transcripts: usize,
    amino_acid: HashMap<String, AminoAcidTranscript>,
    splice: HashMap<String, SpliceTranscript>,
}

const DELTA_SCORE_NAMES: [&str; 4] = [
    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL",
];

fn odd_shape(file: &str, transcript: &str, what: &str) -> Box<dyn std::error::Error> {
    invalid(format!(
        "{file}: transcript {transcript}: {what}; the ClinVar dictionaries are read for the \
         installed shape only"
    ))
}

fn text_list(item: Option<&Item>) -> Option<Vec<Rc<str>>> {
    match item? {
        Item::List(items) => items
            .borrow()
            .iter()
            .map(|item| item.as_str().cloned())
            .collect(),
        _ => None,
    }
}

fn dict_of(item: &Item) -> Option<std::cell::Ref<'_, TextDict>> {
    match item {
        Item::Dict(dict) => Some(dict.borrow()),
        _ => None,
    }
}

impl ClinvarKnownVariants {
    /// Reads the amino-acid and splice pickles (gzip when the name ends in `.gz`), keeping the
    /// transcripts in `transcripts` (the table's `Feature` values as `astype(str)` gives them),
    /// or every transcript for `None`.
    pub(crate) fn load(
        clinvar_aa_dict_pkl: &str,
        clinvar_splice_dict_pkl: &str,
        transcripts: Option<&HashSet<String>>,
    ) -> Result<ClinvarKnownVariants> {
        let wanted = |transcript: &str| transcripts.is_none_or(|kept| kept.contains(transcript));
        let mut known = ClinvarKnownVariants::default();
        // Every amino-acid entry is decoded: PS1_PM5_criteria reads them all.
        known.amino_acid_transcripts = visit_top_level_dict(
            clinvar_aa_dict_pkl,
            |_| true,
            |transcript, value| {
                let keep = wanted(&transcript);
                let decoded =
                    decode_amino_acid_transcript(clinvar_aa_dict_pkl, &transcript, &value, keep)?;
                if keep {
                    known.amino_acid.insert(transcript.to_string(), decoded);
                }
                Ok(())
            },
        )?;
        known.splice_transcripts =
            visit_top_level_dict(clinvar_splice_dict_pkl, wanted, |transcript, value| {
                let decoded =
                    decode_splice_transcript(clinvar_splice_dict_pkl, &transcript, &value)?;
                known.splice.insert(transcript.to_string(), decoded);
                Ok(())
            })?;
        Ok(known)
    }

    /// Whether some record of `transcript` passes the gene summary's test, in the amino-acid
    /// dict (`None` when the transcript is not in it) and in the splice dict.
    pub(crate) fn gene_level_pathogenic(&self, transcript: &str) -> (Option<bool>, Option<bool>) {
        (
            self.amino_acid
                .get(transcript)
                .map(|found| found.gene_level_pathogenic),
            self.splice
                .get(transcript)
                .map(|found| found.gene_level_pathogenic),
        )
    }
}

/// One transcript of the amino-acid dict:
/// `{protein_pos: {hgvsp: {"CLNSIG": [str], "CLNREVSTAT": [str]}}}`. Every transcript is
/// checked (PS1_PM5_criteria reads them all); the facts are gathered only when `keep` is set.
fn decode_amino_acid_transcript(
    file: &str,
    transcript: &str,
    value: &Item,
    keep: bool,
) -> Result<AminoAcidTranscript> {
    let odd = |what: &str| odd_shape(file, transcript, what);
    let positions = dict_of(value).ok_or_else(|| odd("the value is not a dict"))?;
    let mut decoded = AminoAcidTranscript::default();
    for (position, variants) in &positions.entries {
        if keep {
            decoded.positions.insert(position.clone());
        }
        let variants = dict_of(variants).ok_or_else(|| odd("a position maps to a non-dict"))?;
        for (hgvsp, entry) in &variants.entries {
            let entry = dict_of(entry).ok_or_else(|| odd("an HGVSp maps to a non-dict"))?;
            let (Some(significances), Some(statuses)) = (
                text_list(entry.get("CLNSIG")),
                text_list(entry.get("CLNREVSTAT")),
            ) else {
                return Err(odd("an entry lacks CLNSIG or CLNREVSTAT lists of text"));
            };
            if !significances
                .iter()
                .chain(&statuses)
                .all(|text| text.is_ascii())
            {
                return Err(odd("an entry holds non-ASCII CLNSIG or CLNREVSTAT text"));
            }
            let pairs = || significances.iter().zip(&statuses);
            if !keep {
                // Only the HGVSp's type can fail (non-ASCII text), and only for these entries.
                let high_confidence = pairs().any(|(significance, status)| {
                    let stars = review_stars(status);
                    (significance.contains("Pathogenic") && stars >= 2)
                        || (&**significance == "Likely_pathogenic" && stars >= 3)
                });
                if high_confidence {
                    variant_type_of_hgvsp(hgvsp)?;
                }
                continue;
            }
            if pairs().any(|(significance, status)| gene_level_pathogenic(significance, status)) {
                decoded.gene_level_pathogenic = true;
            }
            let high_confidence = pairs().any(|(significance, status)| {
                let stars = review_stars(status);
                (significance.contains("Pathogenic") && stars >= 2)
                    || (&**significance == "Likely_pathogenic" && stars >= 3)
            });
            if high_confidence {
                decoded
                    .same_change
                    .insert((position.clone(), hgvsp.clone()));
                if let Some(kind) = variant_type_of_hgvsp(hgvsp)? {
                    decoded.same_residue.insert((position.clone(), kind));
                }
            }
        }
    }
    Ok(decoded)
}

/// One transcript of the splice dict: a list of records with the eleven fields
/// `chrom, pos, ref, alt, exon, intron, hgvsc, consequence, clinvar_sig, clinvar_review,
/// splice_ai`.
fn decode_splice_transcript(
    file: &str,
    transcript: &str,
    value: &Item,
) -> Result<SpliceTranscript> {
    let odd = |what: &str| odd_shape(file, transcript, what);
    let Item::List(records) = value else {
        return Err(odd("the value is not a list"));
    };
    let mut decoded = SpliceTranscript::default();
    for record in records.borrow().iter() {
        let record = dict_of(record).ok_or_else(|| odd("a record is not a dict"))?;
        let text = |name: &str| -> Result<Rc<str>> {
            match record.get(name) {
                Some(Item::Str(text)) if text.is_ascii() => Ok(text.clone()),
                _ => Err(odd(&format!("a record's {name} is not ASCII text"))),
            }
        };
        let pos = match record.get("pos") {
            Some(Item::Int(number)) => *number,
            _ => return Err(odd("a record's pos is not an int")),
        };
        let splice_ai = record
            .get("splice_ai")
            .and_then(dict_of)
            .ok_or_else(|| odd("a record's splice_ai is not a dict"))?;
        let mut delta_scores: [Option<Rc<str>>; 4] = Default::default();
        for (slot, name) in delta_scores.iter_mut().zip(DELTA_SCORE_NAMES) {
            *slot = match splice_ai.get(name) {
                None => None,
                Some(Item::Str(text)) if text.is_ascii() => Some(text.clone()),
                Some(_) => return Err(odd("a SpliceAI delta score is not ASCII text")),
            };
        }
        let (reference, alternative) = (text("ref")?, text("alt")?);
        let change = (reference.len() as i64 - alternative.len() as i64).abs();
        let (clinvar_sig, clinvar_review) = (text("clinvar_sig")?, text("clinvar_review")?);
        text("consequence")?;
        if record.len() != 11 {
            return Err(odd("a record has other fields than the eleven expected"));
        }
        if gene_level_pathogenic(&clinvar_sig, &clinvar_review) {
            decoded.gene_level_pathogenic = true;
        }
        decoded.records.push(SpliceRecord {
            chrom: text("chrom")?,
            pos,
            length_change: if change == 0 { 1 } else { change },
            hgvsc: text("hgvsc")?,
            exon: text("exon")?,
            intron: text("intron")?,
            clinvar_sig,
            delta_scores,
        });
    }
    Ok(decoded)
}

// ------------------------------------------------------------------ get_variant_type

/// `\w` of Python's `re` for ASCII text.
fn is_word(byte: u8) -> bool {
    byte.is_ascii_alphanumeric() || byte == b'_'
}

/// `re.search(r'p\.\w+\d+' + suffix, text)` for ASCII text: some `p.` is followed by two or
/// more word characters whose last one is a digit, and then by the suffix.
fn search_word_digits_then(text: &[u8], suffix: impl Fn(&[u8]) -> bool) -> bool {
    (0..text.len().saturating_sub(1))
        .filter(|start| &text[*start..*start + 2] == b"p.")
        .any(|start| {
            let body = start + 2;
            let run_end = body
                + text[body..]
                    .iter()
                    .take_while(|byte| is_word(**byte))
                    .count();
            // The suffix starts at `at`; `\w+\d+` covers text[body..at].
            (body + 2..=run_end).any(|at| text[at - 1].is_ascii_digit() && suffix(&text[at..]))
        })
}

/// `re.search(r'p\.[A-Z][a-z]+\d+[A-Z][a-z]+$', text)` for ASCII text (`$` also matches before
/// a final line break).
fn search_missense(text: &[u8]) -> bool {
    let end = if text.ends_with(b"\n") {
        text.len() - 1
    } else {
        text.len()
    };
    (0..end.saturating_sub(1))
        .filter(|start| &text[*start..*start + 2] == b"p.")
        .any(|start| {
            let rest = &text[start + 2..end];
            let mut index = 0;
            let mut run = |test: fn(&u8) -> bool, at_least: usize| -> bool {
                let count = rest[index..].iter().take_while(|byte| test(byte)).count();
                index += count;
                count >= at_least
            };
            run(u8::is_ascii_uppercase, 1)
                && run(u8::is_ascii_lowercase, 1)
                && run(u8::is_ascii_digit, 1)
                && run(u8::is_ascii_uppercase, 1)
                && run(u8::is_ascii_lowercase, 1)
                && index == rest.len()
        })
}

/// `get_variant_type(hgvsp)` (acmg_consequence.py) for a text: `nonsense`, `frameshift`,
/// `delins`, `deletion`, `insertion_type`, `missense` or `None`. Copied until thread T1's
/// `acmg_consequence` module lands. Non-ASCII text is refused (Python's `\w` and `\d` would
/// read letters and digits of other scripts).
pub(crate) fn variant_type_of_hgvsp(hgvsp: &str) -> Result<Option<&'static str>> {
    if !hgvsp.is_ascii() {
        return Err(invalid(format!(
            "HGVSp {hgvsp:?} holds non-ASCII text; Python's \\w and \\d are not reproduced there"
        )));
    }
    let text = hgvsp.rsplit(':').next().unwrap_or("").as_bytes();
    let starts = |prefix: &'static [u8]| move |rest: &[u8]| rest.starts_with(prefix);
    let kind = if search_word_digits_then(text, starts(b"Ter"))
        || search_word_digits_then(text, starts(b"*"))
    {
        Some("nonsense")
    } else if search_word_digits_then(text, starts(b"fs"))
        || search_word_digits_then(text, |rest| {
            rest.len() >= 6
                && matches!(rest[0], b'F' | b'f')
                && rest[1] == b's'
                && matches!(rest[2], b'T' | b't')
                && &rest[3..5] == b"er"
                && rest[5].is_ascii_digit()
        })
    {
        Some("frameshift")
    } else if search_word_digits_then(text, starts(b"delins")) {
        Some("delins")
    } else if search_word_digits_then(text, starts(b"del")) {
        Some("deletion")
    } else if search_word_digits_then(text, starts(b"dup"))
        || search_word_digits_then(text, starts(b"ins"))
    {
        Some("insertion_type")
    } else if search_missense(text) {
        Some("missense")
    } else {
        None
    };
    Ok(kind)
}

// ------------------------------------------------------------------ check_splice_pathogenic

/// A value of the row dict `df.iloc[...].to_dict('records')` hands check_splice_pathogenic.
fn row_value<'a>(table: &'a Table, name: &str, row: usize) -> Option<Cell<'a>> {
    table
        .column_index(name)
        .map(|index| table.columns[index].cell(row))
}

/// Whether that value is Python's `None` (an `object` column can hold it; `Cell::Missing`
/// does not tell it from NaN).
fn python_none(table: &Table, name: &str, row: usize) -> bool {
    match table.column_index(name).map(|index| &table.columns[index]) {
        Some(Column::Object(values)) => matches!(values[row], ObjectValue::None),
        _ => false,
    }
}

fn key_error(name: &str) -> Box<dyn std::error::Error> {
    Box::new(AcmgError::new("KeyError", format!("'{name}'")))
}

/// `abs(value)` of a SpliceAI delta score cell (NaN for an absent column).
fn delta_score(cell: Option<Cell<'_>>, name: &str) -> Result<f64> {
    match cell {
        None | Some(Cell::Missing) => Ok(f64::NAN),
        Some(Cell::Float(value)) => Ok(value.abs()),
        Some(Cell::Int(value)) => Ok((value as f64).abs()),
        Some(Cell::Bool(flag)) => Ok(f64::from(u8::from(flag))),
        Some(Cell::Text(_)) => Err(Box::new(AcmgError::new(
            "TypeError",
            format!("bad operand type for abs(): 'str' ({name})"),
        ))),
    }
}

/// `str(value)`.
fn str_of(cell: Cell<'_>) -> String {
    python_str(cell)
}

/// `len(value)` of `row["ref"]` / `row["alt"]`.
fn text_length(cell: Cell<'_>, name: &str) -> Result<i64> {
    match cell {
        Cell::Text(text) => Ok(text.chars().count() as i64),
        _ => Err(Box::new(AcmgError::new(
            "TypeError",
            format!(
                "object of type '{}' has no len() ({name})",
                python_type(cell)
            ),
        ))),
    }
}

fn python_type(cell: Cell<'_>) -> &'static str {
    match cell {
        Cell::Missing | Cell::Float(_) => "float",
        Cell::Text(_) => "str",
        Cell::Int(_) => "int",
        Cell::Bool(_) => "bool",
    }
}

/// `int(row["pos"])`.
fn int_of(cell: Cell<'_>, name: &str) -> Result<i64> {
    match cell {
        Cell::Int(value) => Ok(value),
        Cell::Bool(flag) => Ok(i64::from(flag)),
        Cell::Float(value) if value.is_finite() && value.abs() < 9.0e18 => Ok(value.trunc() as i64),
        Cell::Missing => Err(Box::new(AcmgError::new(
            "ValueError",
            format!("cannot convert float NaN to integer ({name})"),
        ))),
        _ => Err(invalid(format!(
            "int() of {name} = {cell:?} is not reproduced"
        ))),
    }
}

/// The Python exception a group 2 helper reports as `"<type>: <message>"`; a message of
/// another form is a refusal.
pub(super) fn python_exception(message: &str) -> Box<dyn std::error::Error> {
    let Some((kind, text)) = message.split_once(": ") else {
        return invalid(message.to_owned());
    };
    let kind: &'static str = match kind {
        "ValueError" => "ValueError",
        "TypeError" => "TypeError",
        "IndexError" => "IndexError",
        "KeyError" => "KeyError",
        "AttributeError" => "AttributeError",
        "ZeroDivisionError" => "ZeroDivisionError",
        _ => return invalid(message.to_owned()),
    };
    Box::new(AcmgError::new(kind, text.to_owned()))
}

/// `parse_hgvsc_splice_position(...)["overlapping_canonical_site"]`, False for `None`.
fn overlaps_canonical_site(hgvsc: &str, strand: &str, pos: i64, length: i64) -> Result<bool> {
    match parse_hgvsc_splice_position(Some(hgvsc), strand, pos, length) {
        Ok(found) => Ok(found.is_some_and(|position| position.overlapping_canonical_site)),
        Err(message) => Err(python_exception(&message)),
    }
}

/// What check_splice_pathogenic returns: `False` or one of the PS1 strengths.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum SpliceMatch {
    NoMatch,
    Ps1,
    Ps1Moderate,
    Ps1Supporting,
}

impl SpliceMatch {
    fn python_text(self) -> &'static str {
        match self {
            SpliceMatch::NoMatch => "False",
            SpliceMatch::Ps1 => "PS1",
            SpliceMatch::Ps1Moderate => "PS1_Moderate",
            SpliceMatch::Ps1Supporting => "PS1_Supporting",
        }
    }
}

/// The equality `patho['exon'] == var_exon` of a ClinVar text with a row value.
fn equals_text(record_text: &str, cell: Cell<'_>) -> bool {
    matches!(cell, Cell::Text(text) if text == record_text)
}

/// `check_splice_pathogenic(row, records, pvs1_strength, spliceai_cutoff)` for row `row` of
/// the table and the splice records of its transcript.
fn check_splice_pathogenic(
    table: &Table,
    row: usize,
    records: &[SpliceRecord],
    pvs1_strength: i64,
    spliceai_cutoff: f64,
) -> Result<SpliceMatch> {
    if records.is_empty() {
        return Ok(SpliceMatch::NoMatch);
    }
    let (Some(var_intron), Some(var_exon)) = (
        row_value(table, "INTRON", row),
        row_value(table, "EXON", row),
    ) else {
        return Ok(SpliceMatch::NoMatch);
    };
    if python_none(table, "INTRON", row) || python_none(table, "EXON", row) {
        return Ok(SpliceMatch::NoMatch);
    }
    let matching: Vec<&SpliceRecord> = records
        .iter()
        .filter(|record| {
            equals_text(&record.exon, var_exon) || equals_text(&record.intron, var_intron)
        })
        .collect();
    if matching.is_empty() {
        return Ok(SpliceMatch::NoMatch);
    }
    let deltas = DELTA_SCORE_NAMES
        .iter()
        .map(|name| {
            if python_none(table, name, row) {
                return Err(Box::new(AcmgError::new(
                    "TypeError",
                    format!("bad operand type for abs(): 'NoneType' ({name})"),
                )) as Box<dyn std::error::Error>);
            }
            delta_score(row_value(table, name, row), name)
        })
        .collect::<Result<Vec<f64>>>()?;
    if deltas.iter().all(|delta| *delta < spliceai_cutoff) {
        return Ok(SpliceMatch::NoMatch);
    }
    // `np.argmax`: the first NaN if there is one, else the first largest value.
    let largest = match deltas.iter().position(|delta| delta.is_nan()) {
        Some(index) => index,
        None => {
            let mut best = 0;
            for (index, delta) in deltas.iter().enumerate() {
                if *delta > deltas[best] {
                    best = index;
                }
            }
            best
        }
    };
    let max_delta_score = deltas[largest];
    for record in matching {
        let record_delta = match &record.delta_scores[largest] {
            None => f64::NAN,
            Some(text) => parse_float(text).ok_or_else(|| {
                Box::new(AcmgError::new(
                    "ValueError",
                    format!("could not convert string to float: {text:?}"),
                )) as Box<dyn std::error::Error>
            })?,
        };
        if !(record_delta.abs() <= max_delta_score) {
            continue;
        }
        let cell = |name: &str| row_value(table, name, row).ok_or_else(|| key_error(name));
        let text_of = |name: &str| -> Result<String> {
            let value = cell(name)?;
            Ok(if python_none(table, name, row) {
                "None".to_owned()
            } else {
                str_of(value)
            })
        };
        let vua_hgvsc = text_of("HGVSc")?;
        let ref_length = text_length(cell("ref")?, "ref")?;
        let alt_length = text_length(cell("alt")?, "alt")?;
        let vua_length = match (ref_length - alt_length).abs() {
            0 => 1,
            change => change,
        };
        let strand = text_of("STRAND")?;
        if python_none(table, "pos", row) {
            return Err(Box::new(AcmgError::new(
                "TypeError",
                "int() argument must be a string, a bytes-like object or a real number, not 'NoneType'",
            )));
        }
        let row_pos = int_of(cell("pos")?, "pos")?;
        let vua_overlapping = overlaps_canonical_site(&vua_hgvsc, &strand, row_pos, vua_length)?;
        let patho_overlapping =
            overlaps_canonical_site(&record.hgvsc, &strand, record.pos, record.length_change)?;
        let patho_likely = record.clinvar_sig.contains("ikely");
        return Ok(if vua_overlapping {
            if pvs1_strength >= 4 {
                if patho_overlapping {
                    SpliceMatch::NoMatch
                } else {
                    SpliceMatch::Ps1Supporting
                }
            } else if pvs1_strength >= 1 {
                match (patho_overlapping, patho_likely) {
                    (true, true) => SpliceMatch::NoMatch,
                    (true, false) => SpliceMatch::Ps1,
                    (false, true) => SpliceMatch::Ps1Supporting,
                    (false, false) => SpliceMatch::Ps1Moderate,
                }
            } else {
                SpliceMatch::NoMatch
            }
        } else {
            let same_site = record.pos == row_pos && text_of("chrom")? == *record.chrom;
            match (same_site, patho_likely) {
                (true, true) => SpliceMatch::Ps1Moderate,
                (true, false) => SpliceMatch::Ps1,
                (false, true) => SpliceMatch::Ps1Supporting,
                (false, false) => SpliceMatch::Ps1Moderate,
            }
        });
    }
    Ok(SpliceMatch::NoMatch)
}

// ------------------------------------------------------------------ PS1_PM5_criteria

/// `df[name].astype(str)` for the four columns PS1_PM5_criteria keys on.
fn str_column(table: &Table, name: &str) -> Result<Vec<String>> {
    let column = table.column(name)?;
    if matches!(column, Column::NullableInt(_) | Column::NullableBool(_)) {
        return Err(invalid(format!(
            "{name} is a nullable pandas column; astype(str) is not reproduced there"
        )));
    }
    if let Column::Object(values) = column
        && values
            .iter()
            .any(|value| matches!(value, ObjectValue::None))
    {
        return Err(invalid(format!(
            "{name} holds None; astype(str) gives 'None' there, which is not reproduced"
        )));
    }
    Ok((0..table.rows)
        .map(|row| str_of(column.cell(row)))
        .collect())
}

/// The transcripts a table can look up: `df['Feature'].astype(str)`.
pub(crate) fn table_transcripts(table: &Table) -> Result<HashSet<String>> {
    Ok(str_column(table, "Feature")?.into_iter().collect())
}

/// `PS1_PM5_criteria(df, clinvar_aa_dict, clinvar_splice_dict, ps3_clinvar_patho,
/// pvs1_criteria, threads, spliceai_cutoff)`: `(PS1, PM5)`, int64 strengths.
pub(crate) fn ps1_pm5_criteria(
    table: &Table,
    known: &ClinvarKnownVariants,
    pvs1_criteria: &[i64],
    threads: i64,
    spliceai_cutoff: f64,
) -> Result<StepOutput> {
    let rows = table.rows;
    if pvs1_criteria.len() != rows {
        return Err(invalid(format!(
            "pvs1_criteria has {} values for {rows} rows",
            pvs1_criteria.len()
        )));
    }
    let feature = str_column(table, "Feature")?;
    let raw_pp = str_column(table, "Protein_position")?;
    let hgvsp = str_column(table, "HGVSp")?;
    let consequence = str_column(table, "Consequence")?;
    let no_hgvsp: Vec<bool> = hgvsp
        .iter()
        .map(|text| matches!(text.as_str(), "nan" | "inf" | ""))
        .collect();

    let mut ps1_eligible = vec![false; rows];
    let mut pm5_eligible = vec![false; rows];
    let mut splice_fallback = vec![false; rows];
    for row in 0..rows {
        let found = known.amino_acid.get(&feature[row]);
        let position: Rc<str> = Rc::from(raw_pp[row].as_str());
        let change: Rc<str> = Rc::from(hgvsp[row].as_str());
        ps1_eligible[row] =
            found.is_some_and(|found| found.same_change.contains(&(position.clone(), change)));
        let kind = if no_hgvsp[row] {
            None
        } else {
            variant_type_of_hgvsp(&hgvsp[row])?
        };
        let residue = match (found, kind) {
            (Some(found), Some(kind)) => found.same_residue.contains(&(position.clone(), kind)),
            _ => false,
        };
        // `str.contains('synonymous')` on the str values.
        let synonymous = consequence[row].contains("synonymous");
        pm5_eligible[row] = residue && !synonymous && !ps1_eligible[row];
        let in_clinvar = found.is_some() || known.splice.contains_key(&feature[row]);
        let position_known = found.is_some_and(|found| found.positions.contains(&position));
        splice_fallback[row] = in_clinvar && (no_hgvsp[row] || !position_known);
    }

    let fallback_rows: Vec<usize> = (0..rows).filter(|row| splice_fallback[*row]).collect();
    let mut splice_results = vec![SpliceMatch::NoMatch; rows];
    if !fallback_rows.is_empty() {
        // `chunk_size = max(len(args) // (threads * 4), 1)`, then `mp.Pool(threads)`.
        if threads == 0 {
            return Err(Box::new(AcmgError::new(
                "ZeroDivisionError",
                "integer division or modulo by zero",
            )));
        }
        if threads < 1 {
            return Err(Box::new(AcmgError::new(
                "ValueError",
                "Number of processes must be at least 1",
            )));
        }
        let feature_column = table.column("Feature")?;
        for row in fallback_rows {
            // `clinvar_splice_dict.get(row.get('Feature', ''), [])` with the raw value.
            let records = match feature_column.cell(row) {
                Cell::Text(text) => known
                    .splice
                    .get(text)
                    .map(|found| found.records.as_slice())
                    .unwrap_or(&[]),
                _ => &[],
            };
            splice_results[row] =
                check_splice_pathogenic(table, row, records, pvs1_criteria[row], spliceai_cutoff)?;
        }
    }

    // `df['CLNSIG'].fillna("").str.contains('pathogenic') & (df['CLNREVSTAT'].map(...) >= 2)`.
    let clnsig = table.column("CLNSIG")?;
    let clnsig_pathogenic = (0..rows)
        .map(|row| match clnsig.cell(row) {
            Cell::Text(text) => Ok(text.contains("pathogenic")),
            Cell::Missing => Ok(false),
            _ => Err(invalid(
                "CLNSIG holds values other than text; fillna(\"\").str.contains is not reproduced",
            )),
        })
        .collect::<Result<Vec<bool>>>()?;
    if matches!(clnsig, Column::NullableBool(_) | Column::NullableInt(_)) {
        return Err(invalid(
            "CLNSIG is a nullable pandas column; not reproduced",
        ));
    }
    let high_review = high_confidence_review(table, "CLNREVSTAT")?;

    let mut ps1 = vec![0_i64; rows];
    let mut pm5 = vec![0_i64; rows];
    for row in 0..rows {
        let splice = splice_results[row];
        if ps1_eligible[row] || splice == SpliceMatch::Ps1 {
            ps1[row] = 3;
        } else if splice == SpliceMatch::Ps1Moderate {
            ps1[row] = 2;
        } else if splice == SpliceMatch::Ps1Supporting {
            ps1[row] = 1;
        }
        if clnsig_pathogenic[row] && high_review[row] {
            ps1[row] = 0;
        }
        if pm5_eligible[row] && ps1[row] == 0 {
            pm5[row] = 2;
        }
    }
    Ok(StepOutput {
        criteria: vec![("PS1", Column::Int(ps1)), ("PM5", Column::Int(pm5))],
        ..StepOutput::default()
    })
}

/// `priva acmg-step acmg_ps1_pm5_known_variant <function>`:
///
/// - `PS1_PM5_criteria`, options `clinvar_aa_dict_pkl`, `clinvar_splice_dict_pkl`,
///   `pvs1_criteria` (list of ints), `threads` (default 10), `spliceai_cutoff` (default 0.5):
///   a frame of `PS1`, `PM5` (int64);
/// - `check_splice_pathogenic`, options `clinvar_splice_dict_pkl` (and an amino-acid pickle,
///   read but unused), `pvs1_criteria`, `spliceai_cutoff`: a frame of `result` (text: `False`,
///   `PS1`, `PS1_Moderate` or `PS1_Supporting`), one per row, each row given the records of its
///   `Feature`.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let table = &call.table;
    let options = &call.options;
    let pvs1 =
        harness_options::int_list(options, "pvs1_criteria")?.unwrap_or_else(|| vec![0; table.rows]);
    let cutoff = harness_options::float(options, "spliceai_cutoff")?.unwrap_or(0.5);
    let load = || -> Result<ClinvarKnownVariants> {
        ClinvarKnownVariants::load(
            harness_options::required_text(options, "clinvar_aa_dict_pkl")?,
            harness_options::required_text(options, "clinvar_splice_dict_pkl")?,
            Some(&table_transcripts(table)?),
        )
    };
    match function {
        "PS1_PM5_criteria" => {
            let known = load()?;
            let threads = harness_options::float(options, "threads")?.unwrap_or(10.0) as i64;
            let output = ps1_pm5_criteria(table, &known, &pvs1, threads, cutoff)?;
            Ok(harness_options::criteria_frame(&output, table.rows))
        }
        "check_splice_pathogenic" => {
            let known = load()?;
            if pvs1.len() != table.rows {
                return Err(invalid("pvs1_criteria must give one value per row"));
            }
            let feature = table.column("Feature")?;
            let results = (0..table.rows)
                .map(|row| {
                    let records = match feature.cell(row) {
                        Cell::Text(text) => known
                            .splice
                            .get(text)
                            .map(|found| found.records.as_slice())
                            .unwrap_or(&[]),
                        _ => &[],
                    };
                    check_splice_pathogenic(table, row, records, pvs1[row], cutoff)
                        .map(|found| Some(found.python_text().to_owned()))
                })
                .collect::<Result<Vec<_>>>()?;
            Ok(harness_options::frame(
                vec![("result".to_owned(), Column::Text(results))],
                table.rows,
            ))
        }
        "load_known_variants" => {
            // How much reading the two dictionaries takes: every transcript kept with option
            // `keep_all`, else the table's.
            let keep_all = matches!(
                options.get("keep_all"),
                Some(crate::json_text::Json::Bool(true))
            );
            let transcripts = table_transcripts(table)?;
            let known = ClinvarKnownVariants::load(
                harness_options::required_text(options, "clinvar_aa_dict_pkl")?,
                harness_options::required_text(options, "clinvar_splice_dict_pkl")?,
                (!keep_all).then_some(&transcripts),
            )?;
            let counts = [
                known.amino_acid_transcripts,
                known.splice_transcripts,
                known.amino_acid.len(),
                known.splice.len(),
                known.splice.values().map(|found| found.records.len()).sum(),
            ];
            Ok(harness_options::frame(
                vec![(
                    "count".to_owned(),
                    Column::Int(counts.iter().map(|count| *count as i64).collect()),
                )],
                counts.len(),
            ))
        }
        _ => Err(not_converted("acmg_ps1_pm5_known_variant", function)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn variant_types_follow_the_python_patterns() {
        let kind = |text: &str| variant_type_of_hgvsp(text).unwrap();
        assert_eq!(kind("ENSP1.1:p.Arg123Ter"), Some("nonsense"));
        assert_eq!(kind("p.Arg123*"), Some("nonsense"));
        assert_eq!(kind("p.R1Ter"), Some("nonsense"));
        // `\w+\d+` needs two characters, the last a digit.
        assert_eq!(kind("p.1Ter"), None);
        // The patterns need a digit right before "fs", so VEP's usual frameshift HGVSp
        // (an amino acid between the position and "fs") is no frameshift to Python either.
        assert_eq!(kind("p.Asn3264LeufsTer12"), None);
        assert_eq!(kind("p.Gly123fs"), Some("frameshift"));
        assert_eq!(kind("p.Gly123FsTer5"), Some("frameshift"));
        assert_eq!(kind("p.Gly123FsTer"), None);
        assert_eq!(kind("p.Ala123delinsGly"), Some("delins"));
        assert_eq!(kind("p.Ala123_Gly125del"), Some("deletion"));
        assert_eq!(kind("p.Ala123dup"), Some("insertion_type"));
        assert_eq!(kind("p.Ala123_Gly124insTrp"), Some("insertion_type"));
        assert_eq!(kind("p.Ala123Gly"), Some("missense"));
        assert_eq!(kind("p.Ala123Gly\n"), Some("missense"));
        assert_eq!(kind("p.Ala123Gly\n\n"), None);
        assert_eq!(kind("p.Ala123="), None);
        assert_eq!(kind("p.Leu12%3D"), None);
        assert_eq!(kind("x:p.Ala123Gly:p.Arg5Ter"), Some("nonsense"));
        assert_eq!(kind("p.Ala123Gly p.Arg1Cys"), Some("missense"));
        assert!(variant_type_of_hgvsp("p.Ala١٢٣Gly").is_err());
    }

    #[test]
    fn gene_level_test_tokenizes_clnsig() {
        let two = "criteria_provided,_multiple_submitters,_no_conflicts";
        assert!(gene_level_pathogenic("Pathogenic/Likely_pathogenic", two));
        assert!(gene_level_pathogenic(" likely_PATHOGENIC |x", two));
        assert!(!gene_level_pathogenic(
            "Conflicting_classifications_of_pathogenicity",
            two
        ));
        assert!(!gene_level_pathogenic(
            "Pathogenic",
            "criteria_provided,_single_submitter"
        ));
        // The comma splits off "Pathogenic", as Python's re.split does.
        assert!(gene_level_pathogenic("Pathogenic,_low_penetrance", two));
        assert!(!gene_level_pathogenic("Benign/Likely_benign", two));
    }
}
