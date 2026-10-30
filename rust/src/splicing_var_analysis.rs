//! Rust version of `scripts/splicing_var_analysis.py` (`priva splicing-analysis`).
//!
//! Reads the annotation table, judges every row's effect on splicing from
//! its SpliceVault top events, its SpliceAI delta scores and its HGVSc
//! position, adds nine columns, and replaces the table with the result:
//! `splicing_lof`, `splicing_len_changing`, `splicing_frameshift`,
//! `splicing_span_intol_domain`, `splicing_affected_exons`,
//! `splicing_ten_percent_protein`, `spliceai_canonical_gain_offset`,
//! `spliceai_canonical_loss_offset` and `splicing_last_exon_only`. A column
//! that already exists is replaced where it stands.
//!
//! The table is read and written the way pandas does it (`pandas_table`): each
//! column takes its type from its values, so a value is written back as pandas
//! prints it (`0.10` becomes `0.1`, an integer column with a gap becomes
//! `5.0`) and the type of a cell decides what the row logic does with it, as
//! in Python (`float()` of text, `int()` of a float, text where a number is
//! expected, ...). Every place where Python would raise stops the whole run
//! with exit code 1 and leaves the table untouched, as an exception inside
//! the script's worker pool does.
//!
//! Differences from the Python script:
//! - `splicing_affected_exons` lists the exons in ascending order; the script
//!   joins a set of strings, whose order changes from run to run.
//! - The script's INFO log lines per row are not written.
//! - The pickles are read by the restricted reader (`pickle_reader`), and an
//!   InterPro, transcript-map or intolerant-domain file of an unexpected shape
//!   stops the run when it is loaded rather than when a row first reaches it.
//! - The new table replaces the old one by a rename in the same folder; the
//!   script runs `mv` through a shell, which also mishandles paths with
//!   spaces.
//! - `na_value` is a private copy of the function of that name in
//!   `combine_annotations.py` (another thread converts that script; the two
//!   copies are to be merged later).
use crate::cli::{Result, invalid};
use crate::interpro_mapping_io::{InterproMapping, load_interpro_mapping};
use crate::pandas_table::{Cell, Column, Table, read_table, write_table};
use crate::pickle_reader::{self, Value};
use crate::protein_domain_mapping::{DomainNormalizer, Functionality};
use crate::python_text::{float_repr, parse_float, parse_int, sum_floats};
use crate::splicing_cutoffs::{Cutoffs, cutoffs_for_config};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::path::Path;

const USAGE: &str = "usage: priva splicing-analysis --anno_table ANNO_TABLE --transcript_domain_map PICKLE --intolerant_domains PICKLE --interpro_entry_map_pkl PICKLE [--threads THREADS] [--config CONFIG]";
const LONG_OPTIONS: [&str; 7] = [
    "--anno_table",
    "--transcript_domain_map",
    "--intolerant_domains",
    "--interpro_entry_map_pkl",
    "--threads",
    "--config",
    "--help",
];

/// A failure while judging rows, which stops the run like a Python exception.
type Outcome<T> = std::result::Result<T, String>;

// ------------------------------------------------------------ value rules

/// `na_value` of `combine_annotations.py` (see the module comment).
fn na_value(cell: &Cell) -> bool {
    match cell {
        Cell::Missing => true,
        Cell::Text(text) => [
            "NaN", "nan", "na", "NA", "NAN", "None", "none", "", ".", "-",
        ]
        .contains(text),
        Cell::Float(number) => number.is_nan(),
        Cell::Int(_) | Cell::Bool(_) => false,
    }
}

/// `str(value)` of a table value.
fn cell_text(cell: &Cell) -> String {
    match cell {
        Cell::Missing => "nan".to_owned(),
        Cell::Text(text) => (*text).to_owned(),
        Cell::Int(number) => number.to_string(),
        Cell::Float(number) => float_repr(*number),
        Cell::Bool(flag) => (if *flag { "True" } else { "False" }).to_owned(),
    }
}

/// `float(value)`.
fn cell_float(cell: &Cell) -> Outcome<f64> {
    match cell {
        Cell::Missing => Ok(f64::NAN),
        Cell::Int(number) => Ok(*number as f64),
        Cell::Float(number) => Ok(*number),
        Cell::Bool(flag) => Ok(if *flag { 1.0 } else { 0.0 }),
        Cell::Text(text) => parse_float(text)
            .ok_or_else(|| format!("ValueError: could not convert {text:?} to float")),
    }
}

fn text_int(text: &str) -> Outcome<i64> {
    parse_int(text)
        .ok_or_else(|| format!("ValueError: invalid literal for int() with base 10: {text:?}"))
}

/// `int(value)`.
fn cell_int(cell: &Cell) -> Outcome<i64> {
    match cell {
        Cell::Missing => Err("ValueError: cannot convert float NaN to integer".to_owned()),
        Cell::Int(number) => Ok(*number),
        Cell::Float(number) => {
            if number.is_finite() && number.abs() < 9.0e18 {
                Ok(number.trunc() as i64)
            } else {
                Err("cannot convert this float to an integer".to_owned())
            }
        }
        Cell::Bool(flag) => Ok(i64::from(*flag)),
        Cell::Text(text) => text_int(text),
    }
}

/// An INTRON, EXON or CDS_position value after `None if na_value(x) else x`.
#[derive(Clone, Copy)]
enum Position<'a> {
    Absent,
    Text(&'a str),
    /// A number or bool where the script expects text; only its truth value
    /// is used until something calls a string method on it.
    Other {
        truthy: bool,
    },
}

impl<'a> Position<'a> {
    fn from_cell(cell: Cell<'a>) -> Position<'a> {
        if na_value(&cell) {
            return Position::Absent;
        }
        match cell {
            Cell::Text(text) => Position::Text(text),
            Cell::Int(number) => Position::Other {
                truthy: number != 0,
            },
            Cell::Float(number) => Position::Other {
                truthy: number != 0.0,
            },
            Cell::Bool(flag) => Position::Other { truthy: flag },
            Cell::Missing => Position::Absent,
        }
    }

    fn truthy(&self) -> bool {
        match self {
            Position::Absent => false,
            Position::Text(text) => !text.is_empty(),
            Position::Other { truthy } => *truthy,
        }
    }

    /// The value as text, for a string method; Python raises AttributeError otherwise.
    fn text(&self) -> Outcome<&'a str> {
        match self {
            Position::Text(text) => Ok(*text),
            _ => Err("AttributeError: this position is not text".to_owned()),
        }
    }
}

/// `int(text.split("/")[0])`.
fn first_number(text: &str) -> Outcome<i64> {
    text_int(text.split('/').next().unwrap_or(""))
}

/// `map(int, text.split("/"))` unpacked into two numbers.
fn number_pair(text: &str) -> Outcome<(i64, i64)> {
    let parts: Vec<&str> = text.split('/').collect();
    if parts.len() != 2 {
        return Err(format!(
            "ValueError: {text:?} does not split into two numbers"
        ));
    }
    Ok((text_int(parts[0])?, text_int(parts[1])?))
}

/// The two numbers of INTRON or EXON for the event functions; `None` when absent.
fn pair_of(position: &Position) -> Outcome<Option<(i64, i64)>> {
    match position {
        Position::Absent => Ok(None),
        Position::Text(text) => Ok(Some(number_pair(text)?)),
        Position::Other { .. } => Err("AttributeError: this position is not text".to_owned()),
    }
}

/// `int(cds_pos.split("/")[1])` for the SpliceAI rules.
fn cds_total(cds_pos: &Position) -> Outcome<i64> {
    let text = cds_pos.text()?;
    let part = text
        .split('/')
        .nth(1)
        .ok_or_else(|| format!("IndexError: {text:?} has no total"))?;
    text_int(part)
}

/// The total length from `cds_pos.split("/")` in the cryptic donor function;
/// `lenient_start` is the acceptor function's `try/except ValueError` around
/// the first number of a range.
fn cds_total_of_event(cds_pos: &Position, lenient_start: bool) -> Outcome<Option<i64>> {
    let text = match cds_pos {
        Position::Absent => return Ok(None),
        Position::Text(text) => *text,
        Position::Other { .. } => {
            return Err("AttributeError: this position is not text".to_owned());
        }
    };
    let parts: Vec<&str> = text.split('/').collect();
    if parts.len() != 2 {
        return Err(format!(
            "ValueError: {text:?} does not split into two parts"
        ));
    }
    let (start, total) = (parts[0], parts[1]);
    if start.contains('-') {
        let first = start.split('-').next().unwrap_or("");
        if !lenient_start && !first.is_empty() {
            text_int(first)?;
        }
    } else {
        text_int(start)?;
    }
    Ok(Some(text_int(total)?))
}

/// `abs(offset) / total > 0.1` for a truthy total.
fn exceeds_a_tenth(offset: f64, total: i64) -> Outcome<bool> {
    if total == 0 {
        return Err("ZeroDivisionError: division by zero".to_owned());
    }
    Ok(offset.abs() / (total as f64) > 0.1)
}

// ----------------------------------------------------------------- HGVSc

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum VariantType {
    Unknown,
    DonorRegion,
    AcceptorRegion,
    Spanning,
}

/// What `parse_hgvsc_splice_position` returns; the distances are signed base
/// counts from the variant to the first base of the canonical splice site.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct SplicePosition {
    pub(crate) variant_type: VariantType,
    pub(crate) distance_to_canonical: Option<i128>,
    // The fields below are for the known-variant step (PS1/PM5), which uses
    // this function after its own conversion.
    #[allow(dead_code)]
    pub(crate) exon_boundary_pos: Option<i128>,
    #[allow(dead_code)]
    pub(crate) exon_skipping: bool,
    #[allow(dead_code)]
    pub(crate) overlapping_canonical_site: bool,
}

#[derive(Clone, Copy)]
struct SignedNumber {
    negative: bool,
    digits_start: usize,
    digits_end: usize,
}

fn signed_value(text: &[u8], number: SignedNumber) -> Outcome<i128> {
    let digits = &text[number.digits_start..number.digits_end];
    if digits.len() > 18 {
        return Err("an HGVSc offset has more than 18 digits".to_owned());
    }
    let magnitude: i128 = std::str::from_utf8(digits)
        .ok()
        .and_then(|digits| digits.parse().ok())
        .ok_or_else(|| "an HGVSc offset is not a number".to_owned())?;
    Ok(if number.negative {
        -magnitude
    } else {
        magnitude
    })
}

/// The ways `\d+` can match at `start`, longest first.
fn digit_run_ends(text: &[u8], start: usize) -> Vec<usize> {
    let mut end = start;
    while end < text.len() && text[end].is_ascii_digit() {
        end += 1;
    }
    (start + 1..=end).rev().collect()
}

/// The ways `(\d+)?` can match at `start`, in the order a regex tries them.
fn optional_digits(text: &[u8], start: usize) -> Vec<Option<usize>> {
    let mut ways: Vec<Option<usize>> = digit_run_ends(text, start).into_iter().map(Some).collect();
    ways.push(None);
    ways
}

/// The ways `(?:(\+|-)(\d+))?` can match at `start`, in the order a regex tries them.
fn optional_signed(text: &[u8], start: usize) -> Vec<Option<SignedNumber>> {
    let mut ways = Vec::new();
    if let Some(&sign) = text.get(start) {
        if sign == b'+' || sign == b'-' {
            for end in digit_run_ends(text, start + 1) {
                ways.push(Some(SignedNumber {
                    negative: sign == b'-',
                    digits_start: start + 1,
                    digits_end: end,
                }));
            }
        }
    }
    ways.push(None);
    ways
}

#[derive(Clone, Copy, PartialEq)]
enum RangeKind {
    Deletion,
    Duplication,
    Insertion,
}

fn range_suffix_matches(text: &[u8], position: usize, kind: RangeKind) -> bool {
    let rest = &text[position..];
    match kind {
        RangeKind::Deletion => rest.starts_with(b"del"),
        RangeKind::Duplication => rest.starts_with(b"dup"),
        RangeKind::Insertion => {
            rest.starts_with(b"ins") && rest.get(3).is_some_and(|base| b"ATCG".contains(base))
        }
    }
}

/// `re.search(r"(\d+)?(?:(\+|-)(\d+))?_(\d+)?(?:(\+|-)(\d+))?del", text)` (or
/// `dup`, or `ins[ATCG]+`): the start and end offsets of the first match,
/// a missing offset group counting as 0. The tries follow the order of a
/// backtracking regex engine, so the first match is the one Python finds.
fn search_range(text: &[u8], kind: RangeKind) -> Outcome<Option<(i128, i128)>> {
    for start in 0..=text.len() {
        for first in optional_digits(text, start) {
            let after_first = first.unwrap_or(start);
            for start_offset in optional_signed(text, after_first) {
                let after_start = start_offset.map_or(after_first, |number| number.digits_end);
                if text.get(after_start) != Some(&b'_') {
                    continue;
                }
                let after_underscore = after_start + 1;
                for second in optional_digits(text, after_underscore) {
                    let after_second = second.unwrap_or(after_underscore);
                    for end_offset in optional_signed(text, after_second) {
                        let after_end = end_offset.map_or(after_second, |number| number.digits_end);
                        if range_suffix_matches(text, after_end, kind) {
                            let start_value = match start_offset {
                                Some(number) => signed_value(text, number)?,
                                None => 0,
                            };
                            let end_value = match end_offset {
                                Some(number) => signed_value(text, number)?,
                                None => 0,
                            };
                            return Ok(Some((start_value, end_value)));
                        }
                    }
                }
            }
        }
    }
    Ok(None)
}

#[derive(Clone, Copy, PartialEq)]
enum SingleKind {
    Indel,
    Substitution,
}

/// `re.search(r"(\d+)?(\+|-)(\d+)(del|dup)", text)` or
/// `re.search(r"(\d+)?(\+|-)(\d+)([ATCG]>[ATCG])", text)`: the signed offset
/// of the first match. Digits before the sign do not change which sign wins
/// (the earliest sign with digits and the right ending after it), so only the
/// signs are tried.
fn search_single(text: &[u8], kind: SingleKind) -> Outcome<Option<i128>> {
    for position in 0..text.len() {
        let sign = text[position];
        if sign != b'+' && sign != b'-' {
            continue;
        }
        let mut end = position + 1;
        while end < text.len() && text[end].is_ascii_digit() {
            end += 1;
        }
        if end == position + 1 {
            continue;
        }
        let rest = &text[end..];
        let suffix = match kind {
            SingleKind::Indel => rest.starts_with(b"del") || rest.starts_with(b"dup"),
            SingleKind::Substitution => {
                rest.len() >= 3
                    && b"ATCG".contains(&rest[0])
                    && rest[1] == b'>'
                    && b"ATCG".contains(&rest[2])
            }
        };
        if suffix {
            let number = SignedNumber {
                negative: sign == b'-',
                digits_start: position + 1,
                digits_end: end,
            };
            return Ok(Some(signed_value(text, number)?));
        }
    }
    Ok(None)
}

fn invalid_strand(strand: &str) -> String {
    format!("ValueError: Invalid strand: {strand}")
}

/// `parse_hgvsc_splice_position`: where a variant lies relative to the
/// canonical splice sites, from its HGVSc (for example
/// `ENST00000231790.2:c.116+1G>A`). `None` where the script returns `None`.
/// Every pattern is tried in the script's order and a later one overrides the
/// result of an earlier one, as the script does.
pub(crate) fn parse_hgvsc_splice_position(
    hgvsc: Option<&str>,
    strand: &str,
    var_pos: i64,
    var_length: i64,
) -> Outcome<Option<SplicePosition>> {
    let Some(hgvsc) = hgvsc else {
        return Ok(None);
    };
    if !hgvsc.contains(":c.") {
        return Ok(None);
    }
    if !hgvsc.contains('-') && !hgvsc.contains('+') {
        return Ok(None);
    }
    if hgvsc.contains('*') {
        return Ok(None);
    }
    let coding = hgvsc.split(":c.").nth(1).unwrap_or("");
    let text = coding.as_bytes();
    let var_pos = i128::from(var_pos);
    let mut variant_type = VariantType::Unknown;
    let mut distance: Option<i128> = None;
    let mut boundary: Option<i128> = None;

    // Range deletions, insertions and duplications. The three blocks differ
    // only in their formulas; `formulas` holds them in the order
    // donor/strand 1, donor/strand -1, acceptor/strand 1, acceptor/strand -1.
    type Formula = fn(i128, i128, i128) -> (i128, i128);
    let blocks: [(RangeKind, [Formula; 4]); 3] = [
        (
            RangeKind::Deletion,
            [
                |start, _end, pos| (start - 1 - 1, pos - start + 1),
                |_start, end, pos| (-(end + 1) + 2, pos + end + 1),
                |start, _end, pos| (start - 1 + 2, pos + 1 - start),
                |_start, end, pos| (end + 1 + 1, pos + 1 + end),
            ],
        ),
        (
            RangeKind::Insertion,
            [
                |start, _end, pos| (start - 1, pos - start),
                |_start, end, pos| (2 - end, pos + end),
                |start, _end, pos| (start + 2, pos - start),
                |_start, end, pos| (-1 - end, pos + end),
            ],
        ),
        (
            RangeKind::Duplication,
            [
                |start, _end, pos| (start - 1 - 1, pos - start + 1),
                |_start, end, pos| (-(end + 1) + 2, pos + end + 1),
                |start, _end, pos| (start - 1 + 2, pos - start + 1),
                |_start, end, pos| (-(end + 1) - 1, pos + end + 1),
            ],
        ),
    ];
    for (kind, formulas) in blocks {
        let Some((start, end)) = search_range(text, kind)? else {
            continue;
        };
        if start * end < 0 {
            // Either the HGVSc is odd or the change spans a whole exon.
            return Ok(Some(SplicePosition {
                variant_type: VariantType::Spanning,
                distance_to_canonical: distance,
                exon_boundary_pos: boundary,
                exon_skipping: true,
                overlapping_canonical_site: false,
            }));
        }
        let (donor, acceptor) = (start >= 0 && end >= 0, start <= 0 && end <= 0);
        if donor || acceptor {
            let base = if donor { 0 } else { 2 };
            variant_type = if donor {
                VariantType::DonorRegion
            } else {
                VariantType::AcceptorRegion
            };
            let formula = match strand {
                "1" => formulas[base],
                "-1" => formulas[base + 1],
                _ => return Err(invalid_strand(strand)),
            };
            let (new_distance, new_boundary) = formula(start, end, var_pos);
            distance = Some(new_distance);
            boundary = Some(new_boundary);
        }
    }

    // Single-base insertions and deletions with an intronic offset.
    if let Some(offset) = search_single(text, SingleKind::Indel)? {
        if offset > 0 {
            variant_type = VariantType::DonorRegion;
            match strand {
                "1" => {
                    distance = Some(offset - 1 - 1);
                    boundary = Some(var_pos - offset + 1);
                }
                "-1" => {
                    distance = Some(-(offset + 1) + 2);
                    boundary = Some(var_pos + 1 + offset);
                }
                _ => {}
            }
        } else if offset < 0 {
            variant_type = VariantType::AcceptorRegion;
            match strand {
                "1" => {
                    distance = Some(-(offset + 1));
                    boundary = Some(var_pos + 1 - offset);
                }
                "-1" => {
                    distance = Some(-(offset + 1) - 1);
                    boundary = Some(var_pos - (offset + 1));
                }
                _ => {}
            }
        } else {
            return Err("ValueError: Invalid offset: 0".to_owned());
        }
    }

    // Single-base substitutions with an intronic offset.
    if let Some(offset) = search_single(text, SingleKind::Substitution)? {
        if offset > 0 {
            variant_type = VariantType::DonorRegion;
            match strand {
                "1" => {
                    distance = Some(offset - 1);
                    boundary = Some(var_pos - offset);
                }
                "-1" => {
                    distance = Some(-offset + 2);
                    boundary = Some(var_pos + offset);
                }
                _ => {}
            }
        } else if offset < 0 {
            variant_type = VariantType::AcceptorRegion;
            match strand {
                "1" => {
                    distance = Some(offset + 2);
                    boundary = Some(var_pos - offset);
                }
                "-1" => {
                    distance = Some(-offset - 1);
                    boundary = Some(var_pos + offset);
                }
                _ => {}
            }
        } else {
            return Err("ValueError: Invalid offset: 0".to_owned());
        }
    }

    let Some(distance) = distance else {
        return Ok(None);
    };
    Ok(Some(SplicePosition {
        variant_type,
        distance_to_canonical: Some(distance),
        exon_boundary_pos: boundary,
        exon_skipping: false,
        overlapping_canonical_site: distance + i128::from(var_length) >= 0,
    }))
}

/// `[+-]\d+` found anywhere in the text: the signed numbers, left to right,
/// as `re.findall(r"([+-])(\d+)", text)` gives them (too long to matter: `None`).
fn signed_numbers(coding: &str) -> Vec<Option<i128>> {
    let bytes = coding.as_bytes();
    let mut found = Vec::new();
    let mut index = 0;
    while index < bytes.len() {
        if (bytes[index] == b'+' || bytes[index] == b'-')
            && bytes.get(index + 1).is_some_and(u8::is_ascii_digit)
        {
            let mut end = index + 1;
            while end < bytes.len() && bytes[end].is_ascii_digit() {
                end += 1;
            }
            let value = coding[index + 1..end].parse::<i128>().ok();
            found.push(value.map(|magnitude| {
                if bytes[index] == b'-' {
                    -magnitude
                } else {
                    magnitude
                }
            }));
            index = end;
        } else {
            index += 1;
        }
    }
    found
}

/// `_touches_canonical_splice_site`: the variant changes the canonical
/// +-1/+-2 splice dinucleotide. HGVSc offsets already follow the transcript,
/// so the check does not depend on the strand.
fn touches_canonical_splice_site(hgvsc: Option<&str>, consequence: &str) -> bool {
    let Some(hgvsc) = hgvsc else {
        return false;
    };
    let Some((_, coding)) = hgvsc.split_once(":c.") else {
        return false;
    };
    // 5' UTR (negative cDNA position) and 3' UTR (asterisk) are left to SpliceAI.
    if coding.starts_with('-') || coding.contains('*') {
        return false;
    }
    let is_donor = consequence.contains("splice_donor_variant");
    let hits_dinucleotide = signed_numbers(coding).iter().any(|offset| match offset {
        Some(offset) if is_donor => *offset == 1 || *offset == 2,
        Some(offset) => *offset == -1 || *offset == -2,
        None => false,
    });
    if hits_dinucleotide {
        return true;
    }
    // A range with one exonic and one intronic end always spans +-1/+-2.
    if coding.contains('_') {
        let sides: Vec<&str> = coding.split('_').collect();
        if sides.len() == 2 {
            let intronic: Vec<bool> = sides
                .iter()
                .map(|side| !signed_numbers(side).is_empty())
                .collect();
            if intronic.iter().any(|flag| *flag) && !intronic.iter().all(|flag| *flag) {
                return true;
            }
        }
    }
    false
}

// ------------------------------------------------------------- resources

struct Resources {
    /// Transcript, then exon number as text, then the domain paths.
    exon_domains: HashMap<String, HashMap<String, Vec<String>>>,
    intolerant_domains: HashSet<String>,
    interpro: InterproMapping,
    normalizer: DomainNormalizer,
}

struct ExonFlags {
    intolerant: bool,
    functional: bool,
}

impl Resources {
    /// The domain paths of one exon of a transcript in the map, if both are in it.
    fn domains_of(&self, transcript: Option<&str>, exon: i64) -> Option<&Vec<String>> {
        self.exon_domains.get(transcript?)?.get(&exon.to_string())
    }

    fn is_intolerant(&self, domains: &[String]) -> bool {
        domains
            .iter()
            .any(|domain| self.intolerant_domains.contains(domain))
    }

    /// Whether the exon's domains include an intolerant one or a functional one.
    fn exon_flags(&self, transcript: Option<&str>, exon: i64) -> Outcome<ExonFlags> {
        let Some(domains) = self.domains_of(transcript, exon) else {
            return Ok(ExonFlags {
                intolerant: false,
                functional: false,
            });
        };
        let mut functional = false;
        for path in domains {
            let annotation = path.split_once(':').map_or(path.as_str(), |(_, rest)| rest);
            if self
                .normalizer
                .interpret_functionality(annotation, &self.interpro)?
                == Functionality::Functional
            {
                functional = true;
            }
        }
        Ok(ExonFlags {
            intolerant: self.is_intolerant(domains),
            functional,
        })
    }
}

fn domain_list(value: &Value) -> Result<Vec<String>> {
    let items: Vec<Value> = match value {
        Value::List(items) => items.borrow().clone(),
        Value::Tuple(items) => (**items).clone(),
        Value::Set(set) => set.borrow().items.clone(),
        Value::FrozenSet(set) => set.items.clone(),
        other => {
            return Err(invalid(format!(
                "a domain list in the pickle is a {}",
                other.type_name()
            )));
        }
    };
    items
        .iter()
        .map(|item| {
            item.as_str()
                .map(str::to_owned)
                .ok_or_else(|| invalid("a domain path in the pickle is not text"))
        })
        .collect()
}

/// `{transcript: {"exon_domains": {exon: [domain path, ...]}}}` from the pickle.
fn exon_domain_map(value: &Value) -> Result<HashMap<String, HashMap<String, Vec<String>>>> {
    let Value::Dict(outer) = value else {
        return Err(invalid(format!(
            "the transcript domain map is a {}, not a dict",
            value.type_name()
        )));
    };
    let mut map = HashMap::new();
    for (key, entry) in &outer.borrow().entries {
        // A key that is not text can never equal a transcript ID.
        let Some(transcript) = key.as_str() else {
            continue;
        };
        let Value::Dict(inner) = entry else {
            return Err(invalid(format!(
                "the transcript domain map entry for {transcript} is a {}, not a dict",
                entry.type_name()
            )));
        };
        let mut by_exon: HashMap<String, Vec<String>> = HashMap::new();
        for (name, exons) in &inner.borrow().entries {
            if name.as_str() != Some("exon_domains") {
                continue;
            }
            let Value::Dict(exons) = exons else {
                return Err(invalid(format!(
                    "exon_domains of {transcript} is a {}, not a dict",
                    exons.type_name()
                )));
            };
            for (exon, domains) in &exons.borrow().entries {
                if let Some(exon) = exon.as_str() {
                    by_exon.insert(exon.to_owned(), domain_list(domains)?);
                }
            }
        }
        map.insert(transcript.to_owned(), by_exon);
    }
    Ok(map)
}

fn load_pickle(path: &str) -> Result<Value> {
    // These two files are decompressed by their name, as the script does.
    let bytes = pickle_reader::read_file(Path::new(path), path.ends_with(".gz"))?;
    pickle_reader::load(&bytes)
}

// ----------------------------------------------------------- row analysis

#[derive(Clone, Copy)]
enum Offset {
    Zero,
    Value(f64),
}

struct RowResult {
    lof: bool,
    len_changing: bool,
    frameshift: bool,
    span_intol_domain: bool,
    ten_percent_protein: bool,
    gain_offset: Offset,
    loss_offset: Offset,
    affected_exons: String,
    last_exon_only: bool,
}

struct ColumnIndexes {
    consequence: Option<usize>,
    intron: Option<usize>,
    exon: Option<usize>,
    feature: Option<usize>,
    events: Option<usize>,
    delta: Option<usize>,
    cds_position: Option<usize>,
    hgvsc: Option<usize>,
    pos: Option<usize>,
    strand: Option<usize>,
    dp_ag: Option<usize>,
    dp_al: Option<usize>,
    dp_dg: Option<usize>,
    dp_dl: Option<usize>,
    ds_ag: Option<usize>,
    ds_al: Option<usize>,
    ds_dg: Option<usize>,
    ds_dl: Option<usize>,
}

impl ColumnIndexes {
    fn find(table: &Table) -> ColumnIndexes {
        let at = |name: &str| table.column_index(name);
        ColumnIndexes {
            consequence: at("Consequence"),
            intron: at("INTRON"),
            exon: at("EXON"),
            feature: at("Feature"),
            events: at("SpliceVault_top_events"),
            delta: at("SpliceVault_SpliceAI_delta"),
            cds_position: at("CDS_position"),
            hgvsc: at("HGVSc"),
            pos: at("pos"),
            strand: at("STRAND"),
            dp_ag: at("SpliceAI_pred_DP_AG"),
            dp_al: at("SpliceAI_pred_DP_AL"),
            dp_dg: at("SpliceAI_pred_DP_DG"),
            dp_dl: at("SpliceAI_pred_DP_DL"),
            ds_ag: at("SpliceAI_pred_DS_AG"),
            ds_al: at("SpliceAI_pred_DS_AL"),
            ds_dg: at("SpliceAI_pred_DS_DG"),
            ds_dl: at("SpliceAI_pred_DS_DL"),
        }
    }
}

struct Context<'a> {
    table: &'a Table,
    columns: ColumnIndexes,
    resources: &'a Resources,
    splice_altering_cutoff: f64,
    acceptor_rescue_cutoff: f64,
}

impl<'a> Context<'a> {
    /// `row[name]`: Python raises KeyError for a column the table lacks.
    fn required(&self, name: &str, index: Option<usize>, row: usize) -> Outcome<Cell<'a>> {
        match index {
            Some(index) => Ok(self.table.cell(index, row)),
            None => Err(format!("KeyError: {name:?}")),
        }
    }
}

/// What the SpliceAI part of the interpretation finds.
struct SpliceAiOutcome {
    is_harmful: bool,
    length_changing: bool,
    frame_shift: bool,
    span_intol_domain: bool,
    affected_exons: BTreeSet<i64>,
    ten_percent_protein: bool,
    gain_offset: Offset,
    loss_offset: Offset,
}

struct SpliceAiInput<'a> {
    ds_ag: f64,
    ds_al: f64,
    ds_dg: f64,
    ds_dl: f64,
    dp_ag: f64,
    dp_al: f64,
    dp_dg: f64,
    dp_dl: f64,
    transcript: Option<&'a str>,
    /// `None` when the Consequence value is not text, which makes `"UTR" in
    /// consequence` raise TypeError.
    consequence: Option<&'a str>,
    cds_pos: Position<'a>,
    intron_pos: Position<'a>,
    exon_pos: Position<'a>,
    hgvsc: Option<String>,
    var_pos: i64,
    strand: String,
}

fn contains_utr(consequence: Option<&str>) -> Outcome<bool> {
    consequence
        .map(|text| text.contains("UTR"))
        .ok_or_else(|| "TypeError: Consequence is not text".to_owned())
}

fn checked_next(number: i64) -> Outcome<i64> {
    number
        .checked_add(1)
        .ok_or_else(|| "an exon number overflows".to_owned())
}

/// `SpliceAI_interpretation`.
fn spliceai_interpretation(context: &Context, input: &SpliceAiInput) -> Outcome<SpliceAiOutcome> {
    let cutoff = context.splice_altering_cutoff;
    let mut outcome = SpliceAiOutcome {
        is_harmful: false,
        length_changing: false,
        frame_shift: false,
        span_intol_domain: false,
        affected_exons: BTreeSet::new(),
        ten_percent_protein: false,
        gain_offset: Offset::Zero,
        loss_offset: Offset::Zero,
    };
    let parse =
        || parse_hgvsc_splice_position(input.hgvsc.as_deref(), &input.strand, input.var_pos, 1);
    // The distance when the variant lies in a region of the given type.
    let distance_in = |parsed: &Option<SplicePosition>, wanted: VariantType| -> Option<f64> {
        let parsed = parsed.as_ref()?;
        if parsed.variant_type == wanted {
            parsed.distance_to_canonical.map(|distance| distance as f64)
        } else {
            None
        }
    };
    let frame_and_size = |outcome: &mut SpliceAiOutcome, offset: f64| -> Outcome<()> {
        if offset.abs() < 3.0 {
            outcome.length_changing = false;
        } else if offset.abs() % 3.0 != 0.0 {
            outcome.frame_shift = true;
            outcome.is_harmful = true;
        }
        Ok(())
    };

    // Acceptor gain
    if input.ds_ag > cutoff && !contains_utr(input.consequence)? {
        outcome.length_changing = true;
        let parsed = parse()?;
        let mut offset = 0.0;
        if let Some(distance) = distance_in(&parsed, VariantType::AcceptorRegion) {
            offset = input.dp_ag + distance;
            outcome.gain_offset = Offset::Value(offset);
        }
        frame_and_size(&mut outcome, offset)?;
        if input.cds_pos.truthy() && exceeds_a_tenth(offset, cds_total(&input.cds_pos)?)? {
            outcome.ten_percent_protein = true;
            outcome.is_harmful = true;
        }
    }

    // Donor gain
    if input.ds_dg > cutoff {
        outcome.length_changing = true;
        if input.exon_pos.truthy() {
            outcome
                .affected_exons
                .insert(first_number(input.exon_pos.text()?)?);
        }
        let parsed = parse()?;
        let mut offset = 0.0;
        if let Some(distance) = distance_in(&parsed, VariantType::DonorRegion) {
            offset = input.dp_dg + distance;
            outcome.gain_offset = Offset::Value(offset);
        }
        frame_and_size(&mut outcome, offset)?;
        if input.cds_pos.truthy() && exceeds_a_tenth(offset, cds_total(&input.cds_pos)?)? {
            outcome.ten_percent_protein = true;
            outcome.is_harmful = true;
        }
    }

    // Acceptor loss
    if input.ds_al > cutoff && !contains_utr(input.consequence)? {
        let parsed = parse()?;
        let mut offset = 0.0;
        if let Some(distance) = distance_in(&parsed, VariantType::AcceptorRegion) {
            offset = input.dp_al + distance;
            outcome.loss_offset = Offset::Value(offset);
        }
        if offset.abs() <= 1.0 {
            outcome.length_changing = true;
            if input.ds_ag > context.acceptor_rescue_cutoff {
                // A nearby acceptor site may rescue the lost one.
                offset = 0.0;
                if let Some(distance) = distance_in(&parsed, VariantType::AcceptorRegion) {
                    offset = input.dp_ag + distance;
                    outcome.gain_offset = Offset::Value(offset);
                }
                frame_and_size(&mut outcome, offset)?;
                if offset.abs() >= 30.0 {
                    if input.intron_pos.truthy() {
                        let intron = first_number(input.intron_pos.text()?)?;
                        outcome.affected_exons.insert(checked_next(intron)?);
                    }
                    if input.exon_pos.truthy() {
                        outcome
                            .affected_exons
                            .insert(first_number(input.exon_pos.text()?)?);
                    }
                }
                if input.cds_pos.truthy() && exceeds_a_tenth(offset, cds_total(&input.cds_pos)?)? {
                    outcome.ten_percent_protein = true;
                    outcome.is_harmful = true;
                }
            } else {
                if input.exon_pos.truthy() {
                    outcome
                        .affected_exons
                        .insert(first_number(input.exon_pos.text()?)?);
                }
                if input.intron_pos.truthy() {
                    let intron = first_number(input.intron_pos.text()?)?;
                    outcome.affected_exons.insert(checked_next(intron)?);
                }
            }
        }
    }

    // Donor loss
    if input.ds_dl > cutoff {
        let parsed = parse()?;
        let mut offset = 2.0;
        if let Some(distance) = distance_in(&parsed, VariantType::DonorRegion) {
            offset = input.dp_dl + distance;
            outcome.loss_offset = Offset::Value(offset);
        }
        if offset.abs() <= 1.0 {
            outcome.length_changing = true;
            outcome.is_harmful = true;
            outcome.frame_shift = true;
            outcome.ten_percent_protein = true;
        }
    }

    // Affected exons that hold intolerant or functional domains
    for exon in outcome.affected_exons.clone() {
        let flags = context.resources.exon_flags(input.transcript, exon)?;
        if flags.intolerant || flags.functional {
            outcome.span_intol_domain = true;
            outcome.is_harmful = true;
        }
    }
    Ok(outcome)
}

// SpliceVault events ------------------------------------------------------

#[derive(Default)]
struct EventEffect {
    affected_exons: BTreeSet<i64>,
    large_truncation: bool,
    length_changing: bool,
    is_lof: bool,
    frame_shift: bool,
    loc_intol_domain: bool,
}

struct EventAnalysis {
    affected_exons: BTreeSet<i64>,
    is_harmful: bool,
    length_changing: bool,
    frame_shift: bool,
    span_intol_domain: bool,
    ten_percent_protein: bool,
    fraction: f64,
}

/// Marks the effect when any affected exon holds an intolerant or functional domain.
fn mark_domain_overlaps(
    resources: &Resources,
    transcript: Option<&str>,
    effect: &mut EventEffect,
) -> Outcome<()> {
    for exon in effect.affected_exons.clone() {
        let flags = resources.exon_flags(transcript, exon)?;
        if flags.intolerant || flags.functional {
            effect.is_lof = true;
            effect.loc_intol_domain = true;
        }
    }
    Ok(())
}

fn process_exon_skipping(
    resources: &Resources,
    pos_str: &str,
    transcript: Option<&str>,
) -> Outcome<EventEffect> {
    let mut effect = EventEffect {
        length_changing: true,
        ..EventEffect::default()
    };
    if pos_str.contains('-') {
        let parts: Vec<&str> = pos_str.split('-').collect();
        if parts.len() != 2 {
            return Err(format!(
                "ValueError: {pos_str:?} is not a range of two exons"
            ));
        }
        let (start, end) = (text_int(parts[0])?, text_int(parts[1])?);
        if end.saturating_sub(start) > 10_000_000 {
            return Err("an exon-skipping range is too long".to_owned());
        }
        for exon in start..=end {
            effect.affected_exons.insert(exon);
        }
    } else {
        effect.affected_exons.insert(text_int(pos_str)?);
    }
    mark_domain_overlaps(resources, transcript, &mut effect)?;
    Ok(effect)
}

fn process_cryptic_donor(
    resources: &Resources,
    pos_str: &str,
    transcript: Option<&str>,
    positions: (&Position, &Position, &Position),
) -> Outcome<EventEffect> {
    let (intron_pos, exon_pos, cds_pos) = positions;
    let intron = pair_of(intron_pos)?;
    let exon = pair_of(exon_pos)?;
    let total_cds = cds_total_of_event(cds_pos, false)?;
    let offset = text_int(pos_str.trim_start_matches(['+', '-']))?;
    let magnitude = offset.unsigned_abs();
    let mut effect = EventEffect::default();
    if pos_str.starts_with('+') {
        // A downstream cryptic donor lies in the intron.
        if magnitude > 3 {
            effect.length_changing = true;
        }
        if magnitude % 3 != 0 {
            effect.frame_shift = true;
        }
    } else {
        // An upstream cryptic donor lies in the exon; very close ones are ignored.
        if magnitude > 3 {
            effect.length_changing = true;
            if exon_pos.truthy() {
                if let Some((exon_number, _)) = exon {
                    effect.affected_exons.insert(exon_number);
                }
            } else if intron_pos.truthy() {
                if let Some((intron_number, _)) = intron {
                    effect.affected_exons.insert(intron_number);
                }
            }
        }
        if let Some(total) = total_cds.filter(|total| *total != 0) {
            if (magnitude as f64) / (total as f64) > 0.1 {
                effect.large_truncation = true;
                effect.is_lof = true;
            }
        }
        if magnitude % 3 != 0 {
            effect.frame_shift = true;
        }
    }
    mark_domain_overlaps(resources, transcript, &mut effect)?;
    Ok(effect)
}

fn process_cryptic_acceptor(
    resources: &Resources,
    pos_str: &str,
    transcript: Option<&str>,
    positions: (&Position, &Position, &Position),
) -> Outcome<EventEffect> {
    let (intron_pos, exon_pos, cds_pos) = positions;
    let intron = pair_of(intron_pos)?;
    pair_of(exon_pos)?;
    let total_cds = cds_total_of_event(cds_pos, true)?;
    let offset = text_int(pos_str.trim_start_matches(['+', '-']))?;
    let magnitude = offset.unsigned_abs();
    let mut effect = EventEffect {
        length_changing: true,
        ..EventEffect::default()
    };
    if pos_str.starts_with('-') {
        // An upstream cryptic acceptor lies in the intron.
        if let Some((intron_number, _)) = intron {
            effect.affected_exons.insert(checked_next(intron_number)?);
        }
        if magnitude % 3 != 0 {
            effect.frame_shift = true;
        }
    } else {
        if magnitude % 3 != 0 {
            effect.frame_shift = true;
        }
        if let Some(total) = total_cds.filter(|total| *total != 0) {
            if (magnitude as f64) / (total as f64) > 0.1 {
                effect.large_truncation = true;
                effect.is_lof = true;
            }
        }
    }
    mark_domain_overlaps(resources, transcript, &mut effect)?;
    Ok(effect)
}

/// `analyze_splice_event`: one event such as `Top1:CD:-10:2%:Frameshift`.
fn analyze_splice_event(
    resources: &Resources,
    event: &str,
    transcript: Option<&str>,
    positions: (&Position, &Position, &Position),
) -> Outcome<EventAnalysis> {
    let parts: Vec<&str> = event.split(':').collect();
    if parts.len() != 5 {
        return Err(format!("ValueError: Error splitting event string {event}"));
    }
    let (event_type, pos_str, frequency, frame_impact) = (parts[1], parts[2], parts[3], parts[4]);
    let percent = |text: &str| {
        parse_float(text).ok_or_else(|| format!("ValueError: could not convert {text:?} to float"))
    };
    let fraction = if frequency.contains('%') {
        percent(frequency.trim_end_matches('%'))? / 100.0
    } else if frequency.contains('#') {
        percent(frequency.trim_start_matches('#'))? / 100.0
    } else {
        return Err(format!(
            "ValueError: The frequency format is not recognized: {frequency}, the event string is {event}"
        ));
    };
    let mut ten_percent_protein = false;
    let effect = match event_type {
        "ES" => {
            ten_percent_protein = true;
            process_exon_skipping(resources, pos_str, transcript)?
        }
        "CD" => process_cryptic_donor(resources, pos_str, transcript, positions)?,
        "CA" => process_cryptic_acceptor(resources, pos_str, transcript, positions)?,
        _ => EventEffect::default(),
    };
    // An event is harmful when it touches an intolerant domain or is a frameshift.
    let mut intolerant_affected = false;
    for exon in &effect.affected_exons {
        if resources
            .domains_of(transcript, *exon)
            .is_some_and(|domains| resources.is_intolerant(domains))
        {
            intolerant_affected = true;
        }
    }
    let mut is_harmful = effect.is_lof;
    let mut frame_shift = effect.frame_shift;
    if intolerant_affected {
        is_harmful = true;
    }
    if frame_impact == "Frameshift" {
        is_harmful = true;
        frame_shift = true;
    }
    if effect.large_truncation {
        ten_percent_protein = true;
    }
    Ok(EventAnalysis {
        affected_exons: effect.affected_exons,
        is_harmful,
        length_changing: effect.length_changing,
        frame_shift,
        span_intol_domain: effect.loc_intol_domain || intolerant_affected,
        ten_percent_protein,
        fraction,
    })
}

#[derive(Default)]
struct SpliceVaultOutcome {
    harmful: bool,
    length_changing: bool,
    frame_shift: bool,
    span_intol_domain: bool,
    affected_exons: Vec<i64>,
    ten_percent_protein: bool,
}

fn fraction_sum(analyses: &[EventAnalysis], pick: fn(&EventAnalysis) -> bool) -> f64 {
    let fractions: Vec<f64> = analyses
        .iter()
        .filter(|analysis| pick(analysis))
        .map(|analysis| analysis.fraction)
        .collect();
    sum_floats(&fractions)
}

fn odds_ratio(numerator: f64, denominator: f64) -> f64 {
    if denominator > 0.0 {
        numerator / denominator
    } else {
        20.0
    }
}

/// `SpliceVault_interpretation`: the top-4 SpliceVault events of a variant,
/// for example `Top1:CD:-10:2%:Frameshift&Top2:CD:-4:0.6%:Frameshift`.
fn splicevault_interpretation(
    resources: &Resources,
    events: &Cell,
    spliceai_delta: f64,
    transcript: Option<&str>,
    positions: (&Position, &Position, &Position),
) -> Outcome<SpliceVaultOutcome> {
    let Cell::Text(value) = events else {
        return Ok(SpliceVaultOutcome::default());
    };
    let analyses: Vec<EventAnalysis> = value
        .split('&')
        .map(|event| analyze_splice_event(resources, event, transcript, positions))
        .collect::<Outcome<Vec<_>>>()?;
    let mut exon_fraction: BTreeMap<i64, f64> = BTreeMap::new();
    for analysis in &analyses {
        for exon in &analysis.affected_exons {
            *exon_fraction.entry(*exon).or_insert(0.0) += analysis.fraction;
        }
    }
    let confident = spliceai_delta >= 0.5;
    let mut outcome = SpliceVaultOutcome {
        affected_exons: exon_fraction
            .iter()
            .filter(|(_, fraction)| **fraction > 0.3 && confident)
            .map(|(exon, _)| *exon)
            .collect(),
        ..SpliceVaultOutcome::default()
    };
    // A flag is set when its events carry over 30% of the samples and outweigh
    // the events without it twenty to one.
    let decides =
        |with: f64, without: f64| odds_ratio(with, without) >= 20.0 && with > 0.3 && confident;
    let harmful = fraction_sum(&analyses, |a| a.is_harmful);
    let not_harmful = fraction_sum(&analyses, |a| !a.is_harmful);
    let changing = fraction_sum(&analyses, |a| a.length_changing);
    let not_changing = fraction_sum(&analyses, |a| !a.length_changing);
    let frameshift = fraction_sum(&analyses, |a| a.frame_shift);
    let not_frameshift = fraction_sum(&analyses, |a| !a.frame_shift);
    let domain = fraction_sum(&analyses, |a| a.span_intol_domain);
    let not_domain = fraction_sum(&analyses, |a| !a.span_intol_domain);
    let tenth = fraction_sum(&analyses, |a| a.ten_percent_protein);
    let not_tenth = fraction_sum(&analyses, |a| !a.ten_percent_protein);
    if decides(harmful, not_harmful) {
        outcome.harmful = true;
    }
    if decides(changing, not_changing) {
        outcome.length_changing = true;
    }
    if decides(domain, not_domain) {
        outcome.span_intol_domain = true;
    }
    if decides(frameshift, not_frameshift) {
        outcome.frame_shift = true;
        outcome.harmful = true;
    }
    if decides(tenth, not_tenth) {
        outcome.ten_percent_protein = true;
        outcome.harmful = true;
    }
    Ok(outcome)
}

/// `current_intron == total_exons - 1` for `intron_pos`, whose first number may
/// be a range such as `1-3`; `None` (False) when it does not parse.
fn is_last_intron(intron_pos: &Position, total_exons: i64) -> Outcome<bool> {
    let text = intron_pos.text()?;
    let first = text
        .split('/')
        .next()
        .unwrap_or("")
        .split('-')
        .next()
        .unwrap_or("");
    Ok(parse_int(first).is_some_and(|current| current == total_exons - 1))
}

/// `splicing_altering_per_row` for one row of the table.
fn analyze_row(context: &Context, row: usize) -> Outcome<RowResult> {
    let columns = &context.columns;
    let consequence_cell = columns
        .consequence
        .map(|index| context.table.cell(index, row));
    let consequence = consequence_cell
        .as_ref()
        .map_or_else(String::new, cell_text);

    let intron_pos = Position::from_cell(context.required("INTRON", columns.intron, row)?);
    let exon_pos = Position::from_cell(context.required("EXON", columns.exon, row)?);

    // The number of exons, from EXON ("3/10") or else from INTRON ("3/9" has one more exon).
    let mut total_exons: Option<i64> = None;
    if exon_pos.truthy() {
        total_exons = exon_pos.text()?.split('/').nth(1).and_then(parse_int);
    }
    if total_exons.is_none() && intron_pos.truthy() {
        total_exons = intron_pos
            .text()?
            .split('/')
            .nth(1)
            .and_then(parse_int)
            .and_then(|total| total.checked_add(1));
    }
    let has_total = total_exons.is_some_and(|total| total != 0);

    // Canonical splice site variants are null variants whatever the scores say.
    let hgvsc_for_check = match columns.hgvsc.map(|index| context.table.cell(index, row)) {
        Some(Cell::Text(text)) => Some(text),
        _ => None,
    };
    if (consequence.contains("splice_donor_variant")
        || consequence.contains("splice_acceptor_variant"))
        && touches_canonical_splice_site(hgvsc_for_check, &consequence)
    {
        let mut last_exon_only = false;
        if intron_pos.truthy() && has_total {
            last_exon_only = is_last_intron(&intron_pos, total_exons.unwrap_or(0))?;
        }
        return Ok(RowResult {
            lof: true,
            len_changing: true,
            frameshift: true,
            span_intol_domain: false,
            ten_percent_protein: true,
            gain_offset: Offset::Zero,
            loss_offset: Offset::Zero,
            affected_exons: String::new(),
            last_exon_only,
        });
    }

    let feature = context.required("Feature", columns.feature, row)?;
    let events = context.required("SpliceVault_top_events", columns.events, row)?;
    let spliceai_delta =
        cell_float(&context.required("SpliceVault_SpliceAI_delta", columns.delta, row)?)?;
    let cds_pos =
        Position::from_cell(context.required("CDS_position", columns.cds_position, row)?);
    let transcript = match feature {
        Cell::Text(text) => Some(text),
        _ => None,
    };
    let positions = (&intron_pos, &exon_pos, &cds_pos);

    let splicevault = splicevault_interpretation(
        context.resources,
        &events,
        spliceai_delta,
        transcript,
        positions,
    )?;

    let float_at = |name: &str, index: Option<usize>| -> Outcome<f64> {
        cell_float(&context.required(name, index, row)?)
    };
    let dp_ag = float_at("SpliceAI_pred_DP_AG", columns.dp_ag)?;
    let dp_al = float_at("SpliceAI_pred_DP_AL", columns.dp_al)?;
    let dp_dg = float_at("SpliceAI_pred_DP_DG", columns.dp_dg)?;
    let dp_dl = float_at("SpliceAI_pred_DP_DL", columns.dp_dl)?;
    let ds_ag = float_at("SpliceAI_pred_DS_AG", columns.ds_ag)?;
    let ds_al = float_at("SpliceAI_pred_DS_AL", columns.ds_al)?;
    let ds_dg = float_at("SpliceAI_pred_DS_DG", columns.ds_dg)?;
    let ds_dl = float_at("SpliceAI_pred_DS_DL", columns.ds_dl)?;
    let var_pos = cell_int(&context.required("pos", columns.pos, row)?)?;
    let strand = cell_int(&context.required("STRAND", columns.strand, row)?)?.to_string();
    let hgvsc_cell = context.required("HGVSc", columns.hgvsc, row)?;
    let hgvsc = if na_value(&hgvsc_cell) {
        None
    } else {
        Some(cell_text(&hgvsc_cell))
    };
    let consequence_value = match context.required("Consequence", columns.consequence, row)? {
        Cell::Text(text) => Some(text),
        _ => None,
    };
    let spliceai = spliceai_interpretation(
        context,
        &SpliceAiInput {
            ds_ag,
            ds_al,
            ds_dg,
            ds_dl,
            dp_ag,
            dp_al,
            dp_dg,
            dp_dl,
            transcript,
            consequence: consequence_value,
            cds_pos,
            intron_pos,
            exon_pos,
            hgvsc,
            var_pos,
            strand,
        },
    )?;

    // Affected exons from both SpliceVault and SpliceAI.
    let mut all_affected: BTreeSet<i64> = splicevault.affected_exons.iter().copied().collect();
    all_affected.extend(spliceai.affected_exons.iter().copied());
    let affected_exons = all_affected
        .iter()
        .map(i64::to_string)
        .collect::<Vec<_>>()
        .join(",");

    // NMD escape: the effect reaches only the last exon.
    let mut last_exon_only = false;
    if has_total && !all_affected.is_empty() {
        let numbers: BTreeSet<i64> = all_affected
            .iter()
            .copied()
            .filter(|exon| *exon >= 0)
            .collect();
        if !numbers.is_empty() && numbers == BTreeSet::from([total_exons.unwrap_or(0)]) {
            last_exon_only = true;
        }
    }
    if !last_exon_only && intron_pos.truthy() && has_total {
        last_exon_only = is_last_intron(&intron_pos, total_exons.unwrap_or(0))?;
    }
    if !last_exon_only && exon_pos.truthy() && has_total {
        let text = exon_pos.text()?;
        let first = text
            .split('/')
            .next()
            .unwrap_or("")
            .split('-')
            .next()
            .unwrap_or("");
        last_exon_only = parse_int(first).is_some_and(|current| Some(current) == total_exons);
    }

    Ok(RowResult {
        lof: splicevault.harmful || spliceai.is_harmful,
        len_changing: splicevault.length_changing || spliceai.length_changing,
        frameshift: splicevault.frame_shift || spliceai.frame_shift,
        span_intol_domain: splicevault.span_intol_domain || spliceai.span_intol_domain,
        ten_percent_protein: splicevault.ten_percent_protein || spliceai.ten_percent_protein,
        gain_offset: spliceai.gain_offset,
        loss_offset: spliceai.loss_offset,
        affected_exons,
        last_exon_only,
    })
}

// ------------------------------------------------------------------ driver

struct Options {
    anno_table: String,
    transcript_domain_map: String,
    intolerant_domains: String,
    interpro_entry_map_pkl: String,
    threads: i64,
    config: Option<String>,
}

enum Parsed {
    Help,
    Run(Options),
}

/// argparse's `^-\d+$|^-\d*\.\d+$`.
fn is_negative_number(text: &str) -> bool {
    let Some(rest) = text.strip_prefix('-') else {
        return false;
    };
    let digits = |part: &str| part.bytes().all(|byte| byte.is_ascii_digit());
    match rest.split_once('.') {
        Some((whole, fraction)) => digits(whole) && !fraction.is_empty() && digits(fraction),
        None => !rest.is_empty() && digits(rest),
    }
}

/// argparse takes an argument as an option's value unless it looks like
/// another option; a negative number or a lone `-` is a value.
fn looks_like_value(text: &str) -> bool {
    !text.starts_with('-') || text == "-" || is_negative_number(text) || text.contains(' ')
}

/// argparse accepts any unambiguous prefix of a long option.
fn resolve_option(name: &str) -> Result<&'static str> {
    if let Some(exact) = LONG_OPTIONS.iter().find(|option| **option == name) {
        return Ok(*exact);
    }
    let candidates: Vec<&'static str> = LONG_OPTIONS
        .iter()
        .copied()
        .filter(|option| name.len() > 2 && option.starts_with(name))
        .collect();
    match candidates.len() {
        1 => Ok(candidates[0]),
        0 => Err(invalid(format!("unrecognized argument {name:?}\n{USAGE}"))),
        _ => Err(invalid(format!("ambiguous option {name:?}\n{USAGE}"))),
    }
}

fn parse_options(args: &[String]) -> Result<Parsed> {
    let mut values: HashMap<&'static str, String> = HashMap::new();
    let mut index = 0;
    while index < args.len() {
        let argument = &args[index];
        if argument == "-h" {
            return Ok(Parsed::Help);
        }
        if !argument.starts_with("--") || argument == "--" {
            return Err(invalid(format!(
                "unrecognized arguments: {argument}\n{USAGE}"
            )));
        }
        let (name, inline) = match argument.split_once('=') {
            Some((name, value)) => (name, Some(value.to_owned())),
            None => (argument.as_str(), None),
        };
        let option = resolve_option(name)?;
        if option == "--help" {
            return Ok(Parsed::Help);
        }
        let value = match inline {
            Some(value) => value,
            None => {
                index += 1;
                match args.get(index) {
                    Some(next) if looks_like_value(next) => next.clone(),
                    _ => {
                        return Err(invalid(format!(
                            "argument {option}: expected one argument\n{USAGE}"
                        )));
                    }
                }
            }
        };
        values.insert(option, value);
        index += 1;
    }
    let missing: Vec<&str> = [
        "--anno_table",
        "--transcript_domain_map",
        "--intolerant_domains",
        "--interpro_entry_map_pkl",
    ]
    .into_iter()
    .filter(|option| !values.contains_key(*option))
    .collect();
    if !missing.is_empty() {
        return Err(invalid(format!(
            "the following arguments are required: {}\n{USAGE}",
            missing.join(", ")
        )));
    }
    let threads = match values.remove("--threads") {
        Some(text) => parse_int(&text).ok_or_else(|| {
            invalid(format!(
                "argument --threads: invalid int value: {text:?}\n{USAGE}"
            ))
        })?,
        None => 12,
    };
    let mut take = |option: &str| values.remove(option).unwrap_or_default();
    Ok(Parsed::Run(Options {
        anno_table: take("--anno_table"),
        transcript_domain_map: take("--transcript_domain_map"),
        intolerant_domains: take("--intolerant_domains"),
        interpro_entry_map_pkl: take("--interpro_entry_map_pkl"),
        threads,
        config: values.remove("--config"),
    }))
}

fn judge_rows(context: &Context, threads: usize) -> Outcome<Vec<RowResult>> {
    let rows = context.table.rows;
    let workers = threads.min(rows).max(1);
    let chunk = rows.div_ceil(workers);
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|worker| {
                let first = worker * chunk;
                let last = ((worker + 1) * chunk).min(rows);
                scope.spawn(move || -> Outcome<Vec<RowResult>> {
                    (first..last).map(|row| analyze_row(context, row)).collect()
                })
            })
            .collect();
        let mut results = Vec::with_capacity(rows);
        for handle in handles {
            match handle.join() {
                Ok(part) => results.extend(part?),
                Err(_) => return Err("a worker thread panicked".to_owned()),
            }
        }
        Ok(results)
    })
}

fn bool_column(results: &[RowResult], pick: fn(&RowResult) -> bool) -> Column {
    Column::Bool(results.iter().map(pick).collect())
}

/// A tuple of ints and floats becomes an int column only when every item is an int.
fn offset_column(results: &[RowResult], pick: fn(&RowResult) -> Offset) -> Column {
    let offsets: Vec<Offset> = results.iter().map(pick).collect();
    if offsets
        .iter()
        .any(|offset| matches!(offset, Offset::Value(_)))
    {
        Column::Float(
            offsets
                .iter()
                .map(|offset| match offset {
                    Offset::Zero => 0.0,
                    Offset::Value(value) => *value,
                })
                .collect(),
        )
    } else {
        Column::Int(vec![0; offsets.len()])
    }
}

/// A name that is different on every run, standing in for the script's uuid4.
fn unique_tag() -> String {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |elapsed| elapsed.as_nanos());
    format!("{:x}{:x}", std::process::id(), nanos)
}

/// `splicing_interpretation`.
fn splicing_interpretation(options: &Options, cutoffs: Cutoffs) -> Result<()> {
    let mut table = read_table(Path::new(&options.anno_table))?;
    let intolerant_domains: HashSet<String> =
        domain_list(&load_pickle(&options.intolerant_domains)?)?
            .into_iter()
            .collect();
    let exon_domains = exon_domain_map(&load_pickle(&options.transcript_domain_map)?)?;
    let interpro = load_interpro_mapping(Path::new(&options.interpro_entry_map_pkl))?;
    let resources = Resources {
        exon_domains,
        intolerant_domains,
        interpro,
        normalizer: DomainNormalizer::new()?,
    };
    if options.threads < 1 {
        return Err(invalid(
            "ValueError: Number of processes must be at least 1",
        ));
    }
    eprintln!(
        "Processing {} variants using {} cores, spliceai_cutoff is {}, acceptor_rescue_cutoff is {}",
        table.rows,
        options.threads,
        float_repr(cutoffs.splice_altering),
        float_repr(cutoffs.nearby_acceptor_rescue)
    );
    if table.rows == 0 {
        // zip(*[]) cannot be unpacked into nine names.
        return Err(invalid(
            "ValueError: not enough values to unpack (expected 9, got 0)",
        ));
    }
    let results = {
        let context = Context {
            table: &table,
            columns: ColumnIndexes::find(&table),
            resources: &resources,
            splice_altering_cutoff: cutoffs.splice_altering,
            acceptor_rescue_cutoff: cutoffs.nearby_acceptor_rescue,
        };
        judge_rows(
            &context,
            usize::try_from(options.threads).unwrap_or(usize::MAX),
        )
        .map_err(invalid)?
    };
    table.set_column("splicing_lof", bool_column(&results, |r| r.lof));
    table.set_column(
        "splicing_len_changing",
        bool_column(&results, |r| r.len_changing),
    );
    table.set_column(
        "splicing_frameshift",
        bool_column(&results, |r| r.frameshift),
    );
    table.set_column(
        "splicing_span_intol_domain",
        bool_column(&results, |r| r.span_intol_domain),
    );
    table.set_column(
        "splicing_affected_exons",
        Column::Text(
            results
                .iter()
                .map(|r| Some(r.affected_exons.clone()))
                .collect(),
        ),
    );
    table.set_column(
        "splicing_ten_percent_protein",
        bool_column(&results, |r| r.ten_percent_protein),
    );
    table.set_column(
        "spliceai_canonical_gain_offset",
        offset_column(&results, |r| r.gain_offset),
    );
    table.set_column(
        "spliceai_canonical_loss_offset",
        offset_column(&results, |r| r.loss_offset),
    );
    table.set_column(
        "splicing_last_exon_only",
        bool_column(&results, |r| r.last_exon_only),
    );

    // Write beside the table and rename over it, so the table is never half written.
    let temporary = format!("{}.{}.splicing.anno.txt", options.anno_table, unique_tag());
    write_table(Path::new(&temporary), &table)?;
    std::fs::rename(&temporary, &options.anno_table)?;
    eprintln!("Results saved to {}", options.anno_table);
    Ok(())
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let options = match parse_options(args)? {
        Parsed::Help => {
            println!("{USAGE}");
            return Ok(());
        }
        Parsed::Run(options) => options,
    };
    // Like the Python script, a failure after the arguments exits with code 1.
    let outcome = cutoffs_for_config(options.config.as_deref())
        .and_then(|cutoffs| splicing_interpretation(&options, cutoffs));
    if let Err(error) = outcome {
        eprintln!("priva splicing-analysis: {error}");
        std::process::exit(1);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pandas_table::parse_table;

    const NORM_MAP: &str = include_str!("../../data/InterPro/InterPro_domain_norm_map.json");
    const ROWS: &str = include_str!("../tests/splicing_rows.tsv");
    const EXPECTED: &str = include_str!("../tests/splicing_rows_expected.tsv");

    // (hgvsc, strand, expected) from the Python function with var_pos 1000, var_length 1
    const POSITIONS: &[(
        &str,
        &str,
        Option<(VariantType, Option<i128>, Option<i128>, bool, bool)>,
    )] = &[
        (
            "ENST1:c.116+1G>A",
            "1",
            Some((VariantType::DonorRegion, Some(0), Some(999), false, true)),
        ),
        (
            "ENST1:c.1897-4A>G",
            "1",
            Some((
                VariantType::AcceptorRegion,
                Some(-2),
                Some(1004),
                false,
                false,
            )),
        ),
        (
            "ENST1:c.459+2dup",
            "1",
            Some((VariantType::DonorRegion, Some(0), Some(999), false, true)),
        ),
        (
            "ENST1:c.3354-1del",
            "1",
            Some((
                VariantType::AcceptorRegion,
                Some(0),
                Some(1002),
                false,
                true,
            )),
        ),
        (
            "ENST1:c.3646_3646+3del",
            "1",
            Some((VariantType::DonorRegion, Some(1), Some(998), false, true)),
        ),
        (
            "ENST1:c.207+1_207+2del",
            "1",
            Some((VariantType::DonorRegion, Some(0), Some(999), false, true)),
        ),
        (
            "ENST1:c.44-7_44-5del",
            "1",
            Some((
                VariantType::AcceptorRegion,
                Some(4),
                Some(1006),
                false,
                true,
            )),
        ),
        (
            "ENST1:c.4644+12_4644+13insTA",
            "1",
            Some((VariantType::DonorRegion, Some(11), Some(988), false, true)),
        ),
        (
            "ENST1:c.248-4_248-3insAAGTTTT",
            "1",
            Some((
                VariantType::AcceptorRegion,
                Some(-2),
                Some(1004),
                false,
                false,
            )),
        ),
        (
            "ENST1:c.2746-17_2746-14dup",
            "1",
            Some((
                VariantType::AcceptorRegion,
                Some(13),
                Some(1015),
                false,
                true,
            )),
        ),
        ("ENST1:c.100_101del", "1", None),
        ("ENST1:c.100G>A", "1", None),
        ("ENST1:c.*10+2A>G", "1", None),
        (
            "ENST1:c.-5+2A>G",
            "1",
            Some((VariantType::DonorRegion, Some(1), Some(998), false, true)),
        ),
        (
            "ENST1:c.10-3_10+4del",
            "1",
            Some((VariantType::Spanning, None, None, true, false)),
        ),
        ("ENST1:c.1_2insA", "1", None),
        (
            "ENST1:c.88+1_88+2delinsAT",
            "1",
            Some((VariantType::DonorRegion, Some(0), Some(999), false, true)),
        ),
        (
            "ENST1:c.5+1A>G:c.7-2del",
            "1",
            Some((VariantType::DonorRegion, Some(0), Some(999), false, true)),
        ),
        ("x", "1", None),
        (
            "ENST1:c.7-1_7-1insG",
            "1",
            Some((
                VariantType::AcceptorRegion,
                Some(1),
                Some(1001),
                false,
                true,
            )),
        ),
        (
            "ENST1:c.7_7+1insG",
            "1",
            Some((VariantType::DonorRegion, Some(-1), Some(1000), false, true)),
        ),
        (
            "ENST1:c.12+3_12-3dup",
            "1",
            Some((VariantType::Spanning, None, None, true, false)),
        ),
        (
            "ENST1:c.9-1G>T",
            "1",
            Some((
                VariantType::AcceptorRegion,
                Some(1),
                Some(1001),
                false,
                true,
            )),
        ),
        (
            "ENST1:c.1+9del",
            "1",
            Some((VariantType::DonorRegion, Some(7), Some(992), false, true)),
        ),
        (
            "ENST1:c.007+005A>G",
            "1",
            Some((VariantType::DonorRegion, Some(4), Some(995), false, true)),
        ),
        (
            "ENST1:c.150+10_151-10delinsA",
            "1",
            Some((VariantType::Spanning, None, None, true, false)),
        ),
        (
            "ENST1:c.500-2_500-1del",
            "1",
            Some((
                VariantType::AcceptorRegion,
                Some(0),
                Some(1002),
                false,
                true,
            )),
        ),
        (
            "ENST1:c.116+1G>A",
            "-1",
            Some((VariantType::DonorRegion, Some(1), Some(1001), false, true)),
        ),
        (
            "ENST1:c.1897-4A>G",
            "-1",
            Some((VariantType::AcceptorRegion, Some(3), Some(996), false, true)),
        ),
        (
            "ENST1:c.459+2dup",
            "-1",
            Some((VariantType::DonorRegion, Some(-1), Some(1003), false, true)),
        ),
        (
            "ENST1:c.3354-1del",
            "-1",
            Some((
                VariantType::AcceptorRegion,
                Some(-1),
                Some(1000),
                false,
                true,
            )),
        ),
        (
            "ENST1:c.3646_3646+3del",
            "-1",
            Some((VariantType::DonorRegion, Some(-2), Some(1004), false, false)),
        ),
        (
            "ENST1:c.207+1_207+2del",
            "-1",
            Some((VariantType::DonorRegion, Some(-1), Some(1003), false, true)),
        ),
        (
            "ENST1:c.44-7_44-5del",
            "-1",
            Some((
                VariantType::AcceptorRegion,
                Some(3),
                Some(1004),
                false,
                true,
            )),
        ),
        (
            "ENST1:c.4644+12_4644+13insTA",
            "-1",
            Some((
                VariantType::DonorRegion,
                Some(-11),
                Some(1013),
                false,
                false,
            )),
        ),
        (
            "ENST1:c.248-4_248-3insAAGTTTT",
            "-1",
            Some((VariantType::AcceptorRegion, Some(2), Some(997), false, true)),
        ),
        (
            "ENST1:c.2746-17_2746-14dup",
            "-1",
            Some((
                VariantType::AcceptorRegion,
                Some(12),
                Some(1013),
                false,
                true,
            )),
        ),
        ("ENST1:c.100_101del", "-1", None),
        ("ENST1:c.100G>A", "-1", None),
        ("ENST1:c.*10+2A>G", "-1", None),
        (
            "ENST1:c.-5+2A>G",
            "-1",
            Some((VariantType::DonorRegion, Some(0), Some(1002), false, true)),
        ),
        (
            "ENST1:c.10-3_10+4del",
            "-1",
            Some((VariantType::Spanning, None, None, true, false)),
        ),
        ("ENST1:c.1_2insA", "-1", None),
        (
            "ENST1:c.88+1_88+2delinsAT",
            "-1",
            Some((VariantType::DonorRegion, Some(-1), Some(1003), false, true)),
        ),
        (
            "ENST1:c.5+1A>G:c.7-2del",
            "-1",
            Some((VariantType::DonorRegion, Some(1), Some(1001), false, true)),
        ),
        ("x", "-1", None),
        (
            "ENST1:c.7-1_7-1insG",
            "-1",
            Some((VariantType::AcceptorRegion, Some(0), Some(999), false, true)),
        ),
        (
            "ENST1:c.7_7+1insG",
            "-1",
            Some((VariantType::DonorRegion, Some(1), Some(1001), false, true)),
        ),
        (
            "ENST1:c.12+3_12-3dup",
            "-1",
            Some((VariantType::Spanning, None, None, true, false)),
        ),
        (
            "ENST1:c.9-1G>T",
            "-1",
            Some((VariantType::AcceptorRegion, Some(0), Some(999), false, true)),
        ),
        (
            "ENST1:c.1+9del",
            "-1",
            Some((VariantType::DonorRegion, Some(-8), Some(1010), false, false)),
        ),
        (
            "ENST1:c.007+005A>G",
            "-1",
            Some((VariantType::DonorRegion, Some(-3), Some(1005), false, false)),
        ),
        (
            "ENST1:c.150+10_151-10delinsA",
            "-1",
            Some((VariantType::Spanning, None, None, true, false)),
        ),
        (
            "ENST1:c.500-2_500-1del",
            "-1",
            Some((
                VariantType::AcceptorRegion,
                Some(-1),
                Some(1000),
                false,
                true,
            )),
        ),
    ];
    const POSITION_ERRORS: &[(&str, &str)] =
        &[("ENST1:c.100+0A>G", "1"), ("ENST1:c.100+0A>G", "-1")];
    // strand 0: ENST1:c.10-3_10-1del raises ValueError
    // strand 0: ENST1:c.116+1G>A None
    // strand 0: ENST1:c.3354-1del None
    // (hgvsc, consequence, expected) from _touches_canonical_splice_site
    const TOUCHES: &[(&str, &str, bool)] = &[
        ("ENST1:c.116+1G>A", "splice_donor_variant", true),
        ("ENST1:c.116+1G>A", "splice_acceptor_variant", false),
        (
            "ENST1:c.116+1G>A",
            "splice_donor_variant&splice_acceptor_variant",
            true,
        ),
        ("ENST1:c.116+1G>A", "intron_variant", false),
        ("ENST1:c.116+2G>A", "splice_donor_variant", true),
        ("ENST1:c.116+2G>A", "splice_acceptor_variant", false),
        (
            "ENST1:c.116+2G>A",
            "splice_donor_variant&splice_acceptor_variant",
            true,
        ),
        ("ENST1:c.116+2G>A", "intron_variant", false),
        ("ENST1:c.116+3G>A", "splice_donor_variant", false),
        ("ENST1:c.116+3G>A", "splice_acceptor_variant", false),
        (
            "ENST1:c.116+3G>A",
            "splice_donor_variant&splice_acceptor_variant",
            false,
        ),
        ("ENST1:c.116+3G>A", "intron_variant", false),
        ("ENST1:c.116-1G>A", "splice_donor_variant", false),
        ("ENST1:c.116-1G>A", "splice_acceptor_variant", true),
        (
            "ENST1:c.116-1G>A",
            "splice_donor_variant&splice_acceptor_variant",
            false,
        ),
        ("ENST1:c.116-1G>A", "intron_variant", true),
        ("ENST1:c.116-2G>A", "splice_donor_variant", false),
        ("ENST1:c.116-2G>A", "splice_acceptor_variant", true),
        (
            "ENST1:c.116-2G>A",
            "splice_donor_variant&splice_acceptor_variant",
            false,
        ),
        ("ENST1:c.116-2G>A", "intron_variant", true),
        ("ENST1:c.116-3G>A", "splice_donor_variant", false),
        ("ENST1:c.116-3G>A", "splice_acceptor_variant", false),
        (
            "ENST1:c.116-3G>A",
            "splice_donor_variant&splice_acceptor_variant",
            false,
        ),
        ("ENST1:c.116-3G>A", "intron_variant", false),
        ("ENST1:c.116+5_116+8del", "splice_donor_variant", false),
        ("ENST1:c.116+5_116+8del", "splice_acceptor_variant", false),
        (
            "ENST1:c.116+5_116+8del",
            "splice_donor_variant&splice_acceptor_variant",
            false,
        ),
        ("ENST1:c.116+5_116+8del", "intron_variant", false),
        ("ENST1:c.116-5_116+2del", "splice_donor_variant", true),
        ("ENST1:c.116-5_116+2del", "splice_acceptor_variant", false),
        (
            "ENST1:c.116-5_116+2del",
            "splice_donor_variant&splice_acceptor_variant",
            true,
        ),
        ("ENST1:c.116-5_116+2del", "intron_variant", false),
        ("ENST1:c.100_116+1del", "splice_donor_variant", true),
        ("ENST1:c.100_116+1del", "splice_acceptor_variant", true),
        (
            "ENST1:c.100_116+1del",
            "splice_donor_variant&splice_acceptor_variant",
            true,
        ),
        ("ENST1:c.100_116+1del", "intron_variant", true),
        ("ENST1:c.100_101del", "splice_donor_variant", false),
        ("ENST1:c.100_101del", "splice_acceptor_variant", false),
        (
            "ENST1:c.100_101del",
            "splice_donor_variant&splice_acceptor_variant",
            false,
        ),
        ("ENST1:c.100_101del", "intron_variant", false),
        ("ENST1:c.-5+1A>G", "splice_donor_variant", false),
        ("ENST1:c.-5+1A>G", "splice_acceptor_variant", false),
        (
            "ENST1:c.-5+1A>G",
            "splice_donor_variant&splice_acceptor_variant",
            false,
        ),
        ("ENST1:c.-5+1A>G", "intron_variant", false),
        ("ENST1:c.*5+1A>G", "splice_donor_variant", false),
        ("ENST1:c.*5+1A>G", "splice_acceptor_variant", false),
        (
            "ENST1:c.*5+1A>G",
            "splice_donor_variant&splice_acceptor_variant",
            false,
        ),
        ("ENST1:c.*5+1A>G", "intron_variant", false),
        ("ENST1:c.116+01G>A", "splice_donor_variant", true),
        ("ENST1:c.116+01G>A", "splice_acceptor_variant", false),
        (
            "ENST1:c.116+01G>A",
            "splice_donor_variant&splice_acceptor_variant",
            true,
        ),
        ("ENST1:c.116+01G>A", "intron_variant", false),
        ("ENST1:c.116+1_116+2del", "splice_donor_variant", true),
        ("ENST1:c.116+1_116+2del", "splice_acceptor_variant", false),
        (
            "ENST1:c.116+1_116+2del",
            "splice_donor_variant&splice_acceptor_variant",
            true,
        ),
        ("ENST1:c.116+1_116+2del", "intron_variant", false),
        ("ENST1:c.9-2_10del", "splice_donor_variant", true),
        ("ENST1:c.9-2_10del", "splice_acceptor_variant", true),
        (
            "ENST1:c.9-2_10del",
            "splice_donor_variant&splice_acceptor_variant",
            true,
        ),
        ("ENST1:c.9-2_10del", "intron_variant", true),
        ("ENST1:c.9-2_9-1_10del", "splice_donor_variant", false),
        ("ENST1:c.9-2_9-1_10del", "splice_acceptor_variant", true),
        (
            "ENST1:c.9-2_9-1_10del",
            "splice_donor_variant&splice_acceptor_variant",
            false,
        ),
        ("ENST1:c.9-2_9-1_10del", "intron_variant", true),
        ("ENST1:c.88+0A>G", "splice_donor_variant", false),
        ("ENST1:c.88+0A>G", "splice_acceptor_variant", false),
        (
            "ENST1:c.88+0A>G",
            "splice_donor_variant&splice_acceptor_variant",
            false,
        ),
        ("ENST1:c.88+0A>G", "intron_variant", false),
        ("", "splice_donor_variant", false),
        ("", "splice_acceptor_variant", false),
        ("", "splice_donor_variant&splice_acceptor_variant", false),
        ("", "intron_variant", false),
        ("no colon c.116+1G>A", "splice_donor_variant", false),
        ("no colon c.116+1G>A", "splice_acceptor_variant", false),
        (
            "no colon c.116+1G>A",
            "splice_donor_variant&splice_acceptor_variant",
            false,
        ),
        ("no colon c.116+1G>A", "intron_variant", false),
    ];

    fn strings(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| (*value).to_owned()).collect()
    }

    fn resources() -> Resources {
        let mut exon_domains: HashMap<String, HashMap<String, Vec<String>>> = HashMap::new();
        exon_domains.insert(
            "ENST_A".to_owned(),
            HashMap::from([
                ("2".to_owned(), strings(&["GENE:Pfam:PF1", "GENE:Pfam:PF2"])),
                ("3".to_owned(), strings(&["GENE:Pfam:PF3"])),
                ("5".to_owned(), strings(&["GENE:SMART:SM1"])),
            ]),
        );
        exon_domains.insert(
            "ENST_B".to_owned(),
            HashMap::from([
                ("1".to_owned(), strings(&["GENE:Pfam:PF9"])),
                ("4".to_owned(), strings(&["GENE:Pfam:PF1"])),
            ]),
        );
        Resources {
            exon_domains,
            intolerant_domains: strings(&["GENE:Pfam:PF1", "GENE:SMART:SM1"])
                .into_iter()
                .collect(),
            interpro: InterproMapping::from_entries(Vec::new()),
            normalizer: DomainNormalizer::from_parts(NORM_MAP, HashMap::new()).unwrap(),
        }
    }

    #[test]
    fn hgvsc_positions_match_the_python_function() {
        for (hgvsc, strand, expected) in POSITIONS {
            let found = parse_hgvsc_splice_position(Some(*hgvsc), strand, 1000, 1).unwrap();
            let shown = found.map(|position| {
                (
                    position.variant_type,
                    position.distance_to_canonical,
                    position.exon_boundary_pos,
                    position.exon_skipping,
                    position.overlapping_canonical_site,
                )
            });
            assert_eq!(&shown, expected, "{hgvsc} on strand {strand}");
        }
        for (hgvsc, strand) in POSITION_ERRORS {
            assert!(
                parse_hgvsc_splice_position(Some(*hgvsc), strand, 1000, 1).is_err(),
                "{hgvsc} on strand {strand}"
            );
        }
        // A strand that is neither 1 nor -1 only stops the range patterns.
        assert!(parse_hgvsc_splice_position(Some("ENST1:c.10-3_10-1del"), "0", 1000, 1).is_err());
        assert_eq!(
            parse_hgvsc_splice_position(Some("ENST1:c.116+1G>A"), "0", 1000, 1),
            Ok(None)
        );
        assert_eq!(parse_hgvsc_splice_position(None, "1", 1000, 1), Ok(None));
    }

    #[test]
    fn canonical_splice_dinucleotide_check_matches_the_python_function() {
        for (hgvsc, consequence, expected) in TOUCHES {
            assert_eq!(
                touches_canonical_splice_site(Some(*hgvsc), consequence),
                *expected,
                "{hgvsc} with {consequence}"
            );
        }
        assert!(!touches_canonical_splice_site(None, "splice_donor_variant"));
    }

    #[test]
    fn rows_match_the_python_function() {
        let table = parse_table(ROWS).unwrap();
        let resources = resources();
        let context = Context {
            table: &table,
            columns: ColumnIndexes::find(&table),
            resources: &resources,
            splice_altering_cutoff: 0.5,
            acceptor_rescue_cutoff: 0.25,
        };
        let expected: Vec<&str> = EXPECTED.lines().collect();
        assert_eq!(expected.len(), table.rows);
        let flag = |value: bool| if value { "True" } else { "False" };
        let offset = |value: Offset| match value {
            Offset::Zero => "0".to_owned(),
            Offset::Value(number) => float_repr(number),
        };
        let mut stopped = 0;
        for (row, line) in expected.iter().enumerate() {
            match analyze_row(&context, row) {
                Ok(result) => {
                    let shown = [
                        flag(result.lof).to_owned(),
                        flag(result.len_changing).to_owned(),
                        flag(result.frameshift).to_owned(),
                        flag(result.span_intol_domain).to_owned(),
                        flag(result.ten_percent_protein).to_owned(),
                        offset(result.gain_offset),
                        offset(result.loss_offset),
                        result.affected_exons.clone(),
                        flag(result.last_exon_only).to_owned(),
                    ]
                    .join("\t");
                    assert_eq!(&shown, line, "data row {}", row + 1);
                }
                Err(error) => {
                    stopped += 1;
                    assert!(
                        line.starts_with("ERROR"),
                        "data row {} stopped ({error}) but Python gave {line}",
                        row + 1
                    );
                }
            }
        }
        assert_eq!(
            stopped,
            expected
                .iter()
                .filter(|line| line.starts_with("ERROR"))
                .count()
        );
    }

    #[test]
    fn offset_columns_are_int_only_when_every_offset_is_zero() {
        let result = |gain: Offset| RowResult {
            lof: false,
            len_changing: false,
            frameshift: false,
            span_intol_domain: false,
            ten_percent_protein: false,
            gain_offset: gain,
            loss_offset: Offset::Zero,
            affected_exons: String::new(),
            last_exon_only: false,
        };
        let zeros = [result(Offset::Zero), result(Offset::Zero)];
        assert_eq!(
            offset_column(&zeros, |r| r.gain_offset),
            Column::Int(vec![0, 0])
        );
        let mixed = [result(Offset::Zero), result(Offset::Value(-3.0))];
        assert_eq!(
            offset_column(&mixed, |r| r.gain_offset),
            Column::Float(vec![0.0, -3.0])
        );
    }

    #[test]
    fn options_are_read_like_argparse() {
        let args = |items: &[&str]| -> Vec<String> { strings(items) };
        let required = [
            "--anno_table",
            "t.tsv",
            "--transcript_domain_map",
            "m.pkl",
            "--intolerant_domains",
            "d.pkl.gz",
            "--interpro_entry_map_pkl",
            "i.pkl.gz",
        ];
        let Parsed::Run(options) = parse_options(&args(&required)).unwrap() else {
            panic!("a complete command line runs")
        };
        assert_eq!(options.threads, 12);
        assert_eq!(options.config, None);
        let mut more = required.to_vec();
        more.extend(["--thr", "-1", "--conf=c.yaml"]);
        let Parsed::Run(options) = parse_options(&args(&more)).unwrap() else {
            panic!()
        };
        assert_eq!(options.threads, -1);
        assert_eq!(options.config.as_deref(), Some("c.yaml"));
        assert!(matches!(parse_options(&args(&["-h"])), Ok(Parsed::Help)));
        assert!(matches!(parse_options(&args(&["--he"])), Ok(Parsed::Help)));
        assert!(parse_options(&args(&["--anno_table", "t.tsv"])).is_err());
        assert!(parse_options(&args(&["--bogus", "x"])).is_err());
        assert!(parse_options(&args(&["--an", "t", "--anno_table"])).is_err());
        let mut bad = required.to_vec();
        bad.extend(["--threads", "four"]);
        assert!(parse_options(&args(&bad)).is_err());
        let mut dash = required.to_vec();
        dash.extend(["--threads", "-x"]);
        assert!(parse_options(&args(&dash)).is_err());
        assert!(
            is_negative_number("-5") && is_negative_number("-.5") && is_negative_number("-5.5")
        );
        assert!(!is_negative_number("-") && !is_negative_number("-5.") && !is_negative_number("5"));
    }
}
