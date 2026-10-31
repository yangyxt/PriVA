//! `scripts/acmg_consequence.py` (196 lines): the consequence helpers of the ACMG step that were
//! not yet in Rust (`coerce_pext_value` is `crate::acmg_consequence`):
//!
//! - `vep_consq_interpret` (with `vep_consq_interpret_per_row`): the columns `vep_consq_lof` and
//!   `vep_consq_len_changing` from VEP's `Consequence`, LOFTEE's `LoF` and `NMD`;
//! - `get_variant_type`: the HGVSp variant type (nonsense, frameshift, delins, deletion,
//!   insertion_type, missense), used by PS1/PM5;
//! - `truncate_fraction`: the share of the protein a variant removes or spans, used by PVS1.
//!
//! The Python regular expressions of `get_variant_type` are written out by hand for ASCII text
//! (`\w` is `[A-Za-z0-9_]` and `\d` is `[0-9]` there); a value holding another character is
//! refused, because Python's and Rust's Unicode classes differ. `pd.to_numeric(errors="coerce")`
//! of a text is pandas' float reading (`precise_xstrtod`) with Python's `int()` for integers;
//! integers beyond 2**53 are refused (pandas would mix exact integers with rounded floats).
use super::acmg_pp3_bp4_bp7_insilico::pandas_masks::{python_str, text_contains};
use super::types::{AcmgError, HarnessCall, not_converted};
use crate::cli::{Result, invalid};
use crate::pandas_table::{Cell, Column, Table, float_token};
use crate::python_text::{python_int, python_strip};

/// What `vep_consq_interpret_per_row` reads from one row; `None` is an absent column.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct ConsequenceRow<'a> {
    pub(crate) consequence: Option<Cell<'a>>,
    pub(crate) lof: Option<Cell<'a>>,
    pub(crate) nmd: Option<Cell<'a>>,
}

/// Python's `str(value)` of a row value (`None` for an absent column gives `"None"`).
fn row_str(value: Option<Cell<'_>>) -> String {
    match value {
        None => "None".to_owned(),
        Some(cell) => python_str(cell),
    }
}

/// Python's truth value of a row value (`value or ""`): NaN is true, as are non-empty text,
/// non-zero numbers and True.
fn truthy(value: Option<Cell<'_>>) -> bool {
    match value {
        None => false,
        Some(Cell::Missing) => true,
        Some(Cell::Text(text)) => !text.is_empty(),
        Some(Cell::Int(number)) => number != 0,
        Some(Cell::Float(number)) => number != 0.0,
        Some(Cell::Bool(flag)) => flag,
    }
}

const LOF_CONSEQUENCES: [&str; 4] = [
    "stop_gained",
    "start_lost",
    "transcript_ablation",
    "frameshift_variant",
];

const LENGTH_CHANGING_CONSEQUENCES: [&str; 5] = [
    "inframe_insertion",
    "inframe_deletion",
    "feature_elongation",
    "feature_truncation",
    "stop_lost",
];

/// `vep_consq_interpret_per_row(row)`: `(is_lof, is_length_changing)`.
pub(crate) fn vep_consq_interpret_per_row(row: &ConsequenceRow<'_>) -> (bool, bool) {
    // `"escaping" in str(row.get('NMD', ""))`: an absent column gives "".
    let nmd_escaping = match row.nmd {
        None => false,
        Some(cell) => python_str(cell).contains("escaping"),
    };
    let Some(Cell::Text(consequence)) = row.consequence else {
        return (false, false);
    };
    // `re.split(r"[&,;|]", str(loftee_result or ""))`, stripped and upper-cased tokens.
    let lof_text = if truthy(row.lof) {
        row_str(row.lof)
    } else {
        String::new()
    };
    let mut high_confidence = false;
    let mut low_confidence = false;
    for token in lof_text.split(['&', ',', ';', '|']) {
        let token = python_strip(token);
        if token.is_empty() {
            continue;
        }
        let upper = token.to_uppercase();
        high_confidence |= upper == "HC" || upper == "OS";
        low_confidence |= upper == "LC";
    }
    let consequence_lof = LOF_CONSEQUENCES
        .iter()
        .any(|name| consequence.contains(name));
    let is_lof = high_confidence || (consequence_lof && !low_confidence && !nmd_escaping);
    let is_length_changing = is_lof
        || LENGTH_CHANGING_CONSEQUENCES
            .iter()
            .any(|name| consequence.contains(name))
        || consequence_lof
        || nmd_escaping;
    (is_lof, is_length_changing)
}

/// `vep_consq_interpret(df, threads)`: writes `vep_consq_lof` and `vep_consq_len_changing`
/// (bool columns). `threads` sizes Python's process pool, which refuses fewer than one; an
/// empty table fails where Python unpacks the empty result.
pub(crate) fn vep_consq_interpret(table: &mut Table, threads: i64) -> Result<()> {
    if threads < 1 {
        return Err(Box::new(AcmgError::new(
            "ValueError",
            "Number of processes must be at least 1",
        )));
    }
    if table.rows == 0 {
        return Err(Box::new(AcmgError::new(
            "ValueError",
            "not enough values to unpack (expected 2, got 0)",
        )));
    }
    for name in ["Consequence", "LoF", "NMD"] {
        if let Ok(column) = table.column(name)
            && matches!(
                column,
                Column::NullableInt(_) | Column::NullableBool(_) | Column::Object(_)
            )
        {
            return Err(invalid(format!(
                "refused (not reproduced): {name} is a nullable or mixed pandas column \
                 (Python's None and NaN differ in `value or \"\"`)"
            )));
        }
    }
    let position = |name: &str| table.column_index(name);
    let (consequence, lof, nmd) = (position("Consequence"), position("LoF"), position("NMD"));
    let mut flags_lof = Vec::with_capacity(table.rows);
    let mut flags_length = Vec::with_capacity(table.rows);
    for row in 0..table.rows {
        let cell = |index: Option<usize>| index.map(|index| table.cell(index, row));
        let (is_lof, is_length) = vep_consq_interpret_per_row(&ConsequenceRow {
            consequence: cell(consequence),
            lof: cell(lof),
            nmd: cell(nmd),
        });
        flags_lof.push(is_lof);
        flags_length.push(is_length);
    }
    table.set_column("vep_consq_lof", Column::Bool(flags_lof));
    table.set_column("vep_consq_len_changing", Column::Bool(flags_length));
    Ok(())
}

fn is_word(byte: u8) -> bool {
    byte.is_ascii_alphanumeric() || byte == b'_'
}

/// `re.search(r"p\.\w+\d+" + suffix, text)` for ASCII text, where `suffix(rest)` says whether
/// the rest of the text, right after the digits, starts with the suffix's match.
fn protein_change_then(text: &[u8], suffix: impl Fn(&[u8]) -> bool) -> bool {
    let mut start = 0;
    while let Some(found) = text[start..].windows(2).position(|pair| pair == b"p.") {
        let begin = start + found + 2;
        let mut end = begin;
        while end < text.len() && is_word(text[end]) {
            end += 1;
        }
        // `\w+` takes at least one character and `\d+` at least one digit: the digits end
        // at `stop`, which lies inside the run of word characters.
        for stop in begin + 2..=end {
            if text[stop - 1].is_ascii_digit() && suffix(&text[stop..]) {
                return true;
            }
        }
        start += found + 1;
    }
    false
}

/// `re.search(r'p\.[A-Z][a-z]+\d+[A-Z][a-z]+$', text)` for ASCII text (`$` also matches
/// before a final line break).
fn missense_change(text: &[u8]) -> bool {
    let body = text.strip_suffix(b"\n").unwrap_or(text);
    let mut start = 0;
    while let Some(found) = body[start..].windows(2).position(|pair| pair == b"p.") {
        let rest = &body[start + found + 2..];
        let upper_then_lower = |part: &[u8]| -> usize {
            if part.first().is_some_and(u8::is_ascii_uppercase) {
                let lower = part[1..]
                    .iter()
                    .take_while(|byte| byte.is_ascii_lowercase())
                    .count();
                if lower > 0 { 1 + lower } else { 0 }
            } else {
                0
            }
        };
        let first = upper_then_lower(rest);
        if first > 0 {
            let digits = rest[first..]
                .iter()
                .take_while(|byte| byte.is_ascii_digit())
                .count();
            if digits > 0 {
                let tail = &rest[first + digits..];
                let second = upper_then_lower(tail);
                if second > 0 && second == tail.len() {
                    return true;
                }
            }
        }
        start += found + 1;
    }
    false
}

/// `get_variant_type(hgvsp)` for a text: the part after the last `:` is searched for the
/// nine HGVSp patterns in Python's order. Text with a non-ASCII character is refused.
pub(crate) fn get_variant_type(hgvsp: &str) -> Result<Option<&'static str>> {
    if !hgvsp.is_ascii() {
        return Err(invalid(format!(
            "refused (not reproduced): the HGVSp {hgvsp:?} holds a non-ASCII character, where \
             Python's \\w and \\d differ from Rust's"
        )));
    }
    let text = hgvsp.rsplit(':').next().unwrap_or(hgvsp).as_bytes();
    let literal = |word: &'static [u8]| move |rest: &[u8]| rest.starts_with(word);
    let frameshift_ter = |rest: &[u8]| {
        rest.len() >= 6
            && matches!(rest[0], b'F' | b'f')
            && rest[1] == b's'
            && matches!(rest[2], b'T' | b't')
            && &rest[3..5] == b"er"
            && rest[5].is_ascii_digit()
    };
    let found = if protein_change_then(text, literal(b"Ter"))
        || protein_change_then(text, literal(b"*"))
    {
        Some("nonsense")
    } else if protein_change_then(text, literal(b"fs")) || protein_change_then(text, frameshift_ter)
    {
        Some("frameshift")
    } else if protein_change_then(text, literal(b"delins")) {
        Some("delins")
    } else if protein_change_then(text, literal(b"del")) {
        Some("deletion")
    } else if protein_change_then(text, literal(b"dup"))
        || protein_change_then(text, literal(b"ins"))
    {
        Some("insertion_type")
    } else if missense_change(text) {
        Some("missense")
    } else {
        None
    };
    Ok(found)
}

/// `get_variant_type(value)` for a table value: `None` for anything but text.
pub(crate) fn get_variant_type_of(cell: Cell<'_>) -> Result<Option<&'static str>> {
    match cell {
        Cell::Text(text) => get_variant_type(text),
        _ => Ok(None),
    }
}

/// A number `pd.to_numeric(text, errors="coerce")` gives.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Numeric {
    Int(i64),
    /// NaN where the text is not a number.
    Float(f64),
}

impl Numeric {
    pub(crate) fn as_f64(self) -> f64 {
        match self {
            Numeric::Int(value) => value as f64,
            Numeric::Float(value) => value,
        }
    }
}

/// The largest integer every step of the computation holds exactly as a double.
const EXACT_INTEGER: i128 = 1 << 53;

/// `pd.to_numeric(text, errors="coerce")` of one text (pandas' `maybe_convert_numeric`):
/// pandas' float reading, and for a text without a decimal point or exponent Python's
/// `int()` of it. Integers beyond 2**53 are refused.
pub(crate) fn to_numeric_text(text: &str) -> Result<Numeric> {
    if text.is_empty() {
        return Ok(Numeric::Float(f64::NAN));
    }
    let Some(value) = float_token(text) else {
        return Ok(Numeric::Float(f64::NAN));
    };
    let maybe_int = value.is_finite() && !text.contains(['.', 'e', 'E']);
    if !maybe_int {
        return Ok(Numeric::Float(value));
    }
    match python_int(text) {
        Some(integer) if (-EXACT_INTEGER..=EXACT_INTEGER).contains(&integer) => {
            Ok(Numeric::Int(integer as i64))
        }
        _ => Err(invalid(format!(
            "refused (not reproduced): {text:?} is an integer beyond 2**53 for pd.to_numeric"
        ))),
    }
}

/// `pd.to_numeric(x, errors="coerce")` of `x.split("/")[0].split("-")` parts, and the
/// `abs(end - start + 1)` the span takes.
fn span_of(text: &str) -> Result<f64> {
    let parts: Vec<&str> = text.split('/').next().unwrap_or("").split('-').collect();
    if parts.len() < 2 {
        return Err(Box::new(AcmgError::new(
            "IndexError",
            "list index out of range",
        )));
    }
    let end = to_numeric_text(parts[1])?;
    let start = to_numeric_text(parts[0])?;
    Ok(match (end, start) {
        (Numeric::Int(end), Numeric::Int(start)) => (end - start + 1).abs() as f64,
        (end, start) => (end.as_f64() - start.as_f64() + 1.0).abs(),
    })
}

/// `truncate_fraction(df)`: where `Consequence` holds `stop_gained` or `frameshift`,
/// `1 - position / length`, else `span / length`, from `Protein_position` (`"pos/length"`,
/// `"start-end/length"`). The length is the last field of `.str.split("/", expand=True)`,
/// so a row with fewer `/`-fields than the longest row has none (NaN).
pub(crate) fn truncate_fraction(table: &Table) -> Result<Vec<f64>> {
    if table.rows == 0 {
        return Err(invalid(
            "refused (not reproduced): truncate_fraction of an empty table",
        ));
    }
    let mut found = text_contains(table, "Consequence", None, &["stop_gained", "frameshift"])?;
    let frameshift = found.pop().expect("two patterns");
    let stop_gained = found.pop().expect("two patterns");
    // pandas' `|` of two `str.contains` results is a plain bool Series.
    let nonsense: Vec<bool> = stop_gained.or(&frameshift).plain()?.to_vec();

    let column = table.column("Protein_position")?;
    if matches!(
        column,
        Column::NullableInt(_) | Column::NullableBool(_) | Column::Object(_)
    ) {
        return Err(invalid(
            "refused (not reproduced): Protein_position is a nullable or mixed pandas column",
        ));
    }
    // `df["Protein_position"].astype(str)`.
    let texts: Vec<String> = (0..table.rows)
        .map(|row| python_str(column.cell(row)))
        .collect();
    // `.map(lambda x: ... if "-" in str(x) else 1)` on the column itself.
    let mut spans = Vec::with_capacity(table.rows);
    for (row, text) in texts.iter().enumerate() {
        if !text.contains('-') {
            spans.push(1.0);
            continue;
        }
        match column.cell(row) {
            Cell::Text(value) => spans.push(span_of(value)?),
            other => {
                let kind = match other {
                    Cell::Int(_) => "int",
                    Cell::Bool(_) => "bool",
                    _ => "float",
                };
                return Err(Box::new(AcmgError::new(
                    "AttributeError",
                    format!("'{kind}' object has no attribute 'split'"),
                )));
            }
        }
    }
    // The first `-`-field, then its first `/`-field.
    let positions: Vec<f64> = texts
        .iter()
        .map(|text| {
            let first = text.split('-').next().unwrap_or("");
            to_numeric_text(first.split('/').next().unwrap_or("")).map(Numeric::as_f64)
        })
        .collect::<Result<_>>()?;
    let field_count = texts
        .iter()
        .map(|text| text.split('/').count())
        .max()
        .unwrap_or(1);
    let lengths: Vec<f64> = texts
        .iter()
        .map(|text| match text.split('/').nth(field_count - 1) {
            Some(last) => to_numeric_text(last).map(Numeric::as_f64),
            None => Ok(f64::NAN),
        })
        .collect::<Result<_>>()?;
    Ok((0..table.rows)
        .map(|row| {
            if nonsense[row] {
                1.0 - positions[row] / lengths[row]
            } else {
                spans[row] / lengths[row]
            }
        })
        .collect())
}

/// `priva acmg-step acmg_consequence <function>`:
///
/// - `vep_consq_interpret`, option `threads`: the table with the two columns written;
/// - `get_variant_type`: for each value of the column `hgvsp`, the type as a column
///   `variant_type` (None where Python returns None);
/// - `truncate_fraction`: a frame of the float column `truncate_fraction`.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    match function {
        "vep_consq_interpret" => {
            let threads = match call.options.get("threads") {
                Some(crate::json_text::Json::Int(number)) => {
                    i64::try_from(*number).map_err(|_| invalid("threads must fit in 64 bits"))?
                }
                None => 10,
                Some(_) => return Err(invalid("the option threads is an integer")),
            };
            let mut table = call.table.clone();
            vep_consq_interpret(&mut table, threads)?;
            Ok(table)
        }
        "get_variant_type" => {
            let column = call.table.column("hgvsp")?;
            let types = (0..call.table.rows)
                .map(|row| {
                    get_variant_type_of(column.cell(row)).map(|found| found.map(str::to_owned))
                })
                .collect::<Result<Vec<_>>>()?;
            Ok(Table {
                names: vec!["variant_type".to_owned()],
                columns: vec![crate::pandas_table::Column::Object(
                    types
                        .into_iter()
                        .map(|found| match found {
                            Some(text) => crate::pandas_table::ObjectValue::Text(text),
                            None => crate::pandas_table::ObjectValue::None,
                        })
                        .collect(),
                )],
                rows: call.table.rows,
            })
        }
        "truncate_fraction" => {
            let values = truncate_fraction(&call.table)?;
            Ok(Table {
                names: vec!["truncate_fraction".to_owned()],
                rows: values.len(),
                columns: vec![Column::Float(values)],
            })
        }
        _ => Err(not_converted("acmg_consequence", function)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row<'a>(
        consequence: Option<Cell<'a>>,
        lof: Option<Cell<'a>>,
        nmd: Option<Cell<'a>>,
    ) -> ConsequenceRow<'a> {
        ConsequenceRow {
            consequence,
            lof,
            nmd,
        }
    }

    #[test]
    fn per_row_flags_follow_the_python_rules() {
        let text = |value: &'static str| Some(Cell::Text(value));
        // Not text: both false.
        assert_eq!(
            vep_consq_interpret_per_row(&row(Some(Cell::Missing), text("HC"), None)),
            (false, false)
        );
        assert_eq!(
            vep_consq_interpret_per_row(&row(None, None, None)),
            (false, false)
        );
        // LOFTEE HC or OS alone makes a LoF call, in any case and among other tokens.
        assert_eq!(
            vep_consq_interpret_per_row(&row(text("intron_variant"), text(" os & lc"), None)),
            (true, true)
        );
        // A LoF consequence is vetoed by LC and by NMD escape, and stays length changing.
        assert_eq!(
            vep_consq_interpret_per_row(&row(text("stop_gained"), text("LC"), None)),
            (false, true)
        );
        assert_eq!(
            vep_consq_interpret_per_row(&row(
                text("frameshift_variant"),
                None,
                text("NMD_escaping_variant")
            )),
            (false, true)
        );
        // NaN in LoF is the text "nan" (truthy); 0 and False are "".
        assert_eq!(
            vep_consq_interpret_per_row(&row(text("start_lost"), Some(Cell::Missing), None)),
            (true, true)
        );
        assert_eq!(
            vep_consq_interpret_per_row(&row(text("stop_lost"), Some(Cell::Int(0)), None)),
            (false, true)
        );
        assert_eq!(
            vep_consq_interpret_per_row(&row(text("missense_variant"), None, None)),
            (false, false)
        );
    }

    #[test]
    fn variant_types_follow_the_python_patterns() {
        let kind = |text: &str| get_variant_type(text).unwrap();
        assert_eq!(kind("ENSP00000439902.1:p.Arg123Ter"), Some("nonsense"));
        assert_eq!(kind("p.Arg123*"), Some("nonsense"));
        // The digits must come right before `fs`, so VEP's usual frameshift spelling is no
        // frameshift here (nor anything else), as in Python.
        assert_eq!(kind("p.Asn3264LeufsTer12"), None);
        assert_eq!(kind("ENSP00000439902.1:p.Asn3264LeufsTer12"), None);
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
        assert_eq!(kind("p.A123G"), None);
        // `\w+\d+` needs a word character before the digits.
        assert_eq!(kind("p.1Ter"), None);
        assert_eq!(kind("p.11Ter"), Some("nonsense"));
        assert_eq!(kind("xp.Met1?p.Arg5Ter"), Some("nonsense"));
        assert!(get_variant_type("p.Ala١٢٣Gly").is_err());
        assert_eq!(get_variant_type_of(Cell::Missing).unwrap(), None);
    }

    #[test]
    fn to_numeric_reads_like_pandas() {
        assert_eq!(to_numeric_text("12").unwrap(), Numeric::Int(12));
        assert_eq!(to_numeric_text(" +007 ").unwrap(), Numeric::Int(7));
        assert_eq!(to_numeric_text("1.5").unwrap(), Numeric::Float(1.5));
        assert_eq!(to_numeric_text("1e2").unwrap(), Numeric::Float(100.0));
        assert_eq!(
            to_numeric_text("-inf").unwrap(),
            Numeric::Float(f64::NEG_INFINITY)
        );
        for text in ["", "?", "nan", "1_0", " ", "١٢"] {
            assert!(to_numeric_text(text).unwrap().as_f64().is_nan(), "{text}");
        }
        assert!(to_numeric_text("99999999999999999999").is_err());
    }

    #[test]
    fn truncate_fraction_uses_the_longest_split() {
        let table = Table {
            names: vec!["Consequence".into(), "Protein_position".into()],
            columns: vec![
                Column::Text(vec![
                    Some("stop_gained".into()),
                    Some("inframe_deletion".into()),
                    None,
                    Some("missense_variant".into()),
                ]),
                Column::Text(vec![
                    Some("25/100".into()),
                    Some("10-14/100".into()),
                    None,
                    Some("?-5/100".into()),
                ]),
            ],
            rows: 4,
        };
        let values = truncate_fraction(&table).unwrap();
        assert_eq!(values[0], 1.0 - 25.0 / 100.0);
        assert_eq!(values[1], 5.0 / 100.0);
        // "nan" has one `/`-field where the longest row has two: no length.
        assert!(values[2].is_nan());
        assert!(values[3].is_nan());
    }
}
