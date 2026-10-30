//! A tab-separated table read and written the way pandas 2.2 does it:
//! `pd.read_table(path, low_memory=False)` and
//! `DataFrame.to_csv(path, sep="\t", index=False)`.
//!
//! Reading follows the C parser. Lines that are empty, or hold only spaces,
//! are skipped; `\n`, `\r\n` and a lone `\r` all end a line; a field that
//! starts with a double quote runs to the next lone quote and may hold tabs,
//! line breaks and `""` for a quote. The first line names the columns (an empty
//! name becomes `Unnamed: <position>`, a repeated name gets `.1`, `.2`, ...).
//! Each column then takes the first type that fits all of its values, tried in
//! pandas' order: 64-bit integer, float, bool, text. A value spelled like one of
//! pandas' missing-value markers (`NA`, `N/A`, `null`, `None`, `nan`, an empty
//! field, ...) is missing; an integer column with a missing value becomes a
//! float column, a bool column with one becomes a column of `True`, `False`
//! and missing. Floats are converted with pandas' `precise_xstrtod` (not
//! always the correctly rounded value, so it is ported step by step).
//!
//! Writing follows `to_csv`: missing values are empty, floats are written as
//! Python's `repr`, bools as `True` and `False`, and a field holding a tab,
//! quote or line break is quoted with `""` for each quote.
//!
//! Not reproduced, and refused with an error rather than guessed: a data row
//! with more fields than the header (pandas may then use the first column as
//! the index), an integer that does not fit 64 bits, a quoted field that is
//! never closed, input that is not UTF-8, and compressions other than gzip.
use crate::cli::{Result, invalid};
use crate::python_text::float_repr;
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;
use std::sync::OnceLock;

/// The spellings pandas reads as missing (`STR_NA_VALUES` of pandas 2.2).
const NA_VALUES: [&str; 19] = [
    "", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN",
    "<NA>", "N/A", "NA", "NULL", "NaN", "None", "n/a", "nan", "null",
];

/// One value of a table, as Python sees it in a row of `DataFrame.iterrows()`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Cell<'a> {
    /// A missing value, which Python holds as the float `nan`.
    Missing,
    Text(&'a str),
    Int(i64),
    Float(f64),
    Bool(bool),
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Column {
    Int(Vec<i64>),
    /// NaN marks a missing value.
    Float(Vec<f64>),
    Bool(Vec<bool>),
    /// pandas' object column of `True`, `False` and missing values.
    BoolOrMissing(Vec<Option<bool>>),
    /// pandas' object column of text; `None` is a missing value.
    Text(Vec<Option<String>>),
}

impl Column {
    pub(crate) fn len(&self) -> usize {
        match self {
            Column::Int(values) => values.len(),
            Column::Float(values) => values.len(),
            Column::Bool(values) => values.len(),
            Column::BoolOrMissing(values) => values.len(),
            Column::Text(values) => values.len(),
        }
    }

    pub(crate) fn cell(&self, row: usize) -> Cell<'_> {
        match self {
            Column::Int(values) => Cell::Int(values[row]),
            Column::Float(values) => {
                let value = values[row];
                if value.is_nan() {
                    Cell::Missing
                } else {
                    Cell::Float(value)
                }
            }
            Column::Bool(values) => Cell::Bool(values[row]),
            Column::BoolOrMissing(values) => match values[row] {
                Some(flag) => Cell::Bool(flag),
                None => Cell::Missing,
            },
            Column::Text(values) => match &values[row] {
                Some(text) => Cell::Text(text.as_str()),
                None => Cell::Missing,
            },
        }
    }
}

#[derive(Debug, PartialEq)]
pub(crate) struct Table {
    pub(crate) names: Vec<String>,
    pub(crate) columns: Vec<Column>,
    pub(crate) rows: usize,
}

impl Table {
    pub(crate) fn column_index(&self, name: &str) -> Option<usize> {
        self.names.iter().position(|candidate| candidate == name)
    }

    pub(crate) fn cell(&self, column: usize, row: usize) -> Cell<'_> {
        self.columns[column].cell(row)
    }

    /// `df[name] = values`: replaces a column in place or appends a new one.
    pub(crate) fn set_column(&mut self, name: &str, column: Column) {
        assert_eq!(column.len(), self.rows, "a new column must fill every row");
        match self.column_index(name) {
            Some(index) => self.columns[index] = column,
            None => {
                self.names.push(name.to_owned());
                self.columns.push(column);
            }
        }
    }
}

fn is_na(token: &str) -> bool {
    NA_VALUES.contains(&token)
}

#[derive(Clone, Copy, PartialEq)]
enum State {
    StartRecord,
    StartField,
    InField,
    InQuotedField,
    QuoteInQuotedField,
}

fn end_field(fields: &mut Vec<String>, field: &mut Vec<u8>) -> Result<()> {
    let text = String::from_utf8(std::mem::take(field))
        .map_err(|_| invalid("the table is not UTF-8 text"))?;
    fields.push(text);
    Ok(())
}

fn end_record(
    fields: &mut Vec<String>,
    field: &mut Vec<u8>,
    emit: &mut impl FnMut(Vec<String>) -> Result<()>,
) -> Result<()> {
    end_field(fields, field)?;
    emit(std::mem::take(fields))
}

/// Splits the text into records of fields, as the pandas C tokenizer does, and
/// hands each record to `emit`.
fn tokenize(text: &str, mut emit: impl FnMut(Vec<String>) -> Result<()>) -> Result<()> {
    let bytes = text.as_bytes();
    let mut state = State::StartRecord;
    let mut fields: Vec<String> = Vec::new();
    let mut field: Vec<u8> = Vec::new();
    let mut index = 0;
    while index < bytes.len() {
        let byte = bytes[index];
        let terminator = byte == b'\n' || byte == b'\r';
        match state {
            State::StartRecord => {
                if terminator {
                    index += 1;
                    continue;
                }
                if byte == b' ' {
                    let mut end = index;
                    while end < bytes.len() && bytes[end] == b' ' {
                        end += 1;
                    }
                    if end == bytes.len() || bytes[end] == b'\n' || bytes[end] == b'\r' {
                        index = end;
                        continue;
                    }
                }
                // The same byte is read again as the first byte of a field.
                state = State::StartField;
                continue;
            }
            State::StartField | State::InField | State::QuoteInQuotedField if terminator => {
                end_record(&mut fields, &mut field, &mut emit)?;
                state = State::StartRecord;
                if byte == b'\r' && bytes.get(index + 1) == Some(&b'\n') {
                    index += 1;
                }
            }
            State::StartField => match byte {
                b'"' => state = State::InQuotedField,
                b'\t' => end_field(&mut fields, &mut field)?,
                _ => {
                    field.push(byte);
                    state = State::InField;
                }
            },
            State::InField => {
                if byte == b'\t' {
                    end_field(&mut fields, &mut field)?;
                    state = State::StartField;
                } else {
                    field.push(byte);
                }
            }
            State::InQuotedField => {
                if byte == b'"' {
                    state = State::QuoteInQuotedField;
                } else {
                    field.push(byte);
                }
            }
            State::QuoteInQuotedField => match byte {
                b'"' => {
                    field.push(b'"');
                    state = State::InQuotedField;
                }
                b'\t' => {
                    end_field(&mut fields, &mut field)?;
                    state = State::StartField;
                }
                _ => {
                    field.push(byte);
                    state = State::InField;
                }
            },
        }
        index += 1;
    }
    match state {
        State::StartRecord => Ok(()),
        State::InQuotedField => Err(invalid("the table ends inside a quoted field")),
        State::StartField | State::InField | State::QuoteInQuotedField => {
            end_record(&mut fields, &mut field, &mut emit)
        }
    }
}

/// The column names pandas gives a header line.
fn column_names(raw: Vec<String>) -> Vec<String> {
    let mut header: Vec<String> = Vec::with_capacity(raw.len());
    let mut unnamed: Vec<usize> = Vec::new();
    for (index, name) in raw.into_iter().enumerate() {
        if name.is_empty() {
            header.push(format!("Unnamed: {index}"));
            unnamed.push(index);
        } else {
            header.push(name);
        }
    }
    // Named columns keep their names; repeated and unnamed ones are mangled after them.
    let order: Vec<usize> = (0..header.len())
        .filter(|index| !unnamed.contains(index))
        .chain(unnamed.iter().copied())
        .collect();
    let mut counts: HashMap<String, usize> = HashMap::new();
    for index in order {
        let old = header[index].clone();
        let mut name = old.clone();
        let mut current = counts.get(&name).copied().unwrap_or(0);
        while current > 0 {
            counts.insert(old.clone(), current + 1);
            name = format!("{old}.{current}");
            if header.contains(&name) {
                current += 1;
            } else {
                current = counts.get(&name).copied().unwrap_or(0);
            }
        }
        header[index] = name.clone();
        counts.insert(name, current + 1);
    }
    header
}

fn is_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

enum IntegerToken {
    Value(i64),
    Overflow,
    Invalid,
}

/// pandas' `str_to_int64`: optional blanks, a sign, digits, optional blanks.
fn integer_token(token: &str) -> IntegerToken {
    let bytes = token.as_bytes();
    let mut position = 0;
    while position < bytes.len() && is_space(bytes[position]) {
        position += 1;
    }
    let negative = position < bytes.len() && bytes[position] == b'-';
    if position < bytes.len() && (bytes[position] == b'-' || bytes[position] == b'+') {
        position += 1;
    }
    if position >= bytes.len() || !bytes[position].is_ascii_digit() {
        return IntegerToken::Invalid;
    }
    let mut magnitude: u128 = 0;
    while position < bytes.len() && bytes[position].is_ascii_digit() {
        magnitude = magnitude * 10 + u128::from(bytes[position] - b'0');
        if magnitude > u128::from(u64::MAX) {
            // Keep reading only to tell an invalid token from an overflow.
            magnitude = u128::from(u64::MAX) + 1;
        }
        position += 1;
    }
    while position < bytes.len() && is_space(bytes[position]) {
        position += 1;
    }
    if position != bytes.len() {
        return IntegerToken::Invalid;
    }
    let limit = if negative {
        1_u128 << 63
    } else {
        (1_u128 << 63) - 1
    };
    if magnitude > limit {
        return IntegerToken::Overflow;
    }
    let value = if negative {
        (magnitude as i128).wrapping_neg() as i64
    } else {
        magnitude as i64
    };
    IntegerToken::Value(value)
}

/// 1e0 to 1e308, each the double nearest to the power of ten, as in pandas' table.
fn power_of_ten(exponent: usize) -> f64 {
    static POWERS: OnceLock<Vec<f64>> = OnceLock::new();
    POWERS.get_or_init(|| {
        (0..=308)
            .map(|power| {
                format!("1e{power}")
                    .parse::<f64>()
                    .expect("a power of ten is a float")
            })
            .collect()
    })[exponent]
}

/// pandas' `precise_xstrtod` (`float_precision=None`): up to 17 digits are
/// accumulated in a double and then scaled once by a power of ten. `None` where
/// pandas reports an error, which makes the column text.
fn precise_xstrtod(token: &str) -> Option<f64> {
    let bytes = token.as_bytes();
    let length = bytes.len();
    let mut position = 0;
    while position < length && is_space(bytes[position]) {
        position += 1;
    }
    let mut negative = false;
    if position < length && bytes[position] == b'-' {
        negative = true;
        position += 1;
    } else if position < length && bytes[position] == b'+' {
        position += 1;
    }
    let max_digits = 17;
    let mut number = 0.0_f64;
    let mut exponent: i64 = 0;
    let mut digits = 0;
    let mut decimals: i64 = 0;
    while position < length && bytes[position].is_ascii_digit() {
        if digits < max_digits {
            number = number * 10.0 + f64::from(bytes[position] - b'0');
            digits += 1;
        } else {
            exponent += 1;
        }
        position += 1;
    }
    if position < length && bytes[position] == b'.' {
        position += 1;
        while digits < max_digits && position < length && bytes[position].is_ascii_digit() {
            number = number * 10.0 + f64::from(bytes[position] - b'0');
            position += 1;
            digits += 1;
            decimals += 1;
        }
        if digits >= max_digits {
            while position < length && bytes[position].is_ascii_digit() {
                position += 1;
            }
        }
        exponent -= decimals;
    }
    if digits == 0 {
        return None;
    }
    if negative {
        number = -number;
    }
    if position < length && (bytes[position] == b'e' || bytes[position] == b'E') {
        position += 1;
        let mut exponent_negative = false;
        if position < length && bytes[position] == b'-' {
            exponent_negative = true;
            position += 1;
        } else if position < length && bytes[position] == b'+' {
            position += 1;
        }
        let mut exponent_digits = 0;
        let mut written: i64 = 0;
        while position < length && bytes[position].is_ascii_digit() {
            written = written
                .saturating_mul(10)
                .saturating_add(i64::from(bytes[position] - b'0'));
            exponent_digits += 1;
            position += 1;
        }
        if exponent_digits == 0 {
            // pandas steps back and then finds the leftover text, so the token fails.
            return None;
        }
        exponent = if exponent_negative {
            exponent.saturating_sub(written)
        } else {
            exponent.saturating_add(written)
        };
    }
    if exponent > 308 {
        return None;
    } else if exponent > 0 {
        number *= power_of_ten(exponent as usize);
    } else if exponent < -308 {
        if exponent < -616 {
            number = 0.0;
        } else {
            number /= power_of_ten((-308 - exponent) as usize);
            number /= power_of_ten(308);
        }
    } else {
        number /= power_of_ten((-exponent) as usize);
    }
    if number.is_infinite() {
        return None;
    }
    while position < length && is_space(bytes[position]) {
        position += 1;
    }
    if position != length {
        return None;
    }
    Some(number)
}

/// A float token as pandas reads it: `precise_xstrtod`, or one of the
/// infinity spellings (any case, no blanks).
fn float_token(token: &str) -> Option<f64> {
    if let Some(value) = precise_xstrtod(token) {
        return Some(value);
    }
    let lower = token.to_ascii_lowercase();
    match lower.as_str() {
        "inf" | "+inf" | "infinity" | "+infinity" => Some(f64::INFINITY),
        "-inf" | "-infinity" => Some(f64::NEG_INFINITY),
        _ => None,
    }
}

/// The column type pandas infers for these tokens.
fn infer_column(tokens: Vec<String>) -> Result<Column> {
    if tokens.is_empty() {
        return Ok(Column::Text(Vec::new()));
    }
    let missing: Vec<bool> = tokens.iter().map(|token| is_na(token)).collect();
    let any_missing = missing.iter().any(|flag| *flag);

    // 1. Integers (every value that is not missing).
    let mut integers: Vec<i64> = Vec::with_capacity(tokens.len());
    let mut all_integers = true;
    for (token, absent) in tokens.iter().zip(&missing) {
        if *absent {
            integers.push(0);
            continue;
        }
        match integer_token(token) {
            IntegerToken::Value(value) => integers.push(value),
            IntegerToken::Overflow => {
                return Err(invalid(format!(
                    "the integer {token:?} does not fit 64 bits, which is outside the supported \
                     range of the table reader"
                )));
            }
            IntegerToken::Invalid => {
                all_integers = false;
                break;
            }
        }
    }
    if all_integers {
        if !any_missing {
            return Ok(Column::Int(integers));
        }
        let floats = integers
            .iter()
            .zip(&missing)
            .map(|(value, absent)| if *absent { f64::NAN } else { *value as f64 })
            .collect();
        return Ok(Column::Float(floats));
    }

    // 2. Floats.
    let mut floats: Vec<f64> = Vec::with_capacity(tokens.len());
    let mut all_floats = true;
    for (token, absent) in tokens.iter().zip(&missing) {
        if *absent {
            floats.push(f64::NAN);
            continue;
        }
        match float_token(token) {
            Some(value) => floats.push(value),
            None => {
                all_floats = false;
                break;
            }
        }
    }
    if all_floats {
        return Ok(Column::Float(floats));
    }

    // 3. Booleans.
    let mut flags: Vec<Option<bool>> = Vec::with_capacity(tokens.len());
    let mut all_flags = true;
    for (token, absent) in tokens.iter().zip(&missing) {
        if *absent {
            flags.push(None);
            continue;
        }
        match token.as_str() {
            "True" | "TRUE" | "true" => flags.push(Some(true)),
            "False" | "FALSE" | "false" => flags.push(Some(false)),
            _ => {
                all_flags = false;
                break;
            }
        }
    }
    if all_flags {
        if any_missing {
            return Ok(Column::BoolOrMissing(flags));
        }
        return Ok(Column::Bool(
            flags.into_iter().map(|flag| flag == Some(true)).collect(),
        ));
    }

    // 4. Text.
    Ok(Column::Text(
        tokens
            .into_iter()
            .zip(missing)
            .map(|(token, absent)| if absent { None } else { Some(token) })
            .collect(),
    ))
}

/// Reads a table from text already decoded as UTF-8.
pub(crate) fn parse_table(text: &str) -> Result<Table> {
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    let mut names: Option<Vec<String>> = None;
    let mut tokens: Vec<Vec<String>> = Vec::new();
    let mut rows = 0;
    tokenize(text, |record| {
        match &names {
            None => {
                let header = column_names(record);
                tokens = vec![Vec::new(); header.len()];
                names = Some(header);
            }
            Some(header) => {
                if record.len() > header.len() {
                    return Err(invalid(format!(
                        "data row {} has {} fields but the header has {}",
                        rows + 1,
                        record.len(),
                        header.len()
                    )));
                }
                let width = record.len();
                for (column, field) in tokens.iter_mut().zip(record) {
                    column.push(field);
                }
                for column in tokens.iter_mut().skip(width) {
                    column.push(String::new());
                }
                rows += 1;
            }
        }
        Ok(())
    })?;
    let names = names.ok_or_else(|| invalid("no columns to parse from the table"))?;
    let mut columns = Vec::with_capacity(names.len());
    for column in tokens {
        columns.push(infer_column(column)?);
    }
    Ok(Table {
        names,
        columns,
        rows,
    })
}

/// `pd.read_table(path, low_memory=False)`.
pub(crate) fn read_table(path: &Path) -> Result<Table> {
    let name = path.to_string_lossy().to_lowercase();
    let mut bytes = Vec::new();
    let mut file = File::open(path)?;
    if name.ends_with(".gz") {
        MultiGzDecoder::new(file).read_to_end(&mut bytes)?;
    } else if [".bz2", ".xz", ".zip", ".zst", ".tar"]
        .iter()
        .any(|extension| name.ends_with(*extension))
    {
        return Err(invalid(format!(
            "{}: pandas would decompress this file by its name, which this reader does not do",
            path.display()
        )));
    } else {
        file.read_to_end(&mut bytes)?;
    }
    let text = String::from_utf8(bytes)
        .map_err(|_| invalid(format!("{}: not UTF-8 text", path.display())))?;
    parse_table(&text)
}

fn write_field(out: &mut impl Write, text: &str, only_column: bool) -> std::io::Result<()> {
    if text.is_empty() && only_column {
        // csv.writer quotes a lone empty field so that the row is not a blank line.
        return out.write_all(b"\"\"");
    }
    if text.contains(['\t', '"', '\n', '\r']) {
        out.write_all(b"\"")?;
        out.write_all(text.replace('"', "\"\"").as_bytes())?;
        out.write_all(b"\"")
    } else {
        out.write_all(text.as_bytes())
    }
}

/// `df.to_csv(out, sep="\t", index=False)`.
pub(crate) fn write_table_to(out: &mut impl Write, table: &Table) -> std::io::Result<()> {
    let only_column = table.names.len() == 1;
    for (index, name) in table.names.iter().enumerate() {
        if index > 0 {
            out.write_all(b"\t")?;
        }
        write_field(out, name, only_column)?;
    }
    out.write_all(b"\n")?;
    for row in 0..table.rows {
        for (index, column) in table.columns.iter().enumerate() {
            if index > 0 {
                out.write_all(b"\t")?;
            }
            match column {
                Column::Int(values) => write!(out, "{}", values[row])?,
                Column::Float(values) => {
                    if !values[row].is_nan() {
                        out.write_all(float_repr(values[row]).as_bytes())?;
                    }
                }
                Column::Bool(values) => {
                    let text: &[u8] = if values[row] { b"True" } else { b"False" };
                    out.write_all(text)?;
                }
                Column::BoolOrMissing(values) => match values[row] {
                    Some(true) => out.write_all(b"True")?,
                    Some(false) => out.write_all(b"False")?,
                    None => {}
                },
                Column::Text(values) => {
                    if let Some(text) = &values[row] {
                        write_field(out, text, only_column)?;
                    }
                }
            }
        }
        out.write_all(b"\n")?;
    }
    Ok(())
}

/// Writes the table to a new file at `path`.
pub(crate) fn write_table(path: &Path, table: &Table) -> Result<()> {
    let mut out = BufWriter::new(File::create(path)?);
    write_table_to(&mut out, table)?;
    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn round_trip(text: &str) -> String {
        let table = parse_table(text).unwrap();
        let mut out = Vec::new();
        write_table_to(&mut out, &table).unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn column_types_follow_pandas_order() {
        let table = parse_table(
            "i\tj\tf\tb\tc\tw\n1\t1\t1.5\tTrue\tx\tTRUE\n2\t\t2\tfalse\ty\t\n3\t4\t-0.5\tFALSE\tNA\ttrue\n",
        )
        .unwrap();
        assert_eq!(table.rows, 3);
        assert_eq!(table.columns[0], Column::Int(vec![1, 2, 3]));
        // An integer column with a missing value is a float column.
        assert_eq!(table.columns[1].cell(0), Cell::Float(1.0));
        assert_eq!(table.columns[1].cell(1), Cell::Missing);
        assert_eq!(table.columns[2], Column::Float(vec![1.5, 2.0, -0.5]));
        assert_eq!(table.columns[3], Column::Bool(vec![true, false, false]));
        assert_eq!(
            table.columns[4],
            Column::Text(vec![Some("x".into()), Some("y".into()), None])
        );
        assert_eq!(
            table.columns[5],
            Column::BoolOrMissing(vec![Some(true), None, Some(true)])
        );
    }

    #[test]
    fn missing_spellings_follow_pandas() {
        let missing = [
            "", "#N/A", "N/A", "NA", "NULL", "NaN", "None", "n/a", "nan", "null", "<NA>", "-nan",
        ];
        let kept = ["na", "none", "NAN", ".", "-", "Null", "NONE", "nA"];
        let mut text = String::from("a\tb\n");
        for value in missing.iter().chain(kept.iter()) {
            text.push_str(&format!("{value}\tz\n"));
        }
        let table = parse_table(&text).unwrap();
        let Column::Text(values) = &table.columns[0] else {
            panic!("a column of words is text")
        };
        assert!(values[..missing.len()].iter().all(Option::is_none));
        for (value, expected) in values[missing.len()..].iter().zip(kept) {
            assert_eq!(value.as_deref(), Some(expected));
        }
    }

    #[test]
    fn numbers_are_written_as_python_prints_them() {
        let text = "x\ty\n0.10\t007\n1e-05\t+5\n5.\t 7 \n.5\t-0\n1E3\t12\n";
        // y is an integer column, x a float column.
        assert_eq!(
            round_trip(text),
            "x\ty\n0.1\t7\n1e-05\t5\n5.0\t7\n0.5\t0\n1000.0\t12\n"
        );
    }

    #[test]
    fn float_tokens_follow_pandas_rules() {
        assert_eq!(float_token("1.5"), Some(1.5));
        assert_eq!(float_token(" 2e2 "), Some(200.0));
        assert_eq!(float_token("1e999"), None);
        assert_eq!(float_token("1e"), None);
        assert_eq!(float_token("1,5"), None);
        assert_eq!(float_token("."), None);
        assert_eq!(float_token("1_0"), None);
        assert_eq!(float_token("NAN"), None);
        assert_eq!(float_token("-Infinity"), Some(f64::NEG_INFINITY));
        assert_eq!(float_token("INF"), Some(f64::INFINITY));
        assert_eq!(float_token(" inf"), None);
        // pandas counts a leading zero among the 17 digits it keeps, so the
        // last digit of this value is dropped.
        assert_eq!(float_token("0.30000000000000004"), Some(0.3));
        assert_eq!(float_token("0.1234567"), Some(0.1234567));
        assert!(float_token("1e-320").is_some());
    }

    #[test]
    fn quoted_fields_blank_lines_and_line_ends() {
        let text = "a\tb\r\n\r\n   \r\n\"x\ty\"\t\"q\"\"r\"\r\n\"two\nlines\"\tz\r\nlast\tw";
        assert_eq!(
            round_trip(text),
            "a\tb\n\"x\ty\"\t\"q\"\"r\"\n\"two\nlines\"\tz\nlast\tw\n"
        );
        // A lone carriage return ends a line; a short row is padded with missing values.
        assert_eq!(round_trip("a\tb\rx\ty\rz\r"), "a\tb\nx\ty\nz\t\n");
        // A quote inside an unquoted field is a plain character.
        assert_eq!(round_trip("a\nx\"y\n"), "a\n\"x\"\"y\"\n");
    }

    #[test]
    fn header_names_are_mangled_like_pandas() {
        let table = parse_table("a\ta\t\tb\ta\n1\t2\t3\t4\t5\n").unwrap();
        assert_eq!(table.names, ["a", "a.1", "Unnamed: 2", "b", "a.2"]);
    }

    #[test]
    fn unsupported_shapes_are_refused() {
        assert!(parse_table("a\tb\n1\t2\t3\n").is_err());
        assert!(parse_table("a\n\"open\n").is_err());
        assert!(parse_table("").is_err());
        assert!(parse_table("a\n99999999999999999999\n").is_err());
        // A header alone is a table without rows.
        assert_eq!(parse_table("a\tb\n").unwrap().rows, 0);
    }

    #[test]
    fn a_lone_empty_field_is_quoted_and_set_column_replaces_in_place() {
        let mut table = parse_table("a\n1\n").unwrap();
        table.set_column("b", Column::Text(vec![Some(String::new())]));
        table.set_column("a", Column::Bool(vec![true]));
        let mut out = Vec::new();
        write_table_to(&mut out, &table).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), "a\tb\nTrue\t\n");
        let mut single = parse_table("a\n1\n").unwrap();
        single.set_column("a", Column::Text(vec![Some(String::new())]));
        let mut out = Vec::new();
        write_table_to(&mut out, &single).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), "a\n\"\"\n");
    }
}
