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
//! the index), a column in which the first value that does not parse as a
//! 64-bit integer starts with more digits than a signed 64-bit integer holds
//! (pandas then reads the column as unsigned integers, or keeps the text of
//! every value), a quoted field that is never closed, input that is not
//! UTF-8, and compressions other than gzip.
//!
//! The ACMG step also reads with these `read_table` options (`ReadOptions`):
//! `dtype=str` (every value stays text), `keep_default_na=False` (no value is
//! missing; the empty fields that pad a short row stay empty text), `header=None`
//! with `names=[...]`, `usecols` by position or by name (the file's column order
//! is kept, as pandas keeps it) and `comment="#"` (the rest of a line after an
//! unquoted `#` is dropped; a line that is only a comment is skipped).
//!
//! Besides the column types the reader produces, a column can hold what the
//! step's pandas code makes later: an `object` column that mixes text, numbers,
//! booleans and missing values (`Column::Object`), and the nullable `Int64` and
//! `boolean` columns (`Column::NullableInt`, `Column::NullableBool`). The
//! module `pandas_polars` moves typed columns to polars and back.
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

/// One value of a pandas `object` column that mixes kinds of values.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum ObjectValue {
    /// The float `nan`, which is how pandas marks a missing value it creates.
    Nan,
    /// Python's `None`. It prints and counts as missing like `nan`, but some
    /// pandas calls tell the two apart (`isin([None])` does not match `nan`).
    None,
    Text(String),
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
    /// pandas' object column of any other mix of values.
    Object(Vec<ObjectValue>),
    /// pandas' nullable `Int64` column; `None` is `pd.NA`.
    NullableInt(Vec<Option<i64>>),
    /// pandas' nullable `boolean` column; `None` is `pd.NA`.
    NullableBool(Vec<Option<bool>>),
}

impl Column {
    pub(crate) fn len(&self) -> usize {
        match self {
            Column::Int(values) => values.len(),
            Column::Float(values) => values.len(),
            Column::Bool(values) => values.len(),
            Column::BoolOrMissing(values) => values.len(),
            Column::Text(values) => values.len(),
            Column::Object(values) => values.len(),
            Column::NullableInt(values) => values.len(),
            Column::NullableBool(values) => values.len(),
        }
    }

    /// The value at `row`. `nan`, `None` and `pd.NA` are all `Cell::Missing`;
    /// code that must tell them apart looks at the column itself.
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
            Column::BoolOrMissing(values) | Column::NullableBool(values) => match values[row] {
                Some(flag) => Cell::Bool(flag),
                None => Cell::Missing,
            },
            Column::Text(values) => match &values[row] {
                Some(text) => Cell::Text(text.as_str()),
                None => Cell::Missing,
            },
            Column::Object(values) => match &values[row] {
                ObjectValue::Nan | ObjectValue::None => Cell::Missing,
                ObjectValue::Text(text) => Cell::Text(text.as_str()),
                ObjectValue::Int(number) => Cell::Int(*number),
                ObjectValue::Float(number) if number.is_nan() => Cell::Missing,
                ObjectValue::Float(number) => Cell::Float(*number),
                ObjectValue::Bool(flag) => Cell::Bool(*flag),
            },
            Column::NullableInt(values) => match values[row] {
                Some(number) => Cell::Int(number),
                None => Cell::Missing,
            },
        }
    }

    /// The pandas dtype name (`str(series.dtype)`).
    pub(crate) fn dtype_name(&self) -> &'static str {
        match self {
            Column::Int(_) => "int64",
            Column::Float(_) => "float64",
            Column::Bool(_) => "bool",
            Column::BoolOrMissing(_) | Column::Text(_) | Column::Object(_) => "object",
            Column::NullableInt(_) => "Int64",
            Column::NullableBool(_) => "boolean",
        }
    }

    /// `series.isna()`: NaN, `None` and `pd.NA` are missing.
    pub(crate) fn missing_mask(&self) -> Vec<bool> {
        (0..self.len())
            .map(|row| self.cell(row) == Cell::Missing)
            .collect()
    }

    /// The rows at `rows`, in that order (`series.iloc[rows]`).
    pub(crate) fn take(&self, rows: &[usize]) -> Column {
        fn pick<T: Clone>(values: &[T], rows: &[usize]) -> Vec<T> {
            rows.iter().map(|row| values[*row].clone()).collect()
        }
        match self {
            Column::Int(values) => Column::Int(pick(values, rows)),
            Column::Float(values) => Column::Float(pick(values, rows)),
            Column::Bool(values) => Column::Bool(pick(values, rows)),
            Column::BoolOrMissing(values) => Column::BoolOrMissing(pick(values, rows)),
            Column::Text(values) => Column::Text(pick(values, rows)),
            Column::Object(values) => Column::Object(pick(values, rows)),
            Column::NullableInt(values) => Column::NullableInt(pick(values, rows)),
            Column::NullableBool(values) => Column::NullableBool(pick(values, rows)),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
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

    /// `df[name]`: the first column of that name, or the error pandas raises
    /// (`KeyError: 'name'`), so a missing column stops the step as it does in Python.
    pub(crate) fn column(&self, name: &str) -> Result<&Column> {
        match self.column_index(name) {
            Some(index) => Ok(&self.columns[index]),
            None => Err(invalid(format!("KeyError: '{name}'"))),
        }
    }

    /// `name in df.columns`.
    #[allow(dead_code)]
    pub(crate) fn has_column(&self, name: &str) -> bool {
        self.column_index(name).is_some()
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

    /// The rows at `rows`, in that order, for a filter (`df[mask]`), a sort or
    /// `drop_duplicates`. Row labels are not kept; the step never writes them.
    pub(crate) fn take_rows(&self, rows: &[usize]) -> Table {
        Table {
            names: self.names.clone(),
            columns: self
                .columns
                .iter()
                .map(|column| column.take(rows))
                .collect(),
            rows: rows.len(),
        }
    }

    /// `df[mask]` for a mask without missing values.
    #[allow(dead_code)]
    pub(crate) fn filter_rows(&self, mask: &[bool]) -> Table {
        assert_eq!(mask.len(), self.rows, "a row mask must cover every row");
        let rows: Vec<usize> = (0..self.rows).filter(|row| mask[*row]).collect();
        self.take_rows(&rows)
    }
}

fn is_na(token: &str) -> bool {
    NA_VALUES.contains(&token)
}

/// Where the column names come from (`header=0` or `header=None, names=[...]`).
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) enum Header {
    /// The first line names the columns (pandas' default).
    #[default]
    FirstLine,
    /// `header=None, names=[...]`: every line is data.
    Names(Vec<String>),
}

/// `usecols`: which columns of the file to keep.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum UseColumns {
    /// `usecols=range(n)` or a list of positions in the file.
    Positions(Vec<usize>),
    /// `usecols=[name, ...]`; the columns keep the file's order.
    Names(Vec<String>),
}

/// The `pd.read_table` options the ACMG step uses; the default is
/// `pd.read_table(path, low_memory=False)`.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ReadOptions {
    /// `dtype=str`: no column type inference, every value stays text.
    pub(crate) all_text: bool,
    /// `keep_default_na`: whether pandas' missing-value spellings are missing.
    pub(crate) keep_default_na: bool,
    pub(crate) header: Header,
    pub(crate) use_columns: Option<UseColumns>,
    /// `comment="#"`: the byte that starts a comment.
    pub(crate) comment: Option<u8>,
}

impl Default for ReadOptions {
    fn default() -> Self {
        ReadOptions {
            all_text: false,
            keep_default_na: true,
            header: Header::FirstLine,
            use_columns: None,
            comment: None,
        }
    }
}

#[derive(Clone, Copy, PartialEq)]
enum State {
    StartRecord,
    StartField,
    InField,
    InQuotedField,
    QuoteInQuotedField,
    /// After a comment byte: the rest of the line is dropped.
    InComment,
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
/// hands each record to `emit`. With a `comment` byte, an unquoted comment byte
/// ends the record (a record that has not started is skipped) and the rest of
/// its line is dropped.
fn tokenize(
    text: &str,
    comment: Option<u8>,
    mut emit: impl FnMut(Vec<String>) -> Result<()>,
) -> Result<()> {
    let bytes = text.as_bytes();
    let mut state = State::StartRecord;
    let mut fields: Vec<String> = Vec::new();
    let mut field: Vec<u8> = Vec::new();
    let mut index = 0;
    while index < bytes.len() {
        let byte = bytes[index];
        let terminator = byte == b'\n' || byte == b'\r';
        if comment == Some(byte)
            && matches!(
                state,
                State::StartRecord | State::StartField | State::InField
            )
        {
            if state != State::StartRecord {
                end_record(&mut fields, &mut field, &mut emit)?;
            }
            state = State::InComment;
            index += 1;
            continue;
        }
        match state {
            State::InComment => {
                if terminator {
                    state = State::StartRecord;
                }
            }
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
        State::StartRecord | State::InComment => Ok(()),
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
/// Like the C code it reports an overflow as soon as the digits run past the
/// 64-bit range, before it looks at the characters that follow them.
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
    let limit = if negative {
        1_u128 << 63
    } else {
        (1_u128 << 63) - 1
    };
    let mut magnitude: u128 = 0;
    while position < bytes.len() && bytes[position].is_ascii_digit() {
        magnitude = magnitude * 10 + u128::from(bytes[position] - b'0');
        if magnitude > limit {
            return IntegerToken::Overflow;
        }
        position += 1;
    }
    while position < bytes.len() && is_space(bytes[position]) {
        position += 1;
    }
    if position != bytes.len() {
        return IntegerToken::Invalid;
    }
    let value = if negative {
        (magnitude as i128).wrapping_neg() as i64
    } else {
        magnitude as i64
    };
    IntegerToken::Value(value)
}

enum UnsignedToken {
    /// A sign of `-`: pandas does not look at the rest of the token.
    Negative,
    Value,
    Overflow,
    Invalid,
}

/// pandas' `str_to_uint64`, which it tries on a column after a value ran past
/// the signed 64-bit range.
fn unsigned_token(token: &str) -> UnsignedToken {
    let bytes = token.as_bytes();
    let mut position = 0;
    while position < bytes.len() && is_space(bytes[position]) {
        position += 1;
    }
    if position < bytes.len() && bytes[position] == b'-' {
        return UnsignedToken::Negative;
    }
    if position < bytes.len() && bytes[position] == b'+' {
        position += 1;
    }
    if position >= bytes.len() || !bytes[position].is_ascii_digit() {
        return UnsignedToken::Invalid;
    }
    let mut magnitude: u128 = 0;
    while position < bytes.len() && bytes[position].is_ascii_digit() {
        magnitude = magnitude * 10 + u128::from(bytes[position] - b'0');
        if magnitude > u128::from(u64::MAX) {
            return UnsignedToken::Overflow;
        }
        position += 1;
    }
    while position < bytes.len() && is_space(bytes[position]) {
        position += 1;
    }
    if position != bytes.len() {
        return UnsignedToken::Invalid;
    }
    UnsignedToken::Value
}

/// The first part of a token, for a message.
fn shortened(token: &str) -> String {
    let mut text: String = token.chars().take(40).collect();
    if text.len() < token.len() {
        text.push_str("...");
    }
    format!("{text:?}")
}

/// What pandas does when the first value of a column that is not a 64-bit
/// integer runs past the signed 64-bit range: it tries the whole column as
/// unsigned integers. The column ends up unsigned integers or text (in which
/// case every value keeps its original spelling) unless that first unsigned
/// failure is not an overflow, and then pandas goes on to floats as usual.
/// Only that last outcome is reproduced; the others are refused.
fn refuse_unless_float_follows(
    tokens: &[String],
    missing: &[bool],
    overflowing: &str,
) -> Result<()> {
    for (token, absent) in tokens.iter().zip(missing) {
        if *absent {
            continue;
        }
        match unsigned_token(token) {
            UnsignedToken::Negative | UnsignedToken::Value => {}
            UnsignedToken::Overflow => break,
            UnsignedToken::Invalid => return Ok(()),
        }
    }
    Err(invalid(format!(
        "the value {} lies outside the signed 64-bit integer range, where pandas reads a column \
         as unsigned integers or keeps its text; the table reader does not reproduce that",
        shortened(overflowing)
    )))
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
///
/// Two details of the C code matter for odd exponents. It reads at most 17
/// digits of an exponent (a longer one leaves digits behind, so the token
/// fails), and it adds them up in a 32-bit `int`, so a very long exponent wraps
/// around: `1e999999999999` is read as 0.0 and `1e4294967296` as 1.0.
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
    let mut exponent: i32 = 0;
    let mut digits = 0;
    let mut decimals: i32 = 0;
    while position < length && bytes[position].is_ascii_digit() {
        if digits < max_digits {
            number = number * 10.0 + f64::from(bytes[position] - b'0');
            digits += 1;
        } else {
            exponent = exponent.wrapping_add(1);
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
        exponent = exponent.wrapping_sub(decimals);
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
        let mut written: i32 = 0;
        while exponent_digits < max_digits && position < length && bytes[position].is_ascii_digit()
        {
            written = written
                .wrapping_mul(10)
                .wrapping_add(i32::from(bytes[position] - b'0'));
            exponent_digits += 1;
            position += 1;
        }
        if exponent_digits == 0 {
            // pandas steps back and then finds the leftover text, so the token fails.
            return None;
        }
        exponent = if exponent_negative {
            exponent.wrapping_sub(written)
        } else {
            exponent.wrapping_add(written)
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

/// The column type pandas infers for these tokens; `missing_spellings` says
/// whether pandas' missing-value spellings count as missing (`keep_default_na`).
fn infer_column(tokens: Vec<String>, missing_spellings: bool) -> Result<Column> {
    if tokens.is_empty() {
        return Ok(Column::Text(Vec::new()));
    }
    let missing: Vec<bool> = tokens
        .iter()
        .map(|token| missing_spellings && is_na(token))
        .collect();
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
                refuse_unless_float_follows(&tokens, &missing, token)?;
                all_integers = false;
                break;
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
        // pandas marks a missing integer with the smallest int64, so that value is lost too.
        let floats = integers
            .iter()
            .zip(&missing)
            .map(|(value, absent)| {
                if *absent || *value == i64::MIN {
                    f64::NAN
                } else {
                    *value as f64
                }
            })
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
        // pandas compares with `strcasecmp`, so every capitalisation counts.
        if token.eq_ignore_ascii_case("true") {
            flags.push(Some(true));
        } else if token.eq_ignore_ascii_case("false") {
            flags.push(Some(false));
        } else {
            all_flags = false;
            break;
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
#[allow(dead_code)]
pub(crate) fn parse_table(text: &str) -> Result<Table> {
    parse_table_with(text, &ReadOptions::default())
}

/// Reads a table from text already decoded as UTF-8, with `pd.read_table` options.
pub(crate) fn parse_table_with(text: &str, options: &ReadOptions) -> Result<Table> {
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    let mut names: Option<Vec<String>> = match &options.header {
        Header::FirstLine => None,
        Header::Names(given) => Some(given.clone()),
    };
    // With `usecols` a row may be longer than the columns named; without it
    // pandas would turn the extra leading fields into the row index.
    let mut width: Option<usize> = match (&options.header, &options.use_columns) {
        (Header::Names(_), Some(UseColumns::Positions(positions))) => {
            Some(positions.iter().max().map_or(0, |last| last + 1))
        }
        (Header::Names(given), _) => Some(given.len()),
        (Header::FirstLine, _) => None,
    };
    let mut tokens: Vec<Vec<String>> = vec![Vec::new(); width.unwrap_or(0)];
    let mut rows = 0;
    tokenize(text, options.comment, |record| {
        let Some(expected) = width else {
            let header = column_names(record);
            tokens = vec![Vec::new(); header.len()];
            width = Some(header.len());
            names = Some(header);
            return Ok(());
        };
        let longest = match options.use_columns {
            Some(UseColumns::Positions(_)) if matches!(options.header, Header::Names(_)) => {
                usize::MAX
            }
            _ => expected,
        };
        if record.len() > longest {
            return Err(invalid(format!(
                "data row {} has {} fields but the header has {}",
                rows + 1,
                record.len(),
                expected
            )));
        }
        let filled = record.len().min(expected);
        for (column, field) in tokens.iter_mut().zip(record) {
            column.push(field);
        }
        for column in tokens.iter_mut().skip(filled) {
            column.push(String::new());
        }
        rows += 1;
        Ok(())
    })?;
    let mut names = names.ok_or_else(|| invalid("no columns to parse from the table"))?;
    if let (Header::Names(given), Some(UseColumns::Positions(positions))) =
        (&options.header, &options.use_columns)
    {
        let mut kept = positions.clone();
        kept.sort_unstable();
        kept.dedup();
        if kept.len() != given.len() {
            return Err(invalid(
                "names and usecols must name the same number of columns",
            ));
        }
        // The names go to the kept columns in the file's order.
        let mut kept_tokens = Vec::with_capacity(kept.len());
        for position in &kept {
            kept_tokens.push(std::mem::take(&mut tokens[*position]));
        }
        tokens = kept_tokens;
        names = given.clone();
    } else if let Some(use_columns) = &options.use_columns {
        let kept: Vec<usize> = match use_columns {
            UseColumns::Positions(positions) => {
                if let Some(outside) = positions.iter().find(|position| **position >= names.len()) {
                    return Err(invalid(format!(
                        "usecols position {outside} is outside the {} columns of the table",
                        names.len()
                    )));
                }
                (0..names.len())
                    .filter(|index| positions.contains(index))
                    .collect()
            }
            UseColumns::Names(wanted) => {
                let absent: Vec<&String> =
                    wanted.iter().filter(|name| !names.contains(name)).collect();
                if !absent.is_empty() {
                    return Err(invalid(format!(
                        "ValueError: Usecols do not match columns, columns expected but not found: {absent:?}"
                    )));
                }
                (0..names.len())
                    .filter(|index| wanted.contains(&names[*index]))
                    .collect()
            }
        };
        let mut kept_tokens = Vec::with_capacity(kept.len());
        let mut kept_names = Vec::with_capacity(kept.len());
        for index in kept {
            kept_tokens.push(std::mem::take(&mut tokens[index]));
            kept_names.push(names[index].clone());
        }
        tokens = kept_tokens;
        names = kept_names;
    }
    let mut columns = Vec::with_capacity(names.len());
    for column in tokens {
        columns.push(if options.all_text {
            Column::Text(
                column
                    .into_iter()
                    .map(|token| {
                        if options.keep_default_na && is_na(&token) {
                            None
                        } else {
                            Some(token)
                        }
                    })
                    .collect(),
            )
        } else {
            infer_column(column, options.keep_default_na)?
        });
    }
    Ok(Table {
        names,
        columns,
        rows,
    })
}

/// `pd.read_table(path, low_memory=False)`.
pub(crate) fn read_table(path: &Path) -> Result<Table> {
    read_table_with(path, &ReadOptions::default())
}

/// `pd.read_table(path, low_memory=False, ...)` with the options given.
pub(crate) fn read_table_with(path: &Path, options: &ReadOptions) -> Result<Table> {
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
    parse_table_with(&text, options)
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
        if only_column && table.columns[0].cell(row) == Cell::Missing {
            // A missing value is written as an empty field, which csv.writer
            // quotes when it is the whole row.
            out.write_all(b"\"\"\n")?;
            continue;
        }
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
                Column::Object(values) => match &values[row] {
                    ObjectValue::Nan | ObjectValue::None => {}
                    ObjectValue::Text(text) => write_field(out, text, only_column)?,
                    ObjectValue::Int(number) => write!(out, "{number}")?,
                    // `str()` of a Python float; NaN is written as missing.
                    ObjectValue::Float(number) => {
                        if !number.is_nan() {
                            out.write_all(float_repr(*number).as_bytes())?;
                        }
                    }
                    ObjectValue::Bool(true) => out.write_all(b"True")?,
                    ObjectValue::Bool(false) => out.write_all(b"False")?,
                },
                Column::NullableInt(values) => {
                    if let Some(number) = values[row] {
                        write!(out, "{number}")?;
                    }
                }
                Column::NullableBool(values) => match values[row] {
                    Some(true) => out.write_all(b"True")?,
                    Some(false) => out.write_all(b"False")?,
                    None => {}
                },
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
    fn exponents_follow_the_c_parser() {
        // At most 17 exponent digits are read, leading zeros included.
        assert_eq!(float_token(&format!("1e{}5", "0".repeat(16))), Some(1e5));
        assert_eq!(float_token(&format!("1e{}5", "0".repeat(17))), None);
        // The digits are added up in a 32-bit integer that wraps around.
        assert_eq!(float_token("1e999999999999"), Some(0.0));
        assert_eq!(float_token("1e-999999999999"), None);
        assert_eq!(float_token("1e4294967296"), Some(1.0));
        assert_eq!(float_token("1e4294967297"), Some(10.0));
        assert_eq!(float_token("1e-4294967295"), Some(10.0));
        assert_eq!(float_token("1e2147483647"), None);
        assert_eq!(float_token("1e2147483648"), Some(0.0));
        assert_eq!(float_token("1e-2147483649"), None);
        assert_eq!(float_token("1e11111111111111111"), None);
        assert_eq!(float_token("0e999999"), None);
        assert_eq!(float_token("0e-999999"), Some(0.0));
    }

    fn single_column(values: &[&str]) -> Result<Column> {
        let table = parse_table(&format!("x\n{}\n", values.join("\n")))?;
        Ok(table.columns.into_iter().next().unwrap())
    }

    #[test]
    fn numbers_past_the_signed_range_are_floats_or_refused() {
        let wide = "1".repeat(21);
        let float_text = |values: &[&str]| match single_column(values).unwrap() {
            Column::Float(numbers) => numbers.iter().map(|n| float_repr(*n)).collect::<Vec<_>>(),
            other => panic!("not a float column: {other:?}"),
        };
        // A float with 19 or 20 digits still fits the unsigned range, so pandas reads a float.
        assert_eq!(
            float_text(&["9223372036854775808.5"]),
            ["9.223372036854778e+18"]
        );
        assert_eq!(
            float_text(&["+11111111111111111111.5"]),
            ["1.1111111111111113e+19"]
        );
        assert_eq!(
            float_text(&["-9223372036854775808.5"]),
            ["-9.223372036854778e+18"]
        );
        // After a value that is no integer, a long digit string is only a float.
        assert_eq!(
            float_text(&["1.5", &format!("{}.5", "1".repeat(30))]),
            ["1.5", "1.1111111111111112e+29"]
        );
        assert_eq!(
            float_text(&["1.5", &"1".repeat(30)]),
            ["1.5", "1.1111111111111112e+29"]
        );
        assert_eq!(
            single_column(&["abc", &"1".repeat(30)]).unwrap(),
            Column::Text(vec![Some("abc".into()), Some("1".repeat(30))])
        );
        // Everything else past the signed range is read by pandas as unsigned integers or text.
        for values in [
            vec![format!("{wide}.5")],
            vec![format!("-{}.5", "1".repeat(20))],
            vec!["-9223372036854775809.5".to_owned()],
            vec!["9223372036854775808".to_owned()],
            vec!["18446744073709551616.5".to_owned()],
            vec!["1".repeat(30), "1.5".to_owned()],
            vec!["1".to_owned(), "1".repeat(30)],
            vec![format!("{}.5", "1".repeat(30)), "abc".to_owned()],
        ] {
            let refs: Vec<&str> = values.iter().map(String::as_str).collect();
            let error = single_column(&refs).unwrap_err().to_string();
            assert!(error.contains("64-bit"), "{values:?}: {error}");
        }
    }

    #[test]
    fn booleans_are_matched_without_regard_to_case() {
        assert_eq!(
            single_column(&["tRuE", "fAlSe", "TRUE", "false"]).unwrap(),
            Column::Bool(vec![true, false, true, false])
        );
        assert!(matches!(
            single_column(&["True ", "False"]).unwrap(),
            Column::Text(_)
        ));
        assert!(matches!(
            single_column(&["T", "F"]).unwrap(),
            Column::Text(_)
        ));
    }

    #[test]
    fn floats_that_are_exact_ties_print_with_the_even_digit() {
        assert_eq!(
            round_trip("x\n129.14285278320312\n100.10678100585938\n111111111111111.12\n"),
            "x\n129.14285278320312\n100.10678100585938\n111111111111111.12\n"
        );
    }

    #[test]
    fn the_smallest_int64_is_the_missing_marker_of_an_integer_column_with_gaps() {
        assert_eq!(
            single_column(&["-9223372036854775808", "7"]).unwrap(),
            Column::Int(vec![i64::MIN, 7])
        );
        let Column::Float(numbers) = single_column(&["-9223372036854775808", "7", "NA"]).unwrap()
        else {
            panic!()
        };
        assert!(numbers[0].is_nan() && numbers[2].is_nan());
        assert_eq!(numbers[1], 7.0);
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

    fn texts(column: &Column) -> Vec<Option<&str>> {
        let Column::Text(values) = column else {
            panic!("not a text column: {column:?}")
        };
        values.iter().map(|value| value.as_deref()).collect()
    }

    // The expected values below are what pandas 2.2.3 returned for the same
    // text (group 3 probe, `probe_pandas_semantics.py`).
    #[test]
    fn short_rows_follow_the_missing_value_options() {
        let text = "a\tb\tc\n1\t\n2\tNA\tx\n";
        let all_text = ReadOptions {
            all_text: true,
            keep_default_na: false,
            ..ReadOptions::default()
        };
        let table = parse_table_with(text, &all_text).unwrap();
        assert_eq!(texts(&table.columns[0]), [Some("1"), Some("2")]);
        assert_eq!(texts(&table.columns[1]), [Some(""), Some("NA")]);
        assert_eq!(texts(&table.columns[2]), [Some(""), Some("x")]);
        let with_na = ReadOptions {
            all_text: true,
            ..ReadOptions::default()
        };
        let table = parse_table_with(text, &with_na).unwrap();
        assert_eq!(texts(&table.columns[1]), [None, None]);
        assert_eq!(texts(&table.columns[2]), [None, Some("x")]);
        let typed = ReadOptions {
            keep_default_na: false,
            ..ReadOptions::default()
        };
        let table = parse_table_with(text, &typed).unwrap();
        assert_eq!(table.columns[0], Column::Int(vec![1, 2]));
        assert_eq!(texts(&table.columns[1]), [Some(""), Some("NA")]);
        let numbers = parse_table_with("x\n007\n1.50\n", &all_text).unwrap();
        assert_eq!(texts(&numbers.columns[0]), [Some("007"), Some("1.50")]);
    }

    #[test]
    fn headerless_pedigree_keeps_six_columns() {
        let names: Vec<String> = [
            "#FamilyID",
            "IndividualID",
            "PaternalID",
            "MaternalID",
            "Sex",
            "Phenotype",
        ]
        .map(str::to_owned)
        .to_vec();
        let options = ReadOptions {
            all_text: true,
            keep_default_na: false,
            header: Header::Names(names.clone()),
            use_columns: Some(UseColumns::Positions((0..6).collect())),
            comment: None,
        };
        let table = parse_table_with("F\tI\t0\t0\t1\t2\textra\nF\tJ\t0\t0\t2\n", &options).unwrap();
        assert_eq!(table.names, names);
        assert_eq!(texts(&table.columns[1]), [Some("I"), Some("J")]);
        assert_eq!(texts(&table.columns[5]), [Some("2"), Some("")]);
    }

    #[test]
    fn comments_names_and_usecols_follow_pandas() {
        let gene_list = ReadOptions {
            header: Header::Names(vec!["SYMBOL".into()]),
            comment: Some(b'#'),
            ..ReadOptions::default()
        };
        let table =
            parse_table_with("#comment\nGENE1\nGENE2 # tail\n\n#x\nGENE3\n", &gene_list).unwrap();
        assert_eq!(
            texts(&table.columns[0]),
            [Some("GENE1"), Some("GENE2 "), Some("GENE3")]
        );
        // A second field would become the row index in pandas: refused.
        assert!(parse_table_with("GENE1\tx\n", &gene_list).is_err());
        let by_name = ReadOptions {
            all_text: true,
            use_columns: Some(UseColumns::Names(vec!["a".into(), "c".into()])),
            ..ReadOptions::default()
        };
        let table = parse_table_with("x\tc\tb\ta\n1\t2\t3\t4\n", &by_name).unwrap();
        assert_eq!(table.names, ["c", "a"]);
        assert_eq!(texts(&table.columns[1]), [Some("4")]);
        let error = parse_table_with("x\tc\n1\t2\n", &by_name)
            .unwrap_err()
            .to_string();
        assert!(error.contains("Usecols do not match columns"), "{error}");
        // A quoted comment byte is text.
        let quoted = parse_table_with("a\n\"x#y\"\n", &gene_list).unwrap();
        assert_eq!(texts(&quoted.columns[0]), [Some("a"), Some("x#y")]);
    }

    #[test]
    fn later_column_types_are_written_like_to_csv() {
        let table = Table {
            names: vec!["x".into(), "y".into(), "o".into()],
            columns: vec![
                Column::NullableInt(vec![Some(1), None]),
                Column::NullableBool(vec![Some(true), None]),
                Column::Object(vec![
                    ObjectValue::Float(1e-5),
                    ObjectValue::Text(String::new()),
                ]),
            ],
            rows: 2,
        };
        let mut out = Vec::new();
        write_table_to(&mut out, &table).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "x\ty\to\n1\tTrue\t1e-05\n\t\t\n"
        );
        let single = Table {
            names: vec!["a".into()],
            columns: vec![Column::Float(vec![f64::NAN, 1.0])],
            rows: 2,
        };
        let mut out = Vec::new();
        write_table_to(&mut out, &single).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), "a\n\"\"\n1.0\n");
        assert_eq!(table.columns[0].dtype_name(), "Int64");
        assert_eq!(table.columns[2].missing_mask(), [false, false]);
        assert_eq!(
            table.take_rows(&[1]).columns[0],
            Column::NullableInt(vec![None])
        );
        assert!(
            table
                .column("absent")
                .unwrap_err()
                .to_string()
                .contains("KeyError")
        );
    }
}
