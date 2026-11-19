//! Reading and writing tab-separated tables the way the Python scripts do.
//!
//! * `Table::read` reads like `pandas.read_csv(path, sep="\t", dtype=str,
//!   keep_default_na=False)`: the first record is the header, `"` quotes a
//!   field, blank lines are skipped and every value stays text.
//! * `TableWriter` writes like Python's `csv` module (`QUOTE_MINIMAL`, tab
//!   separator): a field is quoted when it holds a tab, a `"`, or a line-break
//!   character, and quotes inside are doubled.
//! * `python_lines` and `universal_lines` split a file into lines as Python's
//!   text files do with `newline=""` and with the default newline handling.
//! * `python_csv_records` reads like `csv.reader(open(path), delimiter="\t")`.

use crate::hpo_gene_ranking::Result;
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};

/// Opens a file for reading, decompressing it when the name ends in `.gz`.
pub fn open_text(path: &str) -> Result<Box<dyn Read>> {
    let file = File::open(path).map_err(|error| format!("{path}: {error}"))?;
    if path.ends_with(".gz") {
        Ok(Box::new(MultiGzDecoder::new(BufReader::with_capacity(
            1 << 20,
            file,
        ))))
    } else {
        Ok(Box::new(BufReader::with_capacity(1 << 20, file)))
    }
}

/// Reads a whole file as bytes (decompressing `.gz`).
pub fn read_bytes(path: &str) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    open_text(path)?
        .read_to_end(&mut bytes)
        .map_err(|error| format!("{path}: {error}"))?;
    Ok(bytes)
}

/// A table of text values: the wanted columns of every row, in the order asked for.
pub struct Table {
    pub rows: Vec<Vec<String>>,
}

impl Table {
    /// Reads the columns named in `wanted` (in that order). A wanted column
    /// that the file lacks is an error unless it is listed in `optional`,
    /// in which case its values are empty and `has_column` reports false.
    pub fn read(path: &str, wanted: &[&str], optional: &[&str]) -> Result<(Table, Vec<bool>)> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .flexible(true)
            .quoting(true)
            .double_quote(true)
            .terminator(csv::Terminator::CRLF)
            .from_reader(open_text(path)?);
        let mut records = reader.byte_records();
        let header: Vec<String> = match records.next() {
            None => Vec::new(),
            Some(record) => {
                let record = record.map_err(|error| format!("{path}: {error}"))?;
                record
                    .iter()
                    .map(|field| String::from_utf8_lossy(field).into_owned())
                    .collect()
            }
        };
        let mut positions = Vec::with_capacity(wanted.len());
        let mut present = Vec::with_capacity(wanted.len());
        for name in wanted {
            match header.iter().position(|column| column == name) {
                Some(position) => {
                    positions.push(Some(position));
                    present.push(true);
                }
                None if optional.contains(name) => {
                    positions.push(None);
                    present.push(false);
                }
                None => return Err(format!("{path}: no column named {name}")),
            }
        }
        let mut rows = Vec::new();
        for record in records {
            let record = record.map_err(|error| format!("{path}: {error}"))?;
            let row = positions
                .iter()
                .map(|position| match position.and_then(|p| record.get(p)) {
                    Some(field) => String::from_utf8_lossy(field).into_owned(),
                    None => String::new(),
                })
                .collect();
            rows.push(row);
        }
        Ok((Table { rows }, present))
    }

    /// Row numbers per value of one column, in file order.
    pub fn index_by(&self, column: usize) -> HashMap<String, Vec<usize>> {
        let mut index: HashMap<String, Vec<usize>> = HashMap::new();
        for (row_number, row) in self.rows.iter().enumerate() {
            index
                .entry(row[column].clone())
                .or_default()
                .push(row_number);
        }
        index
    }
}

/// A value as the Python scripts hand it to the csv writer.
#[derive(Clone, Debug)]
pub enum Cell {
    Text(String),
    Int(i64),
    Float(f64),
}

impl Cell {
    pub fn render(&self) -> String {
        match self {
            Cell::Text(text) => text.clone(),
            Cell::Int(value) => value.to_string(),
            Cell::Float(value) => crate::hpo_gene_ranking::python_text::float_repr(*value),
        }
    }
}

impl From<&str> for Cell {
    fn from(text: &str) -> Self {
        Cell::Text(text.to_owned())
    }
}

impl From<String> for Cell {
    fn from(text: String) -> Self {
        Cell::Text(text)
    }
}

impl From<&String> for Cell {
    fn from(text: &String) -> Self {
        Cell::Text(text.clone())
    }
}

impl From<i64> for Cell {
    fn from(value: i64) -> Self {
        Cell::Int(value)
    }
}

impl From<usize> for Cell {
    fn from(value: usize) -> Self {
        Cell::Int(value as i64)
    }
}

impl From<f64> for Cell {
    fn from(value: f64) -> Self {
        Cell::Float(value)
    }
}

/// One output row as (field name, value) pairs; fields the header names but
/// the record lacks are written empty, extra fields are ignored (as
/// `csv.DictWriter(..., extrasaction="ignore", restval="")` does).
pub type Record = Vec<(&'static str, Cell)>;

/// Writes rows like Python's `csv.writer(handle, delimiter="\t",
/// lineterminator=TERMINATOR)` with the default `QUOTE_MINIMAL`.
pub struct TableWriter {
    out: BufWriter<File>,
    terminator: &'static str,
    path: String,
}

impl TableWriter {
    pub fn create(path: &str, terminator: &'static str) -> Result<TableWriter> {
        if let Some(parent) = std::path::Path::new(path).parent()
            && !parent.as_os_str().is_empty()
        {
            std::fs::create_dir_all(parent).map_err(|error| format!("{path}: {error}"))?;
        }
        let file = File::create(path).map_err(|error| format!("{path}: {error}"))?;
        Ok(TableWriter {
            out: BufWriter::with_capacity(1 << 20, file),
            terminator,
            path: path.to_owned(),
        })
    }

    fn quoted(&self, field: &str) -> bool {
        field
            .chars()
            .any(|c| c == '\t' || c == '"' || c == '\n' || c == '\r' || self.terminator.contains(c))
    }

    pub fn write_row<S: AsRef<str>>(&mut self, fields: &[S]) -> Result<()> {
        let mut line = String::new();
        if fields.len() == 1 && fields[0].as_ref().is_empty() {
            // Python writes a lone empty field as "" so that the row is not blank.
            line.push_str("\"\"");
        } else {
            for (index, field) in fields.iter().enumerate() {
                if index > 0 {
                    line.push('\t');
                }
                let field = field.as_ref();
                if self.quoted(field) {
                    line.push('"');
                    line.push_str(&field.replace('"', "\"\""));
                    line.push('"');
                } else {
                    line.push_str(field);
                }
            }
        }
        line.push_str(self.terminator);
        self.out
            .write_all(line.as_bytes())
            .map_err(|error| format!("{}: {error}", self.path))
    }

    pub fn write_record(&mut self, fieldnames: &[&str], record: &Record) -> Result<()> {
        let fields: Vec<String> = fieldnames
            .iter()
            .map(|name| {
                record
                    .iter()
                    .find(|(key, _)| key == name)
                    .map(|(_, value)| value.render())
                    .unwrap_or_default()
            })
            .collect();
        self.write_row(&fields)
    }

    pub fn finish(mut self) -> Result<()> {
        self.out
            .flush()
            .map_err(|error| format!("{}: {error}", self.path))
    }
}

/// Writes a header and records, as the scripts' `write_tsv` helpers do.
pub fn write_records(
    path: &str,
    fieldnames: &[&str],
    records: &[Record],
    terminator: &'static str,
) -> Result<()> {
    let mut writer = TableWriter::create(path, terminator)?;
    writer.write_row(fieldnames)?;
    for record in records {
        writer.write_record(fieldnames, record)?;
    }
    writer.finish()
}

/// Calls `handle` with each line of a file as Python's text mode yields them:
/// with `universal` false like `open(path, newline="")` (a line ends after
/// "\n", "\r\n" or a lone "\r" and keeps its ending), with `universal` true
/// like `open(path)` (those endings become "\n"). Invalid UTF-8 becomes
/// U+FFFD, as with `errors="replace"`. The file is read piece by piece.
pub fn for_each_line(
    path: &str,
    universal: bool,
    mut handle: impl FnMut(&str) -> Result<()>,
) -> Result<()> {
    use std::io::BufRead;
    let mut reader = BufReader::with_capacity(1 << 20, open_text(path)?);
    let mut chunk = Vec::new();
    loop {
        chunk.clear();
        let read = reader
            .read_until(b'\n', &mut chunk)
            .map_err(|error| format!("{path}: {error}"))?;
        if read == 0 {
            return Ok(());
        }
        let text = String::from_utf8_lossy(&chunk);
        let bytes = text.as_bytes();
        let mut start = 0;
        let mut index = 0;
        while index < bytes.len() {
            let end = match bytes[index] {
                b'\r' if index + 1 < bytes.len() && bytes[index + 1] == b'\n' => {
                    Some((index, index + 2))
                }
                b'\r' | b'\n' => Some((index, index + 1)),
                _ => None,
            };
            if let Some((body_end, line_end)) = end {
                if universal {
                    handle(&format!("{}\n", &text[start..body_end]))?;
                } else {
                    handle(&text[start..line_end])?;
                }
                start = line_end;
                index = line_end;
            } else {
                index += 1;
            }
        }
        if start < bytes.len() {
            handle(&text[start..])?;
        }
    }
}

/// Lines of a text file opened with `newline=""`: a line ends after "\n",
/// "\r\n" or a lone "\r", and keeps its ending. Invalid UTF-8 becomes U+FFFD
/// (`errors="replace"`).
pub fn python_lines(bytes: &[u8]) -> Vec<String> {
    let text = String::from_utf8_lossy(bytes);
    let mut lines = Vec::new();
    let mut start = 0;
    let raw = text.as_bytes();
    let mut index = 0;
    while index < raw.len() {
        match raw[index] {
            b'\n' => {
                lines.push(text[start..=index].to_owned());
                start = index + 1;
            }
            b'\r' => {
                if index + 1 < raw.len() && raw[index + 1] == b'\n' {
                    index += 1;
                }
                lines.push(text[start..=index].to_owned());
                start = index + 1;
            }
            _ => {}
        }
        index += 1;
    }
    if start < raw.len() {
        lines.push(text[start..].to_owned());
    }
    lines
}

/// Lines of a text file opened with the default newline handling, where
/// "\r\n" and a lone "\r" become "\n"; each line keeps its "\n".
pub fn universal_lines(bytes: &[u8]) -> Vec<String> {
    python_lines(bytes)
        .into_iter()
        .map(|line| {
            if let Some(body) = line.strip_suffix("\r\n") {
                format!("{body}\n")
            } else if let Some(body) = line.strip_suffix('\r') {
                format!("{body}\n")
            } else {
                line
            }
        })
        .collect()
}

/// Records of `csv.reader(open(path, encoding="utf-8", errors="replace"),
/// delimiter="\t")`: the default dialect with a tab separator, reading lines
/// that the default newline handling already turned into "\n" endings.
pub fn python_csv_records(bytes: &[u8]) -> Vec<Vec<String>> {
    #[derive(PartialEq)]
    enum State {
        StartRecord,
        StartField,
        InField,
        InQuotedField,
        QuoteInQuotedField,
        EatCrnl,
    }
    let mut records = Vec::new();
    let mut fields: Vec<String> = Vec::new();
    let mut field = String::new();
    let mut state = State::StartRecord;
    for line in universal_lines(bytes) {
        for c in line.chars() {
            loop {
                match state {
                    State::StartRecord => {
                        if c == '\n' || c == '\r' {
                            state = State::EatCrnl;
                            break;
                        }
                        state = State::StartField;
                        continue;
                    }
                    State::StartField => {
                        if c == '\n' || c == '\r' {
                            fields.push(std::mem::take(&mut field));
                            state = State::EatCrnl;
                        } else if c == '"' {
                            state = State::InQuotedField;
                        } else if c == '\t' {
                            fields.push(std::mem::take(&mut field));
                            state = State::StartField;
                        } else {
                            field.push(c);
                            state = State::InField;
                        }
                        break;
                    }
                    State::InField => {
                        if c == '\n' || c == '\r' {
                            fields.push(std::mem::take(&mut field));
                            state = State::EatCrnl;
                        } else if c == '\t' {
                            fields.push(std::mem::take(&mut field));
                            state = State::StartField;
                        } else {
                            field.push(c);
                        }
                        break;
                    }
                    State::InQuotedField => {
                        if c == '"' {
                            state = State::QuoteInQuotedField;
                        } else {
                            field.push(c);
                        }
                        break;
                    }
                    State::QuoteInQuotedField => {
                        if c == '"' {
                            field.push('"');
                            state = State::InQuotedField;
                        } else if c == '\t' {
                            fields.push(std::mem::take(&mut field));
                            state = State::StartField;
                        } else if c == '\n' || c == '\r' {
                            fields.push(std::mem::take(&mut field));
                            state = State::EatCrnl;
                        } else {
                            field.push(c);
                            state = State::InField;
                        }
                        break;
                    }
                    State::EatCrnl => {
                        if c == '\n' || c == '\r' {
                            break;
                        }
                        // Not reachable: each line ends with its single "\n".
                        state = State::StartRecord;
                        continue;
                    }
                }
            }
            if state == State::EatCrnl {
                records.push(std::mem::take(&mut fields));
                state = State::StartRecord;
            }
        }
        // A record whose quoted field spans lines continues on the next line.
    }
    // End of input: a record without a final line break still counts.
    match state {
        State::StartRecord | State::EatCrnl => {}
        _ => {
            fields.push(std::mem::take(&mut field));
            records.push(std::mem::take(&mut fields));
        }
    }
    records
}

/// Rows of `csv.DictReader`: the first record names the columns; later fully
/// empty records are skipped; a short row leaves later columns empty (Python
/// stores None there, which the csv writer writes as an empty field).
pub fn python_dict_rows(bytes: &[u8]) -> (Vec<String>, Vec<HashMap<String, String>>) {
    let mut records = python_csv_records(bytes).into_iter();
    let fieldnames = records.next().unwrap_or_default();
    let mut rows = Vec::new();
    for record in records {
        if record.is_empty() {
            continue;
        }
        let mut row = HashMap::new();
        for (name, value) in fieldnames.iter().zip(record.iter()) {
            row.insert(name.clone(), value.clone());
        }
        for name in fieldnames.iter().skip(record.len()) {
            row.insert(name.clone(), String::new());
        }
        rows.push(row);
    }
    (fieldnames, rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn python_lines_keep_endings_and_split_lone_carriage_returns() {
        assert_eq!(
            python_lines(b"a\tb\r\nc\rd\ne"),
            vec!["a\tb\r\n", "c\r", "d\n", "e"]
        );
        assert_eq!(universal_lines(b"a\r\nb\rc"), vec!["a\n", "b\n", "c"]);
    }

    #[test]
    fn python_csv_reader_rules() {
        let records = python_csv_records(b"A\tB\n\"x\"y\tz\n\"a\tb\"\t\"c\"\"d\"\n\nlast");
        assert_eq!(
            records,
            vec![
                vec!["A".to_string(), "B".to_string()],
                vec!["xy".to_string(), "z".to_string()],
                vec!["a\tb".to_string(), "c\"d".to_string()],
                vec![],
                vec!["last".to_string()],
            ]
        );
    }
}
