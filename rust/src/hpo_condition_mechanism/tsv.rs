//! `csv.DictReader(handle, delimiter="\t")` over a file opened with
//! `newline=""`: the csv module's default (non-strict) dialect, its 131072
//! character field limit, blank rows skipped, a short row filled with None
//! and the header read as the first record.
use super::*;
use std::collections::HashMap;

const FIELD_LIMIT: usize = 131_072;

/// One row as `DictReader` returns it; `None` is the restval of a short row.
pub(super) type Row = HashMap<String, Option<String>>;

#[derive(Clone, Copy, PartialEq)]
enum State {
    StartRecord,
    StartField,
    InField,
    InQuotedField,
    QuoteInQuotedField,
    EatCrnl,
}

pub(super) struct Reader {
    text: String,
    at: usize,
    pub(super) fieldnames: Option<Vec<String>>,
}

fn csv_error(message: impl Into<String>) -> Box<dyn std::error::Error> {
    error("Error", message)
}

impl Reader {
    /// `_open_text(path)`: gzip when the file name ends in `.gz`, UTF-8.
    pub(super) fn open(path: &Path) -> Result<Self> {
        let bytes = read_maybe_gzip(path, file_name_ends_gz(path))?;
        let text = crate::python_io::decode_utf8(&bytes)
            .map_err(|e| error("UnicodeDecodeError", e))?
            .to_owned();
        Ok(Self {
            text,
            at: 0,
            fieldnames: None,
        })
    }

    /// The next physical line, keeping its `\n`, `\r` or `\r\n` ending.
    fn line(&mut self) -> Option<&str> {
        if self.at >= self.text.len() {
            return None;
        }
        let start = self.at;
        let bytes = self.text.as_bytes();
        let mut end = start;
        while end < bytes.len() {
            match bytes[end] {
                b'\n' => {
                    end += 1;
                    break;
                }
                b'\r' => {
                    end += 1;
                    if bytes.get(end) == Some(&b'\n') {
                        end += 1;
                    }
                    break;
                }
                _ => end += 1,
            }
        }
        self.at = end;
        Some(&self.text[start..end])
    }

    /// `next(csv.reader)`: one record, `[]` for a blank line, `None` at the end.
    fn record(&mut self) -> Result<Option<Vec<String>>> {
        let mut fields = Vec::new();
        let mut field = String::new();
        let mut length = 0usize;
        let mut state = State::StartRecord;
        loop {
            let Some(line) = self.line().map(str::to_owned) else {
                if length != 0 || state == State::InQuotedField {
                    fields.push(field);
                    return Ok(Some(fields));
                }
                return Ok(None);
            };
            let add = |c: char, field: &mut String, length: &mut usize| -> Result<()> {
                if *length >= FIELD_LIMIT {
                    return Err(csv_error(format!(
                        "field larger than field limit ({FIELD_LIMIT})"
                    )));
                }
                field.push(c);
                *length += 1;
                Ok(())
            };
            for c in line.chars().map(Some).chain(std::iter::once(None)) {
                let newline = matches!(c, Some('\n' | '\r'));
                let eol = c.is_none();
                if state == State::StartRecord {
                    if eol {
                        return Ok(Some(fields));
                    }
                    if newline {
                        state = State::EatCrnl;
                        continue;
                    }
                    state = State::StartField;
                }
                match state {
                    State::StartField => {
                        if newline || eol {
                            fields.push(std::mem::take(&mut field));
                            length = 0;
                            state = if eol {
                                State::StartRecord
                            } else {
                                State::EatCrnl
                            };
                        } else if c == Some('"') {
                            state = State::InQuotedField;
                        } else if c == Some('\t') {
                            fields.push(std::mem::take(&mut field));
                            length = 0;
                        } else {
                            add(c.unwrap(), &mut field, &mut length)?;
                            state = State::InField;
                        }
                    }
                    State::InField => {
                        if newline || eol {
                            fields.push(std::mem::take(&mut field));
                            length = 0;
                            state = if eol {
                                State::StartRecord
                            } else {
                                State::EatCrnl
                            };
                        } else if c == Some('\t') {
                            fields.push(std::mem::take(&mut field));
                            length = 0;
                            state = State::StartField;
                        } else {
                            add(c.unwrap(), &mut field, &mut length)?;
                        }
                    }
                    State::InQuotedField => {
                        if eol {
                        } else if c == Some('"') {
                            state = State::QuoteInQuotedField;
                        } else {
                            add(c.unwrap(), &mut field, &mut length)?;
                        }
                    }
                    State::QuoteInQuotedField => {
                        if c == Some('"') {
                            add('"', &mut field, &mut length)?;
                            state = State::InQuotedField;
                        } else if c == Some('\t') {
                            fields.push(std::mem::take(&mut field));
                            length = 0;
                            state = State::StartField;
                        } else if newline || eol {
                            fields.push(std::mem::take(&mut field));
                            length = 0;
                            state = if eol {
                                State::StartRecord
                            } else {
                                State::EatCrnl
                            };
                        } else {
                            add(c.unwrap(), &mut field, &mut length)?;
                            state = State::InField;
                        }
                    }
                    State::EatCrnl => {
                        if eol {
                            state = State::StartRecord;
                        } else if !newline {
                            return Err(csv_error(
                                "new-line character seen in unquoted field - do you need to \
open the file with newline=''?",
                            ));
                        }
                    }
                    State::StartRecord => unreachable!(),
                }
            }
            if state == State::StartRecord {
                return Ok(Some(fields));
            }
        }
    }

    /// The header check `_iter_tsv_rows` makes before its first row.
    pub(super) fn require(&mut self, path: &Path, required: &[&str]) -> Result<()> {
        self.fieldnames = self.record()?;
        let names = self.fieldnames.clone().unwrap_or_default();
        let mut missing: Vec<&str> = required
            .iter()
            .copied()
            .filter(|name| !names.iter().any(|have| have == name))
            .collect();
        missing.sort_unstable();
        if !missing.is_empty() {
            return Err(error(
                "ValueError",
                format!(
                    "{} is missing columns: {}",
                    path.display(),
                    missing.join(", ")
                ),
            ));
        }
        Ok(())
    }

    /// `next(DictReader)` after the header.
    pub(super) fn next_row(&mut self) -> Result<Option<Row>> {
        let names = self.fieldnames.clone().unwrap_or_default();
        loop {
            let Some(values) = self.record()? else {
                return Ok(None);
            };
            if values.is_empty() {
                continue;
            }
            let mut row = Row::new();
            for (name, value) in names.iter().zip(values.iter()) {
                row.insert(name.clone(), Some(value.clone()));
            }
            for name in names.iter().skip(values.len()) {
                row.insert(name.clone(), None);
            }
            return Ok(Some(row));
        }
    }
}

/// `_clean(row.get(name))`.
pub(super) fn cell(row: &Row, name: &str) -> String {
    clean(row.get(name).and_then(Option::as_deref).unwrap_or(""))
}
