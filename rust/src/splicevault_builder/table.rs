//! Streaming, chunk-local pandas C-parser tables for the SpliceVault producer.
//! The common inference and CSV writer stay shared; only this producer owns
//! its implicit index and default low-memory block concatenation rules.
use crate::cli::{Result, invalid};
use crate::pandas_table::{self, Cell, Column, ObjectValue, Table};
use std::collections::HashMap;
use std::io::{BufRead, Write};

#[derive(Clone, Copy, PartialEq)]
enum State {
    Start,
    Field,
    Text,
    Quoted,
    Quote,
}

pub(super) struct Records<R> {
    reader: R,
    first: bool,
    skip_lf: bool,
}

impl<R: BufRead> Records<R> {
    pub(super) fn new(reader: R) -> Self {
        Self {
            reader,
            first: true,
            skip_lf: false,
        }
    }

    fn byte(&mut self) -> Result<Option<u8>> {
        let buffer = self.reader.fill_buf()?;
        if buffer.is_empty() {
            return Ok(None);
        }
        let byte = buffer[0];
        self.reader.consume(1);
        Ok(Some(byte))
    }

    pub(super) fn next(&mut self) -> Result<Option<Vec<String>>> {
        let mut state = State::Start;
        let mut fields = Vec::new();
        let mut field = Vec::new();
        let finish = |fields: &mut Vec<String>, field: &mut Vec<u8>| -> Result<()> {
            let bytes = std::mem::take(field);
            let end = bytes
                .iter()
                .position(|byte| *byte == 0)
                .unwrap_or(bytes.len());
            let text = std::str::from_utf8(&bytes[..end])
                .map_err(|_| invalid("UnicodeDecodeError: input is not UTF-8"))?;
            // The pandas C tokenizer keeps the bytes after NUL internally,
            // but conversion sees the C string only up to its first NUL.
            fields.push(text.to_owned());
            Ok(())
        };
        loop {
            let Some(byte) = self.byte()? else {
                return match state {
                    State::Start => Ok(None),
                    State::Text if fields.is_empty() && field.iter().all(|b| *b == b' ') => {
                        Ok(None)
                    }
                    State::Quoted => Err(invalid("ParserError: EOF inside string")),
                    _ => {
                        finish(&mut fields, &mut field)?;
                        Ok(Some(fields))
                    }
                };
            };
            if self.first {
                self.first = false;
                if byte == 0xef {
                    let second = self.byte()?;
                    let third = self.byte()?;
                    if second == Some(0xbb) && third == Some(0xbf) {
                        continue;
                    }
                    field.push(byte);
                    if let Some(b) = second {
                        field.push(b);
                    }
                    if let Some(b) = third {
                        field.push(b);
                    }
                    state = State::Text;
                    continue;
                }
            }
            if self.skip_lf {
                self.skip_lf = false;
                if byte == b'\n' {
                    continue;
                }
            }
            let end = matches!(byte, b'\r' | b'\n');
            if end && state != State::Quoted {
                self.skip_lf = byte == b'\r';
                if state == State::Start
                    || (fields.is_empty()
                        && field.iter().all(|b| *b == b' ')
                        && state == State::Text)
                {
                    field.clear();
                    state = State::Start;
                    continue;
                }
                finish(&mut fields, &mut field)?;
                return Ok(Some(fields));
            }
            if state == State::Start {
                state = State::Field;
            }
            match state {
                State::Start => unreachable!(),
                State::Field => match byte {
                    b'"' => state = State::Quoted,
                    b'\t' => finish(&mut fields, &mut field)?,
                    _ => {
                        field.push(byte);
                        state = State::Text;
                    }
                },
                State::Text => {
                    if byte == b'\t' {
                        finish(&mut fields, &mut field)?;
                        state = State::Field;
                    } else {
                        field.push(byte);
                    }
                }
                State::Quoted => {
                    if byte == b'"' {
                        state = State::Quote;
                    } else {
                        field.push(byte);
                    }
                }
                State::Quote => match byte {
                    b'"' => {
                        field.push(byte);
                        state = State::Quoted;
                    }
                    b'\t' => {
                        finish(&mut fields, &mut field)?;
                        state = State::Field;
                    }
                    _ => {
                        field.push(byte);
                        state = State::Text;
                    }
                },
            }
        }
    }
}

#[derive(Clone)]
pub(super) struct Values {
    pub(super) column: Column,
    pub(super) unsigned: Option<Vec<u64>>,
}

fn infer(tokens: Vec<String>) -> Result<Values> {
    match pandas_table::infer_column(tokens.clone(), true) {
        Ok(column) => Ok(Values {
            column,
            unsigned: None,
        }),
        Err(error)
            if error
                .to_string()
                .contains("outside the signed 64-bit integer range") =>
        {
            let mut numbers = Vec::new();
            let mut negative = false;
            let mut large = false;
            let mut valid = true;
            for token in tokens.iter().filter(|token| !pandas_table::is_na(token)) {
                let token = token.trim_matches([' ', '\t', '\n', '\r', '\x0b', '\x0c']);
                if token.starts_with('-') {
                    negative = true;
                    continue;
                }
                let digits = token.strip_prefix('+').unwrap_or(token);
                let value = (!digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()))
                    .then(|| digits.parse::<u64>().ok())
                    .flatten();
                if let Some(value) = value {
                    large |= value > i64::MAX as u64;
                    numbers.push(value);
                } else {
                    valid = false;
                    break;
                }
            }
            if valid && large && !negative && !tokens.iter().any(|token| pandas_table::is_na(token))
            {
                let column = Column::Text(numbers.iter().map(|n| Some(n.to_string())).collect());
                return Ok(Values {
                    column,
                    unsigned: Some(numbers),
                });
            }
            // pandas cancels uint64 conversion with NA/negative values and
            // then preserves their literal spelling; true overflow uses normal NA.
            let keep_na = valid && large;
            Ok(Values {
                column: Column::Text(
                    tokens
                        .into_iter()
                        .map(|s| {
                            if !keep_na && pandas_table::is_na(&s) {
                                None
                            } else {
                                Some(s)
                            }
                        })
                        .collect(),
                ),
                unsigned: None,
            })
        }
        Err(error) => Err(error),
    }
}

fn object(column: &Column, row: usize) -> ObjectValue {
    match column.cell(row) {
        Cell::Missing => ObjectValue::Nan,
        Cell::Text(s) => ObjectValue::Text(s.to_owned()),
        Cell::Int(n) => ObjectValue::Int(n),
        Cell::Float(n) => ObjectValue::Float(n),
        Cell::Bool(b) => ObjectValue::Bool(b),
    }
}

fn concatenate(blocks: Vec<Values>) -> Values {
    if blocks.len() == 1 {
        return blocks.into_iter().next().unwrap();
    }
    let all_integer = blocks
        .iter()
        .all(|b| matches!(b.column, Column::Int(_)) && b.unsigned.is_none());
    let all_numeric = blocks
        .iter()
        .all(|b| matches!(b.column, Column::Int(_) | Column::Float(_)) || b.unsigned.is_some());
    let all_unsigned = blocks.iter().all(|b| b.unsigned.is_some());
    let all_bool = blocks.iter().all(|b| matches!(b.column, Column::Bool(_)));
    let all_text = blocks
        .iter()
        .all(|b| matches!(b.column, Column::Text(_)) && b.unsigned.is_none());
    if all_unsigned {
        let values: Vec<u64> = blocks
            .into_iter()
            .flat_map(|b| b.unsigned.unwrap())
            .collect();
        return Values {
            column: Column::Text(values.iter().map(|n| Some(n.to_string())).collect()),
            unsigned: Some(values),
        };
    }
    let column = if all_integer {
        Column::Int(
            blocks
                .iter()
                .flat_map(|b| {
                    (0..b.column.len()).map(|r| match b.column.cell(r) {
                        Cell::Int(n) => n,
                        _ => unreachable!(),
                    })
                })
                .collect(),
        )
    } else if all_numeric {
        Column::Float(
            blocks
                .iter()
                .flat_map(|b| {
                    (0..b.column.len()).map(|r| {
                        if let Some(values) = &b.unsigned {
                            values[r] as f64
                        } else {
                            match b.column.cell(r) {
                                Cell::Int(n) => n as f64,
                                Cell::Float(n) => n,
                                Cell::Missing => f64::NAN,
                                _ => unreachable!(),
                            }
                        }
                    })
                })
                .collect(),
        )
    } else if all_bool {
        Column::Bool(
            blocks
                .iter()
                .flat_map(|b| (0..b.column.len()).map(|r| b.column.cell(r) == Cell::Bool(true)))
                .collect(),
        )
    } else if all_text {
        Column::Text(
            blocks
                .iter()
                .flat_map(|b| {
                    (0..b.column.len()).map(|r| match b.column.cell(r) {
                        Cell::Text(s) => Some(s.to_owned()),
                        _ => None,
                    })
                })
                .collect(),
        )
    } else {
        Column::Object(
            blocks
                .iter()
                .flat_map(|b| (0..b.column.len()).map(|r| object(&b.column, r)))
                .collect(),
        )
    };
    Values {
        column,
        unsigned: None,
    }
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
enum Label {
    Integer(i128),
    Float(u64),
    Text(String),
    Missing,
}

fn label(column: &Values, row: usize) -> Label {
    if let Some(values) = &column.unsigned {
        return Label::Integer(i128::from(values[row]));
    }
    match column.column.cell(row) {
        Cell::Int(n) => Label::Integer(i128::from(n)),
        Cell::Bool(b) => Label::Integer(i128::from(b)),
        Cell::Float(f) if f.is_finite() && f.fract() == 0.0 => Label::Integer(f as i128),
        Cell::Float(f) => Label::Float(f.to_bits()),
        Cell::Text(s) => Label::Text(s.to_owned()),
        Cell::Missing => Label::Missing,
    }
}

pub(super) struct Chunk {
    pub(super) table: Table,
    pub(super) unsigned: HashMap<usize, Vec<u64>>,
    pub(super) same_index: Vec<Vec<usize>>,
    pub(super) duplicate_multi_index: bool,
}

pub(super) struct Reader<R> {
    records: Records<R>,
    names: Vec<String>,
    width: Option<usize>,
    done: bool,
    first: bool,
    size: usize,
}

impl<R: BufRead> Reader<R> {
    pub(super) fn new(reader: R, size: usize) -> Result<Self> {
        let mut records = Records::new(reader);
        let names = pandas_table::column_names(
            records
                .next()?
                .ok_or_else(|| invalid("EmptyDataError: No columns to parse from file"))?,
        );
        Ok(Self {
            records,
            names,
            width: None,
            done: false,
            first: true,
            size,
        })
    }

    pub(super) fn next(&mut self) -> Result<Option<Chunk>> {
        if self.done {
            return Ok(None);
        }
        let mut raw = Vec::new();
        for _ in 0..self.size {
            let Some(mut row) = self.records.next()? else {
                self.done = true;
                break;
            };
            let width = *self.width.get_or_insert(row.len().max(self.names.len()));
            // The pinned pandas C chunk reader does not check extra fields
            // on the first buffered row of a new read. They are discarded,
            // rather than changing the implicit index inferred at the start.
            if row.len() > width && !raw.is_empty() {
                return Err(invalid(
                    "ParserError: Error tokenizing data: too many fields",
                ));
            }
            row.resize(width, String::new());
            raw.push(row);
        }
        if raw.is_empty() && !self.first {
            return Ok(None);
        }
        self.first = false;
        let width = self.width.unwrap_or(self.names.len());
        let leading = width - self.names.len();
        // pandas' C low-memory reader chooses a power of two below 1MiB / width.
        let heuristic = (1 << 20) / width.max(1);
        let mut block_size = 1;
        while block_size * 2 < heuristic {
            block_size *= 2;
        }
        let mut inferred = Vec::new();
        for col in 0..width {
            let values = if raw.is_empty() {
                infer(Vec::new())?
            } else {
                let blocks = raw
                    .chunks(block_size)
                    .map(|rows| infer(rows.iter().map(|r| r[col].clone()).collect()))
                    .collect::<Result<Vec<_>>>()?;
                concatenate(blocks)
            };
            inferred.push(values);
        }
        let mut groups: HashMap<Vec<Label>, Vec<usize>> = HashMap::new();
        if leading > 0 {
            for row in 0..raw.len() {
                groups
                    .entry(inferred[..leading].iter().map(|c| label(c, row)).collect())
                    .or_default()
                    .push(row);
            }
        }
        let same_index = (0..raw.len())
            .map(|row| {
                if leading == 0 {
                    vec![row]
                } else {
                    groups[&inferred[..leading]
                        .iter()
                        .map(|c| label(c, row))
                        .collect::<Vec<_>>()]
                        .clone()
                }
            })
            .collect();
        let duplicate_multi_index = leading > 1 && groups.values().any(|rows| rows.len() > 1);
        let columns: Vec<_> = inferred.into_iter().skip(leading).collect();
        let unsigned = columns
            .iter()
            .enumerate()
            .filter_map(|(i, c)| c.unsigned.clone().map(|v| (i, v)))
            .collect();
        Ok(Some(Chunk {
            table: Table {
                names: self.names.clone(),
                columns: columns.into_iter().map(|c| c.column).collect(),
                rows: raw.len(),
            },
            unsigned,
            same_index,
            duplicate_multi_index,
        }))
    }
}

/// Skip the shared writer's exactly encoded header on later chunks. This
/// also handles quoted column names containing physical line breaks.
pub(super) fn write(out: &mut impl Write, table: &Table, header: bool) -> Result<()> {
    if header {
        pandas_table::write_table_to(out, table)?;
        return Ok(());
    }
    let empty = Table {
        names: table.names.clone(),
        columns: table.columns.iter().map(|c| c.take(&[])).collect(),
        rows: 0,
    };
    let mut encoded = Vec::new();
    pandas_table::write_table_to(&mut encoded, &empty)?;
    struct Skip<'a, W> {
        out: &'a mut W,
        remaining: usize,
    }
    impl<W: Write> Write for Skip<'_, W> {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            let skipped = self.remaining.min(buf.len());
            self.remaining -= skipped;
            if skipped == buf.len() {
                Ok(skipped)
            } else {
                self.out.write(&buf[skipped..]).map(|n| n + skipped)
            }
        }
        fn flush(&mut self) -> std::io::Result<()> {
            self.out.flush()
        }
    }
    pandas_table::write_table_to(
        &mut Skip {
            out,
            remaining: encoded.len(),
        },
        table,
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn quoting_na_and_chunk_local_inference_match_pandas() {
        let mut reader = Reader::new(Cursor::new("a\tb\n007\t\"x\ty\"\ntext\tNA\n"), 1).unwrap();
        let first = reader.next().unwrap().unwrap();
        let second = reader.next().unwrap().unwrap();
        assert_eq!(first.table.cell(0, 0), Cell::Int(7));
        assert_eq!(first.table.cell(1, 0), Cell::Text("x\ty"));
        assert_eq!(second.table.cell(0, 0), Cell::Text("text"));
        assert_eq!(second.table.cell(1, 0), Cell::Missing);
        assert!(reader.next().unwrap().is_none());
    }

    #[test]
    fn first_row_surplus_is_discarded_only_at_the_chunk_boundary() {
        let input = "a\tb\n1\t2\n3\t4\textra\n";
        let mut single = Reader::new(Cursor::new(input), 1).unwrap();
        single.next().unwrap();
        assert_eq!(
            single.next().unwrap().unwrap().table.cell(1, 0),
            Cell::Int(4)
        );
        let mut together = Reader::new(Cursor::new(input), 2).unwrap();
        assert!(together.next().is_err());
    }

    #[test]
    fn duplicate_implicit_labels_retain_scalar_assignment_scope() {
        let mut reader = Reader::new(Cursor::new("a\tb\nx\t1\t2\nx\t3\t4\n"), 2).unwrap();
        let chunk = reader.next().unwrap().unwrap();
        assert_eq!(chunk.same_index, vec![vec![0, 1], vec![0, 1]]);
        assert!(!chunk.duplicate_multi_index);
        let mut reader = Reader::new(Cursor::new("a\tb\nx\ty\t1\t2\nx\ty\t3\t4\n"), 2).unwrap();
        assert!(reader.next().unwrap().unwrap().duplicate_multi_index);
    }

    #[test]
    fn quoted_header_newlines_are_skipped_exactly_when_appending() {
        let mut reader = Reader::new(Cursor::new("\"a\nb\"\tc\n1\t2\n"), 5).unwrap();
        let chunk = reader.next().unwrap().unwrap();
        let mut output = Vec::new();
        write(&mut output, &chunk.table, false).unwrap();
        assert_eq!(output, b"1\t2\n");
    }

    #[test]
    fn final_unquoted_spaces_are_blank_but_quoted_spaces_are_a_row() {
        let mut records = Records::new(Cursor::new("a\n   "));
        assert!(records.next().unwrap().is_some());
        assert!(records.next().unwrap().is_none());
        let mut records = Records::new(Cursor::new("\"   \""));
        assert_eq!(records.next().unwrap(), Some(vec!["   ".into()]));
    }

    #[test]
    fn low_memory_blocks_keep_numeric_objects_before_later_text() {
        let mut input = (0..64)
            .map(|n| format!("c{n}"))
            .collect::<Vec<_>>()
            .join("\t");
        input.push('\n');
        for row in 0..9000 {
            for col in 0..64 {
                if col > 0 {
                    input.push('\t');
                }
                input.push_str(if col == 63 && row >= 8192 {
                    "word"
                } else {
                    "007"
                });
            }
            input.push('\n');
        }
        let mut reader = Reader::new(Cursor::new(input), 9000).unwrap();
        let chunk = reader.next().unwrap().unwrap();
        assert_eq!(chunk.table.cell(63, 0), Cell::Int(7));
        assert_eq!(chunk.table.cell(63, 8192), Cell::Text("word"));
    }
}
