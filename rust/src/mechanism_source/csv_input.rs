//! Python csv.reader(strict=True) record framing around the mature CSV decoder.
//! No pandas type inference applies to this acquisition-only source reader.
use super::*;
use crate::python_lines::Lines;
use std::io::Read;

pub(super) struct Rows<R: Read> {
    lines: Lines<R>,
    first: bool,
    delimiter: u8,
    dosage: bool,
}
impl<R: Read> Rows<R> {
    pub(super) fn new(source: R, kind: &str) -> Self {
        Self {
            lines: Lines::new(source),
            first: true,
            delimiter: if kind == "g2p_ddg2p" { b',' } else { b'\t' },
            dosage: kind.starts_with("clingen_dosage"),
        }
    }
    pub(super) fn next(&mut self) -> Result<Option<Vec<String>>> {
        // 0=start field, 1=unquoted, 2=quoted, 3=closing quote.
        let mut state = 0;
        let mut length = 0;
        let mut text = String::new();
        loop {
            let Some(line) = self.lines.next_csv(&mut self.first)? else {
                if state == 2 {
                    return Err(error("Error", "unexpected end of data"));
                }
                if text.is_empty() {
                    return Ok(None);
                }
                break;
            };
            // The source filters physical lines before csv.reader consumes
            // them, including a comment inside a quoted multiline field.
            if self.dosage && line.starts_with('#') && !line.contains('\t') {
                continue;
            }
            let mut ended = false;
            for c in line.chars() {
                if ended {
                    continue;
                } // only the LF half of a final CRLF
                let newline = matches!(c, '\r' | '\n');
                match state {
                    0 if c == '"' => state = 2,
                    0 | 1 if c == char::from(self.delimiter) => {
                        state = 0;
                        length = 0;
                    }
                    0 | 1 if newline => ended = true,
                    0 | 1 => {
                        state = 1;
                        length += 1;
                    }
                    2 if c == '"' => state = 3,
                    2 => length += 1,
                    3 if c == '"' => {
                        state = 2;
                        length += 1;
                    }
                    3 if c == char::from(self.delimiter) => {
                        state = 0;
                        length = 0;
                    }
                    3 if newline => ended = true,
                    3 => {
                        return Err(error(
                            "Error",
                            format!("'{}' expected after '\"'", char::from(self.delimiter)),
                        ));
                    }
                    _ => unreachable!(),
                }
                if length > 131072 {
                    return Err(error("Error", "field larger than field limit (131072)"));
                }
            }
            text.push_str(&line);
            if ended {
                break;
            }
        }
        if text.chars().all(|c| matches!(c, '\r' | '\n')) {
            return Ok(Some(vec![]));
        }
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .flexible(true)
            .delimiter(self.delimiter)
            .from_reader(text.as_bytes());
        let row = reader
            .records()
            .next()
            .transpose()
            .map_err(|e| error("RuntimeError", e.to_string()))?
            .ok_or_else(|| error("RuntimeError", "CSV decoder lost a nonempty framed record"))?;
        Ok(Some(row.iter().map(str::to_owned).collect()))
    }
}
