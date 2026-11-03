//! Shared Python TextIOWrapper line boundaries, without source-specific policy.
use crate::cli::{Result, invalid};
use std::io::Read;

/// TextIOWrapper's UTF-8 decoding and universal CR/LF/CRLF line boundaries.
/// Decode each input block before returning its lines, so invalid UTF-8 does
/// not silently publish earlier lines from the same decoded block.
pub(crate) struct Lines<R: Read> {
    source: R,
    latin1: bool,
    text: String,
    at: usize,
    pending: Vec<u8>,
    eof: bool,
    skip_lf: bool,
}
impl<R: Read> Lines<R> {
    pub(crate) fn new(source: R) -> Self {
        Self {
            source,
            latin1: false,
            text: String::new(),
            at: 0,
            pending: Vec::new(),
            eof: false,
            skip_lf: false,
        }
    }
    /// The same TextIOWrapper block and universal-newline mechanism with the
    /// source's explicit one-byte Latin-1 decoder. Existing consumers use UTF-8.
    pub(crate) fn latin1(source: R) -> Self {
        let mut lines = Self::new(source);
        lines.latin1 = true;
        lines
    }
    fn fill(&mut self) -> Result<bool> {
        if self.eof {
            return Ok(false);
        }
        let mut block = [0_u8; 8192];
        let count = self.source.read(&mut block)?;
        self.eof = count == 0;
        self.pending.extend_from_slice(&block[..count]);
        if self.latin1 {
            self.text = self.pending.iter().map(|byte| char::from(*byte)).collect();
            self.pending.clear();
            self.at = 0;
            return Ok(!self.eof || !self.text.is_empty());
        }
        let complete = match std::str::from_utf8(&self.pending) {
            Ok(_) => self.pending.len(),
            Err(error) if error.error_len().is_none() && !self.eof => error.valid_up_to(),
            Err(_) => return Err(invalid("UnicodeDecodeError: input is not UTF-8")),
        };
        self.text = std::str::from_utf8(&self.pending[..complete])
            .unwrap()
            .to_owned();
        self.pending.drain(..complete);
        self.at = 0;
        Ok(!self.eof || !self.text.is_empty())
    }
    pub(crate) fn next(&mut self) -> Result<Option<String>> {
        let mut line = String::new();
        loop {
            if self.at == self.text.len() && !self.fill()? {
                return Ok((!line.is_empty()).then_some(line));
            }
            if self.at == self.text.len() {
                continue;
            }
            if self.skip_lf {
                self.skip_lf = false;
                if self.text.as_bytes()[self.at] == b'\n' {
                    self.at += 1;
                    continue;
                }
            }
            let rest = &self.text[self.at..];
            if let Some(end) = rest.bytes().position(|byte| matches!(byte, b'\r' | b'\n')) {
                line.push_str(&rest[..end]);
                self.skip_lf = rest.as_bytes()[end] == b'\r';
                self.at += end + 1;
                return Ok(Some(line));
            }
            line.push_str(rest);
            self.at = self.text.len();
        }
    }

    /// TextIOWrapper.readline(size), retaining the translated newline and
    /// counting Unicode characters rather than UTF-8 bytes. Existing unbounded
    /// and CSV consumers keep their original code and decoding boundaries.
    pub(crate) fn readline(&mut self, limit: usize) -> Result<String> {
        let mut line = String::new();
        let mut count = 0;
        while count < limit {
            if self.at == self.text.len() {
                let filled = self.fill_for_readline()?;
                if !filled {
                    break;
                }
            }
            if self.at == self.text.len() {
                continue;
            }
            let c = self.text[self.at..].chars().next().unwrap();
            self.at += c.len_utf8();
            if self.skip_lf {
                self.skip_lf = false;
                if c == '\n' {
                    continue;
                }
            }
            if matches!(c, '\r' | '\n') {
                // IncrementalNewlineDecoder retains a trailing CR until the
                // next decoded block or EOF. Its next-block decoding failure
                // therefore precedes publication of this bounded line.
                if c == '\r' {
                    while self.at == self.text.len() && !self.eof {
                        self.fill_for_readline()?;
                    }
                }
                self.skip_lf = c == '\r';
                line.push('\n');
                break;
            }
            line.push(c);
            count += 1;
        }
        Ok(line)
    }

    fn fill_for_readline(&mut self) -> Result<bool> {
        self.fill().map_err(|cause| {
            if cause.to_string() == "UnicodeDecodeError: input is not UTF-8" {
                crate::builder_error::category(
                    "UnicodeDecodeError",
                    invalid(crate::python_io::decode_utf8(&self.pending).unwrap_err()),
                )
            } else {
                cause
            }
        })
    }

    /// The `newline="", encoding="utf-8-sig"` form used by csv.reader.
    /// Keep the established universal-newline `next` body unchanged. This
    /// additional consumer reuses its block decoder but retains actual CR/LF
    /// bytes, including CRLF across a decode boundary. `first` belongs to this
    /// consumer and must initially be true; do not mix the two line methods.
    pub(crate) fn next_csv(&mut self, first: &mut bool) -> Result<Option<String>> {
        let mut line = String::new();
        let mut cr = false;
        loop {
            if self.at == self.text.len() {
                let initial = *first;
                *first = false;
                let filled = self.fill().map_err(|cause| {
                    if cause.to_string() == "UnicodeDecodeError: input is not UTF-8" {
                        let bytes = if initial {
                            self.pending
                                .strip_prefix(&[0xef, 0xbb, 0xbf])
                                .unwrap_or(&self.pending)
                        } else {
                            &self.pending
                        };
                        crate::builder_error::category(
                            "UnicodeDecodeError",
                            invalid(crate::python_io::decode_utf8(bytes).unwrap_err()),
                        )
                    } else {
                        cause
                    }
                })?;
                if initial && self.text.starts_with('\u{feff}') {
                    self.at = 3;
                }
                if !filled {
                    return Ok((!line.is_empty()).then_some(line));
                }
                if self.at == self.text.len() {
                    continue;
                }
            }
            if cr {
                if self.text.as_bytes()[self.at] == b'\n' {
                    line.push('\n');
                    self.at += 1;
                }
                return Ok(Some(line));
            }
            let rest = &self.text[self.at..];
            if let Some(end) = rest.bytes().position(|byte| matches!(byte, b'\r' | b'\n')) {
                line.push_str(&rest[..=end]);
                cr = rest.as_bytes()[end] == b'\r';
                self.at += end + 1;
                if !cr {
                    return Ok(Some(line));
                }
            } else {
                line.push_str(rest);
                self.at = self.text.len();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn latin1_keeps_byte_decoding_and_crlf_block_boundaries() {
        let mut bytes = vec![0xff; 8191];
        bytes.extend_from_slice(b"\r\nnext\rfinal\n");
        let mut lines = Lines::latin1(bytes.as_slice());
        assert_eq!(
            lines.readline(usize::MAX).unwrap(),
            format!("{}\n", "ÿ".repeat(8191))
        );
        assert_eq!(lines.readline(usize::MAX).unwrap(), "next\n");
        assert_eq!(lines.readline(usize::MAX).unwrap(), "final\n");
        assert_eq!(lines.readline(usize::MAX).unwrap(), "");
        let mut utf8 = Lines::new(bytes.as_slice());
        assert!(utf8.readline(usize::MAX).is_err());
    }
    #[test]
    fn bounded_cr_waits_for_the_following_decode() {
        let mut bytes = vec![b'a'; 8191];
        bytes.extend_from_slice(b"\r\xff");
        let mut lines = Lines::new(bytes.as_slice());
        assert_eq!(
            lines.readline(1024 * 1024).unwrap_err().to_string(),
            "UnicodeDecodeError: 'utf-8' codec can't decode byte 0xff in position 0: invalid start byte"
        );
        let mut bytes = vec![b'a'; 8191];
        bytes.push(b'\r');
        let mut lines = Lines::new(bytes.as_slice());
        assert_eq!(
            lines.readline(8192).unwrap(),
            format!("{}\n", "a".repeat(8191))
        );
        assert_eq!(lines.readline(8192).unwrap(), "");
    }
    #[test]
    fn bounded_text_readline_preserves_newline_and_unicode_character_limit() {
        let mut lines = Lines::new("éé\r\nx\r\ny".as_bytes());
        assert_eq!(lines.readline(2).unwrap(), "éé");
        assert_eq!(lines.readline(2).unwrap(), "\n");
        assert_eq!(lines.readline(2).unwrap(), "x\n");
        assert_eq!(lines.readline(2).unwrap(), "y");
        assert_eq!(lines.readline(2).unwrap(), "");
        let text = format!("{}\r\nend", "a".repeat(8191));
        let mut lines = Lines::new(text.as_bytes());
        assert_eq!(
            lines.readline(8192).unwrap(),
            format!("{}\n", "a".repeat(8191))
        );
        assert_eq!(lines.readline(8192).unwrap(), "end");
    }
    fn collect(bytes: &[u8]) -> Result<Vec<String>> {
        let mut lines = Lines::new(bytes);
        let mut result = Vec::new();
        while let Some(line) = lines.next()? {
            result.push(line);
        }
        Ok(result)
    }
    #[test]
    fn cr_lf_and_crlf_are_distinct_from_empty_eof() {
        assert_eq!(
            collect(b"a\rb\r\nc\n\nlast").unwrap(),
            ["a", "b", "c", "", "last"]
        );
        assert_eq!(collect(b"\r\n").unwrap(), [""]);
        assert!(collect(b"").unwrap().is_empty());
    }
    #[test]
    fn crlf_can_cross_the_eight_kib_decode_boundary() {
        let text = format!("{}\r\nlast", "a".repeat(8191));
        assert_eq!(
            collect(text.as_bytes()).unwrap(),
            ["a".repeat(8191), "last".into()]
        );
    }
    #[test]
    fn multibyte_character_can_cross_the_eight_kib_decode_boundary() {
        let text = format!("{}α\nlast", "a".repeat(8191));
        assert_eq!(
            collect(text.as_bytes()).unwrap(),
            [format!("{}α", "a".repeat(8191)), "last".into()]
        );
    }
    #[test]
    fn invalid_decoded_block_fails_before_its_first_line_and_partial_eof_fails() {
        let mut lines = Lines::new(&b"first\n\xff"[..]);
        assert_eq!(
            lines.next().unwrap_err().to_string(),
            "UnicodeDecodeError: input is not UTF-8"
        );
        assert!(collect(b"partial\xe2\x82").is_err());
        let mut bytes = vec![b'a'; 8191];
        bytes.push(b'\n');
        bytes.push(0xff);
        let mut lines = Lines::new(bytes.as_slice());
        assert_eq!(lines.next().unwrap(), Some("a".repeat(8191)));
        assert!(lines.next().is_err());
    }
    #[test]
    fn csv_lines_keep_endings_and_strip_only_initial_bom() {
        let mut lines = Lines::new("\u{feff}a\r\nb\rc\n\u{feff}last".as_bytes());
        let mut first = true;
        let mut actual = vec![];
        while let Some(line) = lines.next_csv(&mut first).unwrap() {
            actual.push(line);
        }
        assert_eq!(actual, ["a\r\n", "b\r", "c\n", "\u{feff}last"]);
    }
    #[test]
    fn csv_crlf_and_multibyte_cross_decode_blocks() {
        for offset in [8188, 8191] {
            let text = format!("{}α\r\nlast", "a".repeat(offset));
            let mut lines = Lines::new(text.as_bytes());
            let mut first = true;
            assert_eq!(
                lines.next_csv(&mut first).unwrap(),
                Some(format!("{}α\r\n", "a".repeat(offset)))
            );
            assert_eq!(lines.next_csv(&mut first).unwrap(), Some("last".into()));
        }
    }
}
