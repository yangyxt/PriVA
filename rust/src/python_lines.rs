//! Shared Python TextIOWrapper line boundaries, without source-specific policy.
use crate::cli::{Result, invalid};
use std::io::Read;

/// TextIOWrapper's UTF-8 decoding and universal CR/LF/CRLF line boundaries.
/// Decode each input block before returning its lines, so invalid UTF-8 does
/// not silently publish earlier lines from the same decoded block.
pub(crate) struct Lines<R: Read> {
    source: R,
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
            text: String::new(),
            at: 0,
            pending: Vec::new(),
            eof: false,
            skip_lf: false,
        }
    }
    fn fill(&mut self) -> Result<bool> {
        if self.eof {
            return Ok(false);
        }
        let mut block = [0_u8; 8192];
        let count = self.source.read(&mut block)?;
        self.eof = count == 0;
        self.pending.extend_from_slice(&block[..count]);
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
}

#[cfg(test)]
mod tests {
    use super::*;
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
}
