//! JSON reading and writing that behaves like Python's `json` module for the
//! receipts and manifests the SpliceAI scripts handle.
//!
//! Reading follows `json.loads`: whitespace is space, tab, line feed and
//! carriage return; `NaN`, `Infinity` and `-Infinity` are accepted; control
//! characters inside strings are refused; an object keeps the position of a
//! key that repeats and the last value wins, unless duplicates are refused as
//! the scripts' `_unique_object` hook does. Writing follows `json.dumps`
//! (`", "` and `": "` separators, or one item per line with `indent=2`,
//! `ensure_ascii` on by default, keys in insertion order or sorted).
//! Numbers with a fraction or exponent are read but cannot be written back,
//! because Python's float printing is not reproduced; the receipts hold none.
use std::fmt::Write as _;

#[derive(Debug, Clone)]
pub(crate) enum Json {
    Null,
    Bool(bool),
    Int(i128),
    /// The value and the text it was written as.
    Float(f64, String),
    Str(String),
    Array(Vec<Json>),
    Object(Vec<(String, Json)>),
}

const MAX_DEPTH: usize = 400;

impl Json {
    pub(crate) fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Object(items) => items.iter().find(|(name, _)| name == key).map(|(_, v)| v),
            _ => None,
        }
    }

    pub(crate) fn as_str(&self) -> Option<&str> {
        match self {
            Json::Str(text) => Some(text),
            _ => None,
        }
    }

    pub(crate) fn as_object(&self) -> Option<&Vec<(String, Json)>> {
        match self {
            Json::Object(items) => Some(items),
            _ => None,
        }
    }

    /// Sets a key: an existing key keeps its position, a new one goes last.
    pub(crate) fn set(&mut self, key: &str, value: Json) {
        if let Json::Object(items) = self {
            if let Some(slot) = items.iter_mut().find(|(name, _)| name == key) {
                slot.1 = value;
            } else {
                items.push((key.to_owned(), value));
            }
        }
    }
}

/// Python's `==` on parsed JSON: `True == 1`, `1 == 1.0`, objects ignore order.
pub(crate) fn py_equal(left: &Json, right: &Json) -> bool {
    fn number(value: &Json) -> Option<f64> {
        match value {
            Json::Bool(flag) => Some(if *flag { 1.0 } else { 0.0 }),
            Json::Int(int) => Some(*int as f64),
            Json::Float(float, _) => Some(*float),
            _ => None,
        }
    }
    match (left, right) {
        (Json::Null, Json::Null) => true,
        (Json::Str(a), Json::Str(b)) => a == b,
        (Json::Int(a), Json::Int(b)) => a == b,
        (Json::Array(a), Json::Array(b)) => {
            a.len() == b.len() && a.iter().zip(b.iter()).all(|(x, y)| py_equal(x, y))
        }
        (Json::Object(a), Json::Object(b)) => {
            a.len() == b.len()
                && a.iter().all(|(key, value)| match right.get(key) {
                    Some(other) => py_equal(value, other),
                    None => false,
                })
        }
        _ => match (number(left), number(right)) {
            (Some(a), Some(b)) => a == b,
            _ => false,
        },
    }
}

struct Parser<'a> {
    bytes: &'a [u8],
    text: &'a str,
    pos: usize,
    reject_duplicates: bool,
}

/// Parses a whole document. With `reject_duplicate_keys`, a key repeated in one
/// object is an error naming the key (`Duplicate manifest key: ...`).
pub(crate) fn parse(text: &str, reject_duplicate_keys: bool) -> Result<Json, String> {
    if text.starts_with('\u{feff}') {
        return Err("Unexpected UTF-8 BOM (decode using utf-8-sig)".to_owned());
    }
    let mut parser = Parser {
        bytes: text.as_bytes(),
        text,
        pos: 0,
        reject_duplicates: reject_duplicate_keys,
    };
    parser.skip_whitespace();
    let value = parser.value(0)?;
    parser.skip_whitespace();
    if parser.pos != parser.bytes.len() {
        return Err(format!("Extra data: char {}", parser.pos));
    }
    Ok(value)
}

impl<'a> Parser<'a> {
    fn skip_whitespace(&mut self) {
        while self.pos < self.bytes.len() && matches!(self.bytes[self.pos], b' ' | b'\t' | b'\n' | b'\r')
        {
            self.pos += 1;
        }
    }

    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.pos).copied()
    }

    fn expect_word(&mut self, word: &str, value: Json) -> Result<Json, String> {
        if self.text[self.pos..].starts_with(word) {
            self.pos += word.len();
            Ok(value)
        } else {
            Err(format!("Expecting value: char {}", self.pos))
        }
    }

    fn value(&mut self, depth: usize) -> Result<Json, String> {
        if depth > MAX_DEPTH {
            return Err("JSON nested too deeply".to_owned());
        }
        match self.peek() {
            None => Err(format!("Expecting value: char {}", self.pos)),
            Some(b'"') => Ok(Json::Str(self.string()?)),
            Some(b'{') => self.object(depth),
            Some(b'[') => self.array(depth),
            Some(b'n') => self.expect_word("null", Json::Null),
            Some(b't') => self.expect_word("true", Json::Bool(true)),
            Some(b'f') => self.expect_word("false", Json::Bool(false)),
            Some(b'N') => self.expect_word("NaN", Json::Float(f64::NAN, "NaN".to_owned())),
            Some(b'I') => self.expect_word(
                "Infinity",
                Json::Float(f64::INFINITY, "Infinity".to_owned()),
            ),
            Some(b'-') if self.text[self.pos..].starts_with("-Infinity") => self.expect_word(
                "-Infinity",
                Json::Float(f64::NEG_INFINITY, "-Infinity".to_owned()),
            ),
            Some(b'-') | Some(b'0'..=b'9') => self.number(),
            Some(_) => Err(format!("Expecting value: char {}", self.pos)),
        }
    }

    fn digits(&mut self) -> usize {
        let start = self.pos;
        while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_digit() {
            self.pos += 1;
        }
        self.pos - start
    }

    fn number(&mut self) -> Result<Json, String> {
        let start = self.pos;
        if self.peek() == Some(b'-') {
            self.pos += 1;
        }
        match self.peek() {
            Some(b'0') => self.pos += 1,
            Some(b'1'..=b'9') => {
                self.digits();
            }
            _ => {
                self.pos = start;
                return Err(format!("Expecting value: char {start}"));
            }
        }
        let mut is_float = false;
        if self.peek() == Some(b'.') {
            let save = self.pos;
            self.pos += 1;
            if self.digits() == 0 {
                self.pos = save;
            } else {
                is_float = true;
            }
        }
        if matches!(self.peek(), Some(b'e') | Some(b'E')) {
            let save = self.pos;
            self.pos += 1;
            if matches!(self.peek(), Some(b'+') | Some(b'-')) {
                self.pos += 1;
            }
            if self.digits() == 0 {
                self.pos = save;
            } else {
                is_float = true;
            }
        }
        let literal = &self.text[start..self.pos];
        if is_float {
            let value: f64 = literal
                .parse()
                .map_err(|_| format!("Expecting value: char {start}"))?;
            Ok(Json::Float(value, literal.to_owned()))
        } else {
            let value: i128 = literal
                .parse()
                .map_err(|_| format!("integer too large for this tool: {literal}"))?;
            Ok(Json::Int(value))
        }
    }

    fn hex4(&mut self) -> Result<u32, String> {
        let end = self.pos + 4;
        if end > self.bytes.len() {
            return Err(format!("Invalid \\uXXXX escape: char {}", self.pos));
        }
        let slice = &self.bytes[self.pos..end];
        if !slice.iter().all(|b| b.is_ascii_hexdigit()) {
            return Err(format!("Invalid \\uXXXX escape: char {}", self.pos));
        }
        let digits = std::str::from_utf8(slice).map_err(|_| "Invalid \\uXXXX escape".to_owned())?;
        self.pos = end;
        u32::from_str_radix(digits, 16).map_err(|_| "Invalid \\uXXXX escape".to_owned())
    }

    fn string(&mut self) -> Result<String, String> {
        // self.pos is at the opening quote.
        self.pos += 1;
        let mut out = String::new();
        loop {
            let rest = &self.text[self.pos..];
            let mut chars = rest.chars();
            let c = match chars.next() {
                Some(c) => c,
                None => return Err(format!("Unterminated string starting at: char {}", self.pos)),
            };
            self.pos += c.len_utf8();
            match c {
                '"' => return Ok(out),
                '\\' => {
                    let escape = match self.peek() {
                        Some(b) => b,
                        None => return Err("Unterminated string".to_owned()),
                    };
                    self.pos += 1;
                    match escape {
                        b'"' => out.push('"'),
                        b'\\' => out.push('\\'),
                        b'/' => out.push('/'),
                        b'b' => out.push('\u{8}'),
                        b'f' => out.push('\u{c}'),
                        b'n' => out.push('\n'),
                        b'r' => out.push('\r'),
                        b't' => out.push('\t'),
                        b'u' => {
                            let first = self.hex4()?;
                            let code = if (0xd800..0xdc00).contains(&first) {
                                if self.text[self.pos..].starts_with("\\u") {
                                    self.pos += 2;
                                    let second = self.hex4()?;
                                    if !(0xdc00..0xe000).contains(&second) {
                                        return Err("lone surrogate in a string is not supported".to_owned());
                                    }
                                    0x10000 + ((first - 0xd800) << 10) + (second - 0xdc00)
                                } else {
                                    return Err("lone surrogate in a string is not supported".to_owned());
                                }
                            } else if (0xdc00..0xe000).contains(&first) {
                                return Err("lone surrogate in a string is not supported".to_owned());
                            } else {
                                first
                            };
                            match char::from_u32(code) {
                                Some(c) => out.push(c),
                                None => return Err("Invalid \\uXXXX escape".to_owned()),
                            }
                        }
                        _ => return Err(format!("Invalid \\escape: char {}", self.pos - 2)),
                    }
                }
                c if (c as u32) < 0x20 => {
                    return Err(format!("Invalid control character at: char {}", self.pos - 1));
                }
                c => out.push(c),
            }
        }
    }

    fn array(&mut self, depth: usize) -> Result<Json, String> {
        self.pos += 1;
        let mut items = Vec::new();
        self.skip_whitespace();
        if self.peek() == Some(b']') {
            self.pos += 1;
            return Ok(Json::Array(items));
        }
        loop {
            self.skip_whitespace();
            items.push(self.value(depth + 1)?);
            self.skip_whitespace();
            match self.peek() {
                Some(b',') => self.pos += 1,
                Some(b']') => {
                    self.pos += 1;
                    return Ok(Json::Array(items));
                }
                _ => return Err(format!("Expecting ',' delimiter: char {}", self.pos)),
            }
        }
    }

    fn object(&mut self, depth: usize) -> Result<Json, String> {
        self.pos += 1;
        let mut items: Vec<(String, Json)> = Vec::new();
        self.skip_whitespace();
        if self.peek() == Some(b'}') {
            self.pos += 1;
            return Ok(Json::Object(items));
        }
        loop {
            self.skip_whitespace();
            if self.peek() != Some(b'"') {
                return Err(format!(
                    "Expecting property name enclosed in double quotes: char {}",
                    self.pos
                ));
            }
            let key = self.string()?;
            self.skip_whitespace();
            if self.peek() != Some(b':') {
                return Err(format!("Expecting ':' delimiter: char {}", self.pos));
            }
            self.pos += 1;
            self.skip_whitespace();
            let value = self.value(depth + 1)?;
            if let Some(slot) = items.iter_mut().find(|(name, _)| *name == key) {
                if self.reject_duplicates {
                    return Err(format!("Duplicate manifest key: {key}"));
                }
                slot.1 = value;
            } else {
                items.push((key, value));
            }
            self.skip_whitespace();
            match self.peek() {
                Some(b',') => self.pos += 1,
                Some(b'}') => {
                    self.pos += 1;
                    return Ok(Json::Object(items));
                }
                _ => return Err(format!("Expecting ',' delimiter: char {}", self.pos)),
            }
        }
    }
}

/// A JSON string literal as `json.dumps` writes it.
pub(crate) fn quote(text: &str, ensure_ascii: bool) -> String {
    let mut out = String::with_capacity(text.len() + 2);
    out.push('"');
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{8}' => out.push_str("\\b"),
            '\u{c}' => out.push_str("\\f"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c if ensure_ascii && (c as u32) > 0x7e => {
                let code = c as u32;
                if code > 0xffff {
                    let reduced = code - 0x10000;
                    let _ = write!(
                        out,
                        "\\u{:04x}\\u{:04x}",
                        0xd800 + (reduced >> 10),
                        0xdc00 + (reduced & 0x3ff)
                    );
                } else {
                    let _ = write!(out, "\\u{:04x}", code);
                }
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn write_value(
    value: &Json,
    indent: Option<usize>,
    level: usize,
    sort_keys: bool,
    ensure_ascii: bool,
    out: &mut String,
) -> Result<(), String> {
    match value {
        Json::Null => out.push_str("null"),
        Json::Bool(true) => out.push_str("true"),
        Json::Bool(false) => out.push_str("false"),
        Json::Int(int) => {
            let _ = write!(out, "{int}");
        }
        Json::Float(_, text) => {
            return Err(format!(
                "a JSON number with a fraction or exponent ({text}) cannot be written back by this tool"
            ));
        }
        Json::Str(text) => out.push_str(&quote(text, ensure_ascii)),
        Json::Array(items) => {
            if items.is_empty() {
                out.push_str("[]");
                return Ok(());
            }
            out.push('[');
            for (number, item) in items.iter().enumerate() {
                if number > 0 {
                    out.push_str(if indent.is_some() { "," } else { ", " });
                }
                if let Some(width) = indent {
                    out.push('\n');
                    out.push_str(&" ".repeat(width * (level + 1)));
                }
                write_value(item, indent, level + 1, sort_keys, ensure_ascii, out)?;
            }
            if let Some(width) = indent {
                out.push('\n');
                out.push_str(&" ".repeat(width * level));
            }
            out.push(']');
        }
        Json::Object(items) => {
            if items.is_empty() {
                out.push_str("{}");
                return Ok(());
            }
            let mut order: Vec<&(String, Json)> = items.iter().collect();
            if sort_keys {
                order.sort_by(|a, b| a.0.cmp(&b.0));
            }
            out.push('{');
            for (number, (key, item)) in order.into_iter().enumerate() {
                if number > 0 {
                    out.push_str(if indent.is_some() { "," } else { ", " });
                }
                if let Some(width) = indent {
                    out.push('\n');
                    out.push_str(&" ".repeat(width * (level + 1)));
                }
                out.push_str(&quote(key, ensure_ascii));
                out.push_str(": ");
                write_value(item, indent, level + 1, sort_keys, ensure_ascii, out)?;
            }
            if let Some(width) = indent {
                out.push('\n');
                out.push_str(&" ".repeat(width * level));
            }
            out.push('}');
        }
    }
    Ok(())
}

/// `json.dumps(value)`, optionally with `sort_keys=True`.
pub(crate) fn dumps(value: &Json, sort_keys: bool) -> Result<String, String> {
    let mut out = String::new();
    write_value(value, None, 0, sort_keys, true, &mut out)?;
    Ok(out)
}

/// `json.dumps(value, indent=2)`.
pub(crate) fn dumps_indented(value: &Json) -> Result<String, String> {
    let mut out = String::new();
    write_value(value, Some(2), 0, false, true, &mut out)?;
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn object(items: Vec<(&str, Json)>) -> Json {
        Json::Object(items.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
    }

    #[test]
    fn documents_round_trip_in_python_layout() {
        let value = object(vec![
            ("b", Json::Int(2)),
            ("a", Json::Array(vec![Json::Bool(true), Json::Null, Json::Str("x\"y".into())])),
            ("empty", object(vec![])),
        ]);
        assert_eq!(
            dumps(&value, false).unwrap(),
            r#"{"b": 2, "a": [true, null, "x\"y"], "empty": {}}"#
        );
        assert_eq!(
            dumps(&value, true).unwrap(),
            r#"{"a": [true, null, "x\"y"], "b": 2, "empty": {}}"#
        );
        assert_eq!(
            dumps_indented(&value).unwrap(),
            "{\n  \"b\": 2,\n  \"a\": [\n    true,\n    null,\n    \"x\\\"y\"\n  ],\n  \"empty\": {}\n}"
        );
        let parsed = parse(&dumps_indented(&value).unwrap(), true).unwrap();
        assert!(py_equal(&parsed, &value));
    }

    #[test]
    fn strings_are_escaped_like_json_dumps() {
        assert_eq!(quote("caf\u{e9}", true), "\"caf\\u00e9\"");
        assert_eq!(quote("caf\u{e9}", false), "\"caf\u{e9}\"");
        assert_eq!(quote("\u{1f600}", true), "\"\\ud83d\\ude00\"");
        assert_eq!(quote("a\tb\u{1}\u{7f}", true), "\"a\\tb\\u0001\\u007f\"");
        assert_eq!(quote("\u{7f}", false), "\"\u{7f}\"");
    }

    #[test]
    fn python_json_quirks_are_followed() {
        assert!(matches!(parse("NaN", false), Ok(Json::Float(_, _))));
        assert!(parse("[1,]", false).is_err());
        assert!(parse("{\"a\": 1,}", false).is_err());
        assert!(parse("01", false).is_err());
        assert!(parse("\"a\nb\"", false).is_err());
        assert!(parse("{} x", false).is_err());
        assert!(parse("\u{feff}{}", false).is_err());
        assert!(matches!(parse(" [1, 2.5e3, -0] ", false), Ok(Json::Array(_))));
        assert!(dumps(&parse("[1.5]", false).unwrap(), false).is_err());
    }

    #[test]
    fn repeated_keys_follow_the_hook_or_python_dict_rules() {
        let error = parse("{\"a\": 1, \"b\": {\"c\": 1, \"c\": 2}}", true).unwrap_err();
        assert_eq!(error, "Duplicate manifest key: c");
        let kept = parse("{\"a\": 1, \"b\": 2, \"a\": 3}", false).unwrap();
        assert_eq!(dumps(&kept, false).unwrap(), "{\"a\": 3, \"b\": 2}");
    }

    #[test]
    fn equality_treats_true_as_one_and_ignores_object_order() {
        let a = parse("{\"x\": true, \"y\": [1, 2]}", false).unwrap();
        let b = parse("{\"y\": [1.0, 2], \"x\": 1}", false).unwrap();
        assert!(py_equal(&a, &b));
        let c = parse("{\"x\": true, \"y\": [1, 3]}", false).unwrap();
        assert!(!py_equal(&a, &c));
    }

    #[test]
    fn surrogate_pairs_in_escapes_become_one_character() {
        let value = parse("\"\\ud83d\\ude00\"", false).unwrap();
        assert_eq!(value.as_str(), Some("\u{1f600}"));
        assert!(parse("\"\\ud83d\"", false).is_err());
    }
}
