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
//! Numbers with a fraction or exponent are written the way `float.__repr__`
//! prints them (shortest digits that read back the same, exponent form below
//! 1e-4 and from 1e16 up); `NaN`, `Infinity` and `-Infinity` as `json.dumps` does.
//! Parse errors carry the text Python's `JSONDecodeError` shows:
//! `<message>: line L column C (char N)`, with N counted in characters.
use std::fmt::Write as _;

#[derive(Debug, Clone)]
pub(crate) enum Json {
    Null,
    Bool(bool),
    Int(i128),
    Float(f64),
    Str(String),
    Array(Vec<Json>),
    Object(Vec<(String, Json)>),
}

const MAX_DEPTH: usize = 400;
const LONE_SURROGATE: &str = "lone surrogate in a string is not supported";
const TOO_DEEP: &str = "JSON nested too deeply";

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
            Json::Float(float) => Some(*float),
            _ => None,
        }
    }
    match (left, right) {
        (Json::Null, Json::Null) => true,
        (Json::Str(a), Json::Str(b)) => a == b,
        (Json::Int(a), Json::Int(b)) => a == b,
        // Every NaN `json.loads` returns is the one object `NaN`, and containers
        // compare identical items as equal.
        (Json::Float(a), Json::Float(b)) => a == b || (a.is_nan() && b.is_nan()),
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

/// Whether a parse message is this tool's refusal of something Python reads
/// (an integer beyond 128 bits, a lone surrogate) or fails on in its own way
/// (nesting too deep for the interpreter), as opposed to a `JSONDecodeError`.
pub(crate) fn is_refusal(message: &str) -> bool {
    message.starts_with("integer too large for this tool")
        || message == LONE_SURROGATE
        || message == TOO_DEEP
}

/// A parse failure: Python's message and, when Python shows one, the byte
/// offset (in the UTF-8 text) the error points at.
struct Problem {
    message: String,
    position: Option<usize>,
}

fn problem<T>(message: &str, position: usize) -> Result<T, Problem> {
    Err(Problem {
        message: message.to_owned(),
        position: Some(position),
    })
}

fn plain_problem<T>(message: &str) -> Result<T, Problem> {
    Err(Problem {
        message: message.to_owned(),
        position: None,
    })
}

impl Problem {
    /// `JSONDecodeError.__str__`: line and column count from 1, the offset
    /// counts characters, and the column is the distance from the last line break.
    fn describe(&self, text: &str) -> String {
        let Some(byte) = self.position else {
            return self.message.clone();
        };
        let before = &text[..byte.min(text.len())];
        let offset = before.chars().count();
        let line = before.matches('\n').count() + 1;
        let column = match before.rfind('\n') {
            Some(newline) => before[newline + 1..].chars().count() + 1,
            None => offset + 1,
        };
        format!(
            "{}: line {line} column {column} (char {offset})",
            self.message
        )
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
    parse_document(text, reject_duplicate_keys).map_err(|failure| failure.describe(text))
}

fn parse_document(text: &str, reject_duplicate_keys: bool) -> Result<Json, Problem> {
    if text.starts_with('\u{feff}') {
        return problem("Unexpected UTF-8 BOM (decode using utf-8-sig)", 0);
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
        return problem("Extra data", parser.pos);
    }
    Ok(value)
}

impl<'a> Parser<'a> {
    fn skip_whitespace(&mut self) {
        while self.pos < self.bytes.len()
            && matches!(self.bytes[self.pos], b' ' | b'\t' | b'\n' | b'\r')
        {
            self.pos += 1;
        }
    }

    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.pos).copied()
    }

    fn expect_word(&mut self, word: &str, value: Json) -> Result<Json, Problem> {
        if self.text[self.pos..].starts_with(word) {
            self.pos += word.len();
            Ok(value)
        } else {
            problem("Expecting value", self.pos)
        }
    }

    fn value(&mut self, depth: usize) -> Result<Json, Problem> {
        if depth > MAX_DEPTH {
            return plain_problem(TOO_DEEP);
        }
        match self.peek() {
            None => problem("Expecting value", self.pos),
            Some(b'"') => Ok(Json::Str(self.string()?)),
            Some(b'{') => self.object(depth),
            Some(b'[') => self.array(depth),
            Some(b'n') => self.expect_word("null", Json::Null),
            Some(b't') => self.expect_word("true", Json::Bool(true)),
            Some(b'f') => self.expect_word("false", Json::Bool(false)),
            Some(b'N') => self.expect_word("NaN", Json::Float(f64::NAN)),
            Some(b'I') => self.expect_word("Infinity", Json::Float(f64::INFINITY)),
            Some(b'-') if self.text[self.pos..].starts_with("-Infinity") => {
                self.expect_word("-Infinity", Json::Float(f64::NEG_INFINITY))
            }
            Some(b'-') | Some(b'0'..=b'9') => self.number(),
            Some(_) => problem("Expecting value", self.pos),
        }
    }

    fn digits(&mut self) -> usize {
        let start = self.pos;
        while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_digit() {
            self.pos += 1;
        }
        self.pos - start
    }

    fn number(&mut self) -> Result<Json, Problem> {
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
                return problem("Expecting value", start);
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
            let value: f64 = match literal.parse() {
                Ok(value) => value,
                Err(_) => return problem("Expecting value", start),
            };
            Ok(Json::Float(value))
        } else {
            // Python integers have no size limit; this tool refuses what i128 cannot hold.
            let value: i128 = match literal.parse() {
                Ok(value) => value,
                Err(_) => {
                    return plain_problem(&format!("integer too large for this tool: {literal}"));
                }
            };
            Ok(Json::Int(value))
        }
    }

    /// The four hex digits after a `\\u` at `u_position`. Python points its error at
    /// the `u`, and also refuses digits that end the text (no closing quote can follow).
    fn hex4(&mut self, u_position: usize) -> Result<u32, Problem> {
        let end = self.pos + 4;
        if end >= self.bytes.len() {
            return problem("Invalid \\uXXXX escape", u_position);
        }
        let slice = &self.bytes[self.pos..end];
        if !slice.iter().all(|b| b.is_ascii_hexdigit()) {
            return problem("Invalid \\uXXXX escape", u_position);
        }
        let Ok(digits) = std::str::from_utf8(slice) else {
            return problem("Invalid \\uXXXX escape", u_position);
        };
        self.pos = end;
        match u32::from_str_radix(digits, 16) {
            Ok(code) => Ok(code),
            Err(_) => problem("Invalid \\uXXXX escape", u_position),
        }
    }

    fn string(&mut self) -> Result<String, Problem> {
        // self.pos is at the opening quote; Python reports an unfinished string there.
        let begin = self.pos;
        self.pos += 1;
        let mut out = String::new();
        loop {
            let rest = &self.text[self.pos..];
            let mut chars = rest.chars();
            let c = match chars.next() {
                Some(c) => c,
                None => return problem("Unterminated string starting at", begin),
            };
            self.pos += c.len_utf8();
            match c {
                '"' => return Ok(out),
                '\\' => {
                    let escape = match self.peek() {
                        Some(b) => b,
                        None => return problem("Unterminated string starting at", begin),
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
                            let first = self.hex4(self.pos - 1)?;
                            let code = if (0xd800..0xdc00).contains(&first) {
                                if self.text[self.pos..].starts_with("\\u") {
                                    self.pos += 2;
                                    let second = self.hex4(self.pos - 1)?;
                                    if !(0xdc00..0xe000).contains(&second) {
                                        return plain_problem(LONE_SURROGATE);
                                    }
                                    0x10000 + ((first - 0xd800) << 10) + (second - 0xdc00)
                                } else {
                                    return plain_problem(LONE_SURROGATE);
                                }
                            } else if (0xdc00..0xe000).contains(&first) {
                                return plain_problem(LONE_SURROGATE);
                            } else {
                                first
                            };
                            match char::from_u32(code) {
                                Some(c) => out.push(c),
                                None => return problem("Invalid \\uXXXX escape", self.pos - 1),
                            }
                        }
                        _ => return problem("Invalid \\escape", self.pos - 2),
                    }
                }
                c if (c as u32) < 0x20 => {
                    return problem("Invalid control character at", self.pos - 1);
                }
                c => out.push(c),
            }
        }
    }

    fn array(&mut self, depth: usize) -> Result<Json, Problem> {
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
                _ => return problem("Expecting ',' delimiter", self.pos),
            }
        }
    }

    fn object(&mut self, depth: usize) -> Result<Json, Problem> {
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
                return problem(
                    "Expecting property name enclosed in double quotes",
                    self.pos,
                );
            }
            let key = self.string()?;
            self.skip_whitespace();
            if self.peek() != Some(b':') {
                return problem("Expecting ':' delimiter", self.pos);
            }
            self.pos += 1;
            self.skip_whitespace();
            let value = self.value(depth + 1)?;
            if let Some(slot) = items.iter_mut().find(|(name, _)| *name == key) {
                if self.reject_duplicates {
                    return plain_problem(&format!("Duplicate manifest key: {key}"));
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
                _ => return problem("Expecting ',' delimiter", self.pos),
            }
        }
    }
}

/// `float.__repr__`: the shortest digits that read back as the same number,
/// written plainly for decimal exponents from -4 to 15 and as `1e-05` or
/// `1.5e+16` outside that range, with `.0` after a whole number.
pub(crate) fn python_float_repr(value: f64) -> String {
    if value.is_nan() {
        return "nan".to_owned();
    }
    if value.is_infinite() {
        return if value < 0.0 { "-inf" } else { "inf" }.to_owned();
    }
    let sign = if value.is_sign_negative() { "-" } else { "" };
    // Rust prints the shortest round-trip digits as `d.ddde<exponent>`.
    let scientific = format!("{:e}", value.abs());
    let (mantissa, exponent) = scientific.split_once('e').unwrap_or((&scientific[..], "0"));
    let exponent: i32 = exponent.parse().unwrap_or(0);
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    // Python's `decpt`: the value is 0.<digits> times ten to this power.
    let point = exponent + 1;
    let length = digits.len() as i32;
    let body = if point <= -4 || point > 16 {
        let fraction = if length > 1 {
            format!(".{}", &digits[1..])
        } else {
            String::new()
        };
        let power = point - 1;
        let power_sign = if power < 0 { '-' } else { '+' };
        format!(
            "{}{}e{}{:02}",
            &digits[..1],
            fraction,
            power_sign,
            power.abs()
        )
    } else if point <= 0 {
        format!("0.{}{}", "0".repeat((-point) as usize), digits)
    } else if point >= length {
        format!("{}{}.0", digits, "0".repeat((point - length) as usize))
    } else {
        format!(
            "{}.{}",
            &digits[..point as usize],
            &digits[point as usize..]
        )
    };
    format!("{sign}{body}")
}

/// A float as `json.dumps` writes it: `repr`, except `NaN` and the infinities.
fn float_text(value: f64) -> String {
    if value.is_nan() {
        "NaN".to_owned()
    } else if value.is_infinite() {
        if value < 0.0 { "-Infinity" } else { "Infinity" }.to_owned()
    } else {
        python_float_repr(value)
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
        Json::Float(float) => out.push_str(&float_text(*float)),
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
            (
                "a",
                Json::Array(vec![Json::Bool(true), Json::Null, Json::Str("x\"y".into())]),
            ),
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
        assert!(matches!(parse("NaN", false), Ok(Json::Float(_))));
        assert!(parse("[1,]", false).is_err());
        assert!(parse("{\"a\": 1,}", false).is_err());
        assert!(parse("01", false).is_err());
        assert!(parse("\"a\nb\"", false).is_err());
        assert!(parse("{} x", false).is_err());
        assert!(parse("\u{feff}{}", false).is_err());
        assert!(matches!(
            parse(" [1, 2.5e3, -0] ", false),
            Ok(Json::Array(_))
        ));
        assert_eq!(
            dumps(
                &parse("[1.5, 0.250, 1e-5, 2E+16, -0.0, 7.0]", false).unwrap(),
                false
            )
            .unwrap(),
            "[1.5, 0.25, 1e-05, 2e+16, -0.0, 7.0]"
        );
    }

    #[test]
    fn floats_print_like_python_repr() {
        for (value, text) in [
            (0.0, "0.0"),
            (0.1, "0.1"),
            (0.25, "0.25"),
            (1.0, "1.0"),
            (100.0, "100.0"),
            (123456789.0, "123456789.0"),
            (1e15, "1000000000000000.0"),
            (1e16, "1e+16"),
            (1.5e16, "1.5e+16"),
            (123456789012345680.0, "1.2345678901234568e+17"),
            (0.0001, "0.0001"),
            (0.00001, "1e-05"),
            (0.00012345, "0.00012345"),
            (1.5e-7, "1.5e-07"),
            (1e100, "1e+100"),
            (5e-324, "5e-324"),
            (1.7976931348623157e308, "1.7976931348623157e+308"),
            (0.30000000000000004, "0.30000000000000004"),
            (-2.5, "-2.5"),
        ] {
            assert_eq!(python_float_repr(value), text, "{value:e}");
        }
        assert_eq!(float_text(f64::NAN), "NaN");
        assert_eq!(float_text(f64::NEG_INFINITY), "-Infinity");
        assert_eq!(
            parse("1e400", false).map(|v| dumps(&v, false).unwrap()),
            Ok("Infinity".to_owned())
        );
    }

    #[test]
    fn parse_errors_name_line_column_and_character_like_python() {
        assert_eq!(
            parse("", false).unwrap_err(),
            "Expecting value: line 1 column 1 (char 0)"
        );
        assert_eq!(
            parse("{bad}", false).unwrap_err(),
            "Expecting property name enclosed in double quotes: line 1 column 2 (char 1)"
        );
        assert_eq!(
            parse("{\n  \"a\": \u{e9},\n}", false).unwrap_err(),
            "Expecting value: line 2 column 8 (char 9)"
        );
        assert_eq!(
            parse("\u{feff}{}", false).unwrap_err(),
            "Unexpected UTF-8 BOM (decode using utf-8-sig): line 1 column 1 (char 0)"
        );
        assert_eq!(
            parse("{} x", false).unwrap_err(),
            "Extra data: line 1 column 4 (char 3)"
        );
        assert_eq!(
            parse("[\"abc", false).unwrap_err(),
            "Unterminated string starting at: line 1 column 2 (char 1)"
        );
        assert_eq!(
            parse("\"\\q\"", false).unwrap_err(),
            "Invalid \\escape: line 1 column 2 (char 1)"
        );
        assert_eq!(
            parse("\"\\u12G4\"", false).unwrap_err(),
            "Invalid \\uXXXX escape: line 1 column 3 (char 2)"
        );
        assert_eq!(
            parse("{\"a\" 1}", false).unwrap_err(),
            "Expecting ':' delimiter: line 1 column 6 (char 5)"
        );
        assert_eq!(
            parse("[1 2]", false).unwrap_err(),
            "Expecting ',' delimiter: line 1 column 4 (char 3)"
        );
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
    fn nan_items_are_equal_inside_containers() {
        let a = parse("{\"x\": [NaN, 1.5]}", false).unwrap();
        let b = parse("{\"x\": [NaN, 1.5]}", false).unwrap();
        assert!(py_equal(&a, &b));
    }

    #[test]
    fn surrogate_pairs_in_escapes_become_one_character() {
        let value = parse("\"\\ud83d\\ude00\"", false).unwrap();
        assert_eq!(value.as_str(), Some("\u{1f600}"));
        assert!(parse("\"\\ud83d\"", false).is_err());
    }
}
