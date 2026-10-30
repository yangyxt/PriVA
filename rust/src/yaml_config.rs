//! A reader for the flat YAML configuration files PriVA uses (`config.yaml`).
//!
//! It handles what those files contain and refuses everything else loudly
//! instead of guessing:
//! - a block mapping at column 0 whose keys are plain names (letters, digits,
//!   `_`, `.`, `-`), with `#` comments and blank lines;
//! - values that are plain, single-quoted or double-quoted scalars on one
//!   line, resolved to a null, a boolean, an integer, another typed scalar or
//!   a string the way PyYAML's `safe_load` does;
//! - a value that is a one-line flow list without quotes (`[one, two]`), or a
//!   nested block mapping or block list whose entries are single-line scalars
//!   (kept as text, since callers only need to leave them untouched);
//! - the byte span of each scalar value, so a caller can replace just that
//!   span and keep every comment, quote and line ending around it.
//! Anything else (anchors, tags, block scalars, multi-line scalars, tabs,
//! quoted keys, several documents, deeper nesting) is reported as
//! `Unsupported`; text PyYAML itself would refuse is `Invalid`.
use crate::json_text;

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Scalar {
    Null,
    Bool(bool),
    /// A decimal integer, written as Python prints it.
    Int(String),
    /// Another typed scalar: a float, a timestamp, an integer in another base.
    Other(String),
    Str(String),
}

impl Scalar {
    /// Python's `repr()` of the value `safe_load` would produce.
    pub(crate) fn python_repr(&self) -> String {
        match self {
            Scalar::Null => "None".to_owned(),
            Scalar::Bool(true) => "True".to_owned(),
            Scalar::Bool(false) => "False".to_owned(),
            Scalar::Int(text) | Scalar::Other(text) => text.clone(),
            Scalar::Str(text) => python_str_repr(text),
        }
    }
}

/// Python's `repr()` of a string: single quotes unless the text holds a single
/// quote and no double quote, with backslash escapes for control characters.
pub(crate) fn python_str_repr(text: &str) -> String {
    let quote = if text.contains('\'') && !text.contains('"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::new();
    out.push(quote);
    for c in text.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c == quote => {
                out.push('\\');
                out.push(c);
            }
            c if (c as u32) < 0x20 || c as u32 == 0x7f => {
                out.push_str(&format!("\\x{:02x}", c as u32));
            }
            c => out.push(c),
        }
    }
    out.push(quote);
    out
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum ConfigValue {
    Scalar(Scalar),
    /// A nested block or a flow list, kept as the text it was written as.
    Opaque(String),
}

#[derive(Debug, Clone)]
pub(crate) struct ConfigEntry {
    pub(crate) key: String,
    pub(crate) value: ConfigValue,
    /// Byte span of a scalar value in the text. A null value written as
    /// nothing has an empty span right after the colon.
    pub(crate) value_start: usize,
    pub(crate) value_end: usize,
}

#[derive(Debug, Clone)]
pub(crate) struct ConfigDocument {
    pub(crate) entries: Vec<ConfigEntry>,
}

impl ConfigDocument {
    /// The value `safe_load` returns for a key: the last one written.
    pub(crate) fn get(&self, key: &str) -> Option<&ConfigValue> {
        self.entries
            .iter()
            .rev()
            .find(|entry| entry.key == key)
            .map(|entry| &entry.value)
    }

    /// A scalar value, or `None` when the key is absent or holds a collection.
    pub(crate) fn scalar(&self, key: &str) -> Option<&Scalar> {
        match self.get(key) {
            Some(ConfigValue::Scalar(scalar)) => Some(scalar),
            _ => None,
        }
    }

    pub(crate) fn entry(&self, key: &str) -> Option<&ConfigEntry> {
        self.entries.iter().find(|entry| entry.key == key)
    }

    pub(crate) fn has_duplicate_keys(&self) -> bool {
        self.entries.iter().enumerate().any(|(index, entry)| {
            self.entries[..index]
                .iter()
                .any(|earlier| earlier.key == entry.key)
        })
    }
}

#[derive(Debug)]
pub(crate) enum YamlError {
    /// Text PyYAML would refuse.
    Invalid(String),
    /// Valid YAML this reader does not handle.
    Unsupported(String),
}

impl std::fmt::Display for YamlError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            YamlError::Invalid(message) => write!(formatter, "{message}"),
            YamlError::Unsupported(message) => write!(
                formatter,
                "unsupported YAML construct (this tool reads flat configuration files only): {message}"
            ),
        }
    }
}

fn invalid<T>(message: &str) -> Result<T, YamlError> {
    Err(YamlError::Invalid(message.to_owned()))
}

fn unsupported<T>(message: &str) -> Result<T, YamlError> {
    Err(YamlError::Unsupported(message.to_owned()))
}

/// What a whole file turned out to be.
#[derive(Debug)]
pub(crate) enum Parsed {
    Mapping(ConfigDocument),
    /// No content at all, or only an empty flow collection: `safe_load` gives
    /// `None`, `[]` or `{}`.
    Empty,
    /// A block list: `safe_load` gives a non-empty list.
    List,
}

struct Line<'a> {
    start: usize,
    text: &'a str,
    indent: usize,
}

impl Line<'_> {
    fn is_skippable(&self) -> bool {
        let content = self.text.trim_start_matches(' ');
        content.is_empty() || content.starts_with('#')
    }

    fn is_list_item(&self) -> bool {
        let content = &self.text[self.indent..];
        content == "-" || content.starts_with("- ")
    }
}

fn split_lines(text: &str) -> Vec<Line<'_>> {
    let mut lines = Vec::new();
    let mut start = 0;
    while start < text.len() {
        let end = text[start..]
            .find('\n')
            .map(|offset| start + offset)
            .unwrap_or(text.len());
        let raw = &text[start..end];
        let line_text = raw.strip_suffix('\r').unwrap_or(raw);
        let indent = line_text.len() - line_text.trim_start_matches(' ').len();
        lines.push(Line {
            start,
            text: line_text,
            indent,
        });
        start = end + 1;
    }
    lines
}

fn is_printable(c: char) -> bool {
    let code = c as u32;
    matches!(code, 0x09 | 0x0a | 0x0d | 0x20..=0x7e | 0x85 | 0xa0..=0xd7ff | 0xe000..=0xfffd)
        || code >= 0x10000
}

/// Reads a configuration file's text.
pub(crate) fn parse_config(text: &str) -> Result<Parsed, YamlError> {
    if text.starts_with('\u{feff}') {
        return unsupported("a byte order mark");
    }
    for (index, c) in text.char_indices() {
        if !is_printable(c) {
            return invalid("unacceptable character in the stream");
        }
        if c == '\t' {
            return unsupported("a tab character");
        }
        if matches!(c, '\u{85}' | '\u{2028}' | '\u{2029}') {
            return unsupported("a line break other than LF or CRLF");
        }
        if c == '\r' && !text[index + 1..].starts_with('\n') {
            return unsupported("a bare carriage return");
        }
    }
    let lines = split_lines(text);
    let mut first = 0;
    while first < lines.len() && lines[first].is_skippable() {
        first += 1;
    }
    if first == lines.len() {
        return Ok(Parsed::Empty);
    }
    let head = &lines[first];
    if head.indent == 0 && head.is_list_item() {
        return Ok(Parsed::List);
    }
    if head.indent == 0 {
        let content = head.text.trim_end();
        if content == "[]" || content == "{}" {
            return Ok(Parsed::Empty);
        }
        if content == "---"
            || content.starts_with("--- ")
            || content == "..."
            || content.starts_with('%')
        {
            return unsupported("a document marker or directive");
        }
    }
    let mut entries = Vec::new();
    let mut index = first;
    while index < lines.len() {
        let line = &lines[index];
        if line.is_skippable() {
            index += 1;
            continue;
        }
        if line.indent != 0 {
            return unsupported("a mapping that does not start at column 0");
        }
        let key_len = match key_length(line.text) {
            Some(length) => length,
            None => {
                return if line.text.starts_with('"') || line.text.starts_with('\'') {
                    unsupported("a quoted or complex key")
                } else if line.is_list_item() {
                    invalid("sequence entries are not allowed here")
                } else {
                    unsupported("a line that is not a plain key followed by a colon")
                };
            }
        };
        let key = line.text[..key_len].to_owned();
        let after_colon = line.start + key_len + 1;
        let rest = &line.text[key_len + 1..];
        let spaces = rest.len() - rest.trim_start_matches(' ').len();
        let body = &rest[spaces..];
        let body_offset = after_colon + spaces;
        let mut lookahead = index + 1;
        while lookahead < lines.len() && lines[lookahead].is_skippable() {
            lookahead += 1;
        }
        let next_is_nested = lookahead < lines.len()
            && (lines[lookahead].indent > 0
                || (lines[lookahead].indent == 0 && lines[lookahead].is_list_item()));
        if body.is_empty() || body.starts_with('#') {
            if next_is_nested {
                let (block, end) = nested_block(&lines, index + 1)?;
                entries.push(ConfigEntry {
                    key,
                    value: ConfigValue::Opaque(block),
                    value_start: after_colon,
                    value_end: after_colon,
                });
                index = end;
            } else {
                entries.push(ConfigEntry {
                    key,
                    value: ConfigValue::Scalar(Scalar::Null),
                    value_start: after_colon,
                    value_end: after_colon,
                });
                index += 1;
            }
            continue;
        }
        let (value, length) = inline_value(body)?;
        if next_is_nested {
            if lines[lookahead].indent == 0 {
                return invalid("sequence entries are not allowed here");
            }
            return if body.starts_with('"') || body.starts_with('\'') {
                invalid("expected <block end>, but found a more indented line")
            } else {
                unsupported("a value continued on the next line")
            };
        }
        entries.push(ConfigEntry {
            key,
            value,
            value_start: body_offset,
            value_end: body_offset + length,
        });
        index += 1;
    }
    Ok(Parsed::Mapping(ConfigDocument { entries }))
}

/// Length of a plain key followed by a colon and a space or the line end.
fn key_length(text: &str) -> Option<usize> {
    let bytes = text.as_bytes();
    if bytes.is_empty() || !(bytes[0].is_ascii_alphanumeric() || bytes[0] == b'_') {
        return None;
    }
    let mut length = 1;
    while length < bytes.len()
        && (bytes[length].is_ascii_alphanumeric() || matches!(bytes[length], b'_' | b'.' | b'-'))
    {
        length += 1;
    }
    if length < bytes.len() && bytes[length] == b':' {
        if length + 1 == bytes.len() || bytes[length + 1] == b' ' {
            return Some(length);
        }
    }
    None
}

/// A nested block mapping or list whose entries are single-line scalars.
/// Returns its text and the index of the first line after it.
fn nested_block(lines: &[Line<'_>], from: usize) -> Result<(String, usize), YamlError> {
    let mut index = from;
    while index < lines.len() && lines[index].is_skippable() {
        index += 1;
    }
    let is_list = lines[index].is_list_item();
    let indent = lines[index].indent;
    let mut text = String::new();
    while index < lines.len() {
        let line = &lines[index];
        if line.is_skippable() {
            index += 1;
            continue;
        }
        if line.indent == 0 && !line.is_list_item() {
            break;
        }
        if line.indent != indent || line.is_list_item() != is_list {
            return unsupported("a nested structure that is not one flat block");
        }
        let content = &line.text[indent..];
        let entry = if is_list {
            content[1..].trim_start_matches(' ')
        } else {
            let length = match key_length(content) {
                Some(length) => length,
                None => return unsupported("a nested entry that is not a plain key"),
            };
            content[length + 1..].trim_start_matches(' ')
        };
        if entry.is_empty() || entry.starts_with('#') {
            return unsupported("an empty or deeper nested entry");
        }
        match inline_value(entry)? {
            (ConfigValue::Scalar(_), _) => {}
            _ => return unsupported("a collection nested in a nested block"),
        }
        text.push_str(line.text);
        text.push('\n');
        index += 1;
    }
    Ok((text, index))
}

/// The value written after `key:` on one line, and how many bytes of `body`
/// it covers. What follows must be nothing or a comment.
fn inline_value(body: &str) -> Result<(ConfigValue, usize), YamlError> {
    let first = match body.chars().next() {
        Some(c) => c,
        None => return invalid("missing value"),
    };
    let (value, length) = match first {
        '"' => {
            let (text, length) = scan_double(body)?;
            (ConfigValue::Scalar(Scalar::Str(text)), length)
        }
        '\'' => {
            let (text, length) = scan_single(body)?;
            (ConfigValue::Scalar(Scalar::Str(text)), length)
        }
        '[' | '{' => scan_flow(body)?,
        '|' | '>' => return unsupported("a block scalar"),
        '&' | '*' | '!' => return unsupported("an anchor, alias or tag"),
        '%' | '@' | '`' | ',' | ']' | '}' | '#' => {
            return invalid("found a character that cannot start any token");
        }
        '-' | '?' | ':' => {
            let after = body[1..].chars().next();
            if after.is_none() || after == Some(' ') {
                return invalid("a block indicator is not allowed here");
            }
            plain(body)?
        }
        _ => plain(body)?,
    };
    let tail = body[length..].trim_start_matches(' ');
    if !tail.is_empty() && !tail.starts_with('#') {
        return invalid("expected the end of the line after the value");
    }
    if !tail.is_empty() && body[length..].len() == tail.len() {
        // A comment must be separated from the value by a space.
        return invalid("a comment must be preceded by a space");
    }
    Ok((value, length))
}

fn plain(body: &str) -> Result<(ConfigValue, usize), YamlError> {
    let end = match body.find(" #") {
        Some(position) => position,
        None => body.len(),
    };
    let text = body[..end].trim_end_matches(' ');
    if text.contains(": ") || text.ends_with(':') {
        return invalid("mapping values are not allowed here");
    }
    Ok((ConfigValue::Scalar(resolve_plain(text)), text.len()))
}

fn scan_double(body: &str) -> Result<(String, usize), YamlError> {
    let mut out = String::new();
    let mut chars = body.char_indices().skip(1);
    while let Some((index, c)) = chars.next() {
        match c {
            '"' => return Ok((out, index + 1)),
            '\\' => {
                let escape = match chars.next() {
                    Some((_, escape)) => escape,
                    None => return unsupported("a quoted scalar continued on the next line"),
                };
                match escape {
                    '0' => out.push('\0'),
                    'a' => out.push('\u{7}'),
                    'b' => out.push('\u{8}'),
                    't' | '\t' => out.push('\t'),
                    'n' => out.push('\n'),
                    'v' => out.push('\u{b}'),
                    'f' => out.push('\u{c}'),
                    'r' => out.push('\r'),
                    'e' => out.push('\u{1b}'),
                    ' ' => out.push(' '),
                    '"' => out.push('"'),
                    '/' => out.push('/'),
                    '\\' => out.push('\\'),
                    'N' => out.push('\u{85}'),
                    '_' => out.push('\u{a0}'),
                    'L' => out.push('\u{2028}'),
                    'P' => out.push('\u{2029}'),
                    'x' | 'u' | 'U' => {
                        let count = match escape {
                            'x' => 2,
                            'u' => 4,
                            _ => 8,
                        };
                        let mut code: u32 = 0;
                        for _ in 0..count {
                            let digit = match chars.next() {
                                Some((_, d)) => d.to_digit(16),
                                None => None,
                            };
                            match digit {
                                Some(d) => code = code * 16 + d,
                                None => {
                                    return invalid(
                                        "expected an escape sequence of hexadecimal numbers",
                                    );
                                }
                            }
                        }
                        match char::from_u32(code) {
                            Some(ch) => out.push(ch),
                            None => return invalid("invalid Unicode character escape code"),
                        }
                    }
                    _ => return invalid("found unknown escape character"),
                }
            }
            c => out.push(c),
        }
    }
    unsupported("a quoted scalar continued on the next line")
}

fn scan_single(body: &str) -> Result<(String, usize), YamlError> {
    let mut out = String::new();
    let mut chars = body.char_indices().skip(1).peekable();
    while let Some((index, c)) = chars.next() {
        if c == '\'' {
            if let Some((_, '\'')) = chars.peek() {
                out.push('\'');
                chars.next();
            } else {
                return Ok((out, index + 1));
            }
        } else {
            out.push(c);
        }
    }
    unsupported("a quoted scalar continued on the next line")
}

/// A one-line flow list without quotes or nesting, or an empty `{}`.
fn scan_flow(body: &str) -> Result<(ConfigValue, usize), YamlError> {
    let open = body.chars().next().unwrap_or('[');
    let close = if open == '[' { ']' } else { '}' };
    let end = match body.find(close) {
        Some(position) => position,
        None => {
            return invalid(
                "the flow collection is not closed on its line (multi-line flow collections are not supported)",
            );
        }
    };
    let inner = &body[1..end];
    if inner.contains(['[', '{', '"', '\'', '#']) {
        return unsupported("a nested or quoted flow collection");
    }
    if open == '{' {
        if inner.trim().is_empty() {
            return Ok((ConfigValue::Opaque("{}".to_owned()), end + 1));
        }
        return unsupported("a flow mapping with entries");
    }
    let mut items: Vec<&str> = inner.split(',').collect();
    if items.len() > 1 && items[items.len() - 1].trim().is_empty() {
        items.pop();
    }
    if !(items.len() == 1 && items[0].trim().is_empty()) {
        for item in &items {
            let item = item.trim();
            if item.is_empty() {
                return invalid("expected a flow sequence entry");
            }
            if item.contains(':') || item.starts_with(['&', '*', '!', '|', '>', '%', '@', '`']) {
                return unsupported("a flow sequence entry that is not a plain word");
            }
        }
    }
    Ok((ConfigValue::Opaque(body[..=end].to_owned()), end + 1))
}

fn all_bytes(bytes: &[u8], allowed: fn(u8) -> bool) -> bool {
    bytes.iter().all(|byte| allowed(*byte))
}

fn is_digit_or_underscore(byte: u8) -> bool {
    byte.is_ascii_digit() || byte == b'_'
}

fn without_sign(text: &str) -> &str {
    text.strip_prefix(['+', '-']).unwrap_or(text)
}

fn is_decimal_integer(text: &str) -> bool {
    let rest = without_sign(text).as_bytes();
    rest == b"0"
        || (!rest.is_empty()
            && (b'1'..=b'9').contains(&rest[0])
            && all_bytes(&rest[1..], is_digit_or_underscore))
}

fn is_other_integer(text: &str) -> bool {
    let rest = without_sign(text);
    let bytes = rest.as_bytes();
    if let Some(digits) = rest.strip_prefix("0b") {
        return !digits.is_empty()
            && all_bytes(digits.as_bytes(), |b| matches!(b, b'0' | b'1' | b'_'));
    }
    if let Some(digits) = rest.strip_prefix("0x") {
        return !digits.is_empty()
            && all_bytes(digits.as_bytes(), |b| b.is_ascii_hexdigit() || b == b'_');
    }
    if bytes.len() >= 2 && bytes[0] == b'0' {
        return all_bytes(&bytes[1..], |b| matches!(b, b'0'..=b'7' | b'_'));
    }
    is_sexagesimal_integer(rest)
}

/// `[1-9][0-9_]*(:[0-5]?[0-9])+`
fn is_sexagesimal_integer(text: &str) -> bool {
    let mut parts = text.split(':');
    let head = parts.next().unwrap_or("").as_bytes();
    if head.is_empty()
        || !(b'1'..=b'9').contains(&head[0])
        || !all_bytes(&head[1..], is_digit_or_underscore)
    {
        return false;
    }
    sexagesimal_groups(parts)
}

fn sexagesimal_groups<'a>(parts: impl Iterator<Item = &'a str>) -> bool {
    let mut count = 0;
    for part in parts {
        let bytes = part.as_bytes();
        let fits = match bytes.len() {
            1 => bytes[0].is_ascii_digit(),
            2 => (b'0'..=b'5').contains(&bytes[0]) && bytes[1].is_ascii_digit(),
            _ => false,
        };
        if !fits {
            return false;
        }
        count += 1;
    }
    count > 0
}

fn is_exponent(text: &str) -> bool {
    // [-+][0-9]+ after the e or E
    let bytes = text.as_bytes();
    bytes.len() >= 2 && matches!(bytes[0], b'+' | b'-') && bytes[1..].iter().all(u8::is_ascii_digit)
}

/// `digits? '.' [0-9_]* ([eE][-+][0-9]+)?` with the given rule for the digits before the dot.
fn float_tail(after_dot: &str) -> bool {
    match after_dot.find(['e', 'E']) {
        Some(position) => {
            all_bytes(after_dot[..position].as_bytes(), is_digit_or_underscore)
                && is_exponent(&after_dot[position + 1..])
        }
        None => all_bytes(after_dot.as_bytes(), is_digit_or_underscore),
    }
}

fn is_float(text: &str) -> bool {
    let rest = without_sign(text);
    if rest == ".inf" || rest == ".Inf" || rest == ".INF" {
        return true;
    }
    if text == ".nan" || text == ".NaN" || text == ".NAN" {
        return true;
    }
    // [-+]?[0-9][0-9_]*(:[0-5]?[0-9])+\.[0-9_]*
    if let Some((before, after)) = rest.split_once('.') {
        if before.contains(':') {
            let mut parts = before.split(':');
            let head = parts.next().unwrap_or("").as_bytes();
            if !head.is_empty()
                && head[0].is_ascii_digit()
                && all_bytes(&head[1..], is_digit_or_underscore)
                && sexagesimal_groups(parts)
                && all_bytes(after.as_bytes(), is_digit_or_underscore)
            {
                return true;
            }
        }
        // [-+]?[0-9][0-9_]*\.[0-9_]*([eE][-+][0-9]+)?
        let bytes = before.as_bytes();
        if !bytes.is_empty()
            && bytes[0].is_ascii_digit()
            && all_bytes(&bytes[1..], is_digit_or_underscore)
            && float_tail(after)
        {
            return true;
        }
    }
    // \.[0-9][0-9_]*([eE][-+][0-9]+)?  (no sign)
    if let Some(after) = text.strip_prefix('.') {
        let bytes = after.as_bytes();
        if !bytes.is_empty() && bytes[0].is_ascii_digit() && float_tail(&after[1..]) {
            return true;
        }
    }
    false
}

fn is_timestamp(text: &str) -> bool {
    let bytes = text.as_bytes();
    let digits = |from: usize, count: usize| -> bool {
        bytes.len() >= from + count && bytes[from..from + count].iter().all(u8::is_ascii_digit)
    };
    if !(digits(0, 4) && bytes.get(4) == Some(&b'-')) {
        return false;
    }
    if bytes.len() == 10 && digits(5, 2) && bytes[7] == b'-' && digits(8, 2) {
        return true;
    }
    // YYYY-M(M)-D(D)(T| )H(H):MM:SS(.f)?(zone)?
    let mut position = 5;
    let take_digits = |position: &mut usize, minimum: usize, maximum: usize| -> bool {
        let start = *position;
        while *position < bytes.len()
            && bytes[*position].is_ascii_digit()
            && *position - start < maximum
        {
            *position += 1;
        }
        *position - start >= minimum
    };
    if !take_digits(&mut position, 1, 2) || bytes.get(position) != Some(&b'-') {
        return false;
    }
    position += 1;
    if !take_digits(&mut position, 1, 2) {
        return false;
    }
    match bytes.get(position) {
        Some(b'T') | Some(b't') => position += 1,
        Some(b' ') => {
            while bytes.get(position) == Some(&b' ') {
                position += 1;
            }
        }
        _ => return false,
    }
    if !take_digits(&mut position, 1, 2) || bytes.get(position) != Some(&b':') {
        return false;
    }
    position += 1;
    if !take_digits(&mut position, 2, 2) || bytes.get(position) != Some(&b':') {
        return false;
    }
    position += 1;
    if !take_digits(&mut position, 2, 2) {
        return false;
    }
    if bytes.get(position) == Some(&b'.') {
        position += 1;
        while position < bytes.len() && bytes[position].is_ascii_digit() {
            position += 1;
        }
    }
    let zone = text[position..].trim_start_matches(' ');
    if zone.is_empty() || zone == "Z" {
        return true;
    }
    let zone_bytes = zone.as_bytes();
    if !matches!(zone_bytes[0], b'+' | b'-') {
        return false;
    }
    let (hours, minutes) = match zone[1..].split_once(':') {
        Some((hours, minutes)) => (hours, Some(minutes)),
        None => (&zone[1..], None),
    };
    (hours.len() == 1 || hours.len() == 2)
        && hours.bytes().all(|b| b.is_ascii_digit())
        && minutes.is_none_or(|m| m.len() == 2 && m.bytes().all(|b| b.is_ascii_digit()))
}

/// PyYAML's implicit resolvers for a plain scalar.
pub(crate) fn resolve_plain(text: &str) -> Scalar {
    match text {
        "" | "~" | "null" | "Null" | "NULL" => return Scalar::Null,
        "yes" | "Yes" | "YES" | "true" | "True" | "TRUE" | "on" | "On" | "ON" => {
            return Scalar::Bool(true);
        }
        "no" | "No" | "NO" | "false" | "False" | "FALSE" | "off" | "Off" | "OFF" => {
            return Scalar::Bool(false);
        }
        "<<" | "=" => return Scalar::Other(text.to_owned()),
        _ => {}
    }
    let first = text.as_bytes()[0];
    if !(first.is_ascii_digit() || matches!(first, b'+' | b'-' | b'.')) {
        return Scalar::Str(text.to_owned());
    }
    if is_decimal_integer(text) {
        let digits: String = without_sign(text).chars().filter(|c| *c != '_').collect();
        let negative = text.starts_with('-');
        return match digits.parse::<i128>() {
            Ok(0) => Scalar::Int("0".to_owned()),
            Ok(value) => Scalar::Int(if negative {
                format!("-{value}")
            } else {
                value.to_string()
            }),
            Err(_) => Scalar::Other(text.to_owned()),
        };
    }
    if is_other_integer(text) || is_float(text) || is_timestamp(text) {
        return Scalar::Other(text.to_owned());
    }
    Scalar::Str(text.to_owned())
}

/// A JSON string literal for a value written into a YAML file, as the script
/// writes it (`json.dumps(value, ensure_ascii=False)`).
pub(crate) fn quoted_for_yaml(value: &str) -> String {
    json_text::quote(value, false)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mapping(text: &str) -> ConfigDocument {
        match parse_config(text).unwrap() {
            Parsed::Mapping(document) => document,
            other => panic!("not a mapping: {other:?}"),
        }
    }

    #[test]
    fn scalars_resolve_like_safe_load() {
        let document = mapping(
            "a: hg38   # note\nb: \"x y\"\nc: 'it''s'\nd: ~\ne:\nf: yes\ng: 0.0001\nh: 1_000\ni: 0755\nj: 08\nk: 2024-01-02\nl: -.5\nm: /p/a:b\nn: null\n",
        );
        assert_eq!(document.scalar("a"), Some(&Scalar::Str("hg38".into())));
        assert_eq!(document.scalar("b"), Some(&Scalar::Str("x y".into())));
        assert_eq!(document.scalar("c"), Some(&Scalar::Str("it's".into())));
        assert_eq!(document.scalar("d"), Some(&Scalar::Null));
        assert_eq!(document.scalar("e"), Some(&Scalar::Null));
        assert_eq!(document.scalar("f"), Some(&Scalar::Bool(true)));
        assert!(matches!(document.scalar("g"), Some(Scalar::Other(_))));
        assert_eq!(document.scalar("h"), Some(&Scalar::Int("1000".into())));
        assert!(matches!(document.scalar("i"), Some(Scalar::Other(_))));
        assert_eq!(document.scalar("j"), Some(&Scalar::Str("08".into())));
        assert!(matches!(document.scalar("k"), Some(Scalar::Other(_))));
        assert_eq!(document.scalar("l"), Some(&Scalar::Str("-.5".into())));
        assert_eq!(document.scalar("m"), Some(&Scalar::Str("/p/a:b".into())));
        assert_eq!(document.scalar("n"), Some(&Scalar::Null));
    }

    #[test]
    fn value_spans_cover_exactly_the_scalar() {
        let text = "# top\r\nassembly: hg38\r\nsnv: \"old\" # c\r\nindel:   # keep\r\nz: 1\r\n";
        let document = mapping(text);
        let snv = document.entry("snv").unwrap();
        assert_eq!(&text[snv.value_start..snv.value_end], "\"old\"");
        let assembly = document.entry("assembly").unwrap();
        assert_eq!(&text[assembly.value_start..assembly.value_end], "hg38");
        let indel = document.entry("indel").unwrap();
        assert_eq!(indel.value_start, indel.value_end);
        assert_eq!(
            &text[..indel.value_start],
            "# top\r\nassembly: hg38\r\nsnv: \"old\" # c\r\nindel:"
        );
    }

    #[test]
    fn nested_blocks_and_flow_lists_stay_untouched_text() {
        let text = "a: 1\ncache_folders:\n  vep: VEP\n  gnomad: gnomAD # c\nother: [one, two] # note\nlist:\n- x\n- y\nlast: z\n";
        let document = mapping(text);
        assert_eq!(
            document.get("cache_folders"),
            Some(&ConfigValue::Opaque(
                "  vep: VEP\n  gnomad: gnomAD # c\n".into()
            ))
        );
        assert_eq!(
            document.get("other"),
            Some(&ConfigValue::Opaque("[one, two]".into()))
        );
        assert_eq!(
            document.get("list"),
            Some(&ConfigValue::Opaque("- x\n- y\n".into()))
        );
        assert_eq!(document.scalar("last"), Some(&Scalar::Str("z".into())));
    }

    #[test]
    fn duplicate_keys_are_found_and_the_last_value_wins() {
        let document = mapping("a: 1\nb: 2\na: 3\n");
        assert!(document.has_duplicate_keys());
        assert_eq!(document.scalar("a"), Some(&Scalar::Int("3".into())));
        assert!(!mapping("a: 1\nb: 2\n").has_duplicate_keys());
    }

    #[test]
    fn empty_and_list_documents_are_told_apart() {
        assert!(matches!(parse_config(""), Ok(Parsed::Empty)));
        assert!(matches!(
            parse_config("# only a comment\n\n"),
            Ok(Parsed::Empty)
        ));
        assert!(matches!(parse_config("{}\n"), Ok(Parsed::Empty)));
        assert!(matches!(parse_config("- a\n- b\n"), Ok(Parsed::List)));
    }

    #[test]
    fn broken_or_unsupported_text_fails_loudly() {
        for (text, invalid_text) in [
            ("assembly: [unterminated\n", true),
            ("a: \"x\" y\n", true),
            ("a: b: c\n", true),
            ("a: 'x\n", false),
            ("a: |\n  text\n", false),
            ("a: &anchor x\n", false),
            ("a: x\n  continued\n", false),
            ("a:\tb\n", false),
            ("\"a\": 1\n", false),
            ("a: {b: 1}\n", false),
            ("---\na: 1\n", false),
            ("a: 1\n- b\n", true),
        ] {
            match parse_config(text) {
                Err(YamlError::Invalid(_)) => assert!(invalid_text, "{text:?}"),
                Err(YamlError::Unsupported(_)) => assert!(!invalid_text, "{text:?}"),
                other => panic!("{text:?} was accepted: {other:?}"),
            }
        }
    }

    #[test]
    fn python_repr_of_strings_and_scalars() {
        assert_eq!(python_str_repr("hg37"), "'hg37'");
        assert_eq!(python_str_repr("it's"), "\"it's\"");
        assert_eq!(python_str_repr("a\\b\n"), "'a\\\\b\\n'");
        assert_eq!(Scalar::Null.python_repr(), "None");
        assert_eq!(Scalar::Bool(false).python_repr(), "False");
        assert_eq!(Scalar::Int("38".into()).python_repr(), "38");
    }

    #[test]
    fn double_quoted_escapes_follow_json_and_yaml() {
        let document = mapping("a: \"p\\\\q\\\"r\\n\\u00e9\\x41\"\n");
        assert_eq!(
            document.scalar("a"),
            Some(&Scalar::Str("p\\q\"r\n\u{e9}A".into()))
        );
        assert_eq!(quoted_for_yaml("a\"b\\c"), "\"a\\\"b\\\\c\"");
    }
}
