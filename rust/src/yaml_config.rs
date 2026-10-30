//! A reader for the flat YAML configuration files PriVA uses (`config.yaml`).
//!
//! It handles what those files contain and refuses everything else loudly
//! instead of guessing:
//! - a block mapping at column 0 whose keys are plain names (letters, digits,
//!   `_`, `.`, `-`) or double- or single-quoted names on one line (spaces may
//!   come between a name and its colon), with `#`
//!   comments and blank lines, after at most one byte order mark and one
//!   `---` line (the scanner skips the mark and it takes no column);
//! - values that are plain, single-quoted or double-quoted scalars on one
//!   line, resolved to a null, a boolean, an integer, another typed scalar or
//!   a string the way PyYAML's `safe_load` does;
//! - a value that is a one-line flow list without quotes (`[one, two]`), or a
//!   nested block mapping or block list whose entries are single-line scalars
//!   (kept as text, since callers only need to leave them untouched);
//! - the byte span of each scalar value, so a caller can replace just that
//!   span and keep every comment, quote and line ending around it.
//! Anything else (anchors, aliases, tags, block scalars, multi-line scalars,
//! keys with spaces or other punctuation, deeper nesting) is reported as
//! `Unsupported`: it may be valid YAML that this reader chooses not to guess
//! at, so callers must refuse it loudly. Text PyYAML itself refuses (a tab where
//! a token starts, a second `---` document, a flow list that is not closed)
//! is `Invalid`.
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
            Scalar::Int(text) => text.clone(),
            Scalar::Other(text) => other_repr(text),
            Scalar::Str(text) => python_str_repr(text),
        }
    }

    /// `type(value).__name__` of the value `safe_load` would produce.
    pub(crate) fn python_type_name(&self) -> &'static str {
        match self {
            Scalar::Null => "NoneType",
            Scalar::Bool(_) => "bool",
            Scalar::Int(_) => "int",
            Scalar::Str(_) => "str",
            Scalar::Other(text) => {
                if is_float(text) {
                    "float"
                } else if is_timestamp(text) {
                    if text.len() == 10 { "date" } else { "datetime" }
                } else {
                    "int"
                }
            }
        }
    }
}

/// `repr()` of an integer written in another base or with `:` groups, of a
/// float, or of a date, as PyYAML's constructors build them. A timestamp with a
/// time of day and the merge and value markers (`<<`, `=`) are given as written.
fn other_repr(text: &str) -> String {
    if is_float(text) {
        if let Some(value) = float_value(text) {
            return crate::json_text::python_float_repr(value);
        }
    } else if is_other_integer(text) {
        if let Some(value) = integer_value(text) {
            return value.to_string();
        }
    } else if is_timestamp(text) && text.len() == 10 {
        let number = |part: &str| part.parse::<u32>().ok();
        if let (Some(year), Some(month), Some(day)) = (
            number(&text[..4]),
            number(&text[5..7]),
            number(&text[8..10]),
        ) {
            return format!("datetime.date({year}, {month}, {day})");
        }
    }
    text.to_owned()
}

/// `construct_yaml_float`
fn float_value(text: &str) -> Option<f64> {
    let cleaned: String = text
        .chars()
        .filter(|c| *c != '_')
        .collect::<String>()
        .to_lowercase();
    let (sign, body) = match cleaned.strip_prefix('-') {
        Some(rest) => (-1.0, rest),
        None => (1.0, cleaned.strip_prefix('+').unwrap_or(&cleaned)),
    };
    let value = if body == ".inf" {
        f64::INFINITY
    } else if body == ".nan" {
        f64::NAN
    } else if body.contains(':') {
        let mut total = 0.0;
        let mut base = 1.0;
        for part in body.split(':').rev() {
            total += part.parse::<f64>().ok()? * base;
            base *= 60.0;
        }
        total
    } else {
        body.parse::<f64>().ok()?
    };
    Some(sign * value)
}

/// `construct_yaml_int` for the forms that are not plain decimal numbers.
fn integer_value(text: &str) -> Option<i128> {
    let cleaned: String = text.chars().filter(|c| *c != '_').collect();
    let (sign, body) = match cleaned.strip_prefix('-') {
        Some(rest) => (-1, rest),
        None => (1, cleaned.strip_prefix('+').unwrap_or(&cleaned)),
    };
    let value = if body == "0" {
        0
    } else if let Some(digits) = body.strip_prefix("0b") {
        i128::from_str_radix(digits, 2).ok()?
    } else if let Some(digits) = body.strip_prefix("0x") {
        i128::from_str_radix(digits, 16).ok()?
    } else if body.starts_with('0') {
        i128::from_str_radix(body, 8).ok()?
    } else if body.contains(':') {
        let mut total: i128 = 0;
        let mut base: i128 = 1;
        for part in body.split(':').rev() {
            total += part.parse::<i128>().ok()? * base;
            base *= 60;
        }
        total
    } else {
        body.parse::<i128>().ok()?
    };
    Some(sign * value)
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

impl ConfigValue {
    /// Python's `repr()` of the value `safe_load` would produce: a list or dict
    /// of resolved scalars for a nested block or a flow list.
    pub(crate) fn python_repr(&self) -> String {
        match self {
            ConfigValue::Scalar(scalar) => scalar.python_repr(),
            ConfigValue::Opaque(text) => opaque_repr(text).unwrap_or_else(|| text.clone()),
        }
    }

    /// `type(value).__name__`
    pub(crate) fn python_type_name(&self) -> &'static str {
        match self {
            ConfigValue::Scalar(scalar) => scalar.python_type_name(),
            ConfigValue::Opaque(text) => {
                let first = text.trim_start().chars().next();
                let is_list = text.starts_with('[') || first == Some('-');
                if is_list { "list" } else { "dict" }
            }
        }
    }
}

/// The scalar a one-line nested entry holds.
fn entry_scalar(entry: &str) -> Option<Scalar> {
    match inline_value(entry) {
        Ok((ConfigValue::Scalar(scalar), _)) => Some(scalar),
        _ => None,
    }
}

/// `repr()` of the list or dict an opaque value stands for.
fn opaque_repr(text: &str) -> Option<String> {
    if text == "{}" {
        return Some("{}".to_owned());
    }
    if let Some(inner) = text
        .strip_prefix('[')
        .and_then(|rest| rest.strip_suffix(']'))
    {
        let mut items: Vec<&str> = inner.split(',').collect();
        if items.len() > 1 && items[items.len() - 1].trim().is_empty() {
            items.pop();
        }
        if items.len() == 1 && items[0].trim().is_empty() {
            return Some("[]".to_owned());
        }
        let parts: Vec<String> = items
            .iter()
            .map(|item| resolve_plain(item.trim()).python_repr())
            .collect();
        return Some(format!("[{}]", parts.join(", ")));
    }
    let lines = split_lines(text);
    let mut list: Vec<String> = Vec::new();
    let mut pairs: Vec<(String, String)> = Vec::new();
    for line in lines.iter().filter(|line| !line.is_skippable()) {
        let content = &line.text[line.indent..];
        if line.is_list_item() {
            let entry = content[1..].trim_start_matches(' ');
            list.push(entry_scalar(entry)?.python_repr());
        } else {
            let length = key_length(content)?;
            let key = resolve_plain(content[..length].trim_end_matches(' ')).python_repr();
            let value = entry_scalar(content[length + 1..].trim_start_matches(' '))?.python_repr();
            // A repeated key keeps its first position and its last value.
            match pairs.iter_mut().find(|(name, _)| *name == key) {
                Some(slot) => slot.1 = value,
                None => pairs.push((key, value)),
            }
        }
    }
    if !list.is_empty() {
        return Some(format!("[{}]", list.join(", ")));
    }
    let parts: Vec<String> = pairs.iter().map(|(k, v)| format!("{k}: {v}")).collect();
    Some(format!("{{{}}}", parts.join(", ")))
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
    // PyYAML skips one byte order mark at the very start of the text. The
    // mark counts as a character in the positions it reports but not as a
    // column, so the rest is read as if it were not there and the byte spans
    // are moved back by the mark's length.
    if let Some(body) = text.strip_prefix('\u{feff}') {
        let mark = '\u{feff}'.len_utf8();
        return Ok(match parse_body(body)? {
            Parsed::Mapping(mut document) => {
                for entry in &mut document.entries {
                    entry.value_start += mark;
                    entry.value_end += mark;
                }
                Parsed::Mapping(document)
            }
            other => other,
        });
    }
    parse_body(text)
}

/// A `---` line, alone or followed only by a comment.
fn is_document_start(line: &Line<'_>) -> bool {
    if line.indent != 0 {
        return false;
    }
    let content = line.text.trim_end_matches(' ');
    content == "---"
        || (content.starts_with("--- ") && content[4..].trim_start_matches(' ').starts_with('#'))
}

fn parse_body(text: &str) -> Result<Parsed, YamlError> {
    for (index, c) in text.char_indices() {
        if !is_printable(c) {
            return invalid("unacceptable character in the stream");
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
    // One explicit document start may come before the content; PyYAML reads
    // what follows it as the document (and nothing after it as an empty one).
    if first < lines.len() && is_document_start(&lines[first]) {
        first += 1;
        while first < lines.len() && lines[first].is_skippable() {
            first += 1;
        }
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
        if content == "..." || content.starts_with("--- ") || content.starts_with('%') {
            return unsupported(
                "a document end marker, a document start with content, or a directive",
            );
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
        if is_document_start(line) {
            return invalid("expected a single document in the stream");
        }
        if line.indent != 0 {
            return unsupported("a mapping that does not start at column 0");
        }
        let (key, key_len) = match key_token(line.text)? {
            Some(found) => found,
            None => return Err(not_a_key(line)),
        };
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

/// Why a line that is not a key could not be read.
fn not_a_key(line: &Line<'_>) -> YamlError {
    let content = &line.text[line.indent..];
    if has_stray_tab(content) {
        return tab_error();
    }
    if line.is_list_item() {
        return YamlError::Invalid("sequence entries are not allowed here".to_owned());
    }
    // A plain word with no colon at all cannot be a key or a value here.
    if content.starts_with(|c: char| c.is_ascii_alphanumeric() || c == '_')
        && !content.contains(':')
    {
        return YamlError::Invalid("could not find expected ':'".to_owned());
    }
    YamlError::Unsupported(
        "a line that is not a plain or quoted key followed by a colon".to_owned(),
    )
}

fn tab_error() -> YamlError {
    YamlError::Invalid("found character '\\t' that cannot start any token".to_owned())
}

/// Whether a tab shows up where PyYAML looks for the next token. Tabs inside
/// comments and quoted text are fine, so the search stops at either.
fn has_stray_tab(text: &str) -> bool {
    let mut previous = ' ';
    for c in text.chars() {
        if c == '"' || c == '\'' || (c == '#' && previous == ' ') {
            return false;
        }
        if c == '\t' {
            return true;
        }
        previous = c;
    }
    false
}

/// A PyYAML simple key has to end within 1024 characters of where it starts.
const MAX_KEY_CHARACTERS: usize = 1024;

/// The key at the start of a line and the byte position of the colon after it.
/// A key is a plain name or a quoted scalar, then optional spaces, then a colon
/// followed by a space or the line end.
fn key_token(text: &str) -> Result<Option<(String, usize)>, YamlError> {
    let (name, colon) = if text.starts_with('"') || text.starts_with('\'') {
        let (name, length) = if text.starts_with('"') {
            scan_double(text)?
        } else {
            scan_single(text)?
        };
        let rest = &text[length..];
        let spaces = rest.len() - rest.trim_start_matches(' ').len();
        let after = &rest[spaces..];
        if !(after == ":" || after.starts_with(": ")) {
            return Ok(None);
        }
        (name, length + spaces)
    } else {
        match key_length(text) {
            Some(colon) => (text[..colon].trim_end_matches(' ').to_owned(), colon),
            None => return Ok(None),
        }
    };
    if text[..colon].chars().count() > MAX_KEY_CHARACTERS {
        return Err(YamlError::Invalid("could not find expected ':'".to_owned()));
    }
    Ok(Some((name, colon)))
}

/// Position of the colon after a plain key (letters, digits, `_`, `.`, `-`) and
/// optional spaces, when a space or the line end follows the colon.
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
    while length < bytes.len() && bytes[length] == b' ' {
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
                None => {
                    return Err(if has_stray_tab(content) {
                        tab_error()
                    } else {
                        YamlError::Unsupported("a nested entry that is not a plain key".to_owned())
                    });
                }
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
    if text.contains('\t') {
        return Err(tab_error());
    }
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
    if inner.contains('\t') {
        return Err(tab_error());
    }
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
            ("a:\tb\n", true),
            ("a: b\t\n", true),
            ("\tb: 1\n", true),
            ("a: [x,\ty]\n", true),
            ("a: {b: 1}\n", false),
            ("--- a: 1\n", false),
            ("a: 1\n---\nb: 2\n", true),
            ("---\n---\na: 1\n", true),
            ("a: 1\nword\n", true),
            ("a: 1\nmy key: 2\n", false),
            ("\"a\" b: 1\n", false),
            ("\"a\":b\n", false),
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
    fn a_byte_order_mark_and_a_document_start_are_skipped() {
        let text = "\u{feff}# top\n---\nassembly: hg38\nsnv: old\n";
        let document = mapping(text);
        assert_eq!(
            document.scalar("assembly"),
            Some(&Scalar::Str("hg38".into()))
        );
        // The spans refer to the text with the mark in front.
        let snv = document.entry("snv").unwrap();
        assert_eq!(&text[snv.value_start..snv.value_end], "old");
        let marked = mapping("---   # note\nassembly: hg38\n");
        assert_eq!(marked.scalar("assembly"), Some(&Scalar::Str("hg38".into())));
        assert!(matches!(parse_config("---\n"), Ok(Parsed::Empty)));
        assert!(matches!(parse_config("\u{feff}"), Ok(Parsed::Empty)));
        assert!(matches!(parse_config("---\n- a\n"), Ok(Parsed::List)));
        assert!(matches!(parse_config("---\n[]\n"), Ok(Parsed::Empty)));
        // A mark that is not first is an ordinary character of a plain value.
        let inside = mapping("a: \u{feff}b\n");
        assert_eq!(inside.scalar("a"), Some(&Scalar::Str("\u{feff}b".into())));
    }

    #[test]
    fn quoted_keys_are_names_like_plain_ones() {
        let text = "\"assembly\": hg38\n'it''s': 1 # c\n\"a\\u0062c\":\nplain: 2\n";
        let document = mapping(text);
        assert_eq!(
            document.scalar("assembly"),
            Some(&Scalar::Str("hg38".into()))
        );
        assert_eq!(document.scalar("it's"), Some(&Scalar::Int("1".into())));
        assert_eq!(document.scalar("abc"), Some(&Scalar::Null));
        let entry = document.entry("assembly").unwrap();
        assert_eq!(&text[entry.value_start..entry.value_end], "hg38");
        assert!(mapping("assembly: a\n\"assembly\": b\n").has_duplicate_keys());
        let spaced = mapping("plain : 1\n\"quoted\"  : 2\n");
        assert_eq!(spaced.scalar("plain"), Some(&Scalar::Int("1".into())));
        assert_eq!(spaced.scalar("quoted"), Some(&Scalar::Int("2".into())));
        let long_key = format!("\"{}\": 1\n", "k".repeat(1024));
        assert!(matches!(
            parse_config(&long_key),
            Err(YamlError::Invalid(_))
        ));
        let limit_key = format!("{}: 1\n", "k".repeat(1024));
        assert!(matches!(parse_config(&limit_key), Ok(Parsed::Mapping(_))));
    }

    #[test]
    fn tabs_in_comments_and_quotes_are_fine() {
        let document = mapping("a: \"x\ty\" # note\tmore\n# line\twith tab\nb: c # d\te\n");
        assert_eq!(document.scalar("a"), Some(&Scalar::Str("x\ty".into())));
        assert_eq!(document.scalar("b"), Some(&Scalar::Str("c".into())));
    }

    #[test]
    fn nested_values_have_python_reprs_and_type_names() {
        let document = mapping(
            "a: [hg38, 5]\nb:\n  name: hg38\n  n: 2\nc:\n- x\n- null\nd: []\ne: {}\nf: 1.5\ng: 2024-01-02\nh: 0x1f\n",
        );
        let repr = |key: &str| document.get(key).unwrap().python_repr();
        let kind = |key: &str| document.get(key).unwrap().python_type_name();
        assert_eq!(repr("a"), "['hg38', 5]");
        assert_eq!(repr("b"), "{'name': 'hg38', 'n': 2}");
        assert_eq!(repr("c"), "['x', None]");
        assert_eq!(repr("d"), "[]");
        assert_eq!(repr("e"), "{}");
        assert_eq!(
            [
                kind("a"),
                kind("b"),
                kind("c"),
                kind("d"),
                kind("e"),
                kind("f"),
                kind("g"),
                kind("h")
            ],
            [
                "list", "dict", "list", "list", "dict", "float", "date", "int"
            ]
        );
    }

    #[test]
    fn python_repr_of_strings_and_scalars() {
        assert_eq!(python_str_repr("hg37"), "'hg37'");
        assert_eq!(python_str_repr("it's"), "\"it's\"");
        assert_eq!(python_str_repr("a\\b\n"), "'a\\\\b\\n'");
        assert_eq!(Scalar::Null.python_repr(), "None");
        assert_eq!(Scalar::Bool(false).python_repr(), "False");
        assert_eq!(Scalar::Int("38".into()).python_repr(), "38");
        assert_eq!(resolve_plain("0x26").python_repr(), "38");
        assert_eq!(resolve_plain("0755").python_repr(), "493");
        assert_eq!(resolve_plain("-0b101").python_repr(), "-5");
        assert_eq!(resolve_plain("1:30").python_repr(), "90");
        assert_eq!(resolve_plain("1_0.50").python_repr(), "10.5");
        assert_eq!(resolve_plain("1.0e+2").python_repr(), "100.0");
        assert_eq!(resolve_plain("-.inf").python_repr(), "-inf");
        assert_eq!(resolve_plain(".NaN").python_repr(), "nan");
        assert_eq!(resolve_plain("1:30.5").python_repr(), "90.5");
        assert_eq!(
            resolve_plain("2024-01-02").python_repr(),
            "datetime.date(2024, 1, 2)"
        );
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
