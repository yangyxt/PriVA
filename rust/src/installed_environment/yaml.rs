//! Safe YAML construction plus package-only emission. libyaml owns syntax;
//! ordered nodes preserve PyYAML's duplicate-key and container-alias behavior.
use super::ffi::{Event, Mark};
use super::*;
use crate::yaml_config::{Scalar, resolve_plain};
use std::collections::{HashMap, HashSet, VecDeque};
const TAG: &str = "tag:yaml.org,2002:";
#[derive(Debug)]
struct YamlFailure {
    category: &'static str,
    message: String,
}
impl std::fmt::Display for YamlFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}
impl std::error::Error for YamlFailure {}
fn error(category: &'static str, message: impl Into<String>) -> Box<dyn std::error::Error> {
    Box::new(YamlFailure {
        category,
        message: message.into(),
    })
}
fn named_failure(cause: Box<dyn std::error::Error>, name: &str) -> Box<dyn std::error::Error> {
    match cause.downcast::<YamlFailure>() {
        Ok(cause) => super::error(
            cause.category,
            cause.message.replace("<unicode string>", name),
        ),
        Err(cause) => cause,
    }
}

#[derive(Clone, Debug)]
enum Raw {
    Scalar(String),
    Sequence(Vec<usize>),
    Mapping(Vec<usize>),
}
impl Raw {
    fn kind(&self) -> &'static str {
        match self {
            Self::Scalar(_) => "scalar",
            Self::Sequence(_) => "sequence",
            Self::Mapping(_) => "mapping",
        }
    }
    fn append(&mut self, node: usize) {
        match self {
            Self::Sequence(v) | Self::Mapping(v) => v.push(node),
            _ => unreachable!(),
        }
    }
}
#[derive(Clone, Debug)]
pub(crate) enum Value {
    Null,
    Bool(bool),
    Number(String),
    Other(String),
    Str(String),
    /// Python's pure-Python scanner retains escaped surrogate code points.
    /// This is constructed only by the optional UniqueLoader profile; it is
    /// not coerced to replacement characters or silently combined as UTF-16.
    Codepoints {
        points: Vec<u32>,
        transport: String,
    },
    Sequence(Vec<usize>),
    Mapping(Vec<(usize, usize)>),
}
#[derive(Clone, Debug)]
struct Node {
    raw: Raw,
    tag: String,
    mark: Mark,
    value: Option<Value>,
}
pub(crate) struct Yaml {
    nodes: Vec<Node>,
    pub root: usize,
    unique: bool,
    constructing: HashSet<usize>,
    escaped_surrogates: HashMap<char, u32>,
}
impl Yaml {
    /// The resolved SafeLoader tag distinguishes integer and float values even
    /// when their canonical numeric spelling happens to be identical.
    pub(crate) fn tag(&self, node: usize) -> &str {
        self.nodes[node]
            .tag
            .strip_prefix(TAG)
            .unwrap_or(&self.nodes[node].tag)
    }
    pub fn value(&self, node: usize) -> &Value {
        self.nodes[node].value.as_ref().unwrap()
    }
    pub fn string(&self, node: usize) -> Option<&str> {
        match self.value(node) {
            Value::Str(s) | Value::Codepoints { transport: s, .. } => Some(s),
            _ => None,
        }
    }
    pub fn get(&self, node: usize, key: &str) -> Option<usize> {
        if let Value::Mapping(v) = self.value(node) {
            v.iter()
                .find(|(k, _)| self.string(*k) == Some(key))
                .map(|(_, v)| *v)
        } else {
            None
        }
    }
    pub fn sequence(&self, node: usize) -> Option<&[usize]> {
        if let Value::Sequence(v) = self.value(node) {
            Some(v)
        } else {
            None
        }
    }
    pub fn mapping(&self, node: usize) -> Option<&[(usize, usize)]> {
        if let Value::Mapping(v) = self.value(node) {
            Some(v)
        } else {
            None
        }
    }
    fn add(&mut self, value: Value) -> usize {
        let id = self.nodes.len();
        self.nodes.push(Node {
            raw: Raw::Scalar(String::new()),
            tag: String::new(),
            mark: Mark::default(),
            value: Some(value),
        });
        id
    }
    pub fn package_root(&mut self, name: &str, channels: Option<usize>, dependencies: usize) {
        let name_key = self.add(Value::Str("name".into()));
        let name_value = self.add(Value::Str(name.into()));
        let channels_key = self.add(Value::Str("channels".into()));
        let channels = channels.unwrap_or_else(|| self.add(Value::Sequence(vec![])));
        let deps_key = self.add(Value::Str("dependencies".into()));
        self.root = self.add(Value::Mapping(vec![
            (name_key, name_value),
            (channels_key, channels),
            (deps_key, dependencies),
        ]));
    }
}

fn mark_text(text: &str, mark: Mark) -> String {
    let chars: Vec<char> = text.chars().chain(std::iter::once('\0')).collect();
    let pointer = mark.index.min(chars.len() - 1);
    let stop = |c| matches!(c, '\0' | '\r' | '\n' | '\u{85}' | '\u{2028}' | '\u{2029}');
    let mut start = pointer;
    let mut head = "";
    while start > 0 && !stop(chars[start - 1]) {
        start -= 1;
        if pointer - start > 36 {
            head = " ... ";
            start += 5;
            break;
        }
    }
    let mut end = pointer;
    let mut tail = "";
    while end < chars.len() && !stop(chars[end]) {
        end += 1;
        if end - pointer > 36 {
            tail = " ... ";
            end -= 5;
            break;
        }
    }
    let snippet: String = chars[start..end].iter().collect();
    format!(
        "  in \"<unicode string>\", line {}, column {}:\n    {head}{snippet}{tail}\n{}^",
        mark.line + 1,
        mark.column + 1,
        " ".repeat(4 + pointer - start + head.len())
    )
}
fn marked(
    text: &str,
    category: &'static str,
    context: Option<(&str, Mark)>,
    problem: &str,
    mark: Mark,
) -> Box<dyn std::error::Error> {
    let mut parts = Vec::new();
    if let Some((message, location)) = context {
        parts.push(message.to_owned());
        if location.line != mark.line || location.column != mark.column {
            parts.push(mark_text(text, location));
        }
    }
    parts.push(problem.to_owned());
    parts.push(mark_text(text, mark));
    error(category, parts.join("\n"))
}
fn resolved_tag(text: &str) -> String {
    if text == "<<" {
        return format!("{TAG}merge");
    }
    if text == "=" {
        return format!("{TAG}value");
    }
    let kind = match resolve_plain(text) {
        Scalar::Str(_) if text == "<<" => "merge",
        Scalar::Str(_) if text == "=" => "value",
        Scalar::Str(_) => "str",
        Scalar::Null => "null",
        Scalar::Bool(_) => "bool",
        Scalar::Int(_) => "int",
        scalar @ Scalar::Other(_) => match scalar.python_type_name() {
            "date" | "datetime" => "timestamp",
            "float" => "float",
            _ => "int",
        },
    };
    format!("{TAG}{kind}")
}
fn parser_error(text: &str, mut e: ffi::ParseError) -> Box<dyn std::error::Error> {
    if e.kind == 2 {
        return error(
            "ReaderError",
            format!(
                "unacceptable character #x{:04x}: {}\n  in \"<unicode string>\", position {}",
                e.value, e.problem, e.offset
            ),
        );
    }
    let category = if e.kind == 3 {
        "ScannerError"
    } else {
        "ParserError"
    };
    let chars: Vec<char> = text.chars().collect();
    let token = chars.get(e.problem_mark.index).copied();
    let token_name = match token {
        None => "<stream end>",
        Some('}') => "}",
        Some(']') => "]",
        Some(_) => "<scalar>",
    };
    if e.problem == "did not find expected ',' or ']'" {
        e.problem = format!("expected ',' or ']', but got '{token_name}'");
    } else if e.problem == "did not find expected ',' or '}'" {
        e.problem = format!("expected ',' or '}}', but got '{token_name}'");
    } else if e.problem == "found character that cannot start any token" {
        e.context = "while scanning for the next token".into();
        e.problem = format!(
            "found character {} that cannot start any token",
            repr(&token.unwrap_or('\0').to_string())
        );
        return marked(
            text,
            category,
            None,
            &format!("{}\n{}", e.context, e.problem),
            e.problem_mark,
        );
    } else if e.problem == "found unknown escape character" {
        e.context = "while scanning a double-quoted scalar".into();
        e.problem_mark.index += 1;
        e.problem_mark.column += 1;
        e.problem = format!(
            "found unknown escape character {}",
            repr(&chars.get(e.problem_mark.index).unwrap_or(&'\0').to_string())
        );
    } else if e.problem == "did not find expected key"
        && e.context == "while parsing a block mapping"
    {
        e.problem = format!("expected <block end>, but found '{token_name}'");
    }
    if e.problem == "did not find expected node content" {
        e.problem = format!("expected the node content, but found '{token_name}'");
    }
    let same_context = e.context_mark == e.problem_mark;
    // libyaml advances an implicit stream end to a new line. PyYAML retains
    // the actual final cursor, which is also the snippet's caret position.
    if token.is_none() && !text.ends_with(['\n', '\r', '\u{85}', '\u{2028}', '\u{2029}']) {
        e.problem_mark.index = chars.len();
        e.problem_mark.line = chars.iter().filter(|c| **c == '\n').count();
        e.problem_mark.column = chars.iter().rev().take_while(|c| **c != '\n').count();
    }
    if same_context {
        e.context_mark = e.problem_mark;
    }
    let context = (!e.context.is_empty()).then_some((e.context.as_str(), e.context_mark));
    marked(text, category, context, &e.problem, e.problem_mark)
}

fn load_impl(text: &str, unique: bool) -> Result<Yaml> {
    // PyYAML's Reader checks printable code points before composition.
    for (position, c) in text.chars().enumerate() {
        if !matches!(c as u32,9|10|13|0x20..=0x7e|0x85|0xa0..=0xd7ff|0xe000..=0xfffd|0x10000..=0x10ffff)
        {
            return Err(error(
                "ReaderError",
                format!(
                    "unacceptable character #x{:04x}: special characters are not allowed\n  in \"<unicode string>\", position {position}",
                    c as u32
                ),
            ));
        }
    }
    let (parser_text, escaped_surrogates) = if unique {
        unique_surrogate_input(text)
    } else {
        (text.to_owned(), HashMap::new())
    };
    let mut yaml = Yaml {
        nodes: vec![],
        root: 0,
        unique,
        constructing: HashSet::new(),
        escaped_surrogates,
    };
    let mut stack: Vec<usize> = vec![];
    let mut anchors = HashMap::<String, usize>::new();
    let mut failure = None;
    let mut docs = 0;
    let result = ffi::parse(&parser_text, &mut |event| {
        if event.kind == 3 {
            docs += 1;
            if docs > 1 {
                failure = Some(marked(
                    text,
                    "ComposerError",
                    Some((
                        "expected a single document in the stream",
                        yaml.nodes[yaml.root].mark,
                    )),
                    "but found another document",
                    event.mark,
                ));
                return false;
            }
        }
        if matches!(event.kind, 8 | 10) {
            stack.pop();
            return true;
        }
        if !matches!(event.kind, 5 | 6 | 7 | 9) {
            return true;
        }
        let node = if event.kind == 5 {
            if let Some(node) = anchors.get(&event.anchor) {
                *node
            } else {
                failure = Some(marked(
                    text,
                    "ComposerError",
                    None,
                    &format!("found undefined alias '{}'", event.anchor),
                    event.mark,
                ));
                return false;
            }
        } else {
            if let Some(first) = anchors
                .get(&event.anchor)
                .filter(|_| !event.anchor.is_empty())
            {
                failure = Some(marked(
                    text,
                    "ComposerError",
                    Some((
                        &format!(
                            "found duplicate anchor '{}'; first occurrence",
                            event.anchor
                        ),
                        yaml.nodes[*first].mark,
                    )),
                    "second occurrence",
                    event.mark,
                ));
                return false;
            }
            let raw = match event.kind {
                6 => Raw::Scalar(event.value.clone()),
                7 => Raw::Sequence(vec![]),
                _ => Raw::Mapping(vec![]),
            };
            let tag = if !event.tag.is_empty() && event.tag != "!" {
                event.tag
            } else {
                match &raw {
                    Raw::Scalar(v) if event.implicit => resolved_tag(v),
                    Raw::Scalar(_) => format!("{TAG}str"),
                    Raw::Sequence(_) => format!("{TAG}seq"),
                    Raw::Mapping(_) => format!("{TAG}map"),
                }
            };
            let id = yaml.nodes.len();
            yaml.nodes.push(Node {
                raw,
                tag,
                mark: event.mark,
                value: None,
            });
            if !event.anchor.is_empty() {
                anchors.insert(event.anchor, id);
            }
            id
        };
        if let Some(parent) = stack.last() {
            yaml.nodes[*parent].raw.append(node);
        } else {
            yaml.root = node;
        }
        if matches!(event.kind, 7 | 9) {
            stack.push(node);
        }
        true
    });
    if let Some(error) = failure {
        return Err(error);
    }
    result.map_err(|e| parser_error(text, e))?;
    if yaml.nodes.is_empty() {
        yaml.root = yaml.add(Value::Null);
        return Ok(yaml);
    }
    let mut pending = VecDeque::new();
    yaml.construct(yaml.root, text, &mut pending)?;
    while let Some(node) = pending.pop_front() {
        yaml.fill(node, text, &mut pending)?;
    }
    Ok(yaml)
}

/// libyaml rejects escaped surrogates that PyYAML's Python scanner accepts.
/// Let libyaml identify actual escapes (never rewrite comments/plain/single
/// quoted text), substitute unused BMP scalars without shifting source marks,
/// then restore the original code points when constructing string nodes.
fn unique_surrogate_input(text: &str) -> (String, HashMap<char, u32>) {
    let mut input = text.to_owned();
    let mut replacements = HashMap::new();
    let mut occupied: HashSet<u32> = text.chars().map(u32::from).collect();
    let bytes = text.as_bytes();
    for at in 0..bytes.len().saturating_sub(2) {
        let width = match bytes.get(at..at + 2) {
            Some(b"\\u") => 4,
            Some(b"\\U") => 8,
            Some(b"\\x") => 2,
            _ => continue,
        };
        if let Some(value) = bytes
            .get(at + 2..at + 2 + width)
            .and_then(|v| std::str::from_utf8(v).ok())
            .and_then(|v| u32::from_str_radix(v, 16).ok())
        {
            occupied.insert(value);
        }
    }
    loop {
        let Err(cause) = ffi::parse(&input, &mut |_| true) else {
            break;
        };
        if cause.problem != "found invalid Unicode character escape code" {
            break;
        }
        let Some(at) = input
            .char_indices()
            .nth(cause.problem_mark.index)
            .map(|(i, _)| i)
        else {
            break;
        };
        let width = if at >= 2 && &input.as_bytes()[at - 2..at] == b"\\u" {
            4
        } else if at >= 2 && &input.as_bytes()[at - 2..at] == b"\\U" {
            8
        } else {
            break;
        };
        let Some(code) = input
            .get(at..at + width)
            .and_then(|v| u32::from_str_radix(v, 16).ok())
        else {
            break;
        };
        if !(0xd800..=0xdfff).contains(&code) {
            break;
        }
        let existing = replacements
            .iter()
            .find_map(|(key, value)| (*value == code).then_some(*key));
        let replacement = existing.or_else(|| {
            (0xe000..=0xfffd)
                .chain(0x100..=0xd7ff)
                .find(|v| !occupied.contains(v))
                .and_then(char::from_u32)
        });
        let Some(replacement) = replacement else {
            break;
        };
        occupied.insert(replacement as u32);
        replacements.insert(replacement, code);
        input.replace_range(at..at + width, &format!("{:0width$x}", replacement as u32));
    }
    (input, replacements)
}

pub(crate) fn load(text: &str) -> Result<Yaml> {
    load_impl(text, false).map_err(|e| named_failure(e, "<unicode string>"))
}
/// `installer_config_paths.UniqueLoader`: its custom map constructor is
/// immediate, rejects duplicate Python keys before constructing their values,
/// and deliberately does not apply SafeLoader's merge-key flattening.
/// The ordinary safe-loader entry points keep their preceding behavior.
pub(crate) fn load_unique(text: &str) -> Result<Yaml> {
    load_impl(text, true).map_err(|e| named_failure(e, "<unicode string>"))
}
pub(super) fn load_bytes(bytes: &[u8]) -> Result<Yaml> {
    fn invalid_byte(bytes: &[u8], position: usize, reason: &str) -> Box<dyn std::error::Error> {
        super::error(
            "ReaderError",
            format!(
                "unacceptable character #x{:04x}: {reason}\n  in \"<byte string>\", position {position}",
                bytes[position]
            ),
        )
    }
    let text = if bytes.starts_with(&[0xff, 0xfe]) || bytes.starts_with(&[0xfe, 0xff]) {
        let little = bytes.starts_with(&[0xff, 0xfe]);
        let mut text = String::new();
        let mut at = 0;
        let word = |at| {
            if little {
                u16::from_le_bytes([bytes[at], bytes[at + 1]])
            } else {
                u16::from_be_bytes([bytes[at], bytes[at + 1]])
            }
        };
        while at < bytes.len() {
            if at + 1 >= bytes.len() {
                return Err(invalid_byte(bytes, at, "truncated data"));
            }
            let first = word(at);
            if (0xd800..=0xdbff).contains(&first) {
                if at + 3 >= bytes.len() {
                    return Err(invalid_byte(bytes, at, "unexpected end of data"));
                }
                let second = word(at + 2);
                if !(0xdc00..=0xdfff).contains(&second) {
                    return Err(invalid_byte(bytes, at, "illegal UTF-16 surrogate"));
                }
                text.push(
                    char::from_u32(
                        0x10000 + ((u32::from(first) - 0xd800) << 10) + u32::from(second) - 0xdc00,
                    )
                    .unwrap(),
                );
                at += 4;
            } else if (0xdc00..=0xdfff).contains(&first) {
                return Err(invalid_byte(bytes, at, "illegal encoding"));
            } else {
                text.push(char::from_u32(u32::from(first)).unwrap());
                at += 2;
            }
        }
        text
    } else {
        match std::str::from_utf8(bytes) {
            Ok(text) => text.into(),
            Err(e) => {
                let detail = crate::python_io::decode_utf8(bytes).unwrap_err();
                return Err(invalid_byte(
                    bytes,
                    e.valid_up_to(),
                    detail.rsplit_once(": ").unwrap().1,
                ));
            }
        }
    };
    load_impl(&text, false).map_err(|e| named_failure(e, "<byte string>"))
}

impl Yaml {
    fn scalar(&self, node: usize, text: &str) -> Result<String> {
        self.scalar_inner(node, text, 0)
    }
    fn scalar_inner(&self, node: usize, text: &str, depth: usize) -> Result<String> {
        if depth > 900 {
            return Err(error("RecursionError", "maximum recursion depth exceeded"));
        }
        match &self.nodes[node].raw {
            Raw::Scalar(v) => Ok(v.clone()),
            Raw::Mapping(v) => {
                for pair in v.chunks_exact(2) {
                    if self.nodes[pair[0]].tag == format!("{TAG}value") {
                        return self.scalar_inner(pair[1], text, depth + 1);
                    }
                }
                Err(self.expected(node, text, "scalar"))
            }
            _ => Err(self.expected(node, text, "scalar")),
        }
    }
    fn expected(&self, node: usize, text: &str, kind: &str) -> Box<dyn std::error::Error> {
        marked(
            text,
            "ConstructorError",
            None,
            &format!(
                "expected a {kind} node, but found {}",
                self.nodes[node].raw.kind()
            ),
            self.nodes[node].mark,
        )
    }
    fn construct(&mut self, node: usize, text: &str, pending: &mut VecDeque<usize>) -> Result<()> {
        if self.nodes[node].value.is_some() {
            return Ok(());
        }
        if self.unique && self.nodes[node].tag == format!("{TAG}map") {
            if !self.constructing.insert(node) {
                return Err(marked(
                    text,
                    "ConstructorError",
                    None,
                    "found unconstructable recursive node",
                    self.nodes[node].mark,
                ));
            }
            if self.constructing.len() > 300 {
                return Err(error("RecursionError", "maximum recursion depth exceeded"));
            }
            let Raw::Mapping(children) = self.nodes[node].raw.clone() else {
                return Err(error("TypeError", "cannot unpack non-iterable object"));
            };
            let mut entries: Vec<(usize, usize)> = vec![];
            for pair in children.chunks_exact(2) {
                let (key, value) = (pair[0], pair[1]);
                self.construct(key, text, pending)?;
                let unhashable = match self.value(key) {
                    Value::Sequence(_) => Some("list"),
                    Value::Mapping(_) => Some("dict"),
                    Value::Other(s) if s == "set" => Some("set"),
                    _ => None,
                };
                if let Some(kind) = unhashable {
                    return Err(error("TypeError", format!("unhashable type: '{kind}'")));
                }
                if entries
                    .iter()
                    .any(|(other, _)| self.unique_key_equal(*other, key))
                {
                    return Err(error(
                        "ValueError",
                        format!(
                            "Duplicate configuration key: {}",
                            crate::python_unicode::display(&self.scalar_display(key))
                        ),
                    ));
                }
                self.construct(value, text, pending)?;
                entries.push((key, value));
            }
            self.nodes[node].value = Some(Value::Mapping(entries));
            self.constructing.remove(&node);
            return Ok(());
        }
        let tag = self.nodes[node]
            .tag
            .strip_prefix(TAG)
            .unwrap_or(&self.nodes[node].tag)
            .to_owned();
        let value = match tag.as_str() {
            "seq" => {
                pending.push_back(node);
                Value::Sequence(vec![])
            }
            "map" => {
                pending.push_back(node);
                Value::Mapping(vec![])
            }
            "str" => {
                let value = self.scalar(node, text)?;
                if value
                    .chars()
                    .any(|c| self.escaped_surrogates.contains_key(&c))
                {
                    let points = value
                        .chars()
                        .map(|c| self.escaped_surrogates.get(&c).copied().unwrap_or(c as u32))
                        .collect::<Vec<_>>();
                    Value::Codepoints {
                        transport: crate::python_unicode::encode_points(points.iter().copied()),
                        points,
                    }
                } else {
                    Value::Str(if self.unique {
                        crate::python_unicode::encode(&value)
                    } else {
                        value
                    })
                }
            }
            "null" => {
                self.scalar(node, text)?;
                Value::Null
            }
            "bool" => {
                let v = self.scalar(node, text)?.to_lowercase();
                Value::Bool(match v.as_str() {
                    "yes" | "true" | "on" => true,
                    "no" | "false" | "off" => false,
                    _ => return Err(error("KeyError", repr(&v))),
                })
            }
            "int" => Value::Number(if self.unique {
                unique_integer(&self.scalar(node, text)?)?
            } else {
                integer(&self.scalar(node, text)?)?
            }),
            "float" => Value::Number(float(&self.scalar(node, text)?)?),
            "timestamp" => Value::Other(timestamp(&self.scalar(node, text)?)?),
            "binary" => {
                let value = self.scalar(node, text)?;
                binary(&value).map_err(|m| {
                    marked(text, "ConstructorError", None, &m, self.nodes[node].mark)
                })?;
                Value::Other(format!("binary:{value}"))
            }
            "set" => {
                pending.push_back(node);
                Value::Other("set".into())
            }
            "omap" | "pairs" => {
                pending.push_back(node);
                Value::Sequence(vec![])
            }
            _ => {
                return Err(marked(
                    text,
                    "ConstructorError",
                    None,
                    &format!(
                        "could not determine a constructor for the tag '{}'",
                        self.nodes[node].tag
                    ),
                    self.nodes[node].mark,
                ));
            }
        };
        self.nodes[node].value = Some(value);
        Ok(())
    }
    fn flatten(&mut self, node: usize, text: &str, depth: usize) -> Result<Vec<usize>> {
        if depth > 900 {
            return Err(error("RecursionError", "maximum recursion depth exceeded"));
        }
        let Raw::Mapping(mut entries) = self.nodes[node].raw.clone() else {
            unreachable!()
        };
        let mut merged = vec![];
        let mut index = 0;
        while index < entries.len() {
            let key = entries[index];
            let value = entries[index + 1];
            if self.nodes[key].tag == format!("{TAG}merge") {
                entries.drain(index..index + 2);
                self.nodes[node].raw = Raw::Mapping(entries.clone());
                match self.nodes[value].raw.clone() {
                    Raw::Mapping(_) => merged.extend(self.flatten(value, text, depth + 1)?),
                    Raw::Sequence(v) => {
                        let mut groups = vec![];
                        for id in v {
                            if !matches!(self.nodes[id].raw, Raw::Mapping(_)) {
                                return Err(marked(
                                    text,
                                    "ConstructorError",
                                    Some(("while constructing a mapping", self.nodes[node].mark)),
                                    &format!(
                                        "expected a mapping for merging, but found {}",
                                        self.nodes[id].raw.kind()
                                    ),
                                    self.nodes[id].mark,
                                ));
                            }
                            groups.push(self.flatten(id, text, depth + 1)?);
                        }
                        for group in groups.into_iter().rev() {
                            merged.extend(group);
                        }
                    }
                    _ => {
                        return Err(marked(
                            text,
                            "ConstructorError",
                            Some(("while constructing a mapping", self.nodes[node].mark)),
                            &format!(
                                "expected a mapping or list of mappings for merging, but found {}",
                                self.nodes[value].raw.kind()
                            ),
                            self.nodes[value].mark,
                        ));
                    }
                }
            } else {
                if self.nodes[key].tag == format!("{TAG}value") {
                    self.nodes[key].tag = format!("{TAG}str");
                }
                index += 2;
            }
        }
        merged.extend(entries);
        self.nodes[node].raw = Raw::Mapping(merged.clone());
        Ok(merged)
    }
    fn equal_key(&self, left: usize, right: usize) -> bool {
        if left == right {
            return true;
        }
        match (self.value(left), self.value(right)) {
            (Value::Null, Value::Null) => true,
            (Value::Str(a), Value::Str(b)) => a == b,
            (Value::Codepoints { points: a, .. }, Value::Codepoints { points: b, .. }) => a == b,
            (Value::Bool(a), Value::Bool(b)) => a == b,
            (Value::Number(a), Value::Number(b)) => a == b && a != "nan",
            (Value::Bool(a), Value::Number(b)) | (Value::Number(b), Value::Bool(a)) => {
                b == if *a { "1" } else { "0" }
            }
            (Value::Other(a), Value::Other(b)) => a == b,
            _ => false,
        }
    }
    fn unique_key_equal(&self, left: usize, right: usize) -> bool {
        if self.equal_key(left, right) {
            return true;
        }
        let number = |node| match self.value(node) {
            Value::Bool(v) => Some(if *v { "1" } else { "0" }),
            Value::Number(v) => Some(v.as_str()),
            _ => None,
        };
        match (number(left), number(right)) {
            // SafeConstructor reuses its singleton nan_value, hence Python
            // dict's identity shortcut treats repeated .nan keys as equal.
            (Some("nan"), Some("nan")) => true,
            (Some(a), Some(b)) if self.tag(left) == "float" && self.tag(right) == "float" => a
                .parse::<f64>()
                .ok()
                .zip(b.parse::<f64>().ok())
                .is_some_and(|(a, b)| a == b),
            (Some(a), Some(b)) if self.tag(left) == "float" || self.tag(right) == "float" => {
                let (float, integer) = if self.tag(left) == "float" {
                    (a, b)
                } else {
                    (b, a)
                };
                float.parse::<f64>().ok().is_some_and(|f| {
                    f.is_finite()
                        && f.fract() == 0.0
                        && (if f == 0.0 {
                            "0".into()
                        } else {
                            format!("{f:.0}")
                        }) == integer
                })
            }
            _ => false,
        }
    }
    pub(crate) fn scalar_display(&self, node: usize) -> String {
        match self.value(node) {
            Value::Null => "None".into(),
            Value::Bool(true) => "True".into(),
            Value::Bool(false) => "False".into(),
            Value::Str(s) => s.clone(),
            Value::Codepoints { transport, .. } => transport.clone(),
            Value::Number(s) if self.tag(node) == "float" => s
                .parse::<f64>()
                .map(crate::json_text::python_float_repr)
                .unwrap_or_else(|_| s.clone()),
            Value::Number(s) | Value::Other(s) => s.clone(),
            Value::Sequence(_) => "[]".into(),
            Value::Mapping(_) => "{}".into(),
        }
    }
    fn fill(&mut self, node: usize, text: &str, pending: &mut VecDeque<usize>) -> Result<()> {
        let tag = self.nodes[node].tag.strip_prefix(TAG).unwrap_or("");
        if matches!(tag, "pairs" | "omap") {
            let context = if tag == "omap" {
                "while constructing an ordered map"
            } else {
                "while constructing pairs"
            };
            let original = self.nodes[node].mark;
            let Raw::Sequence(children) = self.nodes[node].raw.clone() else {
                return Err(marked(
                    text,
                    "ConstructorError",
                    Some((context, original)),
                    &format!(
                        "expected a sequence, but found {}",
                        self.nodes[node].raw.kind()
                    ),
                    original,
                ));
            };
            let mut pairs = vec![];
            for child in children {
                let Raw::Mapping(pair) = self.nodes[child].raw.clone() else {
                    return Err(marked(
                        text,
                        "ConstructorError",
                        Some((context, original)),
                        &format!(
                            "expected a mapping of length 1, but found {}",
                            self.nodes[child].raw.kind()
                        ),
                        self.nodes[child].mark,
                    ));
                };
                if pair.len() != 2 {
                    return Err(marked(
                        text,
                        "ConstructorError",
                        Some((context, original)),
                        &format!(
                            "expected a single mapping item, but found {} items",
                            pair.len() / 2
                        ),
                        self.nodes[child].mark,
                    ));
                }
                for id in pair {
                    self.construct(id, text, pending)?;
                }
                // A pair is a Python tuple, hence not an admitted package string.
                pairs.push(self.add(Value::Other("tuple".into())));
            }
            self.nodes[node].value = Some(Value::Sequence(pairs));
            return Ok(());
        }
        if tag == "seq" && !matches!(self.nodes[node].raw, Raw::Sequence(_)) {
            return Err(self.expected(node, text, "sequence"));
        }
        if matches!(tag, "map" | "set") && !matches!(self.nodes[node].raw, Raw::Mapping(_)) {
            return Err(self.expected(node, text, "mapping"));
        }
        match self.nodes[node].raw.clone() {
            Raw::Sequence(children) => {
                for id in &children {
                    self.construct(*id, text, pending)?;
                }
                self.nodes[node].value = Some(Value::Sequence(children));
            }
            Raw::Mapping(_) => {
                let children = self.flatten(node, text, 0)?;
                let mut entries: Vec<(usize, usize)> = vec![];
                for pair in children.chunks_exact(2) {
                    let (key, value) = (pair[0], pair[1]);
                    self.construct(key, text, pending)?;
                    if (matches!(self.value(key), Value::Sequence(_) | Value::Mapping(_))
                        || matches!(self.value(key),Value::Other(k) if k=="set"||k=="pairs"||k=="omap"))
                    {
                        return Err(marked(
                            text,
                            "ConstructorError",
                            Some(("while constructing a mapping", self.nodes[node].mark)),
                            "found unhashable key",
                            self.nodes[key].mark,
                        ));
                    }
                    self.construct(value, text, pending)?;
                    if let Some(entry) = entries.iter_mut().find(|(k, _)| self.equal_key(*k, key)) {
                        entry.1 = value;
                    } else {
                        entries.push((key, value));
                    }
                }
                if matches!(self.value(node), Value::Mapping(_)) {
                    self.nodes[node].value = Some(Value::Mapping(entries));
                }
            }
            _ => unreachable!(),
        }
        Ok(())
    }
    pub fn dump(&self) -> Result<String> {
        fn anchor(
            yaml: &Yaml,
            node: usize,
            seen: &mut HashSet<usize>,
            names: &mut HashMap<usize, String>,
            next: &mut usize,
        ) {
            if !matches!(yaml.value(node), Value::Sequence(_) | Value::Mapping(_)) {
                return;
            }
            if !seen.insert(node) {
                if let std::collections::hash_map::Entry::Vacant(e) = names.entry(node) {
                    *next += 1;
                    e.insert(format!("id{:03}", *next));
                }
                return;
            }
            match yaml.value(node) {
                Value::Sequence(v) => {
                    for id in v {
                        anchor(yaml, *id, seen, names, next)
                    }
                }
                Value::Mapping(v) => {
                    for (k, v) in v {
                        anchor(yaml, *k, seen, names, next);
                        anchor(yaml, *v, seen, names, next)
                    }
                }
                _ => {}
            }
        }
        fn emit(
            yaml: &Yaml,
            node: usize,
            names: &HashMap<usize, String>,
            seen: &mut HashSet<usize>,
            out: &mut Vec<Event>,
        ) {
            let anchor = names.get(&node).cloned().unwrap_or_default();
            if !anchor.is_empty() && !seen.insert(node) {
                out.push(Event {
                    kind: 5,
                    anchor,
                    ..Event::default()
                });
                return;
            }
            match yaml.value(node) {
                Value::Str(v) => out.push(Event {
                    kind: 6,
                    value: v.clone(),
                    tag: format!("{TAG}str"),
                    implicit: resolved_tag(v) == format!("{TAG}str"),
                    ..Event::default()
                }),
                Value::Sequence(v) => {
                    out.push(Event {
                        kind: 7,
                        anchor,
                        ..Event::default()
                    });
                    for id in v {
                        emit(yaml, *id, names, seen, out)
                    }
                    out.push(Event {
                        kind: 8,
                        ..Event::default()
                    });
                }
                Value::Mapping(v) => {
                    out.push(Event {
                        kind: 9,
                        anchor,
                        ..Event::default()
                    });
                    for (k, v) in v {
                        emit(yaml, *k, names, seen, out);
                        emit(yaml, *v, names, seen, out)
                    }
                    out.push(Event {
                        kind: 10,
                        ..Event::default()
                    });
                }
                _ => unreachable!("only validated string packages are emitted"),
            }
        }
        let mut names = HashMap::new();
        anchor(self, self.root, &mut HashSet::new(), &mut names, &mut 0);
        let mut events = vec![
            Event {
                kind: 1,
                ..Event::default()
            },
            Event {
                kind: 3,
                ..Event::default()
            },
        ];
        emit(self, self.root, &names, &mut HashSet::new(), &mut events);
        events.extend([
            Event {
                kind: 4,
                ..Event::default()
            },
            Event {
                kind: 2,
                ..Event::default()
            },
        ]);
        ffi::emit(&events)
    }
}

fn unique_integer(value: &str) -> Result<String> {
    // Run the established constructor first for its ordinary admission/error
    // order, then retain arbitrary-width integer values in the unique-loader
    // profile. No package-export/SafeLoader default is changed.
    let ordinary = integer(value)?;
    let clean = value.replace('_', "");
    let negative = clean.starts_with('-');
    let body = clean.trim_start_matches(['+', '-']);
    let (base, digits) = if let Some(v) = body.strip_prefix("0b") {
        (2, v)
    } else if let Some(v) = body.strip_prefix("0x") {
        (16, v)
    } else if body.starts_with('0') && (body.len() > 1) {
        (8, body)
    } else {
        (10, body)
    };
    if digits.contains(':') || !digits.is_ascii() {
        return Ok(ordinary);
    }
    let mut decimal = vec![0u32];
    for c in digits.chars() {
        let Some(mut carry) = c.to_digit(base) else {
            return Ok(ordinary);
        };
        for digit in &mut decimal {
            let n = *digit * base + carry;
            *digit = n % 10;
            carry = n / 10;
        }
        while carry > 0 {
            decimal.push(carry % 10);
            carry /= 10;
        }
    }
    while decimal.len() > 1 && decimal.last() == Some(&0) {
        decimal.pop();
    }
    let mut text = String::new();
    if negative && decimal != [0] {
        text.push('-');
    }
    text.extend(decimal.iter().rev().map(|n| char::from(b'0' + *n as u8)));
    Ok(text)
}
fn integer(value: &str) -> Result<String> {
    let cleaned = value.replace('_', "");
    let mut body = cleaned.as_str();
    let negative = body.starts_with('-');
    if body.is_empty() {
        return Err(error("IndexError", "string index out of range"));
    }
    if body.starts_with(['+', '-']) {
        body = &body[1..];
    }
    if body.is_empty() {
        return Err(error("IndexError", "string index out of range"));
    }
    let (base, digits) = if body == "0" {
        (10, body)
    } else if let Some(v) = body.strip_prefix("0b") {
        (2, v)
    } else if let Some(v) = body.strip_prefix("0x") {
        (16, v)
    } else if body.starts_with('0') {
        (8, body)
    } else {
        (10, body)
    };
    if base == 10 && digits.contains(':') {
        let mut sum = 0i128;
        for part in digits.split(':') {
            let Some(n) = crate::python_text::python_int(part) else {
                return Err(error(
                    "ValueError",
                    format!("invalid literal for int() with base 10: {}", repr(part)),
                ));
            };
            sum = sum.saturating_mul(60).saturating_add(n);
        }
        return Ok((if negative { -sum } else { sum }).to_string());
    }
    let valid = if base == 10 {
        crate::python_text::python_int_is_valid(digits)
    } else {
        let d = crate::python_text::python_strip(digits);
        !d.is_empty()
            && d.trim_start_matches(['+', '-'])
                .chars()
                .all(|c| c.to_digit(base).is_some())
    };
    if !valid {
        return Err(error(
            "ValueError",
            format!(
                "invalid literal for int() with base {base}: {}",
                repr(digits)
            ),
        ));
    }
    let parsed = if base == 10 {
        crate::python_text::python_int(digits)
    } else {
        i128::from_str_radix(crate::python_text::python_strip(digits), base).ok()
    };
    Ok(parsed
        .and_then(|v| if negative { v.checked_neg() } else { Some(v) })
        .map_or_else(
            || format!("{}{}", if negative { "-" } else { "" }, digits),
            |v| v.to_string(),
        ))
}
fn float(value: &str) -> Result<String> {
    let clean = value.replace('_', "").to_lowercase();
    let mut body = clean.as_str();
    let negative = body.starts_with('-');
    if body.is_empty() {
        return Err(error("IndexError", "string index out of range"));
    }
    if body.starts_with(['+', '-']) {
        body = &body[1..];
    }
    let mut n = match body {
        ".inf" => f64::INFINITY,
        ".nan" => f64::NAN,
        _ => {
            let mut total = 0.0;
            for part in body.split(':') {
                let Some(n) = crate::python_text::parse_float(part) else {
                    return Err(error(
                        "ValueError",
                        format!("could not convert string to float: {}", repr(part)),
                    ));
                };
                total = total * 60.0 + n;
            }
            total
        }
    };
    if negative {
        n = -n;
    }
    Ok(crate::json_text::python_float_repr(n))
}
fn timestamp(value: &str) -> Result<String> {
    // SafeConstructor's timestamp regex is ASCII and its datetime constructor
    // validates the timezone before the calendar/time fields.
    fn fields(value: &str) -> Option<(i32, i32, i32, Option<(i32, i32, i32)>, i32)> {
        fn number(bytes: &[u8], at: &mut usize, min: usize, max: usize) -> Option<i32> {
            let start = *at;
            let mut value = 0;
            while *at < bytes.len() && *at - start < max && bytes[*at].is_ascii_digit() {
                value = value * 10 + i32::from(bytes[*at] - b'0');
                *at += 1;
            }
            (*at - start >= min).then_some(value)
        }
        fn take(bytes: &[u8], at: &mut usize, want: u8) -> Option<()> {
            if bytes.get(*at) == Some(&want) {
                *at += 1;
                Some(())
            } else {
                None
            }
        }
        let bytes = value.as_bytes();
        let mut at = 0;
        let year = number(bytes, &mut at, 4, 4)?;
        take(bytes, &mut at, b'-')?;
        let month = number(bytes, &mut at, 1, 2)?;
        take(bytes, &mut at, b'-')?;
        let day = number(bytes, &mut at, 1, 2)?;
        if at == bytes.len() {
            return Some((year, month, day, None, 0));
        }
        if matches!(bytes.get(at), Some(b'T' | b't')) {
            at += 1;
        } else {
            let old = at;
            while matches!(bytes.get(at), Some(b' ' | b'\t')) {
                at += 1;
            }
            if old == at {
                return None;
            }
        }
        let hour = number(bytes, &mut at, 1, 2)?;
        take(bytes, &mut at, b':')?;
        let minute = number(bytes, &mut at, 2, 2)?;
        take(bytes, &mut at, b':')?;
        let second = number(bytes, &mut at, 2, 2)?;
        if bytes.get(at) == Some(&b'.') {
            at += 1;
            while bytes.get(at).is_some_and(u8::is_ascii_digit) {
                at += 1;
            }
        }
        let space = at;
        while matches!(bytes.get(at), Some(b' ' | b'\t')) {
            at += 1;
        }
        let mut offset = 0;
        if at == bytes.len() {
            if at != space {
                return None;
            }
        } else if bytes.get(at) == Some(&b'Z') {
            at += 1;
        } else {
            let sign = match bytes.get(at) {
                Some(b'+') => 1,
                Some(b'-') => -1,
                _ => return None,
            };
            at += 1;
            let hours = number(bytes, &mut at, 1, 2)?;
            let minutes = if bytes.get(at) == Some(&b':') {
                at += 1;
                number(bytes, &mut at, 2, 2)?
            } else {
                0
            };
            offset = sign * (hours * 3600 + minutes * 60);
        }
        (at == bytes.len()).then_some((year, month, day, Some((hour, minute, second)), offset))
    }
    let (year, month, day, time, offset) = fields(value).ok_or_else(|| {
        error(
            "AttributeError",
            "'NoneType' object has no attribute 'groupdict'",
        )
    })?;
    if offset.abs() >= 86400 {
        let days = offset.div_euclid(86400);
        let seconds = offset.rem_euclid(86400);
        let delta = if seconds == 0 {
            format!("datetime.timedelta(days={days})")
        } else {
            format!("datetime.timedelta(days={days}, seconds={seconds})")
        };
        return Err(error(
            "ValueError",
            format!(
                "offset must be a timedelta strictly between -timedelta(hours=24) and timedelta(hours=24), not {delta}."
            ),
        ));
    }
    if !(1..=9999).contains(&year) {
        return Err(error("ValueError", format!("year {year} is out of range")));
    }
    if !(1..=12).contains(&month) {
        return Err(error("ValueError", "month must be in 1..12"));
    }
    let leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
    let days = [
        31,
        if leap { 29 } else { 28 },
        31,
        30,
        31,
        30,
        31,
        31,
        30,
        31,
        30,
        31,
    ];
    if day < 1 || day > days[month as usize - 1] {
        return Err(error("ValueError", "day is out of range for month"));
    }
    if let Some((hour, minute, second)) = time {
        for (number, max, name) in [
            (hour, 23, "hour"),
            (minute, 59, "minute"),
            (second, 59, "second"),
        ] {
            if number > max {
                return Err(error("ValueError", format!("{name} must be in 0..{max}")));
            }
        }
    }
    Ok(value.into())
}

fn binary(value: &str) -> std::result::Result<(), String> {
    if let Some((position, c)) = value.chars().enumerate().find(|(_, c)| !c.is_ascii()) {
        return Err(format!(
            "failed to convert base64 data into ascii: 'ascii' codec can't encode character {} in position {position}: ordinal not in range(128)",
            if (c as u32) <= 0xff {
                format!("'\\x{:02x}'", c as u32)
            } else if (c as u32) <= 0xffff {
                format!("'\\u{:04x}'", c as u32)
            } else {
                format!("'\\U{:08x}'", c as u32)
            }
        ));
    }
    let data = value
        .bytes()
        .filter(|b| b.is_ascii_alphanumeric() || matches!(b, b'+' | b'/' | b'='))
        .collect::<Vec<_>>();
    let n = data.iter().take_while(|b| **b != b'=').count();
    let pad = data.len() - n;
    if n % 4 == 1 {
        return Err(format!(
            "failed to decode base64 data: Invalid base64-encoded string: number of data characters ({n}) cannot be 1 more than a multiple of 4"
        ));
    }
    if n % 4 != 0 && n % 4 + pad < 4 {
        return Err("failed to decode base64 data: Incorrect padding".into());
    }
    Ok(())
}

#[cfg(test)]
mod unique_loader_tests {
    use super::*;
    #[test]
    fn duplicate_admission_precedes_duplicate_value_construction() {
        let failure = load_unique("a: 1\na: !unknown bad\n")
            .err()
            .unwrap()
            .to_string();
        assert_eq!(failure, "ValueError: Duplicate configuration key: a");
        let ordinary = load("a: 1\na: 2\n").unwrap();
        assert_eq!(
            ordinary.scalar_display(ordinary.get(ordinary.root, "a").unwrap()),
            "2"
        );
    }
    #[test]
    fn immediate_mapping_recursion_differs_from_sequence_aliases() {
        assert!(
            load_unique("a: &a {self: *a}")
                .err()
                .unwrap()
                .to_string()
                .starts_with("ConstructorError: found unconstructable recursive node")
        );
        assert!(load_unique("a: &a [*a]").is_ok());
        assert!(
            load_unique("a: &a {x: 1}\nb: {<<: *a}")
                .err()
                .unwrap()
                .to_string()
                .contains("tag:yaml.org,2002:merge")
        );
    }
    #[test]
    fn numeric_key_equality_does_not_round_integer_keys() {
        assert!(load_unique("18014398509481985: a\n18014398509481984.0: b").is_ok());
        assert!(load_unique("18014398509481984: a\n18014398509481984.0: b").is_err());
        assert!(load_unique(".nan: a\n.nan: b").is_err());
        let y = load_unique("a: 0xffffffffffffffffffffffffffffffffffff").unwrap();
        assert_eq!(
            y.scalar_display(y.get(y.root, "a").unwrap()),
            "22300745198530623141535718272648361505980415"
        );
    }
    #[test]
    fn unique_escaped_codepoints_preserve_python_strings_and_default_loader() {
        let y = load_unique(
            r#"a: "\ud834\udd1e"
b: "\ud800"
c: '\ud800'
d: "\ue000"
e: "\U0000d800"
"#,
        )
        .unwrap();
        assert!(
            matches!(y.value(y.get(y.root, "a").unwrap()), Value::Codepoints { points: v, .. } if v == &[0xd834, 0xdd1e])
        );
        assert!(
            matches!(y.value(y.get(y.root, "b").unwrap()), Value::Codepoints { points: v, .. } if v == &[0xd800])
        );
        assert_eq!(y.string(y.get(y.root, "c").unwrap()), Some(r"\ud800"));
        assert_eq!(y.string(y.get(y.root, "d").unwrap()), Some("\u{e000}"));
        assert!(
            matches!(y.value(y.get(y.root, "e").unwrap()), Value::Codepoints { points: v, .. } if v == &[0xd800])
        );
        assert!(load(r#"a: "\ud800""#).is_err());
        assert!(load_unique(r#"{"\ud800": 1, "\ud800": 2}"#).is_err());
        assert!(load_unique(r#"{"\ud834\udd1e": 1, "𝄞": 2}"#).is_ok());
    }
}
