//! Lossless receipt values for BioMart's explicit JSON publication procedures.
//! The shared reader validates the complete text first. This adapter retains
//! Python integer digits and Unicode code points (including lone surrogates),
//! which an equality-only metadata projection must never publish. Existing
//! readers and writers keep their original representation and behavior.
use super::*;
use std::fmt::Write as _;

#[derive(Clone, Debug)]
enum Node {
    Atom(String),
    Array(Vec<usize>),
    Object(Vec<(Vec<u32>, usize)>),
}

#[derive(Debug)]
pub(super) struct Record {
    nodes: Vec<Node>,
    root: usize,
}

enum Frame {
    Array(Vec<usize>),
    Object(Vec<(Vec<u32>, usize)>, Vec<u32>),
}

fn whitespace(raw: &str, at: &mut usize) {
    while raw.as_bytes().get(*at).is_some_and(u8::is_ascii_whitespace) {
        *at += 1;
    }
}
fn quoted_end(raw: &str, at: &mut usize) -> usize {
    let start = *at;
    *at += 1;
    loop {
        let b = raw.as_bytes()[*at];
        *at += 1;
        if b == b'"' {
            return start;
        }
        if b == b'\\' {
            *at += 1;
        }
    }
}
fn codepoints(raw: &str) -> Vec<u32> {
    let mut result = Vec::new();
    let mut at = 1;
    while at < raw.len() - 1 {
        let c = raw[at..].chars().next().unwrap();
        at += c.len_utf8();
        if c != '\\' {
            result.push(c as u32);
            continue;
        }
        let escape = raw.as_bytes()[at];
        at += 1;
        let mut value = match escape {
            b'"' | b'\\' | b'/' => escape as u32,
            b'b' => 8,
            b'f' => 12,
            b'n' => 10,
            b'r' => 13,
            b't' => 9,
            b'u' => {
                let value = u32::from_str_radix(&raw[at..at + 4], 16).unwrap();
                at += 4;
                value
            }
            _ => unreachable!("the complete document was validated"),
        };
        if (0xd800..0xdc00).contains(&value) && raw[at..].starts_with("\\u") {
            let low = u32::from_str_radix(&raw[at + 2..at + 6], 16).unwrap();
            if (0xdc00..0xe000).contains(&low) {
                value = 0x10000 + ((value - 0xd800) << 10) + low - 0xdc00;
                at += 6;
            }
        }
        result.push(value);
    }
    result
}
fn quote(codes: &[u32]) -> String {
    let mut out = String::from("\"");
    for &c in codes {
        match c {
            34 => out.push_str("\\\""),
            92 => out.push_str("\\\\"),
            8 => out.push_str("\\b"),
            9 => out.push_str("\\t"),
            10 => out.push_str("\\n"),
            12 => out.push_str("\\f"),
            13 => out.push_str("\\r"),
            0x20..=0x7e => out.push(char::from_u32(c).unwrap()),
            0..=0xffff => {
                let _ = write!(out, "\\u{c:04x}");
            }
            _ => {
                let n = c - 0x10000;
                let _ = write!(
                    out,
                    "\\u{:04x}\\u{:04x}",
                    0xd800 + (n >> 10),
                    0xdc00 + (n & 1023)
                );
            }
        }
    }
    out.push('"');
    out
}

impl Record {
    pub(super) fn parse(raw: &str) -> Result<Self> {
        json_text::metadata_field_spans(raw).map_err(json_failure)?;
        let mut nodes = Vec::new();
        let mut frames = Vec::new();
        let mut at = 0;
        let mut pending = None;
        loop {
            if pending.is_none() {
                whitespace(raw, &mut at);
                let opening = raw.as_bytes()[at];
                let node = match opening {
                    b'[' | b'{' => {
                        at += 1;
                        whitespace(raw, &mut at);
                        if raw.as_bytes()[at] == if opening == b'[' { b']' } else { b'}' } {
                            at += 1;
                            if opening == b'[' {
                                Node::Array(Vec::new())
                            } else {
                                Node::Object(Vec::new())
                            }
                        } else {
                            if opening == b'[' {
                                frames.push(Frame::Array(Vec::new()));
                            } else {
                                let start = quoted_end(raw, &mut at);
                                let key = codepoints(&raw[start..at]);
                                whitespace(raw, &mut at);
                                at += 1; // validated colon
                                frames.push(Frame::Object(Vec::new(), key));
                            }
                            continue;
                        }
                    }
                    b'"' => {
                        let start = quoted_end(raw, &mut at);
                        Node::Atom(quote(&codepoints(&raw[start..at])))
                    }
                    _ => {
                        let start = at;
                        while raw
                            .as_bytes()
                            .get(at)
                            .is_some_and(|b| !b.is_ascii_whitespace() && !b",]}".contains(b))
                        {
                            at += 1;
                        }
                        let literal = &raw[start..at];
                        let encoded = if literal.bytes().any(|b| b".eENI".contains(&b)) {
                            // Use the existing qualified Python float formatter.
                            json_text::dumps_with_options(
                                &parse_json(literal)?,
                                json_text::DumpsOptions::default(),
                            )
                            .map_err(json_failure)?
                        } else if literal == "-0" {
                            "0".into()
                        } else {
                            literal.into()
                        };
                        Node::Atom(encoded)
                    }
                };
                nodes.push(node);
                pending = Some(nodes.len() - 1);
            }
            let child = pending.take().unwrap();
            let Some(frame) = frames.last_mut() else {
                return Ok(Self { nodes, root: child });
            };
            match frame {
                Frame::Array(items) => items.push(child),
                Frame::Object(items, key) => {
                    if let Some((_, old)) = items.iter_mut().find(|(k, _)| k == key) {
                        *old = child;
                    } else {
                        items.push((key.clone(), child));
                    }
                }
            }
            whitespace(raw, &mut at);
            let delimiter = raw.as_bytes()[at];
            at += 1;
            if delimiter == b',' {
                if let Frame::Object(_, key) = frame {
                    whitespace(raw, &mut at);
                    let start = quoted_end(raw, &mut at);
                    *key = codepoints(&raw[start..at]);
                    whitespace(raw, &mut at);
                    at += 1;
                }
                continue;
            }
            nodes.push(match frames.pop().unwrap() {
                Frame::Array(items) => Node::Array(items),
                Frame::Object(items, _) => Node::Object(items),
            });
            pending = Some(nodes.len() - 1);
        }
    }

    pub(super) fn portable(mut self) -> Result<Self> {
        let mut items = Vec::new();
        let names = get(facts(), "ftp_specification")
            .as_object()
            .unwrap()
            .iter()
            .map(|(k, _)| k.as_str())
            .chain(["bytes", "sha256", "statistics"]);
        for name in names {
            let key = name.chars().map(u32::from).collect::<Vec<_>>();
            let value = match &self.nodes[self.root] {
                Node::Object(fields) => fields
                    .iter()
                    .find(|(k, _)| k == &key)
                    .map(|(_, v)| *v)
                    .ok_or_else(|| error("KeyError", repr(name)))?,
                Node::Array(_) => {
                    return Err(error(
                        "TypeError",
                        "list indices must be integers or slices, not str",
                    ));
                }
                Node::Atom(v) => {
                    return Err(error(
                        "TypeError",
                        if v.starts_with('"') {
                            "string indices must be integers, not 'str'".into()
                        } else {
                            format!(
                                "'{}' object is not subscriptable",
                                if v == "null" {
                                    "NoneType"
                                } else if matches!(v.as_str(), "true" | "false") {
                                    "bool"
                                } else if v.bytes().any(|b| b".eENI".contains(&b)) {
                                    "float"
                                } else {
                                    "int"
                                }
                            )
                        },
                    ));
                }
            };
            items.push((key, value));
        }
        for (key, value) in [
            ("metadata_kind", "\"priva_bundled_mapping_snapshot\""),
            ("metadata_schema_version", "1"),
        ] {
            self.nodes.push(Node::Atom(value.into()));
            items.push((key.chars().map(u32::from).collect(), self.nodes.len() - 1));
        }
        self.nodes.push(Node::Object(items));
        self.root = self.nodes.len() - 1;
        Ok(self)
    }

    /// Stream json.dump's chunks; callers retain their own open/close/fsync and
    /// publication order. The arena and explicit stack avoid a native stack
    /// overflow while handling source-admitted nested JSON.
    pub(super) fn write(
        &self,
        indent: bool,
        sorted: bool,
        mut output: impl FnMut(&str) -> Result<()>,
    ) -> Result<()> {
        enum Task {
            Node(usize, usize),
            Text(String),
        }
        let mut tasks = vec![Task::Node(self.root, 0)];
        while let Some(task) = tasks.pop() {
            let (index, level) = match task {
                Task::Node(i, l) => (i, l),
                Task::Text(s) => {
                    output(&s)?;
                    continue;
                }
            };
            let (opening, closing, items) = match &self.nodes[index] {
                Node::Atom(s) => {
                    output(s)?;
                    continue;
                }
                Node::Array(items) => (
                    '[',
                    ']',
                    items.iter().map(|i| (None, *i)).collect::<Vec<_>>(),
                ),
                Node::Object(fields) => {
                    let mut ordered = fields.iter().collect::<Vec<_>>();
                    if sorted {
                        ordered.sort_by(|a, b| a.0.cmp(&b.0));
                    }
                    (
                        '{',
                        '}',
                        ordered
                            .into_iter()
                            .map(|(k, i)| (Some(quote(k)), *i))
                            .collect(),
                    )
                }
            };
            if items.is_empty() {
                output(&format!("{opening}{closing}"))?;
                continue;
            }
            output(&opening.to_string())?;
            tasks.push(Task::Text(if indent {
                format!("\n{}{closing}", "  ".repeat(level))
            } else {
                closing.to_string()
            }));
            for (n, (key, child)) in items.into_iter().enumerate().rev() {
                tasks.push(Task::Node(child, level + 1));
                let mut prefix = if n == 0 {
                    String::new()
                } else if indent {
                    ",".into()
                } else {
                    ", ".into()
                };
                if indent {
                    prefix.push('\n');
                    prefix.push_str(&"  ".repeat(level + 1));
                }
                if let Some(key) = key {
                    prefix.push_str(&key);
                    prefix.push_str(": ");
                }
                tasks.push(Task::Text(prefix));
            }
        }
        Ok(())
    }
    pub(super) fn text(&self) -> Result<String> {
        let mut out = String::new();
        self.write(false, false, |s| {
            out.push_str(s);
            Ok(())
        })?;
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn dump(raw: &str) -> String {
        let r = Record::parse(raw).unwrap();
        let mut out = String::new();
        r.write(true, true, |s| {
            out.push_str(s);
            Ok(())
        })
        .unwrap();
        out
    }
    #[test]
    fn lossless_receipts_keep_wide_integers_and_distinct_surrogate_keys() {
        let text = r#"{"\ud800": 1, "\u0000": 2, "wide": 1461501637330902918203684832716283019655932542976, "\ud800": 3}"#;
        assert_eq!(
            dump(text),
            "{\n  \"\\u0000\": 2,\n  \"wide\": 1461501637330902918203684832716283019655932542976,\n  \"\\ud800\": 3\n}"
        );
    }
    #[test]
    fn unicode_sorting_and_scalar_formatting_match_python() {
        assert_eq!(
            dump(r#"["\ud83d\ude00", "\ud800x", "\ud800\ud800", -0, -0.0, 1e-5, 1e999, NaN]"#),
            "[\n  \"\\ud83d\\ude00\",\n  \"\\ud800x\",\n  \"\\ud800\\ud800\",\n  0,\n  -0.0,\n  1e-05,\n  Infinity,\n  NaN\n]"
        );
        assert_eq!(
            dump(r#"{"😀":0,"\ue000":1,"\ud800":2}"#),
            "{\n  \"\\ud800\": 2,\n  \"\\ue000\": 1,\n  \"\\ud83d\\ude00\": 0\n}"
        );
    }
    #[test]
    fn no_projection_value_reaches_publication() {
        let raw = format!("{{\"a\":{},\"b\":\"\\udfff\"}}", "9".repeat(4300));
        let text = Record::parse(&raw).unwrap().text().unwrap();
        assert!(text.contains(&"9".repeat(4300)));
        assert!(text.contains("\\udfff"));
        assert!(
            Record::parse(&format!("[{}]", "9".repeat(4301)))
                .unwrap_err()
                .to_string()
                .starts_with("ValueError:")
        );
    }
}
