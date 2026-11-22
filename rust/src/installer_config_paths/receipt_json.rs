//! Lossless Python-JSON values for the adoption record, which is read, extended
//! and re-emitted. Ignored values cannot be replaced by equality sentinels.
//! The shared complete-text validator retains ownership of syntax admission;
//! ordinary Json reader defaults and scientific cache representations are not
//! widened. Integers remain decimal text and strings retain Unicode code points.
use super::*;
#[derive(Clone, Debug)]
pub(super) enum RecordJson {
    Null,
    Bool(bool),
    Int(String),
    Float(f64),
    Str(Vec<u32>),
    Array(Vec<Self>),
    Object(Vec<(Vec<u32>, Self)>),
}
impl RecordJson {
    pub fn from_json(value: &Json) -> Self {
        match value {
            Json::Null => Self::Null,
            Json::Bool(v) => Self::Bool(*v),
            Json::Int(v) => Self::Int(v.to_string()),
            Json::Float(v) => Self::Float(*v),
            Json::Str(v) => Self::Str(crate::python_unicode::points(v)),
            Json::Array(v) => Self::Array(v.iter().map(Self::from_json).collect()),
            Json::Object(v) => Self::Object(
                v.iter()
                    .map(|(k, v)| (crate::python_unicode::points(k), Self::from_json(v)))
                    .collect(),
            ),
        }
    }
    pub fn parse(text: &str) -> Result<Self> {
        json_text::metadata_field_spans(text).map_err(|e| {
            error(
                if e.starts_with("Exceeds the limit") {
                    "ValueError"
                } else {
                    "JSONDecodeError"
                },
                e,
            )
        })?;
        struct Parser<'a> {
            text: &'a str,
            at: usize,
        }
        impl Parser<'_> {
            fn ws(&mut self) {
                while self
                    .text
                    .as_bytes()
                    .get(self.at)
                    .is_some_and(|b| b" \r\n\t".contains(b))
                {
                    self.at += 1;
                }
            }
            fn string(&mut self) -> Vec<u32> {
                self.at += 1;
                let mut result = vec![];
                while self.text.as_bytes()[self.at] != b'"' {
                    if self.text.as_bytes()[self.at] != b'\\' {
                        let c = self.text[self.at..].chars().next().unwrap();
                        self.at += c.len_utf8();
                        result.push(c as u32);
                        continue;
                    }
                    self.at += 1;
                    let c = self.text.as_bytes()[self.at];
                    self.at += 1;
                    match c {
                        b'u' => {
                            let mut n =
                                u32::from_str_radix(&self.text[self.at..self.at + 4], 16).unwrap();
                            self.at += 4;
                            if (0xd800..=0xdbff).contains(&n)
                                && self.text[self.at..].starts_with("\\u")
                            {
                                let low =
                                    u32::from_str_radix(&self.text[self.at + 2..self.at + 6], 16)
                                        .unwrap();
                                if (0xdc00..=0xdfff).contains(&low) {
                                    self.at += 6;
                                    n = 0x10000 + (n - 0xd800) * 1024 + low - 0xdc00;
                                }
                            }
                            result.push(n);
                        }
                        b'n' => result.push(10),
                        b'r' => result.push(13),
                        b't' => result.push(9),
                        b'b' => result.push(8),
                        b'f' => result.push(12),
                        _ => result.push(c as u32),
                    }
                }
                self.at += 1;
                result
            }
            fn value(&mut self, depth: usize) -> Result<RecordJson> {
                if depth > 990 {
                    return Err(error(
                        "RecursionError",
                        "maximum recursion depth exceeded while decoding a JSON object from a unicode string",
                    ));
                }
                self.ws();
                Ok(match self.text.as_bytes()[self.at] {
                    b'"' => RecordJson::Str(self.string()),
                    b'[' => {
                        self.at += 1;
                        self.ws();
                        let mut values = vec![];
                        while self.text.as_bytes()[self.at] != b']' {
                            values.push(self.value(depth + 1)?);
                            self.ws();
                            if self.text.as_bytes()[self.at] == b',' {
                                self.at += 1;
                                self.ws();
                            } else {
                                break;
                            }
                        }
                        self.at += 1;
                        RecordJson::Array(values)
                    }
                    b'{' => {
                        self.at += 1;
                        self.ws();
                        let mut values: Vec<(Vec<u32>, RecordJson)> = vec![];
                        while self.text.as_bytes()[self.at] != b'}' {
                            let key = self.string();
                            self.ws();
                            self.at += 1;
                            let value = self.value(depth + 1)?;
                            if let Some((_, old)) = values.iter_mut().find(|(k, _)| *k == key) {
                                *old = value;
                            } else {
                                values.push((key, value));
                            }
                            self.ws();
                            if self.text.as_bytes()[self.at] == b',' {
                                self.at += 1;
                                self.ws();
                            } else {
                                break;
                            }
                        }
                        self.at += 1;
                        RecordJson::Object(values)
                    }
                    _ => {
                        let start = self.at;
                        while self
                            .text
                            .as_bytes()
                            .get(self.at)
                            .is_some_and(|b| !b" \t\r\n,]}".contains(b))
                        {
                            self.at += 1;
                        }
                        let token = &self.text[start..self.at];
                        match token {
                            "null" => RecordJson::Null,
                            "true" => RecordJson::Bool(true),
                            "false" => RecordJson::Bool(false),
                            "NaN" => RecordJson::Float(f64::NAN),
                            "Infinity" => RecordJson::Float(f64::INFINITY),
                            "-Infinity" => RecordJson::Float(f64::NEG_INFINITY),
                            _ if token.contains(['.', 'e', 'E']) => {
                                RecordJson::Float(token.parse().unwrap())
                            }
                            _ => RecordJson::Int(if token == "-0" {
                                "0".into()
                            } else {
                                token.into()
                            }),
                        }
                    }
                })
            }
        }
        Parser { text, at: 0 }.value(0)
    }
    pub fn kind(&self) -> &'static str {
        match self {
            Self::Null => "NoneType",
            Self::Bool(_) => "bool",
            Self::Int(_) => "int",
            Self::Float(_) => "float",
            Self::Str(_) => "str",
            Self::Array(_) => "list",
            Self::Object(_) => "dict",
        }
    }
    pub fn get_mut(&mut self, key: &str) -> Result<&mut Self> {
        let name = key.chars().map(u32::from).collect::<Vec<_>>();
        match self {
            Self::Object(v) => v
                .iter_mut()
                .find(|(k, _)| *k == name)
                .map(|(_, v)| v)
                .ok_or_else(|| error("KeyError", repr(key))),
            Self::Array(_) => Err(error(
                "TypeError",
                "list indices must be integers or slices, not str",
            )),
            Self::Str(_) => Err(error(
                "TypeError",
                "string indices must be integers, not 'str'",
            )),
            _ => Err(error(
                "TypeError",
                format!("'{}' object is not subscriptable", self.kind()),
            )),
        }
    }
    fn equivalent(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Null, Self::Null) => true,
            (Self::Bool(a), Self::Bool(b)) => a == b,
            (Self::Int(a), Self::Int(b)) => a == b,
            (Self::Float(a), Self::Float(b)) => a == b || (a.is_nan() && b.is_nan()),
            (Self::Str(a), Self::Str(b)) => a == b,
            (Self::Array(a), Self::Array(b)) => {
                a.len() == b.len() && a.iter().zip(b).all(|(a, b)| a.equivalent(b))
            }
            (Self::Object(a), Self::Object(b)) => {
                a.len() == b.len()
                    && a.iter()
                        .all(|(k, v)| b.iter().any(|(k2, v2)| k == k2 && v.equivalent(v2)))
            }
            (Self::Bool(a), Self::Int(b)) | (Self::Int(b), Self::Bool(a)) => {
                b == if *a { "1" } else { "0" }
            }
            (Self::Bool(a), Self::Float(b)) | (Self::Float(b), Self::Bool(a)) => {
                *b == if *a { 1.0 } else { 0.0 }
            }
            (Self::Int(a), Self::Float(b)) | (Self::Float(b), Self::Int(a)) => {
                b.is_finite()
                    && b.fract() == 0.0
                    && a.parse::<i128>()
                        .ok()
                        .is_some_and(|a| json_text::py_equal(&Json::Int(a), &Json::Float(*b)))
            }
            _ => false,
        }
    }
    pub fn extend_unique(&mut self, values: &[Json]) -> Result<()> {
        let kind = self.kind();
        let Self::Array(old) = self else {
            return Err(error(
                "AttributeError",
                format!("'{kind}' object has no attribute 'extend'"),
            ));
        };
        for value in values {
            let value = Self::from_json(value);
            if !old.iter().any(|v| v.equivalent(&value)) {
                old.push(value);
            }
        }
        Ok(())
    }
    pub fn dumps(&self, ascii: bool) -> Result<String> {
        fn quote(text: &[u32], ascii: bool, out: &mut String) -> Result<()> {
            out.push('"');
            for (at, &c) in text.iter().enumerate() {
                match c {
                    34 => out.push_str("\\\""),
                    92 => out.push_str("\\\\"),
                    8 => out.push_str("\\b"),
                    12 => out.push_str("\\f"),
                    10 => out.push_str("\\n"),
                    13 => out.push_str("\\r"),
                    9 => out.push_str("\\t"),
                    c if c < 32 || (ascii && c >= 127) => {
                        if c <= 0xffff {
                            out.push_str(&format!("\\u{c:04x}"));
                        } else {
                            let v = c - 0x10000;
                            out.push_str(&format!(
                                "\\u{:04x}\\u{:04x}",
                                0xd800 + (v >> 10),
                                0xdc00 + (v & 1023)
                            ));
                        }
                    }
                    c => {
                        let Some(c) = char::from_u32(c) else {
                            let position = out.chars().count();
                            let count = text[at..]
                                .iter()
                                .take_while(|v| (0xd800..=0xdfff).contains(*v))
                                .count();
                            return Err(error(
                                "UnicodeEncodeError",
                                if count == 1 {
                                    format!(
                                        "'utf-8' codec can't encode character '\\u{c:04x}' in position {position}: surrogates not allowed"
                                    )
                                } else {
                                    format!(
                                        "'utf-8' codec can't encode characters in position {position}-{}: surrogates not allowed",
                                        position + count - 1
                                    )
                                },
                            ));
                        };
                        out.push(c);
                    }
                }
            }
            out.push('"');
            Ok(())
        }
        fn value(v: &RecordJson, level: usize, ascii: bool, out: &mut String) -> Result<()> {
            match v {
                RecordJson::Null => out.push_str("null"),
                RecordJson::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
                RecordJson::Int(n) => out.push_str(n),
                RecordJson::Float(f) => out.push_str(
                    &json_text::dumps_with_options(
                        &Json::Float(*f),
                        json_text::DumpsOptions::default(),
                    )
                    .map_err(invalid)?,
                ),
                RecordJson::Str(s) => quote(s, ascii, out)?,
                RecordJson::Array(a) => {
                    out.push('[');
                    for (i, v) in a.iter().enumerate() {
                        out.push_str(if i == 0 { "\n" } else { ",\n" });
                        out.push_str(&" ".repeat((level + 1) * 2));
                        value(v, level + 1, ascii, out)?;
                    }
                    if !a.is_empty() {
                        out.push('\n');
                        out.push_str(&" ".repeat(level * 2));
                    }
                    out.push(']');
                }
                RecordJson::Object(a) => {
                    out.push('{');
                    for (i, (k, v)) in a.iter().enumerate() {
                        out.push_str(if i == 0 { "\n" } else { ",\n" });
                        out.push_str(&" ".repeat((level + 1) * 2));
                        quote(k, ascii, out)?;
                        out.push_str(": ");
                        value(v, level + 1, ascii, out)?;
                    }
                    if !a.is_empty() {
                        out.push('\n');
                        out.push_str(&" ".repeat(level * 2));
                    }
                    out.push('}');
                }
            }
            Ok(())
        }
        let mut out = String::new();
        value(self, 0, ascii, &mut out)?;
        Ok(out)
    }
    pub fn update(&mut self, updates: &Paths) -> Result<()> {
        let Self::Object(values) = self else {
            return Err(error(
                "AttributeError",
                format!("'{}' object has no attribute 'update'", self.kind()),
            ));
        };
        for (k, v) in updates {
            let key = crate::python_unicode::points(k);
            let value = Self::Str(crate::python_unicode::points(v));
            if let Some((_, old)) = values.iter_mut().find(|(k, _)| *k == key) {
                *old = value;
            } else {
                values.push((key, value));
            }
        }
        Ok(())
    }
}
