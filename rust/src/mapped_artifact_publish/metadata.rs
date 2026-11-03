//! Lossless, read-only views of publication metadata. Unknown values are
//! validated in full and never materialized as bounded or replacement values.
use super::{error, require};
use crate::cli::Result;
use crate::json_text::{self, Json};
use crate::python_path::repr;
use std::ops::Range;

#[derive(Clone, Copy)]
pub(super) struct Value<'a>(pub(super) &'a str);
impl<'a> Value<'a> {
    pub(super) fn validate(text: &'a str) -> Result<Self> {
        json_text::metadata_field_spans(text).map_err(failure)?;
        Ok(Self(text.trim()))
    }
    pub(super) fn fields(self) -> Result<Option<Vec<(String, Range<usize>)>>> {
        json_text::metadata_field_spans(self.0).map_err(failure)
    }
    pub(super) fn object(self) -> bool {
        self.0.starts_with('{')
    }
    pub(super) fn get(self, key: &str) -> Result<Self> {
        if !self.object() {
            return Err(error(
                "AttributeError",
                format!("'{}' object has no attribute 'get'", self.kind()),
            ));
        }
        Ok(self.field(key)?.unwrap_or(Value("null")))
    }
    pub(super) fn item(self, key: &str) -> Result<Self> {
        if !self.object() {
            return Err(error(
                "TypeError",
                if self.0.starts_with('[') {
                    "list indices must be integers or slices, not str".to_owned()
                } else if self.0.starts_with('"') {
                    "string indices must be integers, not 'str'".to_owned()
                } else {
                    format!("'{}' object is not subscriptable", self.kind())
                },
            ));
        }
        self.field(key)?.ok_or_else(|| error("KeyError", repr(key)))
    }
    fn field(self, key: &str) -> Result<Option<Self>> {
        Ok(self
            .fields()?
            .and_then(|fields| fields.into_iter().find(|(k, _)| k == key))
            .map(|(_, span)| Value(&self.0[span])))
    }
    pub(super) fn keys(self, expected: &[&str]) -> Result<bool> {
        if let Some(fields) = self.fields()? {
            return Ok(fields.len() == expected.len()
                && fields
                    .iter()
                    .all(|(key, _)| expected.contains(&key.as_str())));
        }
        // set(value) is deliberately not a dictionary-only check in the source.
        if self.0.starts_with('[') {
            let elements = elements(self.0);
            let mut keys = Vec::new();
            let mut unmatched = false;
            for item in elements {
                let value = Value(item);
                if matches!(value.kind(), "list" | "dict") {
                    return Err(error(
                        "TypeError",
                        format!("unhashable type: '{}'", value.kind()),
                    ));
                }
                let mut matched = None;
                for key in expected {
                    if value.equal(&Json::Str((*key).to_owned()))? {
                        matched = Some(*key);
                        break;
                    }
                }
                if let Some(key) = matched {
                    if !keys.contains(&key) {
                        keys.push(key);
                    }
                } else {
                    unmatched = true;
                }
            }
            return Ok(!unmatched && keys.len() == expected.len());
        }
        if self.0.starts_with('"') {
            return Ok(false);
        }
        Err(error(
            "TypeError",
            format!("'{}' object is not iterable", self.kind()),
        ))
    }
    pub(super) fn exact_object(self, keys: &[&str]) -> Result<bool> {
        Ok(self.object() && self.keys(keys)?)
    }
    pub(super) fn equal(self, target: &Json) -> Result<bool> {
        Ok(json_text::py_equal(
            &json_text::parse_metadata_comparison(self.0, target).map_err(failure)?,
            target,
        ))
    }
    pub(super) fn json(self) -> Result<Json> {
        json_text::parse(self.0, false).map_err(failure)
    }
    pub(super) fn truth(self) -> bool {
        match self.0.as_bytes().first() {
            Some(b'n' | b'f') => false,
            Some(b'"') => self.0 != "\"\"",
            Some(b'{' | b'[') => !self.0[1..self.0.len() - 1].trim().is_empty(),
            Some(b't') => true,
            _ => number(self.0).is_none_or(|n| n != 0.0),
        }
    }
    pub(super) fn positive(self) -> Result<bool> {
        match self.kind() {
            "int" => Ok(!self.0.starts_with('-') && self.0.bytes().any(|b| b != b'0')),
            "float" => Ok(number(self.0).is_some_and(|n| n > 0.0)),
            "bool" => Ok(self.0 == "true"),
            other => Err(error(
                "TypeError",
                format!("'>' not supported between instances of '{other}' and 'int'"),
            )),
        }
    }
    pub(super) fn kind(self) -> &'static str {
        match self.0.as_bytes().first() {
            Some(b'n') => "NoneType",
            Some(b't' | b'f') => "bool",
            Some(b'"') => "str",
            Some(b'{') => "dict",
            Some(b'[') => "list",
            _ if self
                .0
                .bytes()
                .any(|b| matches!(b, b'.' | b'e' | b'E' | b'N' | b'I')) =>
            {
                "float"
            }
            _ => "int",
        }
    }
}
fn number(raw: &str) -> Option<f64> {
    match raw {
        "NaN" => Some(f64::NAN),
        "Infinity" => Some(f64::INFINITY),
        "-Infinity" => Some(f64::NEG_INFINITY),
        _ => raw.parse().ok(),
    }
}
pub(super) fn failure(message: String) -> Box<dyn std::error::Error> {
    error(
        if message.starts_with("Exceeds the limit (4300 digits)") {
            "ValueError"
        } else if json_text::is_refusal(&message) {
            "UnsupportedMetadata"
        } else {
            "JSONDecodeError"
        },
        message,
    )
}
/// Delimiters only; the complete document has already passed the JSON reader.
fn elements(raw: &str) -> Vec<&str> {
    let inner = &raw[1..raw.len() - 1];
    if inner.trim().is_empty() {
        return Vec::new();
    }
    let (mut start, mut depth, mut quoted, mut escaped) = (0, 0, false, false);
    let mut result = Vec::new();
    for (index, b) in inner.bytes().enumerate() {
        if quoted {
            if escaped {
                escaped = false;
            } else if b == b'\\' {
                escaped = true;
            } else if b == b'"' {
                quoted = false;
            }
        } else {
            match b {
                b'"' => quoted = true,
                b'[' | b'{' => depth += 1,
                b']' | b'}' => depth -= 1,
                b',' if depth == 0 => {
                    result.push(inner[start..index].trim());
                    start = index + 1;
                }
                _ => {}
            }
        }
    }
    result.push(inner[start..].trim());
    result
}
pub(super) fn require_equal(raw: Value<'_>, expected: &Json, message: &str) -> Result<()> {
    require(raw.equal(expected)?, message)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn set_shape_checks_finish_iteration_before_equality() {
        let keys = ["artifact", "receipt"];
        assert!(
            Value::validate(r#"["artifact","receipt","artifact"]"#)
                .unwrap()
                .keys(&keys)
                .unwrap()
        );
        let error = Value::validate(r#"["unknown",[]]"#)
            .unwrap()
            .keys(&keys)
            .unwrap_err();
        assert_eq!(error.to_string(), "TypeError: unhashable type: 'list'");
    }
    #[test]
    fn ignored_values_remain_validated_without_bounded_projection() {
        let raw = r#" {"selected":1,"ignored":{"wide":1461501637330902918203684832716283019655932542976,"text":"\ud800"}} "#;
        let value = Value::validate(raw).unwrap();
        assert!(value.get("selected").unwrap().equal(&Json::Int(1)).unwrap());
        assert!(value.get("ignored").unwrap().truth());
        assert!(Value::validate(r#"{"selected":1,"ignored":[1,]}"#).is_err());
    }
}
