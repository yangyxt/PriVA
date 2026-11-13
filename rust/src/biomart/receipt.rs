//! Lossless receipt admission. The complete JSON document is validated by the
//! shared reader; the original shape checks inspect raw, validated fields.
//! Equality-only projections never become returned values or rewritten files.
use super::*;
use crate::python_path::is_file;
use std::ops::Range;

pub(crate) struct ArtifactReceipt {
    pub(crate) json: String,
    identity: Range<usize>,
}
impl ArtifactReceipt {
    pub(crate) fn matches(&self, expected: &Json) -> Result<bool> {
        matches_value(&self.json[self.identity.clone()], expected)
    }
}
struct Object<'a> {
    text: &'a str,
    fields: Vec<(String, Range<usize>)>,
}
impl<'a> Object<'a> {
    fn parse(text: &'a str) -> Result<Option<Self>> {
        Ok(json_text::metadata_field_spans(text)
            .map_err(json_failure)?
            .map(|fields| Self { text, fields }))
    }
    fn span(&self, key: &str) -> Option<Range<usize>> {
        self.fields
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.clone())
    }
    fn get(&self, key: &str) -> Option<&'a str> {
        self.span(key).map(|span| &self.text[span])
    }
    fn keys(&self, expected: &[&str]) -> bool {
        self.fields.len() == expected.len()
            && self
                .fields
                .iter()
                .all(|(key, _)| expected.contains(&key.as_str()))
    }
}
fn matches_value(raw: &str, expected: &Json) -> Result<bool> {
    Ok(equal(
        &json_text::parse_metadata_comparison(raw, expected).map_err(json_failure)?,
        expected,
    ))
}
fn integer(raw: &str) -> bool {
    let digits = raw.strip_prefix('-').unwrap_or(raw);
    !digits.is_empty() && digits.bytes().all(|c| c.is_ascii_digit())
}
fn size(raw: &str) -> bool {
    integer(raw) && (!raw.starts_with('-') || raw[1..].bytes().all(|c| c == b'0'))
}
fn digest(raw: &str) -> bool {
    // No unsupported string (e.g. a lone surrogate) can match ASCII hex.
    json_text::parse(raw, false).is_ok_and(|value| hex64(&value))
}

/// Split only a previously validated array. This locates element boundaries;
/// JSON syntax, Unicode escapes and number conversion belong to json_text.
fn elements(raw: &str) -> Option<Vec<&str>> {
    if !raw.starts_with('[') {
        return None;
    }
    let inner = &raw[1..raw.len() - 1];
    if inner.trim().is_empty() {
        return Some(Vec::new());
    }
    let mut result = Vec::new();
    let (mut start, mut depth, mut quoted, mut escaped) = (0, 0, false, false);
    for (index, byte) in inner.bytes().enumerate() {
        if quoted {
            if escaped {
                escaped = false;
            } else if byte == b'\\' {
                escaped = true;
            } else if byte == b'"' {
                quoted = false;
            }
        } else {
            match byte {
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
    Some(result)
}

pub(crate) fn artifact_receipt(artifact: &Path) -> Result<ArtifactReceipt> {
    let receipt = PathBuf::from(format!("{}.nm_to_enst.json", artifact.display()));
    let raw = read_json_text(&receipt)?;
    let recorded = Object::parse(&raw)?
        .ok_or_else(|| error("ValueError", "output receipt must be a JSON object"))?;
    let identity = Object::parse(recorded.get("mapping_identity").unwrap_or("null"))?;
    let valid = identity.as_ref().is_some_and(|identity| {
        identity.keys(&[
            "mapping_sha256",
            "mapping_bytes",
            "provenance_sha256",
            "anchor_inputs",
        ]) && digest(identity.get("mapping_sha256").unwrap())
            && size(identity.get("mapping_bytes").unwrap())
            && identity
                .get("provenance_sha256")
                .is_some_and(|value| value == "null" || digest(value))
            && identity.get("anchor_inputs").unwrap().starts_with('[')
    });
    require(valid, "output receipt has an invalid mapping identity")?;
    for raw_anchor in elements(identity.as_ref().unwrap().get("anchor_inputs").unwrap()).unwrap() {
        let anchor = Object::parse(raw_anchor)?;
        require(
            anchor.as_ref().is_some_and(|anchor| {
                anchor.keys(&["name", "bytes", "sha256"])
                    && anchor
                        .get("name")
                        .is_some_and(|name| name.starts_with('"') && name != "\"\"")
                    && size(anchor.get("bytes").unwrap())
                    && digest(anchor.get("sha256").unwrap())
            }),
            "output receipt has an invalid anchor identity",
        )?;
    }
    // Preserve short-circuit I/O order: the original calls stat twice before
    // hashing, and does not read the artifact when an earlier shape check fails.
    require(
        is_file(artifact)?
            && metadata(artifact)?.len() > 0
            && integer(recorded.get("artifact_bytes").unwrap_or("null"))
            && matches_value(
                recorded.get("artifact_bytes").unwrap(),
                &Json::Int(metadata(artifact)?.len() as i128),
            )?
            && matches_value(
                recorded.get("artifact_sha256").unwrap_or("null"),
                &string(sha256(artifact)?),
            )?,
        "output content changed or a prior writer left an incomplete artifact",
    )?;
    let identity = recorded.span("mapping_identity").unwrap();
    Ok(ArtifactReceipt {
        json: raw,
        identity,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn validated_array_boundaries_keep_quoted_delimiters_and_nested_values() {
        let text = r#"[ {"x":"a,]\\\"","y":[1,{"z":2}]}, "str", null ]"#;
        json_text::parse_metadata_comparison(text, &Json::Array(vec![])).unwrap();
        let items = elements(text).unwrap();
        assert_eq!(items.len(), 3);
        assert!(items[0].starts_with("{\"x\":"));
        assert_eq!(items[1], "\"str\"");
        assert_eq!(items[2], "null");
        assert!(elements("[ \n\t ]").unwrap().is_empty());
        assert!(elements("null").is_none());
    }
    #[test]
    fn exact_integer_types_admit_wide_sizes_without_projecting_values() {
        for value in [
            "0",
            "-0",
            "1",
            "1461501637330902918203684832716283019655932542976",
        ] {
            assert!(size(value));
        }
        for value in ["-1", "true", "1.0", "1e0", "null", "\"1\""] {
            assert!(!size(value));
        }
        assert!(
            !matches_value(
                "1461501637330902918203684832716283019655932542976",
                &Json::Int(1)
            )
            .unwrap()
        );
        assert!(!digest(r#""\ud800""#));
    }
}
