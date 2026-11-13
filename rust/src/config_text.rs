//! Exact text transforms for the two Python installer configuration writers.
//!
//! This is not a YAML parser. Like `installer_config_paths.render_config` and
//! `reference_resources.render_config_values`, it replaces matching top-level
//! lines, preserves recognized trailing comments, and appends missing keys.
//! JSON-form configuration follows each producer's own Unicode policy.
//! File locking, source checks and publication belong to the caller.
//!
//! This shared API is staged before the installer helpers that will call it.
#![allow(dead_code)]

use crate::json_text::{self, DumpsOptions, Json};
use crate::python_text::is_python_space;
use std::fmt;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Writer {
    /// `installer_config_paths.render_config`: sorted update keys; JSON-form
    /// output leaves Unicode unescaped. YAML values still use ASCII escapes.
    Paths,
    /// `reference_resources.render_config_values`: insertion-order updates;
    /// JSON-form output uses ASCII escapes too.
    Reference,
}

#[derive(Debug)]
pub(crate) struct RenderError {
    pub(crate) kind: &'static str,
    pub(crate) message: String,
}

impl fmt::Display for RenderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for RenderError {}

fn error(kind: &'static str, message: String) -> RenderError {
    RenderError { kind, message }
}

/// Return the `\s+#.*` suffix only when it consumes the whole remainder.
fn comment_suffix(text: &str) -> Option<&str> {
    if text.is_empty() {
        return Some("");
    }
    let rest = text.trim_start_matches(is_python_space);
    (rest.len() < text.len() && rest.starts_with('#') && !rest.contains('\n')).then_some(text)
}

/// The Python regex tries double quotes, single quotes, then an unquoted
/// `[^#]*?` scalar. A malformed quote may therefore still use the last branch;
/// do not reinterpret that accepted text as a new YAML validation error.
fn trailing_comment(old: &str) -> &str {
    let old = old.trim_matches(is_python_space);
    let bytes = old.as_bytes();
    if matches!(bytes.first(), Some(b'"' | b'\'')) {
        let quote = bytes[0];
        let mut at = 1;
        while at < bytes.len() {
            if quote == b'"' && bytes[at] == b'\\' {
                at += 2;
            } else if quote == b'\'' && bytes[at] == b'\'' && bytes.get(at + 1) == Some(&b'\'') {
                at += 2;
            } else if bytes[at] == quote {
                if let Some(comment) = comment_suffix(&old[at + 1..]) {
                    return comment;
                }
                break;
            } else {
                at += 1;
            }
        }
    }
    let Some(hash) = old.find('#') else {
        return "";
    };
    let before = &old[..hash];
    let start = before.trim_end_matches(is_python_space).len();
    if start < hash { &old[start..] } else { "" }
}

/// Literal `^KEY[ \t]*:[^\n]*$` with Python's MULTILINE matching. Walking
/// original line starts also preserves its behavior for a key containing LF.
fn line_matches(text: &str, key: &str) -> Vec<std::ops::Range<usize>> {
    let mut matches = Vec::new();
    let mut at = 0;
    loop {
        let remainder = &text[at..];
        let mut next = at;
        if let Some(after_key) = remainder.strip_prefix(key) {
            let after_spaces = after_key.trim_start_matches([' ', '\t']);
            if after_spaces.starts_with(':') {
                let colon = text.len() - after_spaces.len();
                let end = text[colon..]
                    .find('\n')
                    .map_or(text.len(), |offset| colon + offset);
                matches.push(at..end);
                next = end;
            }
        }
        let Some(offset) = text[next..].find('\n') else {
            break;
        };
        at = next + offset + 1;
    }
    matches
}

/// Reproduce either original writer without opening or mutating a file.
/// `updates` has mapping semantics: duplicate keys keep their first position
/// and last value, as Python's input dict and JSON object do.
pub(crate) fn render(
    original: &str,
    updates: &[(String, Json)],
    writer: Writer,
) -> Result<String, RenderError> {
    let mut update_map = Json::Object(Vec::with_capacity(updates.len()));
    for (key, value) in updates {
        update_map.set(key, value.clone());
    }
    let Json::Object(updates) = update_map else {
        unreachable!()
    };
    if original
        .trim_start_matches(is_python_space)
        .starts_with('{')
    {
        let mut config = json_text::parse(original, false)
            .map_err(|message| error("JSONDecodeError", message))?;
        for (key, value) in updates {
            config.set(&key, value);
        }
        let mut result = json_text::dumps_with_options(
            &config,
            DumpsOptions {
                indent: Some(2),
                ensure_ascii: writer == Writer::Reference,
                ..DumpsOptions::default()
            },
        )
        .map_err(|message| error("ValueError", message))?;
        result.push('\n');
        return Ok(result);
    }
    render_lines(original, updates, writer, |value| {
        json_text::dumps_with_options(value, DumpsOptions::default())
            .map_err(|message| error("ValueError", message))
    })
}

/// Additive scalar boundary for callers preserving Python code points that
/// Rust `String` cannot hold. Each value is already a complete validated JSON
/// scalar spelling; the established line/comment/duplicate logic is shared.
pub(crate) fn render_yaml_scalars(
    original: &str,
    values: &[(String, String)],
    writer: Writer,
) -> Result<String, RenderError> {
    let mut updates = Json::Object(vec![]);
    for (key, value) in values {
        updates.set(key, Json::Str(value.clone()));
    }
    let Json::Object(updates) = updates else {
        unreachable!()
    };
    render_lines(original, updates, writer, |v| {
        Ok(v.as_str().unwrap().to_owned())
    })
}
fn render_lines(
    original: &str,
    mut updates: Vec<(String, Json)>,
    writer: Writer,
    mut scalar: impl FnMut(&Json) -> Result<String, RenderError>,
) -> Result<String, RenderError> {
    if writer == Writer::Paths {
        updates.sort_by(|left, right| left.0.cmp(&right.0));
    }
    let mut result = original.to_owned();
    for (key, value) in updates {
        let matches = line_matches(&result, &key);
        if matches.len() > 1 {
            let label = if writer == Writer::Paths {
                "configuration"
            } else {
                "reference configuration"
            };
            return Err(error("ValueError", format!("Duplicate {label} key: {key}")));
        }
        let mut replacement = key.clone();
        replacement.push_str(": ");
        replacement.push_str(&scalar(&value)?);
        if let Some(found) = matches.first() {
            // Python splits on the FIRST colon, even when the literal key
            // itself contains one. Preserve this edge behavior too.
            let old = result[found.clone()].split_once(':').unwrap().1;
            replacement.push_str(trailing_comment(old));
            result.replace_range(found.clone(), &replacement);
        } else {
            result.truncate(result.trim_end_matches('\n').len());
            result.push('\n');
            result.push_str(&replacement);
            result.push('\n');
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn updates(items: &[(&str, &str)]) -> Vec<(String, Json)> {
        items
            .iter()
            .map(|(key, value)| (key.to_string(), Json::Str(value.to_string())))
            .collect()
    }

    #[test]
    fn quoted_hashes_and_escaped_quotes_do_not_become_comments() {
        let original = "a: \"a # b\\\"c\"  # real\nb: 'x''#y'\t# other\n";
        assert_eq!(
            render(original, &updates(&[("b", "B"), ("a", "A")]), Writer::Paths).unwrap(),
            "a: \"A\"  # real\nb: \"B\"\t# other\n"
        );
    }

    #[test]
    fn append_order_and_empty_input_match_each_writer() {
        let changes = updates(&[("z", "1"), ("a", "2")]);
        assert_eq!(
            render("", &changes, Writer::Paths).unwrap(),
            "\na: \"2\"\nz: \"1\"\n"
        );
        assert_eq!(
            render("", &changes, Writer::Reference).unwrap(),
            "\nz: \"1\"\na: \"2\"\n"
        );
    }

    #[test]
    fn json_preserves_key_order_and_selects_original_unicode_policy() {
        let changes = updates(&[("z", "é"), ("a", "λ")]);
        assert_eq!(
            render("{\"z\": 1, \"b\": 2}", &changes, Writer::Paths).unwrap(),
            "{\n  \"z\": \"é\",\n  \"b\": 2,\n  \"a\": \"λ\"\n}\n"
        );
        assert_eq!(
            render("{}", &changes, Writer::Reference).unwrap(),
            "{\n  \"z\": \"\\u00e9\",\n  \"a\": \"\\u03bb\"\n}\n"
        );
    }

    #[test]
    fn duplicate_messages_and_literal_key_matching_are_preserved() {
        let err = render("a: 1\na : 2\n", &updates(&[("a", "x")]), Writer::Paths).unwrap_err();
        assert_eq!(err.kind, "ValueError");
        assert_eq!(err.message, "Duplicate configuration key: a");
        let err = render("a: 1\na\t: 2\n", &updates(&[("a", "x")]), Writer::Reference).unwrap_err();
        assert_eq!(err.message, "Duplicate reference configuration key: a");
        assert_eq!(
            render("ab: 0\na.: 1\n", &updates(&[("a.", "x")]), Writer::Paths).unwrap(),
            "ab: 0\na.: \"x\"\n"
        );
    }

    #[test]
    fn comments_use_python_whitespace_and_preserve_unrelated_line_endings() {
        assert_eq!(
            render(
                "a:\tx\u{1c}# note\r\nb: keep\r\n",
                &updates(&[("a", "x")]),
                Writer::Paths
            )
            .unwrap(),
            "a: \"x\"\u{1c}# note\nb: keep\r\n"
        );
        assert_eq!(trailing_comment("\"unfinished # comment"), " # comment");
        assert_eq!(trailing_comment("value#not_a_comment"), "");
    }

    #[test]
    fn no_updates_preserve_yaml_bytes_but_reformat_json() {
        assert_eq!(
            render("a: 1\r\n\n", &[], Writer::Paths).unwrap(),
            "a: 1\r\n\n"
        );
        assert_eq!(render("{ }", &[], Writer::Paths).unwrap(), "{}\n");
    }
}
