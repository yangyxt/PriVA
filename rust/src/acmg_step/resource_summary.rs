//! `priva acmg-step inspect pickle|json PATH`: what a resource file holds, as
//! counts of Python types by position, read with the crate's own pickle and
//! JSON readers. The same summary from Python (`probe_group3_resources.py` in
//! the group 3 validation folder) must be identical, which shows the readers
//! see the installed caches as Python does. No value is printed.
//!
//! One line per (path, type): `$` is the top object, `.key` and `.value` the
//! keys and values of a dict, `[]` the items of a list or set, `[i]` position
//! `i` of a tuple; the type is Python's `type(value).__name__`.
use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use crate::pickle_reader::{self, Value};
use std::collections::BTreeMap;
use std::io::Read;
use std::path::Path;

type Counts = BTreeMap<(String, &'static str), u64>;

fn pickle_counts(top: &Value) -> Counts {
    let mut counts = Counts::new();
    let mut stack: Vec<(String, Value)> = vec![("$".to_owned(), top.clone())];
    while let Some((path, value)) = stack.pop() {
        *counts.entry((path.clone(), value.type_name())).or_insert(0) += 1;
        match &value {
            Value::Dict(dict) => {
                for (key, item) in &dict.borrow().entries {
                    stack.push((format!("{path}.key"), key.clone()));
                    stack.push((format!("{path}.value"), item.clone()));
                }
            }
            Value::List(items) => {
                for item in items.borrow().iter() {
                    stack.push((format!("{path}[]"), item.clone()));
                }
            }
            Value::Set(set) => {
                for item in &set.borrow().items {
                    stack.push((format!("{path}[]"), item.clone()));
                }
            }
            Value::FrozenSet(set) => {
                for item in &set.items {
                    stack.push((format!("{path}[]"), item.clone()));
                }
            }
            Value::Tuple(items) => {
                for (position, item) in items.iter().enumerate() {
                    stack.push((format!("{path}[{position}]"), item.clone()));
                }
            }
            _ => {}
        }
    }
    counts
}

fn json_type(value: &Json) -> &'static str {
    match value {
        Json::Null => "NoneType",
        Json::Bool(_) => "bool",
        Json::Int(_) => "int",
        Json::Float(_) => "float",
        Json::Str(_) => "str",
        Json::Array(_) => "list",
        Json::Object(_) => "dict",
    }
}

fn json_counts(top: &Json) -> Counts {
    let mut counts = Counts::new();
    let mut stack: Vec<(String, &Json)> = vec![("$".to_owned(), top)];
    while let Some((path, value)) = stack.pop() {
        *counts.entry((path.clone(), json_type(value))).or_insert(0) += 1;
        match value {
            Json::Object(items) => {
                for (_, item) in items {
                    *counts.entry((format!("{path}.key"), "str")).or_insert(0) += 1;
                    stack.push((format!("{path}.value"), item));
                }
            }
            Json::Array(items) => {
                for item in items {
                    stack.push((format!("{path}[]"), item));
                }
            }
            _ => {}
        }
    }
    counts
}

fn lines(counts: &Counts) -> String {
    let mut out = String::new();
    for ((path, name), count) in counts {
        out.push_str(&format!("{path}\t{name}\t{count}\n"));
    }
    out
}

/// The summary of a pickle file (gzip or plain, told apart by its first bytes).
pub(crate) fn pickle_summary(path: &Path) -> Result<String> {
    let bytes = pickle_reader::read_file_by_magic(path)?;
    let top = pickle_reader::load(&bytes)?;
    Ok(lines(&pickle_counts(&top)))
}

/// The summary of a JSON file (gzip or plain, told apart by its first bytes).
pub(crate) fn json_summary(path: &Path) -> Result<String> {
    let mut raw = Vec::new();
    std::fs::File::open(path)?.read_to_end(&mut raw)?;
    if raw.starts_with(&[0x1f, 0x8b]) {
        let mut plain = Vec::new();
        flate2::read::MultiGzDecoder::new(&raw[..]).read_to_end(&mut plain)?;
        raw = plain;
    }
    let text = String::from_utf8(raw).map_err(|_| invalid("the JSON file is not UTF-8"))?;
    let top = json_text::parse(&text, false).map_err(invalid)?;
    Ok(lines(&json_counts(&top)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn json_summary_counts_types_by_position() {
        let top = json_text::parse(
            r#"{"a": [1, 2.5, null], "b": {"c": true}, "a": ["x"]}"#,
            false,
        )
        .unwrap();
        // The repeated key keeps one entry with the last value, as json.loads does.
        assert_eq!(
            lines(&json_counts(&top)),
            "$\tdict\t1\n$.key\tstr\t2\n$.value\tdict\t1\n$.value\tlist\t1\n\
             $.value.key\tstr\t1\n$.value.value\tbool\t1\n$.value[]\tstr\t1\n"
        );
    }
}
