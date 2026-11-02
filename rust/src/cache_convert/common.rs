//! Conversion-only helpers. Runtime archives retain their individual typed schemas.
//!
//! Tagged JSON lines are an inspection/validation transport, never a runtime cache.
//! Floats are hexadecimal IEEE-754 bits, containers keep their type, and dictionaries
//! keep pair order. This lets a separate Python unpickler verify every archived value.
use super::Provenance;
use crate::cli::{Result, invalid};
use crate::json_text::quote;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

/// Exact plain values used only at conversion and comparison boundaries.
#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum ExactValue {
    None,
    Bool(bool),
    Int(i64),
    Float(u64),
    Text(String),
    Tuple(Vec<ExactValue>),
    List(Vec<ExactValue>),
    Dict(Vec<(ExactValue, ExactValue)>),
    Set(Vec<ExactValue>),
    FrozenSet(Vec<ExactValue>),
}

impl ExactValue {
    /// A deterministic tagged JSON value. Sets sort their encodings; free-key
    /// dictionaries and sequences retain order. No float formatting is involved.
    pub(crate) fn encoded(&self) -> String {
        fn sequence(kind: &str, items: &[ExactValue], sorted: bool) -> String {
            let mut values: Vec<String> = items.iter().map(ExactValue::encoded).collect();
            if sorted {
                values.sort();
            }
            format!("[\"{kind}\",[{}]]", values.join(","))
        }
        match self {
            Self::None => "[\"none\"]".to_owned(),
            Self::Bool(value) => format!("[\"bool\",{value}]"),
            Self::Int(value) => format!("[\"int\",\"{value}\"]"),
            Self::Float(bits) => format!("[\"float\",\"{bits:016x}\"]"),
            Self::Text(value) => format!("[\"str\",{}]", quote(value, false)),
            Self::Tuple(items) => sequence("tuple", items, false),
            Self::List(items) => sequence("list", items, false),
            Self::Set(items) => sequence("set", items, true),
            Self::FrozenSet(items) => sequence("frozenset", items, true),
            Self::Dict(items) => {
                let pairs: Vec<String> = items
                    .iter()
                    .map(|(key, value)| format!("[{},{}]", key.encoded(), value.encoded()))
                    .collect();
                format!("[\"dict\",[{}]]", pairs.join(","))
            }
        }
    }
}

/// One top-level key/value per line; schemas with a root set use a None key.
pub(crate) fn write_entry(
    output: &mut dyn Write,
    key: &ExactValue,
    value: &ExactValue,
) -> Result<()> {
    writeln!(output, "[{},{}]", key.encoded(), value.encoded())?;
    Ok(())
}

/// Portable provenance: generated bundled files must not contain host filesystem paths.
/// Complete acquisition provenance and checksums live beside the managed resources.
pub(crate) fn provenance(pickle: &str) -> Result<Vec<Provenance>> {
    let path = Path::new(pickle);
    let source = path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| invalid("the pickle needs a UTF-8 file name"))?;
    Ok(vec![
        Provenance {
            key: "source".to_owned(),
            value: source.to_owned(),
        },
        Provenance {
            key: "source_bytes".to_owned(),
            value: std::fs::metadata(path)?.len().to_string(),
        },
        Provenance {
            key: "converter".to_owned(),
            value: format!("priva cache-convert {}", env!("CARGO_PKG_VERSION")),
        },
    ])
}

/// Stable insertion-ordered text interning; hashing is lookup-only, never output order.
#[derive(Default)]
pub(crate) struct StringPool {
    texts: Vec<String>,
    ids: HashMap<String, u32>,
}

impl StringPool {
    pub(crate) fn intern(&mut self, text: &str) -> Result<u32> {
        if let Some(id) = self.ids.get(text) {
            return Ok(*id);
        }
        let id = u32::try_from(self.texts.len())
            .map_err(|_| invalid("the cache text pool exceeds u32 indices"))?;
        self.texts.push(text.to_owned());
        self.ids.insert(text.to_owned(), id);
        Ok(id)
    }

    pub(crate) fn into_texts(self) -> Vec<String> {
        self.texts
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn comparison_transport_keeps_types_float_bits_and_pair_order() {
        let f = ExactValue::Float(0x7ff8_0000_0000_0042);
        assert_eq!(f.encoded(), "[\"float\",\"7ff8000000000042\"]");
        assert_ne!(
            ExactValue::Float(0).encoded(),
            ExactValue::Float(1 << 63).encoded()
        );
        assert_ne!(
            ExactValue::Int(1).encoded(),
            ExactValue::Bool(true).encoded()
        );
        assert_ne!(
            ExactValue::Tuple(vec![]).encoded(),
            ExactValue::List(vec![]).encoded()
        );
        let a = ExactValue::Text("α\n\"".to_owned());
        let b = ExactValue::Text("b".to_owned());
        let forward = ExactValue::Dict(vec![(a.clone(), b.clone()), (b.clone(), a.clone())]);
        let reverse = ExactValue::Dict(vec![(b.clone(), a.clone()), (a.clone(), b.clone())]);
        assert_ne!(forward.encoded(), reverse.encoded());
        assert_eq!(
            ExactValue::Set(vec![a.clone(), b.clone()]).encoded(),
            ExactValue::Set(vec![b, a]).encoded()
        );
        let mut out = Vec::new();
        write_entry(&mut out, &ExactValue::None, &forward).unwrap();
        assert_eq!(out.last(), Some(&b'\n'));
        crate::json_text::parse(std::str::from_utf8(&out).unwrap(), true).unwrap();
    }

    #[test]
    fn text_pool_reuses_ids_without_reordering() {
        let mut pool = StringPool::default();
        assert_eq!(pool.intern("second").unwrap(), 0);
        assert_eq!(pool.intern("first").unwrap(), 1);
        assert_eq!(pool.intern("second").unwrap(), 0);
        assert_eq!(pool.into_texts(), ["second", "first"]);
    }
}
