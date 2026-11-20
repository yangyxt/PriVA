use super::*;
use crate::installed_environment::yaml::{self, Value as YamlValue, Yaml};
use std::rc::Rc;

#[derive(Clone)]
pub(crate) struct Config {
    yaml: Rc<Yaml>,
    changes: Paths,
    removed: Vec<String>,
}
#[derive(Clone, Copy)]
pub(crate) enum Value<'a> {
    Missing,
    Text(&'a str),
    Node(&'a Yaml, usize),
}
impl Config {
    pub(crate) fn parse(text: &str) -> Result<Self> {
        Ok(Self {
            yaml: Rc::new(yaml::load_unique(text)?),
            changes: vec![],
            removed: vec![],
        })
    }
    pub(super) fn mapping(&self) -> bool {
        self.yaml.mapping(self.yaml.root).is_some()
    }
    pub(super) fn get(&self, key: &str) -> Value<'_> {
        if let Some((_, v)) = self.changes.iter().find(|(k, _)| k == key) {
            return Value::Text(v);
        }
        if self.removed.iter().any(|k| k == key) {
            return Value::Missing;
        }
        self.yaml
            .get(self.yaml.root, key)
            .map_or(Value::Missing, |n| Value::Node(&self.yaml, n))
    }
    pub(super) fn index(&self, key: &str) -> Result<Value<'_>> {
        let v = self.get(key);
        if matches!(v, Value::Missing) {
            Err(error("KeyError", repr(key)))
        } else {
            Ok(v)
        }
    }
    pub(super) fn set(&mut self, key: &str, value: String) {
        set(&mut self.changes, key, value);
    }
    pub(super) fn remove(&mut self, key: &str) {
        self.changes.retain(|(k, _)| k != key);
        self.removed.push(key.into());
    }
    pub(super) fn overlay(&self, values: &Paths) -> Self {
        let mut other = self.clone();
        for (k, v) in values {
            other.set(k, v.clone());
        }
        other
    }
    pub(super) fn keys(&self) -> Vec<String> {
        let mut keys = self
            .yaml
            .mapping(self.yaml.root)
            .into_iter()
            .flatten()
            .filter_map(|(k, _)| self.yaml.string(*k))
            .filter(|k| !self.removed.iter().any(|r| r == k))
            .map(str::to_owned)
            .collect::<Vec<_>>();
        for (k, _) in &self.changes {
            if !keys.contains(k) {
                keys.push(k.clone());
            }
        }
        keys
    }
}
impl<'a> Value<'a> {
    pub(super) fn string(self) -> Option<&'a str> {
        match self {
            Self::Text(s) => Some(s),
            Self::Node(y, n) => y.string(n),
            _ => None,
        }
    }
    pub(super) fn null(self) -> bool {
        matches!(self, Self::Missing)
            || matches!(self,Self::Node(y,n) if matches!(y.value(n),YamlValue::Null))
    }
    pub(super) fn empty(self) -> bool {
        self.null() || matches!(self.string(), Some("" | "null"))
    }
    pub(super) fn kind(self) -> &'static str {
        match self {
            Self::Missing => "NoneType",
            Self::Text(_) => "str",
            Self::Node(y, n) => match y.value(n) {
                YamlValue::Null => "NoneType",
                YamlValue::Bool(_) => "bool",
                YamlValue::Str(_) | YamlValue::Codepoints { .. } => "str",
                YamlValue::Number(_) => {
                    if y.tag(n) == "int" {
                        "int"
                    } else {
                        "float"
                    }
                }
                YamlValue::Sequence(_) => "list",
                YamlValue::Mapping(_) => "dict",
                YamlValue::Other(v) => {
                    if v.starts_with("binary:") {
                        "bytes"
                    } else if v == "set" {
                        "set"
                    } else if v == "tuple" {
                        "tuple"
                    } else if v.contains(['T', ' ']) {
                        "datetime"
                    } else {
                        "date"
                    }
                }
            },
        }
    }
    pub(super) fn truth(self) -> bool {
        match self {
            Self::Missing => false,
            Self::Text(s) => !s.is_empty(),
            Self::Node(y, n) => match y.value(n) {
                YamlValue::Null => false,
                YamlValue::Bool(v) => *v,
                YamlValue::Str(s) => !s.is_empty(),
                YamlValue::Codepoints { points, .. } => !points.is_empty(),
                YamlValue::Number(s) => s.parse::<f64>().map_or(true, |v| v != 0.0),
                YamlValue::Sequence(v) => !v.is_empty(),
                YamlValue::Mapping(v) => !v.is_empty(),
                YamlValue::Other(_) => true,
            },
        }
    }
    pub(super) fn text(self) -> String {
        match self {
            Self::Missing => "None".into(),
            Self::Text(s) => s.into(),
            Self::Node(y, n) => y.scalar_display(n),
        }
    }
    pub(super) fn path(self) -> Result<PathBuf> {
        self.string().map(path_value).ok_or_else(||error("TypeError",format!("argument should be a str or an os.PathLike object where __fspath__ returns a str, not '{}'",self.kind())))
    }
    pub(super) fn mapping(self) -> Option<Vec<(Option<&'a str>, Value<'a>)>> {
        match self {
            Self::Node(y, n) => y.mapping(n).map(|v| {
                v.iter()
                    .map(|(k, v)| (y.string(*k), Value::Node(y, *v)))
                    .collect()
            }),
            _ => None,
        }
    }
    pub(super) fn get(self, key: &str) -> Self {
        match self {
            Self::Node(y, n) => y.get(n, key).map_or(Self::Missing, |v| Self::Node(y, v)),
            _ => Self::Missing,
        }
    }
    pub(super) fn unhashable(self) -> bool {
        matches!(self.kind(), "list" | "dict" | "set")
    }
}
