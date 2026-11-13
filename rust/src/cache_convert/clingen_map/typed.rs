//! Schema 10 version 2 preserves scalar code types emitted by pandas.
//! Version 1's text/NaN layout and NaN normalization remain unchanged.
use super::{Provenance, plain_code};
use crate::acmg_step::pickle_dict_stream::Item;
use crate::cache_convert::common::ExactValue;
use crate::cli::Result;
use crate::rust_cache::Schema;
use rkyv::{Archive, Deserialize, Serialize};

pub(super) const SCHEMA: Schema = Schema {
    name: "clingen_map",
    id: 10,
    version: 2,
};

#[derive(Archive, Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub(crate) enum Code {
    Text(String),
    Missing,
    /// pandas emits signed int64 or uint64 as a Python int. Both fit exactly.
    Integer(i128),
    /// Non-NaN IEEE-754 bits preserve negative zero and infinities.
    Float(u64),
    Boolean(bool),
}

impl Code {
    pub(crate) fn from_item(path: &str, value: &Item) -> Result<Self> {
        Ok(match value {
            Item::Int(value) => Self::Integer(i128::from(*value)),
            Item::Unsigned(value) => Self::Integer(i128::from(*value)),
            Item::Float(value) if !value.is_nan() => Self::Float(value.to_bits()),
            Item::Bool(value) => Self::Boolean(*value),
            _ => plain_code(path, value)?.into(),
        })
    }

    pub(crate) fn is_missing(&self) -> bool {
        matches!(self, Self::Missing)
    }

    pub(super) fn legacy(&self) -> Option<Option<String>> {
        match self {
            Self::Text(value) => Some(Some(value.clone())),
            Self::Missing => Some(None),
            _ => None,
        }
    }

    pub(super) fn exact(self) -> ExactValue {
        match self {
            Self::Text(value) => ExactValue::Text(value),
            Self::Missing => ExactValue::None,
            Self::Integer(value) => ExactValue::WideInt(value),
            Self::Float(bits) => ExactValue::Float(bits),
            Self::Boolean(value) => ExactValue::Bool(value),
        }
    }
}

impl From<Option<String>> for Code {
    fn from(value: Option<String>) -> Self {
        value.map_or(Self::Missing, Self::Text)
    }
}

#[derive(Archive, Serialize, Deserialize)]
pub(super) struct ClingenMap {
    pub(super) provenance: Vec<Provenance>,
    pub(super) variants: Vec<String>,
    pub(super) codes: Vec<Code>,
}

impl ArchivedCode {
    pub(super) fn text(&self) -> Option<&str> {
        match self {
            Self::Text(value) => Some(value.as_str()),
            _ => None,
        }
    }

    pub(super) fn plain(&self) -> Code {
        match self {
            Self::Text(value) => Code::Text(value.as_str().to_owned()),
            Self::Missing => Code::Missing,
            Self::Integer(value) => Code::Integer(value.to_native()),
            Self::Float(bits) => Code::Float(bits.to_native()),
            Self::Boolean(value) => Code::Boolean(*value),
        }
    }

    pub(super) fn valid(&self) -> bool {
        match self {
            Self::Integer(value) => {
                (i128::from(i64::MIN)..=i128::from(u64::MAX)).contains(&value.to_native())
            }
            Self::Float(bits) => !f64::from_bits(bits.to_native()).is_nan(),
            _ => true,
        }
    }
}
