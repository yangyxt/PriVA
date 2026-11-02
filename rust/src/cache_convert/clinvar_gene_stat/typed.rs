//! Schema 5 v2 changes only the stored gene-symbol representation. The complete
//! variant sets, lookup keys and free-key dictionary orders keep their v1 meaning.
use super::{Domain, Transcript, VariantSets};
use crate::cache_convert::Provenance;
use crate::cache_convert::common::ExactValue;
use crate::cli::{Result, invalid};
use crate::pickle_reader::{ClinvarScalar, Value};
use crate::rust_cache::Schema;
use rkyv::{Archive, Deserialize, Serialize};

pub(crate) const SCHEMA: Schema = Schema {
    name: "clinvar_gene_stat",
    id: 5,
    version: 2,
};

/// Scalar kinds observed from `SYMBOL.iloc[0]` after the original collector's
/// pandas casts. A NumPy scalar and a Python float are deliberately distinct.
#[derive(Archive, Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub(crate) enum Symbol {
    Text(String),
    Float(u64),
    PandasNA,
    NumpyI32(i32),
    NumpyF32(u32),
}

impl Symbol {
    pub(super) fn from_pickle(value: &Value) -> Result<Self> {
        Ok(match value {
            Value::Str(text) => Self::Text(text.to_string()),
            Value::Float(value) => Self::Float(value.to_bits()),
            _ => match value.clinvar_scalar() {
                Some(ClinvarScalar::PandasNA) => Self::PandasNA,
                Some(ClinvarScalar::NumpyI32(value)) => Self::NumpyI32(value),
                Some(ClinvarScalar::NumpyF32(bits)) => Self::NumpyF32(bits),
                None => {
                    return Err(invalid(format!(
                        "ClinVar gene symbol has unsupported type {}",
                        value.type_name()
                    )));
                }
            },
        })
    }
    pub(super) fn exact(self) -> ExactValue {
        match self {
            Self::Text(value) => ExactValue::Text(value),
            Self::Float(bits) => ExactValue::Float(bits),
            Self::PandasNA => ExactValue::PandasNA,
            Self::NumpyI32(value) => ExactValue::NumpyI32(value),
            Self::NumpyF32(bits) => ExactValue::NumpyF32(bits),
        }
    }
}

impl ArchivedSymbol {
    pub(super) fn owned(&self) -> Symbol {
        match self {
            Self::Text(value) => Symbol::Text(value.to_string()),
            Self::Float(bits) => Symbol::Float(bits.to_native()),
            Self::PandasNA => Symbol::PandasNA,
            Self::NumpyI32(value) => Symbol::NumpyI32(value.to_native()),
            Self::NumpyF32(bits) => Symbol::NumpyF32(bits.to_native()),
        }
    }
}

#[derive(Archive, Serialize, Deserialize, Debug)]
pub(crate) struct Gene {
    pub(crate) ensembl_id: String,
    pub(crate) gene_symbol: Symbol,
    pub(crate) sets: VariantSets,
    pub(crate) domains: Vec<Domain>,
    pub(crate) transcripts: Vec<Transcript>,
}

#[derive(Archive, Serialize, Deserialize, Debug)]
pub(crate) struct Root {
    pub(crate) provenance: Vec<Provenance>,
    pub(crate) variants: Vec<String>,
    pub(crate) genes: Vec<String>,
    pub(crate) entries: Vec<Gene>,
}
