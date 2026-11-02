//! Python PVAR keys are text or tuples of decoded INFO scalars. Tuple slicing
//! removes an element, while text slicing removes one Unicode character.
//! A repeated NaN-bearing tuple is a new dictionary key for each source record;
//! clones of the same tuple within one record retain Python's identity shortcut.
use super::{Result, exception};
use crate::vcf_variant::{InfoValue, Scalar};
use std::hash::{Hash, Hasher};
use std::sync::Arc;

#[derive(Clone)]
pub(super) enum Key {
    Text(Arc<str>),
    Tuple(Arc<[Scalar]>),
}
impl std::fmt::Debug for Key {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Text(text) => f.debug_tuple("Text").field(text).finish(),
            Self::Tuple(items) => f
                .debug_struct("Tuple")
                .field("elements", &items.len())
                .finish(),
        }
    }
}
impl Key {
    pub fn from_info(value: InfoValue) -> Result<Self> {
        if value.scalar {
            match value.values.into_iter().next().unwrap_or(Scalar::Missing) {
                Scalar::Text(text) => Ok(Self::Text(text.into())),
                _ => Err(exception("TypeError", "PVAR scalar is not subscriptable")),
            }
        } else {
            Ok(Self::Tuple(value.values.into()))
        }
    }
    pub fn residue(&self) -> Self {
        match self {
            Self::Text(text) => {
                let end = text.char_indices().last().map(|(at, _)| at).unwrap_or(0);
                Self::Text(text[..end].into())
            }
            Self::Tuple(items) => Self::Tuple(items[..items.len().saturating_sub(1)].into()),
        }
    }
    pub fn bytes(&self) -> usize {
        match self {
            Self::Text(text) => text.len() + 24,
            Self::Tuple(items) => {
                items
                    .iter()
                    .map(|item| match item {
                        Scalar::Text(text) => text.len() + 32,
                        _ => 32,
                    })
                    .sum::<usize>()
                    + 24
            }
        }
    }
}
fn numeric(value: &Scalar) -> Option<f64> {
    match value {
        Scalar::Int(v) => Some(*v as f64),
        Scalar::Float(v) => Some(*v),
        Scalar::Bool(v) => Some(f64::from(*v)),
        _ => None,
    }
}
impl PartialEq for Key {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Text(a), Self::Text(b)) => a == b,
            (Self::Tuple(a), Self::Tuple(b)) => {
                Arc::ptr_eq(a, b)
                    || (a.len() == b.len()
                        && a.iter().zip(b.iter()).all(|(a, b)| match (a, b) {
                            (Scalar::Missing, Scalar::Missing) => true,
                            (Scalar::Text(a), Scalar::Text(b)) => a == b,
                            _ => matches!((numeric(a), numeric(b)), (Some(a), Some(b)) if a == b),
                        }))
            }
            _ => false,
        }
    }
}
impl Eq for Key {}
impl Hash for Key {
    fn hash<H: Hasher>(&self, state: &mut H) {
        match self {
            Self::Text(text) => {
                0_u8.hash(state);
                text.hash(state);
            }
            Self::Tuple(items) => {
                1_u8.hash(state);
                if items
                    .iter()
                    .any(|v| matches!(v, Scalar::Float(n) if n.is_nan()))
                {
                    // Distinct tuples containing NaN never compare equal. The
                    // pointer is runtime identity only and is never serialized.
                    (items.as_ptr() as usize).hash(state);
                } else {
                    items.len().hash(state);
                    for item in items.iter() {
                        match item {
                            Scalar::Missing => 0_u8.hash(state),
                            Scalar::Text(text) => {
                                1_u8.hash(state);
                                text.hash(state);
                            }
                            _ => {
                                2_u8.hash(state);
                                let n = numeric(item).unwrap();
                                (if n == 0.0 { 0 } else { n.to_bits() }).hash(state);
                            }
                        }
                    }
                }
            }
        }
    }
}
