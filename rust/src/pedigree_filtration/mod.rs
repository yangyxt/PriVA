//! Complete optional pedigree filter: ordered native decisions, mature host I/O.
//!
//! Source contract: pedigree_filtration_per_fam.py SHA256 06d243890da53d172ac1bc4e5d18f4c6897d3eba2b6c0a2d0404564705dd0c5b.
//! 2026-10-27: preserve ALT1-or-ALT2 carrier versus complete ALT1 source,
//! family-wide first parents, lazy sample errors and original output prefixes.
pub mod family;
pub mod genotype;
pub mod stream;
pub mod variant_source;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceError {
    pub class: &'static str,
    pub message: String,
}

impl SourceError {
    pub fn type_error(message: impl Into<String>) -> Self {
        Self {
            class: "TypeError",
            message: message.into(),
        }
    }
}

pub type SourceResult<T> = std::result::Result<T, SourceError>;
