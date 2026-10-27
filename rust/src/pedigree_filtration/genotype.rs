//! Separate literal Python predicates; no ploidy, sex or ALT normalization.
use super::{SourceError, SourceResult};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Allele {
    None,
    Bool(bool),
    /// Canonical arbitrary-precision decimal text avoids an integer-width cast.
    Integer(String),
    String(String),
}

impl Allele {
    fn equals_integer(&self, value: &str) -> bool {
        match self {
            Self::Bool(false) => value == "0",
            Self::Bool(true) => value == "1",
            Self::Integer(text) => text == value,
            _ => false,
        }
    }

    fn python_equal(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Integer(value), Self::Bool(_)) => other.equals_integer(value),
            (Self::Bool(_), Self::Integer(value)) => self.equals_integer(value),
            _ => self == other,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Genotype {
    None,
    Sequence(Vec<Allele>),
    String(String),
    Scalar(&'static str),
}

impl Genotype {
    fn not_iterable(&self) -> SourceError {
        let name = match self {
            Self::None => "NoneType",
            Self::Scalar(name) => name,
            _ => unreachable!("iterable genotype"),
        };
        SourceError::type_error(format!("'{name}' object is not iterable"))
    }

    pub fn homozygous_or_hemizygous(&self) -> SourceResult<bool> {
        match self {
            Self::None | Self::Scalar(_) => Err(self.not_iterable()),
            Self::String(value) => {
                let mut chars = value.chars();
                Ok(chars
                    .next()
                    .is_some_and(|first| chars.all(|next| next == first)))
            }
            Self::Sequence(values) => Ok(values.first().is_some_and(|first| {
                !matches!(first, Allele::None)
                    && !first.equals_integer("0")
                    && values.iter().all(|value| first.python_equal(value))
            })),
        }
    }

    fn contains_integer(&self, value: &str) -> SourceResult<bool> {
        match self {
            Self::None => Err(SourceError::type_error(
                "argument of type 'NoneType' is not iterable",
            )),
            Self::Scalar(name) => Err(SourceError::type_error(format!(
                "argument of type '{name}' is not iterable"
            ))),
            Self::String(_) => Err(SourceError::type_error(
                "'in <string>' requires string as left operand, not int",
            )),
            Self::Sequence(values) => Ok(values.iter().any(|item| item.equals_integer(value))),
        }
    }

    pub fn carrier(&self) -> SourceResult<bool> {
        Ok(self.contains_integer("1")? || self.contains_integer("2")?)
    }

    pub fn exists(&self) -> SourceResult<bool> {
        match self {
            Self::None => Ok(false),
            Self::Scalar(_) => Err(self.not_iterable()),
            Self::String(_) => Ok(true),
            Self::Sequence(values) => Ok(values.iter().all(|value| !matches!(value, Allele::None))),
        }
    }

    pub fn has_variant(&self) -> SourceResult<bool> {
        Ok(self.exists()? && self.contains_integer("1")?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn ints(values: &[&str]) -> Genotype {
        Genotype::Sequence(
            values
                .iter()
                .map(|value| Allele::Integer((*value).into()))
                .collect(),
        )
    }
    #[test]
    fn carrier_and_source_are_different() {
        assert_eq!(ints(&["2"]).carrier(), Ok(true));
        assert_eq!(ints(&["2"]).has_variant(), Ok(false));
        assert_eq!(ints(&["3", "3"]).homozygous_or_hemizygous(), Ok(true));
        assert_eq!(
            Genotype::Sequence(vec![Allele::Integer("1".into()), Allele::None]).carrier(),
            Ok(true)
        );
        assert_eq!(Genotype::Sequence(vec![]).exists(), Ok(true));
        assert_eq!(
            Genotype::Sequence(vec![]).homozygous_or_hemizygous(),
            Ok(false)
        );
    }
    #[test]
    fn python_bool_integer_and_string_equality() {
        assert_eq!(
            Genotype::Sequence(vec![Allele::Bool(true), Allele::Integer("1".into())])
                .homozygous_or_hemizygous(),
            Ok(true)
        );
        assert_eq!(
            Genotype::Sequence(vec![Allele::String("0".into())]).homozygous_or_hemizygous(),
            Ok(true)
        );
        assert_eq!(
            Genotype::Scalar("int").exists().unwrap_err().message,
            "'int' object is not iterable"
        );
    }
}
