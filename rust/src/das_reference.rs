//! DAS reference and whole-protein Family decisions over typed metadata.
//! The two source predicates deliberately differ at equality (>= versus >).
use crate::das_rescue::Entry;

#[derive(Clone, Copy)]
pub(crate) struct Metadata {
    pub(crate) start: i64,
    pub(crate) end: i64,
    pub(crate) protein_length: Option<i64>,
}

#[derive(Debug, Clone, Copy)]
pub(crate) enum Failure {
    ZeroDivision,
    NullLength,
}
impl Failure {
    pub(crate) fn category(self) -> &'static str {
        match self {
            Self::ZeroDivision => "ZeroDivisionError",
            Self::NullLength => "TypeError",
        }
    }
}

/// Correctly rounded Python integer true division for the metadata's signed
/// 64-bit endpoints/length. Subtraction is exact in i128; only the final ratio
/// rounds to binary64. Casting each large integer to f64 first can double-round.
fn coverage(meta: Metadata) -> Result<f64, Failure> {
    let denominator = meta.protein_length.ok_or(Failure::NullLength)?;
    let numerator = i128::from(meta.end) - i128::from(meta.start) + 1;
    integer_ratio(numerator, i128::from(denominator))
}

/// Positive counts use the same single-rounding integer division as metadata.
pub(crate) fn count_ratio(numerator: u64, denominator: u64) -> Result<f64, Failure> {
    integer_ratio(i128::from(numerator), i128::from(denominator))
}

// Inputs are bounded by signed-64-bit endpoints or unsigned-64-bit counts.
// Their magnitudes are at most 2^64; shifted division fits in u128.
fn integer_ratio(numerator: i128, denominator: i128) -> Result<f64, Failure> {
    if denominator == 0 {
        return Err(Failure::ZeroDivision);
    }
    let sign = ((numerator < 0) ^ (denominator < 0)) as u64;
    let numerator = numerator.unsigned_abs();
    let denominator = denominator.unsigned_abs();
    if numerator == 0 {
        return Ok(f64::from_bits(sign << 63));
    }
    let mut exponent = denominator.leading_zeros() as i32 - numerator.leading_zeros() as i32;
    let below = if exponent >= 0 {
        numerator < denominator << exponent
    } else {
        numerator << -exponent < denominator
    };
    if below {
        exponent -= 1;
    }
    let shift = 52 - exponent;
    let (a, b) = if shift >= 0 {
        (numerator << shift, denominator)
    } else {
        (numerator, denominator << -shift)
    };
    let mut mantissa = a / b;
    let remainder = a % b;
    if remainder * 2 > b || (remainder * 2 == b && mantissa & 1 != 0) {
        mantissa += 1;
    }
    if mantissa == 1 << 53 {
        mantissa >>= 1;
        exponent += 1;
    }
    Ok(f64::from_bits(
        (sign << 63) | (((exponent + 1023) as u64) << 52) | (mantissa as u64 & ((1 << 52) - 1)),
    ))
}

pub(crate) fn round_four(value: f64) -> f64 {
    if value.is_finite() {
        format!("{value:.4}").parse().unwrap()
    } else {
        value
    }
}

#[derive(Debug, PartialEq)]
pub(crate) enum Reference {
    Unknown,
    NonFunctional,
    Functional,
}

pub(crate) fn reference(
    entry: Option<&Entry>,
    metadata: Option<Metadata>,
    threshold: f64,
    min_length: i64,
) -> Result<Reference, Failure> {
    let Some(entry) = entry else {
        return Ok(Reference::Unknown);
    };
    if !entry.has_go_molecular_function {
        return Ok(Reference::NonFunctional);
    }
    if ["Family", "Homologous_superfamily", "Repeat"].contains(&entry.kind.as_str()) {
        if let Some(meta) = metadata {
            if meta
                .protein_length
                .is_some_and(|length| length >= min_length)
                && coverage(meta)? >= threshold
            {
                return Ok(Reference::NonFunctional);
            }
        }
    }
    Ok(Reference::Functional)
}

#[derive(Debug)]
pub(crate) struct Family {
    pub(crate) kind: Option<String>,
    pub(crate) coverage: Option<f64>,
    pub(crate) protein_length: Option<i64>,
    pub(crate) removed: bool,
}

pub(crate) fn family(
    entry: Option<&Entry>,
    metadata: Option<Metadata>,
    threshold: f64,
    min_length: i64,
) -> Result<Family, Failure> {
    let empty = Family {
        kind: None,
        coverage: None,
        protein_length: None,
        removed: false,
    };
    let Some(meta) = metadata else {
        return Ok(empty);
    };
    if meta.protein_length == Some(0) {
        return Ok(empty);
    }
    let ratio = coverage(meta)?;
    let kind = entry.map(|entry| entry.kind.clone());
    let removed = kind.as_deref() == Some("Family")
        && ratio > threshold
        && meta.protein_length.unwrap() > min_length;
    Ok(Family {
        kind,
        coverage: Some(round_four(ratio)),
        protein_length: meta.protein_length,
        removed,
    })
}
