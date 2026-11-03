//! DAS high-confidence benign missense counts. Every qualifying annotation
//! occurrence contributes, including repeated CSQs, transcripts and prefixes.
use crate::alpha_domain::vcf::{Context, visit};
use crate::builder_error;
use crate::cli::{Result, invalid};
use crate::vcf_variant::{Definition, InfoValue, Record, Scalar};
use num_bigint::BigInt;
use std::collections::{HashMap, HashSet};
use std::path::Path;

fn error(category: &'static str, message: &str) -> Box<dyn std::error::Error> {
    builder_error::category(category, invalid(message))
}

/// The original argparse type is Python int, including Unicode digits and
/// arbitrary magnitudes within Python 3.12's default decimal-string limit.
pub(crate) fn minimum(text: &str) -> Result<BigInt> {
    use crate::python_text::{decimal_digit_value, python_int_is_valid, python_strip};
    if !python_int_is_valid(text) {
        return Err(invalid("--min-benign requires a signed integer"));
    }
    let digits: String = python_strip(text)
        .chars()
        .filter(|c| *c != '_')
        .map(|c| decimal_digit_value(c).map_or(c, |n| char::from(b'0' + n as u8)))
        .collect();
    digits
        .parse()
        .map_err(|_| invalid("--min-benign requires a signed integer"))
}

fn get(
    record: &Record<'_>,
    key: &str,
    types: &HashMap<String, Definition>,
) -> Result<Option<InfoValue>> {
    if !types.contains_key(key) {
        return Ok(None);
    }
    record.get(key, types).map_err(|cause| {
        let message = cause.to_string();
        if let Some(rest) = message.strip_prefix("ValueError: ") {
            error("ValueError", rest)
        } else {
            cause
        }
    })
}

fn equals_text(value: InfoValue, expected: &[&str]) -> Result<bool> {
    if value.scalar {
        return Ok(
            matches!(value.values.first(), Some(Scalar::Text(text)) if expected.contains(&text.as_str())),
        );
    }
    let words = value
        .values
        .into_iter()
        .map(|value| match value {
            Scalar::Text(text) => Ok(text),
            _ => Err(error("TypeError", "sequence item is not str")),
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(expected.contains(&words.join(",").as_str()))
}

fn truthy(value: &InfoValue) -> bool {
    if !value.scalar {
        return !value.values.is_empty();
    }
    match value.values.first() {
        None | Some(Scalar::Missing) => false,
        Some(Scalar::Text(text)) => !text.is_empty(),
        Some(Scalar::Int(number)) => *number != 0,
        Some(Scalar::Float(number)) => *number != 0.0,
        Some(Scalar::Bool(value)) => *value,
    }
}

pub(crate) fn counts(path: &Path, domains: &HashSet<String>) -> Result<HashMap<String, u64>> {
    counts_for(path, domains, None).map(|(counts, _)| counts)
}

/// An explicit indexed chromosome procedure reuses the identical record
/// selection and accumulation body. The ordinary whole-file reader stays the
/// default; the shard caller owns completeness and ordered aggregation.
pub(crate) fn counts_for(
    path: &Path,
    domains: &HashSet<String>,
    chromosome: Option<&str>,
) -> Result<(HashMap<String, u64>, u64)> {
    let context = Context::with_fields(
        path,
        &["DOMAINS", "Gene", "Consequence", "Feature_type", "Feature"],
    )?;
    let mut counts = HashMap::<String, u64>::new();
    let mut records = 0_u64;
    visit(path, &context, chromosome, |record, types| {
        records = records
            .checked_add(1)
            .ok_or_else(|| error("OverflowError", "VCF record count exceeds u64"))?;
        let Some(significance) = get(record, "CLNSIG", types)? else {
            return Ok(());
        };
        if !equals_text(significance, &["Benign"])? {
            return Ok(());
        }
        let Some(review) = get(record, "CLNREVSTAT", types)? else {
            return Ok(());
        };
        if !equals_text(review, &["practice_guideline", "reviewed_by_expert_panel"])? {
            return Ok(());
        }
        let Some(csq) = get(record, "CSQ", types)? else {
            return Ok(());
        };
        if !truthy(&csq) {
            return Ok(());
        }
        let annotations = if csq.scalar {
            match csq.values.into_iter().next() {
                Some(Scalar::Text(text)) => {
                    text.split(',').map(|s| Scalar::Text(s.into())).collect()
                }
                _ => {
                    return Err(error(
                        "AttributeError",
                        "CSQ scalar has no attribute 'split'",
                    ));
                }
            }
        } else {
            csq.values
        };
        for annotation in annotations {
            let Scalar::Text(text) = annotation else {
                return Err(error(
                    "AttributeError",
                    "CSQ tuple item has no attribute 'split'",
                ));
            };
            let fields: Vec<_> = text.split('|').collect();
            if context.field(&fields, "Feature_type")? != "Transcript"
                || !context.field(&fields, "Feature")?.starts_with("ENST")
            {
                continue;
            }
            if !context
                .field(&fields, "Consequence")?
                .contains("missense_variant")
            {
                continue;
            }
            let annotation = context.field(&fields, "DOMAINS")?;
            let gene = context.field(&fields, "Gene")?;
            if annotation.is_empty() || gene.is_empty() {
                continue;
            }
            for entry in annotation.split('&') {
                let mut parts = entry.split(':');
                let mut key = format!("{gene}:{}", parts.next().unwrap());
                for level in parts {
                    key.push(':');
                    key.push_str(level);
                    if domains.contains(&key) {
                        let count = counts.entry(key.clone()).or_default();
                        *count = count.checked_add(1).ok_or_else(|| {
                            error("OverflowError", "benign annotation count exceeds u64")
                        })?;
                    }
                }
            }
        }
        Ok(())
    })
    .map_err(|cause| {
        // pysam.fetch reports a noninteger position as OSError even though
        // other indexed parsing failures (such as a negative POS) are
        // ValueError. Keep this new procedure's boundary local; do not alter
        // the previously qualified shared Alpha reader or sequential route.
        if chromosome.is_some()
            && cause.to_string() == "a VCF line has a position that is not an integer"
        {
            error("OSError", &cause.to_string())
        } else {
            cause
        }
    })?;
    Ok((counts, records))
}
