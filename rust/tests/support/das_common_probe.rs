//! Ordered binary64 protocol for the independent DAS common-AF comparison.
#[path = "../../src/das_common.rs"]
mod das_common;
#[allow(dead_code)]
#[path = "../../src/das_reference.rs"]
mod das_reference;
#[allow(dead_code)]
#[path = "../../src/das_rescue.rs"]
mod das_rescue;
#[allow(dead_code)]
#[path = "../../src/exon_af_stats.rs"]
mod exon_af_stats;
use std::io::{self, BufRead};

type Error = Box<dyn std::error::Error>;
fn float(word: &str) -> Result<f64, Error> {
    Ok(f64::from_bits(u64::from_str_radix(word, 16)?))
}
fn array(words: &mut std::str::SplitWhitespace<'_>) -> Result<Option<Vec<f64>>, Error> {
    let count: i64 = words.next().ok_or("missing array length")?.parse()?;
    if count == -1 {
        return Ok(None);
    }
    if count < 0 {
        return Err("invalid array length".into());
    }
    Ok(Some(
        (0..count)
            .map(|_| float(words.next().ok_or("missing array element")?))
            .collect::<Result<_, _>>()?,
    ))
}
fn main() -> Result<(), Error> {
    if std::env::args().skip(1).collect::<Vec<_>>() != ["filter"] {
        return Err("usage: probe filter".into());
    }
    for line in io::stdin().lock().lines() {
        let line = line?;
        let mut words = line.split_whitespace();
        let result = match words.next().ok_or("missing procedure")? {
            "counts" => {
                let total = words.next().ok_or("missing total")?.parse()?;
                let common = words.next().ok_or("missing common count")?.parse()?;
                let threshold = float(words.next().ok_or("missing threshold")?)?;
                das_common::residue_counts(total, common, threshold)
            }
            "filter" => {
                let global = words.next().ok_or("missing global flag")? == "1";
                let common = float(words.next().ok_or("missing common AF")?)?;
                let rare = float(words.next().ok_or("missing rare AF")?)?;
                let fraction = float(words.next().ok_or("missing residue fraction")?)?;
                let count: i64 = words.next().ok_or("missing residues length")?.parse()?;
                let residues = if count == -1 {
                    None
                } else if count < 0 {
                    return Err("invalid residues length".into());
                } else {
                    Some(
                        (0..count)
                            .map(|_| array(&mut words)?.ok_or_else(|| "null residue list".into()))
                            .collect::<Result<Vec<_>, Error>>()?,
                    )
                };
                let flat = array(&mut words)?;
                das_common::filter(residues, flat, global, common, rare, fraction)
            }
            _ => return Err("unknown procedure".into()),
        };
        if words.next().is_some() {
            return Err("trailing fields".into());
        }
        println!(
            "{} {} {} {} {}",
            result.fully_common,
            result.total,
            match result.fraction {
                das_common::Fraction::IntegerZero => "int:0".into(),
                das_common::Fraction::Float(value) => format!("float:{:016x}", value.to_bits()),
            },
            u8::from(result.removed),
            u8::from(result.with_data)
        );
    }
    Ok(())
}
