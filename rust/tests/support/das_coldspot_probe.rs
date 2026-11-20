//! Dependency-free probe of the coldspot loop and exact integer division.
#[path = "../../src/das_coldspot.rs"]
mod das_coldspot;
#[path = "../../src/das_reference.rs"]
#[allow(dead_code)]
mod das_reference;
#[path = "../../src/das_rescue.rs"]
#[allow(dead_code)]
mod das_rescue;
#[path = "../../src/exon_af_stats.rs"]
#[allow(dead_code)]
mod exon_af_stats;
use das_coldspot::{Interval, filter};
use std::io::{self, BufRead};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    if std::env::args().nth(1).as_deref() != Some("batch") {
        return Err("usage: das-coldspot-probe batch < REQUESTS".into());
    }
    for line in io::stdin().lock().lines() {
        let line = line?;
        let mut words = line.split_whitespace();
        let mut next = || words.next().ok_or("missing field");
        let kind = match next()? {
            "0" => None,
            "1" => Some("Family"),
            "2" => Some("Homologous_superfamily"),
            "3" => Some("Domain"),
            "4" => Some("family"),
            _ => return Err("invalid kind".into()),
        };
        let metadata = next()? == "1";
        let start = next()?.parse()?;
        let end = next()?.parse()?;
        let present = next()? == "1";
        let count: usize = next()?.parse()?;
        let threshold = f64::from_bits(u64::from_str_radix(next()?, 16)?);
        let mut spots = Vec::new();
        for _ in 0..count {
            spots.push(Interval {
                start: next()?.parse()?,
                end: next()?.parse()?,
            });
        }
        if words.next().is_some() {
            return Err("extra fields".into());
        }
        let summary = filter(
            kind,
            metadata.then_some(Interval { start, end }),
            present.then_some(spots.as_slice()),
            threshold,
        );
        println!(
            "{:016x} {:016x} {} {} {}",
            summary.fraction_of_coldspot.to_bits(),
            summary.fraction_of_domain.to_bits(),
            u8::from(summary.removed),
            u8::from(summary.evaluated),
            u8::from(summary.exempt)
        );
    }
    Ok(())
}
