//! Standalone bit-pattern protocol for independent NumPy comparisons.
#[path = "../../src/exon_af_stats.rs"]
mod exon_af_stats;
use std::io::{self, BufRead};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args != ["calculate"] {
        return Err("usage: exon-af-stats-probe calculate < hexadecimal-f64-lines".into());
    }
    for line in io::stdin().lock().lines() {
        let values = line?
            .split_whitespace()
            .map(|word| u64::from_str_radix(word, 16).map(f64::from_bits))
            .collect::<Result<Vec<_>, _>>()?;
        if values.is_empty() {
            return Err("a sample list must not be empty".into());
        }
        let result = exon_af_stats::calculate(&values);
        println!(
            "{}",
            result
                .iter()
                .map(|value| format!("{:016x}", value.to_bits()))
                .collect::<Vec<_>>()
                .join(" ")
        );
    }
    Ok(())
}
