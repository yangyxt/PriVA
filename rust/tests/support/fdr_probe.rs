//! Exact-bit line protocol for comparison against the pinned SciPy call.
#[path = "../../src/fdr.rs"]
mod fdr;
use std::io::{self, BufRead, Write};

fn main() -> io::Result<()> {
    let mut out = io::BufWriter::new(io::stdout().lock());
    for line in io::stdin().lock().lines() {
        let line = line?;
        let values = line
            .split_whitespace()
            .map(|text| u64::from_str_radix(text, 16).map(f64::from_bits))
            .collect::<Result<Vec<_>, _>>()
            .map_err(io::Error::other)?;
        match fdr::benjamini_hochberg(&values) {
            Ok(adjusted) => {
                write!(out, "ok")?;
                for value in adjusted {
                    write!(out, " {:016x}", value.to_bits())?;
                }
                writeln!(out)?;
            }
            Err(message) => writeln!(out, "ValueError {message}")?,
        }
        out.flush()?;
    }
    Ok(())
}
