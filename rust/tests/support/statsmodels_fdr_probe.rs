//! Named standalone probe for the production statsmodels BH compatibility path.
#[path = "../../src/fdr.rs"]
#[allow(dead_code)]
mod fdr;
use std::io::{self, BufRead, Write};

fn main() -> io::Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args != ["fdr", "--profile", "statsmodels"] {
        return Err(io::Error::other(
            "usage: probe fdr --profile statsmodels < bits",
        ));
    }
    let mut output = io::BufWriter::new(io::stdout().lock());
    for line in io::stdin().lock().lines() {
        let values = line?
            .split_whitespace()
            .map(|text| u64::from_str_radix(text, 16).map(f64::from_bits))
            .collect::<Result<Vec<_>, _>>()
            .map_err(io::Error::other)?;
        let (&alpha, probabilities) = values
            .split_first()
            .ok_or_else(|| io::Error::other("missing alpha bits"))?;
        match fdr::statsmodels_bh(probabilities, alpha) {
            Ok((rejected, adjusted)) => {
                write!(output, "ok")?;
                for (reject, value) in rejected.into_iter().zip(adjusted) {
                    write!(output, " {}:{:016x}", u8::from(reject), value.to_bits())?;
                }
                writeln!(output)?;
            }
            Err(message) => writeln!(output, "ZeroDivisionError {message}")?,
        }
        output.flush()?;
    }
    Ok(())
}
