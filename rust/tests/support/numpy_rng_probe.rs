//! State-and-draw line protocol for the production NumPy random-stream helper.
#[path = "../../src/numpy_rng.rs"]
mod numpy_rng;
use numpy_rng::{Pcg64, seed_sequence};
use std::io::{self, BufRead, Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args != ["compare", "--profile", "numpy-1.26.4"] {
        return Err("usage: probe compare --profile numpy-1.26.4 < commands".into());
    }
    let mut output = io::BufWriter::new(io::stdout().lock());
    let mut generator = Pcg64::new(0);
    for line in io::stdin().lock().lines() {
        let line = line?;
        let parts: Vec<_> = line.split_whitespace().collect();
        let integer = |at: usize| -> Result<u64, Box<dyn std::error::Error>> {
            Ok(parts.get(at).ok_or("missing command field")?.parse()?)
        };
        let result: Result<Vec<u64>, &str> = match parts.first().copied() {
            Some("seed") => {
                generator = Pcg64::new(integer(1)?);
                Ok(vec![])
            }
            Some("words") => {
                let entropy: Vec<u32> = parts[2..]
                    .iter()
                    .map(|value| value.parse())
                    .collect::<Result<_, _>>()?;
                let words = seed_sequence(&entropy, integer(1)? as usize);
                write!(output, "words")?;
                for word in words {
                    write!(output, " {word:016x}")?;
                }
                writeln!(output)?;
                continue;
            }
            Some("raw") => Ok((0..integer(1)?).map(|_| generator.next_u64()).collect()),
            Some("uint32") => Ok((0..integer(1)?)
                .map(|_| generator.next_u32() as u64)
                .collect()),
            Some("bounded") => {
                let maximum = integer(1)?;
                Ok((0..integer(2)?)
                    .map(|_| generator.bounded(maximum))
                    .collect())
            }
            Some("choice") => generator
                .choice(integer(1)? as usize, integer(2)? as usize, integer(3)? != 0)
                .map(|values| values.into_iter().map(|value| value as u64).collect()),
            _ => return Err("unknown command".into()),
        };
        match result {
            Ok(values) => {
                write!(output, "ok")?;
                for value in values {
                    write!(output, " {value:016x}")?;
                }
            }
            Err(message) => write!(output, "ValueError {message}")?,
        }
        let (state, increment, has_uint32, buffer) = generator.state();
        writeln!(
            output,
            "|{state:032x}:{increment:032x}:{}:{buffer:08x}",
            u8::from(has_uint32)
        )?;
    }
    Ok(())
}
