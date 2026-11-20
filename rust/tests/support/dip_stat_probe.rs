//! Line-oriented development probe; no diptest library is linked.
#[path = "../../src/dip_stat.rs"]
mod dip_stat;
use std::io::{self, BufRead, Write};

fn main() {
    let mut output = io::BufWriter::new(io::stdout().lock());
    for line in io::stdin().lock().lines() {
        let line = line.unwrap();
        let values: Vec<_> = line
            .split_whitespace()
            .map(|word| f64::from_bits(u64::from_str_radix(word, 16).unwrap()))
            .collect();
        match dip_stat::dip(&values) {
            Ok(value) => writeln!(output, "ok {:016x}", value.to_bits()).unwrap(),
            Err(_) => writeln!(output, "error RuntimeError").unwrap(),
        }
    }
}
