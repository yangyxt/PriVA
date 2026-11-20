//! Standalone line protocol for the independently installed SciPy oracle.
#[path = "../../src/ks_greater.rs"]
mod ks_greater;
use std::io::{self, BufRead, Write};

fn main() {
    let mut output = io::BufWriter::new(io::stdout().lock());
    for line in io::stdin().lock().lines() {
        let line = line.unwrap();
        let mut words = line.split_whitespace();
        match words.next().unwrap() {
            "binom" => {
                let n = words.next().unwrap().parse().unwrap();
                let k = words.next().unwrap().parse().unwrap();
                writeln!(output, "ok {:016x}", ks_greater::binomial(n, k).to_bits()).unwrap();
            }
            "ks" => {
                let n: usize = words.next().unwrap().parse().unwrap();
                let values: Vec<_> = words
                    .map(|word| f64::from_bits(u64::from_str_radix(word, 16).unwrap()))
                    .collect();
                let value = ks_greater::greater(&values[..n], &values[n..]);
                writeln!(
                    output,
                    "ok {:016x} {:016x} {}",
                    value.statistic.to_bits(),
                    value.pvalue.to_bits(),
                    u8::from(value.fallback)
                )
                .unwrap();
            }
            _ => panic!("unknown probe command"),
        }
    }
}
