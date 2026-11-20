//! Exact decimal integer pairs in; Python-category or binary64 bits out.
#[path = "../../src/das_integer_ratio.rs"]
mod das_integer_ratio;
use num_bigint::BigInt;
use std::io::{self, BufRead};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    if std::env::args().skip(1).collect::<Vec<_>>() != ["divide"] {
        return Err("usage: probe divide < integer-pairs.tsv".into());
    }
    for line in io::stdin().lock().lines() {
        let line = line?;
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 2 {
            return Err("two decimal integer fields required".into());
        }
        let a: BigInt = fields[0].parse()?;
        let b: BigInt = fields[1].parse()?;
        match das_integer_ratio::divide(&a, &b) {
            Ok(value) => println!("ok {:016x}", value.to_bits()),
            Err(error) => println!("err {}", error.category()),
        }
    }
    Ok(())
}
