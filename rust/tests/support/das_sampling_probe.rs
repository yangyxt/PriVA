//! Development-only line protocol for independently executing the Python procedures.
#[path = "../../src/das_sampling.rs"]
mod das_sampling;
#[path = "../../src/dip_stat.rs"]
mod dip_stat;
#[path = "../../src/exon_af_stats.rs"]
mod exon_af_stats;
#[allow(dead_code)]
#[path = "../../src/ks_greater.rs"]
mod ks_greater;
#[path = "../../src/numpy_rng.rs"]
mod numpy_rng;
use std::io::{self, BufRead, Write};

fn sample(words: &mut std::str::SplitWhitespace<'_>) -> Vec<f64> {
    let count: usize = words.next().unwrap().parse().unwrap();
    words
        .take(count)
        .map(|word| f64::from_bits(u64::from_str_radix(word, 16).unwrap()))
        .collect()
}

fn main() {
    let mut output = io::BufWriter::new(io::stdout().lock());
    for line in io::stdin().lock().lines() {
        let line = line.unwrap();
        let mut words = line.split_whitespace();
        let command = words.next().unwrap();
        let seed: u64 = words.next().unwrap().parse().unwrap();
        let mut rng = numpy_rng::Pcg64::new(seed);
        let result = match command {
            "composite" => {
                let target = words.next().unwrap().parse().unwrap();
                let count: usize = words.next().unwrap().parse().unwrap();
                let references: Vec<_> = (0..count).map(|_| sample(&mut words)).collect();
                das_sampling::fixed_composite(&references, target, &mut rng).map(|values| {
                    values
                        .iter()
                        .map(|value| format!("{:016x}", value.to_bits()))
                        .collect::<Vec<_>>()
                        .join(" ")
                })
            }
            "ks" => {
                let cap = words.next().unwrap();
                let cap = if cap == "none" {
                    None
                } else {
                    Some(cap.parse().unwrap())
                };
                let query = sample(&mut words);
                let composite = sample(&mut words);
                das_sampling::ks_sample(&query, &composite, &mut rng, cap).map(|value| {
                    format!(
                        "{:016x} {:016x} {} {} {:016x} {:016x}",
                        value.statistic.to_bits(),
                        value.pvalue.to_bits(),
                        value.n_query,
                        value.n_composite,
                        value.median_query.to_bits(),
                        value.median_composite.to_bits()
                    )
                })
            }
            "dip" => {
                let permutations = words.next().unwrap().parse().unwrap();
                let query = sample(&mut words);
                let composite = sample(&mut words);
                das_sampling::dip_rescue(&query, &composite, permutations, &mut rng).map(|value| {
                    format!(
                        "{:016x} {:016x} {:016x} {:016x} {}",
                        value.query.to_bits(),
                        value.reference_mean.to_bits(),
                        value.reference_maximum.to_bits(),
                        value.pvalue.to_bits(),
                        value.n_query
                    )
                })
            }
            _ => panic!("unknown procedure"),
        };
        assert!(words.next().is_none());
        let value = match result {
            Ok(value) => format!("ok {value}"),
            Err(error) => format!("error {} {}", error.category, error.message),
        };
        let (state, increment, cached, word) = rng.state();
        writeln!(
            output,
            "{value} | {state} {increment} {} {word}",
            u8::from(cached)
        )
        .unwrap();
    }
}
