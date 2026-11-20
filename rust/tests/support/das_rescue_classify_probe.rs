//! Validation-only candidate classification protocol.
#[allow(dead_code)]
#[path = "../../src/das_rescue.rs"]
mod das_rescue;
#[path = "../../src/exon_af_stats.rs"]
mod exon_af_stats;
use std::io::{self, BufRead};

fn text(word: &str) -> Result<String, Box<dyn std::error::Error>> {
    if word == "-" {
        return Ok(String::new());
    }
    if word.len() % 2 != 0 {
        return Err("odd hex string".into());
    }
    let bytes = (0..word.len())
        .step_by(2)
        .map(|at| u8::from_str_radix(&word[at..at + 2], 16))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(String::from_utf8(bytes)?)
}
fn hex(value: &str) -> String {
    if value.is_empty() {
        "-".to_owned()
    } else {
        value
            .as_bytes()
            .iter()
            .map(|value| format!("{value:02x}"))
            .collect()
    }
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    if std::env::args().skip(1).collect::<Vec<_>>() != ["classify"] {
        return Err("usage: das-rescue-classify-probe classify".into());
    }
    for line in io::stdin().lock().lines() {
        let line = line?;
        let words: Vec<_> = line.split_whitespace().collect();
        if words.len() < 8 {
            return Err("short classification request".into());
        }
        let path = text(words[1])?;
        let scores_present = words[2].parse::<u8>()? != 0;
        let entry_present = words[3].parse::<u8>()? != 0;
        let entry = das_rescue::Entry {
            kind: text(words[4])?,
            has_go_molecular_function: words[5].parse::<u8>()? != 0,
        };
        let threshold = f64::from_bits(u64::from_str_radix(words[6], 16)?);
        let count: usize = words[7].parse()?;
        if words.len() != 8 + count {
            return Err("classification value count".into());
        }
        let values = words[8..]
            .iter()
            .map(|word| u64::from_str_radix(word, 16).map(f64::from_bits))
            .collect::<Result<Vec<_>, _>>()?;
        let scores = scores_present.then_some(values.as_slice());
        let first = entry_present.then_some(&entry);
        let selected = match words[0] {
            "go" => das_rescue::go_mf(scores, first, threshold),
            "non" => das_rescue::non_go_mf(&path, scores, first, threshold),
            _ => return Err("unknown classification family".into()),
        };
        let (action, kind) = match selected.action {
            das_rescue::Action::Skip => ("skip", String::new()),
            das_rescue::Action::Direct => ("direct", String::new()),
            das_rescue::Action::Test(kind) => ("test", kind),
        };
        println!("{} {action} {}", u8::from(selected.eligible), hex(&kind));
    }
    Ok(())
}
