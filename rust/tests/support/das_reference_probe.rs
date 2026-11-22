//! Typed metadata protocol for independent original-Python comparison.
#[allow(dead_code)]
#[path = "../../src/das_reference.rs"]
mod das_reference;
#[allow(dead_code)]
#[path = "../../src/das_rescue.rs"]
mod das_rescue;
#[allow(dead_code)]
#[path = "../../src/exon_af_stats.rs"]
mod exon_af_stats;
use std::io::{self, BufRead};

fn hex(value: &str) -> String {
    if value.is_empty() {
        "-".into()
    } else {
        value
            .as_bytes()
            .iter()
            .map(|value| format!("{value:02x}"))
            .collect()
    }
}
fn text(word: &str) -> Result<String, Box<dyn std::error::Error>> {
    if word == "-" {
        return Ok(String::new());
    }
    if word.len() % 2 != 0 {
        return Err("odd encoded text".into());
    }
    Ok(String::from_utf8(
        (0..word.len())
            .step_by(2)
            .map(|at| u8::from_str_radix(&word[at..at + 2], 16))
            .collect::<Result<Vec<_>, _>>()?,
    )?)
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    if std::env::args().skip(1).collect::<Vec<_>>() != ["classify"] {
        return Err("usage: probe classify".into());
    }
    for line in io::stdin().lock().lines() {
        let line = line?;
        let words: Vec<_> = line.split_whitespace().collect();
        if words.len() != 10 {
            return Err("ten metadata fields required".into());
        }
        let entry = das_rescue::Entry {
            kind: text(words[2])?,
            has_go_molecular_function: words[3] == "1",
        };
        let entry = (words[1] == "1").then_some(&entry);
        let metadata = (words[4] != "missing").then_some(das_reference::Metadata {
            start: words[5].parse()?,
            end: words[6].parse()?,
            protein_length: if words[4] == "null" {
                None
            } else {
                Some(words[7].parse()?)
            },
        });
        let threshold = f64::from_bits(u64::from_str_radix(words[8], 16)?);
        let minimum = words[9].parse()?;
        match words[0] {
            "reference" => match das_reference::reference(entry, metadata, threshold, minimum) {
                Ok(value) => println!(
                    "ok {}",
                    match value {
                        das_reference::Reference::Unknown => "unknown",
                        das_reference::Reference::NonFunctional => "non",
                        das_reference::Reference::Functional => "functional",
                    }
                ),
                Err(error) => println!("err {}", error.category()),
            },
            "family" => match das_reference::family(entry, metadata, threshold, minimum) {
                Ok(value) => println!(
                    "ok {} {} {} {}",
                    value.kind.as_deref().map_or("*".into(), hex),
                    value
                        .coverage
                        .map_or("*".into(), |v| format!("{:016x}", v.to_bits())),
                    value.protein_length.map_or("*".into(), |v| v.to_string()),
                    u8::from(value.removed)
                ),
                Err(error) => println!("err {}", error.category()),
            },
            _ => return Err("unknown procedure".into()),
        }
    }
    Ok(())
}
