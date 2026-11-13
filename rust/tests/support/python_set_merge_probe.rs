//! Standalone observer for the explicit CPython string-set copy/update profile.
#[allow(dead_code)]
#[path = "../../src/python_set_order.rs"]
mod python_set_order;
use python_set_order::StringSet;
use std::collections::HashMap;
use std::io::{BufRead, Write};

fn decode(word: &str) -> Result<String, Box<dyn std::error::Error>> {
    if word == "-" {
        return Ok(String::new());
    }
    if !word.is_ascii() || !word.len().is_multiple_of(2) {
        return Err("invalid hex".into());
    }
    let bytes = (0..word.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&word[i..i + 2], 16))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(String::from_utf8(bytes)?)
}
fn encode(text: &str) -> String {
    if text.is_empty() {
        "-".into()
    } else {
        text.as_bytes().iter().map(|b| format!("{b:02x}")).collect()
    }
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 7
        || args[0] != "batch"
        || args[1] != "--input"
        || args[3] != "--output"
        || args[5] != "--seed"
    {
        return Err("usage: probe batch --input LINES --output LINES --seed UINT32".into());
    }
    let seed: u32 = args[6].parse()?;
    let input = std::io::BufReader::new(std::fs::File::open(&args[2])?);
    let mut output = std::io::BufWriter::new(std::fs::File::create(&args[4])?);
    let mut sets: HashMap<String, StringSet> = HashMap::new();
    for line in input.lines() {
        let line = line?;
        let words: Vec<_> = line.split_whitespace().collect();
        match words.as_slice() {
            ["clear"] => sets.clear(),
            ["new", id] => {
                if sets.insert((*id).into(), StringSet::new(seed)).is_some() {
                    return Err("duplicate set".into());
                }
            }
            ["add", id, value] => sets
                .get_mut(*id)
                .ok_or("missing set")?
                .insert(decode(value)?),
            ["copy", dest, source] => {
                let value = sets.get(*source).ok_or("missing set")?.copy();
                if sets.insert((*dest).into(), value).is_some() {
                    return Err("duplicate set".into());
                }
            }
            ["reload", dest, source] => {
                let mut value = StringSet::new(seed);
                for text in sets.get(*source).ok_or("missing set")?.iter() {
                    value.insert(text.into());
                }
                if sets.insert((*dest).into(), value).is_some() {
                    return Err("duplicate set".into());
                }
            }
            ["update", dest, source] if dest != source => {
                let mut value = sets.remove(*dest).ok_or("missing set")?;
                value.update_set(sets.get(*source).ok_or("missing set")?)?;
                sets.insert((*dest).into(), value);
            }
            ["order", label, id] => {
                let values: Vec<_> = sets
                    .get(*id)
                    .ok_or("missing set")?
                    .iter()
                    .map(encode)
                    .collect();
                writeln!(output, "{} {} {}", label, values.len(), values.join(" "))?;
            }
            ["different_seed"] => {
                let mut first = StringSet::new(seed);
                first.insert("left".into());
                let mut other = StringSet::new(seed.wrapping_add(1));
                other.insert("right".into());
                if first.update_set(&other).is_ok() || first.into_order() != ["left"] {
                    return Err("mixed profile was not refused atomically".into());
                }
                writeln!(output, "different_seed refused")?;
            }
            _ => return Err("invalid operation".into()),
        }
    }
    output.flush()?;
    Ok(())
}
