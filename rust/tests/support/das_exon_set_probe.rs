//! Batch interface for exact text-set reconstruction from pickle opcode order.
#[path = "../../src/python_set_order.rs"]
#[allow(dead_code)] // This probe exercises the production string interface only.
mod python_set_order;

use python_set_order::StringSet;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

fn decode(value: &str) -> Result<String, Box<dyn std::error::Error>> {
    if value == "-" {
        return Ok(String::new());
    }
    if !value.len().is_multiple_of(2) || !value.is_ascii() {
        return Err("invalid UTF-8 hex text".into());
    }
    let bytes = (0..value.len())
        .step_by(2)
        .map(|at| u8::from_str_radix(&value[at..at + 2], 16))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(String::from_utf8(bytes)?)
}

fn encode(value: &str) -> String {
    if value.is_empty() {
        "-".to_owned()
    } else {
        value
            .as_bytes()
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect()
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let options: HashMap<&str, &str> = args[usize::from(!args.is_empty())..]
        .chunks_exact(2)
        .map(|pair| (pair[0].as_str(), pair[1].as_str()))
        .collect();
    if args.first().map(String::as_str) != Some("batch")
        || args.len() != 7
        || options.len() != 3
        || !["--input", "--output", "--seed"]
            .iter()
            .all(|key| options.contains_key(key))
    {
        return Err(
            "usage: das-exon-set-probe batch --input LINES --output LINES --seed UINT32".into(),
        );
    }
    let seed = options["--seed"].parse()?;
    let input = BufReader::new(File::open(options["--input"])?);
    let mut output = BufWriter::new(File::create_new(options["--output"])?);
    for line in input.lines() {
        let line = line?;
        let mut words = line.split_whitespace();
        let count: usize = words.next().ok_or("missing item count")?.parse()?;
        let mut set = StringSet::new(seed);
        for _ in 0..count {
            set.insert(decode(words.next().ok_or("too few encoded items")?)?);
        }
        if words.next().is_some() {
            return Err("too many encoded items".into());
        }
        let ordered = set.into_order();
        write!(output, "{}", ordered.len())?;
        for value in ordered {
            write!(output, " {}", encode(&value))?;
        }
        writeln!(output)?;
    }
    output.flush()?;
    Ok(())
}
