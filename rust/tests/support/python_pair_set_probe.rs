//! Independent, dependency-free build of the shared Python set-order module.
#[path = "../../src/python_set_order.rs"]
mod python_set_order;
use python_set_order::{PairSet, StringHash, StringSet};
use std::collections::HashMap;
use std::fs;

fn unhex(value: &str) -> Result<String, Box<dyn std::error::Error>> {
    if value == "-" {
        return Ok(String::new());
    }
    if value.len() % 2 != 0 {
        return Err("odd text encoding".into());
    }
    let bytes = (0..value.len())
        .step_by(2)
        .map(|at| u8::from_str_radix(&value[at..at + 2], 16))
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
            .map(|b| format!("{b:02x}"))
            .collect()
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let pairs = args.first().is_some_and(|s| s == "pairs");
    if pairs {
        args.remove(0);
    }
    let options: HashMap<&str, &str> = args
        .chunks_exact(2)
        .map(|pair| (pair[0].as_str(), pair[1].as_str()))
        .collect();
    if args.len() != 8
        || !["--input", "--output", "--hashes", "--seed"]
            .iter()
            .all(|key| options.contains_key(key))
    {
        return Err(
            "usage: probe [pairs] --input LINES --output LINES --hashes LINES --seed UINT32".into(),
        );
    }
    let seed = options["--seed"].parse()?;
    let hash = StringHash::seeded(seed);
    let mut hashes = String::new();
    let values = if pairs {
        let mut set = PairSet::new(seed);
        for text in fs::read_to_string(options["--input"])?.lines() {
            let Some((first, second)) = text.split_once(' ') else {
                return Err("two encoded texts required".into());
            };
            let pair = (unhex(first)?, unhex(second)?);
            hashes.push_str(&format!("{}\n", hash.pair_hash(&pair.0, &pair.1)));
            set.insert(pair);
        }
        set.into_order()
            .into_iter()
            .map(|(a, b)| format!("{} {}", hex(&a), hex(&b)))
            .collect::<Vec<_>>()
    } else {
        // Preserve the original string oracle's protocol unchanged.
        let mut set = StringSet::new(seed);
        for text in fs::read_to_string(options["--input"])?.lines() {
            set.insert(text.to_owned());
            hashes.push_str(&format!("{}\n", hash.hash(text)));
        }
        set.into_order()
    };
    fs::write(
        options["--output"],
        if values.is_empty() {
            String::new()
        } else {
            values.join("\n") + "\n"
        },
    )?;
    fs::write(options["--hashes"], hashes)?;
    Ok(())
}
