//! Minimal independent build of the production string-set module for its Python differential test.
#[path = "../src/python_set_order.rs"]
#[allow(dead_code)] // This independent probe exercises only the string interface.
mod python_set_order;
use python_set_order::{StringHash, StringSet};
use std::collections::HashMap;
use std::fs;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let options: HashMap<&str, &str> = args
        .chunks_exact(2)
        .map(|pair| (pair[0].as_str(), pair[1].as_str()))
        .collect();
    if args.len() != 8
        || !["--input", "--output", "--hashes", "--seed"]
            .iter()
            .all(|key| options.contains_key(key))
    {
        return Err("usage: python-set-order-probe --input LINES --output LINES --hashes LINES --seed UINT32".into());
    }
    let seed = options["--seed"].parse()?;
    let mut set = StringSet::new(seed);
    let hash = StringHash::seeded(seed);
    let mut hashes = String::new();
    for text in fs::read_to_string(options["--input"])?.lines() {
        set.insert(text.to_owned());
        hashes.push_str(&format!("{}\n", hash.hash(text)));
    }
    let values = set.into_order();
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
