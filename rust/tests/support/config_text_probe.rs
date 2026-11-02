//! Standalone stdin/stdout probe for config_text_compare.py; no production CLI.
//! Build with rustc --edition=2024 --extern ryu=PATH, then run `serve`.
#[path = "../../src/config_text.rs"]
mod config_text;
#[allow(dead_code)]
#[path = "../../src/json_text.rs"]
mod json_text;
#[allow(dead_code)]
#[path = "../../src/python_text.rs"]
mod python_text;

// The shared whitespace table comes from the exact production module above.
// Its unrelated integer helpers import this strip wrapper, which the tested
// config writer does not call. No scientific/scoring module is substituted.
#[allow(dead_code)]
mod acmg_scoring {
    pub(crate) fn python_strip(text: &str) -> &str {
        text.trim_matches(crate::python_text::is_python_space)
    }
}

use json_text::Json;
use std::io::{self, BufRead, Write};

fn request(line: &str) -> Result<Json, String> {
    let request = json_text::parse(line, true)?;
    let original = request
        .get("text")
        .and_then(Json::as_str)
        .ok_or("text missing")?;
    let updates = request
        .get("updates")
        .and_then(Json::as_object)
        .ok_or("updates missing")?;
    let writer = match request.get("mode").and_then(Json::as_str) {
        Some("paths") => config_text::Writer::Paths,
        Some("reference") => config_text::Writer::Reference,
        _ => return Err("mode missing".into()),
    };
    Ok(match config_text::render(original, updates, writer) {
        Ok(output) => Json::Object(vec![
            ("ok".into(), Json::Bool(true)),
            ("output".into(), Json::Str(output)),
        ]),
        Err(error) => Json::Object(vec![
            ("ok".into(), Json::Bool(false)),
            ("kind".into(), Json::Str(error.kind.into())),
            ("error".into(), Json::Str(error.message)),
        ]),
    })
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    if std::env::args().skip(1).collect::<Vec<_>>() != ["serve"] {
        eprintln!("usage: config_text_probe serve < requests.jsonl > responses.jsonl");
        std::process::exit(2);
    }
    let mut out = io::BufWriter::new(io::stdout().lock());
    for line in io::stdin().lock().lines() {
        let response = request(&line?)?;
        writeln!(out, "{}", json_text::dumps(&response, false)?)?;
        out.flush()?;
    }
    Ok(())
}
