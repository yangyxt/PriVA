//! Standalone JSON writer probe for json_text_compare.py; no production CLI.
//! Compile with rustc --edition=2024 and --extern ryu=PATH, then invoke
//! `json_text_probe serve`.
//! `--cfg legacy` uses only the pre-options wrappers for compatibility checks.
#[allow(dead_code)]
#[path = "../../src/json_text.rs"]
mod json_text;

use json_text::Json;
use std::io::{self, BufRead, Write};

fn flag(value: &Json, key: &str, default: bool) -> Result<bool, String> {
    match value.get(key) {
        None => Ok(default),
        Some(Json::Bool(flag)) => Ok(*flag),
        _ => Err(format!("probe option {key} must be boolean")),
    }
}

fn failure(phase: &str, error: String) -> Json {
    Json::Object(vec![
        ("ok".into(), Json::Bool(false)),
        ("phase".into(), Json::Str(phase.into())),
        ("refusal".into(), Json::Bool(json_text::is_refusal(&error))),
        ("error".into(), Json::Str(error)),
    ])
}

fn request(line: &str) -> Result<Json, String> {
    let request = json_text::parse(line, true)?;
    let input = request
        .get("input")
        .and_then(Json::as_str)
        .ok_or("probe input must be text")?;
    let value = match json_text::parse(input, flag(&request, "reject_duplicates", false)?) {
        Ok(value) => value,
        Err(error) => return Ok(failure("parse", error)),
    };
    let mode = request
        .get("mode")
        .and_then(Json::as_str)
        .unwrap_or("options");
    let result = match mode {
        "legacy_default" => json_text::dumps(&value, false),
        "legacy_sorted" => json_text::dumps(&value, true),
        "legacy_indented" => json_text::dumps_indented(&value),
        #[cfg(not(legacy))]
        "options" => {
            let options = request.get("options").ok_or("probe options missing")?;
            let indent = match options.get("indent") {
                None | Some(Json::Null) => None,
                Some(Json::Int(width)) => {
                    Some(usize::try_from(*width).map_err(|_| "probe indent must be nonnegative")?)
                }
                _ => return Err("probe indent must be an integer or null".into()),
            };
            json_text::dumps_with_options(
                &value,
                json_text::DumpsOptions {
                    indent,
                    sort_keys: flag(options, "sort_keys", false)?,
                    ensure_ascii: flag(options, "ensure_ascii", true)?,
                    compact: flag(options, "compact", false)?,
                    allow_nan: flag(options, "allow_nan", true)?,
                },
            )
        }
        _ => return Err("unknown probe mode".into()),
    };
    Ok(match result {
        Ok(output) => Json::Object(vec![
            ("ok".into(), Json::Bool(true)),
            ("output".into(), Json::Str(output)),
        ]),
        Err(error) => failure("dump", error),
    })
}

fn serve() -> Result<(), String> {
    let stdin = io::stdin();
    let mut stdout = io::BufWriter::new(io::stdout().lock());
    for line in stdin.lock().lines() {
        let result = request(&line.map_err(|error| error.to_string())?)?;
        let encoded = json_text::dumps(&result, false)?;
        writeln!(stdout, "{encoded}").map_err(|error| error.to_string())?;
        stdout.flush().map_err(|error| error.to_string())?;
    }
    Ok(())
}

fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args != ["serve"] {
        eprintln!("usage: json_text_probe serve < requests.jsonl > responses.jsonl");
        std::process::exit(2);
    }
    if let Err(error) = serve() {
        eprintln!("{error}");
        std::process::exit(1);
    }
}
