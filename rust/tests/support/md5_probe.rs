//! Standalone probe for md5_compare.py, compiled directly with rustc.
//! Invoke `md5_probe serve` for tab-separated named-argument requests, or call
//! `bytes --input PATH --chunks 1,63,65`, `file --input PATH`,
//! `repeat --pattern HEX --length N --read-chunk N` directly.
#[allow(dead_code)]
#[path = "../../src/md5.rs"]
mod md5;

use std::collections::BTreeMap;
use std::io::{self, BufRead, Read, Write};
use std::path::Path;

struct PatternReader<'a> {
    pattern: &'a [u8],
    length: u64,
    position: u64,
    chunk: usize,
    interrupt_at: Option<u64>,
    interrupted: bool,
    fail_at: Option<u64>,
    calls: u64,
}

impl Read for PatternReader<'_> {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        self.calls += 1;
        if output.is_empty() {
            return Ok(0);
        }
        if !self.interrupted && self.interrupt_at.is_some_and(|at| self.position >= at) {
            self.interrupted = true;
            return Err(io::ErrorKind::Interrupted.into());
        }
        if self.fail_at.is_some_and(|at| self.position >= at) {
            return Err(io::ErrorKind::InvalidData.into());
        }
        let mut count = (self.length - self.position)
            .min(output.len() as u64)
            .min(self.chunk as u64);
        for at in [
            self.interrupt_at.filter(|_| !self.interrupted),
            self.fail_at,
        ]
        .into_iter()
        .flatten()
        {
            count = count.min(at.saturating_sub(self.position));
        }
        let count = count as usize;
        let mut filled = 0;
        while filled < count {
            let offset = (self.position % self.pattern.len() as u64) as usize;
            let take = (self.pattern.len() - offset).min(count - filled);
            output[filled..filled + take].copy_from_slice(&self.pattern[offset..offset + take]);
            self.position += take as u64;
            filled += take;
        }
        Ok(count)
    }
}

fn number(options: &BTreeMap<&str, &str>, key: &str, default: u64) -> Result<u64, String> {
    options.get(key).map_or(Ok(default), |value| {
        value
            .parse()
            .map_err(|_| format!("{key} must be an unsigned integer"))
    })
}

fn digest_result(result: io::Result<String>, name: &str) -> String {
    match result {
        Ok(value) => format!("{{\"ok\":true,\"{name}\":\"{value}\"}}"),
        Err(error) => format!("{{\"ok\":false,\"error_kind\":\"{:?}\"}}", error.kind()),
    }
}

fn reader_result(mut reader: PatternReader<'_>) -> String {
    match md5::md5_reader(&mut reader) {
        Ok(value) => format!(
            "{{\"ok\":true,\"reader\":\"{value}\",\"read_calls\":{},\"consumed\":{}}}",
            reader.calls, reader.position
        ),
        Err(error) => format!(
            "{{\"ok\":false,\"error_kind\":\"{:?}\",\"read_calls\":{},\"consumed\":{}}}",
            error.kind(),
            reader.calls,
            reader.position
        ),
    }
}

fn run(args: &[String]) -> Result<String, String> {
    let operation = args.first().ok_or("missing operation")?;
    if (args.len() - 1) % 2 != 0 {
        return Err("options require named value pairs".into());
    }
    let mut options = BTreeMap::new();
    for pair in args[1..].chunks_exact(2) {
        if ![
            "--input",
            "--chunks",
            "--pattern",
            "--length",
            "--read-chunk",
            "--interrupt-at",
            "--fail-at",
        ]
        .contains(&pair[0].as_str())
            || options.insert(pair[0].as_str(), pair[1].as_str()).is_some()
        {
            return Err("unknown or duplicate option".into());
        }
    }
    if operation == "file" {
        return Ok(digest_result(
            md5::md5_file(Path::new(options.get("--input").ok_or("--input missing")?)),
            "file",
        ));
    }
    let bytes = match operation.as_str() {
        "bytes" => match std::fs::read(options.get("--input").ok_or("--input missing")?) {
            Ok(bytes) => bytes,
            Err(error) => return Ok(digest_result(Err(error), "bytes")),
        },
        "repeat" | "fault" => {
            let hex = options.get("--pattern").ok_or("--pattern missing")?;
            if !hex.is_ascii() || hex.len() % 2 != 0 {
                return Err("pattern must be hex bytes".into());
            }
            (0..hex.len())
                .step_by(2)
                .map(|at| {
                    u8::from_str_radix(&hex[at..at + 2], 16).map_err(|_| "invalid hex byte".into())
                })
                .collect::<Result<Vec<_>, String>>()?
        }
        _ => return Err("unknown operation".into()),
    };
    let length = number(&options, "--length", bytes.len() as u64)?;
    let chunk = usize::try_from(number(&options, "--read-chunk", 1024 * 1024)?)
        .map_err(|_| "read chunk too large")?;
    if chunk == 0 || (bytes.is_empty() && length != 0) {
        return Err("nonempty streams need a pattern and positive read chunk".into());
    }
    let optional = |name: &str| -> Result<Option<u64>, String> {
        if options.contains_key(name) {
            number(&options, name, 0).map(Some)
        } else {
            Ok(None)
        }
    };
    let reader = PatternReader {
        pattern: &bytes,
        length,
        position: 0,
        chunk,
        interrupt_at: optional("--interrupt-at")?,
        interrupted: false,
        fail_at: optional("--fail-at")?,
        calls: 0,
    };
    if operation != "bytes" {
        return Ok(reader_result(reader));
    }
    let chunks = options
        .get("--chunks")
        .unwrap_or(&"64")
        .split(',')
        .map(|value| {
            value
                .parse::<usize>()
                .map_err(|_| "invalid update chunk".to_owned())
        })
        .collect::<Result<Vec<_>, _>>()?;
    if chunks.is_empty() || chunks.contains(&0) {
        return Err("update chunks must be positive".into());
    }
    let mut state = md5::Md5::new();
    let mut offset: usize = 0;
    for count in chunks.into_iter().cycle() {
        state.update(&[]);
        if offset == bytes.len() {
            break;
        }
        let end = bytes.len().min(offset.saturating_add(count));
        state.update(&bytes[offset..end]);
        offset = end;
    }
    let direct = md5::md5_hex(&bytes);
    let streamed = state.finish_hex();
    let read = md5::md5_reader(reader).map_err(|error| error.to_string())?;
    Ok(format!(
        "{{\"ok\":true,\"bytes\":\"{direct}\",\"stream\":\"{streamed}\",\"reader\":\"{read}\"}}"
    ))
}

fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let result = if args == ["serve"] {
        let mut output = io::BufWriter::new(io::stdout().lock());
        let mut result = Ok(());
        for line in io::stdin().lock().lines() {
            result = line
                .map_err(|error| error.to_string())
                .and_then(|line| run(&line.split('\t').map(str::to_owned).collect::<Vec<_>>()))
                .and_then(|text| {
                    writeln!(output, "{text}")
                        .and_then(|_| output.flush())
                        .map_err(|error| error.to_string())
                });
            if result.is_err() {
                break;
            }
        }
        result
    } else {
        run(&args).map(|text| println!("{text}"))
    };
    if let Err(error) = result {
        eprintln!("md5_probe: {error}");
        std::process::exit(2);
    }
}
