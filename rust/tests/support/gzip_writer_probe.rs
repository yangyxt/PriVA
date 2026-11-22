//! Standalone probe, compiled with the pinned flate2 library, for gzip_writer_compare.py.
//! `serve` reads tab-separated requests. Named arguments for `encode` are
//! --input PATH --output PATH --level N --mtime N --name-hex HEX --chunks N,N
//! --flush-every N --sink-chunk N --fault-at N --fault-kind error|zero|interrupt
//! --flush-fault error|interrupt. All output files use exclusive creation.
#[allow(dead_code)]
#[path = "../../src/gzip_writer.rs"]
mod gzip_writer;

use gzip_writer::{GzipOptions, GzipWriter};
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, BufRead, Write};
use std::rc::Rc;

#[derive(Default)]
struct Stats {
    bytes: u64,
    writes: usize,
    flushes: usize,
    fault_fired: bool,
    interruptions: usize,
}

struct Sink {
    file: File,
    stats: Rc<RefCell<Stats>>,
    chunk: usize,
    fault_at: Option<u64>,
    fault_kind: String,
    flush_fault: String,
    flush_fired: bool,
}

impl Write for Sink {
    fn write(&mut self, input: &[u8]) -> io::Result<usize> {
        let mut stats = self.stats.borrow_mut();
        stats.writes += 1;
        let mut count = input.len();
        if self.chunk > 0 {
            count = count.min(self.chunk);
        }
        if let Some(at) = self.fault_at.filter(|_| !stats.fault_fired) {
            if stats.bytes >= at {
                stats.fault_fired = true;
                return match self.fault_kind.as_str() {
                    "zero" => Ok(0),
                    "interrupt" => {
                        stats.interruptions += 1;
                        Err(io::ErrorKind::Interrupted.into())
                    }
                    _ => Err(io::ErrorKind::BrokenPipe.into()),
                };
            }
            count = count.min((at - stats.bytes) as usize);
        }
        let written = self.file.write(&input[..count])?;
        stats.bytes += written as u64;
        Ok(written)
    }

    fn flush(&mut self) -> io::Result<()> {
        let mut stats = self.stats.borrow_mut();
        stats.flushes += 1;
        if !self.flush_fault.is_empty() && !self.flush_fired {
            self.flush_fired = true;
            if self.flush_fault == "interrupt" {
                stats.interruptions += 1;
                return Err(io::ErrorKind::Interrupted.into());
            }
            return Err(io::ErrorKind::PermissionDenied.into());
        }
        self.file.flush()
    }
}

fn integer(options: &BTreeMap<&str, &str>, key: &str, default: u64) -> Result<u64, String> {
    options.get(key).map_or(Ok(default), |value| {
        value.parse().map_err(|_| format!("{key} must be unsigned"))
    })
}

fn report(stats: &Stats, first_error: Option<(&str, io::ErrorKind)>, finish_ok: bool) -> String {
    let (phase, error) = first_error.map_or((String::new(), String::new()), |(phase, kind)| {
        (phase.to_owned(), format!("{kind:?}"))
    });
    format!(
        "{{\"ok\":{},\"error_phase\":\"{phase}\",\"error_kind\":\"{error}\",\"finish_ok\":{finish_ok},\"bytes\":{},\"write_calls\":{},\"flush_calls\":{},\"fault_fired\":{},\"interruptions\":{}}}",
        first_error.is_none() && finish_ok,
        stats.bytes,
        stats.writes,
        stats.flushes,
        stats.fault_fired,
        stats.interruptions
    )
}

fn run(arguments: &[String]) -> Result<String, String> {
    if arguments.first().map(String::as_str) != Some("encode") || (arguments.len() - 1) % 2 != 0 {
        return Err("expected encode and named option/value pairs".into());
    }
    let mut options = BTreeMap::new();
    for pair in arguments[1..].chunks_exact(2) {
        if ![
            "--input",
            "--output",
            "--level",
            "--mtime",
            "--name-hex",
            "--chunks",
            "--flush-every",
            "--sink-chunk",
            "--fault-at",
            "--fault-kind",
            "--flush-fault",
        ]
        .contains(&pair[0].as_str())
            || options.insert(pair[0].as_str(), pair[1].as_str()).is_some()
        {
            return Err("unknown or duplicate option".into());
        }
    }
    let filename = if let Some(hex) = options.get("--name-hex") {
        if !hex.is_ascii() || hex.len() % 2 != 0 {
            return Err("filename must be hex bytes".into());
        }
        Some(
            (0..hex.len())
                .step_by(2)
                .map(|index| {
                    u8::from_str_radix(&hex[index..index + 2], 16)
                        .map_err(|_| "invalid filename hex".to_owned())
                })
                .collect::<Result<Vec<_>, _>>()?,
        )
    } else {
        None
    };
    let metadata = GzipOptions {
        level: u32::try_from(integer(&options, "--level", 9)?).map_err(|_| "level out of range")?,
        mtime: u32::try_from(integer(&options, "--mtime", 0)?).map_err(|_| "mtime out of range")?,
        filename,
    };
    let chunks = options
        .get("--chunks")
        .unwrap_or(&"65536")
        .split(',')
        .map(|value| value.parse::<usize>().map_err(|_| "bad chunk".to_owned()))
        .collect::<Result<Vec<_>, _>>()?;
    if chunks.is_empty() || chunks.contains(&0) {
        return Err("chunks must be positive".into());
    }
    let flush_every = integer(&options, "--flush-every", 0)?;
    let data = std::fs::read(options.get("--input").ok_or("--input missing")?)
        .map_err(|error| error.to_string())?;
    let file = File::create_new(options.get("--output").ok_or("--output missing")?)
        .map_err(|error| error.to_string())?;
    let stats = Rc::new(RefCell::new(Stats::default()));
    let sink = Sink {
        file,
        stats: Rc::clone(&stats),
        chunk: usize::try_from(integer(&options, "--sink-chunk", 0)?)
            .map_err(|_| "sink chunk out of range")?,
        fault_at: if options.contains_key("--fault-at") {
            Some(integer(&options, "--fault-at", 0)?)
        } else {
            None
        },
        fault_kind: options.get("--fault-kind").unwrap_or(&"error").to_string(),
        flush_fault: options.get("--flush-fault").unwrap_or(&"").to_string(),
        flush_fired: false,
    };
    let mut writer = match GzipWriter::new(sink, metadata) {
        Ok(writer) => writer,
        Err(error) => return Ok(report(&stats.borrow(), Some(("new", error.kind())), false)),
    };
    let mut first_error = None;
    let mut offset = 0usize;
    for (index, count) in chunks.iter().cycle().enumerate() {
        // Exercise empty updates too: they must not change the payload.
        if let Err(error) = writer.write_all(&[]) {
            first_error = Some(("write", error.kind()));
            break;
        }
        if offset == data.len() {
            break;
        }
        let end = data.len().min(offset.saturating_add(*count));
        if let Err(error) = writer.write_all(&data[offset..end]) {
            first_error = Some(("write", error.kind()));
            break;
        }
        offset = end;
        if flush_every > 0 && (index as u64 + 1) % flush_every == 0 {
            if let Err(error) = writer.flush() {
                first_error = Some(("flush", error.kind()));
                break;
            }
        }
    }
    let finish_ok = match writer.finish() {
        Ok(_) => true,
        Err(error) => {
            if first_error.is_none() {
                first_error = Some(("finish", error.kind()));
            }
            false
        }
    };
    Ok(report(&stats.borrow(), first_error, finish_ok))
}

fn main() {
    let arguments: Vec<_> = std::env::args().skip(1).collect();
    let result = if arguments == ["serve"] {
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
        run(&arguments).map(|text| println!("{text}"))
    };
    if let Err(error) = result {
        eprintln!("gzip_writer_probe: {error}");
        std::process::exit(2);
    }
}
