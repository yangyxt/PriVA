//! Native counterpart of scripts/indexed_plugin_cache.py.
//!
//! Retain the mature bgzip/tabix/sort commands, their argument and mutation
//! order, publisher coordinate settings, source identities and staged publish.
//! These optional procedures acquire no source and replace no selected data.
use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use crate::python_args::{parse_long_options_converted, require_options};
use crate::python_io::{GzipFailure, PythonGzipReader};
use crate::python_path::{absolute, is_file, metadata, os_error, path_value, repr};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Read};
use std::os::unix::fs::MetadataExt;
use std::os::unix::process::ExitStatusExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

mod primate;
mod temporary;
use temporary::TemporaryDirectory;

const HELP: &str = "usage: priva indexed-plugin-cache {ensure-index|format-primate} --input FILE [OPTIONS]\n\
ensure-index --kind {mavedb,splicevault,primateai}\n\
format-primate --cache-dir DIRECTORY --assembly {hg19,GRCh37,hg38,GRCh38} [--threads N] [--tmp-dir DIRECTORY] [--print-path]\n\
Use selected prepared publisher files. Original installer callers remain unchanged.";
const PROBE: usize = 1024 * 1024;
type Settings = [i32; 6];
type Identity = (u64, u64, u64, i64, i64, i64, i64);

fn error(category: &'static str, message: impl Into<String>) -> Box<dyn std::error::Error> {
    crate::builder_error::category(category, invalid(message))
}
fn require(condition: bool, message: impl Into<String>) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(error("ValueError", message))
    }
}
fn identity(path: &Path) -> Result<Identity> {
    let m = metadata(path)?;
    Ok((
        m.dev(),
        m.ino(),
        m.size(),
        m.mtime(),
        m.mtime_nsec(),
        m.ctime(),
        m.ctime_nsec(),
    ))
}
fn file(path: &Path) -> Result<File> {
    File::open(path).map_err(|e| os_error(e, path, None))
}
fn exists(path: &Path) -> Result<bool> {
    match fs::metadata(path) {
        Ok(_) => Ok(true),
        Err(e)
            if matches!(
                e.raw_os_error(),
                Some(libc::ENOENT | libc::ENOTDIR | libc::EBADF | libc::ELOOP)
            ) =>
        {
            Ok(false)
        }
        Err(e) => Err(os_error(e, path, None)),
    }
}
fn gzip(path: &Path) -> Result<PythonGzipReader<File>> {
    Ok(PythonGzipReader::new(file(path)?))
}
fn gzip_error(cause: Box<dyn std::error::Error>) -> Box<dyn std::error::Error> {
    if let Some(io) = cause.downcast_ref::<io::Error>() {
        let name = match GzipFailure::of(io) {
            Some(GzipFailure::EndOfFile) => Some("EOFError"),
            Some(GzipFailure::Bad(_)) => Some("BadGzipFile"),
            Some(GzipFailure::Zlib(_)) => Some("error"),
            None => None,
        };
        if let Some(name) = name {
            return error(name, cause.to_string());
        }
    }
    cause
}
fn suffix(path: &Path, extra: &str) -> PathBuf {
    PathBuf::from(format!("{}{extra}", path.display()))
}
fn default_settings(kind: &str) -> Settings {
    [
        0,
        1,
        2,
        if kind == "mavedb" { 3 } else { 2 },
        i32::from(b'#'),
        0,
    ]
}
fn object(values: impl IntoIterator<Item = (&'static str, Json)>) -> Json {
    Json::Object(values.into_iter().map(|(k, v)| (k.into(), v)).collect())
}
fn index_result(path: &Path, repaired: bool) -> Json {
    object([
        ("path", Json::Str(path.display().to_string())),
        (
            "index",
            Json::Str(suffix(path, ".tbi").display().to_string()),
        ),
        ("repaired", Json::Bool(repaired)),
    ])
}
fn command_error(argv: &[String], status: std::process::ExitStatus) -> Box<dyn std::error::Error> {
    let command = format!(
        "[{}]",
        argv.iter().map(|s| repr(s)).collect::<Vec<_>>().join(", ")
    );
    let message = if let Some(code) = status.code() {
        format!("Command '{command}' returned non-zero exit status {code}.")
    } else {
        let number = status.signal().unwrap_or(1);
        let signal = crate::python_subprocess::signal_repr(number);
        format!("Command '{command}' died with {signal}.")
    };
    error("CalledProcessError", message)
}
fn command(argv: &[String], capture: bool) -> Result<Option<Vec<u8>>> {
    let mut cmd = Command::new(&argv[0]);
    cmd.args(&argv[1..]);
    if capture {
        let output = cmd
            .output()
            .map_err(|e| os_error(e, Path::new(&argv[0]), None))?;
        if output.status.success() {
            Ok(Some(output.stdout))
        } else {
            Err(command_error(argv, output.status))
        }
    } else {
        let status = cmd
            .status()
            .map_err(|e| os_error(e, Path::new(&argv[0]), None))?;
        if status.success() {
            Ok(None)
        } else {
            Err(command_error(argv, status))
        }
    }
}

/// Decimal coordinates retain Python's arbitrary integer range. They are only
/// compared and adjusted by one before passing their original semantics to
/// tabix; no fixed-width narrowing or scientific rounding is introduced.
#[derive(Clone, Debug, Eq, PartialEq)]
struct Coordinate(String);
impl Ord for Coordinate {
    fn cmp(&self, rhs: &Self) -> std::cmp::Ordering {
        self.0
            .len()
            .cmp(&rhs.0.len())
            .then_with(|| self.0.cmp(&rhs.0))
    }
}
impl PartialOrd for Coordinate {
    fn partial_cmp(&self, rhs: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(rhs))
    }
}
impl Coordinate {
    fn parse(text: &str) -> Option<(bool, Self)> {
        if !crate::python_text::python_int_is_valid(text) {
            return None;
        }
        let stripped = text.trim_matches(crate::python_text::is_python_space);
        let digits: String = stripped
            .chars()
            .filter_map(crate::python_text::decimal_digit_value)
            .map(|v| char::from(b'0' + v as u8))
            .collect();
        let digits = digits.trim_start_matches('0');
        let number = Self(if digits.is_empty() {
            "0".into()
        } else {
            digits.into()
        });
        Some((stripped.starts_with('-') && number.0 != "0", number))
    }
    fn zero() -> Self {
        Self("0".into())
    }
    fn one() -> Self {
        Self("1".into())
    }
    fn next(&self) -> Self {
        let mut bytes = self.0.as_bytes().to_vec();
        for b in bytes.iter_mut().rev() {
            if *b != b'9' {
                *b += 1;
                return Self(String::from_utf8(bytes).unwrap());
            }
            *b = b'0';
        }
        bytes.insert(0, b'1');
        Self(String::from_utf8(bytes).unwrap())
    }
    fn previous(&self) -> Self {
        if self.0 == "0" {
            return self.clone();
        }
        let mut bytes = self.0.as_bytes().to_vec();
        for b in bytes.iter_mut().rev() {
            if *b != b'0' {
                *b -= 1;
                break;
            }
            *b = b'9';
        }
        let text = String::from_utf8(bytes).unwrap();
        Self(text.trim_start_matches('0').to_owned()).max(Self::zero())
    }
}
struct Record {
    chromosome: String,
    start: Coordinate,
    end: Coordinate,
    line: String,
}

fn index_settings(path: &Path) -> Result<Option<Settings>> {
    let read = || -> Result<Vec<u8>> {
        let mut bytes = Vec::new();
        gzip(path)?.take(36).read_to_end(&mut bytes)?;
        Ok(bytes)
    };
    let header = match read() {
        Ok(header) => header,
        Err(cause) => {
            // Python catches OSError/EOFError/ValueError here, but not zlib.error.
            if cause
                .downcast_ref::<io::Error>()
                .and_then(GzipFailure::of)
                .is_some_and(|v| matches!(v, GzipFailure::Zlib(_)))
            {
                return Err(gzip_error(cause));
            }
            return Ok(None);
        }
    };
    if header.len() != 36 || &header[..4] != b"TBI\x01" {
        return Ok(None);
    }
    let ints: Vec<i32> = header[4..]
        .chunks_exact(4)
        .map(|p| i32::from_le_bytes(p.try_into().unwrap()))
        .collect();
    let (n, fmt, seq, begin, end, meta, skip, names) = (
        ints[0], ints[1], ints[2], ints[3], ints[4], ints[5], ints[6], ints[7],
    );
    if n < 1
        || fmt < 0
        || seq.min(begin).min(end) < 1
        || !(0..=127).contains(&meta)
        || skip < 0
        || names < 1
    {
        return Ok(None);
    }
    Ok(Some([fmt, seq, begin, end, meta, skip]))
}
fn first_record(path: &Path, s: Settings) -> Result<Record> {
    let mut lines = crate::python_lines::Lines::new(gzip(path)?);
    for number in 0..4096 {
        let line = lines.readline(PROBE).map_err(gzip_error)?;
        if line.is_empty() {
            break;
        }
        require(
            line.ends_with('\n') || line.chars().count() < PROBE,
            "Plugin table header/record exceeds the bounded probe",
        )?;
        if number < s[5]
            || line
                .trim_matches(crate::python_text::is_python_space)
                .is_empty()
            || (s[4] != 0 && line.starts_with(char::from_u32(s[4] as u32).unwrap()))
        {
            continue;
        }
        let text = line.trim_end_matches(['\r', '\n']);
        let fields: Vec<&str> = text.split('\t').collect();
        let parsed = || -> Option<(&str, (bool, Coordinate), (bool, Coordinate))> {
            Some((
                *fields.get((s[1] - 1) as usize)?,
                Coordinate::parse(fields.get((s[2] - 1) as usize)?)?,
                Coordinate::parse(fields.get((s[3] - 1) as usize)?)?,
            ))
        };
        let (chrom, (negative_start, start), (negative_end, end)) = parsed().ok_or_else(|| {
            error(
                "ValueError",
                "Plugin table has no valid first coordinate record",
            )
        })?;
        let zero_based = s[0] & 0x10000 != 0;
        require(
            !chrom.is_empty()
                && !negative_start
                && !negative_end
                && (zero_based || start >= Coordinate::one())
                && end >= start,
            "Plugin table coordinates disagree with the TBI coordinate settings",
        )?;
        return Ok(Record {
            chromosome: chrom.into(),
            start: if zero_based {
                start.next()
            } else {
                start.clone()
            },
            end: end.max(start.next()),
            line: text.into(),
        });
    }
    Err(error(
        "ValueError",
        "Plugin table contains no bounded first data record",
    ))
}
fn boundary_queries(path: &Path, record: &Record) -> Result<Vec<(String, Vec<u8>)>> {
    let positions = std::collections::BTreeSet::from([
        record.start.previous().max(Coordinate::one()),
        record.start.clone().max(Coordinate::one()),
        record.start.next(),
        record.end.clone().max(Coordinate::one()),
        record.end.next(),
    ]);
    positions
        .into_iter()
        .map(|position| {
            let region = format!("{}:{}-{}", record.chromosome, position.0, position.0);
            let bytes = command(
                &["tabix".into(), path.display().to_string(), region.clone()],
                true,
            )?
            .unwrap();
            Ok((region, bytes))
        })
        .collect()
}
fn readable_queries(path: &Path, record: &Record) -> Result<Option<Vec<(String, Vec<u8>)>>> {
    let run = || -> Result<Vec<(String, Vec<u8>)>> {
        command(
            &["tabix".into(), "-l".into(), path.display().to_string()],
            true,
        )?;
        boundary_queries(path, record)
    };
    match run() {
        Ok(values) => {
            // bytes.splitlines recognizes CR/LF only, including a final CR.
            let found = values.iter().any(|(_, value)| {
                value
                    .split(|b| *b == b'\n' || *b == b'\r')
                    .any(|line| line == record.line.as_bytes())
            });
            Ok(found.then_some(values))
        }
        Err(cause) if cause.to_string().starts_with("CalledProcessError: ") => Ok(None),
        Err(cause) => Err(cause),
    }
}

fn ensure_index(input: &Path, kind: &str) -> Result<Json> {
    let path = PathBuf::from(absolute(input)?);
    require(
        is_file(&path)? && metadata(&path)?.len() != 0,
        format!("Missing or empty selected {kind} table: {}", path.display()),
    )?;
    let before = identity(&path)?;
    let mut prefix = Vec::new();
    file(&path)?.take(18).read_to_end(&mut prefix)?;
    require(
        prefix.len() >= 18
            && &prefix[..4] == b"\x1f\x8b\x08\x04"
            && &prefix[12..16] == b"BC\x02\x00",
        format!("Selected {kind} table is not BGZF: {}", path.display()),
    )?;
    command(
        &["bgzip".into(), "-t".into(), path.display().to_string()],
        false,
    )?;
    let index = suffix(&path, ".tbi");
    let settings = index_settings(&index)?.unwrap_or(default_settings(kind));
    require(
        settings[0] & !0x10000 == 0,
        format!("Selected {kind} TBI has an unsupported non-generic format"),
    )?;
    require(
        settings[1..4] == default_settings(kind)[1..4],
        format!("Selected {kind} TBI columns disagree with the plugin's columns"),
    )?;
    let record = first_record(&path, settings)?;
    let old = if is_file(&index)? {
        readable_queries(&path, &record)?
    } else {
        None
    };
    if old.is_some() {
        let a = metadata(&index)?;
        let b = metadata(&path)?;
        if (a.mtime(), a.mtime_nsec()) >= (b.mtime(), b.mtime_nsec()) {
            require(
                identity(&path)? == before,
                "Selected table changed during index validation",
            )?;
            return Ok(index_result(&path, false));
        }
    }
    let temporary = TemporaryDirectory::new(
        path.parent().unwrap(),
        &format!(".{}.index-", path.file_name().unwrap().to_string_lossy()),
    )?;
    let outcome = (|| {
        let staged = temporary.0.join(path.file_name().unwrap());
        std::os::unix::fs::symlink(&path, &staged)
            .map_err(|e| os_error(e, &path, Some(&staged)))?;
        let mut args = vec![
            "tabix".into(),
            "-s".into(),
            settings[1].to_string(),
            "-b".into(),
            settings[2].to_string(),
            "-e".into(),
            settings[3].to_string(),
            "-c".into(),
            if settings[4] == 0 {
                String::new()
            } else {
                char::from_u32(settings[4] as u32).unwrap().to_string()
            },
            "-S".into(),
            settings[5].to_string(),
        ];
        if settings[0] & 0x10000 != 0 {
            args.push("-0".into());
        }
        args.push(staged.display().to_string());
        command(&args, false)?;
        let new = readable_queries(&staged, &record)?;
        require(
            new.is_some(),
            "Repaired index cannot retrieve the first data record",
        )?;
        require(
            old.is_none() || new == old,
            "Repaired index changed bounded coordinate-boundary query results",
        )?;
        require(
            index_settings(&suffix(&staged, ".tbi"))? == Some(settings),
            "Repaired TBI did not preserve publisher coordinate/header settings",
        )?;
        require(
            identity(&path)? == before,
            "Selected table changed during index repair",
        )?;
        let staged_index = suffix(&staged, ".tbi");
        fs::rename(&staged_index, &index).map_err(|e| os_error(e, &staged_index, Some(&index)))?;
        Ok(index_result(&path, true))
    })();
    temporary.finish(outcome)
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(action) = args.first() else {
        return Err(invalid(HELP));
    };
    if action == "--help" || action == "-h" {
        println!("{HELP}");
        return Ok(());
    }
    let (names, flags): (&[&str], &[&str]) = match action.as_str() {
        "ensure-index" => (&["input", "kind"], &[]),
        "format-primate" => (
            &["input", "cache-dir", "assembly", "threads", "tmp-dir"],
            &["print-path"],
        ),
        _ => return Err(invalid(HELP)),
    };
    let parsed = parse_long_options_converted(&args[1..], names, flags, |name, value| {
        if name == "threads" && Coordinate::parse(value).is_none() {
            return Err(format!(
                "argument --threads: invalid int value: {}",
                repr(value)
            ));
        }
        let choices = match name {
            "kind" => &["mavedb", "splicevault", "primateai"][..],
            "assembly" => &["hg19", "GRCh37", "hg38", "GRCh38"][..],
            _ => &[],
        };
        if !choices.is_empty() && !choices.contains(&value) {
            return Err(format!(
                "argument --{name}: invalid choice: {} (choose from {})",
                repr(value),
                choices
                    .iter()
                    .map(|s| repr(s))
                    .collect::<Vec<_>>()
                    .join(", ")
            ));
        }
        Ok(())
    })
    .map_err(invalid)?;
    match parsed {
        None => println!("{HELP}"),
        Some(occurrences) => {
            let values: std::collections::HashMap<String, String> =
                occurrences.into_iter().collect();
            let required = if action == "ensure-index" {
                &["input", "kind"][..]
            } else {
                &["input", "cache-dir", "assembly"][..]
            };
            require_options(&values, required).map_err(invalid)?;
            let input = path_value(&values["input"]);
            let result = if action == "ensure-index" {
                let kind = &values["kind"];
                if !["mavedb", "splicevault", "primateai"].contains(&kind.as_str()) {
                    return Err(invalid("argument --kind: invalid choice"));
                }
                crate::builder_error::finish(ensure_index(&input, kind))?
            } else {
                let assembly = &values["assembly"];
                if !["hg19", "GRCh37", "hg38", "GRCh38"].contains(&assembly.as_str()) {
                    return Err(invalid("argument --assembly: invalid choice"));
                }
                let threads = values.get("threads").map_or("1", String::as_str);
                let (negative, threads) = Coordinate::parse(threads)
                    .ok_or_else(|| invalid("argument --threads: invalid int value"))?;
                crate::builder_error::finish((|| {
                    require(
                        !negative && threads > Coordinate::zero(),
                        "threads must be positive",
                    )?;
                    primate::format(
                        &input,
                        &path_value(&values["cache-dir"]),
                        assembly,
                        &threads.0,
                        values.get("tmp-dir").map(|v| path_value(v)).as_deref(),
                    )
                })())?
            };
            if values.contains_key("print-path") {
                println!("{}", result.get("path").unwrap().as_str().unwrap());
            } else {
                println!("{}", json_text::dumps(&result, true).map_err(invalid)?);
            }
        }
    }
    Ok(())
}
