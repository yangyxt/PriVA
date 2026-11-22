//! `priva install-task`: the runner of the install steps of
//! `conf/install_steps.json`, which run as Nextflow processes
//! (`workflows/install.nf`, `modules/local/install.nf` and
//! `modules/local/install_steps/`). It replaced `scripts/install_task.sh` with
//! the same commands, files and lineages:
//!
//! ```text
//! priva install-task plan    --registry R --config USER --out DIR --source-root CHECKOUT
//!                            [--adopt-config OLD --adopt all|ID,ID] [--run ID,ID]
//!                            [--pins PINS.json] [--refresh ID,ID] [--releases FILE]...
//! priva install-task preflight --registry R --plan PLAN --source-root CHECKOUT [--stub]
//! priva install-task record  --registry R --plan PLAN --step ID --lineage L --folder F
//!                            --receipt FILE --source-root CHECKOUT [--input RECEIPT]...
//!                            [--key KEY=PATH]... [--stub]
//! priva install-task adopt   --plan PLAN --step ID --lineage L --receipt FILE --source-root CHECKOUT
//! priva install-task publish --registry R --plan PLAN --lineage L --folder F --source-root CHECKOUT
//!                            [--input RECEIPT]... [--pins PINS.json] [--stub]
//! priva install-task release --registry R --config USER --step ID --out FILE [--stub]
//! priva install-task reach   --registry R --step ID --source-root CHECKOUT
//! priva install-task value   --config FILE --key KEY
//! ```
//!
//! release: asks a changing source for its current release (registry
//! `release`: git commit, or the server's ETag and Last-Modified) and writes
//! `{"step": ID, "release": TEXT, "sources": [...]}`. `plan --releases FILE...`
//! puts it in the step's plan entry; a pinned step uses `resource_root/pins.json`
//! instead (`--pins`), unless it is listed in `--refresh`. `publish --pins`
//! records the releases of refreshed pinned steps there.
//!
//! plan: writes DIR/base.yaml (your configuration without any key an install
//! step writes, so each step resolves its own destinations) and DIR/plan.json
//! (per step: code digest, settings, the files you supply, adopted keys).
//!
//! preflight: checks the files you supply for the steps that will run and
//! the data files shipped in the checkout, before any step (see `preflight`).
//!
//! record: the receipt of a step (registry `module`): its process writes into
//! the task's `data/` with the tools' own commands, and this checks the output
//! keys and names them where storeDir moves them, F/data (see `record`).
//! Nextflow (storeDir) skips a step whose F holds a receipt.
//!
//! adopt: writes a receipt naming the files of an existing installation.
//!
//! reach: lists the module, script files, shipped data files, `priva`
//! commands and Rust files a step's code digest covers, and the digest.
//!
//! value: prints one configuration value as the installer reads it (`yq -r
//! '.[KEY] // empty | tostring'`).
//!
//! The plan's settings, file identities and release texts are the values the
//! shell runner wrote; the digests are those of `digest`.
mod content;
mod digest;
mod plan;
mod preflight;
mod publish;
mod record;

use crate::cli::{Result, invalid};
use crate::json_text::{self, DumpsOptions, Json};
use crate::yaml_config::{self, ConfigValue, Parsed, Scalar};
use std::collections::BTreeMap;
use std::ffi::OsStr;
use std::io::Write as _;
use std::os::unix::fs::MetadataExt;
use std::path::{Component, Path, PathBuf};
use std::process::{Command, Stdio};

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(command) = args.first() else {
        return Err(invalid(USAGE));
    };
    let options = Options::parse(&args[1..])?;
    match command.as_str() {
        "plan" => plan::plan(&options),
        "preflight" => preflight::preflight(&options),
        "release" => plan::release(&options),
        "reach" => digest::reach(&options),
        "record" => record::record(&options),
        "adopt" => publish::adopt(&options),
        "publish" => publish::publish(&options),
        "value" => {
            let config = options.path("config")?;
            let key = options.one("key")?;
            println!("{}", yaml_value(&config, key)?);
            Ok(())
        }
        "-h" | "--help" => {
            println!("{USAGE}");
            Ok(())
        }
        other => Err(invalid(format!(
            "install-task: unknown command {other}\n{USAGE}"
        ))),
    }
}

const USAGE: &str = "usage: priva install-task plan|preflight|record|adopt|publish|release|reach|value [options]
  plan    --registry R --config USER --out DIR --source-root CHECKOUT [--adopt-config OLD --adopt all|ID,ID]
          [--run ID,ID] [--pins PINS.json] [--refresh ID,ID] [--releases FILE]...
  preflight --registry R --plan PLAN --source-root CHECKOUT [--stub]
  record  --registry R --plan PLAN --step ID --lineage L --folder F --receipt FILE --source-root CHECKOUT
          [--input RECEIPT]... [--key KEY=PATH]... [--stub]
  adopt   --plan PLAN --step ID --lineage L --receipt FILE --source-root CHECKOUT
  publish --registry R --plan PLAN --lineage L --folder F --source-root CHECKOUT [--input RECEIPT]...
          [--pins PINS.json] [--stub]
  release --registry R --config USER --step ID --out FILE [--stub]
  reach   --registry R --step ID --source-root CHECKOUT
  value   --config FILE --key KEY";

/// Options as `--name VALUE` pairs (repeatable ones keep every value) and
/// `--stub`.
pub(crate) struct Options {
    values: BTreeMap<String, Vec<String>>,
    pub(crate) stub: bool,
}

impl Options {
    fn parse(args: &[String]) -> Result<Options> {
        let mut values: BTreeMap<String, Vec<String>> = BTreeMap::new();
        let mut stub = false;
        let mut at = 0;
        while at < args.len() {
            let arg = &args[at];
            if arg == "--stub" {
                stub = true;
                at += 1;
                continue;
            }
            let Some(name) = arg.strip_prefix("--") else {
                return Err(invalid(format!("install-task: unknown option {arg}")));
            };
            let Some(value) = args.get(at + 1) else {
                return Err(invalid(format!("install-task: {arg} needs a value")));
            };
            values
                .entry(name.to_owned())
                .or_default()
                .push(value.clone());
            at += 2;
        }
        Ok(Options { values, stub })
    }

    pub(crate) fn get(&self, name: &str) -> Option<&str> {
        self.values
            .get(name)
            .and_then(|v| v.last())
            .map(String::as_str)
    }

    pub(crate) fn one(&self, name: &str) -> Result<&str> {
        self.get(name)
            .filter(|value| !value.is_empty())
            .ok_or_else(|| invalid(format!("install-task: --{name} is required")))
    }

    pub(crate) fn path(&self, name: &str) -> Result<PathBuf> {
        self.one(name).map(PathBuf::from)
    }

    pub(crate) fn all(&self, name: &str) -> &[String] {
        self.values.get(name).map(Vec::as_slice).unwrap_or(&[])
    }

    /// The checkout the steps run from: `conf/install_steps.json` is the
    /// step registry, `modules` their processes, `data` the shipped files.
    pub(crate) fn source_root(&self) -> Result<PathBuf> {
        let root = self.path("source-root")?;
        if !root.join("conf").join("install_steps.json").is_file() {
            return Err(invalid(format!(
                "install-task: {} is not a PriVA checkout (no conf/install_steps.json)",
                root.display()
            )));
        }
        Ok(absolute(&root))
    }
}

pub(crate) fn fail(message: impl Into<String>) -> Box<dyn std::error::Error> {
    invalid(format!("install_task: {}", message.into()))
}

// ---------------------------------------------------------------------------
// Registry

pub(crate) fn read_json(path: &Path) -> Result<Json> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| fail(format!("cannot read {}: {e}", path.display())))?;
    json_text::parse(&text, false).map_err(|e| fail(format!("{}: {e}", path.display())))
}

pub(crate) fn steps(registry: &Json) -> Result<&Vec<Json>> {
    match registry.get("steps") {
        Some(Json::Array(steps)) => Ok(steps),
        _ => Err(fail("the registry has no steps")),
    }
}

pub(crate) fn step<'a>(registry: &'a Json, id: &str) -> Result<&'a Json> {
    steps(registry)?
        .iter()
        .find(|s| s.get("id").and_then(Json::as_str) == Some(id))
        .ok_or_else(|| fail(format!("no step {id} in the registry")))
}

/// The strings of a list field (`.field[]?`); other items are skipped.
pub(crate) fn strings(entry: &Json, field: &str) -> Vec<String> {
    match entry.get(field) {
        Some(Json::Array(items)) => items
            .iter()
            .filter_map(|item| item.as_str().map(str::to_owned))
            .collect(),
        _ => Vec::new(),
    }
}

pub(crate) fn text(entry: &Json, field: &str) -> String {
    entry
        .get(field)
        .and_then(Json::as_str)
        .unwrap_or("")
        .to_owned()
}

pub(crate) fn truthy(value: Option<&Json>) -> bool {
    !matches!(value, None | Some(Json::Null) | Some(Json::Bool(false)))
}

pub(crate) fn object(items: Vec<(String, Json)>) -> Json {
    Json::Object(items)
}

pub(crate) fn string(value: impl Into<String>) -> Json {
    Json::Str(value.into())
}

/// jq's output: two-space indentation, Unicode as written.
pub(crate) fn pretty(value: &Json) -> Result<String> {
    let options = DumpsOptions {
        indent: Some(2),
        ensure_ascii: false,
        ..DumpsOptions::default()
    };
    Ok(json_text::dumps_with_options(value, options).map_err(invalid)? + "\n")
}

/// jq -c's output.
pub(crate) fn compact(value: &Json) -> Result<String> {
    let options = DumpsOptions {
        compact: true,
        ensure_ascii: false,
        ..DumpsOptions::default()
    };
    json_text::dumps_with_options(value, options).map_err(invalid)
}

pub(crate) fn write_file(path: &Path, content: &str) -> Result<()> {
    std::fs::write(path, content).map_err(|e| fail(format!("cannot write {}: {e}", path.display())))
}

// ---------------------------------------------------------------------------
// Configuration values

/// A configuration value as the installer reads it: `yq -r '.[KEY] // empty |
/// tostring'`. A missing key, null and false are empty; true is `true`; a list
/// or mapping is its compact JSON.
pub(crate) fn yaml_value(path: &Path, key: &str) -> Result<String> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| fail(format!("cannot read {}: {e}", path.display())))?;
    let document = match yaml_config::parse_config(&content)
        .map_err(|e| fail(format!("{}: {e}", path.display())))?
    {
        Parsed::Mapping(document) => document,
        Parsed::Empty => return Ok(String::new()),
        Parsed::List => return Err(fail(format!("{} is a list, not a mapping", path.display()))),
    };
    // PyYAML keeps the last of repeated keys.
    let Some(entry) = document.entries.iter().rev().find(|e| e.key == key) else {
        return Ok(String::new());
    };
    match &entry.value {
        ConfigValue::Scalar(scalar) => scalar_text(scalar, path, key),
        ConfigValue::Opaque(text) => opaque_json(text).ok_or_else(|| {
            fail(format!(
                "{}: {key} holds a value this reader does not take",
                path.display()
            ))
        }),
    }
}

/// Every top-level key, sorted and once each (`yq -r 'keys[]'`).
pub(crate) fn yaml_keys(path: &Path) -> Result<Vec<String>> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| fail(format!("cannot read {}: {e}", path.display())))?;
    let mut keys: Vec<String> = match yaml_config::parse_config(&content)
        .map_err(|e| fail(format!("{}: {e}", path.display())))?
    {
        Parsed::Mapping(document) => document.entries.into_iter().map(|e| e.key).collect(),
        _ => Vec::new(),
    };
    keys.sort();
    keys.dedup();
    Ok(keys)
}

fn scalar_text(scalar: &Scalar, path: &Path, key: &str) -> Result<String> {
    Ok(match scalar {
        Scalar::Null | Scalar::Bool(false) => String::new(),
        Scalar::Bool(true) => "true".to_owned(),
        Scalar::Int(text) | Scalar::Str(text) => text.clone(),
        Scalar::Other(_) => match scalar.python_type_name() {
            "float" | "int" => scalar.python_repr(),
            _ => {
                return Err(fail(format!(
                    "{}: {key} holds a {} value; write it in quotes",
                    path.display(),
                    scalar.python_type_name()
                )));
            }
        },
    })
}

fn scalar_json(scalar: &Scalar) -> Option<String> {
    Some(match scalar {
        Scalar::Null => "null".to_owned(),
        Scalar::Bool(flag) => flag.to_string(),
        Scalar::Int(text) => text.clone(),
        Scalar::Str(text) => json_text::quote(text, false),
        Scalar::Other(_) => match scalar.python_type_name() {
            "float" | "int" => scalar.python_repr(),
            _ => return None,
        },
    })
}

/// A flow list or a nested mapping as compact JSON.
fn opaque_json(text: &str) -> Option<String> {
    if let Some(inner) = text
        .strip_prefix('[')
        .and_then(|rest| rest.strip_suffix(']'))
    {
        let mut items: Vec<&str> = inner.split(',').collect();
        if items.len() > 1 && items[items.len() - 1].trim().is_empty() {
            items.pop();
        }
        if items.len() == 1 && items[0].trim().is_empty() {
            return Some("[]".to_owned());
        }
        let parts: Option<Vec<String>> = items
            .iter()
            .map(|item| scalar_json(&yaml_config::resolve_plain(item.trim())))
            .collect();
        return Some(format!("[{}]", parts?.join(",")));
    }
    let entries = yaml_config::nested_mapping_entries(text)?;
    let mut pairs: Vec<(String, String)> = Vec::new();
    for (key, value) in entries {
        let value = scalar_json(&value)?;
        match pairs.iter_mut().find(|(name, _)| *name == key) {
            Some(slot) => slot.1 = value,
            None => pairs.push((key, value)),
        }
    }
    let parts: Vec<String> = pairs
        .iter()
        .map(|(k, v)| format!("{}:{v}", json_text::quote(k, false)))
        .collect();
    Some(format!("{{{}}}", parts.join(",")))
}

fn is_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// Whether a line is `KEY[[:space:]]*:...` (sed's `/^KEY[[:space:]]*:/`).
fn line_has_key(line: &[u8], key: &[u8]) -> bool {
    if !line.starts_with(key) {
        return false;
    }
    let rest = &line[key.len()..];
    let spaces = rest.iter().take_while(|b| is_space(**b)).count();
    rest.get(spaces) == Some(&b':')
}

fn lines(content: &[u8]) -> Vec<&[u8]> {
    let mut out: Vec<&[u8]> = content.split(|b| *b == b'\n').collect();
    // A final newline ends the last line instead of starting another.
    if content.ends_with(b"\n") {
        out.pop();
    }
    out
}

fn join_lines(lines: &[Vec<u8>], final_newline: bool) -> Vec<u8> {
    let mut out = Vec::new();
    for (i, line) in lines.iter().enumerate() {
        out.extend_from_slice(line);
        if i + 1 < lines.len() || final_newline {
            out.push(b'\n');
        }
    }
    out
}

fn read_bytes(path: &Path) -> Result<Vec<u8>> {
    std::fs::read(path).map_err(|e| fail(format!("cannot read {}: {e}", path.display())))
}

fn write_bytes(path: &Path, content: &[u8]) -> Result<()> {
    std::fs::write(path, content).map_err(|e| fail(format!("cannot write {}: {e}", path.display())))
}

/// Replaces the value of every `KEY:` line, or appends `KEY: VALUE` (the
/// line's comment goes with the old value).
pub(crate) fn set_key(path: &Path, key: &str, value: &str) -> Result<()> {
    let content = read_bytes(path)?;
    let key_bytes = key.as_bytes();
    if lines(&content)
        .iter()
        .any(|line| line_has_key(line, key_bytes))
    {
        let final_newline = content.ends_with(b"\n");
        let replaced: Vec<Vec<u8>> = lines(&content)
            .into_iter()
            .map(|line| {
                if line_has_key(line, key_bytes) {
                    let colon = line.iter().position(|b| *b == b':').unwrap_or(line.len());
                    let mut out = line[..=colon].to_vec();
                    out.push(b' ');
                    out.extend_from_slice(value.as_bytes());
                    out
                } else {
                    line.to_vec()
                }
            })
            .collect();
        write_bytes(path, &join_lines(&replaced, final_newline))
    } else {
        let mut file = std::fs::OpenOptions::new()
            .append(true)
            .open(path)
            .map_err(|e| fail(format!("cannot write {}: {e}", path.display())))?;
        writeln!(file, "{key}: {value}")
            .map_err(|e| fail(format!("cannot write {}: {e}", path.display())))
    }
}

/// Removes every `KEY:` line.
pub(crate) fn delete_key(path: &Path, key: &str) -> Result<()> {
    let content = read_bytes(path)?;
    let final_newline = content.ends_with(b"\n");
    let all = lines(&content);
    let kept: Vec<Vec<u8>> = all
        .iter()
        .filter(|line| !line_has_key(line, key.as_bytes()))
        .map(|line| line.to_vec())
        .collect();
    if kept.len() == all.len() {
        return Ok(());
    }
    // sed leaves the newline off only when it prints the unterminated last line.
    let last_kept = all
        .last()
        .is_some_and(|line| !line_has_key(line, key.as_bytes()));
    write_bytes(path, &join_lines(&kept, final_newline || !last_kept))
}

// ---------------------------------------------------------------------------
// Files

/// `path size mtime` for a file, `path folder mtime` for a folder (a folder's
/// size is not stable on StorNext: 0, then 1, for an unchanged LOFTEE
/// repository between two launches of job 1540691), or `missing path`.
pub(crate) fn identity(path: &str) -> String {
    match std::fs::metadata(path) {
        Ok(meta) if meta.is_dir() => format!("{path} folder {}", meta.mtime()),
        Ok(meta) => format!("{path} {} {}", meta.size(), meta.mtime()),
        Err(_) => format!("missing {path}"),
    }
}

/// The path made absolute against the working directory, `.` and `..`
/// resolved by name (`realpath -ms`).
pub(crate) fn absolute(path: &Path) -> PathBuf {
    let joined = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir().unwrap_or_default().join(path)
    };
    let mut out = PathBuf::from("/");
    for part in joined.components() {
        match part {
            Component::ParentDir => {
                out.pop();
            }
            Component::Normal(name) => out.push(name),
            _ => {}
        }
    }
    out
}

/// `realpath -m`: symbolic links resolved where the path exists, the rest
/// taken by name.
pub(crate) fn realpath_m(path: &Path) -> PathBuf {
    fn walk(path: &Path, depth: usize) -> PathBuf {
        let joined = if path.is_absolute() {
            path.to_path_buf()
        } else {
            std::env::current_dir().unwrap_or_default().join(path)
        };
        let mut out = PathBuf::from("/");
        for part in joined.components() {
            match part {
                Component::ParentDir => {
                    out.pop();
                }
                Component::Normal(name) => {
                    let next = out.join(name);
                    match std::fs::symlink_metadata(&next) {
                        Ok(meta) if meta.file_type().is_symlink() && depth < 40 => {
                            let target = std::fs::read_link(&next).unwrap_or_default();
                            out = walk(&out.join(target), depth + 1);
                        }
                        _ => out = next,
                    }
                }
                _ => {}
            }
        }
        out
    }
    walk(path, 0)
}

/// Removes a folder a failed attempt left (`chmod -R u+w; rm -rf`).
pub(crate) fn clear(path: &Path) -> Result<()> {
    if std::fs::symlink_metadata(path).is_err() {
        return Ok(());
    }
    set_writable(path, true);
    let removed = if std::fs::symlink_metadata(path).is_ok_and(|m| m.is_dir()) {
        std::fs::remove_dir_all(path)
    } else {
        std::fs::remove_file(path)
    };
    removed.map_err(|e| fail(format!("cannot clear {}: {e}", path.display())))
}

/// `chmod -R u+w` (true) or `chmod -R a-w` (false); symbolic links inside
/// are left alone, as chmod -R does.
pub(crate) fn set_writable(path: &Path, writable: bool) {
    use std::os::unix::fs::PermissionsExt;
    fn apply(path: &Path, writable: bool, top: bool) {
        let meta = if top {
            std::fs::metadata(path)
        } else {
            std::fs::symlink_metadata(path)
        };
        let Ok(meta) = meta else { return };
        if meta.file_type().is_symlink() {
            return;
        }
        let mode = meta.permissions().mode();
        let mode = if writable {
            mode | 0o200
        } else {
            mode & !0o222
        };
        let _ = std::fs::set_permissions(path, std::fs::Permissions::from_mode(mode));
        if meta.is_dir() {
            if let Ok(entries) = std::fs::read_dir(path) {
                for entry in entries.flatten() {
                    apply(&entry.path(), writable, false);
                }
            }
        }
    }
    apply(path, writable, true);
}

// ---------------------------------------------------------------------------
// Time and host

fn civil(days: i64) -> (i64, u32, u32) {
    // Howard Hinnant's days-to-civil.
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    (if m <= 2 { y + 1 } else { y }, m, d)
}

fn now_seconds() -> i64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0)
}

/// Today in UTC, `date -u +%F`.
pub(crate) fn today_utc() -> String {
    let (y, m, d) = civil(now_seconds().div_euclid(86_400));
    format!("{y:04}-{m:02}-{d:02}")
}

/// Local time as `date -Is` prints it: `2026-11-17T18:43:21+00:00`.
pub(crate) fn now_iso() -> String {
    let seconds = now_seconds();
    // SAFETY: localtime_r fills the zeroed tm from a valid time_t.
    let mut tm: libc::tm = unsafe { std::mem::zeroed() };
    let t: libc::time_t = seconds as libc::time_t;
    unsafe { libc::localtime_r(&t, &mut tm) };
    let offset = tm.tm_gmtoff;
    let sign = if offset < 0 { '-' } else { '+' };
    let offset = offset.abs();
    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}{sign}{:02}:{:02}",
        tm.tm_year + 1900,
        tm.tm_mon + 1,
        tm.tm_mday,
        tm.tm_hour,
        tm.tm_min,
        tm.tm_sec,
        offset / 3600,
        (offset % 3600) / 60
    )
}

pub(crate) fn hostname() -> String {
    let mut buffer = [0u8; 256];
    // SAFETY: the buffer is writable for its length.
    let rc = unsafe { libc::gethostname(buffer.as_mut_ptr().cast(), buffer.len()) };
    if rc != 0 {
        return String::new();
    }
    let end = buffer.iter().position(|b| *b == 0).unwrap_or(buffer.len());
    String::from_utf8_lossy(&buffer[..end]).into_owned()
}

// ---------------------------------------------------------------------------
// Programs

/// Runs a program and returns its standard output; `None` when it fails.
pub(crate) fn output_of<I, S>(program: &str, args: I) -> Option<Vec<u8>>
where
    I: IntoIterator<Item = S>,
    S: AsRef<OsStr>,
{
    let output = Command::new(program)
        .args(args)
        .stderr(Stdio::null())
        .output()
        .ok()?;
    output.status.success().then_some(output.stdout)
}

#[cfg(test)]
mod tests;
