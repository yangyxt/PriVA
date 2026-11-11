//! Native `scripts/clinvar_gene_refresh.py`: the three-calendar-month ClinVar
//! gene-statistics refresh, the locked ClinVar annotation deployment with its
//! retained generations, and the unqualified installation-candidate receipts.
//!
//! The eleven subcommands of the script take its options, and `download` is
//! the install workflow's ClinVar download step (native only); the subcommand comes
//! first. Printed JSON, receipts, journals, retained members, configuration
//! rewrites and exit statuses follow the script. A failure that the script
//! reported as `ERROR: ClinVar gene-statistics refresh failed: message` exits 1
//! here as `priva: CATEGORY: message`, the category being the Python exception
//! class. The original scientific stages still run as
//! `bash scripts/install_utils.sh FUNCTION ...`, and `curl` replaces urllib
//! for the NCBI requests (same headers, Range, redirects and HTTP errors).
//!
//! The scientific profile binds the Python runtime the script ran in
//! (`sys.executable`, `sys.version` and the selected `sys.prefix/conda-meta`
//! records), so that its signature, and the refresh workspace it names, stay
//! those of the script and of `vep_chromosome_jobs.py`. That runtime is
//! observed once with `python -c` from PATH, as `alpha-annotation-base
//! motif-runtime` observes its environment; nothing else runs in Python.
//! Differences: temporary file and folder names carry six random characters.
mod deploy;
mod profile;
mod refresh;
mod remote;

use crate::alpha_annotation_base::io::{path as json_path, type_name};
use crate::builder_error::{self, RuntimeError};
use crate::cli::{Result, invalid};
use crate::installed_environment::yaml::{self, Value, Yaml};
use crate::installer_prerequisites::paths::{join, pathlib};
use crate::json_text::{self, DumpsOptions, Json};
use crate::python_args::parse_long_option_occurrences;
use crate::python_path::{os_error, repr, tempfile_parent};
use std::collections::HashMap;
use std::ffi::{CString, OsStr};
use std::fs::{self, File};
use std::io::{self, Read, Write};
use std::os::fd::FromRawFd;
use std::os::unix::ffi::OsStrExt;
use std::os::unix::fs::{MetadataExt, PermissionsExt};
use std::path::{Path, PathBuf};

const RECEIPT_KIND: &str = "priva_clinvar_gene_statistics_refresh_v1";
const CANDIDATE_KIND: &str = "priva_clinvar_annotation_candidate_v1";
const HELP: &str = "usage: priva clinvar-gene-refresh COMMAND [options]\n\
       {plan,install,deploy,retire-owned,retire-shards,check-header} --config CONFIG --source-root ROOT\n\
       plan [--output OUTPUT]; deploy|retire-owned [--source-vcf VCF]\n\
       retire-shards --input-vcf VCF --chromosome CHROM [--chromosome CHROM ...]\n\
       acquire --plan PLAN --output OUTPUT\n\
       download --config CONFIG --source-root ROOT --output FOLDER\n\
       file-date --input ANNOTATED_VCF\n\
       candidate-route --source-root ROOT\n\
       {candidate-begin,candidate-complete} --config CONFIG --source-root ROOT --input-vcf VCF\n\
Refresh ClinVar gene statistics after three calendar months of source-data age; the\n\
outputs and exit statuses of scripts/clinvar_gene_refresh.py. The subcommand comes first.";
const COMMANDS: [&str; 12] = [
    "plan",
    "install",
    "deploy",
    "retire-owned",
    "retire-shards",
    "check-header",
    "acquire",
    "download",
    "file-date",
    "candidate-route",
    "candidate-begin",
    "candidate-complete",
];

fn error(category: &'static str, message: impl Into<String>) -> Box<dyn std::error::Error> {
    builder_error::category(category, invalid(message))
}
fn value_error(message: impl Into<String>) -> Box<dyn std::error::Error> {
    error("ValueError", message)
}
fn require(condition: bool, message: impl Into<String>) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(value_error(message))
    }
}
/// The Python category and `str()` of a failure.
fn exception(failure: &(dyn std::error::Error + 'static)) -> (&'static str, String) {
    match failure.downcast_ref::<RuntimeError>() {
        Some(known) => (known.category_name(), known.message()),
        None => ("RuntimeError", failure.to_string()),
    }
}
fn object<const N: usize>(values: [(&str, Json); N]) -> Json {
    Json::Object(values.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}
fn string(value: impl Into<String>) -> Json {
    Json::Str(value.into())
}

// ---------------------------------------------------------------- Python values

/// `value[key]` on parsed JSON or YAML, with Python's KeyError and TypeError.
fn item<'a>(value: &'a Json, key: &str) -> Result<&'a Json> {
    match value {
        Json::Object(_) => value.get(key).ok_or_else(|| error("KeyError", repr(key))),
        Json::Array(_) => Err(error(
            "TypeError",
            "list indices must be integers or slices, not str",
        )),
        Json::Str(_) => Err(error(
            "TypeError",
            "string indices must be integers, not 'str'",
        )),
        other => Err(error(
            "TypeError",
            format!("'{}' object is not subscriptable", type_name(other)),
        )),
    }
}
/// `value.get(key)` on what must be a dict.
fn lookup<'a>(value: &'a Json, key: &str) -> Result<Option<&'a Json>> {
    match value {
        Json::Object(_) => Ok(value.get(key)),
        other => Err(error(
            "AttributeError",
            format!("'{}' object has no attribute 'get'", type_name(other)),
        )),
    }
}
/// `value.get(key)` with None for an absent key.
fn get_or_null(value: &Json, key: &str) -> Result<Json> {
    Ok(lookup(value, key)?.cloned().unwrap_or(Json::Null))
}
/// `value in (None, "", "null")`
fn empty(value: &Json) -> bool {
    match value {
        Json::Null => true,
        Json::Str(text) => text.is_empty() || text == "null",
        _ => false,
    }
}
/// Python truth of a parsed value.
fn truth(value: &Json) -> bool {
    crate::alpha_annotation_base::io::truth(value)
}
/// `str(value)`
fn py_str(value: &Json) -> String {
    match value {
        Json::Str(text) => text.clone(),
        other => py_repr(other),
    }
}
/// `repr(value)`
fn py_repr(value: &Json) -> String {
    match value {
        Json::Null => "None".into(),
        Json::Bool(flag) => if *flag { "True" } else { "False" }.into(),
        Json::Int(number) => number.to_string(),
        Json::Float(number) => json_text::python_float_repr(*number),
        Json::Str(text) => repr(text),
        Json::Array(items) => format!(
            "[{}]",
            items.iter().map(py_repr).collect::<Vec<_>>().join(", ")
        ),
        Json::Object(items) => format!(
            "{{{}}}",
            items
                .iter()
                .map(|(key, value)| format!("{}: {}", repr(key), py_repr(value)))
                .collect::<Vec<_>>()
                .join(", ")
        ),
    }
}
/// `str(Path(value))`, with pathlib's TypeError for a non-string.
fn path_of(value: &Json) -> Result<String> {
    Ok(json_path(value)?.display().to_string())
}
fn text_of(value: &Json) -> String {
    value.as_str().map_or_else(|| py_str(value), str::to_owned)
}

// ---------------------------------------------------------------- configuration

/// A YAML node as the dict `yaml.load` builds. Values JSON cannot represent
/// (timestamps, binary, sets) keep their tag and spelling for comparisons.
fn yaml_json(doc: &Yaml, node: usize, depth: usize) -> Json {
    let opaque = || {
        string(format!(
            "\u{0}{}:{}",
            doc.tag(node),
            doc.scalar_display(node)
        ))
    };
    if depth > 200 {
        return opaque();
    }
    match doc.value(node) {
        Value::Null => Json::Null,
        Value::Bool(flag) => Json::Bool(*flag),
        Value::Number(text) if doc.tag(node) == "int" => {
            text.parse().map(Json::Int).unwrap_or_else(|_| opaque())
        }
        Value::Number(text) => match text.to_ascii_lowercase().as_str() {
            ".inf" | "+.inf" => Json::Float(f64::INFINITY),
            "-.inf" => Json::Float(f64::NEG_INFINITY),
            ".nan" => Json::Float(f64::NAN),
            lower => lower.parse().map(Json::Float).unwrap_or_else(|_| opaque()),
        },
        Value::Str(text)
        | Value::Codepoints {
            transport: text, ..
        } => string(text.clone()),
        Value::Other(_) => opaque(),
        Value::Sequence(items) => Json::Array(
            items
                .iter()
                .map(|child| yaml_json(doc, *child, depth + 1))
                .collect(),
        ),
        Value::Mapping(pairs) => {
            let mut result = Json::Object(Vec::new());
            for (key, child) in pairs {
                let name = doc.string(*key).map_or_else(
                    || format!("\u{0}{}", doc.scalar_display(*key)),
                    str::to_owned,
                );
                result.set(&name, yaml_json(doc, *child, depth + 1));
            }
            result
        }
    }
}
/// `load_config(path)`: the whole document through installer_config_paths.UniqueLoader.
fn load_config(path: &str) -> Result<Json> {
    let text = read_text(path)?;
    let doc = yaml::load_unique(&text)?;
    require(
        doc.mapping(doc.root).is_some(),
        "The configuration must be a YAML mapping",
    )?;
    Ok(yaml_json(&doc, doc.root, 0))
}
/// `yaml.safe_load(text)` of one configuration value.
fn safe_load(text: &str) -> Result<Json> {
    let doc = yaml::load(text)?;
    Ok(yaml_json(&doc, doc.root, 0))
}
/// `{"hg19": "GRCh37", ...}.get(config.get("assembly"))`
fn assembly_of(value: &Json) -> Result<Option<&'static str>> {
    if matches!(value, Json::Array(_) | Json::Object(_)) {
        return Err(error(
            "TypeError",
            format!("unhashable type: '{}'", type_name(value)),
        ));
    }
    Ok(match value.as_str() {
        Some("hg19" | "GRCh37") => Some("GRCh37"),
        Some("hg38" | "GRCh38") => Some("GRCh38"),
        _ => None,
    })
}
fn ucsc(assembly: &str) -> &'static str {
    if assembly == "GRCh37" { "hg19" } else { "hg38" }
}
/// `re.fullmatch(r"[1-9][0-9]*", value)`
fn positive(value: &str) -> bool {
    value
        .bytes()
        .next()
        .is_some_and(|b| (b'1'..=b'9').contains(&b))
        && value.bytes().all(|b| b.is_ascii_digit())
}
/// `str(config.get("vep_installed_version", ""))`
fn vep_major(config: &Json) -> String {
    config
        .get("vep_installed_version")
        .map_or_else(String::new, py_str)
}

/// `selected_vcf(config, root, ucsc)`
fn selected_vcf(config: &Json, root: &str, ucsc: &str) -> Result<String> {
    let value = get_or_null(config, "clinvar_vcf")?;
    if empty(&value) {
        return Ok(join(root, &format!("clinvar.{ucsc}.vep.vcf.gz")));
    }
    let path = path_of(&value)?;
    Ok(if path.starts_with('/') || text_of(&value).contains('/') {
        path
    } else {
        join(root, &path)
    })
}

// ---------------------------------------------------------------- files

fn read_text(path: &str) -> Result<String> {
    crate::alpha_annotation_base::io::read_text(Path::new(path))
}
fn read_json(path: &str) -> Result<Json> {
    json_text::parse(&read_text(path)?, false).map_err(|message| {
        error(
            if json_text::is_refusal(&message) {
                "UnsupportedMetadata"
            } else {
                "JSONDecodeError"
            },
            message,
        )
    })
}
fn stat(path: &str) -> Result<fs::Metadata> {
    fs::metadata(path).map_err(|e| os_error(e, Path::new(path), None))
}
fn lstat(path: &str) -> Result<fs::Metadata> {
    fs::symlink_metadata(path).map_err(|e| os_error(e, Path::new(path), None))
}
fn nanoseconds(seconds: i64, nanos: i64) -> i128 {
    seconds as i128 * 1_000_000_000 + nanos as i128
}
/// `nonempty(path, label)`
fn nonempty(path: &str, label: &str) -> Result<fs::Metadata> {
    let info = stat(path)?;
    require(info.len() > 0, format!("Missing or empty {label}: {path}"))?;
    Ok(info)
}
fn sha256(path: &str) -> Result<String> {
    crate::sha256::sha256_file(Path::new(path)).map_err(|e| os_error(e, Path::new(path), None))
}
fn md5(path: &str) -> Result<String> {
    crate::md5::md5_file(Path::new(path)).map_err(|e| os_error(e, Path::new(path), None))
}
/// `identity(path, **extra)`
fn identity_with(path: &str, extra: Vec<(&str, Json)>) -> Result<Json> {
    let info = nonempty(path, "completed refresh artifact")?;
    let mut value = Json::Object(extra.into_iter().map(|(k, v)| (k.to_owned(), v)).collect());
    value.set("path", string(path));
    value.set("size", Json::Int(info.len() as i128));
    value.set("sha256", string(sha256(path)?));
    Ok(value)
}
fn identity(path: &str) -> Result<Json> {
    identity_with(path, Vec::new())
}
/// `os.path.samefile`-style comparison through `os.stat`.
fn same_file(left: &str, right: &str) -> Result<bool> {
    let (a, b) = (stat(left)?, stat(right)?);
    Ok(a.dev() == b.dev() && a.ino() == b.ino())
}
/// `tempfile.mkstemp(prefix=prefix, dir=parent)`: an empty 0600 file.
fn mkstemp(parent: &str, prefix: &str) -> Result<(File, String)> {
    let template = tempfile_parent(Path::new(parent))?.join(format!("{prefix}XXXXXX"));
    let mut name = CString::new(template.as_os_str().as_bytes())
        .map_err(|_| value_error("embedded null byte"))?
        .into_bytes_with_nul();
    let fd = unsafe { libc::mkostemp(name.as_mut_ptr().cast(), libc::O_CLOEXEC) };
    let path = PathBuf::from(OsStr::from_bytes(&name[..name.len() - 1]));
    if fd < 0 {
        return Err(os_error(io::Error::last_os_error(), &path, None));
    }
    Ok((unsafe { File::from_raw_fd(fd) }, path.display().to_string()))
}
/// `tempfile.mkdtemp(prefix=prefix, dir=parent)`
fn mkdtemp(parent: &str, prefix: &str) -> Result<String> {
    Ok(
        crate::alpha_annotation_base::io::stage(Path::new(parent), prefix)?
            .display()
            .to_string(),
    )
}
fn mkdirs(path: &str) -> Result<()> {
    fs::create_dir_all(path).map_err(|e| os_error(e, Path::new(path), None))
}
fn rename(source: &str, target: &str) -> Result<()> {
    fs::rename(source, target).map_err(|e| os_error(e, Path::new(source), Some(Path::new(target))))
}
fn unlink(path: &str) -> Result<()> {
    fs::remove_file(path).map_err(|e| os_error(e, Path::new(path), None))
}
fn parent_of(path: &str) -> String {
    crate::installer_prerequisites::paths::parent(path)
}
fn name_of(path: &str) -> &str {
    path.rsplit('/').next().unwrap_or("")
}
/// `atomic_text(path, value)`
fn atomic_text(path: &str, value: &str) -> Result<()> {
    let parent = parent_of(path);
    mkdirs(&parent)?;
    let (mut file, temporary) = mkstemp(&parent, &format!("{}.partial.", name_of(path)))?;
    let result = (|| {
        file.write_all(value.as_bytes())
            .map_err(|e| os_error(e, Path::new(&temporary), None))?;
        drop(file);
        if crate::installer_prerequisites::paths::exists(path)? {
            let mode = stat(path)?.mode() & 0o777;
            fs::set_permissions(&temporary, fs::Permissions::from_mode(mode))
                .map_err(|e| os_error(e, Path::new(&temporary), None))?;
        }
        rename(&temporary, path)
    })();
    if crate::installer_prerequisites::paths::lexists(&temporary) {
        unlink(&temporary)?;
    }
    result
}
fn pretty(value: &Json) -> Result<String> {
    json_text::dumps_with_options(
        value,
        DumpsOptions {
            indent: Some(2),
            sort_keys: true,
            ..DumpsOptions::default()
        },
    )
    .map_err(value_error)
}
/// `write_json(path, value)`
fn write_json(path: &str, value: &Json) -> Result<()> {
    atomic_text(path, &(pretty(value)? + "\n"))
}

// ---------------------------------------------------------------- time

fn civil_from_days(days: i128) -> (i128, u32, u32) {
    let days = days + 719468;
    let era = days.div_euclid(146097);
    let of_era = days - era * 146097;
    let of_cycle = (of_era - of_era / 1460 + of_era / 36524 - of_era / 146096) / 365;
    let of_year = of_era - (365 * of_cycle + of_cycle / 4 - of_cycle / 100);
    let shifted = (5 * of_year + 2) / 153;
    let day = (of_year - (153 * shifted + 2) / 5 + 1) as u32;
    let month = (if shifted < 10 {
        shifted + 3
    } else {
        shifted - 9
    }) as u32;
    let year = of_cycle + era * 400 + if month <= 2 { 1 } else { 0 };
    (year, month, day)
}
fn month_days(year: i128, month: u32) -> u32 {
    match month {
        2 if year % 4 == 0 && (year % 100 != 0 || year % 400 == 0) => 29,
        2 => 28,
        4 | 6 | 9 | 11 => 30,
        _ => 31,
    }
}
/// A `datetime.date` as (year, month, day).
type Date = (i128, u32, u32);
fn date_text((year, month, day): Date) -> String {
    format!("{year:04}-{month:02}-{day:02}")
}
/// `date.fromisoformat(text)` for the YYYY-MM-DD spellings this script writes.
fn parse_date(text: &str) -> Result<Date> {
    let bytes = text.as_bytes();
    let valid = bytes.len() == 10
        && bytes[4] == b'-'
        && bytes[7] == b'-'
        && bytes
            .iter()
            .enumerate()
            .all(|(i, b)| i == 4 || i == 7 || b.is_ascii_digit());
    require(valid, format!("Invalid isoformat string: {}", repr(text)))?;
    let year: i128 = text[..4].parse().unwrap_or_default();
    let month: u32 = text[5..7].parse().unwrap_or_default();
    let day: u32 = text[8..].parse().unwrap_or_default();
    date(year, month, day)
}
/// `datetime.date(year, month, day)`
fn date(year: i128, month: u32, day: u32) -> Result<Date> {
    require(
        (1..=9999).contains(&year),
        format!("year {year} is out of range"),
    )?;
    require((1..=12).contains(&month), "month must be in 1..12")?;
    require(
        (1..=month_days(year, month)).contains(&day),
        "day is out of range for month",
    )?;
    Ok((year, month, day))
}
/// `add_calendar_months(day, months)`: month ends are clamped.
fn add_calendar_months((year, month, day): Date, months: i128) -> Result<Date> {
    let ordinal = year * 12 + month as i128 - 1 + months;
    let (year, month) = (ordinal.div_euclid(12), ordinal.rem_euclid(12) as u32 + 1);
    date(year, month, day.min(month_days(year, month)))
}
/// `datetime.isoformat()` of a UTC instant in whole seconds and microseconds.
fn utc_text(seconds: i128, micros: u32) -> String {
    let (date, second) = (seconds.div_euclid(86400), seconds.rem_euclid(86400));
    let fraction = if micros == 0 {
        String::new()
    } else {
        format!(".{micros:06}")
    };
    format!(
        "{}T{:02}:{:02}:{:02}{fraction}+00:00",
        date_text(civil_from_days(date)),
        second / 3600,
        second / 60 % 60,
        second % 60
    )
}
/// `datetime.now(timezone.utc).isoformat()`
fn utc_now() -> String {
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default();
    utc_text(now.as_secs() as i128, now.subsec_micros())
}
/// `datetime.fromtimestamp(os.stat(path).st_mtime, timezone.utc)` as whole
/// seconds and microseconds: the float timestamp is rounded half to even.
fn float_timestamp(seconds: i64, nanos: i64) -> (i128, u32) {
    let value = seconds as f64 + nanos as f64 * 1e-9;
    let whole = value.trunc();
    let mut fraction = (value - whole) * 1e6;
    let rounded = fraction.round();
    fraction = if (fraction - rounded).abs() == 0.5 {
        2.0 * (fraction / 2.0).round()
    } else {
        rounded
    };
    let (mut whole, mut fraction) = (whole as i128, fraction);
    if fraction >= 1e6 {
        fraction -= 1e6;
        whole += 1;
    } else if fraction < 0.0 {
        fraction += 1e6;
        whole -= 1;
    }
    (whole, fraction as u32)
}

// ---------------------------------------------------------------- command line

fn options(command: &str) -> (Vec<&'static str>, Vec<&'static str>) {
    match command {
        "acquire" => (vec!["plan", "output"], vec!["plan", "output"]),
        "download" => (
            vec!["config", "source-root", "output"],
            vec!["config", "source-root", "output"],
        ),
        "file-date" => (vec!["input"], vec!["input"]),
        "candidate-route" => (vec!["source-root"], vec!["source-root"]),
        "candidate-begin" | "candidate-complete" => (
            vec!["config", "source-root", "input-vcf"],
            vec!["config", "source-root", "input-vcf"],
        ),
        _ => {
            let mut names = vec!["config", "source-root"];
            let required = names.clone();
            match command {
                "plan" => names.push("output"),
                "deploy" | "retire-owned" => names.push("source-vcf"),
                "retire-shards" => names.extend(["input-vcf", "chromosome"]),
                _ => {}
            }
            let mut required = required;
            if command == "retire-shards" {
                required.extend(["input-vcf", "chromosome"]);
            }
            (names, required)
        }
    }
}

fn run(command: &str, values: &HashMap<String, String>, chromosomes: Vec<String>) -> Result<()> {
    let path = |name: &str| values.get(name).map(|value| pathlib(value));
    let config = path("config").unwrap_or_default();
    let root = path("source-root").unwrap_or_default();
    let printed = |value: &Json| -> Result<()> {
        println!("{}", json_text::dumps(value, true).map_err(value_error)?);
        Ok(())
    };
    match command {
        "deploy" => deploy::deploy_cache(&config, &root, path("source-vcf").as_deref()),
        "retire-owned" => deploy::retire_owned_cache(
            &config,
            "existing deployment requires a cold annotation rebuild",
            path("source-vcf").as_deref(),
            None,
            &[],
        )
        .map(|_| ()),
        "retire-shards" => deploy::retire_rebuilt_shards(
            &config,
            &path("input-vcf").unwrap_or_default(),
            &chromosomes,
        )
        .map(|_| ()),
        "check-header" => printed(&deploy::deployment_header(
            &deploy::deployment_context(&config)?,
            None,
        )?),
        "candidate-route" => printed(&profile::candidate_route(&root)?),
        "candidate-begin" => printed(&profile::candidate_begin(
            &config,
            &root,
            &path("input-vcf").unwrap_or_default(),
        )?),
        "candidate-complete" => printed(&profile::candidate_complete(
            &config,
            &root,
            &path("input-vcf").unwrap_or_default(),
        )?),
        "file-date" => {
            println!(
                "{}",
                date_text(remote::local_file_date(&path("input").unwrap_or_default())?)
            );
            Ok(())
        }
        "acquire" => {
            let plan = read_json(&path("plan").unwrap_or_default())?;
            write_json(
                &path("output").unwrap_or_default(),
                &remote::acquire_inputs(&plan)?,
            )
        }
        "download" => printed(&remote::download_release(
            &config,
            &root,
            &path("output").unwrap_or_default(),
        )?),
        "plan" => {
            let value = refresh::plan_refresh(&config, &root)?;
            match path("output") {
                Some(output) => write_json(&output, &value),
                None => {
                    println!("{}", pretty(&value)?);
                    Ok(())
                }
            }
        }
        _ => refresh::install_refresh(&config, &root).map(|_| ()),
    }
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(command) = args.first() else {
        return Err(invalid(format!(
            "the following arguments are required: command\n{HELP}"
        )));
    };
    if matches!(command.as_str(), "-h" | "--help") {
        println!("{HELP}");
        return Ok(());
    }
    if !COMMANDS.contains(&command.as_str()) {
        return Err(invalid(format!(
            "argument command: invalid choice: {} (choose from {})",
            repr(command),
            COMMANDS.map(repr).join(", ")
        )));
    }
    let (names, required) = options(command);
    let Some(occurrences) =
        parse_long_option_occurrences(&args[1..], &names, &[]).map_err(invalid)?
    else {
        println!("{HELP}");
        return Ok(());
    };
    let chromosomes: Vec<String> = occurrences
        .iter()
        .filter(|(name, _)| name == "chromosome")
        .map(|(_, value)| value.clone())
        .collect();
    let values: HashMap<String, String> = occurrences.into_iter().collect();
    crate::python_args::require_options(&values, &required).map_err(invalid)?;
    let result = builder_error::finish(run(command, &values, chromosomes));
    io::stdout().flush().ok();
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Days since 1970-01-01 of a proleptic Gregorian date.
    fn days_from_civil(year: i128, month: u32, day: u32) -> i128 {
        let year = if month <= 2 { year - 1 } else { year };
        let era = year.div_euclid(400);
        let of_era = year - era * 400;
        let month = month as i128;
        let of_year =
            (153 * (if month > 2 { month - 3 } else { month + 9 }) + 2) / 5 + day as i128 - 1;
        let of_cycle = of_era * 365 + of_era / 4 - of_era / 100 + of_year;
        era * 146097 + of_cycle - 719468
    }

    #[test]
    fn calendar_months_clamp_month_ends() {
        assert_eq!(
            add_calendar_months((2026, 11, 30), 3).unwrap(),
            (2027, 2, 28)
        );
        assert_eq!(
            add_calendar_months((2023, 11, 30), 3).unwrap(),
            (2024, 2, 29)
        );
        assert_eq!(
            add_calendar_months((2026, 1, 31), 1).unwrap(),
            (2026, 2, 28)
        );
        assert!(add_calendar_months((9999, 12, 1), 1).is_err());
    }

    #[test]
    fn utc_text_omits_zero_microseconds() {
        assert_eq!(utc_text(0, 0), "1970-01-01T00:00:00+00:00");
        assert_eq!(
            utc_text(1_762_300_800, 5),
            "2025-11-05T00:00:00.000005+00:00"
        );
        assert_eq!(civil_from_days(days_from_civil(2024, 2, 29)), (2024, 2, 29));
    }

    #[test]
    fn float_timestamps_round_half_to_even() {
        // Python: datetime.fromtimestamp(1 + 500e-9, timezone.utc) is 1.000001 s.
        assert_eq!(float_timestamp(1, 500), (1, 1));
        assert_eq!(float_timestamp(1, 999_999_600), (2, 0));
        assert_eq!(
            float_timestamp(1_700_000_000, 123_456_789),
            (1_700_000_000, 123_457)
        );
    }

    #[test]
    fn python_strings_of_values() {
        assert_eq!(py_str(&Json::Int(116)), "116");
        assert_eq!(py_str(&Json::Null), "None");
        assert_eq!(
            py_repr(&Json::Array(vec![string("a"), Json::Bool(true)])),
            "['a', True]"
        );
    }
}
