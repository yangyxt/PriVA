//! Native `scripts/build_gene_nonlof_mechanism_cache.py`: the non-LOF
//! pathogenic-mechanism cache builder with the script's options.
//!
//! Source freshness and receipts, the G2P, ClinGen dosage, PanelApp and
//! GoFCards readers (pandas' string semantics, including the GoFCards
//! workbook that pandas reads through openpyxl), the HGNC mapping, the
//! canonical JSON and its Draft 2020-12 schema check, the build lock, the
//! weekly ClinVar VCV acquisition and `--fetch-sources-only` follow the
//! script. Printed JSON, receipts, the manifest TSV and exit statuses follow
//! the script; an exception the script let escape exits 1 here as
//! `priva: CATEGORY: message`.
//!
//! Differences:
//! - `--source-root CHECKOUT` is required. The defaults the script derives
//!   from its own location come from it, and `build_inputs` keeps the
//!   SHA-256 of the Python builder and `clinvar_vcv.py` from that checkout, so
//!   a cache built by either one is reused by the other.
//! - urllib transfers run through `curl --http1.1` with urllib's headers and
//!   retries; an HTTP status outside 2xx is `HTTP Error CODE: REASON`, other
//!   transfer failures are reported with curl's reason.
//! - The ClinVar XML is read with the same Expat parser ElementTree uses.
//! - Schema patterns are matched by a small regular-expression reader that
//!   covers anchors, classes, `.`, `\d`, `\w`, `\s` and the `*`, `+`, `?`
//!   quantifiers; any other pattern syntax is refused.
mod build;
mod fetch;
mod parse;
mod schema;
mod xlsx;

pub(crate) use build::{eligible, load_gofcards_cache};

use crate::builder_error;
use crate::cli::{Result, invalid};
use crate::json_text::{self, DumpsOptions, Json};
use crate::python_args::{parse_long_options_converted, require_options};
use crate::python_path::{os_error, path_value};
use std::collections::HashMap;
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};
use std::os::unix::fs::MetadataExt;
use std::path::{Path, PathBuf};

const HELP: &str = "usage: priva gene-nonlof-mechanism --source-root CHECKOUT [--cache-dir CACHE_DIR]\n\
       [--output-json OUTPUT_JSON] [--shared-raw-dir SHARED_RAW_DIR]\n\
       [--gofcards-exact-variants GOFCARDS_EXACT_VARIANTS] [--hgnc-table HGNC_TABLE]\n\
       [--output-schema OUTPUT_SCHEMA] [--clinvar-vcv CLINVAR_VCV]\n\
       [--clinvar-min-review-stars {2,3,4}] [--clinvar-max-download-seconds SECONDS]\n\
       [--clinvar-only-refresh] [--fetch-sources-only] [--force]\n\
       [--max-panelapp-panels N] [--timeout TIMEOUT] [--retries RETRIES]\n\
       [--proxy-url PROXY_URL] [--download-tool {auto,urllib,curl}]\n\
       [--validate-only JSON] [--stale-lock-hours HOURS] [--g2p-ddg2p-url URL]\n\
       [--clingen-dosage-grch38-url URL] [--panelapp-panels-url URL]\n\
Build PriVA's portable non-LOF pathogenic-mechanism cache; the options, outputs and\n\
exit statuses of scripts/build_gene_nonlof_mechanism_cache.py.";

const VALUES: &[&str] = &[
    "source-root",
    "cache-dir",
    "output-json",
    "shared-raw-dir",
    "gofcards-exact-variants",
    "hgnc-table",
    "output-schema",
    "clinvar-vcv",
    "clinvar-min-review-stars",
    "clinvar-max-download-seconds",
    "max-panelapp-panels",
    "timeout",
    "retries",
    "proxy-url",
    "download-tool",
    "validate-only",
    "stale-lock-hours",
    "g2p-ddg2p-url",
    "clingen-dosage-grch38-url",
    "panelapp-panels-url",
];
const FLAGS: &[&str] = &["clinvar-only-refresh", "fetch-sources-only", "force"];

const NONLOF_ASSERTIONS_FILENAME: &str = "gene_nonlof_mechanism_curated_assertions.json.gz";
const NONLOF_ASSERTIONS_SCHEMA_FILENAME: &str =
    "gene_nonlof_mechanism_curated_assertions.schema.json";
const NONLOF_SOURCE_MANIFEST_FILENAME: &str = "nonlof_source_manifest.json";
const NONLOF_SOURCE_MANIFEST_TSV_FILENAME: &str = "nonlof_source_manifest.tsv";
const NONLOF_RUN_SUMMARY_FILENAME: &str = "nonlof_run_summary.json";
const USER_AGENT: &str = "PriVA/0.1 (portable non-LOF pathogenic-mechanism cache)";

/// A source as the script's `SourceSpec` records it.
#[derive(Clone, Debug)]
pub(super) struct Spec {
    pub(super) name: &'static str,
    pub(super) url: String,
    pub(super) raw_filename: &'static str,
    pub(super) refresh_days: i64,
    pub(super) source_role: &'static str,
    pub(super) official_update_note: &'static str,
    pub(super) docs_url: &'static str,
    pub(super) parser: &'static str,
}

const G2P_URL: &str = "https://www.ebi.ac.uk/gene2phenotype/api/panel/all/download";
const CLINGEN_URL: &str = "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv";
const GOFCARDS_URL: &str =
    "https://download.genemed.tech/upload/GainFunCards/gofcards_data_download.xlsx";
const PANELAPP_URL: &str = "https://panelapp.genomicsengland.co.uk/api/v1/panels/";
const CLINVAR_URL: &str = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/weekly_release/\
ClinVarVCVRelease_00-latest_weekly.xml.gz";

impl Spec {
    pub(super) fn default_url(&self) -> &'static str {
        match self.name {
            "g2p_ddg2p" => G2P_URL,
            "clingen_dosage_gene_grch38" => CLINGEN_URL,
            "gofcards_gof_variants" => GOFCARDS_URL,
            "panelapp_all_panels" => PANELAPP_URL,
            _ => CLINVAR_URL,
        }
    }
    /// The script's `custom_source`: a selected URL differs from the default.
    pub(super) fn custom(&self) -> bool {
        self.url != self.default_url()
    }
}

pub(super) fn default_specs() -> Vec<Spec> {
    vec![
        Spec {
            name: "g2p_ddg2p",
            url: G2P_URL.into(),
            raw_filename: "AllG2P.csv",
            refresh_days: 30,
            source_role: "gene-disease molecular mechanism, allelic requirement, confidence, PMIDs",
            official_update_note: "G2P was redesigned as a Vue SPA in 2025. The old DDG2P.csv.gz \
static download URL is gone. Use the new API endpoint: /api/panel/{panel}/download (CSV, not \
gzipped). AllG2P covers all panels including DD, Eye, Skin, Cancer, Cardiac.",
            docs_url: "https://www.ebi.ac.uk/gene2phenotype/api/",
            parser: "g2p_ddg2p",
        },
        Spec {
            name: "clingen_dosage_gene_grch38",
            url: CLINGEN_URL.into(),
            raw_filename: "ClinGen_gene_curation_list_GRCh38.tsv",
            refresh_days: 1,
            source_role: "haploinsufficiency and triplosensitivity gene curation",
            official_update_note: "ClinGen says dosage curated gene TSV files are updated daily; \
the FTP index is refreshed nightly.",
            docs_url: "https://search.clinicalgenome.org/kb/downloads",
            parser: "clingen_dosage",
        },
        Spec {
            name: "gofcards_gof_variants",
            url: GOFCARDS_URL.into(),
            raw_filename: "gofcards/gofcards_data_download.xlsx",
            refresh_days: 180,
            source_role: "curated gain-of-function variant database (3161 variants, 579 genes)",
            official_update_note: "GoFCards (NAR 2025, 53:D976-D988, PMID:39578693) is a curated \
GOF variant database. Updates follow publication cycles.",
            docs_url: "https://academic.oup.com/nar/article/53/D1/D976/7907365",
            parser: "gofcards",
        },
    ]
}

pub(super) fn panelapp_spec() -> Spec {
    Spec {
        name: "panelapp_all_panels",
        url: PANELAPP_URL.into(),
        raw_filename: "panelapp_all_panels.json",
        refresh_days: 7,
        source_role: "PanelApp mode of inheritance and mode of pathogenicity",
        official_update_note: "PanelApp is curator/API-backed and does not publish one fixed bulk \
release cycle; panel versions expose version_created timestamps.",
        docs_url: "https://panelapp.genomicsengland.co.uk/api/v1/panels/",
        parser: "panelapp",
    }
}

pub(super) fn clinvar_spec() -> Spec {
    Spec {
        name: "clinvar_vcv_weekly",
        url: CLINVAR_URL.into(),
        raw_filename: "clinvar_vcv/ClinVarVCVRelease_00-latest_weekly.xml.gz",
        refresh_days: 7,
        source_role: "condition-specific germline classifications and submitted observation \
evidence for exact normalized GoFCards GOF alleles",
        official_update_note: "ClinVar generates VCV XML weekly. Check the companion MD5 first \
and download the multi-gigabyte XML only when that checksum changes.",
        docs_url: "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/_README",
        parser: "clinvar_vcv",
    }
}

/// The download settings every acquisition shares.
pub(super) struct Transfer {
    pub(super) timeout: i128,
    pub(super) retries: i128,
    pub(super) proxy_url: String,
    pub(super) download_tool: String,
    pub(super) max_download_seconds: i128,
    pub(super) max_panels: Option<i128>,
}

pub(crate) fn error(
    category: &'static str,
    message: impl Into<String>,
) -> Box<dyn std::error::Error> {
    builder_error::category(category, invalid(message))
}

/// `str(error)` without the category prefix.
pub(super) fn message(failure: &(dyn std::error::Error + 'static)) -> String {
    failure
        .downcast_ref::<builder_error::RuntimeError>()
        .map_or_else(|| failure.to_string(), builder_error::RuntimeError::message)
}

pub(crate) fn object<const N: usize>(values: [(&str, Json); N]) -> Json {
    Json::Object(values.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}
pub(crate) fn text(value: impl Into<String>) -> Json {
    Json::Str(value.into())
}

/// Python's `type(value).__name__` for a parsed JSON value.
pub(crate) fn type_name(value: &Json) -> &'static str {
    match value {
        Json::Null => "NoneType",
        Json::Bool(_) => "bool",
        Json::Int(_) => "int",
        Json::Float(_) => "float",
        Json::Str(_) => "str",
        Json::Array(_) => "list",
        Json::Object(_) => "dict",
    }
}

/// `bool(value)`.
pub(crate) fn truthy(value: &Json) -> bool {
    match value {
        Json::Null => false,
        Json::Bool(flag) => *flag,
        Json::Int(number) => *number != 0,
        Json::Float(number) => *number != 0.0,
        Json::Str(text) => !text.is_empty(),
        Json::Array(items) => !items.is_empty(),
        Json::Object(items) => !items.is_empty(),
    }
}

/// `value.get(key)`: `None` for an absent key, AttributeError for a non-dict.
pub(crate) fn get<'a>(value: &'a Json, key: &str) -> Result<Option<&'a Json>> {
    match value {
        Json::Object(_) => Ok(value.get(key)),
        other => Err(error(
            "AttributeError",
            format!("'{}' object has no attribute 'get'", type_name(other)),
        )),
    }
}

/// `value.get(key, default)` with the default supplied as JSON.
pub(super) fn get_or(value: &Json, key: &str, default: Json) -> Result<Json> {
    Ok(get(value, key)?.cloned().unwrap_or(default))
}

/// `repr(value)` of the Python object `json.load` returns.
pub(crate) fn repr(value: &Json) -> String {
    match value {
        Json::Null => "None".into(),
        Json::Bool(true) => "True".into(),
        Json::Bool(false) => "False".into(),
        Json::Int(number) => number.to_string(),
        Json::Float(number) => crate::python_text::float_repr(*number),
        Json::Str(text) => crate::python_path::repr(text),
        Json::Array(items) => format!(
            "[{}]",
            items.iter().map(repr).collect::<Vec<_>>().join(", ")
        ),
        Json::Object(items) => format!(
            "{{{}}}",
            items
                .iter()
                .map(|(key, item)| format!("{}: {}", crate::python_path::repr(key), repr(item)))
                .collect::<Vec<_>>()
                .join(", ")
        ),
    }
}

/// `str(value)`.
pub(crate) fn py_str(value: &Json) -> String {
    match value {
        Json::Str(text) => text.clone(),
        other => repr(other),
    }
}

/// Python `==` between two parsed values (`1 == 1.0 == True`).
pub(crate) fn equal(left: &Json, right: &Json) -> bool {
    json_text::py_equal(left, right)
}

/// The current UTC time in microseconds since the epoch.
pub(super) fn now_micros() -> i128 {
    let elapsed = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default();
    elapsed.as_secs() as i128 * 1_000_000 + i128::from(elapsed.subsec_micros())
}

fn civil(days: i64) -> (i64, u32, u32) {
    // Howard Hinnant's days-to-civil conversion for the proleptic calendar.
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    (if m <= 2 { y + 1 } else { y }, m, d)
}

pub(super) fn days_from_civil(year: i64, month: u32, day: u32) -> i64 {
    let y = if month <= 2 { year - 1 } else { year };
    let era = y.div_euclid(400);
    let yoe = y - era * 400;
    let m = i64::from(month);
    let doy = (153 * (if m > 2 { m - 3 } else { m + 9 }) + 2) / 5 + i64::from(day) - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146_097 + doe - 719_468
}

/// `datetime.fromtimestamp(micros, utc).isoformat(timespec="seconds")`.
pub(super) fn iso_seconds(micros: i128) -> String {
    let seconds = micros.div_euclid(1_000_000) as i64;
    let (year, month, day) = civil(seconds.div_euclid(86_400));
    let within = seconds.rem_euclid(86_400);
    format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}+00:00",
        within / 3600,
        within / 60 % 60,
        within % 60
    )
}

pub(crate) fn iso_now() -> String {
    iso_seconds(now_micros())
}

/// `st_mtime` as Python computes it, then `fromtimestamp`'s round-half-even
/// microseconds.
pub(super) fn mtime_micros(meta: &fs::Metadata) -> i128 {
    let value = meta.mtime() as f64 + meta.mtime_nsec() as f64 * 1e-9;
    let whole = value.floor();
    let mut fraction = ((value - whole) * 1e6).round_ties_even();
    let mut seconds = whole as i128;
    if fraction >= 1e6 {
        seconds += 1;
        fraction -= 1e6;
    }
    seconds * 1_000_000 + fraction as i128
}

pub(super) fn file_mtime_iso(path: &Path) -> Result<String> {
    Ok(iso_seconds(mtime_micros(&stat(path)?)))
}

/// `datetime.fromisoformat(value)` as microseconds since the epoch (UTC for
/// an aware value) and whether the value carried an offset; `None` where
/// Python raises ValueError.
pub(super) fn parse_iso(value: &str) -> Option<(i128, bool)> {
    let bytes = value.as_bytes();
    let digits = |from: usize, count: usize| -> Option<i64> {
        let part = value.get(from..from + count)?;
        part.bytes()
            .all(|b| b.is_ascii_digit())
            .then(|| part.parse().ok())
            .flatten()
    };
    let (year, month, day, mut at) = if bytes.get(4) == Some(&b'-') {
        (digits(0, 4)?, digits(5, 2)?, digits(8, 2)?, 10)
    } else {
        (digits(0, 4)?, digits(4, 2)?, digits(6, 2)?, 8)
    };
    if bytes.get(at - 3) != Some(&b'-') && bytes.get(4) == Some(&b'-') {
        return None;
    }
    if !(1..=12).contains(&month) || year < 1 {
        return None;
    }
    let leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
    let lengths = [
        31,
        if leap { 29 } else { 28 },
        31,
        30,
        31,
        30,
        31,
        31,
        30,
        31,
        30,
        31,
    ];
    if day < 1 || day > lengths[month as usize - 1] {
        return None;
    }
    let mut micros = days_from_civil(year, month as u32, day as u32) as i128 * 86_400_000_000;
    let mut aware = false;
    if at < bytes.len() {
        at += 1;
        let hour = digits(at, 2)?;
        at += 2;
        let (mut minute, mut second, mut fraction) = (0, 0, 0i128);
        if bytes.get(at) == Some(&b':') {
            minute = digits(at + 1, 2)?;
            at += 3;
            if bytes.get(at) == Some(&b':') {
                second = digits(at + 1, 2)?;
                at += 3;
                if matches!(bytes.get(at), Some(b'.' | b',')) {
                    let start = at + 1;
                    let mut end = start;
                    while bytes.get(end).is_some_and(u8::is_ascii_digit) {
                        end += 1;
                    }
                    if end == start {
                        return None;
                    }
                    let part = &value[start..end.min(start + 6)];
                    fraction = format!("{part:0<6}").parse().ok()?;
                    at = end;
                }
            }
        }
        if hour > 23 || minute > 59 || second > 59 {
            return None;
        }
        micros += ((hour * 3600 + minute * 60 + second) as i128) * 1_000_000 + fraction;
        if at < bytes.len() {
            let sign = match bytes[at] {
                b'Z' if at + 1 == bytes.len() => {
                    return Some((micros, true));
                }
                b'+' => 1,
                b'-' => -1,
                _ => return None,
            };
            let offset_hour = digits(at + 1, 2)?;
            let mut rest = at + 3;
            let mut offset_minute = 0;
            if bytes.get(rest) == Some(&b':') {
                offset_minute = digits(rest + 1, 2)?;
                rest += 3;
            } else if rest < bytes.len() {
                offset_minute = digits(rest, 2)?;
                rest += 2;
            }
            if rest != bytes.len() || offset_hour > 23 || offset_minute > 59 {
                return None;
            }
            micros -= sign * ((offset_hour * 3600 + offset_minute * 60) as i128) * 1_000_000;
            aware = true;
        }
    }
    Some((micros, aware))
}

/// `(utc_now() - checked_at).total_seconds()` for an aware value.
pub(super) fn age_seconds(checked_micros: i128) -> f64 {
    (now_micros() - checked_micros) as f64 / 1e6
}

pub(crate) fn stat(path: &Path) -> Result<fs::Metadata> {
    fs::metadata(path).map_err(|e| os_error(e, path, None))
}

/// `path.exists()`.
pub(super) fn exists(path: &Path) -> bool {
    fs::metadata(path).is_ok()
}

/// `path.exists() and path.stat().st_size > 0`.
pub(super) fn nonempty(path: &Path) -> bool {
    fs::metadata(path).is_ok_and(|m| m.len() > 0)
}

pub(crate) fn sha256_file(path: &Path) -> Result<String> {
    crate::sha256::sha256_file(path).map_err(|e| os_error(e, path, None))
}
pub(super) fn md5_file(path: &Path) -> Result<String> {
    crate::md5::md5_file(path).map_err(|e| os_error(e, path, None))
}

/// `Path.resolve()` (non-strict).
pub(crate) fn resolve(path: &Path) -> Result<PathBuf> {
    crate::installer_config_paths::files::resolve(path)
}

/// `Path.expanduser()` for `~` and `~/...`.
pub(super) fn expanduser(path: &Path) -> Result<PathBuf> {
    let value = path.to_string_lossy();
    if value != "~" && !value.starts_with("~/") {
        return Ok(path.to_owned());
    }
    let home = std::env::var("HOME")
        .map_err(|_| error("RuntimeError", "Could not determine home directory."))?;
    Ok(path_value(&format!("{home}{}", &value[1..])))
}

/// `Path.absolute()`.
pub(super) fn absolute(path: &Path) -> Result<PathBuf> {
    Ok(PathBuf::from(crate::python_path::absolute(path)?))
}

/// `open_text(path)` then `json.load`: gzip when the name ends in `.gz`.
pub(crate) fn load_json(path: &Path) -> Result<Json> {
    let file = File::open(path).map_err(|e| os_error(e, path, None))?;
    let mut bytes = Vec::new();
    if path.to_string_lossy().ends_with(".gz") {
        crate::python_io::PythonGzipReader::new(file)
            .read_to_end(&mut bytes)
            .map_err(|e| match crate::python_io::GzipFailure::of(&e) {
                Some(crate::python_io::GzipFailure::EndOfFile) => error("EOFError", e.to_string()),
                Some(crate::python_io::GzipFailure::Bad(_)) => error("BadGzipFile", e.to_string()),
                Some(crate::python_io::GzipFailure::Zlib(_)) => error("error", e.to_string()),
                None => os_error(e, path, None),
            })?;
    } else {
        File::open(path)
            .and_then(|mut f| f.read_to_end(&mut bytes))
            .map_err(|e| os_error(e, path, None))?;
    }
    let text = crate::python_io::decode_utf8(&bytes)
        .map_err(|e| error("UnicodeDecodeError", e))?
        .to_owned();
    json_text::parse(&text, false).map_err(|e| error("JSONDecodeError", e))
}

/// The script's `read_json`: `{}` for an absent file.
pub(super) fn read_json(path: &Path) -> Result<Json> {
    if !exists(path) {
        return Ok(Json::Object(vec![]));
    }
    load_json(path)
}

/// `json.dumps(value, indent=2, sort_keys=True)`.
pub(super) fn dumps_indented(value: &Json) -> Result<String> {
    json_text::dumps_with_options(
        value,
        DumpsOptions {
            indent: Some(2),
            sort_keys: true,
            ..DumpsOptions::default()
        },
    )
    .map_err(|e| error("ValueError", e))
}

/// `atomic_write_text`: the whole text to a temporary sibling (keeping a
/// `.gz` suffix, gzip with mtime 0 and no name), then one rename.
pub(super) fn atomic_write(path: &Path, body: &str) -> Result<()> {
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    fs::create_dir_all(parent).map_err(|e| os_error(e, parent, None))?;
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let gz = path.to_string_lossy().ends_with(".gz");
    let temporary = if gz {
        path.with_file_name(format!("{}.tmp.gz", &name[..name.len() - 3]))
    } else {
        path.with_file_name(format!("{name}.tmp"))
    };
    let result = (|| -> Result<()> {
        let file = File::create(&temporary).map_err(|e| os_error(e, &temporary, None))?;
        if gz {
            let mut writer = crate::gzip_writer::GzipWriter::new(
                file,
                crate::gzip_writer::GzipOptions::default(),
            )?;
            writer.write_all(body.as_bytes())?;
            writer.finish()?;
        } else {
            let mut file = file;
            file.write_all(body.as_bytes())?;
        }
        fs::rename(&temporary, path).map_err(|e| os_error(e, &temporary, Some(path)))?;
        Ok(())
    })();
    let _ = fs::remove_file(&temporary);
    result
}

/// The script's `write_json`.
pub(super) fn write_json(path: &Path, value: &Json) -> Result<()> {
    atomic_write(path, &format!("{}\n", dumps_indented(value)?))
}

/// The script's `BuildLock`: an exclusive file, removed when stale.
struct BuildLock {
    path: PathBuf,
}

impl BuildLock {
    fn enter(path: &Path, stale_hours: f64) -> Result<Self> {
        if let Ok(meta) = fs::metadata(path) {
            let now = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap_or_default()
                .as_secs_f64();
            let mtime = meta.mtime() as f64 + meta.mtime_nsec() as f64 * 1e-9;
            if (now - mtime) / 3600.0 > stale_hours {
                fs::remove_file(path).map_err(|e| os_error(e, path, None))?;
            } else {
                return Err(error(
                    "SystemExit",
                    format!(
                        "Lock exists: {}. Another cache build may be running.",
                        path.display()
                    ),
                ));
            }
        }
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).map_err(|e| os_error(e, parent, None))?;
        }
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)
            .map_err(|e| os_error(e, path, None))?;
        file.write_all(
            format!("pid={}\nstarted_at={}\n", std::process::id(), iso_now()).as_bytes(),
        )
        .map_err(|e| os_error(e, path, None))?;
        Ok(Self {
            path: path.to_owned(),
        })
    }
}

impl Drop for BuildLock {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

pub(super) struct Options {
    source_root: PathBuf,
    cache_dir: PathBuf,
    output_json: Option<PathBuf>,
    shared_raw_dir: String,
    gofcards_exact_variants: PathBuf,
    hgnc_table: PathBuf,
    output_schema: PathBuf,
    clinvar_vcv: Option<PathBuf>,
    clinvar_only_refresh: bool,
    fetch_sources_only: bool,
    force: bool,
    validate_only: Option<PathBuf>,
    stale_lock_hours: f64,
    transfer: Transfer,
    sources: Vec<Spec>,
    panelapp: Spec,
}

fn int_value(name: &str, value: &str) -> std::result::Result<i128, String> {
    crate::python_text::python_int(value).ok_or_else(|| {
        format!(
            "argument --{name}: invalid int value: {}",
            crate::python_path::repr(value)
        )
    })
}

fn parse_options(args: &[String]) -> Result<Option<Options>> {
    let converted = parse_long_options_converted(args, VALUES, FLAGS, |name, value| {
        match name {
            "clinvar-min-review-stars" => {
                let number = int_value(name, value)?;
                if !(2..=4).contains(&number) {
                    return Err(format!(
                        "argument --{name}: invalid choice: {number} (choose from 2, 3, 4)"
                    ));
                }
            }
            "clinvar-max-download-seconds" | "max-panelapp-panels" | "timeout" | "retries" => {
                int_value(name, value)?;
            }
            "download-tool" if !["auto", "urllib", "curl"].contains(&value) => {
                return Err(format!(
                    "argument --{name}: invalid choice: {} (choose from 'auto', 'urllib', 'curl')",
                    crate::python_path::repr(value)
                ));
            }
            "stale-lock-hours" => {
                crate::python_text::parse_float(value).ok_or_else(|| {
                    format!(
                        "argument --{name}: invalid float value: {}",
                        crate::python_path::repr(value)
                    )
                })?;
            }
            "g2p-ddg2p-url" | "clingen-dosage-grch38-url" | "panelapp-panels-url" => {
                crate::mechanism_source::check_source_url(value).map_err(|_| {
                    format!(
                        "argument --{name}: invalid source_url value: {}",
                        crate::python_path::repr(value)
                    )
                })?;
            }
            _ => {}
        }
        Ok(())
    })
    .map_err(invalid)?;
    let Some(occurrences) = converted else {
        println!("{HELP}");
        return Ok(None);
    };
    let values: HashMap<String, String> = occurrences.into_iter().collect();
    require_options(&values, &["source-root"]).map_err(invalid)?;
    let root = path_value(&values["source-root"]);
    let cache_default = root.join("data/patho_mechanism");
    let path = |name: &str, default: PathBuf| values.get(name).map_or(default, |v| path_value(v));
    let int = |name: &str, default: i128| {
        values
            .get(name)
            .map_or(default, |v| crate::python_text::python_int(v).unwrap())
    };
    let mut sources = default_specs();
    for spec in &mut sources {
        let key = match spec.name {
            "g2p_ddg2p" => "g2p-ddg2p-url",
            "clingen_dosage_gene_grch38" => "clingen-dosage-grch38-url",
            _ => continue,
        };
        if let Some(url) = values.get(key).filter(|url| !url.is_empty()) {
            spec.url = url.clone();
        }
    }
    let mut panelapp = panelapp_spec();
    if let Some(url) = values
        .get("panelapp-panels-url")
        .filter(|url| !url.is_empty())
    {
        panelapp.url = url.clone();
    }
    Ok(Some(Options {
        cache_dir: path("cache-dir", cache_default.clone()),
        output_json: values.get("output-json").map(|v| path_value(v)),
        shared_raw_dir: values
            .get("shared-raw-dir")
            .cloned()
            .unwrap_or_else(|| cache_default.join("raw").display().to_string()),
        gofcards_exact_variants: path(
            "gofcards-exact-variants",
            root.join("data/gofcards/gofcards_exact_gof.json.gz"),
        ),
        hgnc_table: path(
            "hgnc-table",
            root.join("data/hgnc/hgnc_non_alt_loci_set_sha256_890b18743cc8.tsv"),
        ),
        output_schema: path(
            "output-schema",
            cache_default.join(NONLOF_ASSERTIONS_SCHEMA_FILENAME),
        ),
        clinvar_vcv: values.get("clinvar-vcv").map(|v| path_value(v)),
        clinvar_only_refresh: values.contains_key("clinvar-only-refresh"),
        fetch_sources_only: values.contains_key("fetch-sources-only"),
        force: values.contains_key("force"),
        validate_only: values.get("validate-only").map(|v| path_value(v)),
        stale_lock_hours: values
            .get("stale-lock-hours")
            .map_or(24.0, |v| crate::python_text::parse_float(v).unwrap()),
        transfer: Transfer {
            timeout: int("timeout", 120),
            retries: int("retries", 3),
            proxy_url: values.get("proxy-url").cloned().unwrap_or_default(),
            download_tool: values
                .get("download-tool")
                .cloned()
                .unwrap_or_else(|| "auto".into()),
            max_download_seconds: int("clinvar-max-download-seconds", 86_400),
            max_panels: values
                .get("max-panelapp-panels")
                .map(|v| crate::python_text::python_int(v).unwrap()),
        },
        sources,
        panelapp,
        source_root: root,
    }))
}

fn progress(line: &str) {
    eprintln!("{line}");
}

fn print_json(value: &Json) -> Result<()> {
    println!("{}", dumps_indented(value)?);
    std::io::stdout().flush().ok();
    Ok(())
}

fn summary_update(path: &Path, updates: Vec<(&str, Json)>) -> Result<()> {
    let mut summary = read_json(path)?;
    if !matches!(summary, Json::Object(_)) {
        return Err(error(
            "AttributeError",
            format!("'{}' object has no attribute 'update'", type_name(&summary)),
        ));
    }
    for (key, value) in updates {
        summary.set(key, value);
    }
    write_json(path, &summary)
}

fn run(options: Options) -> Result<()> {
    if let Some(target) = &options.validate_only {
        let target = resolve(target)?;
        let canonical = read_json(&target)?;
        if !truthy(&canonical) {
            return Err(error("FileNotFoundError", target.display().to_string()));
        }
        schema::validate_canonical(&canonical, &resolve(&options.output_schema)?)?;
        let genes = match &canonical {
            Json::Object(items) => items.len() as i128 - 1,
            Json::Array(items) => items.len() as i128 - 1,
            Json::Str(text) => text.chars().count() as i128 - 1,
            _ => 0,
        };
        println!(
            "schema_validation=passed genes={genes} json={}",
            target.display()
        );
        return Ok(());
    }
    let selected = match &options.clinvar_vcv {
        Some(path) => Some(fetch::selected_clinvar_vcv_metadata(path)?),
        None => None,
    };
    if let Some(meta) = &selected {
        if !options.fetch_sources_only {
            fetch::require_selected_clinvar_injection(meta, &options.gofcards_exact_variants)?;
        }
    }
    let cache_dir = resolve(&options.cache_dir)?;
    let raw_dir = if crate::python_text::python_strip(&options.shared_raw_dir).is_empty() {
        cache_dir.join("raw")
    } else {
        resolve(&expanduser(&path_value(&options.shared_raw_dir))?)?
    };
    let raw_dir_is_shared = raw_dir != cache_dir.join("raw");
    let manifest_path = cache_dir.join(NONLOF_SOURCE_MANIFEST_FILENAME);
    let manifest_tsv_path = cache_dir.join(NONLOF_SOURCE_MANIFEST_TSV_FILENAME);
    let json_path = match &options.output_json {
        Some(path) => resolve(&expanduser(path)?)?,
        None => cache_dir.join(NONLOF_ASSERTIONS_FILENAME),
    };
    let summary_path = cache_dir.join(NONLOF_RUN_SUMMARY_FILENAME);
    let lock_path = if raw_dir_is_shared {
        raw_dir.join(".build.lock")
    } else {
        cache_dir.join(".build.lock")
    };
    let _lock = BuildLock::enter(&lock_path, options.stale_lock_hours)?;
    let previous_manifest = read_json(&manifest_path)?;
    let previous_sources = get_or(&previous_manifest, "sources", Json::Object(vec![]))?;
    let mut source_meta: Vec<(String, Json)> = Vec::new();
    let all_specs: Vec<Spec> = options
        .sources
        .iter()
        .cloned()
        .chain([options.panelapp.clone()])
        .collect();
    if options.clinvar_only_refresh {
        for spec in &all_specs {
            progress(&format!("carrying_source={}", spec.name));
            source_meta.push((
                spec.name.into(),
                fetch::carry_forward(spec, &raw_dir, &previous_sources)?,
            ));
        }
    } else {
        for spec in &options.sources {
            progress(&format!("checking_source={}", spec.name));
            source_meta.push((
                spec.name.into(),
                fetch::fetch_static_source(
                    spec,
                    &raw_dir,
                    &previous_sources,
                    options.force,
                    &options.transfer,
                )?,
            ));
        }
        progress(&format!("checking_source={}", options.panelapp.name));
        source_meta.push((
            options.panelapp.name.into(),
            fetch::fetch_panelapp(
                &options.panelapp,
                &raw_dir,
                &previous_sources,
                options.force,
                &options.transfer,
            )?,
        ));
    }
    let clinvar = clinvar_spec();
    progress(&format!("checking_source={}", clinvar.name));
    source_meta.push((
        clinvar.name.into(),
        match &selected {
            Some(meta) => meta.clone(),
            None => fetch::fetch_clinvar_vcv(
                &raw_dir,
                &previous_sources,
                options.force,
                &options.transfer,
            )?,
        },
    ));
    let mut expected: Vec<Json> = options
        .sources
        .iter()
        .map(|spec| text(spec.raw_filename))
        .chain([text(options.panelapp.raw_filename)])
        .collect();
    if selected.is_none() {
        expected.push(text(clinvar.raw_filename));
    }
    let built_at = iso_now();
    let mut manifest = object([
        ("schema_version", text("2.0")),
        ("cache_dir", text(cache_dir.display().to_string())),
        ("raw_dir", text(raw_dir.display().to_string())),
        (
            "raw_dir_mode",
            text(if raw_dir_is_shared {
                "shared"
            } else {
                "cache_local"
            }),
        ),
        ("lock_path", text(lock_path.display().to_string())),
        (
            "proxy_mode",
            text(if options.transfer.proxy_url.is_empty() {
                "environment_or_none"
            } else {
                "explicit_proxy_url"
            }),
        ),
        (
            "download_tool",
            text(options.transfer.download_tool.clone()),
        ),
        ("expected_raw_filenames", Json::Array(expected)),
        ("built_at_utc", text(built_at.clone())),
        ("sources", Json::Object(source_meta.clone())),
    ]);
    if let Some(meta) = &selected {
        manifest.set(
            "selected_clinvar_vcv_path",
            get_or(meta, "raw_path", Json::Null)?,
        );
    }
    if options.fetch_sources_only {
        if let Some(previous) = get(&previous_manifest, "build_inputs")? {
            if !matches!(previous, Json::Null) {
                manifest.set("build_inputs", previous.clone());
            }
        }
        write_json(&manifest_path, &manifest)?;
        let rows = fetch::write_source_manifest_tsv(&manifest_tsv_path, &manifest)?;
        summary_update(
            &summary_path,
            vec![
                ("last_checked_at_utc", text(built_at.clone())),
                ("status", text("fetched_sources_only")),
                ("source_manifest_rows", Json::Int(rows as i128)),
            ],
        )?;
        let clinvar_meta = &source_meta.last().unwrap().1;
        print_json(&object([
            ("status", text("fetched_sources_only")),
            (
                "clinvar_vcv_xml",
                get_or(clinvar_meta, "raw_path", text(""))?,
            ),
            (
                "source_manifest_json",
                text(manifest_path.display().to_string()),
            ),
        ]))?;
        return Ok(());
    }
    let gofcards_exact_path = resolve(&options.gofcards_exact_variants)?;
    let hgnc_path = resolve(&options.hgnc_table)?;
    let schema_path = resolve(&options.output_schema)?;
    let scripts = options.source_root.join("scripts");
    let builder = resolve(&scripts.join("build_gene_nonlof_mechanism_cache.py"))?;
    let mut build_inputs = object([
        ("builder_sha256", text(sha256_file(&builder)?)),
        (
            "clinvar_parser_sha256",
            text(sha256_file(&scripts.join("clinvar_vcv.py"))?),
        ),
        ("schema_sha256", text(sha256_file(&schema_path)?)),
        ("hgnc_sha256", text(sha256_file(&hgnc_path)?)),
        (
            "gofcards_exact_sha256",
            text(sha256_file(&gofcards_exact_path)?),
        ),
        (
            "gofcards_exact_path",
            text(gofcards_exact_path.display().to_string()),
        ),
        ("source_sha256", {
            let mut names: Vec<(String, Json)> = source_meta
                .iter()
                .map(|(name, meta)| {
                    Ok((
                        name.clone(),
                        text(py_str(&get_or(meta, "sha256", text(""))?)),
                    ))
                })
                .collect::<Result<_>>()?;
            names.sort_by(|a, b| a.0.cmp(&b.0));
            Json::Object(names)
        }),
    ]);
    if let Some(meta) = &selected {
        build_inputs.set(
            "selected_clinvar_vcv_path",
            get_or(meta, "raw_path", Json::Null)?,
        );
    }
    manifest.set("build_inputs", build_inputs.clone());
    if !options.force
        && exists(&json_path)
        && get(&previous_manifest, "build_inputs")?
            .is_some_and(|previous| equal(previous, &build_inputs))
    {
        schema::validate_canonical(&read_json(&json_path)?, &schema_path)?;
        write_json(&manifest_path, &manifest)?;
        let rows = fetch::write_source_manifest_tsv(&manifest_tsv_path, &manifest)?;
        summary_update(
            &summary_path,
            vec![
                ("last_checked_at_utc", text(built_at.clone())),
                ("status", text("checked_inputs_unchanged")),
                ("source_manifest_rows", Json::Int(rows as i128)),
            ],
        )?;
        print_json(&object([
            ("status", text("checked_inputs_unchanged")),
            (
                "nonlof_assertions_json",
                text(json_path.display().to_string()),
            ),
            ("schema_validation", text("passed")),
        ]))?;
        return Ok(());
    }
    let rows = parse::parse_all_sources(&raw_dir, &all_specs)?;
    let hgnc = parse::load_hgnc_mapping(&hgnc_path)?;
    progress(&format!("hgnc_mapping: {} entries", hgnc.len()));
    let (exact, stats) = build::load_gofcards_exact_records(&gofcards_exact_path, &hgnc)?;
    let mut result = build::assertions_json(&rows, &hgnc, &exact)?;
    let clinvar_annotated = build::integrate(&mut result, &gofcards_exact_path, &stats, &exact)?;
    progress(&format!(
        "clinvar_vcv: attached upstream; variants_with_conditions={clinvar_annotated}"
    ));
    schema::validate_canonical(&result, &schema_path)?;
    progress(&format!(
        "json_schema_validation=passed schema={}",
        schema_path.display()
    ));
    write_json(&json_path, &result)?;
    write_json(&manifest_path, &manifest)?;
    let manifest_rows = fetch::write_source_manifest_tsv(&manifest_tsv_path, &manifest)?;
    let meta = result.get("_meta").unwrap();
    let genes = meta.get("total_genes").cloned().unwrap_or(Json::Int(0));
    let unmapped = match meta.get("unmapped_symbols") {
        Some(Json::Array(items)) => items.len(),
        _ => 0,
    };
    progress(&format!(
        "json_output: {} genes, {unmapped} unmapped to {}",
        py_str(&genes),
        json_path.display()
    ));
    let summary = object([
        ("built_at_utc", text(built_at.clone())),
        ("last_checked_at_utc", text(built_at)),
        ("status", text("rebuilt")),
        ("cache_dir", text(cache_dir.display().to_string())),
        ("total_genes", genes),
        ("unmapped_symbols", Json::Int(unmapped as i128)),
        ("source_manifest_rows", Json::Int(manifest_rows as i128)),
        (
            "clinvar_vcv",
            object([
                ("source", text("inject_clinvar_into_gofcards.py")),
                (
                    "variants_carrying_clinvar_conditions",
                    Json::Int(clinvar_annotated as i128),
                ),
            ]),
        ),
        (
            "outputs",
            object([
                (
                    "source_manifest_json",
                    text(manifest_path.display().to_string()),
                ),
                (
                    "source_manifest_tsv",
                    text(manifest_tsv_path.display().to_string()),
                ),
                (
                    "nonlof_assertions_json",
                    text(json_path.display().to_string()),
                ),
            ]),
        ),
    ]);
    write_json(&summary_path, &summary)?;
    print_json(&summary)
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(options) = parse_options(args)? else {
        return Ok(());
    };
    builder_error::finish(run(options))
}
