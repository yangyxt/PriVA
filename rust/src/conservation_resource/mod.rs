//! Optional port of the Conservation resource guard. BigWig reads remain in
//! Bio::DB::BigFile and the publisher BSD checksum remains `sum -r`.
//! This checks resource identity/integrity, never scientific qualification.
mod process;
mod reference;

use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use crate::python_args::{Parsed, parse_long_options};
use std::fs;
use std::io;
use std::os::unix::fs::MetadataExt;
use std::path::{Path, PathBuf};

const URL: &str = "https://ftp.ensembl.org/pub/release-114/compara/conservation_scores/91_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.GRCh38.bw";
const PUBLISHER_BYTES: i128 = 9_611_134_022;
const PUBLISHER_SUM: [i128; 2] = [30633, 9385874];
const HELP: &str = "usage: priva conservation-resource {validate|publish} --input BIGWIG --assembly ASSEMBLY --reference FASTA --tool-prefix BIN_DIR [--source-url URL] [--publish-to BIGWIG]\n\
Preserve the original resource checks and hard-link publication order.\n\
--tool-prefix is the lexical parent of the selected environment's Python\n\
interpreter (its bin directory). Its perl is preferred, with the original\n\
PATH fallback if absent. The priva executable's directory is never used.\n\
publish requires --publish-to and an explicit nonempty --source-url.\n\
Existing installer and LoFtee callers remain on the Python helper.";

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
fn object(items: impl IntoIterator<Item = (&'static str, Json)>) -> Json {
    Json::Object(items.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}
fn string(value: impl Into<String>) -> Json {
    Json::Str(value.into())
}
fn optional(value: Option<&str>) -> Json {
    value.map_or(Json::Null, string)
}

/// Pathlib collapses `.` and redundant separators but keeps `..` and two
/// leading slashes. In particular absolute() must not resolve a symlink.
fn path_value(value: &str) -> PathBuf {
    let prefix = if value.starts_with("//") && !value.starts_with("///") {
        "//"
    } else if value.starts_with('/') {
        "/"
    } else {
        ""
    };
    let parts = value
        .split('/')
        .filter(|part| !part.is_empty() && *part != ".")
        .collect::<Vec<_>>();
    let text = format!("{prefix}{}", parts.join("/"));
    PathBuf::from(if text.is_empty() {
        ".".to_owned()
    } else {
        text
    })
}
fn absolute(path: &Path) -> Result<String> {
    if path.is_absolute() {
        return Ok(path.display().to_string());
    }
    Ok(path_value(&format!(
        "{}/{}",
        std::env::current_dir()?.display(),
        path.display()
    ))
    .display()
    .to_string())
}
fn repr(text: &str) -> String {
    let quote = if text.contains('\'') && !text.contains('"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::from(quote);
    for c in text.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c == quote => {
                out.push('\\');
                out.push(c);
            }
            c if c.is_control() && (c as u32) <= 255 => {
                out.push_str(&format!("\\x{:02x}", c as u32))
            }
            c if c.is_control() || matches!(c, '\u{2028}' | '\u{2029}') => {
                out.push_str(&format!("\\u{:04x}", c as u32))
            }
            c => out.push(c),
        }
    }
    out.push(quote);
    out
}
fn os_error(
    cause: io::Error,
    path: &Path,
    destination: Option<&Path>,
) -> Box<dyn std::error::Error> {
    let category = match cause.kind() {
        io::ErrorKind::NotFound => "FileNotFoundError",
        io::ErrorKind::PermissionDenied => "PermissionError",
        io::ErrorKind::AlreadyExists => "FileExistsError",
        io::ErrorKind::IsADirectory => "IsADirectoryError",
        io::ErrorKind::NotADirectory => "NotADirectoryError",
        _ => "OSError",
    };
    let message = if let Some(errno) = cause.raw_os_error() {
        let description = cause.to_string();
        let description = description
            .strip_suffix(&format!(" (os error {errno})"))
            .unwrap_or(&description);
        format!(
            "[Errno {errno}] {description}: {}{}",
            repr(&path.display().to_string()),
            destination.map_or(String::new(), |p| format!(
                " -> {}",
                repr(&p.display().to_string())
            ))
        )
    } else {
        cause.to_string()
    };
    error(category, message)
}
fn is_file(path: &Path) -> Result<bool> {
    match fs::metadata(path) {
        Ok(m) => Ok(m.is_file()),
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
fn metadata(path: &Path) -> Result<fs::Metadata> {
    fs::metadata(path).map_err(|e| os_error(e, path, None))
}
fn identity(path: &Path) -> Result<Json> {
    let m = metadata(path)?;
    Ok(object([
        ("device", Json::Int(m.dev() as i128)),
        ("inode", Json::Int(m.ino() as i128)),
        ("bytes", Json::Int(m.len() as i128)),
        (
            "mtime_ns",
            Json::Int(m.mtime() as i128 * 1_000_000_000 + m.mtime_nsec() as i128),
        ),
        (
            "ctime_ns",
            Json::Int(m.ctime() as i128 * 1_000_000_000 + m.ctime_nsec() as i128),
        ),
    ]))
}
fn sha256(path: &Path) -> Result<String> {
    crate::sha256::sha256_file(path).map_err(|e| os_error(e, path, None))
}
fn normalize(name: &str) -> String {
    let value = name.strip_prefix("chr").unwrap_or(name);
    if value == "M" {
        "MT".to_owned()
    } else {
        value.to_owned()
    }
}
fn chr1(assembly: &str) -> &'static str {
    if assembly == "GRCh37" {
        "249250621"
    } else {
        "248956422"
    }
}

fn validate(
    path: &Path,
    assembly: &str,
    reference: &Path,
    source_url: Option<&str>,
    prefix: &Path,
) -> Result<Json> {
    let assembly = match assembly {
        "hg19" | "GRCh37" => "GRCh37",
        "hg38" | "GRCh38" => "GRCh38",
        _ => {
            return Err(error(
                "ValueError",
                "Conservation requires hg19/GRCh37 or hg38/GRCh38",
            ));
        }
    };
    require(
        is_file(path)? && metadata(path)?.len() > 0,
        format!("Missing or empty Conservation BigWig: {}", path.display()),
    )?;
    let before = identity(path)?;
    let known = source_url == Some(URL);
    if known {
        require(
            assembly == "GRCh38",
            "Archived release114 Conservation requires GRCh38",
        )?;
        require(
            matches!(before.get("bytes"), Some(Json::Int(PUBLISHER_BYTES))),
            "Archived release114 Conservation publisher byte count mismatch; resource retained",
        )?;
    }
    let (lengths, reference_record) = reference::lengths(reference, assembly)?;
    let bigwig = process::bigwig(path, prefix)?;
    let chroms = bigwig
        .get("chroms")
        .and_then(Json::as_object)
        .expect("checked chromosome dictionary");
    let mut normalized = std::collections::BTreeMap::new();
    for (name, length) in chroms {
        let key = normalize(name);
        require(
            normalized
                .get(&key)
                .is_none_or(|old| json_text::py_equal(old, length)),
            format!("Conflicting normalized chromosome lengths in BigWig: {name}"),
        )?;
        normalized.insert(key, length.clone());
    }
    let expected: i128 = chr1(assembly).parse().unwrap();
    require(
        normalized
            .get("1")
            .is_some_and(|v| json_text::py_equal(v, &Json::Int(expected))),
        format!(
            "Conservation BigWig chr1 length does not identify {assembly}: {}",
            path.display()
        ),
    )?;
    let shared: Vec<_> = normalized
        .keys()
        .filter(|name| lengths.contains_key(*name))
        .cloned()
        .collect();
    for name in &shared {
        let measured = match &normalized[name] {
            Json::Int(v) => v.to_string(),
            v => json_text::dumps(v, false).map_err(invalid)?,
        };
        require(
            measured == lengths[name],
            format!(
                "Conservation/reference chromosome length mismatch for {name}: {measured} != {}",
                lengths[name]
            ),
        )?;
    }
    let measured_bsd = if known {
        Some(process::bsd_sum(path)?)
    } else {
        None
    };
    if let Some(sum) = &measured_bsd {
        require(
            *sum == PUBLISHER_SUM,
            "Archived release114 Conservation publisher BSD checksum mismatch; resource retained",
        )?;
    }
    let digest = sha256(path)?;
    require(
        json_text::py_equal(&identity(path)?, &before),
        "Conservation resource changed during validation; resource retained",
    )?;
    let field = |key: &str| {
        bigwig
            .get(key)
            .cloned()
            .ok_or_else(|| error("KeyError", repr(key)))
    };
    Ok(object([
        ("path", string(absolute(path)?)),
        ("assembly", string(assembly)),
        ("sha256", string(digest)),
        ("identity", before),
        ("reader", field("reader")?),
        ("reader_module", field("reader_module")?),
        ("chromosome_count", Json::Int(chroms.len() as i128)),
        (
            "shared_reference_chromosomes",
            Json::Array(shared.into_iter().map(string).collect()),
        ),
        (
            "extra_non_reference_chromosomes",
            Json::Array(
                normalized
                    .keys()
                    .filter(|name| !lengths.contains_key(*name))
                    .cloned()
                    .map(string)
                    .collect(),
            ),
        ),
        ("reference", reference_record),
        ("source_url", optional(source_url)),
        (
            "source_kind",
            string(if source_url.is_some_and(|v| !v.is_empty()) {
                "download"
            } else {
                "configured_existing_identity"
            }),
        ),
        ("release", if known { Json::Int(114) } else { Json::Null }),
        (
            "analysis",
            if known {
                string("91_mammals.gerp_conservation_score")
            } else {
                Json::Null
            },
        ),
        (
            "publisher_byte_count",
            if known {
                Json::Int(PUBLISHER_BYTES)
            } else {
                Json::Null
            },
        ),
        (
            "publisher_bsd_sum",
            measured_bsd.map_or(Json::Null, |sum| {
                Json::Array(sum.into_iter().map(Json::Int).collect())
            }),
        ),
        (
            "publisher_cryptographic_checksum_claimed",
            Json::Bool(false),
        ),
        ("scientific_qualified", Json::Bool(false)),
    ]))
}
fn publish(
    path: &Path,
    destination: &Path,
    assembly: &str,
    reference: &Path,
    source_url: Option<&str>,
    prefix: &Path,
) -> Result<Json> {
    require(
        source_url.is_some_and(|url| !url.is_empty()),
        "Download publication requires its explicit source URL",
    )?;
    require(
        absolute(path)? != absolute(destination)?,
        "Download staging and final paths must differ",
    )?;
    let symlink = match fs::symlink_metadata(path) {
        Ok(m) => m.is_symlink(),
        Err(e)
            if matches!(
                e.raw_os_error(),
                Some(libc::ENOENT | libc::ENOTDIR | libc::EBADF | libc::ELOOP)
            ) =>
        {
            false
        }
        Err(e) => return Err(os_error(e, path, None)),
    };
    require(!symlink, "Conservation staging file must not be a symlink")?;
    let mut record = validate(path, assembly, reference, source_url, prefix)?;
    require(
        json_text::py_equal(&identity(path)?, record.get("identity").unwrap()),
        "Conservation staging file changed before publication",
    )?;
    // Keep the original link -> unlink -> final-stat order. A failed unlink
    // leaves BOTH names; neither an existing target nor staging is removed.
    fs::hard_link(path, destination).map_err(|e| os_error(e, path, Some(destination)))?;
    fs::remove_file(path).map_err(|e| os_error(e, path, None))?;
    record.set("path", string(absolute(destination)?));
    record.set("identity", identity(destination)?);
    Ok(record)
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args.is_empty() {
        return Err(invalid(HELP));
    }
    if matches!(args[0].as_str(), "--help" | "-h") {
        println!("{HELP}");
        return Ok(());
    }
    let procedure = args[0].as_str();
    if procedure == "bsd-sum" {
        let Parsed::Values(options) =
            parse_long_options(&args[1..], &["input"]).map_err(invalid)?
        else {
            println!(
                "usage: priva conservation-resource bsd-sum --input FILE\nRun the original BSD checksum procedure using sum -r and LC_ALL=C."
            );
            return Ok(());
        };
        let path = path_value(
            options
                .get("input")
                .ok_or_else(|| invalid("missing --input"))?,
        );
        return crate::builder_error::finish((|| {
            let value = Json::Array(
                process::bsd_sum(&path)?
                    .into_iter()
                    .map(Json::Int)
                    .collect(),
            );
            println!("{}", json_text::dumps(&value, false).map_err(invalid)?);
            Ok(())
        })());
    }
    if !matches!(procedure, "validate" | "publish") {
        return Err(invalid(HELP));
    }
    let Parsed::Values(options) = parse_long_options(
        &args[1..],
        &[
            "input",
            "assembly",
            "reference",
            "tool-prefix",
            "source-url",
            "publish-to",
        ],
    )
    .map_err(invalid)?
    else {
        println!("{HELP}");
        return Ok(());
    };
    let required = |name: &str| {
        options
            .get(name)
            .ok_or_else(|| invalid(format!("missing --{name}")))
    };
    let input = path_value(required("input")?);
    let assembly = required("assembly")?;
    let reference = path_value(required("reference")?);
    let prefix = path_value(required("tool-prefix")?);
    let destination = if procedure == "publish" {
        Some(path_value(required("publish-to")?))
    } else {
        if options.contains_key("publish-to") {
            return Err(invalid("--publish-to requires the publish procedure"));
        }
        None
    };
    let url = options.get("source-url").map(String::as_str);
    crate::builder_error::finish((|| {
        let record = if let Some(target) = destination {
            publish(&input, &target, assembly, &reference, url, &prefix)?
        } else {
            validate(&input, assembly, &reference, url, &prefix)?
        };
        println!("{}", json_text::dumps(&record, true).map_err(invalid)?);
        Ok(())
    })())
}
