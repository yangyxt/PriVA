//! Explicit acquisition-source checks and receipts. Genuine repository bundles
//! retain their existing existence/nonzero-only reuse path, outside this CLI.
mod csv_input;
mod url;
use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use crate::python_args::{parse_long_option_occurrences, require_options};
use crate::python_path::{is_file, metadata, os_error, path_value, tempfile_parent};
use std::collections::HashMap;
use std::ffi::CString;
use std::fs::{self, File};
use std::io::{Read, Write};
use std::os::fd::FromRawFd;
use std::os::unix::ffi::OsStrExt;
use std::path::{Path, PathBuf};
pub(crate) use url::check as check_source_url;
const HELP: &str = "usage: priva mechanism-source validate --input FILE --kind {clingen_erepo,clingen_dosage,clingen_dosage_runtime,g2p_ddg2p,orphadata}\n\
       priva mechanism-source {record-source|check-source} --input FILE --url URL\n\
check-source accepts --receipt JSON and --metadata-only. No downloads or caller activation.";
const KINDS: &[&str] = &[
    "clingen_erepo",
    "clingen_dosage",
    "clingen_dosage_runtime",
    "g2p_ddg2p",
    "orphadata",
];
fn error(category: &'static str, message: impl Into<String>) -> Box<dyn std::error::Error> {
    crate::builder_error::category(category, invalid(message))
}
fn require(ok: bool, message: impl Into<String>) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(error("ValueError", message))
    }
}

fn normalize(text: &str) -> String {
    let lowered = text.to_lowercase();
    let mut result = String::new();
    let mut pending = false;
    for c in lowered.chars() {
        if c.is_ascii_lowercase() || c.is_ascii_digit() {
            if pending && !result.is_empty() {
                result.push('_');
            }
            result.push(c);
            pending = false;
        } else {
            pending = true;
        }
    }
    result
}
fn finite_number(text: &str) -> bool {
    let text = crate::python_text::python_strip(text);
    let converted: String = text
        .chars()
        .map(|c| {
            crate::python_text::decimal_digit_value(c).map_or(c, |n| char::from(b'0' + n as u8))
        })
        .collect();
    crate::python_text::parse_float(&converted).is_some_and(f64::is_finite)
}
fn validate(path: &Path, kind: &str) -> Result<()> {
    require(
        is_file(path)? && metadata(path)?.len() > 0,
        format!(
            "Downloaded {kind} resource is missing or empty: {}",
            path.display()
        ),
    )?;
    let file = File::open(path).map_err(|e| os_error(e, path, None))?;
    if kind == "orphadata" {
        return crate::interpro_mapping_builder::validate_product6(file);
    }
    let mut rows = csv_input::Rows::new(file, kind);
    let header = rows.next()?.unwrap_or_default();
    let exact: &[&str] = match kind {
        "clingen_erepo" => &["ClinVar Variation Id", "Applied Evidence Codes (Met)"],
        "clingen_dosage_runtime" => &[
            "#Gene Symbol",
            "Genomic Location",
            "Haploinsufficiency Score",
            "Haploinsufficiency Description",
            "Haploinsufficiency Disease ID",
            "Triplosensitivity Score",
        ],
        _ => &[],
    };
    require(
        exact.iter().all(|name| header.iter().any(|v| v == name)),
        format!(
            "Downloaded {kind} lacks the original runtime reader's exact columns: {}",
            path.display()
        ),
    )?;
    let groups: &[&[&str]] = match kind {
        "clingen_erepo" => &[&["clinvar_variation_id"], &["applied_evidence_codes_met"]],
        "clingen_dosage" => &[
            &["gene_symbol", "symbol"],
            &["haploinsufficiency_score", "hi_score"],
            &["triplosensitivity_score", "ts_score"],
        ],
        "clingen_dosage_runtime" => &[
            &["gene_symbol"],
            &["haploinsufficiency_score"],
            &["triplosensitivity_score"],
        ],
        "g2p_ddg2p" => &[
            &["gene_symbol"],
            &["disease_name", "disease_label", "disease"],
            &["molecular_mechanism", "mutation_consequence"],
        ],
        _ => unreachable!(),
    };
    let normalized = header.iter().map(|v| normalize(v)).collect::<Vec<_>>();
    require(
        groups
            .iter()
            .all(|group| normalized.iter().any(|v| group.contains(&v.as_str()))),
        format!(
            "Downloaded {kind} has incompatible source columns: {}",
            path.display()
        ),
    )?;
    let numeric = match kind {
        "clingen_erepo" => Some("ClinVar Variation Id"),
        "clingen_dosage_runtime" => Some("Haploinsufficiency Score"),
        _ => None,
    };
    let missing = [
        "", "NaN", "nan", "NA", "N/A", "NULL", "null", "None", "<NA>", "#N/A", "#N/A N/A", "#NA",
        "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "n/a",
    ];
    let mut count = 0;
    while let Some(row) = rows.next()? {
        if row.iter().all(String::is_empty) {
            continue;
        }
        require(
            row.len() == header.len(),
            format!(
                "Downloaded {kind} contains an incomplete source row: {}",
                path.display()
            ),
        )?;
        if let Some(name) = numeric {
            let value = crate::python_text::python_strip(
                &row[header.iter().position(|v| v == name).unwrap()],
            );
            require(
                missing.contains(&value) || finite_number(value),
                format!(
                    "Downloaded {kind} has a nonnumeric {name}: {}",
                    path.display()
                ),
            )?;
        }
        count += 1;
    }
    require(
        count > 0,
        format!(
            "Downloaded {kind} contains no source rows: {}",
            path.display()
        ),
    )
}
fn read_receipt(path: &Path) -> Result<Option<Json>> {
    let bytes = match fs::read(path) {
        Ok(v) => v,
        Err(_) => return Ok(None),
    };
    let text = match crate::python_io::decode_utf8(&bytes) {
        Ok(v) => v,
        Err(_) => return Ok(None),
    };
    match json_text::parse(&text.replace("\r\n", "\n").replace('\r', "\n"), false) {
        Ok(value) => Ok(Some(value)),
        Err(message) if json_text::is_refusal(&message) => Err(error(
            "RuntimeError",
            format!("Unsupported JSON profile: {message}"),
        )),
        Err(_) => Ok(None),
    }
}
fn publish(path: &Path, value: &Json) -> Result<()> {
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let name = path.file_name().unwrap_or_default().to_string_lossy();
    let template = tempfile_parent(parent)?.join(format!(".{name}.XXXXXX"));
    let mut bytes = CString::new(template.as_os_str().as_bytes())
        .map_err(|_| error("ValueError", "embedded null byte"))?
        .into_bytes_with_nul();
    // NamedTemporaryFile(delete=False) creates mode0600. A failed write or
    // rename leaves the sibling, so no cleanup guard or extra fsync is added.
    let fd = unsafe { libc::mkostemp(bytes.as_mut_ptr().cast(), libc::O_CLOEXEC) };
    let temporary = PathBuf::from(std::ffi::OsStr::from_bytes(&bytes[..bytes.len() - 1]));
    if fd < 0 {
        return Err(os_error(std::io::Error::last_os_error(), &temporary, None));
    }
    let mut stream = unsafe { File::from_raw_fd(fd) };
    let text = json_text::dumps_with_options(
        value,
        json_text::DumpsOptions {
            indent: None,
            sort_keys: true,
            ensure_ascii: true,
            compact: false,
            allow_nan: true,
        },
    )
    .map_err(invalid)?
        + "\n";
    stream
        .write_all(text.as_bytes())
        .map_err(|e| os_error(e, &temporary, None))?;
    drop(stream);
    fs::rename(&temporary, path).map_err(|e| os_error(e, &temporary, Some(path)))
}
fn source_receipt(
    path: &Path,
    url: &str,
    record: bool,
    receipt: Option<&Path>,
    metadata_only: bool,
) -> Result<bool> {
    let sidecar = receipt
        .map(Path::to_owned)
        .unwrap_or_else(|| PathBuf::from(format!("{}.source.json", path.display())));
    if metadata_only {
        let Some(previous) = read_receipt(&sidecar)? else {
            return Ok(false);
        };
        let valid = previous.get("kind").and_then(Json::as_str) == Some("priva_clingen_source_v1")
            && previous.get("url").and_then(Json::as_str) == Some(url)
            && previous
                .get("sha256")
                .and_then(Json::as_str)
                .is_some_and(|s| {
                    s.len() == 64
                        && s.bytes()
                            .all(|c| c.is_ascii_digit() || (b'a'..=b'f').contains(&c))
                });
        return Ok(valid);
    }
    let mut source = File::open(path).map_err(|e| os_error(e, path, None))?;
    let mut digest = crate::sha256::Sha256::new();
    let mut block = vec![0; 1024 * 1024];
    loop {
        let count = source
            .read(&mut block)
            .map_err(|e| os_error(e, path, None))?;
        if count == 0 {
            break;
        }
        digest.update(&block[..count]);
    }
    drop(source);
    let mut value = Json::Object(vec![
        ("kind".into(), Json::Str("priva_clingen_source_v1".into())),
        ("url".into(), Json::Str(url.into())),
        ("sha256".into(), Json::Str(digest.finish_hex())),
    ]);
    let previous = read_receipt(&sidecar)?.unwrap_or(Json::Null);
    if record {
        let download = if previous.get("url").and_then(Json::as_str) == Some(url) {
            previous
                .get("download_sha256")
                .unwrap_or(value.get("sha256").unwrap())
                .clone()
        } else {
            value.get("sha256").unwrap().clone()
        };
        value.set("download_sha256", download);
        publish(&sidecar, &value)?;
        Ok(true)
    } else {
        Ok(value.as_object().unwrap().iter().all(|(key, item)| {
            previous
                .get(key)
                .is_some_and(|old| json_text::py_equal(old, item))
        }))
    }
}
pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args.first().is_some_and(|v| v == "--help" || v == "-h") {
        println!("{HELP}");
        return Ok(());
    }
    let command = args.first().ok_or_else(|| invalid(HELP))?;
    let (names, flags, required): (&[&str], &[&str], &[&str]) = match command.as_str() {
        "validate" => (&["input", "kind"], &[], &["input", "kind"]),
        "record-source" => (&["input", "url"], &[], &["input", "url"]),
        "check-source" => (
            &["input", "url", "receipt"],
            &["metadata-only"],
            &["input", "url"],
        ),
        _ => return Err(invalid(HELP)),
    };
    let Some(occurrences) =
        parse_long_option_occurrences(&args[1..], names, flags).map_err(invalid)?
    else {
        println!("{HELP}");
        return Ok(());
    };
    for (name, value) in &occurrences {
        if name == "url" {
            url::check(value).map_err(|_| {
                invalid(format!(
                    "argument --url: invalid source_url value: {}",
                    crate::python_path::repr(value)
                ))
            })?;
        }
        if name == "kind" && !KINDS.contains(&value.as_str()) {
            return Err(invalid(format!(
                "argument --kind: invalid choice: {}",
                crate::python_path::repr(value)
            )));
        }
    }
    let options: HashMap<_, _> = occurrences.into_iter().collect();
    require_options(&options, required).map_err(invalid)?;
    let path = path_value(&options["input"]);
    crate::builder_error::finish(if command == "validate" {
        validate(&path, &options["kind"])
    } else {
        let receipt = options.get("receipt").map(|v| path_value(v));
        source_receipt(
            &path,
            &options["url"],
            command == "record-source",
            receipt.as_deref(),
            options.contains_key("metadata-only"),
        )
        .map(|ok| {
            if !ok {
                std::process::exit(1);
            }
        })
    })
}
