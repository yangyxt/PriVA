//! `priva biomart fetch`: biomart_nm_to_enst.py `fetch` and the acquisition and
//! `install` functions of biomart_ftp_nm_to_enst.py. Requests run through curl
//! (see http.rs); refusals, printed lines, the directory flock, the bounded
//! attempts and the receipts are the scripts'. Receipts name this program as
//! their producer, with the `biomart-fetch` builder identity in place of the
//! script's SHA-256.
use super::http::{self, Response};
use super::publication_json::Record;
use super::*;
use crate::clinvar_contract::resolve_non_strict;
use crate::file_lock::{ExclusiveLock, Wait};
use std::collections::HashMap;
use std::io::{self, Write};
use std::thread::sleep;
use std::time::{Duration, Instant};

const MART: &str = "ENSEMBL_MART_ENSEMBL";
const NATIVE_PRODUCER: &str = "priva biomart fetch";
const USER_AGENT: &str = "PriVA-pinned-BioMart/1";

pub(super) struct Options {
    source: String,
    output: PathBuf,
    raw_output: Option<PathBuf>,
    resource_dir: Option<PathBuf>,
    timeout: f64,
    attempts: i64,
}

impl Options {
    /// The script's argparse admission, which exits 2 on failure.
    pub(super) fn admit(values: &HashMap<String, String>) -> Result<Self> {
        let option = |key: &str| values.get(key).map(String::as_str);
        let selected = option("source").unwrap();
        let sources = array(get(facts(), "sources"));
        if !sources.iter().any(|s| get(s, "source_id").as_str() == Some(selected)) {
            return Err(invalid(format!(
                "argument --source: invalid choice: {} (choose from {})",
                repr(selected),
                sources
                    .iter()
                    .map(|s| repr(text(get(s, "source_id"))))
                    .collect::<Vec<_>>()
                    .join(", ")
            )));
        }
        let timeout = option("timeout").unwrap_or("600");
        let timeout = timeout.parse().map_err(|_| {
            invalid(format!(
                "argument --timeout: invalid float value: {}",
                repr(timeout)
            ))
        })?;
        let attempts = option("attempts").unwrap_or("3");
        let attempts = attempts.parse().map_err(|_| {
            invalid(format!(
                "argument --attempts: invalid int value: {}",
                repr(attempts)
            ))
        })?;
        if selected != FTP && option("raw-output").is_none() {
            return Err(invalid("--raw-output is required for online BioMart sources"));
        }
        Ok(Self {
            source: selected.to_owned(),
            output: path_value(option("output").unwrap()),
            raw_output: option("raw-output").map(path_value),
            resource_dir: option("resource-dir").map(path_value),
            timeout,
            attempts,
        })
    }
}

fn path_value(value: &str) -> PathBuf {
    crate::python_path::path_value(value)
}

fn producer_sha256() -> Result<String> {
    crate::builder_identity::identity("biomart-fetch")
        .ok_or_else(|| error("ValueError", "No biomart-fetch builder identity"))
}

fn slurm_job() -> Option<String> {
    std::env::var("SLURM_JOB_ID").ok()
}

fn lexists(path: &Path) -> bool {
    fs::symlink_metadata(path).is_ok()
}

/// `Path.parent` of an admitted path; a bare file name lives in `.`.
fn parent(path: &Path) -> PathBuf {
    match path.parent() {
        Some(folder) if !folder.as_os_str().is_empty() => folder.to_path_buf(),
        Some(_) => PathBuf::from("."),
        None => path.to_path_buf(),
    }
}

fn make_directory(path: &Path) -> Result<()> {
    fs::create_dir_all(path).map_err(|e| os_error(e, path, None))
}

fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    PathBuf::from(format!("{}{suffix}", path.display()))
}

/// `flock(LOCK_EX | LOCK_NB)` on the output directory's own inode: there is no
/// lock file, and writers in one directory cannot race on its `.part` paths.
fn lock_directory(folder: &Path) -> Result<io::Result<ExclusiveLock>> {
    let handle = File::open(folder).map_err(|e| os_error(e, folder, None))?;
    Ok(ExclusiveLock::acquire(handle, Wait::Nonblocking))
}

pub(super) fn acquire(context: &Context, options: &Options) -> Result<()> {
    require(
        (1..=5).contains(&options.attempts) && !(options.timeout <= 0.0),
        "Use 1-5 bounded attempts and a positive timeout",
    )?;
    let output = &options.output;
    let pinned_name = source_for_basename(&name(output)).is_some();
    let selected = source(&string(&options.source))?;
    if options.source == FTP {
        require(
            !pinned_name || get(selected, "basename").as_str() == Some(&name(output)),
            "FTP116 source conflicts with the selected pinned mapping filename",
        )?;
        require(
            options.resource_dir.is_some() || lexists(output) || lexists(&provenance_path(output)),
            "New FTP acquisition requires --resource-dir outside the checkout",
        )?;
        if let Some(folder) = options.resource_dir.as_ref().filter(|_| !lexists(output)) {
            publication::external_path(context, folder, false)?;
        }
        return install_ftp(context, options, selected);
    }
    let raw_output = options
        .raw_output
        .as_ref()
        .ok_or_else(|| error("ValueError", "Online BioMart acquisition requires its raw_output destination"))?;
    require(
        options.resource_dir.is_none(),
        "resource_dir is specific to the FTP source",
    )?;
    if !lexists(output) {
        publication::external_path(context, raw_output, false)?;
    }
    let folder = parent(output);
    make_directory(&folder)?;
    let _lock = lock_directory(&folder)?.map_err(|cause| {
        if cause.kind() == io::ErrorKind::WouldBlock {
            error(
                "ValueError",
                format!(
                    "Another mapping writer holds {}; retry after it finishes",
                    folder.display()
                ),
            )
        } else {
            os_error(cause, &folder, None)
        }
    })?;
    acquire_service(context, options, selected, raw_output)
}

/// `urlencode(parameters)`.
fn urlencode(parameters: &[(&str, &str)]) -> String {
    let quote = |value: &str| {
        let mut out = String::new();
        for byte in value.bytes() {
            match byte {
                b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'_' | b'.' | b'-' | b'~' => {
                    out.push(char::from(byte));
                }
                b' ' => out.push('+'),
                other => out.push_str(&format!("%{other:02X}")),
            }
        }
        out
    };
    parameters
        .iter()
        .map(|(key, value)| format!("{}={}", quote(key), quote(value)))
        .collect::<Vec<_>>()
        .join("&")
}

/// `_request`: the whole response stays on the selected publisher's host.
fn service_request(
    selected: &Json,
    parameters: &[(&str, &str)],
    timeout: f64,
) -> Result<Response> {
    let endpoint = text(get(selected, "endpoint"));
    let url = format!("{endpoint}?{}", urlencode(parameters));
    let headers = [("User-Agent", USER_AGENT), ("Accept-Encoding", "identity")];
    let response = http::open(&url, &headers, timeout)?;
    require(
        http::hostname(&response.url) == http::hostname(endpoint),
        format!(
            "BioMart redirected away from the explicitly selected {} host",
            text(get(selected, "assembly"))
        ),
    )?;
    Ok(response)
}

/// `response.read(1024 * 1024).decode("utf-8")`.
fn read_text(mut response: Response) -> Result<String> {
    let block = response.read_block()?;
    if block.len() < http::BLOCK {
        response.finish()?;
    }
    crate::python_io::decode_utf8(&block).map_err(|e| error("UnicodeDecodeError", e))
}

/// `source_metadata`: the exact mart and dataset, publisher text kept verbatim.
fn source_metadata(selected: &Json, timeout: f64, attempts: i64) -> Result<Json> {
    let read_both = || -> Result<(String, String)> {
        let registry = service_request(
            selected,
            &[("type", "registry"), ("virtualSchema", "default")],
            timeout,
        )
        .and_then(read_text)?;
        let datasets =
            service_request(selected, &[("type", "datasets"), ("mart", MART)], timeout)
                .and_then(read_text)?;
        Ok((registry, datasets))
    };
    let mut attempt = 1;
    let (registry, datasets) = loop {
        match read_both() {
            Ok(pair) => break pair,
            Err(failure) => {
                let caught = http::is_os_error(failure.as_ref())
                    || http::is_value_error(failure.as_ref());
                if !caught || attempt == attempts {
                    return Err(failure);
                }
                println!(
                    "Metadata request {attempt}/{attempts} failed: {}; retrying the same source",
                    http::message(failure.as_ref())
                );
                sleep(Duration::from_secs((5 * attempt).min(15) as u64));
                attempt += 1;
            }
        }
    };
    verify::parse_metadata(&registry, &datasets, get(selected, "source_id"))
}

/// One export attempt, writing the raw response and then the checked table.
fn export(
    selected: &Json,
    timeout: f64,
    raw_partial: &Path,
    table_partial: &Path,
) -> Result<Json> {
    let query = text(get(facts(), "query_xml"));
    let mut response = service_request(selected, &[("query", query)], timeout)?;
    let mut handle = File::create(raw_partial).map_err(|e| os_error(e, raw_partial, None))?;
    loop {
        let block = response.read_block()?;
        if block.is_empty() {
            break;
        }
        handle
            .write_all(&block)
            .map_err(|e| os_error(e, raw_partial, None))?;
    }
    response.finish()?;
    drop(handle);
    table::validate(raw_partial, true, Some(table_partial))
}

fn identical_source(before: &Json, after: &Json) -> bool {
    ["registry_record", "dataset_record"]
        .iter()
        .all(|key| equal(get(before, key), get(after, key)))
}

fn slurm_value() -> Json {
    slurm_job().map_or(Json::Null, string)
}

/// `_acquire_locked`.
fn acquire_service(
    context: &Context,
    options: &Options,
    selected: &Json,
    raw_output: &Path,
) -> Result<()> {
    let output = &options.output;
    let source_id = &options.source;
    let assembly = text(get(selected, "assembly"));
    require(
        source_for_basename(&name(output)).is_none()
            || get(selected, "basename").as_str() == Some(&name(output)),
        format!(
            "Selected BioMart source conflicts with the requested pinned output filename: {}",
            output.display()
        ),
    )?;
    if exists(output)? || is_link(output)? {
        let record = verify::pinned(output, Some(assembly))?;
        require(
            get(&record, "source_id").as_str() == Some(source_id),
            format!(
                "Existing mapping is not the explicitly requested BioMart source: {}",
                output.display()
            ),
        )?;
        println!(
            "Reusing completed pinned mapping: {} ({} bytes, SHA256 {})",
            output.display(),
            py_str(get(&record, "bytes")),
            py_str(get(&record, "sha256"))
        );
        return Ok(());
    }
    let receipt = provenance_path(output);
    require(
        !exists(&receipt)? && !is_link(&receipt)?,
        format!(
            "Mapping is absent but its provenance receipt exists: {}; retain it and recover explicitly",
            receipt.display()
        ),
    )?;
    require(
        !exists(raw_output)? && !is_link(raw_output)?,
        format!(
            "A published raw export already exists: {}; retain it and select a new output path",
            raw_output.display()
        ),
    )?;
    make_directory(&parent(output))?;
    make_directory(&parent(raw_output))?;
    require(
        raw_output != output.as_path(),
        "The completion-stamped raw response and installed TSV must use different paths",
    )?;
    println!(
        "Source {source_id}; Slurm job {}; receipt {}",
        slurm_job().unwrap_or_else(|| "none".to_owned()),
        receipt.display()
    );
    let producer = producer_sha256()?;
    println!("Producer {NATIVE_PRODUCER}; SHA256 {producer}");
    let metadata_timeout = options.timeout.min(60.0);
    let before = source_metadata(selected, metadata_timeout, options.attempts)?;
    let query = get(facts(), "query_xml");
    identity::atomic_json(
        &with_suffix(raw_output, ".metadata-before.json"),
        &object([
            ("source_id", string(source_id)),
            ("query_xml", query.clone()),
            ("metadata", before.clone()),
            ("producer_script_sha256", string(&producer)),
            ("producer_script_path", string(NATIVE_PRODUCER)),
            ("slurm_job_id", slurm_value()),
            ("retrieval_started_utc", string(publication::utc_now()?)),
        ]),
    )?;
    let raw_partial = with_suffix(raw_output, ".part");
    let table_partial = with_suffix(output, ".part");
    require(
        !is_link(&raw_partial)? && !is_link(&table_partial)?,
        "Refusing to write through a BioMart staging-file symlink",
    )?;
    let mut attempt = 1;
    let statistics = loop {
        println!(
            "Requesting {source_id} whole-dataset export, attempt {attempt}/{}",
            options.attempts
        );
        match export(selected, options.timeout, &raw_partial, &table_partial) {
            Ok(statistics) => break statistics,
            Err(failure) => {
                let caught = http::is_os_error(failure.as_ref())
                    || http::is_incomplete_read(failure.as_ref());
                if !caught || attempt == options.attempts {
                    return Err(failure);
                }
                println!(
                    "Export failed: {}; retrying the same pinned source",
                    http::message(failure.as_ref())
                );
                sleep(Duration::from_secs((10 * attempt).min(30) as u64));
                attempt += 1;
            }
        }
    };
    let after = source_metadata(selected, metadata_timeout, options.attempts)?;
    require(
        identical_source(&before, &after),
        "BioMart source metadata changed during export; no final mapping installed",
    )?;
    let mut record = object([
        ("source_id", string(source_id)),
        ("endpoint", get(selected, "endpoint").clone()),
        ("mart", string(MART)),
        ("mart_release", get(selected, "mart_release").clone()),
        (
            "annotation_base_release",
            get(selected, "annotation_base_release").clone(),
        ),
        ("assembly", get(selected, "assembly").clone()),
        ("dataset", string(DATASET)),
        ("attributes", get(facts(), "attributes").clone()),
        ("header", get(facts(), "headers").clone()),
        ("query_xml", query.clone()),
        ("query_sha256", string(hash(text(query).as_bytes()))),
        ("retrieved_utc", string(publication::utc_now()?)),
        ("producer_script_sha256", string(&producer)),
        ("producer_script_path", string(NATIVE_PRODUCER)),
        ("slurm_job_id", slurm_value()),
        ("metadata_before", before),
        ("metadata_after", after),
        ("completion_stamp_verified", Json::Bool(true)),
        ("statistics", statistics),
    ]);
    let identity = snapshot_identity(&table_partial, get(selected, "source_id"))?;
    record.set("bytes", get(&identity, "bytes").clone());
    record.set("sha256", get(&identity, "sha256").clone());
    record.set("raw_bytes", Json::Int(metadata(&raw_partial)?.len() as i128));
    record.set("raw_sha256", string(sha256(&raw_partial)?));
    fs::rename(&raw_partial, raw_output)
        .map_err(|e| os_error(e, &raw_partial, Some(raw_output)))?;
    fs::rename(&table_partial, output).map_err(|e| os_error(e, &table_partial, Some(output)))?;
    if identity::is_bundled(context, output)? {
        // The full record stays beside the external raw export; source data
        // keeps portable publisher and method facts only.
        let full = json_text::dumps(&record, true).map_err(|e| error("ValueError", e))?;
        publication::runtime(
            context,
            &Record::parse(&full)?,
            &parent(raw_output),
            &name(output),
            false,
        )?;
        record = bundled_service_record(record);
    }
    identity::atomic_json(&receipt, &record)?;
    verify::pinned(output, Some(assembly))?;
    println!(
        "Completed export; provenance receipt closed: {}",
        receipt.display()
    );
    let summary = object(["source_id", "bytes", "sha256", "statistics"].map(|key| {
        (key, get(&record, key).clone())
    }));
    println!(
        "{}",
        json_text::dumps(&summary, true).map_err(|e| error("ValueError", e))?
    );
    Ok(())
}

fn bundled_service_record(record: Json) -> Json {
    let Json::Object(fields) = record else {
        unreachable!("the service record is an object");
    };
    let mut bundled = Json::Object(
        fields
            .into_iter()
            .filter(|(key, _)| {
                !matches!(
                    key.as_str(),
                    "retrieved_utc" | "producer_script_path" | "slurm_job_id"
                )
            })
            .collect(),
    );
    for (key, value) in [
        ("metadata_kind", string("priva_bundled_mapping_snapshot")),
        ("metadata_schema_version", Json::Int(1)),
        ("acquisition_method", string("online_biomart")),
        ("order_policy", string("Preserve BioMart response row order.")),
        (
            "multiplicity_policy",
            string("uniqueRows=0; retain duplicate rows and blank RefSeq fields."),
        ),
        ("scientific_qualified", Json::Bool(false)),
        ("production_selected", Json::Bool(false)),
        ("original_biomart_export_identity", Json::Bool(false)),
        ("original_biomart_response_order_identity", Json::Bool(false)),
        ("online_query_executed", Json::Bool(true)),
    ] {
        bundled.set(key, value);
    }
    bundled
}

/// `int(response.headers["Content-Length"])`.
fn content_length(response: &Response) -> Result<i128> {
    let value = response.header("Content-Length").ok_or_else(|| {
        error(
            "TypeError",
            "int() argument must be a string, a bytes-like object or a real number, not 'NoneType'",
        )
    })?;
    value.parse().map_err(|_| {
        error(
            "ValueError",
            format!("invalid literal for int() with base 10: {}", repr(value)),
        )
    })
}

/// The FTP script's `acquire`: one pinned dump, resumed from its `.part`.
fn acquire_raw(raw: &Path, file: &str, base_url: &str, timeout: f64) -> Result<()> {
    let path = raw.join(file);
    let pin = ftp::named_pin(file)?;
    let total = get_bytes(&pin);
    if exists(&path)? || is_link(&path)? {
        let record = ftp::verify_raw(&path, &pin)?;
        println!("REUSE {file} {}", py_str(get(&record, "bytes")));
        return Ok(());
    }
    let part = raw.join(format!("{file}.part"));
    require(
        !is_link(&part)?,
        format!("Refuse a symlinked transfer partial: {}", part.display()),
    )?;
    let present = exists(&part)?;
    let offset = if present {
        metadata(&part)?.len() as i128
    } else {
        0
    };
    require(
        offset <= total,
        format!("Unexpected partial length: {}", part.display()),
    )?;
    if offset < total {
        let url = format!("{base_url}{file}");
        let range = format!("bytes={offset}-");
        let mut headers = vec![("User-Agent", USER_AGENT), ("Accept-Encoding", "identity")];
        if offset > 0 {
            headers.push(("Range", range.as_str()));
        }
        let mut response = http::open(&url, &headers, timeout)?;
        require(
            response.url == url,
            format!("Unexpected download redirect: {}", response.url),
        )?;
        if offset > 0 {
            let expected = format!("bytes {offset}-{}/{total}", total - 1);
            require(
                response.status == 206
                    && response.header("Content-Range") == Some(expected.as_str()),
                "Server did not honor exact resume range; partial retained",
            )?;
        } else {
            require(response.status == 200, "Full download must return HTTP200")?;
        }
        require(
            content_length(&response)? == total - offset,
            "Publisher transfer length changed",
        )?;
        let mut stream = if present {
            fs::OpenOptions::new().append(true).open(&part)
        } else {
            fs::OpenOptions::new().write(true).create_new(true).open(&part)
        }
        .map_err(|e| os_error(e, &part, None))?;
        let mut received = offset;
        let mut reported = Instant::now();
        loop {
            let block = response.read_block()?;
            if block.is_empty() {
                break;
            }
            received += block.len() as i128;
            require(received <= total, "Download exceeded pinned size")?;
            stream
                .write_all(&block)
                .map_err(|e| os_error(e, &part, None))?;
            if reported.elapsed() >= Duration::from_secs(20) {
                println!("DOWNLOAD {file} {received} / {total}");
                reported = Instant::now();
            }
        }
        response.finish()?;
    }
    let record = ftp::verify_raw(&part, &pin)?;
    fs::rename(&part, &path).map_err(|e| os_error(e, &part, Some(&path)))?;
    println!(
        "ACQUIRED {file} {} {}",
        py_str(get(&record, "bytes")),
        py_str(get(&record, "sha256"))
    );
    Ok(())
}

fn get_bytes(pin: &Json) -> i128 {
    match array(pin).first() {
        Some(Json::Int(bytes)) => *bytes,
        _ => unreachable!("the fixed raw pins hold a byte count"),
    }
}

/// Each publisher BSD sum in CHECKSUMS must equal its fixed pin.
fn check_checksums(raw: &Path, sources: &[Json]) -> Result<()> {
    let path = raw.join("CHECKSUMS");
    let bytes = fs::read(&path).map_err(|e| os_error(e, &path, None))?;
    let document = crate::python_io::decode_utf8(&bytes).map_err(|e| error("UnicodeDecodeError", e))?;
    let mut found: Vec<(String, Json)> = Vec::new();
    for line in document.lines() {
        let cells = line.split_whitespace().collect::<Vec<_>>();
        if cells.len() != 3 || !sources.iter().any(|s| get(s, "name").as_str() == Some(cells[2])) {
            continue;
        }
        require(
            !found.iter().any(|(file, _)| file == cells[2]),
            "Repeated publisher checksum entry",
        )?;
        let number = |cell: &str| -> Result<Json> {
            cell.parse().map(Json::Int).map_err(|_| {
                error(
                    "ValueError",
                    format!("invalid literal for int() with base 10: {}", repr(cell)),
                )
            })
        };
        found.push((
            cells[2].to_owned(),
            Json::Array(vec![number(cells[0])?, number(cells[1])?]),
        ));
    }
    let agrees = sources.iter().all(|pinned| {
        let expected = get(pinned, "publisher_bsd_sum_blocks");
        !truth(expected)
            || found
                .iter()
                .any(|(file, sum)| get(pinned, "name").as_str() == Some(file) && equal(sum, expected))
    });
    require(agrees, "Publisher checksum document differs from fixed FTP pins")
}

fn verify_ftp(output: &Path, assembly: &str) -> Result<()> {
    let record = read_json(&provenance_path(output))?;
    verify::ftp_record(output, record, Some(assembly)).map(|_| ())
}

fn publish_ftp_receipt(path: &Path, record: &Record) -> Result<()> {
    publication::atomic(path, record, true)
}

fn record_of(value: &Json) -> Result<Record> {
    Record::parse(&json_text::dumps(value, true).map_err(|e| error("ValueError", e))?)
}

/// `Path(os.path.abspath(output)) == canonical`, on the unlinked canonical path.
fn is_canonical(context: &Context, output: &Path, selected: &Json) -> Result<bool> {
    let script = context.helper_root.join("scripts/biomart_ftp_nm_to_enst.py");
    let script = resolve_non_strict(&script)?;
    let root = script
        .parent()
        .and_then(Path::parent)
        .ok_or_else(|| error("IndexError", "1"))?;
    let canonical = root.join(text(get(selected, "relative_path")));
    Ok(crate::python_path::tempfile_parent(output)? == canonical
        && !is_link(output)?
        && resolve_non_strict(output)? == canonical)
}

/// `FTP.install`: reuse exact known bytes offline, or acquire the pinned dumps.
fn install_ftp(context: &Context, options: &Options, selected: &Json) -> Result<()> {
    let output = &options.output;
    let assembly = text(get(selected, "assembly"));
    let receipt = provenance_path(output);
    let folder = parent(output);
    make_directory(&folder)?;
    let _lock = lock_directory(&folder)?.map_err(|cause| {
        if cause.kind() == io::ErrorKind::WouldBlock {
            error("BlockingIOError", "[Errno 11] Resource temporarily unavailable")
        } else {
            os_error(cause, &folder, None)
        }
    })?;
    let producer = producer_sha256()?;
    if exists(output)? || is_link(output)? {
        if exists(&receipt)? || is_link(&receipt)? {
            return verify_ftp(output, assembly);
        }
        // The immutable SHA/bytes/statistics pin identifies this reviewed
        // snapshot. Do not invent a new acquisition or service stamp.
        let record = publication::make_provenance(
            output,
            "registered_existing_exact_snapshot",
            NATIVE_PRODUCER,
            &producer,
        )?;
        publish_ftp_receipt(&receipt, &record_of(&record)?)?;
        return verify_ftp(output, assembly);
    }
    require(
        !exists(&receipt)? && !is_link(&receipt)?,
        "Missing FTP table has existing provenance; restore it explicitly",
    )?;
    let root = options
        .resource_dir
        .as_ref()
        .ok_or_else(|| error("ValueError", "New FTP acquisition requires an external resource_dir"))?;
    publication::external_path(context, root, true)?;
    let raw = root.join("raw");
    make_directory(&raw)?;
    let specification = get(facts(), "ftp_specification");
    let base_url = text(get(specification, "source_root"));
    let sources = array(get(specification, "raw_sources"));
    for pinned in sources {
        let file = text(get(pinned, "name"));
        for attempt in 0..options.attempts {
            match acquire_raw(&raw, file, base_url, options.timeout) {
                Ok(()) => break,
                Err(failure) => {
                    if !http::is_os_error(failure.as_ref()) || attempt + 1 == options.attempts {
                        return Err(failure);
                    }
                    sleep(Duration::from_secs((5 * (attempt + 1)).min(15) as u64));
                }
            }
        }
    }
    check_checksums(&raw, sources)?;
    let candidate = with_suffix(output, ".candidate");
    ftp::project_from_dumps(&raw, root, &candidate)?;
    let value = publication::make_provenance(
        &candidate,
        "projected_from_verified_ftp_sources",
        NATIVE_PRODUCER,
        &producer,
    )?;
    let mut record = record_of(&value)?;
    if is_canonical(context, output, selected)? {
        publication::runtime(context, &record, root, &name(output), true)?;
        record = record.portable()?;
    }
    // Same-parent hard linking is atomic and refuses replacement. Preserve
    // staged data and source dumps if verification or publication fails.
    fs::hard_link(&candidate, output).map_err(|e| os_error(e, &candidate, Some(output)))?;
    fs::remove_file(&candidate).map_err(|e| os_error(e, &candidate, None))?;
    publish_ftp_receipt(&receipt, &record)?;
    verify_ftp(output, assembly)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn urlencode_matches_quote_plus() {
        assert_eq!(
            urlencode(&[("type", "a b"), ("q", "<x y=\"1\"/>~-_.")]),
            "type=a+b&q=%3Cx+y%3D%221%22%2F%3E~-_."
        );
    }

    #[test]
    fn a_bare_file_name_is_in_the_current_folder() {
        assert_eq!(parent(Path::new("x.tsv")), PathBuf::from("."));
        assert_eq!(parent(Path::new("a/x.tsv")), PathBuf::from("a"));
    }

    #[test]
    fn the_bundled_record_drops_run_history_and_states_its_kind() {
        let record = object([
            ("retrieved_utc", string("t")),
            ("producer_script_path", string("p")),
            ("slurm_job_id", Json::Null),
            ("producer_script_sha256", string("s")),
        ]);
        let bundled = bundled_service_record(record);
        assert!(bundled.get("retrieved_utc").is_none());
        assert!(bundled.get("producer_script_path").is_none());
        assert!(bundled.get("slurm_job_id").is_none());
        assert_eq!(
            bundled.get("producer_script_sha256").unwrap().as_str(),
            Some("s")
        );
        assert_eq!(
            bundled.get("acquisition_method").unwrap().as_str(),
            Some("online_biomart")
        );
    }

    #[test]
    fn the_native_producer_identity_is_a_digest() {
        assert!(hex64(&string(producer_sha256().unwrap())));
    }
}
