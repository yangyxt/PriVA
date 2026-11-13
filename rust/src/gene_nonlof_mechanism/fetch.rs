//! Source acquisition, freshness receipts and the manifest TSV.
use super::*;
use crate::pandas_table::{Column, ObjectValue, Table};
use std::process::{Command, Stdio};

/// `previous_meta.get(source.name, {})`, emptied when its raw path moved.
fn previous(spec: &Spec, previous_sources: &Json, raw_path: &Path) -> Result<Json> {
    let old = get_or(previous_sources, spec.name, Json::Object(vec![]))?;
    if !get(&old, "raw_path")?
        .is_some_and(|value| equal(value, &text(raw_path.display().to_string())))
    {
        return Ok(Json::Object(vec![]));
    }
    Ok(old)
}

/// The script's `source_identity_matches`.
fn identity_matches(spec: &Spec, metadata: &Json, raw_path: Option<&Path>) -> Result<bool> {
    let previous = get(metadata, "url")?.cloned().unwrap_or(Json::Null);
    let mut matches = equal(&previous, &text(spec.url.clone()))
        || (!truthy(&previous) && spec.url == spec.default_url());
    if let (true, Some(raw_path)) = (matches, raw_path) {
        if let Some(sha) = get(metadata, "sha256")?.filter(|v| truthy(v)) {
            matches =
                crate::python_path::is_file(raw_path)? && equal(&text(sha256_file(raw_path)?), sha);
        }
    }
    Ok(matches)
}

/// The script's `source_is_fresh`.
fn is_fresh(path: &Path, meta: &Json, refresh_days: i64, force: bool) -> Result<bool> {
    if force || !nonempty(path) {
        return Ok(false);
    }
    let checked = match get(meta, "last_checked_at_utc")? {
        Some(value) if truthy(value) => match value {
            Json::Str(stamp) => parse_iso(stamp),
            other => {
                return Err(error(
                    "TypeError",
                    format!(
                        "fromisoformat: argument must be str, not {}",
                        type_name(other)
                    ),
                ));
            }
        },
        _ => None,
    };
    let checked = match checked {
        Some((_, false)) => {
            return Err(error(
                "TypeError",
                "can't subtract offset-naive and offset-aware datetimes",
            ));
        }
        Some((micros, true)) => micros,
        None => mtime_micros(&stat(path)?),
    };
    Ok(age_seconds(checked) / 86_400.0 < refresh_days as f64)
}

/// The script's `update_metadata_for_file`.
fn metadata_for_file(
    spec: &Spec,
    raw_path: &Path,
    status: &str,
    old_sha256: Option<&Json>,
    failure: Option<String>,
) -> Result<Json> {
    let mut metadata = object([
        ("name", text(spec.name)),
        ("url", text(spec.url.clone())),
        ("docs_url", text(spec.docs_url)),
        ("raw_path", text(raw_path.display().to_string())),
        ("source_role", text(spec.source_role)),
        ("official_update_note", text(spec.official_update_note)),
        ("refresh_days", Json::Int(spec.refresh_days.into())),
        ("parser", text(spec.parser)),
        ("last_checked_at_utc", text(iso_now())),
        ("status", text(status)),
        ("error", text(failure.unwrap_or_default())),
    ]);
    if nonempty(raw_path) {
        let sha256 = sha256_file(raw_path)?;
        let changed =
            old_sha256.is_some_and(|old| truthy(old) && !equal(old, &text(sha256.clone())));
        metadata.set("size_bytes", Json::Int(stat(raw_path)?.len().into()));
        metadata.set("sha256", text(sha256));
        metadata.set("md5", text(md5_file(raw_path)?));
        metadata.set("file_mtime_utc", text(file_mtime_iso(raw_path)?));
        metadata.set("content_changed", Json::Bool(changed));
        if matches!(status, "downloaded" | "updated") {
            metadata.set("last_downloaded_at_utc", text(iso_now()));
        }
    }
    Ok(metadata)
}

fn with_name(path: &Path, suffix: &str) -> PathBuf {
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    path.with_file_name(format!("{name}{suffix}"))
}

fn mkdir_parent(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent).map_err(|e| os_error(e, parent, None))?;
    }
    Ok(())
}

/// The script's `download_url_with_curl`, the same curl command.
pub(super) fn curl(
    url: &str,
    out: &Path,
    transfer: &Transfer,
    max_time: Option<i128>,
    resume: bool,
) -> Result<()> {
    mkdir_parent(out)?;
    let timeout = transfer.timeout;
    let mut command = Command::new("curl");
    command.args([
        "-L",
        "--fail",
        "--show-error",
        "--silent",
        "--connect-timeout",
        &timeout.to_string(),
        "--max-time",
        &max_time
            .unwrap_or_else(|| (timeout * 4).max(timeout + 60))
            .to_string(),
        "--retry",
        &(transfer.retries - 1).max(0).to_string(),
        "--user-agent",
        USER_AGENT,
        "--output",
    ]);
    command.arg(out);
    if resume {
        command.args(["--continue-at", "-"]);
    }
    if !transfer.proxy_url.is_empty() {
        command.args(["--proxy", &transfer.proxy_url]);
    }
    command.arg(url);
    let output = command
        .stdin(Stdio::null())
        .output()
        .map_err(|e| os_error(e, Path::new("curl"), None))?;
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stdout = String::from_utf8_lossy(&output.stdout);
        let detail = if stderr.is_empty() { stdout } else { stderr };
        return Err(error(
            "RuntimeError",
            format!(
                "curl failed for {url}: {}",
                crate::python_text::python_strip(&detail)
            ),
        ));
    }
    Ok(())
}

/// One `urlopen` request written to `out`: `Ok(())`, or `str(error)` of the
/// exception urllib would raise.
fn urlopen(
    url: &str,
    out: &Path,
    transfer: &Transfer,
    accept_json: bool,
) -> std::result::Result<(), String> {
    let headers = with_name(out, ".headers");
    let timeout = transfer.timeout.to_string();
    let mut command = Command::new("curl");
    command.args([
        "--http1.1",
        "--silent",
        "--show-error",
        "-L",
        "--max-redirs",
        "10",
        "--connect-timeout",
        &timeout,
        "--speed-limit",
        "1",
        "--speed-time",
        &timeout,
        "--user-agent",
        USER_AGENT,
        "-H",
        "Accept-Encoding: identity",
    ]);
    if accept_json {
        command.args(["-H", "Accept: application/json"]);
    }
    if !transfer.proxy_url.is_empty() {
        command.args(["--proxy", &transfer.proxy_url]);
    }
    command
        .arg("--dump-header")
        .arg(&headers)
        .arg("--output")
        .arg(out)
        .arg(url);
    let output = command
        .stdin(Stdio::null())
        .output()
        .map_err(|e| format!("[Errno {}] {e}", e.raw_os_error().unwrap_or(0)))?;
    let header_text = fs::read_to_string(&headers).unwrap_or_default();
    let _ = fs::remove_file(&headers);
    let status_line = header_text
        .lines()
        .filter(|line| line.starts_with("HTTP/"))
        .last()
        .unwrap_or("")
        .to_owned();
    if !output.status.success() {
        let _ = fs::remove_file(out);
        let reason = String::from_utf8_lossy(&output.stderr).trim().to_owned();
        return Err(match output.status.code() {
            Some(6) => "<urlopen error [Errno -2] Name or service not known>".into(),
            Some(7) => "<urlopen error [Errno 111] Connection refused>".into(),
            Some(28) => "<urlopen error timed out>".into(),
            _ => format!("<urlopen error {reason}>"),
        });
    }
    let mut parts = status_line.splitn(3, ' ');
    parts.next();
    let code = parts.next().unwrap_or("").to_owned();
    if !code.starts_with('2') {
        let _ = fs::remove_file(out);
        let reason = parts.next().unwrap_or("").trim_end().to_owned();
        return Err(format!("HTTP Error {code}: {reason}"));
    }
    Ok(())
}

fn needs_curl(transfer: &Transfer) -> bool {
    transfer.download_tool == "curl"
        || (transfer.download_tool == "auto"
            && ["socks4://", "socks4h://", "socks5://", "socks5h://"]
                .iter()
                .any(|scheme| transfer.proxy_url.to_lowercase().starts_with(scheme)))
}

fn pause(attempt: i128) {
    std::thread::sleep(std::time::Duration::from_secs((2 * attempt).min(30) as u64));
}

/// The script's `download_url`.
pub(super) fn download(url: &str, out: &Path, transfer: &Transfer) -> Result<()> {
    mkdir_parent(out)?;
    if needs_curl(transfer) {
        return curl(url, out, transfer, None, false);
    }
    let mut last = "None".to_owned();
    for attempt in 1..=transfer.retries {
        match urlopen(url, out, transfer, false) {
            Ok(()) => return Ok(()),
            Err(failure) => last = failure,
        }
        if attempt < transfer.retries {
            pause(attempt);
        }
    }
    Err(error(
        "RuntimeError",
        format!("failed to download {url}: {last}"),
    ))
}

fn temporary_path(prefix: &str, suffix: &str) -> PathBuf {
    let folder = std::env::var_os("TMPDIR").map_or_else(|| PathBuf::from("/tmp"), PathBuf::from);
    let tag: String = std::iter::repeat_with(|| {
        let alphabet = b"abcdefghijklmnopqrstuvwxyz0123456789_";
        let mut byte = [0u8; 1];
        File::open("/dev/urandom")
            .and_then(|mut f| f.read_exact(&mut byte))
            .ok();
        alphabet[byte[0] as usize % alphabet.len()] as char
    })
    .take(8)
    .collect();
    folder.join(format!("{prefix}{tag}{suffix}"))
}

/// The script's `fetch_json`.
fn fetch_json(url: &str, transfer: &Transfer) -> Result<Json> {
    let path = temporary_path("panelapp.", ".json");
    let result = (|| -> Result<Json> {
        if needs_curl(transfer) {
            curl(url, &path, transfer, None, false)?;
            let bytes = fs::read(&path).map_err(|e| os_error(e, &path, None))?;
            let body = crate::python_io::decode_utf8(&bytes)
                .map_err(|e| error("UnicodeDecodeError", e))?;
            return json_text::parse(body, false).map_err(|e| error("JSONDecodeError", e));
        }
        let mut last = "None".to_owned();
        for attempt in 1..=transfer.retries {
            match urlopen(url, &path, transfer, true) {
                Ok(()) => {
                    let bytes = fs::read(&path).map_err(|e| os_error(e, &path, None))?;
                    let body = crate::python_io::decode_utf8(&bytes)
                        .map_err(|e| error("UnicodeDecodeError", e))?;
                    match json_text::parse(body, false) {
                        Ok(value) => return Ok(value),
                        Err(failure) => last = failure,
                    }
                }
                Err(failure) => last = failure,
            }
            if attempt < transfer.retries {
                pause(attempt);
            }
        }
        Err(error(
            "RuntimeError",
            format!("failed to fetch JSON {url}: {last}"),
        ))
    })();
    let _ = fs::remove_file(&path);
    result
}

/// The script's `_download_text`.
fn download_text(url: &str, transfer: &Transfer) -> Result<String> {
    let path = temporary_path("clinvar.meta.", "");
    let result = (|| -> Result<String> {
        download(url, &path, transfer)?;
        let bytes = fs::read(&path).map_err(|e| os_error(e, &path, None))?;
        Ok(crate::python_io::decode_utf8(&bytes)
            .map_err(|e| error("UnicodeDecodeError", e))?
            .to_owned())
    })();
    let _ = fs::remove_file(&path);
    result
}

fn exception(failure: &(dyn std::error::Error + 'static)) -> String {
    message(failure)
}

/// The script's `fetch_static_source`.
pub(super) fn fetch_static_source(
    spec: &Spec,
    raw_dir: &Path,
    previous_sources: &Json,
    force: bool,
    transfer: &Transfer,
) -> Result<Json> {
    let raw_path = raw_dir.join(spec.raw_filename);
    let old = previous(spec, previous_sources, &raw_path)?;
    let matches = identity_matches(spec, &old, Some(&raw_path))?;
    let old_sha = get(&old, "sha256")?.cloned();
    if is_fresh(&raw_path, &old, spec.refresh_days, force || !matches)? {
        return metadata_for_file(spec, &raw_path, "skipped_fresh", old_sha.as_ref(), None);
    }
    let temporary = with_name(&raw_path, ".download.tmp");
    let attempt = (|| -> Result<Json> {
        download(&spec.url, &temporary, transfer)?;
        if spec.custom() || !matches {
            crate::mechanism_source::validate_payload(&temporary, spec.parser)?;
        }
        let new_sha = sha256_file(&temporary)?;
        let status = if matches
            && exists(&raw_path)
            && old_sha
                .as_ref()
                .is_some_and(|old| equal(old, &text(new_sha.clone())))
        {
            fs::remove_file(&temporary).map_err(|e| os_error(e, &temporary, None))?;
            "checked_same_checksum"
        } else {
            fs::rename(&temporary, &raw_path)
                .map_err(|e| os_error(e, &temporary, Some(&raw_path)))?;
            if old_sha.as_ref().is_some_and(truthy) {
                "updated"
            } else {
                "downloaded"
            }
        };
        metadata_for_file(spec, &raw_path, status, old_sha.as_ref(), None)
    })();
    match attempt {
        Ok(value) => Ok(value),
        Err(failure) => {
            if exists(&temporary) {
                fs::remove_file(&temporary).map_err(|e| os_error(e, &temporary, None))?;
            }
            let reason = exception(failure.as_ref());
            if spec.custom() || !matches {
                return Err(error(
                    "RuntimeError",
                    format!(
                        "Selected source {} ({}) failed; prior bytes preserved: {reason}",
                        spec.name, spec.url
                    ),
                ));
            }
            let status = if exists(&raw_path) {
                "using_existing_after_download_error"
            } else {
                "missing_after_download_error"
            };
            metadata_for_file(spec, &raw_path, status, old_sha.as_ref(), Some(reason))
        }
    }
}

/// `urllib.parse.urljoin(base, relative)` for the panel detail path.
fn join_url(base: &str, relative: &str) -> String {
    if relative.contains("://") {
        return relative.to_owned();
    }
    if let Some(rest) = relative.strip_prefix('/') {
        let scheme_end = base.find("://").map_or(0, |at| at + 3);
        let host_end = base[scheme_end..]
            .find('/')
            .map_or(base.len(), |at| scheme_end + at);
        return format!("{}/{rest}", &base[..host_end]);
    }
    format!("{base}{relative}")
}

fn schema_failure(message: &str) -> Box<dyn std::error::Error> {
    error("ValueError", message)
}

/// The script's `fetch_panelapp`.
pub(super) fn fetch_panelapp(
    spec: &Spec,
    raw_dir: &Path,
    previous_sources: &Json,
    force: bool,
    transfer: &Transfer,
) -> Result<Json> {
    let raw_path = raw_dir.join(spec.raw_filename);
    let old = previous(spec, previous_sources, &raw_path)?;
    let matches = identity_matches(spec, &old, Some(&raw_path))?;
    let old_sha = get(&old, "sha256")?.cloned();
    if is_fresh(&raw_path, &old, spec.refresh_days, force || !matches)? {
        return metadata_for_file(spec, &raw_path, "skipped_fresh", old_sha.as_ref(), None);
    }
    let checked = spec.custom() || !matches;
    let attempt = (|| -> Result<Json> {
        let mut panels: Vec<Json> = Vec::new();
        let mut next = Some(spec.url.clone());
        while let Some(url) = next.filter(|url| !url.is_empty()) {
            let payload = fetch_json(&url, transfer)?;
            if checked {
                let valid = match &payload {
                    Json::Object(_) => {
                        let results_ok = match payload.get("results") {
                            Some(Json::Array(rows)) => rows.iter().all(|row| {
                                matches!(row, Json::Object(_)) && row.get("id").is_some_and(truthy)
                            }),
                            _ => false,
                        };
                        results_ok
                            && !matches!(payload.get("next"), Some(v) if !matches!(v, Json::Null | Json::Str(_)))
                    }
                    _ => false,
                };
                if !valid {
                    return Err(schema_failure(
                        "Selected PanelApp URL must provide the original results/next list schema",
                    ));
                }
            }
            match get_or(&payload, "results", Json::Array(vec![]))? {
                Json::Array(rows) => panels.extend(rows),
                Json::Object(items) => panels.extend(items.into_iter().map(|(k, _)| Json::Str(k))),
                Json::Str(chars) => panels.extend(chars.chars().map(|c| Json::Str(c.into()))),
                other => {
                    return Err(error(
                        "TypeError",
                        format!("'{}' object is not iterable", type_name(&other)),
                    ));
                }
            }
            next = match get(&payload, "next")? {
                Some(Json::Str(url)) => Some(url.clone()),
                Some(Json::Null) | None => None,
                Some(other) if !truthy(other) => None,
                Some(other) => {
                    return Err(error(
                        "TypeError",
                        format!("'{}' object is not a valid URL", type_name(other)),
                    ));
                }
            };
            if let Some(cap) = transfer.max_panels {
                if panels.len() as i128 >= cap {
                    // panels[:cap], including a negative cap's slice.
                    let keep = if cap >= 0 {
                        cap
                    } else {
                        (panels.len() as i128 + cap).max(0)
                    };
                    panels.truncate(keep as usize);
                    break;
                }
            }
        }
        let base = format!(
            "{}/",
            spec.url
                .split('?')
                .next()
                .unwrap_or("")
                .trim_end_matches('/')
        );
        let mut detailed = Vec::new();
        for (index, panel) in panels.iter().enumerate() {
            let panel_id = get(panel, "id")?.cloned().unwrap_or(Json::Null);
            if !truthy(&panel_id) {
                continue;
            }
            let detail_url = join_url(&base, &format!("{}/?format=json", py_str(&panel_id)));
            let detail = match (|| -> Result<Json> {
                let detail = fetch_json(&detail_url, transfer)?;
                if checked {
                    let valid = match detail.get("genes") {
                        Some(Json::Array(genes)) if matches!(detail, Json::Object(_)) => {
                            genes.iter().all(|gene| {
                                matches!(gene, Json::Object(_))
                                    && !matches!(gene.get("gene_data"), Some(v) if !matches!(v, Json::Null | Json::Object(_)))
                                    && !matches!(gene.get("publications"), Some(v) if !matches!(v, Json::Null | Json::Array(_)))
                            })
                        }
                        _ => false,
                    };
                    if !valid {
                        return Err(schema_failure(
                            "Selected PanelApp detail must retain the original genes/gene_data/publications schema",
                        ));
                    }
                }
                Ok(detail)
            })() {
                Ok(detail) => detail,
                Err(failure) => {
                    if checked {
                        return Err(failure);
                    }
                    let Json::Object(items) = panel else {
                        return Err(error(
                            "TypeError",
                            format!("'{}' object is not a mapping", type_name(panel)),
                        ));
                    };
                    let mut merged = Json::Object(items.clone());
                    merged.set("genes", Json::Array(vec![]));
                    merged.set("detail_fetch_error", text(exception(failure.as_ref())));
                    merged.set("detail_url", text(detail_url.clone()));
                    merged
                }
            };
            detailed.push(detail);
            if (index + 1) % 50 == 0 {
                progress(&format!("panelapp_details_fetched={}", index + 1));
            }
        }
        let temporary = with_name(&raw_path, ".download.tmp");
        let payload = object([
            ("downloaded_at_utc", text(iso_now())),
            ("list_url", text(spec.url.clone())),
            ("panel_count", Json::Int(panels.len() as i128)),
            ("panels", Json::Array(detailed)),
        ]);
        let body = json_text::dumps_with_options(
            &payload,
            DumpsOptions {
                sort_keys: true,
                ..DumpsOptions::default()
            },
        )
        .map_err(|e| error("ValueError", e))?;
        fs::write(&temporary, format!("{body}\n")).map_err(|e| os_error(e, &temporary, None))?;
        let new_sha = sha256_file(&temporary)?;
        let status = if matches
            && exists(&raw_path)
            && old_sha
                .as_ref()
                .is_some_and(|old| equal(old, &text(new_sha.clone())))
        {
            fs::remove_file(&temporary).map_err(|e| os_error(e, &temporary, None))?;
            "checked_same_checksum"
        } else {
            fs::rename(&temporary, &raw_path)
                .map_err(|e| os_error(e, &temporary, Some(&raw_path)))?;
            if old_sha.as_ref().is_some_and(truthy) {
                "updated"
            } else {
                "downloaded"
            }
        };
        metadata_for_file(spec, &raw_path, status, old_sha.as_ref(), None)
    })();
    match attempt {
        Ok(value) => Ok(value),
        Err(failure) => {
            let reason = exception(failure.as_ref());
            if checked {
                return Err(error(
                    "RuntimeError",
                    format!(
                        "Selected source {} ({}) failed; prior bytes preserved: {reason}",
                        spec.name, spec.url
                    ),
                ));
            }
            let status = if exists(&raw_path) {
                "using_existing_after_download_error"
            } else {
                "missing_after_download_error"
            };
            metadata_for_file(spec, &raw_path, status, old_sha.as_ref(), Some(reason))
        }
    }
}

fn write_downloaded_text(path: &Path, body: &str) -> Result<Json> {
    mkdir_parent(path)?;
    let temporary = with_name(path, ".download.tmp");
    fs::write(&temporary, body).map_err(|e| os_error(e, &temporary, None))?;
    fs::rename(&temporary, path).map_err(|e| os_error(e, &temporary, Some(path)))?;
    Ok(object([
        ("path", text(path.display().to_string())),
        ("size_bytes", Json::Int(stat(path)?.len().into())),
        ("sha256", text(sha256_file(path)?)),
        ("md5", text(md5_file(path)?)),
    ]))
}

const README_URL: &str = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/_README";
const XSD_INDEX_URL: &str = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xsd_public/";

/// The script's `_latest_clinvar_xsd_url`.
fn latest_xsd(index: &str) -> Result<(String, String)> {
    let mut best: Option<(Vec<u128>, String)> = None;
    let mut rest = index;
    while let Some(at) = rest.find("ClinVar_VCV_") {
        let tail = &rest[at + "ClinVar_VCV_".len()..];
        rest = &rest[at + 1..];
        let mut end = 0;
        let bytes = tail.as_bytes();
        // (\d+(?:\.\d+)*) followed by ".xsd"
        loop {
            let start = end;
            while end < bytes.len() && bytes[end].is_ascii_digit() {
                end += 1;
            }
            if end == start {
                end = start.saturating_sub(1);
                break;
            }
            if bytes.get(end) == Some(&b'.') && bytes.get(end + 1).is_some_and(u8::is_ascii_digit) {
                end += 1;
                continue;
            }
            break;
        }
        let version = &tail[..end.min(tail.len())];
        if version.is_empty() || !tail[version.len()..].starts_with(".xsd") {
            continue;
        }
        let key: Vec<u128> = version
            .split('.')
            .map(|part| part.parse().unwrap_or(0))
            .collect();
        if best.as_ref().is_none_or(|(old, _)| key > *old) {
            best = Some((key, version.to_owned()));
        }
    }
    let (_, version) = best.ok_or_else(|| {
        error(
            "ValueError",
            "ClinVar XSD index contains no ClinVar_VCV_<version>.xsd",
        )
    })?;
    let url = format!("{XSD_INDEX_URL}ClinVar_VCV_{version}.xsd");
    Ok((version, url))
}

/// The script's `_parse_remote_md5`.
fn remote_md5(body: &str) -> Result<String> {
    let chars: Vec<char> = body.chars().collect();
    let word = |c: char| c.is_alphanumeric() || c == '_';
    for start in 0..chars.len() {
        if start + 32 > chars.len() {
            break;
        }
        if start > 0 && word(chars[start - 1]) {
            continue;
        }
        let candidate = &chars[start..start + 32];
        if candidate.iter().all(char::is_ascii_hexdigit)
            && !chars.get(start + 32).is_some_and(|c| word(*c))
        {
            return Ok(candidate.iter().collect::<String>().to_lowercase());
        }
    }
    let shown: String = body.chars().take(200).collect();
    Err(error(
        "ValueError",
        format!(
            "Could not parse ClinVar companion MD5: {}",
            crate::python_path::repr(&shown)
        ),
    ))
}

/// The script's `fetch_clinvar_vcv`.
pub(super) fn fetch_clinvar_vcv(
    raw_dir: &Path,
    previous_sources: &Json,
    force: bool,
    transfer: &Transfer,
) -> Result<Json> {
    let spec = clinvar_spec();
    let raw_path = raw_dir.join(spec.raw_filename);
    mkdir_parent(&raw_path)?;
    let old = previous(&spec, previous_sources, &raw_path)?;
    let last_check = match get(&old, "last_remote_md5_checked_at_utc")? {
        Some(Json::Str(stamp)) if !stamp.is_empty() => parse_iso(stamp),
        Some(other) if truthy(other) => {
            return Err(error(
                "TypeError",
                format!(
                    "fromisoformat: argument must be str, not {}",
                    type_name(other)
                ),
            ));
        }
        _ => None,
    };
    if !force && nonempty(&raw_path) {
        if let Some((micros, aware)) = last_check {
            if !aware {
                return Err(error(
                    "TypeError",
                    "can't subtract offset-naive and offset-aware datetimes",
                ));
            }
            if age_seconds(micros) < (spec.refresh_days * 86_400) as f64 {
                let mut metadata = old.clone();
                metadata.set("status", text("skipped_fresh"));
                metadata.set("last_build_seen_at_utc", text(iso_now()));
                return Ok(metadata);
            }
        }
    }
    let attempt = (|| -> Result<Json> {
        let md5_text = download_text(&format!("{}.md5", spec.url), transfer)?;
        let md5 = remote_md5(&md5_text)?;
        let readme = download_text(README_URL, transfer)?;
        let index = download_text(XSD_INDEX_URL, transfer)?;
        let (version, xsd_url) = latest_xsd(&index)?;
        let xsd = download_text(&xsd_url, transfer)?;
        let format_dir = raw_path.parent().unwrap().join("format");
        let readme_meta =
            write_downloaded_text(&format_dir.join("NCBI_CLINVAR_XML_README"), &readme)?;
        let xsd_meta =
            write_downloaded_text(&format_dir.join(format!("ClinVar_VCV_{version}.xsd")), &xsd)?;
        let local = if nonempty(&raw_path) {
            md5_file(&raw_path)?
        } else {
            String::new()
        };
        let mut downloaded = false;
        if local != md5 {
            let part = raw_path.with_file_name(format!(
                "{}.{md5}.part",
                raw_path.file_name().unwrap().to_string_lossy()
            ));
            curl(
                &spec.url,
                &part,
                transfer,
                Some(transfer.max_download_seconds),
                true,
            )?;
            let actual = md5_file(&part)?;
            if actual != md5 {
                let _ = fs::remove_file(&part);
                return Err(error(
                    "ValueError",
                    format!("ClinVar VCV MD5 mismatch: expected {md5}, got {actual}"),
                ));
            }
            fs::rename(&part, &raw_path).map_err(|e| os_error(e, &part, Some(&raw_path)))?;
            let stem = format!("{}.", raw_path.file_name().unwrap().to_string_lossy());
            if let Ok(entries) = fs::read_dir(raw_path.parent().unwrap()) {
                for entry in entries.flatten() {
                    let name = entry.file_name().to_string_lossy().into_owned();
                    if name.starts_with(&stem)
                        && name.ends_with(".part")
                        && name.len() > stem.len() + 4
                    {
                        let _ = fs::remove_file(entry.path());
                    }
                }
            }
            downloaded = true;
        }
        let sha256 = match get(&old, "sha256")? {
            Some(value) if !downloaded && truthy(value) => value.clone(),
            _ => text(sha256_file(&raw_path)?),
        };
        let has_old = truthy(&old);
        let changed = downloaded
            && get(&old, "md5")?
                .is_some_and(|value| truthy(value) && !equal(value, &text(md5.clone())));
        let metadata = object([
            ("name", text(spec.name)),
            ("url", text(spec.url.clone())),
            ("md5_url", text(format!("{}.md5", spec.url))),
            ("docs_url", text(spec.docs_url)),
            ("raw_path", text(raw_path.display().to_string())),
            ("source_role", text(spec.source_role)),
            ("official_update_note", text(spec.official_update_note)),
            ("refresh_days", Json::Int(spec.refresh_days.into())),
            ("parser", text(spec.parser)),
            ("last_checked_at_utc", text(iso_now())),
            ("last_remote_md5_checked_at_utc", text(iso_now())),
            (
                "last_downloaded_at_utc",
                if downloaded {
                    text(iso_now())
                } else {
                    get_or(&old, "last_downloaded_at_utc", text(""))?
                },
            ),
            (
                "status",
                text(if downloaded && has_old {
                    "updated"
                } else if downloaded {
                    "downloaded"
                } else {
                    "checked_same_remote_md5"
                }),
            ),
            ("error", text("")),
            ("size_bytes", Json::Int(stat(&raw_path)?.len().into())),
            ("sha256", sha256),
            ("md5", text(md5.clone())),
            ("remote_md5", text(md5.clone())),
            (
                "remote_md5_text",
                text(crate::python_text::python_strip(&md5_text)),
            ),
            ("file_mtime_utc", text(file_mtime_iso(&raw_path)?)),
            ("content_changed", Json::Bool(changed)),
            (
                "format_metadata",
                object([
                    ("readme_url", text(README_URL)),
                    ("readme", readme_meta),
                    ("xsd_index_url", text(XSD_INDEX_URL)),
                    ("xsd_version", text(version)),
                    ("xsd_url", text(xsd_url)),
                    ("xsd", xsd_meta),
                ]),
            ),
        ]);
        let receipt = with_name(&raw_path, ".md5");
        fs::write(&receipt, &md5_text).map_err(|e| os_error(e, &receipt, None))?;
        Ok(metadata)
    })();
    match attempt {
        Ok(value) => Ok(value),
        Err(failure) => {
            let reason = exception(failure.as_ref());
            if nonempty(&raw_path) {
                let mut metadata = old.clone();
                for (key, value) in [
                    ("name", text(spec.name)),
                    ("url", text(spec.url.clone())),
                    ("docs_url", text(spec.docs_url)),
                    ("raw_path", text(raw_path.display().to_string())),
                    ("source_role", text(spec.source_role)),
                    ("official_update_note", text(spec.official_update_note)),
                    ("refresh_days", Json::Int(spec.refresh_days.into())),
                    ("parser", text(spec.parser)),
                    ("last_checked_at_utc", text(iso_now())),
                    ("status", text("using_existing_after_download_error")),
                    ("error", text(reason)),
                    ("size_bytes", Json::Int(stat(&raw_path)?.len().into())),
                    ("file_mtime_utc", text(file_mtime_iso(&raw_path)?)),
                ] {
                    metadata.set(key, value);
                }
                return Ok(metadata);
            }
            Ok(object([
                ("name", text(spec.name)),
                ("url", text(spec.url.clone())),
                ("docs_url", text(spec.docs_url)),
                ("raw_path", text(raw_path.display().to_string())),
                ("source_role", text(spec.source_role)),
                ("official_update_note", text(spec.official_update_note)),
                ("refresh_days", Json::Int(spec.refresh_days.into())),
                ("parser", text(spec.parser)),
                ("last_checked_at_utc", text(iso_now())),
                ("status", text("missing_after_download_error")),
                ("error", text(reason)),
            ]))
        }
    }
}

/// The script's `carry_forward_source_metadata`.
pub(super) fn carry_forward(spec: &Spec, raw_dir: &Path, previous_sources: &Json) -> Result<Json> {
    let raw_path = raw_dir.join(spec.raw_filename);
    let mut old = get_or(previous_sources, spec.name, Json::Object(vec![]))?;
    if !matches!(old, Json::Object(_)) {
        return Err(error(
            "TypeError",
            format!("'{}' object is not iterable", type_name(&old)),
        ));
    }
    if !identity_matches(spec, &old, Some(&raw_path))? {
        return Err(error(
            "ValueError",
            format!(
                "ClinVar-only refresh cannot change selected source URL: {}; run a complete refresh",
                spec.name
            ),
        ));
    }
    for (key, value) in [
        ("name", text(spec.name)),
        ("url", text(spec.url.clone())),
        ("docs_url", text(spec.docs_url)),
        ("raw_path", text(raw_path.display().to_string())),
        ("source_role", text(spec.source_role)),
        ("official_update_note", text(spec.official_update_note)),
        ("refresh_days", Json::Int(spec.refresh_days.into())),
        ("parser", text(spec.parser)),
        ("status", text("not_checked_clinvar_only_refresh")),
        ("last_build_seen_at_utc", text(iso_now())),
    ] {
        old.set(key, value);
    }
    if !nonempty(&raw_path) {
        return Err(error(
            "FileNotFoundError",
            format!(
                "ClinVar-only refresh requires existing source file: {}",
                raw_path.display()
            ),
        ));
    }
    Ok(old)
}

fn local_name(name: &str) -> &str {
    name.rsplit('}').next().unwrap_or(name)
}

/// `str.isdigit()`: decimal digits and the other Unicode digit characters.
pub(crate) fn is_digit(c: char) -> bool {
    crate::python_text::decimal_digit_value(c).is_some()
        || matches!(c,
            '\u{b2}' | '\u{b3}' | '\u{b9}' | '\u{1369}'..='\u{1371}' | '\u{19da}'
            | '\u{2070}' | '\u{2074}'..='\u{2079}' | '\u{2080}'..='\u{2089}'
            | '\u{2460}'..='\u{2468}' | '\u{2474}'..='\u{247c}' | '\u{2488}'..='\u{2490}'
            | '\u{24ea}' | '\u{24f5}'..='\u{24fd}' | '\u{24ff}' | '\u{2776}'..='\u{277e}'
            | '\u{2780}'..='\u{2788}' | '\u{278a}'..='\u{2792}' | '\u{10a40}'..='\u{10a43}'
            | '\u{1f100}'..='\u{1f10a}')
}

/// The script's `selected_clinvar_vcv_metadata`.
pub(super) fn selected_clinvar_vcv_metadata(supplied: &Path) -> Result<Json> {
    let path = absolute(&expanduser(supplied)?)?;
    if !crate::python_path::is_file(&path)? || stat(&path)?.len() == 0 {
        return Err(error(
            "ValueError",
            format!(
                "Missing or empty explicit ClinVar VCV XML: {}",
                path.display()
            ),
        ));
    }
    let identity = |meta: &fs::Metadata| {
        (
            meta.dev(),
            meta.ino(),
            meta.len(),
            meta.mtime(),
            meta.mtime_nsec(),
            meta.ctime(),
            meta.ctime_nsec(),
        )
    };
    let before_meta = stat(&path)?;
    let before = identity(&before_meta);
    let sha256 = sha256_file(&path)?;
    let md5 = md5_file(&path)?;
    let state = std::rc::Rc::new(std::cell::RefCell::new((
        None::<String>,
        None::<String>,
        0u64,
    )));
    let display = path.display().to_string();
    let file = File::open(&path).map_err(|e| os_error(e, &path, None))?;
    let source: Box<dyn Read> = if display.ends_with(".gz") {
        Box::new(crate::python_io::PythonGzipReader::new(file))
    } else {
        Box::new(file)
    };
    let start_state = state.clone();
    let start_display = display.clone();
    let end_state = state.clone();
    let end_display = display.clone();
    let ids_end = std::rc::Rc::new(std::cell::RefCell::new(Vec::<Option<String>>::new()));
    let ids_start = ids_end.clone();
    crate::interpro_mapping_builder::visit_events(
        source,
        Box::new(move |name: &str, attributes: &[(String, String)]| {
            let local = local_name(name).to_owned();
            let mut state = start_state.borrow_mut();
            if state.0.is_none() {
                if local != "ClinVarVariationRelease" {
                    return Err(error(
                        "ValueError",
                        format!("Not a ClinVar VCV XML release: {start_display}; root={local}"),
                    ));
                }
                state.0 = Some(local.clone());
                state.1 = attributes
                    .iter()
                    .find(|(key, _)| key == "ReleaseDate")
                    .map(|(_, value)| value.clone());
            }
            ids_start.borrow_mut().push(if local == "VariationArchive" {
                Some(
                    attributes
                        .iter()
                        .find(|(key, _)| key == "VariationID")
                        .map_or(String::new(), |(_, value)| value.clone()),
                )
            } else {
                None
            });
            Ok(())
        }),
        Box::new(move |_name: &str| {
            let id = ids_end.borrow_mut().pop().flatten();
            if let Some(id) = id {
                if id.is_empty() || !id.chars().all(is_digit) {
                    return Err(error(
                        "ValueError",
                        format!("ClinVar VCV record lacks a numeric VariationID: {end_display}"),
                    ));
                }
                end_state.borrow_mut().2 += 1;
            }
            Ok(())
        }),
    )?;
    let (root, release, records) = state.borrow().clone();
    if records == 0 {
        return Err(error(
            "ValueError",
            format!("ClinVar VCV XML contains no VariationArchive records: {display}"),
        ));
    }
    if identity(&stat(&path)?) != before {
        return Err(error(
            "ValueError",
            format!(
                "Explicit ClinVar VCV XML changed during validation: {display}; wait for its writer"
            ),
        ));
    }
    let spec = clinvar_spec();
    Ok(object([
        ("name", text(spec.name)),
        ("url", Json::Null),
        ("md5_url", Json::Null),
        ("remote_md5", Json::Null),
        ("raw_path", text(display)),
        ("source_role", text(spec.source_role)),
        ("parser", text(spec.parser)),
        ("source_mode", text("explicit_local_input")),
        ("publisher_identity_verified", Json::Bool(false)),
        ("status", text("validated_explicit_local_input")),
        ("last_checked_at_utc", text(iso_now())),
        ("size_bytes", Json::Int(before_meta.len().into())),
        ("sha256", text(sha256)),
        ("md5", text(md5)),
        ("file_mtime_utc", text(file_mtime_iso(&path)?)),
        ("xml_root", root.map_or(Json::Null, text)),
        ("declared_release_date", release.map_or(Json::Null, text)),
        ("variation_archive_records", Json::Int(records.into())),
        ("error", text("")),
    ]))
}

/// The script's `require_selected_clinvar_injection`.
pub(super) fn require_selected_clinvar_injection(meta: &Json, gofcards: &Path) -> Result<()> {
    let cache = build::load_gofcards_cache(gofcards)?;
    let metadata = get_or(&cache, "metadata", Json::Object(vec![]))?;
    let injected = get_or(&metadata, "clinvar", Json::Object(vec![]))?;
    let source = get(&injected, "source")?.cloned().unwrap_or(Json::Null);
    let matches = match &source {
        Json::Str(value) => {
            let resolved = absolute(&expanduser(&path_value(value))?)?;
            meta.get("raw_path")
                .is_some_and(|raw| equal(raw, &text(resolved.display().to_string())))
                && get(&injected, "source_sha256")?
                    .is_some_and(|sha| meta.get("sha256").is_some_and(|own| equal(sha, own)))
        }
        _ => false,
    };
    if !matches {
        return Err(error(
            "ValueError",
            "Explicit ClinVar VCV XML path/SHA256 does not match the injected GoFCards cache; \
run injection with this selected input before the full non-LOF build",
        ));
    }
    Ok(())
}

fn cell(value: Option<&Json>) -> ObjectValue {
    match value {
        None | Some(Json::Null) => ObjectValue::None,
        Some(Json::Str(text)) => ObjectValue::Text(text.clone()),
        Some(Json::Int(number)) => i64::try_from(*number)
            .map_or_else(|_| ObjectValue::Text(number.to_string()), ObjectValue::Int),
        Some(Json::Float(number)) => ObjectValue::Float(*number),
        Some(Json::Bool(flag)) => ObjectValue::Bool(*flag),
        Some(other) => ObjectValue::Text(repr(other)),
    }
}

/// The script's `write_source_manifest_tsv`.
pub(super) fn write_source_manifest_tsv(path: &Path, manifest: &Json) -> Result<usize> {
    let names = [
        "source",
        "status",
        "url",
        "docs_url",
        "raw_path",
        "refresh_days",
        "last_checked_at_utc",
        "last_downloaded_at_utc",
        "file_mtime_utc",
        "size_bytes",
        "sha256",
        "md5",
        "remote_md5",
        "xsd_version",
        "xsd_path",
        "content_changed",
        "error",
        "official_update_note",
    ];
    let mut sources = match get_or(manifest, "sources", Json::Object(vec![]))? {
        Json::Object(items) => items,
        _ => vec![],
    };
    sources.sort_by(|a, b| a.0.cmp(&b.0));
    let mut columns: Vec<Vec<ObjectValue>> = vec![Vec::new(); names.len()];
    for (name, meta) in &sources {
        let field = |key: &str| -> Result<Option<Json>> { Ok(get(meta, key)?.cloned()) };
        let stringified = |key: &str| -> Result<ObjectValue> {
            Ok(ObjectValue::Text(py_str(&field(key)?.unwrap_or(text("")))))
        };
        let plain =
            |key: &str| -> Result<ObjectValue> { Ok(cell(Some(&field(key)?.unwrap_or(text(""))))) };
        let format = get_or(meta, "format_metadata", Json::Object(vec![]))?;
        let xsd = get_or(&format, "xsd", Json::Object(vec![]))?;
        let row = [
            ObjectValue::Text(name.clone()),
            plain("status")?,
            plain("url")?,
            plain("docs_url")?,
            plain("raw_path")?,
            stringified("refresh_days")?,
            plain("last_checked_at_utc")?,
            plain("last_downloaded_at_utc")?,
            plain("file_mtime_utc")?,
            stringified("size_bytes")?,
            plain("sha256")?,
            plain("md5")?,
            plain("remote_md5")?,
            cell(Some(&get_or(&format, "xsd_version", text(""))?)),
            cell(Some(&get_or(&xsd, "path", text(""))?)),
            stringified("content_changed")?,
            plain("error")?,
            plain("official_update_note")?,
        ];
        for (column, value) in columns.iter_mut().zip(row) {
            column.push(value);
        }
    }
    let table = Table {
        names: names.iter().map(|n| (*n).to_owned()).collect(),
        rows: sources.len(),
        columns: columns.into_iter().map(Column::Object).collect(),
    };
    mkdir_parent(path)?;
    let temporary = with_name(path, ".tmp");
    crate::pandas_table::write_table(&temporary, &table)
        .map_err(|e| error("OSError", message(e.as_ref())))?;
    fs::rename(&temporary, path).map_err(|e| os_error(e, &temporary, Some(path)))?;
    Ok(sources.len())
}
