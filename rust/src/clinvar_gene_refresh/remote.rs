//! The bounded ClinVar fileDate checks and the pinned NCBI acquisitions.
//!
//! urlopen's rules are kept with `curl --http1.1`: the same User-Agent and
//! identity encoding, redirects followed, a final status outside 2xx is
//! `HTTPError` "HTTP Error CODE: REASON", and a transfer failure is
//! `URLError`. The compressed header read stops after 1 MiB even when the
//! server ignores Range.
use super::*;
use crate::python_io::{GzipFailure, PythonGzipReader};
use std::io::{BufRead, BufReader};
use std::process::{Child, ChildStdout, Command, Stdio};

const TIMEOUT: &str = "30";
const HEADER_BYTES: usize = 2 * 1024 * 1024;
const COMPRESSED_HEADER_BYTES: usize = 1024 * 1024;
const BOUNDED: &str = "ClinVar compressed header exceeds the bounded metadata limit";

fn url_error(message: impl std::fmt::Display) -> Box<dyn std::error::Error> {
    error("URLError", format!("<urlopen error {message}>"))
}

/// An open HTTP response; its curl process ends when this is dropped.
struct Response {
    child: Child,
    body: BufReader<ChildStdout>,
    headers: Vec<(String, String)>,
}

impl Response {
    fn header(&self, name: &str) -> Json {
        self.headers
            .iter()
            .find(|(key, _)| key.eq_ignore_ascii_case(name))
            .map_or(Json::Null, |(_, value)| string(value.clone()))
    }
    /// Waits for curl after the body was read; a failed transfer is URLError.
    fn finish(mut self) -> Result<()> {
        let output = self.child.wait().map_err(|e| os_error(e, Path::new("curl"), None))?;
        let mut reason = String::new();
        if let Some(mut stderr) = self.child.stderr.take() {
            stderr.read_to_string(&mut reason).ok();
        }
        if output.success() {
            Ok(())
        } else if output.code() == Some(18) {
            Err(error("IncompleteRead", reason.trim().to_owned()))
        } else {
            Err(url_error(reason.trim()))
        }
    }
}

impl Drop for Response {
    fn drop(&mut self) {
        if matches!(self.child.try_wait(), Ok(None)) {
            let _ = self.child.kill();
        }
        let _ = self.child.wait();
    }
}

/// The status line's code and reason, skipping informational replies and
/// the redirects curl followed.
fn read_head(reader: &mut impl BufRead) -> io::Result<Option<(u16, String, Vec<(String, String)>)>> {
    loop {
        let mut line = Vec::new();
        if reader.read_until(b'\n', &mut line)? == 0 {
            return Ok(None);
        }
        let text = String::from_utf8_lossy(&line).trim_end_matches(['\r', '\n']).to_owned();
        let mut words = text.splitn(3, ' ');
        let version = words.next().unwrap_or_default();
        let status = words.next().and_then(|code| code.parse::<u16>().ok());
        let (true, Some(status)) = (version.starts_with("HTTP/"), status) else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("no HTTP status line: {text}"),
            ));
        };
        let reason = words.next().unwrap_or_default().to_owned();
        let mut headers = Vec::new();
        loop {
            line.clear();
            if reader.read_until(b'\n', &mut line)? == 0 {
                return Ok(None);
            }
            let field = String::from_utf8_lossy(&line).trim_end_matches(['\r', '\n']).to_owned();
            if field.is_empty() {
                break;
            }
            if let Some((name, value)) = field.split_once(':') {
                headers.push((name.to_owned(), value.trim_matches([' ', '\t']).to_owned()));
            }
        }
        let location = headers.iter().any(|(name, _)| name.eq_ignore_ascii_case("location"));
        if (100..200).contains(&status)
            || ((300..400).contains(&status) && location && reader.fill_buf()?.starts_with(b"HTTP/"))
        {
            continue;
        }
        return Ok(Some((status, reason, headers)));
    }
}

/// `open_url(url, header_only=...)`
fn open_url(url: &str, header_only: bool) -> Result<Response> {
    let mut command = Command::new("curl");
    command.args([
        "--silent",
        "--show-error",
        "--location",
        "--http1.1",
        "--include",
        "--suppress-connect-headers",
        "--connect-timeout",
        TIMEOUT,
        "--speed-limit",
        "1",
        "--speed-time",
        TIMEOUT,
        "--header",
        "User-Agent: PriVA-ClinVar-refresh/1",
        "--header",
        "Accept-Encoding: identity",
    ]);
    if header_only {
        command.args([
            "--header",
            &format!("Range: bytes=0-{}", COMPRESSED_HEADER_BYTES - 1),
        ]);
    }
    command
        .arg(url)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    let mut child = command
        .spawn()
        .map_err(|e| os_error(e, Path::new("curl"), None))?;
    let mut body = BufReader::new(child.stdout.take().expect("curl output is piped"));
    let head = read_head(&mut body);
    let mut response = Response {
        child,
        body,
        headers: Vec::new(),
    };
    let Some((status, reason, headers)) = head.map_err(url_error)? else {
        let _ = response.child.wait();
        let mut reason = String::new();
        if let Some(mut stderr) = response.child.stderr.take() {
            stderr.read_to_string(&mut reason).ok();
        }
        return Err(url_error(if reason.trim().is_empty() {
            "no HTTP response".to_owned()
        } else {
            reason.trim().to_owned()
        }));
    };
    if !(200..300).contains(&status) {
        return Err(error("HTTPError", format!("HTTP Error {status}: {reason}")));
    }
    response.headers = headers;
    Ok(response)
}

/// A response body that refuses to read past the bounded metadata limit.
struct Bounded<R> {
    source: R,
    total: usize,
}

impl<R: Read> Read for Bounded<R> {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        let remaining = COMPRESSED_HEADER_BYTES - self.total.min(COMPRESSED_HEADER_BYTES);
        let count = buffer.len().min(remaining + 1);
        let read = self.source.read(&mut buffer[..count])?;
        self.total += read;
        if self.total > COMPRESSED_HEADER_BYTES {
            return Err(io::Error::other(BOUNDED));
        }
        Ok(read)
    }
}

/// A failure while reading header lines, as Python raises it.
fn read_failure(cause: io::Error, path: &str) -> Box<dyn std::error::Error> {
    match GzipFailure::of(&cause) {
        Some(GzipFailure::EndOfFile) => error("EOFError", cause.to_string()),
        Some(GzipFailure::Bad(_)) => error("BadGzipFile", cause.to_string()),
        Some(GzipFailure::Zlib(_)) => error("error", cause.to_string()),
        None if cause.to_string() == BOUNDED => value_error(BOUNDED),
        None => os_error(cause, Path::new(path), None),
    }
}

/// `handle.readline(limit)`
fn readline(reader: &mut impl BufRead, limit: usize) -> io::Result<Vec<u8>> {
    let mut line = Vec::new();
    reader.take(limit as u64).read_until(b'\n', &mut line)?;
    Ok(line)
}

/// `datetime.strptime(stamp, "%Y-%m-%d" | "%Y%m%d").date()` for the spellings
/// the fileDate pattern admits: strptime's regex alternatives in their order,
/// with no end anchor, then the date's own range checks.
fn stamp_date(stamp: &str) -> Result<Date> {
    let dashed = stamp.contains('-');
    let format = if dashed { "%Y-%m-%d" } else { "%Y%m%d" };
    let bytes = stamp.as_bytes();
    let digit = |at: usize, low: u8, high: u8| bytes.get(at).is_some_and(|b| (low..=high).contains(b));
    // %Y takes four digits; the dashed form then has its literal '-'.
    let first = if dashed { 5 } else { 4 };
    // %m: 1[0-2]|0[1-9]|[1-9]
    let months = [
        (digit(first, b'1', b'1') && digit(first + 1, b'0', b'2')).then_some(2),
        (digit(first, b'0', b'0') && digit(first + 1, b'1', b'9')).then_some(2),
        digit(first, b'1', b'9').then_some(1),
    ];
    for width in months.into_iter().flatten() {
        let mut at = first + width;
        if dashed {
            if bytes.get(at) != Some(&b'-') {
                continue;
            }
            at += 1;
        }
        // %d: 3[01]|[12]\d|0[1-9]|[1-9]| [1-9]
        let days = [
            (digit(at, b'3', b'3') && digit(at + 1, b'0', b'1')).then_some(2),
            (digit(at, b'1', b'2') && digit(at + 1, b'0', b'9')).then_some(2),
            (digit(at, b'0', b'0') && digit(at + 1, b'1', b'9')).then_some(2),
            digit(at, b'1', b'9').then_some(1),
        ];
        let Some(day_width) = days.into_iter().flatten().next() else {
            continue;
        };
        let end = at + day_width;
        if end != stamp.len() {
            return Err(value_error(format!(
                "unconverted data remains: {}",
                &stamp[end..]
            )));
        }
        let number = |text: &str| text.parse::<u32>().unwrap_or_default();
        return date(
            number(&stamp[..4]) as i128,
            number(&stamp[first..first + width]),
            number(&stamp[at..end]),
        );
    }
    Err(value_error(format!(
        "time data {} does not match format {}",
        repr(stamp),
        repr(format)
    )))
}

/// `##fileDate=(\d{8}|\d{4}-\d{2}-\d{2})\r?\n?` as the whole line.
fn file_date_stamp(line: &[u8]) -> Option<&str> {
    let body = line.strip_prefix(b"##fileDate=")?;
    let body = body.strip_suffix(b"\n").unwrap_or(body);
    let body = body.strip_suffix(b"\r").unwrap_or(body);
    let digits = |part: &[u8]| part.iter().all(u8::is_ascii_digit);
    let valid = match body.len() {
        8 => digits(body),
        10 => body[4] == b'-' && body[7] == b'-' && digits(&body[..4]) && digits(&body[5..7]) && digits(&body[8..]),
        _ => false,
    };
    valid.then(|| std::str::from_utf8(body).unwrap_or_default())
}

/// `header_date(handle)`: the complete bounded header, never variant records.
fn header_date(reader: &mut impl BufRead, path: &str) -> Result<Date> {
    let mut result: Option<Date> = None;
    let mut consumed = 0;
    while consumed < HEADER_BYTES {
        let line = readline(reader, 65537.min(HEADER_BYTES - consumed + 1))
            .map_err(|cause| read_failure(cause, path))?;
        consumed += line.len();
        require(
            line.len() <= 65536 && consumed <= HEADER_BYTES,
            "ClinVar header exceeds the bounded metadata limit",
        )?;
        if line.is_empty() {
            break;
        }
        if line.starts_with(b"##fileDate=") {
            let Some(stamp) = file_date_stamp(&line) else {
                return Err(value_error(
                    "ClinVar fileDate must be a YYYYMMDD or YYYY-MM-DD data date",
                ));
            };
            let found = stamp_date(stamp)?;
            require(
                result.is_none_or(|known| known == found),
                "Conflicting ClinVar fileDate header values",
            )?;
            result = Some(found);
        } else if line.starts_with(b"#CHROM\t") {
            return result.ok_or_else(|| value_error("ClinVar VCF header has no fileDate data date"));
        } else {
            require(
                line.starts_with(b"##"),
                "Expected a ClinVar VCF header before variant records",
            )?;
        }
    }
    Err(value_error(
        "ClinVar metadata ended before a complete VCF header/fileDate",
    ))
}

/// `local_file_date(path)`: plain or gzip-compressed, by its magic bytes.
pub(super) fn local_file_date(path: &str) -> Result<Date> {
    nonempty(path, "selected ClinVar VCF")?;
    let open = || File::open(path).map_err(|e| os_error(e, Path::new(path), None));
    let mut magic = [0; 2];
    let read = open()?
        .read(&mut magic)
        .map_err(|e| os_error(e, Path::new(path), None))?;
    if read == 2 && magic == [0x1f, 0x8b] {
        header_date(&mut BufReader::new(PythonGzipReader::new(open()?)), path)
    } else {
        header_date(&mut BufReader::new(open()?), path)
    }
}

/// `publisher_md5(url)`: the publisher's `.md5` manifest.
pub(super) fn publisher_md5(url: &str) -> Result<String> {
    let manifest = format!("{url}.md5");
    let mut response = open_url(&manifest, false)?;
    let mut text = Vec::new();
    (&mut response.body)
        .take(4097)
        .read_to_end(&mut text)
        .map_err(url_error)?;
    let head = text.get(..32).unwrap_or_default();
    let valid = text.len() <= 4096
        && text.len() >= 32
        && head.iter().all(u8::is_ascii_hexdigit)
        && text
            .get(32)
            .is_none_or(|next| matches!(next, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c));
    require(valid, format!("Invalid publisher MD5 manifest: {manifest}"))?;
    Ok(String::from_utf8_lossy(head).to_lowercase())
}

/// `latest_url_template.format(assembly=assembly)`
fn format_template(template: &str, assembly: &str) -> Result<String> {
    let mut out = String::new();
    let mut rest = template;
    while let Some(at) = rest.find(['{', '}']) {
        out.push_str(&rest[..at]);
        let tail = &rest[at..];
        if tail.starts_with("{{") || tail.starts_with("}}") {
            out.push_str(&tail[..1]);
            rest = &tail[2..];
        } else if let Some(field) = tail.strip_prefix('{') {
            let Some(end) = field.find('}') else {
                return Err(value_error("Single '{' encountered in format string"));
            };
            match &field[..end] {
                "assembly" => out.push_str(assembly),
                "" => {
                    return Err(error(
                        "IndexError",
                        "Replacement index 0 out of range for positional args tuple",
                    ));
                }
                name => return Err(error("KeyError", repr(name))),
            }
            rest = &field[end + 1..];
        } else {
            return Err(value_error("Single '}' encountered in format string"));
        }
    }
    out.push_str(rest);
    Ok(out)
}

/// `latest_metadata(assembly, url_template=...)`: the official fileDate from
/// a bounded header read, with the publisher MD5 stable around it.
pub(super) fn latest_metadata(assembly: &str, template: &str) -> Result<Json> {
    let url = format_template(template, assembly)?;
    let observed = (|| -> Result<(Date, String)> {
        let before = publisher_md5(&url)?;
        let mut response = open_url(&url, true)?;
        let found = header_date(
            &mut BufReader::new(PythonGzipReader::new(Bounded {
                source: &mut response.body,
                total: 0,
            })),
            &url,
        )?;
        drop(response);
        require(
            publisher_md5(&url)? == before,
            "publisher MD5 changed during bounded header detection",
        )?;
        Ok((found, before))
    })();
    let (found, md5) = match observed {
        Ok(value) => value,
        Err(failure) => {
            let (category, message) = exception(failure.as_ref());
            if !matches!(
                category,
                "OSError"
                    | "FileNotFoundError"
                    | "PermissionError"
                    | "URLError"
                    | "HTTPError"
                    | "BadGzipFile"
                    | "EOFError"
                    | "ValueError"
                    | "UnicodeDecodeError"
                    | "error"
            ) {
                return Err(failure);
            }
            return Err(value_error(format!(
                "Cannot determine latest official ClinVar fileDate from {url}: {message}"
            )));
        }
    };
    Ok(object([
        ("url", string(url)),
        ("file_date", string(date_text(found))),
        ("observed_utc", string(utc_now())),
        ("publisher_md5", string(md5)),
        (
            "date_source",
            string("VCF ##fileDate; not download or local cache mtime"),
        ),
    ]))
}

/// `acquire_file(url, destination, expected_md5=...)`: completed external
/// bytes are reused; a download is published only after its checks.
fn acquire_file(url: &str, destination: &str, expected_md5: Option<&str>) -> Result<Json> {
    let sidecar = format!("{destination}.source.json");
    let present = |path: &str| crate::installer_prerequisites::paths::exists(path);
    let recorded = if present(&sidecar)? {
        read_json(&sidecar)?
    } else {
        Json::Object(Vec::new())
    };
    let same_url = lookup(&recorded, "url")?.is_some_and(|value| value.as_str() == Some(url));
    let mut provenance = if same_url {
        get_or_null(&recorded, "http_response")?
    } else {
        Json::Null
    };
    let mut reusable = false;
    if present(destination)? {
        if let Some(expected) = expected_md5 {
            reusable = md5(destination)? == expected;
        } else if present(&sidecar)? {
            let recorded_sha = lookup(&recorded, "sha256")?.cloned();
            reusable = same_url
                && json_text::py_equal(
                    &recorded_sha.unwrap_or(Json::Null),
                    &string(sha256(destination)?),
                );
        }
    }
    if !reusable {
        let parent = parent_of(destination);
        mkdirs(&parent)?;
        let (mut output, temporary) = mkstemp(&parent, &format!("{}.partial.", name_of(destination)))?;
        let result = (|| -> Result<()> {
            let mut response = open_url(url, false)?;
            provenance = object([
                ("last_modified", response.header("Last-Modified")),
                ("etag", response.header("ETag")),
                ("content_length", response.header("Content-Length")),
                ("retrieved_utc", string(utc_now())),
            ]);
            io::copy(&mut response.body, &mut output)
                .map_err(|e| os_error(e, Path::new(&temporary), None))?;
            response.finish()?;
            drop(output);
            nonempty(&temporary, "downloaded ClinVar refresh input")?;
            if let Some(expected) = expected_md5 {
                require(
                    md5(&temporary)? == expected,
                    format!("Downloaded bytes differ from publisher MD5: {url}"),
                )?;
            }
            if destination.ends_with(".gz") {
                let file = File::open(&temporary).map_err(|e| os_error(e, Path::new(&temporary), None))?;
                io::copy(&mut PythonGzipReader::new(file), &mut io::sink())
                    .map_err(|cause| read_failure(cause, &temporary))?;
            }
            rename(&temporary, destination)
        })();
        if crate::installer_prerequisites::paths::lexists(&temporary) {
            unlink(&temporary)?;
        }
        result?;
    }
    let mut result = identity_with(
        destination,
        vec![("url", string(url)), ("md5", string(md5(destination)?))],
    )?;
    if !matches!(provenance, Json::Null) {
        result.set("http_response", provenance);
    }
    if let Some(expected) = expected_md5 {
        result.set("publisher_md5", string(expected));
    }
    write_json(&sidecar, &result)?;
    Ok(result)
}

/// `acquire_inputs(plan)`: the pinned ClinVar VCF (or the explicit complete
/// annotation), variant_summary and MANE, with the release stable around them.
pub(super) fn acquire_inputs(plan: &Json) -> Result<Json> {
    let build = path_of(item(plan, "build_dir")?)?;
    let folder = join(&build, "inputs");
    mkdirs(&folder)?;
    let latest = item(plan, "latest")?;
    let policy = item(plan, "resource_policy")?;
    let file_date = text_of(item(latest, "file_date")?);
    let mut sources = Json::Object(Vec::new());
    if truth(item(plan, "explicit_vcf")?) {
        let path = path_of(item(plan, "input_vcf")?)?;
        sources.set(
            "vcf",
            identity_with(&path, vec![("kind", string("explicit_complete_annotated_vcf"))])?,
        );
    } else {
        let stamp = file_date.replace('-', "");
        // Weekly releases need not have a dated archive URL. Pin the moving
        // endpoint with publisher MD5 and verify it again after acquisition.
        let assembly = text_of(item(plan, "assembly")?);
        let url = format_template(&text_of(item(policy, "latest_url_template")?), &assembly)?;
        let expected = publisher_md5(&url)?;
        require(
            lookup(latest, "publisher_md5")?.and_then(Json::as_str) == Some(expected.as_str()),
            "Official ClinVar changed after latest-data detection; no replacement selected",
        )?;
        let vcf = join(&folder, &format!("clinvar_{stamp}.{assembly}.md5_{expected}.vcf.gz"));
        sources.set("vcf", acquire_file(&url, &vcf, Some(&expected))?);
        sources.set(
            "index",
            acquire_file(&format!("{url}.tbi"), &format!("{vcf}.tbi"), None)?,
        );
    }
    let vcf = path_of(item(item(&sources, "vcf")?, "path")?)?;
    require(
        date_text(local_file_date(&vcf)?) == file_date,
        "Pinned ClinVar input fileDate differs from the detected latest release",
    )?;
    let summary_url = text_of(item(policy, "summary_url")?);
    let mane_url = text_of(item(policy, "mane_url")?);
    let summary_md5 = publisher_md5(&summary_url)?;
    sources.set(
        "summary",
        acquire_file(
            &summary_url,
            &join(&folder, &format!("variant_summary.md5_{summary_md5}.txt.gz")),
            Some(&summary_md5),
        )?,
    );
    let mane_name = mane_url.rsplit('/').next().unwrap_or_default();
    sources.set(
        "mane",
        acquire_file(&mane_url, &join(&folder, mane_name), None)?,
    );
    require(
        publisher_md5(&summary_url)? == summary_md5,
        "ClinVar variant_summary publisher MD5 changed during acquisition; no statistics selected",
    )?;
    let after = latest_metadata(
        &text_of(item(plan, "assembly")?),
        &text_of(item(policy, "latest_url_template")?),
    )?;
    require(
        json_text::py_equal(item(&after, "file_date")?, item(latest, "file_date")?)
            && json_text::py_equal(item(&after, "publisher_md5")?, item(latest, "publisher_md5")?),
        "Official ClinVar changed during acquisition; no partial release was selected",
    )?;
    let release = mane_url
        .split_once("/release_")
        .and_then(|(_, rest)| rest.split_once('/'))
        .filter(|(name, _)| !name.is_empty())
        .map_or_else(|| mane_url.clone(), |(name, _)| name.to_owned());
    let value = object([
        ("kind", string(RECEIPT_KIND)),
        ("sources", sources),
        ("latest", latest.clone()),
        ("acquired_utc", string(utc_now())),
        ("mane_release", string(release)),
        (
            "summary_release_identity",
            string("independent publisher MD5 snapshot; not inferred from filename or VCF date"),
        ),
    ]);
    write_json(&join(&build, "inputs.ready.json"), &value)?;
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn date_of(text: &str) -> std::result::Result<Date, String> {
        header_date(&mut text.as_bytes(), "x").map_err(|failure| exception(failure.as_ref()).1)
    }

    #[test]
    fn file_dates_come_from_the_bounded_header() {
        assert_eq!(date_of("##fileformat=VCFv4.1\n##fileDate=2026-10-27\n#CHROM\tPOS\n"), Ok((2026, 10, 27)));
        assert_eq!(date_of("##fileDate=20261027\r\n#CHROM\tPOS\n"), Ok((2026, 10, 27)));
        assert_eq!(
            date_of("##fileDate=2026-10-27\n##fileDate=20261028\n#CHROM\t\n"),
            Err("Conflicting ClinVar fileDate header values".into())
        );
        assert_eq!(
            date_of("##x\n#CHROM\t\n"),
            Err("ClinVar VCF header has no fileDate data date".into())
        );
        assert_eq!(
            date_of("1\t2\n"),
            Err("Expected a ClinVar VCF header before variant records".into())
        );
        assert_eq!(
            date_of("##fileDate=2026/10/27\n"),
            Err("ClinVar fileDate must be a YYYYMMDD or YYYY-MM-DD data date".into())
        );
    }

    #[test]
    fn stamps_parse_like_strptime() {
        let message = |stamp| stamp_date(stamp).map_err(|failure| exception(failure.as_ref()).1);
        assert_eq!(message("2026-10-27"), Ok((2026, 10, 27)));
        assert_eq!(message("20261027"), Ok((2026, 10, 27)));
        assert_eq!(message("2026-02-30"), Err("day is out of range for month".into()));
        assert_eq!(
            message("2026-13-01"),
            Err("time data '2026-13-01' does not match format '%Y-%m-%d'".into())
        );
        assert_eq!(message("20261301"), Err("unconverted data remains: 1".into()));
        assert_eq!(message("2026-10-32"), Err("unconverted data remains: 2".into()));
        assert_eq!(message("0000-01-01"), Err("year 0 is out of range".into()));
    }

    #[test]
    fn templates_format_like_str_format() {
        assert_eq!(format_template("a/{assembly}/b", "GRCh38").unwrap(), "a/GRCh38/b");
        assert_eq!(format_template("{{x}}{assembly}", "G").unwrap(), "{x}G");
        assert!(format_template("a}{assembly}", "G").is_err());
    }

    #[test]
    fn redirects_are_followed_to_the_final_reply() {
        let mut text = &b"HTTP/1.1 302 Found\r\nLocation: /b\r\n\r\nHTTP/1.1 404 Not Found\r\nETag: x\r\n\r\n"[..];
        let (status, reason, headers) = read_head(&mut text).unwrap().unwrap();
        assert_eq!((status, reason.as_str()), (404, "Not Found"));
        assert_eq!(headers, [("ETag".to_owned(), "x".to_owned())]);
    }
}
