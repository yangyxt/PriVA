//! One `urlopen` request, run with `curl`.
//!
//! As in rust/src/reference_resources/http.rs, the response head arrives on
//! curl's `--include` output and the body streams on after it, so a caller can
//! check the head before it touches a file and an interrupted transfer keeps
//! the bytes already received. curl follows no redirect here: `open` follows
//! them with urllib's rules, because the caller needs the final URL and each
//! hop's refusal. A 4xx/5xx reply is `HTTPError`. A cut body ends like urllib's
//! `read(amt)`: the caller sees end of data, except for a chunked body, which
//! is `IncompleteRead`. Transport failures carry curl's message.
use super::*;
use crate::builder_error::RuntimeError;
use std::io::{self, BufRead, BufReader};
use std::process::{Child, ChildStdout, Command, Stdio};

pub(super) const BLOCK: usize = 1024 * 1024;
const MAX_REDIRECTIONS: usize = 10;
const MAX_REPEATS: u32 = 4;

/// The exception classes under OSError that this acquisition can raise.
const OS_ERRORS: [&str; 14] = [
    "OSError",
    "FileNotFoundError",
    "PermissionError",
    "FileExistsError",
    "IsADirectoryError",
    "NotADirectoryError",
    "BlockingIOError",
    "BrokenPipeError",
    "BadGzipFile",
    "HTTPError",
    "URLError",
    "TimeoutError",
    "ConnectionError",
    "ConnectionResetError",
];
const VALUE_ERRORS: [&str; 4] = [
    "ValueError",
    "UnicodeDecodeError",
    "UnicodeEncodeError",
    "JSONDecodeError",
];

fn category(failure: &(dyn std::error::Error + 'static)) -> &str {
    failure
        .downcast_ref::<RuntimeError>()
        .map_or("", RuntimeError::category_name)
}
pub(super) fn is_os_error(failure: &(dyn std::error::Error + 'static)) -> bool {
    OS_ERRORS.contains(&category(failure))
}
pub(super) fn is_value_error(failure: &(dyn std::error::Error + 'static)) -> bool {
    VALUE_ERRORS.contains(&category(failure))
}
pub(super) fn is_incomplete_read(failure: &(dyn std::error::Error + 'static)) -> bool {
    category(failure) == "IncompleteRead"
}
/// `str()` of the exception.
pub(super) fn message(failure: &(dyn std::error::Error + 'static)) -> String {
    failure
        .downcast_ref::<RuntimeError>()
        .map_or_else(|| failure.to_string(), RuntimeError::message)
}

pub(super) struct Response {
    child: Child,
    body: BufReader<ChildStdout>,
    pub(super) status: u16,
    reason: String,
    pub(super) url: String,
    headers: Vec<(String, String)>,
}

impl Drop for Response {
    /// Closing the response stops whatever curl has not delivered.
    fn drop(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

impl Response {
    /// The first value of a header, as `headers.get(name)`.
    pub(super) fn header(&self, name: &str) -> Option<&str> {
        self.headers
            .iter()
            .find(|(key, _)| key.eq_ignore_ascii_case(name))
            .map(|(_, value)| value.as_str())
    }

    /// `read(1024 * 1024)`: a full block unless the body ends first. Check
    /// `finish` after the empty block, where a transport failure surfaces.
    pub(super) fn read_block(&mut self) -> Result<Vec<u8>> {
        let mut block = Vec::with_capacity(BLOCK);
        io::Read::read_to_end(
            &mut io::Read::take(&mut self.body, BLOCK as u64),
            &mut block,
        )
        .map_err(|e| error("OSError", e.to_string()))?;
        Ok(block)
    }

    /// How the transfer ended, once the body has been read to its end.
    pub(super) fn finish(&mut self) -> Result<()> {
        let status = self
            .child
            .wait()
            .map_err(|e| os_error(e, Path::new("curl"), None))?;
        if status.success() {
            return Ok(());
        }
        let reason = curl_reason(&mut self.child);
        let chunked = self
            .header("Transfer-Encoding")
            .is_some_and(|v| v.to_ascii_lowercase().contains("chunked"));
        // Exit 18 is a body shorter than announced; 56 with this text is a TLS
        // peer that closed without notice, which Python reads as end of data.
        let abrupt = reason.contains("unexpected eof") || reason.ends_with("errno 0");
        let cut = status.code() == Some(18) || (status.code() == Some(56) && abrupt);
        if cut && chunked {
            return Err(error("IncompleteRead", reason));
        }
        if cut {
            return Ok(());
        }
        Err(transport_failure(status.code(), reason))
    }
}

fn curl_reason(child: &mut Child) -> String {
    let mut text = String::new();
    if let Some(stderr) = child.stderr.as_mut() {
        let _ = io::Read::read_to_string(stderr, &mut text);
    }
    crate::python_text::python_strip(&text).to_owned()
}

fn transport_failure(code: Option<i32>, reason: String) -> Box<dyn std::error::Error> {
    match code {
        Some(28) => error("TimeoutError", reason),
        Some(56) if reason.contains("reset") => error("ConnectionResetError", reason),
        _ => error("URLError", format!("<urlopen error {reason}>")),
    }
}

struct Head {
    status: u16,
    reason: String,
    headers: Vec<(String, String)>,
}

/// The response head, past any informational reply. `None` when curl's output
/// ends first. Without --location this is the only final head on the stream.
fn read_head(reader: &mut impl BufRead) -> io::Result<Option<Head>> {
    loop {
        let mut line = Vec::new();
        if reader.read_until(b'\n', &mut line)? == 0 {
            return Ok(None);
        }
        let first = String::from_utf8_lossy(&line).trim_end().to_owned();
        let mut words = first.splitn(3, ' ');
        let version = words.next().unwrap_or("");
        let status = words.next().and_then(|code| code.parse::<u16>().ok());
        let Some(status) = status.filter(|_| version.starts_with("HTTP/")) else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("no HTTP status line: {first}"),
            ));
        };
        let reason = words.next().unwrap_or("").to_owned();
        let mut headers = Vec::new();
        loop {
            line.clear();
            if reader.read_until(b'\n', &mut line)? == 0 {
                return Ok(None);
            }
            let text = String::from_utf8_lossy(&line);
            let field = text.trim_end_matches(['\r', '\n']);
            if field.is_empty() {
                break;
            }
            if let Some((name, value)) = field.split_once(':') {
                headers.push((name.to_owned(), value.trim().to_owned()));
            }
        }
        if status >= 200 {
            return Ok(Some(Head {
                status,
                reason,
                headers,
            }));
        }
    }
}

/// One curl process for one URL, stopped at the end of the response head.
fn request(url: &str, headers: &[(&str, &str)], timeout: f64) -> Result<Response> {
    let mut command = Command::new("curl");
    command.args([
        "--silent",
        "--show-error",
        "--http1.1",
        "--include",
        "--suppress-connect-headers",
    ]);
    if timeout.is_finite() {
        // urlopen's timeout bounds each connection and each wait for data.
        let seconds = timeout.ceil().max(1.0).to_string();
        command.args(["--connect-timeout", &timeout.to_string()]);
        command.args(["--speed-limit", "1", "--speed-time", &seconds]);
    }
    for (name, value) in headers {
        command.args(["--header", &format!("{name}: {value}")]);
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
    match read_head(&mut body) {
        Ok(Some(head)) => Ok(Response {
            child,
            body,
            status: head.status,
            reason: head.reason,
            url: url.to_owned(),
            headers: head.headers,
        }),
        Ok(None) => {
            let status = child
                .wait()
                .map_err(|e| os_error(e, Path::new("curl"), None))?;
            let reason = curl_reason(&mut child);
            Err(transport_failure(status.code(), reason))
        }
        Err(cause) => {
            let _ = child.kill();
            let _ = child.wait();
            Err(error("OSError", cause.to_string()))
        }
    }
}

/// `urljoin(base, location)` for the redirect forms a server sends.
fn join_url(base: &str, location: &str) -> String {
    if location.contains("://") {
        return location.to_owned();
    }
    let (scheme, rest) = base.split_once("://").unwrap_or(("https", base));
    if let Some(network) = location.strip_prefix("//") {
        return format!("{scheme}://{network}");
    }
    let authority = rest.split(['/', '?', '#']).next().unwrap_or(rest);
    if location.starts_with('/') {
        return format!("{scheme}://{authority}{location}");
    }
    let path = rest[authority.len()..].split(['?', '#']).next().unwrap_or("");
    let folder = path.rfind('/').map_or("/", |end| &path[..=end]);
    format!("{scheme}://{authority}{folder}{location}")
}

/// `urlparse(url).hostname`.
pub(super) fn hostname(url: &str) -> Option<String> {
    let (_, rest) = url.split_once("://")?;
    let authority = rest.split(['/', '?', '#']).next().unwrap_or(rest);
    let host = authority.rsplit('@').next().unwrap_or(authority);
    let host = match host.strip_prefix('[') {
        Some(tail) => tail.split(']').next().unwrap_or(tail),
        None => host.split(':').next().unwrap_or(host),
    };
    (!host.is_empty()).then(|| host.to_ascii_lowercase())
}

/// `urlopen(Request(url, headers=headers), timeout=timeout)`.
pub(super) fn open(url: &str, headers: &[(&str, &str)], timeout: f64) -> Result<Response> {
    let mut current = url.to_owned();
    let mut visited: Vec<(String, u32)> = Vec::new();
    loop {
        let response = request(&current, headers, timeout)?;
        let status = response.status;
        if (200..300).contains(&status) {
            return Ok(response);
        }
        let http_error = |text: &str| error("HTTPError", format!("HTTP Error {status}: {text}"));
        let location = response.header("Location").or(response.header("URI"));
        let Some(location) = location.filter(|_| matches!(status, 301 | 302 | 303 | 307 | 308))
        else {
            return Err(http_error(&response.reason));
        };
        if let Some((scheme, _)) = location.split_once("://")
            && !matches!(scheme.to_ascii_lowercase().as_str(), "http" | "https")
        {
            return Err(http_error(&format!(
                "{} - Redirection to url {} is not allowed",
                response.reason,
                repr(location)
            )));
        }
        let next = join_url(&current, location);
        let seen = visited.iter().position(|(url, _)| *url == next);
        if seen.is_some_and(|at| visited[at].1 >= MAX_REPEATS)
            || visited.len() >= MAX_REDIRECTIONS
        {
            return Err(http_error(&format!(
                "The HTTP server returned a redirect error that would lead to an infinite loop.\n\
                 The last 30x error message was:\n{}",
                response.reason
            )));
        }
        match seen {
            Some(at) => visited[at].1 += 1,
            None => visited.push((next.clone(), 1)),
        }
        current = next;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn head(text: &str) -> Option<Head> {
        read_head(&mut text.as_bytes()).unwrap()
    }

    #[test]
    fn the_head_skips_informational_replies_and_keeps_the_reason() {
        let found = head(
            "HTTP/1.1 100 Continue\r\n\r\n\
             HTTP/1.1 503 Service Unavailable\r\nContent-Length: 3\r\n\r\nbody",
        )
        .unwrap();
        assert_eq!(found.status, 503);
        assert_eq!(found.reason, "Service Unavailable");
        assert_eq!(found.headers[0].1, "3");
        assert!(head("HTTP/1.1 200 OK\r\nContent-Length: 3\r\n").is_none());
        assert!(read_head(&mut "garbage\r\n".as_bytes()).is_err());
    }

    #[test]
    fn redirect_locations_join_like_urljoin() {
        let base = "https://a.example/x/y?query=1";
        assert_eq!(join_url(base, "https://b.example/z"), "https://b.example/z");
        assert_eq!(join_url(base, "//b.example/z"), "https://b.example/z");
        assert_eq!(join_url(base, "/z"), "https://a.example/z");
        assert_eq!(join_url(base, "z"), "https://a.example/x/z");
    }

    #[test]
    fn hostnames_ignore_port_userinfo_and_case() {
        assert_eq!(
            hostname("https://u:p@GRCH37.ensembl.org:443/biomart").unwrap(),
            "grch37.ensembl.org"
        );
        assert_eq!(hostname("https://[::1]:8/").unwrap(), "::1");
        assert_eq!(hostname("not a url"), None);
    }

    #[test]
    fn failure_families_follow_the_python_exception_classes() {
        for (failure, os, value) in [
            (error("HTTPError", "HTTP Error 503: x"), true, false),
            (error("BadGzipFile", "x"), true, false),
            (error("ValueError", "x"), false, true),
            (error("UnicodeDecodeError", "x"), false, true),
            (error("IncompleteRead", "x"), false, false),
            (error("EOFError", "x"), false, false),
        ] {
            assert_eq!(is_os_error(failure.as_ref()), os);
            assert_eq!(is_value_error(failure.as_ref()), value);
        }
        assert_eq!(
            message(error("ValueError", "text").as_ref()),
            "text".to_owned()
        );
    }
}
