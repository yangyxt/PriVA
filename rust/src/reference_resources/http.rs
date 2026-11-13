//! The one download of a pinned UCSC archive, run with `curl`.
//!
//! urlopen's rules are kept: redirects are followed; a final status outside
//! 2xx is the script's "UCSC download failed: HTTP <code>"; a 206 reply must
//! carry exactly the requested Content-Range and is appended to the partial;
//! a 200 reply (the server ignored Range) restarts it; any other status is
//! refused before anything is written. The headers are read from curl's
//! `--include` output, so the body streams into the partial as it arrives and
//! an interrupted transfer keeps the bytes already received.
use super::{Result, error, require, value_error};
use crate::installer_prerequisites::paths::read_failure;
use crate::python_path::os_error;
use std::fs::OpenOptions;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::process::{Child, Command, Stdio};

/// urlopen(timeout=60) waits at most a minute for a connection or for data.
const CURL_OPTIONS: [&str; 11] = [
    "--silent",
    "--show-error",
    "--location",
    "--include",
    "--suppress-connect-headers",
    "--connect-timeout",
    "60",
    "--speed-limit",
    "1",
    "--speed-time",
    "60",
];

/// The final response's status and its first Content-Range value.
#[derive(Debug, PartialEq, Eq)]
pub(super) struct Head {
    pub(super) status: u16,
    pub(super) content_range: Option<String>,
}

fn status_code(line: &str) -> Option<u16> {
    let mut words = line.split_ascii_whitespace();
    if !words.next()?.starts_with("HTTP/") {
        return None;
    }
    let code = words.next()?;
    if code.len() != 3 {
        return None;
    }
    code.parse().ok()
}

/// The head of the response curl settled on: informational replies and the
/// redirects it followed are skipped. `None` when the output ends first.
pub(super) fn read_head(reader: &mut impl BufRead) -> io::Result<Option<Head>> {
    loop {
        let mut line = Vec::new();
        if reader.read_until(b'\n', &mut line)? == 0 {
            return Ok(None);
        }
        let text = String::from_utf8_lossy(&line);
        let Some(status) = status_code(&text) else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("no HTTP status line: {}", text.trim_end()),
            ));
        };
        let mut content_range = None;
        let mut location = false;
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
                let value = value.trim_start_matches([' ', '\t']);
                if name.eq_ignore_ascii_case("content-range") && content_range.is_none() {
                    content_range = Some(value.to_owned());
                }
                location |= name.eq_ignore_ascii_case("location");
            }
        }
        if (100..200).contains(&status) {
            continue;
        }
        if (300..400).contains(&status) && location && reader.fill_buf()?.starts_with(b"HTTP/") {
            continue;
        }
        return Ok(Some(Head {
            status,
            content_range,
        }));
    }
}

fn url_error(message: impl std::fmt::Display) -> Box<dyn std::error::Error> {
    error("URLError", format!("<urlopen error {message}>"))
}

/// Fetches `url` from byte `offset` into `partial`, whose complete length is
/// `total`.
pub(super) fn fetch(url: &str, total: u64, offset: u64, partial: &str) -> Result<()> {
    let mut command = Command::new("curl");
    command.args(CURL_OPTIONS);
    if offset > 0 {
        let range = format!("Range: bytes={offset}-");
        command.args(["--header", range.as_str()]);
    }
    command
        .arg(url)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    let spawned = command.spawn();
    let mut child = spawned.map_err(|failure| os_error(failure, Path::new("curl"), None))?;
    let received = transfer(&mut child, url, total, offset, partial);
    if received.is_err() {
        // Stop the transfer; the partial keeps what was already written.
        let _ = child.kill();
    }
    let finished = child.wait_with_output();
    let output = finished.map_err(|failure| os_error(failure, Path::new("curl"), None))?;
    let responded = received?;
    if output.status.success() && responded {
        return Ok(());
    }
    let reason = String::from_utf8_lossy(&output.stderr).trim().to_owned();
    if reason.is_empty() {
        let ending = if responded {
            ""
        } else {
            " before an HTTP response"
        };
        return Err(url_error(format!(
            "curl ended with {}{ending}",
            output.status
        )));
    }
    Err(url_error(reason))
}

/// Reads the response head and streams the body; false without a response.
fn transfer(child: &mut Child, url: &str, total: u64, offset: u64, partial: &str) -> Result<bool> {
    let mut body = BufReader::new(child.stdout.take().expect("curl output is piped"));
    let Some(head) = read_head(&mut body).map_err(url_error)? else {
        return Ok(false);
    };
    if !(200..300).contains(&head.status) {
        return Err(value_error(format!(
            "UCSC download failed: HTTP {}; resumable partial retained: {partial}",
            head.status
        )));
    }
    let append = head.status == 206;
    if append {
        let range = head.content_range.unwrap_or_default();
        require(
            range == format!("bytes {offset}-{}/{total}", total - 1),
            format!("Incorrect UCSC download Content-Range: {range}; partial retained"),
        )?;
    } else {
        require(
            head.status == 200,
            format!("Unexpected UCSC HTTP status: {}", head.status),
        )?;
    }
    eprintln!("Downloading pinned UCSC resource: {url} -> {partial}");
    // A 200 reply ignored Range: restart this staged download.
    let mut file = OpenOptions::new()
        .create(true)
        .write(true)
        .append(append)
        .truncate(!append)
        .open(partial)
        .map_err(|failure| os_error(failure, Path::new(partial), None))?;
    io::copy(&mut body, &mut file).map_err(|failure| read_failure(Box::new(failure), partial))?;
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn head(text: &str) -> Option<Head> {
        read_head(&mut text.as_bytes()).unwrap()
    }

    #[test]
    fn followed_redirects_and_informational_replies_are_skipped() {
        let found = head(
            "HTTP/1.1 302 Found\r\nLocation: /b\r\n\r\nHTTP/1.1 100 Continue\r\n\r\n\
             HTTP/2 206 \r\ncontent-range: bytes 5-9/10\r\nContent-Range: other\r\n\r\nbody",
        );
        let expected = Head {
            status: 206,
            content_range: Some("bytes 5-9/10".to_owned()),
        };
        assert_eq!(found, Some(expected));
    }

    #[test]
    fn a_redirect_curl_did_not_follow_is_final() {
        let found = head("HTTP/1.1 301 Moved\r\nLocation: /b\r\n\r\n<html>");
        assert_eq!(found.map(|head| head.status), Some(301));
        let found = head("HTTP/1.1 404 Not Found\r\n\r\n");
        assert_eq!(found.map(|head| head.status), Some(404));
    }

    #[test]
    fn missing_or_cut_heads_are_not_responses() {
        assert_eq!(head(""), None);
        assert_eq!(head("HTTP/1.1 200 OK\r\nContent-Length: 3\r\n"), None);
        assert!(read_head(&mut "garbage\r\n".as_bytes()).is_err());
    }
}
