//! CPython 3.12 urlsplit checks used by source_url; no network operation.
use super::*;
use crate::python_text::is_python_space;

pub(super) fn check(value: &str) -> Result<()> {
    let cleaned = value
        .trim_start_matches(|c: char| c <= ' ')
        .replace(['\t', '\r', '\n'], "");
    let mut rest = cleaned.as_str();
    let mut scheme = String::new();
    if let Some((head, tail)) = rest.split_once(':') {
        if head.starts_with(|c: char| c.is_ascii_alphabetic())
            && head
                .chars()
                .all(|c| c.is_ascii_alphanumeric() || "+-.".contains(c))
        {
            scheme = head.to_ascii_lowercase();
            rest = tail;
        }
    }
    let netloc = rest
        .strip_prefix("//")
        .map(|tail| tail.split(['/', '?', '#']).next().unwrap())
        .unwrap_or("");
    if netloc.contains('[') != netloc.contains(']') {
        return Err(error("ValueError", "Invalid IPv6 URL"));
    }
    if let Some((_, tail)) = netloc.split_once('[') {
        let host = tail.split(']').next().unwrap();
        if let Some(tail) = host.strip_prefix('v') {
            if !tail.split_once('.').is_some_and(|(v, body)| {
                !v.is_empty()
                    && v.bytes().all(|b| b.is_ascii_hexdigit())
                    && !body.is_empty()
                    && !body.contains('\n')
            }) {
                return Err(error("ValueError", "IPvFuture address is invalid"));
            }
        } else {
            if host.parse::<std::net::Ipv4Addr>().is_ok() {
                return Err(error("ValueError", "An IPv4 address cannot be in brackets"));
            }
            let (address, scope_ok) = match host.split_once('%') {
                Some((a, s)) => (a, !s.is_empty() && !s.contains('%')),
                None => (host, true),
            };
            if !scope_ok || address.parse::<std::net::Ipv6Addr>().is_err() {
                return Err(error(
                    "ValueError",
                    format!(
                        "{} does not appear to be an IPv4 or IPv6 address",
                        crate::python_path::repr(host)
                    ),
                ));
            }
        }
    }
    // Exhaustively derived from CPython's Unicode 15.0 NFKC database: these
    // non-ASCII characters expand to a URL delimiter. Normalization does not
    // combine other characters into any of the five ASCII delimiters.
    if netloc.chars().any(|c| {
        matches!(
            c as u32,
            8263 | 8264
                | 8265
                | 8448
                | 8449
                | 8453
                | 8454
                | 10868
                | 65043
                | 65046
                | 65109
                | 65110
                | 65119
                | 65131
                | 65283
                | 65295
                | 65306
                | 65311
                | 65312
        )
    }) {
        return Err(error(
            "ValueError",
            format!("netloc '{netloc}' contains invalid characters under NFKC normalization"),
        ));
    }
    if !matches!(scheme.as_str(), "http" | "https")
        || netloc.is_empty()
        || value.chars().any(is_python_space)
    {
        return Err(error(
            "ValueError",
            "Source URL must be an absolute HTTP(S) URL",
        ));
    }
    Ok(())
}
