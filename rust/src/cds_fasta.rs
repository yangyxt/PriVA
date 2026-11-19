//! `priva cds-fasta validate --fasta FASTA --assembly ASSEMBLY`: the structure
//! check of an Ensembl CDS FASTA, which `CDS_FASTA_install` ran as an inline
//! Python block (`validate_cds_fasta` in `install_utils.sh` before this
//! command). It reads the file as `AlternativeStartCodonDetector` does, plain
//! or gzip by the `.gz` name, and checks only that:
//!
//! - the first non-blank line is a header (`>` and some text after it);
//! - every header is followed by at least one sequence line;
//! - a header that names `GRCh37` or `GRCh38` (as a whole word) names only the
//!   configured assembly.
//!
//! Lines are split as Python's text mode splits them (`\n`, `\r\n` or `\r`)
//! and stripped of Python's white space. Success prints
//! `Validated CDS FASTA structure: PATH; records=N`; a failure prints the
//! reason and exits with 1. Compared with the Python block by
//! `rust/tests/install_checks_compare.py`.
use crate::cli::{Result, invalid};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufReader, Read};

const USAGE: &str = "usage: priva cds-fasta validate --fasta FASTA --assembly {hg19,GRCh37,hg38,GRCh38}\n\
Check the structure of a CDS FASTA (plain, or gzip when the name ends in .gz) before the installer records it.";

/// `str.isspace()`: Unicode white space and U+001C to U+001F.
fn is_space(c: char) -> bool {
    c.is_whitespace() || ('\u{1c}'..='\u{1f}').contains(&c)
}

/// Python's `\w` for one character of a `str` pattern.
fn is_word(c: char) -> bool {
    c.is_alphanumeric() || c == '_'
}

/// The assemblies a header names: `\bGRCh(?:37|38)\b`.
fn declared(line: &str) -> Vec<&'static str> {
    let mut found = Vec::new();
    let chars: Vec<(usize, char)> = line.char_indices().collect();
    for (at, &(start, _)) in chars.iter().enumerate() {
        if !line[start..].starts_with("GRCh") {
            continue;
        }
        if at > 0 && is_word(chars[at - 1].1) {
            continue;
        }
        let digits = &line[start + 4..];
        let name = if digits.starts_with("37") {
            "GRCh37"
        } else if digits.starts_with("38") {
            "GRCh38"
        } else {
            continue;
        };
        let after = digits[2..].chars().next();
        if after.is_some_and(is_word) {
            continue;
        }
        if !found.contains(&name) {
            found.push(name);
        }
    }
    found
}

/// The text of the file as Python's `open(..., "rt", encoding="utf-8")` gives it.
fn read_text(path: &str) -> std::result::Result<String, String> {
    let file = File::open(path).map_err(|error| format!("{error}: {path}"))?;
    let mut bytes = Vec::new();
    let read = if path.ends_with(".gz") {
        MultiGzDecoder::new(BufReader::new(file)).read_to_end(&mut bytes)
    } else {
        BufReader::new(file).read_to_end(&mut bytes)
    };
    read.map_err(|error| format!("Cannot read CDS FASTA {path}: {error}"))?;
    String::from_utf8(bytes)
        .map_err(|error| format!("CDS FASTA is not UTF-8 text: {path}: {error}"))
}

/// Lines as Python's universal newlines split them.
fn lines(text: &str) -> impl Iterator<Item = &str> {
    text.split('\n').flat_map(|piece| piece.split('\r'))
}

pub(crate) fn validate(path: &str, assembly: &str) -> std::result::Result<usize, String> {
    let assembly = match assembly {
        "hg19" | "GRCh37" => "GRCh37",
        "hg38" | "GRCh38" => "GRCh38",
        other => return Err(format!("Unsupported assembly for the CDS FASTA: {other}")),
    };
    let text = read_text(path)?;
    let mut records = 0usize;
    let mut has_sequence = false;
    for line in lines(&text) {
        let line = line.trim_matches(is_space);
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if records > 0 && !has_sequence {
                return Err(format!("CDS FASTA contains an empty sequence: {path}"));
            }
            if header.trim_matches(is_space).is_empty() {
                return Err(format!("CDS FASTA contains an empty header: {path}"));
            }
            let named = declared(line);
            if !named.is_empty() && named != [assembly] {
                return Err(format!(
                    "CDS FASTA header conflicts with {assembly}: {path}"
                ));
            }
            records += 1;
            has_sequence = false;
        } else {
            if records == 0 {
                return Err(format!("Not a CDS FASTA: {path}"));
            }
            has_sequence = true;
        }
    }
    if records == 0 || !has_sequence {
        return Err(format!("CDS FASTA has no complete sequence: {path}"));
    }
    Ok(records)
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args
        .first()
        .is_some_and(|arg| arg == "--help" || arg == "-h")
    {
        println!("{USAGE}");
        return Ok(());
    }
    let (mut fasta, mut assembly) = (None, None);
    match args.first().map(String::as_str) {
        Some("validate") => {}
        _ => return Err(invalid(USAGE)),
    }
    let mut rest = args[1..].iter();
    while let Some(arg) = rest.next() {
        let slot = match arg.as_str() {
            "--fasta" => &mut fasta,
            "--assembly" => &mut assembly,
            "--help" | "-h" => {
                println!("{USAGE}");
                return Ok(());
            }
            _ => return Err(invalid(format!("unrecognized argument {arg}\n{USAGE}"))),
        };
        *slot = Some(
            rest.next()
                .ok_or_else(|| invalid(format!("{arg} needs a value")))?
                .clone(),
        );
    }
    let (Some(fasta), Some(assembly)) = (fasta, assembly) else {
        return Err(invalid(USAGE));
    };
    match validate(&fasta, &assembly) {
        Ok(records) => {
            println!("Validated CDS FASTA structure: {fasta}; records={records}");
            Ok(())
        }
        Err(message) => {
            eprintln!("priva: {message}");
            std::process::exit(1)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(name: &str, text: &[u8]) -> String {
        let folder = std::env::temp_dir().join(format!("priva-cds-{}", std::process::id()));
        std::fs::create_dir_all(&folder).unwrap();
        let path = folder.join(name);
        std::fs::write(&path, text).unwrap();
        path.to_string_lossy().into_owned()
    }

    #[test]
    fn headers_and_sequences() {
        let good = write(
            "good.fa",
            b">ENST1 cds chromosome:GRCh38:1\nATG\n\n>ENST2\r\nATGC\r\n",
        );
        assert_eq!(validate(&good, "hg38"), Ok(2));
        assert!(
            validate(&good, "hg19")
                .unwrap_err()
                .contains("conflicts with GRCh37")
        );
        let empty = write("empty.fa", b">ENST1\n>ENST2\nATG\n");
        assert!(
            validate(&empty, "GRCh38")
                .unwrap_err()
                .contains("empty sequence")
        );
        let bare = write("bare.fa", b">  \nATG\n");
        assert!(
            validate(&bare, "GRCh38")
                .unwrap_err()
                .contains("empty header")
        );
        let first = write("first.fa", b"ATG\n>ENST1\nATG\n");
        assert!(
            validate(&first, "GRCh38")
                .unwrap_err()
                .contains("Not a CDS FASTA")
        );
        let last = write("last.fa", b">ENST1\nATG\n>ENST2\n");
        assert!(
            validate(&last, "GRCh38")
                .unwrap_err()
                .contains("no complete sequence")
        );
    }

    #[test]
    fn assembly_words() {
        assert_eq!(declared(">x chromosome:GRCh38:1"), vec!["GRCh38"]);
        assert!(declared(">x GRCh380").is_empty());
        assert!(declared(">x aGRCh38").is_empty());
        assert_eq!(declared(">x GRCh37 GRCh38"), vec!["GRCh37", "GRCh38"]);
        // A lone carriage return splits a line, as in Python's text mode.
        let split = write("cr.fa", b">ENST1\rATG\r");
        assert_eq!(validate(&split, "GRCh38"), Ok(1));
    }
}
