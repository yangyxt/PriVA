//! The SpliceAI cutoffs that `splicing_var_analysis.py --config` takes from
//! `spliceai_cutoffs.cutoffs_for_config`.
//!
//! PRIVATE COPY for `priva splicing-analysis`. `spliceai_cutoffs.py` and
//! `spliceai_cache_input.py` (whose `inspect_header` and `canonical_assembly`
//! are reproduced here) belong to another conversion thread; when that work is
//! merged this module should be replaced by theirs.
//!
//! The version of the SpliceAI scores comes from the headers of the two
//! prescore files named in the config (`spliceai_snv_prescore`,
//! `spliceai_indel_prescore`): both must be the same version and agree with
//! the configured `assembly`. SpliceAI v1 uses the cutoffs 0.8 and 0.5,
//! SpliceAI2 uses 0.5 and 0.25. Without a config, or with neither file set,
//! the v1 cutoffs apply.
//!
//! Differences from Python, all failing loudly: the config is read as a flat
//! list of top-level `key: value` lines (nested sections are skipped over, and
//! a document marker, a sequence, flow style, a tab indent or a quoted key at
//! the left margin is refused), and the three keys must hold plain or quoted
//! text or null. Python's check that the file did not change while its header
//! was read is not repeated.
use crate::cli::{Result, invalid};
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read};
use std::path::{Component, Path};

const MAX_LINE: usize = 64 * 1024;
const MAX_HEADER: usize = 1024 * 1024;
const MAX_LINES: usize = 4096;
const VCF_COLUMNS: [&str; 8] = [
    "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
];

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Cutoffs {
    pub(crate) splice_altering: f64,
    pub(crate) nearby_acceptor_rescue: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Format {
    VcfV1,
    VcfSpliceai2,
    Spliceai2Tsv,
}

struct Header {
    format: Format,
    assembly: Option<&'static str>,
}

fn canonical_assembly(assembly: &str) -> Result<&'static str> {
    match assembly.to_lowercase().as_str() {
        "hg19" | "grch37" => Ok("GRCh37"),
        "hg38" | "grch38" => Ok("GRCh38"),
        _ => Err(invalid(format!("Unsupported assembly: {assembly:?}"))),
    }
}

// ---------------------------------------------------------------- config

#[derive(Default)]
struct ConfigValues {
    snv: Option<String>,
    indel: Option<String>,
    assembly: Option<String>,
}

fn unsupported(why: &str) -> Box<dyn std::error::Error> {
    invalid(format!(
        "the config uses YAML this reader does not handle: {why}"
    ))
}

/// The text of a YAML scalar on one line: quoted, plain, or null.
fn scalar(raw: &str, key: &str) -> Result<Option<String>> {
    let raw = raw.trim();
    if let Some(rest) = raw.strip_prefix('"') {
        let mut value = String::new();
        let mut chars = rest.chars();
        loop {
            match chars.next() {
                None => return Err(unsupported("an unclosed double quote")),
                Some('"') => break,
                Some('\\') => match chars.next() {
                    Some('n') => value.push('\n'),
                    Some('t') => value.push('\t'),
                    Some('"') => value.push('"'),
                    Some('\\') => value.push('\\'),
                    Some('/') => value.push('/'),
                    _ => return Err(unsupported("an escape in a double-quoted value")),
                },
                Some(other) => value.push(other),
            }
        }
        let after: String = chars.collect();
        let after = after.trim();
        if !after.is_empty() && !after.starts_with('#') {
            return Err(unsupported("text after a quoted value"));
        }
        return Ok(Some(value));
    }
    if let Some(rest) = raw.strip_prefix('\'') {
        let mut value = String::new();
        let mut chars = rest.chars().peekable();
        loop {
            match chars.next() {
                None => return Err(unsupported("an unclosed single quote")),
                Some('\'') => {
                    if chars.peek() == Some(&'\'') {
                        chars.next();
                        value.push('\'');
                    } else {
                        break;
                    }
                }
                Some(other) => value.push(other),
            }
        }
        let after: String = chars.collect();
        let after = after.trim();
        if !after.is_empty() && !after.starts_with('#') {
            return Err(unsupported("text after a quoted value"));
        }
        return Ok(Some(value));
    }
    // A plain scalar ends at a comment: a hash after a blank.
    let plain = match raw.find(" #") {
        Some(position) => &raw[..position],
        None if raw.starts_with('#') => "",
        None => raw,
    }
    .trim_end();
    if matches!(plain, "" | "~" | "null" | "Null" | "NULL") {
        return Ok(None);
    }
    if plain.starts_with(['[', '{', '|', '>', '&', '*', '!', '%', '@', '`']) {
        return Err(unsupported(&format!("the value of {key}")));
    }
    let lower = plain.to_lowercase();
    if matches!(
        lower.as_str(),
        "true" | "false" | "yes" | "no" | "on" | "off" | "y" | "n"
    ) || (plain
        .bytes()
        .all(|b| b.is_ascii_digit() || b"+-._eE".contains(&b))
        && plain.bytes().any(|b| b.is_ascii_digit()))
    {
        return Err(unsupported(&format!("a non-text value for {key}")));
    }
    Ok(Some(plain.to_owned()))
}

fn read_config(path: Option<&str>) -> Result<ConfigValues> {
    let Some(path) = path.filter(|path| !path.is_empty() && Path::new(path).is_file()) else {
        return Ok(ConfigValues::default());
    };
    let text = fs::read_to_string(path)?;
    let lines: Vec<&str> = text.lines().collect();
    let mut values: HashMap<&str, Option<String>> = HashMap::new();
    for (number, line) in lines.iter().enumerate() {
        let first = line.chars().next();
        match first {
            None | Some('#') => continue,
            Some('\t') => return Err(unsupported("a tab at the start of a line")),
            Some(' ') => continue,
            Some(
                '-' | '?' | '[' | '{' | '&' | '*' | '!' | '|' | '>' | '\'' | '"' | '%' | '@' | '`',
            ) => {
                return Err(unsupported(&format!(
                    "line {} starts with an indicator",
                    number + 1
                )));
            }
            _ => {}
        }
        if line.starts_with("---") || line.starts_with("...") {
            return Err(unsupported("a document marker"));
        }
        let (key, rest) = match line.find(':') {
            Some(position)
                if position + 1 == line.len() || line[position + 1..].starts_with([' ', '\t']) =>
            {
                (&line[..position], &line[position + 1..])
            }
            _ => {
                return Err(unsupported(&format!(
                    "line {} is not a key: value line",
                    number + 1
                )));
            }
        };
        let key = key.trim_end();
        if matches!(
            key,
            "spliceai_snv_prescore" | "spliceai_indel_prescore" | "assembly"
        ) {
            let value = scalar(rest, key)?;
            // A key with no value whose next line is indented holds a nested section.
            if value.is_none() && rest.trim().is_empty() {
                let next = lines[number + 1..].iter().find(|candidate| {
                    !candidate.trim().is_empty() && !candidate.trim_start().starts_with('#')
                });
                if next.is_some_and(|candidate| candidate.starts_with([' ', '\t'])) {
                    return Err(unsupported(&format!("{key} holds a nested section")));
                }
            }
            values.insert(key_name(key), value);
        }
    }
    Ok(ConfigValues {
        snv: values.remove("spliceai_snv_prescore").flatten(),
        indel: values.remove("spliceai_indel_prescore").flatten(),
        assembly: values.remove("assembly").flatten(),
    })
}

fn key_name(key: &str) -> &'static str {
    match key {
        "spliceai_snv_prescore" => "spliceai_snv_prescore",
        "spliceai_indel_prescore" => "spliceai_indel_prescore",
        _ => "assembly",
    }
}

// ---------------------------------------------------------- header reading

fn native_columns() -> Vec<String> {
    let mut columns: Vec<String> = ["chrom", "pos", "ref", "alt", "strand"]
        .iter()
        .map(|name| (*name).to_owned())
        .collect();
    for effect in [
        "donor_gain",
        "donor_loss",
        "acceptor_gain",
        "acceptor_loss",
        "jxn_gain",
        "jxn_loss",
    ] {
        let mut fields = vec!["delta_score", "ref_score", "alt_score"];
        if effect.starts_with("jxn") {
            fields.extend(["donor_dist", "acceptor_dist"]);
        } else {
            fields.push("dist");
        }
        for field in fields {
            for index in 0..3 {
                columns.push(format!("{effect}_{field}_{index}"));
            }
        }
    }
    columns.push("spliceai2_summary_score".to_owned());
    columns
}

fn require(condition: bool, message: &str) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(invalid(message))
    }
}

fn unfinished(part: &str) -> bool {
    let part = part.to_lowercase();
    part.contains(".partial") || part == ".building" || part.starts_with(".building-")
}

fn complete_path(path: &str) -> Result<std::path::PathBuf> {
    let absolute = std::path::absolute(path)?;
    let named = |candidate: &Path| {
        candidate.components().any(|component| match component {
            Component::Normal(part) => unfinished(&part.to_string_lossy()),
            _ => false,
        })
    };
    require(
        !named(&absolute),
        &format!("Incomplete/staging input path: {path}"),
    )?;
    let resolved = fs::canonicalize(path)?;
    require(
        !named(&resolved),
        &format!("Incomplete/staging input target: {path}"),
    )?;
    let metadata = fs::metadata(&resolved)?;
    require(
        metadata.is_file() && metadata.len() > 0,
        &format!("Input is not a nonempty regular file: {path}"),
    )?;
    Ok(resolved)
}

fn assemblies_in(value: &str) -> Result<Vec<&'static str>> {
    let lower = value.to_ascii_lowercase();
    let bytes = lower.as_bytes();
    let mut found = Vec::new();
    let mut index = 0;
    while index < bytes.len() {
        let preceded = index > 0 && bytes[index - 1].is_ascii_alphanumeric();
        let mut matched = 0;
        if !preceded {
            for token in ["grch37", "grch38", "hg19", "hg38"] {
                if bytes[index..].starts_with(token.as_bytes()) {
                    let after = index + token.len();
                    if after == bytes.len() || !bytes[after].is_ascii_alphanumeric() {
                        found.push(canonical_assembly(token)?);
                        matched = token.len();
                        break;
                    }
                }
            }
        }
        index += matched.max(1);
    }
    Ok(found)
}

/// `re.findall(r'(?:^|,)assembly=("[^"\n]*"|[^,>]+)', text, re.I)`.
fn contig_assemblies(text: &[char]) -> Vec<String> {
    let mut found = Vec::new();
    let mut cursor = 0;
    let mut index = 0;
    while index <= text.len() {
        let mut starts = Vec::new();
        if index == 0 {
            starts.push(0);
        }
        if text.get(index) == Some(&',') {
            starts.push(index + 1);
        }
        let mut advanced = false;
        for start in starts {
            if start < cursor {
                continue;
            }
            let word: String = text[start..text.len().min(start + 9)].iter().collect();
            if !word.eq_ignore_ascii_case("assembly=") || word.chars().count() != 9 {
                continue;
            }
            let value_start = start + 9;
            let mut end = value_start;
            if text.get(value_start) == Some(&'"') {
                let mut close = value_start + 1;
                while close < text.len() && text[close] != '"' && text[close] != '\n' {
                    close += 1;
                }
                if close < text.len() && text[close] == '"' {
                    end = close + 1;
                }
            }
            if end == value_start {
                while end < text.len() && text[end] != ',' && text[end] != '>' {
                    end += 1;
                }
            }
            if end > value_start {
                found.push(text[value_start..end].iter().collect());
                cursor = end;
                index = end;
                advanced = true;
                break;
            }
        }
        if !advanced {
            index += 1;
        }
    }
    found
}

fn is_word(character: char) -> bool {
    character.is_alphanumeric() || character == '_'
}

/// `re.search(r"spliceai(?:2|[ _-]+v?2)(?:\b|_)", value, re.I)`.
fn mentions_spliceai2(value: &str) -> bool {
    let lower = value.to_ascii_lowercase();
    let bytes = lower.as_bytes();
    let boundary = |position: usize| -> bool {
        match lower[position..].chars().next() {
            None => true,
            Some(next) => next == '_' || !is_word(next),
        }
    };
    let mut from = 0;
    while let Some(found) = lower[from..].find("spliceai") {
        let start = from + found + 8;
        if bytes.get(start) == Some(&b'2') && boundary(start + 1) {
            return true;
        }
        let mut position = start;
        while position < bytes.len() && matches!(bytes[position], b' ' | b'_' | b'-') {
            position += 1;
        }
        if position > start {
            if bytes.get(position) == Some(&b'v') {
                position += 1;
            }
            if bytes.get(position) == Some(&b'2') && boundary(position + 1) {
                return true;
            }
        }
        from += found + 1;
    }
    false
}

/// `re.fullmatch(r"<first>(?:\.\d+)*", value)` where `first` is one of `firsts`.
fn is_dotted_version(value: &str, firsts: &[char]) -> bool {
    let mut parts = value.split('.');
    let head = parts.next().unwrap_or("");
    let mut head_chars = head.chars();
    let starts_right = head_chars
        .next()
        .is_some_and(|character| firsts.contains(&character))
        && head_chars.next().is_none();
    starts_right
        && parts.all(|part| !part.is_empty() && part.bytes().all(|byte| byte.is_ascii_digit()))
}

fn metadata(line: &str) -> Result<(String, String)> {
    let rest = line
        .strip_prefix("##")
        .ok_or_else(|| invalid("Invalid VCF metadata line"))?;
    let key_length = rest
        .bytes()
        .take_while(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'.' | b'-'))
        .count();
    require(
        key_length > 0
            && rest.as_bytes().get(key_length) == Some(&b'=')
            && rest.len() > key_length + 1,
        "Invalid VCF metadata line",
    )?;
    let key = rest[..key_length].to_owned();
    let value = rest[key_length + 1..].to_owned();
    // Structured metadata must close its brackets and quoted strings.
    if value.starts_with('<') {
        require(value.ends_with('>'), "Truncated structured VCF metadata")?;
        let inner: Vec<char> = value.chars().collect();
        let (mut quoted, mut escaped) = (false, false);
        for &character in &inner[1..inner.len() - 1] {
            if escaped {
                escaped = false;
            } else if character == '\\' && quoted {
                escaped = true;
            } else if character == '"' {
                quoted = !quoted;
            } else if (character == '<' || character == '>') && !quoted {
                return Err(invalid("Invalid structured VCF metadata brackets"));
            }
        }
        require(!quoted && !escaped, "Unclosed VCF metadata quotation")?;
    }
    Ok((key, value))
}

/// One header line of at most `MAX_LINE` bytes, with its line break, like
/// `readline(MAX_LINE + 1)`; an empty result is the end of the input.
fn read_line(reader: &mut dyn BufRead) -> Result<Vec<u8>> {
    let mut line: Vec<u8> = Vec::new();
    loop {
        let available = reader.fill_buf()?;
        if available.is_empty() {
            break;
        }
        let wanted = match available.iter().position(|byte| *byte == b'\n') {
            Some(position) => position + 1,
            None => available.len(),
        };
        let take = wanted.min(MAX_LINE + 1 - line.len());
        line.extend_from_slice(&available[..take]);
        reader.consume(take);
        if line.ends_with(b"\n") || line.len() > MAX_LINE {
            break;
        }
    }
    Ok(line)
}

fn inspect_header(path: &str) -> Result<Header> {
    let resolved = complete_path(path)?;
    let mut magic = [0_u8; 2];
    let read = File::open(&resolved)?.read(&mut magic)?;
    let file = File::open(&resolved)?;
    let mut reader: Box<dyn BufRead> = if read == 2 && magic == [0x1f, 0x8b] {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let mut assemblies: Vec<&'static str> = Vec::new();
    let mut versions: Vec<u8> = Vec::new();
    let mut total = 0;
    let mut result: Option<Header> = None;
    for number in 1..=MAX_LINES {
        let raw = read_line(reader.as_mut())?;
        require(!raw.is_empty(), "Missing complete VCF/SpliceAI2 header")?;
        total += raw.len();
        require(
            raw.len() <= MAX_LINE && total <= MAX_HEADER,
            "Cache header exceeds bounded size limit",
        )?;
        require(
            raw.ends_with(b"\n"),
            "Truncated cache header line (missing newline)",
        )?;
        let mut body = &raw[..raw.len() - 1];
        if body.ends_with(b"\r") {
            body = &body[..body.len() - 1];
        }
        let line = std::str::from_utf8(body)
            .map_err(|_| invalid(format!("Corrupt/truncated cache header: {path}")))?;
        require(
            !line
                .chars()
                .any(|c| ((c as u32) < 32 && c != '\t') || c as u32 == 127),
            "Control byte in cache header",
        )?;
        if number == 1 {
            let columns: Vec<&str> = line.split('\t').collect();
            let native = native_columns();
            if columns.len() == native.len() && columns.iter().zip(&native).all(|(a, b)| *a == b) {
                result = Some(Header {
                    format: Format::Spliceai2Tsv,
                    assembly: Some("GRCh38"),
                });
                break;
            }
            let version_ok = line
                .strip_prefix("##fileformat=VCFv4.")
                .is_some_and(|minor| {
                    minor.len() == 1 && matches!(minor.as_bytes()[0], b'0'..=b'5')
                });
            require(
                version_ok,
                "Unsupported cache header: expected VCF or exact 84-column SpliceAI2 TSV",
            )?;
            continue;
        }
        if line.starts_with("##") {
            let (key, value) = metadata(line)?;
            let lower = key.to_lowercase();
            require(lower != "fileformat", "Duplicate VCF fileformat line")?;
            if lower == "reference" || lower == "assembly" {
                assemblies.extend(assemblies_in(&value)?);
            } else if lower == "contig" {
                let chars: Vec<char> = value.chars().collect();
                let inner: &[char] = if chars.len() >= 2 {
                    &chars[1..chars.len() - 1]
                } else {
                    &[]
                };
                for capture in contig_assemblies(inner) {
                    assemblies.extend(assemblies_in(&capture)?);
                }
            }
            if lower.starts_with("spliceai2") || (lower == "source" && mentions_spliceai2(&value)) {
                versions.push(2);
            }
            if lower == "spliceai2version" {
                require(
                    is_dotted_version(&value, &['2']),
                    "Unsupported explicit SpliceAI2 version",
                )?;
            }
            if lower == "spliceaiversion" || lower == "spliceai_version" {
                require(
                    is_dotted_version(&value, &['1', '2']),
                    "Unsupported explicit SpliceAI version",
                )?;
                versions.push(value.as_bytes()[0] - b'0');
            }
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        require(
            columns.len() >= 8 && columns[..8] == VCF_COLUMNS,
            "Missing valid complete #CHROM VCF header",
        )?;
        let samples_ok = columns.len() == 8
            || (columns.len() >= 10 && columns[8] == "FORMAT" && {
                let samples = &columns[9..];
                let mut unique = samples.to_vec();
                unique.sort_unstable();
                unique.dedup();
                unique.len() == samples.len()
                    && samples.iter().all(|name| {
                        !name.is_empty()
                            && !name
                                .chars()
                                .any(|c| c.is_whitespace() || ('\u{1c}'..='\u{1f}').contains(&c))
                    })
            });
        require(samples_ok, "Invalid VCF sample columns")?;
        assemblies.sort_unstable();
        assemblies.dedup();
        versions.sort_unstable();
        versions.dedup();
        require(
            assemblies.len() <= 1,
            "Conflicting assembly metadata in VCF header",
        )?;
        require(versions.len() <= 1, "Conflicting SpliceAI version metadata")?;
        let version_two = versions.contains(&2);
        require(
            !version_two || assemblies != ["GRCh37"],
            "SpliceAI2 VCF conflicts with GRCh37 assembly",
        )?;
        result = Some(Header {
            format: if version_two {
                Format::VcfSpliceai2
            } else {
                Format::VcfV1
            },
            assembly: assemblies.first().copied().or(if version_two {
                Some("GRCh38")
            } else {
                None
            }),
        });
        break;
    }
    let header = result.ok_or_else(|| invalid("Cache header exceeds bounded line limit"))?;
    // Reading one more byte also finds a truncated gzip footer.
    let mut extra = [0_u8; 1];
    reader.read_exact(&mut extra).or_else(|error| {
        if error.kind() == std::io::ErrorKind::UnexpectedEof {
            Ok(())
        } else {
            Err(error)
        }
    })?;
    Ok(header)
}

// -------------------------------------------------------------- versions

fn score_version(snv: Option<&str>, indel: Option<&str>, assembly: Option<&str>) -> Result<Format> {
    let expected = match assembly {
        None | Some("") | Some("null") => None,
        Some(assembly) => Some(canonical_assembly(assembly)?),
    };
    let mut found: Vec<(&str, Format)> = Vec::new();
    for (kind, value) in [("SNV", snv), ("indel", indel)] {
        let Some(value) = value else { continue };
        if value.is_empty() || value == "null" || !Path::new(value).is_file() {
            continue;
        }
        let header = inspect_header(value)?;
        if let (Some(expected), Some(actual)) = (expected, header.assembly) {
            require(
                actual == expected,
                &format!(
                    "The SpliceAI {kind} score file is {actual}, but the configured assembly is {expected}: {value}"
                ),
            )?;
        }
        found.push((kind, header.format));
    }
    let mut versions: Vec<Format> = Vec::new();
    for (_, format) in &found {
        if !versions.contains(format) {
            versions.push(*format);
        }
    }
    if versions.len() > 1 {
        return Err(invalid(
            "The SpliceAI SNV and indel score files are different versions",
        ));
    }
    let version = versions.first().copied().unwrap_or(Format::VcfV1);
    if version == Format::Spliceai2Tsv {
        return Err(invalid(
            "SpliceAI2 publisher tables must be converted to indexed VCF by the installer first",
        ));
    }
    Ok(version)
}

/// `cutoffs_for_config(config_path)`.
pub(crate) fn cutoffs_for_config(config_path: Option<&str>) -> Result<Cutoffs> {
    let settings = read_config(config_path)?;
    let version = score_version(
        settings.snv.as_deref(),
        settings.indel.as_deref(),
        settings.assembly.as_deref(),
    )?;
    Ok(if version == Format::VcfSpliceai2 {
        Cutoffs {
            splice_altering: 0.5,
            nearby_acceptor_rescue: 0.25,
        }
    } else {
        Cutoffs {
            splice_altering: 0.8,
            nearby_acceptor_rescue: 0.5,
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn temporary(name: &str, bytes: &[u8]) -> String {
        let path =
            std::env::temp_dir().join(format!("priva_cutoffs_{}_{name}", std::process::id()));
        File::create(&path).unwrap().write_all(bytes).unwrap();
        path.to_string_lossy().into_owned()
    }

    #[test]
    fn config_values_are_read_from_top_level_lines() {
        let path = temporary(
            "config.yaml",
            b"# comment\nassembly: hg38   # note\ncache_folders:\n  clinvar: ClinVar\nspliceai_snv_prescore: \"/a b/snv.vcf.gz\"\nspliceai_indel_prescore: ''\nother: [1, 2]\n",
        );
        let values = read_config(Some(&path)).unwrap();
        assert_eq!(values.assembly.as_deref(), Some("hg38"));
        assert_eq!(values.snv.as_deref(), Some("/a b/snv.vcf.gz"));
        assert_eq!(values.indel.as_deref(), Some(""));
        for refused in [
            "spliceai_snv_prescore: 5\n",
            "spliceai_snv_prescore:\n  nested: 1\n",
            "- a\n",
            "---\nassembly: hg38\n",
            "not a mapping line\n",
        ] {
            let path = temporary("refused.yaml", refused.as_bytes());
            assert!(read_config(Some(&path)).is_err(), "{refused}");
        }
        assert!(read_config(None).unwrap().snv.is_none());
        assert!(
            read_config(Some("/definitely/not/here.yaml"))
                .unwrap()
                .snv
                .is_none()
        );
    }

    #[test]
    fn without_score_files_the_original_cutoffs_apply() {
        let cutoffs = cutoffs_for_config(None).unwrap();
        assert_eq!(cutoffs.splice_altering, 0.8);
        assert_eq!(cutoffs.nearby_acceptor_rescue, 0.5);
    }

    #[test]
    fn header_parts_follow_the_python_regular_expressions() {
        assert_eq!(
            assemblies_in("hg19 and GRCh38.p1").unwrap(),
            ["GRCh37", "GRCh38"]
        );
        assert!(assemblies_in("xhg19 hg191").unwrap().is_empty());
        let contig: Vec<char> = "ID=chr1,length=5,assembly=\"GRCh38\"".chars().collect();
        assert_eq!(contig_assemblies(&contig), ["\"GRCh38\""]);
        let plain: Vec<char> = "assembly=hg38,ID=chr1".chars().collect();
        assert_eq!(contig_assemblies(&plain), ["hg38"]);
        assert!(mentions_spliceai2("SpliceAI2 predictions"));
        assert!(mentions_spliceai2("spliceai v2"));
        assert!(mentions_spliceai2("SpliceAI-2_model"));
        assert!(!mentions_spliceai2("spliceai2x"));
        assert!(!mentions_spliceai2("SpliceAI 1.3"));
        assert!(is_dotted_version("2.1.0", &['2']));
        assert!(!is_dotted_version("2.", &['2']));
        assert!(!is_dotted_version("3", &['1', '2']));
        assert_eq!(
            metadata("##source=a=b").unwrap(),
            ("source".to_owned(), "a=b".to_owned())
        );
        assert!(metadata("##source=").is_err());
        assert!(metadata("##INFO=<ID=A,Description=\"x\">").is_ok());
        assert!(metadata("##INFO=<ID=A,Description=\"x>").is_err());
        assert_eq!(native_columns().len(), 84);
    }

    #[test]
    fn score_files_select_the_version() {
        let v1 = temporary(
            "v1.vcf",
            b"##fileformat=VCFv4.2\n##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        );
        let v2 = temporary(
            "v2.vcf",
            b"##fileformat=VCFv4.2\n##source=SpliceAI2 scores\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        );
        let config = |snv: &str, indel: &str, assembly: &str| {
            let text = format!(
                "spliceai_snv_prescore: \"{snv}\"\nspliceai_indel_prescore: \"{indel}\"\nassembly: {assembly}\n"
            );
            temporary("run.yaml", text.as_bytes())
        };
        let cutoffs = |path: &str| cutoffs_for_config(Some(path));
        assert_eq!(
            cutoffs(&config(&v1, &v1, "hg38")).unwrap().splice_altering,
            0.8
        );
        assert_eq!(
            cutoffs(&config(&v2, &v2, "hg38")).unwrap().splice_altering,
            0.5
        );
        assert_eq!(
            cutoffs(&config(&v2, &v2, "hg38"))
                .unwrap()
                .nearby_acceptor_rescue,
            0.25
        );
        assert!(
            cutoffs(&config(&v1, &v2, "hg38")).is_err(),
            "different versions"
        );
        assert!(
            cutoffs(&config(&v1, &v1, "hg19")).is_err(),
            "assembly differs"
        );
        assert!(
            cutoffs(&config(&v2, &v2, "hg19")).is_err(),
            "SpliceAI2 is GRCh38 only"
        );
        assert!(cutoffs(&config(&v2, "", "hg38")).is_ok());
    }
}
