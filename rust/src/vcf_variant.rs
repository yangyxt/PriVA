//! Typed indexed-VCF primitives shared by source-mapped cache builders.
//! Source-specific CSQ interpretation, filtering, reductions and failures remain
//! with each builder. INFO floats preserve htslib's binary32 promotion; records
//! validate unused sample fields and retain the first repeated INFO key.
use crate::cli::{Result, invalid};
use crate::clingen_map_builder::{
    FormatType, check_samples, check_vcf_compression, format_definition,
};
use crate::tabix_index::{TabixIndex, read_variant_index};
use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::path::{Path, PathBuf};

#[derive(Clone, Copy)]
pub(crate) enum InfoType {
    Integer,
    Float,
    Text,
    Flag,
}
#[derive(Clone, Copy)]
pub(crate) struct Definition {
    pub(crate) kind: InfoType,
    pub(crate) scalar: bool,
    pub(crate) genotype_count: bool,
}
pub(crate) struct Header {
    pub(crate) contigs: Vec<String>,
    pub(crate) csq: Option<HashMap<String, usize>>,
    pub(crate) info: HashMap<String, Definition>,
    #[allow(dead_code)] // Used by builders with their own CSQ format rule.
    pub(crate) info_descriptions: HashMap<String, String>,
    pub(crate) formats: HashMap<String, FormatType>,
    pub(crate) samples: usize,
}

pub(crate) fn attributes(line: &str, prefix: &str) -> Option<HashMap<String, String>> {
    let body = line.strip_prefix(prefix)?.strip_suffix('>')?;
    let (mut start, mut quoted, mut escaped) = (0, false, false);
    let mut parts = Vec::new();
    for (at, c) in body.char_indices() {
        if escaped {
            escaped = false;
            continue;
        }
        if c == '\\' && quoted {
            escaped = true;
        } else if c == '"' {
            quoted = !quoted;
        } else if c == ',' && !quoted {
            parts.push(&body[start..at]);
            start = at + 1;
        }
    }
    parts.push(&body[start..]);
    Some(
        parts
            .into_iter()
            .filter_map(|part| part.split_once('='))
            .map(|(k, v)| (k.to_owned(), v.trim_matches('"').to_owned()))
            .collect(),
    )
}

pub(crate) fn header(path: &Path) -> Result<Header> {
    check_vcf_compression(path)?;
    let mut reader = crate::vcf_text::open_reader(path)?;
    read_header(reader.as_mut())
}

/// Parse a header from an already selected VCF text stream. Container handling
/// remains with the caller (the annotation contract also accepts BCF headers).
pub(crate) fn read_header(reader: &mut dyn BufRead) -> Result<Header> {
    let mut line = String::new();
    reader.read_line(&mut line)?;
    if !line.starts_with("##fileformat=VCFv") {
        return Err(invalid("input is not a VCF file"));
    }
    let mut h = Header {
        contigs: vec![],
        csq: None,
        info: HashMap::new(),
        info_descriptions: HashMap::new(),
        formats: HashMap::new(),
        samples: 0,
    };
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            return Err(invalid("VCF column header is absent"));
        }
        let line = line.trim_end_matches(['\r', '\n']);
        if !line.starts_with("##") {
            let fields: Vec<_> = line.split('\t').collect();
            if fields.len() < 8
                || fields[..8]
                    != [
                        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                    ]
            {
                return Err(invalid("invalid VCF column header"));
            }
            if fields.len() > 8 {
                if fields[8] != "FORMAT" {
                    return Err(invalid("invalid VCF FORMAT header"));
                }
                h.samples = fields.len() - 9;
                let names: HashSet<_> = fields[9..].iter().copied().collect();
                if names.len() != h.samples || names.contains("") {
                    return Err(invalid("invalid VCF sample names"));
                }
            }
            break;
        }
        if let Some(attrs) = attributes(line, "##contig=<") {
            if let Some(id) = attrs.get("ID") {
                if !h.contigs.contains(id) {
                    h.contigs.push(id.clone());
                }
            }
        }
        let mut retained = true;
        if let Some(attrs) = attributes(line, "##INFO=<") {
            if let Some(id) = attrs.get("ID") {
                retained = !h.info.contains_key(id);
                if let Some(description) = attrs.get("Description") {
                    h.info_descriptions
                        .entry(id.clone())
                        .or_insert_with(|| description.clone());
                }
                let kind = match attrs.get("Type").map(String::as_str) {
                    Some("Integer") => InfoType::Integer,
                    Some("Float") => InfoType::Float,
                    Some("Flag") => InfoType::Flag,
                    _ => InfoType::Text,
                };
                h.info.entry(id.clone()).or_insert(Definition {
                    kind,
                    scalar: attrs.get("Number").is_some_and(|s| s == "1"),
                    genotype_count: attrs.get("Number").is_some_and(|s| s == "G"),
                });
            }
        }
        if let Some((id, kind)) = format_definition(line.as_bytes()) {
            h.formats.entry(id).or_insert(kind);
        }
        if retained && h.csq.is_none() && line.contains("ID=CSQ") {
            if let Some((_, rest)) = line.split_once("Format: ") {
                let format: String = rest
                    .chars()
                    .take_while(|c| c.is_alphanumeric() || matches!(c, '_' | '|' | '-'))
                    .collect();
                if !format.is_empty() {
                    h.csq = Some(
                        format
                            .split('|')
                            .enumerate()
                            .map(|(i, s)| (s.to_owned(), i))
                            .collect(),
                    );
                }
            }
        }
    }
    Ok(h)
}

#[derive(Clone)]
pub(crate) enum Scalar {
    Missing,
    Text(String),
    Int(i64),
    Float(f64),
    Bool(bool),
}
impl Scalar {
    pub(crate) fn text(&self) -> String {
        match self {
            Self::Missing => "None".into(),
            Self::Text(s) => s.clone(),
            Self::Int(i) => i.to_string(),
            Self::Float(f) => crate::python_text::float_repr(*f),
            Self::Bool(b) => if *b { "True" } else { "False" }.into(),
        }
    }
    pub(crate) fn number(&self) -> Result<f64> {
        match self {
            Self::Int(i) => Ok(*i as f64),
            Self::Float(f) => Ok(*f),
            Self::Bool(b) => Ok(u8::from(*b) as f64),
            _ => Err(invalid("TypeError: AF is not numeric")),
        }
    }
}
pub(crate) struct InfoValue {
    pub(crate) values: Vec<Scalar>,
    pub(crate) scalar: bool,
}
impl InfoValue {
    pub(crate) fn first(&self) -> Result<Scalar> {
        self.values
            .first()
            .cloned()
            .ok_or_else(|| invalid("IndexError: empty AF tuple"))
    }
    pub(crate) fn text(&self) -> String {
        self.values
            .iter()
            .map(Scalar::text)
            .collect::<Vec<_>>()
            .join(",")
    }
}

// htslib INFO integers use a decimal prefix, no leading whitespace. A sign
// alone consumes a character and yields zero; out-of-range int32 is missing.
fn info_integer(s: &str) -> Scalar {
    let sign = usize::from(s.starts_with(['+', '-']));
    let count = s.as_bytes()[sign..]
        .iter()
        .take_while(|b| b.is_ascii_digit())
        .count();
    if sign + count == 0 {
        return Scalar::Missing;
    }
    let value = if count == 0 {
        Some(0)
    } else {
        s[..sign + count].parse::<i64>().ok()
    };
    value
        .filter(|v| (-2147483640..=2147483647).contains(v))
        .map(Scalar::Int)
        .unwrap_or(Scalar::Missing)
}

fn float_prefix(s: &str) -> Option<f64> {
    let s = s.trim_start_matches(|c: char| c.is_ascii_whitespace());
    let sign = usize::from(s.starts_with(['+', '-']));
    let rest = &s[sign..];
    let lower = rest.to_ascii_lowercase();
    if lower.starts_with("nan") {
        return Some(if s.starts_with('-') {
            -f64::NAN
        } else {
            f64::NAN
        });
    }
    if lower.starts_with("inf") {
        return Some(if s.starts_with('-') {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        });
    }
    if lower.starts_with("0x") {
        if let Some(value) = hex_float(&rest[2..]) {
            return Some(if s.starts_with('-') { -value } else { value });
        }
    }
    let bytes = s.as_bytes();
    let mut at = sign;
    let mut digits = 0;
    while at < bytes.len() && bytes[at].is_ascii_digit() {
        at += 1;
        digits += 1;
    }
    if bytes.get(at) == Some(&b'.') {
        at += 1;
        while at < bytes.len() && bytes[at].is_ascii_digit() {
            at += 1;
            digits += 1;
        }
    }
    if digits == 0 {
        return None;
    }
    if matches!(bytes.get(at), Some(b'e' | b'E')) {
        let exp = at;
        at += 1;
        if matches!(bytes.get(at), Some(b'+' | b'-')) {
            at += 1;
        }
        let start = at;
        while at < bytes.len() && bytes[at].is_ascii_digit() {
            at += 1;
        }
        if at == start {
            at = exp;
        }
    }
    s[..at].parse::<f64>().ok()
}

// Exact binary rounding of strtod's hexadecimal syntax, without a C parser.
fn hex_float(s: &str) -> Option<f64> {
    let (mut nibbles, mut point, mut end) = (Vec::new(), None, 0);
    for (at, c) in s.char_indices() {
        if let Some(value) = c.to_digit(16) {
            nibbles.push(value as u8);
            end = at + 1;
        } else if c == '.' && point.is_none() {
            point = Some(nibbles.len());
            end = at + 1;
        } else {
            break;
        }
    }
    if nibbles.is_empty() {
        return None;
    }
    let point = point.unwrap_or(nibbles.len());
    let tail = &s[end..];
    let mut exponent = 0_i64;
    if tail.starts_with(['p', 'P']) {
        let token = &tail[1..];
        let sign = usize::from(token.starts_with(['+', '-']));
        let count = token.as_bytes()[sign..]
            .iter()
            .take_while(|b| b.is_ascii_digit())
            .count();
        if count > 0 {
            exponent = token[..sign + count]
                .parse::<i64>()
                .unwrap_or(if token.starts_with('-') {
                    -1000000
                } else {
                    1000000
                });
        }
    }
    let bits: Vec<bool> = nibbles
        .iter()
        .flat_map(|n| (0..4).rev().map(move |b| n & (1 << b) != 0))
        .collect();
    let Some(first) = bits.iter().position(|b| *b) else {
        return Some(0.0);
    };
    let e = exponent.saturating_add(4 * point as i64 - first as i64 - 1);
    if e > 1023 {
        return Some(f64::INFINITY);
    }
    if e < -1075 {
        return Some(0.0);
    }
    let kept = if e >= -1022 { 53 } else { (e + 1075) as usize };
    let mut mantissa = 0_u64;
    for i in 0..kept {
        mantissa = (mantissa << 1) | u64::from(bits.get(first + i).copied().unwrap_or(false));
    }
    let guard = bits.get(first + kept).copied().unwrap_or(false);
    let sticky = bits.iter().skip(first + kept + 1).any(|b| *b);
    if guard && (sticky || mantissa & 1 != 0) {
        mantissa += 1;
    }
    if e < -1022 {
        return Some(f64::from_bits(mantissa));
    }
    let mut e = e;
    if mantissa == 1 << 53 {
        mantissa >>= 1;
        e += 1;
    }
    if e > 1023 {
        return Some(f64::INFINITY);
    }
    Some(f64::from_bits(
        ((e + 1023) as u64) << 52 | (mantissa & ((1 << 52) - 1)),
    ))
}

pub(crate) fn info_value(raw: Option<&str>, d: Definition) -> Result<InfoValue> {
    if d.genotype_count {
        return Err(invalid(
            "ValueError: INFO Number=G cannot determine genotype count",
        ));
    }
    let mut values = Vec::new();
    if matches!(d.kind, InfoType::Flag) {
        return Ok(InfoValue {
            values: vec![Scalar::Bool(true)],
            scalar: true,
        });
    }
    if raw.is_none() || (raw == Some("") && matches!(d.kind, InfoType::Text)) {
        if matches!(d.kind, InfoType::Flag) {
            values.push(Scalar::Bool(true));
        } else if d.scalar {
            values.push(Scalar::Missing);
        }
    } else {
        let raw = raw.unwrap();
        for token in raw.split(',') {
            values.push(match d.kind {
                InfoType::Integer => info_integer(token),
                InfoType::Float => float_prefix(token)
                    .map(|v| Scalar::Float((v as f32) as f64))
                    .unwrap_or(Scalar::Missing),
                InfoType::Text | InfoType::Flag => {
                    if raw.len() > 1 && token.is_empty() {
                        Scalar::Missing
                    } else {
                        Scalar::Text(token.to_owned())
                    }
                }
            });
        }
    }
    let scalar = (d.scalar && values.len() == 1) || matches!(values.as_slice(), [Scalar::Bool(_)]);
    Ok(InfoValue { values, scalar })
}

/// Decoded VCF columns and first-occurrence INFO tokens. Policy stays in callers.
pub(crate) struct Record<'a> {
    pub(crate) fields: Vec<&'a str>,
    pub(crate) position: i64,
    pub(crate) info: HashMap<&'a str, Option<&'a str>>,
}
impl Record<'_> {
    pub(crate) fn get(
        &self,
        key: &str,
        info_types: &HashMap<String, Definition>,
    ) -> Result<Option<InfoValue>> {
        let d = info_types
            .get(key)
            .copied()
            .ok_or_else(|| invalid("ValueError: unknown INFO header"))?;
        if !self.info.contains_key(key) && matches!(d.kind, InfoType::Flag) {
            return Ok(Some(InfoValue {
                values: vec![Scalar::Bool(false)],
                scalar: true,
            }));
        }
        self.info
            .get(key)
            .map(|raw| info_value(*raw, d))
            .transpose()
    }
}
pub(crate) fn decode_record<'a>(
    line: &'a [u8],
    h: &Header,
    info_types: &mut HashMap<String, Definition>,
) -> Result<Record<'a>> {
    let line = std::str::from_utf8(line)
        .map_err(|_| invalid("VCF record is not UTF-8"))?
        .trim_end_matches(['\r', '\n']);
    let fields: Vec<_> = line.split('\t').collect();
    if h.samples > 0 {
        check_samples(
            &mut fields.iter().skip(5).map(|s| s.as_bytes()),
            h.samples,
            &h.formats,
        )?;
    }
    let pos = fields.get(1).copied().unwrap_or("");
    let position: i64 = if pos.is_empty() {
        0
    } else {
        pos.parse().map_err(|_| invalid("invalid VCF position"))?
    };
    if position < 0 {
        return Err(invalid("invalid VCF position"));
    }
    let mut info = HashMap::new();
    for item in fields.get(7).copied().unwrap_or(".").split(';') {
        let (key, value) = item
            .split_once('=')
            .map(|(k, v)| (k, Some(v)))
            .unwrap_or((item, None));
        if key.is_empty() || key == "." {
            continue;
        }
        info_types.entry(key.to_owned()).or_insert(Definition {
            kind: InfoType::Text,
            scalar: true,
            genotype_count: false,
        });
        info.entry(key).or_insert(value);
    }
    Ok(Record {
        fields,
        position,
        info,
    })
}

fn appended(path: &Path, suffix: &str) -> PathBuf {
    let mut p = path.as_os_str().to_owned();
    p.push(suffix);
    p.into()
}
pub(crate) fn index(vcf: &Path) -> Result<TabixIndex> {
    let csi = appended(vcf, ".csi");
    read_variant_index(&if csi.exists() {
        csi
    } else {
        appended(vcf, ".tbi")
    })
    .map_err(invalid)
}
