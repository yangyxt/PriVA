//! Source-faithful, opt-in `stat_gene_patho_afs.py` producer for schema 6.
//! Chromosome results are merged in header order, never completion order.
//! The partial accumulated-AF snapshot and the original None merge failure
//! are deliberate compatibility behavior, not corrected scientific rules.
use crate::cache_convert::Provenance;
use crate::cache_convert::clinvar_gene_af::{self, ClinvarGeneAfCache, GeneAfEntry, GeneAfValue};
use crate::cli::{Result, invalid};
use crate::clingen_map_builder::{
    FormatType, check_samples, check_vcf_compression, format_definition,
};
use crate::rust_cache::ensure_distinct_paths;
use crate::tabix_index::{TabixIndex, fetch_contig_lines, read_variant_index};
use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

const HELP: &str = "usage: priva clinvar-gene-af build --vcf INDEXED_VCF --output RKYV [--threads N]\n\
       priva clinvar-gene-af chromosome --vcf INDEXED_VCF --chrom CONTIG --output RKYV\n\
       priva clinvar-gene-af merge --input CHR1_RKYV --input CHR2_RKYV ... --output RKYV\n\
Opt-in schema-6 producer. Merge inputs must follow the VCF header's contig order.\n\
Preserves None entries, first equal-AF winner and the original partial accumulated\n\
AF snapshot. A repeated gene with a None value fails as in the Python producer.\n\
An absent CSQ format returns success without creating output, as in the source.";
type Records = Vec<(String, Option<Vec<GeneAfEntry>>)>;

#[derive(Clone, Copy)]
enum InfoType {
    Integer,
    Float,
    Text,
    Flag,
}
#[derive(Clone, Copy)]
struct Definition {
    kind: InfoType,
    scalar: bool,
    genotype_count: bool,
}
struct Header {
    contigs: Vec<String>,
    csq: Option<HashMap<String, usize>>,
    info: HashMap<String, Definition>,
    formats: HashMap<String, FormatType>,
    samples: usize,
}

fn attributes(line: &str, prefix: &str) -> Option<HashMap<String, String>> {
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

fn header(path: &Path) -> Result<Header> {
    check_vcf_compression(path)?;
    let mut reader = crate::vcf_text::open_reader(path)?;
    let mut line = String::new();
    reader.read_line(&mut line)?;
    if !line.starts_with("##fileformat=VCFv") {
        return Err(invalid("input is not a VCF file"));
    }
    let mut h = Header {
        contigs: vec![],
        csq: None,
        info: HashMap::new(),
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
enum Scalar {
    Missing,
    Text(String),
    Int(i64),
    Float(f64),
    Bool(bool),
}
impl Scalar {
    fn text(&self) -> String {
        match self {
            Self::Missing => "None".into(),
            Self::Text(s) => s.clone(),
            Self::Int(i) => i.to_string(),
            Self::Float(f) => crate::python_text::float_repr(*f),
            Self::Bool(b) => if *b { "True" } else { "False" }.into(),
        }
    }
    fn number(&self) -> Result<f64> {
        match self {
            Self::Int(i) => Ok(*i as f64),
            Self::Float(f) => Ok(*f),
            Self::Bool(b) => Ok(u8::from(*b) as f64),
            _ => Err(invalid("TypeError: AF is not numeric")),
        }
    }
    fn cache(&self) -> Result<GeneAfValue> {
        match self {
            Self::Text(s) => Ok(GeneAfValue::Text(s.clone())),
            Self::Int(i) => Ok(GeneAfValue::Int(*i)),
            Self::Float(f) => Ok(GeneAfValue::Float(*f)),
            _ => Err(invalid("schema 6 cannot represent this AF scalar kind")),
        }
    }
    fn add(&self, other: &Self) -> Result<Self> {
        if let (Self::Int(a), Self::Int(b)) = (self, other) {
            return a
                .checked_add(*b)
                .map(Self::Int)
                .ok_or_else(|| invalid("AF integer sum exceeds schema 6"));
        }
        Ok(Self::Float(self.number()? + other.number()?))
    }
}
struct InfoValue {
    values: Vec<Scalar>,
    scalar: bool,
}
impl InfoValue {
    fn first(&self) -> Result<Scalar> {
        self.values
            .first()
            .cloned()
            .ok_or_else(|| invalid("IndexError: empty AF tuple"))
    }
    fn text(&self) -> String {
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

fn info_value(raw: Option<&str>, d: Definition) -> Result<InfoValue> {
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

fn record(
    line: &[u8],
    h: &Header,
    info_types: &mut HashMap<String, Definition>,
    records: &mut Records,
    positions: &mut HashMap<String, usize>,
    accumulated: &mut HashMap<String, Scalar>,
) -> Result<()> {
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
    if !info.contains_key("CLNSIG") || !info.contains_key("CLNREVSTAT") {
        return Ok(());
    }
    let get = |key: &str| -> Result<Option<InfoValue>> {
        let d = info_types
            .get(key)
            .copied()
            .ok_or_else(|| invalid("ValueError: unknown INFO header"))?;
        if !info.contains_key(key) && matches!(d.kind, InfoType::Flag) {
            return Ok(Some(InfoValue {
                values: vec![Scalar::Bool(false)],
                scalar: true,
            }));
        }
        info.get(key).map(|raw| info_value(*raw, d)).transpose()
    };
    let clnsig = get("CLNSIG")?.unwrap().text();
    if !clnsig.contains("athogenic") {
        return Ok(());
    }
    let clnrevstat = get("CLNREVSTAT")?.unwrap().text();
    let lower = clnrevstat.to_lowercase();
    if ![
        "practice_guideline",
        "reviewed_by_expert_panel",
        "criteria_provided,_multiple_submitters,_no_conflicts",
    ]
    .iter()
    .any(|s| lower.contains(s))
    {
        return Ok(());
    }
    let first =
        |key| -> Result<Scalar> { get(key)?.map(|v| v.first()).unwrap_or(Ok(Scalar::Int(0))) };
    let mut af = first("AF_grpmax_joint")?;
    if matches!(af, Scalar::Missing | Scalar::Bool(false))
        || matches!(af, Scalar::Int(0))
        || matches!(af, Scalar::Float(f) if f == 0.0)
    {
        af = first("AF_joint")?;
    }
    if matches!(af, Scalar::Missing) || af.number()? <= 0.0 {
        af = Scalar::Int(0);
    }
    let mut joint = first("AF_joint")?;
    if matches!(joint, Scalar::Missing) || joint.number()? <= 0.0 {
        joint = Scalar::Int(0);
    }
    if !info.contains_key("CSQ") {
        return Ok(());
    }
    let Some(csq) = get("CSQ")? else {
        return Ok(());
    };
    let annotations: Vec<Scalar> = if csq.scalar {
        match &csq.values[0] {
            Scalar::Text(s) => s.chars().map(|c| Scalar::Text(c.to_string())).collect(),
            _ => return Err(invalid("TypeError: CSQ scalar is not iterable")),
        }
    } else {
        csq.values
    };
    let csq_map = h.csq.as_ref().unwrap();
    for annotation in annotations {
        let Scalar::Text(annotation) = annotation else {
            return Err(invalid("AttributeError: CSQ entry has no split method"));
        };
        let annotation: Vec<_> = annotation.split('|').collect();
        let selected = ["Gene", "HGVSp", "IMPACT", "Consequence", "SYMBOL"]
            .map(|name| csq_map.get(name).and_then(|i| annotation.get(*i)).copied());
        let [
            Some(gene),
            Some(hgvsp),
            Some(impact),
            Some(consequence),
            Some(symbol),
        ] = selected
        else {
            continue;
        };
        if !gene.starts_with("ENSG") {
            continue;
        }
        let total = if let Some(previous) = accumulated.get(gene) {
            previous.add(&joint)?
        } else {
            joint.clone()
        };
        accumulated.insert(gene.to_owned(), total.clone());
        let id = *positions.entry(gene.to_owned()).or_insert_with(|| {
            let id = records.len();
            records.push((gene.to_owned(), None));
            id
        });
        let current = match &records[id].1 {
            Some(entries) => entry_af(entries)?,
            None => 0.0,
        };
        if af.number()? > current {
            let alt = fields
                .get(4)
                .filter(|s| **s != ".")
                .ok_or_else(|| invalid("TypeError: record.alts is None"))?;
            let alt = alt
                .split(',')
                .map(|s| if s.is_empty() { "." } else { s })
                .collect::<Vec<_>>()
                .join(",");
            let reference = fields
                .get(3)
                .copied()
                .filter(|s| !s.is_empty())
                .unwrap_or(".");
            let text = |s: &str| GeneAfValue::Text(s.to_owned());
            let values = [
                text(fields[0]),
                GeneAfValue::Int(position),
                text(reference),
                text(&alt),
                af.cache()?,
                text(hgvsp),
                text(impact),
                text(consequence),
                text(&clnsig),
                text(&clnrevstat),
                text(symbol),
                total.cache()?,
            ];
            let names = [
                "chrom",
                "pos",
                "ref",
                "alt",
                "af",
                "hgvsp",
                "impact",
                "consequence",
                "clnsig",
                "clnrevstat",
                "gene_symbol",
                "acc_patho_af",
            ];
            records[id].1 = Some(
                names
                    .into_iter()
                    .zip(values)
                    .map(|(key, value)| GeneAfEntry {
                        key: key.to_owned(),
                        value,
                    })
                    .collect(),
            );
        }
    }
    Ok(())
}

fn entry_af(entries: &[GeneAfEntry]) -> Result<f64> {
    match entries.iter().find(|e| e.key == "af").map(|e| &e.value) {
        Some(GeneAfValue::Int(i)) => Ok(*i as f64),
        Some(GeneAfValue::Float(f)) => Ok(*f),
        _ => Err(invalid("gene record lacks a numeric af")),
    }
}
fn merge(results: Vec<Records>) -> Result<Records> {
    let mut merged: Records = Vec::new();
    let mut positions: HashMap<String, usize> = HashMap::new();
    for result in results {
        for (gene, record) in result {
            if let Some(&id) = positions.get(&gene) {
                let old: &Option<Vec<GeneAfEntry>> = &merged[id].1;
                let (Some(old), Some(new)) = (old.as_ref(), record.as_ref()) else {
                    return Err(invalid(
                        "TypeError: None is not subscriptable during chromosome merge",
                    ));
                };
                if entry_af(old)? < entry_af(new)? {
                    merged[id].1 = record;
                }
            } else {
                positions.insert(gene.clone(), merged.len());
                merged.push((gene, record));
            }
        }
    }
    Ok(merged)
}
fn chromosome(vcf: &Path, index: &TabixIndex, h: &Header, chrom: &str) -> Result<Records> {
    if !h.contigs.iter().any(|name| name == chrom) {
        return Err(invalid("unknown header contig"));
    }
    let (mut records, mut positions, mut accumulated) = (vec![], HashMap::new(), HashMap::new());
    let mut info_types = h.info.clone();
    fetch_contig_lines(vcf, index, chrom, &mut |line| {
        record(
            line,
            h,
            &mut info_types,
            &mut records,
            &mut positions,
            &mut accumulated,
        )
        .map_err(|e| e.to_string())
    })
    .map_err(invalid)?;
    Ok(records)
}
fn appended(path: &Path, suffix: &str) -> PathBuf {
    let mut p = path.as_os_str().to_owned();
    p.push(suffix);
    p.into()
}
fn index(vcf: &Path) -> Result<TabixIndex> {
    let csi = appended(vcf, ".csi");
    read_variant_index(&if csi.exists() {
        csi
    } else {
        appended(vcf, ".tbi")
    })
    .map_err(invalid)
}
fn provenance(inputs: &[PathBuf], mode: &str) -> Result<Vec<Provenance>> {
    let mut result = vec![Provenance {
        key: "builder".into(),
        value: format!("priva clinvar-gene-af {mode} {}", env!("CARGO_PKG_VERSION")),
    }];
    for (i, path) in inputs.iter().enumerate() {
        result.push(Provenance {
            key: format!("source_{i}"),
            value: path
                .file_name()
                .and_then(|s| s.to_str())
                .ok_or_else(|| invalid("input needs a UTF-8 basename"))?
                .to_owned(),
        });
        result.push(Provenance {
            key: format!("source_{i}_bytes"),
            value: std::fs::metadata(path)?.len().to_string(),
        });
    }
    Ok(result)
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args.is_empty() || args.iter().any(|a| a == "--help" || a == "-h") {
        println!("{HELP}");
        return Ok(());
    }
    let mode = &args[0];
    if !["build", "chromosome", "merge"].contains(&mode.as_str()) {
        return Err(invalid(HELP));
    }
    let mut options = HashMap::new();
    let mut inputs = Vec::new();
    let mut at = 1;
    while at < args.len() {
        let name = &args[at];
        let value = args.get(at + 1).ok_or_else(|| invalid(HELP))?;
        if name == "--input" && mode == "merge" {
            inputs.push(PathBuf::from(value));
        } else if ["--vcf", "--output", "--threads", "--chrom"].contains(&name.as_str()) {
            options.insert(name.as_str(), value.as_str());
        } else {
            return Err(invalid(HELP));
        }
        at += 2;
    }
    let required = |key| options.get(key).copied().ok_or_else(|| invalid(HELP));
    let output = Path::new(required("--output")?);
    let records = if mode == "merge" {
        if inputs.is_empty() {
            return Err(invalid("merge requires ordered --input archives"));
        }
        let results = inputs
            .iter()
            .map(|p| {
                ensure_distinct_paths(p, output)?;
                Ok(ClinvarGeneAfCache::open(&p.to_string_lossy())?.owned_entries())
            })
            .collect::<Result<Vec<_>>>()?;
        merge(results)?
    } else {
        let vcf = Path::new(required("--vcf")?);
        let threads: i64 = options
            .get("--threads")
            .map(|s| s.parse())
            .transpose()?
            .unwrap_or_else(|| {
                std::thread::available_parallelism()
                    .map(|n| n.get().min(16) as i64)
                    .unwrap_or(1)
            });
        let mut h = header(vcf)?;
        if h.csq.is_none() {
            eprintln!("Could not extract CSQ format from VCF header");
            return Ok(());
        }
        if threads <= 0 {
            return Err(invalid("number of processes must be at least 1"));
        }
        ensure_distinct_paths(vcf, output)?;
        inputs.push(vcf.to_path_buf());
        let indexed = index(vcf);
        if let Ok(index) = &indexed {
            for name in &index.names {
                if !h.contigs.contains(name) {
                    h.contigs.push(name.clone());
                }
            }
        }
        if h.contigs.is_empty() {
            vec![]
        } else {
            let index = indexed?;
            if mode == "chromosome" {
                chromosome(vcf, &index, &h, required("--chrom")?)?
            } else {
                let next = AtomicUsize::new(0);
                let (send, receive) = std::sync::mpsc::channel();
                std::thread::scope(|scope| {
                    for _ in 0..(threads as usize).min(h.contigs.len()) {
                        let (send, h, index, next) = (send.clone(), &h, &index, &next);
                        scope.spawn(move || {
                            loop {
                                let i = next.fetch_add(1, Ordering::Relaxed);
                                if i >= h.contigs.len() {
                                    break;
                                }
                                let result = chromosome(vcf, index, h, &h.contigs[i])
                                    .map_err(|e| e.to_string());
                                if send.send((i, result)).is_err() {
                                    break;
                                }
                            }
                        });
                    }
                });
                drop(send);
                let mut results = receive.into_iter().collect::<Vec<_>>();
                results.sort_by_key(|(i, _)| *i);
                let results = results
                    .into_iter()
                    .map(|(_, r)| r.map_err(invalid))
                    .collect::<Result<Vec<_>>>()?;
                merge(results)?
            }
        }
    };
    let count = clinvar_gene_af::write_entries(output, records, provenance(&inputs, mode)?)?;
    println!("{{\"genes\":{count},\"schema\":6,\"version\":1}}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn merge_preserves_first_tie_and_none_failure() {
        let value = |af| {
            Some(vec![GeneAfEntry {
                key: "af".into(),
                value: GeneAfValue::Float(af),
            }])
        };
        assert!(
            merge(vec![
                vec![("g".into(), None)],
                vec![("g".into(), value(1.0))]
            ])
            .is_err()
        );
        assert!(
            merge(vec![
                vec![("g".into(), value(1.0))],
                vec![("g".into(), None)]
            ])
            .is_err()
        );
        assert!(merge(vec![vec![("g".into(), None)], vec![("g".into(), None)]]).is_err());
        assert_eq!(
            merge(vec![vec![("a".into(), None)], vec![("b".into(), None)]])
                .unwrap()
                .len(),
            2
        );
    }
    #[test]
    fn info_float_is_widened_binary32_and_keeps_nonfinite() {
        let d = Definition {
            kind: InfoType::Float,
            scalar: false,
            genotype_count: false,
        };
        let v = info_value(Some("0.1,.,NaN,inf"), d).unwrap();
        assert_eq!(
            v.values[0].number().unwrap().to_bits(),
            (0.1_f32 as f64).to_bits()
        );
        assert!(matches!(v.values[1], Scalar::Missing));
        assert!(v.values[2].number().unwrap().is_nan());
        assert_eq!(v.values[3].number().unwrap(), f64::INFINITY);
    }
}
