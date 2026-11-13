//! Header-only validation. Records are deliberately not decoded. BCF carries
//! the same VCF header text; only its container envelope is consumed here.
use super::{error, object, strings};
use crate::cli::Result;
use crate::json_text::Json;
use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufRead, Cursor, Read, Seek, SeekFrom};
use std::path::Path;

const AA: &[&str] = &[
    "Feature",
    "Feature_type",
    "Protein_position",
    "HGVSp",
    "Consequence",
    "EXON",
    "INTRON",
    "HGVSc",
    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL",
];
const GENE: &[&str] = &[
    "SYMBOL",
    "Gene",
    "BIOTYPE",
    "Feature",
    "DOMAINS",
    "Consequence",
    "LoF",
    "NMD",
    "SpliceVault_top_events",
    "SpliceVault_SpliceAI_delta",
    "SpliceAI_pred_DS_DL",
];
const AF: &[&str] = &[
    "Gene",
    "Feature",
    "EXON",
    "HGVSp",
    "IMPACT",
    "Consequence",
    "SYMBOL",
];
const POPS: &[&str] = &[
    "nfe",
    "eas",
    "afr",
    "amr",
    "asj",
    "fin",
    "sas",
    "mid",
    "remaining",
];

pub(crate) fn header_text(path: &Path) -> Result<Vec<u8>> {
    crate::clingen_map_builder::check_vcf_compression(path).map_err(|e| {
        if e.to_string().starts_with("NotImplementedError:") {
            error(
                "NotImplementedError",
                "seek not implemented in files compressed by method 1",
            )
        } else {
            e
        }
    })?;
    let mut file = File::open(path)?;
    let mut magic = [0; 2];
    if file.read(&mut magic)? == 2 && magic == [31, 139] {
        let mut eof = [0; 28];
        if file.metadata()?.len() < 28 {
            return Err(error(
                "OSError",
                "no BGZF EOF marker; file may be truncated",
            ));
        }
        file.seek(SeekFrom::End(-28))?;
        file.read_exact(&mut eof)?;
        if eof != crate::vcf_text::EOF_BLOCK {
            return Err(error(
                "OSError",
                "no BGZF EOF marker; file may be truncated",
            ));
        }
    }
    let mut stream = crate::vcf_text::open_reader(path)?;
    if stream.fill_buf()?.starts_with(b"BCF") {
        let mut envelope = [0; 9];
        stream
            .read_exact(&mut envelope)
            .map_err(|_| error("ValueError", "invalid BCF header"))?;
        if envelope[..5] != *b"BCF\x02\x02" {
            return Err(error("ValueError", "invalid BCF version"));
        }
        let length = u32::from_le_bytes(envelope[5..].try_into().unwrap()) as u64;
        let mut text = Vec::new();
        stream.take(length).read_to_end(&mut text)?;
        if text.len() as u64 != length {
            return Err(error("ValueError", "incomplete BCF header"));
        }
        if text.last() == Some(&0) {
            text.pop();
        }
        return Ok(text);
    }
    let mut text = Vec::new();
    loop {
        let start = text.len();
        if stream.read_until(b'\n', &mut text)? == 0 {
            break;
        }
        if !text[start..].starts_with(b"##") {
            break;
        }
    }
    Ok(text)
}

fn list_repr(values: &[String]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|s| crate::python_path::repr(s))
            .collect::<Vec<_>>()
            .join(", ")
    )
}

fn csq_description(text: &str) -> Option<String> {
    // htslib retains the first named attribute and the first INFO definition.
    // Leading attribute whitespace is ignored; escaped quotes/backslashes are
    // retained in Description. These rules are local to header validation.
    for line in text.lines() {
        let Some(body) = line
            .trim_end_matches('\r')
            .strip_prefix("##INFO=<")
            .and_then(|s| s.strip_suffix('>'))
        else {
            continue;
        };
        let (mut start, mut quoted, mut escaped) = (0, false, false);
        let mut fields = Vec::new();
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
                fields.push(&body[start..at]);
                start = at + 1;
            }
        }
        fields.push(&body[start..]);
        let attributes: Vec<_> = fields
            .iter()
            .filter_map(|f| f.trim_start().split_once('='))
            .collect();
        let get = |name| {
            attributes
                .iter()
                .find(|(key, _)| *key == name)
                .map(|(_, value)| value.trim_matches('"'))
        };
        if get("ID") == Some("CSQ") {
            return get("Description").map(str::to_owned);
        }
    }
    None
}
pub(crate) fn validate(path: &str, stage: &str) -> Result<Json> {
    if !["aa", "gene", "af", "vep", "complete"].contains(&stage) {
        return Err(error(
            "ValueError",
            format!("Unknown ClinVar validation stage: {stage}"),
        ));
    }
    let bytes = header_text(Path::new(path))?;
    let text = crate::python_io::decode_utf8(&bytes).map_err(|e| error("UnicodeDecodeError", e))?;
    let header = crate::vcf_variant::read_header(&mut Cursor::new(&bytes))?;
    // pysam refuses FORMAT without any sample, even though the shared record
    // parser can represent that header. Keep this constructor rule local.
    if text
        .lines()
        .last()
        .is_some_and(|s| s.trim_end_matches('\r').ends_with("\tFORMAT"))
    {
        return Err(error("ValueError", "invalid VCF column header"));
    }
    let info: BTreeSet<_> = header.info.keys().cloned().collect();
    let mut required_info: BTreeSet<String> = ["CSQ", "CLNSIG", "CLNREVSTAT"]
        .into_iter()
        .map(str::to_owned)
        .collect();
    let mut required_csq = BTreeSet::new();
    if ["aa", "vep", "complete"].contains(&stage) {
        required_csq.extend(AA.iter().copied());
    }
    if ["gene", "vep", "complete"].contains(&stage) {
        required_csq.extend(GENE.iter().copied());
        required_info.insert("ALLELEID".into());
    }
    if ["af", "complete"].contains(&stage) {
        required_csq.extend(AF.iter().copied());
        required_info.extend(["AF_joint".into(), "AF_grpmax_joint".into()]);
    }
    if ["gene", "complete"].contains(&stage) {
        required_info.extend(
            [
                "AC_joint",
                "AN_joint",
                "AF_joint",
                "nhomalt_joint",
                "AC_grpmax_joint",
                "AF_grpmax_joint",
                "AN_grpmax_joint",
                "nhomalt_grpmax_joint",
            ]
            .into_iter()
            .map(str::to_owned),
        );
        for pop in POPS {
            for sex in ["XX", "XY"] {
                required_info.insert(format!("AF_joint_{pop}_{sex}"));
            }
            required_info.insert(format!("nhomalt_joint_{pop}_XX"));
        }
    }
    if !info.contains("CSQ") {
        return Err(error("ValueError", "ClinVar annotation lacks INFO/CSQ"));
    }
    let description = csq_description(text)
        .ok_or_else(|| error("TypeError", "argument of type 'NoneType' is not iterable"))?;
    let (_, format) = description.split_once("Format: ").ok_or_else(|| {
        error(
            "ValueError",
            "ClinVar INFO/CSQ header lacks its Format declaration",
        )
    })?;
    let csq: Vec<String> = format.split('|').map(str::to_owned).collect();
    let unique: BTreeSet<_> = csq.iter().map(String::as_str).collect();
    if unique.len() != csq.len() {
        return Err(error(
            "ValueError",
            "ClinVar INFO/CSQ header contains duplicate fields",
        ));
    }
    let missing_info: Vec<_> = required_info.difference(&info).cloned().collect();
    let missing_csq: Vec<_> = required_csq
        .difference(&unique)
        .map(|s| (*s).to_owned())
        .collect();
    if !missing_info.is_empty() || !missing_csq.is_empty() {
        return Err(error(
            "ValueError",
            format!(
                "Incomplete ClinVar {stage} annotation: missing INFO={}; missing CSQ={}",
                list_repr(&missing_info),
                list_repr(&missing_csq)
            ),
        ));
    }
    Ok(object([
        ("path", Json::Str(path.into())),
        ("stage", Json::Str(stage.into())),
        ("success", Json::Bool(true)),
        ("csq_fields", strings(csq)),
        ("info_fields", strings(info)),
        ("contigs", strings(header.contigs)),
    ]))
}
