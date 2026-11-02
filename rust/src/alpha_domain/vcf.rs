//! Source-specific AlphaMissense VCF selection; typed decoding stays shared.
use super::key::Key;
use super::model::{Model, Observation};
use super::{Result, exception, source_error, unsupported};
use crate::vcf_variant::{self, Definition, Header, InfoValue, Record, Scalar};
use std::collections::HashMap;
use std::io::{BufRead, Read};
use std::path::Path;

pub(super) struct Context {
    pub header: Header,
    pub fields: HashMap<&'static str, usize>,
}
impl Context {
    pub fn new(path: &Path, metadata: bool) -> Result<Self> {
        constructor_prefix(path)?;
        let header = vcf_variant::header(path).map_err(|error| {
            if error.to_string() == "failed to fill whole buffer" {
                exception(
                    "ValueError",
                    "invalid file: VCF format cannot be determined",
                )
            } else {
                source_error(error)
            }
        })?;
        let names =
            crate::clinvar_aa_splice_builder::whole_header_names(path).map_err(source_error)?;
        let required: &[&'static str] = if metadata {
            &[
                "DOMAINS",
                "Gene",
                "SYMBOL",
                "Protein_position",
                "Feature_type",
                "Feature",
            ]
        } else {
            &[
                "DOMAINS",
                "Gene",
                "EXON",
                "Protein_position",
                "Feature_type",
                "Feature",
            ]
        };
        let mut fields = HashMap::new();
        for &name in required {
            let at = names
                .iter()
                .position(|s| s == name)
                .ok_or_else(|| exception("ValueError", format!("'{name}' is not in list")))?;
            fields.insert(name, at);
        }
        Ok(Self { header, fields })
    }
    pub fn field<'a>(&self, fields: &[&'a str], name: &str) -> Result<&'a str> {
        fields
            .get(self.fields[name])
            .copied()
            .ok_or_else(|| exception("IndexError", "list index out of range"))
    }
}

fn constructor_prefix(path: &Path) -> Result<()> {
    // Pinned pysam/HTSlib performs format detection before its BGZF EOF and
    // header checks: at most 2048 compressed bytes produce at most 1024 peek
    // bytes. A short gzip prefix with no recognizable VCF prefix is ValueError,
    // while the lone nontext gzip-magic byte fails hts_hopen with OSError.
    // Keep this boundary local; the shared VCF reader's error policy is unchanged.
    let mut prefix = Vec::new();
    std::fs::File::open(path)?
        .take(2048)
        .read_to_end(&mut prefix)?;
    if prefix == [0x1f] {
        return Err(exception(
            "OSError",
            "cannot open unrecognized binary variant input",
        ));
    }
    if prefix.starts_with(&[0x1f, 0x8b]) {
        let mut sniff = Vec::new();
        let _ = flate2::read::MultiGzDecoder::new(prefix.as_slice())
            .take(1024)
            .read_to_end(&mut sniff);
        if !sniff.starts_with(b"##fileformat=VCF") {
            return Err(exception(
                "ValueError",
                "invalid file: VCF format cannot be determined",
            ));
        }
    }
    Ok(())
}

fn get(
    record: &Record<'_>,
    key: &str,
    definitions: &HashMap<String, Definition>,
) -> Result<InfoValue> {
    if !definitions.contains_key(key) {
        return Err(exception(
            "KeyError",
            "Invalid INFO field: ".to_owned() + key,
        ));
    }
    record
        .get(key, definitions)
        .map_err(source_error)?
        .ok_or_else(|| exception("KeyError", "Invalid INFO field: ".to_owned() + key))
}
fn float(value: &Scalar) -> Result<f64> {
    match value {
        Scalar::Text(text) => crate::python_text::parse_float(text)
            .ok_or_else(|| exception("ValueError", "could not convert string to float")),
        Scalar::Missing => Err(exception(
            "TypeError",
            "float() argument must be a string or a real number, not 'NoneType'",
        )),
        _ => value.number().map_err(source_error),
    }
}
fn score(value: InfoValue) -> Result<f64> {
    if !value.scalar {
        return Err(exception(
            "TypeError",
            "float() argument must be a string or a real number, not 'tuple'",
        ));
    }
    float(&value.values.into_iter().next().unwrap_or(Scalar::Missing))
}
fn af(value: InfoValue) -> Result<f64> {
    if value.scalar {
        match value.values.into_iter().next().unwrap_or(Scalar::Missing) {
            Scalar::Text(text) => {
                let first = text
                    .chars()
                    .next()
                    .ok_or_else(|| exception("IndexError", "string index out of range"))?;
                float(&Scalar::Text(first.to_string()))
            }
            _ => Err(exception("TypeError", "AF scalar is not subscriptable")),
        }
    } else {
        float(&value.first().map_err(source_error)?)
    }
}
pub(super) fn csqs(value: InfoValue) -> Result<Vec<String>> {
    if value.scalar {
        match value.values.into_iter().next().unwrap_or(Scalar::Missing) {
            Scalar::Text(text) => Ok(text.split(',').map(str::to_owned).collect()),
            _ => Err(exception(
                "AttributeError",
                "CSQ scalar has no attribute 'split'",
            )),
        }
    } else {
        value
            .values
            .into_iter()
            .map(|value| match value {
                Scalar::Text(text) => Ok(text),
                _ => Err(exception(
                    "AttributeError",
                    "CSQ tuple item has no attribute 'split'",
                )),
            })
            .collect()
    }
}
pub(super) fn position(text: &str) -> Result<Option<i64>> {
    let first = text.split('/').next().unwrap().split('-').next().unwrap();
    integer(first)
}
pub(super) fn integer(text: &str) -> Result<Option<i64>> {
    if !crate::python_text::python_int_is_valid(text) {
        return Ok(None);
    }
    let converted: String = text
        .chars()
        .map(|c| {
            if let Some(digit) = crate::python_text::decimal_digit_value(c) {
                char::from(b'0' + digit as u8)
            } else if !c.is_ascii() && c.is_whitespace() {
                ' '
            } else {
                c
            }
        })
        .collect();
    match crate::python_text::python_int(&converted) {
        Some(value) => i64::try_from(value)
            .map(Some)
            .map_err(|_| unsupported("protein position exceeds signed 64-bit Parquet column")),
        None => Err(unsupported(
            "protein position exceeds signed 64-bit Parquet column",
        )),
    }
}

pub(super) fn record(
    model: &mut Model,
    record: &Record<'_>,
    definitions: &HashMap<String, Definition>,
    context: &Context,
) -> Result<()> {
    let mut run = || -> Result<()> {
        let score = score(get(record, "AM_PATHOGENICITY", definitions)?)?;
        let af = match get(record, "AF_grpmax_joint", definitions).and_then(af) {
            Ok(value) => value,
            Err(error)
                if ["KeyError:", "ValueError:", "IndexError:"]
                    .iter()
                    .any(|kind| error.to_string().starts_with(kind)) =>
            {
                0.0
            }
            Err(error) => return Err(error),
        };
        let variant = Key::from_info(get(record, "PVAR", definitions)?)?;
        let residue = variant.residue();
        // Process one CSQ at a time: later failures must not undo earlier
        // transcript/domain mutations from this same VCF record.
        let value = get(record, "CSQ", definitions)?;
        let annotations = if value.scalar {
            match value.values.into_iter().next().unwrap_or(Scalar::Missing) {
                Scalar::Text(text) => text
                    .split(',')
                    .map(|s| Scalar::Text(s.to_owned()))
                    .collect(),
                _ => {
                    return Err(exception(
                        "AttributeError",
                        "CSQ scalar has no attribute 'split'",
                    ));
                }
            }
        } else {
            value.values
        };
        for annotation in annotations {
            let Scalar::Text(annotation) = annotation else {
                return Err(exception(
                    "AttributeError",
                    "CSQ tuple item has no attribute 'split'",
                ));
            };
            let fields: Vec<_> = annotation.split('|').collect();
            if context.field(&fields, "Feature_type")? != "Transcript"
                || !context.field(&fields, "Feature")?.starts_with("ENST")
            {
                continue;
            }
            let domains = context.field(&fields, "DOMAINS")?;
            let gene = context.field(&fields, "Gene")?;
            let transcript = context.field(&fields, "Feature")?;
            let exon = context.field(&fields, "EXON")?;
            let position = position(context.field(&fields, "Protein_position")?)?;
            if domains.is_empty() || gene.is_empty() {
                continue;
            }
            model.observe(Observation {
                gene,
                domains,
                transcript,
                exon,
                position,
                variant: &variant,
                residue: &residue,
                score,
                af,
            })?;
        }
        Ok(())
    };
    match run() {
        Err(error)
            if error.to_string().starts_with("KeyError:")
                || error.to_string().starts_with("ValueError:") =>
        {
            println!(
                "Error processing record {}:{} - {}",
                record.fields[0], record.position, error
            );
            Ok(())
        }
        other => other,
    }
}

pub(super) fn visit(
    path: &Path,
    context: &Context,
    chrom: Option<&str>,
    mut visit: impl FnMut(&Record<'_>, &HashMap<String, Definition>) -> Result<()>,
) -> Result<()> {
    let mut definitions = context.header.info.clone();
    let mut consume = |line: &[u8]| -> Result<()> {
        let record = vcf_variant::decode_record(line, &context.header, &mut definitions).map_err(
            |error| {
                let message = error.to_string();
                if message.starts_with("VCF record is not UTF-8") {
                    return exception("UnicodeDecodeError", message);
                }
                if message == "invalid VCF position" {
                    let negative = line
                        .split(|b| *b == b'\t')
                        .nth(1)
                        .is_some_and(|s| s.starts_with(b"-"));
                    return exception(
                        if negative && chrom.is_some() {
                            "ValueError"
                        } else {
                            "OSError"
                        },
                        message,
                    );
                }
                if message.starts_with("ValueError:") || message.starts_with("TypeError:") {
                    return source_error(error);
                }
                exception(
                    if chrom.is_some() {
                        "ValueError"
                    } else {
                        "OSError"
                    },
                    message,
                )
            },
        )?;
        visit(&record, &definitions)
    };
    if let Some(chrom) = chrom {
        if !context.header.contigs.iter().any(|s| s == chrom) {
            return Err(exception("ValueError", "invalid contig"));
        }
        let index = vcf_variant::index(path).map_err(|e| exception("ValueError", e.to_string()))?;
        crate::tabix_index::fetch_contig_lines(path, &index, chrom, &mut |line| {
            consume(line).map_err(|e| e.to_string())
        })
        .map_err(|e| source_error(crate::cli::invalid(e)))?;
    } else {
        let mut reader = crate::vcf_text::open_reader(path).map_err(source_error)?;
        let mut line = Vec::new();
        loop {
            line.clear();
            if reader.read_until(b'\n', &mut line)? == 0 {
                break;
            }
            if line.starts_with(b"#") {
                continue;
            }
            consume(&line)?;
        }
    }
    Ok(())
}
