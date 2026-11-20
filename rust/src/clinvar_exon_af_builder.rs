//! Opt-in `stat_exon_patho_afs.py` producer. Keep source clinical filtering,
//! per-record last-exon selection, exon insertion order and ordered AF samples.
//! Chromosome archives carry their raw samples; merge recomputes global stats.
use crate::builder_error;
use crate::cache_convert::Provenance;
use crate::cache_convert::clinvar_exon_af::{self, ClinvarExonAfCache, ExonAf};
use crate::cli::{Result, invalid};
use crate::rust_cache::ensure_distinct_paths;
use crate::tabix_index::{TabixIndex, fetch_contig_lines};
use crate::vcf_variant::{
    Definition, Header, InfoValue, Record, Scalar, decode_record, header as vcf_header, index,
};
use std::collections::HashMap;
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

const HELP: &str = "usage: priva clinvar-exon-af build --vcf INDEXED_VCF --output RKYV [--threads N]\n\
       priva clinvar-exon-af chromosome --vcf INDEXED_VCF --chrom CONTIG --output RKYV\n\
       priva clinvar-exon-af merge --input CHR1_RKYV --input CHR2_RKYV ... --output RKYV\n\
Opt-in schema-4 producer. Supply merge inputs in the VCF header's contig order.\n\
Exon order and every raw AF sample are retained; merge recomputes the four\n\
global statistics. No installer or scientific caller is switched by this command.";
type Entries = Vec<(String, Vec<ExonAf>)>;

fn exception(kind: &'static str, message: impl Into<String>) -> Box<dyn std::error::Error> {
    builder_error::category(kind, invalid(message.into()))
}

// Shared VCF primitives also serve older builders. Interpret their explicit
// exception labels at this source's boundary without changing shared policy.
fn source_error(error: Box<dyn std::error::Error>) -> Box<dyn std::error::Error> {
    let message = error.to_string();
    // The VCF text adapter returns a string for a short BGZF read; pysam's
    // VariantFile constructor exposes that framing failure as OSError.
    if message == "failed to fill whole buffer" {
        return exception("OSError", message);
    }
    for kind in [
        "TypeError",
        "IndexError",
        "ValueError",
        "KeyError",
        "AttributeError",
        "NotImplementedError",
        "OSError",
        "UnicodeDecodeError",
    ] {
        if let Some(rest) = message
            .strip_prefix(kind)
            .and_then(|s| s.strip_prefix(": "))
        {
            return exception(kind, rest);
        }
    }
    error
}

fn header(path: &Path) -> Result<Header> {
    crate::clingen_map_builder::check_vcf_compression(path).map_err(source_error)?;
    let mut file = std::fs::File::open(path)?;
    let mut magic = [0; 2];
    if file.read(&mut magic)? == 2 && magic == [0x1f, 0x8b] {
        // VariantFile's default check_truncation requires the literal EOF
        // block, even when every preceding compressed block is complete.
        let mut eof = [0; 28];
        if file.metadata()?.len() < eof.len() as u64 {
            return Err(exception(
                "OSError",
                "no BGZF EOF marker; file may be truncated",
            ));
        }
        file.seek(SeekFrom::End(-(eof.len() as i64)))?;
        file.read_exact(&mut eof)?;
        if eof != crate::vcf_text::EOF_BLOCK {
            return Err(exception(
                "OSError",
                "no BGZF EOF marker; file may be truncated",
            ));
        }
    }
    vcf_header(path).map_err(source_error)
}

fn csq_fields(h: &Header) -> Result<(usize, usize)> {
    let mut fields = HashMap::new();
    if let Some(description) = h.info_descriptions.get("CSQ") {
        if description.contains("Format:") {
            let format = description
                .split("Format: ")
                .nth(1)
                .ok_or_else(|| exception("IndexError", "list index out of range"))?;
            fields.extend(format.split('|').enumerate().map(|(i, name)| (name, i)));
        }
    }
    if !fields.contains_key("Transcript") || !fields.contains_key("EXON") {
        if let Some(&feature) = fields.get("Feature") {
            fields.insert("Transcript", feature);
        }
        if !fields.contains_key("EXON") {
            if let Some(&exon) = fields.get("exon") {
                fields.insert("EXON", exon);
            }
        }
    }
    match (fields.get("Transcript"), fields.get("EXON")) {
        (Some(&transcript), Some(&exon)) => Ok((transcript, exon)),
        _ => Err(exception(
            "ValueError",
            "Required CSQ fields not found in VCF header",
        )),
    }
}

fn iterable(value: InfoValue) -> Result<Vec<Scalar>> {
    if !value.scalar {
        return Ok(value.values);
    }
    match value.values.into_iter().next().unwrap_or(Scalar::Missing) {
        Scalar::Text(text) => Ok(text.chars().map(|c| Scalar::Text(c.to_string())).collect()),
        _ => Err(exception("TypeError", "INFO scalar is not iterable")),
    }
}

fn join_text(value: InfoValue) -> Result<String> {
    let values = iterable(value)?;
    let text = values
        .iter()
        .map(|v| match v {
            Scalar::Text(text) => Ok(text.as_str()),
            _ => Err(exception(
                "TypeError",
                "clinical INFO join requires string entries",
            )),
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(text.join(","))
}

fn af_value(value: InfoValue) -> Result<Scalar> {
    let first = value.first().map_err(source_error)?;
    if value.scalar {
        return Ok(first);
    }
    // Python casts ONLY tuple entries. Scalar strings stay strings and fail
    // at the later <= 0 comparison; invalid tuple strings skip the record.
    match first {
        Scalar::Text(text) => crate::python_text::parse_float(&text)
            .map(Scalar::Float)
            .ok_or_else(|| exception("ValueError", "AF tuple entry cannot be converted to float")),
        Scalar::Missing => Err(exception(
            "TypeError",
            "float() argument must be a string or a real number, not 'NoneType'",
        )),
        value => Ok(Scalar::Float(value.number().map_err(source_error)?)),
    }
}

#[derive(Default)]
struct Samples {
    entries: Entries,
    transcripts: HashMap<String, usize>,
    exons: Vec<HashMap<String, usize>>,
}
impl Samples {
    fn extend(&mut self, transcript: &str, exon: &str, values: impl IntoIterator<Item = f64>) {
        let t = *self
            .transcripts
            .entry(transcript.to_owned())
            .or_insert_with(|| {
                let t = self.entries.len();
                self.entries.push((transcript.to_owned(), vec![]));
                self.exons.push(HashMap::new());
                t
            });
        let e = *self.exons[t].entry(exon.to_owned()).or_insert_with(|| {
            let entries = &mut self.entries[t].1;
            let e = entries.len();
            entries.push(ExonAf {
                exon_key: exon.to_owned(),
                item0: 0.0,
                item1: 0.0,
                item2: 0.0,
                item3: 0.0,
                item4: vec![],
            });
            e
        });
        self.entries[t].1[e].item4.extend(values);
    }
    fn merge(&mut self, entries: Entries) {
        for (transcript, exons) in entries {
            for exon in exons {
                self.extend(&transcript, &exon.exon_key, exon.item4);
            }
        }
    }
    fn finish(mut self) -> Entries {
        for (_, exons) in &mut self.entries {
            exons.retain(|e| !e.item4.is_empty());
            for exon in exons {
                [exon.item0, exon.item1, exon.item2, exon.item3] =
                    crate::exon_af_stats::calculate(&exon.item4);
            }
        }
        self.entries
    }
}

fn selected_record(
    record: &Record<'_>,
    definitions: &HashMap<String, Definition>,
    fields: (usize, usize),
    samples: &mut Samples,
) -> Result<()> {
    if ["CLNSIG", "CLNREVSTAT", "CSQ"]
        .iter()
        .any(|key| !record.info.contains_key(key))
    {
        return Ok(());
    }
    let get = |key| record.get(key, definitions).map_err(source_error);
    let significance = join_text(get("CLNSIG")?.unwrap())?;
    if !["Pathogenic", "Likely_pathogenic"]
        .iter()
        .any(|s| significance.contains(s))
    {
        return Ok(());
    }
    let review = join_text(get("CLNREVSTAT")?.unwrap())?;
    if ![
        "practice_guideline",
        "reviewed_by_expert_panel",
        "criteria_provided,_multiple_submitters,_no_conflicts",
    ]
    .iter()
    .any(|s| review.contains(s))
    {
        return Ok(());
    }
    let mut af = Scalar::Missing;
    if record.info.contains_key("AF_grpmax_joint") {
        af = af_value(get("AF_grpmax_joint")?.unwrap())?;
    }
    let fallback = matches!(af, Scalar::Missing | Scalar::Int(0) | Scalar::Bool(false))
        || matches!(af, Scalar::Float(v) if v == 0.0);
    if fallback && record.info.contains_key("AF_joint") {
        af = af_value(get("AF_joint")?.unwrap())?;
    }
    let af = if matches!(af, Scalar::Missing) {
        0.0
    } else {
        let number = af.number().map_err(source_error)?;
        if number <= 0.0 { 0.0 } else { number }
    };
    let mut selected = Vec::new();
    let mut positions = HashMap::new();
    for annotation in iterable(get("CSQ")?.unwrap())? {
        let Scalar::Text(annotation) = annotation else {
            return Err(exception("AttributeError", "CSQ entry has no split method"));
        };
        let annotation = annotation.split('|').collect::<Vec<_>>();
        let (Some(&transcript), Some(&exon)) = (annotation.get(fields.0), annotation.get(fields.1))
        else {
            continue;
        };
        if transcript.is_empty() || exon.is_empty() {
            continue;
        }
        let at = *positions.entry(transcript.to_owned()).or_insert_with(|| {
            let at = selected.len();
            selected.push((transcript.to_owned(), String::new()));
            at
        });
        selected[at].1 = exon.to_owned();
    }
    for (transcript, exon) in selected {
        samples.extend(&transcript, &exon, [af]);
    }
    Ok(())
}

fn chromosome(
    vcf: &Path,
    index: &TabixIndex,
    h: &Header,
    fields: (usize, usize),
    chrom: &str,
) -> Result<Entries> {
    if !h.contigs.iter().any(|c| c == chrom) {
        return Err(exception("ValueError", "invalid contig"));
    }
    let mut samples = Samples::default();
    let mut definitions = h.info.clone();
    fetch_contig_lines(vcf, index, chrom, &mut |line| {
        // Decoding is outside the source's KeyError/ValueError skip boundary.
        let record =
            decode_record(line, h, &mut definitions).map_err(|e| format!("ValueError: {e}"))?;
        match selected_record(&record, &definitions, fields, &mut samples) {
            Ok(()) => Ok(()),
            Err(error) => {
                let message = error.to_string();
                if message.starts_with("KeyError: ") || message.starts_with("ValueError: ") {
                    Ok(())
                } else {
                    Err(message)
                }
            }
        }
    })
    .map_err(|message| {
        if message == "a VCF line has a position that is not an integer" {
            exception("OSError", message)
        } else {
            source_error(invalid(message))
        }
    })?;
    Ok(samples.entries)
}

struct Options {
    mode: String,
    output: PathBuf,
    vcf: Option<PathBuf>,
    chrom: Option<String>,
    inputs: Vec<PathBuf>,
    threads: i64,
}
fn execute(options: Options) -> Result<()> {
    let mut samples = Samples::default();
    let inputs;
    if options.mode == "merge" {
        for path in &options.inputs {
            ensure_distinct_paths(path, &options.output)?;
            samples.merge(ClinvarExonAfCache::open(&path.to_string_lossy())?.owned_entries());
        }
        inputs = options.inputs;
    } else {
        let vcf = options.vcf.as_ref().unwrap();
        let mut h = header(vcf).map_err(source_error)?;
        let fields = csq_fields(&h)?;
        if options.threads <= 0 {
            return Err(exception(
                "ValueError",
                "Number of processes must be at least 1",
            ));
        }
        ensure_distinct_paths(vcf, &options.output)?;
        let indexed = index(vcf);
        if let Ok(index) = &indexed {
            for name in &index.names {
                if !h.contigs.contains(name) {
                    h.contigs.push(name.clone());
                }
            }
        }
        if !h.contigs.is_empty() {
            let index = indexed
                .map_err(|e| exception("ValueError", format!("fetch requires an index: {e}")))?;
            if options.mode == "chromosome" {
                samples.merge(chromosome(
                    vcf,
                    &index,
                    &h,
                    fields,
                    options.chrom.as_ref().unwrap(),
                )?);
            } else {
                let next = AtomicUsize::new(0);
                let (send, receive) = std::sync::mpsc::channel();
                std::thread::scope(|scope| {
                    for _ in 0..(options.threads as usize).min(h.contigs.len()) {
                        let (send, h, index, next) = (send.clone(), &h, &index, &next);
                        scope.spawn(move || {
                            loop {
                                let at = next.fetch_add(1, Ordering::Relaxed);
                                if at >= h.contigs.len() {
                                    break;
                                }
                                let result = chromosome(vcf, index, h, fields, &h.contigs[at])
                                    .map_err(|e| e.to_string());
                                if send.send((at, result)).is_err() {
                                    break;
                                }
                            }
                        });
                    }
                });
                drop(send);
                let mut results = receive.into_iter().collect::<Vec<_>>();
                results.sort_by_key(|(at, _)| *at);
                for (_, result) in results {
                    samples.merge(result.map_err(|s| source_error(invalid(s)))?);
                }
            }
        }
        inputs = vec![vcf.clone()];
    }
    let entries = samples.finish();
    let exons: usize = entries.iter().map(|(_, exons)| exons.len()).sum();
    let count: usize = entries
        .iter()
        .flat_map(|(_, exons)| exons)
        .map(|e| e.item4.len())
        .sum();
    let mut provenance = vec![Provenance {
        key: "builder".into(),
        value: format!(
            "priva clinvar-exon-af {} {}",
            options.mode,
            env!("CARGO_PKG_VERSION")
        ),
    }];
    for (i, path) in inputs.iter().enumerate() {
        provenance.push(Provenance {
            key: format!("source_{i}"),
            value: path
                .file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .into_owned(),
        });
        provenance.push(Provenance {
            key: format!("source_{i}_bytes"),
            value: std::fs::metadata(path)?.len().to_string(),
        });
    }
    let transcripts = clinvar_exon_af::write_entries(&options.output, entries, provenance)?;
    println!(
        "{{\"transcripts\":{transcripts},\"exons\":{exons},\"samples\":{count},\"schema\":4,\"version\":1}}"
    );
    Ok(())
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args.iter().any(|s| s == "--help" || s == "-h") {
        println!("{HELP}");
        return Ok(());
    }
    let mode = args
        .first()
        .filter(|s| ["build", "chromosome", "merge"].contains(&s.as_str()))
        .ok_or_else(|| invalid(HELP))?;
    let mut options = HashMap::new();
    let mut inputs = Vec::new();
    let mut at = 1;
    while at < args.len() {
        let name = &args[at];
        let value = args.get(at + 1).ok_or_else(|| invalid(HELP))?;
        let admitted = match mode.as_str() {
            "build" => ["--vcf", "--output", "--threads"].contains(&name.as_str()),
            "chromosome" => ["--vcf", "--output", "--chrom"].contains(&name.as_str()),
            _ => ["--input", "--output"].contains(&name.as_str()),
        };
        if !admitted {
            return Err(invalid(HELP));
        }
        if name == "--input" {
            inputs.push(PathBuf::from(value));
        } else {
            options.insert(name.as_str(), value.as_str());
        }
        at += 2;
    }
    let required = |name| options.get(name).copied().ok_or_else(|| invalid(HELP));
    let output = required("--output")?.into();
    let vcf = if mode != "merge" {
        Some(required("--vcf")?.into())
    } else {
        None
    };
    let chrom = if mode == "chromosome" {
        Some(required("--chrom")?.to_owned())
    } else {
        None
    };
    if mode == "merge" && inputs.is_empty() {
        return Err(invalid("merge requires ordered --input archives"));
    }
    let threads = options
        .get("--threads")
        .map(|s| s.parse::<i64>())
        .transpose()?
        .unwrap_or(4);
    builder_error::finish(execute(Options {
        mode: mode.clone(),
        output,
        vcf,
        chrom,
        inputs,
        threads,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header_with_format(format: &str) -> Header {
        Header {
            contigs: vec![],
            csq: None,
            info: HashMap::new(),
            info_descriptions: HashMap::from([("CSQ".into(), format.into())]),
            formats: HashMap::new(),
            samples: 0,
        }
    }

    #[test]
    fn header_aliases_and_missing_fields_follow_source_order() {
        assert_eq!(
            csq_fields(&header_with_format("Format: Transcript|Feature|exon")).unwrap(),
            (1, 2)
        );
        assert_eq!(
            csq_fields(&header_with_format("Format: Transcript|EXON|Transcript")).unwrap(),
            (2, 1)
        );
        assert!(
            csq_fields(&header_with_format("Format:Transcript|EXON"))
                .unwrap_err()
                .to_string()
                .starts_with("IndexError:")
        );
        assert!(
            csq_fields(&header_with_format("Format: Transcript|EXON tail"))
                .unwrap_err()
                .to_string()
                .starts_with("ValueError:")
        );
    }

    #[test]
    fn merge_reduces_ordered_raw_samples_not_chromosome_statistics() {
        let mut first = Samples::default();
        first.extend("T", "2/4", [0.1, 0.3]);
        first.extend("T", "1/7", [0.7]);
        let mut second = Samples::default();
        second.extend("T", "2/4", [0.9]);
        let mut merged = Samples::default();
        merged.merge(first.finish());
        merged.merge(second.finish());
        let merged = merged.finish();
        assert_eq!(
            merged[0]
                .1
                .iter()
                .map(|e| e.exon_key.as_str())
                .collect::<Vec<_>>(),
            ["2/4", "1/7"]
        );
        assert_eq!(merged[0].1[0].item4, [0.1, 0.3, 0.9]);
        assert_eq!(merged[0].1[0].item0, 0.3);
        assert_eq!(
            merged[0].1[0].item1.to_bits(),
            ((0.1_f64 + 0.3 + 0.9) / 3.0).to_bits()
        );
    }
}
