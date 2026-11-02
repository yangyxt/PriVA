//! Source-faithful, opt-in `stat_gene_patho_afs.py` producer for schema 6.
//! Chromosome results are merged in header order, never completion order.
//! The partial accumulated-AF snapshot and the original None merge failure
//! are deliberate compatibility behavior, not corrected scientific rules.
use crate::cache_convert::Provenance;
use crate::cache_convert::clinvar_gene_af::{
    self, ClinvarGeneAfCache, Entry as GeneAfEntry, Value as GeneAfValue,
};
use crate::cli::{Result, invalid};
use crate::rust_cache::ensure_distinct_paths;
use crate::tabix_index::{TabixIndex, fetch_contig_lines};
use crate::vcf_variant::{Definition, Header, Scalar, decode_record, header, index};
#[cfg(test)]
use crate::vcf_variant::{InfoType, info_value};
use std::collections::HashMap;
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

impl Scalar {
    fn cache(&self) -> Result<GeneAfValue> {
        match self {
            Self::Text(s) => Ok(GeneAfValue::Text(s.clone())),
            Self::Int(i) => Ok(GeneAfValue::Int(*i)),
            Self::Float(f) => Ok(GeneAfValue::Float(*f)),
            Self::Bool(b) => Ok(GeneAfValue::Bool(*b)),
            _ => Err(invalid("schema 6 cannot represent this AF scalar kind")),
        }
    }
    fn add(&self, other: &Self) -> Result<Self> {
        // Python bool is an integer in addition, but the unadded first value
        // remains Boolean in the saved partial accumulated-AF snapshot.
        let integer = |value: &Self| match value {
            Self::Int(i) => Some(*i),
            Self::Bool(b) => Some(i64::from(*b)),
            _ => None,
        };
        if let (Some(a), Some(b)) = (integer(self), integer(other)) {
            return a
                .checked_add(b)
                .map(Self::Int)
                .ok_or_else(|| invalid("AF integer sum exceeds schema 6"));
        }
        Ok(Self::Float(self.number()? + other.number()?))
    }
}

fn record(
    line: &[u8],
    h: &Header,
    info_types: &mut HashMap<String, Definition>,
    records: &mut Records,
    positions: &mut HashMap<String, usize>,
    accumulated: &mut HashMap<String, Scalar>,
) -> Result<()> {
    let decoded = decode_record(line, h, info_types)?;
    let fields = &decoded.fields;
    let position = decoded.position;
    let info = &decoded.info;
    if !info.contains_key("CLNSIG") || !info.contains_key("CLNREVSTAT") {
        return Ok(());
    }
    let get = |key: &str| decoded.get(key, info_types);
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
        Some(GeneAfValue::Bool(b)) => Ok(u8::from(*b) as f64),
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
    let count = clinvar_gene_af::write_records(output, records, provenance(&inputs, mode)?)?;
    let version = ClinvarGeneAfCache::open(&output.to_string_lossy())?.schema_version();
    println!("{{\"genes\":{count},\"schema\":6,\"version\":{version}}}");
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
