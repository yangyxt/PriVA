//! Opt-in port of stat_aachange_clinvar.py. AA and splice mutations share one
//! record traversal; completed chromosome results merge in header order.
//! The source's whole-header CSQ rule, partial SpliceAI dictionaries, duplicate
//! observations and narrow per-record KeyError/ValueError catch are deliberate.
use crate::cache_convert::clinvar_splice_stat::typed::{
    Cache as SpliceCache, Record as SpliceRecord, write_records,
};
use crate::cache_convert::{ClinvarAaCache, PlainTranscript, Provenance, write_clinvar_aa};
use crate::cli::{Result, invalid};
use crate::rust_cache::ensure_distinct_paths;
use crate::tabix_index::{TabixIndex, fetch_contig_lines};
use crate::vcf_variant::{Definition, Header, InfoValue, Scalar, decode_record, header, index};
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

const HELP: &str = "usage: priva clinvar-aa-splice build --vcf INDEXED_VCF --aa-output AA_RKYV [--splice-output SPLICE_RKYV] [--threads N]\n\
       priva clinvar-aa-splice chromosome --vcf INDEXED_VCF --chrom CONTIG --aa-output AA_RKYV --splice-output SPLICE_RKYV\n\
       priva clinvar-aa-splice merge --aa-input AA_CHR1 --splice-input SPLICE_CHR1 ... --aa-output AA_RKYV --splice-output SPLICE_RKYV\n\
Merge pairs must follow the VCF header's contig order. All processing completes\n\
before AA then splice publication. Empty output names suppress that output.\n\
Schema 1 remains v1; schema 2 uses v2 only for source-supported nullable fields\n\
or SpliceAI dictionaries that cannot be represented by frozen v1.";
const CONFIDENCE: [&str; 3] = [
    "practice_guideline",
    "reviewed_by_expert_panel",
    "criteria_provided,_multiple_submitters,_no_conflicts",
];

#[derive(Default)]
struct Pair {
    aa: Vec<(String, PlainTranscript)>,
    splice: Vec<(String, Vec<SpliceRecord>)>,
}
#[derive(Default)]
struct Collector {
    values: Pair,
    aa_keys: HashMap<String, usize>,
    positions: HashMap<(String, String), (usize, usize)>,
    changes: HashMap<(String, String, String), (usize, usize, usize)>,
    splice_keys: HashMap<String, usize>,
}
impl Collector {
    fn aa(&mut self, key: &str, pos: &str, hgvsp: &str, sig: &[String], review: &[String]) {
        let (i, j, k) = *self
            .changes
            .entry((key.into(), pos.into(), hgvsp.into()))
            .or_insert_with(|| {
                let (i, j) = *self
                    .positions
                    .entry((key.into(), pos.into()))
                    .or_insert_with(|| {
                        let i = *self.aa_keys.entry(key.into()).or_insert_with(|| {
                            let i = self.values.aa.len();
                            self.values.aa.push((key.into(), vec![]));
                            i
                        });
                        let j = self.values.aa[i].1.len();
                        self.values.aa[i].1.push((pos.into(), vec![]));
                        (i, j)
                    });
                let k = self.values.aa[i].1[j].1.len();
                self.values.aa[i].1[j]
                    .1
                    .push((hgvsp.into(), vec![], vec![]));
                (i, j, k)
            });
        let change = &mut self.values.aa[i].1[j].1[k];
        change.1.extend_from_slice(sig);
        change.2.extend_from_slice(review);
    }
    fn splice(&mut self, key: &str, records: Vec<SpliceRecord>) {
        let i = *self.splice_keys.entry(key.into()).or_insert_with(|| {
            let i = self.values.splice.len();
            self.values.splice.push((key.into(), vec![]));
            i
        });
        self.values.splice[i].1.extend(records);
    }
    fn merge(&mut self, pair: Pair) {
        for (key, positions) in pair.aa {
            for (pos, changes) in positions {
                for (hgvsp, sig, review) in changes {
                    self.aa(&key, &pos, &hgvsp, &sig, &review);
                }
            }
        }
        for (key, records) in pair.splice {
            self.splice(&key, records);
        }
    }
}

struct Indices {
    fields: [Option<usize>; 8],
    splice: Vec<(String, usize)>,
}
fn indices(vcf: &Path) -> Result<Indices> {
    let owned_names = whole_header_names(vcf)?;
    let names: Vec<_> = owned_names.iter().map(String::as_str).collect();
    indices_from_names(&names)
}

/// The source AA/splice and AlphaMissense builders both split the complete
/// rendered header at its last `Format: `, including its trailing header text.
pub(crate) fn whole_header_names(vcf: &Path) -> Result<Vec<String>> {
    // pysam's constructor checks the exact BGZF EOF marker before any record.
    // Plain text is allowed to reach the original indexed-fetch failure later.
    let mut file = std::fs::File::open(vcf)?;
    let mut magic = [0; 2];
    if file.read_exact(&mut magic).is_ok() && magic == [31, 139] {
        file.seek(SeekFrom::End(-28))?;
        let mut eof = [0; 28];
        file.read_exact(&mut eof)?;
        if eof != crate::vcf_text::EOF_BLOCK {
            return Err(invalid(
                "OSError: no BGZF EOF marker; file may be truncated",
            ));
        }
    }
    let mut reader = crate::vcf_text::open_reader(vcf)?;
    let mut rendered = String::new();
    let mut seen = HashSet::new();
    loop {
        let mut line = String::new();
        if reader.read_line(&mut line)? == 0 {
            return Err(invalid("VCF column header is absent"));
        }
        let line = line.trim_end_matches(['\r', '\n']);
        // htslib retains the first structured definition for each ID.
        let duplicate = [
            "##INFO=<ID=",
            "##FORMAT=<ID=",
            "##FILTER=<ID=",
            "##contig=<ID=",
        ]
        .iter()
        .find_map(|prefix| {
            line.strip_prefix(prefix)
                .map(|s| (*prefix, s.split([',', '>']).next().unwrap()))
        })
        .is_some_and(|(kind, id)| !seen.insert((kind.to_owned(), id.to_owned())));
        if !duplicate {
            rendered.push_str(line);
            rendered.push('\n');
        }
        if line.starts_with("#CHROM\t") {
            break;
        }
    }
    let format = rendered
        .rsplit("Format: ")
        .next()
        .unwrap()
        .trim_matches(['>', '"']);
    Ok(format.split('|').map(str::to_owned).collect())
}

fn indices_from_names(names: &[&str]) -> Result<Indices> {
    let fields = [
        "Protein_position",
        "HGVSp",
        "Consequence",
        "Feature_type",
        "Feature",
        "EXON",
        "INTRON",
        "HGVSc",
    ]
    .map(|name| names.iter().position(|s| *s == name));
    let mut splice: Vec<(String, usize)> = Vec::new();
    for (i, name) in names
        .iter()
        .enumerate()
        .filter(|(_, s)| s.starts_with("SpliceAI_pred_"))
    {
        if let Some((_, index)) = splice.iter_mut().find(|(key, _)| key == name) {
            *index = i;
        } else {
            splice.push(((*name).into(), i));
        }
    }
    Ok(Indices { fields, splice })
}

fn joined(value: InfoValue) -> Result<String> {
    if value.scalar {
        match value.values.first() {
            Some(Scalar::Text(s)) => Ok(s
                .chars()
                .map(|c| c.to_string())
                .collect::<Vec<_>>()
                .join(",")),
            _ => Err(invalid("TypeError: can only join an iterable")),
        }
    } else {
        value
            .values
            .into_iter()
            .map(|value| match value {
                Scalar::Text(s) => Ok(s),
                _ => Err(invalid("TypeError: sequence item is not str")),
            })
            .collect::<Result<Vec<_>>>()
            .map(|v| v.join(","))
    }
}

fn record(
    line: &[u8],
    h: &Header,
    types: &mut HashMap<String, Definition>,
    indices: &Indices,
    out: &mut Collector,
) -> Result<()> {
    // Iterator decoding happens outside Python's record-level exception handler.
    let decoded = decode_record(line, h, types)?;
    let mut work = || -> Result<()> {
        if !decoded.info.contains_key("CLNSIG") || !decoded.info.contains_key("CLNREVSTAT") {
            return Ok(());
        }
        let review = joined(decoded.get("CLNREVSTAT", types)?.unwrap())?;
        let sig = joined(decoded.get("CLNSIG", types)?.unwrap())?;
        if !["Pathogenic", "Likely_pathogenic"]
            .iter()
            .any(|s| sig.contains(s))
            || !CONFIDENCE.iter().any(|s| review.contains(s))
        {
            return Ok(());
        }
        let Some(csq) = decoded.get("CSQ", types)? else {
            return Ok(());
        };
        let annotations = if csq.scalar {
            match csq.values.first() {
                Some(Scalar::Text(s)) => s.split(',').map(|s| Scalar::Text(s.into())).collect(),
                _ => return Err(invalid("AttributeError: CSQ scalar has no split method")),
            }
        } else {
            csq.values
        };
        for annotation in annotations {
            let Scalar::Text(annotation) = annotation else {
                return Err(invalid("AttributeError: CSQ entry has no split method"));
            };
            let fields: Vec<_> = annotation.split('|').collect();
            if indices
                .fields
                .iter()
                .flatten()
                .max()
                .is_some_and(|i| fields.len() <= *i)
            {
                continue;
            }
            let key = fields[indices.fields[4].unwrap_or(fields.len() - 1)];
            if key.is_empty() {
                continue;
            }
            let field = |i: usize| indices.fields[i].map(|index| fields[index]);
            if let (Some(pos), Some(hgvsp)) = (field(0), field(1)) {
                if !pos.is_empty() && !hgvsp.is_empty() {
                    out.aa(
                        key,
                        pos,
                        hgvsp,
                        std::slice::from_ref(&sig),
                        std::slice::from_ref(&review),
                    );
                }
            }
            let (exon, intron, hgvsc, consequence) = (field(5), field(6), field(7), field(2));
            let splice_ai = indices
                .splice
                .iter()
                .filter_map(|(key, i)| {
                    fields
                        .get(*i)
                        .filter(|v| !v.is_empty())
                        .map(|v| (key.clone(), (*v).into()))
                })
                .collect();
            if (exon.is_some_and(|s| !s.is_empty()) || intron.is_some_and(|s| !s.is_empty()))
                && sig.contains("athogenic")
                && CONFIDENCE.contains(&review.as_str())
            {
                let alt = decoded.fields.get(4).copied().unwrap_or(".");
                out.splice(
                    key,
                    vec![SpliceRecord {
                        chrom: decoded.fields[0].into(),
                        pos: decoded.position,
                        reference: decoded
                            .fields
                            .get(3)
                            .copied()
                            .filter(|s| !s.is_empty())
                            .unwrap_or(".")
                            .into(),
                        alternative: if alt == "." {
                            None
                        } else {
                            Some(
                                alt.split(',')
                                    .next()
                                    .filter(|s| !s.is_empty())
                                    .unwrap_or(".")
                                    .into(),
                            )
                        },
                        exon: exon.map(str::to_owned),
                        intron: intron.map(str::to_owned),
                        hgvsc: hgvsc.map(str::to_owned),
                        consequence: consequence.map(str::to_owned),
                        clinvar_sig: sig.clone(),
                        clinvar_review: review.clone(),
                        splice_ai,
                    }],
                );
            }
        }
        Ok(())
    };
    match work() {
        Err(error)
            if error.to_string().starts_with("KeyError:")
                || error.to_string().starts_with("ValueError:") =>
        {
            Ok(())
        }
        other => other,
    }
}

fn chromosome(
    vcf: &Path,
    index: &TabixIndex,
    h: &Header,
    indices: &Indices,
    chrom: &str,
) -> Result<Pair> {
    if !h.contigs.iter().any(|s| s == chrom) {
        return Err(invalid("unknown header contig"));
    }
    let mut out = Collector::default();
    let mut types = h.info.clone();
    fetch_contig_lines(vcf, index, chrom, &mut |line| {
        record(line, h, &mut types, indices, &mut out).map_err(|e| e.to_string())
    })
    .map_err(invalid)?;
    Ok(out.values)
}
fn provenance(inputs: &[PathBuf], mode: &str) -> Result<Vec<Provenance>> {
    let mut result = vec![Provenance {
        key: "builder".into(),
        value: format!(
            "priva clinvar-aa-splice {mode} {}",
            env!("CARGO_PKG_VERSION")
        ),
    }];
    for (i, path) in inputs.iter().enumerate() {
        result.push(Provenance {
            key: format!("source_{i}"),
            value: path
                .file_name()
                .and_then(|s| s.to_str())
                .ok_or_else(|| invalid("input needs a UTF-8 basename"))?
                .into(),
        });
        result.push(Provenance {
            key: format!("source_{i}_bytes"),
            value: std::fs::metadata(path)?.len().to_string(),
        });
    }
    Ok(result)
}

fn protect_input(input: &Path, output: &Path) -> Result<()> {
    // A nonexistent output cannot alias an existing input. Do not validate its
    // parent here: Python writes AA before discovering an invalid splice parent.
    if output.exists() {
        ensure_distinct_paths(input, output)?;
    }
    Ok(())
}

fn separate_outputs(first: &Path, second: &Path) -> Result<()> {
    if first.exists() && second.exists() {
        return ensure_distinct_paths(first, second);
    }
    let resolved = |p: &Path| {
        std::fs::canonicalize(p).ok().or_else(|| {
            let parent = p
                .parent()
                .filter(|p| !p.as_os_str().is_empty())
                .unwrap_or(Path::new("."));
            std::fs::canonicalize(parent)
                .ok()
                .zip(p.file_name())
                .map(|(dir, name)| dir.join(name))
        })
    };
    if first == second
        || resolved(first)
            .zip(resolved(second))
            .is_some_and(|(a, b)| a == b)
    {
        return Err(invalid("AA and splice output paths must be distinct"));
    }
    Ok(())
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args.is_empty() || args.iter().any(|s| s == "--help" || s == "-h") {
        println!("{HELP}");
        return Ok(());
    }
    let mode = &args[0];
    if !["build", "chromosome", "merge"].contains(&mode.as_str()) {
        return Err(invalid(HELP));
    }
    let mut options = HashMap::new();
    let (mut aa_inputs, mut splice_inputs) = (Vec::new(), Vec::new());
    let mut at = 1;
    while at < args.len() {
        let key = args[at].as_str();
        let value = args.get(at + 1).ok_or_else(|| invalid(HELP))?;
        match key {
            "--aa-input" if mode == "merge" => aa_inputs.push(PathBuf::from(value)),
            "--splice-input" if mode == "merge" => splice_inputs.push(PathBuf::from(value)),
            "--vcf" | "--chrom" | "--aa-output" | "--splice-output" | "--threads" => {
                options.insert(key, value.as_str());
            }
            _ => return Err(invalid(HELP)),
        }
        at += 2;
    }
    let required = |key| options.get(key).copied().ok_or_else(|| invalid(HELP));
    let aa_output = required("--aa-output")?;
    let splice_output = options.get("--splice-output").copied().unwrap_or("");
    let outputs: Vec<_> = [aa_output, splice_output]
        .into_iter()
        .filter(|s| !s.is_empty())
        .map(Path::new)
        .collect();
    if outputs.len() == 2 {
        separate_outputs(outputs[0], outputs[1])?;
    }
    let mut inputs = Vec::new();
    let pair = if mode == "merge" {
        if aa_inputs.is_empty() || aa_inputs.len() != splice_inputs.len() {
            return Err(invalid(
                "merge requires equal, ordered AA/splice input pairs",
            ));
        }
        let mut merged = Collector::default();
        for (aa, splice) in aa_inputs.iter().zip(&splice_inputs) {
            for input in [aa, splice] {
                for output in &outputs {
                    protect_input(input, output)?;
                }
                inputs.push(input.clone());
            }
            let aa = ClinvarAaCache::open(&aa.to_string_lossy())?;
            let splice = SpliceCache::open(&splice.to_string_lossy())?;
            merged.merge(Pair {
                aa: aa
                    .keys()
                    .iter()
                    .map(|k| (k.as_str().into(), aa.plain(k.as_str()).unwrap()))
                    .collect(),
                splice: splice
                    .keys()
                    .iter()
                    .map(|k| (k.as_str().into(), splice.plain(k.as_str()).unwrap()))
                    .collect(),
            });
        }
        merged.values
    } else {
        let vcf = Path::new(required("--vcf")?);
        let threads: i64 = options
            .get("--threads")
            .map(|s| {
                if !crate::python_text::python_int_is_valid(s) {
                    return Err(invalid("ValueError: invalid number of processes"));
                }
                let ascii: String = s
                    .chars()
                    .map(|c| {
                        crate::python_text::decimal_digit_value(c)
                            .map(|d| char::from(b'0' + d as u8))
                            .unwrap_or(c)
                    })
                    .collect();
                crate::python_text::python_int(&ascii)
                    .and_then(|i| i64::try_from(i).ok())
                    .ok_or_else(|| invalid("number of processes exceeds supported range"))
            })
            .transpose()?
            .unwrap_or(4);
        let mut h = header(vcf)?;
        let indices = indices(vcf)?;
        let indexed = index(vcf);
        if let Ok(index) = &indexed {
            for name in &index.names {
                if !h.contigs.contains(name) {
                    h.contigs.push(name.clone());
                }
            }
        }
        if threads <= 0 {
            return Err(invalid("number of processes must be at least 1"));
        }
        for output in &outputs {
            protect_input(vcf, output)?;
        }
        inputs.push(vcf.to_path_buf());
        if h.contigs.is_empty() {
            Pair::default()
        } else {
            let index = indexed?;
            if mode == "chromosome" {
                chromosome(vcf, &index, &h, &indices, required("--chrom")?)?
            } else {
                let next = AtomicUsize::new(0);
                let (send, receive) = std::sync::mpsc::channel();
                std::thread::scope(|scope| {
                    for _ in 0..(threads as usize).min(h.contigs.len()) {
                        let (send, h, index, indices, next) =
                            (send.clone(), &h, &index, &indices, &next);
                        scope.spawn(move || {
                            loop {
                                let i = next.fetch_add(1, Ordering::Relaxed);
                                if i >= h.contigs.len() {
                                    break;
                                }
                                let result = chromosome(vcf, index, h, indices, &h.contigs[i])
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
                let mut merged = Collector::default();
                for (_, result) in results {
                    merged.merge(result.map_err(invalid)?);
                }
                merged.values
            }
        }
    };
    let (aa_count, splice_count) = (pair.aa.len(), pair.splice.len());
    // Keep Python's AA-first publication order, including a later splice path
    // failure. Checked native archives use atomic publication individually.
    if !aa_output.is_empty() {
        write_clinvar_aa(pair.aa, provenance(&inputs, mode)?, Path::new(aa_output))?;
    }
    if !splice_output.is_empty() {
        write_records(
            Path::new(splice_output),
            pair.splice,
            provenance(&inputs, mode)?,
        )?;
    }
    println!("{{\"aa_transcripts\":{aa_count},\"splice_transcripts\":{splice_count}}}");
    Ok(())
}
