//! Streaming, annotation-only projection of the published SpliceAI2 top-three schema.
//!
//! This deliberately lossy VEP view does not replace the immutable 84-column source.
//! Gene associations are adapter-derived from caller-verified, exact source-region /
//! GENCODE context, not publisher-supplied per-gene predictions. The deployment driver
//! must bind that context, original tables and executable to their provenance hashes.
//! Multiple context IDs with the same chromosome/symbol/strand retain every gene
//! interval but emit one symbol per native prediction row: the upstream prediction
//! is gene-unconditioned. Opposite-strand symbol contexts remain unsupported. Source
//! prediction rows are never coalesced, including duplicate allele/strand rows.
//! No ACMG threshold, CADD feature or existing scientific module is changed here.
//!
//! Schema audited against Illumina/SpliceAI2 1dcd2a1882a8223095e7529982ae70bcb658340a
//! and precomputed top-three tables at dataset revision
//! cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9. The latter's producer code identity is
//! not established by the former. Directional deltas use sign*(alt-ref), without
//! absolute value or clipping: they are in [-1,1], whereas ref/alt probabilities
//! are in [0,1]. Genomic distance signs are already restored by upstream and must
//! not be reflected again on the minus strand. All selected lexemes are preserved.

use crate::cli::{Result, invalid};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Write};

const CONTEXT_HEADER: &str = "gene_id\tchrom\tgene_start\tgene_end\tstrand\tgene_name";
const USAGE: &str = "usage: priva spliceai2-vcf --genes CONTEXT.tsv --chrom chr1 --kind snv|indel --input PATH_OR_- --output NEW_PATH_OR_- [--no-header]\nInput is uncompressed native 84-column SpliceAI2 TSV; every shard requires its TSV header.\nContext: gene_id,chrom,gene_start,gene_end,strand,gene_name (tab-separated; 1-based closed).\nThe deployment driver must verify exact source-region/GENCODE identities and retain all original fields.\nOutput is an annotation-only, lossy top0 VCF view; no scoring or ACMG activation is implied.\nOutput files must not exist. On error partial output is not a completed artifact.\n--no-header suppresses only the output VCF header for later chromosome shards.";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VariantKind {
    Snv,
    Indel,
}

impl VariantKind {
    fn as_str(self) -> &'static str {
        match self {
            Self::Snv => "snv",
            Self::Indel => "indel",
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct Statistics {
    pub input_rows: u64,
    pub output_records: u64,
}

#[derive(Debug)]
struct Gene {
    id: String,
    start: u64,
    end: u64,
    strand: u8,
    name: String,
}

/// Exact published 84-column schema, including unprojected fields.
pub fn native_columns() -> Vec<String> {
    let mut columns = ["chrom", "pos", "ref", "alt", "strand"]
        .into_iter()
        .map(str::to_owned)
        .collect::<Vec<_>>();
    for effect in ["donor_gain", "donor_loss", "acceptor_gain", "acceptor_loss"] {
        for field in ["delta_score", "ref_score", "alt_score", "dist"] {
            for rank in 0..3 {
                columns.push(format!("{effect}_{field}_{rank}"));
            }
        }
    }
    for effect in ["jxn_gain", "jxn_loss"] {
        for field in [
            "delta_score",
            "ref_score",
            "alt_score",
            "donor_dist",
            "acceptor_dist",
        ] {
            for rank in 0..3 {
                columns.push(format!("{effect}_{field}_{rank}"));
            }
        }
    }
    columns.push("spliceai2_summary_score".to_owned());
    columns
}

fn canonical_chrom(chrom: &str) -> bool {
    matches!(chrom, "chrX" | "chrY")
        || chrom.strip_prefix("chr").is_some_and(|suffix| {
            suffix
                .parse::<u8>()
                .is_ok_and(|n| (1..=22).contains(&n) && suffix == n.to_string())
        })
}

fn positive_coordinate(value: &str) -> Result<u64> {
    if value.starts_with('0') || value.is_empty() || !value.bytes().all(|c| c.is_ascii_digit()) {
        return Err(invalid(format!(
            "expected canonical positive integer coordinate, got {value:?}"
        )));
    }
    Ok(value.parse()?)
}

fn strand(value: &str) -> Result<u8> {
    match value {
        "+" => Ok(b'+'),
        "-" => Ok(b'-'),
        _ => Err(invalid(format!(
            "expected model/gene strand + or -, got {value:?}"
        ))),
    }
}

fn vcf_token(value: &str) -> bool {
    !value.is_empty()
        && value != "."
        && value
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-' | b'.'))
}

fn read_line<R: BufRead>(reader: &mut R, line: &mut String) -> io::Result<bool> {
    line.clear();
    if reader.read_line(line)? == 0 {
        return Ok(false);
    }
    if line.ends_with('\n') {
        line.pop();
        if line.ends_with('\r') {
            line.pop();
        }
    }
    Ok(true)
}

fn genes_for_chrom<G: BufRead>(mut reader: G, chrom: &str) -> Result<Vec<Gene>> {
    let mut line = String::new();
    if !read_line(&mut reader, &mut line)? || line != CONTEXT_HEADER {
        return Err(invalid(format!(
            "gene context requires exact header {CONTEXT_HEADER:?}"
        )));
    }
    let mut ids = HashSet::new();
    let mut symbols = HashMap::new();
    let mut genes = Vec::new();
    let mut line_number = 1;
    while read_line(&mut reader, &mut line)? {
        line_number += 1;
        let fields: Vec<_> = line.split('\t').collect();
        let result = (|| -> Result<()> {
            if fields.len() != 6 {
                return Err(invalid("expected six gene context columns"));
            }
            if !vcf_token(fields[0]) || !vcf_token(fields[5]) || !canonical_chrom(fields[1]) {
                return Err(invalid(
                    "unsupported gene ID, gene symbol or chromosome in context",
                ));
            }
            if !ids.insert(fields[0].to_owned()) {
                return Err(invalid(format!("duplicate gene context ID {}", fields[0])));
            }
            // The official VEP plugin keys matching records by symbol alone.
            // Multiple same-strand IDs refer to one gene-unconditioned prediction
            // at a source POS; retain all intervals and coalesce output labels only.
            // Opposite strands can carry different predictions, so reject those
            // symbol contexts even when their intervals do not overlap.
            // Across chromosomes (including PAR X/Y symbols) there is no collision.
            let gene_strand = strand(fields[4])?;
            if let Some((prior_id, prior_strand)) = symbols.insert(
                (fields[1].to_owned(), fields[5].to_owned()),
                (fields[0].to_owned(), gene_strand),
            ) {
                if prior_strand != gene_strand {
                    return Err(invalid(format!(
                        "ambiguous symbol {} on {} has opposite-strand context IDs {prior_id} and {}",
                        fields[5], fields[1], fields[0]
                    )));
                }
            }
            let start = positive_coordinate(fields[2])?;
            let end = positive_coordinate(fields[3])?;
            if start > end {
                return Err(invalid("gene_start exceeds gene_end"));
            }
            if fields[1] == chrom {
                genes.push(Gene {
                    id: fields[0].to_owned(),
                    start,
                    end,
                    strand: gene_strand,
                    name: fields[5].to_owned(),
                });
            }
            Ok(())
        })();
        result.map_err(|e| invalid(format!("gene context line {line_number}: {e}")))?;
    }
    if genes.is_empty() {
        return Err(invalid(format!(
            "gene context contains no genes for {chrom}"
        )));
    }
    genes.sort_by(|a, b| (a.start, a.end, &a.id).cmp(&(b.start, b.end, &b.id)));
    Ok(genes)
}

fn probability(value: &str, signed: bool) -> Result<f64> {
    if value.is_empty() || value.bytes().any(|b| b.is_ascii_whitespace()) {
        return Err(invalid(format!(
            "missing or whitespace-containing score {value:?}"
        )));
    }
    let number: f64 = value
        .parse()
        .map_err(|_| invalid(format!("invalid score {value:?}")))?;
    let minimum = if signed { -1.0 } else { 0.0 };
    if !number.is_finite() || number < minimum || number > 1.0 {
        return Err(invalid(format!(
            "score {value:?} is not finite in [{minimum},1]"
        )));
    }
    Ok(number)
}

fn distance(value: &str) -> Result<()> {
    let digits = value.strip_prefix('-').unwrap_or(value);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return Err(invalid(format!(
            "distance must be a signed decimal integer, got {value:?}"
        )));
    }
    // Official score.py casts genomic distances to int32 before serializing.
    value
        .parse::<i32>()
        .map_err(|_| invalid(format!("distance outside int32: {value:?}")))?;
    Ok(())
}

fn write_header<W: Write>(writer: &mut W) -> Result<()> {
    writeln!(writer, "##fileformat=VCFv4.2")?;
    writeln!(writer, "##reference=GRCh38")?;
    writeln!(writer, "##source=PriVA_SpliceAI2_projection")?;
    writeln!(writer, "##SpliceAI2Version=2")?;
    writeln!(
        writer,
        "##SpliceAI2Projection=top0_site_deltas_and_genomic_distances_from_native_top3"
    )?;
    writeln!(writer, "##SpliceAI2ProjectionLossy=true")?;
    writeln!(
        writer,
        "##SpliceAI2OriginalFields=84_retained_in_immutable_source_tables"
    )?;
    writeln!(
        writer,
        "##SpliceAI2GeneAssociation=adapter_derived_exact_source_region_GENCODE_context_and_model_strand"
    )?;
    writeln!(writer, "##SpliceAI2PublisherGeneSpecificPrediction=false")?;
    writeln!(
        writer,
        "##SpliceAI2SymbolPolicy=coalesce_same_strand_context_labels_per_source_row_keep_all_gene_ID_intervals"
    )?;
    writeln!(writer, "##SpliceAI2RecommendedThresholds=0.1,0.25,0.5")?;
    writeln!(
        writer,
        "##SpliceAI2Policy=annotation_only_no_ACMG_or_CADD_activation"
    )?;
    for number in 1..=22 {
        writeln!(writer, "##contig=<ID=chr{number}>")?;
    }
    writeln!(writer, "##contig=<ID=chrX>\n##contig=<ID=chrY>")?;
    writeln!(
        writer,
        "##INFO=<ID=SpliceAI,Number=.,Type=String,Description=\"SpliceAI2 top0 projection; adapter-derived gene association, not publisher per-gene prediction. ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL\">"
    )?;
    writeln!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
    Ok(())
}

/// Stream one chromosome shard. The caller must discard partial output on error.
///
/// Memory holds gene intervals plus the current coordinate's allele/strand keys,
/// not all variant rows. POS, alleles and selected score/distance strings are copied
/// unchanged. Every source field is validated, even those absent from this view.
pub fn project<G: BufRead, R: BufRead, W: Write>(
    genes: G,
    mut input: R,
    mut output: W,
    chrom: &str,
    kind: VariantKind,
    with_header: bool,
) -> Result<Statistics> {
    if !canonical_chrom(chrom) {
        return Err(invalid("--chrom must be chr1..chr22, chrX or chrY"));
    }
    let genes = genes_for_chrom(genes, chrom)?;
    let columns = native_columns();
    let mut line = String::new();
    if !read_line(&mut input, &mut line)? || line != columns.join("\t") {
        return Err(invalid(
            "input requires the exact native 84-column SpliceAI2 top-three header",
        ));
    }
    if with_header {
        write_header(&mut output)?;
    }
    let mut stats = Statistics {
        input_rows: 0,
        output_records: 0,
    };
    let mut previous_pos = 0;
    let mut position_keys = HashSet::new();
    let mut next_gene = 0;
    let mut active: Vec<&Gene> = Vec::new();
    while read_line(&mut input, &mut line)? {
        let line_number = stats.input_rows + 2;
        let result = (|| -> Result<()> {
            let fields: Vec<_> = line.split('\t').collect();
            if fields.len() != 84 {
                return Err(invalid(format!(
                    "expected 84 columns, found {}",
                    fields.len()
                )));
            }
            if fields[0] != chrom {
                return Err(invalid(format!(
                    "chromosome {} differs from requested {chrom}",
                    fields[0]
                )));
            }
            let pos = positive_coordinate(fields[1])?;
            let model_strand = strand(fields[4])?;
            if pos < previous_pos {
                return Err(invalid("coordinates are not ascending"));
            }
            if pos != previous_pos {
                position_keys.clear();
            }
            previous_pos = pos;
            for allele in &fields[2..4] {
                if allele.is_empty()
                    || !allele
                        .bytes()
                        .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
                {
                    return Err(invalid(
                        "alleles must be nonempty uppercase A/C/G/T sequences",
                    ));
                }
            }
            if fields[2] == fields[3] {
                return Err(invalid("REF equals ALT"));
            }
            match kind {
                VariantKind::Snv if fields[2].len() != 1 || fields[3].len() != 1 => {
                    return Err(invalid("snv input requires one-base REF and ALT"));
                }
                VariantKind::Indel if fields[2].len() == fields[3].len() => {
                    return Err(invalid("indel input requires unequal REF/ALT lengths"));
                }
                _ => (),
            }
            if !position_keys.insert((fields[2].to_owned(), fields[3].to_owned(), model_strand)) {
                return Err(invalid(
                    "duplicate allele + model strand at this coordinate",
                ));
            }
            for column in 5..83 {
                let name = &columns[column];
                let check = if name.contains("_dist_") {
                    distance(fields[column])
                } else {
                    probability(fields[column], name.contains("_delta_score_")).map(|_| ())
                };
                check.map_err(|e| invalid(format!("{}: {e}", columns[column])))?;
            }
            let summary = probability(fields[83], true)?;
            // Site block offsets: donor gain 5, donor loss 17, acceptor gain 29,
            // acceptor loss 41. Junctions deliberately do not enter the summary.
            let maximum = [5, 17, 29, 41]
                .into_iter()
                .map(|index| {
                    fields[index]
                        .parse::<f64>()
                        .expect("already validated score")
                })
                .fold(f64::NEG_INFINITY, f64::max);
            if summary != maximum {
                return Err(invalid(
                    "spliceai2_summary_score differs from numeric max of the four top0 site deltas",
                ));
            }
            while next_gene < genes.len() && genes[next_gene].start <= pos {
                active.push(&genes[next_gene]);
                next_gene += 1;
            }
            active.retain(|gene| gene.end >= pos);
            let mut matched = BTreeSet::new();
            for gene in &active {
                if gene.strand == model_strand {
                    matched.insert(&gene.name);
                }
            }
            if matched.is_empty() {
                return Err(invalid(
                    "no verified gene interval contains source POS on the model strand",
                ));
            }
            for name in matched {
                writeln!(
                    output,
                    "{}\t{}\t.\t{}\t{}\t.\t.\tSpliceAI={}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
                    fields[0],
                    fields[1],
                    fields[2],
                    fields[3],
                    fields[3],
                    name,
                    fields[29],
                    fields[41],
                    fields[5],
                    fields[17],
                    fields[38],
                    fields[50],
                    fields[14],
                    fields[26]
                )?;
                stats.output_records += 1;
            }
            stats.input_rows += 1;
            Ok(())
        })();
        result.map_err(|e| invalid(format!("native TSV line {line_number}: {e}")))?;
    }
    output.flush()?;
    Ok(stats)
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args == ["--help"] {
        println!("{USAGE}");
        return Ok(());
    }
    let mut values = HashMap::new();
    let mut no_header = false;
    let mut index = 0;
    while index < args.len() {
        let option = args[index].as_str();
        if option == "--no-header" && !no_header {
            no_header = true;
            index += 1;
            continue;
        }
        if !matches!(
            option,
            "--genes" | "--chrom" | "--kind" | "--input" | "--output"
        ) {
            return Err(invalid(format!(
                "unknown or repeated option {option:?}\n{USAGE}"
            )));
        }
        let value = args.get(index + 1).ok_or_else(|| invalid(USAGE))?;
        if value.starts_with("--") || values.insert(option, value.as_str()).is_some() {
            return Err(invalid(format!(
                "missing value or repeated option {option:?}\n{USAGE}"
            )));
        }
        index += 2;
    }
    for option in ["--genes", "--chrom", "--kind", "--input", "--output"] {
        if !values.contains_key(option) {
            return Err(invalid(format!("missing {option}\n{USAGE}")));
        }
    }
    let kind = match values["--kind"] {
        "snv" => VariantKind::Snv,
        "indel" => VariantKind::Indel,
        _ => return Err(invalid("--kind must be snv or indel")),
    };
    let genes = BufReader::new(File::open(values["--genes"])?);
    let input: Box<dyn BufRead> = if values["--input"] == "-" {
        Box::new(BufReader::new(io::stdin()))
    } else {
        Box::new(BufReader::new(File::open(values["--input"])?))
    };
    let output: Box<dyn Write> = if values["--output"] == "-" {
        Box::new(BufWriter::new(io::stdout()))
    } else {
        Box::new(BufWriter::new(
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(values["--output"])?,
        ))
    };
    let stats = project(genes, input, output, values["--chrom"], kind, !no_header)?;
    // chrom and kind have been validated against a fixed ASCII vocabulary.
    eprintln!(
        "{{\"module\":\"spliceai2-vcf\",\"chrom\":\"{}\",\"kind\":\"{}\",\"input_rows\":{},\"output_records\":{}}}",
        values["--chrom"],
        kind.as_str(),
        stats.input_rows,
        stats.output_records
    );
    Ok(())
}
