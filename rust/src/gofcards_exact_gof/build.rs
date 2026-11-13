//! The reference check, the table readers, the HGNC and review rules, and
//! the join and export of `build_cache` and `_nest`.
use super::fasta::Fasta;
use super::fetch::Source;
use super::*;
use crate::gene_nonlof_mechanism::get;
use crate::hpo_condition_mechanism::tsv::Reader;
use crate::pandas_table::Cell;
use std::collections::BTreeMap;

const ELIGIBLE: &str = "eligible";
const QUARANTINE_GENE: &str = "quarantined_gene_discordance";
const QUARANTINE_SIBLING: &str = "quarantined_allele_gene_discordance";
const QUARANTINE_REVIEW_REQUIRED: &str = "quarantined_mechanism_review_required";
const REVIEW_COLUMNS: [&str; 13] = [
    "source_order",
    "source_gene",
    "allele_key",
    "source_mechanism",
    "reviewed_mechanism",
    "review_status",
    "gof_eligibility",
    "reason_code",
    "evidence_summary",
    "pmids",
    "reviewer",
    "reviewed_at",
    "review_version",
];

fn reviewed_quarantine(mechanism: &str) -> &'static str {
    match mechanism {
        "LOF" => "quarantined_reviewed_lof",
        "MIXED" => "quarantined_reviewed_mixed",
        "DOMINANT_NEGATIVE" | "DN" => "quarantined_reviewed_dominant_negative",
        "UNCERTAIN" => "quarantined_reviewed_uncertain",
        "EXCLUDE" => "quarantined_reviewed_exclude",
        _ => QUARANTINE_REVIEW_REQUIRED,
    }
}

/// One GRCh37 record that passed the reference check.
pub(super) struct Core {
    pub(super) vcf_id: String,
    allele_key: String,
    chrom: String,
    contig: String,
    pos: i128,
    reference: String,
    alternate: String,
    status: &'static str,
    source_pos: String,
    source_ref: String,
    source_alt: String,
}

enum Representation {
    Ok {
        kind: &'static str,
        pos: i128,
        reference: String,
        alternate: String,
    },
    Rejected(&'static str, String),
}

fn mismatch(observed: &str, reference: &str) -> Representation {
    Representation::Rejected(
        "ref_mismatch",
        format!("GRCh37 has {observed}, source claims {reference}"),
    )
}

/// `_vcf_representation`, its `try` as the `Err` of the inner closure.
fn vcf_representation(
    fasta: &Fasta,
    chrom: &str,
    start: i128,
    reference: &str,
    alternate: &str,
) -> Representation {
    let width = reference.chars().count().max(1) as i128;
    let decide = || -> std::result::Result<Representation, String> {
        if !reference.is_empty() && !alternate.is_empty() {
            let observed = fasta.slice(chrom, start, width)?;
            if observed != reference {
                return Ok(mismatch(&observed, reference));
            }
            return Ok(Representation::Ok {
                kind: "raw_ref_alt",
                pos: start,
                reference: reference.to_owned(),
                alternate: alternate.to_owned(),
            });
        }
        if !reference.is_empty() {
            let anchor = fasta.slice(chrom, start - 1, 1)?;
            let observed = fasta.slice(chrom, start, width)?;
            if observed != reference {
                return Ok(mismatch(&observed, reference));
            }
            return Ok(Representation::Ok {
                kind: "deletion_padded_ref_match",
                pos: start - 1,
                reference: format!("{anchor}{reference}"),
                alternate: anchor,
            });
        }
        if !alternate.is_empty() {
            let anchor = fasta.slice(chrom, start, 1)?;
            return Ok(Representation::Ok {
                kind: "insertion_padded",
                pos: start,
                alternate: format!("{anchor}{alternate}"),
                reference: anchor,
            });
        }
        Ok(Representation::Rejected(
            "missing_ref_and_alt",
            "source has neither REF nor ALT".to_owned(),
        ))
    };
    decide().unwrap_or_else(|detail| Representation::Rejected("fasta_lookup_error", detail))
}

/// A field as `csv.writer` writes it with a tab delimiter.
fn tsv_field(value: &str) -> String {
    if value.contains(['\t', '"', '\n', '\r']) {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_owned()
    }
}

/// `build_source_vcf`.
pub(super) fn build_source_vcf(
    sources: &[Source],
    hg19_fasta: &Path,
    out_vcf: &Path,
    rejects_tsv: &Path,
) -> Result<Vec<Core>> {
    let fasta = Fasta::open(hg19_fasta)?;
    let mut seen = HashSet::new();
    let unique = sources
        .iter()
        .filter(|source| seen.insert(source.allele_key.as_str()));
    let mut kept: Vec<Core> = Vec::new();
    let mut rejected: Vec<[String; 3]> = Vec::new();
    for (index, row) in unique.enumerate() {
        if row.chrom.is_empty() || !is_digits(&row.start) {
            rejected.push([
                row.allele_key.clone(),
                "missing_coordinate".to_owned(),
                format!(
                    "chrom={} start={}",
                    crate::python_path::repr(&row.chrom),
                    crate::python_path::repr(&row.start)
                ),
            ]);
            continue;
        }
        let start = python_int(&row.start)?;
        match vcf_representation(&fasta, &row.chrom, start, &row.reference, &row.alternate) {
            Representation::Rejected(reason, detail) => {
                rejected.push([row.allele_key.clone(), reason.to_owned(), detail]);
            }
            Representation::Ok {
                kind,
                pos,
                reference,
                alternate,
            } => kept.push(Core {
                vcf_id: format!("gof{index:06}"),
                allele_key: row.allele_key.clone(),
                chrom: row.chrom.clone(),
                contig: fasta.contig_name(&row.chrom).unwrap_or_default(),
                pos,
                reference,
                alternate,
                status: kind,
                source_pos: row.start.clone(),
                source_ref: row.reference.clone(),
                source_alt: row.alternate.clone(),
            }),
        }
    }
    let contigs: std::collections::BTreeSet<&str> =
        kept.iter().map(|record| record.contig.as_str()).collect();
    if let Some(parent) = out_vcf.parent() {
        fs::create_dir_all(parent).map_err(|e| os_error(e, parent, None))?;
    }
    let mut body = String::from("##fileformat=VCFv4.2\n");
    for name in contigs {
        body.push_str(&format!(
            "##contig=<ID={name},length={}>\n",
            fasta.length(name)
        ));
    }
    body.push_str("##INFO=<ID=GOFID,Number=1,Type=String,Description=\"GoFCards record id\">\n");
    body.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    let mut ordered: Vec<&Core> = kept.iter().collect();
    ordered.sort_by(|a, b| (&a.contig, a.pos).cmp(&(&b.contig, b.pos)));
    for record in ordered {
        body.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t.\t.\tGOFID={}\n",
            record.contig,
            record.pos,
            record.vcf_id,
            record.reference,
            record.alternate,
            record.vcf_id
        ));
    }
    fs::write(out_vcf, body).map_err(|e| os_error(e, out_vcf, None))?;

    let mut table = String::new();
    if rejected.is_empty() {
        table.push('\n');
    } else {
        table.push_str("allele_key\treject_reason\tdetail\n");
        for row in &rejected {
            let fields: Vec<String> = row.iter().map(|field| tsv_field(field)).collect();
            table.push_str(&fields.join("\t"));
            table.push('\n');
        }
    }
    fs::write(rejects_tsv, table).map_err(|e| os_error(e, rejects_tsv, None))?;
    log(&format!(
        "reference check: {} alleles kept, {} rejected -> {}",
        kept.len(),
        rejected.len(),
        rejects_tsv.display()
    ));
    Ok(kept)
}

/// One row of `parse_vep`.
pub(super) struct VepRow {
    vcf_id: String,
    assembly: &'static str,
    transcript: String,
    consequence: String,
    canonical: String,
    mane_select: String,
    symbol: String,
    hgvsc: String,
    hgvsp: String,
}

/// `parse_vep`.
pub(super) fn parse_vep(tsv: &Path, assembly: &'static str) -> Result<Vec<VepRow>> {
    let lines = read_lines(tsv)?;
    let kept: String = lines
        .iter()
        .filter(|line| !line.starts_with("##"))
        .map(String::as_str)
        .collect();
    let mut reader = Reader::from_text(kept);
    reader.require(tsv, &[])?;
    let mut rows = Vec::new();
    while let Some(record) = reader.next_row()? {
        let field = |name: &str| {
            record
                .get(name)
                .and_then(Option::as_deref)
                .map_or_else(String::new, clean_text)
        };
        let transcript = field("Feature");
        if transcript.is_empty() {
            continue;
        }
        rows.push(VepRow {
            vcf_id: field("#Uploaded_variation"),
            assembly,
            transcript,
            consequence: field("Consequence"),
            canonical: field("CANONICAL"),
            mane_select: field("MANE_SELECT"),
            symbol: field("SYMBOL"),
            hgvsc: field("HGVSc"),
            hgvsp: field("HGVSp"),
        });
    }
    let records: HashSet<&str> = rows.iter().map(|row| row.vcf_id.as_str()).collect();
    log(&format!(
        "VEP {assembly}: {} transcript rows over {} records",
        rows.len(),
        records.len()
    ));
    Ok(rows)
}

/// CHROM/POS/REF/ALT of one record, keyed by id.
pub(super) struct Coord {
    pub(super) vcf_id: String,
    chrom: String,
    pos: String,
    reference: String,
    alternate: String,
}

/// `read_vcf_coordinates`.
pub(super) fn read_vcf_coordinates(vcf: &Path) -> Result<Vec<Coord>> {
    let mut rows = Vec::new();
    for line in read_lines(vcf)? {
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
        if fields.len() < 5 || fields[2].is_empty() {
            continue;
        }
        rows.push(Coord {
            vcf_id: fields[2].to_owned(),
            chrom: normalize_chrom(&clean_text(fields[0])),
            pos: fields[1].to_owned(),
            reference: fields[3].to_uppercase(),
            alternate: fields[4].to_uppercase(),
        });
    }
    Ok(rows)
}

/// `endpoint_hg38_coordinates`.
pub(super) fn endpoint_hg38_coordinates(
    core: &[&Core],
    summary: &HashMap<String, Json>,
    hg38_fasta: &Path,
) -> Result<(Vec<Coord>, Vec<(&'static str, i128)>)> {
    let fasta = Fasta::open(hg38_fasta)?;
    let mut rows = Vec::new();
    let mut counts: Vec<(&'static str, i128)> = Vec::new();
    let mut count = |name: &'static str| match counts.iter_mut().find(|(key, _)| *key == name) {
        Some((_, value)) => *value += 1,
        None => counts.push((name, 1)),
    };
    for record in core {
        let position = match summary.get(&record.allele_key) {
            Some(data) if crate::gene_nonlof_mechanism::truthy(data) => {
                normalize_int_json(get(data, "hg38_start")?)
            }
            _ => String::new(),
        };
        if !is_digits(&position) {
            count("no_endpoint_position");
            continue;
        }
        let width = record.reference.chars().count() as i128;
        let observed = crate::python_text::python_int(&position)
            .ok_or_else(String::new)
            .and_then(|start| fasta.slice(&record.chrom, start, width));
        match observed {
            Err(_) => count("endpoint_position_off_the_contig"),
            Ok(observed) if observed != record.reference => {
                count("endpoint_position_rejected_by_reference_check");
            }
            Ok(_) => {
                count("endpoint_position_accepted");
                rows.push(Coord {
                    vcf_id: record.vcf_id.clone(),
                    chrom: record.chrom.clone(),
                    pos: position,
                    reference: record.reference.clone(),
                    alternate: record.alternate.clone(),
                });
            }
        }
    }
    Ok((rows, counts))
}

/// `HgncResolver`.
pub(super) struct HgncResolver {
    approved: HashSet<String>,
    aliases: HashMap<String, HashSet<String>>,
    ids: HashMap<String, String>,
}

impl HgncResolver {
    pub(super) fn new(table: &Path) -> Result<Self> {
        if !table.exists() {
            return Err(error(
                "RuntimeError",
                format!(
                    "HGNC complete-set table is required but was not found: {}. \
Install it rather than skipping symbol normalization.",
                    table.display()
                ),
            ));
        }
        let bytes = fs::read(table).map_err(|e| os_error(e, table, None))?;
        let body =
            crate::python_io::decode_utf8(&bytes).map_err(|e| error("UnicodeDecodeError", e))?;
        let mut reader = Reader::from_text(body.to_owned());
        reader.require(table, &[])?;
        let mut resolver = Self {
            approved: HashSet::new(),
            aliases: HashMap::new(),
            ids: HashMap::new(),
        };
        while let Some(row) = reader.next_row()? {
            let field = |name: &str| {
                row.get(name)
                    .and_then(Option::as_deref)
                    .map_or_else(String::new, clean_text)
            };
            let symbol = field("symbol").to_uppercase();
            if symbol.is_empty() {
                continue;
            }
            resolver.approved.insert(symbol.clone());
            let id = field("hgnc_id");
            if !id.is_empty() {
                resolver.ids.insert(symbol.clone(), id);
            }
            for name in ["alias_symbol", "prev_symbol"] {
                for alias in field(name).split(['|', ',']) {
                    let alias = crate::python_text::python_strip(alias).to_uppercase();
                    if !alias.is_empty() {
                        resolver
                            .aliases
                            .entry(alias)
                            .or_default()
                            .insert(symbol.clone());
                    }
                }
            }
        }
        log(&format!(
            "HGNC resolver loaded {} approved symbols from {}",
            resolver.approved.len(),
            table.display()
        ));
        Ok(resolver)
    }

    fn resolve(&self, value: &str) -> String {
        let query = clean_text(value).to_uppercase();
        if query.is_empty() || self.approved.contains(&query) {
            return query;
        }
        match self.aliases.get(&query) {
            Some(matches) if matches.len() == 1 => matches.iter().next().cloned().unwrap_or(query),
            _ => query,
        }
    }

    fn hgnc_id(&self, symbol: &str) -> Json {
        self.ids
            .get(&clean_text(symbol).to_uppercase())
            .map_or(Json::Null, |id| Json::Str(id.clone()))
    }
}

/// The review fields `build_cache` reads.
pub(super) struct Review {
    eligibility: String,
    mechanism: String,
    reason: String,
}

/// `load_mechanism_reviews`.
pub(super) fn load_mechanism_reviews(
    path: Option<&Path>,
) -> Result<HashMap<(String, String), Review>> {
    let Some(path) = path.filter(|path| path.exists()) else {
        log("no mechanism review table supplied; every allele keeps its source GOF claim");
        return Ok(HashMap::new());
    };
    let table = crate::gene_mechanism_builder::table::read(path, b'\t', false)?;
    let missing: Vec<&str> = REVIEW_COLUMNS
        .iter()
        .copied()
        .filter(|name| table.column_index(name).is_none())
        .collect();
    if !missing.is_empty() {
        return Err(error(
            "RuntimeError",
            format!(
                "{} is missing review columns: {}",
                path.display(),
                missing.join(", ")
            ),
        ));
    }
    let index = |name: &str| table.column_index(name).unwrap_or(0);
    let value = |column: usize, row: usize| match table.cell(column, row) {
        Cell::Text(text) => clean_text(text),
        _ => String::new(),
    };
    let [gene, allele, eligibility, mechanism, reason] = [
        index("source_gene"),
        index("allele_key"),
        index("gof_eligibility"),
        index("reviewed_mechanism"),
        index("reason_code"),
    ];
    let mut reviews = HashMap::new();
    for row in 0..table.rows {
        reviews.insert(
            (value(gene, row).to_uppercase(), value(allele, row)),
            Review {
                eligibility: value(eligibility, row),
                mechanism: value(mechanism, row),
                reason: value(reason, row),
            },
        );
    }
    log(&format!(
        "loaded {} reviewed alleles from {}",
        reviews.len(),
        path.display()
    ));
    Ok(reviews)
}

/// `"; ".join` of the distinct non-empty values, in order.
fn join_unique<'a>(values: impl Iterator<Item = &'a str>) -> String {
    let mut seen: Vec<&str> = Vec::new();
    for value in values {
        if !value.is_empty() && !seen.contains(&value) {
            seen.push(value);
        }
    }
    seen.join("; ")
}

pub(super) struct Layers<'a> {
    pub(super) sources: &'a [Source],
    pub(super) core: &'a [Core],
    pub(super) vep: [Vec<VepRow>; 2],
    pub(super) hg19_coords: Vec<Coord>,
    pub(super) hg38_coords: Vec<Coord>,
    pub(super) unmapped: &'a HashSet<String>,
    pub(super) reviews: &'a HashMap<(String, String), Review>,
    pub(super) summary: &'a HashMap<String, Json>,
    pub(super) resolver: &'a HgncResolver,
    pub(super) provenance: &'a [(String, Json)],
    pub(super) output: &'a Path,
}

/// One row of the merged frame.
struct Merged<'a> {
    core: &'a Core,
    hg19: Option<&'a Coord>,
    hg38: Option<&'a Coord>,
    vep: Option<&'a VepRow>,
    symbol: String,
    vep_symbol: String,
    gene_match_status: &'static str,
    eligibility: String,
    reason: String,
    hg38_status: &'static str,
    accession: String,
}

#[derive(Default)]
struct View {
    by_hgvsc: BTreeMap<String, Json>,
    by_hgvsp: BTreeMap<String, Vec<String>>,
}

struct Block {
    genomic: Json,
    transcripts: BTreeMap<String, View>,
}

struct Variant {
    record: Json,
    status: String,
    evidence: usize,
    assemblies: BTreeMap<&'static str, Block>,
}

struct Gene {
    hgnc_id: Json,
    variants: BTreeMap<String, Variant>,
}

fn object(items: Vec<(&str, Json)>) -> Json {
    Json::Object(
        items
            .into_iter()
            .map(|(key, value)| (key.to_owned(), value))
            .collect(),
    )
}

fn optional(value: &str) -> Json {
    if value.is_empty() {
        Json::Null
    } else {
        Json::Str(value.to_owned())
    }
}

fn string(value: &str) -> Json {
    Json::Str(value.to_owned())
}

/// `variant_id_of`, for text already normalized.
fn variant_id(chrom: &str, start: &str, reference: &str, alternate: &str) -> String {
    let allele = |text: &str| {
        let text = normalize_allele(&clean_text(text));
        if text.is_empty() {
            "-".to_owned()
        } else {
            text
        }
    };
    format!(
        "loc_{}:{}:{}->{}_grch37",
        normalize_chrom(&clean_text(chrom)),
        strip_point_zero(start),
        allele(reference),
        allele(alternate)
    )
}

fn coordinate_block(coord: Option<&Coord>, status: &str) -> Result<Block> {
    let (chrom, pos, reference, alternate) = coord.map_or(("", "", "", ""), |c| {
        (
            c.chrom.as_str(),
            c.pos.as_str(),
            c.reference.as_str(),
            c.alternate.as_str(),
        )
    });
    Ok(Block {
        genomic: object(vec![
            ("chrom", string(chrom)),
            ("pos", Json::Int(python_int(pos)?)),
            ("ref", string(reference)),
            ("alt", string(alternate)),
            ("status", string(status)),
        ]),
        transcripts: BTreeMap::new(),
    })
}

/// `_nest` without the metadata block.
fn nest(
    merged: &[Merged],
    evidence: &HashMap<&str, Vec<&Json>>,
    resolver: &HgncResolver,
    unmapped: &HashSet<String>,
) -> Result<BTreeMap<String, Gene>> {
    let mut by_gene: BTreeMap<&str, BTreeMap<&str, Vec<&Merged>>> = BTreeMap::new();
    for row in merged {
        by_gene
            .entry(&row.symbol)
            .or_default()
            .entry(&row.core.allele_key)
            .or_default()
            .push(row);
    }
    let mut genes = BTreeMap::new();
    for (symbol, alleles) in by_gene {
        if symbol.is_empty() {
            continue;
        }
        let gene = genes.entry(symbol.to_owned()).or_insert_with(|| Gene {
            hgnc_id: resolver.hgnc_id(symbol),
            variants: BTreeMap::new(),
        });
        for (source_key, rows) in alleles {
            let first = rows[0];
            let chrom = first.hg19.map_or("", |c| c.chrom.as_str());
            let identifier = variant_id(
                chrom,
                &first.core.source_pos,
                &first.core.source_ref,
                &first.core.source_alt,
            );
            let variant = gene.variants.entry(identifier).or_insert_with(|| {
                let records = evidence.get(source_key).cloned().unwrap_or_default();
                let dash = |text: &str| string(if text.is_empty() { "-" } else { text });
                Variant {
                    status: first.eligibility.clone(),
                    evidence: records.len(),
                    record: object(vec![
                        (
                            "eligibility",
                            object(vec![
                                ("status", string(&first.eligibility)),
                                ("gene_match_status", string(first.gene_match_status)),
                                ("reason", optional(&first.reason)),
                                ("vep_symbol", optional(&first.vep_symbol)),
                            ]),
                        ),
                        (
                            "source",
                            object(vec![
                                ("gofcards_allele_key", string(source_key)),
                                (
                                    "variant_type_label",
                                    string(source_key.split('|').next().unwrap_or("")),
                                ),
                                ("assembly", string("hg19")),
                                ("chrom", string(chrom)),
                                ("start", string(&first.core.source_pos)),
                                ("ref", dash(&first.core.source_ref)),
                                ("alt", dash(&first.core.source_alt)),
                            ]),
                        ),
                        (
                            "liftover_status",
                            string(if unmapped.contains(&first.core.vcf_id) {
                                "unmapped"
                            } else {
                                "mapped"
                            }),
                        ),
                        (
                            "annotations",
                            object(vec![("clinvar_variation_id", optional(&first.accession))]),
                        ),
                        (
                            "evidence",
                            Json::Array(records.iter().map(|item| (*item).clone()).collect()),
                        ),
                    ]),
                    assemblies: BTreeMap::new(),
                }
            });
            if first.hg38_status == "endpoint_ref_match" {
                variant
                    .assemblies
                    .insert("hg38", coordinate_block(first.hg38, first.hg38_status)?);
            }
            let mut by_assembly: BTreeMap<&'static str, Vec<&Merged>> = BTreeMap::new();
            for row in rows.iter().copied() {
                by_assembly
                    .entry(row.vep.map_or("", |vep| vep.assembly))
                    .or_default()
                    .push(row);
            }
            for (assembly, arows) in by_assembly {
                if assembly.is_empty() {
                    continue;
                }
                let head = arows[0];
                let coord = if assembly == "hg19" {
                    head.hg19
                } else {
                    head.hg38
                };
                if coord.is_none_or(|c| c.chrom.is_empty()) {
                    continue;
                }
                if !variant.assemblies.contains_key(assembly) {
                    let status = if assembly == "hg19" {
                        head.core.status
                    } else {
                        head.hg38_status
                    };
                    variant
                        .assemblies
                        .insert(assembly, coordinate_block(coord, status)?);
                }
                let Some(block) = variant.assemblies.get_mut(assembly) else {
                    continue;
                };
                for row in arows {
                    let Some(vep) = row.vep else { continue };
                    let last = |text: &str| text.rsplit(':').next().unwrap_or("").to_owned();
                    let hgvsc = last(&vep.hgvsc);
                    let hgvsp = last(&vep.hgvsp);
                    if hgvsc.is_empty() {
                        continue;
                    }
                    let view = block.transcripts.entry(vep.transcript.clone()).or_default();
                    view.by_hgvsc.insert(
                        hgvsc.clone(),
                        object(vec![
                            ("hgvsp", optional(&hgvsp)),
                            ("consequence", string(&vep.consequence)),
                            ("canonical", Json::Bool(vep.canonical == "YES")),
                            ("mane_select", optional(&vep.mane_select)),
                        ]),
                    );
                    if !hgvsp.is_empty() {
                        let coding = view.by_hgvsp.entry(hgvsp).or_default();
                        if !coding.contains(&hgvsc) {
                            coding.push(hgvsc);
                        }
                    }
                }
            }
        }
    }
    Ok(genes)
}

/// `date.today().isoformat()`.
fn today() -> String {
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |elapsed| elapsed.as_secs() as libc::time_t);
    // SAFETY: localtime_r fills the zeroed struct it is given.
    let local = unsafe {
        let mut local: libc::tm = std::mem::zeroed();
        libc::localtime_r(&now, &mut local);
        local
    };
    format!(
        "{:04}-{:02}-{:02}",
        local.tm_year + 1900,
        local.tm_mon + 1,
        local.tm_mday
    )
}

fn view_json(view: &View) -> Json {
    Json::Object(vec![
        (
            "by_hgvsc".to_owned(),
            Json::Object(
                view.by_hgvsc
                    .iter()
                    .map(|(key, value)| (key.clone(), value.clone()))
                    .collect(),
            ),
        ),
        (
            "by_hgvsp".to_owned(),
            Json::Object(
                view.by_hgvsp
                    .iter()
                    .map(|(key, values)| {
                        (
                            key.clone(),
                            Json::Array(values.iter().map(|value| string(value)).collect()),
                        )
                    })
                    .collect(),
            ),
        ),
    ])
}

fn genes_json(genes: &BTreeMap<String, Gene>) -> Json {
    Json::Object(
        genes
            .iter()
            .map(|(symbol, gene)| {
                let variants = gene.variants.iter().map(|(identifier, variant)| {
                    let assemblies = variant.assemblies.iter().map(|(name, block)| {
                        let transcripts = block
                            .transcripts
                            .iter()
                            .map(|(name, view)| (name.clone(), view_json(view)))
                            .collect();
                        (
                            (*name).to_owned(),
                            object(vec![
                                ("genomic", block.genomic.clone()),
                                ("transcripts", Json::Object(transcripts)),
                            ]),
                        )
                    });
                    (
                        identifier.clone(),
                        object(vec![
                            ("record", variant.record.clone()),
                            ("assemblies", Json::Object(assemblies.collect())),
                        ]),
                    )
                });
                (
                    symbol.clone(),
                    object(vec![
                        ("hgnc_id", gene.hgnc_id.clone()),
                        ("variants", Json::Object(variants.collect())),
                    ]),
                )
            })
            .collect(),
    )
}

/// `build_cache`: the join, the verdicts, the file, and the counts.
pub(super) fn build_cache(layers: Layers) -> Result<Vec<(String, Json)>> {
    let Layers {
        sources,
        core,
        vep,
        hg19_coords,
        hg38_coords,
        unmapped,
        reviews,
        summary,
        resolver,
        provenance,
        output,
    } = layers;
    if vep.iter().all(Vec::is_empty) {
        return Err(error("ValueError", "No objects to concatenate"));
    }
    let mut groups: HashMap<&str, Vec<&Source>> = HashMap::new();
    for source in sources {
        groups.entry(&source.allele_key).or_default().push(source);
    }
    let evidence: HashMap<&str, Vec<&Json>> = groups
        .iter()
        .map(|(key, rows)| (*key, rows.iter().map(|row| &row.evidence).collect()))
        .collect();
    let index = |coords: &[Coord]| {
        let mut map: HashMap<String, Vec<usize>> = HashMap::new();
        for (at, coord) in coords.iter().enumerate() {
            map.entry(coord.vcf_id.clone()).or_default().push(at);
        }
        map
    };
    let hg19_index = index(&hg19_coords);
    let hg38_index = index(&hg38_coords);
    let mut annotation: HashMap<&str, Vec<&VepRow>> = HashMap::new();
    for row in vep.iter().flatten() {
        annotation.entry(&row.vcf_id).or_default().push(row);
    }

    let mut merged: Vec<Merged> = Vec::new();
    for record in core {
        let symbol = groups
            .get(record.allele_key.as_str())
            .map_or_else(String::new, |rows| {
                resolver.resolve(&join_unique(rows.iter().map(|row| row.symbol.as_str())))
            });
        let pick = |map: &HashMap<String, Vec<usize>>| -> Vec<Option<usize>> {
            map.get(&record.vcf_id).map_or_else(
                || vec![None],
                |rows| rows.iter().copied().map(Some).collect(),
            )
        };
        for hg19 in pick(&hg19_index) {
            for hg38 in pick(&hg38_index) {
                let rows: Vec<Option<&VepRow>> =
                    annotation.get(record.vcf_id.as_str()).map_or_else(
                        || vec![None],
                        |rows| rows.iter().copied().map(Some).collect(),
                    );
                for row in rows {
                    if symbol.is_empty() {
                        continue;
                    }
                    merged.push(Merged {
                        core: record,
                        hg19: hg19.map(|at| &hg19_coords[at]),
                        hg38: hg38.map(|at| &hg38_coords[at]),
                        vep: row,
                        symbol: symbol.clone(),
                        vep_symbol: row.map_or_else(String::new, |r| resolver.resolve(&r.symbol)),
                        gene_match_status: "",
                        eligibility: String::new(),
                        reason: String::new(),
                        hg38_status: "",
                        accession: String::new(),
                    });
                }
            }
        }
    }

    let curated_present: HashSet<String> = merged
        .iter()
        .filter(|row| row.vep_symbol == row.symbol)
        .map(|row| row.core.vcf_id.clone())
        .collect();
    let any_symbol: HashSet<String> = merged
        .iter()
        .filter(|row| !row.vep_symbol.is_empty())
        .map(|row| row.core.vcf_id.clone())
        .collect();
    let before = merged.len();
    merged
        .retain(|row| !curated_present.contains(&row.core.vcf_id) || row.vep_symbol == row.symbol);
    log(&format!(
        "dropped {} transcript rows belonging to other genes at the same loci",
        before - merged.len()
    ));
    if merged.is_empty() {
        return Err(error("KeyError", "0"));
    }

    for row in &mut merged {
        row.gene_match_status = if row.vep_symbol.is_empty() {
            "source_gene_only"
        } else if row.vep_symbol == row.symbol {
            "gene_concordant"
        } else {
            "gene_discordant"
        };
        let id = &row.core.vcf_id;
        let (state, reason) = if curated_present.contains(id) || !any_symbol.contains(id) {
            (ELIGIBLE, "")
        } else if row.gene_match_status == "gene_discordant" {
            (QUARANTINE_GENE, "curated_gene_absent_from_locus")
        } else {
            (QUARANTINE_SIBLING, "sibling_of_gene_discordant_assertion")
        };
        row.eligibility = state.to_owned();
        row.reason = reason.to_owned();
        if state == ELIGIBLE {
            let key = (row.symbol.clone(), row.core.allele_key.clone());
            if let Some(review) = reviews
                .get(&key)
                .filter(|review| review.eligibility.to_lowercase() == "quarantined")
            {
                let mechanism = review.mechanism.to_uppercase();
                row.eligibility = reviewed_quarantine(&mechanism).to_owned();
                let reason = if review.reason.is_empty() {
                    "reviewed_non_gof"
                } else {
                    &review.reason
                };
                let named = if mechanism.is_empty() {
                    "unspecified"
                } else {
                    &mechanism
                };
                row.reason = format!("{reason}:{named}");
            }
        }
        let placed = row.hg38.is_some_and(|coord| !coord.chrom.is_empty());
        row.hg38_status = match (unmapped.contains(id), placed) {
            (true, true) => "endpoint_ref_match",
            (true, false) => "liftover_unmapped",
            (false, true) => "lifted_ref_match",
            (false, false) => "no_hg38_coordinate",
        };
        row.accession = match summary.get(&row.core.allele_key) {
            Some(data) if crate::gene_nonlof_mechanism::truthy(data) => {
                clean_json(get(data, "Accession")?)
            }
            _ => String::new(),
        };
    }
    let alleles: HashSet<&str> = merged
        .iter()
        .map(|row| row.core.allele_key.as_str())
        .collect();
    log(&format!(
        "retaining {} transcript annotations across {} source records",
        merged.len(),
        alleles.len()
    ));

    let genes = nest(&merged, &evidence, resolver, unmapped)?;
    let mut metadata = vec![
        ("source".to_owned(), string("GoFCards")),
        ("mechanism".to_owned(), string("GOF")),
        ("derived_on".to_owned(), string(&today())),
        (
            "structure".to_owned(),
            string(
                "genes -> alleles (normalized GRCh37 chrom|pos|ref|alt) -> assemblies -> \
transcripts",
            ),
        ),
    ];
    metadata.extend(provenance.iter().cloned());
    let cache = object(vec![
        ("metadata", Json::Object(metadata)),
        ("genes", genes_json(&genes)),
    ]);
    if let Some(parent) = output.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent).map_err(|e| os_error(e, parent, None))?;
    }
    let body = json_text::dumps_with_options(
        &cache,
        DumpsOptions {
            indent: None,
            sort_keys: true,
            ensure_ascii: false,
            compact: true,
            allow_nan: true,
        },
    )
    .map_err(|e| error("ValueError", e))?;
    fs::write(output, gzip_open_bytes(body.as_bytes(), output)?)
        .map_err(|e| os_error(e, output, None))?;

    let variants: Vec<&Variant> = genes.values().flat_map(|g| g.variants.values()).collect();
    let views: Vec<&View> = variants
        .iter()
        .flat_map(|v| v.assemblies.values())
        .flat_map(|block| block.transcripts.values())
        .collect();
    let eligible = variants.iter().filter(|v| v.status == ELIGIBLE).count();
    let stats: Vec<(&str, usize)> = vec![
        ("genes", genes.len()),
        ("variants", variants.len()),
        ("eligible_variants", eligible),
        ("quarantined_variants", variants.len() - eligible),
        (
            "evidence_entries",
            variants.iter().map(|v| v.evidence).sum(),
        ),
        (
            "variants_on_hg19",
            variants
                .iter()
                .filter(|v| v.assemblies.contains_key("hg19"))
                .count(),
        ),
        (
            "variants_on_hg38",
            variants
                .iter()
                .filter(|v| v.assemblies.contains_key("hg38"))
                .count(),
        ),
        ("transcript_views", views.len()),
        ("hgvsc_keys", views.iter().map(|v| v.by_hgvsc.len()).sum()),
        ("hgvsp_keys", views.iter().map(|v| v.by_hgvsp.len()).sum()),
        ("liftover_unmapped", unmapped.len()),
    ];
    log(&format!(
        "cache written: {}",
        stats
            .iter()
            .map(|(key, value)| format!("{key}={value}"))
            .collect::<Vec<_>>()
            .join("; ")
    ));
    Ok(stats
        .into_iter()
        .map(|(key, value)| (key.to_owned(), Json::Int(value as i128)))
        .collect())
}

/// `gzip.open(path, "wt")` writing `body`: the file name and the time go in
/// the header, and the text wrapper flushes the stream when it closes.
fn gzip_open_bytes(body: &[u8], path: &Path) -> Result<Vec<u8>> {
    let name: String = path
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_default();
    let latin: Option<Vec<u8>> = name.chars().map(|c| u8::try_from(c as u32).ok()).collect();
    let mut stored = latin.unwrap_or_default();
    if stored.ends_with(b".gz") {
        stored.truncate(stored.len() - 3);
    }
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |elapsed| elapsed.as_secs() as u32);
    Ok(crate::python_gzip::compress_named(
        body, 9, true, &stored, now,
    )?)
}
