//! The review workbook of a family: the rows of its top-ranked variants from
//! the final table of the ACMG step, as one Excel sheet in the order clinical
//! reviewers read it. A few cells are condensed from several columns; every
//! other cell is the table's value.
//!
//! - Rows: the table is in rank order, so the rows are read until the rank
//!   passes `top_variants`; every transcript row of those variants is kept.
//! - The first fifteen columns stay in view when the sheet is scrolled:
//!   variant, gene, transcript, consequence, the ACMG result, the patient's
//!   HPO terms the gene matches, the sort index and the rank.
//! - Condensed cells: the highest gnomAD frequency with its allele number and
//!   homozygotes; the strongest SpliceAI score and its effect; where the
//!   inheritance comes from, in words (for an inheritance read from the
//!   gene's constraint, its ClinGen haploinsufficiency score and LOEUF); the
//!   BS1 frequency cutoffs as fractions; the protein domains without the
//!   structure mappings.
//!
//! A column the table lacks gives empty cells, so tables of runs without a
//! pedigree, a gene ranking or PEXT values still give a workbook. The file is
//! written directly as the XML parts of an `.xlsx` archive.
use crate::acmg_step::gene_mechanism_common::{
    DEFAULT_CLINGEN_DOSAGE, DEFAULT_HGNC_TABLE, DEFAULT_LOEUF_TABLE,
};
use crate::acmg_step::gene_mechanism_resources::{GeneMechanismHub, HubPaths};
use crate::acmg_step::types::SymbolResolver;
use crate::cli::{Result, invalid};
use crate::pandas_table::Cell as TableCell;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

const USAGE: &str = "usage: priva review-workbook --table RANKED.tsv --output WORKBOOK.xlsx [--top-variants 1000] [--ped COHORT.ped --family ID] [--gene-ranking RANKED_GENES.tsv] [--expected-incidence 0.0001] [--pext-tissues LIST] [--data-dir DATA]

The rows of the top-ranked variants of a family's final table (the ACMG step's output, in rank order) as one Excel sheet in the reviewers' column order.";

/// How many variants a workbook holds unless told otherwise.
pub(crate) const DEFAULT_TOP_VARIANTS: usize = 1000;
/// The columns that stay in view when the sheet is scrolled sideways.
const FROZEN_COLUMNS: usize = 15;
/// A cell's longest text (Excel's limit is 32,767 characters).
const LONGEST_TEXT: usize = 32_000;
/// The longest list of protein domains shown.
const LONGEST_DOMAINS: usize = 120;
/// Values that stand for "no value" in the table.
const EMPTY: [&str; 5] = ["", "nan", "NaN", "NA", "."];
/// The header fills, one per column group in the order the groups appear.
const FILLS: [&str; 14] = [
    "37474F", "00695C", "4E342E", "283593", "6A1B9A", "AD1457", "1565C0", "2E7D32", "EF6C00",
    "5D4037", "455A64", "00838F", "C62828", "424242",
];
/// The SpliceAI effects, in the order of the table's score columns.
const EFFECTS: [(&str, &str); 4] = [
    ("AG", "acceptor gain"),
    ("AL", "acceptor loss"),
    ("DG", "donor gain"),
    ("DL", "donor loss"),
];

/// What a column of the workbook shows.
#[derive(Clone, Copy)]
enum Kind {
    /// The table's column of this name, unchanged.
    Copy(&'static str),
    Canonical,
    EscapesNmd,
    HpoTerms,
    PextValue,
    PextTissue,
    /// One column per sequenced member.
    Genotypes,
    Depths,
    /// One column per patient (affected member with a phase column).
    Phases,
    InheritanceSource,
    GnomadFrequency,
    GnomadAlleleNumber,
    GnomadHomozygotes,
    IncidenceCutoff,
    GeneCutoff,
    Domains,
    SpliceScore,
    SpliceEffect,
}

/// The sheet's columns in order: group, label, content.
const LAYOUT: &[(&str, &str, Kind)] = &[
    ("Variant and transcript", "Chr", Kind::Copy("chrom")),
    ("Variant and transcript", "Pos", Kind::Copy("pos")),
    ("Variant and transcript", "Ref", Kind::Copy("ref")),
    ("Variant and transcript", "Alt", Kind::Copy("alt")),
    (
        "Variant and transcript",
        "Variant ID",
        Kind::Copy("variant_id"),
    ),
    (
        "Variant and transcript",
        "Gene symbol",
        Kind::Copy("SYMBOL"),
    ),
    (
        "Variant and transcript",
        "Transcript ID",
        Kind::Copy("Feature"),
    ),
    (
        "Variant and transcript",
        "Molecular consequence",
        Kind::Copy("Consequence"),
    ),
    ("ACMG result", "ACMG criteria", Kind::Copy("ACMG_criteria")),
    ("ACMG result", "ACMG class", Kind::Copy("ACMG_class")),
    (
        "ACMG result",
        "ACMG quantified score",
        Kind::Copy("ACMG_quant_score"),
    ),
    (
        "Phenotype match",
        "Patient HPO terms the gene matches",
        Kind::HpoTerms,
    ),
    (
        "Phenotype match",
        "Phenotype gene index (1 to 1.2)",
        Kind::Copy("phenotype_gene_index"),
    ),
    ("Ranking", "Sort index", Kind::Copy("sort_index")),
    ("Ranking", "Rank", Kind::Copy("variant_rank")),
    (
        "Variant and transcript",
        "Canonical transcript",
        Kind::Canonical,
    ),
    ("Variant and transcript", "HGVSc", Kind::Copy("HGVSc")),
    ("Variant and transcript", "HGVSp", Kind::Copy("HGVSp")),
    (
        "Variant and transcript",
        "Escapes nonsense-mediated decay",
        Kind::EscapesNmd,
    ),
    ("Phenotype match", "", Kind::PextValue),
    ("Phenotype match", "", Kind::PextTissue),
    ("Family genotypes", "GT", Kind::Genotypes),
    ("Family genotypes", "AD", Kind::Depths),
    (
        "Inheritance",
        "Phase with the gene's other pathogenic variants",
        Kind::Phases,
    ),
    (
        "Inheritance",
        "Inheritance used for the variant",
        Kind::Copy("variant_inheritance"),
    ),
    (
        "Inheritance",
        "Where that inheritance comes from",
        Kind::InheritanceSource,
    ),
    ("Inheritance", "X-linked", Kind::Copy("variant_x_linked")),
    (
        "Inheritance",
        "Penetrance",
        Kind::Copy("variant_penetrance"),
    ),
    (
        "Inheritance",
        "Gene's inheritance modes in HPO",
        Kind::Copy("HPO_gene_inheritance"),
    ),
    ("ClinVar", "ClinVar significance", Kind::Copy("CLNSIG")),
    ("ClinVar", "ClinVar review status", Kind::Copy("CLNREVSTAT")),
    (
        "Population frequency",
        "gnomAD highest AF",
        Kind::GnomadFrequency,
    ),
    (
        "Population frequency",
        "gnomAD AN at the highest AF",
        Kind::GnomadAlleleNumber,
    ),
    (
        "Population frequency",
        "gnomAD homozygotes at the highest AF",
        Kind::GnomadHomozygotes,
    ),
    (
        "Population frequency",
        "Control cohort AC",
        Kind::Copy("control_AC"),
    ),
    (
        "Population frequency",
        "Control cohort AN",
        Kind::Copy("control_AN"),
    ),
    (
        "Population frequency",
        "Control cohort AF",
        Kind::Copy("control_AF"),
    ),
    (
        "Population frequency",
        "Control cohort homozygotes",
        Kind::Copy("control_nhomalt"),
    ),
    (
        "Population frequency",
        "BS1 disease incidence cutoff",
        Kind::IncidenceCutoff,
    ),
    (
        "Population frequency",
        "BS1 cutoff from the gene's pathogenic variants",
        Kind::GeneCutoff,
    ),
    ("Protein and predictors", "Protein domains", Kind::Domains),
    (
        "Protein and predictors",
        "AlphaMissense pathogenicity",
        Kind::Copy("am_pathogenicity"),
    ),
    (
        "Protein and predictors",
        "AlphaMissense class",
        Kind::Copy("am_class"),
    ),
    (
        "Protein and predictors",
        "PrimateAI",
        Kind::Copy("PrimateAI"),
    ),
    (
        "Protein and predictors",
        "CADD phred",
        Kind::Copy("CADD_phred"),
    ),
    ("Splicing", "SpliceAI2 highest score", Kind::SpliceScore),
    ("Splicing", "SpliceAI2 predicted effect", Kind::SpliceEffect),
    (
        "Splicing",
        "Splicing: loss of function",
        Kind::Copy("splicing_lof"),
    ),
    (
        "Splicing",
        "Splicing: protein length change",
        Kind::Copy("splicing_len_changing"),
    ),
    (
        "Splicing",
        "Splicing: frameshift",
        Kind::Copy("splicing_frameshift"),
    ),
    (
        "Splicing",
        "Splicing: hits an intolerant domain",
        Kind::Copy("splicing_span_intol_domain"),
    ),
    (
        "Splicing",
        "Splicing: affected exons",
        Kind::Copy("splicing_affected_exons"),
    ),
    (
        "Splicing",
        "Splicing: removes 10% of the protein or more",
        Kind::Copy("splicing_ten_percent_protein"),
    ),
    (
        "5' UTR",
        "5' UTR consequence",
        Kind::Copy("5UTR_consequence"),
    ),
    ("5' UTR", "5' UTR: loss of function", Kind::Copy("5UTR_lof")),
    (
        "5' UTR",
        "5' UTR: protein length change",
        Kind::Copy("5UTR_len_changing"),
    ),
    (
        "5' UTR",
        "5' UTR: frameshift",
        Kind::Copy("5UTR_frameshift"),
    ),
    (
        "5' UTR",
        "5' UTR: hits an intolerant domain",
        Kind::Copy("5UTR_span_intol_domain"),
    ),
    (
        "Functional assay",
        "MaveDB supports PS3",
        Kind::Copy("MaveDB_PS3"),
    ),
    (
        "Functional assay",
        "MaveDB supports BS3",
        Kind::Copy("MaveDB_BS3"),
    ),
    (
        "Functional assay",
        "MaveDB score interpretation",
        Kind::Copy("MaveDB_score_interpretation"),
    ),
    (
        "Known pathogenic variants",
        "Gene: highest AF of ClinVar pathogenic variants",
        Kind::Copy("clinvar_patho_gene_max_af"),
    ),
    (
        "Known pathogenic variants",
        "Gene: highest homozygote frequency of ClinVar pathogenic variants",
        Kind::Copy("clinvar_patho_gene_max_hom_freq"),
    ),
    (
        "Known pathogenic variants",
        "Exon: median AF of ClinVar pathogenic variants",
        Kind::Copy("exon_patho_median_af"),
    ),
    (
        "Known pathogenic variants",
        "Exon: highest AF of ClinVar pathogenic variants",
        Kind::Copy("exon_patho_max_af"),
    ),
    (
        "Exact mechanism match",
        "Exact mechanism of this variant",
        Kind::Copy("variant_exact_mechanisms"),
    ),
    (
        "Exact mechanism match",
        "Evidence of the exact mechanism",
        Kind::Copy("variant_exact_mechanism_evidence"),
    ),
    (
        "Exact mechanism match",
        "How the mechanism record matched",
        Kind::Copy("variant_mechanism_match_route"),
    ),
    (
        "Exact mechanism match",
        "Gain-of-function tag",
        Kind::Copy("variant_gof_tag"),
    ),
    (
        "Exact mechanism match",
        "GoFCards accession",
        Kind::Copy("gofcards_accession_id"),
    ),
    (
        "Exact mechanism match",
        "GoFCards variant ID",
        Kind::Copy("gofcards_variant_id"),
    ),
    (
        "Exact mechanism match",
        "How the GoFCards record matched",
        Kind::Copy("gofcards_match_route"),
    ),
    (
        "Ranking",
        "Factor: common in cohort controls",
        Kind::Copy("control_common_index"),
    ),
    (
        "Ranking",
        "Factor: dispensable gene",
        Kind::Copy("dispensable_gene_index"),
    ),
    (
        "Ranking",
        "Factor: genotype fits inheritance and mechanism",
        Kind::Copy("zygosity_inheritance_mechanism_compatibility"),
    ),
    ("Ranking", "Factor: PEXT", Kind::Copy("pext_sort_index")),
];

/// The widths of the columns that are not 13 characters wide.
const WIDTHS: &[(&str, f64)] = &[
    ("Chr", 5.0),
    ("Pos", 10.0),
    ("Ref", 5.0),
    ("Alt", 5.0),
    ("Variant ID", 18.0),
    ("Gene symbol", 10.0),
    ("Transcript ID", 16.0),
    ("Molecular consequence", 20.0),
    ("ACMG criteria", 18.0),
    ("ACMG class", 14.0),
    ("ACMG quantified score", 7.0),
    ("Patient HPO terms the gene matches", 22.0),
    ("Phenotype gene index (1 to 1.2)", 7.0),
    ("Sort index", 7.0),
    ("Rank", 6.0),
    ("HGVSc", 30.0),
    ("HGVSp", 28.0),
    ("Protein domains", 34.0),
    ("ClinVar significance", 22.0),
    ("ClinVar review status", 24.0),
    ("SpliceAI2 predicted effect", 32.0),
    ("Where that inheritance comes from", 44.0),
    ("BS1 disease incidence cutoff", 16.0),
    ("BS1 cutoff from the gene's pathogenic variants", 16.0),
    ("Gene's inheritance modes in HPO", 26.0),
];

/// A cell: nothing, a text, or a number as it is written into the sheet.
#[derive(Clone, Debug, PartialEq)]
enum Value {
    Empty,
    Text(String),
    Number(String),
}

/// What a workbook is made from.
pub(crate) struct Inputs<'a> {
    /// The final table of the ACMG step, in rank order.
    pub(crate) table: &'a Path,
    pub(crate) output: &'a Path,
    pub(crate) top_variants: usize,
    /// The pedigree and the family in it, for the members' roles.
    pub(crate) ped: Option<&'a Path>,
    pub(crate) family: Option<&'a str>,
    /// The family's ranked gene list (`priva hpo-gene-ranking rank-genes`).
    pub(crate) gene_ranking: Option<&'a Path>,
    /// The run's expected disease incidence (`exp_disease_incidence`).
    pub(crate) expected_incidence: f64,
    /// The run's PEXT tissues (`pext_tissues`), comma-separated.
    pub(crate) pext_tissues: &'a str,
    /// The repository's data folder (ClinGen dosage, LOEUF and HGNC tables).
    pub(crate) data_dir: &'a Path,
}

/// What was written.
#[derive(Debug, PartialEq)]
pub(crate) struct Summary {
    pub(crate) rows: usize,
    pub(crate) variants: usize,
    pub(crate) columns: usize,
}

fn is_empty(value: &str) -> bool {
    EMPTY.contains(&value)
}

/// The number a cell holds, if it holds one.
fn number(value: &str) -> Option<f64> {
    if is_empty(value) {
        return None;
    }
    value.parse::<f64>().ok()
}

fn all_digits(text: &str) -> bool {
    !text.is_empty() && text.bytes().all(|byte| byte.is_ascii_digit())
}

/// A number as the sheet holds it.
fn show_number(value: f64) -> String {
    let size = value.abs();
    if value == 0.0 {
        "0".to_owned()
    } else if (1e-4..1e15).contains(&size) {
        format!("{value}")
    } else {
        format!("{value:e}")
    }
}

/// A table value as a cell: an integer or a decimal number becomes a number,
/// everything else stays text.
fn tidy(value: &str) -> Value {
    if is_empty(value) {
        return Value::Empty;
    }
    let digits = value.strip_prefix('-').unwrap_or(value);
    if all_digits(digits) && digits.len() <= 15 && (digits == "0" || !digits.starts_with('0')) {
        return Value::Number(value.to_owned());
    }
    if let Some(whole) = value.strip_suffix(".0") {
        let whole_digits = whole.strip_prefix('-').unwrap_or(whole);
        if all_digits(whole_digits) && value.len() < 12 {
            return Value::Number(whole.to_owned());
        }
    }
    let (mantissa, exponent) = match digits.find(['e', 'E']) {
        Some(at) => (&digits[..at], Some(&digits[at + 1..])),
        None => (digits, None),
    };
    let mantissa_fits = match mantissa.split_once('.') {
        Some((whole, fraction)) => (whole.is_empty() || all_digits(whole)) && all_digits(fraction),
        None => all_digits(mantissa),
    };
    let exponent_fits = exponent
        .is_none_or(|exponent| all_digits(exponent.strip_prefix(['-', '+']).unwrap_or(exponent)));
    if mantissa_fits && exponent_fits {
        if let Ok(parsed) = value.parse::<f64>() {
            if parsed.is_finite() {
                return Value::Number(show_number(parsed));
            }
        }
    }
    Value::Text(value.chars().take(LONGEST_TEXT).collect())
}

/// A frequency as "1/N" with thousands separators; empty for none.
fn fraction(value: Option<f64>) -> String {
    let Some(value) = value.filter(|value| *value > 0.0) else {
        return String::new();
    };
    if value >= 1.0 {
        return "1/1".to_owned();
    }
    let digits = format!("{}", (1.0 / value).round_ties_even() as u64);
    let mut grouped = String::new();
    for (place, digit) in digits.chars().enumerate() {
        if place > 0 && (digits.len() - place) % 3 == 0 {
            grouped.push(',');
        }
        grouped.push(digit);
    }
    format!("1/{grouped}")
}

/// The protein domains without the structure mappings, each once, cut short.
fn domains(value: &str) -> String {
    if is_empty(value) {
        return String::new();
    }
    let mut kept: Vec<&str> = Vec::new();
    for entry in value.split('&') {
        let structure =
            entry.starts_with("PDB-ENSP_mappings") || entry.starts_with("AFDB-ENSP_mappings");
        if !structure && !kept.contains(&entry) {
            kept.push(entry);
        }
    }
    let text = kept.join("; ");
    if text.chars().count() <= LONGEST_DOMAINS {
        return text;
    }
    let cut: String = text.chars().take(LONGEST_DOMAINS - 3).collect();
    format!("{}...", cut.trim_end_matches([';', ' ']))
}

/// The words for the source of a variant's inheritance other than the gene's
/// constraint.
fn inheritance_source(code: &str) -> &str {
    match code {
        "matched_history" => "disease of this gene with the same mechanism as this variant",
        "unresolved_condition_history" => "disease of this gene, mechanism not established",
        "matched_and_unresolved_history" => {
            "diseases of this gene: some with the same mechanism as this variant, some with no established mechanism"
        }
        "mechanism_mismatch" => "no disease of this gene fits this variant's effect",
        other => other,
    }
}

/// A member of the pedigree.
struct Member {
    father: String,
    mother: String,
    affected: bool,
}

/// The member's role in the family, from the proband's parents in the pedigree.
/// `patients` are the affected members, the proband first.
fn role(sample: &str, patients: &[String], members: &HashMap<String, Member>) -> &'static str {
    let Some(proband) = patients.first() else {
        return "member";
    };
    if sample == proband {
        return "proband";
    }
    let Some(parents) = members.get(proband) else {
        return "member";
    };
    if sample == parents.father {
        return "father";
    }
    if sample == parents.mother {
        return "mother";
    }
    match members.get(sample) {
        Some(own)
            if own.father == parents.father
                && own.mother == parents.mother
                && parents.father != "0" =>
        {
            if own.affected || patients.iter().any(|patient| patient == sample) {
                "affected sibling"
            } else {
                "sibling"
            }
        }
        _ => "member",
    }
}

/// The family's members: sample, parents and whether affected (phenotype 2).
fn read_ped(path: &Path, family: Option<&str>) -> Result<HashMap<String, Member>> {
    let mut members = HashMap::new();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let fields: Vec<&str> = line.split_whitespace().collect();
        if line.starts_with('#') || fields.len() < 4 {
            continue;
        }
        if family.is_none_or(|family| family == fields[0]) {
            members.insert(
                fields[1].to_owned(),
                Member {
                    father: fields[2].to_owned(),
                    mother: fields[3].to_owned(),
                    affected: fields.get(5) == Some(&"2"),
                },
            );
        }
    }
    Ok(members)
}

/// The patient's HPO terms each ranked gene matches, by gene symbol.
fn read_gene_ranking(path: &Path) -> Result<HashMap<String, String>> {
    let mut lines = BufReader::new(File::open(path)?).lines();
    let Some(header) = lines.next().transpose()? else {
        return Ok(HashMap::new());
    };
    let names: Vec<&str> = header.split('\t').collect();
    let place = |name: &str| names.iter().position(|field| *field == name);
    let (Some(symbol), Some(matched)) =
        (place("gene_symbol"), place("candidate_generating_hp_names"))
    else {
        return Ok(HashMap::new());
    };
    let supporting = place("support_only_hp_names");
    let mut terms = HashMap::new();
    for line in lines {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        let field = |place: usize| fields.get(place).map_or("", |text| text.trim());
        let mut text = field(matched).to_owned();
        let more = supporting.map_or("", field);
        if !more.is_empty() {
            text.push_str("; supporting: ");
            text.push_str(more);
        }
        terms.insert(field(symbol).to_owned(), text);
    }
    Ok(terms)
}

/// The gene resources behind an inheritance read from the gene's constraint.
struct Constraint {
    hub: GeneMechanismHub,
    shown: HashMap<String, String>,
}

impl Constraint {
    fn new(data_dir: &Path) -> Result<Self> {
        let path = |relative: &str| data_dir.join(relative).to_string_lossy().into_owned();
        let hub = GeneMechanismHub::new(HubPaths {
            condition_cache: String::new(),
            clingen_dosage: path(DEFAULT_CLINGEN_DOSAGE),
            loeuf_table: path(DEFAULT_LOEUF_TABLE),
            hgnc_table: path(DEFAULT_HGNC_TABLE),
        })?;
        Ok(Self {
            hub,
            shown: HashMap::new(),
        })
    }

    /// The gene's ClinGen haploinsufficiency score and LOEUF, in words.
    fn of(&mut self, gene: &str) -> Result<String> {
        if let Some(text) = self.shown.get(gene) {
            return Ok(text.clone());
        }
        let symbol = self.hub.try_resolve_symbol(TableCell::Text(gene));
        let clingen = self.hub.load_clingen()?;
        let loeuf = self.hub.load_loeuf()?;
        let dosage = match clingen.get(&symbol) {
            Some(record) => match record.haploinsufficiency_score {
                Some(score) => format!(
                    "ClinGen haploinsufficiency score {score} ({})",
                    record.haploinsufficiency_description.to_lowercase()
                ),
                None => "ClinGen haploinsufficiency not scored".to_owned(),
            },
            None => "no ClinGen haploinsufficiency record".to_owned(),
        };
        let constraint = match loeuf.get(&symbol).copied().filter(|value| !value.is_nan()) {
            Some(value) => format!("LOEUF {value:.2}"),
            None => "no LOEUF".to_owned(),
        };
        let text = format!("gene constraint: {dosage}; {constraint}");
        self.shown.insert(gene.to_owned(), text.clone());
        Ok(text)
    }
}

/// One column of the sheet.
struct Column {
    group: &'static str,
    label: String,
    kind: Kind,
    /// The member a genotype, depth or phase column belongs to.
    sample: Option<String>,
}

/// The table's header and what the cells are made from.
struct Sheet {
    places: HashMap<String, usize>,
    columns: Vec<Column>,
    terms: HashMap<String, String>,
    constraint: Constraint,
    incidence: String,
    /// The PEXT column shown and its tissue, when the table has none of its own.
    pext: Option<(String, String)>,
}

impl Sheet {
    fn cell<'a>(&self, row: &[&'a str], name: &str) -> &'a str {
        self.places
            .get(name)
            .and_then(|place| row.get(*place))
            .copied()
            .unwrap_or("")
    }

    /// The highest gnomAD frequency with its allele number and homozygotes.
    fn gnomad(&self, row: &[&str]) -> [Value; 3] {
        let whole = number(self.cell(row, "gnomAD_joint_AF"));
        let highest = number(self.cell(row, "gnomAD_joint_AF_max"));
        match (whole, highest) {
            (None, None) => [Value::Empty, Value::Empty, Value::Empty],
            (whole, Some(highest)) if whole.is_none_or(|whole| highest >= whole) => [
                Value::Number(show_number(highest)),
                tidy(self.cell(row, "gnomAD_joint_AN_max")),
                tidy(self.cell(row, "gnomAD_nhomalt_max")),
            ],
            (whole, _) => {
                let counts = [
                    number(self.cell(row, "gnomAD_nhomalt_XX")),
                    number(self.cell(row, "gnomAD_nhomalt_XY")),
                ];
                let homozygotes = if counts.iter().any(Option::is_some) {
                    Value::Number(format!("{}", counts.iter().flatten().sum::<f64>() as i64))
                } else {
                    Value::Empty
                };
                [
                    Value::Number(show_number(whole.unwrap_or(f64::NAN))),
                    tidy(self.cell(row, "gnomAD_joint_AN")),
                    homozygotes,
                ]
            }
        }
    }

    /// The strongest SpliceAI score and its effect in words.
    fn splice(&self, row: &[&str]) -> (Value, String) {
        let mut best: Option<(f64, &str, &str)> = None;
        for (code, effect) in EFFECTS {
            let score = number(self.cell(row, &format!("SpliceAI_pred_DS_{code}")));
            if let Some(score) = score {
                if best.is_none_or(|(top, _, _)| score > top) {
                    best = Some((
                        score,
                        effect,
                        self.cell(row, &format!("SpliceAI_pred_DP_{code}")),
                    ));
                }
            }
        }
        match best {
            None => (Value::Empty, "no score".to_owned()),
            Some((score, _, _)) if score < 0.1 => (
                Value::Number(show_number(score)),
                "none (below 0.1)".to_owned(),
            ),
            Some((score, effect, distance)) => {
                let place = match number(distance) {
                    Some(distance) => format!(", {:+} bases from the variant", distance as i64),
                    None => String::new(),
                };
                (
                    Value::Number(show_number(score)),
                    format!("{effect}{place}"),
                )
            }
        }
    }

    fn value(&mut self, column: usize, row: &[&str]) -> Result<Value> {
        let kind = self.columns[column].kind;
        let text = |text: String| tidy(&text);
        Ok(match kind {
            Kind::Copy(name) => tidy(self.cell(row, name)),
            Kind::Canonical => {
                let canonical = self.cell(row, "CANONICAL") == "YES";
                Value::Text(if canonical { "yes" } else { "no" }.to_owned())
            }
            Kind::EscapesNmd => {
                if self.cell(row, "NMD").contains("NMD_escaping") {
                    Value::Text("yes".to_owned())
                } else {
                    Value::Empty
                }
            }
            Kind::HpoTerms => text(
                self.terms
                    .get(self.cell(row, "SYMBOL"))
                    .cloned()
                    .unwrap_or_default(),
            ),
            Kind::PextValue => match &self.pext {
                Some((column, _)) => tidy(self.cell(row, column)),
                None => tidy(self.cell(row, "PEXT_FamilyTissues")),
            },
            Kind::PextTissue => match &self.pext {
                Some((column, tissue)) if !is_empty(self.cell(row, column)) => {
                    Value::Text(tissue.clone())
                }
                Some(_) => Value::Empty,
                None => tidy(self.cell(row, "PEXT_FamilyTissues_top")),
            },
            Kind::Genotypes => {
                let sample = self.columns[column].sample.as_deref().unwrap_or("");
                tidy(self.cell(row, sample))
            }
            Kind::Depths => {
                let sample = self.columns[column].sample.as_deref().unwrap_or("");
                tidy(self.cell(row, &format!("{sample}_AD")))
            }
            Kind::Phases => match self.columns[column].sample.as_deref() {
                Some(patient) => tidy(self.cell(row, &format!("{patient}_phase_patho"))),
                None => Value::Empty,
            },
            Kind::InheritanceSource => {
                let code = self.cell(row, "variant_inheritance_basis");
                if code == "gene_constraint" {
                    Value::Text(self.constraint.of(self.cell(row, "SYMBOL"))?)
                } else if is_empty(code) {
                    Value::Empty
                } else {
                    Value::Text(inheritance_source(code).to_owned())
                }
            }
            Kind::GnomadFrequency => self.gnomad(row)[0].clone(),
            Kind::GnomadAlleleNumber => self.gnomad(row)[1].clone(),
            Kind::GnomadHomozygotes => self.gnomad(row)[2].clone(),
            // The cutoffs only: which of them a BS1 came from depends on the
            // inheritance (the dominant tests use the gene's own cutoff when it
            // has one), so the cells do not say where BS1 was applied.
            Kind::IncidenceCutoff => Value::Text(self.incidence.clone()),
            Kind::GeneCutoff => {
                let frequency = number(self.cell(row, "clinvar_patho_gene_max_af"));
                if frequency.is_some_and(|frequency| frequency > 0.0) {
                    Value::Text(fraction(frequency))
                } else {
                    Value::Empty
                }
            }
            Kind::Domains => text(domains(self.cell(row, "DOMAINS"))),
            Kind::SpliceScore => self.splice(row).0,
            Kind::SpliceEffect => Value::Text(self.splice(row).1),
        })
    }
}

/// The sheet's columns for a table with this header.
fn plan(header: &[&str], inputs: &Inputs<'_>) -> Result<Sheet> {
    let places: HashMap<String, usize> = header
        .iter()
        .enumerate()
        .map(|(place, name)| ((*name).to_owned(), place))
        .collect();
    let members = match inputs.ped {
        Some(ped) => read_ped(ped, inputs.family)?,
        None => HashMap::new(),
    };
    // The members: the pedigree's samples the table has, else the columns
    // that come with allele depths.
    let samples: Vec<&str> = header
        .iter()
        .copied()
        .filter(|name| {
            if members.is_empty() {
                places.contains_key(&format!("{name}_AD"))
            } else {
                members.contains_key(*name)
            }
        })
        .collect();
    // The patients: the members with a phase column, the first one the proband;
    // without such a column, the pedigree's affected members.
    let phased: Vec<String> = header
        .iter()
        .filter_map(|name| name.strip_suffix("_phase_patho"))
        .map(str::to_owned)
        .collect();
    let patients: Vec<String> = if phased.is_empty() {
        samples
            .iter()
            .filter(|sample| members.get(**sample).is_some_and(|member| member.affected))
            .map(|sample| (*sample).to_owned())
            .collect()
    } else {
        phased.clone()
    };
    // The PEXT value the ranking used: the first configured tissue, else the
    // mean over all tissues. A table with the family tissues' own column keeps it.
    let family_tissues = places.contains_key("PEXT_FamilyTissues");
    let pext = if family_tissues {
        None
    } else {
        let first = inputs.pext_tissues.split(',').next().unwrap_or("").trim();
        let column = format!("PEXT_{}", first.replace('-', "_"));
        if !first.is_empty() && first != "null" && places.contains_key(&column) {
            Some((column, first.to_owned()))
        } else {
            Some(("PEXT_MEAN".to_owned(), "mean over all tissues".to_owned()))
        }
    };
    let mut columns = Vec::new();
    for (group, label, kind) in LAYOUT {
        let mut push = |label: String, sample: Option<String>| {
            columns.push(Column {
                group,
                label,
                kind: *kind,
                sample,
            })
        };
        match kind {
            Kind::Genotypes | Kind::Depths => {
                for sample in &samples {
                    let role = role(sample, &patients, &members);
                    push(
                        format!("{label} {sample} ({role})"),
                        Some((*sample).to_owned()),
                    );
                }
            }
            // One phase column per patient; a single one needs no name.
            Kind::Phases => match phased.as_slice() {
                [] => push((*label).to_owned(), None),
                [one] => push((*label).to_owned(), Some(one.clone())),
                several => {
                    for patient in several {
                        push(format!("{label} ({patient})"), Some(patient.clone()));
                    }
                }
            },
            Kind::PextValue if family_tissues => {
                push("PEXT in the phenotype's tissues".to_owned(), None)
            }
            Kind::PextValue => push("PEXT in the relevant tissue".to_owned(), None),
            Kind::PextTissue if family_tissues => {
                push("Phenotype tissue of that PEXT value".to_owned(), None)
            }
            Kind::PextTissue => push("Tissue of that PEXT value".to_owned(), None),
            _ => push((*label).to_owned(), None),
        }
    }
    let terms = match inputs.gene_ranking {
        Some(path) => read_gene_ranking(path)?,
        None => HashMap::new(),
    };
    Ok(Sheet {
        places,
        columns,
        terms,
        constraint: Constraint::new(inputs.data_dir)?,
        incidence: fraction(Some(inputs.expected_incidence)),
        pext,
    })
}

/// The letters of a column (0 is A).
fn letters(mut column: usize) -> String {
    let mut name = Vec::new();
    loop {
        name.push(b'A' + (column % 26) as u8);
        if column < 26 {
            break;
        }
        column = column / 26 - 1;
    }
    name.reverse();
    String::from_utf8(name).expect("letters")
}

/// A text inside the sheet's XML: escaped, without the characters XML refuses.
fn escaped(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for letter in text.chars() {
        match letter {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\t' | '\n' | '\r' => out.push(letter),
            control if (control as u32) < 0x20 => {}
            '\u{FFFE}' | '\u{FFFF}' => {}
            other => out.push(other),
        }
    }
    out
}

const CONTENT_TYPES: &str = r#"<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/><Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/></Types>"#;
const ROOT_RELATIONS: &str = r#"<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>"#;
const WORKBOOK: &str = r#"<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Review_table" sheetId="1" r:id="rId1"/></sheets></workbook>"#;
const WORKBOOK_RELATIONS: &str = r#"<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/><Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>"#;

/// The styles: Arial 9 for every cell; a bold white header on each group's fill.
fn styles() -> String {
    let mut fills = String::new();
    let mut formats = String::new();
    for (place, fill) in FILLS.iter().enumerate() {
        fills.push_str(&format!(
            r#"<fill><patternFill patternType="solid"><fgColor rgb="FF{fill}"/><bgColor indexed="64"/></patternFill></fill>"#
        ));
        formats.push_str(&format!(
            r#"<xf numFmtId="0" fontId="1" fillId="{}" borderId="0" xfId="0" applyFont="1" applyFill="1" applyAlignment="1"><alignment vertical="center" wrapText="1"/></xf>"#,
            place + 2
        ));
    }
    format!(
        r#"<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><fonts count="2"><font><sz val="9"/><name val="Arial"/><family val="2"/></font><font><b/><sz val="9"/><color rgb="FFFFFFFF"/><name val="Arial"/><family val="2"/></font></fonts><fills count="{}"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill>{fills}</fills><borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders><cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs><cellXfs count="{}"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>{formats}</cellXfs><cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles></styleSheet>"#,
        FILLS.len() + 2,
        FILLS.len() + 1
    )
}

fn width(label: &str) -> f64 {
    if let Some((_, width)) = WIDTHS.iter().find(|(name, _)| *name == label) {
        return *width;
    }
    if label.starts_with("AD ") { 20.0 } else { 13.0 }
}

/// Writes the workbook of `inputs.table`.
pub(crate) fn write(inputs: &Inputs<'_>) -> Result<Summary> {
    let mut lines = BufReader::new(File::open(inputs.table)?).lines();
    let header_line = lines.next().transpose()?.unwrap_or_default();
    let header: Vec<&str> = if header_line.is_empty() {
        Vec::new()
    } else {
        header_line.split('\t').collect()
    };
    let mut sheet = plan(&header, inputs)?;
    let rank_place = sheet.places.get("variant_rank").copied();
    let id_place = sheet.places.get("variant_id").copied();

    let staged = PathBuf::from(format!("{}.partial", inputs.output.display()));
    let mut archive = zip::ZipWriter::new(BufWriter::new(File::create(&staged)?));
    let options = zip::write::FileOptions::default()
        .compression_method(zip::CompressionMethod::Deflated)
        // A fixed time, so the same table gives the same bytes.
        .last_modified_time(zip::DateTime::default());
    for (name, text) in [
        ("[Content_Types].xml", CONTENT_TYPES.to_owned()),
        ("_rels/.rels", ROOT_RELATIONS.to_owned()),
        ("xl/workbook.xml", WORKBOOK.to_owned()),
        ("xl/_rels/workbook.xml.rels", WORKBOOK_RELATIONS.to_owned()),
        ("xl/styles.xml", styles()),
    ] {
        archive.start_file(name, options)?;
        archive.write_all(text.as_bytes())?;
    }
    archive.start_file("xl/worksheets/sheet1.xml", options.large_file(true))?;
    let count = sheet.columns.len();
    let frozen = FROZEN_COLUMNS.min(count.saturating_sub(1));
    write!(
        archive,
        r#"<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetViews><sheetView workbookViewId="0"><pane xSplit="{frozen}" ySplit="1" topLeftCell="{}2" activePane="bottomRight" state="frozen"/></sheetView></sheetViews><sheetFormatPr defaultRowHeight="12"/><cols>"#,
        letters(frozen)
    )?;
    for (place, column) in sheet.columns.iter().enumerate() {
        write!(
            archive,
            r#"<col min="{0}" max="{0}" width="{1}" customWidth="1"/>"#,
            place + 1,
            width(&column.label)
        )?;
    }
    write!(
        archive,
        r#"</cols><sheetData><row r="1" ht="66" customHeight="1">"#
    )?;
    let mut groups: Vec<&str> = Vec::new();
    for (place, column) in sheet.columns.iter().enumerate() {
        if !groups.contains(&column.group) {
            groups.push(column.group);
        }
        let fill = groups
            .iter()
            .position(|group| *group == column.group)
            .unwrap_or(0)
            % FILLS.len();
        write!(
            archive,
            r#"<c r="{}1" s="{}" t="inlineStr"><is><t xml:space="preserve">{}</t></is></c>"#,
            letters(place),
            fill + 1,
            escaped(&column.label)
        )?;
    }
    write!(archive, "</row>")?;

    let (mut rows, mut variants, mut last_variant) = (0usize, 0usize, String::new());
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
    for line in lines {
        let line = line?;
        let row: Vec<&str> = line.split('\t').collect();
        if let Some(rank) = rank_place.and_then(|place| row.get(place)) {
            let rank = rank.parse::<f64>().map_err(|_| {
                invalid(format!(
                    "{}: variant_rank is not a number in row {}",
                    inputs.table.display(),
                    rows + 2
                ))
            })?;
            if rank as i64 > inputs.top_variants as i64 {
                break;
            }
        }
        let variant = id_place
            .and_then(|place| row.get(place))
            .copied()
            .unwrap_or("");
        if variant != last_variant {
            if seen.insert(variant.to_owned()) {
                variants += 1;
            }
            last_variant = variant.to_owned();
        }
        rows += 1;
        write!(archive, r#"<row r="{}">"#, rows + 1)?;
        for column in 0..count {
            match sheet.value(column, &row)? {
                Value::Empty => {}
                Value::Number(number) => write!(
                    archive,
                    r#"<c r="{}{}"><v>{number}</v></c>"#,
                    letters(column),
                    rows + 1
                )?,
                Value::Text(text) => write!(
                    archive,
                    r#"<c r="{}{}" t="inlineStr"><is><t xml:space="preserve">{}</t></is></c>"#,
                    letters(column),
                    rows + 1,
                    escaped(&text)
                )?,
            }
        }
        write!(archive, "</row>")?;
    }
    write!(
        archive,
        r#"</sheetData><autoFilter ref="A1:{}{}"/></worksheet>"#,
        letters(count.saturating_sub(1)),
        rows + 1
    )?;
    archive.finish()?.flush()?;
    fs::rename(&staged, inputs.output)?;
    Ok(Summary {
        rows,
        variants,
        columns: count,
    })
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args.is_empty() || args[0] == "--help" || args[0] == "-h" {
        println!("{USAGE}");
        return Ok(());
    }
    let mut options: HashMap<&str, &str> = HashMap::new();
    let mut rest = args.iter();
    while let Some(option) = rest.next() {
        let name = option
            .strip_prefix("--")
            .ok_or_else(|| invalid(format!("unexpected argument {option}\n{USAGE}")))?;
        let value = rest
            .next()
            .ok_or_else(|| invalid(format!("{option} needs a value\n{USAGE}")))?;
        options.insert(name, value);
    }
    let known = [
        "table",
        "output",
        "top-variants",
        "ped",
        "family",
        "gene-ranking",
        "expected-incidence",
        "pext-tissues",
        "data-dir",
    ];
    if let Some(name) = options.keys().find(|name| !known.contains(name)) {
        return Err(invalid(format!("unknown option --{name}\n{USAGE}")));
    }
    let need = |name: &str| {
        options
            .get(name)
            .map(|value| Path::new(*value))
            .ok_or_else(|| invalid(format!("--{name} is required\n{USAGE}")))
    };
    let data_dir = match options.get("data-dir") {
        Some(folder) => PathBuf::from(folder),
        None => crate::protein_domain_mapping::data_directory(),
    };
    let summary = write(&Inputs {
        table: need("table")?,
        output: need("output")?,
        top_variants: match options.get("top-variants") {
            Some(count) => count.parse()?,
            None => DEFAULT_TOP_VARIANTS,
        },
        ped: options.get("ped").map(|path| Path::new(*path)),
        family: options.get("family").copied(),
        gene_ranking: options.get("gene-ranking").map(|path| Path::new(*path)),
        expected_incidence: match options.get("expected-incidence") {
            Some(value) => value.parse()?,
            None => 0.0001,
        },
        pext_tissues: options.get("pext-tissues").copied().unwrap_or(""),
        data_dir: &data_dir,
    })?;
    println!(
        "review workbook: {} rows of {} variants, {} columns",
        summary.rows, summary.variants, summary.columns
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    #[test]
    fn table_values_become_numbers_or_stay_text() {
        let number = |text: &str| Value::Number(text.to_owned());
        assert_eq!(tidy("12"), number("12"));
        assert_eq!(tidy("-3"), number("-3"));
        assert_eq!(tidy("1.0"), number("1"));
        assert_eq!(tidy("0.25"), number("0.25"));
        assert_eq!(tidy("6.19587e-07"), number("6.19587e-7"));
        assert_eq!(tidy(".5"), number("0.5"));
        // A number with leading zeros still reads as a number.
        assert_eq!(tidy("007"), number("7"));
        for empty in EMPTY {
            assert_eq!(tidy(empty), Value::Empty);
        }
        for text in ["0/1", "chr1", "5.", "1e", "True", "ENST00000335137", "-"] {
            assert_eq!(tidy(text), Value::Text(text.to_owned()), "{text}");
        }
    }

    #[test]
    fn frequencies_read_as_one_in_n() {
        assert_eq!(fraction(Some(0.0001)), "1/10,000");
        assert_eq!(fraction(Some(0.00023)), "1/4,348");
        assert_eq!(fraction(Some(1.0e-6)), "1/1,000,000");
        assert_eq!(fraction(Some(0.02)), "1/50");
        assert_eq!(fraction(Some(1.5)), "1/1");
        assert_eq!(fraction(Some(0.0)), "");
        assert_eq!(fraction(None), "");
    }

    #[test]
    fn domains_lose_the_structure_mappings_and_are_cut() {
        assert_eq!(
            domains(
                "Gene3D:3.90.550.10&PDB-ENSP_mappings:7d72.E&Pfam:PF00483&Pfam:PF00483&AFDB-ENSP_mappings:AF-Q9Y5P6-F1"
            ),
            "Gene3D:3.90.550.10; Pfam:PF00483"
        );
        let long = (0..40)
            .map(|place| format!("Pfam:PF{place:05}"))
            .collect::<Vec<_>>()
            .join("&");
        let shown = domains(&long);
        assert_eq!(shown.chars().count(), 120);
        assert!(shown.ends_with("..."));
        assert_eq!(domains("nan"), "");
    }

    #[test]
    fn columns_are_lettered_as_in_excel() {
        assert_eq!(letters(0), "A");
        assert_eq!(letters(25), "Z");
        assert_eq!(letters(26), "AA");
        assert_eq!(letters(79), "CB");
        assert_eq!(letters(701), "ZZ");
        assert_eq!(letters(702), "AAA");
    }

    #[test]
    fn roles_come_from_the_probands_parents() {
        let members: HashMap<String, Member> = [
            ("kid", "dad", "mum", true),
            ("sis", "dad", "mum", true),
            ("bro", "dad", "mum", false),
            ("dad", "0", "0", false),
            ("mum", "0", "0", false),
            ("aunt", "0", "0", false),
        ]
        .into_iter()
        .map(|(sample, father, mother, affected)| {
            (
                sample.to_owned(),
                Member {
                    father: father.to_owned(),
                    mother: mother.to_owned(),
                    affected,
                },
            )
        })
        .collect();
        let patients = ["kid".to_owned(), "sis".to_owned()];
        let of = |sample: &str| role(sample, &patients, &members);
        assert_eq!(
            [
                of("kid"),
                of("dad"),
                of("mum"),
                of("sis"),
                of("bro"),
                of("aunt")
            ],
            [
                "proband",
                "father",
                "mother",
                "affected sibling",
                "sibling",
                "member"
            ]
        );
        assert_eq!(role("kid", &[], &members), "member");
    }

    /// A data folder with the three gene tables the constraint text reads.
    fn data_folder(root: &Path) -> PathBuf {
        let data = root.join("data");
        for (relative, text) in [
            (
                DEFAULT_CLINGEN_DOSAGE,
                "#Gene Symbol\tGene ID\tcytoBand\tGenomic Location\tHaploinsufficiency Score\tHaploinsufficiency Description\tTriplosensitivity Score\tTriplosensitivity Description\nG1\t1\t1p1\tchr1:1-2\t3\tSufficient evidence for dosage pathogenicity\t0\tNo evidence available\n",
            ),
            (
                DEFAULT_HGNC_TABLE,
                "hgnc_id\tsymbol\tname\tlocus_group\tlocus_type\tstatus\tlocation\tlocation_sortable\talias_symbol\talias_name\tprev_symbol\tprev_name\nHGNC:1\tG1\tgene one\tprotein-coding gene\tgene with protein product\tApproved\t1p1\t01p1\t\t\t\t\nHGNC:2\tG2\tgene two\tprotein-coding gene\tgene with protein product\tApproved\t1p1\t01p1\t\t\t\t\n",
            ),
        ] {
            let path = data.join(relative);
            fs::create_dir_all(path.parent().unwrap()).unwrap();
            fs::write(path, text).unwrap();
        }
        let loeuf = data.join(DEFAULT_LOEUF_TABLE);
        fs::create_dir_all(loeuf.parent().unwrap()).unwrap();
        let mut encoder = flate2::write::GzEncoder::new(
            File::create(&loeuf).unwrap(),
            flate2::Compression::default(),
        );
        encoder
            .write_all(b"#gene\ttranscript\tcanonical\toe_lof_upper\nG1\tT1\ttrue\t0.21\nG2\tT2\ttrue\t1.4\n")
            .unwrap();
        encoder.finish().unwrap();
        data
    }

    #[test]
    fn the_workbook_holds_the_top_variants_in_the_reviewers_order() {
        let root =
            std::env::temp_dir().join(format!("priva_review_workbook_{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&root).unwrap();
        let data = data_folder(&root);
        let header = [
            "chrom",
            "pos",
            "ref",
            "alt",
            "variant_id",
            "SYMBOL",
            "Feature",
            "Consequence",
            "CANONICAL",
            "NMD",
            "ACMG_criteria",
            "ACMG_class",
            "ACMG_quant_score",
            "sort_index",
            "variant_rank",
            "kid",
            "mum",
            "kid_AD",
            "mum_AD",
            "kid_phase_patho",
            "variant_inheritance_basis",
            "gnomAD_joint_AF",
            "gnomAD_joint_AF_max",
            "gnomAD_joint_AN",
            "gnomAD_joint_AN_max",
            "gnomAD_nhomalt_max",
            "gnomAD_nhomalt_XX",
            "gnomAD_nhomalt_XY",
            "clinvar_patho_gene_max_af",
            "audit_BS1_disease_incidence_gate",
            "audit_BS1_gene_pathogenic_AF_gate",
            "DOMAINS",
            "SpliceAI_pred_DS_AG",
            "SpliceAI_pred_DS_AL",
            "SpliceAI_pred_DS_DG",
            "SpliceAI_pred_DS_DL",
            "SpliceAI_pred_DP_AG",
            "SpliceAI_pred_DP_AL",
            "SpliceAI_pred_DP_DG",
            "SpliceAI_pred_DP_DL",
            "PEXT_MEAN",
            "PEXT_Muscle_Skeletal",
        ];
        let rows = [
            // Rank 1, two transcripts: a constrained gene without disease records.
            "chr1\t100\tA\tAT\tchr1:100:A-AT\tG1\tT1\tframeshift_variant\tYES\tNMD_escaping_variant\tPVS1;PM2\tLikely Pathogenic\t0.95\t1.14\t1\t0/1\t0/0\t10,9\t20,0\ttrans\tgene_constraint\t1e-06\t5e-05\t1000\t200\t0\t0\t0\t0.00023\t0\t1\tPfam:PF1&PDB-ENSP_mappings:1abc.A\t0.01\t0.75\tnan\t0.02\t5\t6\t7\t8\t0.9\t0.8",
            "chr1\t100\tA\tAT\tchr1:100:A-AT\tG2\tT2\tintron_variant\t\t\tPM2\tUncertain Significance\t0.5\t0.5\t1\t0/1\t0/0\t10,9\t20,0\t\tmatched_history\t2e-05\t\t1000\t\t\t1\t2\t\t1\t0\t\t\t\t\t\t\t\t\t\t\t",
            // Rank 2 is kept, rank 3 is beyond the two variants asked for.
            "chr2\t5\tC\tT\tchr2:5:C-T\tG2\tT2\tmissense_variant\tYES\t\tPM2\tUncertain Significance\t0.4\t0.4\t2\t1/1\t0/1\t0,30\t11,12\t\tmechanism_mismatch\t\t\t\t\t\t\t\t\t0\t0\t\t0.02\t0.03\t0.01\t0.0\t1\t2\t3\t4\t\t",
            "chr3\t7\tG\tA\tchr3:7:G-A\tG2\tT2\tmissense_variant\tYES\t\t\tBenign\t0.1\t0.1\t3\t0/1\t0/1\t5,5\t5,5\t\t\t\t\t\t\t\t\t\t\t0\t0\t\t\t\t\t\t\t\t\t\t\t",
        ];
        let table = root.join("family.tsv");
        fs::write(
            &table,
            format!("{}\n{}\n", header.join("\t"), rows.join("\n")),
        )
        .unwrap();
        let ped = root.join("family.ped");
        fs::write(&ped, "F1 kid 0 mum 1 2\nF1 mum 0 0 2 1\nF2 other 0 0 1 1\n").unwrap();
        let ranking = root.join("ranking.tsv");
        fs::write(
            &ranking,
            "sample_id\tgene_symbol\tcandidate_generating_hp_names\tsupport_only_hp_names\nkid\tG1\tMuscular dystrophy\tScoliosis\n",
        )
        .unwrap();
        let output = root.join("family.review.xlsx");
        let summary = write(&Inputs {
            table: &table,
            output: &output,
            top_variants: 2,
            ped: Some(&ped),
            family: Some("F1"),
            gene_ranking: Some(&ranking),
            expected_incidence: 0.0001,
            pext_tissues: "Muscle-Skeletal,Spleen",
            data_dir: &data,
        })
        .unwrap();
        assert_eq!(
            summary,
            Summary {
                rows: 3,
                variants: 2,
                columns: LAYOUT.len() + 2
            }
        );

        let mut archive = zip::ZipArchive::new(File::open(&output).unwrap()).unwrap();
        let mut xml = String::new();
        archive
            .by_name("xl/worksheets/sheet1.xml")
            .unwrap()
            .read_to_string(&mut xml)
            .unwrap();
        let has = |part: &str| assert!(xml.contains(part), "the sheet lacks {part}");
        // The header: fixed columns first, one genotype and depth column per member.
        has(r#"<c r="A1" s="1" t="inlineStr"><is><t xml:space="preserve">Chr</t></is></c>"#);
        has(r#"<t xml:space="preserve">GT kid (proband)</t>"#);
        has(r#"<t xml:space="preserve">AD mum (mother)</t>"#);
        has(r#"<t xml:space="preserve">Tissue of that PEXT value</t>"#);
        has(r#"<pane xSplit="15" ySplit="1" topLeftCell="P2""#);
        // Row 2: numbers as numbers, the condensed cells in words.
        has(r#"<c r="B2"><v>100</v></c>"#);
        has(r#"<c r="K2"><v>0.95</v></c>"#);
        has(r#"<t xml:space="preserve">Muscular dystrophy; supporting: Scoliosis</t>"#);
        has(
            r#"<t xml:space="preserve">gene constraint: ClinGen haploinsufficiency score 3 (sufficient evidence for dosage pathogenicity); LOEUF 0.21</t>"#,
        );
        has(r#"<t xml:space="preserve">1/10,000</t>"#);
        has(r#"<t xml:space="preserve">1/4,348</t>"#);
        has(r#"<t xml:space="preserve">acceptor loss, +6 bases from the variant</t>"#);
        has(r#"<t xml:space="preserve">Pfam:PF1</t>"#);
        has(r#"<t xml:space="preserve">Muscle-Skeletal</t>"#);
        // The highest frequency: the population's in row 2, whole gnomAD's in row 3
        // (homozygotes of both sexes added).
        has(r#"<v>5e-5</v>"#);
        has(r#"<v>2e-5</v>"#);
        has(r#"<v>3</v>"#);
        has(
            r#"<t xml:space="preserve">disease of this gene with the same mechanism as this variant</t>"#,
        );
        assert!(!xml.contains("BS1 applied"));
        // Row 4 (rank 2): a score below 0.1 names no effect; rank 3 is left out.
        has(r#"<t xml:space="preserve">none (below 0.1)</t>"#);
        has(r#"<t xml:space="preserve">no disease of this gene fits this variant's effect</t>"#);
        assert!(!xml.contains("chr3:7:G-A"));
        has(r#"<autoFilter ref="A1:"#);
        assert!(xml.ends_with("4\"/></worksheet>"));

        // The same table gives the same bytes.
        let first = fs::read(&output).unwrap();
        write(&Inputs {
            table: &table,
            output: &output,
            top_variants: 2,
            ped: Some(&ped),
            family: Some("F1"),
            gene_ranking: Some(&ranking),
            expected_incidence: 0.0001,
            pext_tissues: "Muscle-Skeletal,Spleen",
            data_dir: &data,
        })
        .unwrap();
        assert_eq!(first, fs::read(&output).unwrap());
        let _ = fs::remove_dir_all(&root);
    }
}
