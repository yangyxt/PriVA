//! `prepare-gene-evidence`: `scripts/prepare_gene_evidence_resources.py`,
//! which turns the downloaded gene-evidence files into simple term-to-gene
//! tables and a `manifest.json` of what it found.

use crate::hpo_gene_ranking::python_text::{path_join, path_text, strip};
use crate::hpo_gene_ranking::tsv::{TableWriter, for_each_line, python_dict_rows, read_bytes};
use crate::hpo_gene_ranking::xml_tree::{Element, read_top_level};
use crate::hpo_gene_ranking::{Options, Result};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::Path;

const HELP: &str = "\
hpo-gene-ranking prepare-gene-evidence: term-to-gene tables from the raw downloads

  --raw-dir DIR   downloaded files (data/gene_evidence/raw)
  --out-dir DIR   prepared tables and manifest.json (data/gene_evidence/prepared)";

/// A JSON value written as `json.dumps(value, indent=2, sort_keys=True)` writes it.
enum Json {
    Text(String),
    Int(i64),
    List(Vec<Json>),
    Object(BTreeMap<String, Json>),
}

fn json_string(text: &str) -> String {
    let mut out = String::from("\"");
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{8}' => out.push_str("\\b"),
            '\u{c}' => out.push_str("\\f"),
            c if (c as u32) < 0x20 || (c as u32) > 0x7e => {
                let mut units = [0u16; 2];
                for unit in c.encode_utf16(&mut units) {
                    out.push_str(&format!("\\u{unit:04x}"));
                }
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

impl Json {
    fn write(&self, level: usize, out: &mut String) {
        let indent = |level: usize| "  ".repeat(level);
        match self {
            Json::Text(text) => out.push_str(&json_string(text)),
            Json::Int(value) => out.push_str(&value.to_string()),
            Json::List(items) if items.is_empty() => out.push_str("[]"),
            Json::Object(items) if items.is_empty() => out.push_str("{}"),
            Json::List(items) => {
                out.push_str("[\n");
                for (index, item) in items.iter().enumerate() {
                    if index > 0 {
                        out.push_str(",\n");
                    }
                    out.push_str(&indent(level + 1));
                    item.write(level + 1, out);
                }
                out.push('\n');
                out.push_str(&indent(level));
                out.push(']');
            }
            Json::Object(items) => {
                out.push_str("{\n");
                for (index, (key, item)) in items.iter().enumerate() {
                    if index > 0 {
                        out.push_str(",\n");
                    }
                    out.push_str(&indent(level + 1));
                    out.push_str(&json_string(key));
                    out.push_str(": ");
                    item.write(level + 1, out);
                }
                out.push('\n');
                out.push_str(&indent(level));
                out.push('}');
            }
        }
    }
}

fn object(items: Vec<(&str, Json)>) -> Json {
    Json::Object(
        items
            .into_iter()
            .map(|(key, value)| (key.to_owned(), value))
            .collect(),
    )
}

fn missing(file: &str) -> Json {
    object(vec![
        ("status", Json::Text("missing".to_owned())),
        ("file", Json::Text(file.to_owned())),
    ])
}

/// `read_tsv`: lines of `open(path, newline="")` without their "\n", split
/// on tabs; lines starting with `comment` are skipped.
fn each_tsv_row(
    path: &str,
    comment: Option<&str>,
    mut handle: impl FnMut(Vec<&str>) -> Result<()>,
) -> Result<()> {
    for_each_line(path, false, |line| {
        if comment.is_some_and(|comment| line.starts_with(comment)) {
            return Ok(());
        }
        let line = line.trim_end_matches('\n');
        handle(line.split('\t').collect())
    })
}

struct Writer {
    writer: TableWriter,
    count: i64,
}

impl Writer {
    fn new(path: &str, header: &[&str]) -> Result<Writer> {
        let mut writer = TableWriter::create(path, "\n")?;
        writer.write_row(header)?;
        Ok(Writer { writer, count: 0 })
    }

    fn row<S: AsRef<str>>(&mut self, fields: &[S]) -> Result<()> {
        self.count += 1;
        self.writer.write_row(fields)
    }

    fn finish(self) -> Result<i64> {
        self.writer.finish()?;
        Ok(self.count)
    }
}

fn prepare_hpo(raw: &str, out: &str) -> Result<Json> {
    let src = path_join(raw, "phenotype_to_genes.txt");
    if !Path::new(&src).exists() {
        return Ok(missing(&src));
    }
    let mut writer = Writer::new(
        &path_join(out, "hp_to_genes.tsv"),
        &["hp_id", "hp_name", "ncbi_gene_id", "gene_symbol", "source"],
    )?;
    let mut genes_by_hp: HashMap<String, HashSet<String>> = HashMap::new();
    each_tsv_row(&src, Some("#"), |fields| {
        if fields.len() < 4 || fields[0] == "hpo_id" {
            return Ok(());
        }
        genes_by_hp
            .entry(fields[0].to_owned())
            .or_default()
            .insert(fields[3].to_owned());
        writer.row(&[
            fields[0],
            fields[1],
            fields[2],
            fields[3],
            "HPO_phenotype_to_genes",
        ])
    })?;
    let rows = writer.finish()?;
    let genes: HashSet<&String> = genes_by_hp.values().flatten().collect();
    Ok(object(vec![
        ("status", Json::Text("ok".to_owned())),
        ("rows", Json::Int(rows)),
        ("unique_hp_terms", Json::Int(genes_by_hp.len() as i64)),
        ("unique_gene_symbols", Json::Int(genes.len() as i64)),
        ("file", Json::Text(src)),
    ]))
}

fn prepare_goa(raw: &str, out: &str) -> Result<Json> {
    let src = path_join(raw, "goa_human.gaf");
    if !Path::new(&src).exists() {
        return Ok(missing(&src));
    }
    let mut writer = Writer::new(
        &path_join(out, "go_to_genes.tsv"),
        &[
            "go_id",
            "aspect",
            "gene_symbol",
            "db",
            "db_object_id",
            "db_object_name",
            "db_object_type",
            "evidence_code",
            "reference",
            "with_or_from",
            "taxon",
            "source",
        ],
    )?;
    let mut genes_by_go: HashMap<String, HashSet<String>> = HashMap::new();
    each_tsv_row(&src, Some("!"), |f| {
        if f.len() < 15 {
            return Ok(());
        }
        if f[3].split('|').any(|part| part == "NOT") {
            return Ok(());
        }
        genes_by_go
            .entry(f[4].to_owned())
            .or_default()
            .insert(f[2].to_owned());
        writer.row(&[
            f[4],
            f[8],
            f[2],
            f[0],
            f[1],
            f[9],
            f[11],
            f[6],
            f[5],
            f[7],
            f[12],
            "GOA_human_gaf",
        ])
    })?;
    let rows = writer.finish()?;
    let genes: HashSet<&String> = genes_by_go.values().flatten().collect();
    Ok(object(vec![
        ("status", Json::Text("ok".to_owned())),
        ("rows", Json::Int(rows)),
        ("unique_go_terms", Json::Int(genes_by_go.len() as i64)),
        ("unique_gene_symbols", Json::Int(genes.len() as i64)),
        ("file", Json::Text(src)),
    ]))
}

fn prepare_gene_info(raw: &str, out: &str) -> Result<Json> {
    let src = path_join(raw, "Homo_sapiens.gene_info");
    if !Path::new(&src).exists() {
        return Ok(missing(&src));
    }
    let mut writer = Writer::new(
        &path_join(out, "human_gene_info.tsv"),
        &[
            "ncbi_gene_id",
            "gene_symbol",
            "synonyms",
            "db_xrefs",
            "chromosome",
            "map_location",
            "description",
            "gene_type",
            "nomenclature_symbol",
            "nomenclature_name",
            "tax_id",
        ],
    )?;
    each_tsv_row(&src, Some("#"), |f| {
        if f.len() < 12 {
            return Ok(());
        }
        writer.row(&[
            f[1], f[2], f[4], f[5], f[6], f[7], f[8], f[9], f[10], f[11], f[0],
        ])
    })?;
    let rows = writer.finish()?;
    Ok(object(vec![
        ("status", Json::Text("ok".to_owned())),
        ("rows", Json::Int(rows)),
        ("file", Json::Text(src)),
    ]))
}

fn prepare_reactome(raw: &str, out: &str) -> Result<Json> {
    let src = path_join(raw, "NCBI2Reactome_All_Levels.txt");
    if !Path::new(&src).exists() {
        return Ok(missing(&src));
    }
    let mut writer = Writer::new(
        &path_join(out, "reactome_to_ncbi_genes.tsv"),
        &[
            "reactome_id",
            "pathway_name",
            "ncbi_gene_id",
            "evidence_code",
            "url",
            "species",
            "source",
        ],
    )?;
    let mut genes_by_pathway: HashMap<String, HashSet<String>> = HashMap::new();
    each_tsv_row(&src, None, |f| {
        if f.len() < 6 || f[5] != "Homo sapiens" {
            return Ok(());
        }
        genes_by_pathway
            .entry(f[1].to_owned())
            .or_default()
            .insert(f[0].to_owned());
        writer.row(&[
            f[1],
            f[3],
            f[0],
            f[4],
            f[2],
            f[5],
            "Reactome_NCBI2Reactome_All_Levels",
        ])
    })?;
    let rows = writer.finish()?;
    let genes: HashSet<&String> = genes_by_pathway.values().flatten().collect();
    Ok(object(vec![
        ("status", Json::Text("ok".to_owned())),
        ("rows", Json::Int(rows)),
        (
            "unique_reactome_pathways",
            Json::Int(genes_by_pathway.len() as i64),
        ),
        ("unique_ncbi_gene_ids", Json::Int(genes.len() as i64)),
        ("file", Json::Text(src)),
    ]))
}

/// The accession of `^(UniProtKB|UniProtKB/Swiss-Prot|UniProtKB/TrEMBL):([A-Z0-9]+)`.
fn uniprot_accession(xref: &str) -> Option<&str> {
    for prefix in ["UniProtKB:", "UniProtKB/Swiss-Prot:", "UniProtKB/TrEMBL:"] {
        if let Some(rest) = xref.strip_prefix(prefix) {
            let end = rest
                .find(|c: char| !(c.is_ascii_uppercase() || c.is_ascii_digit()))
                .unwrap_or(rest.len());
            if end > 0 {
                return Some(&rest[..end]);
            }
        }
    }
    None
}

/// The number of `^NCBIGene:(\d+)`.
fn ncbi_gene(xref: &str) -> Option<&str> {
    let rest = xref.strip_prefix("NCBIGene:")?;
    let end = rest
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(rest.len());
    (end > 0).then(|| &rest[..end])
}

fn prepare_pro(raw: &str, out: &str) -> Result<Json> {
    let src = path_join(raw, "pr.obo");
    if !Path::new(&src).exists() {
        return Ok(missing(&src));
    }
    let mut writer = Writer::new(
        &path_join(out, "pr_to_gene_or_protein.tsv"),
        &[
            "pr_id",
            "pr_name",
            "xref_db",
            "uniprot_accession",
            "ncbi_gene_id",
            "source",
        ],
    )?;
    #[derive(Default)]
    struct Block {
        keys: usize,
        ids: Vec<String>,
        names: Vec<String>,
        xrefs: Vec<String>,
    }
    let mut pr_terms: HashSet<String> = HashSet::new();
    let mut emit = |block: &Block, writer: &mut Writer| -> Result<()> {
        let Some(pr_id) = block.ids.first().filter(|id| id.starts_with("PR:")) else {
            return Ok(());
        };
        let name = block.names.first().map(String::as_str).unwrap_or("");
        for xref in &block.xrefs {
            let value = xref
                .split_once(" ! ")
                .map_or(xref.as_str(), |(value, _)| value);
            if let Some(accession) = uniprot_accession(value) {
                pr_terms.insert(pr_id.clone());
                writer.row(&[
                    pr_id.as_str(),
                    name,
                    "UniProtKB",
                    accession,
                    "",
                    "PRO_pr_obo_xref",
                ])?;
            } else if let Some(gene) = ncbi_gene(value) {
                pr_terms.insert(pr_id.clone());
                writer.row(&[
                    pr_id.as_str(),
                    name,
                    "NCBIGene",
                    "",
                    gene,
                    "PRO_pr_obo_xref",
                ])?;
            }
        }
        Ok(())
    };
    let mut block: Option<Block> = None;
    for_each_line(&src, true, |raw_line| {
        let line = raw_line.trim_end_matches('\n');
        if line == "[Term]" {
            if let Some(done) = block.take().filter(|b| b.keys > 0) {
                emit(&done, &mut writer)?;
            }
            block = Some(Block::default());
            return Ok(());
        }
        if line.starts_with('[') {
            if let Some(done) = block.take().filter(|b| b.keys > 0) {
                emit(&done, &mut writer)?;
            }
            return Ok(());
        }
        let Some(current) = block.as_mut() else {
            return Ok(());
        };
        let Some((key, value)) = line.split_once(": ") else {
            return Ok(());
        };
        current.keys += 1;
        match key {
            "id" => current.ids.push(value.to_owned()),
            "name" => current.names.push(value.to_owned()),
            "xref" => current.xrefs.push(value.to_owned()),
            _ => {}
        }
        Ok(())
    })?;
    if let Some(done) = block.take().filter(|b| b.keys > 0) {
        emit(&done, &mut writer)?;
    }
    let rows = writer.finish()?;
    Ok(object(vec![
        ("status", Json::Text("ok".to_owned())),
        ("rows", Json::Int(rows)),
        ("unique_pr_terms", Json::Int(pr_terms.len() as i64)),
        ("file", Json::Text(src)),
    ]))
}

fn prepare_mgi(raw: &str, out: &str) -> Result<Json> {
    let src = path_join(raw, "MGI_GenePheno.rpt");
    let homology = path_join(raw, "HOM_MouseHumanSequence.rpt");
    if !Path::new(&src).exists() {
        return Ok(missing(&src));
    }
    struct Mouse {
        class_key: String,
        symbol: String,
        entrez: String,
    }
    struct Human {
        symbol: String,
        entrez: String,
        hgnc: String,
        omim: String,
    }
    let mut mouse_by_mgi: HashMap<String, Mouse> = HashMap::new();
    let mut human_by_class: HashMap<String, Vec<Human>> = HashMap::new();
    let homology_present = Path::new(&homology).exists();
    if homology_present {
        each_tsv_row(&homology, None, |f| {
            if f.len() < 8 || f[0] == "DB Class Key" {
                return Ok(());
            }
            if f[2] == "10090" && !f[5].is_empty() {
                mouse_by_mgi.insert(
                    f[5].to_owned(),
                    Mouse {
                        class_key: f[0].to_owned(),
                        symbol: f[3].to_owned(),
                        entrez: f[4].to_owned(),
                    },
                );
            } else if f[2] == "9606" {
                human_by_class
                    .entry(f[0].to_owned())
                    .or_default()
                    .push(Human {
                        symbol: f[3].to_owned(),
                        entrez: f[4].to_owned(),
                        hgnc: f[6].to_owned(),
                        omim: f[7].to_owned(),
                    });
            }
            Ok(())
        })?;
    }
    let mut mouse_writer = Writer::new(
        &path_join(out, "mp_to_mouse_genes.tsv"),
        &[
            "mp_id",
            "mgi_marker_id",
            "mouse_gene_symbol",
            "mouse_entrez_gene_id",
            "allelic_composition",
            "allele_symbols",
            "allele_ids",
            "genetic_background",
            "pubmed_ids",
            "mgi_genotype_id",
            "homology_class_key",
            "source",
        ],
    )?;
    let mut human_writer = Writer::new(
        &path_join(out, "mp_to_human_orthologs.tsv"),
        &[
            "mp_id",
            "mouse_mgi_marker_id",
            "mouse_gene_symbol",
            "mouse_entrez_gene_id",
            "human_gene_symbol",
            "human_entrez_gene_id",
            "hgnc_id",
            "omim_gene_id",
            "homology_class_key",
            "pubmed_ids",
            "mgi_genotype_id",
            "source",
        ],
    )?;
    let mut mouse_genes_by_mp: HashMap<String, HashSet<String>> = HashMap::new();
    let mut human_genes_by_mp: HashMap<String, HashSet<String>> = HashMap::new();
    let no_mouse = Mouse {
        class_key: String::new(),
        symbol: String::new(),
        entrez: String::new(),
    };
    each_tsv_row(&src, None, |f| {
        if f.len() < 8 {
            return Ok(());
        }
        let (allelic, symbols, allele_ids, background) = (f[0], f[1], f[2], f[3]);
        let (mp_id, pubmed, markers, genotype) = (f[4], f[5], f[6], f[7]);
        if mp_id.is_empty() {
            return Ok(());
        }
        for marker in markers
            .split('|')
            .map(strip)
            .filter(|marker| !marker.is_empty())
        {
            let mouse = mouse_by_mgi.get(marker).unwrap_or(&no_mouse);
            if !mouse.symbol.is_empty() {
                mouse_genes_by_mp
                    .entry(mp_id.to_owned())
                    .or_default()
                    .insert(mouse.symbol.clone());
            }
            mouse_writer.row(&[
                mp_id,
                marker,
                mouse.symbol.as_str(),
                mouse.entrez.as_str(),
                allelic,
                symbols,
                allele_ids,
                background,
                pubmed,
                genotype,
                mouse.class_key.as_str(),
                "MGI_GenePheno",
            ])?;
            for human in human_by_class.get(&mouse.class_key).into_iter().flatten() {
                if !human.symbol.is_empty() {
                    human_genes_by_mp
                        .entry(mp_id.to_owned())
                        .or_default()
                        .insert(human.symbol.clone());
                }
                human_writer.row(&[
                    mp_id,
                    marker,
                    mouse.symbol.as_str(),
                    mouse.entrez.as_str(),
                    human.symbol.as_str(),
                    human.entrez.as_str(),
                    human.hgnc.as_str(),
                    human.omim.as_str(),
                    mouse.class_key.as_str(),
                    pubmed,
                    genotype,
                    "MGI_GenePheno+HOM_MouseHumanSequence",
                ])?;
            }
        }
        Ok(())
    })?;
    let mouse_rows = mouse_writer.finish()?;
    let human_rows = human_writer.finish()?;
    let mouse_genes: HashSet<&String> = mouse_genes_by_mp.values().flatten().collect();
    let human_genes: HashSet<&String> = human_genes_by_mp.values().flatten().collect();
    Ok(object(vec![
        ("status", Json::Text("ok".to_owned())),
        ("mouse_rows", Json::Int(mouse_rows)),
        ("human_ortholog_rows", Json::Int(human_rows)),
        (
            "unique_mp_terms_mouse",
            Json::Int(mouse_genes_by_mp.len() as i64),
        ),
        (
            "unique_mouse_gene_symbols",
            Json::Int(mouse_genes.len() as i64),
        ),
        (
            "unique_mp_terms_human",
            Json::Int(human_genes_by_mp.len() as i64),
        ),
        (
            "unique_human_gene_symbols",
            Json::Int(human_genes.len() as i64),
        ),
        ("file", Json::Text(src)),
        ("homology_file", Json::Text(homology)),
        (
            "homology_status",
            Json::Text(if homology_present { "ok" } else { "missing" }.to_owned()),
        ),
    ]))
}

fn prepare_uniprot(raw: &str, out: &str) -> Result<Json> {
    let src = path_join(raw, "uniprot_human_reviewed.tsv");
    if !Path::new(&src).exists() {
        return Ok(missing(&src));
    }
    let (_, rows) = python_dict_rows(&read_bytes(&src)?);
    let mut writer = Writer::new(
        &path_join(out, "uniprot_human_reviewed_evidence.tsv"),
        &[
            "uniprot_accession",
            "primary_gene_symbol",
            "gene_names",
            "protein_names",
            "function_text",
            "subcellular_location_text",
            "tissue_specificity_text",
            "disease_text",
            "ncbi_gene_ids",
        ],
    )?;
    let columns = [
        "Entry",
        "Gene Names (primary)",
        "Gene Names",
        "Protein names",
        "Function [CC]",
        "Subcellular location [CC]",
        "Tissue specificity",
        "Involvement in disease",
        "GeneID",
    ];
    for row in &rows {
        let fields: Vec<&str> = columns
            .iter()
            .map(|column| row.get(*column).map(String::as_str).unwrap_or(""))
            .collect();
        writer.row(&fields)?;
    }
    let count = writer.finish()?;
    Ok(object(vec![
        ("status", Json::Text("ok".to_owned())),
        ("rows", Json::Int(count)),
        ("file", Json::Text(src)),
    ]))
}

fn prepare_hpa(raw: &str, out: &str) -> Result<Json> {
    let mut summaries: Vec<(&str, Json)> = Vec::new();
    for filename in [
        "rna_tissue_consensus.tsv",
        "rna_single_cell_type.tsv",
        "rna_immune_cell_schmiedel.tsv",
    ] {
        let src = path_join(raw, filename);
        if !Path::new(&src).exists() {
            summaries.push((filename, missing(&src)));
            continue;
        }
        let (fieldnames, rows) = python_dict_rows(&read_bytes(&src)?);
        let header: Vec<&str> = fieldnames.iter().map(String::as_str).collect();
        let mut writer = Writer::new(&path_join(out, filename), &header)?;
        for row in &rows {
            let fields: Vec<&str> = fieldnames
                .iter()
                .map(|column| row.get(column).map(String::as_str).unwrap_or(""))
                .collect();
            writer.row(&fields)?;
        }
        let count = writer.finish()?;
        summaries.push((
            filename,
            object(vec![
                ("status", Json::Text("ok".to_owned())),
                ("rows", Json::Int(count)),
                (
                    "columns",
                    Json::List(fieldnames.into_iter().map(Json::Text).collect()),
                ),
                ("file", Json::Text(src)),
            ]),
        ));
    }
    Ok(object(summaries))
}

fn stripped_text(element: &Element) -> String {
    strip(&element.text).to_owned()
}

fn prepare_orphadata(raw: &str, out: &str) -> Result<Json> {
    let src = path_join(raw, "orphadata_en_product6.xml");
    if !Path::new(&src).exists() {
        return Ok(missing(&src));
    }
    let mut writer = Writer::new(
        &path_join(out, "orpha_disease_to_genes.tsv"),
        &[
            "orpha_disease_id",
            "disease_name",
            "orpha_gene_id",
            "gene_symbol",
            "gene_name",
            "association_type",
            "association_status",
            "external_refs",
            "source",
        ],
    )?;
    let mut diseases: HashSet<String> = HashSet::new();
    let mut genes: HashSet<String> = HashSet::new();
    let mut handle_top = |top: Element| -> Result<()> {
        let mut disorders = Vec::new();
        top.iter_tag("Disorder", &mut disorders);
        for disorder in disorders {
            let mut disorder_orpha = String::new();
            let mut disorder_name = String::new();
            for child in &disorder.children {
                if child.local == "OrphaCode" {
                    disorder_orpha = stripped_text(child);
                } else if child.local == "Name" {
                    disorder_name = stripped_text(child);
                }
            }
            let mut associations = Vec::new();
            disorder.iter_tag("DisorderGeneAssociation", &mut associations);
            for association in associations {
                let Some(gene) = association.find("Gene") else {
                    continue;
                };
                let mut gene_orpha = String::new();
                let mut gene_symbol = String::new();
                let mut gene_name = String::new();
                let mut references: Vec<String> = Vec::new();
                for child in &gene.children {
                    match child.local.as_str() {
                        "OrphaCode" => gene_orpha = stripped_text(child),
                        "Symbol" => gene_symbol = stripped_text(child),
                        "Name" => gene_name = stripped_text(child),
                        "ExternalReferenceList" => {
                            let mut found = Vec::new();
                            child.iter_tag("ExternalReference", &mut found);
                            for reference in found {
                                let source = strip(&reference.findtext("Source")).to_owned();
                                let value = strip(&reference.findtext("Reference")).to_owned();
                                if !source.is_empty() && !value.is_empty() {
                                    references.push(format!("{source}:{value}"));
                                }
                            }
                        }
                        _ => {}
                    }
                }
                let association_type =
                    strip(&association.findtext("DisorderGeneAssociationType/Name")).to_owned();
                let association_status =
                    strip(&association.findtext("DisorderGeneAssociationStatus/Name")).to_owned();
                let disease_id = if disorder_orpha.is_empty() {
                    String::new()
                } else {
                    format!("ORPHA:{disorder_orpha}")
                };
                let gene_id = if gene_orpha.is_empty() {
                    String::new()
                } else {
                    format!("ORPHA:{gene_orpha}")
                };
                if !disease_id.is_empty() {
                    diseases.insert(disease_id.clone());
                }
                if !gene_symbol.is_empty() {
                    genes.insert(gene_symbol.clone());
                }
                let joined_references = references.join("|");
                writer.row(&[
                    disease_id.as_str(),
                    disorder_name.as_str(),
                    gene_id.as_str(),
                    gene_symbol.as_str(),
                    gene_name.as_str(),
                    association_type.as_str(),
                    association_status.as_str(),
                    joined_references.as_str(),
                    "Orphadata_product6",
                ])?;
            }
        }
        Ok(())
    };
    read_top_level(&src, &mut handle_top)?;
    let rows = writer.finish()?;
    Ok(object(vec![
        ("status", Json::Text("ok".to_owned())),
        ("rows", Json::Int(rows)),
        ("unique_orpha_diseases", Json::Int(diseases.len() as i64)),
        ("unique_gene_symbols", Json::Int(genes.len() as i64)),
        ("file", Json::Text(src)),
    ]))
}

pub fn run(arguments: &[String]) -> Result<()> {
    let Some(options) = Options::parse(arguments, &["raw-dir", "out-dir"], &[])? else {
        println!("{HELP}");
        return Ok(());
    };
    let raw = path_text(options.value("raw-dir").unwrap_or("data/gene_evidence/raw"));
    let out = path_text(
        options
            .value("out-dir")
            .unwrap_or("data/gene_evidence/prepared"),
    );
    std::fs::create_dir_all(&out).map_err(|error| format!("{out}: {error}"))?;
    let resources = object(vec![
        ("hpo", prepare_hpo(&raw, &out)?),
        ("goa", prepare_goa(&raw, &out)?),
        ("gene_info", prepare_gene_info(&raw, &out)?),
        ("reactome", prepare_reactome(&raw, &out)?),
        ("pro", prepare_pro(&raw, &out)?),
        ("mgi_mp", prepare_mgi(&raw, &out)?),
        ("uniprot", prepare_uniprot(&raw, &out)?),
        ("hpa_expression", prepare_hpa(&raw, &out)?),
        ("orphadata", prepare_orphadata(&raw, &out)?),
    ]);
    let manifest = object(vec![
        ("raw_dir", Json::Text(raw.clone())),
        ("out_dir", Json::Text(out.clone())),
        ("resources", resources),
    ]);
    let mut text = String::new();
    manifest.write(0, &mut text);
    let manifest_path = path_join(&out, "manifest.json");
    std::fs::write(&manifest_path, text).map_err(|error| format!("{manifest_path}: {error}"))?;
    println!("wrote={manifest_path}");
    Ok(())
}
