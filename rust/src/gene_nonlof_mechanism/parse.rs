//! The script's evidence readers and the HGNC mapping, on pandas' string tables.
use super::*;
use crate::gene_mechanism_builder::{mechanism_confidence, norm_key, normalize_mechanism};
use crate::pandas_table::{Column, Table};
use crate::python_text::{is_python_space, python_strip};
use std::collections::BTreeSet;

pub(super) const UNIFIED: [&str; 21] = [
    "gene_symbol",
    "mechanism",
    "assertion_level",
    "source",
    "source_record_id",
    "disease",
    "inheritance",
    "confidence",
    "disease_confidence",
    "pmids",
    "patho_mode_raw",
    "assembly",
    "chromosome",
    "position",
    "ref",
    "alt",
    "hgvsc",
    "hgvsp",
    "consequence",
    "allele_key",
    "notes",
];

/// One `make_unified_row` record, in `UNIFIED` order.
#[derive(Clone, Debug)]
pub(super) struct Row(pub(super) Vec<String>);

impl Row {
    pub(super) fn get(&self, name: &str) -> &str {
        &self.0[UNIFIED.iter().position(|n| *n == name).unwrap()]
    }
}

const PANELAPP_NON_LOF_LABEL: &str = "PANELAPP_GREEN_NON_LOF_PATHO_HISTORY";
const KEPT_LABELS: [&str; 6] = [
    "GOF",
    "ACTIVATING",
    "DOMINANT_NEGATIVE",
    "TRIPLOSENSITIVITY",
    "INCREASED_DOSAGE",
    PANELAPP_NON_LOF_LABEL,
];

/// The script's `make_unified_row`.
fn unified(fields: &[(&str, String)], mechanism: &[String]) -> Row {
    let mut values = vec![String::new(); UNIFIED.len()];
    values[2] = "gene_level".into();
    for (name, value) in fields {
        values[UNIFIED.iter().position(|n| n == name).unwrap()] = value.clone();
    }
    values[0] = python_strip(&values[0]).to_owned();
    let labels: BTreeSet<&str> = mechanism.iter().map(String::as_str).collect();
    values[1] = labels.into_iter().collect::<Vec<_>>().join("|");
    Row(values)
}

/// The script's `normalize_mechanism(raw, source)`.
fn mechanism(raw: &str, source: &str) -> Vec<String> {
    if source == "PanelApp" {
        let text = raw.to_lowercase().replace(['_', '-'], " ");
        if text.contains("do not cause")
            && ["loss of function", "loss function", "lof"]
                .iter()
                .any(|term| text.contains(term))
        {
            return vec![PANELAPP_NON_LOF_LABEL.into()];
        }
        return vec![];
    }
    normalize_mechanism(raw).into_iter().map(str::to_owned).collect()
}

/// A pandas string table as text cells (`fillna("")`).
struct Frame {
    names: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl Frame {
    fn from_table(table: Table) -> Frame {
        let mut rows = vec![Vec::with_capacity(table.names.len()); table.rows];
        for column in &table.columns {
            for (row, value) in rows.iter_mut().enumerate() {
                value.push(match column {
                    Column::Text(values) => values[row].clone().unwrap_or_default(),
                    _ => String::new(),
                });
            }
        }
        Frame {
            names: table.names,
            rows,
        }
    }
    /// The script's `resolve_column`: the last column of each normalized name.
    fn resolve(&self, candidates: &[&str]) -> Option<usize> {
        let mut map: HashMap<String, usize> = HashMap::new();
        for (index, name) in self.names.iter().enumerate() {
            map.insert(norm_key(name), index);
        }
        candidates.iter().find_map(|c| map.get(&norm_key(c)).copied())
    }
    /// `row.get(name, "")` by label (the first column carrying it).
    fn by_name<'a>(&self, row: &'a [String], name: &str) -> &'a str {
        self.names
            .iter()
            .position(|n| n == name)
            .map_or("", |index| row[index].as_str())
    }
}

/// `_val(row, col)`.
fn val(row: &[String], column: Option<usize>) -> String {
    column.map_or(String::new(), |index| python_strip(&row[index]).to_owned())
}

fn parse_g2p(path: &Path) -> Result<Vec<Row>> {
    if !exists(path) {
        return Ok(vec![]);
    }
    let frame = Frame::from_table(crate::gene_mechanism_builder::table::read(path, b',', true)?);
    let gene = frame.resolve(&["gene symbol", "gene_symbol"]);
    let disease = frame.resolve(&["disease name", "disease_label", "disease"]);
    let consequence = frame.resolve(&[
        "molecular mechanism",
        "mutation consequence",
        "mutation_consequence",
    ]);
    let confidence = frame.resolve(&["confidence", "confidence category", "g2p relation label"]);
    let inheritance = frame.resolve(&["allelic requirement", "allelic_requirement"]);
    let pmids = frame.resolve(&["pmids", "publications"]);
    let hgnc = frame.resolve(&["hgnc id", "hgnc_id"]);
    let mut rows = Vec::new();
    for row in &frame.rows {
        let symbol = val(row, gene);
        if symbol.is_empty() {
            continue;
        }
        let raw = val(row, consequence);
        let level = val(row, confidence);
        rows.push(unified(
            &[
                ("gene_symbol", symbol),
                ("source", "G2P_DDG2P".into()),
                ("source_record_id", val(row, hgnc)),
                ("disease", val(row, disease)),
                ("inheritance", val(row, inheritance)),
                ("confidence", mechanism_confidence(&level).into()),
                ("disease_confidence", level),
                ("pmids", val(row, pmids)),
                ("patho_mode_raw", raw.clone()),
            ],
            &mechanism(&raw, "G2P"),
        ));
    }
    Ok(rows)
}

fn parse_panelapp(path: &Path) -> Result<Vec<Row>> {
    if !exists(path) {
        return Ok(vec![]);
    }
    let payload = read_json(path)?;
    let panels = get_or(&payload, "panels", Json::Array(vec![]))?;
    let panels: Vec<Json> = iterate(&panels)?;
    let mut rows = Vec::new();
    for panel in &panels {
        let name = py_str(&get_or(panel, "name", text(""))?);
        let id = py_str(&get_or(panel, "id", text(""))?);
        let version = py_str(&get_or(panel, "version", text(""))?);
        let genes = get_or(panel, "genes", Json::Array(vec![]))?;
        let genes = if truthy(&genes) { genes } else { Json::Array(vec![]) };
        let Json::Array(genes) = genes else {
            continue;
        };
        for entry in &genes {
            let data = get(entry, "gene_data")?.cloned().unwrap_or(Json::Null);
            let data = if truthy(&data) { data } else { Json::Object(vec![]) };
            let mut symbol = Json::Null;
            for (source, key) in [
                (&data, "gene_symbol"),
                (&data, "hgnc_symbol"),
                (entry, "entity_name"),
                (entry, "gene_symbol"),
            ] {
                let value = get(source, key)?.cloned().unwrap_or(Json::Null);
                if truthy(&value) {
                    symbol = value;
                    break;
                }
            }
            let gene = if truthy(&symbol) { py_str(&symbol) } else { String::new() };
            let gene = python_strip(&gene).to_owned();
            if gene.is_empty() {
                continue;
            }
            let first = |keys: &[&str]| -> Result<String> {
                for key in keys {
                    let value = get(entry, key)?.cloned().unwrap_or(Json::Null);
                    if truthy(&value) {
                        return Ok(py_str(&value));
                    }
                }
                Ok(String::new())
            };
            let mop = python_strip(&first(&[
                "mode_of_pathogenicity",
                "mode_of_pathogenicity_description",
            ])?)
            .to_owned();
            let moi = python_strip(&first(&["mode_of_inheritance"])?).to_owned();
            let publications = get(entry, "publications")?.cloned().unwrap_or(Json::Null);
            let publications = if truthy(&publications) {
                iterate(&publications)?
            } else {
                vec![]
            };
            let pmids = publications
                .iter()
                .filter(|item| truthy(item))
                .map(py_str)
                .collect::<Vec<_>>()
                .join("|");
            let confidence = first(&["confidence_level", "penetrance"])?;
            let labels = mechanism(&mop, "PanelApp");
            if python_strip(&confidence) != "3" || !labels.iter().any(|l| l == PANELAPP_NON_LOF_LABEL)
            {
                continue;
            }
            rows.push(unified(
                &[
                    ("gene_symbol", gene),
                    ("source", "PanelApp".into()),
                    ("source_record_id", format!("panel:{id};version:{version}")),
                    ("disease", name.clone()),
                    ("inheritance", moi),
                    ("confidence", mechanism_confidence(&confidence).into()),
                    ("disease_confidence", confidence.clone()),
                    ("pmids", pmids),
                    ("patho_mode_raw", mop),
                ],
                &labels,
            ));
        }
    }
    Ok(rows)
}

/// `for item in value` over a parsed JSON value.
fn iterate(value: &Json) -> Result<Vec<Json>> {
    match value {
        Json::Array(items) => Ok(items.clone()),
        Json::Object(items) => Ok(items.iter().map(|(k, _)| text(k.clone())).collect()),
        Json::Str(chars) => Ok(chars.chars().map(|c| text(c.to_string())).collect()),
        other => Err(error(
            "TypeError",
            format!("'{}' object is not iterable", type_name(other)),
        )),
    }
}

fn parse_clingen(path: &Path) -> Result<Vec<Row>> {
    if !exists(path) {
        return Ok(vec![]);
    }
    let bytes = fs::read(path).map_err(|e| os_error(e, path, None))?;
    let body = crate::python_io::decode_utf8(&bytes).map_err(|e| error("UnicodeDecodeError", e))?;
    let body = body.strip_prefix('\u{feff}').unwrap_or(body);
    let mut skip = 0;
    for line in universal_lines(body) {
        if !line.starts_with('#') || line.contains('\t') {
            break;
        }
        skip += 1;
    }
    let mut frame = Frame::from_table(crate::gene_mechanism_builder::table::read_skipping(
        path, b'\t', true, skip,
    )?);
    if frame.names.first().is_some_and(|name| name.starts_with('#')) {
        frame.names[0] = python_strip(frame.names[0].trim_start_matches('#')).to_owned();
    }
    let gene = frame.resolve(&["Gene Symbol", "#Gene Symbol", "gene_symbol", "symbol"]);
    let disease = frame.resolve(&["Disease Name", "Disease", "disease_label"]);
    let hi_score = frame.resolve(&["Haploinsufficiency Score", "HI Score"]);
    let hi_desc = frame.resolve(&[
        "Haploinsufficiency Description",
        "HI Evidence Strength",
        "Haploinsufficiency",
    ]);
    let ts_score = frame.resolve(&["Triplosensitivity Score", "TS Score"]);
    let ts_desc = frame.resolve(&[
        "Triplosensitivity Description",
        "TS Evidence Strength",
        "Triplosensitivity",
    ]);
    let gene_id = frame.resolve(&["Gene ID", "HGNC ID"]);
    let pmid_columns: Vec<usize> = (0..frame.names.len())
        .filter(|index| norm_key(&frame.names[*index]).contains("pmid"))
        .collect();
    let mut rows = Vec::new();
    for row in &frame.rows {
        let symbol = val(row, gene);
        if symbol.is_empty() {
            continue;
        }
        let disease = val(row, disease);
        let (hs, hd, ts, td) = (
            val(row, hi_score),
            val(row, hi_desc),
            val(row, ts_score),
            val(row, ts_desc),
        );
        let values: Vec<String> = pmid_columns
            .iter()
            .map(|index| python_strip(&row[*index]).to_owned())
            .filter(|value| !value.is_empty())
            .collect();
        let pmids = if values.is_empty() {
            String::new()
        } else {
            let joined = values.join("|");
            let unique: BTreeSet<&str> = joined.split('|').collect();
            unique.into_iter().collect::<Vec<_>>().join("|")
        };
        let id = val(row, gene_id);
        for (score, description, label, prefix) in [
            (&hs, &hd, "HAPLOINSUFFICIENCY", "HI_score"),
            (&ts, &td, "TRIPLOSENSITIVITY", "TS_score"),
        ] {
            if score.is_empty() && description.is_empty() {
                continue;
            }
            let labels = if ["3", "2"].contains(&python_strip(score)) {
                vec![label.to_owned()]
            } else {
                vec![]
            };
            let raw = format!("{prefix}={score}; {description}");
            rows.push(unified(
                &[
                    ("gene_symbol", symbol.clone()),
                    ("source", "ClinGen_Dosage".into()),
                    ("source_record_id", id.clone()),
                    ("disease", disease.clone()),
                    ("confidence", mechanism_confidence(description).into()),
                    ("disease_confidence", description.clone()),
                    ("pmids", pmids.clone()),
                    ("patho_mode_raw", raw.trim_matches([';', ' ']).to_owned()),
                ],
                &labels,
            ));
        }
    }
    Ok(rows)
}

/// `for line in open(path)` with universal newlines.
fn universal_lines(body: &str) -> impl Iterator<Item = &str> {
    let mut rest = body;
    std::iter::from_fn(move || {
        if rest.is_empty() {
            return None;
        }
        let end = rest.find(['\r', '\n']).map_or(rest.len(), |at| {
            if rest[at..].starts_with("\r\n") { at + 2 } else { at + 1 }
        });
        let (line, tail) = rest.split_at(end);
        rest = tail;
        Some(line)
    })
}

fn placeholder(value: &str) -> String {
    if ["_", "-", "nan", "NA", "N/A"].contains(&python_strip(value)) {
        String::new()
    } else {
        value.to_owned()
    }
}

/// `re.sub(r"^chr", "", text, flags=re.I)`.
pub(super) fn without_chr(value: &str) -> String {
    match value.get(..3) {
        Some(lead) if lead.eq_ignore_ascii_case("chr") => value[3..].to_owned(),
        _ => value.to_owned(),
    }
}

fn allele_key(frame: &Frame, row: &[String]) -> String {
    let reference = python_strip(&placeholder(frame.by_name(row, "ref"))).to_uppercase();
    let alternate = python_strip(&placeholder(frame.by_name(row, "alt"))).to_uppercase();
    let kind = if reference.chars().count() == 1 && alternate.chars().count() == 1 {
        "SNV"
    } else {
        "Indel"
    };
    let chrom = without_chr(python_strip(frame.by_name(row, "chr")));
    let integer = |name: &str| {
        let value = python_strip(frame.by_name(row, name));
        value.strip_suffix(".0").unwrap_or(value).to_owned()
    };
    [
        kind.to_owned(),
        chrom,
        integer("hg19start"),
        integer("hg19end"),
        reference,
        alternate,
    ]
    .join("|")
}

fn parse_gofcards(path: &Path) -> Result<Vec<Row>> {
    if !exists(path) {
        return Ok(vec![]);
    }
    let (names, cells) = super::xlsx::read_sheet(path, "total3161")?;
    let frame = Frame { names, rows: cells };
    let mut rows = Vec::new();
    for row in &frame.rows {
        let get = |name: &str| frame.by_name(row, name).to_owned();
        let gene = python_strip(&get("genesymbol")).to_owned();
        if gene.is_empty() {
            continue;
        }
        let pscore_raw = python_strip(&get("Pscore")).to_owned();
        let pscore = crate::python_text::parse_float(&pscore_raw).unwrap_or(0.0);
        let confidence = if pscore >= 5.0 {
            "high"
        } else if pscore >= 3.0 {
            "moderate"
        } else {
            "limited"
        };
        let function = placeholder(&get("Function"));
        let disorder = placeholder(&get("Disorder involved"));
        let mut pmid = python_strip(&get("PMID")).to_owned();
        if pmid == "nan" {
            pmid.clear();
        }
        let chrom_raw = python_strip(&get("chr")).to_owned();
        let chrom = if !chrom_raw.is_empty() && !chrom_raw.starts_with("chr") {
            format!("chr{chrom_raw}")
        } else {
            chrom_raw
        };
        let position = python_strip(&get("hg19start")).to_owned();
        let reference = python_strip(&get("ref")).to_owned();
        let alternate = python_strip(&get("alt")).to_owned();
        let transcript = python_strip(&get("transcript")).to_owned();
        let hgvsc = if !transcript.is_empty() && !chrom.is_empty() && !position.is_empty() {
            format!("{transcript}:{chrom}:g.{position}{reference}>{alternate}")
        } else {
            String::new()
        };
        let consequence = if reference.chars().count() == 1
            && alternate.chars().count() == 1
            && !["-", ""].contains(&reference.as_str())
            && !["-", ""].contains(&alternate.as_str())
        {
            "missense_variant"
        } else {
            "indel"
        };
        let pathway = placeholder(&get("Pathways proteins involved"));
        let animal = placeholder(&get("Animal model"));
        let cell = placeholder(&get("Cell model"));
        let notes: Vec<String> = [("pathway", pathway), ("animal", animal), ("cell", cell)]
            .into_iter()
            .filter(|(_, value)| !value.is_empty())
            .map(|(name, value)| format!("{name}:{value}"))
            .collect();
        rows.push(unified(
            &[
                ("gene_symbol", gene),
                ("assertion_level", "variant_level".into()),
                ("source", "GoFCards".into()),
                ("source_record_id", get("Order numbe")),
                ("disease", disorder),
                ("confidence", confidence.into()),
                (
                    "disease_confidence",
                    if pscore_raw.is_empty() {
                        String::new()
                    } else {
                        format!("Pscore={pscore_raw}")
                    },
                ),
                ("pmids", pmid),
                (
                    "patho_mode_raw",
                    if function.is_empty() {
                        "curated GOF variant".into()
                    } else {
                        function
                    },
                ),
                ("assembly", "hg19".into()),
                ("chromosome", chrom),
                ("position", position),
                ("ref", reference),
                ("alt", alternate),
                ("hgvsc", hgvsc),
                ("consequence", consequence.into()),
                ("allele_key", allele_key(&frame, row)),
                ("notes", notes.join("; ")),
            ],
            &["GOF".to_owned()],
        ));
    }
    Ok(rows)
}

/// The script's `parse_all_sources`.
pub(super) fn parse_all_sources(raw_dir: &Path, specs: &[Spec]) -> Result<Vec<Row>> {
    let mut unified_rows = Vec::new();
    for spec in specs {
        let path = raw_dir.join(spec.raw_filename);
        let parsed = match spec.parser {
            "g2p_ddg2p" => parse_g2p(&path),
            "panelapp" => parse_panelapp(&path),
            "clingen_dosage" => parse_clingen(&path),
            "gofcards" => parse_gofcards(&path),
            _ => continue,
        };
        match parsed {
            Ok(rows) => unified_rows.extend(rows),
            Err(failure) => {
                if spec.custom() {
                    return Err(error(
                        "ValueError",
                        format!(
                            "Incompatible selected source {} ({}): {}",
                            spec.name,
                            spec.url,
                            message(failure.as_ref())
                        ),
                    ));
                }
                progress(&format!(
                    "parse_error source={}: {}",
                    spec.name,
                    message(failure.as_ref())
                ));
            }
        }
    }
    if unified_rows.is_empty() {
        return Ok(unified_rows);
    }
    let total = unified_rows.len();
    let kept: Vec<Row> = unified_rows
        .into_iter()
        .filter(|row| {
            let value = row.get("mechanism");
            !python_strip(value).is_empty()
                && value.split('|').any(|label| KEPT_LABELS.contains(&label))
        })
        .collect();
    if total != kept.len() {
        progress(&format!(
            "filter: dropped {} rows without curated non-LOF mechanism",
            total - kept.len()
        ));
    }
    Ok(kept)
}

#[derive(Clone, Debug)]
pub(super) struct Hgnc {
    pub(super) hgnc_id: String,
    pub(super) symbol: String,
    pub(super) entrez_id: String,
    pub(super) ensembl_id: String,
}

/// The script's `load_hgnc_mapping`.
pub(super) fn load_hgnc_mapping(path: &Path) -> Result<HashMap<String, Hgnc>> {
    if !exists(path) {
        return Ok(HashMap::new());
    }
    let frame = Frame::from_table(crate::gene_mechanism_builder::table::read(path, b'\t', false)?);
    let mut columns: HashMap<String, usize> = HashMap::new();
    for (index, name) in frame.names.iter().enumerate() {
        columns.insert(name.to_lowercase(), index);
    }
    let find = |names: &[&str], required: bool| -> Result<Option<usize>> {
        for name in names {
            if let Some(index) = columns.get(&name.to_lowercase()) {
                return Ok(Some(*index));
            }
        }
        if required {
            return Err(error(
                "ValueError",
                format!(
                    "{} lacks required HGNC column; expected one of ({})",
                    path.display(),
                    names
                        .iter()
                        .map(|n| crate::python_path::repr(n))
                        .collect::<Vec<_>>()
                        .join(", ")
                ),
            ));
        }
        Ok(None)
    };
    let hgnc_column = find(&["hgnc_id", "HGNC ID"], true)?;
    let symbol_column = find(&["symbol", "Approved symbol"], true)?;
    let alias_column = find(&["alias_symbol", "Alias symbols"], false)?;
    let previous_column = find(&["prev_symbol", "Previous symbols"], false)?;
    let entrez_column = find(&["entrez_id", "NCBI Gene ID(supplied by NCBI)"], false)?;
    let ensembl_column = find(&["ensembl_gene_id", "Ensembl ID(supplied by Ensembl)"], false)?;
    let cell = |row: &[String], column: Option<usize>| column.map_or("", |i| row[i].as_str()).to_owned();
    let mut lookup: HashMap<String, Hgnc> = HashMap::new();
    for row in &frame.rows {
        let entrez = python_strip(&cell(row, entrez_column)).to_owned();
        let entry = Hgnc {
            hgnc_id: python_strip(&cell(row, hgnc_column)).to_owned(),
            symbol: python_strip(&cell(row, symbol_column)).to_owned(),
            entrez_id: entrez.strip_suffix(".0").unwrap_or(&entrez).to_owned(),
            ensembl_id: python_strip(&cell(row, ensembl_column)).to_owned(),
        };
        let symbol = entry.symbol.to_uppercase();
        if !symbol.is_empty() {
            lookup.insert(symbol, entry.clone());
        }
        for column in [alias_column, previous_column] {
            for name in cell(row, column).split(['|', ',']) {
                let name = python_strip(name).to_uppercase();
                if !name.is_empty() && !lookup.contains_key(&name) {
                    lookup.insert(name, entry.clone());
                }
            }
        }
    }
    Ok(lookup)
}

/// The script's `_split_pmids`.
pub(super) fn split_pmids(raw: &str) -> Vec<String> {
    if python_strip(raw).is_empty() {
        return vec![];
    }
    python_strip(raw)
        .split(|c: char| c == '|' || c == ';' || c == ',' || is_python_space(c))
        .map(python_strip)
        .filter(|part| {
            let digits = part.replace("PMID:", "");
            !part.is_empty() && !digits.is_empty() && digits.chars().all(super::fetch::is_digit)
        })
        .map(str::to_owned)
        .collect()
}
