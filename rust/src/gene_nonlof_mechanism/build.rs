//! The GoFCards cache reader and the canonical non-LOF assertion JSON.
use super::parse::{Hgnc, Row, split_pmids, without_chr};
use super::*;
use crate::python_text::python_strip;
use std::collections::{BTreeMap, BTreeSet};

/// `clinvar_vcv.load_gofcards_cache`.
pub(super) fn load_gofcards_cache(path: &Path) -> Result<Json> {
    let cache = load_json(path)?;
    let has_genes = match &cache {
        Json::Object(_) => cache.get("genes").is_some(),
        Json::Array(items) => items.iter().any(|item| equal(item, &text("genes"))),
        Json::Str(body) => body.contains("genes"),
        other => {
            return Err(error(
                "TypeError",
                format!("argument of type '{}' is not iterable", type_name(other)),
            ));
        }
    };
    if !has_genes {
        return Err(error(
            "ValueError",
            format!("GoFCards cache has no 'genes' block: {}", path.display()),
        ));
    }
    let meta = get(&cache, "metadata")?.cloned().unwrap_or(Json::Null);
    let meta = if truthy(&meta) {
        meta
    } else {
        Json::Object(vec![])
    };
    let source = get(&meta, "source")?.cloned().unwrap_or(Json::Null);
    let mechanism = get(&meta, "mechanism")?.cloned().unwrap_or(Json::Null);
    if !equal(&source, &text("GoFCards")) || !equal(&mechanism, &text("GOF")) {
        return Err(error(
            "ValueError",
            format!(
                "GoFCards cache does not declare GoFCards/GOF: {}",
                path.display()
            ),
        ));
    }
    Ok(cache)
}

fn items(value: &Json) -> Result<Vec<(String, Json)>> {
    if !truthy(value) {
        return Ok(vec![]);
    }
    match value {
        Json::Object(entries) => Ok(entries.clone()),
        other => Err(error(
            "AttributeError",
            format!("'{}' object has no attribute 'items'", type_name(other)),
        )),
    }
}

/// `clinvar_vcv.gofcards_variant_is_eligible`.
pub(super) fn eligible(variant: &Json) -> Result<bool> {
    let record = get(variant, "record")?.cloned().unwrap_or(Json::Null);
    let record = if truthy(&record) {
        record
    } else {
        Json::Object(vec![])
    };
    let eligibility = get(&record, "eligibility")?.cloned().unwrap_or(Json::Null);
    let eligibility = if truthy(&eligibility) {
        eligibility
    } else {
        Json::Object(vec![])
    };
    Ok(get(&eligibility, "status")?.is_some_and(|status| equal(status, &text("eligible"))))
}

fn quarantine_status(variant: &Json) -> Result<&'static str> {
    let record = get(variant, "record")?.cloned().unwrap_or(Json::Null);
    let record = if truthy(&record) {
        record
    } else {
        Json::Object(vec![])
    };
    let eligibility = get(&record, "eligibility")?.cloned().unwrap_or(Json::Null);
    let eligibility = if truthy(&eligibility) {
        eligibility
    } else {
        Json::Object(vec![])
    };
    let status = py_str(&get_or(&eligibility, "status", text(""))?);
    let status = python_strip(&status).to_lowercase();
    Ok(
        if status.starts_with("quarantined_reviewed_")
            || status.starts_with("quarantined_mechanism_")
        {
            "quarantined_upstream_mechanism_review"
        } else {
            "quarantined_upstream_gene_discordance"
        },
    )
}

fn resolved_hgnc_id(symbol: &str, hgnc: &HashMap<String, Hgnc>) -> String {
    let cleaned = python_strip(symbol).to_uppercase();
    if cleaned.is_empty() {
        return String::new();
    }
    hgnc.get(&cleaned).map_or_else(
        || format!("SYMBOL:{cleaned}"),
        |entry| entry.hgnc_id.clone(),
    )
}

pub(super) type Exact = HashMap<(String, String), Json>;

/// The script's `load_gofcards_exact_records`.
pub(super) fn load_gofcards_exact_records(
    path: &Path,
    hgnc: &HashMap<String, Hgnc>,
) -> Result<(Exact, Json)> {
    let cache = load_gofcards_cache(path)?;
    let mut index = Exact::new();
    let mut genes = BTreeSet::new();
    let (mut eligible_count, mut quarantined) = (0i128, 0i128);
    for (symbol, gene) in items(&get(&cache, "genes")?.cloned().unwrap_or(Json::Null))? {
        let variants = get(&gene, "variants")?.cloned().unwrap_or(Json::Null);
        for (variant_id, variant) in items(&variants)? {
            if symbol.is_empty() || variant_id.is_empty() {
                return Err(error(
                    "ValueError",
                    "Every GoFCards variant requires a gene symbol and an identifier",
                ));
            }
            if eligible(&variant)? {
                eligible_count += 1;
            } else {
                quarantined += 1;
            }
            genes.insert(python_strip(&symbol).to_uppercase());
            let Json::Object(fields) = &variant else {
                return Err(error(
                    "TypeError",
                    format!("'{}' object is not a mapping", type_name(&variant)),
                ));
            };
            let mut merged = object([
                ("source", text("GoFCards")),
                ("mechanism", text("GOF")),
                ("variant_id", text(variant_id.clone())),
            ]);
            for (key, value) in fields {
                merged.set(key, value.clone());
            }
            index.insert((resolved_hgnc_id(&symbol, hgnc), variant_id), merged);
        }
    }
    let stats = object([
        ("gofcards_genes", Json::Int(genes.len() as i128)),
        ("gofcards_variants", Json::Int(eligible_count + quarantined)),
        ("gofcards_eligible_variants", Json::Int(eligible_count)),
        ("gofcards_quarantined_variants", Json::Int(quarantined)),
    ]);
    Ok((index, stats))
}

fn clean_text(value: &str) -> String {
    let text = python_strip(value);
    if ["nan", "none", "_", ".", "-"].contains(&text.to_lowercase().as_str()) {
        String::new()
    } else {
        text.to_owned()
    }
}

/// `clinvar_vcv.variant_id_of` for text fields.
pub(super) fn variant_id_of(chrom: &str, start: &str, reference: &str, alternate: &str) -> String {
    let allele = |value: &str| {
        let text = clean_text(value);
        let text = if ["-", "."].contains(&text.as_str()) {
            String::new()
        } else {
            text.to_uppercase()
        };
        if text.is_empty() {
            "-".to_owned()
        } else {
            text
        }
    };
    let start = python_strip(start);
    let start = start.strip_suffix(".0").unwrap_or(start);
    format!(
        "loc_{}:{start}:{}->{}_grch37",
        without_chr(&clean_text(chrom)),
        allele(reference),
        allele(alternate)
    )
}

/// `_nonempty`: drop `None`, `""` and `[]`.
fn nonempty(fields: Vec<(&str, Json)>) -> Json {
    Json::Object(
        fields
            .into_iter()
            .filter(|(_, value)| match value {
                Json::Null => false,
                Json::Str(text) => !text.is_empty(),
                Json::Array(items) => !items.is_empty(),
                _ => true,
            })
            .map(|(key, value)| (key.to_owned(), value))
            .collect(),
    )
}

fn retain_nonempty(value: Json) -> Json {
    match value {
        Json::Object(fields) => nonempty(
            fields
                .iter()
                .map(|(k, v)| (k.as_str(), v.clone()))
                .collect(),
        ),
        other => other,
    }
}

fn pmids(row: &Row) -> Json {
    Json::Array(
        split_pmids(row.get("pmids"))
            .into_iter()
            .map(Json::Str)
            .collect(),
    )
}

fn assertion(source: &str, row: &Row) -> Option<Json> {
    let field = |name: &str| text(row.get(name));
    Some(match source {
        "G2P_DDG2P" => nonempty(vec![
            ("mechanism", field("mechanism")),
            ("mechanism_raw", field("patho_mode_raw")),
            ("disease", field("disease")),
            ("inheritance", field("inheritance")),
            ("confidence", field("confidence")),
            ("pmids", pmids(row)),
        ]),
        "PanelApp" => {
            let record = row.get("source_record_id");
            let panel_id = record.split_once("panel:").map_or("", |(_, rest)| {
                let between = rest.split("panel:").next().unwrap_or("");
                between.split(';').next().unwrap_or("")
            });
            nonempty(vec![
                ("mechanism", field("mechanism")),
                ("panel", field("disease")),
                ("panel_id", text(panel_id)),
                ("inheritance", field("inheritance")),
                ("confidence", field("confidence")),
            ])
        }
        "ClinGen_Dosage" => {
            let mut score = "";
            for part in row.get("patho_mode_raw").split(';') {
                let part = python_strip(part);
                if part.starts_with("TS_score=") || part.starts_with("HI_score=") {
                    score = python_strip(part.split('=').nth(1).unwrap_or(""));
                    break;
                }
            }
            nonempty(vec![
                ("mechanism", field("mechanism")),
                ("score", text(score)),
                ("score_description", field("disease_confidence")),
                ("pmids", pmids(row)),
            ])
        }
        "GoFCards" => {
            let confidence = row.get("disease_confidence");
            let pscore = confidence
                .strip_prefix("Pscore=")
                .and_then(|rest| {
                    crate::python_text::parse_float(rest.split('=').next().unwrap_or(""))
                })
                .map_or(Json::Null, Json::Float);
            let hgvsc = row.get("hgvsc");
            let transcript = if !hgvsc.is_empty() && hgvsc.contains(':') {
                hgvsc.split(':').next().unwrap_or("")
            } else {
                ""
            };
            let raw = row.get("patho_mode_raw");
            let mut value = nonempty(vec![
                ("mechanism", field("mechanism")),
                ("source_record_id", field("source_record_id")),
                ("allele_key", field("allele_key")),
                ("disease", field("disease")),
                ("chr", field("chromosome")),
                ("pos", field("position")),
                ("ref", field("ref")),
                ("alt", field("alt")),
                ("transcript", text(transcript)),
                ("consequence", field("consequence")),
                ("pscore", pscore),
                ("pmids", pmids(row)),
                (
                    "function",
                    text(if raw != "curated GOF variant" {
                        raw
                    } else {
                        ""
                    }),
                ),
                ("pathway", text("")),
                ("animal_model", text("")),
                ("cell_model", text("")),
            ]);
            for part in row.get("notes").split(';') {
                let part = python_strip(part);
                if let Some(rest) = part.strip_prefix("pathway:") {
                    value.set("pathway", text(python_strip(rest)));
                } else if let Some(rest) = part.strip_prefix("animal:") {
                    if python_strip(rest).to_uppercase() == "Y" {
                        value.set("animal_model", Json::Bool(true));
                    }
                } else if let Some(rest) = part.strip_prefix("cell:") {
                    if python_strip(rest).to_uppercase() == "Y" {
                        value.set("cell_model", Json::Bool(true));
                    }
                }
            }
            value
        }
        _ => return None,
    })
}

struct Gene {
    ids: Json,
    mechanisms: BTreeSet<String>,
    gene_level: Vec<Json>,
    variant_level: Vec<Json>,
}

fn strings(values: &[&str]) -> Json {
    Json::Array(values.iter().map(|v| text(*v)).collect())
}

fn meta(genes: usize, unmapped: Vec<String>) -> Json {
    object([
        ("version", text("2.0")),
        ("built_at", text(iso_now())),
        ("total_genes", Json::Int(genes as i128)),
        (
            "unmapped_symbols",
            Json::Array(unmapped.into_iter().map(Json::Str).collect()),
        ),
        (
            "sources",
            object([
                (
                    "G2P_DDG2P",
                    object([
                        ("level", text("gene_level")),
                        (
                            "keys",
                            strings(&[
                                "mechanism",
                                "mechanism_raw",
                                "disease",
                                "inheritance",
                                "confidence",
                                "pmids",
                            ]),
                        ),
                        (
                            "mechanism_raw",
                            text("original G2P molecular mechanism text"),
                        ),
                        (
                            "confidence_values",
                            text(
                                "definitive | strong | moderate | limited | conflicting_or_refuted",
                            ),
                        ),
                    ]),
                ),
                (
                    "GoFCards",
                    object([
                        ("level", text("variant_level")),
                        (
                            "keys",
                            strings(&[
                                "mechanism",
                                "source_record_id",
                                "allele_key",
                                "exact_normalization_status",
                                "exact_cache_gene_symbols",
                                "exact_normalized_variants",
                                "disease",
                                "chr",
                                "pos",
                                "ref",
                                "alt",
                                "transcript",
                                "consequence",
                                "pscore",
                                "pmids",
                                "function",
                                "pathway",
                                "animal_model",
                                "cell_model",
                            ]),
                        ),
                        (
                            "exact_normalized_variant_shape",
                            text(
                                "gene -> variants -> <loc_...._grch37> -> {record, assemblies -> hg19|hg38 -> {genomic, transcripts}}",
                            ),
                        ),
                        (
                            "raw_public_allele_fields",
                            object([
                                ("assembly", text("hg19")),
                                ("keys", strings(&["chr", "pos", "ref", "alt"])),
                                (
                                    "note",
                                    text(
                                        "Applies only to the legacy top-level GoFCards source fields; exact_normalized_variants identify hg19 and hg38 fields explicitly",
                                    ),
                                ),
                            ]),
                        ),
                        (
                            "exact_match_methods",
                            strings(&[
                                "normalized gene symbol plus exact HGVSp",
                                "normalized gene symbol plus exact hg19/hg38 genomic allele",
                            ]),
                        ),
                        (
                            "pscore",
                            text(
                                "GoFCards pathogenicity score (float); >=5 high, >=3 moderate, else limited",
                            ),
                        ),
                        ("animal_model", text("true if animal model evidence exists")),
                        ("cell_model", text("true if cell model evidence exists")),
                    ]),
                ),
                (
                    "PanelApp",
                    object([
                        ("level", text("gene_level")),
                        (
                            "keys",
                            strings(&[
                                "mechanism",
                                "panel",
                                "panel_id",
                                "inheritance",
                                "confidence",
                            ]),
                        ),
                        (
                            "note",
                            text(
                                "non-LOF history marker; LOF variants do NOT cause phenotype in this panel",
                            ),
                        ),
                    ]),
                ),
                (
                    "ClinGen_Dosage",
                    object([
                        ("level", text("gene_level")),
                        (
                            "keys",
                            strings(&["mechanism", "score", "score_description", "pmids"]),
                        ),
                        (
                            "score",
                            text("ClinGen dosage score; 3=sufficient, 2=emerging evidence"),
                        ),
                    ]),
                ),
                (
                    "ClinVar_VCV",
                    object([
                        ("level", text("variant_level")),
                        (
                            "scope",
                            text(
                                "condition-specific germline RCV classifications and linked SCV observations for exact normalized GoFCards alleles",
                            ),
                        ),
                        ("minimum_review_stars", Json::Int(2)),
                        (
                            "eligible_review_statuses",
                            strings(&[
                                "criteria provided, multiple submitters, no conflicts",
                                "reviewed by expert panel",
                                "practice guideline",
                            ]),
                        ),
                        (
                            "mechanism_caveat",
                            text(
                                "ClinVar does not assign GOF here; GOF context comes from the exact matched GoFCards allele.",
                            ),
                        ),
                        (
                            "inheritance_caveat",
                            text(
                                "Observed zygosity and submitted mode of inheritance are retained as evidence and do not set allelic_requirement.",
                            ),
                        ),
                    ]),
                ),
            ]),
        ),
        (
            "mechanism_values",
            strings(&[
                "GOF",
                "DOMINANT_NEGATIVE",
                "TRIPLOSENSITIVITY",
                "PANELAPP_GREEN_NON_LOF_PATHO_HISTORY",
            ]),
        ),
    ])
}

/// The script's `build_nonlof_assertions_json`.
pub(super) fn assertions_json(
    rows: &[Row],
    hgnc: &HashMap<String, Hgnc>,
    exact: &Exact,
) -> Result<Json> {
    let mut genes: BTreeMap<String, Gene> = BTreeMap::new();
    let mut unmapped = Vec::new();
    for row in rows {
        let symbol = python_strip(row.get("gene_symbol")).to_owned();
        if symbol.is_empty() {
            continue;
        }
        let (id, ids) = match hgnc.get(&symbol.to_uppercase()) {
            None => {
                unmapped.push(symbol.clone());
                (
                    format!("SYMBOL:{symbol}"),
                    nonempty(vec![("symbol", text(symbol.clone()))]),
                )
            }
            Some(entry) => (
                entry.hgnc_id.clone(),
                nonempty(vec![
                    ("symbol", text(entry.symbol.clone())),
                    ("entrez_id", text(entry.entrez_id.clone())),
                    ("ensembl_id", text(entry.ensembl_id.clone())),
                ]),
            ),
        };
        let gene = genes.entry(id.clone()).or_insert_with(|| Gene {
            ids,
            mechanisms: BTreeSet::new(),
            gene_level: vec![],
            variant_level: vec![],
        });
        let source = row.get("source");
        let Some(mut value) = assertion(source, row) else {
            continue;
        };
        if source == "GoFCards" {
            let variant_id = variant_id_of(
                row.get("chromosome"),
                row.get("position"),
                row.get("ref"),
                row.get("alt"),
            );
            match exact.get(&(id.clone(), variant_id)) {
                None => value.set(
                    "exact_normalization_status",
                    text("unmatched_public_source_allele"),
                ),
                Some(variant) if eligible(variant)? => {
                    value.set(
                        "exact_normalization_status",
                        text("matched_gene_concordant"),
                    );
                    value.set(
                        "exact_normalized_variants",
                        Json::Array(vec![variant.clone()]),
                    );
                }
                Some(variant) => {
                    value.set(
                        "exact_normalization_status",
                        text(quarantine_status(variant)?),
                    );
                    value.set(
                        "quarantined_exact_normalized_variants",
                        Json::Array(vec![variant.clone()]),
                    );
                }
            }
            value = retain_nonempty(value);
        }
        for label in row.get("mechanism").split('|') {
            let label = python_strip(label);
            if !label.is_empty() {
                gene.mechanisms.insert(label.to_owned());
            }
        }
        let entry = Json::Object(vec![(source.to_owned(), value)]);
        if source == "GoFCards" {
            gene.variant_level.push(entry);
        } else {
            gene.gene_level.push(entry);
        }
    }
    let unmapped: BTreeSet<String> = unmapped.into_iter().collect();
    let mut result = vec![(
        "_meta".to_owned(),
        meta(genes.len(), unmapped.into_iter().collect()),
    )];
    for (id, gene) in genes {
        let mut value = gene.ids;
        value.set(
            "mechanisms",
            Json::Array(gene.mechanisms.into_iter().map(Json::Str).collect()),
        );
        if !gene.gene_level.is_empty() {
            value.set("gene_level", Json::Array(gene.gene_level));
        }
        if !gene.variant_level.is_empty() {
            value.set("variant_level", Json::Array(gene.variant_level));
        }
        result.push((id, value));
    }
    Ok(Json::Object(result))
}

fn meta_mut(result: &mut Json) -> &mut Json {
    match result {
        Json::Object(items) => &mut items[0].1,
        _ => unreachable!("the canonical result is an object"),
    }
}

/// The GoFCards and ClinVar integration summaries; returns the number of
/// variants that carry ClinVar conditions.
pub(super) fn integrate(
    result: &mut Json,
    path: &Path,
    stats: &Json,
    exact: &Exact,
) -> Result<usize> {
    let mut assertions = Vec::new();
    if let Json::Object(items) = &*result {
        for (id, gene) in items {
            if id == "_meta" {
                continue;
            }
            if let Some(Json::Array(entries)) = gene.get("variant_level") {
                for entry in entries {
                    if let Some(value) = entry.get("GoFCards") {
                        assertions.push(value.clone());
                    }
                }
            }
        }
    }
    let mut counts: BTreeMap<String, i128> = BTreeMap::new();
    let mut emitted = 0;
    let mut quarantined = 0;
    for value in &assertions {
        *counts
            .entry(py_str(
                value
                    .get("exact_normalization_status")
                    .unwrap_or(&Json::Null),
            ))
            .or_default() += 1;
        if let Some(Json::Array(items)) = value.get("exact_normalized_variants") {
            emitted += items.len();
        }
        if let Some(Json::Array(items)) = value.get("quarantined_exact_normalized_variants") {
            quarantined += items.len();
        }
    }
    let digest = sha256_file(path)?;
    let mut integration = object([
        ("exact_cache_path", text(path.display().to_string())),
        ("exact_cache_sha256", text(digest.clone())),
    ]);
    if let Json::Object(fields) = stats {
        for (key, value) in fields {
            integration.set(key, value.clone());
        }
    }
    integration.set("public_assertions", Json::Int(assertions.len() as i128));
    integration.set(
        "emitted_exact_normalized_records",
        Json::Int(emitted as i128),
    );
    integration.set(
        "emitted_quarantined_records",
        Json::Int(quarantined as i128),
    );
    integration.set(
        "assertion_status_counts",
        Json::Object(counts.into_iter().map(|(k, v)| (k, Json::Int(v))).collect()),
    );
    let annotated = exact
        .values()
        .filter(|variant| variant.get("clinvar").is_some_and(truthy))
        .count();
    let meta = meta_mut(result);
    meta.set("gofcards_exact_integration", integration);
    meta.set(
        "clinvar_vcv_integration",
        object([
            ("source", text("inject_clinvar_into_gofcards.py")),
            ("gofcards_exact_sha256", text(digest)),
            (
                "variants_carrying_clinvar_conditions",
                Json::Int(annotated as i128),
            ),
        ]),
    );
    Ok(annotated)
}
