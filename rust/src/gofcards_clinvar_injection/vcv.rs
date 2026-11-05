//! `stream_parse_clinvar_vcv` and the helpers beneath it. A match keeps what
//! the injection and the audit read from the script's match object: its genes,
//! its GoFCards records, the variation's names and HGVS expressions, and the
//! condition assertions in full.
use super::{Lookup, VariantRow, error};
use crate::cli::Result;
use crate::gene_nonlof_mechanism::{object, text};
use crate::interpro_mapping_builder::{Element, visit_elements_by_local_name};
use crate::json_text::Json;
use crate::python_io::PythonGzipReader;
use crate::python_path::os_error;
use crate::python_text::{is_python_space, python_int, python_strip};
use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::rc::Rc;

type Node = Rc<RefCell<Element>>;

fn is_named(node: &Node, name: &str) -> bool {
    node.borrow().name.rsplit('}').next() == Some(name)
}

/// `_child` for an element that exists.
fn child(element: &Element, name: &str) -> Option<Node> {
    element
        .children
        .iter()
        .find(|item| is_named(item, name))
        .cloned()
}

/// `_child`.
fn child_of(node: Option<&Node>, name: &str) -> Option<Node> {
    child(&node?.borrow(), name)
}

/// `_children`.
fn children(node: Option<&Node>, name: &str) -> Vec<Node> {
    node.map_or_else(Vec::new, |node| {
        node.borrow()
            .children
            .iter()
            .filter(|item| is_named(item, name))
            .cloned()
            .collect()
    })
}

/// `_text`.
fn text_of(node: Option<&Node>) -> String {
    node.map_or_else(String::new, |node| {
        python_strip(&node.borrow().text).to_owned()
    })
}

fn collect_named(node: &Node, name: &str, found: &mut Vec<Node>) {
    if is_named(node, name) {
        found.push(node.clone());
    }
    for item in &node.borrow().children {
        collect_named(item, name, found);
    }
}

/// `_descendants`: the element itself comes first, as in `Element.iter()`.
fn descendants(node: Option<&Node>, name: &str) -> Vec<Node> {
    let mut found = Vec::new();
    if let Some(node) = node {
        collect_named(node, name, &mut found);
    }
    found
}

/// `element.get(key)`.
fn attribute(node: &Node, key: &str) -> Option<String> {
    node.borrow().attribute(key).map(str::to_owned)
}

/// `element.get(key, "")`.
fn attr(node: &Node, key: &str) -> String {
    attribute(node, key).unwrap_or_default()
}

/// `element.get(key, "") if element is not None else ""`.
fn attr_of(node: Option<&Node>, key: &str) -> String {
    node.map_or_else(String::new, |node| attr(node, key))
}

/// `_nonempty`.
fn nonempty<const N: usize>(values: [(&str, Json); N]) -> Json {
    Json::Object(
        values
            .into_iter()
            .filter(|(_, value)| match value {
                Json::Null => false,
                Json::Str(item) => !item.is_empty(),
                Json::Array(items) => !items.is_empty(),
                Json::Object(items) => !items.is_empty(),
                _ => true,
            })
            .map(|(key, value)| (key.to_owned(), value))
            .collect(),
    )
}

/// `_int_or_text`. A number too large for 128 bits stays text.
fn int_or_text(value: Option<&str>) -> Json {
    match value {
        None | Some("") => Json::Null,
        Some(value) => python_int(value).map_or_else(|| text(value), Json::Int),
    }
}

/// `re.sub(r"\s+", " ", value)`.
fn collapse_space(value: &str) -> String {
    let mut result = String::with_capacity(value.len());
    let mut in_space = false;
    for c in value.chars() {
        if is_python_space(c) {
            if !in_space {
                result.push(' ');
            }
            in_space = true;
        } else {
            result.push(c);
            in_space = false;
        }
    }
    result
}

fn folded(value: &str) -> String {
    python_strip(value).to_lowercase()
}

/// `review_stars`.
fn review_stars(status: &str) -> i128 {
    match collapse_space(&folded(status)).as_str() {
        "criteria provided, single submitter" | "criteria provided, conflicting classifications" => {
            1
        }
        "criteria provided, multiple submitters, no conflicts" => 2,
        "reviewed by expert panel" => 3,
        "practice guideline" => 4,
        _ => 0,
    }
}

/// `normalize_vcf_key`.
pub(super) fn normalize_vcf_key(chromosome: &str, position: &str, reference: &str, alternate: &str) -> String {
    let stripped = python_strip(chromosome);
    let mut chrom: String = if stripped.to_lowercase().starts_with("chr") {
        stripped.chars().skip(3).collect()
    } else {
        stripped.to_owned()
    };
    if chrom == "MT" {
        chrom = "M".to_owned();
    }
    [
        chrom.as_str(),
        python_strip(position),
        &python_strip(reference).to_uppercase(),
        &python_strip(alternate).to_uppercase(),
    ]
    .join("|")
}

/// The builds and keys of `_allele_locations`: one entry per assembly, the
/// last usable location holding the first one's place.
fn allele_keys(allele: &Node) -> Vec<(&'static str, &'static str, String)> {
    let mut keys: Vec<(&'static str, &'static str, String)> = Vec::new();
    let location = child(&allele.borrow(), "Location");
    for sequence in descendants(location.as_ref(), "SequenceLocation") {
        let assembly = attr(&sequence, "Assembly");
        let (assembly, build) = match python_strip(&assembly) {
            "GRCh37" => ("GRCh37", "hg19"),
            "GRCh38" => ("GRCh38", "hg38"),
            _ => continue,
        };
        let value = |key: &str| python_strip(&attr(&sequence, key)).to_owned();
        let (chrom, position) = (value("Chr"), value("positionVCF"));
        let (reference, alternate) = (value("referenceAlleleVCF"), value("alternateAlleleVCF"));
        if chrom.is_empty() || position.is_empty() || reference.is_empty() || alternate.is_empty() {
            continue;
        }
        let key = normalize_vcf_key(&chrom, &position, &reference, &alternate);
        match keys.iter_mut().find(|entry| entry.0 == assembly) {
            Some(entry) => entry.2 = key,
            None => keys.push((assembly, build, key)),
        }
    }
    keys
}

/// The symbols `_make_match` reads from `_allele_genes`.
fn allele_symbols(allele: &Node) -> BTreeSet<String> {
    let list = child(&allele.borrow(), "GeneList");
    children(list.as_ref(), "Gene")
        .iter()
        .map(|gene| attr(gene, "Symbol"))
        .filter(|symbol| !symbol.is_empty())
        .collect()
}

/// The expressions of `_allele_hgvs`.
fn allele_hgvs(allele: &Node) -> Vec<String> {
    let list = child(&allele.borrow(), "HGVSlist");
    children(list.as_ref(), "HGVS")
        .iter()
        .filter_map(|hgvs| {
            let hgvs = hgvs.borrow();
            let nucleotide = child(&hgvs, "NucleotideExpression");
            let protein = child(&hgvs, "ProteinExpression");
            let expression = text_of(child_of(nucleotide.as_ref(), "Expression").as_ref());
            let expression = if expression.is_empty() {
                text_of(child_of(protein.as_ref(), "Expression").as_ref())
            } else {
                expression
            };
            (!expression.is_empty()).then_some(expression)
        })
        .collect()
}

/// `_variant_alleles`.
fn variant_alleles(classified: &Node) -> Vec<(Node, &'static str, &'static str)> {
    let record = classified.borrow();
    if let Some(direct) = child(&record, "SimpleAllele") {
        return vec![(direct, "simple_allele", "allele")];
    }
    if let Some(haplotype) = child(&record, "Haplotype") {
        return children(Some(&haplotype), "SimpleAllele")
            .into_iter()
            .map(|allele| (allele, "haplotype", "haplotype_component"))
            .collect();
    }
    let Some(genotype) = child(&record, "Genotype") else {
        return vec![];
    };
    let mut rows = Vec::new();
    for item in &genotype.borrow().children {
        if is_named(item, "SimpleAllele") {
            rows.push((item.clone(), "genotype", "genotype_component"));
        } else if is_named(item, "Haplotype") {
            for allele in children(Some(item), "SimpleAllele") {
                rows.push((allele, "genotype", "genotype_haplotype_component"));
            }
        }
    }
    rows
}

fn add_identifier_tokens(tokens: &mut HashSet<String>, database: &str, identifier: &str) {
    let (database, identifier) = (folded(database), folded(identifier));
    if !identifier.is_empty() {
        if !database.is_empty() {
            tokens.insert(format!("{database}:{identifier}"));
        }
        tokens.insert(identifier);
    }
}

/// `_trait_tokens`.
fn trait_tokens(assertion: &Node) -> HashSet<String> {
    let mut tokens = HashSet::new();
    let trait_set = child(&assertion.borrow(), "TraitSet");
    for item in descendants(trait_set.as_ref(), "Trait") {
        for xref in children(Some(&item), "XRef") {
            add_identifier_tokens(&mut tokens, &attr(&xref, "DB"), &attr(&xref, "ID"));
        }
        let name = child(&item.borrow(), "Name");
        for preferred in descendants(name.as_ref(), "ElementValue") {
            let value = collapse_space(&text_of(Some(&preferred)).to_lowercase());
            if !value.is_empty() {
                tokens.insert(format!("name:{value}"));
            }
        }
    }
    tokens
}

/// One `_trait_mapping_by_assertion` row, with its `_mapping_tokens` input.
struct TraitMapping {
    row: Json,
    cui: String,
    name: String,
}

/// `_trait_mapping_by_assertion`.
fn trait_mappings(classified: &Node) -> HashMap<String, Vec<TraitMapping>> {
    let mut mapping: HashMap<String, Vec<TraitMapping>> = HashMap::new();
    let list = child(&classified.borrow(), "TraitMappingList");
    for item in children(list.as_ref(), "TraitMapping") {
        let assertion_id = attr(&item, "ClinicalAssertionID");
        if assertion_id.is_empty() {
            continue;
        }
        let medgen = child(&item.borrow(), "MedGen");
        let (name, cui) = (attr_of(medgen.as_ref(), "Name"), attr_of(medgen.as_ref(), "CUI"));
        let row = nonempty([
            ("clinical_assertion_id", text(assertion_id.clone())),
            ("trait_type", text(attr(&item, "TraitType"))),
            ("mapping_type", text(attr(&item, "MappingType"))),
            ("mapping_value", text(attr(&item, "MappingValue"))),
            ("mapping_ref", text(attr(&item, "MappingRef"))),
            ("medgen_name", text(name.clone())),
            ("medgen_cui", text(cui.clone())),
        ]);
        mapping
            .entry(assertion_id)
            .or_default()
            .push(TraitMapping { row, cui, name });
    }
    mapping
}

/// `_mapping_tokens`.
fn mapping_tokens(rows: &[TraitMapping]) -> HashSet<String> {
    let mut tokens = HashSet::new();
    for row in rows {
        let cui = folded(&row.cui);
        let name = collapse_space(&folded(&row.name));
        if !cui.is_empty() {
            tokens.insert(format!("medgen:{cui}"));
            tokens.insert(cui);
        }
        if !name.is_empty() {
            tokens.insert(format!("name:{name}"));
        }
    }
    tokens
}

/// `_attribute_values`.
fn attribute_values(assertion: &Node, wanted_type: &str) -> Json {
    let mut values = BTreeSet::new();
    for attribute_set in children(Some(assertion), "AttributeSet") {
        if let Some(item) = child(&attribute_set.borrow(), "Attribute") {
            if attribute(&item, "Type").as_deref() == Some(wanted_type) {
                let value = text_of(Some(&item));
                if !value.is_empty() {
                    values.insert(value);
                }
            }
        }
    }
    Json::Array(values.into_iter().map(Json::Str).collect())
}

fn text_field(node: Option<&Node>, name: &str) -> Json {
    text(text_of(child_of(node, name).as_ref()))
}

fn int_field(node: Option<&Node>, name: &str) -> Json {
    int_or_text(Some(&text_of(child_of(node, name).as_ref())))
}

/// `_cooccurrence`.
fn cooccurrence(item: &Node) -> Json {
    let alleles = children(Some(item), "AlleleDescSet")
        .iter()
        .map(|allele| {
            let significance = child(&allele.borrow(), "ClinicalSignificance");
            let description = text_of(child_of(significance.as_ref(), "Description").as_ref());
            let description = if description.is_empty() {
                text_of(significance.as_ref())
            } else {
                description
            };
            nonempty([
                ("name", text_field(Some(allele), "Name")),
                (
                    "relative_orientation",
                    text_field(Some(allele), "RelativeOrientation"),
                ),
                ("zygosity", text_field(Some(allele), "Zygosity")),
                ("clinical_significance", text(description)),
            ])
        })
        .collect();
    nonempty([
        ("asserted_variant_zygosity", text_field(Some(item), "Zygosity")),
        ("count", int_field(Some(item), "Count")),
        ("cooccurring_alleles", Json::Array(alleles)),
    ])
}

/// `_observation`.
fn observation(item: &Node) -> Json {
    let sample = child(&item.borrow(), "Sample");
    let sample = sample.as_ref();
    let data = children(Some(item), "ObservedData")
        .iter()
        .filter_map(|observed| child(&observed.borrow(), "Attribute"))
        .map(|value| {
            nonempty([
                ("type", text(attr(&value, "Type"))),
                ("value", text(text_of(Some(&value)))),
            ])
        })
        .collect();
    let sets = children(Some(item), "Co-occurrenceSet")
        .iter()
        .map(cooccurrence)
        .collect();
    nonempty([
        ("origin", text_field(sample, "Origin")),
        ("affected_status", text_field(sample, "AffectedStatus")),
        ("sex", text_field(sample, "Sex")),
        ("number_tested", int_field(sample, "NumberTested")),
        ("number_males", int_field(sample, "NumberMales")),
        ("number_females", int_field(sample, "NumberFemales")),
        ("number_chromosomes_tested", int_field(sample, "NumberChrTested")),
        ("observed_data", Json::Array(data)),
        ("cooccurrence_sets", Json::Array(sets)),
    ])
}

/// `_clinical_assertion`.
fn clinical_assertion(assertion: &Node, mappings: &[TraitMapping], linkage_method: &str) -> Json {
    let (accession, classification, observed_in) = {
        let item = assertion.borrow();
        (
            child(&item, "ClinVarAccession"),
            child(&item, "Classification"),
            child(&item, "ObservedInList"),
        )
    };
    let (accession, classification) = (accession.as_ref(), classification.as_ref());
    let germline = child_of(classification, "GermlineClassification");
    let observations: Vec<Json> = children(observed_in.as_ref(), "ObservedIn")
        .iter()
        .map(observation)
        .collect();
    let mut zygosity_counts: BTreeMap<String, i128> = BTreeMap::new();
    for item in &observations {
        let Some(Json::Array(sets)) = item.get("cooccurrence_sets") else {
            continue;
        };
        for set in sets {
            if let Some(Json::Str(value)) = set.get("asserted_variant_zygosity") {
                let count = match set.get("count") {
                    Some(Json::Int(number)) => *number,
                    _ => 1,
                };
                *zygosity_counts.entry(value.clone()).or_default() += count;
            }
        }
    }
    let contributes = folded_exact(&attr(assertion, "ContributesToAggregateClassification"));
    nonempty([
        ("clinical_assertion_id", text(attr(assertion, "ID"))),
        ("scv_accession", text(attr_of(accession, "Accession"))),
        (
            "scv_version",
            accession.map_or(Json::Null, |item| {
                int_or_text(attribute(item, "Version").as_deref())
            }),
        ),
        ("record_status", text_field(Some(assertion), "RecordStatus")),
        (
            "contributes_to_aggregate_classification",
            Json::Bool(contributes == "true"),
        ),
        (
            "submitter",
            nonempty([
                ("name", text(attr_of(accession, "SubmitterName"))),
                ("organization_id", text(attr_of(accession, "OrgID"))),
                (
                    "organization_category",
                    text(attr_of(accession, "OrganizationCategory")),
                ),
            ]),
        ),
        (
            "classification",
            nonempty([
                ("germline_significance", text(text_of(germline.as_ref()))),
                ("review_status", text_field(classification, "ReviewStatus")),
                (
                    "date_last_evaluated",
                    text(attr_of(classification, "DateLastEvaluated")),
                ),
            ]),
        ),
        (
            "submitted_mode_of_inheritance",
            attribute_values(assertion, "ModeOfInheritance"),
        ),
        ("penetrance", attribute_values(assertion, "Penetrance")),
        ("trait_linkage_method", text(linkage_method)),
        (
            "trait_mappings",
            Json::Array(mappings.iter().map(|item| item.row.clone()).collect()),
        ),
        (
            "observed_zygosity_counts",
            Json::Object(
                zygosity_counts
                    .into_iter()
                    .map(|(key, count)| (key, Json::Int(count)))
                    .collect(),
            ),
        ),
        ("observations", Json::Array(observations)),
    ])
}

/// `value.lower()`.
fn folded_exact(value: &str) -> String {
    value.to_lowercase()
}

/// A dict of `_eligible_conditions` before `_attach_scvs` finishes it.
struct Condition {
    fields: Json,
    rcv_accession: String,
    tokens: HashSet<String>,
    scvs: Vec<Json>,
}

/// `_eligible_conditions`.
fn eligible_conditions(classified: &Node, min_review_stars: i128) -> Vec<Condition> {
    let mut eligible = Vec::new();
    let list = child(&classified.borrow(), "RCVList");
    for rcv in children(list.as_ref(), "RCVAccession") {
        let classifications = child(&rcv.borrow(), "RCVClassifications");
        let Some(germline) = child_of(classifications.as_ref(), "GermlineClassification") else {
            continue;
        };
        let status = text_of(child(&germline.borrow(), "ReviewStatus").as_ref());
        let stars = review_stars(&status);
        if stars < min_review_stars {
            continue;
        }
        let description = child(&germline.borrow(), "Description");
        let condition_list = child(&rcv.borrow(), "ClassifiedConditionList");
        let mut tokens = HashSet::new();
        let mut conditions = Vec::new();
        for condition in children(condition_list.as_ref(), "ClassifiedCondition") {
            let (name, database) = (text_of(Some(&condition)), attr(&condition, "DB"));
            let identifier = attr(&condition, "ID");
            add_identifier_tokens(&mut tokens, &database, &identifier);
            let token = collapse_space(&folded(&name));
            if !token.is_empty() {
                tokens.insert(format!("name:{token}"));
            }
            conditions.push(nonempty([
                ("name", text(name)),
                ("database", text(database)),
                ("id", text(identifier)),
            ]));
        }
        let rcv_accession = attr(&rcv, "Accession");
        let fields = object([
            ("rcv_accession", text(rcv_accession.clone())),
            ("rcv_version", int_or_text(attribute(&rcv, "Version").as_deref())),
            ("title", text(attr(&rcv, "Title"))),
            ("trait_set_id", text(attr_of(condition_list.as_ref(), "TraitSetID"))),
            ("conditions", Json::Array(conditions)),
            (
                "germline_classification",
                nonempty([
                    ("clinical_significance", text(text_of(description.as_ref()))),
                    ("review_status", text(status)),
                    ("review_stars", Json::Int(stars)),
                    (
                        "date_last_evaluated",
                        text(attr_of(description.as_ref(), "DateLastEvaluated")),
                    ),
                    (
                        "submission_count",
                        int_or_text(
                            description
                                .as_ref()
                                .and_then(|item| attribute(item, "SubmissionCount"))
                                .as_deref(),
                        ),
                    ),
                ]),
            ),
        ]);
        eligible.push(Condition {
            fields,
            rcv_accession,
            tokens,
            scvs: vec![],
        });
    }
    eligible
}

/// `str(scv.get("clinical_assertion_id") or scv.get("scv_accession") or "")`.
fn scv_key(scv: &Json) -> String {
    ["clinical_assertion_id", "scv_accession"]
        .iter()
        .find_map(|key| match scv.get(key) {
            Some(Json::Str(value)) if !value.is_empty() => Some(value.clone()),
            _ => None,
        })
        .unwrap_or_default()
}

/// `_attach_scvs`, returning the finished condition dicts.
fn attach_scvs(classified: &Node, mut conditions: Vec<Condition>) -> Vec<Json> {
    let mapping_by_id = trait_mappings(classified);
    let list = child(&classified.borrow(), "ClinicalAssertionList");
    let assertions = children(list.as_ref(), "ClinicalAssertion");
    let mut linkage: Vec<(i128, i128)> = Vec::new();
    for condition in &mut conditions {
        let mut unlinked = 0;
        for assertion in &assertions {
            if text_of(child(&assertion.borrow(), "RecordStatus").as_ref()) != "current" {
                continue;
            }
            if folded_exact(&attr(assertion, "ContributesToAggregateClassification")) != "true" {
                continue;
            }
            let classification = child(&assertion.borrow(), "Classification");
            if child_of(classification.as_ref(), "GermlineClassification").is_none() {
                continue;
            }
            let mapping_rows = mapping_by_id
                .get(&attr(assertion, "ID"))
                .map_or(&[][..], Vec::as_slice);
            let method = if !condition.tokens.is_disjoint(&mapping_tokens(mapping_rows)) {
                "trait_mapping_clinical_assertion_id"
            } else if !condition.tokens.is_disjoint(&trait_tokens(assertion)) {
                "scv_trait_identifier_or_name"
            } else {
                unlinked += 1;
                continue;
            };
            condition
                .scvs
                .push(clinical_assertion(assertion, mapping_rows, method));
        }
        linkage.push((condition.scvs.len() as i128, unlinked));
    }
    let mut scv_to_rcvs: HashMap<String, BTreeSet<String>> = HashMap::new();
    for condition in &conditions {
        for scv in &condition.scvs {
            let key = scv_key(scv);
            if !key.is_empty() && !condition.rcv_accession.is_empty() {
                scv_to_rcvs
                    .entry(key)
                    .or_default()
                    .insert(condition.rcv_accession.clone());
            }
        }
    }
    conditions
        .into_iter()
        .zip(linkage)
        .map(|(mut condition, (matched, unlinked))| {
            for scv in &mut condition.scvs {
                let linked = scv_to_rcvs.get(&scv_key(scv));
                if let Some(linked) = linked.filter(|linked| linked.len() > 1) {
                    scv.set("trait_linkage_ambiguous_across_eligible_rcvs", Json::Bool(true));
                    scv.set(
                        "linked_eligible_rcv_accessions",
                        Json::Array(linked.iter().cloned().map(Json::Str).collect()),
                    );
                }
            }
            condition.fields.set("matched_scvs", Json::Array(condition.scvs));
            condition.fields.set(
                "scv_linkage",
                object([
                    (
                        "method",
                        text("condition identifiers/names plus TraitMapping ClinicalAssertionID"),
                    ),
                    (
                        "matched_current_contributing_germline_scvs",
                        Json::Int(matched),
                    ),
                    (
                        "unlinked_current_contributing_germline_scvs",
                        Json::Int(unlinked),
                    ),
                ]),
            );
            condition.fields
        })
        .collect()
}

/// A GoFCards row of one allele: symbol, variant identifier and route.
type Row = (String, String, &'static str);

/// `_match_rows_for_allele`.
fn match_rows(allele: &Node, exact_lookup: &Lookup) -> Vec<Row> {
    let mut rows: Vec<Row> = Vec::new();
    let mut seen: HashSet<(String, String)> = HashSet::new();
    for (_, build, key) in allele_keys(allele) {
        let Some(found) = exact_lookup.get(&(build.to_owned(), key)) else {
            continue;
        };
        for (symbol, variant_id) in found {
            if seen.insert((symbol.clone(), variant_id.clone())) {
                rows.push((symbol.clone(), variant_id.clone(), "genomic_coordinates"));
            }
        }
    }
    rows
}

/// What the injection and the audit read from one `_make_match` result.
pub(super) struct Match {
    pub(super) gene_symbols: Vec<String>,
    /// The `gofcards_variant_id` and `link_route` of each retained record.
    pub(super) records: Vec<(String, &'static str)>,
    pub(super) vcv_accession: String,
    pub(super) variation_id: String,
    pub(super) variation_name: String,
    pub(super) classification_scope: &'static str,
    pub(super) component_context: &'static str,
    pub(super) hgvs: Vec<String>,
    pub(super) conditions: Vec<Json>,
}

/// `_make_match`.
fn make_match(
    archive: &Element,
    allele: &Node,
    scope: &'static str,
    context: &'static str,
    rows: Vec<Row>,
    conditions: &[Json],
) -> Match {
    let clinvar_symbols = allele_symbols(allele);
    let matched_symbols: BTreeSet<String> = rows
        .iter()
        .filter(|row| !row.0.is_empty())
        .map(|row| row.0.clone())
        .collect();
    let concordant: Vec<String> = clinvar_symbols
        .intersection(&matched_symbols)
        .cloned()
        .collect();
    let by_clinvar_gene = !clinvar_symbols.is_empty();
    let gene_symbols = if by_clinvar_gene {
        concordant.clone()
    } else {
        matched_symbols.into_iter().collect()
    };
    let records = rows
        .into_iter()
        .filter(|row| !by_clinvar_gene || concordant.contains(&row.0))
        .map(|(_, variant_id, route)| (variant_id, route))
        .collect();
    let value = |key: &str| archive.attribute(key).unwrap_or_default().to_owned();
    Match {
        gene_symbols,
        records,
        vcv_accession: value("Accession"),
        variation_id: value("VariationID"),
        variation_name: value("VariationName"),
        classification_scope: scope,
        component_context: context,
        hgvs: allele_hgvs(allele),
        conditions: conditions.to_vec(),
    }
}

/// The loop of `stream_parse_clinvar_vcv`.
struct Scan {
    exact_lookup: Lookup,
    variation_id_lookup: HashMap<String, Vec<VariantRow>>,
    min_review_stars: i128,
    matches: Vec<Match>,
    stats: BTreeMap<String, i128>,
}

impl Scan {
    fn add(&mut self, key: &str, count: i128) {
        match self.stats.get_mut(key) {
            Some(value) => *value += count,
            None => {
                self.stats.insert(key.to_owned(), count);
            }
        }
    }

    fn count(&self, key: &str) -> i128 {
        self.stats.get(key).copied().unwrap_or(0)
    }

    fn visit(&mut self, archive: &Element) -> Result<()> {
        self.add("variation_archives_seen", 1);
        let seen = self.count("variation_archives_seen");
        if seen % 100_000 == 0 {
            eprintln!(
                "clinvar_vcv_progress: variation_archives={seen}, exact_archives={}, retained_matches={}",
                self.count("exact_variant_archives"),
                self.count("retained_component_matches"),
            );
        }
        let status = text_of(child(archive, "RecordStatus").as_ref());
        let classified = child(archive, "ClassifiedRecord");
        if status != "current" {
            self.add("skipped_noncurrent", 1);
            return Ok(());
        }
        let Some(classified) = classified else {
            self.add("skipped_included_record", 1);
            return Ok(());
        };
        let identifier = archive.attribute("VariationID").unwrap_or_default();
        let by_variation_id: &[VariantRow] = self
            .variation_id_lookup
            .get(python_strip(identifier))
            .map_or(&[], Vec::as_slice);
        let named = !by_variation_id.is_empty();
        let mut candidates = Vec::new();
        for (allele, scope, context) in variant_alleles(&classified) {
            let mut rows = match_rows(&allele, &self.exact_lookup);
            if named {
                let known: HashSet<(String, String)> = rows
                    .iter()
                    .map(|row| (row.0.clone(), row.1.clone()))
                    .collect();
                for (symbol, variant_id) in by_variation_id {
                    if !known.contains(&(symbol.clone(), variant_id.clone())) {
                        rows.push((symbol.clone(), variant_id.clone(), "gofcards_variation_id"));
                    }
                }
            }
            if !rows.is_empty() {
                candidates.push((allele, scope, context, rows));
            }
        }
        if named {
            self.add("archives_matched_by_gofcards_variation_id", 1);
        }
        if candidates.is_empty() {
            return Ok(());
        }
        self.add("exact_variant_archives", 1);
        let conditions = eligible_conditions(&classified, self.min_review_stars);
        if conditions.is_empty() {
            self.add("exact_archives_without_eligible_rcv", 1);
            return Ok(());
        }
        let conditions = attach_scvs(&classified, conditions);
        for (allele, scope, context, rows) in candidates {
            let found = make_match(archive, &allele, scope, context, rows, &conditions);
            if found.gene_symbols.is_empty() {
                self.add("skipped_gene_discordant_component_matches", 1);
            } else {
                self.matches.push(found);
                self.add("retained_component_matches", 1);
                self.add(&format!("retained_{scope}"), 1);
            }
        }
        self.add("eligible_rcv_conditions", conditions.len() as i128);
        Ok(())
    }
}

/// `stream_parse_clinvar_vcv`: the matches and the sorted statistics.
pub(super) fn stream_parse_clinvar_vcv(
    xml_path: &Path,
    exact_lookup: Lookup,
    variation_id_lookup: HashMap<String, Vec<VariantRow>>,
    min_review_stars: i128,
) -> Result<(Vec<Match>, BTreeMap<String, i128>)> {
    let file = File::open(xml_path).map_err(|e| os_error(e, xml_path, None))?;
    let source: Box<dyn Read> = if xml_path.extension().is_some_and(|suffix| suffix == "gz") {
        Box::new(PythonGzipReader::new(file))
    } else {
        Box::new(file)
    };
    let scan = Rc::new(RefCell::new(Scan {
        exact_lookup,
        variation_id_lookup,
        min_review_stars,
        matches: vec![],
        stats: BTreeMap::new(),
    }));
    let visiting = scan.clone();
    visit_elements_by_local_name(source, "VariationArchive", move |archive| {
        visiting.borrow_mut().visit(archive)
    })?;
    let Ok(scan) = Rc::try_unwrap(scan) else {
        return Err(error("RuntimeError", "the XML reader kept its visitor"));
    };
    let mut scan = scan.into_inner();
    let assignments = scan
        .matches
        .iter()
        .map(|found| found.gene_symbols.len() as i128)
        .sum();
    scan.stats
        .insert("retained_gene_assignments".to_owned(), assignments);
    Ok((scan.matches, scan.stats))
}
