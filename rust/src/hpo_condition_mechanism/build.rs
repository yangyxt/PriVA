//! The build, in the script's order: the HPO frame and its MONDO merge, the
//! G2P/Orphadata mechanisms and penetrance, ClinGen haploinsufficiency, the
//! GoFCards exact alleles, the deduced recessive LOF and the scope gate.
use super::tsv::{Reader, Row, cell};
use super::*;
use crate::acmg_step::gene_mechanism_conditions::normalize_inheritance;
use crate::hpo_penetrance::{Values, normalize_penetrance_values, penetrance_status};
use std::collections::{BTreeMap, BTreeSet, HashSet};

const SCHEMA_VERSION: &str = "1.2";

const HPO_INHERITANCE_TERMS: &[(&str, &str)] = &[
    ("HP:0000006", "autosomal_dominant"),
    ("HP:0000007", "autosomal_recessive"),
    ("HP:0001417", "x_linked"),
    ("HP:0001419", "x_linked_recessive"),
    ("HP:0001423", "x_linked_dominant"),
    ("HP:0001426", "non_mendelian"),
    ("HP:0001427", "mitochondrial"),
    ("HP:0001450", "y_linked"),
    ("HP:0010982", "polygenic"),
    ("HP:0010983", "oligogenic"),
    ("HP:0010984", "digenic"),
    ("HP:0012275", "autosomal_dominant_maternal_imprinting"),
    ("HP:0034341", "pseudoautosomal_recessive"),
];

const HPO_ONSET_TERMS: &[(&str, &str)] = &[
    ("HP:0030674", "antenatal"),
    ("HP:0011460", "embryonal"),
    ("HP:0011461", "fetal"),
    ("HP:0003577", "congenital"),
    ("HP:0003623", "neonatal"),
    ("HP:0003593", "infantile"),
    ("HP:0011463", "childhood"),
    ("HP:0003621", "juvenile"),
    ("HP:0410280", "pediatric"),
    ("HP:0003581", "adult"),
    ("HP:0011462", "young_adult"),
    ("HP:0003596", "middle_age"),
    ("HP:0003584", "late"),
    ("HP:0034857", "variable_age"),
    ("HP:0003587", "insidious"),
];

const HPO_REQUIRED_COLUMNS: &[&str] = &[
    "gene_symbol",
    "disease_id",
    "hpo_id",
    "frequency",
    "evidence",
    "reference",
    "mondo_id",
    "mondo_name",
    "disease_scope",
    "priva_scope",
    "scope_evidence",
    "scope_reference",
    "scope_review_status",
];

const CONDITION_MECHANISM_SOURCES: &[&str] = &["G2P_DDG2P", "Orphadata"];
pub(super) const CANONICAL_MECHANISMS: &[&str] = &["LOF", "GOF", "DOMINANT_NEGATIVE"];
const MECHANISM_REQUIRED_COLUMNS: &[&str] = &[
    "gene_symbol",
    "source",
    "source_record_id",
    "source_condition_id",
    "mondo_id",
    "disease_scope",
    "priva_scope",
    "scope_review_status",
    "disease_label",
    "inheritance",
    "penetrance_raw",
    "penetrance_hpo_ids",
    "normalized_penetrance",
    "patho_mode_raw",
    "normalized_mechanisms",
    "mechanism_confidence",
    "disease_confidence",
    "pmids",
    "evidence_url",
];
const CLINGEN_REQUIRED_COLUMNS: &[&str] = &[
    "#Gene Symbol",
    "Genomic Location",
    "Haploinsufficiency Score",
    "Haploinsufficiency Description",
    "Haploinsufficiency Disease ID",
];

/// Most restrictive first; a merged condition takes the most cautious decision.
const SCOPE_PRECEDENCE: [&str; 4] = ["exclude", "review", "include", ""];

fn lookup(table: &[(&str, &'static str)], term: &str) -> Option<&'static str> {
    table.iter().find(|(id, _)| *id == term).map(|(_, v)| *v)
}

fn upper(text: &str) -> String {
    text.to_uppercase()
}

fn strings(values: &[String]) -> Json {
    Json::Array(values.iter().map(|v| Json::Str(v.clone())).collect())
}

fn obj(fields: Vec<(&str, Json)>) -> Json {
    Json::Object(fields.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}

fn st(value: &str) -> Json {
    Json::Str(value.to_owned())
}

/// A text field of an evidence record this builder wrote.
fn field<'a>(record: &'a Json, name: &str) -> &'a str {
    record.get(name).and_then(Json::as_str).unwrap_or("")
}

/// `";".join(record.get(name, []) or [])` for a list of texts.
fn joined(record: &Json, name: &str) -> String {
    match record.get(name) {
        Some(Json::Array(items)) => items
            .iter()
            .map(|item| item.as_str().unwrap_or("").to_owned())
            .collect::<Vec<_>>()
            .join(";"),
        _ => String::new(),
    }
}

fn push_unique(target: &mut Vec<String>, value: &str) {
    if !target.iter().any(|v| v == value) {
        target.push(value.to_owned());
    }
}

/// `_split_multi`.
fn split_multi(value: &str) -> Vec<String> {
    let mut parts = Vec::new();
    for part in clean(value).split(';') {
        let part = crate::hpo_penetrance::py_strip(part);
        if !part.is_empty() && part != "-" {
            push_unique(&mut parts, part);
        }
    }
    parts
}

/// `_penetrance_evidence_key`.
pub(super) fn penetrance_key(evidence: &Json) -> Vec<String> {
    vec![
        clean(field(evidence, "source")),
        clean(field(evidence, "source_record_id")),
        joined(evidence, "condition_identifiers"),
        joined(evidence, "raw_values"),
        clean(field(evidence, "normalized_value")),
        clean(field(evidence, "reference")),
    ]
}

/// `_mechanism_evidence_key`.
fn mechanism_key(evidence: &Json) -> Vec<String> {
    vec![
        field(evidence, "source").to_owned(),
        field(evidence, "source_record_id").to_owned(),
        joined(evidence, "condition_identifiers"),
        field(evidence, "mechanism").to_owned(),
        field(evidence, "allelic_requirement").to_owned(),
    ]
}

fn hpo_assertion_key(evidence: &Json) -> Vec<String> {
    ["hpo_id", "reference", "evidence", "frequency"]
        .iter()
        .map(|name| field(evidence, name).to_owned())
        .collect()
}

#[derive(Default)]
struct Axis {
    values: Vec<String>,
    assertions: Vec<Json>,
    keys: HashSet<Vec<String>>,
}

struct Scope {
    decision: String,
    category: String,
    review_status: String,
    review_statuses: Vec<String>,
    evidence: Vec<String>,
    references: Vec<String>,
}

#[derive(Default)]
struct Block {
    allelic_requirements: Vec<String>,
    evidence: Vec<Json>,
    variants: BTreeMap<String, Json>,
}

struct Condition {
    label: String,
    /// Namespace to identifiers, each list in insertion order until sorted.
    identifiers: BTreeMap<String, Vec<String>>,
    scope: Scope,
    inheritance: Axis,
    penetrance: Axis,
    onset: Axis,
    mechanisms: BTreeMap<String, Block>,
    hpo_assertion_count: i128,
}

#[derive(Default)]
struct Gene {
    conditions: BTreeMap<String, Condition>,
    unmapped_mechanisms: Vec<Json>,
    unmapped_penetrance: Vec<Json>,
    unmapped_variants: BTreeMap<String, Json>,
    excluded: i128,
}

type Genes = BTreeMap<String, Gene>;

/// `_add_identifier`.
fn add_identifier(target: &mut BTreeMap<String, Vec<String>>, identifier: &str) {
    let cleaned = upper(&clean(identifier));
    let Some((namespace, _)) = cleaned.split_once(':') else {
        return;
    };
    push_unique(target.entry(namespace.to_owned()).or_default(), &cleaned);
}

fn new_condition(row: &Row) -> Condition {
    let mut identifiers = BTreeMap::new();
    add_identifier(&mut identifiers, &cell(row, "disease_id"));
    add_identifier(&mut identifiers, &cell(row, "mondo_id"));
    Condition {
        label: cell(row, "mondo_name"),
        identifiers,
        scope: Scope {
            decision: cell(row, "priva_scope"),
            category: cell(row, "disease_scope"),
            review_status: cell(row, "scope_review_status"),
            review_statuses: vec![],
            evidence: vec![],
            references: vec![],
        },
        inheritance: Axis::default(),
        penetrance: Axis::default(),
        onset: Axis::default(),
        mechanisms: BTreeMap::new(),
        hpo_assertion_count: 0,
    }
}

/// `_merge_scope`.
fn merge_scope(scope: &mut Scope, row: &Row) -> Result<()> {
    for (target, source) in [("decision", "priva_scope"), ("category", "disease_scope")] {
        let incoming = cell(row, source);
        let slot = if target == "decision" {
            &mut scope.decision
        } else {
            &mut scope.category
        };
        let current = clean(slot);
        if !current.is_empty() && !incoming.is_empty() && current != incoming {
            return Err(error(
                "ValueError",
                format!(
                    "Conflicting HPO condition scope after canonical identity merge: \
{target}={} versus {}",
                    crate::python_path::repr(&current),
                    crate::python_path::repr(&incoming)
                ),
            ));
        }
        if !incoming.is_empty() {
            *slot = incoming;
        }
    }
    let review_status = cell(row, "scope_review_status");
    if !review_status.is_empty() {
        push_unique(&mut scope.review_statuses, &review_status);
    }
    if scope.review_statuses.iter().any(|s| s == "manually_confirmed") {
        scope.review_status = "manually_confirmed".into();
    } else if scope.review_statuses.len() == 1 {
        scope.review_status = scope.review_statuses[0].clone();
    } else if scope.review_statuses.len() > 1 {
        scope.review_status = "mixed".into();
    }
    for value in split_multi(&cell(row, "scope_evidence")) {
        push_unique(&mut scope.evidence, &value);
    }
    for value in split_multi(&cell(row, "scope_reference")) {
        push_unique(&mut scope.references, &value);
    }
    Ok(())
}

/// `_hpo_penetrance_evidence_record`.
fn hpo_penetrance_record(row: &Row, normalized_value: &str) -> Json {
    let hpo_id = upper(&cell(row, "hpo_id"));
    let disease = upper(&cell(row, "disease_id"));
    let mut identifiers = Vec::new();
    for value in [disease.clone(), upper(&cell(row, "mondo_id"))] {
        if !value.is_empty() {
            push_unique(&mut identifiers, &value);
        }
    }
    let reference = cell(row, "reference");
    let record_id = [disease, hpo_id.clone(), reference.clone()]
        .into_iter()
        .filter(|v| !v.is_empty() && v != "-")
        .collect::<Vec<_>>()
        .join("|");
    let ids = if hpo_id.is_empty() { vec![] } else { vec![hpo_id] };
    obj(vec![
        ("source", st("HPO")),
        ("source_record_id", st(&record_id)),
        ("condition_identifiers", strings(&identifiers)),
        ("raw_values", strings(&ids)),
        ("normalized_value", st(normalized_value)),
        ("hpo_ids", strings(&ids)),
        ("frequency", st(&or_dash(cell(row, "frequency")))),
        ("evidence", st(&or_dash(cell(row, "evidence")))),
        ("reference", st(&or_dash(reference))),
        ("evidence_url", st("")),
    ])
}

fn or_dash(value: String) -> String {
    if value.is_empty() { "-".into() } else { value }
}

/// `_evidence_record`.
fn hpo_evidence_record(row: &Row) -> Json {
    obj(vec![
        ("hpo_id", st(&upper(&cell(row, "hpo_id")))),
        ("frequency", st(&or_dash(cell(row, "frequency")))),
        ("evidence", st(&or_dash(cell(row, "evidence")))),
        ("reference", st(&or_dash(cell(row, "reference")))),
    ])
}

/// `_append_hpo_axis_assertion`.
fn append_axis(axis: &mut Axis, penetrance: bool, value: &str, row: &Row) {
    push_unique(&mut axis.values, value);
    let (evidence, key) = if penetrance {
        let evidence = hpo_penetrance_record(row, value);
        let key = penetrance_key(&evidence);
        (evidence, key)
    } else {
        let evidence = hpo_evidence_record(row);
        let key = ["hpo_id", "frequency", "evidence", "reference"]
            .iter()
            .map(|n| field(&evidence, n).to_owned())
            .collect();
        (evidence, key)
    };
    if axis.keys.insert(key) {
        axis.assertions.push(evidence);
    }
}

fn scope_rank(decision: &str) -> Result<usize> {
    SCOPE_PRECEDENCE
        .iter()
        .position(|v| *v == decision)
        .ok_or_else(|| {
            error(
                "ValueError",
                format!("{} is not in tuple", crate::python_path::repr(decision)),
            )
        })
}

fn union_axis(target: &mut Axis, source: Axis) {
    for value in source.values {
        push_unique(&mut target.values, &value);
    }
    for assertion in source.assertions {
        if !target.assertions.iter().any(|a| equal(a, &assertion)) {
            target.assertions.push(assertion);
        }
    }
}

/// `_merge_conditions_by_disease`.
fn merge_conditions(conditions: BTreeMap<String, Condition>) -> Result<BTreeMap<String, Condition>> {
    let mut groups: BTreeMap<String, Vec<(String, Condition)>> = BTreeMap::new();
    for (key, condition) in conditions {
        let group = condition
            .identifiers
            .get("MONDO")
            .and_then(|mondo| mondo.first())
            .map_or_else(|| key.clone(), |first| clean(first));
        groups.entry(group).or_default().push((key, condition));
    }
    let mut merged = BTreeMap::new();
    for (_, mut members) in groups {
        if members.len() == 1 {
            let (key, condition) = members.pop().unwrap();
            merged.insert(key, condition);
            continue;
        }
        let rank = |key: &str| -> (usize, String) {
            let prefix = upper(key.split(':').next().unwrap_or(""));
            let order = match prefix.as_str() {
                "OMIM" => 0,
                "ORPHA" => 1,
                "MONDO" => 2,
                _ => 3,
            };
            (order, key.to_owned())
        };
        members.sort_by(|a, b| rank(&a.0).cmp(&rank(&b.0)));
        let mut members = members.into_iter();
        let (primary_key, mut primary) = members.next().unwrap();
        for (_, other) in members {
            for (namespace, values) in &other.identifiers {
                let target = primary.identifiers.entry(namespace.clone()).or_default();
                for value in values {
                    push_unique(target, value);
                }
            }
            union_axis(&mut primary.inheritance, other.inheritance);
            union_axis(&mut primary.penetrance, other.penetrance);
            union_axis(&mut primary.onset, other.onset);
            for (mechanism, block) in other.mechanisms {
                match primary.mechanisms.get_mut(&mechanism) {
                    None => {
                        primary.mechanisms.insert(mechanism, block);
                    }
                    Some(existing) => {
                        for requirement in &block.allelic_requirements {
                            push_unique(&mut existing.allelic_requirements, requirement);
                        }
                        existing.evidence.extend(block.evidence);
                        existing.variants.extend(block.variants);
                    }
                }
            }
            primary.hpo_assertion_count += other.hpo_assertion_count;
            if clean(&primary.label).is_empty() {
                primary.label = other.label.clone();
            }
            let other_decision = clean(&other.scope.decision);
            if scope_rank(&other_decision)? < scope_rank(&clean(&primary.scope.decision))? {
                primary.scope.decision = other_decision;
                primary.scope.category = clean(&other.scope.category);
                primary.scope.review_status = clean(&other.scope.review_status);
            }
        }
        merged.insert(primary_key, primary);
    }
    Ok(merged)
}

/// `build_hpo_gene_condition_frame`.
fn hpo_frame(path: &Path) -> Result<Genes> {
    let mut reader = Reader::open(path)?;
    reader.require(path, HPO_REQUIRED_COLUMNS)?;
    let mut genes = Genes::new();
    while let Some(row) = reader.next_row()? {
        let symbol = upper(&cell(&row, "gene_symbol"));
        let key = {
            let disease = upper(&cell(&row, "disease_id"));
            if disease.is_empty() {
                upper(&cell(&row, "mondo_id"))
            } else {
                disease
            }
        };
        if symbol.is_empty() || key.is_empty() {
            continue;
        }
        let gene = genes.entry(symbol).or_default();
        let condition = gene
            .conditions
            .entry(key)
            .or_insert_with(|| new_condition(&row));
        let name = cell(&row, "mondo_name");
        if condition.label.is_empty() && !name.is_empty() {
            condition.label = name;
        }
        add_identifier(&mut condition.identifiers, &cell(&row, "disease_id"));
        add_identifier(&mut condition.identifiers, &cell(&row, "mondo_id"));
        merge_scope(&mut condition.scope, &row)?;
        condition.hpo_assertion_count += 1;
        let hpo_id = upper(&cell(&row, "hpo_id"));
        if let Some(mode) = lookup(HPO_INHERITANCE_TERMS, &hpo_id) {
            append_axis(&mut condition.inheritance, false, mode, &row);
        }
        if let Some(status) = penetrance_status(&hpo_id) {
            append_axis(&mut condition.penetrance, true, status, &row);
        }
        if let Some(term) = lookup(HPO_ONSET_TERMS, &hpo_id) {
            append_axis(&mut condition.onset, false, term, &row);
        }
    }
    for gene in genes.values_mut() {
        gene.conditions = merge_conditions(std::mem::take(&mut gene.conditions))?;
        for condition in gene.conditions.values_mut() {
            for values in condition.identifiers.values_mut() {
                values.sort();
            }
            for (axis, penetrance) in [
                (&mut condition.inheritance, false),
                (&mut condition.penetrance, true),
                (&mut condition.onset, false),
            ] {
                axis.keys.clear();
                axis.values.sort();
                if penetrance {
                    axis.assertions.sort_by_key(penetrance_key);
                } else {
                    axis.assertions.sort_by_key(hpo_assertion_key);
                }
            }
        }
    }
    Ok(genes)
}

/// `_condition_alias_index`: only identifiers naming one condition per gene.
fn alias_index(genes: &Genes) -> BTreeMap<String, BTreeMap<String, String>> {
    let mut index = BTreeMap::new();
    for (symbol, gene) in genes {
        let mut aliases: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
        for (key, condition) in &gene.conditions {
            aliases.entry(upper(key)).or_default().insert(key.clone());
            for identifiers in condition.identifiers.values() {
                for identifier in identifiers {
                    aliases
                        .entry(upper(identifier))
                        .or_default()
                        .insert(key.clone());
                }
            }
        }
        index.insert(
            symbol.clone(),
            aliases
                .into_iter()
                .filter(|(_, keys)| keys.len() == 1)
                .map(|(alias, keys)| (alias, keys.into_iter().next().unwrap()))
                .collect(),
        );
    }
    index
}


fn matched_keys(
    aliases: &BTreeMap<String, BTreeMap<String, String>>,
    symbol: &str,
    identifiers: &[String],
) -> BTreeSet<String> {
    identifiers
        .iter()
        .filter_map(|identifier| aliases.get(symbol).and_then(|a| a.get(identifier)))
        .cloned()
        .collect()
}

/// `_condition_identifiers_from_row`.
fn row_identifiers(row: &Row) -> Vec<String> {
    let mut identifiers = Vec::new();
    for value in [
        upper(&cell(row, "source_condition_id")),
        upper(&cell(row, "mondo_id")),
    ] {
        if !value.is_empty() {
            push_unique(&mut identifiers, &value);
        }
    }
    identifiers
}

fn empty_scope() -> Json {
    obj(vec![
        ("decision", st("")),
        ("category", st("")),
        ("review_status", st("")),
    ])
}

/// `_mechanism_evidence_record`.
fn mechanism_record(row: &Row, mechanism: &str) -> Json {
    obj(vec![
        ("source", st(&cell(row, "source"))),
        ("source_record_id", st(&cell(row, "source_record_id"))),
        ("condition_identifiers", strings(&row_identifiers(row))),
        ("condition_label", st(&cell(row, "disease_label"))),
        ("mechanism", st(mechanism)),
        ("mechanism_raw", st(&cell(row, "patho_mode_raw"))),
        ("allelic_requirement", st(&cell(row, "inheritance"))),
        ("mechanism_confidence", st(&cell(row, "mechanism_confidence"))),
        ("disease_confidence", st(&cell(row, "disease_confidence"))),
        (
            "source_scope",
            obj(vec![
                ("decision", st(&cell(row, "priva_scope"))),
                ("category", st(&cell(row, "disease_scope"))),
                ("review_status", st(&cell(row, "scope_review_status"))),
            ]),
        ),
        ("assertion_basis", st("curated")),
        ("pmids", strings(&split_multi(&cell(row, "pmids")))),
        ("evidence_url", st(&cell(row, "evidence_url"))),
    ])
}

/// `_source_penetrance_evidence_record`.
fn source_penetrance_record(row: &Row) -> Result<Option<Json>> {
    let raw_status = cell(row, "normalized_penetrance").to_lowercase();
    if !raw_status.is_empty() && raw_status != "complete" && raw_status != "incomplete" {
        let id = cell(row, "source_record_id");
        return Err(error(
            "ValueError",
            format!(
                "Unsupported normalized_penetrance {} for {}",
                crate::python_path::repr(&raw_status),
                if id.is_empty() { "<unknown>" } else { &id }
            ),
        ));
    }
    let normalized = normalize_penetrance_values(Values::List(&[raw_status.as_str()]));
    if normalized.is_empty() {
        return Ok(None);
    }
    let hpo_ids: Vec<String> = split_multi(&cell(row, "penetrance_hpo_ids"))
        .into_iter()
        .map(|v| upper(&v))
        .filter(|v| v.starts_with("HP:"))
        .collect();
    let mut raw_values = Vec::new();
    for value in split_multi(&cell(row, "penetrance_raw"))
        .iter()
        .chain(hpo_ids.iter())
    {
        push_unique(&mut raw_values, value);
    }
    Ok(Some(obj(vec![
        ("source", st(&cell(row, "source"))),
        ("source_record_id", st(&cell(row, "source_record_id"))),
        ("condition_identifiers", strings(&row_identifiers(row))),
        ("raw_values", strings(&raw_values)),
        ("normalized_value", st(normalized)),
        ("hpo_ids", strings(&hpo_ids)),
        ("frequency", st("")),
        ("evidence", st("")),
        ("reference", st("")),
        ("evidence_url", st(&cell(row, "evidence_url"))),
    ])))
}

/// `_attach_penetrance_to_condition`.
fn attach_penetrance(condition: &mut Condition, evidence: Json) {
    let status = clean(field(&evidence, "normalized_value")).to_lowercase();
    if !status.is_empty() {
        push_unique(&mut condition.penetrance.values, &status);
    }
    let key = penetrance_key(&evidence);
    if !condition
        .penetrance
        .assertions
        .iter()
        .any(|item| penetrance_key(item) == key)
    {
        condition.penetrance.assertions.push(evidence);
    }
}

/// `_attach_mechanism_to_condition`.
fn attach_mechanism(condition: &mut Condition, evidence: Json) {
    let mechanism = field(&evidence, "mechanism").to_owned();
    let block = condition.mechanisms.entry(mechanism).or_default();
    let requirement = field(&evidence, "allelic_requirement");
    if !requirement.is_empty() {
        push_unique(&mut block.allelic_requirements, requirement);
    }
    let key = mechanism_key(&evidence);
    if !block.evidence.iter().any(|item| mechanism_key(item) == key) {
        block.evidence.push(evidence);
    }
}

fn with_reason(evidence: &Json, reason: &str) -> Json {
    let mut record = evidence.clone();
    record.set("reason", st(reason));
    record
}

fn counter(stats: BTreeMap<String, i128>) -> Json {
    Json::Object(stats.into_iter().map(|(k, v)| (k, Json::Int(v))).collect())
}

fn ordered(stats: Vec<(&str, i128)>) -> Json {
    Json::Object(
        stats
            .into_iter()
            .map(|(k, v)| (k.to_owned(), Json::Int(v)))
            .collect(),
    )
}

/// `attach_condition_mechanisms`.
fn attach_condition_mechanisms(genes: &mut Genes, path: &Path) -> Result<Json> {
    let aliases = alias_index(genes);
    let (mut rows, mut mechanisms, mut matched, mut unmapped) = (0, 0, 0, 0);
    let (mut penetrance, mut penetrance_matched, mut penetrance_unmapped) = (0, 0, 0);
    let mut seen_mechanisms: BTreeMap<String, HashSet<Vec<String>>> = BTreeMap::new();
    let mut seen_penetrance: BTreeMap<String, HashSet<Vec<String>>> = BTreeMap::new();
    let mut reader = Reader::open(path)?;
    reader.require(path, MECHANISM_REQUIRED_COLUMNS)?;
    while let Some(row) = reader.next_row()? {
        if !CONDITION_MECHANISM_SOURCES.contains(&cell(&row, "source").as_str()) {
            continue;
        }
        rows += 1;
        let symbol = upper(&cell(&row, "gene_symbol"));
        if symbol.is_empty() {
            continue;
        }
        let gene = genes.entry(symbol.clone()).or_default();
        let keys = matched_keys(&aliases, &symbol, &row_identifiers(&row));
        let reason = if keys.len() > 1 {
            "conflicting_condition_identifiers"
        } else {
            "no_exact_hpo_condition_identifier"
        };
        if let Some(evidence) = source_penetrance_record(&row)? {
            penetrance += 1;
            if keys.len() == 1 {
                let key = keys.iter().next().unwrap();
                attach_penetrance(gene.conditions.get_mut(key).unwrap(), evidence);
                penetrance_matched += 1;
            } else if seen_penetrance
                .entry(symbol.clone())
                .or_default()
                .insert(penetrance_key(&evidence))
            {
                gene.unmapped_penetrance.push(with_reason(&evidence, reason));
                penetrance_unmapped += 1;
            }
        }
        let tokens: BTreeSet<String> = regex_split(&cell(&row, "normalized_mechanisms"))
            .into_iter()
            .map(|t| upper(crate::hpo_penetrance::py_strip(&t)))
            .filter(|t| !t.is_empty())
            .collect();
        for mechanism in tokens
            .iter()
            .filter(|t| CANONICAL_MECHANISMS.contains(&t.as_str()))
        {
            mechanisms += 1;
            let evidence = mechanism_record(&row, mechanism);
            if keys.len() == 1 {
                let condition = gene.conditions.get_mut(keys.iter().next().unwrap()).unwrap();
                let label = field(&evidence, "condition_label").to_owned();
                attach_mechanism(condition, evidence);
                if condition.label.is_empty() && !label.is_empty() {
                    condition.label = label;
                }
                matched += 1;
                continue;
            }
            let record = with_reason(&evidence, reason);
            if seen_mechanisms
                .entry(symbol.clone())
                .or_default()
                .insert(mechanism_key(&record))
            {
                gene.unmapped_mechanisms.push(record);
                unmapped += 1;
            }
        }
    }
    for gene in genes.values_mut() {
        for condition in gene.conditions.values_mut() {
            for block in condition.mechanisms.values_mut() {
                block.allelic_requirements.sort();
                block.evidence.sort_by_key(mechanism_key);
            }
            condition.penetrance.values.sort();
            condition.penetrance.assertions.sort_by_key(penetrance_key);
        }
        gene.unmapped_mechanisms.sort_by_key(mechanism_key);
        gene.unmapped_penetrance.sort_by_key(penetrance_key);
    }
    Ok(ordered(vec![
        ("source_rows", rows),
        ("mechanism_records", mechanisms),
        ("mechanisms_matched", matched),
        ("mechanisms_unmapped", unmapped),
        ("penetrance_records", penetrance),
        ("penetrance_matched", penetrance_matched),
        ("penetrance_unmapped", penetrance_unmapped),
    ]))
}

/// `re.split(r"[|;,]", text)`.
fn regex_split(text: &str) -> Vec<String> {
    text.split(['|', ';', ',']).map(str::to_owned).collect()
}

/// `_condition_has_explicit_inheritance`.
fn has_explicit_inheritance(condition: &Condition) -> bool {
    !condition.inheritance.values.is_empty()
        || condition.mechanisms.values().any(|block| {
            block
                .evidence
                .iter()
                .any(|e| !clean(field(e, "allelic_requirement")).is_empty())
        })
}

/// `_clingen_hi3_allelic_requirement`.
fn hi3_requirement(row: &Row, condition: &Condition) -> &'static str {
    if cell(row, "Haploinsufficiency Score") != "3" {
        return "";
    }
    let location = cell(row, "Genomic Location");
    let mut chromosome = location
        .split(':')
        .next()
        .unwrap_or("")
        .to_lowercase();
    if let Some(rest) = chromosome.strip_prefix("chr") {
        chromosome = rest.to_owned();
    }
    if !chromosome.is_empty() && chromosome.bytes().all(|b| b.is_ascii_digit()) {
        let digits = chromosome.trim_start_matches('0');
        if digits.len() <= 2 && (1..=22).contains(&digits.parse::<u32>().unwrap_or(0)) {
            return "monoallelic_autosomal";
        }
        return "";
    }
    match chromosome.as_str() {
        "x" if has_explicit_inheritance(condition) => "",
        "x" => "monoallelic_X_heterozygous",
        "y" => "monoallelic_Y_hemizygous",
        _ => "",
    }
}

/// `attach_clingen_haploinsufficiency`.
fn attach_clingen(genes: &mut Genes, path: &Path) -> Result<Json> {
    let mut stats: BTreeMap<String, i128> = BTreeMap::new();
    let count = |stats: &mut BTreeMap<String, i128>, key: &str| {
        *stats.entry(key.to_owned()).or_default() += 1;
    };
    let mut reader = Reader::open(path)?;
    reader.require(path, CLINGEN_REQUIRED_COLUMNS)?;
    while let Some(row) = reader.next_row()? {
        let score = cell(&row, "Haploinsufficiency Score");
        let strength = match score.as_str() {
            "3" => "sufficient_evidence",
            "2" => "some_evidence",
            "1" => "little_evidence",
            _ => "",
        };
        let disease = upper(&cell(&row, "Haploinsufficiency Disease ID"));
        let symbol = upper(&cell(&row, "#Gene Symbol"));
        if strength.is_empty() || disease.is_empty() || symbol.is_empty() {
            count(&mut stats, "rows_without_a_usable_score_or_disease");
            continue;
        }
        let Some(gene) = genes.get_mut(&symbol) else {
            count(&mut stats, "gene_absent_from_the_hpo_frame");
            continue;
        };
        let matched: Vec<String> = gene
            .conditions
            .iter()
            .filter(|(_, c)| {
                c.identifiers
                    .get("MONDO")
                    .is_some_and(|ids| ids.contains(&disease))
            })
            .map(|(k, _)| k.clone())
            .collect();
        if matched.is_empty() {
            count(&mut stats, "disease_absent_from_the_gene");
            continue;
        }
        let description = cell(&row, "Haploinsufficiency Description");
        let pmids: Vec<String> = (1..=6)
            .map(|i| cell(&row, &format!("Haploinsufficiency PMID{i}")))
            .filter(|v| !v.is_empty())
            .collect();
        for key in matched {
            let condition = gene.conditions.get_mut(&key).unwrap();
            let requirement = hi3_requirement(&row, condition);
            let had_lof = condition.mechanisms.contains_key("LOF");
            attach_mechanism(
                condition,
                obj(vec![
                    ("source", st("ClinGen_haploinsufficiency")),
                    ("source_record_id", st(&disease)),
                    ("condition_identifiers", strings(std::slice::from_ref(&disease))),
                    ("condition_label", st(&description)),
                    ("mechanism", st("LOF")),
                    ("mechanism_raw", st(&description)),
                    ("allelic_requirement", st(requirement)),
                    ("mechanism_confidence", st(strength)),
                    ("disease_confidence", st("")),
                    ("assertion_basis", st("curated")),
                    ("source_scope", empty_scope()),
                    ("pmids", strings(&pmids)),
                    (
                        "evidence_url",
                        st("https://search.clinicalgenome.org/kb/gene-dosage"),
                    ),
                ]),
            );
            count(&mut stats, "conditions_matched");
            if score == "3" {
                let name = if requirement.is_empty() {
                    "allelic_state_deferred_to_condition_evidence"
                } else {
                    requirement
                };
                count(&mut stats, &format!("hi3_{name}"));
            }
            count(
                &mut stats,
                if had_lof {
                    "conditions_that_already_had_lof"
                } else {
                    "conditions_that_gained_lof"
                },
            );
        }
    }
    Ok(counter(stats))
}

// ---------------------------------------------------------------- GoFCards

/// `value or {}` then `.items()`.
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

/// Iterating a parsed value as Python does.
fn iterate(value: &Json) -> Result<Vec<Json>> {
    match value {
        Json::Array(items) => Ok(items.clone()),
        Json::Object(entries) => Ok(entries.iter().map(|(k, _)| st(k)).collect()),
        Json::Str(text) => Ok(text.chars().map(|c| Json::Str(c.to_string())).collect()),
        other => Err(error(
            "TypeError",
            format!("'{}' object is not iterable", type_name(other)),
        )),
    }
}

/// `value.get(key)`, `None` when absent.
fn pyget(value: &Json, key: &str) -> Result<Json> {
    Ok(get(value, key)?.cloned().unwrap_or(Json::Null))
}

/// `value or default`.
fn or(value: Json, default: Json) -> Json {
    if truthy(&value) { value } else { default }
}

fn empty_object() -> Json {
    Json::Object(vec![])
}

/// `_normalize_clinvar_condition_identifier`.
fn normalize_clinvar_identifier(namespace: &Json, value: &Json) -> String {
    let database: String = upper(&clean_json(namespace))
        .chars()
        .filter(|c| c.is_ascii_uppercase() || c.is_ascii_digit())
        .collect();
    let identifier = upper(&clean_json(value));
    if database.is_empty() || identifier.is_empty() {
        return String::new();
    }
    let prefix = match database.as_str() {
        "OMIM" => "OMIM",
        "MONDO" => "MONDO",
        "ORPHA" | "ORPHANET" => "ORPHA",
        _ => return String::new(),
    };
    let mut identifier = identifier.as_str();
    if let Some((incoming, rest)) = identifier.split_once(':') {
        let incoming = if incoming == "ORPHANET" || incoming == "ORPHA" {
            "ORPHA"
        } else {
            incoming
        };
        if incoming != prefix {
            return String::new();
        }
        identifier = rest;
    }
    format!("{prefix}:{identifier}")
}

/// `_clinvar_condition_identities`.
fn clinvar_identities(assertion: &Json) -> Result<Vec<String>> {
    let mut identities = Vec::new();
    let conditions = or(get(assertion, "conditions")?.cloned().unwrap_or(Json::Array(vec![])), Json::Array(vec![]));
    for condition in iterate(&conditions)? {
        if !matches!(condition, Json::Object(_)) {
            continue;
        }
        let identifier = normalize_clinvar_identifier(
            &pyget(&condition, "database")?,
            &pyget(&condition, "id")?,
        );
        if !identifier.is_empty() {
            push_unique(&mut identities, &identifier);
        }
    }
    let scvs = or(get(assertion, "matched_scvs")?.cloned().unwrap_or(Json::Array(vec![])), Json::Array(vec![]));
    for scv in iterate(&scvs)? {
        if !matches!(scv, Json::Object(_)) {
            continue;
        }
        let mappings = or(get(&scv, "trait_mappings")?.cloned().unwrap_or(Json::Array(vec![])), Json::Array(vec![]));
        for mapping in iterate(&mappings)? {
            if !matches!(mapping, Json::Object(_)) {
                continue;
            }
            let identifier = normalize_clinvar_identifier(
                &pyget(&mapping, "mapping_ref")?,
                &pyget(&mapping, "mapping_value")?,
            );
            if !identifier.is_empty() {
                push_unique(&mut identities, &identifier);
            }
        }
    }
    Ok(identities)
}

/// `clinvar_vcv.iter_gofcards_variants`.
fn gofcards_variants(cache: &Json) -> Result<Vec<(String, String, Json)>> {
    let mut variants = Vec::new();
    for (symbol, gene) in items(&pyget(cache, "genes")?)? {
        for (variant_id, variant) in items(&pyget(&gene, "variants")?)? {
            variants.push((symbol.clone(), variant_id, variant));
        }
    }
    Ok(variants)
}

/// `_clinvar_links_by_variant`.
fn clinvar_links(cache: &Json) -> Result<BTreeMap<(String, String), Vec<Json>>> {
    let mut index: BTreeMap<(String, String), Vec<Json>> = BTreeMap::new();
    for (symbol, variant_id, variant) in gofcards_variants(cache)? {
        let mut blocks = Vec::new();
        let clinvar = pyget(&variant, "clinvar")?;
        if truthy(&clinvar) {
            blocks.push(clinvar);
        }
        blocks.extend(iterate(&or(
            pyget(&variant, "clinvar_additional")?,
            Json::Array(vec![]),
        ))?);
        let key = (upper(crate::hpo_penetrance::py_strip(&symbol)), variant_id);
        for block in blocks {
            let vcv = clean_json(&pyget(&block, "vcv_accession")?);
            let hgvs: Vec<Json> = iterate(&or(pyget(&block, "hgvs")?, Json::Array(vec![])))?
                .into_iter()
                .filter(|h| !clean_json(h).is_empty())
                .collect();
            for assertion in iterate(&or(
                pyget(&block, "condition_assertions")?,
                Json::Array(vec![]),
            ))? {
                let conditions = or(pyget(&assertion, "conditions")?, Json::Array(vec![]));
                let classification = or(pyget(&assertion, "germline_classification")?, empty_object());
                let rcv = clean_json(&pyget(&assertion, "rcv_accession")?);
                let significance = clean_json(&pyget(&classification, "clinical_significance")?);
                let stars = pyget(&classification, "review_stars")?;
                let identities = clinvar_identities(&assertion)?;
                let mut names = Vec::new();
                for condition in iterate(&conditions)? {
                    let name = clean_json(&pyget(&condition, "name")?);
                    if !name.is_empty() {
                        names.push(name);
                    }
                }
                index.entry(key.clone()).or_default().push(obj(vec![
                    ("vcv_accession", st(&vcv)),
                    ("hgvs", Json::Array(hgvs.clone())),
                    ("rcv_accession", st(&rcv)),
                    ("clinical_significance", st(&significance)),
                    ("review_stars", stars),
                    ("condition_identifiers", strings(&identities)),
                    ("condition_names", strings(&names)),
                ]));
            }
        }
    }
    Ok(index)
}

/// `_gofcards_handover`.
fn handover(symbol: &str, variant_id: &str, variant: &Json) -> Result<Json> {
    let mut views = Vec::new();
    let mut assemblies = items(&pyget(variant, "assemblies")?)?;
    assemblies.sort_by(|a, b| a.0.cmp(&b.0));
    for (assembly, block) in assemblies {
        let mut transcripts = items(&pyget(&block, "transcripts")?)?;
        transcripts.sort_by(|a, b| a.0.cmp(&b.0));
        for (transcript, view) in transcripts {
            let mut by_hgvsc = items(&pyget(&view, "by_hgvsc")?)?;
            by_hgvsc.sort_by(|a, b| a.0.cmp(&b.0));
            for (hgvsc, detail) in by_hgvsc {
                let detail = or(detail, empty_object());
                views.push(obj(vec![
                    ("assembly", st(&assembly)),
                    ("transcript", st(&transcript)),
                    ("hgvsc", st(&hgvsc)),
                    ("hgvsp", st(&clean_json(&pyget(&detail, "hgvsp")?))),
                ]));
            }
        }
    }
    Ok(obj(vec![
        ("mechanism", st("GOF")),
        ("symbol", st(symbol)),
        ("gofcards_variant_id", st(variant_id)),
        ("transcripts", Json::Array(views)),
        ("clinvar_links", Json::Array(vec![])),
    ]))
}

/// `_public_variant_record`.
fn public_variant(variant: &Json, status: &str, key: &str, links: &[Json], reason: &str) -> Json {
    let mut record = variant.clone();
    let mut link = obj(vec![("status", st(status)), ("condition_key", st(key))]);
    if !reason.is_empty() {
        link.set("reason", st(reason));
    }
    record.set("condition_link", link);
    record.set("clinvar_links", Json::Array(links.to_vec()));
    record
}

fn link_strings(link: &Json, name: &str) -> Vec<String> {
    match link.get(name) {
        Some(Json::Array(items)) => items
            .iter()
            .map(|i| i.as_str().unwrap_or("").to_owned())
            .collect(),
        _ => vec![],
    }
}

/// `_condition_link_mechanism_evidence`.
fn link_evidence(variant_key: &str, links: &[Json]) -> Json {
    let identifiers: BTreeSet<String> = links
        .iter()
        .flat_map(|l| link_strings(l, "condition_identifiers"))
        .collect();
    let vcvs: BTreeSet<String> = links
        .iter()
        .map(|l| field(l, "vcv_accession").to_owned())
        .filter(|v| !v.is_empty())
        .collect();
    let record_id = if vcvs.is_empty() {
        variant_key.to_owned()
    } else {
        vcvs.into_iter().collect::<Vec<_>>().join(";")
    };
    obj(vec![
        ("source", st("GoFCards_exact+ClinVar_VCV")),
        ("source_record_id", st(&record_id)),
        (
            "condition_identifiers",
            strings(&identifiers.into_iter().collect::<Vec<_>>()),
        ),
        ("condition_label", st("")),
        ("mechanism", st("GOF")),
        ("mechanism_raw", st("gain of function")),
        ("allelic_requirement", st("")),
        ("assertion_basis", st("curated")),
        ("mechanism_confidence", st("exact_variant")),
        ("disease_confidence", st("ClinVar_germline_assertion")),
        ("source_scope", empty_scope()),
        ("pmids", Json::Array(vec![])),
        ("evidence_url", st("https://www.ncbi.nlm.nih.gov/clinvar/")),
    ])
}

/// `attach_gofcards_variants`.
fn attach_gofcards(genes: &mut Genes, path: &Path) -> Result<Json> {
    let cache = shared::load_gofcards_cache(path)?;
    let links_by_variant = clinvar_links(&cache)?;
    let aliases = alias_index(genes);
    let mut eligible = Vec::new();
    let mut quarantined = 0i128;
    for (symbol, variant_id, variant) in gofcards_variants(&cache)? {
        if shared::eligible(&variant)? {
            eligible.push((symbol, variant_id, variant));
        } else {
            quarantined += 1;
        }
    }
    let eligible_count = eligible.len() as i128;
    let (mut linked, mut link_count, mut unmapped) = (0i128, 0i128, 0i128);
    eligible.sort_by(|a, b| (upper(&a.0), &a.1).cmp(&(upper(&b.0), &b.1)));
    for (source_symbol, variant_id, source_variant) in eligible {
        let symbol = upper(&source_symbol);
        let cleaned = clean(&variant_id);
        if cleaned.is_empty() {
            return Err(error("ValueError", "GoFCards variant has no identifier"));
        }
        let variant_key = format!("GOFCARDS:{cleaned}");
        let variant = handover(&source_symbol, &variant_id, &source_variant)?;
        let gene = genes.entry(symbol.clone()).or_default();
        // A repeated (VCV, RCV, identifiers) keeps its first place and last value.
        let mut links: Vec<(Vec<String>, Json)> = Vec::new();
        for link in links_by_variant
            .get(&(symbol.clone(), variant_id.clone()))
            .into_iter()
            .flatten()
        {
            let key = vec![
                field(link, "vcv_accession").to_owned(),
                field(link, "rcv_accession").to_owned(),
                link_strings(link, "condition_identifiers").join(";"),
            ];
            match links.iter_mut().find(|(k, _)| *k == key) {
                Some(slot) => slot.1 = link.clone(),
                None => links.push((key, link.clone())),
            }
        }
        let links: Vec<Json> = links.into_iter().map(|(_, l)| l).collect();
        let mut by_condition: BTreeMap<String, Vec<Json>> = BTreeMap::new();
        let mut ambiguous = false;
        for link in &links {
            let keys = matched_keys(&aliases, &symbol, &link_strings(link, "condition_identifiers"));
            if keys.len() == 1 {
                by_condition
                    .entry(keys.into_iter().next().unwrap())
                    .or_default()
                    .push(link.clone());
            } else if keys.len() > 1 {
                ambiguous = true;
            }
        }
        if !by_condition.is_empty() {
            linked += 1;
            for (key, condition_links) in by_condition {
                let condition = gene.conditions.get_mut(&key).unwrap();
                attach_mechanism(condition, link_evidence(&variant_key, &condition_links));
                condition
                    .mechanisms
                    .get_mut("GOF")
                    .unwrap()
                    .variants
                    .insert(
                        variant_key.clone(),
                        public_variant(&variant, "exact", &key, &condition_links, ""),
                    );
                link_count += 1;
            }
            continue;
        }
        let reason = if ambiguous {
            "ambiguous_clinvar_condition_identifiers"
        } else if !links.is_empty() {
            "no_exact_hpo_condition_identifier"
        } else {
            "no_exact_clinvar_condition_link"
        };
        gene.unmapped_variants.insert(
            variant_key,
            public_variant(&variant, "unresolved", "", &links, reason),
        );
        unmapped += 1;
    }
    for gene in genes.values_mut() {
        for condition in gene.conditions.values_mut() {
            for block in condition.mechanisms.values_mut() {
                block.evidence.sort_by_key(mechanism_key);
            }
        }
    }
    Ok(ordered(vec![
        ("source_variants", eligible_count + quarantined),
        ("eligible_variants", eligible_count),
        ("quarantined_variants", quarantined),
        ("condition_linked_variants", linked),
        ("condition_variant_links", link_count),
        ("unmapped_variants", unmapped),
    ]))
}

/// `deduce_mechanisms_from_inheritance`.
fn deduce(genes: &mut Genes) -> Json {
    let mut stats: BTreeMap<String, i128> = BTreeMap::new();
    for gene in genes.values_mut() {
        for (key, condition) in gene.conditions.iter_mut() {
            if condition.scope.decision != "include" {
                continue;
            }
            if !condition.mechanisms.is_empty() {
                *stats.entry("already_had_a_mechanism".into()).or_default() += 1;
                continue;
            }
            let (inheritance, _) = normalize_inheritance("", &condition.inheritance.values);
            if inheritance != "recessive" && inheritance != "x_linked_recessive" {
                let name = if inheritance.is_empty() {
                    "no_inheritance"
                } else {
                    inheritance.as_str()
                };
                *stats.entry(format!("left_unresolved_{name}")).or_default() += 1;
                continue;
            }
            let label = clean(&condition.label);
            attach_mechanism(
                condition,
                obj(vec![
                    ("source", st("deduced_from_inheritance")),
                    ("source_record_id", st(key)),
                    ("condition_identifiers", strings(std::slice::from_ref(key))),
                    ("condition_label", st(&label)),
                    ("mechanism", st("LOF")),
                    (
                        "mechanism_raw",
                        st(&format!("{inheritance} inheritance implies loss of function")),
                    ),
                    ("allelic_requirement", st("")),
                    ("mechanism_confidence", st("")),
                    ("disease_confidence", st("")),
                    ("assertion_basis", st("deduced")),
                    ("source_scope", empty_scope()),
                    ("pmids", Json::Array(vec![])),
                    ("evidence_url", st("")),
                ]),
            );
            *stats
                .entry(format!("{inheritance}_conditions_given_lof"))
                .or_default() += 1;
        }
    }
    counter(stats)
}

/// `apply_condition_scope_gate`.
fn scope_gate(genes: &mut Genes) -> Json {
    let (mut examined, mut excluded, mut include, mut review, mut unscoped) = (0, 0, 0, 0, 0);
    let (mut genes_excluding, mut genes_emptied) = (0, 0);
    for gene in genes.values_mut() {
        let mut count = 0;
        let conditions = std::mem::take(&mut gene.conditions);
        for (key, condition) in conditions {
            examined += 1;
            let decision = clean(&condition.scope.decision).to_lowercase();
            if decision == "exclude" {
                count += 1;
                excluded += 1;
                continue;
            }
            match decision.as_str() {
                "include" => include += 1,
                "review" => review += 1,
                _ => unscoped += 1,
            }
            gene.conditions.insert(key, condition);
        }
        gene.excluded = count;
        if count > 0 {
            genes_excluding += 1;
            if gene.conditions.is_empty() {
                genes_emptied += 1;
            }
        }
    }
    ordered(vec![
        ("conditions_examined", examined),
        ("conditions_excluded", excluded),
        ("conditions_retained_include", include),
        ("conditions_retained_review", review),
        ("conditions_retained_unscoped", unscoped),
        ("genes_with_excluded_conditions", genes_excluding),
        ("genes_left_without_conditions", genes_emptied),
    ])
}

fn sorted_counts(counts: BTreeMap<String, i128>) -> Json {
    counter(counts)
}

/// `_refresh_gene_summary` over the final gene.
fn summary(gene: &Gene) -> Json {
    let conditions: Vec<&Condition> = gene.conditions.values().collect();
    let mut mechanisms: BTreeSet<String> = BTreeSet::new();
    let mut condition_counts: BTreeMap<String, i128> = BTreeMap::new();
    for condition in &conditions {
        for mechanism in condition.mechanisms.keys() {
            mechanisms.insert(mechanism.clone());
            *condition_counts.entry(mechanism.clone()).or_default() += 1;
        }
    }
    let mut unmapped_counts: BTreeMap<String, i128> = BTreeMap::new();
    for record in &gene.unmapped_mechanisms {
        let mechanism = upper(&clean(field(record, "mechanism")));
        if !mechanism.is_empty() {
            mechanisms.insert(mechanism.clone());
            *unmapped_counts.entry(mechanism).or_default() += 1;
        }
    }
    let mut variant_counts: BTreeMap<String, i128> = BTreeMap::new();
    for variant in gene.unmapped_variants.values() {
        let mechanism = upper(&clean(field(variant, "mechanism")));
        if !mechanism.is_empty() {
            mechanisms.insert(mechanism.clone());
            *variant_counts.entry(mechanism).or_default() += 1;
        }
    }
    let count = |test: &dyn Fn(&Condition) -> bool| -> Json {
        Json::Int(conditions.iter().filter(|c| test(c)).count() as i128)
    };
    obj(vec![
        ("condition_count", Json::Int(conditions.len() as i128)),
        (
            "conditions_with_inheritance",
            count(&|c| !c.inheritance.values.is_empty()),
        ),
        (
            "conditions_with_penetrance",
            count(&|c| !c.penetrance.values.is_empty()),
        ),
        ("conditions_with_onset", count(&|c| !c.onset.values.is_empty())),
        (
            "pathogenic_mechanisms",
            strings(&mechanisms.into_iter().collect::<Vec<_>>()),
        ),
        (
            "condition_resolved_mechanism_counts",
            sorted_counts(condition_counts),
        ),
        (
            "unmapped_mechanism_count",
            Json::Int(unmapped_counts.values().sum()),
        ),
        ("unmapped_mechanism_counts", sorted_counts(unmapped_counts)),
        (
            "unmapped_penetrance_count",
            Json::Int(gene.unmapped_penetrance.len() as i128),
        ),
        (
            "unmapped_variant_count",
            Json::Int(gene.unmapped_variants.len() as i128),
        ),
        ("unmapped_variant_mechanism_counts", sorted_counts(variant_counts)),
    ])
}

fn axis_json(axis: &Axis, values_key: &str) -> Json {
    obj(vec![
        (values_key, strings(&axis.values)),
        ("assertions", Json::Array(axis.assertions.clone())),
    ])
}

fn condition_json(condition: &Condition) -> Json {
    let scope = &condition.scope;
    obj(vec![
        ("label", st(&condition.label)),
        (
            "identifiers",
            Json::Object(
                condition
                    .identifiers
                    .iter()
                    .map(|(ns, ids)| (ns.clone(), strings(ids)))
                    .collect(),
            ),
        ),
        (
            "priva_scope",
            obj(vec![
                ("decision", st(&scope.decision)),
                ("category", st(&scope.category)),
                ("review_status", st(&scope.review_status)),
                ("review_statuses", strings(&scope.review_statuses)),
                ("evidence", strings(&scope.evidence)),
                ("references", strings(&scope.references)),
            ]),
        ),
        ("inheritance", axis_json(&condition.inheritance, "modes")),
        ("penetrance", axis_json(&condition.penetrance, "statuses")),
        ("onset", axis_json(&condition.onset, "terms")),
        (
            "pathogenic_mechanisms",
            Json::Object(
                condition
                    .mechanisms
                    .iter()
                    .map(|(name, block)| {
                        (
                            name.clone(),
                            obj(vec![
                                ("allelic_requirements", strings(&block.allelic_requirements)),
                                ("evidence", Json::Array(block.evidence.clone())),
                                (
                                    "variants",
                                    Json::Object(
                                        block
                                            .variants
                                            .iter()
                                            .map(|(k, v)| (k.clone(), v.clone()))
                                            .collect(),
                                    ),
                                ),
                            ]),
                        )
                    })
                    .collect(),
            ),
        ),
        ("hpo_assertion_count", Json::Int(condition.hpo_assertion_count)),
    ])
}

fn gene_json(gene: &Gene) -> Json {
    obj(vec![
        ("summary", summary(gene)),
        (
            "conditions",
            Json::Object(
                gene.conditions
                    .iter()
                    .map(|(k, c)| (k.clone(), condition_json(c)))
                    .collect(),
            ),
        ),
        (
            "unmapped_evidence",
            obj(vec![
                ("mechanisms", Json::Array(gene.unmapped_mechanisms.clone())),
                ("penetrance", Json::Array(gene.unmapped_penetrance.clone())),
                (
                    "variants",
                    Json::Object(
                        gene.unmapped_variants
                            .iter()
                            .map(|(k, v)| (k.clone(), v.clone()))
                            .collect(),
                    ),
                ),
            ]),
        ),
        (
            "scope_gate",
            obj(vec![("excluded_condition_count", Json::Int(gene.excluded))]),
        ),
    ])
}

pub(super) struct Inputs {
    pub(super) builder: PathBuf,
    pub(super) hpo_assertions: PathBuf,
    pub(super) mechanism_evidence: PathBuf,
    pub(super) clingen_dosage: PathBuf,
    pub(super) gofcards_variants: PathBuf,
    pub(super) hpo_release: String,
    pub(super) mondo_release: String,
}

/// `build_cache_payload`.
pub(super) fn build_cache_payload(inputs: &Inputs) -> Result<Json> {
    let mut genes = hpo_frame(&inputs.hpo_assertions)?;
    let mechanism_stats = attach_condition_mechanisms(&mut genes, &inputs.mechanism_evidence)?;
    let clingen_stats = attach_clingen(&mut genes, &inputs.clingen_dosage)?;
    let variant_stats = attach_gofcards(&mut genes, &inputs.gofcards_variants)?;
    let deduced_stats = deduce(&mut genes);
    let gate_stats = scope_gate(&mut genes);
    let genes_json = Json::Object(
        genes
            .iter()
            .map(|(symbol, gene)| (symbol.clone(), gene_json(gene)))
            .collect(),
    );
    let built_at = shared::iso_now();
    let sources = obj(vec![
        ("builder_script", source_metadata(&inputs.builder)?),
        ("hpo_assertions", source_metadata(&inputs.hpo_assertions)?),
        ("mechanism_evidence", source_metadata(&inputs.mechanism_evidence)?),
        ("clingen_dosage", source_metadata(&inputs.clingen_dosage)?),
        ("gofcards_variants", source_metadata(&inputs.gofcards_variants)?),
    ]);
    let counts = super::validate::cache_counts(&genes_json)?;
    let payload = obj(vec![
        (
            "_meta",
            obj(vec![
                ("schema_version", st(SCHEMA_VERSION)),
                ("built_at", st(&built_at)),
                (
                    "releases",
                    obj(vec![
                        ("HPO", st(&clean(&inputs.hpo_release))),
                        ("MONDO", st(&clean(&inputs.mondo_release))),
                    ]),
                ),
                ("sources", sources),
                (
                    "build_statistics",
                    obj(vec![
                        ("mechanisms", mechanism_stats),
                        ("clingen_haploinsufficiency", clingen_stats),
                        ("variants", variant_stats),
                        ("deduced_mechanisms", deduced_stats),
                        ("scope_gate", gate_stats),
                    ]),
                ),
                ("counts", counts),
            ]),
        ),
        ("genes", genes_json),
    ]);
    super::validate::validate_payload(&payload)?;
    Ok(payload)
}
