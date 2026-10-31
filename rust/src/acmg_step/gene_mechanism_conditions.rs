//! `scripts/gene_mechanism_conditions.py` (1,060 lines): the condition cache's assertions,
//! inheritance and penetrance normalisation.
//!
//! Thread T2 (gene-condition hub). Converted, the functions on the step's path:
//! `condition_cache_context`, `condition_cache_mechanism_assertions` (with
//! `_hpo_allelic_requirements` and `_condition_allelic_requirements`), `normalize_inheritance`,
//! `normalize_inheritance_histories`, `gene_has_lof_mechanism_history`,
//! `gene_inheritance_from_constraint` and `normalize_penetrance`.
//!
//! Not converted (they serve the hub's diagnostic summaries and the group 4 builders):
//! `_hpo_inheritance_flags`, `build_hpo_condition_index`, `condition_cache_mechanism_entries`,
//! `enrich_condition_mechanism_assertion`, `gene_inheritance_consensus`.
//!
//! The cache is read as JSON (`json_text`), so its values keep Python's kinds. Where Python would
//! raise on a value of the wrong kind (iterating `None`, calling `.get` on a list) the same
//! exception is raised here; where Python would read a text's characters or a dict's keys as if
//! they were list items, the value is refused.
use super::gene_mechanism_common::{
    CANONICAL_MECHANISMS, CONDITION_MECHANISM_SOURCES, HPO_CACHE_INHERITANCE_LABELS,
    UNRESOLVED_MECHANISM, clean, clean_json, clean_text, deduplicate_by, json_column_value,
    json_type_name, python_error, refused, result_table, show_json, show_outcome, show_text,
    truthy,
};
use super::types::{HarnessCall, not_converted};
use crate::cli::{Result, invalid};
use crate::hpo_penetrance::{Values, normalize_penetrance_values};
use crate::json_text::Json;
use crate::pandas_table::{Cell, Column, Table};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::rc::Rc;

/// Inheritance values that are germline but not single-gene.
pub(crate) const NON_MENDELIAN_INHERITANCE: [&str; 4] =
    ["non_mendelian", "polygenic", "digenic", "oligogenic"];
const X_LINKED_INHERITANCE: [&str; 3] = [
    "x_linked_recessive",
    "x_linked_dominant",
    "x_linked_unspecified",
];

// ------------------------------------------------------------------- Python access to JSON

/// `value.get(key)`; on anything but a dict it is Python's AttributeError.
fn py_get<'a>(value: &'a Json, key: &str) -> Result<Option<&'a Json>> {
    match value {
        Json::Object(_) => Ok(value.get(key)),
        other => Err(python_error(
            "AttributeError",
            format!("'{}' object has no attribute 'get'", json_type_name(other)),
        )),
    }
}

/// `value.get(key, {})` kept only when it is a dict (`x if isinstance(x, dict) else {}`).
fn dict_or_none<'a>(value: Option<&'a Json>) -> Option<&'a Json> {
    value.filter(|value| matches!(value, Json::Object(_)))
}

/// `for item in value` where the items are used as values: a list gives its
/// items, an empty text or dict nothing; the characters of a text or the keys
/// of a dict are refused; `None`, a bool or a number is Python's TypeError.
fn iterate(value: &Json) -> Result<&[Json]> {
    match value {
        Json::Array(items) => Ok(items),
        Json::Str(text) if text.is_empty() => Ok(&[]),
        Json::Object(items) if items.is_empty() => Ok(&[]),
        Json::Str(_) | Json::Object(_) => Err(refused(format!(
            "iterating a {} where the cache holds a list",
            json_type_name(value)
        ))),
        other => Err(python_error(
            "TypeError",
            format!("'{}' object is not iterable", json_type_name(other)),
        )),
    }
}

/// `for item in (value or [])`.
fn iterate_or_empty(value: Option<&Json>) -> Result<&[Json]> {
    match value {
        Some(value) if truthy(value) => iterate(value),
        _ => Ok(&[]),
    }
}

/// `for item in (value or []) if isinstance(item, dict)`: the characters of a
/// text and the keys of a dict are never dicts, so they give nothing.
fn dicts_or_empty(value: Option<&Json>) -> Result<Vec<&Json>> {
    match value {
        Some(Json::Str(_)) | Some(Json::Object(_)) => Ok(Vec::new()),
        _ => Ok(iterate_or_empty(value)?
            .iter()
            .filter(|item| matches!(item, Json::Object(_)))
            .collect()),
    }
}

/// The items of a dict sorted by key (`sorted(d.items())`; keys are unique).
fn sorted_items(items: &[(String, Json)]) -> Vec<(&String, &Json)> {
    let mut sorted: Vec<(&String, &Json)> = items.iter().map(|(key, value)| (key, value)).collect();
    sorted.sort_by(|a, b| a.0.cmp(b.0));
    sorted
}

/// `list(value or [])` of a JSON value.
fn list_or_empty(value: Option<&Json>) -> Result<Vec<Json>> {
    Ok(iterate_or_empty(value)?.to_vec())
}

// --------------------------------------------------------------------------- the context

/// One HPO assertion as the context keeps it (`hpo_id`, `frequency`, `evidence`, `reference`).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct HpoAssertion {
    pub(crate) hpo_id: String,
    pub(crate) frequency: String,
    pub(crate) evidence: String,
    pub(crate) reference: String,
}

impl HpoAssertion {
    fn from_fields(hpo_id: String, record: &Json) -> Result<HpoAssertion> {
        Ok(HpoAssertion {
            hpo_id,
            frequency: clean_json(record.get("frequency"))?,
            evidence: clean_json(record.get("evidence"))?,
            reference: clean_json(record.get("reference"))?,
        })
    }

    fn to_json(&self) -> Json {
        json_object(vec![
            ("hpo_id", Json::Str(self.hpo_id.clone())),
            ("frequency", Json::Str(self.frequency.clone())),
            ("evidence", Json::Str(self.evidence.clone())),
            ("reference", Json::Str(self.reference.clone())),
        ])
    }
}

/// What `condition_cache_context` returns for an included condition (the
/// fields `hpo_match_status` and `priva_scope` are always
/// `matched_gene_and_condition` and `include`).
#[derive(Clone, Debug)]
pub(crate) struct ConditionContext {
    pub(crate) hpo_disease_id: String,
    pub(crate) mondo_id: String,
    pub(crate) disease_scope: String,
    pub(crate) scope_review_status: String,
    pub(crate) hpo_inheritance_modes: Vec<String>,
    pub(crate) penetrance_statuses: Vec<String>,
    pub(crate) penetrance_assertions: Vec<Json>,
    pub(crate) penetrance_hpo_ids: Vec<String>,
    pub(crate) onset_hpo_ids: Vec<String>,
    pub(crate) hpo_assertions: Vec<HpoAssertion>,
    pub(crate) hpo_assertion_count: Json,
}

fn json_object(fields: Vec<(&str, Json)>) -> Json {
    Json::Object(
        fields
            .into_iter()
            .map(|(key, value)| (key.to_owned(), value))
            .collect(),
    )
}

fn json_texts(values: &[String]) -> Json {
    Json::Array(values.iter().cloned().map(Json::Str).collect())
}

impl ConditionContext {
    fn fields(&self) -> Vec<(&'static str, Json)> {
        vec![
            (
                "hpo_match_status",
                Json::Str("matched_gene_and_condition".to_owned()),
            ),
            ("hpo_disease_id", Json::Str(self.hpo_disease_id.clone())),
            ("mondo_id", Json::Str(self.mondo_id.clone())),
            ("disease_scope", Json::Str(self.disease_scope.clone())),
            ("priva_scope", Json::Str("include".to_owned())),
            (
                "scope_review_status",
                Json::Str(self.scope_review_status.clone()),
            ),
            (
                "hpo_inheritance_modes",
                json_texts(&self.hpo_inheritance_modes),
            ),
            ("penetrance_statuses", json_texts(&self.penetrance_statuses)),
            (
                "penetrance_assertions",
                Json::Array(self.penetrance_assertions.clone()),
            ),
            ("penetrance_hpo_ids", json_texts(&self.penetrance_hpo_ids)),
            ("onset_hpo_ids", json_texts(&self.onset_hpo_ids)),
            (
                "hpo_assertions",
                Json::Array(
                    self.hpo_assertions
                        .iter()
                        .map(HpoAssertion::to_json)
                        .collect(),
                ),
            ),
            ("hpo_assertion_count", self.hpo_assertion_count.clone()),
        ]
    }
}

/// `HPO_CACHE_INHERITANCE_LABELS.get(mode, mode)`.
fn inheritance_label(mode: String) -> String {
    match HPO_CACHE_INHERITANCE_LABELS
        .iter()
        .find(|(key, _)| *key == mode)
    {
        Some((_, label)) => (*label).to_owned(),
        None => mode,
    }
}

/// Keeps the first of each HPO assertion.
fn push_new(seen: &mut HashSet<HpoAssertion>, out: &mut Vec<HpoAssertion>, item: HpoAssertion) {
    if seen.insert(item.clone()) {
        out.push(item);
    }
}

/// `condition_cache_context(condition_key, condition)`: `None` for Python's
/// empty dict (a condition outside the germline `include` scope). `condition`
/// is a dict.
pub(crate) fn condition_cache_context(
    condition_key: &str,
    condition: &Json,
) -> Result<Option<ConditionContext>> {
    let Some(scope) = dict_or_none(condition.get("priva_scope")) else {
        return Ok(None);
    };
    if clean_json(scope.get("decision"))?.to_lowercase() != "include" {
        return Ok(None);
    }
    let identifiers = dict_or_none(condition.get("identifiers"));
    let mondo_id = match identifiers.and_then(|identifiers| identifiers.get("MONDO")) {
        Some(Json::Array(ids)) if !ids.is_empty() => clean_json(Some(&ids[0]))?,
        _ => String::new(),
    };
    let inheritance = dict_or_none(condition.get("inheritance"));
    let mut hpo_inheritance_modes = Vec::new();
    if let Some(modes) = inheritance.and_then(|inheritance| inheritance.get("modes")) {
        for mode in iterate(modes)? {
            let cleaned = clean_json(Some(mode))?;
            if !cleaned.is_empty() {
                hpo_inheritance_modes.push(inheritance_label(cleaned));
            }
        }
    }
    let penetrance = dict_or_none(condition.get("penetrance"));
    let onset = dict_or_none(condition.get("onset"));

    let mut hpo_assertions = Vec::new();
    let mut seen = HashSet::new();
    for axis in [inheritance, onset] {
        let records = dicts_or_empty(axis.and_then(|axis| axis.get("assertions")))?;
        for record in records {
            let item = HpoAssertion::from_fields(clean_json(record.get("hpo_id"))?, record)?;
            push_new(&mut seen, &mut hpo_assertions, item);
        }
    }
    let penetrance_assertions: Vec<Json> =
        dicts_or_empty(penetrance.and_then(|penetrance| penetrance.get("assertions")))?
            .into_iter()
            .cloned()
            .collect();
    let mut penetrance_hpo_ids: Vec<String> = Vec::new();
    for assertion in &penetrance_assertions {
        if clean_json(assertion.get("source"))?.to_uppercase() != "HPO" {
            continue;
        }
        for hpo_id in iterate_or_empty(assertion.get("hpo_ids"))? {
            let cleaned = clean_json(Some(hpo_id))?;
            if !cleaned.is_empty() && !penetrance_hpo_ids.contains(&cleaned) {
                penetrance_hpo_ids.push(cleaned);
            }
        }
    }
    for assertion in &penetrance_assertions {
        if clean_json(assertion.get("source"))?.to_uppercase() != "HPO" {
            continue;
        }
        for hpo_id in iterate_or_empty(assertion.get("hpo_ids"))? {
            // Unlike the list above, an empty HPO id is kept here, as in Python.
            let item = HpoAssertion::from_fields(clean_json(Some(hpo_id))?, assertion)?;
            push_new(&mut seen, &mut hpo_assertions, item);
        }
    }

    let hpo_disease_id = clean_text(condition_key);
    let disease_scope = clean_json(scope.get("category"))?;
    let scope_review_status = clean_json(scope.get("review_status"))?;
    let mut penetrance_statuses = Vec::new();
    for status in iterate_or_empty(penetrance.and_then(|penetrance| penetrance.get("statuses")))? {
        let cleaned = clean_json(Some(status))?;
        if !cleaned.is_empty() {
            penetrance_statuses.push(cleaned.to_lowercase());
        }
    }
    let mut onset_hpo_ids = Vec::new();
    for item in dicts_or_empty(onset.and_then(|onset| onset.get("assertions")))? {
        let cleaned = clean_json(item.get("hpo_id"))?;
        if !cleaned.is_empty() {
            onset_hpo_ids.push(cleaned);
        }
    }
    Ok(Some(ConditionContext {
        hpo_disease_id,
        mondo_id,
        disease_scope,
        scope_review_status,
        hpo_inheritance_modes,
        penetrance_statuses,
        penetrance_assertions,
        penetrance_hpo_ids,
        onset_hpo_ids,
        hpo_assertions,
        hpo_assertion_count: condition
            .get("hpo_assertion_count")
            .cloned()
            .unwrap_or(Json::Int(0)),
    }))
}

// ------------------------------------------------------------------------ the assertions

/// One record of `condition_cache_mechanism_assertions`: the condition's
/// context and the assertion's own fields.
#[derive(Clone, Debug)]
pub(crate) struct ConditionAssertion {
    pub(crate) context: Rc<ConditionContext>,
    pub(crate) source: String,
    pub(crate) source_record_id: String,
    pub(crate) source_condition_id: String,
    pub(crate) disease: String,
    pub(crate) mechanism: String,
    pub(crate) assertion_basis: String,
    pub(crate) mechanism_raw: String,
    pub(crate) allelic_requirement: String,
    pub(crate) confidence: String,
    pub(crate) mechanism_confidence: String,
    pub(crate) disease_confidence: String,
    pub(crate) pmids: Vec<Json>,
}

impl ConditionAssertion {
    /// The Python dict, for the harness (`json.dumps(..., sort_keys=True)`).
    pub(crate) fn to_json(&self) -> Json {
        let mut fields = self.context.fields();
        fields.extend([
            ("source", Json::Str(self.source.clone())),
            ("source_record_id", Json::Str(self.source_record_id.clone())),
            (
                "source_condition_id",
                Json::Str(self.source_condition_id.clone()),
            ),
            ("disease", Json::Str(self.disease.clone())),
            ("mechanism", Json::Str(self.mechanism.clone())),
            ("assertion_basis", Json::Str(self.assertion_basis.clone())),
            ("mechanism_raw", Json::Str(self.mechanism_raw.clone())),
            (
                "allelic_requirement",
                Json::Str(self.allelic_requirement.clone()),
            ),
            ("confidence", Json::Str(self.confidence.clone())),
            (
                "mechanism_confidence",
                Json::Str(self.mechanism_confidence.clone()),
            ),
            (
                "disease_confidence",
                Json::Str(self.disease_confidence.clone()),
            ),
            ("pmids", Json::Array(self.pmids.clone())),
        ]);
        json_object(fields)
    }
}

/// `_hpo_allelic_requirements`: HPO inheritance labels (joined by `;`) as the
/// compact allelic requirements, sorted as `sorted(set)` gives them.
fn hpo_allelic_requirements(inheritance_modes: &str) -> BTreeSet<String> {
    let modes: HashSet<String> = clean_text(inheritance_modes)
        .split(';')
        .map(clean_text)
        .filter(|mode| !mode.is_empty())
        .collect();
    let has = |label: &str| modes.contains(label);
    let mut requirements = BTreeSet::new();
    let mut add = |requirement: &str| {
        requirements.insert(requirement.to_owned());
    };
    if has("Autosomal recessive inheritance") {
        add("biallelic_autosomal");
    }
    if has("Pseudoautosomal recessive inheritance") {
        add("recessive");
    }
    if has("X-linked recessive inheritance") {
        add("x_linked_recessive");
    }
    if has("X-linked dominant inheritance") {
        add("x_linked_dominant");
    }
    if has("X-linked inheritance") {
        add("x_linked_unspecified");
    }
    if has("Autosomal dominant inheritance")
        || has("Autosomal dominant inheritance with maternal imprinting")
    {
        add("monoallelic_autosomal");
    }
    if has("Y-linked inheritance") {
        add("monoallelic_Y_hemizygous");
    }
    if has("Mitochondrial inheritance") {
        add("mitochondrial");
    }
    requirements
}

/// `(value or {}).values()`: a non-dict that is true is Python's AttributeError.
fn dict_values_or_empty(value: Option<&Json>) -> Result<Vec<&Json>> {
    match value {
        Some(Json::Object(items)) => Ok(items.iter().map(|(_, value)| value).collect()),
        Some(value) if truthy(value) => Err(python_error(
            "AttributeError",
            format!(
                "'{}' object has no attribute 'values'",
                json_type_name(value)
            ),
        )),
        _ => Ok(Vec::new()),
    }
}

/// `_condition_allelic_requirements(condition, hpo_inheritance_modes)`.
fn condition_allelic_requirements(
    condition: &Json,
    inheritance_modes: &str,
) -> Result<BTreeSet<String>> {
    let mut requirements = hpo_allelic_requirements(inheritance_modes);
    for block in dict_values_or_empty(condition.get("pathogenic_mechanisms"))? {
        if !matches!(block, Json::Object(_)) {
            continue;
        }
        for evidence in dicts_or_empty(block.get("evidence"))? {
            let requirement = clean_json(evidence.get("allelic_requirement"))?;
            if !requirement.is_empty() {
                requirements.insert(requirement);
            }
        }
    }
    Ok(requirements)
}

/// `condition_cache_mechanism_assertions(gene_record)`: every included
/// condition's mechanism assertions (one per allelic requirement), or one
/// `UNRESOLVED` history per requirement for a condition without a usable
/// mechanism; repeats dropped. `gene_record` is the gene's dict.
pub(crate) fn condition_cache_mechanism_assertions(
    gene_record: &Json,
) -> Result<Vec<ConditionAssertion>> {
    let mut assertions: Vec<ConditionAssertion> = Vec::new();
    let Some(Json::Object(conditions)) = py_get(gene_record, "conditions")? else {
        return Ok(assertions);
    };
    for (condition_key, condition) in sorted_items(conditions) {
        if !matches!(condition, Json::Object(_)) {
            continue;
        }
        let Some(context) = condition_cache_context(condition_key, condition)? else {
            continue;
        };
        let context = Rc::new(context);
        let joined_modes = context.hpo_inheritance_modes.join(";");
        let mut emitted = false;
        if let Some(Json::Object(blocks)) = condition.get("pathogenic_mechanisms") {
            for (mechanism_key, block) in sorted_items(blocks) {
                let mechanism = clean_text(mechanism_key).to_uppercase();
                if !CANONICAL_MECHANISMS.contains(&mechanism.as_str())
                    || !matches!(block, Json::Object(_))
                {
                    continue;
                }
                for evidence in dicts_or_empty(block.get("evidence"))? {
                    let source = clean_json(evidence.get("source"))?;
                    if !CONDITION_MECHANISM_SOURCES.contains(&source.as_str()) {
                        continue;
                    }
                    let identifiers: &[Json] = match evidence.get("condition_identifiers") {
                        Some(Json::Array(items)) => items,
                        _ => &[],
                    };
                    // `next(...)` stops at the first identifier that is not a MONDO id.
                    let mut source_condition_id = None;
                    for identifier in identifiers {
                        let cleaned = clean_json(Some(identifier))?;
                        if !cleaned.is_empty() && !cleaned.to_uppercase().starts_with("MONDO:") {
                            source_condition_id = Some(cleaned);
                            break;
                        }
                    }
                    let source_condition_id =
                        source_condition_id.unwrap_or_else(|| clean_text(condition_key));
                    let requirement = clean_json(evidence.get("allelic_requirement"))?;
                    let mut requirements: Vec<String> = if requirement.is_empty() {
                        condition_allelic_requirements(condition, &joined_modes)?
                            .into_iter()
                            .collect()
                    } else {
                        vec![requirement]
                    };
                    if requirements.is_empty() {
                        requirements.push(String::new());
                    }
                    for effective_requirement in requirements {
                        let source_record_id = clean_json(evidence.get("source_record_id"))?;
                        let mut disease = clean_json(evidence.get("condition_label"))?;
                        if disease.is_empty() {
                            disease = clean_json(condition.get("label"))?;
                        }
                        let assertion_basis = clean_json(evidence.get("assertion_basis"))?;
                        let mechanism_raw = clean_json(evidence.get("mechanism_raw"))?;
                        let confidence = clean_json(evidence.get("mechanism_confidence"))?;
                        let disease_confidence = clean_json(evidence.get("disease_confidence"))?;
                        let pmids = list_or_empty(evidence.get("pmids"))?;
                        assertions.push(ConditionAssertion {
                            context: Rc::clone(&context),
                            source: source.clone(),
                            source_record_id,
                            source_condition_id: source_condition_id.clone(),
                            disease,
                            mechanism: mechanism.clone(),
                            assertion_basis,
                            mechanism_raw,
                            allelic_requirement: effective_requirement,
                            mechanism_confidence: confidence.clone(),
                            confidence,
                            disease_confidence,
                            pmids,
                        });
                        emitted = true;
                    }
                }
            }
        }
        if emitted {
            continue;
        }
        let mut requirements: Vec<String> =
            condition_allelic_requirements(condition, &joined_modes)?
                .into_iter()
                .collect();
        if requirements.is_empty() {
            requirements.push(String::new());
        }
        for requirement in requirements {
            assertions.push(ConditionAssertion {
                context: Rc::clone(&context),
                source: "HPO_condition_history".to_owned(),
                source_record_id: clean_text(condition_key),
                source_condition_id: clean_text(condition_key),
                disease: clean_json(condition.get("label"))?,
                mechanism: UNRESOLVED_MECHANISM.to_owned(),
                assertion_basis: "condition_without_explicit_mechanism".to_owned(),
                mechanism_raw: String::new(),
                allelic_requirement: requirement,
                confidence: "condition_without_explicit_mechanism".to_owned(),
                mechanism_confidence: String::new(),
                disease_confidence: String::new(),
                pmids: Vec::new(),
            });
        }
    }
    // ASSERTION_IDENTITY_FIELDS, each through `_clean`.
    Ok(deduplicate_by(assertions, |a| {
        [
            &a.source,
            &a.source_record_id,
            &a.source_condition_id,
            &a.context.mondo_id,
            &a.disease,
            &a.mechanism,
            &a.allelic_requirement,
            &a.confidence,
        ]
        .map(|field| clean_text(field))
    }))
}

// -------------------------------------------------------------------- inheritance

/// `_INHERITANCE_LABEL_TO_KEY.get(mode, mode)` of a lower-cased mode.
fn inheritance_key(lowered: String) -> String {
    match HPO_CACHE_INHERITANCE_LABELS
        .iter()
        .find(|(_, label)| label.to_lowercase() == lowered)
    {
        Some((key, _)) => (*key).to_owned(),
        None => lowered,
    }
}

/// `normalize_inheritance(allelic_requirement, hpo_inheritance_modes)`: one
/// inheritance value and whether it is X-linked. `modes` are the items of the
/// Python list (a single text is a one-item list).
pub(crate) fn normalize_inheritance(allelic_requirement: &str, modes: &[String]) -> (String, bool) {
    let requirement = clean_text(allelic_requirement).to_lowercase();
    let modes: Vec<String> = modes
        .iter()
        .map(|mode| clean_text(mode))
        .filter(|mode| !mode.is_empty())
        .map(|mode| inheritance_key(mode.to_lowercase()))
        .collect();
    let x_linked = X_LINKED_INHERITANCE.contains(&requirement.as_str())
        || requirement.starts_with("monoallelic_x")
        || modes.iter().any(|mode| mode.starts_with("x_linked"));
    let mendelian: Vec<&String> = modes
        .iter()
        .filter(|mode| !NON_MENDELIAN_INHERITANCE.contains(&mode.as_str()))
        .collect();
    if mendelian.is_empty() {
        if let Some(mode) = modes
            .iter()
            .find(|mode| NON_MENDELIAN_INHERITANCE.contains(&mode.as_str()))
        {
            return (mode.clone(), x_linked);
        }
    }
    let answer = |value: &str, x: bool| (value.to_owned(), x);
    let r = requirement.as_str();
    if r == "recessive" || r.starts_with("biallelic_") {
        return answer("recessive", x_linked);
    }
    if r == "x_linked_recessive" {
        return answer("x_linked_recessive", true);
    }
    if r == "x_linked_dominant" || r == "monoallelic_x_heterozygous" {
        return answer("x_linked_dominant", true);
    }
    if r == "x_linked_unspecified" || r.starts_with("monoallelic_x") {
        return answer("x_linked_unspecified", true);
    }
    if r == "dominant" {
        return answer("dominant", x_linked);
    }
    if r == "mitochondrial" {
        return answer("mitochondrial", false);
    }
    if r == "y_linked" || r.starts_with("monoallelic_y") {
        return answer("y_linked", false);
    }
    if r.starts_with("monoallelic_") {
        return answer("dominant", x_linked);
    }
    for mode in mendelian {
        match mode.as_str() {
            "autosomal_recessive" | "pseudoautosomal_recessive" => {
                return answer("recessive", x_linked);
            }
            "x_linked_recessive" => return answer("x_linked_recessive", true),
            "x_linked_dominant" => return answer("x_linked_dominant", true),
            "x_linked" => return answer("x_linked_unspecified", true),
            "autosomal_dominant" | "autosomal_dominant_maternal_imprinting" => {
                return answer("dominant", x_linked);
            }
            "mitochondrial" => return answer("mitochondrial", false),
            "y_linked" => return answer("y_linked", false),
            _ => {}
        }
    }
    (String::new(), x_linked)
}

/// `normalize_inheritance_histories(allelic_requirement, hpo_inheritance_modes)`:
/// every inheritance category the requirement and each mode state, repeats
/// dropped, or one fallback pair.
pub(crate) fn normalize_inheritance_histories(
    allelic_requirement: &str,
    modes: &[String],
) -> Vec<(String, bool)> {
    let mut normalized = Vec::new();
    if !clean_text(allelic_requirement).is_empty() {
        normalized.push(normalize_inheritance(allelic_requirement, &[]));
    }
    for mode in modes {
        if !clean_text(mode).is_empty() {
            normalized.push(normalize_inheritance("", std::slice::from_ref(mode)));
        }
    }
    let mut values: Vec<(String, bool)> = Vec::new();
    for value in normalized {
        if !value.0.is_empty() && !values.contains(&value) {
            values.push(value);
        }
    }
    if !values.is_empty() {
        return values;
    }
    let fallback = normalize_inheritance(allelic_requirement, modes);
    if !fallback.0.is_empty() || fallback.1 {
        vec![fallback]
    } else {
        vec![(String::new(), false)]
    }
}

/// `gene_has_lof_mechanism_history(gene_record)`: whether an included
/// condition of the gene has a curated LOF assertion (ClinGen's only with
/// sufficient evidence). Python reads the dicts without type checks, so the
/// errors it raises on other kinds of values are raised here too.
pub(crate) fn gene_has_lof_mechanism_history(gene_record: &Json) -> Result<bool> {
    let missing_get = |value: &Json| {
        python_error(
            "AttributeError",
            format!("'{}' object has no attribute 'get'", json_type_name(value)),
        )
    };
    // `(x or {}).get(key)`
    let or_empty_get = |value: Option<&Json>, key: &str| -> Result<Option<Json>> {
        match value {
            Some(value) if truthy(value) => Ok(py_get(value, key)?.cloned()),
            _ => Ok(None),
        }
    };
    for condition in dict_values_or_empty(py_get(gene_record, "conditions")?)? {
        let decision = or_empty_get(py_get(condition, "priva_scope")?, "decision")?;
        if !matches!(&decision, Some(Json::Str(text)) if text == "include") {
            continue;
        }
        let Some(lof) = or_empty_get(py_get(condition, "pathogenic_mechanisms")?, "LOF")? else {
            continue;
        };
        if !truthy(&lof) {
            continue;
        }
        let evidence_list = py_get(&lof, "evidence")?;
        let items: &[Json] = match evidence_list {
            Some(value) if truthy(value) => match value {
                Json::Array(items) => items,
                // The characters of a text and the keys of a dict are texts.
                Json::Str(_) | Json::Object(_) => {
                    return Err(missing_get(&Json::Str(String::new())));
                }
                other => {
                    return Err(python_error(
                        "TypeError",
                        format!("'{}' object is not iterable", json_type_name(other)),
                    ));
                }
            },
            _ => &[],
        };
        for evidence in items {
            let basis = py_get(evidence, "assertion_basis")?;
            if !matches!(basis, Some(Json::Str(text)) if text == "curated") {
                continue;
            }
            if clean_json(py_get(evidence, "source")?)? == "ClinGen_haploinsufficiency"
                && clean_json(py_get(evidence, "mechanism_confidence")?)? != "sufficient_evidence"
            {
                continue;
            }
            return Ok(true);
        }
    }
    Ok(false)
}

/// `gene_inheritance_from_constraint(symbol, clingen=..., loeuf=..., chromosome=...)`:
/// the inheritance read from the gene's ClinGen haploinsufficiency score and
/// LOEUF when no condition states one. `haploinsufficiency_score` is the
/// gene's ClinGen score (`None` without a record or a score), `loeuf` its LOEUF.
pub(crate) fn gene_inheritance_from_constraint(
    haploinsufficiency_score: Option<i128>,
    loeuf: Option<f64>,
    chromosome: Option<Cell<'_>>,
) -> &'static str {
    let lowered = clean(chromosome).to_lowercase();
    let normalized = lowered.strip_prefix("chr").unwrap_or(&lowered);
    let haploinsufficient = haploinsufficiency_score == Some(3);
    if normalized == "x" {
        return if haploinsufficient {
            "x_linked_dominant"
        } else {
            "x_linked_unspecified"
        };
    }
    if normalized == "y" {
        return if haploinsufficient { "y_linked" } else { "" };
    }
    let constrained = loeuf.is_some_and(|value| !value.is_nan() && value < 0.35);
    if haploinsufficient || constrained {
        "dominant"
    } else {
        "recessive"
    }
}

/// `normalize_penetrance(values)` of a list of texts.
pub(crate) fn normalize_penetrance(values: &[String]) -> &'static str {
    let values: Vec<&str> = values.iter().map(String::as_str).collect();
    normalize_penetrance_values(Values::List(&values))
}

// ------------------------------------------------------------------------------- harness

/// The modes argument the capture script passes (JSON): a text is one mode,
/// `None` or another false value none, a list its items.
fn modes_from_json(value: &Json) -> Result<Vec<String>> {
    match value {
        Json::Str(text) => Ok(vec![text.clone()]),
        value if !truthy(value) => Ok(Vec::new()),
        Json::Array(items) => items.iter().map(|item| clean_json(Some(item))).collect(),
        Json::Object(_) => Err(refused("modes given as a dict")),
        other => Err(python_error(
            "TypeError",
            format!("'{}' object is not iterable", json_type_name(other)),
        )),
    }
}

fn show_pair((inheritance, x_linked): &(String, bool)) -> Json {
    Json::Array(vec![Json::Str(inheritance.clone()), Json::Bool(*x_linked)])
}

/// The ClinGen scores and LOEUF values a `gene_inheritance_from_constraint` case gives.
fn constraint_options(
    options: &Json,
) -> Result<(HashMap<String, Option<i128>>, HashMap<String, f64>)> {
    let mut clingen = HashMap::new();
    if let Some(Json::Object(items)) = options.get("clingen") {
        for (symbol, score) in items {
            let score = match score {
                Json::Null => None,
                Json::Int(number) => Some(*number),
                _ => return Err(invalid("ClinGen scores must be integers or null")),
            };
            clingen.insert(symbol.clone(), score);
        }
    }
    let mut loeuf = HashMap::new();
    if let Some(Json::Object(items)) = options.get("loeuf") {
        for (symbol, value) in items {
            let value = match value {
                Json::Float(number) => *number,
                Json::Int(number) => *number as f64,
                _ => return Err(invalid("LOEUF values must be numbers")),
            };
            loeuf.insert(symbol.clone(), value);
        }
    }
    Ok((clingen, loeuf))
}

/// `priva acmg-step gene_mechanism_conditions <function>`: each row's inputs
/// (`gene_record`, `condition_key`/`condition`, `requirement`/`modes`,
/// `symbol`/`chromosome`, `values`; JSON as text) given to the function, the
/// result shown as the capture script shows Python's. For
/// `condition_cache_mechanism_assertions` the output has one row per assertion
/// (`row`, `assertion`), or one row with the exception.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    const FUNCTIONS: [&str; 7] = [
        "condition_cache_mechanism_assertions",
        "gene_has_lof_mechanism_history",
        "condition_cache_context",
        "normalize_inheritance",
        "normalize_inheritance_histories",
        "gene_inheritance_from_constraint",
        "normalize_penetrance",
    ];
    if !FUNCTIONS.contains(&function) {
        return Err(not_converted("gene_mechanism_conditions", function));
    }
    let table = &call.table;
    if function == "condition_cache_mechanism_assertions" {
        let mut rows = Vec::new();
        let mut shown = Vec::new();
        for row in 0..table.rows {
            let record = json_column_value(table, "gene_record", row)?;
            let outcome = condition_cache_mechanism_assertions(&record).and_then(|assertions| {
                assertions
                    .iter()
                    .map(|assertion| show_json(&assertion.to_json()))
                    .collect::<Result<Vec<String>>>()
            });
            match outcome {
                Ok(texts) => {
                    for text in texts {
                        rows.push(row as i64);
                        shown.push(Some(text));
                    }
                }
                Err(error) => {
                    rows.push(row as i64);
                    shown.push(Some(show_outcome(Err(error))?));
                }
            }
        }
        return Ok(Table {
            names: vec!["row".to_owned(), "assertion".to_owned()],
            rows: rows.len(),
            columns: vec![Column::Int(rows), Column::Text(shown)],
        });
    }
    let constraint = if function == "gene_inheritance_from_constraint" {
        Some(constraint_options(&call.options)?)
    } else {
        None
    };
    let mut results = Vec::with_capacity(table.rows);
    for row in 0..table.rows {
        let outcome: Result<String> = match function {
            "gene_has_lof_mechanism_history" => (|| {
                let record = json_column_value(table, "gene_record", row)?;
                Ok(format!(
                    "b:{}",
                    super::gene_mechanism_common::python_bool(gene_has_lof_mechanism_history(
                        &record
                    )?)
                ))
            })(),
            "condition_cache_context" => (|| {
                let key = clean(Some(table.column("condition_key")?.cell(row)));
                let condition = json_column_value(table, "condition", row)?;
                if !matches!(condition, Json::Object(_)) {
                    return Err(invalid("condition must be a JSON object"));
                }
                match condition_cache_context(&key, &condition)? {
                    Some(context) => show_json(&json_object(context.fields())),
                    None => Ok("j:{}".to_owned()),
                }
            })(),
            "normalize_inheritance" | "normalize_inheritance_histories" => (|| {
                let requirement = clean(Some(table.column("requirement")?.cell(row)));
                let modes = modes_from_json(&json_column_value(table, "modes", row)?)?;
                if function == "normalize_inheritance" {
                    show_json(&show_pair(&normalize_inheritance(&requirement, &modes)))
                } else {
                    show_json(&Json::Array(
                        normalize_inheritance_histories(&requirement, &modes)
                            .iter()
                            .map(show_pair)
                            .collect(),
                    ))
                }
            })(),
            "gene_inheritance_from_constraint" => (|| {
                let (clingen, loeuf) = constraint.as_ref().expect("read above");
                let symbol = match table.column("symbol")?.cell(row) {
                    Cell::Text(text) => text.to_owned(),
                    _ => return Err(invalid("symbol must be text")),
                };
                let chromosome = table.column("chromosome")?.cell(row);
                Ok(show_text(gene_inheritance_from_constraint(
                    clingen.get(&symbol).copied().flatten(),
                    loeuf.get(&symbol).copied(),
                    Some(chromosome),
                )))
            })(),
            "normalize_penetrance" => (|| {
                let values = json_column_value(table, "values", row)?;
                let result = match &values {
                    Json::Str(text) => normalize_penetrance_values(Values::Text(text)),
                    Json::Array(items) => {
                        let texts: Vec<String> = items
                            .iter()
                            .map(|item| match item {
                                Json::Str(text) => Ok(text.clone()),
                                _ => Err(invalid("penetrance values must be texts")),
                            })
                            .collect::<Result<_>>()?;
                        normalize_penetrance(&texts)
                    }
                    _ => return Err(invalid("values must be a text or a list")),
                };
                Ok(show_text(result))
            })(),
            _ => return Err(not_converted("gene_mechanism_conditions", function)),
        };
        results.push(show_outcome(outcome)?);
    }
    Ok(result_table(results))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::json_text::parse;

    fn texts(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| (*value).to_owned()).collect()
    }

    #[test]
    fn inheritance_vocabularies_meet() {
        assert_eq!(
            normalize_inheritance("biallelic_autosomal", &[]),
            ("recessive".to_owned(), false)
        );
        assert_eq!(
            normalize_inheritance("monoallelic_X_hemizygous", &[]),
            ("x_linked_unspecified".to_owned(), true)
        );
        assert_eq!(
            normalize_inheritance("", &texts(&["X-linked dominant inheritance"])),
            ("x_linked_dominant".to_owned(), true)
        );
        assert_eq!(
            normalize_inheritance("", &texts(&["polygenic", "Autosomal dominant inheritance"])),
            ("dominant".to_owned(), false)
        );
        assert_eq!(
            normalize_inheritance("", &texts(&["Polygenic inheritance"])),
            ("polygenic".to_owned(), false)
        );
        assert_eq!(
            normalize_inheritance_histories(
                "",
                &texts(&[
                    "Autosomal dominant inheritance",
                    "polygenic",
                    "autosomal_dominant"
                ])
            ),
            vec![
                ("dominant".to_owned(), false),
                ("polygenic".to_owned(), false)
            ]
        );
        assert_eq!(
            normalize_inheritance_histories("monoallelic_x_unknown", &[]),
            vec![("x_linked_unspecified".to_owned(), true)]
        );
        assert_eq!(
            normalize_inheritance_histories("nan", &texts(&["odd"])),
            vec![(String::new(), false)]
        );
    }

    #[test]
    fn constraint_fallback() {
        assert_eq!(
            gene_inheritance_from_constraint(Some(3), None, Some(Cell::Text("chrX"))),
            "x_linked_dominant"
        );
        assert_eq!(
            gene_inheritance_from_constraint(None, Some(0.2), Some(Cell::Text("chrY"))),
            ""
        );
        assert_eq!(
            gene_inheritance_from_constraint(Some(30), Some(0.2), Some(Cell::Int(7))),
            "dominant"
        );
        assert_eq!(
            gene_inheritance_from_constraint(None, Some(f64::NAN), None),
            "recessive"
        );
    }

    #[test]
    fn assertions_from_a_cache_record() {
        let gene = parse(
            r#"{"conditions": {
                "OMIM:2": {"label": "B", "priva_scope": {"decision": "Include", "category": "c"},
                           "identifiers": {"MONDO": ["MONDO:1"]},
                           "inheritance": {"modes": ["autosomal_recessive", "x_linked"],
                                           "assertions": [{"hpo_id": "HP:0000007"}, {"hpo_id": "HP:0000007"}]},
                           "penetrance": {"statuses": [" Incomplete "],
                                          "assertions": [{"source": "hpo", "hpo_ids": ["HP:0003829", ""]}]},
                           "pathogenic_mechanisms": {
                               "LOF": {"evidence": [
                                   {"source": "G2P_DDG2P", "allelic_requirement": "",
                                    "condition_identifiers": ["MONDO:9", "OMIM:22"],
                                    "assertion_basis": "curated", "pmids": ["1"]},
                                   {"source": "Unknown"}]},
                               "TOXIC": {"evidence": [{"source": "Orphadata"}]}}},
                "OMIM:1": {"label": "A", "priva_scope": {"decision": "include"},
                           "inheritance": {"modes": []}},
                "OMIM:3": {"priva_scope": {"decision": "review"}}}}"#,
            false,
        )
        .unwrap();
        let assertions = condition_cache_mechanism_assertions(&gene).unwrap();
        let summary: Vec<(String, String, String, String)> = assertions
            .iter()
            .map(|a| {
                (
                    a.context.hpo_disease_id.clone(),
                    a.mechanism.clone(),
                    a.allelic_requirement.clone(),
                    a.source_condition_id.clone(),
                )
            })
            .collect();
        assert_eq!(
            summary,
            vec![
                (
                    "OMIM:1".into(),
                    "UNRESOLVED".into(),
                    "".into(),
                    "OMIM:1".into()
                ),
                (
                    "OMIM:2".into(),
                    "LOF".into(),
                    "biallelic_autosomal".into(),
                    "OMIM:22".into()
                ),
                (
                    "OMIM:2".into(),
                    "LOF".into(),
                    "x_linked_unspecified".into(),
                    "OMIM:22".into()
                ),
            ]
        );
        let context = &assertions[1].context;
        assert_eq!(context.mondo_id, "MONDO:1");
        assert_eq!(
            context.hpo_inheritance_modes,
            ["Autosomal recessive inheritance", "X-linked inheritance"]
        );
        assert_eq!(context.penetrance_statuses, ["incomplete"]);
        assert_eq!(context.penetrance_hpo_ids, ["HP:0003829"]);
        assert_eq!(context.hpo_assertions.len(), 3);
        assert!(gene_has_lof_mechanism_history(&gene).is_ok_and(|found| !found));
    }

    #[test]
    fn python_errors_on_odd_kinds() {
        let gene = parse(
            r#"{"conditions": {"A": {"priva_scope": {"decision": "include"}, "inheritance": {"modes": null}}}}"#,
            false,
        )
        .unwrap();
        assert_eq!(
            condition_cache_mechanism_assertions(&gene)
                .unwrap_err()
                .to_string(),
            "TypeError: 'NoneType' object is not iterable"
        );
        let gene = parse(
            r#"{"conditions": {"A": {"priva_scope": {"decision": "include"},
                "pathogenic_mechanisms": {"LOF": {"evidence": [{"assertion_basis": "curated",
                "source": "ClinGen_haploinsufficiency", "mechanism_confidence": "sufficient_evidence"}]}}},
                "B": []}}"#,
            false,
        )
        .unwrap();
        assert!(gene_has_lof_mechanism_history(&gene).unwrap());
        let gene = parse(r#"{"conditions": {"B": [1]}}"#, false).unwrap();
        assert_eq!(
            gene_has_lof_mechanism_history(&gene)
                .unwrap_err()
                .to_string(),
            "AttributeError: 'list' object has no attribute 'get'"
        );
        assert_eq!(
            normalize_penetrance(&texts(&["incomplete", "complete"])),
            "incomplete"
        );
    }
}
