//! `validate_cache_payload`, `_cache_counts`, `load_and_validate_cache` and
//! `validate_cache_provenance` over the parsed JSON, with the exceptions the
//! script raises for a value of another shape.
use super::build::CANONICAL_MECHANISMS;
use super::*;

const SCHEMA_VERSION: &str = "1.2";
const PENETRANCE_ASSERTION_KEYS: [&str; 10] = [
    "source",
    "source_record_id",
    "condition_identifiers",
    "raw_values",
    "normalized_value",
    "hpo_ids",
    "frequency",
    "evidence",
    "reference",
    "evidence_url",
];

fn invalid_cache(message: String) -> Box<dyn std::error::Error> {
    error("ValueError", message)
}

/// `value[key]` for a text key.
fn index<'a>(value: &'a Json, key: &str) -> Result<&'a Json> {
    match value {
        Json::Object(_) => value
            .get(key)
            .ok_or_else(|| error("KeyError", crate::python_path::repr(key))),
        Json::Array(_) => Err(error(
            "TypeError",
            "list indices must be integers or slices, not str",
        )),
        Json::Str(_) => Err(error("TypeError", "string indices must be integers, not 'str'")),
        other => Err(error(
            "TypeError",
            format!("'{}' object is not subscriptable", type_name(other)),
        )),
    }
}

/// `len(value)`.
fn length(value: &Json) -> Result<i128> {
    match value {
        Json::Array(items) => Ok(items.len() as i128),
        Json::Object(items) => Ok(items.len() as i128),
        Json::Str(text) => Ok(text.chars().count() as i128),
        other => Err(error(
            "TypeError",
            format!("object of type '{}' has no len()", type_name(other)),
        )),
    }
}

/// `value.get(key, default)`.
fn get_or(value: &Json, key: &str, default: Json) -> Result<Json> {
    Ok(get(value, key)?.cloned().unwrap_or(default))
}

/// `value.items()` / `value.values()`.
fn entries<'a>(value: &'a Json, method: &str) -> Result<&'a Vec<(String, Json)>> {
    match value {
        Json::Object(items) => Ok(items),
        other => Err(error(
            "AttributeError",
            format!("'{}' object has no attribute '{method}'", type_name(other)),
        )),
    }
}

/// Iterating a value as Python does.
fn iterate(value: &Json) -> Result<Vec<Json>> {
    match value {
        Json::Array(items) => Ok(items.clone()),
        Json::Object(items) => Ok(items.iter().map(|(k, _)| Json::Str(k.clone())).collect()),
        Json::Str(text) => Ok(text.chars().map(|c| Json::Str(c.to_string())).collect()),
        other => Err(error(
            "TypeError",
            format!("'{}' object is not iterable", type_name(other)),
        )),
    }
}

/// Adding a value to a set: a list or dict is unhashable.
fn hashable(value: &Json) -> Result<()> {
    match value {
        Json::Array(_) | Json::Object(_) => Err(error(
            "TypeError",
            format!("unhashable type: '{}'", type_name(value)),
        )),
        _ => Ok(()),
    }
}

fn is_text(value: &Json, text: &str) -> bool {
    matches!(value, Json::Str(have) if have == text)
}

fn sum_over(values: impl Iterator<Item = Result<i128>>) -> Result<Json> {
    let mut total = 0;
    for value in values {
        total += value?;
    }
    Ok(Json::Int(total))
}

/// `_cache_counts(genes)` for a mapping the payload check has accepted.
pub(super) fn cache_counts(genes: &Json) -> Result<Json> {
    let genes = entries(genes, "values")?;
    let mut conditions = Vec::new();
    for (_, gene) in genes {
        for (_, condition) in entries(index(gene, "conditions")?, "values")? {
            conditions.push(condition);
        }
    }
    let mut resolved = 0;
    for condition in &conditions {
        resolved += i128::from(truthy(index(condition, "pathogenic_mechanisms")?));
    }
    let mut variants = 0;
    for condition in &conditions {
        for (_, block) in entries(index(condition, "pathogenic_mechanisms")?, "values")? {
            variants += length(index(block, "variants")?)?;
        }
    }
    let with = |axis: &str, values: &str| -> Result<Json> {
        sum_over(
            conditions
                .iter()
                .map(|c| Ok(i128::from(truthy(index(index(c, axis)?, values)?)))),
        )
    };
    Ok(Json::Object(vec![
        ("genes".into(), Json::Int(genes.len() as i128)),
        (
            "genes_with_hpo_conditions".into(),
            sum_over(
                genes
                    .iter()
                    .map(|(_, g)| Ok(i128::from(truthy(index(g, "conditions")?)))),
            )?,
        ),
        (
            "excluded_conditions".into(),
            sum_over(genes.iter().map(|(_, g)| {
                let gate = get(g, "scope_gate")?.cloned().unwrap_or(Json::Null);
                let gate = if truthy(&gate) { gate } else { Json::Object(vec![]) };
                match get_or(&gate, "excluded_condition_count", Json::Int(0))? {
                    Json::Int(n) => Ok(n),
                    Json::Bool(b) => Ok(i128::from(b)),
                    other => Err(error(
                        "TypeError",
                        format!(
                            "int() argument must be a string, a bytes-like object or a real number, not '{}'",
                            type_name(&other)
                        ),
                    )),
                }
            }))?,
        ),
        ("conditions".into(), Json::Int(conditions.len() as i128)),
        ("conditions_with_inheritance".into(), with("inheritance", "modes")?),
        ("conditions_with_penetrance".into(), with("penetrance", "statuses")?),
        (
            "conditions_with_pathogenic_mechanisms".into(),
            Json::Int(resolved),
        ),
        ("condition_variant_links".into(), Json::Int(variants)),
        (
            "unmapped_mechanisms".into(),
            sum_over(
                genes
                    .iter()
                    .map(|(_, g)| length(index(index(g, "unmapped_evidence")?, "mechanisms")?)),
            )?,
        ),
        (
            "unmapped_penetrance_assertions".into(),
            sum_over(genes.iter().map(|(_, g)| {
                length(&get_or(
                    index(g, "unmapped_evidence")?,
                    "penetrance",
                    Json::Array(vec![]),
                )?)
            }))?,
        ),
        (
            "unmapped_variants".into(),
            sum_over(
                genes
                    .iter()
                    .map(|(_, g)| length(index(index(g, "unmapped_evidence")?, "variants")?)),
            )?,
        ),
    ]))
}

fn validate_condition(symbol: &str, key: &str, condition: &Json) -> Result<()> {
    let scope = get(condition, "priva_scope")?.cloned().unwrap_or(Json::Null);
    if matches!(scope, Json::Object(_))
        && clean_json(&get(&scope, "decision")?.cloned().unwrap_or(Json::Null)).to_lowercase()
            == "exclude"
    {
        return Err(invalid_cache(format!(
            "{symbol}/{key} survived the final condition scope gate"
        )));
    }
    let identifiers = get_or(condition, "identifiers", Json::Object(vec![]))?;
    let mut present = false;
    for (_, values) in entries(&identifiers, "values")? {
        for identifier in iterate(values)? {
            hashable(&identifier)?;
            present |= is_text(&identifier, key);
        }
    }
    if !present {
        return Err(invalid_cache(format!(
            "{symbol}/{key} canonical key is absent from identifiers"
        )));
    }
    let penetrance = get(condition, "penetrance")?.cloned().unwrap_or(Json::Null);
    if !matches!(penetrance, Json::Object(_)) {
        return Err(invalid_cache(format!(
            "{symbol}/{key} lacks penetrance axis"
        )));
    }
    let statuses = get(&penetrance, "statuses")?.cloned().unwrap_or(Json::Null);
    let assertions = get(&penetrance, "assertions")?.cloned().unwrap_or(Json::Null);
    let Json::Array(statuses) = statuses else {
        return Err(invalid_cache(format!(
            "{symbol}/{key} has unsupported penetrance statuses"
        )));
    };
    for status in &statuses {
        hashable(status)?;
    }
    if !statuses
        .iter()
        .all(|s| is_text(s, "complete") || is_text(s, "incomplete"))
    {
        return Err(invalid_cache(format!(
            "{symbol}/{key} has unsupported penetrance statuses"
        )));
    }
    let Json::Array(assertions) = assertions else {
        return Err(invalid_cache(format!(
            "{symbol}/{key} penetrance assertions must be a list"
        )));
    };
    for assertion in &assertions {
        let schema_ok = match assertion {
            Json::Object(items) => {
                items.len() == PENETRANCE_ASSERTION_KEYS.len()
                    && PENETRANCE_ASSERTION_KEYS
                        .iter()
                        .all(|name| items.iter().any(|(k, _)| k == name))
            }
            _ => false,
        };
        if !schema_ok {
            return Err(invalid_cache(format!(
                "{symbol}/{key} has an invalid penetrance assertion schema"
            )));
        }
        let normalized = clean_json(&get_or(assertion, "normalized_value", Json::Null)?).to_lowercase();
        if normalized != "complete" && normalized != "incomplete" {
            return Err(invalid_cache(format!(
                "{symbol}/{key} has invalid normalized penetrance"
            )));
        }
        if !statuses.iter().any(|s| is_text(s, &normalized)) {
            return Err(invalid_cache(format!(
                "{symbol}/{key} penetrance assertion is absent from statuses"
            )));
        }
    }
    let mechanisms = get_or(condition, "pathogenic_mechanisms", Json::Object(vec![]))?;
    for (mechanism, block) in entries(&mechanisms, "items")? {
        if !CANONICAL_MECHANISMS.contains(&mechanism.as_str()) {
            return Err(invalid_cache(format!(
                "{symbol}/{key} has unsupported mechanism {mechanism}"
            )));
        }
        let variants = get_or(block, "variants", Json::Object(vec![]))?;
        for (variant_key, variant) in entries(&variants, "items")? {
            let link = get_or(variant, "condition_link", Json::Object(vec![]))?;
            let status = get_or(&link, "status", Json::Null)?;
            if !is_text(&status, "exact")
                || !is_text(&get_or(&link, "condition_key", Json::Null)?, key)
            {
                return Err(invalid_cache(format!(
                    "{symbol}/{key}/{variant_key} lacks an exact link"
                )));
            }
        }
    }
    Ok(())
}

/// `validate_cache_payload`: the counts when the payload is valid.
pub(super) fn validate_payload(payload: &Json) -> Result<Json> {
    if !matches!(payload, Json::Object(_)) {
        return Err(invalid_cache("Cache payload must be a JSON object".into()));
    }
    let meta = payload.get("_meta").cloned().unwrap_or(Json::Null);
    let genes = payload.get("genes").cloned().unwrap_or(Json::Null);
    if !matches!(meta, Json::Object(_))
        || !is_text(meta.get("schema_version").unwrap_or(&Json::Null), SCHEMA_VERSION)
    {
        return Err(invalid_cache(format!(
            "Cache schema_version must be {SCHEMA_VERSION}"
        )));
    }
    let Json::Object(gene_items) = &genes else {
        return Err(invalid_cache("Cache must contain at least one gene".into()));
    };
    if gene_items.is_empty() {
        return Err(invalid_cache("Cache must contain at least one gene".into()));
    }
    for (symbol, gene) in gene_items {
        if symbol.is_empty() || !matches!(gene, Json::Object(_)) {
            return Err(invalid_cache("Cache contains an invalid gene record".into()));
        }
        let conditions = gene.get("conditions").cloned().unwrap_or(Json::Null);
        let unmapped = gene.get("unmapped_evidence").cloned().unwrap_or(Json::Null);
        if !matches!(conditions, Json::Object(_)) || !matches!(unmapped, Json::Object(_)) {
            return Err(invalid_cache(format!(
                "{symbol} is missing condition or unmapped sections"
            )));
        }
        let gate = gene.get("scope_gate").cloned().unwrap_or(Json::Null);
        if !matches!(gate, Json::Object(_)) {
            return Err(invalid_cache(format!(
                "{symbol} is missing its scope gate record"
            )));
        }
        if !matches!(gate.get("excluded_condition_count"), Some(Json::Int(n)) if *n >= 0) {
            return Err(invalid_cache(format!(
                "{symbol} has an invalid excluded condition count"
            )));
        }
        for (key, condition) in entries(&conditions, "items")? {
            validate_condition(symbol, key, condition)?;
        }
        let variants = get_or(&unmapped, "variants", Json::Object(vec![]))?;
        for (variant_key, variant) in entries(&variants, "items")? {
            let link = get_or(variant, "condition_link", Json::Object(vec![]))?;
            if is_text(&get_or(&link, "status", Json::Null)?, "exact") {
                return Err(invalid_cache(format!(
                    "{symbol}/{variant_key} is exact but stored as unmapped"
                )));
            }
        }
        for assertion in iterate(&get_or(&unmapped, "penetrance", Json::Array(vec![]))?)? {
            let complete = match &assertion {
                Json::Object(items) => PENETRANCE_ASSERTION_KEYS
                    .iter()
                    .all(|name| items.iter().any(|(k, _)| k == name)),
                _ => false,
            };
            if !complete {
                return Err(invalid_cache(format!(
                    "{symbol} has an invalid unmapped penetrance assertion"
                )));
            }
        }
    }
    let actual = cache_counts(&genes)?;
    let recorded = meta.get("counts").cloned().unwrap_or(Json::Null);
    if !equal(&recorded, &actual) {
        return Err(invalid_cache(format!(
            "Cache count mismatch: recorded={}, actual={}",
            repr(&recorded),
            repr(&actual)
        )));
    }
    let statistics = meta.get("build_statistics").cloned().unwrap_or(Json::Null);
    let statistics = if truthy(&statistics) {
        statistics
    } else {
        Json::Object(vec![])
    };
    let scope_stats = get(&statistics, "scope_gate")?.cloned().unwrap_or(Json::Null);
    if !matches!(scope_stats, Json::Object(_)) {
        return Err(invalid_cache(
            "Cache metadata lacks final condition scope-gate statistics".into(),
        ));
    }
    let recorded_excluded = scope_stats
        .get("conditions_excluded")
        .cloned()
        .unwrap_or(Json::Null);
    let actual_excluded = actual.get("excluded_conditions").cloned().unwrap_or(Json::Null);
    if !equal(&recorded_excluded, &actual_excluded) {
        return Err(invalid_cache(format!(
            "Cache scope-gate count mismatch: recorded={}, actual={}",
            py_str(&recorded_excluded),
            py_str(&actual_excluded)
        )));
    }
    Ok(actual)
}

/// `load_and_validate_cache`.
pub(super) fn load_and_validate(path: &Path) -> Result<Json> {
    validate_payload(&shared::load_json(path)?)
}

/// `validate_cache_provenance`.
pub(super) fn validate_provenance(
    path: &Path,
    sources: &[(&str, PathBuf)],
    hpo_release: &str,
    mondo_release: &str,
) -> Result<Json> {
    let payload = shared::load_json(path)?;
    let counts = validate_payload(&payload)?;
    let meta = payload.get("_meta").cloned().unwrap_or(Json::Null);
    let recorded_sources = get_or(&meta, "sources", Json::Object(vec![]))?;
    let mut mismatches = Vec::new();
    for (name, source) in sources {
        let actual = source_metadata(source)?;
        let recorded = get(&recorded_sources, name)?.cloned().unwrap_or(Json::Null);
        if !matches!(recorded, Json::Object(_)) {
            mismatches.push(format!("{name}: missing provenance"));
            continue;
        }
        let changed: Vec<&str> = ["path", "size_bytes", "sha256"]
            .into_iter()
            .filter(|field| {
                !equal(
                    recorded.get(field).unwrap_or(&Json::Null),
                    actual.get(field).unwrap_or(&Json::Null),
                )
            })
            .collect();
        if !changed.is_empty() {
            mismatches.push(format!("{name}: changed {}", changed.join(",")));
        }
    }
    let releases = get_or(&meta, "releases", Json::Object(vec![]))?;
    for (name, expected) in [("HPO", clean(hpo_release)), ("MONDO", clean(mondo_release))] {
        let recorded = get(&releases, name)?.cloned().unwrap_or(Json::Null);
        if !is_text(&recorded, &expected) {
            mismatches.push(format!(
                "{name} release: recorded={}, configured={}",
                repr(&recorded),
                crate::python_path::repr(&expected)
            ));
        }
    }
    if !mismatches.is_empty() {
        return Err(invalid_cache(format!(
            "Cache provenance mismatch: {}",
            mismatches.join("; ")
        )));
    }
    Ok(counts)
}
