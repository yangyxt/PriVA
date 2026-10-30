//! Rust version of `scripts/hpo_penetrance.py`, the binary penetrance
//! vocabulary shared by PriVA's observation-based ACMG criteria.
//!
//! A condition is `complete` (explicit complete or high penetrance) or
//! `incomplete` (explicit incomplete, moderate or low penetrance, or late
//! onset; incomplete wins over complete); an empty string means no usable
//! penetrance assertion. Age-dependent penetrance and variable expressivity are
//! deliberately not mapped.
//!
//! Text is handled as Python does: `strip`, `lower` and `\s` use Python's
//! whitespace, which also counts the separators U+001C to U+001F.

/// HPO terms with a binary penetrance status: complete penetrance, high
/// penetrance, explicit incomplete, moderate and low penetrance, late onset.
pub const PENETRANCE_STATUS_BY_TERM: [(&str, &str); 6] = [
    ("HP:0034950", "complete"),
    ("HP:4000158", "complete"),
    ("HP:0003829", "incomplete"),
    ("HP:4000159", "incomplete"),
    ("HP:4000160", "incomplete"),
    ("HP:0003584", "incomplete"),
];

/// A value the Python functions accept: one text, or a list of values.
#[derive(Clone, Copy, Debug)]
pub enum Values<'a> {
    Text(&'a str),
    List(&'a [&'a str]),
}

pub fn penetrance_status(term: &str) -> Option<&'static str> {
    PENETRANCE_STATUS_BY_TERM
        .iter()
        .find(|(id, _)| *id == term)
        .map(|(_, status)| *status)
}

/// Python's `str.isspace` for one character.
pub fn py_is_space(c: char) -> bool {
    c.is_whitespace() || ('\u{1c}'..='\u{1f}').contains(&c)
}

pub fn py_strip(text: &str) -> &str {
    text.trim_matches(py_is_space)
}

/// The script's `_as_values`: text is split at `;` and `,`; list items are kept
/// whole; both are stripped and empty ones dropped.
pub fn as_values(values: Values) -> Vec<String> {
    let parts: Vec<&str> = match values {
        Values::Text(text) => text.split([';', ',']).collect(),
        Values::List(items) => items.to_vec(),
    };
    parts
        .into_iter()
        .map(py_strip)
        .filter(|part| !part.is_empty())
        .map(str::to_owned)
        .collect()
}

/// Recognised penetrance HPO IDs, upper-cased and sorted.
pub fn recognized_penetrance_hpo_ids(values: Values) -> Vec<String> {
    let mut ids: Vec<String> = as_values(values)
        .into_iter()
        .map(|value| value.to_uppercase())
        .filter(|value| penetrance_status(value).is_some())
        .collect();
    ids.sort();
    ids.dedup();
    ids
}

fn resolve_statuses<S: AsRef<str>>(statuses: &[S]) -> &'static str {
    let normalized: Vec<String> = statuses
        .iter()
        .map(|status| py_strip(status.as_ref()).to_lowercase())
        .filter(|status| !status.is_empty())
        .collect();
    if normalized.iter().any(|status| status == "incomplete") {
        "incomplete"
    } else if normalized
        .iter()
        .any(|status| status == "complete" || status == "high")
    {
        "complete"
    } else {
        ""
    }
}

pub fn normalize_penetrance_hpo_ids(values: Values) -> &'static str {
    let statuses: Vec<&str> = recognized_penetrance_hpo_ids(values)
        .iter()
        .filter_map(|id| penetrance_status(id))
        .collect();
    resolve_statuses(&statuses)
}

/// Controlled source text such as G2P cross-cutting modifiers.
pub fn normalize_penetrance_text(values: Values) -> &'static str {
    let joined = as_values(values)
        .join(" ; ")
        .to_lowercase()
        .replace('_', "-");
    let mut text = String::with_capacity(joined.len());
    let mut in_space = false;
    for c in joined.chars() {
        if py_is_space(c) {
            if !in_space {
                text.push(' ');
            }
            in_space = true;
        } else {
            text.push(c);
            in_space = false;
        }
    }
    let incomplete = [
        "incomplete penetrance",
        "moderate penetrance",
        "low penetrance",
        "late onset",
    ];
    if incomplete.iter().any(|phrase| text.contains(phrase)) {
        "incomplete"
    } else if text.contains("complete penetrance") || text.contains("high penetrance") {
        "complete"
    } else {
        ""
    }
}

/// Source assertions combined with incomplete-over-complete precedence.
pub fn normalize_penetrance_evidence(hpo_ids: Values, raw_values: Values) -> &'static str {
    resolve_statuses(&[
        normalize_penetrance_hpo_ids(hpo_ids),
        normalize_penetrance_text(raw_values),
    ])
}

/// Cache statuses, HPO IDs or controlled raw values; `unknown` is no assertion.
pub fn normalize_penetrance_values(values: Values) -> &'static str {
    let raw = as_values(values);
    let statuses: Vec<String> = raw.iter().map(|value| value.to_lowercase()).collect();
    let direct = resolve_statuses(&statuses);
    let raw: Vec<&str> = raw.iter().map(String::as_str).collect();
    resolve_statuses(&[
        direct,
        normalize_penetrance_hpo_ids(Values::List(&raw)),
        normalize_penetrance_text(Values::List(&raw)),
    ])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn values_split_and_strip_like_python() {
        assert_eq!(
            as_values(Values::Text(" HP:0003829; hp:0034950 ,,\u{1f}")),
            ["HP:0003829", "hp:0034950"]
        );
        assert_eq!(as_values(Values::List(&[" a;b ", "", " "])), ["a;b"]);
    }

    #[test]
    fn hpo_ids_resolve_with_incomplete_first() {
        assert_eq!(
            recognized_penetrance_hpo_ids(Values::Text("hp:4000158;HP:0003584;HP:0003831")),
            ["HP:0003584", "HP:4000158"]
        );
        assert_eq!(
            normalize_penetrance_hpo_ids(Values::Text("HP:0034950;HP:0003829")),
            "incomplete"
        );
        assert_eq!(
            normalize_penetrance_hpo_ids(Values::Text("HP:4000158")),
            "complete"
        );
        assert_eq!(normalize_penetrance_hpo_ids(Values::Text("HP:0003831")), "");
    }

    #[test]
    fn text_and_combined_values() {
        assert_eq!(
            normalize_penetrance_text(Values::Text("Typified by LOW\tpenetrance")),
            "incomplete"
        );
        assert_eq!(
            normalize_penetrance_text(Values::List(&["Low  Penetrance"])),
            "incomplete"
        );
        assert_eq!(normalize_penetrance_text(Values::List(&["late_onset"])), "");
        assert_eq!(
            normalize_penetrance_text(Values::List(&["high penetrance"])),
            "complete"
        );
        assert_eq!(
            normalize_penetrance_evidence(
                Values::Text("HP:4000158"),
                Values::List(&["incomplete penetrance"])
            ),
            "incomplete"
        );
        assert_eq!(
            normalize_penetrance_values(Values::Text("High;unknown")),
            "complete"
        );
        assert_eq!(normalize_penetrance_values(Values::Text("unknown")), "");
        assert_eq!(
            normalize_penetrance_values(Values::Text("hp:0003829")),
            "incomplete"
        );
    }
}
