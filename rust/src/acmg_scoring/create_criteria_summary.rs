//! Corresponds to acmg_scoring.create_criteria_summary, source lines 117–152.
use super::{SourceError, fix_clingen_term::fix_clingen_term, python_strip};

pub fn create_criteria_summary(
    ordered_criteria: &[(&str, i64)],
    clingen_evidence: Option<&str>,
) -> Result<Option<String>, SourceError> {
    let mut active_criteria = Vec::new();
    // Finish PriVA's summary BEFORE examining a ClinGen replacement. Invalid
    // source strengths must still fail even when curated evidence is available.
    for &(name, strength) in ordered_criteria {
        if strength <= 0 {
            continue;
        }
        let suffix = match strength {
            1 => Some("Supporting"),
            2 => Some("Moderate"),
            3 => Some("Strong"),
            4 => Some("Very Strong"),
            5 => Some("Standalone"),
            _ => None,
        };
        let label = match name.chars().nth(1) {
            Some('P') => Some("Supporting"),
            Some('M') => Some("Moderate"),
            Some('S') => Some("Strong"),
            Some('V') => Some("Very Strong"),
            Some('A') => Some("Standalone"),
            Some(_) => None,
            None => {
                return Err(SourceError(
                    "IndexError",
                    "criterion name has no second character".into(),
                ));
            }
        };
        if suffix == label {
            active_criteria.push(name.to_owned());
        } else {
            let suffix = suffix.ok_or_else(|| {
                SourceError(
                    "TypeError",
                    "cannot concatenate a missing strength suffix".into(),
                )
            })?;
            active_criteria.push(format!("{name}_{suffix}"));
        }
    }
    if let Some(evidence) = clingen_evidence.filter(|value| !python_strip(value).is_empty()) {
        active_criteria = evidence
            .split(',')
            .map(python_strip)
            .filter(|c| !c.is_empty())
            .map(fix_clingen_term)
            .collect();
    }
    Ok(if active_criteria.is_empty() {
        None
    } else {
        Some(active_criteria.join(";"))
    })
}
