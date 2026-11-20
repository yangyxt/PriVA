//! Corresponds to acmg_scoring.translate_strength_level, source lines 88–114.
use super::SourceError;

pub fn translate_strength_level(
    criteria_name: &str,
    criteria_list_str: Option<&str>,
) -> Result<i64, SourceError> {
    let Some(summary) = criteria_list_str else {
        return Ok(0);
    };
    if !summary.contains(criteria_name) {
        return Ok(0);
    }
    let criteria: Vec<&str> = summary.split(';').collect();
    if criteria.contains(&criteria_name) {
        return match criteria_name.chars().nth(1) {
            Some('V') => Ok(4),
            Some('S') => Ok(3),
            Some('M') => Ok(2),
            Some('P') => Ok(1),
            Some('A') => Ok(5),
            Some(value) => Err(SourceError("KeyError", value.to_string())),
            None => Err(SourceError(
                "IndexError",
                "criterion name has no second character".into(),
            )),
        };
    }
    // Keep substring check → first startswith match → split()[1] ordering.
    let matching = criteria
        .iter()
        .find(|c| c.starts_with(criteria_name))
        .ok_or_else(|| SourceError("IndexError", "no prefix-matching criterion".into()))?;
    let suffix = matching
        .split('_')
        .nth(1)
        .ok_or_else(|| SourceError("IndexError", "criterion has no strength suffix".into()))?;
    match suffix {
        "Supporting" => Ok(1),
        "Moderate" => Ok(2),
        "Strong" => Ok(3),
        "Very Strong" | "Very" => Ok(4),
        "Stand Alone" | "Stand" | "Standalone" => Ok(5),
        value => Err(SourceError("KeyError", value.into())),
    }
}
