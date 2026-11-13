//! Corresponds to acmg_scoring.fix_clingen_term, source lines 60–85.
pub fn fix_clingen_term(clingen_str: &str) -> String {
    if !clingen_str.contains('_') {
        return clingen_str.to_owned();
    }
    let mut parts = clingen_str.split('_');
    let prefix = parts.next().unwrap();
    // The reference deliberately takes ONLY split()[1], discarding later parts.
    let suffix = parts.next().unwrap();
    let suffix = if suffix.starts_with("Very") {
        "Very Strong"
    } else if suffix.starts_with("Stand") {
        "Stand Alone"
    } else {
        suffix
    };
    format!("{prefix}_{suffix}")
}
