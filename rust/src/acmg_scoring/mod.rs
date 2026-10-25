//! Frozen source: scripts/acmg_scoring.py at 09223ef5a8b0ceb4c6d61055c6eefe62b51ba903.
//! SHA256: 9111a76e3fc94a064e14d97b535bb1dd2132f7c61836f8e8892477d4c5e69646.
//!
//! Scalar functions retain their source boundaries. This does not assign evidence,
//! load pickle resources, calculate ranking penalties or qualify a family run.
pub mod calculate_posterior_probability;
pub mod create_criteria_summary;
pub mod fix_clingen_term;
pub mod ranking_setup;
pub(crate) mod ranking_setup_transport;
pub(crate) mod score_table;
pub mod terminal_order;
pub(crate) mod terminal_order_transport;
pub mod translate_strength_level;

pub const CRITERIA_ORDER: [&str; 28] = [
    "PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4", "PM5", "PM6", "PP1", "PP2",
    "PP3", "PP4", "PP5", "BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP5",
    "BP6", "BP7",
];

/// An original Python failure category, retained rather than silently repaired.
#[derive(Debug, PartialEq)]
pub struct SourceError(pub &'static str, pub String);

impl std::fmt::Display for SourceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.0, self.1)
    }
}
impl std::error::Error for SourceError {}

// Python str.strip also recognizes the four information separators. Rust's
// Unicode White_Space predicate alone omits these characters.
pub(crate) fn python_strip(value: &str) -> &str {
    value.trim_matches(|c: char| c.is_whitespace() || ('\u{001c}'..='\u{001f}').contains(&c))
}
