//! Complete frequency unit, opt-in and unregistered until separately qualified.
//! Original nine functions: ef2e70115ddfe7087b4d24002ce4af95c4f4b91fbf38a06e8a52275eecca4fb3.
//! This module evaluates PriVA table/rule operations. SciPy sf/cdf and the
//! original gzip/pickle codecs remain mature host dependencies; no tail
//! probability implementation, default switch or original-source edit lives here.
pub mod values;
pub mod columns;
pub mod binomial;
pub mod rules;
pub mod transport;
pub use values::{Cell, Value, Result};
pub use columns::calculate;

pub const PROCEDURES: [&str; 10] = [
    "_gene_max_patho_hom_freq", "control_false_neg_rate",
    "control_false_positive_rate", "gnomAD_rare_AF", "_rare_af_binomial",
    "_rare_hom_binomial", "PM2_criteria", "BS1_criteria", "BA1_criteria",
    "ACMG_criteria_assign",
];

/// Source-site tokens are checked independently of the untrusted wire operands.
/// Exact expression-to-operation correspondence is additionally source-pinned
/// by the host lowering manifest; it is not inferred from output agreement.
pub fn valid_site(site: &str, operation: &str) -> bool {
    if operation == "quit" { return site == "close@0"; }
    let Some((name, line)) = site.rsplit_once('@') else { return false; };
    let Ok(line) = line.parse::<usize>() else { return false; };
    let range = match name {
        "_gene_max_patho_hom_freq" => (43,65),
        "control_false_neg_rate" => (68,163),
        "control_false_positive_rate" => (166,253),
        "gnomAD_rare_AF" => (256,269),
        "_rare_af_binomial" => (272,310),
        "_rare_hom_binomial" => (313,343),
        "PM2_criteria" => (346,413),
        "BS1_criteria" | "_series_or_empty" | "_numeric_series" |
        "_frequency_above" | "_hom_bs1" => (416,776),
        "BA1_criteria" => (779,799),
        "_variant_condition_masks" | "has_token" => (267,340),
        "ACMG_criteria_assign" => (391,443),
        _ => return false,
    };
    (range.0..=range.1).contains(&line)
}
