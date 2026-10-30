//! `scripts/gene_mechanism_variants.py` (642 lines): `infer_query_variant_effect` and the exact
//! curated-allele matcher.
//!
//! Thread T3 (exact-allele matching and the scoring tail). Converted: `infer_query_variant_effect`,
//! which the hub (thread T2) calls for every row. Still to convert: ExactVariantMatcher
//! (match_curated_nonlof_variant, match_gofcards_variant_gof, match_gofcards_variant_gof_by_genomic
//! and loaders).
//!
//! The text helpers of `gene_mechanism_common.py` this module needs (`_clean`, `_bool_value`,
//! `_split_annotation_tokens`, `_bounded_mechanism_score`) are kept here as private functions,
//! because `gene_mechanism_common` belongs to thread T2; they can move there at integration.
use super::types::{AcmgError, EffectInputs, HarnessCall, QueryEffect, not_converted};
use crate::cli::{Result, invalid};
use crate::pandas_table::{Cell, Column, Table};
use crate::python_text::{float_repr, parse_float, python_strip};

/// `priva acmg-step gene_mechanism_variants <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    match function {
        "infer_query_variant_effect" => infer_frame(&call.table),
        _ => Err(not_converted("gene_mechanism_variants", function)),
    }
}

/// `pd.DataFrame([infer_query_variant_effect(row) for _, row in df.iterrows()])`.
fn infer_frame(table: &Table) -> Result<Table> {
    let column = |name: &str| table.column_index(name);
    let positions = [
        column("LoF"),
        column("Consequence"),
        column("NMD"),
        column("LoF_filter"),
        column("splicing_frameshift"),
        column("vep_consq_lof"),
        column("splicing_lof"),
        column("5UTR_lof"),
        column("variant_lof_score"),
        column("variant_gof_score"),
        column("variant_dn_score"),
        column("variant_gof_tag"),
    ];
    let mut effects = Vec::with_capacity(table.rows);
    for row in 0..table.rows {
        let cell = |index: usize| positions[index].map(|position| table.cell(position, row));
        let inputs = EffectInputs {
            lof: cell(0),
            consequence: cell(1),
            nmd: cell(2),
            lof_filter: cell(3),
            splicing_frameshift: cell(4),
            vep_consq_lof: cell(5),
            splicing_lof: cell(6),
            utr5_lof: cell(7),
            variant_lof_score: cell(8),
            variant_gof_score: cell(9),
            variant_dn_score: cell(10),
            variant_gof_tag: cell(11),
        };
        effects.push(infer_query_variant_effect(&inputs)?);
    }
    let texts = |pick: fn(&QueryEffect) -> &String| {
        Column::Text(effects.iter().map(|e| Some(pick(e).clone())).collect())
    };
    let numbers = |pick: fn(&QueryEffect) -> i64| Column::Int(effects.iter().map(pick).collect());
    let names = [
        "variant_effect",
        "variant_effect_evidence",
        "variant_effect_suppressed_evidence",
        "variant_lof_score",
        "variant_gof_score",
        "variant_dn_score",
        "variant_mechanism_exclusive",
        "variant_exact_mechanisms",
    ];
    // An empty list gives pandas a frame without columns.
    if effects.is_empty() {
        return Ok(Table {
            names: Vec::new(),
            columns: Vec::new(),
            rows: 0,
        });
    }
    let columns = vec![
        texts(|e| &e.variant_effect),
        texts(|e| &e.variant_effect_evidence),
        texts(|e| &e.variant_effect_suppressed_evidence),
        numbers(|e| e.variant_lof_score),
        numbers(|e| e.variant_gof_score),
        numbers(|e| e.variant_dn_score),
        Column::Bool(
            effects
                .iter()
                .map(|e| e.variant_mechanism_exclusive)
                .collect(),
        ),
        texts(|e| &e.variant_exact_mechanisms),
    ];
    Ok(Table {
        names: names.iter().map(|name| (*name).to_owned()).collect(),
        columns,
        rows: effects.len(),
    })
}

/// Python's `str(value)` of a table value (`nan` for a missing one).
fn python_str(cell: Cell<'_>) -> String {
    match cell {
        Cell::Missing => "nan".to_owned(),
        Cell::Text(text) => text.to_owned(),
        Cell::Int(number) => number.to_string(),
        Cell::Float(number) => float_repr(number),
        Cell::Bool(flag) => (if flag { "True" } else { "False" }).to_owned(),
    }
}

/// `gene_mechanism_common._clean`: the stripped text, or `""` for a missing
/// value and the spellings `nan`, `none` and `<na>` in any case.
pub(crate) fn clean(cell: Option<Cell<'_>>) -> String {
    let Some(cell) = cell else {
        return String::new();
    };
    clean_text(&python_str(cell))
}

/// `_clean` of a text.
pub(crate) fn clean_text(text: &str) -> String {
    let text = python_strip(text);
    // Only ASCII letters lower-case to these spellings, so an ASCII comparison is Python's.
    if ["", "nan", "none", "<na>"]
        .iter()
        .any(|word| text.eq_ignore_ascii_case(word))
    {
        String::new()
    } else {
        text.to_owned()
    }
}

/// `gene_mechanism_common._bool_value`: True, a non-zero integer, or the texts
/// `1`, `true`, `t`, `yes`, `y` in any case. A float (other than NaN, which is
/// false) goes through its text, which is never one of those spellings.
fn bool_value(cell: Option<Cell<'_>>) -> bool {
    match cell {
        None | Some(Cell::Missing) => false,
        Some(Cell::Bool(flag)) => flag,
        Some(Cell::Int(number)) => number != 0,
        Some(other) => {
            let text = clean(Some(other));
            ["1", "true", "t", "yes", "y"]
                .iter()
                .any(|word| text.eq_ignore_ascii_case(word))
        }
    }
}

/// `gene_mechanism_common._split_annotation_tokens`: the upper-cased tokens
/// between `&`, `,`, `;` and `|`, as a set (here a list without repeats).
fn split_annotation_tokens(cell: Option<Cell<'_>>) -> Vec<String> {
    let mut tokens: Vec<String> = Vec::new();
    for token in clean(cell).split(['&', ',', ';', '|']) {
        let token = python_strip(token);
        if token.is_empty() {
            continue;
        }
        let token = token.to_uppercase();
        if !tokens.contains(&token) {
            tokens.push(token);
        }
    }
    tokens
}

/// `gene_mechanism_common._bounded_mechanism_score`: `int(float(_clean(value)))`
/// when it is 0, 1 or 2, else 0; text Python cannot read as a float gives 0.
/// An infinite value raises OverflowError in Python, and here.
fn bounded_mechanism_score(cell: Option<Cell<'_>>) -> Result<i64> {
    let overflow = || -> Box<dyn std::error::Error> {
        Box::new(AcmgError::new(
            "OverflowError",
            "cannot convert float infinity to integer",
        ))
    };
    let number = match cell {
        None | Some(Cell::Missing) | Some(Cell::Bool(_)) => return Ok(0),
        Some(Cell::Int(number)) => return Ok(if (0..=2).contains(&number) { number } else { 0 }),
        Some(Cell::Float(number)) => number,
        Some(Cell::Text(text)) => {
            let text = clean_text(text);
            match parse_float(&text) {
                Some(number) => number,
                // Python's float() also reads digits of other scripts: not reproduced.
                None if !text.is_ascii() => {
                    return Err(invalid(format!(
                        "refused (not reproduced): float() of the non-ASCII text {text:?}"
                    )));
                }
                None => return Ok(0),
            }
        }
    };
    if number.is_nan() {
        return Ok(0);
    }
    if number.is_infinite() {
        return Err(overflow());
    }
    let truncated = number.trunc();
    Ok(
        if truncated == 0.0 || truncated == 1.0 || truncated == 2.0 {
            truncated as i64
        } else {
            0
        },
    )
}

/// The curated mechanisms a score of 2 can establish (`CURATED_EXACT_MECHANISMS`), sorted.
const CURATED_EXACT_MECHANISMS: [&str; 2] = ["DOMINANT_NEGATIVE", "GOF"];

/// `infer_query_variant_effect(row)` for one table row. The hub (T2) calls it
/// for every row.
pub(crate) fn infer_query_variant_effect(row: &EffectInputs<'_>) -> Result<QueryEffect> {
    let mut predicted: Vec<&'static str> = Vec::new();
    let loftee = split_annotation_tokens(row.lof);
    if loftee.iter().any(|token| token == "HC") {
        predicted.push("LOFTEE_HC");
    }
    if loftee.iter().any(|token| token == "OS") {
        predicted.push("LOFTEE_OS");
    }
    let consequence = clean(row.consequence).to_lowercase();
    let nmd = clean(row.nmd).to_lowercase();
    let lof_filter = clean(row.lof_filter).to_uppercase();
    let truncating = consequence.contains("stop_gained") || consequence.contains("frameshift");
    let escapes_nmd = nmd.contains("escaping") || lof_filter.contains("END_TRUNC");
    let splice_frameshift = bool_value(row.splicing_frameshift);
    if (truncating || splice_frameshift) && !escapes_nmd {
        predicted.push("NMD_PREDICTED_LOF");
    } else if truncating || splice_frameshift {
        predicted.push("NMD_ESCAPING_TRUNCATION");
    }
    if bool_value(row.vep_consq_lof) {
        predicted.push("VEP_LOF");
    }
    if bool_value(row.splicing_lof) {
        predicted.push("PRIVA_SPLICE_LOF");
    }
    if bool_value(row.utr5_lof) {
        predicted.push("PRIVA_5UTR_LOF");
    }

    // Scores in the order LOF, GOF, DOMINANT_NEGATIVE.
    let mut scores = [
        bounded_mechanism_score(row.variant_lof_score)?,
        bounded_mechanism_score(row.variant_gof_score)?,
        bounded_mechanism_score(row.variant_dn_score)?,
    ];
    if split_annotation_tokens(row.variant_gof_tag)
        .iter()
        .any(|token| token == "GOF")
    {
        scores[1] = 2;
    }
    // Sorted, as `sorted(curated_mechanisms)` gives them.
    let curated: Vec<&str> = CURATED_EXACT_MECHANISMS
        .into_iter()
        .filter(|mechanism| scores[if *mechanism == "GOF" { 1 } else { 2 }] == 2)
        .collect();
    let canonical = |mechanisms: &[&str]| -> Vec<String> {
        mechanisms
            .iter()
            .map(|mechanism| format!("CANONICAL_EXACT_{mechanism}"))
            .collect()
    };
    let predicted: Vec<String> = predicted.into_iter().map(str::to_owned).collect();
    let (effect, evidence, suppressed, exact): (String, Vec<String>, Vec<String>, Vec<&str>) =
        if predicted.iter().any(|token| token == "NMD_PREDICTED_LOF") {
            scores = [2, 0, 0];
            (
                "exact_known_LOF".to_owned(),
                predicted,
                canonical(&curated),
                vec!["LOF"],
            )
        } else if !curated.is_empty() {
            scores = [
                0,
                if curated.contains(&"GOF") { 2 } else { 0 },
                if curated.contains(&"DOMINANT_NEGATIVE") {
                    2
                } else {
                    0
                },
            ];
            (
                format!("exact_known_{}", curated.join("+")),
                canonical(&curated),
                predicted,
                curated.clone(),
            )
        } else if !predicted.is_empty() {
            let floor = if predicted.iter().any(|token| token == "LOFTEE_HC") {
                2
            } else {
                1
            };
            scores[0] = scores[0].max(floor);
            (
                "predicted_LOF_high_confidence".to_owned(),
                predicted,
                Vec::new(),
                Vec::new(),
            )
        } else {
            scores = scores.map(|score| score.max(1));
            ("uncertain".to_owned(), Vec::new(), Vec::new(), Vec::new())
        };
    Ok(QueryEffect {
        variant_effect: effect,
        variant_effect_evidence: evidence.join(";"),
        variant_effect_suppressed_evidence: suppressed.join(";"),
        variant_lof_score: scores[0],
        variant_gof_score: scores[1],
        variant_dn_score: scores[2],
        variant_mechanism_exclusive: !exact.is_empty(),
        // `exact` is already sorted (LOF alone, or DOMINANT_NEGATIVE before GOF).
        variant_exact_mechanisms: exact.join(";"),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn effect(inputs: EffectInputs<'_>) -> QueryEffect {
        infer_query_variant_effect(&inputs).unwrap()
    }

    #[test]
    fn decay_lof_outranks_a_curated_allele() {
        let call = effect(EffectInputs {
            consequence: Some(Cell::Text("Stop_Gained&splice_region_variant")),
            variant_gof_score: Some(Cell::Int(2)),
            variant_dn_score: Some(Cell::Text(" 2.7 ")),
            lof: Some(Cell::Text("HC")),
            ..EffectInputs::default()
        });
        assert_eq!(call.variant_effect, "exact_known_LOF");
        assert_eq!(call.variant_effect_evidence, "LOFTEE_HC;NMD_PREDICTED_LOF");
        assert_eq!(
            call.variant_effect_suppressed_evidence,
            "CANONICAL_EXACT_DOMINANT_NEGATIVE;CANONICAL_EXACT_GOF"
        );
        assert_eq!(
            (
                call.variant_lof_score,
                call.variant_gof_score,
                call.variant_dn_score
            ),
            (2, 0, 0)
        );
        assert!(call.variant_mechanism_exclusive);
        assert_eq!(call.variant_exact_mechanisms, "LOF");
    }

    #[test]
    fn curated_alleles_and_the_gof_tag() {
        let call = effect(EffectInputs {
            variant_gof_tag: Some(Cell::Text("x;gof")),
            vep_consq_lof: Some(Cell::Bool(true)),
            consequence: Some(Cell::Text("frameshift_variant")),
            nmd: Some(Cell::Text("NMD_escaping_variant")),
            ..EffectInputs::default()
        });
        assert_eq!(call.variant_effect, "exact_known_GOF");
        assert_eq!(call.variant_effect_evidence, "CANONICAL_EXACT_GOF");
        assert_eq!(
            call.variant_effect_suppressed_evidence,
            "NMD_ESCAPING_TRUNCATION;VEP_LOF"
        );
        assert_eq!(
            (
                call.variant_lof_score,
                call.variant_gof_score,
                call.variant_dn_score
            ),
            (0, 2, 0)
        );
        assert_eq!(call.variant_exact_mechanisms, "GOF");
    }

    #[test]
    fn predicted_and_uncertain_calls() {
        let call = effect(EffectInputs {
            splicing_lof: Some(Cell::Text(" Yes ")),
            variant_lof_score: Some(Cell::Float(1.0)),
            ..EffectInputs::default()
        });
        assert_eq!(call.variant_effect, "predicted_LOF_high_confidence");
        assert_eq!(call.variant_lof_score, 1);
        assert!(!call.variant_mechanism_exclusive);
        // A float 1.0 is not a true flag in _bool_value; NaN gives 0.
        let call = effect(EffectInputs {
            utr5_lof: Some(Cell::Float(1.0)),
            variant_gof_score: Some(Cell::Missing),
            variant_dn_score: Some(Cell::Text("nan")),
            ..EffectInputs::default()
        });
        assert_eq!(call.variant_effect, "uncertain");
        assert_eq!(
            (
                call.variant_lof_score,
                call.variant_gof_score,
                call.variant_dn_score
            ),
            (1, 1, 1)
        );
        assert_eq!(call.variant_effect_evidence, "");
    }

    #[test]
    fn scores_follow_int_of_float() {
        assert_eq!(bounded_mechanism_score(Some(Cell::Text("1_0"))).unwrap(), 0);
        assert_eq!(bounded_mechanism_score(Some(Cell::Text("2.9"))).unwrap(), 2);
        assert_eq!(
            bounded_mechanism_score(Some(Cell::Text("-0.5"))).unwrap(),
            0
        );
        assert_eq!(bounded_mechanism_score(Some(Cell::Text("abc"))).unwrap(), 0);
        assert_eq!(
            bounded_mechanism_score(Some(Cell::Text("-nan"))).unwrap(),
            0
        );
        assert_eq!(bounded_mechanism_score(Some(Cell::Int(3))).unwrap(), 0);
        assert_eq!(bounded_mechanism_score(Some(Cell::Bool(true))).unwrap(), 0);
        assert_eq!(
            bounded_mechanism_score(Some(Cell::Text("Infinity")))
                .unwrap_err()
                .to_string(),
            "OverflowError: cannot convert float infinity to integer"
        );
        assert!(bounded_mechanism_score(Some(Cell::Text("\u{ff12}"))).is_err());
    }
}
