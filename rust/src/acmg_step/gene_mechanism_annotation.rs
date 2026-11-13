//! `scripts/gene_mechanism_annotation.py` (411 lines): the per-row chain that writes the hub's 16
//! columns.
//!
//! Thread T2 (gene-condition hub). Converted: `annotate_gene_mechanism_categories` with
//! `select_condition_histories_for_variant`, `_classify_variant_applicability`,
//! `_mechanism_profile_tag` and `_compact_profile_tags`. `_compact_inheritance` is not called
//! by anything. The variant's own effect comes from
//! `gene_mechanism_variants::infer_query_variant_effect` (thread T3).
//!
//! Python walks the table with `DataFrame.iterrows()`, which hands each row over in one common
//! type: a table without `object` columns whose columns are all numbers gives every integer as a
//! float. That is reproduced; a table whose rows pandas would hand over in a nullable type (no
//! `object` column, an `Int64` or `boolean` one) is refused.
use super::gene_mechanism_common::{
    UNRESOLVED_MECHANISM, VARIANT_MECHANISM_OUTPUT_COLUMNS, VARIANT_MECHANISM_SCORE_KEYS, clean,
    clean_text, deduplicate_by, exact_mechanisms_from_effect, python_error, python_text_repr,
    refused,
};
use super::gene_mechanism_conditions::{
    ConditionAssertion, NON_MENDELIAN_INHERITANCE, condition_cache_mechanism_assertions,
    gene_has_lof_mechanism_history, gene_inheritance_from_constraint,
    normalize_inheritance_histories, normalize_penetrance,
};
use super::gene_mechanism_resources::{GeneMechanismHub, HubPaths, harness_paths};
use super::gene_mechanism_variants::infer_query_variant_effect;
use super::types::{EffectInputs, HarnessCall, SymbolResolver, not_converted};
use crate::cli::{Result, invalid};
use crate::json_text::Json;
use crate::pandas_table::{Cell, Column, Table};
use std::collections::HashMap;
use std::rc::Rc;

/// The columns `infer_query_variant_effect` reads, in the order of `EffectInputs`.
const EFFECT_COLUMNS: [&str; 12] = [
    "LoF",
    "Consequence",
    "NMD",
    "LoF_filter",
    "splicing_frameshift",
    "vep_consq_lof",
    "splicing_lof",
    "5UTR_lof",
    "variant_lof_score",
    "variant_gof_score",
    "variant_dn_score",
    "variant_gof_tag",
];

/// One condition history a variant's mechanism reaches.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct History {
    pub(crate) mechanism: String,
    pub(crate) condition: String,
    pub(crate) inheritance: String,
    pub(crate) x_linked: bool,
    pub(crate) penetrance: &'static str,
}

/// `select_condition_histories_for_variant(assertions, variant_effect=...)`:
/// the histories whose mechanism the variant's effect allows (exact mechanisms,
/// LOF for a predicted loss, all three for an uncertain effect), and every
/// `UNRESOLVED` history; one per normalised inheritance.
pub(crate) fn select_condition_histories_for_variant(
    assertions: &[ConditionAssertion],
    variant_effect: &str,
) -> Vec<History> {
    let effect = clean_text(variant_effect);
    let exact = exact_mechanisms_from_effect(&effect);
    let allowed: Vec<String> = if !exact.is_empty() {
        exact
    } else if effect == "predicted_LOF_high_confidence" {
        vec!["LOF".to_owned()]
    } else if effect == "uncertain" {
        VARIANT_MECHANISM_SCORE_KEYS
            .iter()
            .map(|key| (*key).to_owned())
            .collect()
    } else {
        Vec::new()
    };
    let mut histories = Vec::new();
    for assertion in assertions {
        let mechanism = clean_text(&assertion.mechanism).to_uppercase();
        if mechanism != UNRESOLVED_MECHANISM && !allowed.contains(&mechanism) {
            continue;
        }
        let context = &assertion.context;
        for (inheritance, x_linked) in normalize_inheritance_histories(
            &assertion.allelic_requirement,
            &context.hpo_inheritance_modes,
        ) {
            let condition = if context.hpo_disease_id.is_empty() {
                &assertion.source_condition_id
            } else {
                &context.hpo_disease_id
            };
            let penetrance_values = if context.penetrance_statuses.is_empty() {
                &context.penetrance_hpo_ids
            } else {
                &context.penetrance_statuses
            };
            histories.push(History {
                mechanism: mechanism.clone(),
                condition: clean_text(condition),
                inheritance,
                x_linked,
                penetrance: normalize_penetrance(penetrance_values),
            });
        }
    }
    histories
}

/// `_mechanism_profile_tag(history)`: `<inheritance>_<MECHANISM>` (DN for
/// dominant negative), the inheritance alone for an unresolved mechanism, and
/// nothing for a non-Mendelian inheritance.
fn mechanism_profile_tag(history: &History) -> String {
    let inheritance = clean_text(&history.inheritance);
    if NON_MENDELIAN_INHERITANCE.contains(&inheritance.as_str()) {
        return String::new();
    }
    let mut mechanism = clean_text(&history.mechanism).to_uppercase();
    if mechanism.is_empty() {
        mechanism = "UNRESOLVED".to_owned();
    }
    if mechanism == "UNRESOLVED" {
        return if inheritance.is_empty() {
            "uncertain".to_owned()
        } else {
            inheritance
        };
    }
    if mechanism == "DOMINANT_NEGATIVE" {
        mechanism = "DN".to_owned();
    }
    if inheritance.is_empty() {
        mechanism
    } else {
        format!("{inheritance}_{mechanism}")
    }
}

/// `_compact_profile_tags(histories)`: the distinct tags, without a bare
/// inheritance or a bare mechanism that a qualified tag already names, and
/// without `uncertain`.
fn compact_profile_tags(histories: &[History]) -> Vec<String> {
    let mut tags: Vec<String> = Vec::new();
    for history in histories {
        let tag = mechanism_profile_tag(history);
        if !tag.is_empty() && !tags.contains(&tag) {
            tags.push(tag);
        }
    }
    const INHERITANCE_VALUES: [&str; 7] = [
        "x_linked_recessive",
        "x_linked_dominant",
        "x_linked_unspecified",
        "mitochondrial",
        "recessive",
        "dominant",
        "y_linked",
    ];
    let split = |tag: &str| -> (String, String) {
        for inheritance in INHERITANCE_VALUES {
            if let Some(mechanism) = tag.strip_prefix(&format!("{inheritance}_")) {
                return (inheritance.to_owned(), mechanism.to_owned());
            }
        }
        (String::new(), String::new())
    };
    let mut inheritance_with_mechanism: Vec<String> = Vec::new();
    let mut qualified_mechanisms: Vec<String> = Vec::new();
    for tag in &tags {
        let (inheritance, mechanism) = split(tag);
        if !mechanism.is_empty() {
            inheritance_with_mechanism.push(inheritance);
            qualified_mechanisms.push(mechanism);
        }
    }
    tags.into_iter()
        .filter(|tag| {
            !inheritance_with_mechanism.contains(tag)
                && !qualified_mechanisms.contains(tag)
                && tag != "uncertain"
        })
        .collect()
}

/// What `_classify_variant_applicability` returns, each a `;`-joined text.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Applicability {
    pub(crate) plausible: String,
    pub(crate) applicable: String,
    pub(crate) uncertain: String,
}

/// `_classify_variant_applicability(histories, scores)`: a history whose
/// mechanism the variant scores 2 is applicable, any other uncertain.
/// `scores` are the variant's LOF, GOF and DOMINANT_NEGATIVE scores.
fn classify_variant_applicability(histories: &[History], scores: [i64; 3]) -> Applicability {
    let mut applicable: Vec<String> = Vec::new();
    let mut uncertain: Vec<String> = Vec::new();
    for history in histories {
        let mechanism = clean_text(&history.mechanism).to_uppercase();
        let score = VARIANT_MECHANISM_SCORE_KEYS
            .iter()
            .position(|key| *key == mechanism)
            .map_or(0, |index| scores[index]);
        let established = mechanism != UNRESOLVED_MECHANISM && score == 2;
        let tag = mechanism_profile_tag(history);
        if tag.is_empty() {
            continue;
        }
        let target = if established {
            &mut applicable
        } else {
            &mut uncertain
        };
        if !target.contains(&tag) {
            target.push(tag);
        }
    }
    Applicability {
        plausible: compact_profile_tags(histories).join(";"),
        applicable: applicable.join(";"),
        uncertain: uncertain.join(";"),
    }
}

/// `dict.fromkeys(values)` joined by `;`.
fn join_distinct<'a>(values: impl Iterator<Item = &'a str>) -> String {
    let mut seen: Vec<&str> = Vec::new();
    for value in values {
        if !seen.contains(&value) {
            seen.push(value);
        }
    }
    seen.join(";")
}

/// How `DataFrame.iterrows()` hands a row's values over.
#[derive(Clone, Copy, Debug, PartialEq)]
enum RowValues {
    /// An `object` row (any `object` column, or booleans mixed with numbers),
    /// or a row of one numpy type: the values as stored.
    AsStored,
    /// Integer and float columns only: numpy's common type is float64.
    IntegersAsFloats,
}

fn row_values(table: &Table) -> Result<RowValues> {
    let all = |test: fn(&Column) -> bool| table.columns.iter().all(test);
    let any = |test: fn(&Column) -> bool| table.columns.iter().any(test);
    if any(|column| {
        matches!(
            column,
            Column::Text(_) | Column::Object(_) | Column::BoolOrMissing(_)
        )
    }) || all(|column| matches!(column, Column::Int(_)))
        || all(|column| matches!(column, Column::Float(_)))
        || all(|column| matches!(column, Column::Bool(_)))
    {
        return Ok(RowValues::AsStored);
    }
    if all(|column| matches!(column, Column::Int(_) | Column::Float(_))) {
        return Ok(RowValues::IntegersAsFloats);
    }
    if any(|column| matches!(column, Column::Bool(_)))
        && all(|column| matches!(column, Column::Int(_) | Column::Float(_) | Column::Bool(_)))
    {
        return Ok(RowValues::AsStored);
    }
    Err(refused(
        "a table without object columns whose rows pandas hands over in a nullable type",
    ))
}

/// `row.get(name)` of the row `iterrows()` gives (`None` for an absent column).
fn row_cell<'a>(
    table: &'a Table,
    kind: RowValues,
    column: Option<usize>,
    row: usize,
) -> Option<Cell<'a>> {
    let cell = table.cell(column?, row);
    Some(match (kind, cell) {
        (RowValues::IntegersAsFloats, Cell::Int(number)) => Cell::Float(number as f64),
        _ => cell,
    })
}

/// The column of a name the chain reads; a repeated name (pandas would hand a
/// Series over) is refused.
fn single_column(table: &Table, name: &str) -> Result<Option<usize>> {
    if table
        .names
        .iter()
        .filter(|candidate| *candidate == name)
        .count()
        > 1
    {
        return Err(refused(format!("the column name {name:?} repeats")));
    }
    Ok(table.column_index(name))
}

/// `annotate_gene_mechanism_categories(df, condition_cache=..., symbol_col=...,
/// output_col=..., use_hgnc_package=False, clingen_dosage=..., loeuf_table=...,
/// hgnc_table=...)`: the table with `output_col` and the 15
/// `VARIANT_MECHANISM_OUTPUT_COLUMNS` set (new columns at the end, existing
/// ones replaced in place).
pub(crate) fn annotate_gene_mechanism_categories(
    table: &Table,
    paths: HubPaths,
    symbol_col: &str,
    output_col: &str,
) -> Result<Table> {
    if !table.has_column(symbol_col) {
        let message = format!("missing symbol column: {symbol_col}");
        return Err(python_error("KeyError", python_text_repr(&message)?));
    }
    let mut hub = GeneMechanismHub::new(paths)?;
    let condition_records = hub.load_condition_cache()?;
    let rows = table.rows;
    // Only rows are read, so an empty table needs none of these checks.
    let (kind, symbol_position, chrom_position, effect_positions) = if rows == 0 {
        (RowValues::AsStored, None, None, Vec::new())
    } else {
        (
            row_values(table)?,
            single_column(table, symbol_col)?,
            single_column(table, "chrom")?,
            EFFECT_COLUMNS
                .iter()
                .map(|name| single_column(table, name))
                .collect::<Result<Vec<Option<usize>>>>()?,
        )
    };
    let empty_record = Json::Object(Vec::new());
    let mut assertion_cache: HashMap<String, Rc<Vec<ConditionAssertion>>> = HashMap::new();
    let mut lof_history_cache: HashMap<String, bool> = HashMap::new();

    let mut plausible = Vec::with_capacity(rows);
    let mut texts: HashMap<&str, Vec<Option<String>>> = HashMap::new();
    let mut scores: [Vec<i64>; 3] = [Vec::new(), Vec::new(), Vec::new()];
    let mut exclusive = Vec::with_capacity(rows);
    for row in 0..rows {
        let gene = row_cell(table, kind, symbol_position, row).expect("the symbol column exists");
        let symbol = hub.resolve_symbol(gene)?;
        let gene_record = condition_records.gene(&symbol).unwrap_or(&empty_record);
        if !assertion_cache.contains_key(&symbol) {
            let assertions = condition_cache_mechanism_assertions(gene_record)?;
            assertion_cache.insert(symbol.clone(), Rc::new(assertions));
        }
        let assertions = Rc::clone(&assertion_cache[&symbol]);
        if !lof_history_cache.contains_key(&symbol) {
            let history = gene_has_lof_mechanism_history(gene_record)?;
            lof_history_cache.insert(symbol.clone(), history);
        }

        let cell = |index: usize| row_cell(table, kind, effect_positions[index], row);
        let effect = infer_query_variant_effect(&EffectInputs {
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
        })?;
        let histories = deduplicate_by(
            select_condition_histories_for_variant(&assertions, &effect.variant_effect),
            |history| {
                (
                    clean_text(&history.mechanism),
                    clean_text(&history.condition),
                    clean_text(&history.inheritance),
                    history.x_linked,
                    clean_text(history.penetrance),
                )
            },
        );
        let effect_scores = [
            effect.variant_lof_score,
            effect.variant_gof_score,
            effect.variant_dn_score,
        ];
        let applicability = classify_variant_applicability(&histories, effect_scores);

        let chromosome = row_cell(table, kind, chrom_position, row);
        let lowered = clean(chromosome).to_lowercase();
        let row_chromosome = lowered.strip_prefix("chr").unwrap_or(&lowered);
        let mut matched_inheritance: Vec<String> = Vec::new();
        for history in &histories {
            if !history.inheritance.is_empty()
                && !matched_inheritance.contains(&history.inheritance)
            {
                matched_inheritance.push(history.inheritance.clone());
            }
        }
        let mut matched_x_linked =
            row_chromosome == "x" || histories.iter().any(|history| history.x_linked);
        let basis;
        if !matched_inheritance.is_empty() {
            let mechanisms: Vec<String> = histories
                .iter()
                .map(|history| clean_text(&history.mechanism).to_uppercase())
                .collect();
            let known = mechanisms
                .iter()
                .any(|mechanism| mechanism != UNRESOLVED_MECHANISM);
            let unresolved = mechanisms
                .iter()
                .any(|mechanism| mechanism == UNRESOLVED_MECHANISM);
            basis = if known && unresolved {
                "matched_and_unresolved_history"
            } else if unresolved {
                "unresolved_condition_history"
            } else {
                "matched_history"
            };
        } else if !histories.is_empty() || effect.variant_lof_score >= 1 {
            // A condition matched without an inheritance, or none matched for a
            // plausibly loss-of-function variant: the gene's constraint fills in.
            let clingen = hub.load_clingen()?;
            let loeuf = hub.load_loeuf()?;
            let inferred = gene_inheritance_from_constraint(
                clingen
                    .get(&symbol)
                    .and_then(|record| record.haploinsufficiency_score),
                loeuf.get(&symbol).copied(),
                chromosome,
            );
            if !inferred.is_empty() {
                matched_inheritance = vec![inferred.to_owned()];
            }
            matched_x_linked = inferred.starts_with("x_");
            basis = "gene_constraint";
        } else {
            matched_x_linked = row_chromosome == "x";
            basis = "mechanism_mismatch";
        }

        let mut push = |name: &'static str, value: String| {
            texts.entry(name).or_default().push(Some(value));
        };
        plausible.push(Some(applicability.plausible));
        push("variant_effect", effect.variant_effect.clone());
        for (index, score) in effect_scores.iter().enumerate() {
            scores[index].push(*score);
        }
        exclusive.push(effect.variant_mechanism_exclusive);
        push(
            "variant_exact_mechanisms",
            effect.variant_exact_mechanisms.clone(),
        );
        push("variant_mechanism_applicable", applicability.applicable);
        push("variant_mechanism_uncertain", applicability.uncertain);
        push(
            "variant_condition_ids",
            join_distinct(
                histories
                    .iter()
                    .map(|history| history.condition.as_str())
                    .filter(|condition| !condition.is_empty()),
            ),
        );
        let entries: Vec<String> = histories
            .iter()
            .filter(|history| !history.condition.is_empty())
            .map(|history| {
                let mechanism = if history.mechanism == "DOMINANT_NEGATIVE" {
                    "DN"
                } else {
                    history.mechanism.as_str()
                };
                let inheritance = if history.inheritance.is_empty() {
                    "unknown"
                } else {
                    history.inheritance.as_str()
                };
                [
                    history.condition.as_str(),
                    mechanism,
                    inheritance,
                    history.penetrance,
                ]
                .join("|")
            })
            .collect();
        push(
            "variant_condition_histories",
            join_distinct(entries.iter().map(String::as_str)),
        );
        push("variant_inheritance", matched_inheritance.join(";"));
        push("variant_inheritance_basis", basis.to_owned());
        push(
            "variant_x_linked",
            (if matched_x_linked { "true" } else { "false" }).to_owned(),
        );
        push(
            "gene_lof_mechanism_history",
            (if lof_history_cache[&symbol] {
                "true"
            } else {
                "false"
            })
            .to_owned(),
        );
        push(
            "variant_penetrance",
            join_distinct(
                histories
                    .iter()
                    .map(|history| history.penetrance)
                    .filter(|penetrance| !penetrance.is_empty()),
            ),
        );
    }

    let mut out = table.clone();
    // `df[name] = []` on an empty table makes a float64 column.
    let column_of = |values: Column| {
        if rows == 0 {
            Column::Float(Vec::new())
        } else {
            values
        }
    };
    out.set_column(output_col, column_of(Column::Text(plausible)));
    let [lof, gof, dn] = scores;
    let mut scores = [Some(lof), Some(gof), Some(dn)];
    let mut exclusive = Some(exclusive);
    for name in VARIANT_MECHANISM_OUTPUT_COLUMNS {
        let values = match name {
            "variant_lof_score" => Column::Int(scores[0].take().unwrap_or_default()),
            "variant_gof_score" => Column::Int(scores[1].take().unwrap_or_default()),
            "variant_dn_score" => Column::Int(scores[2].take().unwrap_or_default()),
            "variant_mechanism_exclusive" => Column::Bool(exclusive.take().unwrap_or_default()),
            _ => Column::Text(texts.remove(name).unwrap_or_default()),
        };
        out.set_column(name, column_of(values));
    }
    Ok(out)
}

/// `priva acmg-step gene_mechanism_annotation annotate_gene_mechanism_categories`:
/// options `condition_cache`, `hgnc_table`, `clingen_dosage`, `loeuf_table`
/// (absent ones are the repository's), `symbol_col` and `output_col`.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    match function {
        "annotate_gene_mechanism_categories" => {
            let text = |name: &str, default: &str| -> Result<String> {
                match call.options.get(name) {
                    None => Ok(default.to_owned()),
                    Some(Json::Str(text)) => Ok(text.clone()),
                    Some(_) => Err(invalid(format!("the option {name} must be text"))),
                }
            };
            let symbol_col = text("symbol_col", "SYMBOL")?;
            let output_col = text("output_col", "var_plausible_patho_mechs")?;
            // The paths only; the hub itself is built inside, after the column check.
            let paths = harness_paths(&call.options)?;
            annotate_gene_mechanism_categories(&call.table, paths, &symbol_col, &output_col)
        }
        _ => Err(not_converted("gene_mechanism_annotation", function)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn history(mechanism: &str, inheritance: &str) -> History {
        History {
            mechanism: mechanism.to_owned(),
            condition: "OMIM:1".to_owned(),
            inheritance: inheritance.to_owned(),
            x_linked: false,
            penetrance: "",
        }
    }

    #[test]
    fn profile_tags_and_applicability() {
        let histories = vec![
            history("LOF", "dominant"),
            history("UNRESOLVED", "dominant"),
            history("UNRESOLVED", ""),
            history("DOMINANT_NEGATIVE", "x_linked_recessive"),
            history("GOF", "polygenic"),
            history("LOF", ""),
        ];
        assert_eq!(
            mechanism_profile_tag(&histories[3]),
            "x_linked_recessive_DN"
        );
        assert_eq!(mechanism_profile_tag(&histories[2]), "uncertain");
        assert_eq!(mechanism_profile_tag(&histories[4]), "");
        let applicability = classify_variant_applicability(&histories, [2, 1, 1]);
        assert_eq!(
            applicability.plausible,
            "dominant_LOF;x_linked_recessive_DN"
        );
        assert_eq!(applicability.applicable, "dominant_LOF;LOF");
        assert_eq!(
            applicability.uncertain,
            "dominant;uncertain;x_linked_recessive_DN"
        );
    }

    #[test]
    fn rows_of_numbers_only_become_floats() {
        let table = Table {
            names: vec!["a".to_owned(), "b".to_owned()],
            columns: vec![Column::Int(vec![1]), Column::Float(vec![2.5])],
            rows: 1,
        };
        assert_eq!(row_values(&table).unwrap(), RowValues::IntegersAsFloats);
        assert_eq!(
            row_cell(&table, RowValues::IntegersAsFloats, Some(0), 0),
            Some(Cell::Float(1.0))
        );
        let table = Table {
            names: vec!["a".to_owned(), "b".to_owned()],
            columns: vec![Column::Int(vec![1]), Column::Bool(vec![true])],
            rows: 1,
        };
        assert_eq!(row_values(&table).unwrap(), RowValues::AsStored);
        let table = Table {
            names: vec!["a".to_owned()],
            columns: vec![Column::NullableInt(vec![Some(1)])],
            rows: 1,
        };
        assert!(row_values(&table).is_err());
    }
}
