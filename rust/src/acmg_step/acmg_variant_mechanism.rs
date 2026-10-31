//! `scripts/acmg_variant_mechanism.py` (591 lines): the exact-allele columns and the adapters to
//! the existing mask modules.
//!
//! Thread T3 (exact-allele matching and the scoring tail). Converted: the two mask adapters
//! (`_variant_condition_masks`, `_variant_mechanism_masks`), which hand the table's columns to the
//! earlier mask modules in `crate::acmg_variant_mechanism` with their pandas dtypes, and
//! `annotate_exact_nonlof_variants` with its helpers (`_clean_text`, `_row_assembly_for_gof_lookup`,
//! `_merge_semicolon_values`, `_ensure_exact_nonlof_columns`, `_record_exact_nonlof_match`).
use super::gene_mechanism_common::{clean, clean_text};
use super::gene_mechanism_resources::harness_resolver;
use super::gene_mechanism_variants::{CuratedMatch, ExactVariantMatcher, text_option};
use super::types::{AcmgError, HarnessCall, NamedMasks, SymbolResolver, not_converted};
use crate::acmg_scoring::SourceError;
use crate::acmg_variant_mechanism::variant_mechanism_masks::numeric_text;
use crate::acmg_variant_mechanism::variant_mechanism_masks::{
    Cell as MaskCell, Dtype as MaskDtype, Input, InputColumn, Masks,
};
use crate::acmg_variant_mechanism::{
    variant_condition_masks as condition, variant_mechanism_masks as mechanism,
};
use crate::cli::{Result, invalid};
use crate::pandas_table::Cell;
use crate::pandas_table::{Column, ObjectValue, Table};
use std::collections::HashMap;

/// The columns `_variant_condition_masks` reads.
const CONDITION_COLUMNS: [&str; 2] = ["variant_inheritance", "variant_penetrance"];
/// The columns `_variant_mechanism_masks` reads (the scores only when present).
const MECHANISM_COLUMNS: [&str; 5] = [
    "var_plausible_patho_mechs",
    "variant_effect",
    "variant_lof_score",
    "variant_gof_score",
    "variant_dn_score",
];

/// `priva acmg-step acmg_variant_mechanism <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    match function {
        "_variant_condition_masks" => Ok(masks_table(
            variant_condition_masks(&call.table)?,
            call.table.rows,
        )),
        "_variant_mechanism_masks" => Ok(masks_table(
            variant_mechanism_masks(&call.table)?,
            call.table.rows,
        )),
        "annotate_exact_nonlof_variants" => {
            let resolver = harness_resolver(&call.options)?;
            let mut table = call.table.clone();
            annotate_exact_nonlof_variants(
                &mut table,
                None,
                text_option(&call.options, "mechanism_json")?,
                &resolver,
            )?;
            Ok(table)
        }
        _ => Err(not_converted("acmg_variant_mechanism", function)),
    }
}

/// The masks as a frame of bool columns, as `pd.DataFrame(masks)` gives.
fn masks_table(masks: NamedMasks, rows: usize) -> Table {
    let (names, columns) = masks
        .masks
        .into_iter()
        .map(|(name, values)| (name, Column::Bool(values)))
        .unzip();
    Table {
        names,
        columns,
        rows,
    }
}

/// `_variant_condition_masks(df)`: the 13 condition masks from
/// `variant_inheritance` and `variant_penetrance`, by their Python names. Used
/// by PVS1 (T4), PM2/BS1 and BS2/BP5 (T5) and PP1/BP2/PM3 (T7).
pub(crate) fn variant_condition_masks(table: &Table) -> Result<NamedMasks> {
    let input = mask_input(table, &CONDITION_COLUMNS)?;
    let masks = condition::variant_condition_masks(&input).map_err(source_error)?;
    named_masks(&masks, &condition::NAMES)
}

/// `_variant_mechanism_masks(df)`: the 27 mechanism masks, by their Python names.
pub(crate) fn variant_mechanism_masks(table: &Table) -> Result<NamedMasks> {
    let input = mask_input(table, &MECHANISM_COLUMNS)?;
    let masks = mechanism::variant_mechanism_masks(&input).map_err(source_error)?;
    named_masks(&masks, &mechanism::NAMES)
}

/// The columns a mask builder reads, in the typed cells of the earlier mask modules.
/// pandas' reader gives each column one name, so a name is never repeated here.
fn mask_input(table: &Table, names: &[&str]) -> Result<Input<usize>> {
    let mut columns = Vec::new();
    for name in names {
        if let Ok(column) = table.column(name) {
            columns.push(mask_column(name, column));
        }
    }
    Ok(Input {
        columns,
        index: (0..table.rows).collect(),
        index_names: vec![None],
    })
}

/// One column with its pandas dtype: numpy columns, `object` columns (NaN and
/// `None` kept apart as pandas holds them) and the nullable `Int64`/`boolean`.
fn mask_column(name: &str, column: &Column) -> InputColumn {
    let nan = || MaskCell::Float(f64::NAN);
    let (dtype, values) = match column {
        Column::Int(values) => (
            MaskDtype::Numpy,
            values
                .iter()
                .map(|value| MaskCell::Integer(*value))
                .collect(),
        ),
        Column::Float(values) => (
            MaskDtype::Numpy,
            values.iter().map(|value| MaskCell::Float(*value)).collect(),
        ),
        Column::Bool(values) => (
            MaskDtype::Numpy,
            values.iter().map(|value| MaskCell::Bool(*value)).collect(),
        ),
        Column::BoolOrMissing(values) => (
            MaskDtype::Object,
            values
                .iter()
                .map(|value| value.map_or_else(nan, MaskCell::Bool))
                .collect(),
        ),
        Column::Text(values) => (
            MaskDtype::Object,
            values
                .iter()
                .map(|value| match value {
                    Some(text) => MaskCell::Text(text.clone()),
                    None => nan(),
                })
                .collect(),
        ),
        Column::Object(values) => (
            MaskDtype::Object,
            values
                .iter()
                .map(|value| match value {
                    ObjectValue::Nan => nan(),
                    ObjectValue::None => MaskCell::Null,
                    ObjectValue::Text(text) => MaskCell::Text(text.clone()),
                    ObjectValue::Int(number) => MaskCell::Integer(*number),
                    ObjectValue::Float(number) => MaskCell::Float(*number),
                    ObjectValue::Bool(flag) => MaskCell::Bool(*flag),
                })
                .collect(),
        ),
        Column::NullableInt(values) => (
            MaskDtype::Nullable("Int64"),
            values
                .iter()
                .map(|value| value.map_or(MaskCell::Null, MaskCell::Integer))
                .collect(),
        ),
        Column::NullableBool(values) => (
            MaskDtype::Nullable("boolean"),
            values
                .iter()
                .map(|value| value.map_or(MaskCell::Null, MaskCell::Bool))
                .collect(),
        ),
    };
    InputColumn {
        name: name.to_owned(),
        dtype,
        values,
    }
}

/// A Python exception of the mask modules stays that exception; their
/// "retain Python" answers become a refusal naming the value.
fn source_error(error: SourceError) -> Box<dyn std::error::Error> {
    if error.0 == "UnfinishedParity" {
        invalid(format!("refused (not reproduced): {}", error.1))
    } else {
        Box::new(AcmgError::new(error.0, error.1))
    }
}

/// The polars frame of the mask modules as named boolean vectors, in the Python dict order.
fn named_masks<I>(masks: &Masks<I>, names: &[&str]) -> Result<NamedMasks> {
    let mut out = Vec::with_capacity(names.len());
    for name in names {
        let values = masks
            .frame
            .column(name)
            .and_then(|column| column.bool().cloned())
            .map_err(|error| invalid(format!("mask {name}: {error}")))?;
        if values.null_count() > 0 {
            return Err(invalid(format!("mask {name} holds missing values")));
        }
        out.push(((*name).to_owned(), values.into_no_null_iter().collect()));
    }
    Ok(NamedMasks { masks: out })
}

/// The 11 columns `annotate_exact_nonlof_variants` makes sure of (`EXACT_NONLOF_OUTPUT_COLS`).
const EXACT_NONLOF_OUTPUT_COLS: [&str; 11] = [
    "variant_lof_score",
    "variant_gof_score",
    "variant_dn_score",
    "variant_mechanism_exclusive",
    "variant_exact_mechanisms",
    "variant_exact_mechanism_evidence",
    "variant_mechanism_match_route",
    "variant_gof_tag",
    "gofcards_accession_id",
    "gofcards_variant_id",
    "gofcards_match_route",
];

/// `_merge_semicolon_values(*values)`: the distinct tokens between `;` and `,`, in order.
fn merge_semicolon_values(values: &[&str]) -> String {
    let mut merged: Vec<String> = Vec::new();
    for value in values {
        for token in clean_text(value).replace(',', ";").split(';') {
            let token = crate::python_text::python_strip(token);
            if !token.is_empty() && !merged.iter().any(|seen| seen == token) {
                merged.push(token.to_owned());
            }
        }
    }
    merged.join(";")
}

/// `pd.to_numeric(column, errors="coerce").fillna(0).astype(int)` for the column
/// types the step can hold here; booleans in an object or nullable column are refused.
fn integer_scores(name: &str, column: &Column) -> Result<Vec<i64>> {
    let not_finite = || -> Box<dyn std::error::Error> {
        Box::new(AcmgError::new(
            "IntCastingNaNError",
            "Cannot convert non-finite values (NA or inf) to integer",
        ))
    };
    // Past 2**53 pandas' choice between an integer and a float column decides the value.
    let float_to_int = |number: f64| -> Result<i64> {
        if number.is_nan() {
            return Ok(0);
        }
        if number.is_infinite() {
            return Err(not_finite());
        }
        if number.abs() >= 9_007_199_254_740_992.0 {
            return Err(invalid(format!(
                "refused (not reproduced): {name} holds a number past 2**53"
            )));
        }
        Ok(number.trunc() as i64)
    };
    let refuse = || {
        invalid(format!(
            "refused (not reproduced): to_numeric of {name} ({} with booleans)",
            column.dtype_name()
        ))
    };
    match column {
        Column::Int(values) => Ok(values.clone()),
        Column::Float(values) => values.iter().map(|value| float_to_int(*value)).collect(),
        Column::Bool(values) => Ok(values.iter().map(|value| i64::from(*value)).collect()),
        Column::NullableInt(values) => Ok(values.iter().map(|value| value.unwrap_or(0)).collect()),
        Column::Text(values) => values
            .iter()
            .map(|value| match value {
                Some(text) => float_to_int(numeric_text(text)),
                None => Ok(0),
            })
            .collect(),
        Column::Object(values) => values
            .iter()
            .map(|value| match value {
                ObjectValue::Nan | ObjectValue::None => Ok(0),
                ObjectValue::Int(number) => Ok(*number),
                ObjectValue::Float(number) => float_to_int(*number),
                ObjectValue::Text(text) => float_to_int(numeric_text(text)),
                ObjectValue::Bool(_) => Err(refuse()),
            })
            .collect(),
        Column::BoolOrMissing(_) | Column::NullableBool(_) => Err(refuse()),
    }
}

/// `column.fillna("").astype(str)`; a nullable column holding a missing value
/// raises pandas' TypeError.
fn text_values(column: &Column) -> Result<Vec<String>> {
    let nullable_missing = match column {
        Column::NullableInt(values) => values.iter().any(Option::is_none).then_some("Int64"),
        Column::NullableBool(values) => values.iter().any(Option::is_none).then_some("boolean"),
        _ => None,
    };
    if let Some(dtype) = nullable_missing {
        return Err(Box::new(AcmgError::new(
            "TypeError",
            format!("Invalid value '' for dtype {dtype}"),
        )));
    }
    Ok((0..column.len())
        .map(|row| match column.cell(row) {
            Cell::Missing => String::new(),
            Cell::Text(text) => text.to_owned(),
            Cell::Int(number) => number.to_string(),
            Cell::Float(number) => crate::python_text::float_repr(number),
            Cell::Bool(flag) => (if flag { "True" } else { "False" }).to_owned(),
        })
        .collect())
}

/// `column.fillna(False).astype(str).str.lower().isin({"1", "true", "yes"})`.
fn exclusive_flags(column: &Column) -> Vec<bool> {
    (0..column.len())
        .map(|row| match column.cell(row) {
            Cell::Missing | Cell::Float(_) => false,
            Cell::Bool(flag) => flag,
            Cell::Int(number) => number == 1,
            Cell::Text(text) => ["1", "true", "yes"].contains(&text.to_lowercase().as_str()),
        })
        .collect()
}

/// The 11 columns after `_ensure_exact_nonlof_columns`, held while matches are recorded.
struct ExactColumns {
    scores: [Vec<i64>; 3],
    exclusive: Vec<bool>,
    /// The 7 text columns, in `EXACT_NONLOF_OUTPUT_COLS` order.
    texts: Vec<Vec<String>>,
}

impl ExactColumns {
    fn text(&mut self, name: &str) -> &mut Vec<String> {
        let index = EXACT_NONLOF_OUTPUT_COLS[4..]
            .iter()
            .position(|candidate| *candidate == name)
            .expect("one of the text columns");
        &mut self.texts[index]
    }

    /// `_record_exact_nonlof_match(df, rows, match)`.
    fn record(&mut self, rows: &[usize], found: &CuratedMatch) {
        if !found.mechanism_exclusive() {
            return;
        }
        let route = clean_text(&found.match_route);
        let evidence = found.evidence_json();
        for (index, score) in found.scores.iter().enumerate() {
            for row in rows {
                self.scores[index][*row] = self.scores[index][*row].max(*score);
            }
        }
        let mechanisms = found.mechanisms.join(";");
        let gof = found.mechanisms.iter().any(|mechanism| mechanism == "GOF");
        let accession = clean_text(&found.gofcards_accession_id);
        let variant = clean_text(&found.gofcards_variant_id);
        for row in rows {
            let row = *row;
            self.exclusive[row] = true;
            let merge = |values: &mut Vec<String>, addition: &str| {
                values[row] = merge_semicolon_values(&[&values[row], addition]);
            };
            merge(self.text("variant_exact_mechanisms"), &mechanisms);
            self.text("variant_exact_mechanism_evidence")[row] = evidence.clone();
            if !route.is_empty() {
                merge(self.text("variant_mechanism_match_route"), &route);
            }
            if !gof {
                continue;
            }
            merge(self.text("variant_gof_tag"), "GOF");
            if !accession.is_empty() {
                merge(self.text("gofcards_accession_id"), &accession);
            }
            if !variant.is_empty() {
                merge(self.text("gofcards_variant_id"), &variant);
            }
            if !route.is_empty() {
                merge(self.text("gofcards_match_route"), &route);
            }
        }
    }
}

/// `_ensure_exact_nonlof_columns(df)`: the 11 columns exist with their types
/// (an absent one is added at the end, as pandas does).
fn ensure_exact_nonlof_columns(table: &mut Table) -> Result<ExactColumns> {
    let mut scores: [Vec<i64>; 3] = Default::default();
    for (index, name) in EXACT_NONLOF_OUTPUT_COLS[..3].iter().enumerate() {
        scores[index] = match table.column(name) {
            Ok(column) => integer_scores(name, column)?,
            Err(_) => vec![0; table.rows],
        };
        table.set_column(name, Column::Int(scores[index].clone()));
    }
    let exclusive = match table.column("variant_mechanism_exclusive") {
        Ok(column) => exclusive_flags(column),
        Err(_) => vec![false; table.rows],
    };
    table.set_column(
        "variant_mechanism_exclusive",
        Column::Bool(exclusive.clone()),
    );
    let mut texts = Vec::new();
    for name in &EXACT_NONLOF_OUTPUT_COLS[4..] {
        let values = match table.column(name) {
            Ok(column) => text_values(column)?,
            Err(_) => vec![String::new(); table.rows],
        };
        table.set_column(name, text_column(&values));
        texts.push(values);
    }
    let mut columns = ExactColumns {
        scores,
        exclusive,
        texts,
    };
    if let Ok(legacy) = table.column("gofcards_accession") {
        let merged: Vec<String> = (0..table.rows)
            .map(|row| {
                let legacy = clean(Some(legacy.cell(row)));
                merge_semicolon_values(&[&columns.texts[4][row], &legacy])
            })
            .collect();
        // pandas makes a column from an empty list float64.
        let merged_column = if merged.is_empty() {
            Column::Float(Vec::new())
        } else {
            text_column(&merged)
        };
        table.set_column("gofcards_accession_id", merged_column);
        *columns.text("gofcards_accession_id") = merged;
    }
    Ok(columns)
}

fn text_column(values: &[String]) -> Column {
    Column::Text(values.iter().map(|value| Some(value.clone())).collect())
}

/// `_row_assembly_for_gof_lookup(row)`: the first of `assembly`, `genome_build`,
/// `reference_build`, `build` with a value, else `hg19`.
fn row_assembly(table: &Table, row: usize) -> String {
    for name in ["assembly", "genome_build", "reference_build", "build"] {
        if let Ok(column) = table.column(name) {
            let value = clean(Some(column.cell(row)));
            if !value.is_empty() {
                return value;
            }
        }
    }
    "hg19".to_owned()
}

/// A lookup error Python catches (`except Exception: continue`) is skipped;
/// a refusal stops the step.
fn skip_python_error(outcome: Result<CuratedMatch>) -> Result<Option<CuratedMatch>> {
    match outcome {
        Ok(found) => Ok(Some(found)),
        Err(error) if error.downcast_ref::<AcmgError>().is_some() => Ok(None),
        Err(error) => Err(error),
    }
}

/// `annotate_exact_nonlof_variants(df, row_mask, mechanism_json=...)`: the exact
/// curated GoF/DN matches of the non-LOF mechanism JSON, by HGVSp and then by
/// genomic allele, written into the 11 exact-allele columns. `row_mask` stands for
/// an already boolean mask; `resolver` is the hub's HGNC resolver. The table is
/// changed in place, as the Python function changes its frame.
pub(crate) fn annotate_exact_nonlof_variants(
    table: &mut Table,
    row_mask: Option<&[bool]>,
    mechanism_json: &str,
    resolver: &dyn SymbolResolver,
) -> Result<()> {
    if !table.has_column("SYMBOL") {
        return Ok(());
    }
    let mut columns = ensure_exact_nonlof_columns(table)?;
    let symbols = table.column("SYMBOL")?;
    let mask: Vec<bool> = (0..table.rows)
        .map(|row| {
            row_mask.is_none_or(|mask| mask[row]) && !clean(Some(symbols.cell(row))).is_empty()
        })
        .collect();
    if !mask.iter().any(|keep| *keep) {
        return Ok(());
    }
    let mut matcher = ExactVariantMatcher::new(mechanism_json, resolver);

    if let Ok(hgvsp) = table.column("HGVSp") {
        let mut groups: Vec<((String, String), Vec<usize>)> = Vec::new();
        let mut place: HashMap<(String, String), usize> = HashMap::new();
        for row in (0..table.rows).filter(|row| mask[*row]) {
            let protein = clean(Some(hgvsp.cell(row)));
            if protein.is_empty() {
                continue;
            }
            let key = (clean(Some(symbols.cell(row))), protein);
            match place.get(&key) {
                Some(index) => groups[*index].1.push(row),
                None => {
                    place.insert(key.clone(), groups.len());
                    groups.push((key, vec![row]));
                }
            }
        }
        for ((symbol, protein), rows) in &groups {
            let outcome = matcher.match_curated_nonlof_variant(
                Cell::Text(symbol),
                Some(Cell::Text(protein)),
                None,
                [None; 4],
                "hg38",
                "auto",
            );
            if let Some(found) = skip_python_error(outcome)? {
                columns.record(rows, &found);
            }
        }
    }

    let genomic_names = ["chrom", "pos", "ref", "alt"];
    if genomic_names.iter().all(|name| table.has_column(name)) {
        let genomic: Vec<&Column> = genomic_names
            .iter()
            .map(|name| table.column(name))
            .collect::<Result<_>>()?;
        // Rows with equal values form one pandas group; a group's lookup is the
        // same for each of its rows, so rows are grouped by the values' texts.
        let mut groups: Vec<(usize, String, Vec<usize>)> = Vec::new();
        let mut place: HashMap<Vec<String>, usize> = HashMap::new();
        for row in 0..table.rows {
            if !mask[row] || columns.exclusive[row] {
                continue;
            }
            if genomic
                .iter()
                .any(|column| clean(Some(column.cell(row))).is_empty())
            {
                continue;
            }
            let assembly = row_assembly(table, row);
            let mut key: Vec<String> = std::iter::once(symbols)
                .chain(genomic.iter().copied())
                .map(|column| format!("{:?}", column.cell(row)))
                .collect();
            key.push(assembly.clone());
            match place.get(&key) {
                Some(index) => groups[*index].2.push(row),
                None => {
                    place.insert(key, groups.len());
                    groups.push((row, assembly, vec![row]));
                }
            }
        }
        for (first, assembly, rows) in &groups {
            let cell = |index: usize| Some(genomic[index].cell(*first));
            let outcome = matcher.match_curated_nonlof_variant(
                symbols.cell(*first),
                None,
                None,
                [cell(0), cell(1), cell(2), cell(3)],
                assembly,
                "auto",
            );
            if let Some(found) = skip_python_error(outcome)? {
                columns.record(rows, &found);
            }
        }
    }

    for (index, name) in EXACT_NONLOF_OUTPUT_COLS[..3].iter().enumerate() {
        table.set_column(name, Column::Int(columns.scores[index].clone()));
    }
    table.set_column(
        "variant_mechanism_exclusive",
        Column::Bool(columns.exclusive.clone()),
    );
    for (index, name) in EXACT_NONLOF_OUTPUT_COLS[4..].iter().enumerate() {
        table.set_column(name, text_column(&columns.texts[index]));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table(columns: Vec<(&str, Column)>) -> Table {
        let rows = columns.first().map_or(0, |(_, column)| column.len());
        let (names, columns) = columns
            .into_iter()
            .map(|(name, column)| (name.to_owned(), column))
            .unzip();
        Table {
            names,
            columns,
            rows,
        }
    }

    fn text(values: &[Option<&str>]) -> Column {
        Column::Text(values.iter().map(|v| v.map(str::to_owned)).collect())
    }

    #[test]
    fn condition_masks_follow_the_semicolon_tokens() {
        let input = table(vec![
            (
                "variant_inheritance",
                text(&[
                    Some("Recessive;dominant"),
                    None,
                    Some("x_linked_dominant"),
                    Some("polygenic;non_mendelian"),
                ]),
            ),
            (
                "variant_penetrance",
                text(&[
                    Some("incomplete"),
                    Some("complete"),
                    None,
                    Some("incomplete;complete"),
                ]),
            ),
        ]);
        let masks = variant_condition_masks(&input).unwrap();
        assert_eq!(masks.masks.len(), 13);
        assert_eq!(masks.masks[0].0, "has_recessive");
        assert_eq!(
            masks.get("has_recessive").unwrap(),
            &[true, false, false, false]
        );
        assert_eq!(
            masks.get("has_dominant_like").unwrap(),
            &[true, false, true, false]
        );
        assert_eq!(
            masks.get("has_non_monogenic").unwrap(),
            &[false, false, false, true]
        );
        assert_eq!(
            masks.get("has_non_mendelian").unwrap(),
            &[false, false, false, true]
        );
        assert_eq!(
            masks.get("has_complete_penetrance").unwrap(),
            &[false, true, false, true]
        );
    }

    #[test]
    fn a_missing_condition_column_is_the_python_key_error() {
        let input = table(vec![("variant_penetrance", text(&[Some("complete")]))]);
        let error = variant_condition_masks(&input).unwrap_err().to_string();
        assert_eq!(
            error,
            "KeyError: 'variant-level condition annotations are required; missing columns: variant_inheritance'"
        );
    }

    #[test]
    fn mechanism_masks_read_profile_effect_and_scores() {
        let input = table(vec![
            (
                "var_plausible_patho_mechs",
                text(&[Some("dominant_GOF;recessive_LOF"), Some(" "), None]),
            ),
            (
                "variant_effect",
                text(&[
                    Some("exact_known_GOF"),
                    Some("uncertain"),
                    Some("predicted_LOF_high_confidence"),
                ]),
            ),
            ("variant_lof_score", Column::Int(vec![0, 1, 2])),
            ("variant_dn_score", Column::Float(vec![2.0, f64::NAN, 0.0])),
        ]);
        let masks = variant_mechanism_masks(&input).unwrap();
        assert_eq!(masks.masks.len(), 27);
        assert_eq!(masks.get("modern_profile").unwrap(), &[true, false, false]);
        assert_eq!(masks.get("is_exact_gof").unwrap(), &[true, false, false]);
        assert_eq!(masks.get("is_exact_dn").unwrap(), &[true, false, false]);
        assert_eq!(masks.get("is_uncertain").unwrap(), &[false, true, false]);
        assert_eq!(
            masks.get("is_predicted_lof").unwrap(),
            &[false, false, true]
        );
        assert_eq!(
            masks.get("has_dom_gof_history").unwrap(),
            &[true, false, false]
        );
        assert_eq!(
            masks.get("has_recessive_compatible").unwrap(),
            &[true, false, false]
        );
    }

    #[test]
    fn nullable_columns_keep_the_pandas_errors() {
        let input = table(vec![
            (
                "var_plausible_patho_mechs",
                Column::NullableInt(vec![Some(1), None]),
            ),
            ("variant_effect", text(&[None, None])),
        ]);
        let error = variant_mechanism_masks(&input).unwrap_err().to_string();
        assert_eq!(error, "TypeError: Invalid value '' for dtype Int64");
    }

    #[test]
    fn semicolon_values_merge_in_order() {
        assert_eq!(
            merge_semicolon_values(&["A;B", "b,A ; C", "nan", ""]),
            "A;B;b;C"
        );
        assert_eq!(merge_semicolon_values(&["", " "]), "");
    }

    #[test]
    fn exact_columns_are_ensured_with_their_types() {
        let mut input = table(vec![
            ("SYMBOL", text(&[None, Some("nan")])),
            ("variant_gof_score", text(&[Some(" 2 "), Some("x")])),
            ("variant_mechanism_exclusive", text(&[Some("Yes"), None])),
            ("variant_gof_tag", Column::Float(vec![1.5, f64::NAN])),
            ("gofcards_accession", text(&[Some("1,2"), None])),
        ]);
        struct Never;
        impl SymbolResolver for Never {
            fn resolve_symbol(&self, _: Cell<'_>) -> Result<String> {
                panic!("no symbol is looked up")
            }
            fn try_resolve_symbol(&self, _: Cell<'_>) -> String {
                panic!("no symbol is looked up")
            }
        }
        annotate_exact_nonlof_variants(&mut input, None, "/no/such.json", &Never).unwrap();
        assert_eq!(
            input.column("variant_gof_score").unwrap(),
            &Column::Int(vec![2, 0])
        );
        assert_eq!(
            input.column("variant_lof_score").unwrap(),
            &Column::Int(vec![0, 0])
        );
        assert_eq!(
            input.column("variant_mechanism_exclusive").unwrap(),
            &Column::Bool(vec![true, false])
        );
        assert_eq!(
            input.column("variant_gof_tag").unwrap(),
            &text(&[Some("1.5"), Some("")])
        );
        assert_eq!(
            input.column("gofcards_accession_id").unwrap(),
            &text(&[Some("1;2"), Some("")])
        );
        assert_eq!(input.names.len(), 13);
    }
}
