//! `scripts/acmg_variant_mechanism.py` (591 lines): the exact-allele columns and the adapters to
//! the existing mask modules.
//!
//! Thread T3 (exact-allele matching and the scoring tail). Converted: the two mask adapters
//! (`_variant_condition_masks`, `_variant_mechanism_masks`), which hand the table's columns to the
//! earlier mask modules in `crate::acmg_variant_mechanism` with their pandas dtypes. Still to
//! convert: _clean_text, _row_assembly_for_gof_lookup, _merge_semicolon_values,
//! _ensure_exact_nonlof_columns, _record_exact_nonlof_match, annotate_exact_nonlof_variants.
use super::types::{AcmgError, HarnessCall, NamedMasks, not_converted};
use crate::acmg_scoring::SourceError;
use crate::acmg_variant_mechanism::variant_mechanism_masks::{
    Cell as MaskCell, Dtype as MaskDtype, Input, InputColumn, Masks,
};
use crate::acmg_variant_mechanism::{
    variant_condition_masks as condition, variant_mechanism_masks as mechanism,
};
use crate::cli::{Result, invalid};
use crate::pandas_table::{Column, ObjectValue, Table};

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
}
