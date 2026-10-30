//! pandas columns as polars Series and back, and the polars forms of the
//! pandas calls the ACMG step uses most, with pandas' rules for missing values.
//!
//! A column moves to polars only when polars can hold it without changing
//! what the pandas code computes: `int64`, `float64`, `bool`, the nullable
//! `Int64` and `boolean`, an `object` column of text and missing values and an
//! `object` column of booleans and missing values. A missing value becomes a
//! polars null. An `object` column that mixes kinds of values stays in
//! `pandas_table` and is handled by own code (`to_series` refuses it).
//!
//! What differs between the two libraries, and how the functions below keep
//! pandas' result (each point checked against pandas 2.2.3 and polars 0.51.0):
//!
//! - polars orders NaN above every number and counts `NaN == NaN` as true;
//!   pandas compares NaN false (and `!=` true). Floats therefore go to polars
//!   with NaN as null, and the comparisons fill nulls the way pandas answers.
//!   Arithmetic in polars can make NaN again (`0.0 / 0.0`): pass such a result
//!   through `nan_to_null` before comparing it.
//! - polars divides integers with floor division and gives null for a zero
//!   divisor; pandas `/` is true division. Cast integers to `Float64` first.
//! - `isin`: pandas answers false for a missing value unless the list holds
//!   NaN; polars answers null. `isin_text` and `isin_numbers` fill accordingly.
//! - `str.contains`: both leave missing values missing; pandas uses Python's
//!   `re`, where `$` also matches before a final line break; `str_contains`
//!   refuses a pattern with `$` when a value holds a line break.
//! - `np.where` on a mask with missing values takes the first choice (NaN is
//!   true in Python) while polars `zip_with` takes the second; `where_mask`
//!   refuses masks with missing values.
//! - `groupby(sort=False)` and polars `group_by_stable` both list groups in
//!   order of first appearance; pandas drops groups whose key is missing
//!   unless `dropna=False`, polars keeps them (`first_appearance_groups`).
//! - `drop_duplicates(keep="first")` and polars `unique_stable` agree, both
//!   counting `-0.0` equal to `0.0` and missing values equal to each other.
//! - `sort_values` on several keys is a stable sort with missing values last;
//!   polars with `maintain_order` and `nulls_last` gives the same order. A
//!   single-key `sort_values` uses an unstable quicksort in pandas, so ties
//!   there need own care.
// The group 3 threads are the callers of most of these functions.
#![allow(dead_code)]
use crate::cli::{Result, invalid};
use crate::pandas_table::{Column as PandasColumn, ObjectValue, Table};
use polars_core::prelude::*;
use polars_ops::chunked_array::strings::StringNameSpaceImpl;

/// The pandas column type a polars result goes back to.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum PandasDtype {
    Int64,
    Float64,
    Bool,
    /// `object` column of `True`, `False` and NaN.
    ObjectBool,
    /// `object` column of text and NaN.
    ObjectText,
    /// pandas' nullable `Int64`.
    NullableInt64,
    /// pandas' nullable `boolean`.
    NullableBoolean,
}

fn polars_error(error: PolarsError) -> Box<dyn std::error::Error> {
    invalid(format!("polars: {error}"))
}

/// The column as a polars Series named `name`; missing values become nulls.
pub(crate) fn to_series(column: &PandasColumn, name: &str) -> Result<Series> {
    Ok(match column {
        PandasColumn::Int(values) => Series::new(name.into(), values),
        PandasColumn::Float(values) => Series::new(
            name.into(),
            values
                .iter()
                .map(|value| (!value.is_nan()).then_some(*value))
                .collect::<Vec<Option<f64>>>(),
        ),
        PandasColumn::Bool(values) => Series::new(name.into(), values),
        PandasColumn::BoolOrMissing(values) | PandasColumn::NullableBool(values) => {
            Series::new(name.into(), values)
        }
        PandasColumn::NullableInt(values) => Series::new(name.into(), values),
        PandasColumn::Text(values) => Series::new(
            name.into(),
            values
                .iter()
                .map(|value| value.as_deref())
                .collect::<Vec<Option<&str>>>(),
        ),
        PandasColumn::Object(_) => {
            return Err(invalid(format!(
                "column {name} is a pandas object column that mixes kinds of values; \
                 it has no polars equivalent and needs own code"
            )));
        }
    })
}

/// A polars Series back as a pandas column of type `dtype`. Nulls become NaN
/// (`Float64`, the object types) or `pd.NA` (the nullable types); a null where
/// the type cannot hold a missing value (`Int64`, `Bool`) is an error.
pub(crate) fn from_series(series: &Series, dtype: PandasDtype) -> Result<PandasColumn> {
    let name = series.name().to_string();
    let no_nulls = |what: &str| -> Result<()> {
        if series.null_count() > 0 {
            return Err(invalid(format!(
                "column {name} has missing values, which a pandas {what} column cannot hold"
            )));
        }
        Ok(())
    };
    Ok(match dtype {
        PandasDtype::Int64 => {
            no_nulls("int64")?;
            let values = series.strict_cast(&DataType::Int64).map_err(polars_error)?;
            PandasColumn::Int(
                values
                    .i64()
                    .map_err(polars_error)?
                    .into_no_null_iter()
                    .collect(),
            )
        }
        PandasDtype::Float64 => {
            let values = series
                .strict_cast(&DataType::Float64)
                .map_err(polars_error)?;
            PandasColumn::Float(
                values
                    .f64()
                    .map_err(polars_error)?
                    .into_iter()
                    .map(|value| value.unwrap_or(f64::NAN))
                    .collect(),
            )
        }
        PandasDtype::Bool => {
            no_nulls("bool")?;
            PandasColumn::Bool(
                series
                    .bool()
                    .map_err(polars_error)?
                    .into_no_null_iter()
                    .collect(),
            )
        }
        PandasDtype::ObjectBool => {
            PandasColumn::BoolOrMissing(series.bool().map_err(polars_error)?.into_iter().collect())
        }
        PandasDtype::NullableBoolean => {
            PandasColumn::NullableBool(series.bool().map_err(polars_error)?.into_iter().collect())
        }
        PandasDtype::NullableInt64 => {
            let values = series.strict_cast(&DataType::Int64).map_err(polars_error)?;
            PandasColumn::NullableInt(values.i64().map_err(polars_error)?.into_iter().collect())
        }
        PandasDtype::ObjectText => PandasColumn::Text(
            series
                .str()
                .map_err(polars_error)?
                .into_iter()
                .map(|value| value.map(str::to_owned))
                .collect(),
        ),
    })
}

/// The values of a mask that has no missing values.
pub(crate) fn mask_values(mask: &BooleanChunked) -> Result<Vec<bool>> {
    if mask.null_count() > 0 {
        return Err(invalid(format!(
            "mask {} has missing values where pandas would have a plain boolean",
            mask.name()
        )));
    }
    Ok(mask.into_no_null_iter().collect())
}

/// Float results of polars arithmetic with NaN turned into null, as pandas
/// counts NaN as missing.
pub(crate) fn nan_to_null(series: &Series) -> Result<Series> {
    let values = series
        .strict_cast(&DataType::Float64)
        .map_err(polars_error)?;
    Ok(Series::new(
        series.name().clone(),
        values
            .f64()
            .map_err(polars_error)?
            .into_iter()
            .map(|value| value.filter(|number| !number.is_nan()))
            .collect::<Vec<Option<f64>>>(),
    ))
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Comparison {
    Equal,
    NotEqual,
    Greater,
    GreaterOrEqual,
    Less,
    LessOrEqual,
}

/// The answer pandas gives where a value is missing: false, except for `!=`.
fn fill_missing_comparison(mask: BooleanChunked, comparison: Comparison) -> Result<BooleanChunked> {
    mask.fill_null_with_values(comparison == Comparison::NotEqual)
        .map_err(polars_error)
}

/// `series <op> number` for a numeric column, as pandas answers it (a plain
/// boolean mask). Integers are compared as doubles, as numpy compares an
/// int64 column with a float.
pub(crate) fn compare_number(
    series: &Series,
    comparison: Comparison,
    number: f64,
) -> Result<BooleanChunked> {
    let values = nan_to_null(series)?;
    let mask = match comparison {
        Comparison::Equal => values.equal(number),
        Comparison::NotEqual => values.not_equal(number),
        Comparison::Greater => values.gt(number),
        Comparison::GreaterOrEqual => values.gt_eq(number),
        Comparison::Less => values.lt(number),
        Comparison::LessOrEqual => values.lt_eq(number),
    }
    .map_err(polars_error)?;
    fill_missing_comparison(mask, comparison)
}

/// `left <op> right` for two numeric columns of the same length, as pandas answers it.
pub(crate) fn compare_numbers(
    left: &Series,
    comparison: Comparison,
    right: &Series,
) -> Result<BooleanChunked> {
    let left = nan_to_null(left)?;
    let right = nan_to_null(right)?;
    let mask = match comparison {
        Comparison::Equal => left.equal(&right),
        Comparison::NotEqual => left.not_equal(&right),
        Comparison::Greater => left.gt(&right),
        Comparison::GreaterOrEqual => left.gt_eq(&right),
        Comparison::Less => left.lt(&right),
        Comparison::LessOrEqual => left.lt_eq(&right),
    }
    .map_err(polars_error)?;
    fill_missing_comparison(mask, comparison)
}

/// `series == text` / `series != text` for a text column (a missing value is
/// never equal to the text).
pub(crate) fn compare_text(
    series: &Series,
    comparison: Comparison,
    text: &str,
) -> Result<BooleanChunked> {
    let values = series.str().map_err(polars_error)?;
    let mask = match comparison {
        Comparison::Equal => values.equal(text),
        Comparison::NotEqual => values.not_equal(text),
        _ => return Err(invalid("text columns are only compared for equality here")),
    };
    fill_missing_comparison(mask, comparison)
}

fn as_list(values: Series) -> Series {
    Series::new("values".into(), &[values])
}

/// `series.isin(values)` for a text column and a list of texts.
pub(crate) fn isin_text(series: &Series, values: &[&str]) -> Result<BooleanChunked> {
    let wanted = as_list(Series::new("values".into(), values));
    polars_ops::series::is_in(series, &wanted, false)
        .and_then(|mask| mask.fill_null_with_values(false))
        .map_err(polars_error)
}

/// `series.isin(values)` for a numeric column; NaN in `values` matches the
/// missing values, as it does in pandas.
pub(crate) fn isin_numbers(series: &Series, values: &[f64]) -> Result<BooleanChunked> {
    let column = nan_to_null(series)?;
    let numbers: Vec<f64> = values.iter().copied().filter(|v| !v.is_nan()).collect();
    let wanted = as_list(Series::new("values".into(), numbers));
    let mask = polars_ops::series::is_in(&column, &wanted, false)
        .and_then(|mask| mask.fill_null_with_values(false))
        .map_err(polars_error)?;
    if values.iter().any(|value| value.is_nan()) {
        Ok(&mask | &column.is_null())
    } else {
        Ok(mask)
    }
}

/// `series.str.contains(pattern, regex=regex)` for a text column. The result
/// keeps missing values missing (pandas' `na=None`); fill them as the Python
/// code does (`na=False`, `.fillna(False)`).
pub(crate) fn str_contains(series: &Series, pattern: &str, regex: bool) -> Result<BooleanChunked> {
    let values = series.str().map_err(polars_error)?;
    if regex
        && pattern.contains('$')
        && values.into_iter().flatten().any(|text| text.contains('\n'))
    {
        return Err(invalid(format!(
            "pattern {pattern:?} on {}: Python's `$` also matches before a final line \
             break, which the polars pattern does not",
            series.name()
        )));
    }
    if regex {
        values.contains(pattern, true)
    } else {
        values.contains_literal(pattern)
    }
    .map_err(polars_error)
}

/// `np.where(mask, then, otherwise)` for a mask without missing values.
pub(crate) fn where_mask(
    mask: &BooleanChunked,
    then: &Series,
    otherwise: &Series,
) -> Result<Series> {
    if mask.null_count() > 0 {
        return Err(invalid(
            "np.where on a mask with missing values picks the first choice there; \
             fill the mask first",
        ));
    }
    then.zip_with(mask, otherwise).map_err(polars_error)
}

/// `series.astype(bool)` for a float column: NaN counts as true, as in Python.
pub(crate) fn float_astype_bool(series: &Series) -> Result<BooleanChunked> {
    let values = nan_to_null(series)?;
    let nonzero = values.not_equal(0.0).map_err(polars_error)?;
    nonzero.fill_null_with_values(true).map_err(polars_error)
}

/// `series.astype(int)` for a float column: truncates toward zero, and fails
/// on a missing or infinite value as pandas does. Values outside the 64-bit
/// range are refused (numpy gives an undefined number there).
pub(crate) fn float_astype_int(series: &Series) -> Result<Vec<i64>> {
    let values = series
        .strict_cast(&DataType::Float64)
        .map_err(polars_error)?;
    let mut out = Vec::with_capacity(values.len());
    for value in values.f64().map_err(polars_error)?.into_iter() {
        match value {
            Some(number) if number.is_finite() => {
                let truncated = number.trunc();
                if !(-9.223_372_036_854_776e18..9.223_372_036_854_776e18).contains(&truncated) {
                    return Err(invalid(format!(
                        "{number} in {} is outside the 64-bit integer range",
                        series.name()
                    )));
                }
                out.push(truncated as i64);
            }
            _ => {
                return Err(invalid(format!(
                    "IntCastingNaNError: Cannot convert non-finite values (NA or inf) to integer ({})",
                    series.name()
                )));
            }
        }
    }
    Ok(out)
}

/// Builds a polars frame of the named columns plus a row-number column.
fn keyed_frame(table: &Table, keys: &[&str]) -> Result<(DataFrame, Vec<String>)> {
    let mut columns: Vec<Column> = Vec::with_capacity(keys.len() + 1);
    let mut names = Vec::with_capacity(keys.len());
    for (position, key) in keys.iter().enumerate() {
        // Own names, so that a key that repeats or clashes stays one column each.
        let name = format!("key_{position}");
        let series = match table.column(key)? {
            PandasColumn::Object(values) => object_key_series(&name, values),
            column => to_series(column, &name)?,
        };
        columns.push(series.into_column());
        names.push(name);
    }
    let rows: Vec<IdxSize> = (0..table.rows as IdxSize).collect();
    columns.push(Series::new("row_number".into(), rows).into_column());
    Ok((DataFrame::new(columns).map_err(polars_error)?, names))
}

/// An object column as text keys that are equal exactly when pandas' hashing
/// counts the values equal: every missing value is one null, numbers and
/// booleans compare by value (`1 == 1.0 == True`, `-0.0 == 0.0`), text by text.
fn object_key_series(name: &str, values: &[ObjectValue]) -> Series {
    let number = |value: f64| {
        let value = if value == 0.0 { 0.0 } else { value };
        Some(format!("n{:016x}", value.to_bits()))
    };
    let keys: Vec<Option<String>> = values
        .iter()
        .map(|value| match value {
            ObjectValue::Nan | ObjectValue::None => None,
            ObjectValue::Float(x) if x.is_nan() => None,
            ObjectValue::Float(x) => number(*x),
            ObjectValue::Int(x) => number(*x as f64),
            ObjectValue::Bool(flag) => number(if *flag { 1.0 } else { 0.0 }),
            ObjectValue::Text(text) => Some(format!("s{text}")),
        })
        .collect();
    Series::new(name.into(), keys)
}

fn row_numbers(frame: &DataFrame) -> Result<Vec<usize>> {
    Ok(frame
        .column("row_number")
        .map_err(polars_error)?
        .as_materialized_series()
        .idx()
        .map_err(polars_error)?
        .into_no_null_iter()
        .map(|row| row as usize)
        .collect())
}

/// `df.drop_duplicates(subset=keys, keep="first")` (all columns when `keys`
/// is `None`): the rows kept, in their original order.
pub(crate) fn drop_duplicates_first(table: &Table, keys: Option<&[&str]>) -> Result<Vec<usize>> {
    let all: Vec<&str> = table.names.iter().map(String::as_str).collect();
    let mut distinct = all.clone();
    distinct.sort_unstable();
    distinct.dedup();
    if keys.is_none() && distinct.len() != all.len() {
        return Err(invalid(
            "drop_duplicates over a table with repeated column names is not reproduced",
        ));
    }
    let (frame, names) = keyed_frame(table, keys.unwrap_or(&all[..]))?;
    let unique = frame
        .unique_stable(Some(names.as_slice()), UniqueKeepStrategy::First, None)
        .map_err(polars_error)?;
    row_numbers(&unique)
}

/// `df.groupby(keys, sort=False, dropna=dropna)`: the rows of each group, the
/// groups in order of first appearance.
pub(crate) fn first_appearance_groups(
    table: &Table,
    keys: &[&str],
    dropna: bool,
) -> Result<Vec<Vec<usize>>> {
    let (frame, names) = keyed_frame(table, keys)?;
    let grouped = frame
        .group_by_stable(names.iter().map(String::as_str))
        .map_err(polars_error)?;
    let key_missing: Vec<Vec<bool>> = keys
        .iter()
        .map(|key| table.column(key).map(PandasColumn::missing_mask))
        .collect::<Result<_>>()?;
    let mut groups = Vec::new();
    for group in grouped.get_groups().iter() {
        let rows: Vec<usize> = match group {
            GroupsIndicator::Idx((_, all)) => all.iter().map(|row| *row as usize).collect(),
            GroupsIndicator::Slice([first, len]) => {
                (first..first + len).map(|row| row as usize).collect()
            }
        };
        if dropna && key_missing.iter().any(|missing| missing[rows[0]]) {
            continue;
        }
        groups.push(rows);
    }
    Ok(groups)
}

/// `df.sort_values(keys, ascending=[...])` with several keys: a stable sort,
/// missing values last. Returns the rows in their new order.
pub(crate) fn stable_sort(table: &Table, keys: &[&str], descending: &[bool]) -> Result<Vec<usize>> {
    if keys.len() != descending.len() {
        return Err(invalid("one sort direction per key"));
    }
    if keys
        .iter()
        .any(|key| matches!(table.column(key), Ok(PandasColumn::Object(_))))
    {
        return Err(invalid(
            "sorting by an object column that mixes kinds of values is not reproduced",
        ));
    }
    let (frame, names) = keyed_frame(table, keys)?;
    let sorted = frame
        .sort(
            names.iter().map(String::as_str).collect::<Vec<&str>>(),
            SortMultipleOptions::default()
                .with_order_descending_multi(descending.to_vec())
                .with_nulls_last(true)
                .with_maintain_order(true),
        )
        .map_err(polars_error)?;
    row_numbers(&sorted)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn values(mask: &BooleanChunked) -> Vec<Option<bool>> {
        mask.into_iter().collect()
    }

    fn table(columns: Vec<(&str, PandasColumn)>) -> Table {
        let rows = columns[0].1.len();
        Table {
            names: columns.iter().map(|(name, _)| name.to_string()).collect(),
            columns: columns.into_iter().map(|(_, column)| column).collect(),
            rows,
        }
    }

    fn text(values: &[Option<&str>]) -> PandasColumn {
        PandasColumn::Text(values.iter().map(|v| v.map(str::to_owned)).collect())
    }

    // The expected values below are what pandas 2.2.3 printed for the same input
    // (group 3 probe, `probe_pandas_semantics.py`).
    #[test]
    fn comparisons_answer_missing_values_like_pandas() {
        let column = PandasColumn::Float(vec![1.0, f64::NAN, 0.5, -0.0, 0.0]);
        let series = to_series(&column, "x").unwrap();
        let some = |flags: &[bool]| flags.iter().map(|f| Some(*f)).collect::<Vec<_>>();
        assert_eq!(
            values(&compare_number(&series, Comparison::Greater, 0.6).unwrap()),
            some(&[true, false, false, false, false])
        );
        assert_eq!(
            values(&compare_number(&series, Comparison::Equal, f64::NAN).unwrap()),
            some(&[false, false, false, false, false])
        );
        assert_eq!(
            values(&compare_number(&series, Comparison::NotEqual, 0.5).unwrap()),
            some(&[true, true, false, true, true])
        );
        assert_eq!(
            values(&compare_number(&series, Comparison::Equal, 0.0).unwrap()),
            some(&[false, false, false, true, true])
        );
        let zero = to_series(&PandasColumn::Float(vec![0.0, 1.0, 0.0, 2.0, 0.0]), "z").unwrap();
        let quotient = nan_to_null(&(&series / &zero).unwrap()).unwrap();
        // pandas: [inf, nan, inf, -0.0, nan]
        let quotient = from_series(&quotient, PandasDtype::Float64).unwrap();
        let PandasColumn::Float(numbers) = quotient else {
            panic!()
        };
        assert_eq!(numbers[0], f64::INFINITY);
        assert!(numbers[1].is_nan() && numbers[4].is_nan());
        assert!(numbers[3] == 0.0 && numbers[3].is_sign_negative());
        let left = to_series(&PandasColumn::Float(vec![1.0, f64::NAN, 2.0]), "l").unwrap();
        let right = to_series(&PandasColumn::Int(vec![1, 1, 3]), "r").unwrap();
        assert_eq!(
            values(&compare_numbers(&left, Comparison::NotEqual, &right).unwrap()),
            some(&[false, true, true])
        );
    }

    #[test]
    fn isin_matches_pandas() {
        let column = text(&[
            Some("Pathogenic"),
            None,
            Some("Likely_pathogenic"),
            Some("benign"),
            Some(""),
        ]);
        let series = to_series(&column, "t").unwrap();
        assert_eq!(
            mask_values(&isin_text(&series, &["benign", "Pathogenic"]).unwrap()).unwrap(),
            [true, false, false, true, false]
        );
        let numbers = to_series(
            &PandasColumn::Float(vec![1.0, f64::NAN, 0.5, -0.0, 0.0]),
            "x",
        )
        .unwrap();
        assert_eq!(
            mask_values(&isin_numbers(&numbers, &[0.0]).unwrap()).unwrap(),
            [false, false, false, true, true]
        );
        assert_eq!(
            mask_values(&isin_numbers(&numbers, &[f64::NAN]).unwrap()).unwrap(),
            [false, true, false, false, false]
        );
    }

    #[test]
    fn text_tests_match_pandas() {
        let column = text(&[
            Some("Pathogenic"),
            None,
            Some("Likely_pathogenic"),
            Some("benign"),
            Some(""),
        ]);
        let series = to_series(&column, "t").unwrap();
        // pandas: [True, nan, True, False, False] (object), so missing stays missing.
        assert_eq!(
            values(&str_contains(&series, "athogenic", true).unwrap()),
            [Some(true), None, Some(true), Some(false), Some(false)]
        );
        assert_eq!(
            values(&str_contains(&series, "^p", true).unwrap()),
            [Some(false), None, Some(false), Some(false), Some(false)]
        );
        assert!(str_contains(&series, "(", true).is_err());
        let newline = to_series(&text(&[Some("ab\n")]), "n").unwrap();
        assert!(str_contains(&newline, "b$", true).is_err());
        let lowered = series.str().unwrap().to_lowercase();
        assert_eq!(
            lowered.into_iter().collect::<Vec<_>>(),
            [
                Some("pathogenic"),
                None,
                Some("likely_pathogenic"),
                Some("benign"),
                Some("")
            ]
        );
    }

    #[test]
    fn where_mask_and_casts_follow_pandas() {
        let mask = BooleanChunked::new("m".into(), &[true, false, true]);
        let then = Series::new("a".into(), &[1.0, 2.0, 3.0]);
        let otherwise = Series::new("b".into(), &[9.0, 8.0, 7.0]);
        let picked = where_mask(&mask, &then, &otherwise).unwrap();
        assert_eq!(
            picked
                .f64()
                .unwrap()
                .into_no_null_iter()
                .collect::<Vec<_>>(),
            [1.0, 8.0, 3.0]
        );
        let with_missing = BooleanChunked::new("m".into(), &[Some(true), None, Some(true)]);
        assert!(where_mask(&with_missing, &then, &otherwise).is_err());
        let floats = to_series(
            &PandasColumn::Float(vec![1.0, f64::NAN, 0.5, -0.0, 0.0]),
            "x",
        )
        .unwrap();
        assert_eq!(
            mask_values(&float_astype_bool(&floats).unwrap()).unwrap(),
            [true, true, true, false, false]
        );
        assert!(float_astype_int(&floats).is_err());
        let fractions = to_series(&PandasColumn::Float(vec![1.9, -1.9]), "f").unwrap();
        assert_eq!(float_astype_int(&fractions).unwrap(), [1, -1]);
        assert!(
            float_astype_int(&to_series(&PandasColumn::Float(vec![1e30]), "f").unwrap()).is_err()
        );
    }

    #[test]
    fn columns_survive_the_trip_to_polars_and_back() {
        let columns = [
            (PandasColumn::Int(vec![1, -2]), PandasDtype::Int64),
            (
                PandasColumn::Float(vec![1.5, f64::NAN]),
                PandasDtype::Float64,
            ),
            (PandasColumn::Bool(vec![true, false]), PandasDtype::Bool),
            (
                PandasColumn::BoolOrMissing(vec![Some(true), None]),
                PandasDtype::ObjectBool,
            ),
            (text(&[Some("a"), None]), PandasDtype::ObjectText),
            (
                PandasColumn::NullableInt(vec![Some(3), None]),
                PandasDtype::NullableInt64,
            ),
            (
                PandasColumn::NullableBool(vec![None, Some(false)]),
                PandasDtype::NullableBoolean,
            ),
        ];
        for (column, dtype) in columns {
            let back = from_series(&to_series(&column, "c").unwrap(), dtype).unwrap();
            match (&column, &back) {
                (PandasColumn::Float(a), PandasColumn::Float(b)) => {
                    assert_eq!(a[0], b[0]);
                    assert!(b[1].is_nan());
                }
                _ => assert_eq!(column, back),
            }
        }
        assert!(to_series(&PandasColumn::Object(vec![ObjectValue::Int(1)]), "o").is_err());
        let with_null = to_series(&PandasColumn::NullableInt(vec![Some(1), None]), "c").unwrap();
        assert!(from_series(&with_null, PandasDtype::Int64).is_err());
    }

    #[test]
    fn kleene_logic_matches_the_nullable_boolean_dtype() {
        // pandas "boolean": & gives [NA, False, False, NA], | gives [True, NA, NA, True].
        let left = to_series(
            &PandasColumn::NullableBool(vec![Some(true), None, Some(false), None]),
            "a",
        )
        .unwrap();
        let right = to_series(
            &PandasColumn::NullableBool(vec![None, Some(false), None, Some(true)]),
            "b",
        )
        .unwrap();
        let (left, right) = (left.bool().unwrap(), right.bool().unwrap());
        assert_eq!(
            values(&(left & right)),
            [None, Some(false), Some(false), None]
        );
        assert_eq!(
            values(&(left | right)),
            [Some(true), None, None, Some(true)]
        );
    }

    #[test]
    fn duplicates_groups_and_sorts_follow_pandas() {
        // pandas drop_duplicates keeps rows [0, 2, 4]: -0.0 equals 0.0 and None equals NaN.
        let frame = table(vec![
            (
                "a",
                PandasColumn::Float(vec![0.0, -0.0, f64::NAN, f64::NAN, 1.0]),
            ),
            (
                "b",
                PandasColumn::Object(vec![
                    ObjectValue::Text("x".into()),
                    ObjectValue::Text("x".into()),
                    ObjectValue::None,
                    ObjectValue::Nan,
                    ObjectValue::Text("x".into()),
                ]),
            ),
        ]);
        assert_eq!(drop_duplicates_first(&frame, None).unwrap(), [0, 2, 4]);
        assert_eq!(
            drop_duplicates_first(&frame, Some(&["a"])).unwrap(),
            [0, 2, 4]
        );
        // groupby("g", sort=False): b, a, c; with dropna=False the missing key is third.
        let keys = table(vec![(
            "g",
            text(&[Some("b"), Some("a"), None, Some("b"), Some("c"), None]),
        )]);
        assert_eq!(
            first_appearance_groups(&keys, &["g"], true).unwrap(),
            vec![vec![0, 3], vec![1], vec![4]]
        );
        assert_eq!(
            first_appearance_groups(&keys, &["g"], false).unwrap(),
            vec![vec![0, 3], vec![1], vec![2, 5], vec![4]]
        );
        // sort_values(["s", "t"], ascending=[False, True]) gives [2, 0, 3, 5, 6, 1, 4].
        let sort_frame = table(vec![
            (
                "s",
                PandasColumn::Float(vec![1.0, f64::NAN, 2.0, 1.0, f64::NAN, -0.0, 0.0]),
            ),
            ("t", PandasColumn::Int(vec![0; 7])),
        ]);
        assert_eq!(
            stable_sort(&sort_frame, &["s", "t"], &[true, false]).unwrap(),
            [2, 0, 3, 5, 6, 1, 4]
        );
        assert_eq!(
            stable_sort(&sort_frame, &["s", "t"], &[false, false]).unwrap(),
            [5, 6, 0, 3, 2, 1, 4]
        );
    }
}
