//! `scripts/acmg_pp3_bp4_bp7_insilico.py` (139 lines): PP3, BP4 and BP7 from the in-silico
//! predictors.
//!
//! `pp3_bp4_criteria` is `PP3_BP4_criteria`, `bp7_criteria` is `BP7_criteria`. Both return
//! int64 strengths (0 or 1), as the Python arrays are.
//!
//! The submodule `pandas_masks` holds what the four T6 criterion modules share: the boolean
//! masks the Python code builds with `&`, `|`, `~`, `np.logical_not` and `fillna`, with pandas'
//! rules for an `object` mask that holds NaN; `str.contains` on a text column; the numeric
//! columns read through `pd.to_numeric(errors="coerce")`; ClinVar's review-status map.
//!
//! Column types the Python code would treat in another way than the real tables need (an
//! `object` column mixing kinds of values, the nullable `Int64`/`boolean` columns, a text column
//! where a number is expected) are refused with an error naming the column.
use super::types::{HarnessCall, StepOutput, not_converted};
use crate::cli::{Result, invalid};
use crate::pandas_table::{Column, Table};
use pandas_masks::{
    Mask, bool_input, bool_input_fillna_false, high_confidence_review, numeric_column, strengths,
    text_contains,
};

/// The default `high_confidence_status` of PP3_BP4_criteria (and PS3_BS3_criteria):
/// ClinVar review statuses of two stars and more. Every value is at least 2, so
/// `map(...) >= 2` is "the status is one of these keys".
pub(super) const HIGH_CONFIDENCE_STATUS: [&str; 3] = [
    "practice_guideline",
    "reviewed_by_expert_panel",
    "criteria_provided,_multiple_submitters,_no_conflicts",
];

/// The SpliceAI cutoff PP3_BP4_criteria uses when the caller gives none.
pub(crate) const DEFAULT_SPLICEAI_NO_IMPACT_CUTOFF: f64 = 0.1;

/// `(?:^|[&,;|])(?:HC|OS)(?:[&,;|]|$)` searched in a text: a field of `HC` or `OS`
/// between the separators `& , ; |`. Python's `$` also matches before a final line
/// break, so the last field may be followed by one `\n`.
fn loftee_high_confidence(text: &str) -> bool {
    let fields: Vec<&str> = text.split(['&', ',', ';', '|']).collect();
    let last = fields.len() - 1;
    fields.iter().enumerate().any(|(index, field)| {
        matches!(*field, "HC" | "OS") || (index == last && matches!(*field, "HC\n" | "OS\n"))
    })
}

/// `df.get("LoF", pd.Series("", index=df.index)).fillna("").astype(str)`, then the LOFTEE
/// pattern: true where a LOFTEE high-confidence (or `OS`) call is listed.
fn loftee_mask(table: &Table) -> Result<Vec<bool>> {
    let Some(index) = table.column_index("LoF") else {
        return Ok(vec![false; table.rows]);
    };
    let column = &table.columns[index];
    let text_of = pandas_masks::python_str_filled(column, "LoF", "")?;
    Ok(text_of
        .iter()
        .map(|text| loftee_high_confidence(text))
        .collect())
}

/// `PP3_BP4_criteria(df, pvs1_criteria, spliceai_no_impact_cutoff=...)`: `(PP3, BP4)`.
pub(crate) fn pp3_bp4_criteria(
    table: &Table,
    pvs1_criteria: Option<&[i64]>,
    spliceai_no_impact_cutoff: f64,
) -> Result<StepOutput> {
    let rows = table.rows;
    if let Some(pvs1) = pvs1_criteria
        && pvs1.len() != rows
    {
        return Err(invalid(format!(
            "pvs1_criteria has {} values for {rows} rows",
            pvs1.len()
        )));
    }
    let mut flags = text_contains(
        table,
        "Consequence",
        None,
        &[
            "missense_variant",
            "stop",
            "frameshift",
            "splice",
            "5_prime_UTR_variant",
        ],
    )?
    .into_iter();
    let mut next_flag = || flags.next().expect("five patterns");
    let (missense, stop, frameshift, splice, five_utr) = (
        next_flag(),
        next_flag(),
        next_flag(),
        next_flag(),
        next_flag(),
    );
    // `contains(kind) & np.logical_not(contains('stop')) & np.logical_not(contains('frameshift'))`.
    let missense_variant = missense
        .and(&stop.logical_not())
        .and(&frameshift.logical_not());
    let splice_variant = splice
        .and(&stop.logical_not())
        .and(&frameshift.logical_not());
    let five_utr_variant = five_utr.fillna(false);
    let missense_variant = missense_variant.fillna(false);
    let splice_variant = splice_variant.fillna(false);
    let loftee_splice_lof = splice_variant.and(&Mask::Plain(loftee_mask(table)?));

    let primateai = fill_missing(numeric_column(table, "PrimateAI")?, 0.0);
    let am_pathogenicity = fill_missing(numeric_column(table, "am_pathogenicity")?, 0.0);
    let cadd_phred = numeric_column(table, "CADD_phred")?;
    let cadd_reg_phred = numeric_column(table, "CADD_reg_phred")?;
    let cadd_phred_high = Mask::Plain(cadd_phred.iter().map(|value| *value >= 20.0).collect());
    let cadd_bp4_phred_cutoff = 15.0;
    let cadd_phred_low = Mask::Plain(
        cadd_phred
            .iter()
            .map(|value| !value.is_nan() && *value < cadd_bp4_phred_cutoff)
            .collect(),
    );
    let cadd_reg_phred_not_high = Mask::Plain(
        cadd_reg_phred
            .iter()
            .map(|value| !(*value >= cadd_bp4_phred_cutoff))
            .collect(),
    );
    let splice_computational_lof =
        bool_input_fillna_false(table, "splicing_lof")?.or(&loftee_splice_lof);
    let spliceai_columns = [
        "SpliceAI_pred_DS_AG",
        "SpliceAI_pred_DS_AL",
        "SpliceAI_pred_DS_DG",
        "SpliceAI_pred_DS_DL",
    ];
    let mut spliceai_max = vec![f64::NAN; rows];
    for name in spliceai_columns {
        for (best, value) in spliceai_max.iter_mut().zip(numeric_column(table, name)?) {
            // `max(axis=1)` skips missing values; all missing stays NaN.
            if !value.is_nan() && (best.is_nan() || value > *best) {
                *best = value;
            }
        }
    }
    let spliceai_no_impact = Mask::Plain(
        spliceai_max
            .iter()
            .map(|value| !value.is_nan() && *value < spliceai_no_impact_cutoff)
            .collect(),
    );

    let not_splice = splice_variant.logical_not();
    let not_missense = missense_variant.logical_not();
    let not_five_utr = five_utr_variant.logical_not();
    let primateai_high = Mask::Plain(primateai.iter().map(|value| *value > 0.8).collect());
    let am_class_pathogenic =
        text_contains(table, "am_class", Some(""), &["pathogenic"])?.remove(0);
    // Columns in the order Python reads them (a missing one raises the same KeyError).
    let vep_consq_lof = bool_input(table, "vep_consq_lof")?;
    let utr_lof = bool_input(table, "5UTR_lof")?;
    let mut clnsig = text_contains(table, "CLNSIG", Some(""), &["enign", "athogenic"])?;
    let high_confidence = Mask::Plain(high_confidence_review(table, "CLNREVSTAT")?);
    let pp3 = primateai_high
        .and(&missense_variant)
        .or(&cadd_phred_high
            .and(&not_splice)
            .and(&not_missense)
            .and(&not_five_utr))
        .or(&am_class_pathogenic.and(&missense_variant))
        .or(&vep_consq_lof
            .and(&not_splice)
            .and(&not_missense)
            .and(&not_five_utr))
        .or(&splice_computational_lof
            .or(&cadd_phred_high)
            .and(&splice_variant)
            .or(&utr_lof.and(&five_utr_variant)));
    let clinvar_patho = clnsig.pop().expect("two patterns").and(&high_confidence);
    let clinvar_benign = clnsig.pop().expect("two patterns").and(&high_confidence);
    let mut pp3 = pp3.and(&clinvar_benign.invert()?);
    if let Some(pvs1) = pvs1_criteria {
        pp3 = pp3.and(&Mask::Plain(pvs1.iter().map(|value| *value < 3).collect()));
    }

    let missense_benign = Mask::Plain(primateai.iter().map(|value| *value < 0.8).collect())
        .and(&Mask::Plain(
            am_pathogenicity
                .iter()
                .map(|value| *value < 0.564)
                .collect(),
        ))
        .and(&missense_variant);
    let splice_benign = splice_computational_lof
        .logical_not()
        .and(&splice_variant)
        .and(&cadd_phred_low.or(&spliceai_no_impact));
    let utr_benign = bool_input_fillna_false(table, "5UTR_lof")?
        .logical_not()
        .and(&five_utr_variant);
    let other_benign = bool_input_fillna_false(table, "vep_consq_lof")?
        .logical_not()
        .and(&cadd_phred_low)
        .and(&cadd_reg_phred_not_high)
        .and(&not_splice)
        .and(&not_missense)
        .and(&not_five_utr);
    let mut bp4 = missense_benign
        .or(&splice_benign)
        .or(&utr_benign)
        .or(&other_benign)
        .and(&clinvar_patho.invert()?)
        .and(&pp3.invert()?);
    if let Some(pvs1) = pvs1_criteria {
        bp4 = bp4.and(&Mask::Plain(pvs1.iter().map(|value| *value < 3).collect()));
    }
    Ok(StepOutput {
        criteria: vec![
            ("PP3", Column::Int(strengths(&pp3, 1)?)),
            ("BP4", Column::Int(strengths(&bp4, 1)?)),
        ],
        ..StepOutput::default()
    })
}

fn fill_missing(mut values: Vec<f64>, fill: f64) -> Vec<f64> {
    for value in &mut values {
        if value.is_nan() {
            *value = fill;
        }
    }
    values
}

/// `BP7_criteria(df)`: synonymous, no splicing-altering consequence, GERP++ at most 5.
pub(crate) fn bp7_criteria(table: &Table) -> Result<StepOutput> {
    let consequence = table.column("Consequence")?;
    let length_changing = table.column("splicing_len_changing")?;
    let splicing_lof = table.column("splicing_lof")?;
    let conservation = table.column("Conservation")?;
    if !matches!(conservation, Column::Float(_) | Column::Int(_)) {
        return Err(invalid(format!(
            "Conservation is a pandas {} column; `<= 5` is reproduced for number columns only",
            conservation.dtype_name()
        )));
    }
    let checked = |column: &Column, name: &str| -> Result<()> {
        if matches!(column, Column::NullableBool(_) | Column::NullableInt(_)) {
            return Err(invalid(format!(
                "{name} is a nullable pandas column; `== False` there gives pd.NA, which BP7 does not expect"
            )));
        }
        Ok(())
    };
    checked(length_changing, "splicing_len_changing")?;
    checked(splicing_lof, "splicing_lof")?;
    let values: Vec<i64> = (0..table.rows)
        .map(|row| {
            let synonymous =
                pandas_masks::cell_equals_text(consequence.cell(row), "synonymous_variant");
            let no_splicing_altering = pandas_masks::cell_equals_false(length_changing.cell(row))
                || pandas_masks::cell_is_missing(splicing_lof.cell(row))
                || pandas_masks::cell_equals_false(splicing_lof.cell(row))
                || pandas_masks::cell_is_missing(length_changing.cell(row));
            let not_conserved = match conservation.cell(row) {
                crate::pandas_table::Cell::Float(value) => value <= 5.0,
                crate::pandas_table::Cell::Int(value) => (value as f64) <= 5.0,
                _ => false,
            };
            i64::from(synonymous && no_splicing_altering && not_conserved)
        })
        .collect();
    Ok(StepOutput {
        criteria: vec![("BP7", Column::Int(values))],
        ..StepOutput::default()
    })
}

/// `priva acmg-step acmg_pp3_bp4_bp7_insilico <function>`: runs one converted function on a frame.
///
/// - `PP3_BP4_criteria`, options `pvs1_criteria` (list of ints or null) and
///   `spliceai_no_impact_cutoff` (default 0.1): a frame of `PP3`, `BP4` (int64);
/// - `BP7_criteria`: a frame of `BP7` (int64).
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let output = match function {
        "PP3_BP4_criteria" => {
            let pvs1 = harness_options::int_list(&call.options, "pvs1_criteria")?;
            let cutoff = harness_options::float(&call.options, "spliceai_no_impact_cutoff")?
                .unwrap_or(DEFAULT_SPLICEAI_NO_IMPACT_CUTOFF);
            pp3_bp4_criteria(&call.table, pvs1.as_deref(), cutoff)?
        }
        "BP7_criteria" => bp7_criteria(&call.table)?,
        _ => return Err(not_converted("acmg_pp3_bp4_bp7_insilico", function)),
    };
    Ok(harness_options::criteria_frame(&output, call.table.rows))
}

/// The `--options` JSON of the T6 harness entries.
pub(super) mod harness_options {
    use super::super::types::StepOutput;
    use crate::cli::{Result, invalid};
    use crate::json_text::Json;
    use crate::pandas_table::{Column, Table};

    fn entry<'a>(options: &'a Json, key: &str) -> Option<&'a Json> {
        options
            .get(key)
            .filter(|value| !matches!(value, Json::Null))
    }

    pub(in super::super) fn int_list(options: &Json, key: &str) -> Result<Option<Vec<i64>>> {
        let Some(value) = entry(options, key) else {
            return Ok(None);
        };
        let Json::Array(items) = value else {
            return Err(invalid(format!("option {key} must be a list of integers")));
        };
        items
            .iter()
            .map(|item| match item {
                Json::Int(number) => i64::try_from(*number)
                    .map_err(|_| invalid(format!("option {key}: {number} beyond 64 bits"))),
                _ => Err(invalid(format!("option {key} must be a list of integers"))),
            })
            .collect::<Result<Vec<i64>>>()
            .map(Some)
    }

    pub(in super::super) fn bool_list(options: &Json, key: &str) -> Result<Option<Vec<bool>>> {
        let Some(value) = entry(options, key) else {
            return Ok(None);
        };
        let Json::Array(items) = value else {
            return Err(invalid(format!("option {key} must be a list of booleans")));
        };
        items
            .iter()
            .map(|item| match item {
                Json::Bool(flag) => Ok(*flag),
                _ => Err(invalid(format!("option {key} must be a list of booleans"))),
            })
            .collect::<Result<Vec<bool>>>()
            .map(Some)
    }

    pub(in super::super) fn float(options: &Json, key: &str) -> Result<Option<f64>> {
        match entry(options, key) {
            None => Ok(None),
            Some(Json::Float(value)) => Ok(Some(*value)),
            Some(Json::Int(value)) => Ok(Some(*value as f64)),
            Some(_) => Err(invalid(format!("option {key} must be a number"))),
        }
    }

    pub(in super::super) fn text<'a>(options: &'a Json, key: &str) -> Result<Option<&'a str>> {
        match entry(options, key) {
            None => Ok(None),
            Some(Json::Str(value)) => Ok(Some(value.as_str())),
            Some(_) => Err(invalid(format!("option {key} must be a text"))),
        }
    }

    pub(in super::super) fn required_text<'a>(options: &'a Json, key: &str) -> Result<&'a str> {
        text(options, key)?.ok_or_else(|| invalid(format!("option {key} is required")))
    }

    /// The criteria of `output` as a frame, one column per criterion, in order.
    pub(in super::super) fn criteria_frame(output: &StepOutput, rows: usize) -> Table {
        Table {
            names: output
                .criteria
                .iter()
                .map(|(name, _)| (*name).to_owned())
                .collect(),
            columns: output
                .criteria
                .iter()
                .map(|(_, values)| values.clone())
                .collect(),
            rows,
        }
    }

    /// A frame of named columns of equal length.
    pub(in super::super) fn frame(columns: Vec<(String, Column)>, rows: usize) -> Table {
        let (names, columns) = columns.into_iter().unzip();
        Table {
            names,
            columns,
            rows,
        }
    }
}

/// The pandas boolean algebra and column readers the T6 criteria share.
pub(super) mod pandas_masks {
    use crate::cli::{Result, invalid};
    use crate::pandas_polars;
    use crate::pandas_table::{Cell, Column, ObjectValue, Table};
    use crate::python_text::float_repr;

    /// A boolean pandas Series as the criteria build it: a plain `bool` Series, or an
    /// `object` Series of True, False and NaN (`None`), which `str.contains` makes on a
    /// column with missing values and which the boolean input columns are when they hold
    /// a missing value.
    #[derive(Clone, Debug, PartialEq)]
    pub(in super::super) enum Mask {
        Plain(Vec<bool>),
        Object(Vec<Option<bool>>),
    }

    impl Mask {
        /// What `str.contains` returns: `bool` without missing values, else `object`.
        pub(in super::super) fn from_options(values: Vec<Option<bool>>) -> Mask {
            if values.iter().all(Option::is_some) {
                Mask::Plain(
                    values
                        .into_iter()
                        .map(|value| value.unwrap_or(false))
                        .collect(),
                )
            } else {
                Mask::Object(values)
            }
        }

        pub(in super::super) fn len(&self) -> usize {
            match self {
                Mask::Plain(values) => values.len(),
                Mask::Object(values) => values.len(),
            }
        }

        fn value(&self, row: usize) -> Option<bool> {
            match self {
                Mask::Plain(values) => Some(values[row]),
                Mask::Object(values) => values[row],
            }
        }

        /// pandas' `logical_op` for `self <op> right` (or `ndarray <op> self`, which pandas
        /// turns into the same call with `self` first): the right operand's NaN become
        /// False; where `self` holds NaN, Python's `nan & x` raises TypeError, pandas
        /// keeps NaN there and finally turns it into False. The result is a plain mask.
        fn combine(&self, right: &Mask, op: fn(bool, bool) -> bool) -> Mask {
            assert_eq!(self.len(), right.len(), "masks of one table");
            Mask::Plain(
                (0..self.len())
                    .map(|row| match self.value(row) {
                        None => false,
                        Some(left) => op(left, right.value(row).unwrap_or(false)),
                    })
                    .collect(),
            )
        }

        /// `self & right`.
        pub(in super::super) fn and(&self, right: &Mask) -> Mask {
            self.combine(right, |left, right| left && right)
        }

        /// `self | right` (NaN in `self` gives False, also against True).
        pub(in super::super) fn or(&self, right: &Mask) -> Mask {
            self.combine(right, |left, right| left || right)
        }

        /// `np.logical_not(self)`: `not value`, so NaN (truthy) gives False; an object
        /// Series stays an object Series.
        pub(in super::super) fn logical_not(&self) -> Mask {
            match self {
                Mask::Plain(values) => Mask::Plain(values.iter().map(|value| !value).collect()),
                Mask::Object(values) => Mask::Object(
                    values
                        .iter()
                        .map(|value| Some(!value.unwrap_or(true)))
                        .collect(),
                ),
            }
        }

        /// `~self`. On an object Series Python's `~True` is `-2` (and `~nan` raises), so
        /// only a plain mask is inverted here.
        pub(in super::super) fn invert(&self) -> Result<Mask> {
            match self {
                Mask::Plain(values) => Ok(Mask::Plain(values.iter().map(|value| !value).collect())),
                Mask::Object(_) => Err(invalid(
                    "`~` of an object boolean Series gives -2/-1 (or raises on NaN) in pandas; \
                     not reproduced",
                )),
            }
        }

        /// `self.fillna(value)`: an object Series with a missing value becomes `bool`
        /// (pandas 2.2 downcasts after filling); one without stays `object`.
        pub(in super::super) fn fillna(&self, value: bool) -> Mask {
            match self {
                Mask::Object(values) if values.iter().any(Option::is_none) => {
                    Mask::Plain(values.iter().map(|item| item.unwrap_or(value)).collect())
                }
                other => other.clone(),
            }
        }

        /// The values of a mask numpy can index with (`array[mask] = k`).
        pub(in super::super) fn plain(&self) -> Result<&[bool]> {
            match self {
                Mask::Plain(values) => Ok(values),
                Mask::Object(_) => Err(invalid(
                    "numpy cannot index with an object boolean Series (IndexError); not reproduced",
                )),
            }
        }
    }

    /// `array = np.zeros(len(df), dtype=int); array[mask] = strength`.
    pub(in super::super) fn strengths(mask: &Mask, strength: i64) -> Result<Vec<i64>> {
        Ok(mask
            .plain()?
            .iter()
            .map(|flag| if *flag { strength } else { 0 })
            .collect())
    }

    fn attribute_error() -> Box<dyn std::error::Error> {
        Box::new(super::super::types::AcmgError::new(
            "AttributeError",
            "Can only use .str accessor with string values!",
        ))
    }

    /// The text values of `df[name]` (optionally after `.fillna(fill)`) as the `.str`
    /// accessor sees them; `None` is a missing value. A number column raises pandas'
    /// AttributeError, except a float column of only missing values that `fillna` turns
    /// into text.
    pub(in super::super) fn text_values(
        table: &Table,
        name: &str,
        fill: Option<&str>,
    ) -> Result<Vec<Option<String>>> {
        let column = table.column(name)?;
        match (column, fill) {
            (Column::Text(values), None) => Ok(values.clone()),
            (Column::Text(values), Some(fill)) => Ok(values
                .iter()
                .map(|value| Some(value.clone().unwrap_or_else(|| fill.to_owned())))
                .collect()),
            (Column::Float(values), Some(fill)) if values.iter().any(|value| value.is_nan()) => {
                if values.iter().all(|value| value.is_nan()) {
                    Ok(vec![Some(fill.to_owned()); values.len()])
                } else {
                    Err(invalid(format!(
                        "{name} is a float column with numbers and missing values; \
                         fillna({fill:?}) then .str is not reproduced"
                    )))
                }
            }
            (Column::Int(_) | Column::Float(_) | Column::Bool(_), _) => Err(attribute_error()),
            (other, _) => Err(invalid(format!(
                "{name} is a pandas {} column of other values than text; .str is not reproduced",
                other.dtype_name()
            ))),
        }
    }

    /// `df[name]` (optionally `.fillna(fill)`) `.str.contains(pattern)` for each literal
    /// pattern (the step's patterns are plain words, the same as a regular expression).
    pub(in super::super) fn text_contains(
        table: &Table,
        name: &str,
        fill: Option<&str>,
        patterns: &[&str],
    ) -> Result<Vec<Mask>> {
        let values = text_values(table, name, fill)?;
        let series = pandas_polars::to_series(&Column::Text(values), name)?;
        patterns
            .iter()
            .map(|pattern| {
                assert!(
                    pattern
                        .bytes()
                        .all(|byte| byte.is_ascii_alphanumeric() || byte == b'_'),
                    "only literal words are searched here"
                );
                let found = pandas_polars::str_contains(&series, pattern, false)?;
                Ok(Mask::from_options(found.into_iter().collect()))
            })
            .collect()
    }

    /// `str(value)` of each cell after `.fillna(fill)` (`astype(str)`).
    pub(in super::super) fn python_str_filled(
        column: &Column,
        name: &str,
        fill: &str,
    ) -> Result<Vec<String>> {
        if matches!(column, Column::NullableInt(_) | Column::NullableBool(_)) {
            return Err(invalid(format!(
                "{name} is a nullable pandas column; fillna({fill:?}) is not reproduced there"
            )));
        }
        Ok((0..column.len())
            .map(|row| match column.cell(row) {
                Cell::Missing => fill.to_owned(),
                cell => python_str(cell),
            })
            .collect())
    }

    /// Python's `str()` of a value that is not missing.
    pub(in super::super) fn python_str(cell: Cell<'_>) -> String {
        match cell {
            Cell::Missing => "nan".to_owned(),
            Cell::Text(text) => text.to_owned(),
            Cell::Int(number) => number.to_string(),
            Cell::Float(number) => float_repr(number),
            Cell::Bool(flag) => (if flag { "True" } else { "False" }).to_owned(),
        }
    }

    /// `pd.to_numeric(df[name], errors="coerce")` for a number column, as doubles with NaN
    /// for missing values (numpy compares an int64 column with a float as doubles).
    pub(in super::super) fn numeric_column(table: &Table, name: &str) -> Result<Vec<f64>> {
        match table.column(name)? {
            Column::Float(values) => Ok(values.clone()),
            Column::Int(values) => Ok(values.iter().map(|value| *value as f64).collect()),
            other => Err(invalid(format!(
                "{name} is a pandas {} column; pd.to_numeric(errors=\"coerce\") is reproduced \
                 for number columns only",
                other.dtype_name()
            ))),
        }
    }

    /// `df[name].map(high_confidence_status) >= 2`: the review status is one of the
    /// three two-star-or-better texts.
    pub(in super::super) fn high_confidence_review(table: &Table, name: &str) -> Result<Vec<bool>> {
        let column = table.column(name)?;
        if matches!(column, Column::NullableInt(_) | Column::NullableBool(_)) {
            return Err(invalid(format!(
                "{name} is a nullable pandas column; map(...) >= 2 is not reproduced there"
            )));
        }
        Ok((0..table.rows)
            .map(|row| match column.cell(row) {
                Cell::Text(text) => super::HIGH_CONFIDENCE_STATUS.contains(&text),
                _ => false,
            })
            .collect())
    }

    /// A boolean input column used as it is (`df[name] & ...`): `bool`, or an `object`
    /// column of True, False and NaN. A column of only missing values (`object` or
    /// `float64`) acts as NaN everywhere: on the left of `&`/`|` pandas gives False, on
    /// the right NaN counts as False, and `fillna(False)` makes all False (an `object`
    /// Series of False, whose `~` values `-1` count as True where the step uses them).
    pub(in super::super) fn bool_input(table: &Table, name: &str) -> Result<Mask> {
        match table.column(name)? {
            Column::Bool(values) => Ok(Mask::Plain(values.clone())),
            Column::BoolOrMissing(values) => Ok(Mask::Object(values.clone())),
            Column::Text(values) if values.iter().all(Option::is_none) => {
                Ok(Mask::Object(vec![None; values.len()]))
            }
            Column::Float(values) if values.iter().all(|value| value.is_nan()) => {
                Ok(Mask::Object(vec![None; values.len()]))
            }
            Column::Object(values) => Ok(Mask::Object(object_booleans(values, name)?)),
            other => Err(invalid(format!(
                "{name} is a pandas {} column; the boolean algebra on it is reproduced for \
                 boolean columns only",
                other.dtype_name()
            ))),
        }
    }

    /// `df[name].fillna(False)` of a boolean input column.
    pub(in super::super) fn bool_input_fillna_false(table: &Table, name: &str) -> Result<Mask> {
        Ok(bool_input(table, name)?.fillna(false))
    }

    fn object_booleans(values: &[ObjectValue], name: &str) -> Result<Vec<Option<bool>>> {
        values
            .iter()
            .map(|value| match value {
                ObjectValue::Bool(flag) => Ok(Some(*flag)),
                ObjectValue::Nan => Ok(None),
                ObjectValue::Float(number) if number.is_nan() => Ok(None),
                _ => Err(invalid(format!(
                    "{name} mixes booleans with other values; not reproduced"
                ))),
            })
            .collect()
    }

    /// `pd.isna(value)`.
    pub(in super::super) fn cell_is_missing(cell: Cell<'_>) -> bool {
        cell == Cell::Missing
    }

    /// `value == False` in Python (numpy for number columns): `False`, `0` and `0.0`.
    pub(in super::super) fn cell_equals_false(cell: Cell<'_>) -> bool {
        match cell {
            Cell::Bool(flag) => !flag,
            Cell::Int(number) => number == 0,
            Cell::Float(number) => number == 0.0,
            Cell::Text(_) | Cell::Missing => false,
        }
    }

    /// `value == text`.
    pub(in super::super) fn cell_equals_text(cell: Cell<'_>, text: &str) -> bool {
        matches!(cell, Cell::Text(value) if value == text)
    }
}

#[cfg(test)]
mod tests {
    use super::pandas_masks::Mask;
    use super::*;

    fn text(values: &[Option<&str>]) -> Column {
        Column::Text(
            values
                .iter()
                .map(|value| value.map(str::to_owned))
                .collect(),
        )
    }

    #[test]
    fn loftee_pattern_follows_python_re() {
        assert!(loftee_high_confidence("HC"));
        assert!(loftee_high_confidence("LC&HC"));
        assert!(loftee_high_confidence("OS,LC"));
        assert!(loftee_high_confidence("LC|HC\n"));
        assert!(!loftee_high_confidence("HC\n&LC"));
        assert!(!loftee_high_confidence("HCX"));
        assert!(!loftee_high_confidence(""));
        assert!(!loftee_high_confidence("1.0"));
    }

    #[test]
    fn object_masks_follow_pandas_logical_op() {
        let object = Mask::Object(vec![Some(true), None, Some(false), None]);
        let plain = Mask::Plain(vec![true, true, true, false]);
        // NaN on the left gives False, also for `|`; on the right it counts as False.
        assert_eq!(
            object.or(&plain),
            Mask::Plain(vec![true, false, true, false])
        );
        assert_eq!(
            plain.or(&object),
            Mask::Plain(vec![true, true, true, false])
        );
        assert_eq!(
            object.and(&plain),
            Mask::Plain(vec![true, false, false, false])
        );
        assert_eq!(
            object.logical_not(),
            Mask::Object(vec![Some(false), Some(false), Some(true), Some(false)])
        );
        assert_eq!(
            object.fillna(false),
            Mask::Plain(vec![true, false, false, false])
        );
        let no_missing = Mask::Object(vec![Some(true), Some(false)]);
        assert_eq!(no_missing.fillna(false), no_missing);
        assert!(no_missing.invert().is_err());
    }

    #[test]
    fn bp7_reads_python_equality() {
        let table = Table {
            names: vec![
                "Consequence".into(),
                "splicing_len_changing".into(),
                "splicing_lof".into(),
                "Conservation".into(),
            ],
            columns: vec![
                text(&[Some("synonymous_variant"), Some("synonymous_variant"), None]),
                Column::BoolOrMissing(vec![Some(true), None, Some(false)]),
                Column::Float(vec![1.0, 1.0, 0.0]),
                Column::Float(vec![5.0, f64::NAN, 1.0]),
            ],
            rows: 3,
        };
        let output = bp7_criteria(&table).unwrap();
        assert_eq!(output.criterion("BP7"), Some(&Column::Int(vec![0, 0, 0])));
    }
}
