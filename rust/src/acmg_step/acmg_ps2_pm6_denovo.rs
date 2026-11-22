//! `scripts/acmg_ps2_pm6_denovo.py` (215 lines): PS2 and PM6 from the family genotypes.
//!
//! `identify_fam_members` finds the proband (the family's first row with
//! phenotype 2 after `pd.to_numeric`), and the father and mother named on the
//! proband's row; `determine_denovo_per_row` compares the parents' genotype
//! texts with the substring tests of the script (`"1" in gt`, `"." in gt`);
//! `PS2_PM6_criteria` judges each distinct `chrom_pos_ref_alt` once (its first
//! row) and gives every row of that key the same answer: PS2 = 3, PM6 = 2.
//!
//! Reproduced on purpose: the proband's own genotype is never tested; the
//! proband's sex is read from the pedigree as text (`"1"` is male), unlike the
//! numeric conversion `identify_fam_members` makes; on chrY and chrM (and on
//! chrX for a male proband) only the transmitting parent is read, and a missing
//! genotype there gives neither criterion (the Python function returns `None`
//! there); PM6 needs
//! `gnomAD_joint_AF` equal to 0 after `fillna(0)`.
//!
//! The PED table is the step's (`dtype=str, keep_default_na=False`): every
//! value is text. `pd.to_numeric(errors="coerce")` of the phenotype is
//! reproduced for plain decimal integers and decimal or exponent floats;
//! other text holding a digit (`1_0`, `0x1`) and infinities are refused, as
//! are non-integral phenotypes (pandas raises when casting them to `Int64`).
use super::frame_file::read_frame;
use super::types::{AcmgError, HarnessCall, StepOutput};
use crate::cli::{Result, invalid};
use crate::json_text::Json;
use crate::pandas_table::{Cell, Column, ObjectValue, Table};
use crate::python_text::float_repr;
use std::collections::HashMap;
use std::path::Path;

/// A value after `pd.to_numeric(errors="coerce")`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum PedNumber {
    Missing,
    Int(i64),
    Float(f64),
}

impl PedNumber {
    /// `value == number` as pandas compares a numeric cell (NaN is never equal).
    pub(crate) fn equals(self, number: i64) -> bool {
        match self {
            PedNumber::Missing => false,
            PedNumber::Int(value) => value == number,
            PedNumber::Float(value) => value == number as f64,
        }
    }
}

fn refused(message: impl Into<String>) -> Box<dyn std::error::Error> {
    invalid(format!("refused: {}", message.into()))
}

/// Whether `text` is a decimal number `pd.to_numeric` reads: an optional sign,
/// digits with an optional fraction (or a fraction alone), an optional exponent.
fn decimal_number(text: &str) -> Option<bool> {
    let bytes = text.as_bytes();
    let mut index = usize::from(matches!(bytes.first(), Some(b'+' | b'-')));
    let digits_before = bytes[index..]
        .iter()
        .take_while(|byte| byte.is_ascii_digit())
        .count();
    index += digits_before;
    let mut digits_after = 0;
    let mut integer = true;
    if bytes.get(index) == Some(&b'.') {
        integer = false;
        index += 1;
        digits_after = bytes[index..]
            .iter()
            .take_while(|byte| byte.is_ascii_digit())
            .count();
        index += digits_after;
    }
    if digits_before + digits_after == 0 {
        return None;
    }
    if matches!(bytes.get(index), Some(b'e' | b'E')) {
        integer = false;
        index += 1;
        index += usize::from(matches!(bytes.get(index), Some(b'+' | b'-')));
        let exponent = bytes[index..]
            .iter()
            .take_while(|byte| byte.is_ascii_digit())
            .count();
        if exponent == 0 {
            return None;
        }
        index += exponent;
    }
    (index == bytes.len()).then_some(integer)
}

/// `pd.to_numeric(value, errors="coerce")` of one PED value (probed on pandas
/// 2.2.3: surrounding whitespace is ignored, `+2`, `007`, `2.`, `.5` and `1e0`
/// are numbers, text without digits is missing).
pub(crate) fn ped_to_numeric(cell: Cell<'_>) -> Result<PedNumber> {
    let text = match cell {
        Cell::Missing => return Ok(PedNumber::Missing),
        Cell::Int(value) => return Ok(PedNumber::Int(value)),
        Cell::Float(value) => return Ok(PedNumber::Float(value)),
        Cell::Bool(_) => return Err(refused("a boolean pedigree value")),
        Cell::Text(text) => text,
    };
    let trimmed = text.trim_matches(|c: char| c.is_ascii_whitespace());
    match decimal_number(trimmed) {
        Some(true) => {
            let digits = trimmed.trim_start_matches(['+', '-']);
            if digits.len() > 18 {
                return Err(refused(format!(
                    "pedigree value {text:?} has more digits than this tool converts"
                )));
            }
            let value: i64 = trimmed
                .trim_start_matches('+')
                .parse()
                .map_err(|_| refused(format!("pedigree value {text:?} is not converted")))?;
            Ok(PedNumber::Int(value))
        }
        Some(false) => {
            let significant = trimmed
                .split(['e', 'E'])
                .next()
                .unwrap_or("")
                .bytes()
                .filter(u8::is_ascii_digit)
                .count();
            if significant > 15 {
                return Err(refused(format!(
                    "pedigree value {text:?} has more digits than this tool converts"
                )));
            }
            let value: f64 = trimmed
                .parse()
                .map_err(|_| refused(format!("pedigree value {text:?} is not converted")))?;
            if !value.is_finite() {
                return Err(refused(format!("pedigree value {text:?} is not finite")));
            }
            Ok(PedNumber::Float(value))
        }
        None => {
            let lower = trimmed.to_ascii_lowercase();
            let unsigned = lower.trim_start_matches(['+', '-']);
            if unsigned.starts_with("inf") {
                Err(refused(format!(
                    "pedigree value {text:?} reads as an infinity"
                )))
            } else if text.bytes().any(|byte| byte.is_ascii_digit()) {
                Err(refused(format!(
                    "pedigree value {text:?} holds digits but is not a plain number"
                )))
            } else {
                Ok(PedNumber::Missing)
            }
        }
    }
}

/// `pd.to_numeric(..., errors="coerce").astype("Int64")`: `None` is `pd.NA`.
pub(crate) fn ped_to_nullable_int(cell: Cell<'_>) -> Result<Option<i64>> {
    match ped_to_numeric(cell)? {
        PedNumber::Missing => Ok(None),
        PedNumber::Int(value) => Ok(Some(value)),
        PedNumber::Float(value) if value.is_nan() => Ok(None),
        PedNumber::Float(value) if value.fract() == 0.0 && value.abs() < 9.0e15 => {
            Ok(Some(value as i64))
        }
        PedNumber::Float(value) => Err(refused(format!(
            "phenotype {} is not an integer (pandas cannot cast it to Int64)",
            float_repr(value)
        ))),
    }
}

/// The text of a cell the way Python's `str()` writes it; `None` for a
/// missing value.
pub(crate) fn python_str(cell: Cell<'_>) -> Option<String> {
    match cell {
        Cell::Missing => None,
        Cell::Text(text) => Some(text.to_owned()),
        Cell::Int(value) => Some(value.to_string()),
        Cell::Float(value) => Some(float_repr(value)),
        Cell::Bool(flag) => Some(if flag { "True" } else { "False" }.to_owned()),
    }
}

/// The pedigree rows of one family (`ped_df['#FamilyID'] == fam_name`).
pub(crate) fn family_rows(ped: &Table, fam_name: &str) -> Result<Vec<usize>> {
    let families = ped.column("#FamilyID")?;
    Ok((0..ped.rows)
        .filter(|row| families.cell(*row) == Cell::Text(fam_name))
        .collect())
}

/// What `identify_fam_members` returns: `(proband, 2), (father, father_pheno),
/// (mother, mother_pheno), sib_info`. A phenotype is `None` when Python has
/// `None` (the parent is not in the family) and `Some(None)` for `pd.NA`.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct FamilyMembers {
    pub(crate) proband: String,
    pub(crate) father: Option<String>,
    pub(crate) father_pheno: Option<Option<i64>>,
    pub(crate) mother: Option<String>,
    pub(crate) mother_pheno: Option<Option<i64>>,
    pub(crate) siblings: Vec<(String, Option<i64>)>,
}

/// `str(value)` of a pedigree ID, as `.astype(str)` gives it.
fn id_text(cell: Cell<'_>) -> String {
    python_str(cell).unwrap_or_else(|| "nan".to_owned())
}

/// The parent named on the proband's row, `None` for the script's missing values.
fn parent(cell: Option<Cell<'_>>) -> Option<String> {
    match cell? {
        Cell::Missing => None,
        Cell::Text(text) if ["0", ".", "", "NA", "NaN", "nan", "-9"].contains(&text) => None,
        Cell::Int(0 | -9) => None,
        Cell::Float(value) if value == 0.0 || value == -9.0 => None,
        Cell::Bool(false) => None,
        other => python_str(other),
    }
}

/// `identify_fam_members(ped_df, fam_name)`.
pub(crate) fn identify_fam_members(ped: &Table, fam_name: &str) -> Result<FamilyMembers> {
    let rows = family_rows(ped, fam_name)?;
    if rows.is_empty() {
        return Err(AcmgError::new(
            "ValueError",
            format!("Family {fam_name} not found in pedigree table"),
        )
        .into());
    }
    let phenotype_column = ped.column("Phenotype")?;
    let phenotypes: Vec<Option<i64>> = rows
        .iter()
        .map(|row| ped_to_nullable_int(phenotype_column.cell(*row)))
        .collect::<Result<_>>()?;
    ped.column("Sex")?;
    let ids = ped.column("IndividualID")?;
    let patients: Vec<usize> = (0..rows.len())
        .filter(|position| phenotypes[*position] == Some(2))
        .collect();
    let Some(&first) = patients.first() else {
        return Err(AcmgError::new(
            "ValueError",
            format!("No affected individual (Phenotype==2) found for family {fam_name}"),
        )
        .into());
    };
    let proband = match ids.cell(rows[first]) {
        Cell::Text(text) => text.to_owned(),
        other => {
            return Err(refused(format!(
                "the proband's IndividualID {other:?} is not text"
            )));
        }
    };
    // The first family row of that ID, which may come before the first patient row.
    let proband_row = *rows
        .iter()
        .find(|row| ids.cell(**row) == Cell::Text(proband.as_str()))
        .unwrap_or(&rows[first]);
    let field = |name: &str| {
        ped.column_index(name)
            .map(|column| ped.cell(column, proband_row))
    };
    let father = parent(field("PaternalID"));
    let mother = parent(field("MaternalID"));
    let id_texts: Vec<String> = rows.iter().map(|row| id_text(ids.cell(*row))).collect();
    let phenotype_of = |member: &Option<String>| -> Option<Option<i64>> {
        let member = member.as_deref()?;
        id_texts
            .iter()
            .position(|id| id == member)
            .map(|position| phenotypes[position])
    };
    let father_pheno = phenotype_of(&father);
    let mother_pheno = phenotype_of(&mother);
    let mut siblings: Vec<(String, Option<i64>)> = Vec::new();
    for (position, id) in id_texts.iter().enumerate() {
        let known = [&father, &mother]
            .iter()
            .any(|member| member.as_deref() == Some(id.as_str()))
            || *id == proband;
        if known {
            continue;
        }
        match siblings.iter_mut().find(|(sibling, _)| sibling == id) {
            Some(slot) => slot.1 = phenotypes[position],
            None => siblings.push((id.clone(), phenotypes[position])),
        }
    }
    Ok(FamilyMembers {
        proband,
        father,
        father_pheno,
        mother,
        mother_pheno,
        siblings,
    })
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Denovo {
    Ps2,
    Pm6,
    /// `False`, or `None` where the Python function falls off its branches.
    No,
}

/// A genotype as `determine_denovo_per_row` sees it: the cell's text, `'.'`
/// for a value that is not text, `'.'` for a sample without a column.
fn genotype_text<'a>(table: &'a Table, sample: Option<&str>, row: usize) -> Result<&'a str> {
    let Some(column) = sample.and_then(|sample| table.column_index(sample)) else {
        return Ok(".");
    };
    match &table.columns[column] {
        Column::Object(values) if values.iter().any(|value| *value == ObjectValue::None) => {
            Err(refused(format!(
                "sample column {:?} holds None, on which Python's substring tests fail",
                table.names[column]
            )))
        }
        values => Ok(match values.cell(row) {
            Cell::Text(text) => text,
            _ => ".",
        }),
    }
}

/// `determine_denovo_per_row` for one variant.
fn denovo(chrom: &str, male: bool, father: &str, mother: &str, af_absent: bool) -> Denovo {
    let carries = |gt: &str| gt.contains('1');
    let called = |gt: &str| !gt.contains('.');
    let both_parents = || {
        if carries(father) || carries(mother) {
            Denovo::No
        } else if called(father) && called(mother) {
            Denovo::Ps2
        } else if af_absent && (called(father) || called(mother)) {
            Denovo::Pm6
        } else {
            Denovo::No
        }
    };
    let one_parent = |gt: &str| {
        if carries(gt) {
            Denovo::No
        } else if called(gt) {
            Denovo::Ps2
        } else {
            Denovo::No
        }
    };
    match (male, chrom) {
        (true, "chrX") | (true, "chrM") | (false, "chrM") => one_parent(mother),
        (_, "chrY") => one_parent(father),
        _ => both_parents(),
    }
}

/// The text of a key column for `df['chrom'] + '_' + df['pos'].astype(str) + ...`.
fn key_part(table: &Table, name: &str, allow_int: bool) -> Result<Vec<String>> {
    let column = table.column(name)?;
    (0..table.rows)
        .map(|row| match column.cell(row) {
            Cell::Text(text) if matches!(column, Column::Text(_)) => Ok(text.to_owned()),
            Cell::Int(value) if allow_int && matches!(column, Column::Int(_)) => {
                Ok(value.to_string())
            }
            _ => Err(refused(format!(
                "column {name} of type {} with a missing or non-text value builds no variant key here",
                column.dtype_name()
            ))),
        })
        .collect()
}

/// `PS2_PM6_criteria(df, ped_df, fam_name)`: PS2 (3) and PM6 (2) per row.
pub(crate) fn ps2_pm6_criteria(table: &Table, ped: &Table, fam_name: &str) -> Result<StepOutput> {
    let members = identify_fam_members(ped, fam_name)?;
    let family = family_rows(ped, fam_name)?;
    let samples = [
        Some(members.proband.as_str()),
        members.father.as_deref(),
        members.mother.as_deref(),
    ];
    const OWN: [&str; 6] = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "gnomAD_joint_AF",
        "variant_key",
    ];
    if table.has_column("variant_key") {
        return Err(refused(
            "the table already has a variant_key column, which the Python function overwrites and drops",
        ));
    }
    if samples.iter().flatten().any(|sample| OWN.contains(sample)) {
        return Err(refused(
            "a family member is named like a column the function reads",
        ));
    }
    // The key columns in the order Python evaluates the sum (a KeyError or a
    // TypeError of one column comes before the next column is read).
    let chrom = key_part(table, "chrom", false)?;
    let pos = key_part(table, "pos", true)?;
    let reference = key_part(table, "ref", false)?;
    let alt = key_part(table, "alt", false)?;
    let keys: Vec<String> = (0..table.rows)
        .map(|row| {
            format!(
                "{}_{}_{}_{}",
                chrom[row], pos[row], reference[row], alt[row]
            )
        })
        .collect();
    let frequency = table.column("gnomAD_joint_AF")?;
    if table.rows > 0 && !matches!(frequency, Column::Float(_) | Column::Int(_)) {
        return Err(refused(format!(
            "gnomAD_joint_AF of type {} is not compared with 0 here",
            frequency.dtype_name()
        )));
    }
    let proband_sex = {
        let ids = ped.column("IndividualID")?;
        let sexes = ped.column("Sex")?;
        family
            .iter()
            .find(|row| ids.cell(**row) == Cell::Text(members.proband.as_str()))
            .map(|row| sexes.cell(*row))
            .ok_or_else(|| refused("the proband's pedigree row is not found again"))?
    };
    let male = matches!(
        proband_sex,
        Cell::Text("1") | Cell::Int(1) | Cell::Bool(true)
    ) || proband_sex == Cell::Float(1.0);
    let mut answers: HashMap<&str, Denovo> = HashMap::new();
    let mut results = Vec::with_capacity(table.rows);
    for (row, key) in keys.iter().enumerate() {
        let answer = match answers.get(key.as_str()) {
            Some(answer) => *answer,
            None => {
                let af_absent = match frequency.cell(row) {
                    Cell::Missing => true,
                    Cell::Float(value) => value == 0.0,
                    Cell::Int(value) => value == 0,
                    _ => false,
                };
                let father = genotype_text(table, members.father.as_deref(), row)?;
                let mother = genotype_text(table, members.mother.as_deref(), row)?;
                genotype_text(table, Some(members.proband.as_str()), row)?;
                let answer = denovo(&chrom[row], male, father, mother, af_absent);
                answers.insert(key.as_str(), answer);
                answer
            }
        };
        results.push(answer);
    }
    let ps2 = results
        .iter()
        .map(|answer| if *answer == Denovo::Ps2 { 3 } else { 0 })
        .collect();
    let pm6 = results
        .iter()
        .map(|answer| if *answer == Denovo::Pm6 { 2 } else { 0 })
        .collect();
    Ok(StepOutput {
        criteria: vec![("PS2", Column::Int(ps2)), ("PM6", Column::Int(pm6))],
        ..StepOutput::default()
    })
}

/// A text option of a harness call.
pub(crate) fn harness_text<'a>(call: &'a HarnessCall, name: &str) -> Result<&'a str> {
    call.options
        .get(name)
        .and_then(Json::as_str)
        .ok_or_else(|| invalid(format!("option {name} (text) is required")))
}

/// The pedigree frame a harness call names in `ped_frame`.
pub(crate) fn harness_ped(call: &HarnessCall) -> Result<Table> {
    read_frame(Path::new(harness_text(call, "ped_frame")?))
}

fn optional_text(value: Option<&str>) -> Option<String> {
    value.map(str::to_owned)
}

fn phenotype_text(value: Option<Option<i64>>) -> Option<String> {
    match value {
        None => Some("None".to_owned()),
        Some(None) => Some("<NA>".to_owned()),
        Some(Some(number)) => Some(number.to_string()),
    }
}

/// `priva acmg-step acmg_ps2_pm6_denovo <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    match function {
        "identify_fam_members" => {
            let ped = harness_ped(call)?;
            let members = identify_fam_members(&ped, harness_text(call, "fam_name")?)?;
            let mut roles = vec![
                (
                    Some("proband".to_owned()),
                    Some(members.proband.clone()),
                    Some("2".to_owned()),
                ),
                (
                    Some("father".to_owned()),
                    optional_text(members.father.as_deref()),
                    phenotype_text(members.father_pheno),
                ),
                (
                    Some("mother".to_owned()),
                    optional_text(members.mother.as_deref()),
                    phenotype_text(members.mother_pheno),
                ),
            ];
            for (id, phenotype) in &members.siblings {
                roles.push((
                    Some("sibling".to_owned()),
                    Some(id.clone()),
                    phenotype_text(Some(*phenotype)),
                ));
            }
            Ok(Table {
                names: vec!["role".into(), "id".into(), "phenotype".into()],
                rows: roles.len(),
                columns: vec![
                    Column::Text(roles.iter().map(|role| role.0.clone()).collect()),
                    Column::Object(
                        roles
                            .iter()
                            .map(|role| match &role.1 {
                                Some(id) => ObjectValue::Text(id.clone()),
                                None => ObjectValue::None,
                            })
                            .collect(),
                    ),
                    Column::Text(roles.iter().map(|role| role.2.clone()).collect()),
                ],
            })
        }
        "PS2_PM6_criteria" => {
            let ped = harness_ped(call)?;
            let output = ps2_pm6_criteria(&call.table, &ped, harness_text(call, "fam_name")?)?;
            criteria_frame(&output, call.table.rows)
        }
        _ => Err(super::types::not_converted("acmg_ps2_pm6_denovo", function)),
    }
}

/// The criteria of a step output as a frame, one column per criterion.
pub(crate) fn criteria_frame(output: &StepOutput, rows: usize) -> Result<Table> {
    Ok(Table {
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
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn text_table(columns: &[(&str, &[&str])]) -> Table {
        Table {
            names: columns.iter().map(|(name, _)| (*name).to_owned()).collect(),
            rows: columns[0].1.len(),
            columns: columns
                .iter()
                .map(|(_, values)| {
                    Column::Text(
                        values
                            .iter()
                            .map(|value| Some((*value).to_owned()))
                            .collect(),
                    )
                })
                .collect(),
        }
    }

    fn trio() -> Table {
        text_table(&[
            ("#FamilyID", &["F", "F", "F", "G"]),
            ("IndividualID", &["kid", "dad", "mum", "x"]),
            ("PaternalID", &["dad", "0", "0", "0"]),
            ("MaternalID", &["mum", "0", "0", "0"]),
            ("Sex", &["1", "1", "2", "2"]),
            ("Phenotype", &["2", "1", "1", "2"]),
        ])
    }

    #[test]
    fn to_numeric_follows_the_pandas_probe() {
        let cases = [
            ("", PedNumber::Missing),
            (" 2", PedNumber::Int(2)),
            ("+2", PedNumber::Int(2)),
            ("007", PedNumber::Int(7)),
            ("2.", PedNumber::Float(2.0)),
            (".5", PedNumber::Float(0.5)),
            ("1e0", PedNumber::Float(1.0)),
            ("x", PedNumber::Missing),
            ("NA", PedNumber::Missing),
            ("nan", PedNumber::Missing),
            ("true", PedNumber::Missing),
            ("\u{0662}", PedNumber::Missing),
        ];
        for (text, expected) in cases {
            assert_eq!(
                ped_to_numeric(Cell::Text(text)).unwrap(),
                expected,
                "{text:?}"
            );
        }
        for text in [
            "inf",
            "-Infinity",
            "1_0",
            "0x1",
            "1,0",
            "99999999999999999999",
        ] {
            assert!(ped_to_numeric(Cell::Text(text)).is_err(), "{text:?}");
        }
        assert_eq!(ped_to_nullable_int(Cell::Text("2.0")).unwrap(), Some(2));
        assert!(ped_to_nullable_int(Cell::Text("2.5")).is_err());
    }

    #[test]
    fn family_members_come_from_the_probands_row() {
        let members = identify_fam_members(&trio(), "F").unwrap();
        assert_eq!(members.proband, "kid");
        assert_eq!(members.father.as_deref(), Some("dad"));
        assert_eq!(members.father_pheno, Some(Some(1)));
        assert!(members.siblings.is_empty());
        let error = identify_fam_members(&trio(), "H").unwrap_err().to_string();
        assert_eq!(error, "ValueError: Family H not found in pedigree table");
    }

    #[test]
    fn denovo_rules_follow_the_script() {
        assert_eq!(denovo("chr1", true, "0/0", "0/0", false), Denovo::Ps2);
        assert_eq!(denovo("chr1", true, "0/1", "0/0", false), Denovo::No);
        assert_eq!(denovo("chr1", false, "./.", "0/0", true), Denovo::Pm6);
        assert_eq!(denovo("chr1", false, "./.", "0/0", false), Denovo::No);
        assert_eq!(denovo("chrX", true, "0/1", "0/0", false), Denovo::Ps2);
        assert_eq!(denovo("chrX", false, "0/1", "0/0", false), Denovo::No);
        assert_eq!(denovo("chrY", false, "0", "1", false), Denovo::Ps2);
        assert_eq!(denovo("chrM", true, "0", ".", true), Denovo::No);
        assert_eq!(denovo("chrM", true, ".", "0", true), Denovo::Ps2);
    }

    #[test]
    fn criteria_judge_each_key_once() {
        let mut table = text_table(&[
            ("chrom", &["chr1", "chr1", "chr2"]),
            ("ref", &["A", "A", "C"]),
            ("alt", &["G", "G", "T"]),
            ("kid", &["0/1", "0/1", "0/1"]),
            ("dad", &["0/0", "0/1", "./."]),
            ("mum", &["0/0", "0/0", "0/0"]),
        ]);
        table.set_column("pos", Column::Int(vec![5, 5, 9]));
        table.set_column(
            "gnomAD_joint_AF",
            Column::Float(vec![f64::NAN, 0.1, f64::NAN]),
        );
        let output = ps2_pm6_criteria(&table, &trio(), "F").unwrap();
        assert_eq!(output.criterion("PS2"), Some(&Column::Int(vec![3, 3, 0])));
        assert_eq!(output.criterion("PM6"), Some(&Column::Int(vec![0, 0, 2])));
    }
}
