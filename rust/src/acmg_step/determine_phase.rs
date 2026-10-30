//! `scripts/determine_phase.py` (521 lines): cis and trans calls for pairs of variants in a gene.
//!
//! `determine_cis_trans_relationships` writes `is_pathogenic` into the table,
//! finds the patients of the pedigree (phenotype `"2"` among the table's
//! columns), and, when some patient has a parent ID or a phased genotype in
//! the first 50 called values, pairs each non-pathogenic variant of a patient
//! with every pathogenic variant of the same gene: the phase is that of the
//! last pathogenic variant in `variant_id` order. The per-gene results are
//! merged in the genes' first-appearance order (owner's decision 2026-10-30,
//! as main's Python does with `pool.imap`), so a variant analysed under two
//! genes takes the later gene's phase. Each patient gets a column
//! `<patient>_phase_patho` (`cis`, `trans`, `both`, `unknown`, or missing; a
//! float column of NaN when the patient has no phase at all, as
//! `Series.map({})` gives), and the counts of patients with a cis (or trans)
//! call, `both` counting for both.
//!
//! Reproduced on purpose: `Genotype.parse` takes the second character as the
//! delimiter (a one-character genotype raises IndexError, a third allele
//! raises ValueError); on chrY, chrM (any chromosome name holding `Y` or `M`)
//! and on a male's `X` the genotype is haploid; parent IDs other than `0` and
//! missing (`.`, `-`, empty text) count as parents for the short cut; the
//! memo of the Python loop never hits; every `Variant` reports its sex as
//! `female` in error messages.
//!
//! Refused (exit 1): a patient on two pedigree rows (Python fails on the
//! ambiguous truth value of a Series), genotype cells that are not text, and
//! the columns the grouping and parsing need with types the step never has.
use super::acmg_ps2_pm6_denovo::{harness_ped, python_str};
use super::types::{AcmgError, HarnessCall};
use crate::cli::{Result, invalid};
use crate::json_text::Json;
use crate::pandas_polars::first_appearance_groups;
use crate::pandas_table::{Cell, Column, ObjectValue, Table};
use std::collections::HashMap;

fn refused(message: impl Into<String>) -> Box<dyn std::error::Error> {
    invalid(format!("refused: {}", message.into()))
}

fn value_error(message: String) -> Box<dyn std::error::Error> {
    AcmgError::new("ValueError", message).into()
}

/// `Genotype` of the script.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Genotype {
    pub(crate) allele1: Option<String>,
    pub(crate) allele2: Option<String>,
    pub(crate) is_phased: bool,
    pub(crate) is_hemizygous: bool,
    pub(crate) chromosome: String,
}

/// Python's `repr` of a text, for the texts this tool prints back.
fn text_repr(text: &str) -> Result<String> {
    if text
        .chars()
        .all(|c| (' '..='~').contains(&c) && c != '\'' && c != '\\')
    {
        Ok(format!("'{text}'"))
    } else {
        Err(refused(format!(
            "the genotype text {text:?} would need Python's repr escapes"
        )))
    }
}

fn option_repr(value: &Option<String>) -> Result<String> {
    match value {
        Some(text) => text_repr(text),
        None => Ok("None".to_owned()),
    }
}

fn bool_repr(flag: bool) -> &'static str {
    if flag { "True" } else { "False" }
}

impl Genotype {
    fn new(
        allele1: Option<&str>,
        allele2: Option<&str>,
        is_phased: bool,
        is_hemizygous: bool,
        chromosome: &str,
    ) -> Genotype {
        Genotype {
            allele1: allele1.map(str::to_owned),
            allele2: allele2.map(str::to_owned),
            is_phased,
            is_hemizygous,
            chromosome: chromosome.to_owned(),
        }
    }

    /// `Genotype.parse(gt_str, patient_gender, chromosome)`.
    pub(crate) fn parse(gt: &str, gender: &str, chromosome: &str) -> Result<Genotype> {
        if gt.is_empty() || [".", "./.", ".|."].contains(&gt) {
            return Ok(Genotype::new(None, None, false, false, chromosome));
        }
        let Some(delimiter) = gt.chars().nth(1) else {
            return Err(AcmgError::new("IndexError", "string index out of range").into());
        };
        let invalid_genotype = || {
            value_error(format!(
                "Invalid genotype: {gt} for variant in patient (gender {gender}) at chromosome {chromosome}"
            ))
        };
        let alleles: Vec<&str> = gt.split(delimiter).collect();
        if (gender == "male" && chromosome.contains('X'))
            || chromosome.contains('Y')
            || chromosome.contains('M')
        {
            if alleles.len() == 1 && alleles[0] != "." {
                return Ok(Genotype::new(
                    Some(alleles[0]),
                    Some(alleles[0]),
                    true,
                    true,
                    chromosome,
                ));
            }
            if alleles.len() > 2 {
                return Err(value_error(
                    "too many values to unpack (expected 2)".to_owned(),
                ));
            }
            let (allele1, allele2) = (alleles[0], alleles[1]);
            return if allele1 == "." {
                Ok(Genotype::new(
                    Some(allele2),
                    Some(allele2),
                    true,
                    true,
                    chromosome,
                ))
            } else if allele2 == "." || allele1 == allele2 {
                Ok(Genotype::new(
                    Some(allele1),
                    Some(allele1),
                    true,
                    true,
                    chromosome,
                ))
            } else if alleles.contains(&"1") {
                Ok(Genotype::new(Some("1"), Some("1"), true, true, chromosome))
            } else {
                Err(invalid_genotype())
            };
        }
        if alleles.len() != 2 {
            return Err(invalid_genotype());
        }
        fn allele(text: &str) -> Option<&str> {
            (text != ".").then_some(text)
        }
        Ok(Genotype::new(
            allele(alleles[0]),
            allele(alleles[1]),
            delimiter == '|',
            false,
            chromosome,
        ))
    }

    pub(crate) fn has_alt(&self) -> bool {
        self.allele1.as_deref() == Some("1") || self.allele2.as_deref() == Some("1")
    }

    pub(crate) fn full_genotype(&self) -> bool {
        if self.is_hemizygous {
            self.allele1.is_some() || self.allele2.is_some()
        } else {
            self.allele1.is_some() && self.allele2.is_some()
        }
    }

    /// The dataclass `repr`, which one error message prints.
    fn repr(&self) -> Result<String> {
        Ok(format!(
            "Genotype(allele1={}, allele2={}, is_phased={}, is_hemizygous={}, chromosome={})",
            option_repr(&self.allele1)?,
            option_repr(&self.allele2)?,
            bool_repr(self.is_phased),
            bool_repr(self.is_hemizygous),
            text_repr(&self.chromosome)?
        ))
    }
}

/// The parts of the script's `Variant` the phase rules read.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Variant {
    pub(crate) variant_id: String,
    pub(crate) chromosome: String,
    pub(crate) patient_gt: Genotype,
    pub(crate) father_gt: Option<Genotype>,
    pub(crate) mother_gt: Option<Genotype>,
}

impl Variant {
    /// `self.determine_phase_with(other)` for two variants of one patient.
    pub(crate) fn determine_phase_with(&self, other: &Variant) -> Result<&'static str> {
        if self.chromosome != other.chromosome
            || !self.patient_gt.has_alt()
            || !other.patient_gt.has_alt()
        {
            return Ok("unknown");
        }
        if self.patient_gt.is_hemizygous && other.patient_gt.is_hemizygous {
            return Ok("cis");
        }
        if self.patient_gt.is_phased && other.patient_gt.is_phased {
            return self.phase_from_haplotypes(other);
        }
        Ok(self.phase_with_parents(other))
    }

    fn phase_from_haplotypes(&self, other: &Variant) -> Result<&'static str> {
        let homozygous = |gt: &Genotype| {
            gt.allele1.as_deref() == Some("1") && gt.allele2.as_deref() == Some("1")
        };
        if homozygous(&self.patient_gt) || homozygous(&other.patient_gt) {
            return Ok("both");
        }
        let (mine, theirs) = (&self.patient_gt, &other.patient_gt);
        if mine.allele1 == theirs.allele1 || mine.allele2 == theirs.allele2 {
            Ok("cis")
        } else if mine.allele1 == theirs.allele2 || mine.allele2 == theirs.allele1 {
            Ok("trans")
        } else {
            Err(value_error(format!(
                "Invalid genotype: {} for variant in patient (gender female) at chromosome {}, please try to convert the genotype to bi-alleleic first",
                mine.repr()?,
                self.chromosome
            )))
        }
    }

    fn phase_with_parents(&self, other: &Variant) -> &'static str {
        let (Some(self_father), Some(self_mother), Some(other_father), Some(other_mother)) = (
            &self.father_gt,
            &self.mother_gt,
            &other.father_gt,
            &other.mother_gt,
        ) else {
            return "unknown";
        };
        if !(self_father.full_genotype()
            && self_mother.full_genotype()
            && other_father.full_genotype()
            && other_mother.full_genotype())
        {
            return "unknown";
        }
        let (self_from_father, self_from_mother) = (self_father.has_alt(), self_mother.has_alt());
        let (other_from_father, other_from_mother) =
            (other_father.has_alt(), other_mother.has_alt());
        if !self_from_father && !self_from_mother || !other_from_father && !other_from_mother {
            return "unknown";
        }
        match (
            self_from_father,
            self_from_mother,
            other_from_father,
            other_from_mother,
        ) {
            (true, false, true, false) | (false, true, false, true) => "cis",
            (true, false, false, true) | (false, true, true, false) => "trans",
            _ => "unknown",
        }
    }
}

/// One entry of `build_patient_info`.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct PatientInfo {
    pub(crate) patient_id: String,
    pub(crate) gender: &'static str,
    pub(crate) father_id: Option<String>,
    pub(crate) mother_id: Option<String>,
}

/// `str(value) if pd.notna(value) and value not in [0, '0'] else None`.
fn parent_id(cell: Option<Cell<'_>>) -> Option<String> {
    match cell? {
        Cell::Missing | Cell::Text("0") | Cell::Int(0) | Cell::Bool(false) => None,
        Cell::Float(value) if value == 0.0 => None,
        other => python_str(other),
    }
}

/// `build_patient_info(pedigree_df, variants_df)`.
pub(crate) fn build_patient_info(ped: &Table, table: &Table) -> Result<Vec<PatientInfo>> {
    let ids = ped.column("IndividualID")?;
    let rows: Vec<usize> = (0..ped.rows)
        .filter(|row| matches!(ids.cell(*row), Cell::Text(id) if table.has_column(id)))
        .collect();
    let phenotypes = ped.column("Phenotype")?;
    let mut patients: Vec<&str> = Vec::new();
    for row in &rows {
        let affected = match phenotypes.cell(*row) {
            Cell::Text(text) => text == "2",
            Cell::Int(value) => value == 2,
            Cell::Float(value) => value == 2.0,
            _ => false,
        };
        if let (true, Cell::Text(id)) = (affected, ids.cell(*row))
            && !patients.contains(&id)
        {
            patients.push(id);
        }
    }
    let field = |name: &str, row: usize| ped.column_index(name).map(|column| ped.cell(column, row));
    let mut out = Vec::with_capacity(patients.len());
    for patient in patients {
        let own: Vec<usize> = rows
            .iter()
            .copied()
            .filter(|row| ids.cell(*row) == Cell::Text(patient))
            .collect();
        if own.len() != 1 {
            return Err(value_error(
                "The truth value of a Series is ambiguous. Use a.empty, a.bool(), a.item(), a.any() or a.all()."
                    .to_owned(),
            ));
        }
        let row = own[0];
        let gender = match field("Sex", row) {
            Some(Cell::Text("1") | Cell::Int(1) | Cell::Bool(true)) => "male",
            Some(Cell::Float(value)) if value == 1.0 => "male",
            Some(Cell::Text("2") | Cell::Int(2)) => "female",
            Some(Cell::Float(value)) if value == 2.0 => "female",
            _ => gender_from_chry(table, patient)?,
        };
        out.push(PatientInfo {
            patient_id: patient.to_owned(),
            gender,
            father_id: parent_id(field("PaternalID", row)),
            mother_id: parent_id(field("MaternalID", row)),
        });
    }
    Ok(out)
}

/// The fallback for an unknown sex: male when a genotype holding `1` lies on chrY.
fn gender_from_chry(table: &Table, patient: &str) -> Result<&'static str> {
    let genotypes = table.column(patient)?;
    let Column::Text(values) = genotypes else {
        return Err(refused(format!(
            "the sex of a patient is unknown and the sample column is of type {}, where Python's .str fails",
            genotypes.dtype_name()
        )));
    };
    if values.iter().any(Option::is_none) {
        return Err(value_error(
            "Cannot mask with non-boolean array containing NA / NaN values".to_owned(),
        ));
    }
    let chrom = table.column("chrom")?;
    let on_y = values.iter().enumerate().any(|(row, value)| {
        value.as_deref().is_some_and(|gt| gt.contains('1')) && chrom.cell(row) == Cell::Text("chrY")
    });
    Ok(if on_y { "male" } else { "female" })
}

/// `_can_resolve_phase(variants_df, patient_info)`.
fn can_resolve_phase(table: &Table, patients: &[PatientInfo]) -> Result<bool> {
    if patients
        .iter()
        .any(|info| info.father_id.is_some() || info.mother_id.is_some())
    {
        return Ok(true);
    }
    for info in patients {
        let column = table.column(&info.patient_id)?;
        let phased = (0..table.rows)
            .filter_map(|row| python_str(column.cell(row)))
            .take(50)
            .any(|text| text.contains('|'));
        if phased {
            return Ok(true);
        }
    }
    Ok(false)
}

/// A genotype cell: `None` when missing, the text otherwise.
fn genotype_cell<'a>(column: &'a Column, row: usize, what: &str) -> Result<Option<&'a str>> {
    match column.cell(row) {
        Cell::Missing => Ok(None),
        Cell::Text(text) => Ok(Some(text)),
        other => Err(refused(format!(
            "the {what} genotype {other:?} is not text"
        ))),
    }
}

/// The patients' variants of one gene, split into pathogenic and not
/// (`create_variants_from_gene_df`).
struct PatientVariants {
    pathogenic: Vec<Variant>,
    non_pathogenic: Vec<Variant>,
}

fn create_variants(
    table: &Table,
    gene_rows: &[usize],
    pathogenic: &[bool],
    patients: &[PatientInfo],
) -> Result<Vec<PatientVariants>> {
    let ids = table.column("variant_id")?;
    let Column::Text(id_values) = ids else {
        return Err(refused(format!(
            "variant_id of type {} is not grouped here",
            ids.dtype_name()
        )));
    };
    // groupby('variant_id'): sorted keys, missing keys dropped, rows in table order.
    let mut groups: Vec<(&str, Vec<usize>)> = Vec::new();
    let mut position: HashMap<&str, usize> = HashMap::new();
    for row in gene_rows {
        let Some(id) = id_values[*row].as_deref() else {
            continue;
        };
        match position.get(id) {
            Some(index) => groups[*index].1.push(*row),
            None => {
                position.insert(id, groups.len());
                groups.push((id, vec![*row]));
            }
        }
    }
    groups.sort_by(|left, right| left.0.cmp(right.0));
    let mut out: Vec<PatientVariants> = patients
        .iter()
        .map(|_| PatientVariants {
            pathogenic: Vec::new(),
            non_pathogenic: Vec::new(),
        })
        .collect();
    for (variant_id, rows) in groups {
        let is_pathogenic = rows.iter().any(|row| pathogenic[*row]);
        let first = rows[0];
        // first_row['chrom'], ['pos'], ['ref'], ['alt'] in the script's order.
        let chromosome = match table.column("chrom")?.cell(first) {
            Cell::Text(text) => text.to_owned(),
            Cell::Int(value) => value.to_string(),
            other => {
                return Err(refused(format!("chrom value {other:?} is not text")));
            }
        };
        match table.column("pos")?.cell(first) {
            Cell::Int(_) => {}
            Cell::Float(value) if value.is_finite() => {}
            Cell::Missing => {
                return Err(value_error(
                    "cannot convert float NaN to integer".to_owned(),
                ));
            }
            other => return Err(refused(format!("pos value {other:?} is not an integer"))),
        }
        table.column("ref")?;
        table.column("alt")?;
        for (index, info) in patients.iter().enumerate() {
            let column = table.column(&info.patient_id)?;
            let Some(text) = genotype_cell(column, first, "patient")? else {
                continue;
            };
            let patient_gt = Genotype::parse(text, info.gender, &chromosome)?;
            if !patient_gt.has_alt() {
                continue;
            }
            let parent = |id: &Option<String>, gender: &str| -> Result<Option<Genotype>> {
                let Some(id) = id.as_deref().filter(|id| !id.is_empty()) else {
                    return Ok(None);
                };
                let Some(column) = table.column_index(id) else {
                    return Ok(None);
                };
                match genotype_cell(&table.columns[column], first, "parent")? {
                    Some(text) => Genotype::parse(text, gender, &chromosome).map(Some),
                    None => Ok(None),
                }
            };
            let father_gt = parent(&info.father_id, "male")?;
            let mother_gt = parent(&info.mother_id, "female")?;
            let variant = Variant {
                variant_id: variant_id.to_owned(),
                chromosome: chromosome.clone(),
                patient_gt,
                father_gt,
                mother_gt,
            };
            if is_pathogenic {
                out[index].pathogenic.push(variant);
            } else {
                out[index].non_pathogenic.push(variant);
            }
        }
    }
    Ok(out)
}

/// What `determine_cis_trans_relationships` returns, with the columns it
/// writes into the table (`is_pathogenic`, then one phase column per patient).
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct PhaseResult {
    pub(crate) in_cis_pathogenic: Vec<i64>,
    pub(crate) in_trans_pathogenic: Vec<i64>,
    pub(crate) columns: Vec<(String, Column)>,
}

/// `determine_cis_trans_relationships(variants_df, pathogenic_array, pedigree_df)`.
pub(crate) fn determine_cis_trans_relationships(
    table: &Table,
    pathogenic: &[bool],
    ped: &Table,
) -> Result<PhaseResult> {
    let mut table = table.clone();
    table.set_column("is_pathogenic", Column::Bool(pathogenic.to_vec()));
    let mut columns = vec![(
        "is_pathogenic".to_owned(),
        Column::Bool(pathogenic.to_vec()),
    )];
    let patients = build_patient_info(ped, &table)?;
    let rows = table.rows;
    if !can_resolve_phase(&table, &patients)? {
        return Ok(PhaseResult {
            in_cis_pathogenic: vec![0; rows],
            in_trans_pathogenic: vec![0; rows],
            columns,
        });
    }
    if matches!(table.column("Gene")?, Column::Object(_)) {
        return Err(refused("an object Gene column mixing kinds of values"));
    }
    let mut phases: Vec<HashMap<String, &'static str>> = vec![HashMap::new(); patients.len()];
    for gene_rows in first_appearance_groups(&table, &["Gene"], true)? {
        let variants = create_variants(&table, &gene_rows, pathogenic, &patients)?;
        for (index, found) in variants.iter().enumerate() {
            let mut gene_phase: HashMap<&str, &'static str> = HashMap::new();
            for pvar in &found.pathogenic {
                for upvar in &found.non_pathogenic {
                    gene_phase.insert(&upvar.variant_id, upvar.determine_phase_with(pvar)?);
                }
            }
            for (variant_id, phase) in gene_phase {
                phases[index].insert(variant_id.to_owned(), phase);
            }
        }
    }
    let ids = table.column("variant_id")?;
    let mut in_cis = vec![0i64; rows];
    let mut in_trans = vec![0i64; rows];
    for (info, found) in patients.iter().zip(&phases) {
        let values: Vec<Option<&'static str>> = (0..rows)
            .map(|row| match ids.cell(row) {
                Cell::Text(id) => found.get(id).copied(),
                _ => None,
            })
            .collect();
        for (row, value) in values.iter().enumerate() {
            if matches!(value, Some("cis" | "both")) {
                in_cis[row] += 1;
            }
            if matches!(value, Some("trans" | "both")) {
                in_trans[row] += 1;
            }
        }
        let column = if found.is_empty() {
            Column::Float(vec![f64::NAN; rows])
        } else {
            Column::Text(
                values
                    .iter()
                    .map(|value| value.map(str::to_owned))
                    .collect(),
            )
        };
        columns.push((format!("{}_phase_patho", info.patient_id), column));
    }
    Ok(PhaseResult {
        in_cis_pathogenic: in_cis,
        in_trans_pathogenic: in_trans,
        columns,
    })
}

/// A table with the given columns written in (`df[name] = values`), and the
/// columns among them in the table's order, as the harness reports them.
pub(crate) fn written_columns(
    table: &Table,
    columns: &[(String, Column)],
) -> Vec<(String, Column)> {
    let mut written = table.clone();
    for (name, column) in columns {
        written.set_column(name, column.clone());
    }
    written
        .names
        .iter()
        .zip(&written.columns)
        .filter(|(name, _)| columns.iter().any(|(written, _)| written == *name))
        .map(|(name, column)| (name.clone(), column.clone()))
        .collect()
}

/// A JSON list of booleans or integers as a mask (nonzero is true).
pub(crate) fn json_mask(value: Option<&Json>, name: &str) -> Result<Vec<bool>> {
    match value {
        Some(Json::Array(items)) => items
            .iter()
            .map(|item| match item {
                Json::Bool(flag) => Ok(*flag),
                Json::Int(number) => Ok(*number != 0),
                _ => Err(invalid(format!(
                    "option {name} must hold booleans or integers"
                ))),
            })
            .collect(),
        _ => Err(invalid(format!("option {name} (a list) is required"))),
    }
}

fn text_cell(value: &Option<String>) -> ObjectValue {
    match value {
        Some(text) => ObjectValue::Text(text.clone()),
        None => ObjectValue::None,
    }
}

/// `priva acmg-step determine_phase <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let table = &call.table;
    match function {
        "Genotype.parse" => {
            // One parse per row (columns gt, gender, chromosome); an error is
            // reported in the row, as the Python capture does.
            let (gt, gender, chrom) = (
                table.column("gt")?,
                table.column("gender")?,
                table.column("chromosome")?,
            );
            let mut cells: Vec<Vec<ObjectValue>> = vec![Vec::new(); 5];
            for row in 0..table.rows {
                let text = |column: &Column| match column.cell(row) {
                    Cell::Text(text) => Ok(text.to_owned()),
                    _ => Err(invalid("frame cells must be text")),
                };
                let parsed = Genotype::parse(&text(gt)?, &text(gender)?, &text(chrom)?);
                let row_cells = match parsed {
                    Ok(found) => [
                        text_cell(&found.allele1),
                        text_cell(&found.allele2),
                        ObjectValue::Bool(found.is_phased),
                        ObjectValue::Bool(found.is_hemizygous),
                        ObjectValue::None,
                    ],
                    Err(error) => [
                        ObjectValue::None,
                        ObjectValue::None,
                        ObjectValue::None,
                        ObjectValue::None,
                        ObjectValue::Text(error.to_string()),
                    ],
                };
                for (column, value) in cells.iter_mut().zip(row_cells) {
                    column.push(value);
                }
            }
            Ok(Table {
                names: ["allele1", "allele2", "is_phased", "is_hemizygous", "error"]
                    .map(str::to_owned)
                    .to_vec(),
                columns: cells.into_iter().map(Column::Object).collect(),
                rows: table.rows,
            })
        }
        "determine_phase_with" => {
            // Columns chromosome, gender, then p1, f1, m1 (self) and p2, f2, m2
            // (other): genotype texts or None.
            let text = |name: &str, row: usize| -> Result<Option<String>> {
                Ok(match table.column(name)?.cell(row) {
                    Cell::Text(text) => Some(text.to_owned()),
                    _ => None,
                })
            };
            let mut out = Vec::with_capacity(table.rows);
            for row in 0..table.rows {
                let chromosome = text("chromosome", row)?.unwrap_or_default();
                let gender = text("gender", row)?.unwrap_or_default();
                let make = |p: &str, f: &str, m: &str, id: &str| -> Result<Variant> {
                    let parse = |name: &str, gender: &str| -> Result<Option<Genotype>> {
                        text(name, row)?
                            .map(|gt| Genotype::parse(&gt, gender, &chromosome))
                            .transpose()
                    };
                    Ok(Variant {
                        variant_id: id.to_owned(),
                        chromosome: chromosome.clone(),
                        patient_gt: parse(p, &gender)?
                            .ok_or_else(|| invalid("p1/p2 are required"))?,
                        father_gt: parse(f, "male")?,
                        mother_gt: parse(m, "female")?,
                    })
                };
                let result = make("p1", "f1", "m1", "v1").and_then(|first| {
                    let second = make("p2", "f2", "m2", "v2")?;
                    first.determine_phase_with(&second).map(str::to_owned)
                });
                out.push(Some(match result {
                    Ok(phase) => phase,
                    Err(error) => format!("error {error}"),
                }));
            }
            Ok(Table {
                names: vec!["phase".to_owned()],
                columns: vec![Column::Text(out)],
                rows: table.rows,
            })
        }
        "build_patient_info" => {
            let ped = harness_ped(call)?;
            let patients = build_patient_info(&ped, table)?;
            let column = |pick: &dyn Fn(&PatientInfo) -> Option<String>| {
                Column::Object(patients.iter().map(|info| text_cell(&pick(info))).collect())
            };
            Ok(Table {
                names: ["patient", "gender", "father_id", "mother_id"]
                    .map(str::to_owned)
                    .to_vec(),
                columns: vec![
                    column(&|info| Some(info.patient_id.clone())),
                    column(&|info| Some(info.gender.to_owned())),
                    column(&|info| info.father_id.clone()),
                    column(&|info| info.mother_id.clone()),
                ],
                rows: patients.len(),
            })
        }
        "determine_cis_trans_relationships" => {
            let ped = harness_ped(call)?;
            let pathogenic = json_mask(call.options.get("pathogenic"), "pathogenic")?;
            let result = determine_cis_trans_relationships(table, &pathogenic, &ped)?;
            let mut names = vec!["in_cis".to_owned(), "in_trans".to_owned()];
            let mut columns = vec![
                Column::Int(result.in_cis_pathogenic),
                Column::Int(result.in_trans_pathogenic),
            ];
            for (name, column) in written_columns(table, &result.columns) {
                names.push(name);
                columns.push(column);
            }
            Ok(Table {
                names,
                columns,
                rows: table.rows,
            })
        }
        _ => Err(super::types::not_converted("determine_phase", function)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gt(text: &str, gender: &str, chromosome: &str) -> Result<Genotype> {
        Genotype::parse(text, gender, chromosome)
    }

    #[test]
    fn genotype_grammar_follows_the_script() {
        let het = gt("0/1", "female", "chr1").unwrap();
        assert_eq!(
            (
                het.allele1.as_deref(),
                het.allele2.as_deref(),
                het.is_phased
            ),
            (Some("0"), Some("1"), false)
        );
        assert!(gt("0|1", "male", "chr1").unwrap().is_phased);
        let male_x = gt("0/1", "male", "chrX").unwrap();
        assert!(male_x.is_hemizygous && male_x.allele1.as_deref() == Some("1"));
        assert_eq!(
            gt("./1", "female", "chrY").unwrap().allele1.as_deref(),
            Some("1")
        );
        assert!(
            gt("0/2", "female", "chrM")
                .unwrap_err()
                .to_string()
                .starts_with("ValueError: Invalid genotype: 0/2")
        );
        assert_eq!(
            gt("1", "female", "chr1").unwrap_err().to_string(),
            "IndexError: string index out of range"
        );
        assert_eq!(
            gt("0/1/1", "male", "chrX").unwrap_err().to_string(),
            "ValueError: too many values to unpack (expected 2)"
        );
        assert!(gt("0/1/1", "male", "chr1").is_err());
        assert_eq!(gt("", "male", "chr1").unwrap().allele1, None);
        // The second character is the delimiter, whatever it is.
        let odd = gt("10/1", "female", "chr1").unwrap();
        assert_eq!(
            (odd.allele1.as_deref(), odd.allele2.as_deref()),
            (Some("1"), Some("/1"))
        );
    }

    fn variant(id: &str, patient: &str, father: Option<&str>, mother: Option<&str>) -> Variant {
        Variant {
            variant_id: id.to_owned(),
            chromosome: "chr1".to_owned(),
            patient_gt: gt(patient, "female", "chr1").unwrap(),
            father_gt: father.map(|text| gt(text, "male", "chr1").unwrap()),
            mother_gt: mother.map(|text| gt(text, "female", "chr1").unwrap()),
        }
    }

    #[test]
    fn phase_rules_follow_the_script() {
        let from_father = variant("a", "0/1", Some("0/1"), Some("0/0"));
        let from_mother = variant("b", "0/1", Some("0/0"), Some("1/1"));
        let either = variant("c", "0/1", Some("0/1"), Some("0/1"));
        assert_eq!(
            from_father.determine_phase_with(&from_mother).unwrap(),
            "trans"
        );
        assert_eq!(
            from_father.determine_phase_with(&from_father).unwrap(),
            "cis"
        );
        assert_eq!(
            either.determine_phase_with(&from_father).unwrap(),
            "unknown"
        );
        let phased_a = variant("d", "0|1", None, None);
        let phased_b = variant("e", "1|0", None, None);
        let hom = variant("f", "1|1", None, None);
        assert_eq!(phased_a.determine_phase_with(&phased_b).unwrap(), "trans");
        assert_eq!(phased_a.determine_phase_with(&hom).unwrap(), "both");
        assert_eq!(
            phased_a.determine_phase_with(&from_father).unwrap(),
            "unknown"
        );
    }

    fn table(columns: Vec<(&str, Column)>) -> Table {
        Table {
            rows: columns[0].1.len(),
            names: columns.iter().map(|(name, _)| (*name).to_owned()).collect(),
            columns: columns.into_iter().map(|(_, column)| column).collect(),
        }
    }

    fn texts(values: &[&str]) -> Column {
        Column::Text(
            values
                .iter()
                .map(|value| Some((*value).to_owned()))
                .collect(),
        )
    }

    #[test]
    fn a_variant_under_two_genes_takes_the_later_genes_phase() {
        let variants = table(vec![
            ("variant_id", texts(&["v1", "v2", "v1", "v3"])),
            ("Gene", texts(&["A", "A", "B", "B"])),
            ("chrom", texts(&["chr1", "chr1", "chr1", "chr1"])),
            ("pos", Column::Int(vec![1, 2, 1, 3])),
            ("ref", texts(&["A", "C", "A", "G"])),
            ("alt", texts(&["G", "T", "G", "A"])),
            ("kid", texts(&["0/1", "0/1", "0/1", "0/1"])),
            ("dad", texts(&["0/1", "0/1", "0/1", "0/0"])),
            ("mum", texts(&["0/0", "0/0", "0/0", "0/1"])),
        ]);
        let ped = table(vec![
            ("#FamilyID", texts(&["F", "F", "F"])),
            ("IndividualID", texts(&["kid", "dad", "mum"])),
            ("PaternalID", texts(&["dad", "0", "0"])),
            ("MaternalID", texts(&["mum", "0", "0"])),
            ("Sex", texts(&["2", "1", "2"])),
            ("Phenotype", texts(&["2", "1", "1"])),
        ]);
        let result =
            determine_cis_trans_relationships(&variants, &[false, true, false, true], &ped)
                .unwrap();
        // v1 is cis with v2 in gene A, trans with v3 in gene B; B comes later.
        assert_eq!(result.columns[1].0, "kid_phase_patho");
        assert_eq!(
            result.columns[1].1,
            Column::Text(vec![Some("trans".into()), None, Some("trans".into()), None])
        );
        assert_eq!(result.in_trans_pathogenic, vec![1, 0, 1, 0]);
        assert_eq!(result.in_cis_pathogenic, vec![0, 0, 0, 0]);
    }
}
