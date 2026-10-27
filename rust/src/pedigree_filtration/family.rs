//! Ordered Polars selections, distinct from the ranking pedigree_context module.
use crate::cli::{Result, invalid};
use polars_core::prelude::*;

pub const COLUMNS: [&str; 6] = [
    "#FamilyID",
    "IndividualID",
    "PaternalID",
    "MaternalID",
    "Sex",
    "Phenotype",
];

pub struct FamilyPlan {
    pub positions: Vec<usize>,
    pub patients: Vec<String>,
    pub controls: Vec<String>,
    pub proband: Option<String>,
    pub sex: Option<String>,
}

pub fn select(rows: &[Vec<String>], target: Option<&str>) -> Result<FamilyPlan> {
    if rows.iter().any(|row| row.len() != 6) {
        return Err(invalid(
            "pedigree rows require exactly six parsed string cells",
        ));
    }
    let mut columns = Vec::with_capacity(7);
    for (position, name) in COLUMNS.iter().enumerate() {
        let values: Vec<&str> = rows.iter().map(|row| row[position].as_str()).collect();
        columns.push(Column::new((*name).into(), values));
    }
    columns.push(Column::new(
        "source_row".into(),
        (0..rows.len() as u64).collect::<Vec<_>>(),
    ));
    let frame = DataFrame::new(columns)?;
    let mask = match target {
        Some(value) => frame.column("#FamilyID")?.str()?.equal(value),
        None => BooleanChunked::full("family_mask".into(), false, rows.len()),
    };
    let selected = frame.filter(&mask)?;
    let positions = selected
        .column("source_row")?
        .u64()?
        .into_no_null_iter()
        .map(|value| value as usize)
        .collect();
    let patients_frame = selected.filter(&selected.column("Phenotype")?.str()?.equal("2"))?;
    let controls_frame = selected.filter(&selected.column("Phenotype")?.str()?.equal("1"))?;
    let patients: Vec<String> = patients_frame
        .column("IndividualID")?
        .str()?
        .into_no_null_iter()
        .map(str::to_owned)
        .collect();
    let controls = controls_frame
        .column("IndividualID")?
        .str()?
        .into_no_null_iter()
        .map(str::to_owned)
        .collect();
    let proband = patients.first().cloned();
    let sex = if let Some(name) = proband.as_deref() {
        let first = selected.filter(&selected.column("IndividualID")?.str()?.equal(name))?;
        first.column("Sex")?.str()?.get(0).map(str::to_owned)
    } else {
        None
    };
    Ok(FamilyPlan {
        positions,
        patients,
        controls,
        proband,
        sex,
    })
}

pub fn first_parent(values: &[Option<String>]) -> Option<(usize, String)> {
    // dropna().drop_duplicates().tolist()[0]: first surviving unique is the
    // first surviving item; retaining its source ordinal protects that order.
    values
        .iter()
        .enumerate()
        .find_map(|(position, value)| value.clone().map(|text| (position, text)))
}
