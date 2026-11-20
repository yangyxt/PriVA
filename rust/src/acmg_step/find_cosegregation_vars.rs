//! `scripts/find_cosegregation_vars.py` (228 lines): PP1 co-segregation counts from a multi-family
//! VCF.
//!
//! The pedigree is read as the script reads it (`pd.read_csv(sep="\t",
//! header=None, names=[six names])`, pandas' type inference and missing
//! values; a header line in the file is a data row), sex and phenotype go
//! through `pd.to_numeric(errors="coerce").astype("Int64")`, and families are
//! taken in sorted `FamilyID` order. For each VCF record, every family must
//! have all its affected members carrying the first ALT allele (allele `1`)
//! and no unaffected member homozygous or hemizygous; the counts are those of
//! `find_cosegregating_variants`, kept per `(chrom, pos, ref, alts[0])` (a
//! later record of the same key replaces an earlier one).
//!
//! The VCF is read the way pysam gives `record.samples[name]['GT']`: allele
//! indexes, `.` and alleles beyond the record's ALT count missing (`parse_gt`
//! of `pedigree_filter`). A pedigree member who is not a VCF sample raises
//! pysam's KeyError when first looked up; a record with ALT `.` raises the
//! TypeError of `record.alts[0]`.
//!
//! Refused (exit 1): numeric individual IDs (pysam would take them as sample
//! positions), a record without GT whose samples are looked up, a GT that
//! htslib would not parse (an empty GT stops pysam with "truncated file"),
//! records with a sample count other than the header's, and a VCF not starting
//! with `##fileformat=VCF`.
use super::acmg_ps2_pm6_denovo::{harness_text, ped_to_nullable_int};
use super::types::{AcmgError, HarnessCall};
use crate::acmg_scoring::python_strip;
use crate::cli::{Result, invalid};
use crate::pandas_table::{Cell, Column, Header, ReadOptions, Table, read_table_with};
use crate::pedigree_filter::{next_line, parse_gt};
use crate::vcf_text::open_reader;
use std::collections::HashMap;
use std::path::Path;

fn refused(message: impl Into<String>) -> Box<dyn std::error::Error> {
    invalid(format!("refused: {}", message.into()))
}

/// The counts `find_cosegregating_variants` keeps for a variant, in the
/// order of the Python dict.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct SegregationCounts {
    pub(crate) affected_segregated: i64,
    pub(crate) unaffected_segregated: i64,
    pub(crate) affected_homozygous: i64,
    pub(crate) unaffected_nonhomo: i64,
    pub(crate) male_affected_segregated: i64,
    pub(crate) male_unaffected_segregated: i64,
    pub(crate) female_affected_segregated: i64,
    pub(crate) female_unaffected_segregated: i64,
}

impl SegregationCounts {
    pub(crate) const NAMES: [&'static str; 8] = [
        "Affected_segregated_inds",
        "Unaffected_segregated_inds",
        "Affected_homozygous_inds",
        "Unaffected_nonhomo_inds",
        "Male_affected_segregated_inds",
        "Male_unaffected_segregated_inds",
        "Female_affected_segregated_inds",
        "Female_unaffected_segregated_inds",
    ];

    pub(crate) fn values(&self) -> [i64; 8] {
        [
            self.affected_segregated,
            self.unaffected_segregated,
            self.affected_homozygous,
            self.unaffected_nonhomo,
            self.male_affected_segregated,
            self.male_unaffected_segregated,
            self.female_affected_segregated,
            self.female_unaffected_segregated,
        ]
    }
}

/// The variant key `(chrom, pos, ref, alts[0])`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct VariantKey {
    pub(crate) chrom: String,
    pub(crate) pos: i64,
    pub(crate) reference: String,
    pub(crate) alt: String,
}

/// What `check_variant_presence` returns.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Presence {
    None,
    False,
    True,
    Heterozygous,
    Homozygous,
    Hemizygous,
}

impl Presence {
    fn truthy(self) -> bool {
        !matches!(self, Presence::None | Presence::False)
    }
}

/// `_normalized_sex(gender)` of an `Int64` value.
fn normalized_sex(sex: Option<i64>) -> Option<i64> {
    sex.filter(|value| *value == 1 || *value == 2)
}

/// `_normalized_chromosome(chrom)`.
fn normalized_chromosome(chrom: &str) -> String {
    let upper = python_strip(chrom).to_uppercase();
    upper.strip_prefix("CHR").unwrap_or(&upper).to_owned()
}

/// `check_variant_presence(parse_gt(gt), gender, chrom)`.
fn check_variant_presence(
    genotype: &[Option<u32>],
    sex: Option<i64>,
    chromosome: &str,
) -> Presence {
    let alleles: Vec<u32> = genotype.iter().flatten().copied().collect();
    if alleles.is_empty() {
        return Presence::None;
    }
    let count = alleles.iter().filter(|allele| **allele == 1).count();
    match (chromosome, sex) {
        ("Y", Some(1)) if count >= 1 => Presence::True,
        ("Y", Some(1) | Some(2)) => Presence::False,
        ("Y", _) => Presence::None,
        ("X", Some(1)) if count >= 1 => Presence::Hemizygous,
        ("X", Some(1)) => Presence::False,
        ("X", None) => Presence::None,
        _ if count > 1 => Presence::Homozygous,
        _ if count == 1 => Presence::Heterozygous,
        _ => Presence::False,
    }
}

/// One family of the pedigree: `(IndividualID, Sex)` of its affected and
/// unaffected members, and its parent IDs (text IDs only: a numeric parent
/// ID never equals a text individual ID).
#[derive(Clone, Debug, Default)]
struct Family {
    affected: Vec<(String, Option<i64>)>,
    unaffected: Vec<(String, Option<i64>)>,
    parents: Vec<String>,
}

/// The first parent ID of a family that is not missing and not `0`, `"0"`, `"."` or `"-"`.
fn family_parent(column: &Column, rows: &[usize]) -> Option<Option<String>> {
    rows.iter().find_map(|row| match column.cell(*row) {
        Cell::Missing | Cell::Text("0" | "." | "-") | Cell::Int(0) | Cell::Bool(false) => None,
        Cell::Float(value) if value == 0.0 => None,
        Cell::Text(text) => Some(Some(text.to_owned())),
        _ => Some(None),
    })
}

/// The pedigree as `find_cosegregating_variants` reads it, families in `groupby` order.
fn read_families(ped_file: &Path, mode: &str) -> Result<Vec<Family>> {
    let names = [
        "FamilyID",
        "IndividualID",
        "PaternalID",
        "MaternalID",
        "Sex",
        "Phenotype",
    ];
    let table = read_table_with(
        ped_file,
        &ReadOptions {
            header: Header::Names(names.map(str::to_owned).to_vec()),
            ..ReadOptions::default()
        },
    )?;
    if !["dominant", "recessive", "both"].contains(&mode) {
        return Err(AcmgError::new(
            "ValueError",
            format!("Unsupported cosegregation mode: {mode}"),
        )
        .into());
    }
    let sexes: Vec<Option<i64>> = (0..table.rows)
        .map(|row| ped_to_nullable_int(table.cell(4, row)))
        .collect::<Result<_>>()?;
    let phenotypes: Vec<Option<i64>> = (0..table.rows)
        .map(|row| ped_to_nullable_int(table.cell(5, row)))
        .collect::<Result<_>>()?;
    // groupby('FamilyID'): sorted keys, missing keys dropped.
    let family_ids = &table.columns[0];
    let mut keyed: Vec<(FamilyKey, Vec<usize>)> = Vec::new();
    for row in 0..table.rows {
        let key = match family_ids.cell(row) {
            Cell::Missing => continue,
            Cell::Text(text) => FamilyKey::Text(text.to_owned()),
            Cell::Int(value) => FamilyKey::Number(value as f64),
            Cell::Float(value) => FamilyKey::Number(value),
            Cell::Bool(_) => return Err(refused("boolean family IDs")),
        };
        match keyed.iter_mut().find(|(existing, _)| *existing == key) {
            Some((_, rows)) => rows.push(row),
            None => keyed.push((key, vec![row])),
        }
    }
    if matches!(family_ids, Column::Int(values) if values.iter().any(|value| value.unsigned_abs() > 1 << 53))
    {
        return Err(refused("family IDs beyond 2**53"));
    }
    keyed.sort_by(|left, right| left.0.order(&right.0));
    let ids = &table.columns[1];
    let member = |row: usize| -> Result<String> {
        match (ids, ids.cell(row)) {
            (Column::Text(_), Cell::Text(text)) => Ok(text.to_owned()),
            (_, cell) => Err(refused(format!(
                "individual ID {cell:?} (column type {}) is not text; pysam would look a number up by sample position",
                ids.dtype_name()
            ))),
        }
    };
    let mut families = Vec::with_capacity(keyed.len());
    for (_, rows) in keyed {
        let mut family = Family::default();
        for row in &rows {
            match phenotypes[*row] {
                Some(2) => family.affected.push((member(*row)?, sexes[*row])),
                Some(1) => family.unaffected.push((member(*row)?, sexes[*row])),
                _ => {}
            }
        }
        for column in [&table.columns[2], &table.columns[3]] {
            if let Some(Some(parent)) = family_parent(column, &rows) {
                family.parents.push(parent);
            }
        }
        families.push(family);
    }
    Ok(families)
}

#[derive(Clone, Debug, PartialEq)]
enum FamilyKey {
    Text(String),
    Number(f64),
}

impl FamilyKey {
    fn order(&self, other: &FamilyKey) -> std::cmp::Ordering {
        match (self, other) {
            (FamilyKey::Text(left), FamilyKey::Text(right)) => left.cmp(right),
            (FamilyKey::Number(left), FamilyKey::Number(right)) => left.total_cmp(right),
            _ => std::cmp::Ordering::Equal,
        }
    }
}

/// One VCF record's fields the function reads.
struct Record<'a> {
    chrom: &'a str,
    pos: i64,
    reference: &'a str,
    alts: Option<&'a str>,
    fields: Vec<&'a str>,
    gt_index: Option<usize>,
    alleles: u32,
}

impl Record<'_> {
    /// `record.samples[name]['GT']`.
    fn genotype(&self, column: usize, line: usize) -> Result<Vec<Option<u32>>> {
        let Some(index) = self.gt_index else {
            return Err(refused(format!(
                "line {line}: a pedigree member is looked up in a record without GT"
            )));
        };
        let value = self.fields[column].split(':').nth(index).unwrap_or(".");
        parse_gt(value, self.alleles)
    }
}

/// `find_cosegregating_variants(vcf_file, ped_file, mode)`.
pub(crate) fn find_cosegregating_variants(
    vcf_file: &Path,
    ped_file: &Path,
    mode: &str,
) -> Result<Vec<(VariantKey, SegregationCounts)>> {
    let families = read_families(ped_file, mode)?;
    let mut reader = open_reader(vcf_file)?;
    let (mut buffer, mut number) = (Vec::new(), 0usize);
    let mut first = true;
    let samples: HashMap<String, usize> = loop {
        let Some(line) = next_line(reader.as_mut(), &mut buffer, &mut number)? else {
            return Err(refused("the VCF has no #CHROM line"));
        };
        if first && !line.starts_with("##fileformat=VCF") {
            return Err(refused("the file does not start with ##fileformat=VCF"));
        }
        first = false;
        if line.starts_with("#CHROM") {
            let names: Vec<&str> = line.split('\t').skip(9).collect();
            let mut samples = HashMap::new();
            for (index, name) in names.iter().enumerate() {
                if samples.insert((*name).to_owned(), index + 9).is_some() {
                    return Err(refused(format!(
                        "sample {name:?} appears twice in the header"
                    )));
                }
            }
            break samples;
        }
        if !line.starts_with("##") {
            return Err(refused(format!("line {number}: a record before #CHROM")));
        }
    };
    let width = 9 + samples.len();
    let mut kept: Vec<(VariantKey, SegregationCounts)> = Vec::new();
    let mut position: HashMap<VariantKey, usize> = HashMap::new();
    while let Some(line) = next_line(reader.as_mut(), &mut buffer, &mut number)? {
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != width && !(samples.is_empty() && fields.len() == 8) {
            return Err(refused(format!(
                "line {number}: {} columns where the header has {width}",
                fields.len()
            )));
        }
        let pos: i64 = fields[1]
            .parse()
            .map_err(|_| refused(format!("line {number}: POS is not an integer")))?;
        let gt_index = fields
            .get(8)
            .and_then(|format| format.split(':').position(|key| key == "GT"));
        let record = Record {
            chrom: fields[0],
            pos,
            reference: fields[3],
            alts: (fields[4] != ".").then_some(fields[4]),
            alleles: match fields[4] {
                "." => 1,
                alts => 1 + alts.split(',').count() as u32,
            },
            gt_index,
            fields: fields.clone(),
        };
        if let Some(index) = gt_index {
            for column in 9..fields.len() {
                let value = fields[column].split(':').nth(index).unwrap_or(".");
                if value.is_empty() {
                    return Err(refused(format!(
                        "line {number}: an empty GT, on which pysam stops with \"truncated file\""
                    )));
                }
                parse_gt(value, record.alleles)
                    .map_err(|error| refused(format!("line {number}: {error}")))?;
            }
        }
        let Some(alts) = record.alts else {
            return Err(
                AcmgError::new("TypeError", "'NoneType' object is not subscriptable").into(),
            );
        };
        let key = VariantKey {
            chrom: record.chrom.to_owned(),
            pos: record.pos,
            reference: record.reference.to_owned(),
            alt: alts.split(',').next().unwrap_or_default().to_owned(),
        };
        if let Some(counts) = segregate(&record, &families, &samples, number)? {
            match position.get(&key) {
                Some(index) => kept[*index].1 = counts,
                None => {
                    position.insert(key.clone(), kept.len());
                    kept.push((key, counts));
                }
            }
        }
    }
    Ok(kept)
}

/// The counts of one record, or `None` when some family does not segregate.
fn segregate(
    record: &Record<'_>,
    families: &[Family],
    samples: &HashMap<String, usize>,
    line: usize,
) -> Result<Option<SegregationCounts>> {
    let chromosome = normalized_chromosome(record.chrom);
    let sex_linked = chromosome == "X" || chromosome == "Y";
    let genotype = |id: &str| -> Result<Vec<Option<u32>>> {
        match samples.get(id) {
            Some(column) => record.genotype(*column, line),
            None => Err(AcmgError::new("KeyError", format!("'invalid sample name: {id}'")).into()),
        }
    };
    let mut counts = SegregationCounts::default();
    for family in families {
        let mut affected: Vec<Presence> = Vec::new();
        let (mut male_affected, mut female_affected) = (0, 0);
        for (id, sex) in &family.affected {
            let sex = normalized_sex(*sex);
            if sex_linked && sex.is_none() {
                continue;
            }
            let presence = check_variant_presence(&genotype(id)?, sex, &chromosome);
            if !presence.truthy() {
                return Ok(None);
            }
            affected.push(presence);
            match sex {
                Some(1) => male_affected += 1,
                Some(2) => female_affected += 1,
                _ => {}
            }
        }
        let mut controls: Vec<Presence> = Vec::new();
        let (mut nonparental, mut male_nonparental, mut female_nonparental) = (0, 0, 0);
        for (id, sex) in &family.unaffected {
            let sex = normalized_sex(*sex);
            if sex_linked && sex.is_none() {
                continue;
            }
            let presence = check_variant_presence(&genotype(id)?, sex, &chromosome);
            if matches!(presence, Presence::Homozygous | Presence::Hemizygous) {
                return Ok(None);
            }
            controls.push(presence);
            if !family.parents.iter().any(|parent| parent == id) {
                nonparental += 1;
                match sex {
                    Some(1) => male_nonparental += 1,
                    Some(2) => female_nonparental += 1,
                    _ => {}
                }
            }
        }
        if !controls.iter().any(|presence| presence.truthy()) {
            counts.affected_segregated += affected.len() as i64;
            counts.unaffected_segregated += nonparental;
            counts.male_affected_segregated += male_affected;
            counts.male_unaffected_segregated += male_nonparental;
            counts.female_affected_segregated += female_affected;
            counts.female_unaffected_segregated += female_nonparental;
        } else if affected.contains(&Presence::Heterozygous) {
            return Ok(None);
        } else {
            counts.affected_homozygous += affected.len() as i64;
            counts.unaffected_nonhomo += nonparental;
        }
    }
    Ok(Some(counts))
}

/// `priva acmg-step find_cosegregation_vars <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    match function {
        "find_cosegregating_variants" => {
            let found = find_cosegregating_variants(
                Path::new(harness_text(call, "vcf")?),
                Path::new(harness_text(call, "ped")?),
                harness_text(call, "mode")?,
            )?;
            let mut names = vec!["chrom", "pos", "ref", "alt"];
            names.extend(SegregationCounts::NAMES);
            let mut columns = vec![
                Column::Text(
                    found
                        .iter()
                        .map(|(key, _)| Some(key.chrom.clone()))
                        .collect(),
                ),
                Column::Int(found.iter().map(|(key, _)| key.pos).collect()),
                Column::Text(
                    found
                        .iter()
                        .map(|(key, _)| Some(key.reference.clone()))
                        .collect(),
                ),
                Column::Text(found.iter().map(|(key, _)| Some(key.alt.clone())).collect()),
            ];
            for index in 0..8 {
                columns.push(Column::Int(
                    found
                        .iter()
                        .map(|(_, counts)| counts.values()[index])
                        .collect(),
                ));
            }
            Ok(Table {
                names: names.into_iter().map(str::to_owned).collect(),
                columns,
                rows: found.len(),
            })
        }
        _ => Err(super::types::not_converted(
            "find_cosegregation_vars",
            function,
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn presence_follows_the_python_tests() {
        assert_eq!(
            check_variant_presence(&[Some(1)], Some(1), "X"),
            Presence::Hemizygous
        );
        assert_eq!(
            check_variant_presence(&[Some(1), Some(1)], Some(1), "X"),
            Presence::Hemizygous
        );
        assert_eq!(
            check_variant_presence(&[Some(0), Some(1)], Some(2), "X"),
            Presence::Heterozygous
        );
        assert_eq!(
            check_variant_presence(&[Some(1), Some(1)], Some(2), "X"),
            Presence::Homozygous
        );
        assert_eq!(
            check_variant_presence(&[Some(0), Some(1)], None, "X"),
            Presence::None
        );
        assert_eq!(
            check_variant_presence(&[None, None], Some(1), "1"),
            Presence::None
        );
        assert_eq!(
            check_variant_presence(&[Some(2), Some(2)], Some(1), "1"),
            Presence::False
        );
        assert_eq!(normalized_chromosome(" chrx "), "X");
    }

    #[test]
    fn counts_known_x_sexes_and_excludes_unknown() {
        // tests/test_find_cosegregation_vars.py, second test.
        let folder = std::env::temp_dir().join(format!("priva_coseg_{}", std::process::id()));
        fs::create_dir_all(&folder).unwrap();
        let ped = folder.join("family.ped");
        fs::write(
            &ped,
            "F1\tMALE\t0\t0\t1\t2\nF1\tFEMALE\t0\t0\t2\t2\nF1\tUNKNOWN\t0\t0\t3\t2\n",
        )
        .unwrap();
        let vcf = folder.join("family.vcf");
        fs::write(
            &vcf,
            "##fileformat=VCFv4.2\n##contig=<ID=X>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMALE\tFEMALE\tUNKNOWN\n\
             X\t100\t.\tA\tG\t.\tPASS\t.\tGT\t1\t0/1\t0/1\n",
        )
        .unwrap();
        let found = find_cosegregating_variants(&vcf, &ped, "both").unwrap();
        fs::remove_dir_all(&folder).unwrap();
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].0.pos, 100);
        assert_eq!(found[0].1.values(), [2, 0, 0, 0, 1, 0, 1, 0]);
    }
}
