//! Rust version of `scripts/pedigree_filtration_per_fam.py`.
//!
//! Keeps the records of a family VCF that at least one patient (pedigree
//! phenotype 2) carries, with allele 1 or 2 in GT, and that no control
//! (phenotype 1) has as a homozygous or hemizygous non-reference genotype.
//! Each kept record gets INFO VARIANT_SOURCE: PATERNAL, MATERNAL, DENOVO,
//! INHERITED or UNKNOWN.
//!
//! Behaviour reproduced from the Python script:
//! - the pedigree is read as pandas does with `dtype=str, keep_default_na=False`:
//!   every line is a row (a `#FamilyID` header line too), only the first six
//!   tab-separated fields count and missing trailing fields are missing values;
//! - the proband is the family's first patient row; the father and the mother
//!   are the first PaternalID and MaternalID of any family row that is not a
//!   missing value (`0`, `.`, empty, `NA`, `nan`, `NaN`, `-`; `-9` stays an ID),
//!   not necessarily the proband's own parents;
//! - the proband carries a record only with allele 1, so allele 2 alone gives
//!   UNKNOWN, while patients and controls are judged on alleles 1 and 2;
//! - a family missing from the pedigree logs an error and exits 0 without
//!   writing the output; a family without a patient, or a patient sample in a
//!   record without GT, stops the run with exit code 1;
//! - pysam answers a lookup of a sample name the VCF lacks with a KeyError, so
//!   a family member who is not a VCF sample stops the run when the script
//!   first looks them up: a patient at the first record, a control at the
//!   first record a patient carries, the proband's father or mother at the
//!   first kept record the proband carries with allele 1. The script checks
//!   patients and controls in set order, which changes from run to run; they
//!   are checked here in pedigree order;
//! - a VCF whose header already defines INFO VARIANT_SOURCE (the script's own
//!   output) is refused before the output is opened: pysam will not add the
//!   definition twice.
//!
//! Records are written the way pysam writes them (see `vcf_text`): the PASS
//! filter line second in the header, the VARIANT_SOURCE definition appended
//! last, INFO END synchronised. Other columns are copied unchanged, which
//! equals htslib's own formatting for bcftools-written input such as PriVA's
//! family VCFs. When a record stops the run, the output keeps the header and
//! the records before it, as pysam closes the file.
use crate::cli::{Result, invalid};
use crate::vcf_text::{
    TextSink, create_writer, info_value, open_reader, pysam_header, pysam_sync_end,
};
use std::collections::HashMap;
use std::error::Error;
use std::fs;
use std::io::BufRead;
use std::path::Path;

const USAGE: &str = "usage: priva pedigree-filter -v VCF -p PED_TABLE -f TARGET_FAM -o OUTPUT_VCF";
const SOURCE_LINE: &str =
    "##INFO=<ID=VARIANT_SOURCE,Number=1,Type=String,Description=\"Source of the variant\">";
/// Parent IDs the Python script replaces with a missing value.
const MISSING_PARENT: [&str; 7] = ["0", ".", "", "NA", "nan", "NaN", "-"];

/// Allele indexes of one sample's GT; `None` is a missing allele (`.`).
type Genotype = Vec<Option<u32>>;

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args == ["--help"] || args == ["-h"] {
        println!("{USAGE}");
        return Ok(());
    }
    let mut options: HashMap<&str, String> = HashMap::new();
    let mut index = 0;
    while index < args.len() {
        let (option, inline) = match args[index].split_once('=') {
            Some((name, value)) if name.starts_with("--") => (name, Some(value.to_owned())),
            _ => (args[index].as_str(), None),
        };
        let name = match option {
            "-v" | "--vcf" => "vcf",
            "-p" | "--ped_table" => "ped_table",
            "-f" | "--target_fam" => "target_fam",
            "-o" | "--output" => "output",
            _ => return Err(invalid(format!("unknown option {option:?}\n{USAGE}"))),
        };
        let value = match inline {
            Some(value) => value,
            None => {
                index += 1;
                args.get(index)
                    .cloned()
                    .ok_or_else(|| invalid(format!("{option} needs a value\n{USAGE}")))?
            }
        };
        options.insert(name, value);
        index += 1;
    }
    let mut take = |name: &str| {
        options
            .remove(name)
            .ok_or_else(|| invalid(format!("--{name} is required\n{USAGE}")))
    };
    let (vcf, ped, family, output) = (
        take("vcf")?,
        take("ped_table")?,
        take("target_fam")?,
        take("output")?,
    );
    // Like the Python script, a failure after the arguments exits with code 1.
    if let Err(error) = pedigree_filter(
        Path::new(&vcf),
        Path::new(&ped),
        &family,
        Path::new(&output),
    ) {
        eprintln!("priva pedigree-filter: {error}");
        std::process::exit(1);
    }
    Ok(())
}

/// The family members the filter needs, chosen as the Python script does.
#[derive(Debug, PartialEq)]
struct Family {
    proband: String,
    proband_sex: String,
    patients: Vec<String>,
    controls: Vec<String>,
    father: Option<String>,
    mother: Option<String>,
}

/// Pedigree rows as pandas reads them: six fields, missing trailing ones `None`.
fn read_pedigree(text: &str) -> Result<Vec<[Option<String>; 6]>> {
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    let mut rows = Vec::new();
    for (number, line) in text.split('\n').enumerate() {
        let line = line.strip_suffix('\r').unwrap_or(line);
        if line.is_empty() {
            continue;
        }
        if line.contains('"') {
            return Err(invalid(format!(
                "pedigree line {}: quoted fields are outside the validated contract",
                number + 1
            )));
        }
        let mut row: [Option<String>; 6] = Default::default();
        for (cell, value) in row.iter_mut().zip(line.split('\t')) {
            *cell = Some(value.to_owned());
        }
        rows.push(row);
    }
    Ok(rows)
}

/// `None` when the pedigree has no row of the family.
fn family_members(rows: &[[Option<String>; 6]], target: &str) -> Result<Option<Family>> {
    let rows: Vec<_> = rows
        .iter()
        .filter(|row| row[0].as_deref() == Some(target))
        .collect();
    if rows.is_empty() {
        return Ok(None);
    }
    let with_phenotype = |phenotype: &str| -> Vec<String> {
        rows.iter()
            .filter(|row| row[5].as_deref() == Some(phenotype))
            .filter_map(|row| row[1].clone())
            .collect()
    };
    let patients = with_phenotype("2");
    let proband = patients
        .first()
        .cloned()
        .ok_or_else(|| invalid(format!("family {target} has no patient (phenotype 2)")))?;
    let proband_sex = rows
        .iter()
        .find(|row| row[1].as_deref() == Some(proband.as_str()))
        .map(|row| row[4].clone().unwrap_or_else(|| "nan".to_owned()))
        .unwrap_or_default();
    let parent = |column: usize| {
        rows.iter()
            .filter_map(|row| row[column].clone())
            .find(|id| !MISSING_PARENT.contains(&id.as_str()))
    };
    Ok(Some(Family {
        proband,
        proband_sex,
        patients,
        controls: with_phenotype("1"),
        father: parent(2),
        mother: parent(3),
    }))
}

/// A family member's sample column, or the name pysam does not find.
enum Sample {
    Column(usize),
    Absent(String),
}

/// The family as found in the VCF header, and whether the header defines END.
struct Roles {
    patients: Vec<Sample>,
    controls: Vec<Sample>,
    /// `None` is a member the script never looks up: no parent, or an empty ID.
    proband: Option<Sample>,
    father: Option<Sample>,
    mother: Option<Sample>,
    proband_sex: String,
    end_defined: bool,
}

fn pedigree_filter(vcf: &Path, ped: &Path, target: &str, output: &Path) -> Result<()> {
    let rows = read_pedigree(&fs::read_to_string(ped)?)?;
    let Some(family) = family_members(&rows, target)? else {
        eprintln!("ERROR: Family {target} not found in pedigree file.");
        return Ok(());
    };
    eprintln!(
        "INFO: The proband is {} and the proband sex is {}; patients {:?}, controls {:?}, father {:?}, mother {:?}",
        family.proband,
        family.proband_sex,
        family.patients,
        family.controls,
        family.father,
        family.mother
    );
    let mut reader = open_reader(vcf)?;
    let (mut buffer, mut number, mut header) = (Vec::new(), 0, Vec::new());
    let samples_line = loop {
        match next_line(reader.as_mut(), &mut buffer, &mut number)? {
            Some(text) if text.starts_with("#CHROM") => break text.to_owned(),
            Some(text) if text.starts_with('#') => header.push(text.to_owned()),
            Some(_) => return Err(invalid(format!("line {number}: a record before #CHROM"))),
            None => return Err(invalid("the VCF has no #CHROM line")),
        }
    };
    if header
        .iter()
        .any(|line: &String| line.starts_with("##INFO=<ID=VARIANT_SOURCE,"))
    {
        return Err(invalid(
            "the VCF header already defines INFO VARIANT_SOURCE, which pysam will not add twice",
        ));
    }
    header.push(SOURCE_LINE.to_owned());
    let columns: HashMap<&str, usize> = samples_line
        .split('\t')
        .skip(9)
        .enumerate()
        .map(|(index, sample)| (sample, index + 9))
        .collect();
    let sample = |id: &str| match columns.get(id) {
        Some(&column) => Sample::Column(column),
        None => Sample::Absent(id.to_owned()),
    };
    let unique = |ids: &[String]| {
        let mut seen: Vec<&String> = Vec::new();
        for id in ids {
            if !seen.contains(&id) {
                seen.push(id);
            }
        }
        seen.into_iter().map(|id| sample(id)).collect::<Vec<_>>()
    };
    // The script's _get_genotype skips an empty ID instead of looking it up.
    let person = |id: Option<&str>| id.filter(|id| !id.is_empty()).map(sample);
    let roles = Roles {
        patients: unique(&family.patients),
        controls: unique(&family.controls),
        proband: person(Some(family.proband.as_str())),
        father: person(family.father.as_deref()),
        mother: person(family.mother.as_deref()),
        proband_sex: family.proband_sex.clone(),
        end_defined: header
            .iter()
            .any(|line| line.starts_with("##INFO=<ID=END,")),
    };
    header.push(samples_line.clone());
    let mut writer = create_writer(output)?;
    for line in pysam_header(&header) {
        writer.write_line(line.as_bytes())?;
        writer.write_line(b"\n")?;
    }
    // Like pysam's `with` block, close the output even when a record stops the run.
    let outcome = filter_records(reader.as_mut(), writer.as_mut(), &roles, number);
    let closed = writer.finish();
    outcome?;
    closed?;
    eprintln!("INFO: Filtered VCF written to {}", output.display());
    Ok(())
}

/// The next line without its line ending, or `None` at the end of the input.
fn next_line<'a>(
    reader: &mut dyn BufRead,
    buffer: &'a mut Vec<u8>,
    number: &mut usize,
) -> Result<Option<&'a str>> {
    buffer.clear();
    if reader.read_until(b'\n', buffer)? == 0 {
        return Ok(None);
    }
    *number += 1;
    if buffer.last() == Some(&b'\n') {
        buffer.pop();
    }
    let buffer: &'a Vec<u8> = buffer;
    let line = std::str::from_utf8(buffer)
        .map_err(|_| invalid(format!("line {number}: not UTF-8 text")))?;
    Ok(Some(line.strip_suffix('\r').unwrap_or(line)))
}

fn absent(line: usize, id: &str) -> Box<dyn Error> {
    invalid(format!(
        "line {line}: pedigree member {id:?} is not a sample in the VCF (pysam: invalid sample name)"
    ))
}

fn filter_records(
    reader: &mut dyn BufRead,
    writer: &mut dyn TextSink,
    roles: &Roles,
    mut number: usize,
) -> Result<()> {
    let mut buffer = Vec::new();
    while let Some(text) = next_line(reader, &mut buffer, &mut number)? {
        if text.is_empty() {
            continue;
        }
        if text.starts_with('#') {
            return Err(invalid(format!(
                "line {number}: a header line after #CHROM"
            )));
        }
        let mut fields: Vec<&str> = text.split('\t').collect();
        if fields.len() < 8 {
            return Err(invalid(format!("line {number}: fewer than 8 VCF columns")));
        }
        let line = number;
        let alleles_in_record = match fields[4] {
            "." => 1,
            alts => 1 + alts.split(',').count() as u32,
        };
        let gt_index = fields
            .get(8)
            .and_then(|format| format.split(':').position(|key| key == "GT"));
        let genotype = |column: usize| -> Result<Genotype> {
            let index = gt_index.ok_or_else(|| {
                invalid(format!(
                    "line {line}: a family sample has no GT in this record"
                ))
            })?;
            let value = fields
                .get(column)
                .and_then(|sample| sample.split(':').nth(index))
                .unwrap_or(".");
            parse_gt(value, alleles_in_record)
                .map_err(|error| invalid(format!("line {line}: {error}")))
        };
        // record.samples[sample]['GT']: an absent sample fails, then a missing GT.
        let required = |sample: &Sample| -> Result<Genotype> {
            match sample {
                Sample::Column(column) => genotype(*column),
                Sample::Absent(id) => Err(absent(line, id)),
            }
        };
        // _get_genotype reads GT with .get: a missing GT gives no genotype.
        let optional = |sample: &Option<Sample>| -> Result<Option<Genotype>> {
            match sample {
                None => Ok(None),
                Some(Sample::Absent(id)) => Err(absent(line, id)),
                Some(Sample::Column(column)) if gt_index.is_some() => genotype(*column).map(Some),
                Some(Sample::Column(_)) => Ok(None),
            }
        };
        let mut carried = false;
        for patient in &roles.patients {
            if is_variant_carrier(&required(patient)?) {
                carried = true;
                break;
            }
        }
        if !carried {
            continue;
        }
        let mut control_homozygous = false;
        for control in &roles.controls {
            if is_homozygous_or_hemizygous(&required(control)?) {
                control_homozygous = true;
                break;
            }
        }
        if control_homozygous {
            continue;
        }
        // The parents are looked up only once the proband carries allele 1.
        let proband = optional(&roles.proband)?;
        let source = if has_variant(proband.as_ref()) {
            let father = optional(&roles.father)?;
            let mother = optional(&roles.mother)?;
            variant_source(
                fields[0],
                &roles.proband_sex,
                proband.as_ref(),
                father.as_ref(),
                mother.as_ref(),
            )
        } else {
            "UNKNOWN"
        };
        let pos: i64 = fields[1]
            .parse()
            .map_err(|_| invalid(format!("line {line}: POS is not an integer")))?;
        let alts: Vec<&str> = fields[4].split(',').collect();
        let info = pysam_sync_end(
            pos,
            fields[3],
            &alts,
            &set_info(fields[7], "VARIANT_SOURCE", source),
        );
        if !roles.end_defined && info_value(&info, "END").is_some() {
            // pysam fails here too (OSError) after writing the header.
            return Err(invalid(format!(
                "line {line}: the record needs INFO END, which the header does not define"
            )));
        }
        fields[7] = &info;
        writer.write_line(fields.join("\t").as_bytes())?;
        writer.write_line(b"\n")?;
    }
    Ok(())
}

/// GT text as pysam's allele tuple: an allele beyond the record's alleles is
/// missing, like `.`.
fn parse_gt(value: &str, alleles_in_record: u32) -> Result<Genotype> {
    if value.is_empty() || value == "." {
        return Ok(vec![None]);
    }
    value
        .split(['/', '|'])
        .map(|allele| match allele {
            "." => Ok(None),
            _ => allele
                .parse::<u32>()
                .map(|allele| (allele < alleles_in_record).then_some(allele))
                .map_err(|_| invalid(format!("GT {value:?} is not a genotype"))),
        })
        .collect()
}

fn is_variant_carrier(genotype: &[Option<u32>]) -> bool {
    genotype.contains(&Some(1)) || genotype.contains(&Some(2))
}

fn is_homozygous_or_hemizygous(genotype: &[Option<u32>]) -> bool {
    match genotype.first() {
        Some(Some(first)) => *first != 0 && genotype.iter().all(|allele| *allele == Some(*first)),
        _ => false,
    }
}

fn genotype_exists(genotype: Option<&Genotype>) -> bool {
    genotype.is_some_and(|genotype| genotype.iter().all(Option::is_some))
}

fn has_variant(genotype: Option<&Genotype>) -> bool {
    genotype_exists(genotype) && genotype.is_some_and(|genotype| genotype.contains(&Some(1)))
}

fn variant_source(
    chrom: &str,
    proband_sex: &str,
    proband: Option<&Genotype>,
    father: Option<&Genotype>,
    mother: Option<&Genotype>,
) -> &'static str {
    if !has_variant(proband) {
        return "UNKNOWN";
    }
    let autosomal = || {
        if !genotype_exists(father) || !genotype_exists(mother) {
            return if has_variant(father) || has_variant(mother) {
                "INHERITED"
            } else {
                "UNKNOWN"
            };
        }
        match (has_variant(father), has_variant(mother)) {
            (false, false) => "DENOVO",
            (true, false) => "PATERNAL",
            (false, true) => "MATERNAL",
            (true, true) => "INHERITED",
        }
    };
    match chrom.replace("chr", "").to_uppercase().as_str() {
        "X" => match proband_sex {
            "1" if !genotype_exists(mother) => "UNKNOWN",
            "1" if has_variant(mother) => "MATERNAL",
            "1" => "DENOVO",
            "2" => autosomal(),
            _ => "UNKNOWN",
        },
        "Y" => match proband_sex {
            "1" if !genotype_exists(father) => "UNKNOWN",
            "1" if has_variant(father) => "PATERNAL",
            "1" => "DENOVO",
            _ => "UNKNOWN",
        },
        _ => autosomal(),
    }
}

/// INFO after pysam's `record.info[key] = value`: an existing key is updated
/// in place, a new one is appended.
fn set_info(info: &str, key: &str, value: &str) -> String {
    let item = format!("{key}={value}");
    if info == "." {
        return item;
    }
    let mut placed = false;
    let mut items: Vec<&str> = info
        .split(';')
        .map(|existing| {
            let name = existing.split_once('=').map_or(existing, |(name, _)| name);
            if name == key && !placed {
                placed = true;
                item.as_str()
            } else {
                existing
            }
        })
        .collect();
    if !placed {
        items.push(&item);
    }
    items.join(";")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gt(text: &str) -> Genotype {
        parse_gt(text, 4).unwrap()
    }

    #[test]
    fn genotypes_parse_like_pysam_allele_tuples() {
        assert_eq!(gt("0/1"), vec![Some(0), Some(1)]);
        assert_eq!(gt("1|2"), vec![Some(1), Some(2)]);
        assert_eq!(gt("./."), vec![None, None]);
        assert_eq!(gt("."), vec![None]);
        assert_eq!(gt("1"), vec![Some(1)]);
        assert_eq!(gt("0/."), vec![Some(0), None]);
        assert!(parse_gt("A/1", 4).is_err());
        assert_eq!(parse_gt("0/1", 1).unwrap(), vec![Some(0), None]);
    }

    #[test]
    fn carriers_count_alleles_one_and_two_only() {
        assert!(is_variant_carrier(&gt("0/1")));
        assert!(is_variant_carrier(&gt("0/2")));
        assert!(!is_variant_carrier(&gt("0/3")));
        assert!(!is_variant_carrier(&gt("./.")));
        assert!(is_homozygous_or_hemizygous(&gt("1/1")));
        assert!(is_homozygous_or_hemizygous(&gt("2|2")));
        assert!(is_homozygous_or_hemizygous(&gt("1")));
        assert!(is_homozygous_or_hemizygous(&gt("3/3")));
        assert!(!is_homozygous_or_hemizygous(&gt("0/0")));
        assert!(!is_homozygous_or_hemizygous(&gt("./.")));
        assert!(!is_homozygous_or_hemizygous(&gt("1/2")));
    }

    #[test]
    fn variant_sources_follow_the_python_rules() {
        let (het, ref_, miss, alt2) = (gt("0/1"), gt("0/0"), gt("./."), gt("0/2"));
        let source = variant_source;
        assert_eq!(
            source("chr1", "1", Some(&het), Some(&ref_), Some(&ref_)),
            "DENOVO"
        );
        assert_eq!(
            source("chr1", "1", Some(&het), Some(&het), Some(&ref_)),
            "PATERNAL"
        );
        assert_eq!(
            source("chr1", "1", Some(&het), Some(&ref_), Some(&het)),
            "MATERNAL"
        );
        assert_eq!(
            source("chr1", "1", Some(&het), Some(&het), Some(&het)),
            "INHERITED"
        );
        assert_eq!(
            source("chr1", "1", Some(&het), Some(&het), Some(&miss)),
            "INHERITED"
        );
        assert_eq!(
            source("chr1", "1", Some(&het), None, Some(&ref_)),
            "UNKNOWN"
        );
        assert_eq!(
            source("chr1", "1", Some(&alt2), Some(&ref_), Some(&ref_)),
            "UNKNOWN"
        );
        assert_eq!(source("chrX", "1", Some(&het), None, Some(&ref_)), "DENOVO");
        assert_eq!(source("X", "1", Some(&het), None, Some(&het)), "MATERNAL");
        assert_eq!(source("chrX", "1", Some(&het), None, None), "UNKNOWN");
        assert_eq!(
            source("chrX", "2", Some(&het), Some(&ref_), Some(&ref_)),
            "DENOVO"
        );
        assert_eq!(
            source("chrX", "nan", Some(&het), Some(&ref_), Some(&ref_)),
            "UNKNOWN"
        );
        assert_eq!(
            source("chrY", "1", Some(&gt("1")), Some(&gt("1")), None),
            "PATERNAL"
        );
        assert_eq!(
            source("chrY", "1", Some(&gt("1")), Some(&gt("0")), None),
            "DENOVO"
        );
        assert_eq!(
            source("chrY", "2", Some(&gt("1")), Some(&gt("1")), None),
            "UNKNOWN"
        );
    }

    #[test]
    fn family_members_are_chosen_like_the_python_script() {
        let rows = read_pedigree(
            "#FamilyID\tIndividualID\tPaternalID\tMaternalID\tSex\tPhenotype\r\n\
             F1\tsib\t-9\t0\t2\t1\n\
             F1\tkid\tdad\tmom\t1\t2\textra\n\
             \n\
             F1\tdad\t0\t0\t1\t1\n\
             F1\tmom\tNA\t.\t2\t2\n\
             F2\tother\n",
        )
        .unwrap();
        assert_eq!(
            family_members(&rows, "F1").unwrap().unwrap(),
            Family {
                proband: "kid".into(),
                proband_sex: "1".into(),
                patients: vec!["kid".into(), "mom".into()],
                controls: vec!["sib".into(), "dad".into()],
                father: Some("-9".into()),
                mother: Some("mom".into()),
            }
        );
        assert_eq!(family_members(&rows, "F3").unwrap(), None);
        assert!(family_members(&rows, "F2").is_err());
    }

    #[test]
    fn info_updates_in_place_or_appends() {
        assert_eq!(
            set_info(".", "VARIANT_SOURCE", "DENOVO"),
            "VARIANT_SOURCE=DENOVO"
        );
        assert_eq!(
            set_info("DP=3", "VARIANT_SOURCE", "DENOVO"),
            "DP=3;VARIANT_SOURCE=DENOVO"
        );
        assert_eq!(
            set_info("VARIANT_SOURCE=UNKNOWN;DP=3", "VARIANT_SOURCE", "DENOVO"),
            "VARIANT_SOURCE=DENOVO;DP=3"
        );
    }
}
