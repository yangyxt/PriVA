//! Joins the INFO fields that separate annotation steps added to copies of
//! the same preprocessed VCF, as if the steps had run one after another on
//! one file: each step's new header lines and INFO fields follow the
//! preprocessed ones, in step order.
//!
//! Every step must keep every record, in order, and change nothing but INFO.
//! The merge refuses a record whose position, alleles, QUAL, FILTER or
//! genotype columns differ between the files, and a step that changed an
//! INFO field the preprocessed VCF already had.
use super::{refuse, Result};
use std::collections::HashSet;
use std::fs;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::Command;

type Lines = io::Lines<Box<dyn BufRead>>;

fn open(path: &Path) -> Result<Lines> {
    let file = fs::File::open(path).map_err(|e| crate::python_path::os_error(e, path, None))?;
    let reader: Box<dyn BufRead> = Box::new(BufReader::with_capacity(
        1 << 20,
        flate2::read::MultiGzDecoder::new(file),
    ));
    Ok(reader.lines())
}

/// The `##` lines and the `#CHROM` line.
fn header(lines: &mut Lines, name: &str) -> Result<(Vec<String>, String)> {
    let mut meta = Vec::new();
    for line in lines.by_ref() {
        let line = line?;
        if line.starts_with("##") {
            meta.push(line);
        } else if line.starts_with("#CHROM") {
            return Ok((meta, line));
        } else {
            break;
        }
    }
    refuse(format!("{name} VCF has no #CHROM header line"))
}

fn info_entries(field: &str) -> Vec<&str> {
    if field == "." {
        Vec::new()
    } else {
        field.split(';').collect()
    }
}

fn key(entry: &str) -> &str {
    entry.split_once('=').map_or(entry, |(key, _)| key)
}

/// One merged record line from the preprocessed record and each step's.
fn merge_record(base: &str, steps: &[(&str, String)], number: usize) -> Result<String> {
    let fields: Vec<&str> = base.split('\t').collect();
    if fields.len() < 8 {
        return refuse(format!("Preprocessed VCF record {number} has fewer than 8 columns"));
    }
    let mut info: Vec<String> = info_entries(fields[7]).iter().map(|e| e.to_string()).collect();
    let base_info: HashSet<String> = info.iter().cloned().collect();
    let base_keys: HashSet<String> = info.iter().map(|e| key(e).to_owned()).collect();
    let mut added: HashSet<String> = HashSet::new();
    for (name, line) in steps {
        let other: Vec<&str> = line.split('\t').collect();
        let same_outside_info = other.len() == fields.len()
            && (0..7).all(|i| other[i] == fields[i])
            && (8..fields.len()).all(|i| other[i] == fields[i]);
        if !same_outside_info {
            return refuse(format!(
                "The {name} step changed record {number} outside INFO, or kept a different record"
            ));
        }
        for entry in info_entries(other[7]) {
            let name_of = key(entry);
            if base_keys.contains(name_of) {
                if !base_info.contains(entry) {
                    return refuse(format!(
                        "The {name} step changed the existing INFO field {name_of} at record {number}"
                    ));
                }
                continue;
            }
            if !added.insert(name_of.to_owned()) {
                return refuse(format!(
                    "Two annotation steps added the INFO field {name_of} at record {number}"
                ));
            }
            info.push(entry.to_owned());
        }
    }
    let joined = if info.is_empty() {
        ".".to_owned()
    } else {
        info.join(";")
    };
    let mut out: Vec<&str> = fields.clone();
    out[7] = &joined;
    Ok(out.join("\t"))
}

/// Writes `output` (BGZF with a tabix index) from the preprocessed VCF and
/// the step VCFs, in step order. `output` may be the preprocessed file.
pub(crate) fn merge(base: &Path, steps: &[(&str, &Path)], output: &Path) -> Result<()> {
    let mut base_lines = open(base)?;
    let (base_meta, chrom) = header(&mut base_lines, "Preprocessed")?;
    let mut meta = base_meta.clone();
    let mut seen: HashSet<String> = base_meta.into_iter().collect();
    let mut readers = Vec::new();
    for (name, path) in steps {
        let mut lines = open(path)?;
        let (step_meta, step_chrom) = header(&mut lines, name)?;
        if step_chrom != chrom {
            return refuse(format!("The {name} step VCF has different sample columns"));
        }
        for line in step_meta {
            if seen.insert(line.clone()) {
                meta.push(line);
            }
        }
        readers.push((*name, lines));
    }
    let plain = output.with_extension("merge.vcf");
    {
        let file =
            fs::File::create(&plain).map_err(|e| crate::python_path::os_error(e, &plain, None))?;
        let mut writer = BufWriter::with_capacity(1 << 20, file);
        for line in &meta {
            writeln!(writer, "{line}")?;
        }
        writeln!(writer, "{chrom}")?;
        let mut number = 0usize;
        for line in base_lines {
            let line = line?;
            number += 1;
            let mut step_lines = Vec::with_capacity(readers.len());
            for (name, lines) in readers.iter_mut() {
                match lines.next() {
                    Some(step) => step_lines.push((*name, step?)),
                    None => {
                        return refuse(format!(
                            "The {name} step VCF ended before record {number}"
                        ));
                    }
                }
            }
            writeln!(writer, "{}", merge_record(&line, &step_lines, number)?)?;
        }
        for (name, lines) in readers.iter_mut() {
            if lines.next().is_some() {
                return refuse(format!("The {name} step VCF has more records than the preprocessed VCF"));
            }
        }
        writer.flush()?;
    }
    let compressed = plain.with_extension("vcf.gz");
    let status = Command::new("bgzip").arg("-f").arg(&plain).status()?;
    if !status.success() {
        return refuse("bgzip failed on the merged VCF");
    }
    fs::rename(&compressed, output)
        .map_err(|e| crate::python_path::os_error(e, &compressed, Some(output)))?;
    let status = Command::new("tabix").args(["-f", "-p", "vcf"]).arg(output).status()?;
    if !status.success() {
        return refuse("tabix failed on the merged VCF");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn steps_add_their_fields_in_order_after_the_preprocessed_ones() {
        let base = "chr1\t10\t.\tA\tG\t50\tPASS\tDP=3\tGT\t0/1";
        let steps = [
            ("gnomAD", "chr1\t10\t.\tA\tG\t50\tPASS\tDP=3;AF_joint=0.1\tGT\t0/1".to_owned()),
            ("ClinVar", "chr1\t10\t.\tA\tG\t50\tPASS\tDP=3;CLNSIG=Benign\tGT\t0/1".to_owned()),
            ("VEP", "chr1\t10\t.\tA\tG\t50\tPASS\tDP=3;CSQ=G|x\tGT\t0/1".to_owned()),
        ];
        assert_eq!(
            merge_record(base, &steps, 1).unwrap(),
            "chr1\t10\t.\tA\tG\t50\tPASS\tDP=3;AF_joint=0.1;CLNSIG=Benign;CSQ=G|x\tGT\t0/1"
        );
        let empty = "chr1\t10\t.\tA\tG\t50\tPASS\t.";
        let one = [("VEP", "chr1\t10\t.\tA\tG\t50\tPASS\tCSQ=G".to_owned())];
        assert_eq!(merge_record(empty, &one, 1).unwrap(), "chr1\t10\t.\tA\tG\t50\tPASS\tCSQ=G");
    }

    #[test]
    fn changes_outside_new_info_fields_are_refused() {
        let base = "chr1\t10\t.\tA\tG\t50\tPASS\tDP=3";
        let message = |line: &str| match merge_record(base, &[("step", line.to_owned())], 4) {
            Err(error) => error.to_string(),
            Ok(_) => String::new(),
        };
        assert!(message("chr1\t11\t.\tA\tG\t50\tPASS\tDP=3").contains("outside INFO"));
        assert!(message("chr1\t10\t.\tA\tG\t50\tLowQ\tDP=3").contains("outside INFO"));
        assert!(message("chr1\t10\t.\tA\tG\t50\tPASS\tDP=4").contains("existing INFO field DP"));
        let twice = [
            ("one", "chr1\t10\t.\tA\tG\t50\tPASS\tDP=3;X=1".to_owned()),
            ("two", "chr1\t10\t.\tA\tG\t50\tPASS\tDP=3;X=2".to_owned()),
        ];
        assert!(merge_record(base, &twice, 1).is_err());
    }
}
