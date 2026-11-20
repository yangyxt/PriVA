//! Rust version of `scripts/generate_varID_for_vcf.py`.
//!
//! encode: write `chrom:pos:ref:alt` into each record's ID column, or
//! `chrom:pos-end:ref:alt` for a record with a symbolic `<...>` allele. The
//! Python script also asks for INFO SVTYPE plus END, but pysam never reports
//! END as an INFO key, so in practice only symbolic alleles count; their end is
//! pysam's `record.stop` (INFO END when present, else the end of REF). As in
//! the Python script, a multi-allelic record keeps the ID of its last
//! alternative allele.
//! decode: rewrite such an ID as `chrom:pos-end:ref->alt` (`end` is `None`
//! when the ID has no end coordinate).
//! Both modes write records the way pysam does: the PASS filter line second in
//! the header and INFO END synchronised (see `vcf_text`); like pysam, a record
//! that needs END stops the run when the header does not define END. Other columns are
//! copied unchanged, which equals htslib's own formatting for bcftools-written
//! input such as PriVA's normalised VCFs.
use crate::cli::{Result, invalid};
use crate::vcf_text::{create_writer, info_value, open_reader, pysam_header, pysam_sync_end};
use std::io::BufRead;
use std::path::Path;

const USAGE: &str = "usage: priva variant-ids -m {encode,decode} -i INPUT_VCF [-o OUTPUT_VCF]";

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args == ["--help"] || args == ["-h"] {
        println!("{USAGE}");
        return Ok(());
    }
    let (mut mode, mut input, mut output) = (None, None, String::new());
    let mut index = 0;
    while index < args.len() {
        let (option, inline) = match args[index].split_once('=') {
            Some((name, value)) if name.starts_with("--") => (name, Some(value.to_owned())),
            _ => (args[index].as_str(), None),
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
        match option {
            "-m" | "--mode" => mode = Some(value),
            "-i" | "--input_vcf" => input = Some(value),
            "-o" | "--output_vcf" => output = value,
            _ => return Err(invalid(format!("unknown option {option:?}\n{USAGE}"))),
        }
        index += 1;
    }
    let mode = mode.ok_or_else(|| invalid(format!("-m/--mode is required\n{USAGE}")))?;
    let input = input.ok_or_else(|| invalid(format!("-i/--input_vcf is required\n{USAGE}")))?;
    let (suffix, rewrite): (&str, fn(&str, &[&str]) -> Result<String>) = match mode.as_str() {
        "encode" => (".uniqID.vcf", encode_id),
        "decode" => (".decodeID.vcf", decode_id),
        _ => return Err(invalid(format!("mode must be encode or decode\n{USAGE}"))),
    };
    // The Python script derives the default name with str.replace (every ".vcf").
    let output = if output.is_empty() {
        input.replace(".vcf", suffix)
    } else {
        output
    };
    rewrite_ids(Path::new(&input), Path::new(&output), rewrite)
}

fn rewrite_ids(
    input: &Path,
    output: &Path,
    rewrite: fn(&str, &[&str]) -> Result<String>,
) -> Result<()> {
    let mut reader = open_reader(input)?;
    let mut writer = create_writer(output)?;
    let (mut buffer, mut number, mut header) = (Vec::new(), 0, Vec::new());
    let mut end_defined = false;
    loop {
        buffer.clear();
        if reader.read_until(b'\n', &mut buffer)? == 0 {
            break;
        }
        number += 1;
        if buffer.last() == Some(&b'\n') {
            buffer.pop();
        }
        let line = std::str::from_utf8(&buffer)
            .map_err(|_| invalid(format!("line {number}: not UTF-8 text")))?;
        let text = line.strip_suffix('\r').unwrap_or(line);
        if text.starts_with('#') {
            header.push(text.to_owned());
            if text.starts_with("#CHROM") {
                end_defined = header
                    .iter()
                    .any(|line| line.starts_with("##INFO=<ID=END,"));
                // pysam writes the header as htslib holds it.
                for line in pysam_header(&header) {
                    writer.write_line(line.as_bytes())?;
                    writer.write_line(b"\n")?;
                }
            }
            continue;
        }
        if text.is_empty() {
            continue;
        }
        let mut fields: Vec<&str> = text.split('\t').collect();
        if fields.len() < 8 {
            return Err(invalid(format!("line {number}: fewer than 8 VCF columns")));
        }
        let id =
            rewrite(text, &fields).map_err(|error| invalid(format!("line {number}: {error}")))?;
        let pos: i64 = fields[1]
            .parse()
            .map_err(|_| invalid(format!("line {number}: POS is not an integer")))?;
        let alts: Vec<&str> = fields[4].split(',').collect();
        let info = pysam_sync_end(pos, fields[3], &alts, fields[7]);
        if !end_defined && info_value(&info, "END").is_some() {
            // pysam fails here too (OSError) after writing the header.
            return Err(invalid(format!(
                "line {number}: the record needs INFO END, which the header does not define"
            )));
        }
        fields[2] = &id;
        fields[7] = &info;
        writer.write_line(fields.join("\t").as_bytes())?;
        writer.write_line(b"\n")?;
    }
    writer.finish()?;
    Ok(())
}

fn encode_id(_line: &str, fields: &[&str]) -> Result<String> {
    let (chrom, pos, reference, alts, info) =
        (fields[0], fields[1], fields[3], fields[4], fields[7]);
    // pysam reports a missing ALT as None, and iterating it fails.
    if alts == "." {
        return Err(invalid("record has no alternative allele"));
    }
    let alts: Vec<&str> = alts.split(',').collect();
    let pos: i64 = pos.parse().map_err(|_| invalid("POS is not an integer"))?;
    let structural = alts
        .iter()
        .any(|alt| alt.starts_with('<') && alt.ends_with('>'));
    let end = if structural {
        Some(match info_value(info, "END") {
            Some(Some(value)) => value
                .parse::<i64>()
                .map_err(|_| invalid("INFO END is not an integer"))?,
            _ => pos + reference.len() as i64 - 1,
        })
    } else {
        None
    };
    let alt = alts[alts.len() - 1];
    Ok(match end {
        Some(end) => format!("{chrom}:{pos}-{end}:{reference}:{alt}"),
        None => format!("{chrom}:{pos}:{reference}:{alt}"),
    })
}

fn decode_id(_line: &str, fields: &[&str]) -> Result<String> {
    if fields[2] == "." {
        return Err(invalid("record has no ID to decode"));
    }
    let parts: Vec<&str> = fields[2].split(':').collect();
    if parts.len() < 4 {
        return Err(invalid("ID is not chrom:pos:ref:alt"));
    }
    let integer = |text: &str| {
        text.trim()
            .parse::<i64>()
            .map_err(|_| invalid("ID position is not an integer"))
    };
    let (pos, end) = match parts[1].split_once('-') {
        Some((pos, end)) if !end.contains('-') => (integer(pos)?, integer(end)?.to_string()),
        Some(_) => return Err(invalid("ID position has more than one '-'")),
        None => (integer(parts[1])?, "None".to_owned()),
    };
    Ok(format!(
        "{}:{pos}-{end}:{}->{}",
        parts[0], parts[2], parts[3]
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fields(line: &str) -> Vec<&str> {
        line.split('\t').collect()
    }

    #[test]
    fn encodes_small_variants_with_the_last_alternative_allele() {
        let line = "chr1\t100\t.\tA\tG,T\t50\tPASS\tDP=3";
        assert_eq!(encode_id(line, &fields(line)).unwrap(), "chr1:100:A:T");
    }

    #[test]
    fn encodes_structural_variants_with_their_end() {
        let line = "chr2\t10\t.\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=500";
        assert_eq!(
            encode_id(line, &fields(line)).unwrap(),
            "chr2:10-500:N:<DEL>"
        );
        let line = "chr2\t10\t.\tACG\t<INV>\t.\t.\t.";
        assert_eq!(
            encode_id(line, &fields(line)).unwrap(),
            "chr2:10-12:ACG:<INV>"
        );
        // Sequence alleles never count, even with SVTYPE and END (pysam hides END).
        let line = "chr2\t10\t.\tACGT\tA\t.\t.\tSVTYPE=DEL;END=13";
        assert_eq!(encode_id(line, &fields(line)).unwrap(), "chr2:10:ACGT:A");
    }

    #[test]
    fn decodes_both_id_forms() {
        let line = "chr1\t100\tchr1:100:A:T\tA\tT\t.\t.\t.";
        assert_eq!(
            decode_id(line, &fields(line)).unwrap(),
            "chr1:100-None:A->T"
        );
        let line = "chr2\t10\tchr2:10-500:N:<DEL>\tN\t<DEL>\t.\t.\t.";
        assert_eq!(
            decode_id(line, &fields(line)).unwrap(),
            "chr2:10-500:N-><DEL>"
        );
    }

    #[test]
    fn missing_alleles_and_ids_fail_like_the_python_script() {
        let line = "chr1\t100\t.\tA\t.\t.\t.\t.";
        assert!(encode_id(line, &fields(line)).is_err());
        assert!(decode_id(line, &fields(line)).is_err());
    }
}
