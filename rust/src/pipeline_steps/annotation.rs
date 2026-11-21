//! What the annotation steps share from `common_bash_utils.sh`: telling the
//! assembly of a VCF or FASTA, and renaming chromosomes with bcftools.
use super::tools::{Context, Result, args, is_file, log};
use std::path::Path;

/// `check_vcf_contig_size VCF chr1`: the lengths of the `##contig` lines
/// whose text holds `ID=chr1` (chr10 to chr19 too), each the field after a
/// `length` field when the line is split at `,` and `=` (so a length that
/// ends the line keeps its `>`), one per line.
fn chr1_sizes(header: &str) -> String {
    let mut sizes = Vec::new();
    for line in header.lines() {
        if !line.starts_with("##contig") || !line.contains("ID=chr1") {
            continue;
        }
        let fields: Vec<&str> = line.split([',', '=']).collect();
        for (index, field) in fields.iter().enumerate() {
            if *field == "length" {
                if let Some(next) = fields.get(index + 1) {
                    sizes.push(*next);
                }
            }
        }
    }
    sizes.join("\n")
}

/// The `##contig` header lines whose third `=` field matches NAME.
fn has_contig(header: &str, name: &str) -> bool {
    header.lines().any(|line| {
        let fields: Vec<&str> = line.split('=').collect();
        fields[0].contains("##contig") && fields.get(2).is_some_and(|f| f.contains(name))
    })
}

/// `check_vcf_assembly_version VCF`: hg19 or hg38 from the chr1 length in
/// the header, else from a contig only one of them has; empty when it
/// cannot tell or the VCF is not valid.
pub(crate) fn vcf_assembly(context: &Context, vcf: &str) -> String {
    if !context.vcf_valid(vcf, 1, None) {
        log(&format!(
            "Cant locate the required contig name in the VCF header line to tell the assembly version of this VCF file {vcf}"
        ));
        return String::new();
    }
    let header = context.header(vcf);
    match chr1_sizes(&header).as_str() {
        "249250621" => "hg19".into(),
        "248956422" => "hg38".into(),
        _ if has_contig(&header, "chr7_gl000195_random") => "hg19".into(),
        _ if has_contig(&header, "chr6_GL000256v2_alt") => "hg38".into(),
        _ => String::new(),
    }
}

/// `check_fasta_assembly_version FASTA`: hg19 or hg38 from the chr1 (or 1)
/// length in the `.fai`, else from `hg19` or `hg38` in the file name.
pub(crate) fn fasta_assembly(fasta: &str) -> String {
    let index = format!("{fasta}.fai");
    let mut size = String::new();
    if is_file(&index) {
        if let Ok(text) = std::fs::read_to_string(&index) {
            for line in text.lines() {
                let mut fields = line.split('\t');
                if matches!(fields.next(), Some("chr1" | "1")) {
                    size = fields.next().unwrap_or("").to_owned();
                    break;
                }
            }
        }
    }
    match size.as_str() {
        "249250621" => return "hg19".into(),
        "248956422" => return "hg38".into(),
        _ => {}
    }
    let name = Path::new(fasta)
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let first = ["hg19", "hg38"]
        .iter()
        .filter_map(|tag| name.find(tag).map(|at| (at, *tag)))
        .min();
    first.map(|(_, tag)| tag.to_owned()).unwrap_or_default()
}

/// `liftover_from_GRCh_to_ucsc` and `liftover_from_ucsc_to_GRCh`:
/// `bcftools annotate --rename-chrs MAP VCF -Ou | bcftools sort -Oz -o OUT`,
/// then `tabix`. Whether it worked.
pub(crate) fn rename_chromosomes(
    context: &Context,
    vcf: &str,
    map: &str,
    output: &str,
) -> Result<bool> {
    let renamed = context.pipeline(
        &[
            (
                "bcftools",
                args(&["annotate", "--rename-chrs", map, vcf, "-Ou"]),
            ),
            ("bcftools", args(&["sort", "-Oz", "-o", output])),
        ],
        false,
    )?;
    Ok(renamed && context.succeeds("tabix", &args(&["-f", "-p", "vcf", output])))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chr1_length_is_read_as_the_shell_read_it() {
        let header = "##contig=<ID=chr1,length=248956422>\n##contig=<ID=chr10,length=133797422>\n";
        assert_eq!(chr1_sizes(header), "248956422>\n133797422>");
        let gatk = "##contig=<ID=chr1,length=249250621,assembly=hg19>\n";
        assert_eq!(chr1_sizes(gatk), "249250621");
        assert!(has_contig(
            "##contig=<ID=chr6_GL000256v2_alt,length=4929269>",
            "chr6_GL000256v2_alt"
        ));
    }

    #[test]
    fn fasta_assembly_falls_back_to_the_name() {
        assert_eq!(fasta_assembly("/missing/ucsc.hg19.fasta"), "hg19");
        assert_eq!(fasta_assembly("/missing/GRCh38.fa"), "");
    }
}
