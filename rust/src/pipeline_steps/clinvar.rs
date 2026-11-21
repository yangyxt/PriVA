//! The ClinVar step without `anno_clinvar_data` of annotation_vcf.sh:
//! `bcftools annotate -a clinvar_vcf -c ... | bcftools sort`, then `tabix`,
//! and the VCF replaced by its annotated copy.
use super::tools::{
    Context, Result, args, log, move_with_index, option, random_tag, refuse, replace_first, words,
};

const TAGS: &str = "CLNDN,CLNHGVS,CLNREVSTAT,CLNSIG,GENEINFO,CLNCSQ";
const COLUMNS: &str =
    ".INFO/CLNDN,.INFO/CLNHGVS,.INFO/CLNREVSTAT,.INFO/CLNSIG,.INFO/GENEINFO,.INFO/CLNCSQ:=INFO/CSQ";

/// `anno_clinvar_data VCF clinvar_vcf threads`.
pub(crate) fn run(context: &Context, input: &str) -> Result<()> {
    let clinvar = context.value("clinvar_vcf");
    let configured = context.value("threads");
    let threads = if configured.is_empty() {
        "1"
    } else {
        configured.as_str()
    };
    let output = replace_first(input, ".vcf", &format!(".{}.vcf", random_tag()));
    if !input.ends_with(".vcf.gz") {
        return refuse(format!(
            "The input vcf {input} is not a gzipped vcf file. Quit with Error!"
        ));
    }
    if context.has_info_tags(input, TAGS) {
        log(&format!(
            "The input vcf {input} already contains the INFO tags {TAGS}. We do not need to add them again"
        ));
        return Ok(());
    }
    let mut annotate = vec!["annotate".to_owned()];
    annotate.extend(option("--threads", threads));
    annotate.extend(option("-a", &clinvar));
    annotate.extend(args(&["-c", COLUMNS, "-Ou"]));
    annotate.extend(words(input));
    let done = context.pipeline(
        &[
            ("bcftools", annotate),
            ("bcftools", args(&["sort", "-Oz", "-o", &output])),
        ],
        false,
    )? && context.succeeds("tabix", &args(&["-f", "-p", "vcf", &output]))
        && move_with_index(&output, input).is_ok()
        && context.has_info_tags(input, TAGS);
    if !done {
        return refuse(format!(
            "Failed to add ClinVar annotation on {input}. Quit now"
        ));
    }
    Ok(())
}
