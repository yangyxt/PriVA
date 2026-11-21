//! The control cohort step without `anno_control_vcf_allele` of
//! annotation_vcf.sh: the control VCF's AC, AN, AF and AC_Hom (counted with
//! `bcftools +fill-tags` when it lacks them) copied as control_AC,
//! control_AN, control_AF and control_nhomalt (`bcftools annotate`). Without
//! a control VCF the counts come from the cohort's own unaffected samples of
//! the PED file (`bcftools view -s | bcftools +fill-tags`); when there are
//! none, `control_vcf_optional: true` keeps the VCF as it is.
use super::tools::{
    Context, Result, args, log, move_with_index, quoted_path, random_tag, refuse, remove,
    replace_first,
};
use std::collections::BTreeSet;
use std::fs;

const CONTROL_TAGS: &str = "control_AC,control_AN,control_AF,control_nhomalt";
const SOURCE_TAGS: &str = "AC,AN,AF,AC_Hom";

/// `mktemp --tmpdir=$TMPDIR NAME.XXXXXX.vcf.gz`.
fn temporary(context: &Context, name: &str) -> String {
    context
        .temp_dir()
        .join(format!("{name}.{}.vcf.gz", random_tag()))
        .display()
        .to_string()
}

/// `compute_control_stats_from_cohort VCF PED THREADS`: a VCF of the PED's
/// unaffected samples (phenotype 1) found in the cohort, with fill-tags'
/// counts; None when there is none.
fn cohort_controls(
    context: &Context,
    vcf: &str,
    ped: &str,
    threads: &str,
) -> Result<Option<String>> {
    if let Some(error) = quoted_path(ped, "file", "ped_file") {
        log(&format!("Error: {error}"));
        log("ERROR: ped_file not available for cohort-based control stat fallback");
        return Ok(None);
    }
    if let Some(error) = quoted_path(vcf, "file", "input_vcf") {
        log(&format!("Error: {error}"));
        return Ok(None);
    }
    let mut controls = BTreeSet::new();
    for line in fs::read_to_string(ped)?.lines() {
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.get(5) == Some(&"1") {
            controls.insert(fields.get(1).copied().unwrap_or("").to_owned());
        }
    }
    if controls.is_empty() {
        log(&format!(
            "WARNING: No control samples (Phenotype=1) found in PED file: {ped}"
        ));
        return Ok(None);
    }
    let samples: BTreeSet<String> = context
        .output("bcftools", &["query", "-l", vcf])
        .lines()
        .map(str::to_owned)
        .collect();
    let found: Vec<String> = controls.intersection(&samples).cloned().collect();
    if found.is_empty() || found.iter().all(String::is_empty) {
        log(&format!(
            "WARNING: No control samples from PED file found in cohort VCF: {vcf}"
        ));
        return Ok(None);
    }
    let list = found.join(",");
    log(&format!(
        "Found {} control samples in cohort VCF: {list}",
        found.len()
    ));
    let output = temporary(context, "cohort_control_stats");
    let mut view = vec!["view".to_owned()];
    view.extend(args(&[
        "--threads",
        threads,
        "--force-samples",
        "-s",
        &list,
        vcf,
    ]));
    let filled = context.pipeline(
        &[
            ("bcftools", view),
            (
                "bcftools",
                args(&[
                    "+fill-tags",
                    "--threads",
                    threads,
                    "-Oz",
                    "-o",
                    &output,
                    "--",
                    "-t",
                    SOURCE_TAGS,
                ]),
            ),
        ],
        false,
    )? && context.succeeds("tabix", &args(&["-p", "vcf", "-f", &output]));
    if !filled {
        log("ERROR: Failed to compute control stats from cohort VCF");
        remove(&[&output, &format!("{output}.tbi")]);
        return Ok(None);
    }
    log(&format!(
        "Cohort-derived control stats written to: {output}"
    ));
    Ok(Some(output))
}

/// `anno_control_vcf_allele VCF control_vcf ped_file threads
/// control_vcf_optional`.
pub(crate) fn run(context: &Context, input: &str) -> Result<()> {
    let mut control = context.value("control_vcf");
    let ped = context.value("ped_file");
    let configured = context.value("threads");
    let threads = if configured.is_empty() {
        "1"
    } else {
        configured.as_str()
    };
    let optional = context.value("control_vcf_optional");
    let optional = if optional.is_empty() {
        "false"
    } else {
        optional.as_str()
    };
    let output = replace_first(input, ".vcf.gz", &format!(".{}.vcf.gz", random_tag()));
    if let Some(error) = quoted_path(input, "file", "input_vcf") {
        return refuse(format!("ERROR: input_vcf is not a valid file ({error})"));
    }
    let mut cohort = false;
    if let Some(error) = quoted_path(&control, "file", "control_vcf") {
        log(&format!("Error: {error}"));
        log(
            "WARNING: control_vcf is not provided or not a valid file. Attempting fallback using control samples from cohort VCF.",
        );
        match cohort_controls(context, input, &ped, threads)? {
            Some(path) => {
                control = path;
                cohort = true;
                log(&format!("Using cohort-derived control VCF: {control}"));
            }
            None if optional == "true" => {
                log(
                    "WARNING: Could not derive control stats from cohort VCF either. Skipping control allele annotation (control_vcf_optional=true).",
                );
                return Ok(());
            }
            None => {
                return refuse(
                    "ERROR: Could not derive control stats from cohort VCF either. Aborting control allele annotation.",
                );
            }
        }
    }
    if !input.ends_with(".vcf.gz") {
        return refuse("ERROR: input_vcf must be bgzipped (.vcf.gz)");
    }
    if !control.ends_with(".vcf.gz") {
        return refuse("ERROR: control_vcf must be bgzipped (.vcf.gz)");
    }
    if !context.vcf_valid(input, 1, None) || !context.vcf_valid(&control, 1, None) {
        return refuse(format!("{input} or {control} is not a valid VCF"));
    }
    if context.has_info_tags(input, CONTROL_TAGS) {
        log(&format!(
            "Input VCF {input} already contains control allele info tags: {CONTROL_TAGS}. Skipping."
        ));
        return Ok(());
    }
    let processed = if context.has_info_tags(&control, SOURCE_TAGS) {
        log(&format!(
            "Control VCF {control} already has required tags ({SOURCE_TAGS})."
        ));
        control.clone()
    } else {
        log(&format!(
            "Control VCF {control} missing one or more required tags ({SOURCE_TAGS}). Calculating..."
        ));
        let filled = temporary(context, "control_processed");
        if !std::path::Path::new(&format!("{control}.tbi")).is_file()
            && !context.succeeds("tabix", &args(&["-p", "vcf", "-f", &control]))
        {
            return refuse(format!("ERROR: Failed to index control VCF {control}"));
        }
        let fill = args(&[
            "+fill-tags",
            &control,
            "-Oz",
            "-o",
            &filled,
            "--",
            "-t",
            SOURCE_TAGS,
        ]);
        if !context.succeeds("bcftools", &fill)
            || !context.succeeds("tabix", &args(&["-p", "vcf", "-f", &filled]))
        {
            return refuse(format!(
                "ERROR: Failed to calculate tags for control VCF {control}"
            ));
        }
        filled
    };
    let cleanup = |processed: &str| {
        if processed != control {
            remove(&[processed, &format!("{processed}.tbi")]);
        }
    };
    log(&format!(
        "Annotating {input} with control allele info from {processed}"
    ));
    let annotated = context.succeeds(
        "bcftools",
        &args(&[
            "annotate",
            "-a",
            &processed,
            "-c",
            "INFO/control_AC:=INFO/AC,INFO/control_AN:=INFO/AN,INFO/control_AF:=INFO/AF,INFO/control_nhomalt:=INFO/AC_Hom",
            "-Oz",
            "-o",
            &output,
            input,
        ]),
    ) && context.succeeds("tabix", &args(&["-p", "vcf", "-f", &output]));
    if !annotated {
        cleanup(&processed);
        remove(&[&output, &format!("{output}.tbi")]);
        return refuse(format!(
            "ERROR: Failed to annotate {input} with control allele info."
        ));
    }
    if move_with_index(&output, input).is_err() {
        cleanup(&processed);
        remove(&[&output, &format!("{output}.tbi")]);
        return refuse(format!(
            "ERROR: Failed to move temporary output {output} to {input}"
        ));
    }
    cleanup(&processed);
    log(&format!(
        "Annotation with control cohort allele frequencies complete for {input}."
    ));
    if cohort {
        remove(&[&control, &format!("{control}.tbi")]);
    }
    Ok(())
}
