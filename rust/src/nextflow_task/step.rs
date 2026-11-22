//! `priva task step`: one pipeline step outside Nextflow, in the current
//! folder, with CONFIG as its configuration and ROOT as the checkout whose
//! `data/` it reads. It verifies no contract; it is for checking a step
//! against the shell function it replaced and for looking into one family
//! or cohort by hand. The priva program is PRIVA_BINARY, else this program.
use super::{Options, Result};
use crate::pipeline_steps::{cadd, clinvar, control, family, gnomad, preprocess, tools, vep};
use std::path::PathBuf;

/// The steps `--step` takes.
pub(super) const STEPS: &[&str] = &[
    "preprocess",
    "gnomad",
    "clinvar",
    "control",
    "vep",
    "vep-worker",
    "cadd",
    "filtration",
    "prioritization",
    "combine",
    "splicing",
    "utr",
    "acmg",
];

/// The options each step needs besides --step, --root and --config.
pub(super) fn required(step: &str) -> &'static [&'static str] {
    match step {
        "preprocess" | "cadd" => &["vcf", "output"],
        "gnomad" | "clinvar" | "control" | "vep" | "filtration" => &["vcf"],
        "prioritization" | "combine" => &["vcf", "cadd"],
        "vep-worker" => &["vcf", "output", "stats", "forks"],
        _ => &["table"],
    }
}

pub(super) fn run(options: &Options) -> Result<()> {
    let root = PathBuf::from(options.value("root"));
    let config_path = PathBuf::from(options.value("config"));
    let settings = tools::Settings::read(&config_path)?;
    let binary = match std::env::var("PRIVA_BINARY") {
        Ok(path) if !path.is_empty() => path,
        _ => std::env::current_exe()?.display().to_string(),
    };
    let mut env = vec![("PRIVA_BINARY".to_owned(), binary.clone())];
    if let Ok(temp) = std::env::var("TMPDIR") {
        env.push(("TMPDIR".to_owned(), temp));
    }
    let context = tools::Context {
        root: &root,
        binary: &binary,
        config_path: &config_path,
        settings: &settings,
        env: &env,
    };
    let family_id = options.get("family").filter(|v| !v.is_empty());
    let value = |name: &str| options.value(name);
    match options.value("step") {
        // The annotation steps change VCF in place, but for these two.
        "preprocess" => preprocess::run(&context, value("vcf"), value("output")),
        "cadd" => cadd::run(&context, value("vcf"), value("output")),
        "gnomad" => gnomad::run(&context, value("vcf")),
        "clinvar" => clinvar::run(&context, value("vcf")),
        "control" => control::run(&context, value("vcf")),
        "vep" => vep::run(&context, value("vcf")),
        // anno_VEP_cache_worker: VEP on one chromosome of the ClinVar release
        // (install steps D1_chrN).
        "vep-worker" => vep::worker(
            &context,
            value("vcf"),
            &vep::Worker {
                output: value("output"),
                stats: value("stats"),
                forks: value("forks"),
            },
        ),
        "filtration" => family::filtration(&context, value("vcf"), family_id),
        "prioritization" => family::prioritize(&context, value("vcf"), value("cadd"), family_id),
        "combine" => family::combine(&context, value("vcf"), value("cadd")),
        "splicing" => family::splicing(&context, value("table")),
        "utr" => family::utr(&context, value("table")),
        _ => family::acmg(&context, value("table"), family_id),
    }
}
