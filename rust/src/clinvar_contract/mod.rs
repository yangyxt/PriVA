//! Opt-in procedures corresponding to scripts/clinvar_annotation_contract.py.
//! Header validation, native cache shapes and completion metadata are separate
//! operations. The original Python source asset is explicit in pair operations:
//! its content identity must not be replaced by the native executable's hash.
mod cache;
pub(crate) mod header;
mod legacy;
mod metadata;
pub(crate) use metadata::resolve as resolve_non_strict;

use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use crate::python_args::{Parsed, parse_long_options};
use std::collections::HashMap;

fn error(category: &'static str, message: impl Into<String>) -> Box<dyn std::error::Error> {
    crate::builder_error::category(category, invalid(message))
}
fn object<const N: usize>(values: [(&str, Json); N]) -> Json {
    Json::Object(values.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}
fn strings(values: impl IntoIterator<Item = String>) -> Json {
    Json::Array(values.into_iter().map(Json::Str).collect())
}
fn required<'a>(options: &'a HashMap<String, String>, name: &str) -> Result<&'a str> {
    options
        .get(name)
        .map(String::as_str)
        .ok_or_else(|| invalid(format!("missing --{name}")))
}
fn print(value: &Json) -> Result<()> {
    println!("{}", json_text::dumps(value, true).map_err(invalid)?);
    Ok(())
}
const USAGE: &str = "priva clinvar-contract header --vcf VCF_OR_BCF [--stage aa|gene|af|vep|complete]\n\
priva clinvar-contract cache --input NATIVE_ARCHIVE --kind aa|splice|gene-af|exon-af|gene\n\
priva clinvar-contract legacy-cache --input PYTHON_PICKLE --kind aa|splice|gene-af|exon-af|gene\n\
priva clinvar-contract af-paths-check --vcf VCF --config YAML --gene GENE --exon EXON\n\
priva clinvar-contract pair-metadata --input PATH\n\
priva clinvar-contract pair-inputs --vcf VCF --producer PYTHON_ASSET --contract-source PYTHON_CONTRACT\n\
priva clinvar-contract pair-capture|pair-check|pair-complete --vcf VCF --producer PYTHON_ASSET \
--contract-source PYTHON_CONTRACT --aa AA --splice SPLICE --receipt JSON [--captured JSON] [--config YAML]\n\
Header, archive and pickle shape, output path and AA/splice pair procedures.";

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(command) = args.first().map(String::as_str) else {
        return Err(invalid(USAGE));
    };
    if matches!(command, "--help" | "-h") {
        println!("{USAGE}");
        return Ok(());
    }
    let names: &[&str] = match command {
        "header" => &["vcf", "stage"],
        "cache" | "legacy-cache" => &["input", "kind"],
        "af-paths-check" => &["vcf", "config", "gene", "exon"],
        "pair-metadata" => &["input"],
        "pair-inputs" => &["vcf", "producer", "contract-source"],
        "pair-capture" | "pair-check" | "pair-complete" => &[
            "vcf",
            "producer",
            "contract-source",
            "aa",
            "splice",
            "receipt",
            "captured",
            "config",
        ],
        _ => return Err(invalid(USAGE)),
    };
    let Parsed::Values(options) = parse_long_options(&args[1..], names).map_err(invalid)? else {
        println!("{USAGE}");
        return Ok(());
    };
    // Admission is exit 2; source procedure failures after admission are exit 1.
    for name in names
        .iter()
        .filter(|n| !matches!(**n, "stage" | "captured" | "config"))
    {
        required(&options, name)?;
    }
    if command == "af-paths-check" {
        required(&options, "config")?;
    }
    if command == "header" {
        let stage = options.get("stage").map_or("complete", String::as_str);
        if !["aa", "gene", "af", "vep", "complete"].contains(&stage) {
            return Err(invalid("unknown validation stage"));
        }
    }
    if matches!(command, "cache" | "legacy-cache")
        && !["aa", "splice", "gene-af", "exon-af", "gene"].contains(&required(&options, "kind")?)
    {
        return Err(invalid("unknown cache kind"));
    }
    crate::builder_error::finish(execute(command, &options))
}

fn execute(command: &str, options: &HashMap<String, String>) -> Result<()> {
    match command {
        "header" => print(&header::validate(
            required(options, "vcf")?,
            options.get("stage").map_or("complete", String::as_str),
        )?),
        "cache" => print(&cache::validate(
            required(options, "input")?,
            required(options, "kind")?,
        )?),
        "legacy-cache" => print(&legacy::validate(
            required(options, "input")?,
            required(options, "kind")?,
        )?),
        "pair-metadata" => print(&metadata::file(required(options, "input")?)?),
        "pair-inputs" => print(&metadata::inputs(
            required(options, "vcf")?,
            required(options, "producer")?,
            required(options, "contract-source")?,
        )?),
        "af-paths-check" => {
            let vcf = required(options, "vcf")?;
            metadata::distinct(
                &[
                    vcf.to_owned(),
                    format!("{vcf}.tbi"),
                    format!("{vcf}.csi"),
                    required(options, "config")?.to_owned(),
                    required(options, "gene")?.to_owned(),
                    required(options, "exon")?.to_owned(),
                ],
                &[],
            )
        }
        _ => metadata::pair(command, options),
    }
}
