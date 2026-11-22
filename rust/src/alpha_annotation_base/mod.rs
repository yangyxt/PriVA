//! Optional procedures corresponding to alphamissense_annotation_base.py.
//! Original installer callers stay unchanged. VEP is queried only for runtime
//! provenance; scientific annotation and model generation are never invoked.
pub(crate) mod header;
pub(crate) mod io;
mod motif;
mod profile;
mod publication;

use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use crate::python_args::parse_long_option_occurrences;
use crate::python_path::path_value;
use std::collections::HashMap;
use std::path::Path;

const SCHEMA: &str = "priva-alphamissense-annotation-base-v1";
const MOTIF_SCHEMA: &str = "priva-alphamissense-motif-input-v1";
fn error(category: &'static str, message: impl Into<String>) -> Box<dyn std::error::Error> {
    crate::builder_error::category(category, invalid(message))
}
fn require(condition: bool, message: impl Into<String>) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(error("ValueError", message))
    }
}
fn object<const N: usize>(values: [(&str, Json); N]) -> Json {
    Json::Object(values.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}
const HELP: &str = "usage: priva alpha-annotation-base {path,capture,check,recover,publish,restore,adopt,runtime-associated,motif-binding,motif-check,adopt-motif,motif-runtime}\n\
--config YAML --source-root CHECKOUT [--base VCF] [--python PYTHON] [--helper-root CHECKOUT]\n\
capture [--output JSON]; publish --staged-converted VCF --staged-annotation VCF --expected-inputs JSON; adopt --lineage JSON\n\
motif-check --artifact FILE; adopt-motif --artifact FILE --lineage JSON --mapping TSV [--input FILE ...]\n\
Opt-in only. motif-runtime observes the selected Python environment; adopt-motif retains the original read-only BioMart contract dependency.";
fn option<'a>(options: &'a HashMap<String, String>, name: &str) -> Result<&'a str> {
    options
        .get(name)
        .map(String::as_str)
        .ok_or_else(|| invalid(format!("missing --{name}")))
}
pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(command) = args.first().map(String::as_str) else {
        return Err(invalid(HELP));
    };
    if matches!(command, "--help" | "-h") {
        println!("{HELP}");
        return Ok(());
    }
    let mut names = vec!["config", "source-root", "base", "python", "helper-root"];
    match command {
        "path" | "check" | "recover" | "restore" | "runtime-associated" | "motif-binding"
        | "motif-runtime" => {}
        "capture" => names.push("output"),
        "publish" => names.extend(["staged-converted", "staged-annotation", "expected-inputs"]),
        "adopt" => names.push("lineage"),
        "motif-check" => names.push("artifact"),
        "adopt-motif" => names.extend(["artifact", "lineage", "mapping", "input"]),
        _ => return Err(invalid(HELP)),
    }
    let Some(occurrences) =
        parse_long_option_occurrences(&args[1..], &names, &[]).map_err(invalid)?
    else {
        println!("{HELP}");
        return Ok(());
    };
    let inputs = occurrences
        .iter()
        .filter(|(k, _)| k == "input")
        .map(|(_, v)| path_value(v))
        .collect::<Vec<_>>();
    let options: HashMap<_, _> = occurrences.into_iter().collect();
    if command != "motif-runtime" {
        option(&options, "config")?;
        option(&options, "source-root")?;
    }
    for name in match command {
        "publish" => vec!["staged-converted", "staged-annotation", "expected-inputs"],
        "adopt" => vec!["lineage"],
        "motif-check" => vec!["artifact"],
        "adopt-motif" => vec!["artifact", "lineage", "mapping"],
        _ => vec![],
    } {
        option(&options, name)?;
    }
    crate::builder_error::finish((|| {
        let python = options.get("python").map_or("python", String::as_str);
        if command == "motif-runtime" {
            println!(
                "{}",
                json_text::dumps(&profile::motif_runtime(python)?, true).map_err(invalid)?
            );
            return Ok(());
        }
        let config = path_value(option(&options, "config")?);
        let source = path_value(option(&options, "source-root")?);
        let base = options.get("base").map(|p| path_value(p));
        let helper = options
            .get("helper-root")
            .map(|p| path_value(p))
            .unwrap_or_else(|| source.clone());
        let context = publication::Context {
            config: &config,
            source: &source,
            base: base.as_deref(),
            python,
            helper: &helper,
        };
        match command {
            "path" => println!("{}", io::text(io::get(&context.settings()?.paths, "base")?)),
            "capture" => {
                let value = context.capture()?;
                if let Some(output) = options.get("output") {
                    io::atomic_json(&path_value(output), &value)?;
                } else {
                    println!("{}", json_text::dumps(&value, true).map_err(invalid)?);
                }
            }
            "check" => {
                context.check(false)?;
            }
            "recover" => context.recover()?,
            "restore" => context.restore()?,
            "publish" => context.publish(
                &path_value(option(&options, "staged-converted")?),
                &path_value(option(&options, "staged-annotation")?),
                &io::read_json(&path_value(option(&options, "expected-inputs")?))?,
                Json::Null,
            )?,
            "adopt" => context.adopt(&path_value(option(&options, "lineage")?))?,
            "runtime-associated" => {
                let record = context.check_fields(false, &["native_runtime"])?;
                let paths = context.settings()?.paths;
                publication::validate_runtime(
                    io::optional(&record, "native_runtime")?,
                    &paths,
                    io::get(&record, "science")?,
                    io::get(&record, "base_sha256")?,
                )?;
            }
            "motif-binding" => println!(
                "{}",
                json_text::dumps(&motif::binding(&context)?, true).map_err(invalid)?
            ),
            "motif-check" => motif::check(&context, Path::new(option(&options, "artifact")?))?,
            "adopt-motif" => motif::adopt(
                &context,
                Path::new(option(&options, "artifact")?),
                Path::new(option(&options, "lineage")?),
                Path::new(option(&options, "mapping")?),
                &inputs,
            )?,
            _ => unreachable!(),
        }
        Ok(())
    })())
}

/// Original motif binding for native publication; runtime observation remains
/// in the explicitly selected Python environment, never a scientific rerun.
pub(crate) fn motif_binding(
    config: &Path,
    source: &Path,
    base: Option<&Path>,
    python: &str,
) -> Result<Json> {
    motif::binding(&publication::Context {
        config,
        source,
        base,
        python,
        helper: source,
    })
}
