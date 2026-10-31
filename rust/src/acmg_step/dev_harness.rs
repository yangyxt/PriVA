//! `priva acmg-step`: runs one piece of the converted ACMG step on files the
//! Python capture script wrote, so each thread can compare its module with the
//! Python function before the whole step exists. Not part of the pipeline.
//!
//! ```text
//! priva acmg-step table read --input TSV --kwargs JSON --output FRAME
//! priva acmg-step table write --input FRAME --output TSV
//! priva acmg-step table drop-duplicates --input FRAME --output FRAME [--subset JSON]
//! priva acmg-step inspect {pickle|json} PATH
//! priva acmg-step MODULE FUNCTION --input FRAME --output FRAME [--options JSON]
//! ```
//!
//! `--kwargs` is the JSON of the keyword arguments given to `pd.read_table`
//! (`dtype`, `keep_default_na`, `header`, `names`, `usecols`, `comment`,
//! `low_memory`, `sep`). A failure prints the error and exits with status 1;
//! a bad command line exits with status 2.
use super::frame_file::{read_frame, write_frame};
use super::types::{HarnessCall, HarnessEntry};
use super::*;
use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use crate::pandas_polars;
use crate::pandas_table::{Header, ReadOptions, Table, UseColumns, read_table_with, write_table};
use std::path::Path;

const USAGE: &str = "usage: priva acmg-step table read --input TSV --kwargs JSON --output FRAME
       priva acmg-step table write --input FRAME --output TSV
       priva acmg-step table drop-duplicates --input FRAME --output FRAME [--subset JSON]
       priva acmg-step inspect {pickle|json} PATH
       priva acmg-step MODULE FUNCTION --input FRAME --output FRAME [--options JSON]";

/// The modules `priva acmg-step MODULE FUNCTION` reaches, by Python script name.
const MODULES: [(&str, HarnessEntry); 26] = [
    ("acmg_consequence", acmg_consequence::run_harness),
    ("stats", stats::run_harness),
    ("acmg_criteria_assign", acmg_criteria_assign::run_harness),
    ("acmg_inheritance", acmg_inheritance::run_harness),
    ("clinvar_vcv", clinvar_vcv::run_harness),
    (
        "gene_mechanism_annotation",
        gene_mechanism_annotation::run_harness,
    ),
    ("gene_mechanism_common", gene_mechanism_common::run_harness),
    (
        "gene_mechanism_conditions",
        gene_mechanism_conditions::run_harness,
    ),
    ("gene_mechanism_hub", gene_mechanism_hub::run_harness),
    (
        "gene_mechanism_resources",
        gene_mechanism_resources::run_harness,
    ),
    ("acmg_scoring", acmg_scoring::run_harness),
    (
        "acmg_variant_mechanism",
        acmg_variant_mechanism::run_harness,
    ),
    (
        "gene_mechanism_variants",
        gene_mechanism_variants::run_harness,
    ),
    (
        "acmg_ps1_pm5_known_variant",
        acmg_ps1_pm5_known_variant::run_harness,
    ),
    (
        "acmg_pvs1_null_variant",
        acmg_pvs1_null_variant::run_harness,
    ),
    (
        "alternative_start_codon",
        alternative_start_codon::run_harness,
    ),
    (
        "acmg_bs2_bp5_observation",
        acmg_bs2_bp5_observation::run_harness,
    ),
    (
        "acmg_pm2_bs1_ba1_frequency",
        acmg_pm2_bs1_ba1_frequency::run_harness,
    ),
    (
        "acmg_pm1_pp2_bp1_missense",
        acmg_pm1_pp2_bp1_missense::run_harness,
    ),
    (
        "acmg_pm4_bp3_protein_length",
        acmg_pm4_bp3_protein_length::run_harness,
    ),
    (
        "acmg_pp3_bp4_bp7_insilico",
        acmg_pp3_bp4_bp7_insilico::run_harness,
    ),
    (
        "acmg_ps3_bs3_pp5_bp6_reported",
        acmg_ps3_bs3_pp5_bp6_reported::run_harness,
    ),
    ("mavedb_interpreter", mavedb_interpreter::run_harness),
    (
        "acmg_pp1_bs4_bp2_pm3_family",
        acmg_pp1_bs4_bp2_pm3_family::run_harness,
    ),
    ("acmg_ps2_pm6_denovo", acmg_ps2_pm6_denovo::run_harness),
    ("determine_phase", determine_phase::run_harness),
    (
        "find_cosegregation_vars",
        find_cosegregation_vars::run_harness,
    ),
];

fn usage_error(message: &str) -> ! {
    eprintln!("{USAGE}");
    eprintln!("priva acmg-step: error: {message}");
    std::process::exit(2);
}

/// `--name value` pairs after the positional words.
fn named_options(args: &[String], allowed: &[&str]) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut index = 0;
    while index < args.len() {
        let name = args[index].trim_start_matches("--");
        if !args[index].starts_with("--") || !allowed.contains(&name) {
            usage_error(&format!("unexpected argument {}", args[index]));
        }
        let Some(value) = args.get(index + 1) else {
            usage_error(&format!("argument {} expects a value", args[index]));
        };
        out.push((name.to_owned(), value.clone()));
        index += 2;
    }
    out
}

fn option<'a>(options: &'a [(String, String)], name: &str) -> Option<&'a str> {
    options
        .iter()
        .rev()
        .find(|(key, _)| key == name)
        .map(|(_, value)| value.as_str())
}

fn required<'a>(options: &'a [(String, String)], name: &str) -> &'a str {
    option(options, name).unwrap_or_else(|| usage_error(&format!("--{name} is required")))
}

fn parse_json(text: &str) -> Result<Json> {
    json_text::parse(text, true).map_err(invalid)
}

fn text_list(value: &Json, what: &str) -> Result<Vec<String>> {
    match value {
        Json::Array(items) => items
            .iter()
            .map(|item| {
                item.as_str()
                    .map(str::to_owned)
                    .ok_or_else(|| invalid(format!("{what} must be a list of texts")))
            })
            .collect(),
        _ => Err(invalid(format!("{what} must be a list"))),
    }
}

/// The `pd.read_table` keyword arguments as `ReadOptions`.
pub(crate) fn read_options(kwargs: &Json) -> Result<ReadOptions> {
    let mut options = ReadOptions::default();
    let items = kwargs
        .as_object()
        .ok_or_else(|| invalid("--kwargs must be a JSON object"))?;
    let mut names: Option<Vec<String>> = None;
    let mut header_none = false;
    for (key, value) in items {
        match (key.as_str(), value) {
            ("dtype", Json::Str(dtype)) if dtype == "str" => options.all_text = true,
            ("keep_default_na", Json::Bool(flag)) => options.keep_default_na = *flag,
            ("header", Json::Null) => header_none = true,
            ("header", Json::Int(0)) => {}
            ("header", Json::Str(header)) if header == "infer" => {}
            ("names", value) => names = Some(text_list(value, "names")?),
            ("usecols", Json::Array(items)) => {
                options.use_columns =
                    Some(if items.iter().all(|item| matches!(item, Json::Int(_))) {
                        UseColumns::Positions(
                            items
                                .iter()
                                .map(|item| match item {
                                    Json::Int(position) if *position >= 0 => Ok(*position as usize),
                                    _ => Err(invalid("usecols positions must not be negative")),
                                })
                                .collect::<Result<_>>()?,
                        )
                    } else {
                        UseColumns::Names(text_list(value, "usecols")?)
                    })
            }
            ("comment", Json::Str(comment)) if comment.len() == 1 => {
                options.comment = Some(comment.as_bytes()[0])
            }
            ("low_memory", Json::Bool(_)) => {}
            ("sep", Json::Str(separator)) if separator == "\t" => {}
            _ => {
                return Err(invalid(format!(
                    "read_table option {key}={} is not reproduced",
                    json_text::dumps(value, false).unwrap_or_default()
                )));
            }
        }
    }
    match (header_none, names) {
        (true, Some(names)) => options.header = Header::Names(names),
        (false, None) => {}
        _ => return Err(invalid("header=None goes with names=[...] here")),
    }
    Ok(options)
}

fn table_command(args: &[String]) -> Result<()> {
    let Some(action) = args.first() else {
        usage_error("table needs read, write or drop-duplicates");
    };
    match action.as_str() {
        "read" => {
            let options = named_options(&args[1..], &["input", "kwargs", "output"]);
            let kwargs = parse_json(required(&options, "kwargs"))?;
            let table = read_table_with(
                Path::new(required(&options, "input")),
                &read_options(&kwargs)?,
            )?;
            write_frame(Path::new(required(&options, "output")), &table)
        }
        "write" => {
            let options = named_options(&args[1..], &["input", "output"]);
            let table = read_frame(Path::new(required(&options, "input")))?;
            write_table(Path::new(required(&options, "output")), &table)
        }
        "drop-duplicates" => {
            let options = named_options(&args[1..], &["input", "output", "subset"]);
            let table = read_frame(Path::new(required(&options, "input")))?;
            let subset = match option(&options, "subset") {
                Some(text) => Some(text_list(&parse_json(text)?, "--subset")?),
                None => None,
            };
            let keys: Option<Vec<&str>> = subset
                .as_ref()
                .map(|names| names.iter().map(String::as_str).collect());
            let rows = pandas_polars::drop_duplicates_first(&table, keys.as_deref())?;
            write_frame(
                Path::new(required(&options, "output")),
                &table.take_rows(&rows),
            )
        }
        other => usage_error(&format!("unknown table action {other}")),
    }
}

fn inspect_command(args: &[String]) -> Result<()> {
    match args {
        [kind, path] if kind == "pickle" => {
            print!("{}", resource_summary::pickle_summary(Path::new(path))?);
            Ok(())
        }
        [kind, path] if kind == "json" => {
            print!("{}", resource_summary::json_summary(Path::new(path))?);
            Ok(())
        }
        _ => usage_error("inspect needs pickle or json and a path"),
    }
}

fn module_command(module: &str, args: &[String]) -> Result<()> {
    let Some((_, entry)) = MODULES.iter().find(|(name, _)| *name == module) else {
        usage_error(&format!("unknown module {module}"));
    };
    let Some(function) = args.first() else {
        usage_error("a module needs a function name");
    };
    let options = named_options(&args[1..], &["input", "output", "options"]);
    let call = HarnessCall {
        table: read_frame(Path::new(required(&options, "input")))?,
        options: match option(&options, "options") {
            Some(text) => parse_json(text)?,
            None => Json::Object(Vec::new()),
        },
    };
    let result: Table = entry(function, &call)?;
    write_frame(Path::new(required(&options, "output")), &result)
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let outcome = match args.first().map(String::as_str) {
        None | Some("-h") | Some("--help") => {
            println!("{USAGE}");
            return Ok(());
        }
        Some("table") => table_command(&args[1..]),
        Some("inspect") => inspect_command(&args[1..]),
        Some(module) => module_command(module, &args[1..]),
    };
    if let Err(error) = outcome {
        eprintln!("priva acmg-step: {error}");
        std::process::exit(1);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_table_keywords_become_read_options() {
        let kwargs = parse_json(
            r##"{"header": null, "names": ["#FamilyID", "IndividualID"], "usecols": [0, 1],
                "dtype": "str", "keep_default_na": false, "low_memory": false}"##,
        )
        .unwrap();
        let options = read_options(&kwargs).unwrap();
        assert!(options.all_text && !options.keep_default_na);
        assert_eq!(
            options.header,
            Header::Names(vec!["#FamilyID".into(), "IndividualID".into()])
        );
        assert_eq!(options.use_columns, Some(UseColumns::Positions(vec![0, 1])));
        let by_name =
            read_options(&parse_json(r##"{"usecols": ["a", "c"], "comment": "#"}"##).unwrap())
                .unwrap();
        assert_eq!(
            by_name.use_columns,
            Some(UseColumns::Names(vec!["a".into(), "c".into()]))
        );
        assert_eq!(by_name.comment, Some(b'#'));
        assert!(read_options(&parse_json(r#"{"dtype": "float"}"#).unwrap()).is_err());
        assert!(read_options(&parse_json(r#"{"header": null}"#).unwrap()).is_err());
    }

    #[test]
    fn every_module_has_a_harness_entry() {
        let call = HarnessCall {
            table: Table {
                names: Vec::new(),
                columns: Vec::new(),
                rows: 0,
            },
            options: Json::Object(Vec::new()),
        };
        for (name, entry) in MODULES {
            let error = entry("anything", &call).unwrap_err().to_string();
            assert!(error.contains(name), "{name}: {error}");
        }
    }
}
