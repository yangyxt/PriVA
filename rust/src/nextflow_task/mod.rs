//! `priva task` and `priva task-boundary`: the Nextflow task contract.
//!
//! This reproduces the former `bin/priva_task.py` and `bin/priva_task_boundary.py`. It
//! checks inputs and prepared resources, gives each task private copies of
//! the files the shell steps may change, and runs the unchanged shell
//! workflow scripts. It implements no variant scoring, filtering, criteria or
//! sorting.
//!
//!   priva task manifest --config CONFIG --legacy-root ROOT --output JSON
//!   priva task launch|validate|run ...
//!   priva task-boundary check|exec --task-dir DIR [-- COMMAND...]
mod boundary;
pub(crate) mod contract;
pub(crate) mod hub_lineage;
mod info_merge;
mod table_join;
mod tasks;
#[cfg(test)]
mod tests;

use std::collections::BTreeMap;
use std::error::Error;
use std::fmt;

/// A refusal reported as `PriVA compatibility contract failed: ...`.
#[derive(Debug)]
pub(crate) struct ContractError(pub(crate) String);

impl fmt::Display for ContractError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl Error for ContractError {}

pub(crate) type Result<T> = std::result::Result<T, Box<dyn Error>>;

pub(crate) fn refuse<T>(message: impl Into<String>) -> Result<T> {
    Err(Box::new(ContractError(message.into())))
}

/// Parsed `--name value` options and `--flag` switches, as argparse reads
/// them for these subcommands (`--name=value` is accepted too).
pub(crate) struct Options {
    values: BTreeMap<String, String>,
    flags: Vec<String>,
}

impl Options {
    fn parse(
        args: &[String],
        valued: &[&str],
        switches: &[&str],
    ) -> std::result::Result<Self, String> {
        let mut values = BTreeMap::new();
        let mut flags = Vec::new();
        let mut index = 0;
        while index < args.len() {
            let argument = &args[index];
            let Some(name) = argument.strip_prefix("--") else {
                return Err(format!("unrecognized arguments: {argument}"));
            };
            let (name, inline) = match name.split_once('=') {
                Some((name, value)) => (name, Some(value.to_owned())),
                None => (name, None),
            };
            if switches.contains(&name) && inline.is_none() {
                flags.push(name.to_owned());
                index += 1;
                continue;
            }
            if !valued.contains(&name) {
                return Err(format!("unrecognized arguments: {argument}"));
            }
            let value = match inline {
                Some(value) => value,
                None => {
                    index += 1;
                    match args.get(index) {
                        Some(value) if !value.starts_with("--") => value.clone(),
                        _ => return Err(format!("argument --{name}: expected one argument")),
                    }
                }
            };
            values.insert(name.to_owned(), value);
            index += 1;
        }
        Ok(Self { values, flags })
    }

    fn require(&self, names: &[&str]) -> std::result::Result<(), String> {
        let missing: Vec<String> = names
            .iter()
            .filter(|name| !self.values.contains_key(**name))
            .map(|name| format!("--{name}"))
            .collect();
        if missing.is_empty() {
            Ok(())
        } else {
            Err(format!(
                "the following arguments are required: {}",
                missing.join(", ")
            ))
        }
    }

    pub(crate) fn get(&self, name: &str) -> Option<&str> {
        self.values.get(name).map(String::as_str)
    }

    pub(crate) fn value(&self, name: &str) -> &str {
        self.get(name).unwrap_or_default()
    }

    pub(crate) fn flag(&self, name: &str) -> bool {
        self.flags.iter().any(|flag| flag == name)
    }

    fn choice(&self, name: &str, choices: &[&str]) -> std::result::Result<(), String> {
        match self.get(name) {
            Some(value) if !choices.contains(&value) => Err(format!(
                "argument --{name}: invalid choice: '{value}' (choose from {})",
                choices
                    .iter()
                    .map(|c| format!("'{c}'"))
                    .collect::<Vec<_>>()
                    .join(", ")
            )),
            _ => Ok(()),
        }
    }

    pub(crate) fn int(&self, name: &str) -> std::result::Result<i64, String> {
        let value = self.value(name);
        value
            .trim()
            .parse::<i64>()
            .map_err(|_| format!("argument --{name}: invalid int value: '{value}'"))
    }
}

fn usage_error(command: &str, message: &str) -> ! {
    eprintln!("usage: priva {command} --help\npriva {command}: error: {message}");
    std::process::exit(2)
}

// --ped is optional (no pedigree: one sample-set run without family
// filtering); --family selects one PED family, else every family runs.
const VALIDATE_REQUIRED: &[&str] = &[
    "vcf",
    "vcf-index",
    "config",
    "resources",
    "base",
    "assembly",
    "mode",
    "legacy-root",
    "source-digest",
    "software-id",
    "vcf-hub-mode",
    "cadd-hub-mode",
];
const LAUNCH_REQUIRED: &[&str] = &[
    "vcf",
    "vcf-index",
    "config",
    "resources",
    "assembly",
    "mode",
    "legacy-root",
    "source-digest",
    "software-id",
    "vcf-hub-mode",
    "cadd-hub-mode",
    "outdir",
    "annotation-cpus",
    "family-cpus",
];

/// Task phases. `annotation` runs the whole annotation script; the
/// `annotation-*` steps run one of its functions each, and `annotation-merge`
/// joins their INFO fields. `prioritization` runs the whole prioritization
/// script; `prioritization-combine`, `-splicing`, `-utr` and `-acmg` run one of
/// its steps each, and `prioritization-join` joins the splicing and 5' UTR
/// copies of the family table.
pub(crate) const PHASES: &[&str] = &[
    "annotation",
    "annotation-preprocess",
    "annotation-gnomad",
    "annotation-clinvar",
    "annotation-control",
    "annotation-vep",
    "annotation-cadd",
    "annotation-merge",
    "filtration",
    "prioritization",
    "prioritization-combine",
    "prioritization-splicing",
    "prioritization-utr",
    "prioritization-join",
    "prioritization-acmg",
];

fn parse_task(args: &[String]) -> std::result::Result<(String, Options), String> {
    let Some(command) = args.first() else {
        return Err(
            "the following arguments are required: command (manifest, validate, launch, run)"
                .to_owned(),
        );
    };
    let rest = &args[1..];
    let options = match command.as_str() {
        "manifest" => {
            let options = Options::parse(
                rest,
                &[
                    "config",
                    "legacy-root",
                    "mode",
                    "vcf-hub-mode",
                    "cadd-hub-mode",
                    "output",
                ],
                &["checksum"],
            )?;
            options.require(&["config", "legacy-root", "output"])?;
            options.choice("mode", &["raw", "annotated"])?;
            options.choice("vcf-hub-mode", &["off", "qualified-copy"])?;
            options.choice("cadd-hub-mode", &["off", "copy"])?;
            options
        }
        "validate" => {
            let mut valued = VALIDATE_REQUIRED.to_vec();
            valued.extend(["cadd", "task-boundary", "ped", "family", "resources-lock"]);
            let options = Options::parse(rest, &valued, &["stub"])?;
            options.require(VALIDATE_REQUIRED)?;
            options.choice("task-boundary", &["off", "bwrap"])?;
            options
        }
        "launch" => {
            let mut valued = LAUNCH_REQUIRED.to_vec();
            valued.extend(["cadd", "task-boundary", "ped", "family"]);
            let options = Options::parse(rest, &valued, &["stub"])?;
            options.require(LAUNCH_REQUIRED)?;
            options.int("annotation-cpus")?;
            options.int("family-cpus")?;
            options.choice("task-boundary", &["off", "bwrap"])?;
            options
        }
        "run" => {
            let options = Options::parse(
                rest,
                &[
                    "phase",
                    "bundle",
                    "vcf",
                    "vcf-index",
                    "cadd",
                    "cpus",
                    "family",
                    "gnomad-vcf",
                    "clinvar-vcf",
                    "control-vcf",
                    "vep-vcf",
                    "table",
                    "splicing-table",
                    "utr-table",
                ],
                &["stub"],
            )?;
            options.require(&["phase", "bundle", "cpus"])?;
            options.choice("phase", PHASES)?;
            // The steps after prioritization-combine read the family table.
            if matches!(
                options.value("phase"),
                "prioritization-splicing"
                    | "prioritization-utr"
                    | "prioritization-join"
                    | "prioritization-acmg"
            ) {
                options.require(&["table"])?;
            } else {
                options.require(&["vcf", "vcf-index"])?;
            }
            options.int("cpus")?;
            options
        }
        other => {
            return Err(format!(
                "argument command: invalid choice: '{other}' (choose from 'manifest', 'validate', 'launch', 'run')"
            ));
        }
    };
    Ok((command.clone(), options))
}

/// `priva task COMMAND ...`: exit 2 with the adapter's message on refusal.
pub(crate) fn dispatch(args: &[String]) -> crate::cli::Result<()> {
    if args.first().map(String::as_str) == Some("--help") {
        println!(
            "usage: priva task {{manifest,validate,launch,run}} ...\n\nNon-scientific Nextflow task contract: validate inputs and resources, isolate\nwrites, and run the unchanged PriVA shell workflow scripts."
        );
        return Ok(());
    }
    let (command, options) = parse_task(args).unwrap_or_else(|m| usage_error("task", &m));
    let outcome = match command.as_str() {
        "manifest" => tasks::build_manifest(&options),
        "validate" => tasks::validate(&options),
        "launch" => tasks::launch(&options),
        _ => tasks::run(&options),
    };
    if let Err(error) = outcome {
        eprintln!("PriVA compatibility contract failed: {error}");
        std::process::exit(2);
    }
    Ok(())
}

/// `priva task-boundary check|exec`: exit 125 on refusal, as before.
pub(crate) fn dispatch_boundary(args: &[String]) -> crate::cli::Result<()> {
    std::process::exit(boundary::main(args))
}
