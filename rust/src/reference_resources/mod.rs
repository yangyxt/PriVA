//! Native `scripts/reference_resources.py`: the pinned UCSC hg19 and hg38
//! FASTAs and the forward hg19-to-hg38 chain, planned, checked and prepared
//! without replacing a supplied input.
//!
//! `sources`, `preflight` and `install` keep the script's options, printed
//! lines, written files (indexes, source receipts, lock files and the
//! configuration rewrite) and exit statuses. A failure that main() printed as
//! `ERROR: message` exits 1 here as `priva: CATEGORY: message`, CATEGORY
//! being the Python exception class; failures the script did not catch keep
//! their class too. Differences:
//! - `--data-dir` has no default (the script used its own checkout's data
//!   folder); `preflight` and `install` require it.
//! - The command and CONFIG come before the options.
//! - The configuration is read by `yaml_config`, which refuses anchors, tags,
//!   flow (JSON) documents and other unsupported constructs with
//!   NotImplementedError.
//! - Downloads run `curl` with the pinned URL, the same Range request, status
//!   and Content-Range checks and the same partial-file handling. A failed
//!   connection or transfer is reported with curl's own message as URLError.
mod files;
mod http;
mod install;
mod plan;
mod verify;

pub(crate) use plan::{Scope, resource_plan};
pub(crate) use verify::preflight;

use crate::builder_error;
use crate::cli::{Result, invalid};
use crate::config_text::{self, Writer};
use crate::installer_prerequisites::Config;
use crate::installer_prerequisites::paths::pathlib;
use crate::json_text::{self, DumpsOptions, Json};
use crate::python_args::{parse_long_option_occurrences, require_options};
use crate::python_path::repr;
use std::collections::HashMap;
use std::error::Error;

const HELP: &str = "usage: priva reference-resources {sources,preflight,install} [CONFIG]\n\
       --data-dir DATA [--scope {selected,both,gofcards}] [--use-cache-root] [--write-config]\n\
The checks, written files, printed lines and exit statuses of scripts/reference_resources.py.\n\
The command and CONFIG come first. preflight and install need --data-dir (the script's\n\
default was its own checkout's data folder); sources prints the pinned UCSC sources.";
const COMMANDS: [&str; 3] = ["sources", "preflight", "install"];
const SCOPES: [&str; 3] = ["selected", "both", "gofcards"];
const NAMES: [&str; 2] = ["data-dir", "scope"];
const FLAGS: [&str; 2] = ["use-cache-root", "write-config"];

fn error(category: &'static str, message: impl Into<String>) -> Box<dyn Error> {
    builder_error::category(category, invalid(message))
}
fn value_error(message: impl Into<String>) -> Box<dyn Error> {
    error("ValueError", message)
}
fn require(condition: bool, message: impl Into<String>) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(value_error(message))
    }
}

/// One pinned UCSC publication. Publisher manifests and HTTP lengths were
/// verified 2026-10-16; they identify the archives the previous installer
/// used, including the original hg19 chrM.
struct Source {
    url: &'static str,
    manifest: &'static str,
    md5: &'static str,
    bytes: u64,
    identity: &'static str,
}

static SOURCES: [(&str, Source); 3] = [
    (
        "hg19",
        Source {
            url: "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
            manifest: "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/md5sum.txt",
            md5: "806c02398f5ac5da8ffd6da2d1d5d1a9",
            bytes: 948731419,
            identity: "UCSC hg19 original bigZips FASTA; original chrM; not latest/analysisSet",
        },
    ),
    (
        "hg38",
        Source {
            url: "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
            manifest: "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt",
            md5: "1c9dcaddfa41027f17cd8f7a82c7293b",
            bytes: 983659424,
            identity: "UCSC hg38 original bigZips FASTA; not latest/analysisSet",
        },
    ),
    (
        "hg19_to_hg38",
        Source {
            url: concat!(
                "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/",
                "hg19ToHg38.over.chain.gz"
            ),
            manifest: "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/md5sum.txt",
            md5: "35887f73fe5e2231656504d1f6430900",
            bytes: 227698,
            identity: "UCSC forward hg19-to-hg38 over.chain; target=hg19, query=hg38",
        },
    ),
];

/// `SOURCES[name]` for a planned resource's assembly or chain name.
fn source(name: &str) -> &'static Source {
    &SOURCES
        .iter()
        .find(|(key, _)| *key == name)
        .expect("every planned resource names a pinned source")
        .1
}

impl Source {
    /// The source as the script's dictionary, for receipts and `sources`.
    fn json(&self) -> Json {
        let text = |value: &str| Json::Str(value.to_owned());
        Json::Object(vec![
            ("url".to_owned(), text(self.url)),
            ("manifest".to_owned(), text(self.manifest)),
            ("md5".to_owned(), text(self.md5)),
            ("bytes".to_owned(), Json::Int(i128::from(self.bytes))),
            ("identity".to_owned(), text(self.identity)),
        ])
    }
}

/// `json.dumps(value, indent=2, sort_keys=True)`
fn dumps_sorted(value: &Json) -> Result<String> {
    json_text::dumps_with_options(
        value,
        DumpsOptions {
            indent: Some(2),
            sort_keys: true,
            ..DumpsOptions::default()
        },
    )
    .map_err(|message| error("ValueError", message))
}

struct Options {
    install: bool,
    config: Option<String>,
    data_dir: String,
    scope: Scope,
    use_cache_root: bool,
    write_config: bool,
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(command) = args.first() else {
        return Err(invalid(format!(
            "the following arguments are required: command\n{HELP}"
        )));
    };
    if matches!(command.as_str(), "-h" | "--help") {
        println!("{HELP}");
        return Ok(());
    }
    if !COMMANDS.contains(&command.as_str()) {
        return Err(invalid(format!(
            "argument command: invalid choice: {} (choose from {})",
            repr(command),
            COMMANDS.map(repr).join(", ")
        )));
    }
    let mut rest = &args[1..];
    let mut config = None;
    if let Some(first) = rest.first()
        && (!first.starts_with('-') || first == "-")
    {
        config = Some(pathlib(first));
        rest = &rest[1..];
    }
    let Some(occurrences) = parse_long_option_occurrences(rest, &NAMES, &FLAGS).map_err(invalid)?
    else {
        println!("{HELP}");
        return Ok(());
    };
    let mut values = HashMap::new();
    let mut flags = Vec::new();
    for (name, value) in occurrences {
        if FLAGS.contains(&name.as_str()) {
            flags.push(name);
        } else {
            values.insert(name, value);
        }
    }
    let given = values.remove("scope");
    let scope = match given.as_deref().unwrap_or("gofcards") {
        "selected" => Scope::Selected,
        "both" => Scope::Both,
        "gofcards" => Scope::Gofcards,
        _ => {
            return Err(invalid(format!(
                "argument --scope: invalid choice: {} (choose from {})",
                repr(given.as_deref().unwrap_or_default()),
                SCOPES.map(repr).join(", ")
            )));
        }
    };
    if command == "sources" {
        let sources = Json::Object(
            SOURCES
                .iter()
                .map(|(name, source)| ((*name).to_owned(), source.json()))
                .collect(),
        );
        println!("{}", dumps_sorted(&sources)?);
        return Ok(());
    }
    require_options(&values, &["data-dir"]).map_err(invalid)?;
    let flag = |name: &str| flags.iter().any(|seen| seen == name);
    let options = Options {
        install: command == "install",
        config,
        data_dir: pathlib(&values["data-dir"]),
        scope,
        use_cache_root: flag("use-cache-root"),
        write_config: flag("write-config"),
    };
    builder_error::finish(run(&options))
}

fn run(options: &Options) -> Result<()> {
    let Some(config) = &options.config else {
        return Err(value_error("Supply a configuration file"));
    };
    require(
        !options.write_config || options.install,
        "--write-config requires install",
    )?;
    let config_path = plan::checked_path(config)?;
    let original = files::read_text(&config_path)?;
    let (resources, updates) = resource_plan(
        &Config::load(&original, false)?,
        &options.data_dir,
        options.scope,
        options.use_cache_root,
    )?;
    let updated = if options.write_config {
        render_config_values(&original, &updates)?
    } else {
        original.clone()
    };
    preflight(&resources)?;
    if options.install {
        require(
            install::samtools_on_path(),
            "samtools is required to install/index references",
        )?;
        for resource in &resources {
            install::install_resource(resource)?;
        }
        if options.write_config {
            install::write_config_values(&config_path, &original, &updated)?;
        }
        for (key, value) in &updates {
            println!("{key}\t{value}");
        }
    }
    Ok(())
}

/// `render_config_values(original, updates)`: quoted path scalars; unrelated
/// YAML text, ordering and comments preserved.
fn render_config_values(original: &str, updates: &[(String, String)]) -> Result<String> {
    let values = updates
        .iter()
        .map(|(key, value)| (key.clone(), Json::Str(value.clone())))
        .collect::<Vec<_>>();
    config_text::render(original, &values, Writer::Reference)
        .map_err(|failure| error(failure.kind, failure.message))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder_error::RuntimeError;

    fn args(items: &[&str]) -> Vec<String> {
        items.iter().map(|item| (*item).to_owned()).collect()
    }

    #[test]
    fn sources_serialize_as_the_sorted_python_dictionary() {
        let text = dumps_sorted(&source("hg19_to_hg38").json()).unwrap();
        let expected = concat!(
            "{\n",
            "  \"bytes\": 227698,\n",
            "  \"identity\": \"UCSC forward hg19-to-hg38 over.chain; target=hg19, query=hg38\",\n",
            "  \"manifest\": ",
            "\"https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/md5sum.txt\",\n",
            "  \"md5\": \"35887f73fe5e2231656504d1f6430900\",\n",
            "  \"url\": \"https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/",
            "hg19ToHg38.over.chain.gz\"\n",
            "}"
        );
        assert_eq!(text, expected);
        assert_eq!(source("hg38").bytes, 983659424);
    }

    #[test]
    fn admission_errors_exit_two_and_runtime_checks_exit_one() {
        for bad in [
            args(&[]),
            args(&["fetch"]),
            args(&["install", "c.yaml", "--scope", "all", "--data-dir", "d"]),
            args(&["install", "c.yaml", "--write=1", "--data-dir", "d"]),
            args(&["install", "c.yaml", "extra", "--data-dir", "d"]),
            args(&["install", "c.yaml"]),
        ] {
            let failure = dispatch(&bad).unwrap_err();
            assert!(!failure.is::<RuntimeError>(), "{failure}");
        }
        let failure = dispatch(&args(&["install", "--data-dir", "d"])).unwrap_err();
        assert_eq!(
            failure.to_string(),
            "ValueError: Supply a configuration file"
        );
        let options = args(&["preflight", "c.yaml", "--data-dir", "d", "--write"]);
        let failure = dispatch(&options).unwrap_err();
        assert_eq!(
            failure.to_string(),
            "ValueError: --write-config requires install"
        );
    }

    #[test]
    fn updates_are_rendered_in_insertion_order() {
        let updates = vec![
            ("ref_genome".to_owned(), "/r/b.fa".to_owned()),
            ("gofcards_hg19_fasta".to_owned(), "/r/a.fa".to_owned()),
        ];
        assert_eq!(
            render_config_values("ref_genome: null  # keep\n", &updates).unwrap(),
            "ref_genome: \"/r/b.fa\"  # keep\ngofcards_hg19_fasta: \"/r/a.fa\"\n"
        );
        let duplicated = "ref_genome: a\nref_genome: b\n";
        let failure = render_config_values(duplicated, &updates).unwrap_err();
        assert_eq!(
            failure.to_string(),
            "ValueError: Duplicate reference configuration key: ref_genome"
        );
    }
}
