//! Native installer input checks, from `scripts/installer_prerequisites.py`.
//!
//! The stages read metadata, small indexes and headers only. They never create,
//! index or download a resource. Only `--stage manual-inputs --prompt` at a
//! terminal writes the configuration, through the native `render_config` of
//! installer_config_paths.py. Printed paths, exit statuses and the Python
//! exception category and message of each failure follow the script.
//! `--allow-pinned-reference-defaults` plans and checks the references with
//! the native `priva reference-resources` procedures, as the script did with
//! reference_resources.py. The script's default `--source-root` is its own
//! checkout; this command has no such location, so the option is required.
pub(crate) mod checks;
mod mapping;
pub(crate) mod paths;

use crate::builder_error::{self, RuntimeError};
use crate::cli::{Result, invalid};
use crate::python_args::{parse_long_option_occurrences, require_options};
use crate::python_path::repr;
use crate::yaml_config::{self, ConfigValue, Parsed, Scalar, YamlError};
use paths::pathlib;
use std::collections::HashMap;
use std::error::Error;

const HELP: &str = "usage: priva installer-prerequisites CONFIG --source-root CHECKOUT\n\
       [--stage STAGE] [--prompt] [--print-mapping] [--allow-missing-bundled-mapping]\n\
       [--allow-missing-bundled-clinvar] [--allow-pinned-reference-defaults]\n\
       [--allow-missing-merged-caches] [--vep-version RELEASE]\n\
STAGE: main (default), manual-inputs, snapshots, mapping, mapping-target, mapping-source,\n\
mapping-workspace, clinvar, gofcards, gofcards-main.\n\
The checks, printed paths and exit statuses of scripts/installer_prerequisites.py.\n\
CONFIG comes first; --source-root names the PriVA checkout.";
const STAGES: [&str; 10] = [
    "main",
    "manual-inputs",
    "snapshots",
    "mapping",
    "mapping-target",
    "mapping-source",
    "mapping-workspace",
    "clinvar",
    "gofcards",
    "gofcards-main",
];
const NAMES: [&str; 3] = ["source-root", "stage", "vep-version"];
const FLAGS: [&str; 6] = [
    "prompt",
    "print-mapping",
    "allow-missing-bundled-mapping",
    "allow-missing-bundled-clinvar",
    "allow-pinned-reference-defaults",
    "allow-missing-merged-caches",
];

fn error(category: &'static str, message: impl Into<String>) -> Box<dyn Error> {
    builder_error::category(category, invalid(message))
}
fn value_error(message: impl Into<String>) -> Box<dyn Error> {
    error("ValueError", message)
}
fn not_native(message: impl Into<String>) -> Box<dyn Error> {
    error("NotImplementedError", message)
}
fn require(condition: bool, message: impl Into<String>) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(value_error(message))
    }
}
/// The Python category and message of a failure raised in this module.
fn exception(failure: &(dyn Error + 'static)) -> Option<(&'static str, String)> {
    failure
        .downcast_ref::<RuntimeError>()
        .map(|known| (known.category_name(), known.message()))
}

/// A value `yaml.safe_load` gives for one configuration key.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Value {
    None,
    Bool(bool),
    /// A decimal integer as Python prints it.
    Int(String),
    /// A float, an integer in another notation, a date or a timestamp, as written.
    Other(String),
    Str(String),
    /// A nested block or a flow collection, as written.
    Collection(String),
}

static NONE: Value = Value::None;

impl Value {
    fn of(value: &ConfigValue) -> Result<Self> {
        match value {
            ConfigValue::Scalar(scalar) => Self::scalar(scalar),
            ConfigValue::Opaque(text) => Ok(Self::Collection(text.clone())),
        }
    }
    fn scalar(scalar: &Scalar) -> Result<Self> {
        Ok(match scalar {
            Scalar::Null => Self::None,
            Scalar::Bool(value) => Self::Bool(*value),
            Scalar::Int(text) => Self::Int(text.clone()),
            // safe_load has no constructor for the merge and value tags.
            Scalar::Other(text) if text == "<<" || text == "=" => {
                return Err(error(
                    "YAMLError",
                    format!("could not determine a constructor for the plain scalar {text}"),
                ));
            }
            Scalar::Other(text) => Self::Other(text.clone()),
            Scalar::Str(text) => Self::Str(text.clone()),
        })
    }
    /// `value in (None, "", "null")`
    pub(crate) fn is_unset(&self) -> bool {
        match self {
            Self::None => true,
            Self::Str(text) => text.is_empty() || text == "null",
            _ => false,
        }
    }
    pub(super) fn text(&self) -> Option<&str> {
        match self {
            Self::Str(text) => Some(text),
            _ => None,
        }
    }
    /// `type(value).__name__`
    pub(super) fn type_name(&self) -> &'static str {
        match self {
            Self::None => "NoneType",
            Self::Bool(_) => "bool",
            Self::Int(_) => "int",
            Self::Other(text) => Scalar::Other(text.clone()).python_type_name(),
            Self::Str(_) => "str",
            Self::Collection(text) => ConfigValue::Opaque(text.clone()).python_type_name(),
        }
    }
    /// `bool(value)`
    pub(super) fn truthy(&self) -> bool {
        match self {
            Self::None => false,
            Self::Bool(value) => *value,
            Self::Int(text) => text != "0",
            Self::Other(text) => {
                let scalar = Scalar::Other(text.clone());
                !(matches!(scalar.python_type_name(), "int" | "float")
                    && matches!(scalar.python_repr().as_str(), "0" | "0.0" | "-0.0"))
            }
            Self::Str(text) => !text.is_empty(),
            Self::Collection(text) => !matches!(
                ConfigValue::Opaque(text.clone()).python_repr().as_str(),
                "[]" | "{}"
            ),
        }
    }
    /// `str(value)`
    pub(super) fn py_str(&self) -> Result<String> {
        Ok(match self {
            Self::None => "None".to_owned(),
            Self::Bool(value) => (if *value { "True" } else { "False" }).to_owned(),
            Self::Int(text) | Self::Str(text) => text.clone(),
            Self::Other(text) => {
                let scalar = Scalar::Other(text.clone());
                match scalar.python_type_name() {
                    // The YAML date form is exactly date.isoformat().
                    "date" => text.clone(),
                    "datetime" => {
                        return Err(not_native(format!(
                            "a YAML timestamp is not used as configuration text: {text}"
                        )));
                    }
                    _ => scalar.python_repr(),
                }
            }
            Self::Collection(text) => ConfigValue::Opaque(text.clone()).python_repr(),
        })
    }
}

/// The top-level mapping: each key once, at its first position, with its last value.
pub(crate) struct Config {
    entries: Vec<(String, Value)>,
}

impl Config {
    /// `config.get(key)`
    pub(crate) fn get(&self, key: &str) -> &Value {
        self.entries
            .iter()
            .find(|(name, _)| name == key)
            .map_or(&NONE, |(_, value)| value)
    }
    pub(super) fn entries(&self) -> &[(String, Value)] {
        &self.entries
    }
    fn set(&mut self, key: &str, value: Value) {
        match self.entries.iter_mut().find(|(name, _)| name == key) {
            Some(entry) => entry.1 = value,
            None => self.entries.push((key.to_owned(), value)),
        }
    }
    /// `yaml.safe_load` of configuration text; `unique` adds the duplicate-key
    /// refusal of installer_config_paths.UniqueLoader, nested mappings included.
    pub(crate) fn load(text: &str, unique: bool) -> Result<Self> {
        let document = match yaml_config::parse_config(text) {
            Ok(Parsed::Mapping(document)) => document,
            Ok(Parsed::Empty) if empty_mapping(text) => {
                return Ok(Self {
                    entries: Vec::new(),
                });
            }
            Ok(Parsed::Empty | Parsed::List) => {
                return Err(value_error("The configuration must be a YAML mapping"));
            }
            Err(failure @ YamlError::Invalid(_)) => {
                return Err(error("YAMLError", failure.to_string()));
            }
            Err(failure @ YamlError::Unsupported(_)) => {
                return Err(not_native(failure.to_string()));
            }
        };
        let mut config = Self {
            entries: Vec::new(),
        };
        for entry in &document.entries {
            if unique {
                require(
                    config.entries.iter().all(|(name, _)| *name != entry.key),
                    format!("Duplicate configuration key: {}", entry.key),
                )?;
                if let ConfigValue::Opaque(text) = &entry.value
                    && let Some(nested) = yaml_config::nested_mapping_entries(text)
                {
                    for (index, (key, _)) in nested.iter().enumerate() {
                        require(
                            nested[..index].iter().all(|(earlier, _)| earlier != key),
                            format!("Duplicate configuration key: {key}"),
                        )?;
                    }
                }
            }
            config.set(&entry.key, Value::of(&entry.value)?);
        }
        Ok(config)
    }
}

/// Whether text the reader calls empty is `{}` (a mapping) rather than nothing or `[]`.
fn empty_mapping(text: &str) -> bool {
    let body = text.strip_prefix('\u{feff}').unwrap_or(text);
    let mut lines = body
        .split('\n')
        .map(|line| line.strip_suffix('\r').unwrap_or(line))
        .filter(|line| {
            let content = line.trim_start_matches(' ');
            !content.is_empty() && !content.starts_with('#')
        });
    let mut first = lines.next();
    if first.is_some_and(|line| {
        let content = line.trim_end_matches(' ');
        content == "---"
            || (content.starts_with("--- ")
                && content[4..].trim_start_matches(' ').starts_with('#'))
    }) {
        first = lines.next();
    }
    first.is_some_and(|line| line.trim_end() == "{}")
}

struct Options {
    stage: String,
    source_root: String,
    vep_version: Option<String>,
    flags: Vec<String>,
}

impl Options {
    fn flag(&self, name: &str) -> bool {
        self.flags.iter().any(|flag| flag == name)
    }
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args
        .first()
        .is_some_and(|arg| matches!(arg.as_str(), "-h" | "--help"))
    {
        println!("{HELP}");
        return Ok(());
    }
    let config = match args.first() {
        Some(arg) if !arg.starts_with('-') || arg == "-" => arg,
        _ => {
            return Err(invalid(format!(
                "the following arguments are required: config (it comes first)\n{HELP}"
            )));
        }
    };
    let Some(occurrences) =
        parse_long_option_occurrences(&args[1..], &NAMES, &FLAGS).map_err(invalid)?
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
    require_options(&values, &["source-root"]).map_err(invalid)?;
    let stage = values.remove("stage").unwrap_or_else(|| "main".to_owned());
    if !STAGES.contains(&stage.as_str()) {
        return Err(invalid(format!(
            "argument --stage: invalid choice: {} (choose from {})",
            repr(&stage),
            STAGES.map(repr).join(", ")
        )));
    }
    let options = Options {
        stage,
        source_root: pathlib(&values["source-root"]),
        vep_version: values.remove("vep-version"),
        flags,
    };
    builder_error::finish(run(&pathlib(config), &options))
}

fn run(config_path: &str, options: &Options) -> Result<()> {
    let stage = options.stage.as_str();
    let gofcards = matches!(stage, "gofcards" | "gofcards-main");
    require(
        !options.flag("prompt") || stage == "manual-inputs",
        "--prompt is restricted to --stage manual-inputs",
    )?;
    require(
        !options.flag("allow-missing-bundled-clinvar") || stage == "main",
        "--allow-missing-bundled-clinvar is restricted to main preflight",
    )?;
    require(
        !options.flag("allow-pinned-reference-defaults") || gofcards,
        "--allow-pinned-reference-defaults is restricted to GoFCards stages",
    )?;
    require(
        !options.flag("allow-missing-merged-caches") || gofcards,
        "--allow-missing-merged-caches is restricted to GoFCards stages",
    )?;
    if stage == "manual-inputs" {
        return checks::manual_inputs(config_path, options.flag("prompt"));
    }
    let mut config = Config::load(&paths::read_text(config_path)?, false)?;
    if let Some(release) = &options.vep_version {
        config.set("vep_installed_version", Value::Str(release.clone()));
    }
    let root = options.source_root.as_str();
    match stage {
        "main" => checks::validate_main(
            &config,
            root,
            options.flag("allow-missing-bundled-mapping"),
            options.flag("allow-missing-bundled-clinvar"),
        ),
        "snapshots" => checks::validate_bundled_snapshots(&config, root),
        "mapping" => {
            let path = mapping::validate_mapping(&config, root)?;
            if options.flag("print-mapping") {
                println!("{path}");
            }
            Ok(())
        }
        "mapping-target" => {
            println!("{}", mapping::nm_to_enst_path(&config, root)?);
            Ok(())
        }
        "mapping-source" => {
            println!("{}", mapping::mapping_fetch_source(&config, root)?);
            Ok(())
        }
        "mapping-workspace" => {
            println!("{}", mapping::mapping_workspace(&config, root)?);
            Ok(())
        }
        "clinvar" => checks::validate_selected_clinvar_stats(&config, root),
        _ => checks::validate_gofcards(
            &config,
            root,
            stage == "gofcards-main",
            options.flag("allow-pinned-reference-defaults"),
            options.flag("allow-missing-merged-caches"),
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn load(text: &str) -> Result<Config> {
        Config::load(text, false)
    }

    #[test]
    fn values_follow_python_truth_and_text() {
        let config = load(
            "a: 0\nb: 0.0\nc: ''\nd: '0'\ne: []\nf: [x]\ng: 2024-01-02\nh: yes\ni: ~\nj: 0x10\n",
        )
        .unwrap();
        for (key, truth) in [
            ("a", false),
            ("b", false),
            ("c", false),
            ("d", true),
            ("e", false),
            ("f", true),
            ("g", true),
            ("h", true),
            ("i", false),
            ("j", true),
        ] {
            assert_eq!(config.get(key).truthy(), truth, "{key}");
        }
        assert_eq!(config.get("g").py_str().unwrap(), "2024-01-02");
        assert_eq!(config.get("h").py_str().unwrap(), "True");
        assert_eq!(config.get("j").py_str().unwrap(), "16");
        assert_eq!(config.get("f").py_str().unwrap(), "['x']");
        assert_eq!(config.get("i").py_str().unwrap(), "None");
        assert_eq!(config.get("absent"), &Value::None);
        assert!(config.get("c").is_unset() && !config.get("d").is_unset());
        assert!(Value::Str("null".into()).is_unset());
    }

    #[test]
    fn duplicate_keys_keep_first_position_and_last_value_unless_unique() {
        let text = "a: 1\nb: 2\na: 3\n";
        let config = load(text).unwrap();
        assert_eq!(config.entries().len(), 2);
        assert_eq!(
            config.entries()[0],
            ("a".to_owned(), Value::Int("3".into()))
        );
        let failure = Config::load(text, true).err().unwrap();
        assert_eq!(
            exception(failure.as_ref()).unwrap(),
            ("ValueError", "Duplicate configuration key: a".to_owned())
        );
        let nested = Config::load("cache_folders:\n  x: 1\n  x: 2\n", true).err();
        assert!(nested.is_some());
    }

    #[test]
    fn empty_documents_are_mappings_only_when_written_as_braces() {
        assert!(load("{}\n").unwrap().entries().is_empty());
        assert!(load("---\n{}\n").unwrap().entries().is_empty());
        for text in ["", "# only a comment\n", "[]\n", "- a\n"] {
            let failure = load(text).err().unwrap();
            assert_eq!(
                exception(failure.as_ref()).unwrap().1,
                "The configuration must be a YAML mapping"
            );
        }
        let failure = load("a: &anchor x\n").err().unwrap();
        assert_eq!(
            exception(failure.as_ref()).unwrap().0,
            "NotImplementedError"
        );
        let failure = load("a: <<\n").err().unwrap();
        assert_eq!(exception(failure.as_ref()).unwrap().0, "YAMLError");
    }

    #[test]
    fn admission_keeps_argparse_exit_two_and_the_restrictions_exit_one() {
        let args = |items: &[&str]| items.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>();
        for bad in [
            args(&[]),
            args(&["--stage", "main"]),
            args(&["config.yaml"]),
            args(&["config.yaml", "--source-root", "x", "--stage", "other"]),
            args(&["config.yaml", "--source-root", "x", "--prompt=yes"]),
        ] {
            let failure = dispatch(&bad).err().unwrap();
            assert!(!failure.is::<RuntimeError>(), "{failure}");
        }
        let failure = dispatch(&args(&[
            "missing-config.yaml",
            "--source-root",
            "x",
            "--stage",
            "mapping",
            "--prompt",
        ]))
        .err()
        .unwrap();
        assert_eq!(
            exception(failure.as_ref()).unwrap(),
            (
                "ValueError",
                "--prompt is restricted to --stage manual-inputs".to_owned()
            )
        );
    }
}
