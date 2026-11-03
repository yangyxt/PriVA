use super::*;
use crate::python_args::{parse_long_option_occurrences, require_options};
use crate::python_path::path_value;
use std::collections::HashMap;

const HELP: &str = "usage: priva mapped-artifact {run,recover,adopt-native,context}\n\
--artifact FILE --mapping TSV [--input FILE ...] --cache-root DIR --output-dir DIR\n\
--assembly NAME --config YAML --producer SCRIPT --kind {motifs,pm1} --helper-root CHECKOUT\n\
[--annotation-base VCF] [--native-lineage JSON] [--python PYTHON]\n\
run --expected-identity JSON -- PRODUCER_ARGUMENTS\n\
adopt-native --source FILE --lineage JSON_FILE\n\
Optional procedures only. The scientific producer is executed as a literal argument vector.\n\
context prints the original read-only publication binding; no installed caller is changed.";

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(command) = args.first().map(String::as_str) else {
        return Err(invalid(HELP));
    };
    if matches!(command, "-h" | "--help") {
        println!("{HELP}");
        return Ok(());
    }
    if !matches!(command, "run" | "recover" | "adopt-native" | "context") {
        return Err(invalid(HELP));
    }
    let mut names = vec![
        "artifact",
        "mapping",
        "input",
        "cache-root",
        "output-dir",
        "assembly",
        "config",
        "producer",
        "kind",
        "annotation-base",
        "native-lineage",
        "helper-root",
        "python",
    ];
    if command == "run" {
        names.push("expected-identity");
    }
    if command == "adopt-native" {
        names.extend(["source", "lineage"]);
    }
    let arguments = &args[1..];
    let mut split = arguments.len();
    if command == "run" {
        let mut at = 0;
        while at < arguments.len() {
            let token = &arguments[at];
            if token == "--" || !token.starts_with('-') || token == "-" {
                split = at;
                break;
            }
            let bare = token.trim_start_matches('-').split('=').next().unwrap();
            if !token.contains('=')
                && names.iter().any(|n| n.starts_with(bare))
                && !matches!(token.as_str(), "--help" | "-h")
            {
                at += 1;
            }
            at += 1;
        }
    }
    let Some(occurrences) =
        parse_long_option_occurrences(&arguments[..split], &names, &[]).map_err(invalid)?
    else {
        println!("{HELP}");
        return Ok(());
    };
    let inputs = occurrences
        .iter()
        .filter(|(k, _)| k == "input")
        .map(|(_, v)| path_value(v))
        .collect();
    let options: HashMap<_, _> = occurrences.into_iter().collect();
    let mut required = vec![
        "artifact",
        "mapping",
        "cache-root",
        "output-dir",
        "assembly",
        "config",
        "producer",
        "kind",
        "helper-root",
    ];
    if command == "run" {
        required.push("expected-identity");
    }
    if command == "adopt-native" {
        required.extend(["source", "lineage"]);
    }
    require_options(&options, &required).map_err(invalid)?;
    let select = |name: &str| options.get(name).expect("required option");
    if !matches!(select("kind").as_str(), "motifs" | "pm1") {
        return Err(invalid(format!(
            "argument --kind: invalid choice: '{}' (choose from 'motifs', 'pm1')",
            select("kind")
        )));
    }
    crate::builder_error::finish((|| {
        let helper = path_value(select("helper-root"));
        let context = Context {
            artifact: path_value(select("artifact")),
            mapping: path_value(select("mapping")),
            inputs,
            cache_root: path_value(select("cache-root")),
            output_dir: path_value(select("output-dir")),
            assembly: select("assembly").clone(),
            config: path_value(select("config")),
            producer: path_value(select("producer")),
            kind: select("kind").clone(),
            annotation_base: options.get("annotation-base").map(|v| path_value(v)),
            native_lineage: options.get("native-lineage").map(|v| path_value(v)),
            python: options
                .get("python")
                .cloned()
                .unwrap_or_else(|| "python".to_owned()),
            biomart: biomart::Context::new(&helper)?,
            helper,
        };
        match command {
            "context" => println!(
                "{}",
                crate::json_text::dumps(&context.expected()?, true).map_err(invalid)?
            ),
            "recover" => {
                context.recover()?;
            }
            "run" => context.run(select("expected-identity"), &arguments[split..])?,
            "adopt-native" => {
                let encoded = read_text(&path_value(select("lineage")))?;
                context.adopt(&path_value(select("source")), Value::validate(&encoded)?)?;
            }
            _ => unreachable!(),
        }
        Ok(())
    })())
}
