//! Native `scripts/gnomad_cache.py`: `members`, `validate` and `default-path`
//! over the selected indexed gnomAD chromosome family. Read-only: no download,
//! index creation or byte change.
//!
//! The checks are those of `installer_config_paths::gnomad` (its source
//! adoption uses the same procedures). Printed lines and exit statuses follow
//! the script; a failure that main() printed as `ERROR: message` exits 1 here
//! as `priva: CATEGORY: message`. Differences: the subcommand comes first, and
//! headers and the first indexed record are read by the native BGZF and
//! tabix/CSI readers instead of htslib through pysam.
use crate::builder_error;
use crate::cli::{Result, invalid};
use crate::installer_config_paths::gnomad::{default_path, family_paths, validate_selected};
use crate::installer_prerequisites::paths::pathlib;
use crate::json_text::{self, Json};
use crate::python_args::{parse_long_option_occurrences, require_options};
use crate::python_path::repr;
use std::collections::HashMap;
use std::path::Path;

const HELP: &str = "usage: priva gnomad-cache {members,validate,default-path}\n\
       members --chr-x CHR_X\n\
       validate --chr-x CHR_X --assembly ASSEMBLY [--reference REFERENCE]\n\
       default-path --cache-dir CACHE_DIR --assembly ASSEMBLY\n\
Inspect the exact indexed gnomAD chromosome family consumed by PriVA; the output and exit\n\
statuses of scripts/gnomad_cache.py. The subcommand comes first.";
const COMMANDS: [&str; 3] = ["members", "validate", "default-path"];

fn options(command: &str) -> (&'static [&'static str], &'static [&'static str]) {
    match command {
        "members" => (&["chr-x"], &["chr-x"]),
        "validate" => (&["chr-x", "assembly", "reference"], &["chr-x", "assembly"]),
        _ => (&["cache-dir", "assembly"], &["cache-dir", "assembly"]),
    }
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
    let (names, required) = options(command);
    let Some(occurrences) =
        parse_long_option_occurrences(&args[1..], names, &[]).map_err(invalid)?
    else {
        println!("{HELP}");
        return Ok(());
    };
    let values: HashMap<String, String> = occurrences.into_iter().collect();
    require_options(&values, required).map_err(invalid)?;
    let path = |name: &str| values.get(name).map(|value| pathlib(value));
    builder_error::finish(run(
        command,
        path("chr-x").as_deref(),
        values.get("assembly").map(String::as_str),
        path("reference").as_deref(),
        path("cache-dir").as_deref(),
    ))
}

fn run(
    command: &str,
    chr_x: Option<&str>,
    assembly: Option<&str>,
    reference: Option<&str>,
    cache_dir: Option<&str>,
) -> Result<()> {
    let printed = match command {
        "members" => {
            let members = family_paths(Path::new(chr_x.unwrap_or_default()))?
                .into_iter()
                .map(|(name, path)| (name, Json::Str(path.display().to_string())))
                .collect();
            json_text::dumps(&Json::Object(members), false).map_err(invalid)?
        }
        "default-path" => default_path(
            Path::new(cache_dir.unwrap_or_default()),
            assembly.unwrap_or_default(),
        )?
        .display()
        .to_string(),
        _ => {
            let result = validate_selected(
                Path::new(chr_x.unwrap_or_default()),
                assembly.unwrap_or_default(),
                reference,
            )?;
            json_text::dumps(&result, false).map_err(invalid)?
        }
    };
    println!("{printed}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(items: &[&str]) -> Vec<String> {
        items.iter().map(|item| (*item).to_owned()).collect()
    }

    #[test]
    fn admission_errors_exit_two_and_checks_exit_one() {
        for bad in [
            args(&[]),
            args(&["inspect"]),
            args(&["members"]),
            args(&["validate", "--chr-x", "/g/a.chrX.vcf.bgz"]),
            args(&[
                "default-path",
                "--cache-dir",
                "/g",
                "--assembly",
                "hg38",
                "--reference",
                "/r",
            ]),
        ] {
            let failure = dispatch(&bad).unwrap_err();
            assert!(!failure.is::<builder_error::RuntimeError>(), "{failure}");
        }
        let failure = dispatch(&args(&[
            "default-path",
            "--cache-dir",
            "/g",
            "--assembly",
            "hg39",
        ]))
        .unwrap_err();
        assert_eq!(
            failure.to_string(),
            "ValueError: Select hg19/GRCh37 or hg38/GRCh38"
        );
        let failure = dispatch(&args(&["members", "--chr-x", "/g/a.chr1.vcf.bgz"])).unwrap_err();
        assert!(failure.to_string().starts_with(
            "ValueError: Selected gnomAD filename must retain the caller's .chrX.vcf.gz"
        ));
    }
}
