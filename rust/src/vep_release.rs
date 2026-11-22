//! `priva vep-release pins --packages FILE`: the exact package pins that keep
//! every other Conda package of the PriVA environment fixed while
//! `conda_install_vep` changes its VEP release with the package manager's own
//! `install --file PINS ensembl-vep=RELEASE.*`. The inline Python block it
//! replaces read the same file.
//!
//! FILE holds the manager's `list --json` output (a list of package records,
//! or an object with a `packages` list). It is replaced by one
//! `NAME=VERSION=BUILD` line per package, sorted by name, leaving out
//! `ensembl-vep` itself and packages installed by pip (channel or platform
//! `pypi`, build `pypi_0`). The build is `build_string`, or `build` when that
//! is empty. A name, version or build outside `[A-Za-z0-9_.+!-]+`, a package
//! listed twice, or nothing left to pin stops with
//! `ERROR: Cannot preserve installed Conda packages while changing VEP: ...`
//! and exit code 1, before the environment is touched.
use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use std::collections::BTreeMap;
use std::fs;

const USAGE: &str = "usage: priva vep-release pins --packages FILE\n\
Replace the package manager's `list --json` output in FILE with exact NAME=VERSION=BUILD pins of every package except ensembl-vep and pip packages.";

/// `str()` of a JSON value as Python prints it.
fn python_text(value: &Json) -> String {
    match value {
        Json::Null => "None".to_owned(),
        Json::Bool(true) => "True".to_owned(),
        Json::Bool(false) => "False".to_owned(),
        Json::Str(text) => text.clone(),
        other => json_text::dumps(other, false).unwrap_or_default(),
    }
}

/// Python truthiness of a JSON value.
fn truthy(value: Option<&Json>) -> bool {
    match value {
        None | Some(Json::Null) | Some(Json::Bool(false)) => false,
        Some(Json::Int(n)) => *n != 0,
        Some(Json::Float(f)) => *f != 0.0,
        Some(Json::Str(s)) => !s.is_empty(),
        Some(Json::Array(items)) => !items.is_empty(),
        Some(Json::Object(items)) => !items.is_empty(),
        Some(Json::Bool(true)) => true,
    }
}

fn pinnable(value: Option<&Json>) -> bool {
    matches!(value, Some(Json::Str(text)) if !text.is_empty()
        && text.bytes().all(|b| b.is_ascii_alphanumeric() || b"_.+!-".contains(&b)))
}

pub(crate) fn pins(listing: &str) -> std::result::Result<String, String> {
    let mut records = json_text::parse(listing, false)?;
    if let Json::Object(_) = records {
        records = records.get("packages").cloned().unwrap_or(Json::Null);
    }
    let records = match records {
        Json::Array(items) if !items.is_empty() => items,
        _ => return Err("package list must be a nonempty JSON array".to_owned()),
    };
    let mut pins = BTreeMap::new();
    for record in &records {
        if !matches!(record, Json::Object(_)) {
            return Err("installed package records must be JSON objects".to_owned());
        }
        let name = record.get("name").ok_or("'name'")?;
        let build = if truthy(record.get("build_string")) {
            record.get("build_string")
        } else {
            record.get("build")
        };
        let is = |value: Option<&Json>, wanted: &str| value.and_then(Json::as_str) == Some(wanted);
        if is(record.get("channel"), "pypi")
            || is(record.get("platform"), "pypi")
            || is(build, "pypi_0")
            || is(Some(name), "ensembl-vep")
        {
            continue;
        }
        let version = record.get("version").ok_or("'version'")?;
        let shown = python_text(name);
        if !(pinnable(Some(name)) && pinnable(Some(version)) && pinnable(build)) {
            return Err(format!("cannot pin installed package {shown}"));
        }
        if pins.contains_key(&shown) {
            return Err(format!("duplicate installed package {shown}"));
        }
        let line = format!(
            "{shown}={}={}",
            python_text(version),
            build.map(python_text).unwrap_or_default()
        );
        pins.insert(shown, line);
    }
    if pins.is_empty() {
        return Err("no installed Conda dependencies could be pinned".to_owned());
    }
    let mut out = pins.into_values().collect::<Vec<_>>().join("\n");
    out.push('\n');
    Ok(out)
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let path = match args {
        [mode, option, path] if mode == "pins" && option == "--packages" => path,
        [flag, ..] if flag == "--help" || flag == "-h" => {
            println!("{USAGE}");
            return Ok(());
        }
        _ => return Err(invalid(USAGE)),
    };
    let outcome = fs::read_to_string(path)
        .map_err(|error| format!("{error}: {path}"))
        .and_then(|listing| pins(&listing))
        .and_then(|text| fs::write(path, text).map_err(|error| format!("{error}: {path}")));
    if let Err(message) = outcome {
        eprintln!("ERROR: Cannot preserve installed Conda packages while changing VEP: {message}");
        std::process::exit(1);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::pins;

    #[test]
    fn exact_pins_without_vep_or_pip() {
        let listing = r#"[
            {"name": "zlib", "version": "1.3.1", "build_string": "hb9d3cd8_2", "channel": "conda-forge"},
            {"name": "ensembl-vep", "version": "113.2", "build_string": "pl5321h2a3209d_0"},
            {"name": "diptest", "version": "0.8.1", "build_string": "pypi_0", "channel": "pypi"},
            {"name": "bcftools", "version": "1.21", "build": "h3a4d415_1", "build_string": ""}
        ]"#;
        assert_eq!(
            pins(listing).unwrap(),
            "bcftools=1.21=h3a4d415_1\nzlib=1.3.1=hb9d3cd8_2\n"
        );
        assert_eq!(
            pins(r#"{"packages": [{"name": "a b", "version": "1", "build": "0"}]}"#),
            Err("cannot pin installed package a b".to_owned())
        );
        assert_eq!(
            pins("[]"),
            Err("package list must be a nonempty JSON array".to_owned())
        );
        assert_eq!(pins(r#"[{"version": "1"}]"#), Err("'name'".to_owned()));
        assert_eq!(
            pins(r#"[{"name": "ensembl-vep", "version": "1", "build": "0"}]"#),
            Err("no installed Conda dependencies could be pinned".to_owned())
        );
    }
}
