//! preflight: priva's checks of what an installation needs, run once after
//! the plan and before any step (`INSTALL_PREFLIGHT` in
//! `modules/local/install.nf`), so a missing publisher file stops the run at
//! once instead of after hours of downloads, as `main_install` stops before
//! its first resource:
//!
//! - each file you supply (registry `user_inputs`) for a step that will run:
//!   a set path must be an absolute path to a completed, nonempty, readable
//!   file; the SpliceAI and PrimateAI publisher files must be set, since no
//!   step can download them. Steps taken from an older installation
//!   (`--adopt`) are left out.
//! - the data files shipped in the checkout, with this program's own
//!   `installer-prerequisites --stage snapshots`.
//!
//! Both read the plan's configuration (`base_config`). With `--stub` only the
//! shipped files are checked, since a stub run supplies no publisher files.
//!
//! The file rules are those of `installer-prerequisites --stage
//! manual-inputs`, kept here rather than shared: that module is reached by the
//! code digest of 37 install steps (all of the ClinVar VEP chain among them),
//! and no step runs `install-task`, so this check can change without
//! rebuilding any resource.
use super::{Options, read_json, steps, strings, text, yaml_value};
use crate::cli::{Result, invalid};
use crate::json_text::Json;
use std::fs::File;
use std::path::{Component, Path};
use std::process::Command;

/// Publisher files no step can download (`MANUAL_INPUT_KEYS` of
/// `installer-prerequisites`).
const PUBLISHER_FILES: [&str; 3] = [
    "primateai_prescore",
    "spliceai_snv_prescore",
    "spliceai_indel_prescore",
];

/// The steps that will run and the files they need from you, in registry order.
fn running(registry: &Json, plan: &Json) -> Result<(Vec<String>, Vec<String>)> {
    let planned = plan
        .get("steps")
        .ok_or_else(|| invalid("install-task preflight: the plan has no steps"))?;
    let mut ids = Vec::new();
    let mut keys: Vec<String> = Vec::new();
    for entry in steps(registry)? {
        let id = text(entry, "id");
        let Some(step) = planned.get(&id) else {
            continue;
        };
        if !matches!(step.get("adopt"), None | Some(Json::Null)) {
            continue;
        }
        for key in strings(entry, "user_inputs") {
            if !keys.contains(&key) {
                keys.push(key);
            }
        }
        ids.push(id);
    }
    Ok((ids, keys))
}

/// A path part that names an unfinished or staging copy.
fn unfinished(part: &str) -> bool {
    let part = part.to_lowercase();
    part.contains(".partial")
        || part == ".building"
        || part.starts_with(".building-")
        || part.starts_with(".download.")
        || [".part", ".crdownload", ".download", ".copying", ".tmp"]
            .iter()
            .any(|ending| part.ends_with(ending))
}

fn any_unfinished(path: &Path) -> bool {
    path.components().any(|part| match part {
        Component::Normal(name) => unfinished(&name.to_string_lossy()),
        _ => false,
    })
}

/// One file you supply, as the configuration names it.
fn supplied_file(key: &str, value: &str) -> std::result::Result<(), String> {
    let path = Path::new(value);
    if !path.is_absolute() {
        return Err(format!(
            "{key} must be an absolute path, since each install step runs in its own folder: {value}"
        ));
    }
    let usable = std::fs::metadata(path).is_ok_and(|m| m.is_file() && m.len() > 0)
        && File::open(path).is_ok();
    if !usable {
        return Err(format!("Missing, empty, or unreadable {key}: {value}"));
    }
    // A final-looking link must not conceal an unfinished source.
    let resolved = std::fs::canonicalize(path).map_err(|error| format!("{error}: {value}"))?;
    if any_unfinished(path) || any_unfinished(&resolved) {
        return Err(format!(
            "Incomplete/staging manual input {key}: {value}; select the completed original file"
        ));
    }
    Ok(())
}

/// The files you supply for the steps that will run.
fn supplied_files(config: &Path, keys: &[String]) -> Result<()> {
    let mut missing = Vec::new();
    for key in keys {
        let value = yaml_value(config, key)?;
        if matches!(value.trim(), "" | "null") {
            if PUBLISHER_FILES.contains(&key.as_str()) {
                missing.push(key.as_str());
            }
            continue;
        }
        supplied_file(key, &value).map_err(invalid)?;
    }
    if missing.is_empty() {
        return Ok(());
    }
    Err(invalid(format!(
        "Manual publisher inputs must be prepared before resource installation: {}; \
         set these paths in the install configuration. No install step was started",
        missing.join(", ")
    )))
}

pub(crate) fn preflight(options: &Options) -> Result<()> {
    let registry = read_json(&options.path("registry")?)?;
    let plan = read_json(&options.path("plan")?)?;
    let root = options.source_root()?;
    let (ids, keys) = running(&registry, &plan)?;
    let config = text(&plan, "base_config");
    if config.is_empty() {
        return Err(invalid(
            "install-task preflight: the plan names no base_config",
        ));
    }
    let config = Path::new(&config);
    if !options.stub {
        supplied_files(config, &keys)
            .map_err(|failure| invalid(format!("install-task preflight: {failure}")))?;
    }
    let program = std::env::current_exe()
        .map_err(|error| invalid(format!("install-task preflight: {error}")))?;
    let status = Command::new(program)
        .arg("installer-prerequisites")
        .arg(config)
        .arg("--source-root")
        .arg(&root)
        .args(["--stage", "snapshots"])
        .status()
        .map_err(|error| invalid(format!("install-task preflight: {error}")))?;
    if !status.success() {
        return Err(invalid(
            "install-task preflight: data files shipped in the checkout are missing or empty (see above)",
        ));
    }
    let supplied = if options.stub {
        "files you supply not checked (stub)".to_owned()
    } else if keys.is_empty() {
        "no files you supply".to_owned()
    } else {
        format!("files you supply checked ({})", keys.join(", "))
    };
    println!(
        "preflight: {} steps to run; {supplied}; shipped data files present",
        ids.len()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{running, supplied_file, unfinished};
    use crate::json_text::parse;

    #[test]
    fn adopted_steps_need_no_files_from_you() {
        let registry = parse(
            r#"{"steps": [
                {"id": "A1", "user_inputs": []},
                {"id": "B4", "user_inputs": ["spliceai_snv_prescore", "spliceai_indel_prescore"]},
                {"id": "B7", "user_inputs": ["primateai_prescore"]},
                {"id": "F5", "user_inputs": ["gofcards_hg19_fasta", "spliceai_snv_prescore"]}
            ]}"#,
            false,
        )
        .unwrap();
        let plan = parse(
            r#"{"steps": {"A1": {"adopt": null}, "B4": {"adopt": null},
                "B7": {"adopt": {"identity": "x"}}, "F5": {}}}"#,
            false,
        )
        .unwrap();
        let (ids, keys) = running(&registry, &plan).unwrap();
        assert_eq!(ids, ["A1", "B4", "F5"]);
        assert_eq!(
            keys,
            [
                "spliceai_snv_prescore",
                "spliceai_indel_prescore",
                "gofcards_hg19_fasta"
            ]
        );
    }

    #[test]
    fn a_supplied_file_is_absolute_complete_and_nonempty() {
        let folder = std::env::temp_dir().join(format!("priva-preflight-{}", std::process::id()));
        std::fs::create_dir_all(&folder).unwrap();
        let good = folder.join("snv.vcf.gz");
        std::fs::write(&good, b"x").unwrap();
        let partial = folder.join("snv.vcf.gz.partial");
        std::fs::write(&partial, b"x").unwrap();
        let empty = folder.join("empty.vcf.gz");
        std::fs::write(&empty, b"").unwrap();
        let link = folder.join("final.vcf.gz");
        std::os::unix::fs::symlink(&partial, &link).unwrap();
        let text = |path: &std::path::Path| path.display().to_string();
        assert_eq!(supplied_file("k", &text(&good)), Ok(()));
        assert!(
            supplied_file("k", "relative.vcf.gz")
                .unwrap_err()
                .contains("absolute")
        );
        assert!(
            supplied_file("k", &text(&empty))
                .unwrap_err()
                .starts_with("Missing, empty")
        );
        assert!(
            supplied_file("k", &text(&partial))
                .unwrap_err()
                .starts_with("Incomplete")
        );
        assert!(
            supplied_file("k", &text(&link))
                .unwrap_err()
                .starts_with("Incomplete")
        );
        std::fs::remove_dir_all(&folder).unwrap();
        assert!(unfinished("X.CRDOWNLOAD") && unfinished(".building-1") && !unfinished("a.vcf.gz"));
    }
}
