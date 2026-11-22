//! `due --config CONFIG --source-root ROOT`: whether install step D1 keeps
//! the selected ClinVar gene statistics or builds new ones from the release
//! its download step (D0, `clinvar_release_inputs`) acquired, by the rule of
//! `plan`: an explicitly selected file is kept; the bundled statistics are
//! kept until that release's fileDate reaches the bundled file's
//! modification date plus `threshold_months` calendar months
//! (cache_resources.json, `clinvar_gene_statistics.refresh`), or while the
//! registry disables the automatic refresh; missing or empty bundled
//! statistics are rebuilt.
//!
//! It prints `{"action": "preserve" | "rebuild", "reason", "selected",
//! "file_date"}`. The plan's managed-root parts (the receipt of an earlier
//! refresh, the workspace and its scientific profile) are left out: the step
//! writes into its own folder, and its lineage stands for what it was built
//! from.
use super::policy::resource_policy;
use super::*;
use crate::installer_prerequisites::paths::{exists, is_file};
use crate::python_path::absolute;

pub(super) fn due(config_path: &str, source_root: &str) -> Result<Json> {
    let config = load_config(config_path)?;
    let Some(assembly) = assembly_of(&get_or_null(&config, "assembly")?)? else {
        return Err(value_error(
            "ClinVar refresh requires hg19/GRCh37 or hg38/GRCh38",
        ));
    };
    let baseline = join(source_root, "data/ClinVar/clinvar_2star_stats.pkl.gz");
    let selected = clinvar_stats_path(&config, source_root)?;
    let decision = |action: &str, reason: String, file_date: Json| {
        object([
            ("action", string(action)),
            ("reason", string(reason)),
            ("selected", string(selected.clone())),
            ("file_date", file_date),
        ])
    };
    if absolute(Path::new(&selected))? != absolute(Path::new(&baseline))? {
        nonempty(&selected, "selected ClinVar gene statistics")?;
        return Ok(decision(
            "preserve",
            "explicitly selected statistics are retained".into(),
            Json::Null,
        ));
    }
    let policy = resource_policy(source_root)?;
    let state = baseline_metadata(&baseline)?;
    let missing = !state.get("size").is_some_and(truth);
    let months = match item(&policy, "threshold_months")? {
        Json::Int(months) => *months,
        _ => 0,
    };
    if !truth(item(&policy, "automatic")?) && !missing {
        return Ok(decision(
            "preserve",
            "automatic ClinVar gene-statistics refresh is disabled in the resource registry".into(),
            Json::Null,
        ));
    }
    let release = read_json(&path_of(item(&config, "clinvar_release_inputs")?)?)?;
    require(
        json_text::py_equal(&get_or_null(&release, "assembly")?, &string(assembly)),
        "The ClinVar download step fetched another assembly",
    )?;
    let file_date = item(item(&release, "latest")?, "file_date")?.clone();
    if !missing {
        let stamp = text_of(item(&state, "mtime_utc")?);
        let boundary = add_calendar_months(parse_date(&stamp[..10])?, months)?;
        if parse_date(&text_of(&file_date))? < boundary {
            return Ok(decision(
                "preserve",
                format!("latest ClinVar data date is before the {months}-calendar-month boundary"),
                file_date,
            ));
        }
    }
    let reason = if missing {
        "missing/empty bundled ClinVar gene statistics require managed reconstruction".to_owned()
    } else {
        format!("latest ClinVar data date reaches the {months}-calendar-month boundary")
    };
    Ok(decision("rebuild", reason, file_date))
}

/// `clinvar_stats_path(config, source_root)` of installer_prerequisites.py.
fn clinvar_stats_path(config: &Json, source_root: &str) -> Result<String> {
    let configured = get_or_null(config, "clinvar_gene_stat")?;
    let selected = if truth(&configured) {
        py_str(&configured)
    } else {
        "clinvar_2star_stats.pkl.gz".to_owned()
    };
    let path = pathlib(&selected);
    Ok(if path.starts_with('/') || selected.contains('/') {
        path
    } else {
        join(&join(source_root, "data/ClinVar"), &path)
    })
}

/// `baseline_metadata(path)`: existence, size and mtime only.
fn baseline_metadata(path: &str) -> Result<Json> {
    let present = exists(path)?;
    let mut value = object([("path", string(path)), ("exists", Json::Bool(present))]);
    if !present {
        return Ok(value);
    }
    require(
        is_file(path)?,
        format!("Bundled ClinVar gene statistics is not a file: {path}"),
    )?;
    let info = stat(path)?;
    let (seconds, micros) = float_timestamp(info.mtime(), info.mtime_nsec());
    value.set("size", Json::Int(info.len() as i128));
    value.set(
        "mtime_ns",
        Json::Int(nanoseconds(info.mtime(), info.mtime_nsec())),
    );
    value.set("mtime_utc", string(utc_text(seconds, micros)));
    Ok(value)
}
