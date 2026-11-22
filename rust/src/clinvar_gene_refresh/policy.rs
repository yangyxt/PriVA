//! The ClinVar sources and refresh rule of the resource registry
//! (`cache_resources.json`, `clinvar_gene_statistics.refresh` and the
//! ClinVar, variant_summary and MANE sources).
use super::*;

const NCBI: &str = "https://ftp.ncbi.nlm.nih.gov";

/// Whether a failure is one of the classes `except (OSError, ValueError,
/// KeyError, TypeError)` catches.
fn policy_failure(category: &str) -> bool {
    matches!(
        category,
        "OSError"
            | "FileNotFoundError"
            | "PermissionError"
            | "IsADirectoryError"
            | "NotADirectoryError"
            | "ValueError"
            | "UnicodeDecodeError"
            | "JSONDecodeError"
            | "KeyError"
            | "TypeError"
    )
}

/// `re.findall(r"\{([^}]+)\}", text)`
fn braced_names(text: &str) -> Vec<&str> {
    let mut names = Vec::new();
    let mut at = 0;
    while let Some(offset) = text[at..].find('{') {
        let start = at + offset;
        match text[start + 1..].find('}') {
            Some(length) if length > 0 => {
                names.push(&text[start + 1..start + 1 + length]);
                at = start + 2 + length;
            }
            _ => at = start + 1,
        }
    }
    names
}

/// `resource_policy(source_root)`: the active source and refresh choices.
pub(super) fn resource_policy(source_root: &str) -> Result<Json> {
    let path = join(source_root, "cache_resources.json");
    let loaded = (|| -> Result<Json> {
        let value = read_json(&path)?;
        if !matches!(value, Json::Object(_))
            || !lookup(&value, "schema_version")?
                .is_some_and(|version| json_text::py_equal(version, &Json::Int(1)))
        {
            return Err(value_error("expected schema_version 1"));
        }
        let resources = item(&value, "resources")?;
        let refresh = item(item(resources, "clinvar_gene_statistics")?, "refresh")?;
        let months = item(refresh, "threshold_months")?;
        let automatic = item(refresh, "automatic")?;
        if !matches!(automatic, Json::Bool(_))
            || lookup(refresh, "policy")?.and_then(Json::as_str)
                != Some("latest_filedate_vs_bundled_mtime")
        {
            return Err(value_error(
                "expected boolean automatic and policy latest_filedate_vs_bundled_mtime",
            ));
        }
        let latest = item(
            item(item(resources, "clinvar_vcf")?, "source")?,
            "latest_url_template",
        )?;
        let summary = item(
            item(item(resources, "clinvar_variant_summary")?, "source")?,
            "url",
        )?;
        let mane = item(item(item(resources, "mane_select")?, "source")?, "url")?;
        if !matches!(months, Json::Int(number) if *number > 0) {
            return Err(value_error("threshold_months must be a positive integer"));
        }
        let template_valid = latest.as_str().is_some_and(|text| {
            text.contains("{assembly}")
                && braced_names(text) == ["assembly"]
                && text.ends_with("/clinvar.vcf.gz")
        });
        if !template_valid {
            return Err(value_error(
                "latest_url_template must contain {assembly} and end in /clinvar.vcf.gz",
            ));
        }
        for url in [latest, summary, mane] {
            if !url
                .as_str()
                .is_some_and(|text| text.starts_with(&format!("{NCBI}/")))
            {
                return Err(value_error(
                    "ClinVar/MANE source URLs must use official NCBI HTTPS",
                ));
            }
        }
        Ok(object([
            ("threshold_months", months.clone()),
            ("automatic", automatic.clone()),
            ("latest_url_template", latest.clone()),
            ("summary_url", summary.clone()),
            ("mane_url", mane.clone()),
        ]))
    })();
    let fields = match loaded {
        Ok(fields) => fields,
        Err(failure) => {
            let (category, message) = exception(failure.as_ref());
            if !policy_failure(category) {
                return Err(failure);
            }
            return Err(value_error(format!(
                "Cannot load ClinVar refresh policy from {path}: {message}"
            )));
        }
    };
    let mut policy = object([
        ("registry", string(path.clone())),
        ("registry_sha256", string(sha256(&path)?)),
    ]);
    if let Json::Object(items) = fields {
        for (key, value) in items {
            policy.set(&key, value);
        }
    }
    Ok(policy)
}
