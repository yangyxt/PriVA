//! `priva mondo-scope validate --registry TSV[.gz] --hpo-assertions TSV[.gz]`:
//! the check of the disease-scope registry and the scoped HPO assertions that
//! `mondo_hpo_scope_install` ran as an inline Python block
//! (`validate_mondo_hpo_scope_cache` before this command; the shell function
//! still checks the MONDO OBO header lines). Each table, read as gzip when its
//! name ends in `.gz`, must have the columns PriVA reads and at least one row.
//! Success prints `disease_scope_rows=N; scoped_hpo_assertions=M`; a failure
//! prints the reason and exits with 1. Compared with the Python block by
//! `rust/tests/install_checks_compare.py`.
use crate::cli::{Result, invalid};
use crate::tsv_dict;
use flate2::read::MultiGzDecoder;
use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufReader, Read};

const USAGE: &str = "usage: priva mondo-scope validate --registry TSV[.gz] --hpo-assertions TSV[.gz]\n\
Check the disease-scope registry and the scoped HPO assertions before the installer records them.";

const REGISTRY: [&str; 6] = [
    "disease_id",
    "mondo_id",
    "disease_scope",
    "priva_scope",
    "scope_evidence",
    "scope_review_status",
];
const ASSERTIONS: [&str; 11] = [
    "gene_symbol",
    "disease_id",
    "hpo_id",
    "frequency",
    "evidence",
    "reference",
    "mondo_id",
    "disease_scope",
    "priva_scope",
    "scope_evidence",
    "scope_review_status",
];

/// The number of data rows of a table that has every `required` column.
fn inspect(path: &str, required: &[&str]) -> std::result::Result<usize, String> {
    let file = File::open(path).map_err(|error| format!("{error}: {path}"))?;
    let source: Box<dyn Read> = if path.ends_with(".gz") {
        Box::new(MultiGzDecoder::new(BufReader::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let rows = tsv_dict::open(source)?;
    let missing: BTreeSet<&str> = required
        .iter()
        .copied()
        .filter(|name| !rows.header.iter().any(|field| field == name))
        .collect();
    if !missing.is_empty() {
        let names: Vec<&str> = missing.into_iter().collect();
        return Err(format!(
            "{path} missing required columns: {}",
            names.join(", ")
        ));
    }
    let mut count = 0;
    for row in rows {
        row.map_err(|error| format!("{path}: {error}"))?;
        count += 1;
    }
    if count == 0 {
        return Err(format!("{path} contains no data rows"));
    }
    Ok(count)
}

pub(crate) fn validate(registry: &str, assertions: &str) -> std::result::Result<String, String> {
    let registry_count = inspect(registry, &REGISTRY)?;
    let assertion_count = inspect(assertions, &ASSERTIONS)?;
    Ok(format!(
        "disease_scope_rows={registry_count}; scoped_hpo_assertions={assertion_count}"
    ))
}

/// `ARGS` after `validate`.
pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let (mut registry, mut assertions) = (None, None);
    let mut rest = args.iter();
    while let Some(arg) = rest.next() {
        let slot = match arg.as_str() {
            "--registry" => &mut registry,
            "--hpo-assertions" => &mut assertions,
            "--help" | "-h" => {
                println!("{USAGE}");
                return Ok(());
            }
            _ => return Err(invalid(format!("unrecognized argument {arg}\n{USAGE}"))),
        };
        *slot = Some(
            rest.next()
                .ok_or_else(|| invalid(format!("{arg} needs a value")))?
                .clone(),
        );
    }
    let (Some(registry), Some(assertions)) = (registry, assertions) else {
        return Err(invalid(USAGE));
    };
    match validate(&registry, &assertions) {
        Ok(summary) => {
            println!("{summary}");
            Ok(())
        }
        Err(message) => {
            eprintln!("priva: {message}");
            std::process::exit(1)
        }
    }
}
