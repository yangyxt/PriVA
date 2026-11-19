//! `priva mechanism-evidence validate --evidence TSV`: the check of the
//! DDG2P/G2P and Orphadata mechanism evidence table that
//! `gene_pathogenic_mechanism_cache_install` ran as an inline Python block
//! (`validate_gene_pathogenic_mechanism_evidence` before this command).
//!
//! The table is read as `csv.DictReader` reads a tab-separated file (blank
//! lines skipped, a short row's missing fields absent). It must have the 18
//! columns PriVA reads, at least one row, only the sources `G2P_DDG2P` and
//! `Orphadata`, and at least one each of: a G2P row, a G2P loss-of-function
//! row (`LOF` in `normalized_mechanisms` or "loss of function" in
//! `patho_mode_raw`), a G2P row with a penetrance tied to a condition, an
//! Orphadata row and an Orphadata row mapped to MONDO. A penetrance other than
//! `complete` or `incomplete`, or one without a condition identifier, stops
//! the check at its row. Success prints the counts as the Python block did; a
//! failure prints the reason and exits with 1. Compared with the Python block
//! by `rust/tests/install_checks_compare.py`.
use crate::cli::{Result, invalid};
use crate::tsv_dict;
use std::collections::BTreeSet;

const USAGE: &str = "usage: priva mechanism-evidence validate --evidence TSV\n\
Check the gene condition-mechanism evidence table before the installer records it.";

const REQUIRED: [&str; 18] = [
    "gene_symbol",
    "source",
    "source_record_id",
    "source_condition_id",
    "mondo_id",
    "disease_scope",
    "priva_scope",
    "scope_review_status",
    "disease_label",
    "inheritance",
    "penetrance_raw",
    "penetrance_hpo_ids",
    "normalized_penetrance",
    "patho_mode_raw",
    "normalized_mechanisms",
    "mechanism_confidence",
    "disease_confidence",
    "pmids",
];

/// `str.isspace()`: Unicode white space and U+001C to U+001F.
fn is_space(c: char) -> bool {
    c.is_whitespace() || ('\u{1c}'..='\u{1f}').contains(&c)
}

/// Python's `repr()` of a `str`.
fn repr(text: &str) -> String {
    let quote = if text.contains('\'') && !text.contains('"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::new();
    out.push(quote);
    for c in text.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c == quote => {
                out.push('\\');
                out.push(c);
            }
            c if (c as u32) < 0x20 || c as u32 == 0x7f => {
                out.push_str(&format!("\\x{:02x}", c as u32))
            }
            c => out.push(c),
        }
    }
    out.push(quote);
    out
}

pub(crate) fn validate(path: &str) -> std::result::Result<String, String> {
    let file = std::fs::File::open(path).map_err(|error| format!("{error}: {path}"))?;
    let rows = tsv_dict::open(std::io::BufReader::new(file))?;
    let missing: BTreeSet<&str> = REQUIRED
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
    let (source_at, mechanisms_at, raw_at, penetrance_at, condition_at, mondo_at) = (
        rows.column("source"),
        rows.column("normalized_mechanisms"),
        rows.column("patho_mode_raw"),
        rows.column("normalized_penetrance"),
        rows.column("source_condition_id"),
        rows.column("mondo_id"),
    );
    let (mut row_count, mut g2p_rows, mut g2p_lof_rows, mut g2p_penetrance_rows) = (0, 0, 0, 0);
    let (mut orphadata_rows, mut orphadata_mondo_rows) = (0, 0);
    let mut unexpected = BTreeSet::new();
    for (number, row) in rows.enumerate() {
        let row = row.map_err(|error| format!("{path}: {error}"))?;
        row_count += 1;
        let field = |at: Option<usize>| at.and_then(|at| row.get(at)).map(String::as_str);
        let source = field(source_at).unwrap_or("");
        if source != "G2P_DDG2P" && source != "Orphadata" {
            unexpected.insert(if source.is_empty() {
                "<blank>".to_owned()
            } else {
                source.to_owned()
            });
        }
        if source == "G2P_DDG2P" {
            g2p_rows += 1;
            let short = || format!("{path} row {} has fewer fields than its header", number + 2);
            let normalized = field(mechanisms_at).ok_or_else(short)?.to_uppercase();
            let raw = field(raw_at).ok_or_else(short)?.to_lowercase();
            if normalized.contains("LOF") || raw.contains("loss of function") {
                g2p_lof_rows += 1;
            }
            let penetrance = field(penetrance_at)
                .ok_or_else(short)?
                .trim_matches(is_space)
                .to_lowercase();
            if !penetrance.is_empty() {
                if penetrance != "complete" && penetrance != "incomplete" {
                    return Err(format!(
                        "{path} has unsupported penetrance {}",
                        repr(&penetrance)
                    ));
                }
                let named = |at| field(at).is_some_and(|value| !value.is_empty());
                if !(named(condition_at) || named(mondo_at)) {
                    return Err(format!(
                        "{path} has penetrance without a condition identifier"
                    ));
                }
                g2p_penetrance_rows += 1;
            }
        }
        if source == "Orphadata" {
            orphadata_rows += 1;
            if field(mondo_at).is_some_and(|value| !value.is_empty()) {
                orphadata_mondo_rows += 1;
            }
        }
    }
    if row_count == 0 {
        return Err(format!("{path} contains no evidence rows"));
    }
    if !unexpected.is_empty() {
        let names: Vec<String> = unexpected.into_iter().collect();
        return Err(format!(
            "{path} contains non-condition sources: {}",
            names.join(", ")
        ));
    }
    for (count, what) in [
        (g2p_rows, "G2P_DDG2P rows"),
        (g2p_lof_rows, "DDG2P/G2P LoF rows"),
        (g2p_penetrance_rows, "condition-linked G2P penetrance rows"),
        (orphadata_rows, "explicit Orphadata mechanism rows"),
        (orphadata_mondo_rows, "MONDO-mapped Orphadata rows"),
    ] {
        if count == 0 {
            return Err(format!("{path} contains no {what}"));
        }
    }
    Ok(format!(
        "evidence_rows={row_count}; g2p_rows={g2p_rows}; g2p_lof_rows={g2p_lof_rows}; \
g2p_penetrance_rows={g2p_penetrance_rows}; orphadata_rows={orphadata_rows}; \
orphadata_mondo_rows={orphadata_mondo_rows}"
    ))
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args
        .first()
        .is_some_and(|arg| arg == "--help" || arg == "-h")
    {
        println!("{USAGE}");
        return Ok(());
    }
    let evidence = match args {
        [mode, option, path] if mode == "validate" && option == "--evidence" => path,
        _ => return Err(invalid(USAGE)),
    };
    match validate(evidence) {
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

#[cfg(test)]
mod tests {
    use super::*;

    fn table(rows: &[&str]) -> String {
        let mut text = REQUIRED.join("\t");
        text.push('\n');
        for row in rows {
            text.push_str(row);
            text.push('\n');
        }
        let folder = std::env::temp_dir().join(format!("priva-evidence-{}", std::process::id()));
        std::fs::create_dir_all(&folder).unwrap();
        let path = folder.join(format!("{}.tsv", rows.len()));
        std::fs::write(&path, text).unwrap();
        path.to_string_lossy().into_owned()
    }

    fn row(
        source: &str,
        condition: &str,
        mondo: &str,
        penetrance: &str,
        raw: &str,
        mechanisms: &str,
    ) -> String {
        let mut fields = vec![""; 18];
        fields[0] = "GENE";
        fields[1] = source;
        fields[3] = condition;
        fields[4] = mondo;
        fields[12] = penetrance;
        fields[13] = raw;
        fields[14] = mechanisms;
        fields.join("\t")
    }

    #[test]
    fn complete_table_and_failures() {
        let g2p = row(
            "G2P_DDG2P",
            "OMIM:1",
            "",
            "Complete ",
            "loss of function",
            "",
        );
        let orpha = row("Orphadata", "", "MONDO:1", "", "", "GOF");
        let good = table(&[&g2p, "", &orpha]);
        assert_eq!(
            validate(&good).unwrap(),
            "evidence_rows=2; g2p_rows=1; g2p_lof_rows=1; g2p_penetrance_rows=1; orphadata_rows=1; orphadata_mondo_rows=1"
        );
        let odd = row("G2P_DDG2P", "OMIM:1", "", "partial", "", "LOF");
        assert!(
            validate(&table(&[&odd]))
                .unwrap_err()
                .ends_with("has unsupported penetrance 'partial'")
        );
        let other = row("ClinGen", "", "", "", "", "");
        assert!(
            validate(&table(&[&g2p, &orpha, &other]))
                .unwrap_err()
                .ends_with("non-condition sources: ClinGen")
        );
        assert!(
            validate(&table(&[&orpha]))
                .unwrap_err()
                .ends_with("contains no G2P_DDG2P rows")
        );
        assert!(
            validate(&table(&[]))
                .unwrap_err()
                .ends_with("contains no evidence rows")
        );
        assert_eq!(repr("it's"), "\"it's\"");
    }
}
