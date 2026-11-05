//! Optional prepared G2P/Orphadata condition-evidence procedures.
//! Source order, filtering, default-source catches and publication match the
//! Python builder. These commands do not acquire sources or activate callers.
pub(crate) mod table;
mod word;
pub(crate) use word::is_word;
use crate::builder_table::{strings, text};
use crate::cli::{Result, invalid};
use crate::hpo_penetrance::{self, Values};
use crate::interpro_mapping_builder::Element;
use crate::pandas_table::{Column, Table};
use crate::python_args::parse_long_options_converted;
use crate::python_path::{os_error, path_value};
use crate::python_text::{is_python_space, python_strip};
use std::cell::RefCell;
use std::collections::{BTreeSet, HashMap};
use std::fs::{self, File};
use std::path::Path;
use std::rc::Rc;

const G2P_URL: &str = "https://www.ebi.ac.uk/gene2phenotype/api/panel/all/download";
const ORPHA_URL: &str = "https://www.orphadata.com/data/xml/en_product6.xml";
const HELP: &str = "usage: priva gene-mechanism parse --kind {g2p,orphadata} --input FILE --output TSV\n\
       priva gene-mechanism prepared --raw-dir DIR --output TSV [--disease-scope-registry TSV]\n\
       priva gene-mechanism attach --input UNIFIED_TSV --disease-scope-registry TSV --output TSV\n\
prepared also accepts --g2p-ddg2p-url and --orphadata-gene-disease-url to retain explicit source identity/failure policy.\n\
Inputs must already be prepared. Acquisition and installer routes remain separate.";

const UNIFIED: &[&str] = &[
    "gene_symbol",
    "mechanism",
    "assertion_level",
    "source",
    "source_record_id",
    "source_condition_id",
    "mondo_id",
    "disease_scope",
    "priva_scope",
    "scope_review_status",
    "disease",
    "inheritance",
    "penetrance_raw",
    "penetrance_hpo_ids",
    "normalized_penetrance",
    "phenotype_terms",
    "confidence",
    "disease_confidence",
    "pmids",
    "patho_mode_raw",
    "assembly",
    "chromosome",
    "position",
    "ref",
    "alt",
    "hgvsc",
    "hgvsp",
    "consequence",
    "notes",
];
const EVIDENCE: &[&str] = &[
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
    "evidence_url",
    "panel_or_submitter",
    "phenotype_terms",
    "notes",
];
fn error(kind: &'static str, message: impl Into<String>) -> Box<dyn std::error::Error> {
    crate::builder_error::category(kind, invalid(message))
}
fn exists(path: &Path) -> Result<bool> {
    match fs::metadata(path) {
        Ok(_) => Ok(true),
        Err(e)
            if matches!(
                e.raw_os_error(),
                Some(libc::ENOENT | libc::ENOTDIR | libc::EBADF | libc::ELOOP)
            ) =>
        {
            Ok(false)
        }
        Err(e) => Err(os_error(e, path, None)),
    }
}
pub(crate) fn norm_key(value: &str) -> String {
    let mut out = String::new();
    let mut separator = false;
    for c in value.to_lowercase().chars() {
        if c.is_ascii_lowercase() || c.is_ascii_digit() {
            if separator && !out.is_empty() {
                out.push('_');
            }
            out.push(c);
            separator = false;
        } else {
            separator = true;
        }
    }
    out
}
fn word(c: char) -> bool {
    word::is_word(c)
}
fn dominant_negative(text: &str) -> bool {
    text.match_indices("dominant").any(|(at, _)| {
        if text[..at].chars().next_back().is_some_and(word) {
            return false;
        }
        let rest = &text[at + 8..];
        let trimmed = rest.trim_start_matches(is_python_space);
        rest.len() != trimmed.len()
            && trimmed
                .strip_prefix("negative")
                .is_some_and(|tail| !tail.chars().next().is_some_and(word))
    })
}
pub(crate) fn normalize_mechanism(raw: &str) -> Vec<&'static str> {
    // Rust's Unicode version can assign new capital letters that CPython 3.12
    // still treats as non-word characters (for example U+A7CB). Lowering such
    // a letter to an older word character would change Python's regex boundary.
    // Only lower the pinned Letter/Number set for these ASCII term searches.
    let mut lowered = String::new();
    for character in raw.chars() {
        if word(character) {
            lowered.extend(character.to_lowercase());
        } else {
            lowered.push(character);
        }
    }
    let value = lowered.replace(['_', '-'], " ");
    let mut labels = BTreeSet::new();
    if dominant_negative(&value) {
        labels.insert("DOMINANT_NEGATIVE");
    }
    for (terms, label) in [
        (
            &["gain of function", "gain function", "activating"][..],
            "GOF",
        ),
        (
            &["activated", "activation", "increased activity"][..],
            "ACTIVATING",
        ),
        (
            &[
                "loss of function",
                "loss function",
                "lof",
                "absent gene product",
                "absence of gene product",
            ][..],
            "LOF",
        ),
        (
            &["reduced function", "decreased function", "partial loss"][..],
            "HYPMORPHIC_LOF",
        ),
        (
            &[
                "haploinsufficiency",
                "haplo insufficient",
                "decreased gene dosage",
                "decreased dosage",
            ][..],
            "HAPLOINSUFFICIENCY",
        ),
        (
            &["triplosensitivity", "triplo sensitivity"][..],
            "TRIPLOSENSITIVITY",
        ),
        (
            &["increased gene dosage", "increased dosage"][..],
            "INCREASED_DOSAGE",
        ),
        (
            &["altered gene product", "altered protein"][..],
            "ALTERED_PRODUCT",
        ),
        (
            &["all missense", "missense", "in frame", "inframe"][..],
            "MISSENSE_OR_INFRAME",
        ),
    ] {
        if terms.iter().any(|term| value.contains(term)) {
            labels.insert(label);
        }
    }
    if labels.is_empty() && !python_strip(raw).is_empty() {
        labels.insert("UNMAPPED_MECHANISM");
    }
    labels.into_iter().collect()
}
/// `mechanism_confidence_from_text`, shared with the non-LOF builder.
pub(crate) fn mechanism_confidence(value: &str) -> &'static str {
    confidence(value)
}
fn confidence(value: &str) -> &'static str {
    match python_strip(value) {
        "3" => return "high",
        "2" => return "moderate",
        "1" => return "limited",
        _ => {}
    }
    let value = value.to_lowercase();
    for (terms, label) in [
        (&["definitive", "sufficient", "green", "strong"][..], "high"),
        (&["moderate", "emerging"][..], "moderate"),
        (&["limited", "amber", "supportive"][..], "limited"),
        (
            &["refuted", "disputed", "red", "unlikely"][..],
            "conflicting_or_refuted",
        ),
    ] {
        if terms.iter().any(|term| value.contains(term)) {
            return label;
        }
    }
    ""
}
fn row(fields: &[(&str, String)]) -> Vec<String> {
    let mut result = vec![String::new(); UNIFIED.len()];
    result[2] = "gene_level".into();
    for (name, value) in fields {
        result[UNIFIED.iter().position(|n| n == name).unwrap()] = value.clone();
    }
    result[0] = python_strip(&result[0]).to_owned();
    result
}
fn parse_g2p(path: &Path) -> Result<Table> {
    if !exists(path)? {
        return Ok(strings(UNIFIED, vec![]));
    }
    let input = table::read(path, b',', true)?;
    let by_name = input
        .names
        .iter()
        .enumerate()
        .map(|(i, name)| (norm_key(name), i))
        .collect::<HashMap<_, _>>();
    let resolve = |names: &[&str]| {
        names
            .iter()
            .find_map(|name| by_name.get(&norm_key(name)).copied())
    };
    let columns = [
        resolve(&["gene symbol", "gene_symbol"]),
        resolve(&["disease name", "disease_label", "disease"]),
        resolve(&[
            "molecular mechanism",
            "mutation consequence",
            "mutation_consequence",
        ]),
        resolve(&["confidence", "confidence category", "g2p relation label"]),
        resolve(&["allelic requirement", "allelic_requirement"]),
        resolve(&["cross cutting modifier", "cross_cutting_modifier"]),
        resolve(&["phenotypes", "phenotype_terms"]),
        resolve(&["pmids", "publications"]),
        resolve(&["g2p id", "g2p_id"]),
        resolve(&["disease mim", "disease_mim"]),
        resolve(&["disease MONDO", "disease_mondo"]),
    ];
    let mut rows = Vec::new();
    for at in 0..input.rows {
        let get = |index: usize| {
            columns[index].map_or("", |col| python_strip(text(&input, &input.names[col], at)))
        };
        let gene = get(0);
        if gene.is_empty() {
            continue;
        }
        let ids = hpo_penetrance::recognized_penetrance_hpo_ids(Values::Text(get(6)));
        let refs = ids.iter().map(String::as_str).collect::<Vec<_>>();
        rows.push(row(&[
            ("gene_symbol", gene.into()),
            ("mechanism", normalize_mechanism(get(2)).join("|")),
            ("source", "G2P_DDG2P".into()),
            ("source_record_id", get(8).into()),
            (
                "source_condition_id",
                if get(9).is_empty() {
                    String::new()
                } else {
                    format!("OMIM:{}", get(9))
                },
            ),
            ("mondo_id", get(10).into()),
            ("disease", get(1).into()),
            ("inheritance", get(4).into()),
            ("penetrance_raw", get(5).into()),
            ("penetrance_hpo_ids", ids.join(";")),
            (
                "normalized_penetrance",
                hpo_penetrance::normalize_penetrance_evidence(
                    Values::List(&refs),
                    Values::Text(get(5)),
                )
                .into(),
            ),
            ("phenotype_terms", get(6).into()),
            ("confidence", confidence(get(3)).into()),
            ("disease_confidence", get(3).into()),
            ("pmids", get(7).into()),
            ("patho_mode_raw", get(2).into()),
        ]));
    }
    Ok(strings(UNIFIED, rows))
}

fn nested_text(element: &Element, parent: &str, child: &str) -> String {
    element
        .child(parent)
        .map_or_else(String::new, |node| node.child_text(child))
}
fn pmids(value: &str) -> String {
    let mut values = Vec::new();
    let mut start = None;
    for (at, c) in value
        .char_indices()
        .chain(std::iter::once((value.len(), '\0')))
    {
        if crate::python_text::decimal_digit_value(c).is_some() {
            start.get_or_insert(at);
        } else if let Some(begin) = start.take() {
            if value[at..]
                .trim_start_matches(is_python_space)
                .starts_with("[PMID]")
            {
                let id = &value[begin..at];
                if !values.contains(&id) {
                    values.push(id);
                }
            }
        }
    }
    values.join(";")
}
fn disorder(element: &Element, rows: &mut Vec<Vec<String>>) {
    let code = element.child_text("OrphaCode");
    let disease = element.child_text("Name");
    if code.is_empty() {
        return;
    }
    let Some(associations) = element.child("DisorderGeneAssociationList") else {
        return;
    };
    for association in associations
        .children
        .iter()
        .map(|c| c.borrow())
        .filter(|c| c.name == "DisorderGeneAssociation")
    {
        let status = nested_text(&association, "DisorderGeneAssociationStatus", "Name");
        let relationship = nested_text(&association, "DisorderGeneAssociationType", "Name");
        if status != "Assessed"
            || !relationship
                .to_lowercase()
                .contains("disease-causing germline")
        {
            continue;
        }
        let mechanisms = normalize_mechanism(&relationship)
            .into_iter()
            .filter(|m| ["GOF", "LOF"].contains(m))
            .collect::<Vec<_>>();
        if mechanisms.is_empty() {
            continue;
        }
        let gene = nested_text(&association, "Gene", "Symbol");
        if gene.is_empty() {
            continue;
        }
        let condition = format!("ORPHA:{code}");
        rows.push(row(&[
            ("gene_symbol", gene.clone()),
            ("mechanism", mechanisms.join("|")),
            ("source", "Orphadata".into()),
            (
                "source_record_id",
                format!("{condition}|{gene}|{}", norm_key(&relationship)),
            ),
            ("source_condition_id", condition),
            ("disease", disease.clone()),
            ("confidence", "high".into()),
            ("disease_confidence", status.clone()),
            (
                "pmids",
                pmids(&association.child_text("SourceOfValidation")),
            ),
            ("patho_mode_raw", relationship),
            ("notes", format!("association_status:{status}")),
        ]));
    }
}
fn parse_orphadata(path: &Path) -> Result<Table> {
    if !exists(path)? {
        return Ok(strings(UNIFIED, vec![]));
    }
    let file = File::open(path).map_err(|e| os_error(e, path, None))?;
    if file.metadata()?.is_dir() {
        return Err(os_error(
            std::io::Error::from_raw_os_error(libc::EISDIR),
            path,
            None,
        ));
    }
    let rows = Rc::new(RefCell::new(Vec::new()));
    let collected = rows.clone();
    crate::interpro_mapping_builder::visit_selected_elements(file, "Disorder", move |element| {
        disorder(element, &mut collected.borrow_mut())
    })?;
    let rows = std::mem::take(&mut *rows.borrow_mut());
    Ok(strings(UNIFIED, rows))
}

fn attach(mut output: Table, registry: Option<&Path>) -> Result<Table> {
    const FIELDS: &[&str] = &[
        "mondo_id",
        "disease_scope",
        "priva_scope",
        "scope_review_status",
    ];
    for name in FIELDS {
        if !output.has_column(name) {
            output.set_column(name, Column::Text(vec![Some(String::new()); output.rows]));
        }
    }
    let Some(path) = registry else {
        return Ok(output);
    };
    if output.rows == 0 || !exists(path)? {
        return Ok(output);
    }
    let registry = table::read(path, b'\t', false)?;
    let mut missing = std::iter::once("disease_id")
        .chain(FIELDS.iter().copied())
        .filter(|name| !registry.has_column(name))
        .collect::<Vec<_>>();
    missing.sort_unstable();
    if !missing.is_empty() {
        return Err(error(
            "ValueError",
            format!(
                "disease scope registry missing columns: {}",
                missing.join(", ")
            ),
        ));
    }
    if !output.has_column("source_condition_id") {
        return Err(error("KeyError", "'source_condition_id'"));
    }
    let retained = crate::pandas_polars::drop_duplicates_first(&registry, Some(&["disease_id"]))?;
    let registry = registry.take_rows(&retained);
    let left = strings(
        &["disease_id"],
        (0..output.rows)
            .map(|i| vec![text(&output, "source_condition_id", i).into()])
            .collect(),
    );
    let joins = crate::pandas_polars::left_join_text_rows(&left, &registry, &["disease_id"], true)?;
    let mut conflicts = 0;
    for name in FIELDS {
        let mut values = Vec::with_capacity(output.rows);
        for &(at, right) in &joins {
            let old = text(&output, name, at);
            let mapped = right.map_or("", |i| text(&registry, name, i));
            if *name == "mondo_id" && !old.is_empty() && !mapped.is_empty() && old != mapped {
                conflicts += 1;
            }
            values.push(Some(if old.is_empty() { mapped } else { old }.to_owned()));
        }
        output.set_column(name, Column::Text(values));
    }
    if conflicts > 0 {
        eprintln!(
            "condition_identity: retained source MONDO ID for {conflicts} registry conflicts"
        );
    }
    Ok(output)
}
fn evidence(unified: &Table, g2p_url: &str, orpha_url: &str) -> Table {
    let mut rows = Vec::new();
    for at in 0..unified.rows {
        let source = text(unified, "source", at);
        let url = match source {
            "G2P_DDG2P" => g2p_url,
            "Orphadata" => orpha_url,
            _ => continue,
        };
        rows.push(
            EVIDENCE
                .iter()
                .map(|name| match *name {
                    "disease_label" => text(unified, "disease", at).into(),
                    "normalized_mechanisms" => text(unified, "mechanism", at).replace('|', ";"),
                    "mechanism_confidence" => text(unified, "confidence", at).into(),
                    "evidence_url" => url.into(),
                    "panel_or_submitter" => String::new(),
                    name => text(unified, name, at).into(),
                })
                .collect(),
        );
    }
    strings(EVIDENCE, rows)
}
fn exception_message(cause: &(dyn std::error::Error + 'static)) -> String {
    // These parser procedures return typed Python categories. Python's local
    // catch formats str(exc), without the category or traceback envelope.
    let full = cause.to_string();
    if cause.is::<crate::builder_error::RuntimeError>() {
        full.split_once(": ")
            .map_or(full.clone(), |(_, message)| message.to_owned())
    } else {
        full
    }
}
fn prepared(raw: &Path, registry: Option<&Path>, g2p_url: &str, orpha_url: &str) -> Result<Table> {
    let mut rows = Vec::new();
    for (name, url, default, filename, parse) in [
        (
            "g2p_ddg2p",
            g2p_url,
            G2P_URL,
            "AllG2P.csv",
            parse_g2p as fn(&Path) -> Result<Table>,
        ),
        (
            "orphadata_gene_disease",
            orpha_url,
            ORPHA_URL,
            "orphadata/en_product6.xml",
            parse_orphadata as fn(&Path) -> Result<Table>,
        ),
    ] {
        match parse(&raw.join(filename)) {
            Ok(frame) => {
                for at in 0..frame.rows {
                    rows.push(
                        UNIFIED
                            .iter()
                            .map(|name| text(&frame, name, at).into())
                            .collect(),
                    );
                }
            }
            Err(cause) => {
                let message = exception_message(cause.as_ref());
                if url != default {
                    return Err(error(
                        "ValueError",
                        format!("Incompatible selected source {name} ({url}): {message}"),
                    ));
                }
                eprintln!("parse_error source={name}: {message}");
            }
        }
    }
    attach(strings(UNIFIED, rows), registry)
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(command) = args.first().map(String::as_str) else {
        return Err(invalid(HELP));
    };
    if ["-h", "--help"].contains(&command) {
        println!("{HELP}");
        return Ok(());
    }
    let (names, required): (&[&str], &[&str]) = match command {
        "parse" => (&["kind", "input", "output"], &["kind", "input", "output"]),
        "prepared" => (
            &[
                "raw-dir",
                "output",
                "disease-scope-registry",
                "g2p-ddg2p-url",
                "orphadata-gene-disease-url",
            ],
            &["raw-dir", "output"],
        ),
        "attach" => (
            &["input", "output", "disease-scope-registry"],
            &["input", "output", "disease-scope-registry"],
        ),
        _ => return Err(invalid(HELP)),
    };
    let Some(occurrences) = parse_long_options_converted(&args[1..], names, &[], |name, value| {
        if name.ends_with("-url") {
            crate::mechanism_source::check_source_url(value).map_err(|_| {
                format!(
                    "argument --{name}: invalid source_url value: {}",
                    crate::python_path::repr(value)
                )
            })?;
        }
        if name == "kind" && !["g2p", "orphadata"].contains(&value) {
            return Err("invalid --kind".into());
        }
        Ok(())
    })
    .map_err(invalid)?
    else {
        println!("{HELP}");
        return Ok(());
    };
    let values: HashMap<_, _> = occurrences.into_iter().collect();
    for name in required {
        if !values.contains_key(*name) {
            return Err(invalid(format!("required option --{name}")));
        }
    }
    if command == "parse" && !["g2p", "orphadata"].contains(&values["kind"].as_str()) {
        return Err(invalid("invalid --kind"));
    }
    crate::builder_error::finish((|| {
        let output = path_value(&values["output"]);
        let registry = values.get("disease-scope-registry").map(|v| path_value(v));
        let data = match command {
            "parse" => {
                let input = path_value(&values["input"]);
                if values["kind"] == "g2p" {
                    parse_g2p(&input)?
                } else {
                    parse_orphadata(&input)?
                }
            }
            "attach" => attach(
                table::read(&path_value(&values["input"]), b'\t', false)?,
                registry.as_deref(),
            )?,
            "prepared" => {
                let g2p = values.get("g2p-ddg2p-url").map_or(G2P_URL, String::as_str);
                let orpha = values
                    .get("orphadata-gene-disease-url")
                    .map_or(ORPHA_URL, String::as_str);
                let data = prepared(
                    &path_value(&values["raw-dir"]),
                    registry.as_deref(),
                    g2p,
                    orpha,
                )?;
                evidence(&data, g2p, orpha)
            }
            _ => unreachable!(),
        };
        table::publish(&output, &data)?;
        Ok(())
    })())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mechanism_vocabulary_retains_source_order_and_spelling() {
        assert_eq!(
            normalize_mechanism("dominant\u{1c}negative; partial loss; activating; missense"),
            [
                "DOMINANT_NEGATIVE",
                "GOF",
                "HYPMORPHIC_LOF",
                "MISSENSE_OR_INFRAME"
            ]
        );
        assert_eq!(
            normalize_mechanism("dominant negativex"),
            ["UNMAPPED_MECHANISM"]
        );
        assert_eq!(normalize_mechanism(""), Vec::<&str>::new());
        assert_eq!(confidence("not sufficient, disputed"), "high");
    }
    #[test]
    fn orphadata_pmids_preserve_unicode_digits_first_order() {
        assert_eq!(
            pmids("123 [PMID] ١٢٣\u{2003}[PMID] 123[PMID] 456 [OTHER]"),
            "123;١٢٣"
        );
    }
}
