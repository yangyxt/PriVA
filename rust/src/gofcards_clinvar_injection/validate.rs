//! `scripts/validate_gofcards_cache.py`, the structure check behind the
//! installer's `validate_gofcards_exact_gof_cache`; the main preflight, which
//! runs before priva_binary is recorded, keeps the script.
//! `priva gofcards-cache-validate CACHE [KEY=VALUE ...]` checks the cache and
//! any recorded provenance; a value left empty is skipped.
use super::{entries, or_empty};
use crate::builder_error;
use crate::cli::{Result, invalid};
use crate::gene_nonlof_mechanism::{equal, error, get, load_json, py_str, repr, truthy, type_name};
use crate::json_text::Json;
use crate::python_path::repr as text_repr;
use crate::python_text::{decimal_digit_value, python_strip};
use std::path::Path;

const HELP: &str = "usage: priva gofcards-cache-validate CACHE [KEY=VALUE ...]\n\
Check a GoFCards cache and the provenance it records; a KEY may be dotted, and an\n\
empty VALUE is skipped. Prints the cache's counts.";

const ALLOWED: &[&str] = &[
    "eligible",
    "quarantined_gene_discordance",
    "quarantined_allele_gene_discordance",
    "quarantined_missing_source_gene",
    "quarantined_reviewed_lof",
    "quarantined_reviewed_mixed",
    "quarantined_reviewed_dominant_negative",
    "quarantined_reviewed_uncertain",
    "quarantined_reviewed_exclude",
    "quarantined_mechanism_review_required",
];

fn stop(message: String) -> Box<dyn std::error::Error> {
    error("SystemExit", message)
}

/// `^loc_[0-9XYMT]+:\d+:[ACGTN-]+->[ACGTN-]+_grch37$`.
fn well_formed(variant_id: &str) -> bool {
    let body = variant_id.strip_suffix('\n').unwrap_or(variant_id);
    let Some(body) = body
        .strip_prefix("loc_")
        .and_then(|body| body.strip_suffix("_grch37"))
    else {
        return false;
    };
    let mut parts = body.splitn(3, ':');
    let (Some(chrom), Some(position), Some(alleles)) = (parts.next(), parts.next(), parts.next())
    else {
        return false;
    };
    let allele = |text: &str| {
        !text.is_empty()
            && text
                .bytes()
                .all(|byte| matches!(byte, b'A' | b'C' | b'G' | b'T' | b'N' | b'-'))
    };
    let chrom_ok = !chrom.is_empty()
        && chrom
            .bytes()
            .all(|byte| byte.is_ascii_digit() || matches!(byte, b'X' | b'Y' | b'M' | b'T'));
    let position_ok =
        !position.is_empty() && position.chars().all(|c| decimal_digit_value(c).is_some());
    let Some((reference, alternate)) = alleles.split_once("->") else {
        return false;
    };
    chrom_ok && position_ok && allele(reference) && allele(alternate)
}

/// `len(value)`.
fn length(value: &Json) -> Result<usize> {
    match value {
        Json::Array(items) => Ok(items.len()),
        Json::Object(items) => Ok(items.len()),
        Json::Str(body) => Ok(body.chars().count()),
        other => Err(error(
            "TypeError",
            format!("object of type '{}' has no len()", type_name(other)),
        )),
    }
}

/// `key in cache`.
fn contains(cache: &Json, key: &str) -> Result<bool> {
    match cache {
        Json::Object(_) => Ok(cache.get(key).is_some()),
        Json::Array(items) => Ok(items.iter().any(|item| equal(item, &Json::Str(key.into())))),
        Json::Str(body) => Ok(body.contains(key)),
        other => Err(error(
            "TypeError",
            format!("argument of type '{}' is not iterable", type_name(other)),
        )),
    }
}

/// `item in container` for the coding keys of one transcript view.
fn holds(container: &Json, item: &Json) -> Result<bool> {
    match (container, item) {
        (Json::Object(_), Json::Str(key)) => Ok(container.get(key).is_some()),
        (Json::Object(_), Json::Array(_) | Json::Object(_)) => Err(error(
            "TypeError",
            format!("unhashable type: '{}'", type_name(item)),
        )),
        (Json::Object(_), _) => Ok(false),
        (Json::Array(items), _) => Ok(items.iter().any(|entry| equal(entry, item))),
        (Json::Str(body), Json::Str(key)) => Ok(body.contains(key.as_str())),
        (Json::Str(_), other) => Err(error(
            "TypeError",
            format!(
                "'in <string>' requires string as left operand, not {}",
                type_name(other)
            ),
        )),
        (other, _) => Err(error(
            "TypeError",
            format!("argument of type '{}' is not iterable", type_name(other)),
        )),
    }
}

#[derive(Default)]
struct Totals {
    genes: usize,
    variants: usize,
    eligible: usize,
    quarantined: usize,
    with_protein: usize,
    with_hg38: usize,
    hgvsc_keys: usize,
    hgvsp_keys: usize,
    evidence: usize,
}

fn check_variant(
    path: &str,
    symbol: &str,
    variant_id: &str,
    variant: &Json,
    totals: &mut Totals,
) -> Result<()> {
    let empty = Json::Object(vec![]);
    totals.variants += 1;
    if !well_formed(variant_id) {
        return Err(stop(format!(
            "{path}: {symbol}/{variant_id} is not a well-formed variant identifier"
        )));
    }
    let record = or_empty(get(variant, "record")?, &empty);
    let eligibility = or_empty(get(record, "eligibility")?, &empty);
    let state = get(eligibility, "status")?
        .cloned()
        .unwrap_or_else(|| Json::Str(String::new()));
    if !matches!(&state, Json::Str(value) if ALLOWED.contains(&value.as_str())) {
        return Err(stop(format!(
            "{path}: {symbol}/{variant_id} has unknown eligibility {}",
            repr(&state)
        )));
    }
    let assemblies = entries(get(variant, "assemblies")?)?;
    if assemblies.is_empty() {
        return Err(stop(format!(
            "{path}: {symbol}/{variant_id} is present on no assembly"
        )));
    }
    let mut protein = false;
    for (assembly, block) in assemblies {
        if assembly != "hg19" && assembly != "hg38" {
            return Err(stop(format!(
                "{path}: {symbol}/{variant_id} has unknown assembly {}",
                text_repr(assembly)
            )));
        }
        let coords = or_empty(get(block, "genomic")?, &empty);
        let chrom = get(coords, "chrom")?.is_some_and(truthy);
        let position = matches!(get(coords, "pos")?, Some(Json::Int(_) | Json::Bool(_)));
        if !chrom || !position {
            return Err(stop(format!(
                "{path}: {symbol}/{variant_id}/{assembly} has no usable coordinates"
            )));
        }
        if !get(coords, "ref")?.is_some_and(truthy) || !get(coords, "alt")?.is_some_and(truthy) {
            return Err(stop(format!(
                "{path}: {symbol}/{variant_id}/{assembly} has no reference or alternate"
            )));
        }
        for (_, view) in entries(get(block, "transcripts")?)? {
            let by_c = or_empty(get(view, "by_hgvsc")?, &empty);
            let by_p = or_empty(get(view, "by_hgvsp")?, &empty);
            totals.hgvsc_keys += length(by_c)?;
            totals.hgvsp_keys += length(by_p)?;
            protein |= truthy(by_p);
            for (hgvsp, coding_keys) in entries(Some(by_p))? {
                let keys = match coding_keys {
                    Json::Array(keys) if !keys.is_empty() => keys,
                    _ => {
                        return Err(stop(format!(
                            "{path}: {symbol}/{variant_id} protein key {} maps to no coding change",
                            text_repr(hgvsp)
                        )));
                    }
                };
                for key in keys {
                    if !holds(by_c, key)? {
                        return Err(stop(format!(
                            "{path}: {symbol}/{variant_id} protein key {} points at unknown {}",
                            text_repr(hgvsp),
                            repr(key)
                        )));
                    }
                }
            }
        }
    }
    if assemblies.iter().any(|(assembly, _)| assembly == "hg38") {
        totals.with_hg38 += 1;
    }
    if protein {
        totals.with_protein += 1;
    }
    if matches!(&state, Json::Str(value) if value == "eligible") {
        totals.eligible += 1;
    } else {
        totals.quarantined += 1;
        if !get(eligibility, "reason")?.is_some_and(truthy) {
            return Err(stop(format!(
                "{path}: {symbol}/{variant_id} is quarantined without a reason"
            )));
        }
    }
    if let Some(evidence) = get(record, "evidence")?.filter(|value| truthy(value)) {
        totals.evidence += length(evidence)?;
    }
    Ok(())
}

/// `cache[key]` after `key in cache` held.
fn item<'a>(cache: &'a Json, key: &str) -> Result<&'a Json> {
    match cache {
        Json::Object(_) => cache
            .get(key)
            .ok_or_else(|| error("KeyError", text_repr(key))),
        Json::Array(_) => Err(error(
            "TypeError",
            "list indices must be integers or slices, not str",
        )),
        _ => Err(error(
            "TypeError",
            "string indices must be integers, not 'str'",
        )),
    }
}

fn validate(path: &str, expected: &[String]) -> Result<()> {
    let cache = load_json(Path::new(path))?;
    for key in ["metadata", "genes"] {
        if !contains(&cache, key)? {
            return Err(stop(format!(
                "{path} has no top-level {} block",
                text_repr(key)
            )));
        }
    }
    let meta = item(&cache, "metadata")?;
    let declares = |key: &str, want: &str| -> Result<bool> {
        Ok(get(meta, key)?.is_some_and(|value| equal(value, &Json::Str(want.into()))))
    };
    if !declares("source", "GoFCards")? {
        return Err(stop(format!(
            "{path} metadata does not declare GoFCards as the source"
        )));
    }
    if !declares("mechanism", "GOF")? {
        return Err(stop(format!(
            "{path} metadata does not declare GOF as the mechanism"
        )));
    }
    for pair in expected {
        let (key, want) = pair.split_once('=').unwrap_or((pair, ""));
        if want.is_empty() {
            continue;
        }
        let mut have = Some(meta);
        for part in key.split('.') {
            have = match have {
                Some(value @ Json::Object(_)) => value.get(part),
                _ => None,
            };
        }
        let have = have.unwrap_or(&Json::Null);
        if py_str(have) != want {
            return Err(stop(format!(
                "{path} was built with {key}={}, but {} is configured now",
                repr(have),
                text_repr(want)
            )));
        }
    }
    let mut totals = Totals::default();
    let genes = item(&cache, "genes")?;
    let Json::Object(genes) = genes else {
        return Err(error(
            "AttributeError",
            format!("'{}' object has no attribute 'items'", type_name(genes)),
        ));
    };
    for (symbol, gene) in genes {
        totals.genes += 1;
        if python_strip(symbol).is_empty() {
            return Err(stop(format!("{path} has an empty gene symbol key")));
        }
        for (variant_id, variant) in entries(get(gene, "variants")?)? {
            check_variant(path, symbol, variant_id, variant, &mut totals)?;
        }
    }
    if totals.variants == 0 {
        return Err(stop(format!("{path} contains no variants")));
    }
    if totals.eligible == 0 {
        return Err(stop(format!("{path} contains no runtime-eligible variants")));
    }
    if totals.with_protein == 0 {
        return Err(stop(format!(
            "{path} contains no variant with a protein-change key"
        )));
    }
    println!(
        "gofcards_genes={}; variants={}; eligible={}; quarantined={}; with_protein_key={}; \
         on_hg38={}; hgvsc_keys={}; hgvsp_keys={}; evidence_entries={}",
        totals.genes,
        totals.variants,
        totals.eligible,
        totals.quarantined,
        totals.with_protein,
        totals.with_hg38,
        totals.hgvsc_keys,
        totals.hgvsp_keys,
        totals.evidence
    );
    Ok(())
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    match args.first().map(String::as_str) {
        Some("-h" | "--help") => {
            println!("{HELP}");
            Ok(())
        }
        Some(path) => builder_error::finish(validate(path, &args[1..])),
        None => Err(invalid("the following arguments are required: CACHE")),
    }
}
