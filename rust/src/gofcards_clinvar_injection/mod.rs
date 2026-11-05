//! Native `scripts/inject_clinvar_into_gofcards.py`: nests ClinVar VCV
//! condition evidence under each matched GoFCards variant, with the script's
//! options, and `scripts/validate_gofcards_cache.py`.
//!
//! The streamed VCV reader and its matching rules from
//! `scripts/clinvar_vcv.py`, the injection, the atomic gzip output, the audit
//! table, the statistics, the `[clinvar-inject]` log lines and exit statuses
//! follow the script. An exception the script let escape exits 1 here as
//! `priva: CATEGORY: message`.
//!
//! Differences:
//! - `--source-root CHECKOUT` is required: the cache records the SHA-256 of
//!   the checkout's Python injector as `clinvar.injector_sha256`, the value
//!   the installer compares, so a cache from either one is reused by the other.
//! - The `confirmed_by_*` counts of a variant linked by both routes are added
//!   in sorted order; the script's order follows its string hash.
//! - An integer field too large for 128 bits stays text.
//! - A VCV root element itself named `VariationArchive` keeps its attributes
//!   when a nested archive ends; the script's `root.clear()` would remove them.
mod validate;
mod vcv;

use crate::builder_error;
use crate::cli::{Result, invalid};
use crate::gene_nonlof_mechanism::{
    eligible, equal, error, get, load_gofcards_cache, object, py_str, sha256_file, text, truthy,
    type_name,
};
use crate::json_text::{self, DumpsOptions, Json};
use crate::python_args::{parse_long_options_converted, require_options};
use crate::python_path::{os_error, path_value, repr};
use crate::python_text::{python_int, python_strip};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::{self, File};
use std::io::Write;
use std::os::unix::fs::MetadataExt;
use std::path::{Path, PathBuf};
use vcv::Match;

pub(crate) use validate::dispatch as dispatch_validate;

const HELP: &str = "usage: priva gofcards-clinvar-injection --source-root CHECKOUT\n\
       --gofcards-cache GOFCARDS_CACHE --clinvar-vcv CLINVAR_VCV\n\
       [--clinvar-vcv-sha256 CLINVAR_VCV_SHA256] --out-json OUT_JSON\n\
       [--min-review-stars {2,3,4}] [--audit-tsv AUDIT_TSV] [--stats-json STATS_JSON]\n\
Attach ClinVar VCV condition evidence to each GoFCards variant; the options, outputs\n\
and exit statuses of scripts/inject_clinvar_into_gofcards.py.";

const VALUES: &[&str] = &[
    "source-root",
    "gofcards-cache",
    "clinvar-vcv",
    "clinvar-vcv-sha256",
    "out-json",
    "min-review-stars",
    "audit-tsv",
    "stats-json",
];

const AUDIT_COLUMNS: &str = "gene_symbol\tvcv_accession\tvariation_id\tvariation_name\t\
classification_scope\tmatched_component_context\trcv_accession\tcondition\t\
clinical_significance\treview_status\treview_stars\tmatched_scv_count\n";

/// The symbol and identifier of an indexed variant.
type VariantRow = (String, String);
/// `(build, normalized VCF allele)` to the variants on it.
type Lookup = HashMap<(String, String), Vec<VariantRow>>;

fn log(message: &str) {
    eprintln!("[clinvar-inject] {message}");
}

struct Options {
    source_root: PathBuf,
    gofcards_cache: PathBuf,
    clinvar_vcv: PathBuf,
    clinvar_vcv_sha256: Option<String>,
    out_json: PathBuf,
    min_review_stars: i128,
    audit_tsv: Option<PathBuf>,
    stats_json: Option<PathBuf>,
}

fn parse_options(args: &[String]) -> Result<Option<Options>> {
    let converted = parse_long_options_converted(args, VALUES, &[], |name, value| {
        if name == "min-review-stars" {
            let number = python_int(value)
                .ok_or_else(|| format!("argument --{name}: invalid int value: {}", repr(value)))?;
            if !(2..=4).contains(&number) {
                return Err(format!(
                    "argument --{name}: invalid choice: {number} (choose from 2, 3, 4)"
                ));
            }
        }
        Ok(())
    })
    .map_err(invalid)?;
    let Some(occurrences) = converted else {
        println!("{HELP}");
        return Ok(None);
    };
    let values: HashMap<String, String> = occurrences.into_iter().collect();
    require_options(
        &values,
        &["source-root", "gofcards-cache", "clinvar-vcv", "out-json"],
    )
    .map_err(invalid)?;
    let path = |name: &str| path_value(&values[name]);
    let optional = |name: &str| values.get(name).map(|value| path_value(value));
    Ok(Some(Options {
        source_root: path("source-root"),
        gofcards_cache: path("gofcards-cache"),
        clinvar_vcv: path("clinvar-vcv"),
        clinvar_vcv_sha256: values.get("clinvar-vcv-sha256").cloned(),
        out_json: path("out-json"),
        min_review_stars: values
            .get("min-review-stars")
            .and_then(|value| python_int(value))
            .unwrap_or(2),
        audit_tsv: optional("audit-tsv"),
        stats_json: optional("stats-json"),
    }))
}

/// `value or {}` for a looked-up value.
fn or_empty<'a>(value: Option<&'a Json>, empty: &'a Json) -> &'a Json {
    value.filter(|value| truthy(value)).unwrap_or(empty)
}

/// `(value or {}).items()`.
fn entries<'a>(value: Option<&'a Json>) -> Result<&'a [(String, Json)]> {
    match value.filter(|value| truthy(value)) {
        None => Ok(&[]),
        Some(Json::Object(items)) => Ok(items),
        Some(other) => Err(error(
            "AttributeError",
            format!("'{}' object has no attribute 'items'", type_name(other)),
        )),
    }
}

/// `iter_gofcards_variants`.
fn walk_variants(
    cache: &Json,
    mut visit: impl FnMut(&str, &str, &Json) -> Result<()>,
) -> Result<()> {
    for (symbol, gene) in entries(get(cache, "genes")?)? {
        for (variant_id, variant) in entries(get(gene, "variants")?)? {
            visit(symbol, variant_id, variant)?;
        }
    }
    Ok(())
}

/// `gofcards_variation_id`.
fn gofcards_variation_id(variant: &Json) -> Result<String> {
    let empty = Json::Object(vec![]);
    let record = or_empty(get(variant, "record")?, &empty);
    let annotations = or_empty(get(record, "annotations")?, &empty);
    Ok(match get(annotations, "clinvar_variation_id")? {
        Some(value) if truthy(value) => python_strip(&py_str(value)).to_owned(),
        _ => String::new(),
    })
}

/// `gofcards_genomic_index_key`.
fn genomic_keys(variant: &Json) -> Result<Vec<(String, String)>> {
    let empty = Json::Object(vec![]);
    let mut keys = Vec::new();
    for (assembly, block) in entries(get(variant, "assemblies")?)? {
        let coords = or_empty(get(block, "genomic")?, &empty);
        let chrom = get(coords, "chrom")?.filter(|value| truthy(value));
        let position = get(coords, "pos")?.filter(|value| !matches!(value, Json::Null));
        let (Some(chrom), Some(position)) = (chrom, position) else {
            continue;
        };
        let allele = |key: &str| get(coords, key).map(|value| value.map(py_str).unwrap_or_default());
        let joined = format!(
            "{}|{}|{}|{}",
            py_str(chrom),
            py_str(position),
            allele("ref")?,
            allele("alt")?
        );
        let parts: Vec<&str> = joined.split('|').collect();
        let [chrom, position, reference, alternate] = parts[..] else {
            return Err(error(
                "TypeError",
                format!(
                    "normalize_vcf_key() takes 4 positional arguments but {} were given",
                    parts.len()
                ),
            ));
        };
        keys.push((
            assembly.clone(),
            vcv::normalize_vcf_key(chrom, position, reference, alternate),
        ));
    }
    Ok(keys)
}

/// `index_gofcards_variants` over the genomic keys.
fn genomic_lookup(cache: &Json) -> Result<Lookup> {
    let mut index = Lookup::new();
    walk_variants(cache, |symbol, variant_id, variant| {
        if eligible(variant)? {
            for key in genomic_keys(variant)? {
                index
                    .entry(key)
                    .or_default()
                    .push((symbol.to_owned(), variant_id.to_owned()));
            }
        }
        Ok(())
    })?;
    Ok(index)
}

/// `index_gofcards_variants` over GoFCards' own ClinVar VariationIDs.
fn variation_id_lookup(cache: &Json) -> Result<HashMap<String, Vec<VariantRow>>> {
    let mut index: HashMap<String, Vec<VariantRow>> = HashMap::new();
    walk_variants(cache, |symbol, variant_id, variant| {
        if eligible(variant)? {
            let identifier = gofcards_variation_id(variant)?;
            if !identifier.is_empty() {
                index
                    .entry(identifier)
                    .or_default()
                    .push((symbol.to_owned(), variant_id.to_owned()));
            }
        }
        Ok(())
    })?;
    Ok(index)
}

/// A `collections.Counter` that keeps first-insertion order.
#[derive(Default)]
struct Counts(Vec<(String, i128)>);

impl Counts {
    fn add(&mut self, key: &str, count: i128) {
        match self.0.iter_mut().find(|(name, _)| name == key) {
            Some((_, value)) => *value += count,
            None => self.0.push((key.to_owned(), count)),
        }
    }
}

/// `clinvar_block`.
fn clinvar_block(found: &Match) -> Json {
    let mut hgvs: Vec<Json> = Vec::new();
    for expression in &found.hgvs {
        let expression = text(expression.clone());
        if !hgvs.iter().any(|item| equal(item, &expression)) {
            hgvs.push(expression);
        }
    }
    object([
        ("vcv_accession", text(python_strip(&found.vcv_accession))),
        ("variation_id", text(python_strip(&found.variation_id))),
        ("variation_name", text(python_strip(&found.variation_name))),
        ("hgvs", Json::Array(hgvs)),
        ("condition_assertions", Json::Array(found.conditions.clone())),
    ])
}

fn value_mut<'a>(value: &'a mut Json, key: &str) -> Option<&'a mut Json> {
    match value {
        Json::Object(items) => items
            .iter_mut()
            .find(|(name, _)| name == key)
            .map(|(_, item)| item),
        _ => None,
    }
}

fn variant_mut<'a>(cache: &'a mut Json, symbol: &str, variant_id: &str) -> Option<&'a mut Json> {
    let gene = value_mut(value_mut(cache, "genes")?, symbol)?;
    value_mut(value_mut(gene, "variants")?, variant_id)
}

/// The `clinvar_additional` branch of `inject`.
fn add_additional(variant: &mut Json, block: &Json, accession: &Json, counts: &mut Counts) -> Result<()> {
    if variant.get("clinvar_additional").is_none() {
        variant.set("clinvar_additional", Json::Array(vec![]));
    }
    match value_mut(variant, "clinvar_additional") {
        Some(Json::Array(extra)) => {
            for entry in extra.iter() {
                if get(entry, "vcv_accession")?.is_some_and(|seen| equal(seen, accession)) {
                    return Ok(());
                }
            }
            extra.push(block.clone());
            counts.add("additional_clinvar_records", 1);
            Ok(())
        }
        other => Err(error(
            "TypeError",
            format!(
                "clinvar_additional is a {}, not a list",
                other.map_or("NoneType", |value| type_name(value))
            ),
        )),
    }
}

/// `inject`.
fn inject(cache: &mut Json, matches: &[Match]) -> Result<Counts> {
    // The script's dict by variant identifier: the last gene holding it wins.
    let mut owners: HashMap<String, String> = HashMap::new();
    walk_variants(cache, |symbol, variant_id, _| {
        owners.insert(variant_id.to_owned(), symbol.to_owned());
        Ok(())
    })?;
    let mut counts = Counts::default();
    for found in matches {
        let block = clinvar_block(found);
        let accession = text(python_strip(&found.vcv_accession));
        let variation_id = python_strip(&found.variation_id);
        for (variant_id, route) in &found.records {
            let variant_id = python_strip(variant_id);
            let variant = owners
                .get(variant_id)
                .and_then(|symbol| variant_mut(cache, symbol, variant_id));
            let Some(variant) = variant else {
                counts.add("clinvar_matches_without_a_variant", 1);
                continue;
            };
            let existing = match get(variant, "clinvar")? {
                None | Some(Json::Null) => None,
                Some(existing) => Some(get(existing, "vcv_accession")?.cloned().unwrap_or(Json::Null)),
            };
            match existing {
                None => {
                    let mut routes = BTreeSet::from([*route]);
                    if !variation_id.is_empty() && gofcards_variation_id(variant)? == variation_id {
                        routes.insert("gofcards_variation_id");
                    }
                    let mut annotated = block.clone();
                    annotated.set(
                        "link_routes",
                        Json::Array(routes.iter().map(|route| text(*route)).collect()),
                    );
                    variant.set("clinvar", annotated);
                    counts.add("variants_annotated", 1);
                    for route in &routes {
                        counts.add(&format!("confirmed_by_{route}"), 1);
                    }
                    if routes.len() > 1 {
                        counts.add("confirmed_by_both_routes", 1);
                    }
                    counts.add("condition_assertions", found.conditions.len() as i128);
                }
                Some(seen) if !equal(&seen, &accession) => {
                    add_additional(variant, &block, &accession, &mut counts)?;
                }
                Some(_) => {}
            }
        }
    }
    Ok(counts)
}

type Identity = (u64, u64, u64, i64, i64, i64, i64);

fn source_identity(path: &Path) -> Result<Identity> {
    let meta = fs::metadata(path).map_err(|e| os_error(e, path, None))?;
    Ok((
        meta.dev(),
        meta.ino(),
        meta.size(),
        meta.mtime(),
        meta.mtime_nsec(),
        meta.ctime(),
        meta.ctime_nsec(),
    ))
}

/// `open_text(tmp, "wt")` and `atomic_write_text`.
fn write_atomic(path: &Path, body: &str) -> Result<()> {
    let name = path
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_default();
    let compressed = path.to_string_lossy().ends_with(".gz");
    let temporary = path.with_file_name(match name.strip_suffix(".gz") {
        Some(stem) => format!("{stem}.tmp.gz"),
        None => format!("{name}.tmp"),
    });
    if let Some(parent) = path.parent().filter(|parent| !parent.as_os_str().is_empty()) {
        fs::create_dir_all(parent).map_err(|e| os_error(e, parent, None))?;
    }
    let bytes = if compressed {
        crate::python_gzip::compress(body.as_bytes(), 9, true)
            .map_err(|e| error("OSError", e.to_string()))?
    } else {
        body.as_bytes().to_vec()
    };
    let written = fs::write(&temporary, bytes)
        .map_err(|e| os_error(e, &temporary, None))
        .and_then(|()| {
            fs::rename(&temporary, path).map_err(|e| os_error(e, &temporary, Some(path)))
        });
    if written.is_err() {
        let _ = fs::remove_file(&temporary);
    }
    written
}

fn field(value: &Json, key: &str) -> String {
    value.get(key).map(py_str).unwrap_or_default()
}

/// The rows of `flatten_match_audit`, in the script's column order.
fn audit_rows(matches: &[Match]) -> Vec<String> {
    let mut rows = Vec::new();
    for found in matches {
        for symbol in &found.gene_symbols {
            for condition in &found.conditions {
                let empty = Json::Object(vec![]);
                let classification = condition.get("germline_classification").unwrap_or(&empty);
                let names: Vec<String> = match condition.get("conditions") {
                    Some(Json::Array(items)) => items
                        .iter()
                        .map(|item| field(item, "name"))
                        .filter(|name| !name.is_empty())
                        .collect(),
                    _ => vec![],
                };
                let scvs = match condition.get("matched_scvs") {
                    Some(Json::Array(items)) => items.len(),
                    _ => 0,
                };
                rows.push(
                    [
                        symbol.as_str(),
                        &found.vcv_accession,
                        &found.variation_id,
                        &found.variation_name,
                        found.classification_scope,
                        found.component_context,
                        &field(condition, "rcv_accession"),
                        &names.join("|"),
                        &field(classification, "clinical_significance"),
                        &field(classification, "review_status"),
                        &field(classification, "review_stars"),
                        &scvs.to_string(),
                    ]
                    .join("\t"),
                );
            }
        }
    }
    rows
}

fn write_audit(path: &Path, matches: &[Match]) -> Result<()> {
    if let Some(parent) = path.parent().filter(|parent| !parent.as_os_str().is_empty()) {
        fs::create_dir_all(parent).map_err(|e| os_error(e, parent, None))?;
    }
    let rows = audit_rows(matches);
    let mut file = File::create(path).map_err(|e| os_error(e, path, None))?;
    let mut body = AUDIT_COLUMNS.to_owned();
    for row in &rows {
        body.push_str(row);
        body.push('\n');
    }
    file.write_all(body.as_bytes())
        .map_err(|e| os_error(e, path, None))?;
    log(&format!("wrote {} audit rows to {}", rows.len(), path.display()));
    Ok(())
}

fn run(options: Options) -> Result<()> {
    let source = &options.clinvar_vcv;
    let expected = options.clinvar_vcv_sha256.as_deref();
    let mut source_before = None;
    if let Some(expected) = expected {
        let before = source_identity(source)?;
        if sha256_file(source)? != expected || source_identity(source)? != before {
            return Err(error(
                "SystemExit",
                format!(
                    "Explicit ClinVar VCV XML checksum/identity mismatch before injection: {}",
                    source.display()
                ),
            ));
        }
        source_before = Some(before);
    }

    let mut cache = load_gofcards_cache(&options.gofcards_cache)?;
    let lookup = genomic_lookup(&cache)?;
    let by_variation_id = variation_id_lookup(&cache)?;
    log(&format!(
        "indexed {} GoFCards ClinVar variation ids and {} genomic keys from eligible variants",
        by_variation_id.len(),
        lookup.len()
    ));

    let (matches, parse_stats) =
        vcv::stream_parse_clinvar_vcv(source, lookup, by_variation_id, options.min_review_stars)?;
    if let (Some(expected), Some(before)) = (expected, source_before) {
        if source_identity(source)? != before
            || sha256_file(source)? != expected
            || source_identity(source)? != before
        {
            return Err(error(
                "SystemExit",
                format!(
                    "Explicit ClinVar VCV XML changed during injection; output not published: {}",
                    source.display()
                ),
            ));
        }
    }
    log(&format!("ClinVar returned {} matched records", matches.len()));

    let counts = inject(&mut cache, &matches)?;
    let injector = options.source_root.join("scripts/inject_clinvar_into_gofcards.py");
    let mut provenance = object([
        ("source", text(source.display().to_string())),
        ("min_review_stars", Json::Int(options.min_review_stars)),
        ("injector_sha256", text(sha256_file(&injector)?)),
        (
            "source_cache_sha256",
            text(sha256_file(&options.gofcards_cache)?),
        ),
    ]);
    let stats = stats_object(&parse_stats, &counts);
    if let Json::Object(items) = &stats {
        for (key, value) in items {
            provenance.set(key, value.clone());
        }
    }
    if let Some(expected) = expected {
        provenance.set("source_sha256", text(expected));
    }
    if let Some(metadata) = value_mut(&mut cache, "metadata") {
        metadata.set("clinvar", provenance);
    }

    let body = json_text::dumps_with_options(
        &cache,
        DumpsOptions {
            indent: None,
            sort_keys: true,
            ensure_ascii: false,
            compact: true,
            allow_nan: true,
        },
    )
    .map_err(|e| error("ValueError", e))?;
    write_atomic(&options.out_json, &body)?;
    let injected: Vec<String> = counts
        .0
        .iter()
        .map(|(key, value)| format!("{key}={value}"))
        .collect();
    log(&format!("injected: {}", injected.join("; ")));

    if let Some(path) = &options.audit_tsv {
        write_audit(path, &matches)?;
    }
    if let Some(path) = &options.stats_json {
        let body = json_text::dumps_with_options(
            &stats,
            DumpsOptions {
                indent: Some(2),
                sort_keys: false,
                ensure_ascii: true,
                compact: false,
                allow_nan: true,
            },
        )
        .map_err(|e| error("ValueError", e))?;
        fs::write(path, body).map_err(|e| os_error(e, path, None))?;
    }
    Ok(())
}

/// `{**parse_stats, **counts}`.
fn stats_object(parse_stats: &BTreeMap<String, i128>, counts: &Counts) -> Json {
    Json::Object(
        parse_stats
            .iter()
            .map(|(key, value)| (key.clone(), *value))
            .chain(counts.0.iter().cloned())
            .map(|(key, value)| (key, Json::Int(value)))
            .collect(),
    )
}

/// The checks before the script's `parser.error`, which exits 2.
fn require_source(path: &Path) -> Result<()> {
    match fs::metadata(path) {
        Ok(meta) if meta.len() > 0 => Ok(()),
        _ => Err(error(
            "SystemExit",
            format!("ClinVar VCV is required and was not found: {}", path.display()),
        )),
    }
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(options) = parse_options(args)? else {
        return Ok(());
    };
    builder_error::finish(require_source(&options.clinvar_vcv))?;
    if let Some(digest) = &options.clinvar_vcv_sha256 {
        let lower_hex = |c: u8| c.is_ascii_digit() || (b'a'..=b'f').contains(&c);
        if digest.len() != 64 || !digest.bytes().all(lower_hex) {
            return Err(invalid(
                "--clinvar-vcv-sha256 must be a lowercase SHA256 digest",
            ));
        }
    }
    builder_error::finish(run(options))
}
