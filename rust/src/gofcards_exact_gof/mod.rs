//! Native `scripts/build_gofcards_exact_gof_cache.py`: the GoFCards exact
//! gain-of-function cache builder with the script's options.
//!
//! The workbook download, the summary cache and its retries, the GRCh37
//! reference check, the bcftools, VEP and CrossMap runs, the endpoint
//! recovery of unplaced variants, the HGNC and review rules and the export
//! follow the script. Intermediate files, the rejects table, the stats file,
//! the `[gofcards]` log lines and exit statuses follow it too; the gzip output
//! is the same zlib stream, with the time of writing in its header as
//! `gzip.open` leaves it. An exception the script let escape exits 1 here as
//! `priva: CATEGORY: message`.
//!
//! Differences:
//! - Requests run through `curl`, so a transport failure carries curl's message.
//! - `--summary-endpoint URL` replaces the script's constant; tests use it.
//! - A workbook with two columns that map to one source field is refused, and
//!   so is a compressed FASTA.
mod build;
mod fasta;
mod fetch;

use crate::builder_error;
use crate::cli::{Result, invalid};
use crate::gene_nonlof_mechanism::error;
use crate::json_text::{self, DumpsOptions, Json};
use crate::python_args::parse_long_options_converted;
use crate::python_path::{os_error, path_value};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const HELP: &str = "usage: priva gofcards-exact-gof --workdir WORKDIR --out-json OUT_JSON\n\
       --hg19-fasta HG19_FASTA --hg38-fasta HG38_FASTA --chain CHAIN\n\
       --vep-cache-dir VEP_CACHE_DIR [--vep-cache-version VEP_CACHE_VERSION]\n\
       --hgnc-table HGNC_TABLE [--mechanism-review-tsv MECHANISM_REVIEW_TSV]\n\
       [--public-excel-url PUBLIC_EXCEL_URL] [--api-timeout API_TIMEOUT]\n\
       [--stats-json STATS_JSON] [--summary-endpoint URL]\n\
Build the PriVA GoFCards exact gain-of-function evidence cache; the options, outputs\n\
and exit statuses of scripts/build_gofcards_exact_gof_cache.py.";

const VALUES: &[&str] = &[
    "workdir",
    "out-json",
    "out-tsv",
    "hg19-fasta",
    "hg38-fasta",
    "chain",
    "vep-cache-dir",
    "vep-cache-version",
    "hgnc-table",
    "mechanism-review-tsv",
    "public-excel-url",
    "api-timeout",
    "stats-json",
    "summary-endpoint",
];

const PUBLIC_EXCEL_URL: &str =
    "https://download.genemed.tech/upload/GainFunCards/gofcards_data_download.xlsx";
const SUMMARY_ENDPOINT: &str =
    "https://java.genemed.tech/admin-api/backend/data/hg19/variantLevel/summary";
const VEP_FIELDS: &str = "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,\
Consequence,HGVSc,HGVSp,MANE_SELECT,MANE_PLUS_CLINICAL,CANONICAL,SYMBOL";

fn log(message: &str) {
    eprintln!("[gofcards] {message}");
}

/// `clean_text` for text.
fn clean_text(value: &str) -> String {
    let text = crate::python_text::python_strip(value);
    if matches!(
        text.to_lowercase().as_str(),
        "nan" | "none" | "_" | "." | "-"
    ) {
        String::new()
    } else {
        text.to_owned()
    }
}

/// `clean_text` for a parsed JSON value.
fn clean_json(value: Option<&Json>) -> String {
    match value {
        None | Some(Json::Null) => String::new(),
        Some(Json::Float(number)) if number.is_nan() => String::new(),
        Some(other) => clean_text(&crate::gene_nonlof_mechanism::py_str(other)),
    }
}

/// `normalize_chrom` for text `clean_text` has seen.
fn normalize_chrom(clean: &str) -> String {
    match clean.get(..3) {
        Some(prefix) if prefix.eq_ignore_ascii_case("chr") => clean[3..].to_owned(),
        _ => clean.to_owned(),
    }
}

/// `normalize_allele` for text `clean_text` has seen.
fn normalize_allele(clean: &str) -> String {
    clean.to_uppercase()
}

/// The tail of `normalize_int`, after `str(value)`.
fn strip_point_zero(value: &str) -> String {
    let text = crate::python_text::python_strip(value);
    text.strip_suffix(".0").unwrap_or(text).to_owned()
}

/// `str(int(value))` for a whole float.
fn whole_float(value: f64) -> String {
    if value == 0.0 {
        "0".to_owned()
    } else {
        format!("{value:.0}")
    }
}

/// `normalize_int` for a parsed JSON value.
fn normalize_int_json(value: Option<&Json>) -> String {
    match value {
        None | Some(Json::Null) => String::new(),
        Some(Json::Float(number)) if number.is_nan() => String::new(),
        Some(Json::Float(number)) if number.is_finite() && number.fract() == 0.0 => {
            whole_float(*number)
        }
        Some(other) => strip_point_zero(&crate::gene_nonlof_mechanism::py_str(other)),
    }
}

/// `str.isdigit()`.
fn is_digits(text: &str) -> bool {
    !text.is_empty() && text.chars().all(crate::gene_nonlof_mechanism::is_digit)
}

/// `int(text)`.
fn python_int(text: &str) -> Result<i128> {
    crate::python_text::python_int(text).ok_or_else(|| {
        error(
            "ValueError",
            format!(
                "invalid literal for int() with base 10: {}",
                crate::python_path::repr(text)
            ),
        )
    })
}

/// `work / name`.
fn child(work: &Path, name: &str) -> PathBuf {
    path_value(&work.join(name).to_string_lossy())
}

/// The text of a file opened with `open(path, "r", encoding="utf-8")`.
fn read_text(path: &Path) -> Result<String> {
    let bytes = fs::read(path).map_err(|e| os_error(e, path, None))?;
    let text =
        crate::python_io::decode_utf8(&bytes).map_err(|e| error("UnicodeDecodeError", e))?;
    Ok(text.replace("\r\n", "\n").replace('\r', "\n"))
}

/// `list(handle)` of the same file, each line keeping its `\n`.
fn read_lines(path: &Path) -> Result<Vec<String>> {
    Ok(read_text(path)?
        .split_inclusive('\n')
        .map(str::to_owned)
        .collect())
}

/// `run_command`.
fn run_command(command: &[String], log_path: &Path) -> Result<()> {
    if let Some(parent) = log_path.parent() {
        fs::create_dir_all(parent).map_err(|e| os_error(e, parent, None))?;
    }
    let shown = command.join(" ");
    log(&format!("$ {shown}"));
    let file = File::create(log_path).map_err(|e| os_error(e, log_path, None))?;
    let duplicate = file.try_clone().map_err(|e| os_error(e, log_path, None))?;
    let status = Command::new(&command[0])
        .args(&command[1..])
        .stdin(Stdio::inherit())
        .stdout(file)
        .stderr(duplicate)
        .status()
        .map_err(|e| os_error(e, Path::new(&command[0]), None))?;
    if !status.success() {
        use std::os::unix::process::ExitStatusExt;
        let code = status
            .code()
            .unwrap_or_else(|| -status.signal().unwrap_or(0));
        return Err(error(
            "RuntimeError",
            format!(
                "Command failed ({code}); see {}: {shown}",
                log_path.display()
            ),
        ));
    }
    Ok(())
}

/// `require_executable`.
fn require_executable(name: &str) -> Result<String> {
    crate::spliceai2_prepare::which(name).ok_or_else(|| {
        error(
            "RuntimeError",
            format!(
                "Required executable {} is not on PATH. Install it into the PriVA \
environment rather than working around it.",
                crate::python_path::repr(name)
            ),
        )
    })
}

fn text(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

/// `normalize_vcf`.
fn normalize_vcf(input: &Path, fasta: &Path, output: &Path, logs: &Path) -> Result<()> {
    require_executable("bcftools")?;
    let stem = output
        .file_stem()
        .map(|stem| stem.to_string_lossy().into_owned())
        .unwrap_or_default();
    run_command(
        &[
            "bcftools".to_owned(),
            "norm".to_owned(),
            "-f".to_owned(),
            text(fasta),
            "-c".to_owned(),
            "e".to_owned(),
            "-m-any".to_owned(),
            "-Ov".to_owned(),
            "-o".to_owned(),
            text(output),
            text(input),
        ],
        &child(logs, &format!("bcftools_norm.{stem}.log")),
    )
}

/// `run_vep`.
fn run_vep(
    input: &Path,
    assembly: &str,
    options: &Options,
    fasta: &Path,
    output: &Path,
    logs: &Path,
) -> Result<()> {
    let vep = require_executable("vep")?;
    let mut command: Vec<String> = vec![
        vep,
        "--input_file".into(),
        text(input),
        "--output_file".into(),
        text(output),
        "--format".into(),
        "vcf".into(),
        "--tab".into(),
        "--force_overwrite".into(),
        "--assembly".into(),
        assembly.into(),
        "--fasta".into(),
        text(fasta),
        "--symbol".into(),
        "--hgvs".into(),
        "--mane".into(),
        "--canonical".into(),
        "--transcript_version".into(),
        "--fields".into(),
        VEP_FIELDS.into(),
        "--cache".into(),
        "--offline".into(),
        "--dir_cache".into(),
        text(&options.vep_cache_dir),
    ];
    if !options.vep_cache_version.is_empty() {
        command.push("--cache_version".into());
        command.push(options.vep_cache_version.clone());
    }
    command.push("--merged".into());
    run_command(&command, &child(logs, &format!("vep.{assembly}.log")))
}

/// `liftover_vcf`: the records CrossMap could not place.
fn liftover_vcf(
    input: &Path,
    options: &Options,
    output: &Path,
    logs: &Path,
) -> Result<HashSet<String>> {
    require_executable("CrossMap")?;
    run_command(
        &[
            "CrossMap".to_owned(),
            "vcf".to_owned(),
            text(&options.chain),
            text(input),
            text(&options.hg38_fasta),
            text(output),
        ],
        &child(logs, "crossmap.log"),
    )?;
    let mut unmapped = HashSet::new();
    let unmap_path = PathBuf::from(format!("{}.unmap", output.display()));
    if unmap_path.exists() {
        for line in read_lines(&unmap_path)? {
            if line.starts_with('#') {
                continue;
            }
            let body = line.trim_end_matches('\n');
            if let Some(id) = body.split('\t').nth(2).filter(|id| !id.is_empty()) {
                unmapped.insert(id.to_owned());
            }
        }
    }
    log(&format!(
        "liftover: {} records could not be placed on GRCh38",
        unmapped.len()
    ));
    Ok(unmapped)
}

struct Options {
    workdir: PathBuf,
    out_json: PathBuf,
    hg19_fasta: PathBuf,
    hg38_fasta: PathBuf,
    chain: PathBuf,
    vep_cache_dir: PathBuf,
    vep_cache_version: String,
    hgnc_table: PathBuf,
    mechanism_review: Option<PathBuf>,
    public_excel_url: String,
    summary_endpoint: String,
    api_timeout: i128,
    stats_json: Option<PathBuf>,
}

fn parse_options(args: &[String]) -> Result<Option<Options>> {
    let converted = parse_long_options_converted(args, VALUES, &[], |name, value| {
        if name == "api-timeout" && crate::python_text::python_int(value).is_none() {
            return Err(format!(
                "argument --api-timeout: invalid int value: {}",
                crate::python_path::repr(value)
            ));
        }
        Ok(())
    })
    .map_err(invalid)?;
    let Some(occurrences) = converted else {
        println!("{HELP}");
        return Ok(None);
    };
    let mut values: HashMap<String, String> = HashMap::new();
    for (name, value) in occurrences {
        let name = if name == "out-tsv" {
            "out-json".to_owned()
        } else {
            name
        };
        values.insert(name, value);
    }
    let missing: Vec<&str> = [
        ("workdir", "--workdir"),
        ("out-json", "--out-json/--out-tsv"),
        ("hg19-fasta", "--hg19-fasta"),
        ("hg38-fasta", "--hg38-fasta"),
        ("chain", "--chain"),
        ("vep-cache-dir", "--vep-cache-dir"),
        ("hgnc-table", "--hgnc-table"),
    ]
    .into_iter()
    .filter(|(name, _)| !values.contains_key(*name))
    .map(|(_, shown)| shown)
    .collect();
    if !missing.is_empty() {
        return Err(invalid(format!(
            "the following arguments are required: {}",
            missing.join(", ")
        )));
    }
    let path = |name: &str| path_value(&values[name]);
    let optional = |name: &str| values.get(name).map(|value| path_value(value));
    let text_or = |name: &str, default: &str| {
        values
            .get(name)
            .cloned()
            .unwrap_or_else(|| default.to_owned())
    };
    Ok(Some(Options {
        workdir: path("workdir"),
        out_json: path("out-json"),
        hg19_fasta: path("hg19-fasta"),
        hg38_fasta: path("hg38-fasta"),
        chain: path("chain"),
        vep_cache_dir: path("vep-cache-dir"),
        vep_cache_version: text_or("vep-cache-version", ""),
        hgnc_table: path("hgnc-table"),
        mechanism_review: optional("mechanism-review-tsv"),
        public_excel_url: text_or("public-excel-url", PUBLIC_EXCEL_URL),
        summary_endpoint: text_or("summary-endpoint", SUMMARY_ENDPOINT),
        api_timeout: values
            .get("api-timeout")
            .and_then(|value| crate::python_text::python_int(value))
            .unwrap_or(300),
        stats_json: optional("stats-json"),
    }))
}

fn run(options: Options) -> Result<()> {
    let work = options.workdir.clone();
    let logs = child(&work, "logs");
    fs::create_dir_all(&logs).map_err(|e| os_error(e, &logs, None))?;
    let at = |name: &str| child(&work, name);

    // Step 1
    let (sources, mut provenance) =
        fetch::fetch_sources(&work, &options.public_excel_url, options.api_timeout)?;
    let summary = fetch::fetch_summary_annotations(
        &sources,
        &at("gofcards_summary_cache.jsonl"),
        &options.summary_endpoint,
        options.api_timeout,
        &work,
    )?;

    // Step 2
    let core = build::build_source_vcf(
        &sources,
        &options.hg19_fasta,
        &at("gofcards.hg19.vcf"),
        &at("gofcards_reference_rejects.tsv"),
    )?;

    // Step 3
    let norm19 = at("gofcards.hg19.norm.vcf");
    normalize_vcf(&at("gofcards.hg19.vcf"), &options.hg19_fasta, &norm19, &logs)?;

    // Step 4
    let vep19_path = at("gofcards.hg19.vep.tsv");
    run_vep(
        &norm19,
        "GRCh37",
        &options,
        &options.hg19_fasta,
        &vep19_path,
        &logs,
    )?;
    let vep19 = build::parse_vep(&vep19_path, "hg19")?;

    // Step 5
    let lifted = at("gofcards.hg38.vcf");
    let unmapped = liftover_vcf(&norm19, &options, &lifted, &logs)?;
    let norm38 = at("gofcards.hg38.norm.vcf");
    normalize_vcf(&lifted, &options.hg38_fasta, &norm38, &logs)?;

    // Step 6
    let vep38_path = at("gofcards.hg38.vep.tsv");
    run_vep(
        &norm38,
        "GRCh38",
        &options,
        &options.hg38_fasta,
        &vep38_path,
        &logs,
    )?;
    let vep38 = build::parse_vep(&vep38_path, "hg38")?;

    // Step 7
    provenance.push((
        "vep_cache_dir".to_owned(),
        Json::Str(text(&options.vep_cache_dir)),
    ));
    provenance.push((
        "vep_cache_version".to_owned(),
        Json::Str(options.vep_cache_version.clone()),
    ));
    provenance.push(("liftover_chain".to_owned(), Json::Str(text(&options.chain))));
    let lifted_coords = build::read_vcf_coordinates(&norm38)?;
    let placed: HashSet<&str> = lifted_coords.iter().map(|c| c.vcf_id.as_str()).collect();
    if core.is_empty() {
        return Err(error("KeyError", "'vcf_id'"));
    }
    let unplaced: Vec<&build::Core> = core
        .iter()
        .filter(|record| !placed.contains(record.vcf_id.as_str()))
        .collect();
    let (recovered, counts) =
        build::endpoint_hg38_coordinates(&unplaced, &summary, &options.hg38_fasta)?;
    log(&format!(
        "GRCh38 coordinates: {} from the liftover, {} recovered from the GoFCards endpoint",
        lifted_coords.len(),
        recovered.len()
    ));
    for (key, value) in counts {
        provenance.push((format!("hg38_{key}"), Json::Int(value)));
    }
    let mut hg38_coords = lifted_coords;
    hg38_coords.extend(recovered);
    let hg19_coords = build::read_vcf_coordinates(&norm19)?;
    let reviews = build::load_mechanism_reviews(options.mechanism_review.as_deref())?;
    let resolver = build::HgncResolver::new(&options.hgnc_table)?;
    let stats = build::build_cache(build::Layers {
        sources: &sources,
        core: &core,
        vep: [vep19, vep38],
        hg19_coords,
        hg38_coords,
        unmapped: &unmapped,
        reviews: &reviews,
        summary: &summary,
        resolver: &resolver,
        provenance: &provenance,
        output: &options.out_json,
    })?;
    if let Some(path) = &options.stats_json {
        let mut all = provenance;
        all.extend(stats);
        let body = json_text::dumps_with_options(
            &Json::Object(all),
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

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(options) = parse_options(args)? else {
        return Ok(());
    };
    builder_error::finish(run(options))
}
