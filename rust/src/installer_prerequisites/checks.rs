//! The manual-input, main, bundled-snapshot, ClinVar statistics and GoFCards
//! stages of the installer prerequisite checks.
use super::mapping::{mapping_fetch_source, nm_to_enst_path, validate_mapping};
use super::paths::{
    TextLines, abspath, checked, checked_text, exists, is_dir, is_file, is_relative_to, is_symlink,
    join, lexists, nonempty_readable, open, parts, pathlib, read_failure, read_text, readable,
    readable_text, resolve, suffix,
};
use super::{Config, Result, Value, error, exception, require, value_error};
use crate::pickle_reader::{self, Value as Pickled};
use crate::python_io::{PythonGzipReader, resolve_strict};
use crate::python_path::{absolute, metadata, os_error, repr};
use crate::python_text::{decimal_digit_value, is_python_space};
use crate::reference_resources;
use crate::yaml_config;
use std::collections::HashMap;
use std::io::Read;
use std::path::Path;
use std::rc::Rc;

const WORKFLOW_ONLY_KEYS: [&str; 9] = [
    "input_vcf",
    "ped_file",
    "output_dir",
    "hub_file_dir",
    "alt_disease_vcf",
    "control_vcf",
    "pp1_vcf",
    "pp1_ped",
    "conda_env_yaml",
];
const MANUAL_INPUT_KEYS: [&str; 3] = [
    "primateai_prescore",
    "spliceai_snv_prescore",
    "spliceai_indel_prescore",
];
/// `str.splitlines()` boundaries.
const LINE_BREAKS: [char; 10] = [
    '\n', '\r', '\u{b}', '\u{c}', '\u{1c}', '\u{1d}', '\u{1e}', '\u{85}', '\u{2028}', '\u{2029}',
];
/// The exception classes `validate_clinvar_stats` turns into its own ValueError:
/// OSError, ValueError, EOFError and pickle.UnpicklingError, with subclasses.
const REUSE_FAILURES: [&str; 11] = [
    "OSError",
    "FileNotFoundError",
    "PermissionError",
    "IsADirectoryError",
    "NotADirectoryError",
    "FileExistsError",
    "BadGzipFile",
    "ValueError",
    "UnicodeDecodeError",
    "EOFError",
    "UnpicklingError",
];
const BUNDLED_CLINVAR: &str = "data/ClinVar/clinvar_2star_stats.pkl.gz";
const REGENERATION: &str = "Separate candidate regeneration requires explicit ClinVar \
    variant_summary and MANE inputs: run ClinVar_Gene_stat_regenerate CONFIG SUMMARY MANE \
    OUTPUT with a distinct output; this does not replace the baseline.";

/// `CHR1_LENGTH[assembly]`
fn chr1_length(assembly: &str) -> &'static str {
    if assembly == "hg19" {
        "249250621"
    } else {
        "248956422"
    }
}

/// `{"hg19": "hg19", "GRCh37": "hg19", "hg38": "hg38", "GRCh38": "hg38"}.get(value)`
pub(crate) fn assembly_alias(value: &Value) -> Result<Option<&'static str>> {
    Ok(match value {
        Value::Str(text) => match text.as_str() {
            "hg19" | "GRCh37" => Some("hg19"),
            "hg38" | "GRCh38" => Some("hg38"),
            _ => None,
        },
        Value::Collection(_) => {
            return Err(error(
                "TypeError",
                format!("unhashable type: '{}'", value.type_name()),
            ));
        }
        _ => None,
    })
}

/// `int(text)` in base 10 as a canonical decimal, or Python's ValueError.
pub(crate) fn python_int(text: &str) -> Result<String> {
    // Python converts other scripts' digits and white space to ASCII first.
    let converted: String = text
        .chars()
        .map(|c| match decimal_digit_value(c) {
            _ if c.is_ascii() => c,
            Some(digit) => char::from(b'0' + digit as u8),
            None if is_python_space(c) => ' ',
            None => c,
        })
        .collect();
    let body = converted.trim_matches([' ', '\t', '\n', '\u{b}', '\u{c}', '\r']);
    let (negative, digits) = match body.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, body.strip_prefix('+').unwrap_or(body)),
    };
    if digits.is_empty()
        || digits.starts_with('_')
        || digits.ends_with('_')
        || digits.contains("__")
        || !digits.bytes().all(|b| b.is_ascii_digit() || b == b'_')
    {
        // CPython formats the literal with %.200R.
        let shown: String = repr(text).chars().take(200).collect();
        return Err(value_error(format!(
            "invalid literal for int() with base 10: {shown}"
        )));
    }
    let digits: String = digits.chars().filter(|c| *c != '_').collect();
    if digits.len() > 4300 {
        return Err(value_error(format!(
            "Exceeds the limit (4300 digits) for integer string conversion: value has {} \
             digits; use sys.set_int_max_str_digits() to increase the limit",
            digits.len()
        )));
    }
    let trimmed = digits.trim_start_matches('0');
    Ok(match (trimmed.is_empty(), negative) {
        (true, _) => "0".to_owned(),
        (false, true) => format!("-{trimmed}"),
        (false, false) => trimmed.to_owned(),
    })
}

/// `str.splitlines()`
pub(crate) fn splitlines(text: &str) -> Vec<&str> {
    let mut lines = Vec::new();
    let mut start = 0;
    let mut characters = text.char_indices().peekable();
    while let Some((index, c)) = characters.next() {
        if LINE_BREAKS.contains(&c) {
            lines.push(&text[start..index]);
            start = index + c.len_utf8();
            if c == '\r' && characters.peek().is_some_and(|(_, next)| *next == '\n') {
                characters.next();
                start += 1;
            }
        }
    }
    if start < text.len() {
        lines.push(&text[start..]);
    }
    lines
}

/// `line.split(None, 1)`
fn split_once_on_space(line: &str) -> Vec<&str> {
    let rest = line.trim_start_matches(is_python_space);
    match rest.find(is_python_space) {
        _ if rest.is_empty() => Vec::new(),
        None => vec![rest],
        Some(index) => {
            let tail = rest[index..].trim_start_matches(is_python_space);
            if tail.is_empty() {
                vec![&rest[..index]]
            } else {
                vec![&rest[..index], tail]
            }
        }
    }
}

/// Whether a path part names an unfinished or staging copy.
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

/// `validate_manual_file(value, key)`: readiness only; no payload is read.
fn validate_manual_file(value: &Value, key: &str) -> Result<()> {
    let path = readable(value, key)?;
    let text = value.text().unwrap_or_default();
    // A final-looking link must not conceal an explicitly unfinished source.
    let logical = abspath(text)?;
    let resolved = resolve(&path)?;
    if parts(&logical)
        .into_iter()
        .chain(parts(&resolved))
        .any(unfinished)
    {
        return Err(value_error(format!(
            "Incomplete/staging manual input {key}: {text}; select the completed original file"
        )));
    }
    Ok(())
}

/// `manual_inputs(config_file, prompt=prompt)`. Prepared inputs never write the
/// configuration. At a terminal, every unset path is collected, then all are
/// checked before one atomic configuration rewrite.
pub(super) fn manual_inputs(config_path: &str, prompt: bool) -> Result<()> {
    let destination = resolve_strict(Path::new(config_path))
        .map_err(|(failure, named)| os_error(failure, &named, None))?;
    let destination = destination.display().to_string();
    let original = read_text(&destination)?;
    let mut config = Config::load(&original, true)?;
    let mut missing = Vec::new();
    for key in MANUAL_INPUT_KEYS {
        let value = config.get(key);
        let unset = match value {
            Value::None => true,
            Value::Str(text) => matches!(text.trim_matches(is_python_space), "" | "null"),
            _ => false,
        };
        if unset {
            missing.push(key);
        } else {
            validate_manual_file(value, key)?;
        }
    }
    if missing.is_empty() {
        return Ok(());
    }
    // SAFETY: isatty only inspects the process's standard input descriptor.
    if !prompt || unsafe { libc::isatty(0) } != 1 {
        return Err(value_error(format!(
            "Manual publisher inputs must be prepared before resource installation: {}; \
             set these paths in config.yaml, or use --stage manual-inputs --prompt from an \
             interactive terminal. No files were acquired or configuration changed",
            missing.join(", ")
        )));
    }
    eprintln!(
        "Select completed manual publisher files; all paths are checked before the \
         configuration is saved."
    );
    let mut updates = Vec::new();
    for key in &missing {
        let Some(answer) = input_line(&format!("{key}: "))? else {
            return Err(value_error(
                "Manual input collection ended; configuration is unchanged",
            ));
        };
        let answer = answer.trim_matches(is_python_space).to_owned();
        updates.push((key.to_string(), answer));
    }
    // Collect all missing paths, then validate them as one transaction. A later
    // invalid answer must not leave earlier answers saved in the configuration.
    for (key, value) in &updates {
        config.set(key, Value::Str(value.clone()));
    }
    for key in MANUAL_INPUT_KEYS {
        validate_manual_file(config.get(key), key)?;
    }
    let updated = crate::installer_config_paths::render_text(&original, &updates)?;
    // Refuse an edit the shared renderer cannot make without changing another
    // value (a hand-edited construct such as a quoted duplicate key).
    let reloaded = Config::load(&updated, true)?;
    let (before, after) = (config.entries(), reloaded.entries());
    require(
        before.len() == after.len() && before.iter().all(|entry| after.contains(entry)),
        "Manual input paths could not be updated without changing other configuration values",
    )?;
    crate::installer_config_paths::replace_text_synced(
        Path::new(&destination),
        &original,
        &updated,
        "Configuration changed during manual input collection; retry after its writer closes",
    )
}

/// `input(prompt)`: CPython writes the prompt to standard error when standard
/// input and output are both terminals (no readline module is loaded), and to
/// standard output otherwise. `None` is the end of input before any text.
fn input_line(prompt: &str) -> Result<Option<String>> {
    use std::io::{BufRead, Write};
    let os = |failure: std::io::Error| error("OSError", failure.to_string());
    // SAFETY: isatty only inspects the process's standard output descriptor.
    if unsafe { libc::isatty(1) } == 1 {
        let mut stream = std::io::stderr();
        stream.write_all(prompt.as_bytes()).map_err(os)?;
        stream.flush().map_err(os)?;
    } else {
        let mut stream = std::io::stdout();
        stream.write_all(prompt.as_bytes()).map_err(os)?;
        stream.flush().map_err(os)?;
    }
    let mut line = Vec::new();
    let mut input = std::io::stdin().lock();
    input.read_until(b'\n', &mut line).map_err(os)?;
    if line.is_empty() {
        return Ok(None);
    }
    if line.last() == Some(&b'\n') {
        line.pop();
    }
    match crate::python_io::decode_utf8(&line) {
        Ok(text) => Ok(Some(text.to_owned())),
        Err(message) => Err(error("UnicodeDecodeError", message)),
    }
}

/// `validate_main(config, source_root, ...)`
pub(super) fn validate_main(
    config: &Config,
    root: &str,
    allow_missing_mapping: bool,
    allow_missing_clinvar: bool,
) -> Result<()> {
    // Reject old-host values before opening any of their paths.
    let problems: Vec<String> = config
        .entries()
        .iter()
        .filter(|(key, _)| !WORKFLOW_ONLY_KEYS.contains(&key.as_str()))
        .filter_map(|(key, value)| match value {
            Value::Str(text) if text.starts_with("/paedyl01/") => Some(format!(
                "Unresolved old-host path {key}: {text}; supply its actual local path"
            )),
            _ => None,
        })
        .collect();
    require(problems.is_empty(), problems.join("\n"))?;
    let mapping = nm_to_enst_path(config, root)?;
    if allow_missing_mapping && config.get("nm_to_enst_tsv").is_unset() && !exists(&mapping)? {
        // Reject orphan provenance before deferring acquisition.
        mapping_fetch_source(config, root)?;
    } else {
        validate_mapping(config, root)?;
    }
    // A builder output during a fresh installation; only its path is checked.
    let selected_hpo = config.get("hpo_assertions");
    if selected_hpo.is_unset() {
        let expected = join(root, "data/hpo/genes_to_phenotype.assertions.tsv.gz");
        checked_text(&expected, "hpo_assertions", 0)?;
    } else {
        checked(selected_hpo, "hpo_assertions")?;
    }
    let baseline = clinvar_stats_path(config, root)?;
    let bundled = join(root, BUNDLED_CLINVAR);
    // Defer only the genuine default: an explicit external input or a
    // compatibility link must not be silently replaced with a new dataset.
    let missing_default = baseline == bundled
        && resolve(&baseline)? == absolute(Path::new(&bundled))?
        && (!exists(&baseline)? || (is_file(&baseline)? && file_size(&baseline)? == 0));
    if !(allow_missing_clinvar && missing_default) {
        validate_selected_clinvar_stats(config, root)?;
    }
    Ok(())
}

fn file_size(path: &str) -> Result<u64> {
    Ok(metadata(Path::new(path))?.len())
}

/// `validate_bundled_snapshots(config, source_root)`: existence and size only.
pub(super) fn validate_bundled_snapshots(config: &Config, root: &str) -> Result<()> {
    let assembly = assembly_alias(config.get("assembly"))?
        .ok_or_else(|| value_error("Select hg19/GRCh37 or hg38/GRCh38"))?;
    let configured = config.get("repeat_region_file_name");
    let repeat_name = if configured.truthy() {
        configured.clone()
    } else {
        Value::Str(format!("ucsc.{assembly}.fasta.repeat.bed.gz"))
    };
    let repeat_name = match &repeat_name {
        Value::Str(text)
            if !matches!(text.as_str(), "" | "." | "..")
                && !text.contains(['/', '\\', '\r', '\n']) =>
        {
            text.clone()
        }
        _ => {
            return Err(value_error(
                "repeat_region_file_name must be one filename beneath data/repeats, \
                 not an absolute or relative path",
            ));
        }
    };
    let rmc = format!("data/gnomAD_RMC/gnomad_v2.1.1_rmc_flat.{assembly}.tsv");
    let mut paths = vec![
        join(root, &format!("data/HCSeeker/HC_spots.{assembly}.tsv")),
        join(root, &rmc),
        join(root, "data/alphamissense/alphamissense_mean_score.tsv"),
        join(root, "data/vep/csq_interpretations.tsv"),
        join(&join(root, "data/repeats"), &repeat_name),
    ];
    for (key, relative) in [
        (
            "hgnc_table",
            "data/hgnc/hgnc_non_alt_loci_set_sha256_890b18743cc8.tsv",
        ),
        (
            "dispensable_gene_list",
            "data/dispensable_genes/dispensable_gene_list.txt",
        ),
        (
            "gofcards_mechanism_review_tsv",
            "data/gofcards/gofcards_mechanism_reviews.tsv",
        ),
        (
            "mondo_disease_scope_overrides",
            "data/mondo/disease_scope_overrides.tsv",
        ),
        (
            "gene_nonlof_mechanism_schema",
            "data/patho_mechanism/gene_nonlof_mechanism_curated_assertions.schema.json",
        ),
    ] {
        let selected = config.get(key);
        let path = if selected.is_unset() {
            join(root, relative)
        } else {
            join(root, &selected.py_str()?)
        };
        // An explicit review is a build-only input: GoFCards requires it before
        // any acquisition when rebuilding, not for reuse of a current cache.
        if key == "gofcards_mechanism_review_tsv"
            && !selected.is_unset()
            && abspath(&path)? != abspath(&join(root, relative))?
            && !exists(&path)?
            && !is_symlink(&path)?
        {
            continue;
        }
        paths.push(path);
    }
    for path in &paths {
        readable_text(path, "bundled snapshot or explicitly selected replacement")?;
    }
    Ok(())
}

/// `clinvar_stats_path(config, source_root)`: a bare legacy name is
/// bundled-relative; absolute paths and relative paths with folders are kept.
fn clinvar_stats_path(config: &Config, root: &str) -> Result<String> {
    let configured = config.get("clinvar_gene_stat");
    let selected = if configured.truthy() {
        configured.py_str()?
    } else {
        "clinvar_2star_stats.pkl.gz".to_owned()
    };
    let path = pathlib(&selected);
    Ok(if path.starts_with('/') || selected.contains('/') {
        path
    } else {
        join(&join(root, "data/ClinVar"), &path)
    })
}

/// `validate_selected_clinvar_stats(config, source_root)`: bundled original data
/// need only exist; supplied or external bytes keep their checks.
pub(super) fn validate_selected_clinvar_stats(config: &Config, root: &str) -> Result<()> {
    let path = clinvar_stats_path(config, root)?;
    let bundled = join(root, BUNDLED_CLINVAR);
    let source_data = join(&resolve(root)?, "data/ClinVar");
    if path == bundled && is_relative_to(&resolve(&path)?, &source_data) {
        return readable_text(&path, "bundled clinvar_gene_stat")
            .map(|_| ())
            .map_err(|failure| match exception(failure.as_ref()) {
                Some(("ValueError", message)) => value_error(format!(
                    "{message}. Restore the nonempty original bundled artifact for baseline \
                     reuse. {REGENERATION}"
                )),
                _ => failure,
            });
    }
    validate_clinvar_stats(&path)
}

/// `validate_clinvar_stats(path)`
fn validate_clinvar_stats(path: &str) -> Result<()> {
    let checked = if suffix(path) == ".rkyv" {
        validate_native_clinvar_stats(path)
    } else {
        validate_clinvar_pickle(path)
    };
    checked.map_err(|failure| match exception(failure.as_ref()) {
        Some((category, message)) if REUSE_FAILURES.contains(&category) => value_error(format!(
            "Cannot reuse bundled ClinVar gene statistics {path}: {message}. Restore a \
             validated original artifact for baseline reuse. {REGENERATION}"
        )),
        _ => failure,
    })
}

/// The PRIVARKV header of a gene-statistics archive (rust/CACHE_FORMATS.md 2).
fn validate_native_clinvar_stats(path: &str) -> Result<()> {
    readable_text(path, "clinvar_gene_stat archive")?;
    let mut header = Vec::with_capacity(32);
    open(path)?
        .take(32)
        .read_to_end(&mut header)
        .map_err(|failure| os_error(failure, Path::new(path), None))?;
    require(
        header.len() == 32 && header.starts_with(b"PRIVARKV"),
        "not a PRIVARKV archive",
    )?;
    let word = |at: usize| u32::from_le_bytes(header[at..at + 4].try_into().unwrap());
    let (container, schema) = (word(8), word(12));
    require(
        (container, schema) == (1, 5),
        format!("container {container}, schema {schema} is not a ClinVar gene-statistics archive"),
    )?;
    let length = u64::from_le_bytes(header[24..32].try_into().unwrap());
    require(
        32 + u128::from(length) == u128::from(metadata(Path::new(path))?.len()),
        "archive length differs from its header",
    )
}

/// A gzip pickle of gene statistics: a nonempty dict of dicts holding the five
/// fields. The restricted reader reads plain data only, so a pickle that needs
/// a reducer (sets in protocol 2 or 3) is refused rather than evaluated.
fn validate_clinvar_pickle(path: &str) -> Result<()> {
    readable_text(path, "bundled clinvar_gene_stat")?;
    let mut bytes = Vec::new();
    PythonGzipReader::new(open(path)?)
        .read_to_end(&mut bytes)
        .map_err(|failure| read_failure(Box::new(failure), path))?;
    let stats = pickle_reader::load(&bytes)
        .map_err(|failure| error("UnpicklingError", failure.to_string()))?;
    let fields = [
        "ensembl_id",
        "clinvar_pathogenic",
        "large_aachange",
        "small_aachange",
        "putative_nmd_variants",
    ];
    let complete = match &stats {
        Pickled::Dict(genes) => {
            let genes = genes.borrow();
            !genes.entries.is_empty()
                && genes.entries.iter().all(|(_, gene)| match gene {
                    Pickled::Dict(gene) => fields.iter().all(|field| {
                        let key = Pickled::Str(Rc::from(*field));
                        gene.borrow().get(&key).is_ok_and(|found| found.is_some())
                    }),
                    _ => false,
                })
        }
        _ => false,
    };
    require(complete, "missing gene statistics fields")
}

/// `validate_fasta(value, assembly)`: the FASTA and its .fai identify the
/// assembly's chromosome 1 length.
fn validate_fasta(value: &Value, assembly: &str) -> Result<()> {
    let key = format!("gofcards_{assembly}_fasta");
    let fasta = readable(value, &key)?;
    let index = readable_text(
        &format!("{fasta}.fai"),
        &format!("{key} index (run samtools faidx)"),
    )?;
    let mut lines = TextLines::new(open(&index)?, &index);
    let mut lengths = HashMap::new();
    while let Some(line) = lines.readline()? {
        let line = line.trim_end_matches(is_python_space);
        let row: Vec<&str> = line.split('\t').collect();
        if row.len() >= 2 {
            lengths.insert(row[0].to_owned(), python_int(row[1])?);
        }
    }
    let observed = lengths.get("chr1").or_else(|| lengths.get("1"));
    require(
        observed.map(String::as_str) == Some(chr1_length(assembly)),
        format!(
            "{key} .fai does not identify the required {assembly} chromosome 1 length: {}",
            value.text().unwrap_or_default()
        ),
    )
}

/// `line.split()` of the first line starting with "chain ", or no words.
fn chain_header<R: Read>(mut lines: TextLines<R>) -> Result<Vec<String>> {
    while let Some(line) = lines.readline()? {
        if line.starts_with("chain ") {
            return Ok(line
                .split(is_python_space)
                .filter(|word| !word.is_empty())
                .map(str::to_owned)
                .collect());
        }
    }
    Ok(Vec::new())
}

/// `validate_forward_chain(value)`: the first chain header is UCSC's forward
/// hg19-to-hg38 chromosome 1 chain.
fn validate_forward_chain(value: &Value) -> Result<()> {
    let path = readable(value, "gofcards_hg19_to_hg38_chain")?;
    let source = open(&path)?;
    // UCSC's target is the source coordinate system consumed by CrossMap.
    let header = if suffix(&path) == ".gz" {
        chain_header(TextLines::new(PythonGzipReader::new(source), &path))?
    } else {
        chain_header(TextLines::new(source, &path))?
    };
    let unprefixed = |index: usize| {
        let word = header[index].as_str();
        word.strip_prefix("chr").unwrap_or(word)
    };
    let mut forward = header.len() == 13 && unprefixed(2) == "1" && unprefixed(7) == "1";
    if forward {
        // Both lengths are converted before they are compared, as in the tuple.
        let source_length = python_int(&header[3])?;
        let target_length = python_int(&header[8])?;
        forward = source_length == chr1_length("hg19") && target_length == chr1_length("hg38");
    }
    require(
        forward,
        format!(
            "gofcards_hg19_to_hg38_chain must be the forward UCSC hg19-to-hg38 chain: {}",
            value.text().unwrap_or_default()
        ),
    )
}

/// `validate_merged_cache(directory, version, assembly, allow_missing=...)`
fn validate_merged_cache(
    directory: &Value,
    version: &str,
    assembly: &str,
    allow_missing: bool,
) -> Result<()> {
    require(
        crate::conservation_resource::python_isdigit(version),
        "Set vep_installed_version to the selected numeric VEP release",
    )?;
    let root = checked(directory, "vep_cache_dir")?;
    let merged = join(&root, "homo_sapiens_merged");
    let selected = join(&merged, &format!("{version}_{assembly}"));
    let cache = checked_text(&selected, "GoFCards merged VEP cache", 0)?;
    for parent in [pathlib(directory.text().unwrap_or_default()), merged] {
        require(
            !lexists(&parent) || is_dir(&parent)?,
            format!("GoFCards merged VEP cache parent is not a directory: {parent}"),
        )?;
    }
    if allow_missing && !lexists(&selected) {
        return Ok(());
    }
    for relative in ["info.txt", "chr_synonyms.txt", "1/all_vars.gz"] {
        readable_text(
            &join(&cache, relative),
            &format!("GoFCards {assembly} merged VEP {version} cache ({relative})"),
        )?;
    }
    let mut indexes = Vec::new();
    for name in ["1/all_vars.gz.csi", "1/all_vars.gz.tbi"] {
        let index = checked_text(&join(&cache, name), "GoFCards merged VEP index", 0)?;
        indexes.push(index);
    }
    let mut indexed = false;
    for index in &indexes {
        if nonempty_readable(index)? {
            indexed = true;
            break;
        }
    }
    require(
        indexed,
        format!(
            "GoFCards {assembly} merged VEP {version} cache lacks a converted variant index: \
             {cache}"
        ),
    )?;
    let information = read_text(&join(&cache, "info.txt"))?;
    let mut recorded = HashMap::new();
    for line in splitlines(&information) {
        if let [name, value] = split_once_on_space(line)[..] {
            recorded.insert(name, value.trim_matches(is_python_space));
        }
    }
    let matching = [
        ("assembly", assembly),
        ("species", "homo_sapiens"),
        ("var_type", "tabix"),
    ]
    .iter()
    .all(|(name, expected)| recorded.get(name) == Some(expected));
    require(
        matching,
        format!(
            "GoFCards {assembly} merged VEP {version} cache has mismatched \
             assembly/species/converted metadata: {cache}"
        ),
    )
}

/// `str(Path(config["cache_root"]) / (config.get("cache_folders") or {}).get(...))`
fn managed_vep_cache(config: &Config) -> Result<String> {
    let root = match config.get("cache_root") {
        Value::Str(text) => pathlib(text),
        other => {
            return Err(error(
                "TypeError",
                // CPython 3.12 pathlib.Path wording.
                format!(
                    "argument should be a str or an os.PathLike object where __fspath__ returns \
                     a str, not '{}'",
                    other.type_name()
                ),
            ));
        }
    };
    let folders = config.get("cache_folders");
    let mut folder = Value::Str("VEP_cache".to_owned());
    if folders.truthy() {
        let entries = match folders {
            Value::Collection(text) => yaml_config::nested_mapping_entries(text),
            _ => None,
        };
        let Some(entries) = entries else {
            return Err(error(
                "AttributeError",
                format!("'{}' object has no attribute 'get'", folders.type_name()),
            ));
        };
        if let Some((_, value)) = entries.iter().rev().find(|(key, _)| key == "vep_cache") {
            folder = Value::scalar(value)?;
        }
    }
    match folder {
        Value::Str(name) => Ok(join(&root, &name)),
        other => Err(error(
            "TypeError",
            format!(
                "unsupported operand type(s) for /: 'PosixPath' and '{}'",
                other.type_name()
            ),
        )),
    }
}

/// `value == other` for the configuration values compared here.
fn same(value: &Value, other: &Value) -> bool {
    match (value, other) {
        (Value::None, Value::None) => true,
        (Value::Str(a), Value::Str(b)) => a == b,
        _ => false,
    }
}

/// `validate_gofcards(config, main_preflight=..., ...)`
pub(super) fn validate_gofcards(
    config: &Config,
    root: &str,
    main_preflight: bool,
    allow_pinned_reference_defaults: bool,
    allow_missing_merged_caches: bool,
) -> Result<()> {
    let selected = assembly_alias(config.get("assembly"))?;
    if allow_pinned_reference_defaults {
        // This resolves and checks every explicit reference without creating
        // directories, indexes or downloads. Only pinned missing defaults wait
        // for preparation; an original DRAGEN path is never a UCSC destination.
        let (resources, _) = reference_resources::resource_plan(
            config,
            &join(root, "data"),
            reference_resources::Scope::Gofcards,
            true,
        )?;
        reference_resources::preflight(&resources)?;
    } else {
        for assembly in ["hg19", "hg38"] {
            let key = format!("gofcards_{assembly}_fasta");
            let value = config.get(&key);
            if main_preflight
                && selected == Some(assembly)
                && same(value, config.get("ref_genome"))
            {
                // The selected reference may still be installed by the main run.
                let fasta = checked(value, &key)?;
                let index = checked_text(&format!("{fasta}.fai"), &format!("{key} index"), 0)?;
                if is_file(&fasta)? && is_file(&index)? {
                    validate_fasta(value, assembly)?;
                }
            } else {
                validate_fasta(value, assembly)?;
            }
        }
        require(
            !same(
                config.get("gofcards_hg19_fasta"),
                config.get("gofcards_hg38_fasta"),
            ),
            "GoFCards hg19 and hg38 FASTAs must be distinct assembly-correct inputs",
        )?;
        validate_forward_chain(config.get("gofcards_hg19_to_hg38_chain"))?;
    }
    let directory = if main_preflight && config.get("cache_root").truthy() {
        Value::Str(managed_vep_cache(config)?)
    } else {
        config.get("vep_cache_dir").clone()
    };
    let release = config.get("vep_installed_version");
    let version = if release.truthy() {
        release.py_str()?
    } else {
        String::new()
    };
    for assembly in ["GRCh37", "GRCh38"] {
        let smoke = std::env::var_os("PRIVA_SMOKE_SKIP_VEP").is_some_and(|value| value == "1");
        validate_merged_cache(
            &directory,
            &version,
            assembly,
            allow_missing_merged_caches && !smoke,
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn python_int_reads_what_int_reads() {
        assert_eq!(python_int(" 0248956422 ").unwrap(), "248956422");
        assert_eq!(python_int("2_4").unwrap(), "24");
        assert_eq!(python_int("-0").unwrap(), "0");
        assert_eq!(python_int("\u{ff11}\u{ff12}").unwrap(), "12");
        for bad in ["", "1__2", "_1", "x", "\u{b2}", "1\u{1c}"] {
            let message = python_int(bad).err().unwrap().to_string();
            assert!(message.contains("invalid literal for int()"), "{bad}");
        }
        let message = python_int(&"9".repeat(4301)).err().unwrap().to_string();
        assert!(message.contains("Exceeds the limit (4300 digits)"));
    }

    #[test]
    fn line_and_word_splitting_follow_str_methods() {
        assert_eq!(splitlines("a\r\nb\rc\u{2028}d\n"), ["a", "b", "c", "d"]);
        assert_eq!(splitlines("a\n\nb"), ["a", "", "b"]);
        assert_eq!(
            split_once_on_space("  key  value with space "),
            ["key", "value with space "]
        );
        assert_eq!(split_once_on_space("key   "), ["key"]);
        assert!(split_once_on_space("   ").is_empty());
    }

    #[test]
    fn staging_names_are_recognized_case_insensitively() {
        for part in [
            "x.PARTIAL",
            ".building",
            ".building-1",
            ".download.x",
            "a.part",
            "a.TMP",
        ] {
            assert!(unfinished(part), "{part}");
        }
        for part in ["complete.tsv.gz", "building", "partial"] {
            assert!(!unfinished(part), "{part}");
        }
    }

    #[test]
    fn unhashable_assembly_is_a_type_error() {
        let failure = assembly_alias(&Value::Collection("[a]".into()))
            .err()
            .unwrap();
        assert_eq!(failure.to_string(), "TypeError: unhashable type: 'list'");
        assert_eq!(assembly_alias(&Value::Int("19".into())).unwrap(), None);
        assert_eq!(
            assembly_alias(&Value::Str("GRCh38".into())).unwrap(),
            Some("hg38")
        );
    }
}
