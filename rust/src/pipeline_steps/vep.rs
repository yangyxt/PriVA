//! The VEP step without `anno_VEP_data` of annotation_vcf.sh: the checks of
//! every cache and plugin file, priva's SpliceAI version check, then `vep`
//! with the function's command line (the plugins and their files, MaveDB on
//! GRCh38, the pext tracks when `pext_enabled`), then `bcftools sort` and
//! `tabix` of its output in place of the VCF.
//!
//! The function skipped a VCF that already had CSQ and was newer than the
//! script; here, newer than the priva program that runs the step.
//!
//! [`worker`] is `anno_VEP_cache_worker` of the same script, which install
//! steps D1_chr1 to D1_chrOther ran on each chromosome of the ClinVar
//! release: the same checks and command line (one VEP command for the
//! pipeline and for ClinVar), with the output, statistics file and forks it
//! is given, without the SpliceAI version check, the skip and the sort.
use super::annotation::vcf_assembly;
use super::tools::{
    Context, Result, args, is_file, log, newer, quoted_path, random_tag, refuse, replace_vcf, words,
};
use std::os::unix::fs::MetadataExt;
use std::path::{Component, Path};

const MAVEDB_COLUMNS: &str =
    "urn:score:nt:pro:se:sd:df:epsilon:vtype:high_conf:pvalue:ci95_lower:ci95_upper";

/// `$(basename ${VALUE})` with VALUE unquoted: one word is its last path
/// part, two words drop the second as a suffix, none or more is empty
/// (basename's usage error).
fn basename(value: &str) -> String {
    let words = words(value);
    let last = |path: &str| {
        let trimmed = path.trim_end_matches('/');
        if trimmed.is_empty() {
            return if path.is_empty() {
                String::new()
            } else {
                "/".to_owned()
            };
        }
        trimmed.rsplit('/').next().unwrap_or(trimmed).to_owned()
    };
    match words.as_slice() {
        [one] => last(one),
        [one, suffix] => {
            let name = last(one);
            match name.strip_suffix(suffix.as_str()) {
                Some(stem) if !stem.is_empty() => stem.to_owned(),
                _ => name,
            }
        }
        _ => String::new(),
    }
}

/// Python's `int(TEXT)` printed back: white space around it, an optional
/// sign, digits with single underscores between them; None where int()
/// raises.
fn python_int(text: &str) -> Option<String> {
    let trimmed = text.trim();
    let (negative, digits) = match trimmed.as_bytes().first() {
        Some(b'-') => (true, &trimmed[1..]),
        Some(b'+') => (false, &trimmed[1..]),
        _ => (false, trimmed),
    };
    let valid = !digits.is_empty()
        && !digits.starts_with('_')
        && !digits.ends_with('_')
        && !digits.contains("__")
        && digits.chars().all(|c| c.is_ascii_digit() || c == '_');
    if !valid {
        return None;
    }
    let number = digits.replace('_', "");
    let number = number.trim_start_matches('0');
    Some(match (negative, number.is_empty()) {
        (_, true) => "0".to_owned(),
        (true, false) => format!("-{number}"),
        (false, false) => number.to_owned(),
    })
}

/// `get_pext_tissues_array`: the comma-separated names, each with its white
/// space trimmed and collapsed as `xargs` printed it, empty ones left out.
fn tissues(value: &str) -> Vec<String> {
    if value.is_empty() || value == "null" {
        return Vec::new();
    }
    value
        .split(',')
        .map(|tissue| tissue.split_whitespace().collect::<Vec<_>>().join(" "))
        .filter(|tissue| !tissue.is_empty())
        .collect()
}

/// What `anno_VEP_cache_worker` is given besides the configuration and the
/// input VCF: VEP's plain VCF output, its statistics file and its forks.
pub(crate) struct Worker<'a> {
    pub(crate) output: &'a str,
    pub(crate) stats: &'a str,
    pub(crate) forks: &'a str,
}

/// `anno_VEP_data --input_vcf VCF --config CONFIG --assembly ...
/// --threads ...`, every other setting from the configuration.
pub(crate) fn run(context: &Context, input: &str) -> Result<()> {
    annotate(context, input, None)
}

/// `anno_VEP_cache_worker --config CONFIG --input-vcf VCF --output-vcf OUTPUT
/// --stats-file STATS --forks N`: its checks of the paths and forks, then
/// VEP on VCF into OUTPUT, which it leaves unsorted.
pub(crate) fn worker(context: &Context, input: &str, given: &Worker) -> Result<()> {
    check_worker(context, input, given)?;
    annotate(context, input, Some(given))
}

/// The worker's own checks (its PY_VALIDATE): positive forks, within
/// SLURM_CPUS_PER_TASK with one more child and the parent when it is set;
/// absolute paths; a regular readable configuration and input; a fresh plain
/// .vcf output and another statistics file, in existing folders of this user
/// reached without links.
fn check_worker(context: &Context, input: &str, given: &Worker) -> Result<()> {
    let positive = |text: &str| {
        text.bytes()
            .next()
            .is_some_and(|b| (b'1'..=b'9').contains(&b))
            && text.bytes().all(|b| b.is_ascii_digit())
    };
    if !positive(given.forks) {
        return refuse("VEP worker forks must be a positive integer");
    }
    if let Some(cpus) = std::env::var_os("SLURM_CPUS_PER_TASK") {
        let cpus = cpus.to_string_lossy().into_owned();
        // Python's int(): any number of digits.
        let number = |text: &str| text.parse::<num_bigint::BigUint>().ok();
        let fits = positive(&cpus)
            && match (number(given.forks), number(&cpus)) {
                (Some(forks), Some(cpus)) => forks + 2u32 <= cpus,
                _ => false,
            };
        if !fits {
            return refuse(
                "VEP worker requires forks + one additional child + parent CPU within SLURM_CPUS_PER_TASK",
            );
        }
    }
    let config = context.config_text();
    let paths = [config.as_str(), input, given.output, given.stats];
    if paths.iter().any(|path| !Path::new(path).is_absolute()) {
        return refuse("VEP worker paths must be absolute");
    }
    for path in [config.as_str(), input] {
        let regular = std::fs::symlink_metadata(path)
            .map(|m| m.file_type().is_file())
            .map_err(|e| format!("{path}: {e}"))?;
        let readable = std::ffi::CString::new(path)
            .map(|c| unsafe { libc::access(c.as_ptr(), libc::R_OK) } == 0)
            .unwrap_or(false);
        if !regular || !readable {
            return refuse(format!(
                "VEP worker requires a readable regular input/config: {path}"
            ));
        }
    }
    let output = Path::new(given.output);
    if output.extension().is_none_or(|e| e != "vcf") || given.output == given.stats {
        return refuse("VEP worker needs a plain .vcf output and a distinct stats path");
    }
    let uid = unsafe { libc::getuid() };
    for path in [given.output, given.stats] {
        let path = Path::new(path);
        let parent = path.parent().unwrap_or(Path::new("/"));
        let owned = !path.components().any(|c| c == Component::ParentDir)
            && std::fs::canonicalize(parent).map_err(|e| format!("{}: {e}", parent.display()))?
                == parent
            && parent.is_dir()
            && std::fs::metadata(parent).is_ok_and(|m| m.uid() == uid);
        if !owned {
            return refuse(
                "VEP worker output parents must be existing owned directories without symlink aliases",
            );
        }
        if std::fs::symlink_metadata(path).is_ok() {
            return refuse(format!(
                "VEP worker output/stats must be fresh: {}",
                path.display()
            ));
        }
    }
    Ok(())
}

/// The body both share (`_priva_anno_VEP_data`).
fn annotate(context: &Context, input: &str, worker: Option<&Worker>) -> Result<()> {
    let v = |key: &str| context.value(key);
    let config = context.config_text();
    let mut assembly = v("assembly");
    let reference = v("ref_genome");
    let cache = v("vep_cache_dir");
    let plugins = v("vep_plugins_dir");
    let plugins_cache = v("vep_plugins_cachedir");
    // The wrapper always passed --threads, so argparse's int() read it; the
    // worker passed its forks.
    let configured = worker.map_or_else(|| v("threads"), |w| w.forks.to_owned());
    let Some(threads) = python_int(&configured) else {
        return refuse(format!(
            "argument --threads: invalid int value: '{configured}'"
        ));
    };
    let utr = v("utr_annotator_file");
    let loeuf = context.data(&format!("loeuf/{}", basename(&v("loeuf_prescore"))));
    let alphamissense = v("alphamissense_prescore");
    let spliceai_snv = v("spliceai_snv_prescore");
    let spliceai_indel = v("spliceai_indel_prescore");
    let primateai = v("primateai_prescore");
    let conservation = v("conservation_file");
    let splicevault = v("splicevault_prescore");
    let loftee = v("loftee_repo");
    let ancestor = v("human_ancestor_fasta");
    let loftee_conservation = v("loftee_conservation_file");
    let gerp = v("gerp_bigwig");
    let pext_enabled = match v("pext_enabled") {
        text if text.is_empty() => "false".to_owned(),
        text => text,
    };
    let pext_mean = v("pext_mean_bw");
    let pext_dir = v("pext_tissue_bw_dir");
    let pext_tissues = v("pext_tissues");

    // priva_check_spliceai_versions: both score files the same version (the
    // run of an analysis only).
    if worker.is_none()
        && !context.priva_succeeds(
            "spliceai-cutoffs",
            &[&config, &spliceai_snv, &spliceai_indel],
        )
    {
        return refuse("The SpliceAI score files cannot be used");
    }

    let mut errors = Vec::new();
    for (path, kind, name) in [
        (input, "file", "input_vcf"),
        (reference.as_str(), "file", "ref_genome"),
        (cache.as_str(), "dir", "vep_cache_dir"),
        (plugins.as_str(), "dir", "vep_plugins_dir"),
        (plugins_cache.as_str(), "dir", "vep_plugins_cachedir"),
        (utr.as_str(), "file", "utr_annotator_file"),
        (loeuf.as_str(), "file", "loeuf_prescore"),
        (alphamissense.as_str(), "file", "alphamissense_prescore"),
        (spliceai_snv.as_str(), "file", "spliceai_snv_prescore"),
        (spliceai_indel.as_str(), "file", "spliceai_indel_prescore"),
        (primateai.as_str(), "file", "primateai_prescore"),
        (conservation.as_str(), "file", "conservation_file"),
        (splicevault.as_str(), "file", "splicevault_prescore"),
        (loftee.as_str(), "dir", "loftee_repo"),
        (ancestor.as_str(), "file", "human_ancestor_fasta"),
        (
            loftee_conservation.as_str(),
            "file",
            "loftee_conservation_file",
        ),
    ] {
        errors.extend(quoted_path(path, kind, name));
    }
    if assembly.is_empty() {
        assembly = vcf_assembly(context, input);
    }
    // extract_assembly_from_fasta set the assembly only inside its own
    // subshell, so an assembly still unknown stayed empty.
    match assembly.as_str() {
        "hg19" => assembly = "GRCh37".into(),
        "hg38" => assembly = "GRCh38".into(),
        _ => {}
    }
    let (mut gerp_argument, mut conservation_argument, mut mavedb) =
        (String::new(), String::new(), Vec::new());
    if assembly == "GRCh37" {
        conservation_argument = format!(",conservation_file:{loftee_conservation}");
    } else if assembly == "GRCh38" {
        gerp_argument = format!(",gerp_bigwig:{gerp}");
        errors.extend(quoted_path(&gerp, "file", "gerp_bigwig"));
        let mavedb_file = v("mavedb_file");
        errors.extend(quoted_path(&mavedb_file, "file", "mavedb_file"));
        mavedb = words(&format!(
            "-plugin MaveDB,file={mavedb_file},cols={MAVEDB_COLUMNS}"
        ));
    }
    let mut pext = String::new();
    if pext_enabled == "true" {
        if !pext_mean.is_empty() && is_file(&pext_mean) {
            pext = format!("--custom {pext_mean},PEXT_MEAN,bigwig,overlap,0");
            log(&format!("Adding mean pext annotation from: {pext_mean}"));
        } else {
            log(&format!(
                "WARNING: pext_enabled=true but pext_mean_bw not found: {pext_mean}"
            ));
        }
        if !pext_tissues.is_empty() && pext_tissues != "null" && Path::new(&pext_dir).is_dir() {
            let folder = format!("{pext_dir}/tissues");
            if Path::new(&folder).is_dir() {
                for tissue in tissues(&pext_tissues) {
                    let track = format!("{folder}/{tissue}.bw");
                    if is_file(&track) {
                        let tag = format!("PEXT_{}", tissue.replace('-', "_"));
                        let custom = format!("--custom {track},{tag},bigwig,overlap,0");
                        pext = if pext.is_empty() {
                            custom
                        } else {
                            format!("{pext} {custom}")
                        };
                        log(&format!("Adding tissue pext annotation: {tissue} -> {tag}"));
                    } else {
                        log(&format!("WARNING: Tissue pext file not found: {track}"));
                    }
                }
            } else {
                log(&format!(
                    "WARNING: Tissue pext directory not found: {folder}"
                ));
            }
        }
    }
    if !errors.is_empty() {
        for error in &errors {
            log(&format!("Error: {error}"));
        }
        return refuse(format!("VEP cannot run: {}", errors.join("; ")));
    }
    if worker.is_none() && context.has_info_tags(input, "CSQ") && newer(input, context.binary) {
        log(&format!(
            "The input vcf {input} is already annotated by VEP. Skip this step"
        ));
        return Ok(());
    }

    let (output, stats) = match worker {
        Some(given) => (given.output.to_owned(), given.stats.to_owned()),
        None => (
            replace_vcf(input, &format!(".{}.vcf", random_tag())),
            replace_vcf(input, ".vep.stats.html"),
        ),
    };
    let mut command = args(&[
        "-i",
        input,
        "--format",
        "vcf",
        "--verbose",
        "--vcf",
        "--species",
        "homo_sapiens",
        "--use_transcript_ref",
    ]);
    command.push("--assembly".into());
    command.extend(words(&assembly));
    command.extend(args(&[
        "--cache",
        "--offline",
        "--merged",
        "--domains",
        "--hgvs",
        "--numbers",
        "--symbol",
        "--canonical",
        "--total_length",
        "--variant_class",
        "--gene_phenotype",
        "--stats_file",
        &stats,
        "--fork",
        &threads,
        "--buffer_size",
        "10000",
    ]));
    command.push("--fasta".into());
    command.extend(words(&reference));
    command.push("--dir_cache".into());
    command.extend(words(&cache));
    command.extend(["--dir_plugins".to_owned(), format!("{plugins} {loftee}")]);
    for plugin in [
        format!("UTRAnnotator,file={utr}"),
        format!("LOEUF,file={loeuf},match_by=transcript"),
        format!("AlphaMissense,file={alphamissense}"),
        format!(
            "LoF,loftee_path:{loftee},human_ancestor_fa:{ancestor}{conservation_argument}{gerp_argument}"
        ),
        format!("SpliceAI,snv={spliceai_snv},indel={spliceai_indel},cutoff=0.5"),
        format!("PrimateAI,{primateai}"),
        format!("SpliceVault,file={splicevault}"),
        format!("Conservation,{conservation},MAX"),
    ] {
        command.push("-plugin".into());
        command.extend(words(&plugin));
    }
    command.extend(args(&["-plugin", "NMD"]));
    command.extend(mavedb);
    command.extend(words(&pext));
    command.extend(args(&["--force_overwrite", "-o", &output]));
    if !context.succeeds("vep", &command) {
        return refuse(format!("VEP failed on {input}"));
    }
    if worker.is_some() {
        return Ok(());
    }
    let sorted = replace_vcf(input, ".vcf.gz");
    if !context.succeeds("bcftools", &args(&["sort", "-Oz", "-o", &sorted, &output]))
        || !context.succeeds("tabix", &args(&["-f", "-p", "vcf", &sorted]))
    {
        return refuse(format!("Failed to sort the VEP output {output}"));
    }
    let _ = std::fs::remove_file(&output);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basename_takes_the_words_as_the_shell_passed_them() {
        assert_eq!(basename("/a/b/loeuf.tsv.gz"), "loeuf.tsv.gz");
        assert_eq!(basename("null"), "null");
        assert_eq!(basename(""), "");
        assert_eq!(basename("/a/b.tsv .tsv"), "b");
        assert_eq!(basename("/a/b/"), "b");
    }

    #[test]
    fn threads_are_read_as_int_read_them() {
        for (text, number) in [
            ("8", Some("8")),
            (" +08 ", Some("8")),
            ("1_000", Some("1000")),
            ("-0", Some("0")),
            ("", None),
            ("null", None),
            ("1__0", None),
            ("4.0", None),
        ] {
            assert_eq!(python_int(text).as_deref(), number, "{text}");
        }
    }

    #[test]
    fn tissue_names_are_trimmed() {
        assert_eq!(
            tissues(" Brain_Cortex , Liver,,"),
            ["Brain_Cortex", "Liver"]
        );
        assert!(tissues("null").is_empty());
    }
}
