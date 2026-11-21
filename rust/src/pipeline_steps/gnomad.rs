//! The gnomAD step without `anno_agg_gnomAD_data` of annotation_vcf.sh: the
//! VCF split per chromosome (`bcftools view -r | bcftools sort`), each piece
//! annotated from that chromosome's gnomAD v4.1 joint VCF (`bcftools
//! annotate -a -c` over the piece's span), chrM kept without values, then
//! the pieces joined again (`bcftools concat | bcftools sort`).
//!
//! GNU parallel ran the pieces; here the same commands run on as many
//! threads as it had jobs at once: half the step's threads for the split
//! (at most SLURM_CPUS_PER_TASK), the step's threads for the annotation.
//! The split stops at the first failed piece, as `--halt soon,fail=1` did;
//! the annotation runs every piece, then fails if one did. Not written any
//! more: the job logs, which only GNU parallel needed. The split pieces are
//! removed with the annotated ones (the function left them in the folder).
//! The ClinVar recovery hook (`PRIVA_CLINVAR_RECOVERY_CONFIG`) belongs to the
//! installer's ClinVar chain, which still runs the shell function.
use super::tools::{
    Context, Result, args, is_file, log, not_empty, older, random_tag, refuse, replace_first,
    replace_vcf,
};
use std::collections::BTreeSet;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::Stdio;
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, Ordering};

const POPULATIONS: [&str; 9] = [
    "nfe",
    "eas",
    "afr",
    "amr",
    "asj",
    "fin",
    "sas",
    "mid",
    "remaining",
];

/// The INFO fields copied from gnomAD, as `-c` takes them.
fn columns() -> String {
    let mut fields = String::from(
        ".INFO/AC_joint,.INFO/AN_joint,.INFO/AF_joint,.INFO/nhomalt_joint\
         ,.INFO/AC_grpmax_joint,.INFO/AF_grpmax_joint,.INFO/AN_grpmax_joint,.INFO/nhomalt_grpmax_joint\
         ,.INFO/AC_joint_XX,.INFO/AN_joint_XX,.INFO/AC_joint_XY,.INFO/AN_joint_XY,.INFO/AF_joint_XX\
         ,.INFO/AF_joint_XY,.INFO/nhomalt_joint_XX,.INFO/nhomalt_joint_XY",
    );
    for population in POPULATIONS {
        fields.push_str(&format!(
            ",.INFO/AF_joint_{population}_XX,.INFO/AF_joint_{population}_XY"
        ));
    }
    for population in POPULATIONS {
        fields.push_str(&format!(",.INFO/nhomalt_joint_{population}_XX"));
    }
    fields
}

/// `[[ CHR =~ chr[0-9XY]+ ]]`: `chr` and a digit, X or Y anywhere in it.
fn supported(chromosome: &str) -> bool {
    chromosome.match_indices("chr").any(|(at, _)| {
        chromosome[at + 3..]
            .chars()
            .next()
            .is_some_and(|c| c.is_ascii_digit() || c == 'X' || c == 'Y')
    })
}

/// `${path/*.chrX.vcf/}`: what follows the last `.chrX.vcf`.
fn after_last_chrx(path: &str) -> String {
    match path.rfind(".chrX.vcf") {
        Some(at) => path[at + ".chrX.vcf".len()..].to_owned(),
        None => String::new(),
    }
}

/// `${path/.chrX.vcf*gz/}`: from the first `.chrX.vcf` to the last `gz`
/// taken out.
fn without_chrx(path: &str) -> String {
    let Some(at) = path.find(".chrX.vcf") else {
        return path.to_owned();
    };
    match path
        .rfind("gz")
        .filter(|end| *end >= at + ".chrX.vcf".len())
    {
        Some(end) => format!("{}{}", &path[..at], &path[end + 2..]),
        None => path.to_owned(),
    }
}

/// `mktemp -d DIR/.NAME.XXXXXX`.
fn stage(beside: &str, name: &str) -> Result<PathBuf> {
    let folder = Path::new(beside)
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let path = folder.join(format!(".{name}.{}", random_tag()));
    fs::create_dir(&path)?;
    Ok(path)
}

fn file_name(path: &str) -> String {
    Path::new(path)
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default()
}

/// A publication marker must be a regular file when it is there at all.
fn marker_ok(path: &str) -> bool {
    match fs::symlink_metadata(path) {
        Ok(meta) => meta.file_type().is_file(),
        Err(_) => true,
    }
}

/// `anno_agg_gnomAD_index --input-vcf VCF [--force]`: a `.tbi` not older
/// than the VCF, made through a link in a folder of its own and moved in.
fn index(context: &Context, vcf: &str, force: bool) -> Result<()> {
    let tbi = format!("{vcf}.tbi");
    if !not_empty(vcf) || Path::new(&tbi).is_dir() {
        return refuse(format!("Cannot index {vcf}"));
    }
    if !force && is_file(&tbi) && not_empty(&tbi) && !older(&tbi, vcf) {
        return Ok(());
    }
    let folder = stage(vcf, "gnomad-index")?;
    let canonical = fs::canonicalize(vcf)?;
    let alias = folder.join(file_name(vcf)).display().to_string();
    std::os::unix::fs::symlink(&canonical, &alias)?;
    if !context.succeeds("tabix", &args(&["-f", "-p", "vcf", &alias]))
        || !context
            .command("tabix")
            .args(["-l", &alias])
            .stdout(Stdio::null())
            .status()
            .is_ok_and(|s| s.success())
    {
        return refuse(format!("Failed to index {vcf}"));
    }
    fs::rename(format!("{alias}.tbi"), &tbi)?;
    fs::remove_file(&alias)?;
    fs::remove_dir(&folder)?;
    Ok(())
}

/// `mv STAGED OUT && mv -T STAGED.tbi OUT.tbi && rmdir STAGE`.
fn publish(staged: &str, output: &str, folder: &Path) -> Result<()> {
    fs::rename(staged, output)?;
    fs::rename(format!("{staged}.tbi"), format!("{output}.tbi"))?;
    fs::remove_dir(folder)?;
    Ok(())
}

/// `anno_agg_gnomAD_split --input-vcf VCF --chromosome CHR`: the piece
/// `${VCF/.vcf*/.CHR.vcf.gz}` with its index.
fn split(context: &Context, vcf: &str, chromosome: &str) -> Result<()> {
    if !is_file(vcf) || !supported(chromosome) {
        return refuse(format!("Cannot split {vcf} at {chromosome}"));
    }
    let piece = replace_vcf(vcf, &format!(".{chromosome}.vcf.gz"));
    let tbi = format!("{piece}.tbi");
    if Path::new(&tbi).is_dir() {
        return refuse(format!("{tbi} is a folder"));
    }
    if !is_file(&piece) || older(&piece, vcf) {
        let folder = stage(vcf, "gnomad-split")?;
        let staged = folder.join(file_name(&piece)).display().to_string();
        let sorted = context.pipeline(
            &[
                ("bcftools", args(&["view", "-r", chromosome, "-Ou", vcf])),
                ("bcftools", args(&["sort", "-Oz", "-o", &staged])),
            ],
            true,
        )?;
        if !sorted
            || !context.succeeds("tabix", &args(&["-f", "-p", "vcf", &staged]))
            || !context.vcf_valid(&staged, 1, None)
        {
            return refuse(format!("Failed to split {vcf} at {chromosome}"));
        }
        publish(&staged, &piece, &folder)?;
    }
    if (!is_file(&tbi) || older(&tbi, &piece))
        && !context.succeeds("tabix", &args(&["-f", "-p", "vcf", &piece]))
    {
        return refuse(format!("Failed to index {piece}"));
    }
    if !context.vcf_valid(&piece, 1, None) {
        return refuse(format!("{piece} is not valid"));
    }
    Ok(())
}

/// The `-r CONTIG:FIRST-LAST` span of a piece's records when they are on one
/// contig, as the function's awk read `bcftools query -f '%CHROM\t%POS\n'`.
fn span(context: &Context, vcf: &str) -> Result<String> {
    let mut child = context
        .command("bcftools")
        .args(["query", "-f", "%CHROM\\t%POS\\n", vcf])
        .stdout(Stdio::piped())
        .spawn()?;
    let reader = BufReader::new(child.stdout.take().expect("piped output"));
    let mut contig: Option<String> = None;
    let mut mixed = false;
    let mut first: Option<(f64, String)> = None;
    let mut last: Option<(f64, String)> = None;
    for line in reader.lines() {
        let line = line?;
        let (chrom, pos) = line.split_once('\t').unwrap_or((&line, ""));
        let number = pos.trim().parse::<f64>().unwrap_or(0.0);
        match &contig {
            None => {
                contig = Some(chrom.to_owned());
                first = Some((number, pos.to_owned()));
                last = Some((number, pos.to_owned()));
            }
            Some(name) => {
                if name != chrom {
                    mixed = true;
                }
                if first.as_ref().is_some_and(|(value, _)| number < *value) {
                    first = Some((number, pos.to_owned()));
                }
                if last.as_ref().is_some_and(|(value, _)| number > *value) {
                    last = Some((number, pos.to_owned()));
                }
            }
        }
    }
    if !child.wait()?.success() {
        return refuse(format!("bcftools query failed on {vcf}"));
    }
    Ok(match (contig, first, last) {
        (Some(contig), Some((_, first)), Some((_, last))) if !mixed => {
            format!("{contig}:{first}-{last}")
        }
        _ => String::new(),
    })
}

/// `anno_agg_gnomAD_chromosome --input-vcf PIECE --gnomad-vcf GNOMAD
/// --output-vcf OUT --columns FIELDS`.
fn annotate(
    context: &Context,
    piece: &str,
    gnomad: &str,
    output: &str,
    fields: &str,
) -> Result<()> {
    if output == piece || output == gnomad {
        return refuse(format!("{output} cannot be its own input"));
    }
    let pending = format!("{output}.gnomAD.publish.pending");
    if Path::new(&format!("{output}.tbi")).is_dir() || !marker_ok(&pending) {
        return refuse(format!("{output} cannot be written"));
    }
    if is_file(output) && Path::new(&pending).exists() {
        index(context, output, true)?;
        if !context.vcf_valid(output, 1, None) {
            return refuse(format!("{output} is not valid"));
        }
        fs::remove_file(&pending)?;
    }
    if is_file(output) && !older(output, piece) && !older(output, gnomad) {
        index(context, output, false)?;
        return if context.vcf_valid(output, 1, None) {
            Ok(())
        } else {
            refuse(format!("{output} is not valid"))
        };
    }
    let folder = stage(output, "gnomad-annotation")?;
    let staged = folder.join(file_name(output)).display().to_string();
    let region = span(context, piece)?;
    let mut command = args(&["annotate", "-a", gnomad, "-c", fields]);
    if !region.is_empty() {
        command.extend(args(&["-r", &region]));
    }
    command.extend(args(&["-Oz", "-o", &staged, piece]));
    if !context.succeeds("bcftools", &command)
        || !context.succeeds("tabix", &args(&["-f", "-p", "vcf", &staged]))
        || !context.vcf_valid(&staged, 1, None)
    {
        return refuse(format!("Failed to annotate {piece} with {gnomad}"));
    }
    let marker = folder.join("publication.pending");
    fs::write(&marker, format!("{output}\n"))?;
    fs::rename(&marker, &pending)?;
    fs::rename(&staged, output)?;
    fs::rename(format!("{staged}.tbi"), format!("{output}.tbi"))?;
    fs::remove_file(&pending)?;
    fs::remove_dir(&folder)?;
    Ok(())
}

/// `anno_agg_gnomAD_passthrough --input-vcf VCF --chromosome CHR
/// --output-vcf OUT`: a chromosome gnomAD has no data for, as it is.
fn passthrough(context: &Context, vcf: &str, chromosome: &str, output: &str) -> Result<()> {
    if Path::new(&format!("{output}.tbi")).is_dir() {
        return refuse(format!("{output}.tbi is a folder"));
    }
    let folder = stage(output, "gnomad-passthrough")?;
    let staged = folder.join(file_name(output)).display().to_string();
    if !context.succeeds(
        "bcftools",
        &args(&["view", "-r", chromosome, "-Oz", "-o", &staged, vcf]),
    ) || !context.succeeds("tabix", &args(&["-f", "-p", "vcf", &staged]))
        || !context.vcf_valid(&staged, 1, None)
    {
        return refuse(format!("Failed to keep {chromosome} of {vcf}"));
    }
    publish(&staged, output, &folder)
}

/// Runs JOB for each chromosome on WORKERS threads. With HALT no new job
/// starts after one failed; without it every job runs. The first failure.
fn each<F>(chromosomes: &[String], workers: usize, halt: bool, job: F) -> Result<()>
where
    F: Fn(&str) -> Result<()> + Sync,
{
    let next = Mutex::new(0usize);
    let failed = AtomicBool::new(false);
    let errors: Mutex<Vec<String>> = Mutex::new(Vec::new());
    std::thread::scope(|scope| {
        for _ in 0..workers.max(1).min(chromosomes.len().max(1)) {
            scope.spawn(|| {
                loop {
                    if halt && failed.load(Ordering::SeqCst) {
                        return;
                    }
                    let at = {
                        let mut next = next.lock().unwrap();
                        let at = *next;
                        *next += 1;
                        at
                    };
                    let Some(chromosome) = chromosomes.get(at) else {
                        return;
                    };
                    if let Err(error) = job(chromosome) {
                        failed.store(true, Ordering::SeqCst);
                        errors
                            .lock()
                            .unwrap()
                            .push(format!("{chromosome}: {error}"));
                    }
                }
            });
        }
    });
    let errors = errors.into_inner().unwrap();
    match errors.first() {
        Some(_) => refuse(errors.join("; ")),
        None => Ok(()),
    }
}

/// The VCF's chromosomes in file order, as `bcftools query -f '%CHROM\n' |
/// uniq` gave them.
fn chromosomes(context: &Context, vcf: &str) -> Result<Vec<String>> {
    let mut child = context
        .command("bcftools")
        .args(["query", "-f", "%CHROM\\n", vcf])
        .stdout(Stdio::piped())
        .spawn()?;
    let reader = BufReader::new(child.stdout.take().expect("piped output"));
    let mut out: Vec<String> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if out.last() != Some(&line) {
            out.push(line);
        }
    }
    if !child.wait()?.success() {
        return refuse(format!("bcftools query failed on {vcf}"));
    }
    Ok(out)
}

fn positive(text: &str) -> Option<usize> {
    let ok = !text.is_empty() && !text.starts_with('0') && text.chars().all(|c| c.is_ascii_digit());
    ok.then(|| text.parse().ok()).flatten()
}

/// `anno_agg_gnomAD_data VCF threads assembly gnomad_vcf_chrX`: the VCF
/// replaced by its annotated copy.
pub(crate) fn run(context: &Context, input: &str) -> Result<()> {
    let threads = context.value("threads");
    let gnomad_x = context.value("gnomad_vcf_chrX");
    let output = replace_first(input, ".vcf", &format!(".{}.vcf", random_tag()));
    if !input.ends_with(".vcf.gz") {
        return refuse(format!(
            "The input vcf {input} is not a gzipped vcf file. Quit with Error!"
        ));
    }
    let pending = format!("{input}.gnomAD.publish.pending");
    if !marker_ok(&pending) {
        return refuse(format!("{pending} is not a regular file"));
    }
    if Path::new(&pending).exists() {
        index(context, input, true)?;
        if !context.vcf_valid(input, 1, None) {
            return refuse(format!("{input} is not valid"));
        }
        fs::remove_file(&pending)?;
    } else {
        index(context, input, false)?;
    }
    let fields = columns();
    if !context.vcf_valid(input, 1, None) {
        return refuse(format!("The input vcf {input} is not valid. Quit now"));
    }
    if context.has_info_tags(input, &fields.replace(".INFO/", "")) {
        log(&format!(
            "The input vcf {input} already contains the gnomAD INFO tags. We do not need to add them again"
        ));
        return Ok(());
    }
    let listed = chromosomes(context, input)?;
    let (mut valid, mut kept) = (Vec::new(), Vec::new());
    for chromosome in listed {
        if supported(&chromosome) {
            valid.push(chromosome);
        } else if chromosome == "chrM" || chromosome == "chrMT" {
            log(&format!(
                "WARNING!!! gnomAD has no {chromosome} data. Its records are kept without gnomAD values"
            ));
            kept.push(chromosome);
        } else {
            log(&format!(
                "WARNING!!! The chromosome {chromosome} is not supported by gnomAD now. Currently only support chr1 to chrY. Skip this chromosome"
            ));
        }
    }
    if valid.is_empty() {
        return refuse("No supported chromosomes are available for gnomAD annotation.");
    }
    let Some(jobs) = positive(&threads) else {
        return refuse(format!("threads must be a positive number, not {threads}"));
    };
    let mut split_cpus = jobs;
    if let Ok(slurm) = std::env::var("SLURM_CPUS_PER_TASK") {
        let Some(slurm) = positive(&slurm) else {
            return refuse(format!(
                "SLURM_CPUS_PER_TASK must be a positive number, not {slurm}"
            ));
        };
        split_cpus = split_cpus.min(slurm);
    }
    each(&valid, (split_cpus / 2).max(1), true, |chromosome| {
        split(context, input, chromosome)
    })?;

    let stem = replace_vcf(input, "");
    let gnomad_stem = without_chrx(&gnomad_x);
    let suffix = after_last_chrx(&gnomad_x);
    each(&valid, jobs, false, |chromosome| {
        annotate(
            context,
            &format!("{stem}.{chromosome}.vcf.gz"),
            &format!("{gnomad_stem}.{chromosome}.vcf{suffix}"),
            &format!("{stem}.{chromosome}.gnomAD.vcf.gz"),
            &fields,
        )
    })?;
    for chromosome in &kept {
        passthrough(
            context,
            input,
            chromosome,
            &format!("{stem}.{chromosome}.gnomAD.vcf.gz"),
        )?;
    }
    // The shards, in the order the shell's glob listed them.
    let mut shards = BTreeSet::new();
    for chromosome in valid.iter().chain(&kept) {
        let shard = format!("{stem}.{chromosome}.gnomAD.vcf.gz");
        if !context.vcf_valid(&shard, 1, None) {
            return refuse(format!("{shard} is not valid"));
        }
        shards.insert(shard);
    }
    let mut concat = args(&["concat", "-Ou"]);
    concat.extend(shards.iter().cloned());
    let joined = context.pipeline(
        &[
            ("bcftools", concat),
            ("bcftools", args(&["sort", "-Oz", "-o", &output])),
        ],
        true,
    )?;
    if !joined
        || !context.succeeds("tabix", &args(&["-f", "-p", "vcf", &output]))
        || !context.vcf_valid(&output, 1, None)
    {
        return refuse(format!("Failed to join the gnomAD shards of {input}"));
    }
    fs::write(format!("{output}.publish.pending"), format!("{input}\n"))?;
    fs::rename(format!("{output}.publish.pending"), &pending)?;
    fs::rename(&output, input)?;
    fs::rename(format!("{output}.tbi"), format!("{input}.tbi"))?;
    fs::remove_file(&pending)?;
    for shard in &shards {
        let _ = fs::remove_file(shard);
        let _ = fs::remove_file(format!("{shard}.tbi"));
    }
    for chromosome in &valid {
        let piece = format!("{stem}.{chromosome}.vcf.gz");
        let _ = fs::remove_file(&piece);
        let _ = fs::remove_file(format!("{piece}.tbi"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gnomad_file_names_follow_the_shell_patterns() {
        let x = "/g/gnomad.joint.v4.1.sites.chrX.vcf.bgz";
        assert_eq!(without_chrx(x), "/g/gnomad.joint.v4.1.sites");
        assert_eq!(after_last_chrx(x), ".bgz");
        assert_eq!(without_chrx("/g/a.chrX.vcf.gz"), "/g/a");
        assert_eq!(after_last_chrx("/g/a.chrX.vcf.gz"), ".gz");
    }

    #[test]
    fn supported_chromosomes_match_the_unanchored_pattern() {
        assert!(supported("chr1"));
        assert!(supported("chrX"));
        assert!(supported("chr1_KI270706v1_random"));
        assert!(!supported("chrM"));
        assert!(!supported("chrUn_KI270742v1"));
        assert!(!supported("1"));
    }

    #[test]
    fn the_field_list_is_the_functions() {
        let fields = columns();
        assert!(fields.starts_with(".INFO/AC_joint,.INFO/AN_joint,"));
        assert!(
            fields.contains(",.INFO/AF_joint_nfe_XX,.INFO/AF_joint_nfe_XY,.INFO/AF_joint_eas_XX")
        );
        assert!(fields.ends_with(",.INFO/nhomalt_joint_mid_XX,.INFO/nhomalt_joint_remaining_XX"));
        assert_eq!(fields.split(',').count(), 16 + 18 + 9);
    }
}
