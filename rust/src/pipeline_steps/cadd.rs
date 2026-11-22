//! The CADD step without `Calculate_CADD` of annotation_vcf.sh: a Conda
//! frontend for CADD's Snakemake found as `cadd_conda.py check` found it,
//! GRC chromosome names (`bcftools annotate --rename-chrs | bcftools sort`),
//! the sites without MT records and genotypes (`bcftools view -G -t ^MT`),
//! then upstream `CADD.sh` with the function's options and environment,
//! `gunzip` and the table without its `##` lines.
use super::annotation::rename_chromosomes;
use super::tools::{Context, Result, args, is_file, log, option, refuse, replace_vcf, words};
use std::fs;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// A version as packaging's Version reads it ([v][N!]N(.N)*[pre][post][dev]
/// [+local], any case, pre, post and dev with or without a separator and
/// number), reduced to what a comparison with a plain release needs: the
/// epoch, the release without trailing zeros, and whether the version sorts
/// below that release (a pre-release, or a development release that is no
/// post-release). None when packaging would refuse it.
fn release(version: &str) -> Option<(u64, Vec<u64>, bool)> {
    fn digits(text: &str) -> usize {
        text.len() - text.trim_start_matches(|c: char| c.is_ascii_digit()).len()
    }
    /// A label from `labels` at the start of `text`, after an optional
    /// separator, then its optional number; what follows.
    fn segment<'a>(text: &'a str, labels: &[&str]) -> Option<&'a str> {
        let after = text.strip_prefix(['-', '_', '.']).unwrap_or(text);
        let label = labels.iter().find(|label| after.starts_with(*label))?;
        let after = &after[label.len()..];
        let number = after.strip_prefix(['-', '_', '.']).unwrap_or(after);
        Some(match digits(number) {
            0 => after,
            n => &number[n..],
        })
    }
    let lower = version.trim().to_ascii_lowercase();
    let text = lower.strip_prefix('v').unwrap_or(&lower);
    let text = match text.split_once('+') {
        Some((public, local)) => {
            let valid = local
                .split(['-', '_', '.'])
                .all(|part| !part.is_empty() && part.chars().all(|c| c.is_ascii_alphanumeric()));
            if !valid {
                return None;
            }
            public
        }
        None => text,
    };
    let (epoch, mut rest) = match text.split_once('!') {
        Some((epoch, rest)) if digits(epoch) == epoch.len() && !epoch.is_empty() => {
            (epoch.parse().ok()?, rest)
        }
        Some(_) => return None,
        None => (0, text),
    };
    let mut numbers = Vec::new();
    loop {
        let n = digits(rest);
        if n == 0 {
            return None;
        }
        numbers.push(rest[..n].parse().ok()?);
        rest = &rest[n..];
        match rest.strip_prefix('.') {
            Some(after) if digits(after) > 0 => rest = after,
            _ => break,
        }
    }
    let mut pre = false;
    if let Some(after) = segment(
        rest,
        &["alpha", "beta", "preview", "pre", "rc", "a", "b", "c"],
    ) {
        (pre, rest) = (true, after);
    }
    let mut post = false;
    if let Some(after) = segment(rest, &["post", "rev", "r"]) {
        (post, rest) = (true, after);
    } else if let Some(number) = rest.strip_prefix('-').filter(|n| digits(n) > 0) {
        (post, rest) = (true, &number[digits(number)..]);
    }
    let mut dev = false;
    if let Some(after) = segment(rest, &["dev"]) {
        (dev, rest) = (true, after);
    }
    if !rest.is_empty() {
        return None;
    }
    while numbers.len() > 1 && numbers.last() == Some(&0) {
        numbers.pop();
    }
    Some((epoch, numbers, pre || (dev && !post)))
}

/// Snakemake 8.25.2 needs Conda 24.7.1 or newer.
fn new_enough(version: &str) -> bool {
    let Some((epoch, numbers, below)) = release(version) else {
        return false;
    };
    let minimum = vec![24, 7, 1];
    epoch > 0 || numbers > minimum || (numbers == minimum && !below)
}

fn json_of(program: &Path, args: &[&str]) -> std::result::Result<crate::json_text::Json, String> {
    let out = Command::new(program)
        .args(args)
        .stderr(Stdio::null())
        .output()
        .map_err(|e| e.to_string())?;
    if !out.status.success() {
        return Err(format!("{} {} failed", program.display(), args.join(" ")));
    }
    let text = String::from_utf8_lossy(&out.stdout);
    let value = crate::json_text::parse(&text, false).map_err(|e| e.to_string())?;
    match value {
        crate::json_text::Json::Object(_) => Ok(value),
        _ => Err(format!(
            "expected an object from {} {}",
            program.display(),
            args.join(" ")
        )),
    }
}

/// `inspect_conda`: version, platform, a prefix with bin/activate and a
/// readable channel priority. (executable, prefix, version).
fn inspect(executable: &Path) -> std::result::Result<(PathBuf, String, String), String> {
    let executable = std::path::absolute(executable).map_err(|e| e.to_string())?;
    let out = Command::new(&executable)
        .arg("--version")
        .stderr(Stdio::null())
        .output()
        .map_err(|e| e.to_string())?;
    if !out.status.success() {
        return Err(format!("{} --version failed", executable.display()));
    }
    let version = String::from_utf8_lossy(&out.stdout).trim().to_owned();
    let number = version
        .strip_prefix("conda ")
        .filter(|n| !n.contains(char::is_whitespace));
    if !number.is_some_and(new_enough) {
        return Err(format!(
            "{}: CADD requires Conda >= 24.7.1; found '{version}'",
            executable.display()
        ));
    }
    let info = json_of(&executable, &["info", "--json"])?;
    if !info
        .get("platform")
        .and_then(crate::json_text::Json::as_str)
        .is_some_and(|p| !p.is_empty())
    {
        return Err(format!(
            "{}: Conda info lacks its platform",
            executable.display()
        ));
    }
    let prefix = info
        .get("conda_prefix")
        .and_then(crate::json_text::Json::as_str)
        .unwrap_or("")
        .to_owned();
    if !Path::new(&prefix).is_absolute() || !Path::new(&prefix).join("bin/activate").is_file() {
        return Err(format!(
            "{}: reported Conda prefix lacks bin/activate",
            executable.display()
        ));
    }
    let config = json_of(
        &executable,
        &["config", "--get", "channel_priority", "--json"],
    )?;
    let failed = matches!(
        config.get("success"),
        Some(crate::json_text::Json::Bool(false))
    );
    if failed || !matches!(config.get("get"), Some(crate::json_text::Json::Object(_))) {
        return Err(format!(
            "{}: cannot read Conda channel priority",
            executable.display()
        ));
    }
    Ok((executable, prefix, version))
}

/// `shutil.which(NAME, path=PATH)`.
fn which(name: &str, path: &str) -> Option<PathBuf> {
    use std::os::unix::fs::PermissionsExt;
    path.split(':').find_map(|folder| {
        let folder = if folder.is_empty() { "." } else { folder };
        let candidate = Path::new(folder).join(name);
        let executable = fs::metadata(&candidate)
            .is_ok_and(|m| m.is_file() && m.permissions().mode() & 0o111 != 0);
        executable.then_some(candidate)
    })
}

fn same_file(left: Option<&Path>, right: &Path) -> bool {
    use std::os::unix::fs::MetadataExt;
    let Some(left) = left else {
        return false;
    };
    match (fs::metadata(left), fs::metadata(right)) {
        (Ok(a), Ok(b)) => a.dev() == b.dev() && a.ino() == b.ino(),
        _ => false,
    }
}

/// `child_path`: PATH with the selected Conda reachable and the scientific
/// Python unchanged.
fn child_path(executable: &Path, prefix: &str) -> std::result::Result<String, String> {
    let original = std::env::var("PATH").unwrap_or_default();
    if same_file(which("conda", &original).as_deref(), executable) {
        return Ok(original);
    }
    let folder = executable
        .parent()
        .unwrap_or(Path::new("/"))
        .display()
        .to_string();
    let mut path = format!("{original}:{folder}");
    if !same_file(which("conda", &path).as_deref(), executable) {
        let condabin = Path::new(prefix).join("condabin");
        let (_, alternate_prefix, _) = inspect(&condabin.join("conda"))?;
        if alternate_prefix != prefix || condabin.join("python").exists() {
            return Err(
                "cannot expose the selected Conda without changing scientific Python".into(),
            );
        }
        path = format!("{}:{original}", condabin.display());
    }
    if which("python", &path) != which("python", &original) {
        return Err("the Conda frontend must not replace the active scientific Python".into());
    }
    Ok(path)
}

/// `cadd_conda.py check --frontend-prefix PREFIX --print-path`: the PATH
/// CADD.sh runs with.
fn conda_path(frontend: &Path) -> Result<String> {
    if !frontend.is_absolute() || frontend.file_name().is_none_or(|n| n != "conda_frontend") {
        return refuse("ERROR: use the absolute cadd_parent_dir/conda_frontend prefix");
    }
    let original = std::env::var("PATH").unwrap_or_default();
    let candidates = [
        which("conda", &original).map(|p| p.display().to_string()),
        std::env::var("CONDA_EXE").ok(),
        Some(frontend.join("bin/conda").display().to_string()),
    ];
    let mut seen = Vec::new();
    let mut failures = Vec::new();
    for candidate in candidates.into_iter().flatten() {
        if candidate.is_empty() || seen.contains(&candidate) {
            continue;
        }
        seen.push(candidate.clone());
        if !Path::new(&candidate).exists() {
            continue;
        }
        match inspect(Path::new(&candidate))
            .and_then(|(exe, prefix, version)| Ok((child_path(&exe, &prefix)?, exe, version)))
        {
            Ok((path, executable, version)) => {
                log(&format!("Using {version}: {}", executable.display()));
                return Ok(path);
            }
            Err(error) => failures.push(format!("{candidate}: {error}")),
        }
    }
    if !failures.is_empty() {
        log(&format!(
            "Unusable Conda candidates: {}",
            failures.join("; ")
        ));
    }
    refuse(
        "ERROR: CADD.sh needs a conda executable (>= 24.7.1); install step E1 (`nextflow run main.nf --action install`) sets one up",
    )
}

/// `${path%/*}`.
fn parent_text(path: &str) -> String {
    match path.rfind('/') {
        Some(at) => path[..at].to_owned(),
        None => path.to_owned(),
    }
}

/// `Calculate_CADD VCF CONFIG OUTPUT`: OUTPUT is the CADD table.
pub(crate) fn run(context: &Context, input: &str, output: &str) -> Result<()> {
    let threads = context.value("threads");
    let base = context.value("cadd_base_dir");
    let assembly = context.value("assembly");
    let temp = context.value("tmp_dir");
    let script = format!("{base}/CADD.sh");
    let runnable = {
        use std::os::unix::fs::PermissionsExt;
        fs::metadata(&script).is_ok_and(|m| m.permissions().mode() & 0o111 != 0)
    };
    if !runnable {
        return refuse(format!(
            "CADD.sh is missing from cadd_base_dir {base}; install step E1 (`nextflow run main.nf --action install`) installs it"
        ));
    }
    let path = conda_path(Path::new(&format!("{}/conda_frontend", parent_text(&base))))?;
    let first = context.first_line("bcftools", &["view", "-H", input]);
    if is_file(input) && first.is_empty() {
        let table = output.strip_suffix(".gz").unwrap_or(output);
        fs::write(table, "")?;
        log(&format!(
            "{input} has no records; wrote an empty CADD table {table}"
        ));
        return Ok(());
    }
    if !context.vcf_valid(input, 1, None) {
        return refuse(format!(
            "The input vcf file {input} is not valid. Please check the file"
        ));
    }
    let compressed = if output.ends_with(".gz") {
        output.to_owned()
    } else {
        format!("{output}.gz")
    };
    let table = compressed.replacen(".gz", "", 1);
    let mut grc = replace_vcf(input, ".nochr.vcf.gz");
    if grc == input {
        grc = match input.find(".bcf") {
            Some(at) => format!("{}.nochr.vcf.gz", &input[..at]),
            None => input.to_owned(),
        };
    }
    if !rename_chromosomes(
        context,
        input,
        &context.data("liftover/ucsc_to_GRC.contig.map.tsv"),
        &grc,
    )? {
        return refuse(format!(
            "Failed to convert chromosome names from UCSC to GRCh style for {input}. Quit now"
        ));
    }
    let sites = format!(
        "{}.cadd_input.vcf.gz",
        grc.strip_suffix(".vcf.gz").unwrap_or(&grc)
    );
    let mut view = vec!["view".to_owned()];
    view.extend(option("--threads", &threads));
    view.extend(args(&["-G", "-t", "^MT", "-Oz", "-o", &sites, &grc]));
    if !context.succeeds("bcftools", &view) {
        return refuse(format!(
            "Failed to write the CADD input {sites} without MT records and genotypes. Quit now"
        ));
    }
    let genome = match assembly.as_str() {
        "GRCh37" | "hg19" => "GRCh37",
        "GRCh38" | "hg38" => "GRCh38",
        _ => {
            return refuse(format!(
                "The assembly {assembly} is not supported. Please check the config file {}",
                context.config_text()
            ));
        }
    };
    let exon_tables = format!("{base}/data/annotations/{genome}_v1.7/mmsplice");
    let exon_tables = if Path::new(&exon_tables).is_dir() {
        exon_tables
    } else {
        String::new()
    };
    let mut command = vec!["-c".to_owned()];
    command.extend(words(&threads));
    command.extend(args(&[
        "-a",
        "-p",
        "-m",
        "-d",
        "-g",
        genome,
        "-o",
        &compressed,
        &sites,
    ]));
    log(&format!("Running: {script} {}", command.join(" ")));
    let mut cadd = context.command(&script);
    cadd.args(&command)
        .env("PATH", &path)
        .env("TMPDIR", &temp)
        .env("PRIVA_CADD_MMSPLICE_EXON_CACHE", &exon_tables)
        .env("PRIVA_CADD_PRESCORE_BIN", context.binary);
    let ran = cadd.status().is_ok_and(|s| s.success())
        && context.succeeds("gunzip", &args(&["-f", &compressed]));
    if !ran {
        return refuse(format!("Failed to run CADD on {input}. Quit now"));
    }
    // mawk -F '\t' '$1 !~ /^##/'
    let staged = format!("{table}.tmp");
    {
        let reader = BufReader::new(fs::File::open(&table)?);
        let mut writer = BufWriter::new(fs::File::create(&staged)?);
        for line in reader.split(b'\n') {
            let line = line?;
            if !line.starts_with(b"##") {
                writer.write_all(&line)?;
                writer.write_all(b"\n")?;
            }
        }
        writer.flush()?;
    }
    fs::rename(&staged, &table)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conda_versions_compare_as_packaging_did() {
        assert!(new_enough("24.7.1"));
        assert!(new_enough("24.9.2"));
        assert!(new_enough("25.1"));
        assert!(!new_enough("24.7.0"));
        assert!(!new_enough("24.7.1rc1"));
        assert!(new_enough("24.7.1.post1"));
        assert!(new_enough("24.7.1-1"));
        assert!(new_enough("24.7.1.post1.dev0"));
        assert!(!new_enough("24.7.1.dev3"));
        assert!(new_enough("V24.7.1+local.7"));
        assert!(new_enough("1!1.0"));
        assert!(new_enough("24.7.1.0"));
        assert!(!new_enough("24.7.1+"));
        assert!(!new_enough("24.7.1x"));
        assert!(!new_enough("x"));
    }

    #[test]
    fn the_parent_is_taken_as_the_shell_took_it() {
        assert_eq!(parent_text("/c/CADD-scripts"), "/c");
        assert_eq!(parent_text("CADD"), "CADD");
    }
}
