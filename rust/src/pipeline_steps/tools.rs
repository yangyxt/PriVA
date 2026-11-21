//! What every pipeline step shares: its settings, the environment its tools
//! run in, and the checks the shell scripts made with
//! `common_bash_utils.sh` (`read_yaml`, `check_vcf_validity`,
//! `check_path`, `check_vcf_infotags`, `check_table_column`).
//!
//! The steps' hub lineages count this file (rust_reach.rs from each step's
//! own file), so it names nothing of the task contract.
use crate::yaml_config::{ConfigDocument, ConfigValue, Parsed, Scalar};
use std::collections::BTreeSet;
use std::error::Error;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

pub(crate) type Result<T> = std::result::Result<T, Box<dyn Error>>;

pub(crate) fn refuse<T>(message: impl Into<String>) -> Result<T> {
    Err(message.into().into())
}

pub(crate) fn log(message: &str) {
    eprintln!("[priva task] {message}");
}

/// A step's configuration: the task's private copy.
pub(crate) struct Settings {
    doc: ConfigDocument,
}

impl Settings {
    pub(crate) fn read(path: &Path) -> Result<Settings> {
        let text = fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
        match crate::yaml_config::parse_config(&text) {
            Ok(Parsed::Mapping(doc)) => Ok(Settings { doc }),
            Ok(_) => refuse(format!("{} is not a YAML mapping", path.display())),
            Err(error) => Err(error.to_string().into()),
        }
    }

    #[cfg(test)]
    pub(crate) fn parse(text: &str) -> Settings {
        match crate::yaml_config::parse_config(text).unwrap() {
            Parsed::Mapping(doc) => Settings { doc },
            _ => panic!("not a mapping"),
        }
    }

    /// `read_yaml FILE KEY` of common_bash_utils.sh: yq's JSON of the value
    /// with the quotes of a string taken off. A missing key or null is
    /// `null`, a boolean `true` or `false`.
    pub(crate) fn value(&self, key: &str) -> String {
        match self.doc.get(key) {
            None | Some(ConfigValue::Scalar(Scalar::Null)) => "null".to_owned(),
            Some(ConfigValue::Scalar(Scalar::Bool(flag))) => flag.to_string(),
            Some(ConfigValue::Scalar(Scalar::Str(text))) => text.clone(),
            Some(ConfigValue::Scalar(Scalar::Int(text))) => text.clone(),
            Some(ConfigValue::Scalar(other)) => other.python_repr(),
            Some(ConfigValue::Opaque(text)) => text.clone(),
        }
    }
}

/// What a step needs: the source checkout (its `data/` folder), the pinned
/// priva program, the private configuration and the environment the task
/// gives its commands.
pub(crate) struct Context<'a> {
    pub(crate) root: &'a Path,
    pub(crate) binary: &'a str,
    pub(crate) config_path: &'a Path,
    pub(crate) settings: &'a Settings,
    pub(crate) env: &'a [(String, String)],
}

impl Context<'_> {
    /// `read_yaml CONFIG KEY`.
    pub(crate) fn value(&self, key: &str) -> String {
        self.settings.value(key)
    }

    /// A file shipped in the checkout's `data/` folder.
    pub(crate) fn data(&self, relative: &str) -> String {
        self.root.join("data").join(relative).display().to_string()
    }

    pub(crate) fn config_text(&self) -> String {
        self.config_path.display().to_string()
    }

    /// The temporary folder of the task (TMPDIR), else the system's.
    pub(crate) fn temp_dir(&self) -> PathBuf {
        self.env
            .iter()
            .find(|(k, _)| k == "TMPDIR")
            .map(|(_, v)| PathBuf::from(v))
            .unwrap_or_else(std::env::temp_dir)
    }

    pub(crate) fn command(&self, program: &str) -> Command {
        let mut command = Command::new(program);
        for (key, value) in self.env {
            command.env(key, value);
        }
        // priva_resolve_binary pointed the program at this checkout.
        command.env("PRIVA_DATA_DIR", self.root.join("data"));
        command.env("PRIVA_SCRIPT_DIR", self.root.join("scripts"));
        command
    }

    /// Runs PROGRAM ARGS and fails with its exit status.
    pub(crate) fn run(&self, program: &str, args: &[String]) -> Result<()> {
        log(&format!("Running: {program} {}", args.join(" ")));
        let status = self.command(program).args(args).status()?;
        if !status.success() {
            return refuse(format!(
                "{program} {} failed ({status})",
                args.first().map(String::as_str).unwrap_or("")
            ));
        }
        Ok(())
    }

    /// `"${PRIVA_BIN}" COMMAND ARGS`.
    pub(crate) fn priva(&self, command: &str, args: Vec<String>) -> Result<()> {
        let mut all = vec![command.to_owned()];
        all.extend(args);
        self.run(self.binary, &all)
    }

    /// `"${PRIVA_BIN}" COMMAND ARGS >&2`: whether it succeeded, its output
    /// sent to the log.
    pub(crate) fn priva_succeeds(&self, command: &str, args: &[&str]) -> bool {
        log(&format!("Running: priva {command} {}", args.join(" ")));
        self.command(self.binary)
            .arg(command)
            .args(args)
            .stdout(Stdio::from(std::io::stderr()))
            .status()
            .is_ok_and(|status| status.success())
    }

    /// Runs PROGRAM ARGS; whether it succeeded.
    pub(crate) fn succeeds(&self, program: &str, args: &[String]) -> bool {
        log(&format!("Running: {program} {}", args.join(" ")));
        self.command(program)
            .args(args)
            .status()
            .is_ok_and(|status| status.success())
    }

    /// `A | B | ...`: the last command's status, as the scripts had it
    /// without pipefail, or every command's with `pipefail`.
    pub(crate) fn pipeline(&self, stages: &[(&str, Vec<String>)], pipefail: bool) -> Result<bool> {
        let text: Vec<String> = stages
            .iter()
            .map(|(program, args)| format!("{program} {}", args.join(" ")))
            .collect();
        log(&format!("Running: {}", text.join(" | ")));
        let mut children = Vec::new();
        let mut previous: Option<std::process::ChildStdout> = None;
        for (index, (program, args)) in stages.iter().enumerate() {
            let mut command = self.command(program);
            command.args(args);
            if let Some(input) = previous.take() {
                command.stdin(input);
            }
            if index + 1 < stages.len() {
                command.stdout(Stdio::piped());
            }
            let mut child = command.spawn()?;
            previous = child.stdout.take();
            children.push(child);
        }
        let mut statuses = Vec::new();
        for mut child in children {
            statuses.push(child.wait()?.success());
        }
        Ok(if pipefail {
            statuses.iter().all(|ok| *ok)
        } else {
            *statuses.last().unwrap_or(&false)
        })
    }

    /// The first line a command prints (`COMMAND | head -1`), without its
    /// newline; the command is stopped after it.
    pub(crate) fn first_line(&self, program: &str, args: &[&str]) -> String {
        let Ok(mut child) = self
            .command(program)
            .args(args)
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
        else {
            return String::new();
        };
        let mut line = String::new();
        if let Some(out) = child.stdout.take() {
            let _ = BufReader::new(out).read_line(&mut line);
        }
        let _ = child.kill();
        let _ = child.wait();
        line.trim_end_matches('\n').to_owned()
    }

    /// The output of a command, empty when it cannot run.
    pub(crate) fn output(&self, program: &str, args: &[&str]) -> String {
        self.command(program)
            .args(args)
            .stderr(Stdio::null())
            .output()
            .map(|out| String::from_utf8_lossy(&out.stdout).into_owned())
            .unwrap_or_default()
    }

    /// `count_vcf_records VCF`: `bcftools query -f '\n' VCF | wc -l`.
    pub(crate) fn count_records(&self, vcf: &str) -> usize {
        self.output("bcftools", &["query", "-f", "\\n", vcf])
            .bytes()
            .filter(|b| *b == b'\n')
            .count()
    }

    /// `bcftools view -h VCF`.
    pub(crate) fn header(&self, vcf: &str) -> String {
        self.output("bcftools", &["view", "-h", vcf])
    }

    /// `check_vcf_index VCF`: a current `.tbi` (or `.csi`), made or renewed
    /// with bcftools when it is missing or older than the VCF.
    pub(crate) fn vcf_index(&self, vcf: &str) -> bool {
        let ok = |args: &[&str]| {
            self.command("bcftools")
                .args(args)
                .status()
                .is_ok_and(|s| s.success())
        };
        if vcf.ends_with(".vcf.gz") || vcf.ends_with(".vcf.bgz") {
            let tbi = format!("{vcf}.tbi");
            let csi = format!("{vcf}.csi");
            if Path::new(&tbi).is_file() {
                if newer(&tbi, vcf) {
                    let _ = fs::remove_file(&csi);
                    return true;
                }
                let renewed = ok(&["index", "-t", "-f", vcf]);
                if renewed {
                    let _ = fs::remove_file(&csi);
                }
                renewed
            } else if Path::new(&csi).is_file() {
                !older(&csi, vcf) || ok(&["index", "-f", vcf])
            } else {
                ok(&["index", "-t", "-f", vcf])
            }
        } else if vcf.ends_with(".vcf") {
            true
        } else if vcf.ends_with(".bcf") {
            let csi = format!("{vcf}.csi");
            (Path::new(&csi).is_file() && newer(&csi, vcf)) || ok(&["index", "-f", vcf])
        } else {
            false
        }
    }

    /// `check_vcf_validity VCF [LINES [SAMPLES]]`: a compressed, plain or
    /// BCF file with its index, at least LINES records (1 when not given)
    /// and, when SAMPLES (comma-separated) are given, exactly those samples.
    /// The script's `bcftools view | head` format check could not fail
    /// without pipefail, so it is not repeated. Not repeated either: its
    /// fallback to a plain copy beside a compressed file (a task's folder
    /// never holds one), and the compressed copy it wrote beside a plain
    /// file to list its samples.
    pub(crate) fn vcf_valid(&self, vcf: &str, lines: usize, samples: Option<&str>) -> bool {
        if !Path::new(vcf).is_file() {
            return false;
        }
        let known = vcf.ends_with("gz") || vcf.ends_with(".vcf") || vcf.ends_with(".bcf");
        if !known || !self.vcf_index(vcf) || self.count_records(vcf) < lines {
            return false;
        }
        match samples {
            None | Some("") => true,
            Some(expected) => {
                let sorted = |items: Vec<String>| {
                    let set: BTreeSet<String> =
                        items.into_iter().filter(|s| !s.is_empty()).collect();
                    set.into_iter().collect::<Vec<_>>().join(",")
                };
                let actual = sorted(
                    self.output("bcftools", &["query", "-l", vcf])
                        .lines()
                        .map(str::to_owned)
                        .collect(),
                );
                actual == sorted(expected.split(',').map(str::to_owned).collect())
            }
        }
    }

    /// `check_vcf_infotags VCF TAG,TAG,...`: every tag has an INFO header
    /// line (`##INFO=<ID=TAG,`).
    pub(crate) fn has_info_tags(&self, vcf: &str, tags: &str) -> bool {
        if !is_file(vcf) {
            return false;
        }
        let wanted: Vec<&str> = tags.split(',').filter(|t| !t.is_empty()).collect();
        if wanted.is_empty() {
            return false;
        }
        let defined: BTreeSet<String> = self
            .header(vcf)
            .lines()
            .filter_map(|line| line.strip_prefix("##INFO=<ID="))
            .map(|rest| rest.split(',').next().unwrap_or("").to_owned())
            .collect();
        wanted.iter().all(|tag| defined.contains(*tag))
    }
}

/// A value the scripts expanded without quotes: its words.
pub(crate) fn words(value: &str) -> Vec<String> {
    value.split_whitespace().map(str::to_owned).collect()
}

/// `--OPTION VALUE` with VALUE unquoted.
pub(crate) fn option(name: &str, value: &str) -> Vec<String> {
    let mut out = vec![name.to_owned()];
    out.extend(words(value));
    out
}

/// Owned arguments.
pub(crate) fn args(items: &[&str]) -> Vec<String> {
    items.iter().map(|item| (*item).to_owned()).collect()
}

/// `${path/.vcf*/SUFFIX}`: everything from the first `.vcf` replaced.
pub(crate) fn replace_vcf(path: &str, suffix: &str) -> String {
    match path.find(".vcf") {
        Some(at) => format!("{}{suffix}", &path[..at]),
        None => path.to_owned(),
    }
}

/// `${path/.vcf/SUFFIX}`: the first `.vcf` replaced.
pub(crate) fn replace_first(path: &str, pattern: &str, with: &str) -> String {
    path.replacen(pattern, with, 1)
}

fn mtime(path: &str) -> Option<std::time::SystemTime> {
    fs::metadata(path).and_then(|m| m.modified()).ok()
}

/// `[[ A -nt B ]]`.
pub(crate) fn newer(a: &str, b: &str) -> bool {
    match (mtime(a), mtime(b)) {
        (Some(x), Some(y)) => x > y,
        (Some(_), None) => true,
        _ => false,
    }
}

/// `[[ A -ot B ]]`.
pub(crate) fn older(a: &str, b: &str) -> bool {
    match (mtime(a), mtime(b)) {
        (Some(x), Some(y)) => x < y,
        (None, Some(_)) => true,
        _ => false,
    }
}

pub(crate) fn is_file(path: &str) -> bool {
    Path::new(path).is_file()
}

/// `[[ -s PATH ]]`.
pub(crate) fn not_empty(path: &str) -> bool {
    fs::metadata(path).is_ok_and(|m| m.len() > 0)
}

/// `check_table_column TABLE COLUMN`: the header names COLUMN.
pub(crate) fn has_column(table: &str, column: &str) -> bool {
    let Ok(file) = fs::File::open(table) else {
        return false;
    };
    let mut header = String::new();
    if BufReader::new(file).read_line(&mut header).is_err() {
        return false;
    }
    header
        .trim_end_matches(['\n', '\r'])
        .split('\t')
        .any(|name| name == column)
}

pub(crate) fn line_count(path: &str) -> usize {
    fs::read(path)
        .map(|bytes| bytes.iter().filter(|b| **b == b'\n').count())
        .unwrap_or(0)
}

/// `grep PATTERN` for a pattern whose only special character is `.` (any
/// character); a pattern with other special characters is looked for as it
/// is.
pub(crate) fn grep_matches(line: &str, pattern: &str) -> bool {
    if pattern.contains(['*', '[', '^', '$', '\\']) {
        return line.contains(pattern);
    }
    let needle: Vec<char> = pattern.chars().collect();
    let hay: Vec<char> = line.chars().collect();
    hay.len() >= needle.len()
        && (0..=hay.len() - needle.len()).any(|start| {
            needle
                .iter()
                .zip(&hay[start..])
                .all(|(p, h)| *p == '.' || p == h)
        })
}

/// `check_path PATH TYPE NAME` with PATH unquoted: the shell passed PATH's
/// words and then `TYPE NAME`, so an empty value made TYPE the path and NAME
/// the type (which passes), and a value of two words made its second word
/// the type. The error message, when it fails.
pub(crate) fn unquoted_path(value: &str, kind: &str, name: &str) -> Option<String> {
    let mut arguments = words(value);
    arguments.extend([kind.to_owned(), name.to_owned()]);
    let label = arguments.get(2).map_or("", String::as_str);
    quoted_path(&arguments[0], &arguments[1], label)
}

/// `check_path "PATH" TYPE NAME`.
pub(crate) fn quoted_path(path: &str, kind: &str, name: &str) -> Option<String> {
    if path.is_empty() {
        return Some(format!("{name} not specified"));
    }
    match kind {
        "file" if !is_file(path) => Some(format!("{name} file not found: {path}")),
        "dir" if !Path::new(path).is_dir() => Some(format!("{name} directory not found: {path}")),
        _ => None,
    }
}

/// A name for a temporary file, as `randomID` gave.
pub(crate) fn random_tag() -> String {
    use std::sync::atomic::{AtomicU64, Ordering};
    use std::time::{SystemTime, UNIX_EPOCH};
    static COUNT: AtomicU64 = AtomicU64::new(0);
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos());
    let count = COUNT.fetch_add(1, Ordering::Relaxed);
    format!("{:x}{:x}{count:x}", std::process::id(), nanos)
}

/// `mv A B` for a file and its `.tbi`.
pub(crate) fn move_with_index(from: &str, to: &str) -> Result<()> {
    fs::rename(from, to)?;
    fs::rename(format!("{from}.tbi"), format!("{to}.tbi"))?;
    Ok(())
}

/// Removes files, as `rm -f`.
pub(crate) fn remove(paths: &[&str]) {
    for path in paths {
        let _ = fs::remove_file(path);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vcf_suffix_replacement_follows_the_shell_pattern() {
        assert_eq!(replace_vcf("/r/a.anno.vcf.gz", ".tsv"), "/r/a.anno.tsv");
        assert_eq!(
            replace_vcf("a.anno.F1.vcf.gz", ".filtered.vcf.gz"),
            "a.anno.F1.filtered.vcf.gz"
        );
        assert_eq!(replace_vcf("table.tsv", ".x"), "table.tsv");
        assert_eq!(replace_first("a.vcf.gz", ".vcf", ".x.vcf"), "a.x.vcf.gz");
    }

    #[test]
    fn unquoted_values_split_into_words() {
        assert_eq!(option("--x", ""), vec!["--x"]);
        assert_eq!(option("--x", "a, b"), vec!["--x", "a,", "b"]);
        assert_eq!(unquoted_path("", "file", "k"), None);
        assert_eq!(
            quoted_path("", "file", "k"),
            Some("k not specified".to_owned())
        );
    }

    #[test]
    fn unquoted_paths_are_checked_as_the_shell_passed_their_words() {
        // Two words: the second became check_path's type, so nothing failed.
        assert_eq!(unquoted_path("/missing/a /missing/b", "file", "k"), None);
        assert_eq!(
            unquoted_path("/missing/a", "file", "k"),
            Some("k file not found: /missing/a".to_owned())
        );
        assert_eq!(
            unquoted_path("null", "dir", "n"),
            Some("n directory not found: null".to_owned())
        );
    }

    #[test]
    fn grep_dot_matches_any_character() {
        let line = "##bcftools_viewCommand=view -e AF_grpmax_joint > 0.01 -Ou x.vcf.gz";
        assert!(grep_matches(line, "AF_grpmax_joint > 0.01"));
        assert!(grep_matches(
            "AF_grpmax_joint > 0x01",
            "AF_grpmax_joint > 0.01"
        ));
        assert!(!grep_matches(
            "AF_grpmax_joint > 0.05",
            "AF_grpmax_joint > 0.01"
        ));
        assert!(!grep_matches("short", "AF_grpmax_joint > 0.01"));
    }

    #[test]
    fn settings_render_values_as_yq_did() {
        let settings = Settings::parse("a: \"x\"\nb: null\nc: 0.0001\nd: false\ne: 8\n");
        let got: Vec<String> = ["a", "b", "c", "d", "e", "missing"]
            .iter()
            .map(|k| settings.value(k))
            .collect();
        assert_eq!(got, ["x", "null", "0.0001", "false", "8", "null"]);
    }
}
