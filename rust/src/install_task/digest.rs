//! A step's code digest: its registry entry, the text of the installer
//! functions it reaches, the shared shell files' lines outside functions, and
//! the content of the script files and of the shipped data files it names.
//!
//! The installer functions are still shell, so their text is what bash prints
//! for them (`declare -f`) after sourcing `scripts/install_utils.sh`, as the
//! shell runner read it. Two details keep the digests the shell runner wrote:
//!
//! - The shell runner defined its own functions before it sourced the
//!   installer, so a step whose functions use the word `identity` (B4, C6, C8,
//!   C10, C11, D1, D6 and D9 when this was ported) counted the runner's
//!   `identity` function. That text is defined here the same way
//!   (`RUNNER_IDENTITY`); none of the runner's other function names was reached.
//! - The functions, scripts and data files are ordered by the system's `sort`
//!   in the environment's locale, as the shell runner ordered them.
use super::{Options, Result, compact, fail, read_json, step, strings, system_sort};
use crate::json_text::Json;
use crate::sha256::{Sha256, sha256_file};
use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// The shell runner's `identity` function, as `scripts/install_task.sh`
/// defined it (see the module comment).
const RUNNER_IDENTITY: &str = r#"identity() { # path size mtime for a file, path "folder" mtime for a folder, or "missing"
    # A folder's size is not stable on StorNext (0, then 1, for an unchanged
    # LOFTEE repository between two launches of job 1540691).
    local p=$1
    if [[ -d $p ]]; then stat -L -c '%n folder %Y' -- "$p"
    elif [[ -e $p ]]; then stat -L -c '%n %s %Y' -- "$p"; else echo "missing $p"; fi
}"#;

/// The shell files whose functions every step shares. A step counts the
/// functions it reaches, not these whole files; their lines outside any
/// function (comments and blank lines left out) count for every step.
pub(crate) const SHELL_LIBRARIES: [&str; 3] = [
    "install_utils.sh",
    "common_bash_utils.sh",
    "annotation_vcf.sh",
];

/// The installer's shell functions: name to the text `declare -f` prints.
pub(crate) struct Installer {
    scripts: PathBuf,
    root: PathBuf,
    functions: HashMap<String, String>,
}

impl Installer {
    pub(crate) fn load(root: &Path) -> Result<Installer> {
        let scripts = root.join("scripts");
        let program = format!(
            "set -o pipefail\n{RUNNER_IDENTITY}\nsource \"$1/install_utils.sh\" > /dev/null 2>&1 || exit 97\n\
             while read -r _ _ f; do printf '%s\\0' \"$f\"; declare -f \"$f\"; printf '\\0'; done < <(declare -F)\n"
        );
        let output = Command::new("bash")
            .arg("-c")
            .arg(program)
            .arg("install-task")
            .arg(&scripts)
            .stdin(Stdio::null())
            .stderr(Stdio::inherit())
            .output()
            .map_err(|e| fail(format!("cannot run bash: {e}")))?;
        if !output.status.success() {
            return Err(fail("cannot source install_utils.sh"));
        }
        let mut functions = HashMap::new();
        let mut parts = output.stdout.split(|b| *b == 0);
        while let (Some(name), Some(body)) = (parts.next(), parts.next()) {
            if name.is_empty() {
                break;
            }
            functions.insert(
                String::from_utf8_lossy(name).into_owned(),
                String::from_utf8_lossy(body).into_owned(),
            );
        }
        Ok(Installer {
            scripts,
            root: root.to_path_buf(),
            functions,
        })
    }

    pub(crate) fn has(&self, name: &str) -> bool {
        self.functions.contains_key(name)
    }

    fn text(&self, name: &str) -> &str {
        self.functions.get(name).map(String::as_str).unwrap_or("")
    }

    /// The installer functions a step reaches, sorted: from its calls and
    /// its `runs` (installer functions that a program of the step runs
    /// itself), through every function those name.
    pub(crate) fn step_functions(&self, entry: &Json) -> Result<Vec<String>> {
        let runs = strings(entry, "runs");
        for name in &runs {
            if !self.has(name) {
                return Err(fail(format!(
                    "step {} runs unknown installer function {name}",
                    super::text(entry, "id")
                )));
            }
        }
        let mut queue: VecDeque<String> = VecDeque::new();
        if let Some(Json::Array(calls)) = entry.get("calls") {
            for call in calls {
                if let Json::Array(words) = call {
                    if let Some(first) = words.first().and_then(Json::as_str) {
                        queue.push_back(first.to_owned());
                    }
                }
            }
        }
        queue.extend(runs);
        let mut seen: HashSet<String> = HashSet::new();
        let mut found = Vec::new();
        while let Some(name) = queue.pop_front() {
            if !seen.insert(name.clone()) {
                continue;
            }
            if self.has(&name) {
                found.push(name.clone());
            }
            for word in called_words(self.text(&name)) {
                if self.has(&word) && !seen.contains(&word) {
                    queue.push_back(word);
                }
            }
        }
        system_sort(&found, false)
    }

    /// The script files those functions name: by a script file name, or as
    /// SCRIPT_DIR/NAME, and the Python files of the scripts folder that such
    /// a Python script imports (and those import). The shared shell files are
    /// counted by function instead.
    pub(crate) fn step_scripts(&self, names: &[String]) -> Result<Vec<String>> {
        let text: String = names.iter().map(|n| self.text(n)).collect();
        let named: Vec<String> = script_names(&text);
        let named = system_sort(&named, true)?;
        let mut queue: VecDeque<String> = named
            .into_iter()
            .filter(|w| !SHELL_LIBRARIES.contains(&w.as_str()) && self.scripts.join(w).is_file())
            .collect();
        let mut seen: HashSet<String> = HashSet::new();
        let mut found = Vec::new();
        while let Some(name) = queue.pop_front() {
            if !seen.insert(name.clone()) {
                continue;
            }
            found.push(name.clone());
            if !name.ends_with(".py") {
                continue;
            }
            let content = std::fs::read(self.scripts.join(&name)).unwrap_or_default();
            let modules: BTreeSet<String> = python_imports(&content).into_iter().collect();
            let modules = system_sort(&modules.into_iter().collect::<Vec<_>>(), true)?;
            for module in modules {
                let file = format!("{module}.py");
                if self.scripts.join(&file).is_file() && !seen.contains(&file) {
                    queue.push_back(file);
                }
            }
        }
        system_sort(&found, true)
    }

    /// The data files shipped with the repository that those functions name
    /// as DATA_DIR/PATH, with `${assembly}` read as hg19 and as hg38. A file
    /// counts when the source tree has it and, in a git checkout, git tracks
    /// it (an installer run from a checkout may have left its own results in
    /// the data folder).
    pub(crate) fn step_data(&self, names: &[String]) -> Result<Vec<String>> {
        let text: String = names.iter().map(|n| self.text(n)).collect();
        let named = system_sort(&data_names(&text), true)?;
        let git = self.root.join(".git").exists();
        let mut found = Vec::new();
        for name in named {
            for assembly in ["hg19", "hg38"] {
                let path = name.replace("${assembly}", assembly);
                if path.contains('$') || !self.root.join("data").join(&path).is_file() {
                    continue;
                }
                if git {
                    let tracked = Command::new("git")
                        .arg("-C")
                        .arg(&self.root)
                        .args(["ls-files", "--error-unmatch", "--"])
                        .arg(format!("data/{path}"))
                        .stdout(Stdio::null())
                        .stderr(Stdio::null())
                        .status()
                        .is_ok_and(|s| s.success());
                    if !tracked {
                        continue;
                    }
                }
                found.push(path);
            }
        }
        system_sort(&found, true)
    }

    /// The lines of the shared shell files outside any function.
    fn toplevel(&self) -> Result<String> {
        let mut out = String::new();
        for name in SHELL_LIBRARIES {
            let path = self.scripts.join(name);
            let content = std::fs::read(&path)
                .map_err(|e| fail(format!("cannot read {}: {e}", path.display())))?;
            out.push_str(&toplevel_lines(&String::from_utf8_lossy(&content)));
        }
        Ok(out)
    }

    pub(crate) fn digest(&self, entry: &Json) -> Result<String> {
        let names = self.step_functions(entry)?;
        let mut hash = Sha256::new();
        hash.update(compact(entry)?.as_bytes());
        hash.update(b"\n");
        for name in &names {
            hash.update(self.text(name).as_bytes());
        }
        hash.update(b"== outside functions\n");
        hash.update(self.toplevel()?.as_bytes());
        for script in self.step_scripts(&names)? {
            let sum = sha256_file(&self.scripts.join(&script))?;
            hash.update(format!("{script} {sum}  -\n").as_bytes());
        }
        for data in self.step_data(&names)? {
            let sum = sha256_file(&self.root.join("data").join(&data))?;
            hash.update(format!("data {data} {sum}  -\n").as_bytes());
        }
        Ok(hash.finish_hex())
    }
}

fn is_word_start(b: u8) -> bool {
    b.is_ascii_alphabetic() || b == b'_'
}

fn is_word(b: u8) -> bool {
    b.is_ascii_alphanumeric() || b == b'_'
}

fn is_space(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// Whether a line is a message: `^[[:space:]]*(log|echo|printf)([[:space:]]|$)`.
fn is_message(line: &[u8]) -> bool {
    let start = line.iter().take_while(|b| is_space(**b)).count();
    let rest = &line[start..];
    ["log", "echo", "printf"].iter().any(|command| {
        rest.starts_with(command.as_bytes()) && rest.get(command.len()).is_none_or(|b| is_space(*b))
    })
}

/// The words of a function that may name a function it runs. In a message
/// (a log, echo or printf line) only a word where a command can start counts:
/// after `$(`, a backquote, `|`, `&` or `;`. Read whole, the message "run
/// main_install" in priva_resolve_binary tied almost every step to the whole
/// installer.
pub(crate) fn called_words(text: &str) -> BTreeSet<String> {
    let mut words = BTreeSet::new();
    for line in text.as_bytes().split(|b| *b == b'\n') {
        if is_message(line) {
            let mut at = 0;
            while at < line.len() {
                let after = if line[at..].starts_with(b"$(") {
                    Some(at + 2)
                } else if matches!(line[at], b'`' | b'|' | b'&' | b';') {
                    Some(at + 1)
                } else {
                    None
                };
                let matched = after.and_then(|mut end| {
                    while end < line.len() && is_space(line[end]) {
                        end += 1;
                    }
                    if end < line.len() && is_word_start(line[end]) {
                        let start = end;
                        while end < line.len() && is_word(line[end]) {
                            end += 1;
                        }
                        Some((start, end))
                    } else {
                        None
                    }
                });
                match matched {
                    Some((start, end)) => {
                        words.insert(String::from_utf8_lossy(&line[start..end]).into_owned());
                        at = end;
                    }
                    None => at += 1,
                }
            }
        } else {
            for word in line.split(|b| !is_word(*b)) {
                if word.first().is_some_and(|b| is_word_start(*b)) {
                    words.insert(String::from_utf8_lossy(word).into_owned());
                }
            }
        }
    }
    words
}

const SCRIPT_EXTENSIONS: [&str; 5] = ["py", "sh", "R", "pl", "bash"];

fn is_name_char(b: u8) -> bool {
    b.is_ascii_alphanumeric() || matches!(b, b'_' | b'.' | b'-')
}

/// Leftmost-longest matches of `prefix\}?/[class]+` at `at`, as grep -o reads them.
fn prefixed_at(text: &[u8], at: usize, prefix: &[u8], class: fn(u8) -> bool) -> Option<usize> {
    if !text[at..].starts_with(prefix) {
        return None;
    }
    let mut end = at + prefix.len();
    if text.get(end) == Some(&b'}') {
        end += 1;
    }
    if text.get(end) != Some(&b'/') {
        return None;
    }
    end += 1;
    let start = end;
    while end < text.len() && class(text[end]) {
        end += 1;
    }
    (end > start).then_some(end)
}

/// The script names in a text: `grep -oE
/// '[A-Za-z0-9_.-]+\.(py|sh|R|pl|bash)|SCRIPT_DIR\}?/[A-Za-z0-9_./-]+'`, the
/// `SCRIPT_DIR/` prefix removed.
pub(crate) fn script_names(text: &str) -> Vec<String> {
    let text = text.as_bytes();
    let mut found = Vec::new();
    let mut at = 0;
    while at < text.len() {
        let mut best: Option<usize> = None;
        // A file name: the longest prefix of the run of name characters that
        // ends in .EXT with at least one character before the dot.
        let mut run_end = at;
        while run_end < text.len() && is_name_char(text[run_end]) {
            run_end += 1;
        }
        for end in (at + 1..=run_end).rev() {
            let candidate = &text[at..end];
            if SCRIPT_EXTENSIONS.iter().any(|ext| {
                candidate.len() > ext.len() + 1
                    && candidate.ends_with(ext.as_bytes())
                    && candidate[candidate.len() - ext.len() - 1] == b'.'
            }) {
                best = Some(end);
                break;
            }
        }
        let path_class =
            |b: u8| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'.' | b'/' | b'-');
        if let Some(end) = prefixed_at(text, at, b"SCRIPT_DIR", path_class) {
            best = Some(best.map_or(end, |b| b.max(end)));
        }
        match best {
            Some(end) => {
                let word = &text[at..end];
                let word = strip_prefix(word, b"SCRIPT_DIR");
                found.push(String::from_utf8_lossy(word).into_owned());
                at = end;
            }
            None => at += 1,
        }
    }
    found
}

fn strip_prefix<'a>(word: &'a [u8], prefix: &[u8]) -> &'a [u8] {
    if let Some(rest) = word.strip_prefix(prefix) {
        let rest = rest.strip_prefix(b"}").unwrap_or(rest);
        if let Some(rest) = rest.strip_prefix(b"/") {
            return rest;
        }
    }
    word
}

/// The data paths in a text: `grep -oE 'DATA_DIR\}?/[A-Za-z0-9_./${}-]+'`,
/// the `DATA_DIR/` prefix removed.
pub(crate) fn data_names(text: &str) -> Vec<String> {
    let text = text.as_bytes();
    let class = |b: u8| {
        b.is_ascii_alphanumeric() || matches!(b, b'_' | b'.' | b'/' | b'$' | b'{' | b'}' | b'-')
    };
    let mut found = Vec::new();
    let mut at = 0;
    while at < text.len() {
        match prefixed_at(text, at, b"DATA_DIR", class) {
            Some(end) => {
                found.push(
                    String::from_utf8_lossy(strip_prefix(&text[at..end], b"DATA_DIR")).into_owned(),
                );
                at = end;
            }
            None => at += 1,
        }
    }
    found
}

/// The modules a Python file imports at the start of a line:
/// `^[[:space:]]*(from|import)[[:space:]]+NAME`.
fn python_imports(content: &[u8]) -> Vec<String> {
    let mut found = Vec::new();
    for line in content.split(|b| *b == b'\n') {
        let start = line.iter().take_while(|b| is_space(**b)).count();
        let rest = &line[start..];
        let Some(rest) = rest
            .strip_prefix(b"from")
            .or_else(|| rest.strip_prefix(b"import"))
        else {
            continue;
        };
        let spaces = rest.iter().take_while(|b| is_space(**b)).count();
        if spaces == 0 {
            continue;
        }
        let rest = &rest[spaces..];
        if !rest.first().is_some_and(|b| is_word_start(*b)) {
            continue;
        }
        let end = rest.iter().position(|b| !is_word(*b)).unwrap_or(rest.len());
        found.push(String::from_utf8_lossy(&rest[..end]).into_owned());
    }
    found
}

/// The lines of a shell file outside any function, comments and blank lines
/// left out. A function starts at `function NAME` and ends at its last `}` or
/// `)` in column 0 before the next function (main_install has a `(` body).
pub(crate) fn toplevel_lines(content: &str) -> String {
    fn skippable(line: &str) -> bool {
        let rest = line.trim_start_matches(|c: char| c.is_ascii() && is_space(c as u8));
        rest.is_empty() || rest.starts_with('#')
    }
    let mut out = String::new();
    let mut pending: Vec<&str> = Vec::new();
    let mut inside = false;
    let flush = |pending: &mut Vec<&str>, out: &mut String| {
        for line in pending.drain(..) {
            if !skippable(line) {
                out.push_str(line);
                out.push('\n');
            }
        }
    };
    let mut all: Vec<&str> = content.split('\n').collect();
    if content.ends_with('\n') {
        all.pop();
    }
    for line in all {
        let starts_function = line
            .strip_prefix("function ")
            .is_some_and(|rest| rest.as_bytes().first().is_some_and(|b| is_word_start(*b)));
        if starts_function {
            flush(&mut pending, &mut out);
            inside = true;
        } else if inside && (line.starts_with('}') || line.starts_with(')')) {
            pending.clear();
        } else if inside {
            pending.push(line);
        } else if !skippable(line) {
            out.push_str(line);
            out.push('\n');
        }
    }
    flush(&mut pending, &mut out);
    out
}

/// reach: what a step's code digest covers, to see why a step ran again.
pub(crate) fn reach(options: &Options) -> Result<()> {
    let registry = read_json(&options.path("registry")?)?;
    let id = options.one("step")?;
    let installer = Installer::load(&options.source_root()?)?;
    let entry = step(&registry, id)?;
    let names = installer.step_functions(entry)?;
    for name in &names {
        println!("function {name}");
    }
    for script in installer.step_scripts(&names)? {
        println!("script {script}");
    }
    for data in installer.step_data(&names)? {
        println!("data {data}");
    }
    println!("digest {}", installer.digest(entry)?);
    Ok(())
}
