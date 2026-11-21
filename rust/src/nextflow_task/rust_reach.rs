//! The Rust files behind a `priva` command or an annotation step, for the
//! lineages of the install steps and the annotation hub folders: code that
//! runs `priva COMMAND` depends on the Rust code of that command, and an
//! annotation step on its file under pipeline_steps/ and what that names.
//!
//! The files of a command are found from the checkout's rust/src, file by
//! file: the file the command table (`run` in cli.rs) hands the command to,
//! every file of the crate that one names (`crate::`, `super::` and `self::`
//! paths, with `{a, b}` groups), and those they name in turn. A child module
//! (`mod NAME;`) counts when its parent file uses it (`NAME::`), or when it
//! implements a trait for a type defined outside itself, which acts without
//! the child being named. Files a reached file embeds (`include_str!`,
//! `include_bytes!`) count by content.
//!
//! Left out: comment lines, blank lines and test modules (`#[cfg(test)] mod`),
//! and the command table itself, which names every command; cli.rs counts
//! only above `run`. The reach is by file, not by function: a file counts
//! whole when anything in it is named.
use super::Result;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Component, Path, PathBuf};

const TABLE: &str = "cli.rs";
const TABLE_START: &str = "pub(crate) fn run(";

fn is_name(c: char) -> bool {
    c.is_ascii_alphanumeric() || c == '_'
}

/// The lines of a Rust file that count: no comment lines, no blank lines, no
/// test modules.
fn counted(content: &str) -> String {
    let lines: Vec<&str> = content.lines().collect();
    let mut out = Vec::new();
    let mut at = 0;
    while at < lines.len() {
        let line = lines[at];
        let trimmed = line.trim_start();
        let module = lines.get(at + 1).map(|next| next.trim());
        if trimmed == "#[cfg(test)]" && module.is_some_and(|next| next.starts_with("mod ")) {
            let close = format!("{}}}", &line[..line.len() - trimmed.len()]);
            let one_line = module.is_some_and(|next| next.ends_with(';') || next.ends_with('}'));
            at += 2;
            if !one_line {
                while at < lines.len() && lines[at].trim_end() != close {
                    at += 1;
                }
                at += 1;
            }
            continue;
        }
        if !trimmed.is_empty() && !trimmed.starts_with("//") {
            out.push(line.trim_end());
        }
        at += 1;
    }
    out.join("\n")
}

/// Reads one module path, or one `{group}` of them, at `at`; adds every full
/// path to `found` and returns the position after it.
fn path_tree<'a>(
    text: &'a str,
    mut at: usize,
    mut prefix: Vec<&'a str>,
    found: &mut Vec<Vec<&'a str>>,
) -> usize {
    loop {
        if text[at..].starts_with('{') {
            at += 1;
            loop {
                let rest = &text[at..];
                at += rest.len() - rest.trim_start_matches([' ', '\n', '\t', ',']).len();
                match text[at..].chars().next() {
                    None => return at,
                    Some('}') => return at + 1,
                    Some(_) => {
                        let next = path_tree(text, at, prefix.clone(), found);
                        if next == at {
                            return at;
                        }
                        at = next;
                    }
                }
            }
        }
        let rest = &text[at..];
        let end = rest.find(|c: char| !is_name(c)).unwrap_or(rest.len());
        if end == 0 {
            found.push(prefix);
            return at + usize::from(rest.starts_with('*'));
        }
        prefix.push(&rest[..end]);
        at += end;
        if text[at..].starts_with("::") {
            at += 2;
            continue;
        }
        if let Some(alias) = text[at..].strip_prefix(" as ") {
            at += 4 + alias.find(|c: char| !is_name(c)).unwrap_or(alias.len());
        }
        found.push(prefix);
        return at;
    }
}

/// The module paths a text names from `crate::`, `super::` or `self::`, each
/// starting with that word, with groups expanded.
fn module_paths(text: &str) -> Vec<Vec<&str>> {
    let mut found = Vec::new();
    for head in ["crate", "super", "self"] {
        let mut from = 0;
        while let Some(offset) = text[from..].find(head) {
            let start = from + offset;
            from = start + head.len();
            let joined = text[..start]
                .chars()
                .next_back()
                .is_some_and(|c| is_name(c) || c == ':');
            if !joined && text[from..].starts_with("::") {
                path_tree(text, from + 2, vec![head], &mut found);
            }
        }
    }
    found
}

/// The file of the longest leading part of a module path that is a module
/// file, relative to rust/src: `a/b.rs` or `a/b/mod.rs`.
fn module_file(src: &Path, parts: &[&str]) -> Option<String> {
    (1..=parts.len()).rev().find_map(|count| {
        let base = parts[..count].join("/");
        [format!("{base}.rs"), format!("{base}/mod.rs")]
            .into_iter()
            .find(|file| src.join(file).is_file())
    })
}

/// The module path of a file relative to rust/src.
fn own_module(file: &str) -> Vec<&str> {
    let mut parts: Vec<&str> = file.trim_end_matches(".rs").split('/').collect();
    if parts.last() == Some(&"mod") {
        parts.pop();
    }
    parts
}

/// The child modules a file declares with `mod NAME;`.
fn declared_children(text: &str) -> Vec<&str> {
    let mut children = Vec::new();
    for line in text.lines() {
        let mut rest = line.trim();
        if let Some(after) = rest.strip_prefix("pub") {
            rest = match after.strip_prefix('(') {
                Some(scope) => scope.split_once(')').map_or("", |(_, tail)| tail),
                None => after,
            }
            .trim_start();
        }
        let name = rest
            .strip_prefix("mod ")
            .and_then(|tail| tail.strip_suffix(';'));
        if let Some(name) = name.filter(|name| !name.is_empty() && name.chars().all(is_name)) {
            children.push(name);
        }
    }
    children
}

/// Whether a text uses `NAME::` as the start of a path.
fn uses_child(text: &str, name: &str) -> bool {
    let pattern = format!("{name}::");
    let mut from = 0;
    while let Some(offset) = text[from..].find(&pattern) {
        let start = from + offset;
        from = start + pattern.len();
        let joined = text[..start]
            .chars()
            .next_back()
            .is_some_and(|c| is_name(c) || c == ':');
        if !joined {
            return true;
        }
    }
    false
}

/// Whether a file's text implements a trait for a type it does not define
/// itself (`impl Trait for super::Thing`, `for u8`, `for Imported`). Such an
/// implementation acts wherever that type is used. One for the file's own
/// `struct`, `enum`, `union` or `type` needs the file to be named first.
fn acts_unnamed(text: &str) -> bool {
    let defined: BTreeSet<&str> = text
        .lines()
        .filter_map(|line| {
            let mut rest = line.trim_start();
            if let Some(after) = rest.strip_prefix("pub") {
                rest = match after.strip_prefix('(') {
                    Some(scope) => scope.split_once(')').map_or("", |(_, tail)| tail),
                    None => after,
                }
                .trim_start();
            }
            let name = ["struct ", "enum ", "union ", "type "]
                .iter()
                .find_map(|keyword| rest.strip_prefix(keyword))?;
            Some(&name[..name.find(|c: char| !is_name(c)).unwrap_or(name.len())])
        })
        .collect();
    text.lines().any(|line| {
        let line = line.trim_start();
        let line = line.strip_prefix("unsafe ").unwrap_or(line);
        let target = match line.strip_prefix("impl") {
            Some(rest) => rest.rsplit_once(" for ").map(|(_, target)| target),
            // A long header continues on its own line: `    for Type<T>`.
            None => line
                .strip_prefix("for ")
                .filter(|rest| !rest.contains(" in ")),
        };
        let Some(target) = target else {
            return false;
        };
        let end = target.find(|c: char| !is_name(c)).unwrap_or(target.len());
        target[end..].starts_with("::") || !defined.contains(&target[..end])
    })
}

/// The files of the checkout a text embeds, relative to the checkout.
fn embedded(text: &str, file: &Path, root: &Path, out: &mut BTreeSet<String>) {
    for opener in ["include_str!(", "include_bytes!("] {
        for piece in text.split(opener).skip(1) {
            let Some(inner) = piece.trim_start().strip_prefix('"') else {
                continue;
            };
            let Some((name, _)) = inner.split_once('"') else {
                continue;
            };
            let mut path = PathBuf::new();
            for part in file.with_file_name(name).components() {
                match part {
                    Component::ParentDir => {
                        path.pop();
                    }
                    Component::CurDir => {}
                    other => path.push(other),
                }
            }
            if let Ok(relative) = path.strip_prefix(root) {
                if path.is_file() {
                    out.insert(relative.to_string_lossy().into_owned());
                }
            }
        }
    }
}

/// {command: its file relative to rust/src}, read from the command table.
pub(crate) fn commands(root: &Path) -> BTreeMap<String, String> {
    let src = root.join("rust/src");
    let mut out = BTreeMap::new();
    let Ok(table) = std::fs::read_to_string(src.join(TABLE)) else {
        return out;
    };
    let mut pieces = table.split("Some(\"");
    pieces.next();
    for piece in pieces {
        let Some((name, body)) = piece.split_once('"') else {
            continue;
        };
        let command = name.starts_with(|c: char| c.is_ascii_lowercase())
            && name
                .chars()
                .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '-');
        if !command {
            continue;
        }
        let file = module_paths(body)
            .into_iter()
            .filter(|path| path[0] == "crate")
            .find_map(|path| module_file(&src, &path[1..]));
        if let Some(file) = file {
            out.entry(name.to_owned()).or_insert(file);
        }
    }
    out
}

/// The `priva` commands a shell line runs: the word after the program
/// (`"${PRIVA_BIN}" COMMAND`), and the command a tool is handed the program
/// for (`PRIVA_CADD_PRESCORE_BIN=...` is `cadd-prescore`). Only names of the
/// command table count; the caller checks that.
pub(crate) fn called_commands(line: &str) -> Vec<String> {
    let mut found = Vec::new();
    let mut from = 0;
    while let Some(offset) = line[from..].find("PRIVA_") {
        let start = from + offset;
        from = start + "PRIVA_".len();
        if line[..start].chars().next_back().is_some_and(is_name) {
            continue;
        }
        let rest = &line[from..];
        let end = rest.find(|c: char| !is_name(c)).unwrap_or(rest.len());
        let (variable, after) = rest.split_at(end);
        if variable == "BIN" {
            let word = after.trim_start_matches(['}', '"', '\'', ' ']);
            let end = word
                .find(|c: char| !(c.is_ascii_lowercase() || c.is_ascii_digit() || c == '-'))
                .unwrap_or(word.len());
            if word.starts_with(|c: char| c.is_ascii_lowercase()) {
                found.push(word[..end].to_owned());
            }
        } else if let Some(tool) = variable.strip_suffix("_BIN") {
            if after.starts_with('=') {
                found.push(tool.to_ascii_lowercase().replace('_', "-"));
            }
        }
    }
    found
}

/// {file relative to the checkout: SHA-256} for the Rust files behind the
/// command whose file is `entry`, and the files they embed. A Rust file is
/// hashed over its counted lines, an embedded file over its content.
pub(crate) fn command_files(root: &Path, entry: &str) -> Result<BTreeMap<String, String>> {
    let (texts, embedded_files) = reach(root, entry)?;
    let mut out = BTreeMap::new();
    for (file, text) in texts {
        let digest = crate::sha256::sha256_hex(text.as_bytes());
        out.insert(file, digest);
    }
    for file in embedded_files {
        let digest = crate::sha256::sha256_file(&root.join(&file))?;
        out.entry(file).or_insert(digest);
    }
    Ok(out)
}

/// {file relative to the checkout: its counted lines} for the Rust files
/// `entry` reaches, the ones command_files hashes.
pub(crate) fn reached_texts(root: &Path, entry: &str) -> Result<BTreeMap<String, String>> {
    Ok(reach(root, entry)?.0)
}

/// The counted text of each Rust file `entry` reaches, and the files they
/// embed, both relative to the checkout.
fn reach(root: &Path, entry: &str) -> Result<(BTreeMap<String, String>, BTreeSet<String>)> {
    let src = root.join("rust/src");
    let mut texts: BTreeMap<String, String> = BTreeMap::new();
    let mut embedded_files = BTreeSet::new();
    let mut pending = vec![entry.to_owned()];
    while let Some(file) = pending.pop() {
        if texts.contains_key(&file) {
            continue;
        }
        let content = std::fs::read_to_string(src.join(&file))
            .map_err(|error| format!("{}: {error}", src.join(&file).display()))?;
        let mut text = counted(&content);
        if file == TABLE {
            let table = text
                .lines()
                .position(|line| line.starts_with(TABLE_START))
                .ok_or_else(|| format!("rust/src/{TABLE} has no line `{TABLE_START}`"))?;
            text = text.lines().take(table).collect::<Vec<_>>().join("\n");
        }
        let me = own_module(&file);
        for path in module_paths(&text) {
            let levels = path.iter().take_while(|part| **part == "super").count();
            let rest = &path[levels.max(1)..];
            // Inside an inline module `super` is this file's own module, so
            // both readings count.
            let bases: Vec<&[&str]> = match path[0] {
                "crate" => vec![&me[..0]],
                "self" => vec![&me[..]],
                _ => (levels.saturating_sub(1)..=levels)
                    .filter_map(|up| me.len().checked_sub(up).map(|keep| &me[..keep]))
                    .collect(),
            };
            for base in bases {
                let full: Vec<&str> = base.iter().chain(rest).copied().collect();
                pending.extend(module_file(&src, &full));
            }
        }
        for child in declared_children(&text) {
            let mut full = me.clone();
            full.push(child);
            let base = full.join("/");
            let Some(child_file) = [format!("{base}.rs"), format!("{base}/mod.rs")]
                .into_iter()
                .find(|candidate| src.join(candidate).is_file())
            else {
                continue;
            };
            let unnamed = || {
                std::fs::read_to_string(src.join(&child_file))
                    .is_ok_and(|content| acts_unnamed(&counted(&content)))
            };
            if uses_child(&text, child) || unnamed() {
                pending.push(child_file);
            }
        }
        embedded(&text, &src.join(&file), root, &mut embedded_files);
        texts.insert(file, text);
    }
    let texts = texts
        .into_iter()
        .map(|(file, text)| (format!("rust/src/{file}"), text))
        .collect();
    Ok((texts, embedded_files))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn checkout(name: &str, files: &[(&str, &str)]) -> PathBuf {
        let root = std::env::temp_dir().join(format!(
            "priva_rust_reach_{name}_{}_{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        let _ = std::fs::remove_dir_all(&root);
        for (file, content) in files {
            let path = root.join(file);
            std::fs::create_dir_all(path.parent().unwrap()).unwrap();
            std::fs::write(path, content).unwrap();
        }
        root
    }

    const CLI: &str = "\
use crate::shared::Outcome;
pub(crate) fn invalid() {}
pub(crate) fn run() -> Outcome {
    if first() == Some(\"one-thing\") {
        return crate::one::dispatch();
    }
    if first() == Some(\"two\") {
        return crate::pack::two::dispatch();
    }
    crate::other::fallback()
}
";

    const ONE: &str = "\
use crate::cli::invalid;
use crate::{pack::thing, shared::{self as common, Outcome}};
// crate::other is only named in this comment
const MAP: &str = include_str!(\"../../data/map.txt\");
pub(crate) fn dispatch() -> Outcome {
    thing()
}

#[cfg(test)]
mod tests {
    use crate::other;
    fn check() {
        other::fallback();
    }
}
";

    fn files() -> Vec<(&'static str, &'static str)> {
        vec![
            ("rust/src/cli.rs", CLI),
            ("rust/src/one.rs", ONE),
            ("rust/src/shared.rs", "pub(crate) struct Outcome;\n"),
            ("rust/src/other.rs", "pub(crate) fn fallback() {}\n"),
            (
                "rust/src/pack/mod.rs",
                "mod own;\nmod quiet;\nmod shown;\npub(crate) mod two;\nmod used;\n\
                 pub(crate) fn thing() {\n    used::go()\n}\n",
            ),
            ("rust/src/pack/used.rs", "pub(super) fn go() {}\n"),
            ("rust/src/pack/quiet.rs", "pub(super) fn never() {}\n"),
            (
                "rust/src/pack/own.rs",
                "pub(crate) struct Own<T>(T);\nimpl<T> std::fmt::Debug for Own<T> {\n}\n",
            ),
            (
                "rust/src/pack/shown.rs",
                "impl std::fmt::Display for super::Thing {\n}\n",
            ),
            (
                "rust/src/pack/two.rs",
                "pub(crate) fn dispatch() {\n    super::thing()\n}\n",
            ),
            ("data/map.txt", "a\tb\n"),
        ]
    }

    fn names(root: &Path, entry: &str) -> Vec<String> {
        command_files(root, entry).unwrap().into_keys().collect()
    }

    #[test]
    fn a_command_counts_the_files_its_file_names_in_turn() {
        let root = checkout("names", &files());
        let table = commands(&root);
        assert_eq!(table["one-thing"], "one.rs");
        assert_eq!(table["two"], "pack/two.rs");
        assert_eq!(
            names(&root, "one.rs"),
            [
                "data/map.txt",
                "rust/src/cli.rs",
                "rust/src/one.rs",
                "rust/src/pack/mod.rs",
                "rust/src/pack/shown.rs",
                "rust/src/pack/used.rs",
                "rust/src/shared.rs",
            ]
        );
        assert_eq!(
            names(&root, "pack/two.rs"),
            [
                "rust/src/pack/mod.rs",
                "rust/src/pack/shown.rs",
                "rust/src/pack/two.rs",
                "rust/src/pack/used.rs",
            ]
        );
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn tests_comments_and_the_command_table_do_not_count() {
        let root = checkout("stable", &files());
        let before = command_files(&root, "one.rs").unwrap();
        let write = |file: &str, content: String| {
            std::fs::write(root.join(file), content).unwrap();
        };
        write(
            "rust/src/cli.rs",
            CLI.replace("    crate::other", "    // a new command\n    crate::other"),
        );
        write(
            "rust/src/one.rs",
            ONE.replace("other::fallback();", "other::fallback(); // twice"),
        );
        write(
            "rust/src/other.rs",
            "pub(crate) fn fallback() -> u8 { 1 }\n".into(),
        );
        write(
            "rust/src/pack/quiet.rs",
            "pub(super) fn never() -> u8 { 2 }\n".into(),
        );
        assert_eq!(command_files(&root, "one.rs").unwrap(), before);
        // A line above the command table, an embedded file and a used child do.
        for (file, content) in [
            (
                "rust/src/cli.rs",
                CLI.replace("invalid() {}", "invalid() { }"),
            ),
            ("data/map.txt", "a\tc\n".to_owned()),
            (
                "rust/src/pack/used.rs",
                "pub(super) fn go() -> u8 { 3 }\n".to_owned(),
            ),
        ] {
            let kept = std::fs::read_to_string(root.join(file)).unwrap();
            write(file, content);
            let changed: Vec<String> = command_files(&root, "one.rs")
                .unwrap()
                .into_iter()
                .filter(|(name, digest)| before[name] != *digest)
                .map(|(name, _)| name)
                .collect();
            assert_eq!(changed, [file.to_owned()]);
            write(file, kept);
        }
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn a_trait_for_a_type_defined_elsewhere_acts_unnamed() {
        for (text, expected) in [
            (
                "pub struct Own;\nimpl std::fmt::Display for Own {\n}",
                false,
            ),
            (
                "pub(crate) enum Kind<T> {\n}\nimpl<T> From<T> for Kind<T> {\n}",
                false,
            ),
            ("impl std::fmt::Display for super::Thing {\n}", true),
            (
                "use super::Thing;\nimpl std::fmt::Display for Thing {\n}",
                true,
            ),
            ("struct Own;\nimpl Marker for u8 {}", true),
            ("struct Own;\nimpl Marker for &Own {}", true),
            ("struct Own;\nimpl Own {\n    fn before(&self) {}\n}", false),
            (
                "struct Own;\nimpl<T> LongTrait<T>\n    for Own\n{\n}",
                false,
            ),
            ("impl<T> LongTrait<T>\n    for Other<T>\n{\n}", true),
            (
                "struct Own;\nfn f() {\n    for item in items {\n    }\n}",
                false,
            ),
            ("unsafe impl Send for super::Thing {}", true),
        ] {
            assert_eq!(acts_unnamed(text), expected, "{text}");
        }
    }

    #[test]
    fn module_paths_expand_groups() {
        let text = "use crate::{a, b::{c, d as e}, f::*};\nlet x = self::g::h(super::super::i);";
        let mut paths = module_paths(text);
        paths.sort();
        assert_eq!(
            paths,
            [
                vec!["crate", "a"],
                vec!["crate", "b", "c"],
                vec!["crate", "b", "d"],
                vec!["crate", "f"],
                vec!["self", "g", "h"],
                vec!["super", "super", "i"],
            ]
        );
        assert!(module_paths("myself::a, Self::b, x::crate::c").is_empty());
    }

    #[test]
    fn shell_lines_name_the_commands_they_run() {
        for (line, expected) in [
            ("    \"${PRIVA_BIN}\" variant-ids \\", vec!["variant-ids"]),
            ("    $PRIVA_BIN cadd-hub match --hub x", vec!["cadd-hub"]),
            (
                "    \"${PRIVA_BIN}\" clinvar-gene-refresh \"candidate-${phase}\" \\",
                vec!["clinvar-gene-refresh"],
            ),
            (
                "    PRIVA_CADD_PRESCORE_BIN=${prescore_bin} \\",
                vec!["cadd-prescore"],
            ),
            ("    PRIVA_BIN=${candidate}", vec![]),
            ("        prescore_bin=${PRIVA_BIN}", vec![]),
            ("    \"${PRIVA_BIN}\" --version", vec![]),
            ("    if [[ -n ${PRIVA_BINARY:-} ]]; then", vec![]),
            ("    MY_PRIVA_BIN x", vec![]),
        ] {
            assert_eq!(called_commands(line), expected, "{line}");
        }
    }
}
