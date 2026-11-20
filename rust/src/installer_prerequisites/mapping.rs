//! NM-to-ENST mapping selection and checks (`nm_to_enst_path`, `validate_mapping`,
//! `mapping_fetch_source`, `mapping_workspace`). Pinned BioMart facts and the
//! pinned-provenance verification come from the native BioMart module.
use super::paths::{
    abspath, checked_text, exists, is_relative_to, is_symlink, join, lexists, name, open, pathlib,
    read_failure, readable_text, resolve,
};
use super::{Config, Result, Value, error, require, value_error};
use crate::biomart::{self, PinnedSource};
use crate::json_text::Json;
use crate::python_lines::Lines;
use crate::python_text::is_python_space;
use std::io::Read;
use std::path::Path;

/// `_mapping_default_source(config)`
fn default_source(config: &Config) -> Result<PinnedSource> {
    let release = match config.get("vep_installed_version") {
        value if value.is_unset() => None,
        Value::Bool(flag) => Some(Json::Bool(*flag)),
        value => Some(Json::Str(value.py_str()?)),
    };
    biomart::pinned_source_for_assembly(&config.get("assembly").py_str()?, release.as_ref())
}

/// `_mapping_default_path(config, source_root, selected)`: pinned mappings
/// belong to source data, never beneath cache_root/data.
fn default_path(root: &str, selected: &PinnedSource) -> Result<String> {
    abspath(&join(root, &selected.relative_path))
}

/// `config.get("assembly") in named["aliases"]`
fn compatible(named: &PinnedSource, assembly: &Value) -> bool {
    assembly
        .text()
        .is_some_and(|text| named.aliases.iter().any(|alias| alias == text))
}

/// The checked `<mapping>.provenance.json` sidecar path.
fn provenance(path: &str) -> Result<String> {
    checked_text(
        &format!("{path}.provenance.json"),
        "nm_to_enst_tsv provenance",
        0,
    )
}

/// `nm_to_enst_path(config, source_root)`: explicit inputs are kept; a new
/// default follows the assembly and selected VEP release.
pub(super) fn nm_to_enst_path(config: &Config, root: &str) -> Result<String> {
    let value = match config.get("nm_to_enst_tsv") {
        selected if selected.is_unset() => default_path(root, &default_source(config)?)?,
        Value::Str(text) if text.starts_with('/') => text.clone(),
        Value::Str(text) => join(root, text),
        _ => return Err(value_error("nm_to_enst_tsv must be a path string")),
    };
    // Keep the logical file name and adjacent provenance of a DataLad payload
    // resolved into .git/annex; the path policy still checks every link.
    checked_text(&value, "nm_to_enst_tsv", 0)?;
    let path = abspath(&value)?;
    let assembly = config.get("assembly");
    if let Some(named) = biomart::pinned_source_for_basename(name(&path))
        && !compatible(&named, assembly)
    {
        return Err(value_error(format!(
            "Pinned BioMart mapping {} is not compatible with {}; \
             supply an explicitly assembly-compatible nm_to_enst_tsv",
            name(&path),
            assembly.py_str()?
        )));
    }
    Ok(path)
}

/// `mapping_fetch_source(config, source_root)`: acquisition is authorized only
/// for the missing, automatically planned destination.
pub(super) fn mapping_fetch_source(config: &Config, root: &str) -> Result<String> {
    let target = nm_to_enst_path(config, root)?;
    let missing = |what: &str| {
        value_error(format!(
            "Explicit nm_to_enst_tsv is missing: {target}; restore that selected {what}"
        ))
    };
    if is_symlink(&target)? && !exists(&target)? {
        return Err(missing("linked input"));
    }
    if !config.get("nm_to_enst_tsv").is_unset() {
        return Err(missing("input"));
    }
    let selected = default_source(config)?;
    if target != default_path(root, &selected)? {
        return Err(missing("input"));
    }
    if lexists(&target) {
        return Err(value_error(format!(
            "Bundled nm_to_enst_tsv already exists; validate or recover it explicitly: {target}"
        )));
    }
    let metadata = provenance(&target)?;
    if lexists(&metadata) {
        return Err(value_error(format!(
            "Missing bundled nm_to_enst_tsv has existing provenance: {metadata}; \
             retain it and recover explicitly"
        )));
    }
    Ok(selected.source_id)
}

/// `mapping_workspace(config, source_root)`: the raw acquisition folder outside
/// the checkout, planned without creating anything.
pub(super) fn mapping_workspace(config: &Config, root: &str) -> Result<String> {
    let mut selected = config.get("cache_root");
    if selected.is_unset() {
        selected = config.get("tmp_dir");
    }
    let base = match selected {
        Value::Str(text) if !selected.is_unset() => pathlib(text),
        _ => {
            return Err(value_error(
                "Set cache_root or tmp_dir for an external BioMart acquisition workspace",
            ));
        }
    };
    let base = if base.starts_with('/') {
        base
    } else {
        join(root, &base)
    };
    let source = default_source(config)?;
    let workspace = abspath(&join(&join(&base, "BioMart"), &source.source_id))?;
    require(
        !is_relative_to(&resolve(&workspace)?, &resolve(root)?),
        "BioMart acquisition workspace must be outside the checkout",
    )?;
    checked_text(&workspace, "BioMart acquisition workspace", 0)?;
    Ok(workspace)
}

/// `validate_mapping(config, source_root)`
pub(super) fn validate_mapping(config: &Config, root: &str) -> Result<String> {
    let path = nm_to_enst_path(config, root)?;
    validate_nm_to_enst(&path, config.get("assembly"), root)?;
    Ok(path)
}

/// `BIOMART.is_bundled_mapping(mapping, source_root)` for an absolute path.
fn is_bundled_mapping(logical: &str, root: &str) -> Result<bool> {
    let basename = name(logical);
    let relative = format!("data/hgnc/{basename}");
    Ok(biomart::pinned_source_for_basename(basename).is_some()
        && logical == abspath(&join(root, &relative))?
        && !is_symlink(logical)?
        && resolve(logical)? == join(&resolve(root)?, &relative))
}

fn strip(text: &str) -> &str {
    text.trim_matches(is_python_space)
}

/// `validate_nm_to_enst(value, assembly, source_root)`
fn validate_nm_to_enst(value: &str, assembly: &Value, root: &str) -> Result<()> {
    // Check the policy before opening the chosen path or its sidecar.
    readable_text(value, "nm_to_enst_tsv")?;
    let path = abspath(value)?;
    let named = biomart::pinned_source_for_basename(name(&path));
    if let Some(named) = &named
        && *assembly != Value::None
        && !compatible(named, assembly)
    {
        return Err(value_error(format!(
            "Pinned BioMart mapping {} is not compatible with {}",
            name(&path),
            assembly.py_str()?
        )));
    }
    if is_bundled_mapping(&path, root)? {
        // Canonical bundled snapshots need only exist and be nonempty.
        return Ok(());
    }
    let metadata = provenance(&path)?;
    if exists(&metadata)? || is_symlink(&metadata)? || named.is_some() {
        let selected = match assembly {
            Value::None => None,
            value => Some(value.py_str()?),
        };
        return biomart::verify_pinned_mapping(Path::new(&path), selected.as_deref());
    }
    // Both unchanged scientific readers use the original columns 2 and 5.
    let mut rows = CsvRows::new(open(&path)?, &path);
    let header = rows.next()?.unwrap_or_default();
    let layout = header.len() >= 5
        && strip(&header[1]) == "Transcript stable ID"
        && strip(&header[4]) == "RefSeq mRNA ID";
    if !layout {
        return Err(value_error(format!(
            "nm_to_enst_tsv must retain the original five-column Ensembl/RefSeq layout: {value}"
        )));
    }
    while let Some(row) = rows.next()? {
        if row.len() >= 5 && row[1].starts_with("ENST") && row[4].starts_with("NM_") {
            return Ok(());
        }
    }
    Err(value_error(format!(
        "nm_to_enst_tsv contains no usable original NM-to-ENST rows: {value}"
    )))
}

/// Python's default csv field size limit, in characters.
const FIELD_LIMIT: usize = 128 * 1024;

#[derive(Clone, Copy, PartialEq)]
enum State {
    StartRecord,
    StartField,
    InField,
    InQuotedField,
    QuoteInQuotedField,
    EatCrnl,
}

/// One row being read by the `_csv` reader state machine (excel dialect, tab
/// delimiter, `strict=False`).
struct Row {
    fields: Vec<String>,
    field: String,
    length: usize,
    state: State,
}

impl Row {
    fn save(&mut self) {
        self.fields.push(std::mem::take(&mut self.field));
        self.length = 0;
    }
    fn add(&mut self, c: char) -> Result<()> {
        if self.length >= FIELD_LIMIT {
            return Err(error(
                "Error",
                format!("field larger than field limit ({FIELD_LIMIT})"),
            ));
        }
        self.field.push(c);
        self.length += 1;
        Ok(())
    }
    /// `parse_process_char`; `None` is the end-of-line marker after each line.
    fn step(&mut self, c: Option<char>) -> Result<()> {
        let end = matches!(c, None | Some('\n' | '\r'));
        if self.state == State::StartRecord {
            match c {
                // An empty line is the row [].
                None => return Ok(()),
                Some(_) if end => {
                    self.state = State::EatCrnl;
                    return Ok(());
                }
                Some(_) => self.state = State::StartField,
            }
        }
        if end && !matches!(self.state, State::InQuotedField | State::EatCrnl) {
            self.save();
            self.state = if c.is_none() {
                State::StartRecord
            } else {
                State::EatCrnl
            };
            return Ok(());
        }
        let Some(c) = c else {
            if self.state == State::EatCrnl {
                self.state = State::StartRecord;
            }
            return Ok(());
        };
        match (self.state, c) {
            (State::StartField, '"') => self.state = State::InQuotedField,
            (State::InQuotedField, '"') => self.state = State::QuoteInQuotedField,
            (State::QuoteInQuotedField, '"') => {
                self.add('"')?;
                self.state = State::InQuotedField;
            }
            (State::StartField | State::InField | State::QuoteInQuotedField, '\t') => {
                self.save();
                self.state = State::StartField;
            }
            (State::StartField | State::InField | State::QuoteInQuotedField, c) => {
                self.add(c)?;
                self.state = State::InField;
            }
            (State::InQuotedField, c) => self.add(c)?,
            (State::EatCrnl, '\n' | '\r') => {}
            (State::EatCrnl | State::StartRecord, _) => {
                return Err(error(
                    "Error",
                    "new-line character seen in unquoted field - \
                     do you need to open the file with newline=''?",
                ));
            }
        }
        Ok(())
    }
}

/// `csv.reader(open(path, newline="", encoding="utf-8-sig"), delimiter="\t")`
struct CsvRows<R: Read> {
    lines: Lines<R>,
    first: bool,
    path: String,
}

impl<R: Read> CsvRows<R> {
    fn new(source: R, path: &str) -> Self {
        Self {
            lines: Lines::new(source),
            first: true,
            path: path.to_owned(),
        }
    }
    fn next(&mut self) -> Result<Option<Vec<String>>> {
        let mut row = Row {
            fields: Vec::new(),
            field: String::new(),
            length: 0,
            state: State::StartRecord,
        };
        loop {
            let line = self
                .lines
                .next_csv(&mut self.first)
                .map_err(|failure| read_failure(failure, &self.path))?;
            let Some(line) = line else {
                // Input ends inside a quoted field or after an unsaved field.
                if row.length != 0 || row.state == State::InQuotedField {
                    row.save();
                    return Ok(Some(row.fields));
                }
                return Ok(None);
            };
            for c in line.chars() {
                row.step(Some(c))?;
            }
            row.step(None)?;
            if row.state == State::StartRecord {
                return Ok(Some(row.fields));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rows(text: &str) -> Result<Vec<Vec<String>>> {
        let mut reader = CsvRows::new(text.as_bytes(), "fixture.tsv");
        let mut rows = Vec::new();
        while let Some(row) = reader.next()? {
            rows.push(row);
        }
        Ok(rows)
    }

    #[test]
    fn csv_rows_follow_the_python_reader() {
        let text = "\u{feff}a\tb\r\n\nq\"x\t\"t\tu\"\"v\"w\n\"open\nline\"\tz\rlast";
        let parsed = rows(text).unwrap();
        let expected: Vec<Vec<&str>> = vec![
            vec!["a", "b"],
            vec![],
            vec!["q\"x", "t\tu\"vw"],
            vec!["open\nline", "z"],
            vec!["last"],
        ];
        assert_eq!(parsed, expected);
        assert_eq!(rows("\"unterminated\n").unwrap(), [["unterminated\n"]]);
        assert_eq!(rows("a\t\n").unwrap(), [["a", ""]]);
        let long = "x".repeat(FIELD_LIMIT + 1);
        let failure = rows(&long).err().unwrap();
        assert_eq!(
            failure.to_string(),
            "Error: field larger than field limit (131072)"
        );
    }
}
