//! Checked Rust-native InterPro mapping access for the domain and splicing rules.
//!
//! The complete schema-9 archive stays memory-mapped. A key lookup lazily copies only the
//! first four texts retained by those rules, preserving entry order. The producer's GO terms
//! are tuples inside a list: Python's list-only GO-term rule skips them, as this view does.
//! Full six-tuples, GO triples and abstracts remain available through the typed cache/export.
use crate::cache_convert::interpro_mapping::InterproCache;
use crate::cli::{Result, invalid};
#[cfg(test)]
use crate::pickle_reader::Value;
#[cfg(test)]
use std::collections::HashMap;
use std::path::Path;
use std::sync::OnceLock;

const USAGE: &str = "usage: priva interpro-mapping inspect --input INPUT";

#[cfg_attr(not(test), allow(dead_code))]
pub(crate) enum GoTerm {
    /// Not a list of three or more items, so no rule applies to it.
    Skipped,
    /// A list `[go_id, category, description, ...]`; `None` is a non-text item.
    Annotated {
        category: Option<String>,
        description: Option<String>,
    },
}

#[cfg_attr(not(test), allow(dead_code))]
pub(crate) enum GoTerms {
    NotAList,
    List(Vec<GoTerm>),
}

pub(crate) struct UsableEntry {
    pub(crate) id: String,
    pub(crate) kind: String,
    pub(crate) short_name: String,
    pub(crate) full_name: String,
    pub(crate) go_terms: GoTerms,
}

#[cfg_attr(not(test), allow(dead_code))]
pub(crate) enum InterproEntry {
    Usable(UsableEntry),
    /// Python's `is_functional_domain` raises ValueError for this entry.
    Unusable,
}

enum MappingSource {
    Mapped {
        cache: InterproCache,
        entries: Vec<OnceLock<Vec<InterproEntry>>>,
    },
    #[cfg(test)]
    Fixture(HashMap<String, HashMap<String, Vec<InterproEntry>>>),
}

pub(crate) struct InterproMapping {
    source: MappingSource,
}

impl InterproMapping {
    pub(crate) fn get(&self, database: &str, accession: &str) -> Option<&Vec<InterproEntry>> {
        match &self.source {
            MappingSource::Mapped { cache, entries } => {
                let at = cache.position(database, accession)?;
                Some(entries[at].get_or_init(|| {
                    cache
                        .get(database, accession)
                        .unwrap()
                        .iter()
                        .map(|entry| {
                            InterproEntry::Usable(UsableEntry {
                                id: cache.text(&entry.id).to_owned(),
                                kind: cache.text(&entry.kind).to_owned(),
                                short_name: cache.text(&entry.short_name).to_owned(),
                                full_name: cache.text(&entry.name).to_owned(),
                                go_terms: GoTerms::List(
                                    entry.go_terms.iter().map(|_| GoTerm::Skipped).collect(),
                                ),
                            })
                        })
                        .collect()
                }))
            }
            #[cfg(test)]
            MappingSource::Fixture(entries) => entries.get(database)?.get(accession),
        }
    }

    #[cfg(test)]
    pub(crate) fn from_entries(rows: Vec<(&str, &str, Vec<InterproEntry>)>) -> Self {
        let mut entries: HashMap<String, HashMap<String, Vec<InterproEntry>>> = HashMap::new();
        for (database, accession, list) in rows {
            entries
                .entry(database.to_owned())
                .or_default()
                .insert(accession.to_owned(), list);
        }
        Self {
            source: MappingSource::Fixture(entries),
        }
    }
}

#[cfg(test)]
fn go_term_from_value(value: &Value) -> GoTerm {
    match value {
        Value::List(items) => {
            let items = items.borrow();
            if items.len() >= 3 {
                GoTerm::Annotated {
                    category: items[1].as_str().map(str::to_owned),
                    description: items[2].as_str().map(str::to_owned),
                }
            } else {
                GoTerm::Skipped
            }
        }
        _ => GoTerm::Skipped,
    }
}

#[cfg(test)]
fn entry_from_value(value: &Value) -> InterproEntry {
    let Value::Tuple(fields) = value else {
        return InterproEntry::Unusable;
    };
    if fields.len() < 5 {
        return InterproEntry::Unusable;
    }
    let (Some(id), Some(kind), Some(short_name), Some(full_name)) = (
        fields[0].as_str(),
        fields[1].as_str(),
        fields[2].as_str(),
        fields[3].as_str(),
    ) else {
        return InterproEntry::Unusable;
    };
    let go_terms = match &fields[4] {
        Value::List(items) => {
            GoTerms::List(items.borrow().iter().map(go_term_from_value).collect())
        }
        _ => GoTerms::NotAList,
    };
    InterproEntry::Usable(UsableEntry {
        id: id.to_owned(),
        kind: kind.to_owned(),
        short_name: short_name.to_owned(),
        full_name: full_name.to_owned(),
        go_terms,
    })
}

/// Maps and validates the complete schema-9 archive; Python pickles are refused.
pub(crate) fn load_interpro_mapping(path: &Path) -> Result<InterproMapping> {
    let path = path
        .to_str()
        .ok_or_else(|| invalid("the InterPro cache path is not UTF-8"))?;
    let cache = InterproCache::open(path)?;
    let entries = (0..cache.len()).map(|_| OnceLock::new()).collect();
    Ok(InterproMapping {
        source: MappingSource::Mapped { cache, entries },
    })
}

/// Python's `json.dumps` of a string (ASCII only, as `ensure_ascii` writes it).
fn json_string(text: &str) -> String {
    let mut out = String::from("\"");
    for character in text.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{8}' => out.push_str("\\b"),
            '\u{c}' => out.push_str("\\f"),
            c if (c as u32) < 0x20 || (c as u32) > 0x7e => {
                let mut units = [0_u16; 2];
                for unit in c.encode_utf16(&mut units).iter() {
                    out.push_str(&format!("\\u{unit:04x}"));
                }
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// `str(Path(argument))`: repeated slashes and `.` parts are dropped.
fn path_text(argument: &str) -> String {
    let absolute = argument.starts_with('/');
    let parts: Vec<&str> = argument
        .split('/')
        .filter(|part| !part.is_empty() && *part != ".")
        .collect();
    let joined = parts.join("/");
    if absolute {
        // POSIX keeps exactly two leading slashes, and one for any other count.
        let leading = argument.len() - argument.trim_start_matches('/').len();
        format!("{}{joined}", if leading == 2 { "//" } else { "/" })
    } else if joined.is_empty() {
        ".".to_owned()
    } else {
        joined
    }
}

/// The metadata line of `inspect`: the pickle's type and, for a sized value, its length.
pub(crate) fn inspect_mapping(path: &str) -> Result<String> {
    let entries = InterproCache::open(path)?.len();
    Ok(format!(
        "{{\"entries\": {entries}, \"input\": {}, \"value_type\": {}}}",
        json_string(&path_text(path)),
        json_string("dict")
    ))
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args == ["--help"]
        || args == ["-h"]
        || args == ["inspect", "--help"]
        || args == ["inspect", "-h"]
    {
        println!("{USAGE}");
        return Ok(());
    }
    if args.first().map(String::as_str) != Some("inspect") {
        return Err(invalid(format!("the command must be inspect\n{USAGE}")));
    }
    let mut input: Option<String> = None;
    let mut index = 1;
    while index < args.len() {
        let (option, inline) = match args[index].split_once('=') {
            Some((name, value)) if name.starts_with("--") => (name, Some(value.to_owned())),
            _ => (args[index].as_str(), None),
        };
        // argparse accepts any unambiguous prefix of a long option.
        if option.len() <= 2 || !"--input".starts_with(option) {
            return Err(invalid(format!("unknown argument {option:?}\n{USAGE}")));
        }
        input = Some(match inline {
            Some(value) => value,
            None => {
                index += 1;
                args.get(index)
                    .cloned()
                    .ok_or_else(|| invalid(format!("{option} needs a value\n{USAGE}")))?
            }
        });
        index += 1;
    }
    let input = input.ok_or_else(|| invalid(format!("--input is required\n{USAGE}")))?;
    match inspect_mapping(&input) {
        Ok(line) => {
            println!("{line}");
            Ok(())
        }
        Err(error) => {
            eprintln!("ERROR: {error}");
            std::process::exit(1)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn json_strings_are_ascii_like_python() {
        assert_eq!(json_string("a/b"), "\"a/b\"");
        assert_eq!(json_string("q\"\\\n\t"), "\"q\\\"\\\\\\n\\t\"");
        assert_eq!(json_string("é"), "\"\\u00e9\"");
        assert_eq!(json_string("\u{1f600}"), "\"\\ud83d\\ude00\"");
        assert_eq!(json_string("\u{1}\u{7f}"), "\"\\u0001\\u007f\"");
    }

    #[test]
    fn paths_print_like_pathlib() {
        assert_eq!(path_text("/a//b/./c/"), "/a/b/c");
        assert_eq!(path_text("./x.pkl"), "x.pkl");
        assert_eq!(path_text("."), ".");
        assert_eq!(path_text("//a"), "//a");
        assert_eq!(path_text("///a"), "/a");
        assert_eq!(path_text("/"), "/");
    }

    #[test]
    fn entries_keep_the_types_the_rules_look_at() {
        use std::cell::RefCell;
        use std::rc::Rc;
        let text = |value: &str| Value::Str(Rc::from(value));
        let go_list = Value::List(Rc::new(RefCell::new(vec![
            // A tuple, as the installer writes GO terms: no rule applies.
            Value::Tuple(Rc::new(vec![
                text("GO:1"),
                text("molecular_function"),
                text("binding"),
            ])),
            Value::List(Rc::new(RefCell::new(vec![
                text("GO:2"),
                text("biological_process"),
                text("signaling"),
            ]))),
        ])));
        let entry = Value::Tuple(Rc::new(vec![
            text("IPR000001"),
            text("Domain"),
            text("Short"),
            text("Long name"),
            go_list,
            text("abstract"),
        ]));
        let InterproEntry::Usable(usable) = entry_from_value(&entry) else {
            panic!("a six-item tuple of text and a list is usable")
        };
        assert_eq!(usable.id, "IPR000001");
        let GoTerms::List(terms) = &usable.go_terms else {
            panic!()
        };
        assert!(matches!(terms[0], GoTerm::Skipped));
        assert!(matches!(
            &terms[1],
            GoTerm::Annotated { category: Some(category), .. } if category == "biological_process"
        ));
        let short = Value::Tuple(Rc::new(vec![text("a"), text("b")]));
        assert!(matches!(entry_from_value(&short), InterproEntry::Unusable));
        assert!(matches!(
            entry_from_value(&text("x")),
            InterproEntry::Unusable
        ));
    }

    #[test]
    fn the_bundled_mapping_is_read() {
        let path = crate::protein_domain_mapping::data_directory()
            .join("InterPro")
            .join("Interpro_entry_mapping.v1.rkyv");
        let line = inspect_mapping(&path.to_string_lossy()).unwrap();
        assert!(
            line.starts_with("{\"entries\": 60301, \"input\": "),
            "{line}"
        );
        assert!(line.ends_with("\"value_type\": \"dict\"}"), "{line}");
        let mapping = load_interpro_mapping(&path).unwrap();
        let MappingSource::Mapped { entries, .. } = &mapping.source else {
            panic!("native map")
        };
        assert!(entries.iter().all(|entry| entry.get().is_none()));
        let entries = mapping.get("PFAM", "PF00001").unwrap();
        let InterproEntry::Usable(first) = &entries[0] else {
            panic!("the first PF00001 entry is usable")
        };
        assert_eq!(first.id, "IPR000276");
        assert_eq!(first.short_name, "GPCR_Rhodpsn");
        assert!(matches!(first.go_terms, GoTerms::List(_)));
        assert!(
            matches!(&first.go_terms, GoTerms::List(terms) if terms.iter().all(|term| matches!(term, GoTerm::Skipped)))
        );
        let MappingSource::Mapped { entries, .. } = &mapping.source else {
            unreachable!()
        };
        assert_eq!(
            entries.iter().filter(|entry| entry.get().is_some()).count(),
            1
        );
    }
}
