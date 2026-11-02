//! Direct schema-9 producer for protein_domain_mapping.parse_interpro_xml_for_mapping.
//!
//! This procedure performs only the source's full-mapping XML pass. It does not
//! replace the separate VCF-format inventory/report workflow in that script.
//! The original mature Expat backend owns XML syntax, encodings and namespaces;
//! this module owns only the source's element selection and cache semantics.
use crate::builder_error::{category, finish};
use crate::cache_convert::Provenance;
use crate::cache_convert::interpro_mapping::{PlainEntry, PlainMapping, write_produced};
use crate::cli::{Result, invalid};
use crate::python_args::{Parsed, parse_long_options};
use crate::python_io::PythonGzipReader;
use crate::python_text::python_strip;
use std::cell::{Ref, RefCell};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::rc::Rc;
mod encodings;
mod xml;

/// Reuse the mature XML backend for the acquisition helper's independent
/// Product 6 shape check. It does not construct or publish an InterPro mapping.
pub(crate) fn validate_product6(source: impl Read) -> Result<()> {
    xml::validate_product6(source)
}

const HELP: &str = "usage: priva interpro-mapping build --xml XML.gz --output CACHE.rkyv\n\
Build the full member-signature mapping directly from an explicit gzipped InterPro XML.\n\
No download, VCF-format report, installer caller change or pickle output.";

#[derive(Default)]
struct Element {
    name: String,
    attributes: BTreeMap<String, String>,
    text: String,
    tail: String,
    children: Vec<Rc<RefCell<Element>>>,
}

impl Element {
    fn child(&self, name: &str) -> Option<Ref<'_, Element>> {
        self.children.iter().find_map(|child| {
            let child = child.borrow();
            (child.name == name).then_some(child)
        })
    }

    fn child_text(&self, name: &str) -> String {
        self.child(name)
            .map_or_else(String::new, |child| python_strip(&child.text).to_owned())
    }

    fn append_text(&mut self, text: &str) {
        // ElementTree concatenates text/CDATA/entity callbacks until an actual
        // child boundary. Comments and processing instructions add no child.
        match self.children.last_mut() {
            Some(child) => child.borrow_mut().tail.push_str(text),
            None => self.text.push_str(text),
        }
    }

    fn itertext(&self, parts: &mut Vec<String>) {
        if !self.text.is_empty() {
            parts.push(self.text.clone());
        }
        for child in &self.children {
            let child = child.borrow();
            child.itertext(parts);
            if !child.tail.is_empty() {
                parts.push(child.tail.clone());
            }
        }
    }

    fn clear(&mut self) {
        self.attributes.clear();
        self.text.clear();
        self.tail.clear();
        self.children.clear();
    }
}

type Mapping = BTreeMap<(String, String), Vec<PlainEntry>>;

fn record(element: &Element, result: &mut Mapping) {
    let get = |key: &str| element.attributes.get(key);
    let or = |key: &str, default: &str| {
        get(key)
            .filter(|value| !value.is_empty())
            .map_or(default, String::as_str)
            .to_owned()
    };
    let mut go_terms = Vec::new();
    if let Some(classes) = element.child("class_list") {
        for class in classes
            .children
            .iter()
            .map(|c| c.borrow())
            .filter(|c| c.name == "classification")
        {
            if class.attributes.get("class_type").map(String::as_str) == Some("GO") {
                if let Some(id) = class.attributes.get("id").filter(|id| !id.is_empty()) {
                    go_terms.push((
                        id.clone(),
                        class.child_text("category"),
                        class.child_text("description"),
                    ));
                }
            }
        }
    }
    let mut text_parts = Vec::new();
    if let Some(abstract_element) = element.child("abstract") {
        abstract_element.itertext(&mut text_parts);
    }
    let entry = PlainEntry {
        id: or("id", "UnknownIPR"),
        kind: or("type", "UnknownType"),
        short_name: or("short_name", "UnknownShortName"),
        name: element.child_text("name"),
        go_terms,
        // XML 1.0 excludes U+001C..U+001F, the only Python whitespace
        // characters absent from Rust's Unicode White_Space property.
        abstract_text: text_parts
            .join(" ")
            .split_whitespace()
            .collect::<Vec<_>>()
            .join(" "),
    };
    // Presence, not truthiness: an explicitly empty id still publishes entries
    // whose id field is UnknownIPR, exactly as the source does.
    if get("id").is_some() {
        if let Some(members) = element.child("member_list") {
            for member in members
                .children
                .iter()
                .map(|m| m.borrow())
                .filter(|m| m.name == "db_xref")
            {
                if let (Some(db), Some(key)) = (
                    member.attributes.get("db").filter(|x| !x.is_empty()),
                    member.attributes.get("dbkey").filter(|x| !x.is_empty()),
                ) {
                    result
                        .entry((db.to_uppercase(), key.clone()))
                        .or_default()
                        .push(entry.clone());
                }
            }
        }
    }
}

fn parse(source: impl Read) -> Result<PlainMapping> {
    let result = xml::parse(source)?;
    if result.is_empty() {
        return Err(category(
            "ValueError",
            invalid("No usable InterPro member mappings were generated"),
        ));
    }
    Ok(result.into_iter().collect())
}

fn build(xml: &Path, output: &Path) -> Result<()> {
    let hash = crate::sha256::sha256_file(xml)?;
    let mapping = parse(PythonGzipReader::new(File::open(xml)?))?;
    if crate::sha256::sha256_file(xml)? != hash {
        return Err(invalid(
            "the InterPro XML changed during parsing; no cache published",
        ));
    }
    crate::rust_cache::ensure_distinct_paths(xml, output)?;
    let provenance = vec![
        Provenance {
            key: "producer".into(),
            value: "protein_domain_mapping.parse_interpro_xml_for_mapping".into(),
        },
        Provenance {
            key: "source_sha256".into(),
            value: hash,
        },
    ];
    write_produced(output, &mapping, provenance)?;
    eprintln!("XML backend: {}", xml::version());
    println!(
        "{{\"cache\":\"interpro_mapping_pickle\",\"keys\":{},\"records\":{}}}",
        mapping.len(),
        mapping
            .iter()
            .map(|(_, values)| values.len())
            .sum::<usize>()
    );
    Ok(())
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Parsed::Values(values) = parse_long_options(args, &["xml", "output"]).map_err(invalid)?
    else {
        println!("{HELP}");
        return Ok(());
    };
    let required = |key: &str| {
        values
            .get(key)
            .ok_or_else(|| invalid(format!("missing --{key}")))
    };
    let xml = required("xml")?;
    let output = required("output")?;
    finish(build(Path::new(xml), Path::new(output)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn direct_children_order_duplicates_and_defaults() {
        let entries = parse(br#"<root><interpro id=""><name>A<b>B</b>C</name><name>ignored</name><member_list><db_xref db="pfam" dbkey="P"/><db_xref db="pfam" dbkey="P"/></member_list><abstract>x<b>y</b>z<!--ignored-->q</abstract></interpro></root>"#.as_slice()).unwrap();
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].0, ("PFAM".into(), "P".into()));
        assert_eq!(entries[0].1.len(), 2);
        assert_eq!(entries[0].1[0].id, "UnknownIPR");
        assert_eq!(entries[0].1[0].name, "A");
        assert_eq!(entries[0].1[0].abstract_text, "x y zq");
    }

    #[test]
    fn namespaces_are_not_matched_by_local_name() {
        assert!(parse(br#"<root xmlns="u"><interpro id="I"><member_list><db_xref db="p" dbkey="q"/></member_list></interpro></root>"#.as_slice()).is_err());
    }

    #[test]
    fn missing_id_is_not_an_empty_id() {
        assert!(
            parse(
                br#"<interpro><member_list><db_xref db="p" dbkey="q"/></member_list></interpro>"#
                    .as_slice()
            )
            .is_err()
        );
    }
}
