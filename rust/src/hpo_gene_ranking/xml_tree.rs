//! A small element tree read with quick-xml, enough for hp.owl (RDF/XML) and
//! Orphadata's product 6 XML: element names resolved to namespace and local
//! name, attributes with entity references resolved, and the text before an
//! element's first child (ElementTree's `.text`). Top-level elements (the
//! children of the root) are handed over one at a time, so a large file never
//! sits in memory as one tree.

use crate::hpo_gene_ranking::Result;
use quick_xml::XmlVersion;
use quick_xml::escape::resolve_xml_entity;
use quick_xml::events::Event;
use std::collections::HashMap;

pub const XML_NAMESPACE: &str = "http://www.w3.org/XML/1998/namespace";

#[derive(Debug, Default)]
pub struct Attribute {
    pub namespace: String,
    pub local: String,
    pub value: String,
}

#[derive(Debug, Default)]
pub struct Element {
    pub namespace: String,
    pub local: String,
    pub attributes: Vec<Attribute>,
    pub children: Vec<Element>,
    /// Text before the first child element.
    pub text: String,
    started_child: bool,
}

impl Element {
    pub fn is(&self, namespace: &str, local: &str) -> bool {
        self.namespace == namespace && self.local == local
    }

    pub fn attribute(&self, namespace: &str, local: &str) -> Option<&str> {
        self.attributes
            .iter()
            .find(|attribute| attribute.namespace == namespace && attribute.local == local)
            .map(|attribute| attribute.value.as_str())
    }

    /// `element.find(tag)`: the first child with that tag.
    pub fn find(&self, tag: &str) -> Option<&Element> {
        self.children.iter().find(|child| child.local == tag)
    }

    /// `element.findtext("a/b", default="")`: the text of the first b child
    /// of an a child, or "" when there is none.
    pub fn findtext(&self, path: &str) -> String {
        let mut current: Vec<&Element> = vec![self];
        for step in path.split('/') {
            current = current
                .iter()
                .flat_map(|element| element.children.iter().filter(|child| child.local == step))
                .collect();
        }
        current
            .first()
            .map(|element| element.text.clone())
            .unwrap_or_default()
    }

    /// `element.iter(tag)`: the element itself and every descendant with that
    /// tag, in document order.
    pub fn iter_tag<'a>(&'a self, tag: &str, found: &mut Vec<&'a Element>) {
        if self.local == tag {
            found.push(self);
        }
        for child in &self.children {
            child.iter_tag(tag, found);
        }
    }
}

/// The text of an XML file: Latin-1 when its declaration says so, else UTF-8.
fn decode(bytes: &[u8], path: &str) -> Result<String> {
    let bytes = bytes.strip_prefix(b"\xEF\xBB\xBF").unwrap_or(bytes);
    let head: String = bytes.iter().take(200).map(|&byte| byte as char).collect();
    let declared = head
        .split_once("encoding=")
        .map(|(_, rest)| {
            let quote = rest.chars().next().unwrap_or('"');
            rest[1..]
                .split(quote)
                .next()
                .unwrap_or("")
                .to_ascii_lowercase()
        })
        .unwrap_or_default();
    if matches!(
        declared.as_str(),
        "iso-8859-1" | "iso8859-1" | "latin-1" | "latin1"
    ) {
        Ok(bytes.iter().map(|&byte| byte as char).collect())
    } else {
        String::from_utf8(bytes.to_vec()).map_err(|error| format!("{path}: not UTF-8: {error}"))
    }
}

/// `<!ENTITY name "value">` declarations of a DOCTYPE.
fn doctype_entities(doctype: &str, entities: &mut HashMap<String, String>) {
    let mut rest = doctype;
    while let Some(position) = rest.find("<!ENTITY") {
        rest = &rest[position + "<!ENTITY".len()..];
        let body = rest.trim_start();
        let name: String = body.chars().take_while(|c| !c.is_whitespace()).collect();
        let after = body[name.len()..].trim_start();
        if let Some(quote) = after.chars().next().filter(|c| *c == '"' || *c == '\'')
            && let Some(end) = after[1..].find(quote)
        {
            entities.insert(name, after[1..1 + end].to_owned());
        }
    }
}

/// The value of a predefined or DOCTYPE-declared entity.
fn resolve_entity<'a>(name: &str, entities: &'a HashMap<String, String>) -> Option<&'a str> {
    if let Some(value) = resolve_xml_entity(name) {
        return Some(value);
    }
    entities.get(name).map(String::as_str)
}

fn split_name(name: &str) -> (&str, &str) {
    match name.split_once(':') {
        Some((prefix, local)) => (prefix, local),
        None => ("", name),
    }
}

/// Reads an XML file, calls `handle` with each child element of the root
/// (with its whole subtree) and returns the root without children.
pub fn read_top_level(
    path: &str,
    mut handle: impl FnMut(Element) -> Result<()>,
) -> Result<Element> {
    let bytes = crate::hpo_gene_ranking::tsv::read_bytes(path)?;
    let text = decode(&bytes, path)?;
    drop(bytes);
    let mut reader = quick_xml::Reader::from_str(&text);
    reader.config_mut().expand_empty_elements = true;
    let mut entities: HashMap<String, String> = HashMap::new();
    let mut scopes: Vec<HashMap<String, String>> = vec![HashMap::from([(
        "xml".to_owned(),
        XML_NAMESPACE.to_owned(),
    )])];
    let mut stack: Vec<Element> = Vec::new();
    let mut root: Option<Element> = None;
    let failure = |what: String| format!("{path}: {what}");
    loop {
        let event = reader
            .read_event()
            .map_err(|error| failure(error.to_string()))?;
        match event {
            Event::DocType(doctype) => doctype_entities(&doctype, &mut entities),
            Event::Start(start) => {
                let mut declared: HashMap<String, String> = HashMap::new();
                let mut raw: Vec<(String, String)> = Vec::new();
                for attribute in start.attributes() {
                    let attribute = attribute.map_err(|error| failure(error.to_string()))?;
                    let key = AsRef::<str>::as_ref(&attribute.key).to_owned();
                    let value = attribute
                        .normalized_value_with(XmlVersion::Implicit1_0, 128, |name| {
                            resolve_entity(name, &entities)
                        })
                        .map_err(|error| failure(error.to_string()))?
                        .into_owned();
                    if key == "xmlns" {
                        declared.insert(String::new(), value);
                    } else if let Some(prefix) = key.strip_prefix("xmlns:") {
                        declared.insert(prefix.to_owned(), value);
                    } else {
                        raw.push((key, value));
                    }
                }
                let mut scope = scopes.last().cloned().unwrap_or_default();
                scope.extend(declared);
                let qualified = start.name();
                let name = AsRef::<str>::as_ref(&qualified).to_owned();
                let (prefix, local) = split_name(&name);
                let namespace = scope.get(prefix).cloned().unwrap_or_default();
                let attributes = raw
                    .into_iter()
                    .map(|(key, value)| {
                        let (prefix, local) = split_name(&key);
                        // Attributes without a prefix have no namespace.
                        let namespace = if prefix.is_empty() {
                            String::new()
                        } else {
                            scope.get(prefix).cloned().unwrap_or_default()
                        };
                        Attribute {
                            namespace,
                            local: local.to_owned(),
                            value,
                        }
                    })
                    .collect();
                scopes.push(scope);
                if let Some(parent) = stack.last_mut() {
                    parent.started_child = true;
                }
                stack.push(Element {
                    namespace,
                    local: local.to_owned(),
                    attributes,
                    ..Element::default()
                });
            }
            Event::End(_) => {
                scopes.pop();
                let element = stack
                    .pop()
                    .ok_or_else(|| failure("unbalanced end tag".to_owned()))?;
                match stack.len() {
                    0 => root = Some(element),
                    1 => handle(element)?,
                    _ => stack.last_mut().expect("a parent").children.push(element),
                }
            }
            Event::Text(content) => {
                if let Some(top) = stack.last_mut()
                    && !top.started_child
                {
                    top.text.push_str(&content.xml10_content());
                }
            }
            Event::CData(content) => {
                if let Some(top) = stack.last_mut()
                    && !top.started_child
                {
                    top.text.push_str(&content.xml10_content());
                }
            }
            Event::GeneralRef(reference) => {
                if let Some(top) = stack.last_mut()
                    && !top.started_child
                {
                    match reference
                        .resolve_char_ref()
                        .map_err(|error| failure(error.to_string()))?
                    {
                        Some(character) => top.text.push(character),
                        None => {
                            let name: &str = &reference;
                            match resolve_entity(name, &entities) {
                                Some(value) => top.text.push_str(value),
                                None => return Err(failure(format!("unknown entity &{name};"))),
                            }
                        }
                    }
                }
            }
            Event::Eof => break,
            Event::Empty(_) | Event::Comment(_) | Event::Decl(_) | Event::PI(_) => {}
        }
    }
    root.ok_or_else(|| failure("no root element".to_owned()))
}
