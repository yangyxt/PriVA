//! `build-hpo-graph`: the nodes and edges tables that
//! `scripts/build_hpo_owl_graph_tool.py --hpo-only` writes from hp.owl, read
//! directly from the RDF/XML file instead of through owlready2.
//!
//! Left out, because they are owlready2 renderings that no later step reads:
//! the graph-tool `.gt` file, the annotation, logical-axiom and prefix
//! sidecar tables, and the node columns `logical_definition_text`,
//! `annotation_property_count` and `annotation_value_count`.

use crate::hpo_gene_ranking::python_text::with_suffix;
use crate::hpo_gene_ranking::tsv::{Cell, Record, write_records};
use crate::hpo_gene_ranking::xml_tree::{Element, read_top_level};
use crate::hpo_gene_ranking::{Options, Result};
use std::collections::{BTreeMap, HashMap};

const RDF: &str = "http://www.w3.org/1999/02/22-rdf-syntax-ns#";
const RDFS: &str = "http://www.w3.org/2000/01/rdf-schema#";
const OWL: &str = "http://www.w3.org/2002/07/owl#";
const OBO: &str = "http://purl.obolibrary.org/obo/";
const OBO_IN_OWL: &str = "http://www.geneontology.org/formats/oboInOwl#";

const HELP: &str = "\
hpo-gene-ranking build-hpo-graph: HPO nodes and edges tables from hp.owl

  --hpo-owl PATH          hp.owl (RDF/XML)
  --output-prefix PREFIX  writes PREFIX.nodes.tsv and PREFIX.edges.tsv
  --hpo-only              accepted for compatibility; only HP classes are written";

/// `curie_from_iri`.
pub fn curie_from_iri(iri: &str) -> String {
    if let Some(tail) = iri.strip_prefix(OBO) {
        return match tail.split_once('_') {
            Some((prefix, local)) => format!("{prefix}:{local}"),
            None => tail.to_owned(),
        };
    }
    let known = [
        (RDFS, "rdfs"),
        (OWL, "owl"),
        (RDF, "rdf"),
        (OBO_IN_OWL, "oboInOwl"),
        ("http://purl.org/dc/elements/1.1/", "dc"),
        ("http://purl.org/dc/terms/", "dcterms"),
        ("http://www.w3.org/2004/02/skos/core#", "skos"),
    ];
    for (prefix, label) in known {
        if let Some(rest) = iri.strip_prefix(prefix) {
            return format!("{label}:{rest}");
        }
    }
    // owlready2's entity name: the part after the last '#', else after the last '/'.
    match iri.rfind('#') {
        Some(position) => iri[position + 1..].to_owned(),
        None => iri.rsplit('/').next().unwrap_or(iri).to_owned(),
    }
}

/// A class expression of an RDF/XML axiom, as far as the tables need it.
enum Expression {
    Named(String),
    /// A blank node written elsewhere in the file (`rdf:nodeID`), resolved after reading.
    NodeRef(String),
    Restriction(Option<Box<Expression>>),
    And(Vec<Expression>),
    Or(Vec<Expression>),
    Other,
}

/// One annotation value: a literal text or another entity.
enum Value {
    Literal { text: String, datatype: String },
    Resource(String),
}

#[derive(Default)]
struct ClassData {
    annotations: Vec<(String, Value)>,
    /// Named superclasses in document order.
    named_parents: Vec<String>,
    /// Anonymous superclass expressions.
    superclass_expressions: Vec<Expression>,
    /// Equivalent classes, named or anonymous.
    equivalents: Vec<Expression>,
}

fn element_iri(element: &Element) -> String {
    format!("{}{}", element.namespace, element.local)
}

/// An expression given by an `rdf:resource` attribute or the element's first child.
fn expression_of(holder: &Element) -> Expression {
    if let Some(iri) = holder.attribute(RDF, "resource") {
        return Expression::Named(iri.to_owned());
    }
    if let Some(node) = holder.attribute(RDF, "nodeID") {
        return Expression::NodeRef(node.to_owned());
    }
    match holder.children.first() {
        Some(child) => parse_expression(child),
        None => Expression::Other,
    }
}

fn parse_expression(element: &Element) -> Expression {
    if let Some(iri) = element.attribute(RDF, "about") {
        return Expression::Named(iri.to_owned());
    }
    // A reference to a blank node, not its definition (which has children).
    if element.children.is_empty()
        && let Some(node) = element.attribute(RDF, "nodeID")
    {
        return Expression::NodeRef(node.to_owned());
    }
    if element.is(OWL, "Restriction") {
        let mut filler = None;
        for child in &element.children {
            if child.namespace != OWL {
                continue;
            }
            match child.local.as_str() {
                "someValuesFrom" | "allValuesFrom" | "hasValue" | "onClass" | "onDataRange" => {
                    filler = Some(Box::new(expression_of(child)));
                }
                _ => {}
            }
        }
        return Expression::Restriction(filler);
    }
    if element.is(OWL, "Class") || element.is(RDF, "Description") {
        for child in &element.children {
            if child.is(OWL, "intersectionOf") {
                return Expression::And(child.children.iter().map(parse_expression).collect());
            }
            if child.is(OWL, "unionOf") {
                return Expression::Or(child.children.iter().map(parse_expression).collect());
            }
        }
    }
    Expression::Other
}

/// `iter_restrictions`: the named filler of every restriction, recursing into
/// restriction fillers and the parts of intersections and unions.
fn restriction_fillers(expression: &Expression, fillers: &mut Vec<String>) {
    match expression {
        Expression::Restriction(filler) => {
            if let Some(filler) = filler {
                if let Expression::Named(iri) = filler.as_ref() {
                    fillers.push(iri.clone());
                }
                restriction_fillers(filler, fillers);
            }
        }
        Expression::And(parts) | Expression::Or(parts) => {
            for part in parts {
                restriction_fillers(part, fillers);
            }
        }
        Expression::Named(_) | Expression::NodeRef(_) | Expression::Other => {}
    }
}

/// Replaces every blank-node reference by the expression its definition holds.
fn resolve(
    expression: Expression,
    blank_nodes: &HashMap<String, Element>,
    depth: usize,
) -> Expression {
    match expression {
        Expression::NodeRef(node) => match blank_nodes.get(&node) {
            Some(element) if depth < 64 => {
                resolve(parse_expression(element), blank_nodes, depth + 1)
            }
            _ => Expression::Other,
        },
        Expression::Restriction(Some(filler)) => {
            Expression::Restriction(Some(Box::new(resolve(*filler, blank_nodes, depth + 1))))
        }
        Expression::And(parts) => Expression::And(
            parts
                .into_iter()
                .map(|part| resolve(part, blank_nodes, depth + 1))
                .collect(),
        ),
        Expression::Or(parts) => Expression::Or(
            parts
                .into_iter()
                .map(|part| resolve(part, blank_nodes, depth + 1))
                .collect(),
        ),
        other => other,
    }
}

/// Reads hp.owl: data of every class and the first label of every entity.
fn read_owl(path: &str) -> Result<(BTreeMap<String, ClassData>, HashMap<String, String>)> {
    let mut classes: BTreeMap<String, ClassData> = BTreeMap::new();
    let mut labels: HashMap<String, String> = HashMap::new();
    let mut blank_nodes: HashMap<String, Element> = HashMap::new();
    read_top_level(path, |element| {
        let Some(about) = element.attribute(RDF, "about").map(str::to_owned) else {
            if let Some(node) = element.attribute(RDF, "nodeID").map(str::to_owned) {
                blank_nodes.insert(node, element);
            }
            return Ok(());
        };
        if element.is(OWL, "Axiom") {
            return Ok(());
        }
        for child in &element.children {
            if child.is(RDFS, "label") && child.attribute(RDF, "resource").is_none() {
                labels
                    .entry(about.clone())
                    .or_insert_with(|| child.text.clone());
            }
        }
        if !element.is(OWL, "Class") {
            return Ok(());
        }
        let data = classes.entry(about.clone()).or_default();
        for child in &element.children {
            if child.is(RDFS, "subClassOf") {
                match expression_of(child) {
                    Expression::Named(iri) => data.named_parents.push(iri),
                    expression => data.superclass_expressions.push(expression),
                }
            } else if child.is(OWL, "equivalentClass") {
                data.equivalents.push(expression_of(child));
            } else if !child.is(RDF, "type") && !child.is(OWL, "disjointWith") {
                let property = element_iri(child);
                let value = match child.attribute(RDF, "resource") {
                    Some(iri) => Value::Resource(iri.to_owned()),
                    None => Value::Literal {
                        text: child.text.clone(),
                        datatype: child.attribute(RDF, "datatype").unwrap_or("").to_owned(),
                    },
                };
                data.annotations.push((property, value));
            }
        }
        Ok(())
    })?;
    for data in classes.values_mut() {
        let superclasses = std::mem::take(&mut data.superclass_expressions);
        for expression in superclasses {
            match resolve(expression, &blank_nodes, 0) {
                Expression::Named(iri) => data.named_parents.push(iri),
                resolved => data.superclass_expressions.push(resolved),
            }
        }
        let equivalents = std::mem::take(&mut data.equivalents);
        data.equivalents = equivalents
            .into_iter()
            .map(|expression| resolve(expression, &blank_nodes, 0))
            .collect();
    }
    Ok((classes, labels))
}

/// `value_to_text`: an entity as "CURIE label", a literal as Python's `str`
/// of owlready2's value (booleans become True/False).
fn value_text(value: &Value, labels: &HashMap<String, String>) -> String {
    match value {
        Value::Resource(iri) => {
            let curie = curie_from_iri(iri);
            let label = labels.get(iri).map(String::as_str).unwrap_or("");
            crate::hpo_gene_ranking::python_text::strip(&format!("{curie} {label}")).to_owned()
        }
        Value::Literal { text, datatype } => {
            if datatype == "http://www.w3.org/2001/XMLSchema#boolean" {
                match text.trim() {
                    "true" | "1" => "True".to_owned(),
                    "false" | "0" => "False".to_owned(),
                    other => other.to_owned(),
                }
            } else {
                text.clone()
            }
        }
    }
}

pub fn run(arguments: &[String]) -> Result<()> {
    let Some(options) = Options::parse(
        arguments,
        &["hpo-owl", "output-prefix"],
        &[
            "hpo-only",
            "export-graphml",
            "no-annotation-json",
            "no-logical-json",
        ],
    )?
    else {
        println!("{HELP}");
        return Ok(());
    };
    let owl_path = options.required("hpo-owl")?;
    let prefix = options
        .value("output-prefix")
        .unwrap_or("results/hpo_owl_graph_tool/hp_owl");
    let (classes, labels) = read_owl(owl_path)?;

    // HP classes, sorted by CURIE (they all share the prefix "HP").
    let mut hp: Vec<(String, &String, &ClassData)> = classes
        .iter()
        .map(|(iri, data)| (curie_from_iri(iri), iri, data))
        .filter(|(curie, _, _)| curie.starts_with("HP:"))
        .collect();
    hp.sort_by(|a, b| a.0.cmp(&b.0));
    let indexed: std::collections::HashSet<&str> =
        hp.iter().map(|(curie, _, _)| curie.as_str()).collect();

    let core = |data: &ClassData, property: &str| -> String {
        data.annotations
            .iter()
            .filter(|(name, _)| name == property)
            .map(|(_, value)| value_text(value, &labels))
            .collect::<Vec<_>>()
            .join(" | ")
    };
    let mut child_counts: HashMap<String, i64> = HashMap::new();
    let mut parents_of: HashMap<&str, Vec<String>> = HashMap::new();
    for (curie, _, data) in &hp {
        let parents: Vec<String> = data
            .named_parents
            .iter()
            .map(|iri| curie_from_iri(iri))
            .filter(|parent| indexed.contains(parent.as_str()))
            .collect();
        for parent in &parents {
            *child_counts.entry(parent.clone()).or_default() += 1;
        }
        parents_of.insert(curie, parents);
    }
    let names: HashMap<&str, String> = hp
        .iter()
        .map(|(curie, iri, data)| {
            let label = core(data, &format!("{RDFS}label"));
            let label = if label.is_empty() {
                labels.get(*iri).cloned().unwrap_or_default()
            } else {
                label
            };
            (curie.as_str(), label)
        })
        .collect();

    let mut node_records: Vec<Record> = Vec::new();
    let mut edge_records: Vec<Record> = Vec::new();
    for (curie, iri, data) in &hp {
        // Fillers of equivalent-class axioms, then of top-level restriction superclasses.
        let mut fillers = Vec::new();
        for expression in &data.equivalents {
            restriction_fillers(expression, &mut fillers);
        }
        for expression in &data.superclass_expressions {
            if matches!(expression, Expression::Restriction(_)) {
                restriction_fillers(expression, &mut fillers);
            }
        }
        let mut by_prefix: BTreeMap<String, BTreeMap<String, String>> = BTreeMap::new();
        for filler in &fillers {
            let filler_curie = curie_from_iri(filler);
            let Some((filler_prefix, _)) = filler_curie.split_once(':') else {
                continue;
            };
            let label = labels.get(filler).map(String::as_str).unwrap_or("");
            let text =
                crate::hpo_gene_ranking::python_text::strip(&format!("{filler_curie} {label}"))
                    .to_owned();
            by_prefix
                .entry(filler_prefix.to_owned())
                .or_default()
                .insert(filler_curie.clone(), text);
        }
        let linked = |prefix: &str| -> String {
            let mut values: Vec<&String> = by_prefix
                .get(prefix)
                .map(|terms| terms.values().collect())
                .unwrap_or_default();
            values.sort();
            values
                .iter()
                .map(|value| value.as_str())
                .collect::<Vec<_>>()
                .join(" | ")
        };
        let mut other: Vec<&String> = by_prefix
            .iter()
            .filter(|(prefix, _)| !["GO", "UBERON", "PATO", "CL"].contains(&prefix.as_str()))
            .flat_map(|(_, terms)| terms.values())
            .collect();
        other.sort();
        let linked_count: usize = by_prefix.values().map(BTreeMap::len).sum();
        let restriction_superclasses = data
            .superclass_expressions
            .iter()
            .filter(|expression| matches!(expression, Expression::Restriction(_)))
            .count();
        let obsolete = core(data, &format!("{OWL}deprecated")).to_lowercase() == "true";
        let parents = &parents_of[curie.as_str()];
        node_records.push(vec![
            ("hp_id", Cell::from(curie.as_str())),
            ("iri", Cell::from(iri.as_str())),
            ("name", Cell::from(names[curie.as_str()].as_str())),
            (
                "definition",
                Cell::from(core(data, &format!("{OBO}IAO_0000115"))),
            ),
            ("comment", Cell::from(core(data, &format!("{RDFS}comment")))),
            (
                "alternative_ids",
                Cell::from(core(data, &format!("{OBO_IN_OWL}hasAlternativeId"))),
            ),
            (
                "exact_synonyms",
                Cell::from(core(data, &format!("{OBO_IN_OWL}hasExactSynonym"))),
            ),
            (
                "related_synonyms",
                Cell::from(core(data, &format!("{OBO_IN_OWL}hasRelatedSynonym"))),
            ),
            (
                "broad_synonyms",
                Cell::from(core(data, &format!("{OBO_IN_OWL}hasBroadSynonym"))),
            ),
            (
                "narrow_synonyms",
                Cell::from(core(data, &format!("{OBO_IN_OWL}hasNarrowSynonym"))),
            ),
            (
                "xrefs",
                Cell::from(core(data, &format!("{OBO_IN_OWL}hasDbXref"))),
            ),
            (
                "subsets",
                Cell::from(core(data, &format!("{OBO_IN_OWL}inSubset"))),
            ),
            (
                "replaced_by",
                Cell::from(core(data, &format!("{OBO}IAO_0100001"))),
            ),
            (
                "consider",
                Cell::from(core(data, &format!("{OBO_IN_OWL}consider"))),
            ),
            ("linked_go_terms", Cell::from(linked("GO"))),
            ("linked_uberon_terms", Cell::from(linked("UBERON"))),
            ("linked_pato_terms", Cell::from(linked("PATO"))),
            ("linked_cl_terms", Cell::from(linked("CL"))),
            (
                "linked_other_ontology_terms",
                Cell::from(
                    other
                        .iter()
                        .map(|value| value.as_str())
                        .collect::<Vec<_>>()
                        .join(" | "),
                ),
            ),
            (
                "obsolete",
                Cell::from(if obsolete { "true" } else { "false" }),
            ),
            (
                "has_logical_definition",
                Cell::from(if data.equivalents.is_empty() {
                    "false"
                } else {
                    "true"
                }),
            ),
            ("parent_count", Cell::from(parents.len())),
            (
                "child_count",
                Cell::Int(child_counts.get(curie.as_str()).copied().unwrap_or(0)),
            ),
            (
                "logical_definition_count",
                Cell::from(data.equivalents.len() + restriction_superclasses),
            ),
            ("linked_ontology_count", Cell::from(linked_count)),
        ]);
        for parent in parents {
            edge_records.push(vec![
                ("child_hp_id", Cell::from(curie.as_str())),
                ("child_name", Cell::from(names[curie.as_str()].as_str())),
                ("parent_hp_id", Cell::from(parent.as_str())),
                ("parent_name", Cell::from(names[parent.as_str()].as_str())),
                ("edge_class", Cell::from("is_a")),
                ("owl_axiom", Cell::from("SubClassOf")),
                (
                    "predicate_iri",
                    Cell::from("http://www.w3.org/2000/01/rdf-schema#subClassOf"),
                ),
                ("edge_source", Cell::from("owl_named_subclass")),
            ]);
        }
    }
    const NODE_FIELDS: [&str; 25] = [
        "hp_id",
        "iri",
        "name",
        "definition",
        "comment",
        "alternative_ids",
        "exact_synonyms",
        "related_synonyms",
        "broad_synonyms",
        "narrow_synonyms",
        "xrefs",
        "subsets",
        "replaced_by",
        "consider",
        "linked_go_terms",
        "linked_uberon_terms",
        "linked_pato_terms",
        "linked_cl_terms",
        "linked_other_ontology_terms",
        "obsolete",
        "has_logical_definition",
        "parent_count",
        "child_count",
        "logical_definition_count",
        "linked_ontology_count",
    ];
    const EDGE_FIELDS: [&str; 8] = [
        "child_hp_id",
        "child_name",
        "parent_hp_id",
        "parent_name",
        "edge_class",
        "owl_axiom",
        "predicate_iri",
        "edge_source",
    ];
    // The Python script writes these two tables with csv's default "\r\n" line ends.
    if !node_records.is_empty() {
        write_records(
            &with_suffix(prefix, ".nodes.tsv"),
            &NODE_FIELDS,
            &node_records,
            "\r\n",
        )?;
    }
    if !edge_records.is_empty() {
        write_records(
            &with_suffix(prefix, ".edges.tsv"),
            &EDGE_FIELDS,
            &edge_records,
            "\r\n",
        )?;
    }
    println!("owl={owl_path}");
    println!("classes={}", hp.len());
    println!("named_subclass_edges={}", edge_records.len());
    Ok(())
}
