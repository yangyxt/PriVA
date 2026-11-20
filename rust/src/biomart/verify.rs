use super::*;
use crate::python_path::is_file;
use std::cell::RefCell;
use std::rc::Rc;

#[derive(Clone)]
pub(super) struct Node {
    pub(super) name: String,
    pub(super) attributes: Json,
    pub(super) children: Vec<Node>,
    pub(super) text: String,
}
impl Node {
    fn copy(value: &crate::interpro_mapping_builder::Element) -> Self {
        Self {
            name: value.name.clone(),
            attributes: Json::Object(
                value
                    .ordered_attributes
                    .iter()
                    .map(|(k, v)| (k.clone(), string(v)))
                    .collect(),
            ),
            children: value
                .children
                .iter()
                .map(|n| Self::copy(&n.borrow()))
                .collect(),
            text: value.text.clone(),
        }
    }
    pub(super) fn walk<'a>(&'a self, name: &str, result: &mut Vec<&'a Node>) {
        if self.name == name {
            result.push(self)
        }
        for child in &self.children {
            child.walk(name, result)
        }
    }
}
pub(super) fn document(raw: &[u8], unicode_input: bool) -> Result<Node> {
    let result = Rc::new(RefCell::new(None));
    let sink = result.clone();
    crate::interpro_mapping_builder::visit_document(raw, unicode_input, move |element| {
        *sink.borrow_mut() = Some(Node::copy(element));
    })?;
    Rc::try_unwrap(result)
        .ok()
        .and_then(|v| v.into_inner())
        .ok_or_else(|| error("RuntimeError", "XML document produced no root"))
}
pub(super) fn parse_metadata(registry: &str, datasets: &str, source_id: &Json) -> Result<Json> {
    let selected = source(source_id)?;
    let root = document(registry.as_bytes(), true)?;
    let mut nodes = Vec::new();
    root.walk("MartURLLocation", &mut nodes);
    let locations = Json::Array(
        nodes
            .into_iter()
            .filter(|node| get(&node.attributes, "name").as_str() == Some("ENSEMBL_MART_ENSEMBL"))
            .map(|n| n.attributes.clone())
            .collect(),
    );
    let release = py_str(get(selected, "mart_release"));
    let host = match get(selected, "endpoint").as_str() {
        Some("https://grch37.ensembl.org/biomart/martservice") => "grch37.ensembl.org",
        Some("https://oct2024.archive.ensembl.org/biomart/martservice") => {
            "oct2024.archive.ensembl.org"
        }
        _ => return Err(error("KeyError", repr("endpoint"))),
    };
    let required = object([
        ("database", string(format!("ensembl_mart_{release}"))),
        ("displayName", string(format!("Ensembl Genes {release}"))),
        ("host", string(host)),
        ("path", string("/biomart/martservice")),
        ("serverVirtualSchema", string("default")),
    ]);
    let list = array(&locations);
    require(
        list.len() == 1
            && required
                .as_object()
                .unwrap()
                .iter()
                .all(|(k, v)| equal(get(&list[0], k), v)),
        format!(
            "Pinned BioMart source changed; expected {}, observed {}",
            py_repr(&required),
            py_repr(&locations)
        ),
    )?;
    // str.splitlines includes Unicode line boundaries; a terminal separator
    // does not add another line (empty trailing fields are otherwise retained).
    let lines = datasets.replace("\r\n", "\n");
    let records = Json::Array(
        lines
            .split([
                '\n', '\r', '\u{b}', '\u{c}', '\u{1c}', '\u{1d}', '\u{1e}', '\u{85}', '\u{2028}',
                '\u{2029}',
            ])
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .filter(|v| v.len() >= 5 && v[1] == DATASET)
            .map(|v| Json::Array(v.into_iter().map(string).collect()))
            .collect(),
    );
    let list = array(&records);
    let assembly = text(get(selected, "assembly"));
    require(
        list.len() == 1
            && array(&list[0])[2].as_str() == Some(&format!("Human genes ({assembly})"))
            && array(&list[0])[4].as_str() == Some(assembly),
        format!(
            "Pinned BioMart dataset/assembly changed: {}",
            py_repr(&records)
        ),
    )?;
    Ok(object([
        ("registry_record", array(&locations)[0].clone()),
        ("dataset_record", list[0].clone()),
        ("registry_xml", string(registry)),
        ("datasets_tsv", string(datasets)),
    ]))
}
pub(super) fn pinned(mapping: &Path, assembly: Option<&str>) -> Result<Json> {
    let record = read_json(&provenance_path(mapping))?;
    require(
        matches!(record, Json::Object(_)),
        format!(
            "Mapping provenance must be a JSON object: {}",
            mapping.display()
        ),
    )?;
    let selected = source(get(&record, "source_id"))?;
    require(
        source_for_basename(&name(mapping)).is_none()
            || get(selected, "basename").as_str() == Some(&name(mapping)),
        format!(
            "Known BioMart filename conflicts with its pinned source provenance: {}",
            mapping.display()
        ),
    )?;
    if get(selected, "source_id").as_str() == Some(FTP) {
        return ftp_record(mapping, record, assembly);
    }
    let mut expected = object([
        ("source_id", get(selected, "source_id").clone()),
        ("endpoint", get(selected, "endpoint").clone()),
        ("mart", string("ENSEMBL_MART_ENSEMBL")),
        ("mart_release", get(selected, "mart_release").clone()),
        (
            "annotation_base_release",
            get(selected, "annotation_base_release").clone(),
        ),
        ("assembly", get(selected, "assembly").clone()),
        ("dataset", string(DATASET)),
        ("query_xml", get(facts(), "query_xml").clone()),
    ]);
    expected.set("attributes", get(facts(), "attributes").clone());
    expected.set("header", get(facts(), "headers").clone());
    require(
        expected
            .as_object()
            .unwrap()
            .iter()
            .all(|(k, v)| equal(get(&record, k), v)),
        format!(
            "Mapping provenance does not match the explicitly pinned {} source: {}",
            text(get(selected, "source_id")),
            mapping.display()
        ),
    )?;
    if let Some(assembly) = assembly {
        require(
            array(get(selected, "aliases"))
                .iter()
                .any(|a| a.as_str() == Some(assembly)),
            format!(
                "Pinned {} BioMart mapping cannot be selected implicitly for {assembly}; supply an assembly-compatible mapping",
                text(get(selected, "assembly"))
            ),
        )?;
    }
    require(
        truth(get(&record, "completion_stamp_verified")),
        format!(
            "Mapping provenance lacks a completed-export receipt: {}",
            mapping.display()
        ),
    )?;
    let before = get(&record, "metadata_before");
    let after = get(&record, "metadata_after");
    for metadata in [before, after] {
        require(
            matches!(metadata, Json::Object(_))
                && ["registry_xml", "datasets_tsv"]
                    .iter()
                    .all(|k| get(metadata, k).as_str().is_some()),
            format!(
                "Mapping provenance metadata must contain publisher XML/TSV strings: {}",
                mapping.display()
            ),
        )?;
        let parsed = parse_metadata(
            text(get(metadata, "registry_xml")),
            text(get(metadata, "datasets_tsv")),
            get(selected, "source_id"),
        )?;
        require(
            ["registry_record", "dataset_record"]
                .iter()
                .all(|k| equal(get(metadata, k), get(&parsed, k))),
            format!(
                "Mapping provenance metadata is inconsistent: {}",
                mapping.display()
            ),
        )?;
    }
    require(
        ["registry_record", "dataset_record"]
            .iter()
            .all(|k| equal(get(before, k), get(after, k))),
        format!(
            "Mapping source changed during acquisition: {}",
            mapping.display()
        ),
    )?;
    require(
        get(&record, "query_sha256").as_str()
            == Some(&hash(text(get(facts(), "query_xml")).as_bytes())),
        format!("Mapping query identity mismatch: {}", mapping.display()),
    )?;
    let actual = snapshot_identity(mapping, get(selected, "source_id"))?;
    require(
        equal(
            &actual,
            &object([
                ("bytes", get(&record, "bytes").clone()),
                ("sha256", get(&record, "sha256").clone()),
            ]),
        ),
        format!(
            "Pinned BioMart mapping size/SHA256 mismatch: {}; retain it and recover explicitly",
            mapping.display()
        ),
    )?;
    let observed = table::validate(mapping, false, None)?;
    let stats = get(&record, "statistics");
    require(
        exact_keys(
            stats,
            &observed
                .as_object()
                .unwrap()
                .iter()
                .map(|(k, _)| k.as_str())
                .collect::<Vec<_>>(),
        ),
        format!(
            "Mapping provenance is missing full-export statistics: {}",
            mapping.display()
        ),
    )?;
    require(
        stats
            .as_object()
            .unwrap()
            .iter()
            .all(|(k, v)| k == "completion_stamp_verified" || equal(get(&observed, k), v)),
        format!(
            "Pinned BioMart mapping statistics mismatch: {}",
            mapping.display()
        ),
    )?;
    Ok(record)
}
pub(super) fn ftp_table_identity(mapping: &Path) -> Result<Json> {
    let selected = source(&string(FTP))?;
    require(
        is_file(mapping)?
            && equal(
                &Json::Int(metadata(mapping)?.len() as i128),
                get(selected, "expected_bytes"),
            ),
        "Pinned FTP table size differs; existing bytes retained",
    )?;
    require(
        sha256(mapping)? == text(get(selected, "expected_sha256")),
        "Pinned FTP table SHA256 differs; existing bytes retained",
    )?;
    let observed = table::validate(mapping, false, None)?;
    require(
        equal(&observed, get(facts(), "ftp_statistics")),
        "Pinned FTP table schema/order/multiplicity statistics differ",
    )?;
    require(
        sha256(mapping)? == text(get(selected, "expected_sha256")),
        "FTP mapping changed during verification",
    )?;
    Ok(observed)
}
pub(super) fn ftp_record(mapping: &Path, record: Json, assembly: Option<&str>) -> Result<Json> {
    let selected = source(&string(FTP))?;
    if let Some(assembly) = assembly {
        require(
            array(get(selected, "aliases"))
                .iter()
                .any(|v| v.as_str() == Some(assembly)),
            format!("FTP GRCh38 mapping cannot be used for {assembly}"),
        )?;
    }
    require(
        matches!(record, Json::Object(_)),
        "FTP provenance must be a JSON object",
    )?;
    for key in [
        "scientific_qualified",
        "production_selected",
        "online_query_executed",
        "original_biomart_export_identity",
        "original_biomart_response_order_identity",
    ] {
        require(
            matches!(get(&record, key), Json::Bool(false)),
            "FTP resource provenance must retain explicit unqualified/offline status",
        )?;
    }
    require(
        get(facts(), "ftp_specification")
            .as_object()
            .unwrap()
            .iter()
            .all(|(k, v)| equal(get(&record, k), v)),
        "FTP source/release/assembly/method provenance differs",
    )?;
    require(
        ![
            "query_xml",
            "metadata_before",
            "metadata_after",
            "registry_xml",
            "completion_stamp_verified",
        ]
        .iter()
        .any(|k| record.get(k).is_some()),
        "FTP provenance must not claim an online BioMart response or completion stamp",
    )?;
    require(
        equal(get(&record, "bytes"), get(selected, "expected_bytes"))
            && equal(get(&record, "sha256"), get(selected, "expected_sha256")),
        "FTP provenance table identity differs",
    )?;
    if get(&record, "metadata_kind").as_str() == Some("priva_bundled_mapping_snapshot") {
        require(
            matches!(get(&record, "metadata_schema_version"), Json::Int(1)),
            "Unknown bundled FTP metadata schema",
        )?;
    } else {
        require(
            matches!(
                get(&record, "operation").as_str(),
                Some("registered_existing_exact_snapshot" | "projected_from_verified_ftp_sources")
            ),
            "Unknown FTP resource operation",
        )?;
        require(
            hex64(get(&record, "registration_producer_sha256")),
            "Missing FTP registration producer identity",
        )?;
    }
    require(
        equal(get(&record, "statistics"), &ftp_table_identity(mapping)?),
        "FTP provenance statistics differ",
    )?;
    Ok(record)
}

pub(super) fn parse_xml(blob: &[u8], dataset_row: &Json) -> Result<Json> {
    let row = match dataset_row {
        Json::Array(v) => v,
        _ => return Err(error("TypeError", "BioMart dataset row must be a list")),
    };
    let identity = row.get(1).is_some_and(|v| v.as_str() == Some(DATASET))
        && row
            .get(2)
            .is_some_and(|v| v.as_str() == Some("Human genes (GRCh38.p14)"));
    // Preserve Python's short-circuit before the fixed index access.
    let same = if identity {
        row.get(6)
            .ok_or_else(|| error("IndexError", "list index out of range"))?
            .as_str()
            == Some("GRCh38.p14")
    } else {
        false
    };
    require(same, "Publisher dataset/assembly metadata changed")?;
    let root = document(blob, false)?;
    require(
        root.name == "DatasetConfig"
            && get(&root.attributes, "dataset").as_str() == Some(DATASET)
            && get(&root.attributes, "version").as_str() == Some("GRCh38.p14")
            && equal(get(&root.attributes, "datasetID"), &row[0]),
        "Human XML identity differs from dataset metadata",
    )?;
    require(
        equal(get(&root.attributes, "displayName"), &row[2])
            && get(&root.attributes, "interfaces").as_str() == Some("default"),
        "Human dataset display/interface changed",
    )?;
    let attrs = [
        (
            "ensembl_gene_id",
            "gene_id_1020_key",
            "main",
            "stable_id_1023",
        ),
        (
            "ensembl_transcript_id",
            "transcript_id_1064_key",
            "main",
            "stable_id_1066",
        ),
        (
            "hgnc_symbol",
            "gene_id_1020_key",
            "hsapiens_gene_ensembl__ox_hgnc__dm",
            "display_label_1074",
        ),
        (
            "uniprotswissprot",
            "translation_id_1068_key",
            "hsapiens_gene_ensembl__ox_uniprotswissprot__dm",
            "dbprimary_acc_1074",
        ),
        (
            "refseq_mrna",
            "transcript_id_1064_key",
            "hsapiens_gene_ensembl__ox_refseq_mrna__dm",
            "dbprimary_acc_1074",
        ),
    ];
    let mut nodes = Vec::new();
    root.walk("AttributeDescription", &mut nodes);
    let mut records = Json::Object(vec![]);
    for (i, (name, key, table, field)) in attrs.iter().enumerate() {
        let selected = nodes
            .iter()
            .filter(|n| get(&n.attributes, "internalName").as_str() == Some(name))
            .collect::<Vec<_>>();
        require(
            selected.len() == 1,
            format!("Requested attribute is not unique: {name}"),
        )?;
        let node = selected[0];
        let expected = [
            string(*key),
            string(*table),
            string(*field),
            array(get(facts(), "headers"))[i].clone(),
        ];
        require(
            ["key", "tableConstraint", "field", "displayName"]
                .iter()
                .zip(expected)
                .all(|(k, v)| equal(get(&node.attributes, k), &v)),
            format!("Requested attribute definition changed: {name}"),
        )?;
        records.set(name, node.attributes.clone());
    }
    let children = |name: &str| {
        Json::Array(
            root.children
                .iter()
                .filter(|n| n.name == name)
                .map(|n| {
                    if n.text.is_empty() {
                        Json::Null
                    } else {
                        string(&n.text)
                    }
                })
                .collect(),
        )
    };
    let mains = children("MainTable");
    let keys = children("Key");
    require(
        equal(
            &mains,
            &Json::Array(
                ["gene", "transcript", "translation"]
                    .map(|level| string(format!("{DATASET}__{level}__main")))
                    .into(),
            ),
        ),
        "Unexpected main table hierarchy",
    )?;
    require(
        equal(
            &keys,
            &Json::Array(
                [
                    "gene_id_1020_key",
                    "transcript_id_1064_key",
                    "translation_id_1068_key",
                ]
                .map(string)
                .into(),
            ),
        ),
        "Unexpected main table key hierarchy",
    )?;
    Ok(object([
        ("dataset_row", dataset_row.clone()),
        ("root", root.attributes),
        ("attributes", records),
        ("main_tables", mains),
        ("main_keys", keys),
        ("xml_sha256", string(hash(blob))),
    ]))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn document_keeps_attribute_order_and_direct_text() {
        let root = document(
            b"<r z='1' a='2'><child n='x'>before<i/>after</child></r>",
            false,
        )
        .unwrap();
        assert_eq!(py_repr(&root.attributes), "{'z': '1', 'a': '2'}");
        assert_eq!(root.children[0].text, "before");
    }
    #[test]
    fn unicode_feed_ignores_declared_byte_encoding() {
        let raw = "<?xml version='1.0' encoding='utf-16'?><r a='α'/>";
        let root = document(raw.as_bytes(), true).unwrap();
        assert_eq!(get(&root.attributes, "a").as_str(), Some("α"));
        assert!(document(raw.as_bytes(), false).is_err());
    }
    #[test]
    fn complete_document_errors_precede_metadata_validation() {
        let error = parse_metadata("<r/>junk", "", &string(GRCH37))
            .err()
            .unwrap()
            .to_string();
        assert!(
            error.starts_with("ParseError: junk after document element:"),
            "{error}"
        );
    }
    #[test]
    fn namespace_expansion_does_not_match_unqualified_registry_tag() {
        let root = document(b"<r xmlns:p='u'><p:MartURLLocation/></r>", true).unwrap();
        assert_eq!(root.children[0].name, "{u}MartURLLocation");
    }
}
