//! Opt-in equivalent of `scripts/build_mondo_disease_scope.py`.
//! OBO xrefs/ancestry, inheritance precedence and manual review overrides are
//! preserved. Building a registry publishes it before optional annotation, as
//! the Python procedure does; the two outputs are not one transaction.
use crate::builder_table::{self as table, Publication};
use crate::cli::{Result, invalid};
use crate::hpo_assertions::COLUMNS as BASE_COLUMNS;
use crate::json_text::Json;
use crate::pandas_polars;
use crate::pandas_table::{Column, Table};
use crate::python_args::{self, Parsed};
use crate::python_text::python_strip;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs;
use std::io;
use std::path::Path;

const HEREDITARY: &str = "MONDO:0003847";
const NEOPLASM: &str = "MONDO:0005070";
const MENDELIAN: [&str; 9] = [
    "HP:0000006",
    "HP:0000007",
    "HP:0001417",
    "HP:0001419",
    "HP:0001423",
    "HP:0001427",
    "HP:0001450",
    "HP:0012275",
    "HP:0034341",
];
const NON_MONOGENIC: [&str; 4] = ["HP:0001426", "HP:0010982", "HP:0010983", "HP:0010984"];
const SCOPE_COLUMNS: [&str; 7] = [
    "mondo_id",
    "mondo_name",
    "disease_scope",
    "priva_scope",
    "scope_evidence",
    "scope_reference",
    "scope_review_status",
];
const REGISTRY_COLUMNS: [&str; 16] = [
    "disease_id",
    "mondo_id",
    "mondo_name",
    "mondo_version",
    "hpo_mendelian_inheritance",
    "hpo_non_monogenic_inheritance",
    "mondo_hereditary",
    "mondo_neoplasm",
    "mondo_predisposition",
    "mondo_omim_susceptibility",
    "disease_scope",
    "priva_scope",
    "scope_evidence",
    "scope_reference",
    "scope_review_status",
    "scope_note",
];
const DISEASE_SCOPES: [&str; 9] = [
    "mendelian_non_neoplastic",
    "mendelian_cancer_predisposition",
    "mendelian_susceptibility_uncertain",
    "complex_or_non_monogenic",
    "complex_susceptibility_candidate",
    "somatic_oncogenic",
    "neoplastic_uncertain",
    "conflicting_inheritance",
    "uncertain",
];

#[derive(Default)]
struct Term {
    id: String,
    name: String,
    parents: Vec<String>,
    xrefs: Vec<String>,
    subsets: Vec<String>,
    obsolete: bool,
}

fn finish(raw: Option<Term>, terms: &mut BTreeMap<String, Term>) {
    if let Some(term) = raw.filter(|term| term.id.starts_with("MONDO:")) {
        terms.insert(term.id.clone(), term);
    }
}

fn first_word(value: &str) -> Result<String> {
    value
        .split(|ch: char| ch.is_whitespace() || matches!(ch, '\u{1c}'..='\u{1f}'))
        .find(|word| !word.is_empty())
        .map(str::to_owned)
        .ok_or_else(|| {
            crate::builder_error::category("IndexError", invalid("list index out of range"))
        })
}

fn parse_obo(path: &Path) -> Result<(String, BTreeMap<String, Term>)> {
    let contents = String::from_utf8(fs::read(path)?)
        .map_err(|error| crate::builder_error::category("UnicodeDecodeError", error.into()))?;
    let mut version = String::new();
    let mut terms = BTreeMap::new();
    let mut raw: Option<Term> = None;
    // Python open(..., encoding='utf-8') accepts CR, LF and CRLF. Empty lines
    // between CR and LF have no OBO meaning in this deliberately narrow parser.
    for line in contents.split(['\r', '\n']) {
        if raw.is_none() {
            if let Some(value) = line.strip_prefix("data-version:") {
                version = python_strip(value).to_owned();
            }
        }
        if line == "[Term]" {
            finish(raw.take(), &mut terms);
            raw = Some(Term::default());
            continue;
        }
        if line.starts_with('[') {
            finish(raw.take(), &mut terms);
            continue;
        }
        let Some(term) = &mut raw else {
            continue;
        };
        if let Some(value) = line.strip_prefix("id: ") {
            term.id = python_strip(value).to_owned();
        } else if let Some(value) = line.strip_prefix("name: ") {
            term.name = python_strip(value).to_owned();
        } else if let Some(value) = line.strip_prefix("is_a: ") {
            term.parents.push(first_word(value)?);
        } else if let Some(value) = line.strip_prefix("subset: ") {
            term.subsets.push(first_word(value)?);
        } else if line == "is_obsolete: true" {
            term.obsolete = true;
        } else if let Some(value) = line.strip_prefix("xref: ") {
            let value = first_word(value)?;
            if !line.contains('{') || line.contains("MONDO:equivalentTo") {
                if value.starts_with("OMIM:") || value.starts_with("ORPHA:") {
                    term.xrefs.push(value);
                } else if let Some(value) = value.strip_prefix("Orphanet:") {
                    term.xrefs.push(format!("ORPHA:{value}"));
                }
            }
        }
    }
    finish(raw, &mut terms);
    if version.is_empty() {
        return Err(invalid(format!(
            "MONDO OBO has no data-version header: {}",
            path.display()
        )));
    }
    for root in [HEREDITARY, NEOPLASM] {
        if !terms.contains_key(root) {
            return Err(invalid(format!(
                "MONDO OBO is missing required term {root}: {}",
                path.display()
            )));
        }
    }
    Ok((version, terms))
}

fn ancestors(
    id: &str,
    terms: &BTreeMap<String, Term>,
    cache: &mut HashMap<String, HashSet<String>>,
) -> HashSet<String> {
    if let Some(found) = cache.get(id) {
        return found.clone();
    }
    let mut seen = HashSet::new();
    let mut stack = terms
        .get(id)
        .map(|term| term.parents.clone())
        .unwrap_or_default();
    while let Some(parent) = stack.pop() {
        if !seen.insert(parent.clone()) {
            continue;
        }
        if let Some(known) = cache.get(&parent) {
            seen.extend(known.iter().cloned());
        } else if let Some(term) = terms.get(&parent) {
            stack.extend(term.parents.iter().cloned());
        }
    }
    cache.insert(id.to_owned(), seen.clone());
    seen
}

fn ordered_join<'a>(values: impl IntoIterator<Item = &'a str>) -> String {
    let mut seen = HashSet::new();
    values
        .into_iter()
        .filter(|value| !value.is_empty() && seen.insert(*value))
        .collect::<Vec<_>>()
        .join(";")
}

fn require_columns(input: &Table, names: &[&str], label: &str) -> Result<()> {
    let mut missing: Vec<&str> = names
        .iter()
        .copied()
        .filter(|name| !input.has_column(name))
        .collect();
    missing.sort_unstable();
    if missing.is_empty() {
        Ok(())
    } else {
        Err(invalid(format!(
            "{label} are missing columns: {}",
            table::repr_list(&missing)
        )))
    }
}

fn load_overrides(path: Option<&Path>) -> Result<HashMap<String, HashMap<String, String>>> {
    let Some(path) = path else {
        return Ok(HashMap::new());
    };
    let metadata = match fs::metadata(path) {
        Ok(metadata) => metadata,
        Err(error)
            if matches!(
                error.kind(),
                io::ErrorKind::NotFound | io::ErrorKind::NotADirectory
            ) =>
        {
            return Ok(HashMap::new());
        }
        // pathlib in the pinned Python also treats EBADF and ELOOP as absent.
        // Rust 1.90's FilesystemLoop ErrorKind is not stable yet; use the OS's
        // errno values at this small compatibility boundary, without new FFI.
        #[cfg(any(target_os = "linux", target_os = "android"))]
        Err(error) if matches!(error.raw_os_error(), Some(9 | 40)) => return Ok(HashMap::new()),
        #[cfg(any(
            target_os = "macos",
            target_os = "ios",
            target_os = "freebsd",
            target_os = "openbsd",
            target_os = "netbsd",
            target_os = "dragonfly"
        ))]
        Err(error) if matches!(error.raw_os_error(), Some(9 | 62)) => return Ok(HashMap::new()),
        #[cfg(windows)]
        Err(error) if matches!(error.raw_os_error(), Some(21 | 123 | 1921)) => {
            return Ok(HashMap::new());
        }
        Err(error) => return Err(error.into()),
    };
    if metadata.len() == 0 {
        return Ok(HashMap::new());
    }
    let overrides = table::read(path, true)?;
    require_columns(
        &overrides,
        &[
            "disease_id",
            "disease_scope",
            "priva_scope",
            "scope_review_status",
        ],
        "Disease-scope overrides",
    )?;
    let unique = pandas_polars::drop_duplicates_first(&overrides, Some(&["disease_id"]))?;
    if unique.len() != overrides.rows {
        let unique: HashSet<usize> = unique.into_iter().collect();
        let duplicates: Vec<&str> = (0..overrides.rows)
            .filter(|row| !unique.contains(row))
            .take(5)
            .map(|row| table::text(&overrides, "disease_id", row))
            .collect();
        return Err(invalid(format!(
            "Duplicate disease-scope overrides: {}",
            table::repr_list(&duplicates)
        )));
    }
    let mut result = HashMap::new();
    for row in 0..overrides.rows {
        let id = python_strip(table::text(&overrides, "disease_id", row));
        if id.is_empty() {
            continue;
        }
        let disease = python_strip(table::text(&overrides, "disease_scope", row));
        let priva = python_strip(table::text(&overrides, "priva_scope", row));
        if !DISEASE_SCOPES.contains(&disease) {
            return Err(invalid(format!(
                "Invalid disease_scope for {id}: {disease}"
            )));
        }
        if !["include", "exclude", "review"].contains(&priva) {
            return Err(invalid(format!("Invalid priva_scope for {id}: {priva}")));
        }
        // Validation is on raw IDs before stripping; normalized collisions keep
        // the last row, exactly as the reference dict assignment does.
        result.insert(
            id.to_owned(),
            overrides
                .names
                .iter()
                .map(|name| {
                    (
                        name.clone(),
                        python_strip(table::text(&overrides, name, row)).to_owned(),
                    )
                })
                .collect(),
        );
    }
    Ok(result)
}

fn automatic(flags: [bool; 6]) -> [&'static str; 4] {
    let [
        mendelian,
        non_monogenic,
        hereditary,
        neoplasm,
        predisposition,
        susceptibility,
    ] = flags;
    let inherited = mendelian || hereditary;
    if mendelian && non_monogenic {
        [
            "conflicting_inheritance",
            "review",
            "needs_review",
            "HPO contains both Mendelian and non-monogenic inheritance assertions",
        ]
    } else if non_monogenic {
        [
            "complex_or_non_monogenic",
            "exclude",
            "auto_supported",
            "Explicit HPO non-monogenic inheritance",
        ]
    } else if inherited && neoplasm {
        [
            "mendelian_cancer_predisposition",
            "include",
            "auto_supported",
            "Inherited evidence and MONDO neoplasm ancestry",
        ]
    } else if inherited && (predisposition || susceptibility) {
        [
            "mendelian_susceptibility_uncertain",
            "review",
            "needs_review",
            "Inherited susceptibility/predisposition lacks explicit MONDO neoplasm ancestry",
        ]
    } else if inherited {
        [
            "mendelian_non_neoplastic",
            "include",
            "auto_supported",
            "Inherited evidence without MONDO neoplasm ancestry",
        ]
    } else if neoplasm {
        [
            "neoplastic_uncertain",
            "review",
            "needs_review",
            "MONDO neoplasm ancestry does not establish somatic-only etiology",
        ]
    } else if predisposition || susceptibility {
        [
            "complex_susceptibility_candidate",
            "review",
            "needs_review",
            "Susceptibility/predisposition lacks explicit inherited evidence",
        ]
    } else {
        [
            "uncertain",
            "review",
            "needs_review",
            "No explicit inherited, non-monogenic, or neoplastic scope evidence",
        ]
    }
}

pub(crate) fn build(
    obo_path: &Path,
    assertion_path: &Path,
    override_path: Option<&Path>,
) -> Result<(Table, Json)> {
    let (version, terms) = parse_obo(obo_path)?;
    let assertions = table::read(assertion_path, true)?;
    require_columns(&assertions, &BASE_COLUMNS, "HPO assertions")?;
    let assertions = assertions.take_rows(&pandas_polars::stable_sort(
        &assertions,
        &["disease_id"],
        &[false],
    )?);
    let groups = pandas_polars::first_appearance_groups(&assertions, &["disease_id"], true)?;
    let mut xrefs: HashMap<&str, BTreeSet<&str>> = HashMap::new();
    for term in terms.values().filter(|term| !term.obsolete) {
        for xref in &term.xrefs {
            xrefs.entry(xref).or_default().insert(&term.id);
        }
    }
    let overrides = load_overrides(override_path)?;
    let mut cache = HashMap::new();
    let mut records = Vec::with_capacity(groups.len());
    let mut mapped_diseases = 0;
    let mut override_diseases = 0;
    for group in &groups {
        let disease = table::text(&assertions, "disease_id", group[0]);
        let hpos: HashSet<&str> = group
            .iter()
            .map(|&row| table::text(&assertions, "hpo_id", row))
            .collect();
        let inheritance_refs = ordered_join(
            group
                .iter()
                .filter(|&&row| {
                    let hpo = table::text(&assertions, "hpo_id", row);
                    MENDELIAN.contains(&hpo) || NON_MONOGENIC.contains(&hpo)
                })
                .map(|&row| table::text(&assertions, "reference", row))
                .filter(|value| *value != "-"),
        );
        let ids: Vec<&str> = xrefs
            .get(disease)
            .map(|values| values.iter().copied().collect())
            .unwrap_or_default();
        if !ids.is_empty() {
            mapped_diseases += 1;
        }
        let mapped: Vec<&Term> = ids.iter().map(|id| &terms[*id]).collect();
        let mut all_ancestors = HashSet::new();
        for id in &ids {
            all_ancestors.insert((*id).to_owned());
            all_ancestors.extend(ancestors(id, &terms, &mut cache));
        }
        let subsets: HashSet<&str> = mapped
            .iter()
            .flat_map(|term| term.subsets.iter().map(String::as_str))
            .collect();
        let flags = [
            MENDELIAN.iter().any(|hpo| hpos.contains(hpo)),
            NON_MONOGENIC.iter().any(|hpo| hpos.contains(hpo)),
            all_ancestors.contains(HEREDITARY),
            all_ancestors.contains(NEOPLASM),
            subsets.contains("predisposition"),
            subsets.contains("omim_susceptibility"),
        ];
        let [scope, priva, review, note] = automatic(flags);
        let evidence_names = [
            "HPO_mendelian_inheritance",
            "HPO_non_monogenic_inheritance",
            "MONDO_ancestor:MONDO:0003847",
            "MONDO_ancestor:MONDO:0005070",
            "MONDO_subset:predisposition",
            "MONDO_subset:omim_susceptibility",
        ];
        let evidence = ordered_join(
            evidence_names
                .iter()
                .enumerate()
                .filter(|(index, _)| flags[*index])
                .map(|(_, value)| *value),
        );
        let release_ref = format!("MONDO_RELEASE:{version}");
        let mut record = vec![
            disease.to_owned(),
            ordered_join(ids.iter().copied()),
            ordered_join(mapped.iter().map(|term| term.name.as_str())),
            version.clone(),
        ];
        record.extend(flags.iter().map(|flag| u8::from(*flag).to_string()));
        record.extend([
            scope.to_owned(),
            priva.to_owned(),
            if evidence.is_empty() {
                "none".to_owned()
            } else {
                evidence
            },
            ordered_join([release_ref.as_str(), inheritance_refs.as_str()]),
            review.to_owned(),
            note.to_owned(),
        ]);
        if let Some(overrides) = overrides.get(disease) {
            override_diseases += 1;
            for name in [
                "disease_scope",
                "priva_scope",
                "scope_review_status",
                "scope_evidence",
                "scope_reference",
                "scope_note",
            ] {
                if let Some(value) = overrides.get(name).filter(|value| !value.is_empty()) {
                    record[REGISTRY_COLUMNS
                        .iter()
                        .position(|column| *column == name)
                        .unwrap()] = value.clone();
                }
            }
            if overrides.get("scope_evidence").is_none_or(String::is_empty) {
                record[12] = ordered_join([record[12].as_str(), "manual_override"]);
            }
            if overrides
                .get("scope_review_status")
                .is_none_or(String::is_empty)
            {
                record[14] = "manually_confirmed".to_owned();
            }
        }
        records.push(record);
    }
    let mut registry = table::strings(&REGISTRY_COLUMNS, records);
    for name in &REGISTRY_COLUMNS[4..10] {
        let flags = (0..registry.rows)
            .map(|row| i64::from(table::text(&registry, name, row) == "1"))
            .collect();
        registry.set_column(name, Column::Int(flags));
    }
    let mut summary = match table::json_counts([
        ("assertion_rows", assertions.rows),
        ("diseases", registry.rows),
        ("mondo_terms", terms.len()),
        ("mapped_diseases", mapped_diseases),
        ("unmapped_diseases", registry.rows - mapped_diseases),
        ("override_diseases", override_diseases),
        ("loaded_overrides", overrides.len()),
    ]) {
        Json::Object(fields) => fields,
        _ => unreachable!(),
    };
    summary.push(("mondo_version".to_owned(), Json::Str(version)));
    for (name, key) in [
        ("disease_scope", "scope_counts"),
        ("priva_scope", "priva_counts"),
        ("scope_review_status", "review_counts"),
    ] {
        let mut counts = BTreeMap::new();
        for row in 0..registry.rows {
            *counts
                .entry(table::text(&registry, name, row))
                .or_insert(0i128) += 1;
        }
        summary.push((
            key.to_owned(),
            Json::Object(
                counts
                    .into_iter()
                    .map(|(name, count)| (name.to_owned(), Json::Int(count)))
                    .collect(),
            ),
        ));
    }
    Ok((registry, Json::Object(summary)))
}

pub(crate) fn annotate(assertions: &Table, registry: &Table) -> Result<Table> {
    let names: Vec<&str> = assertions
        .names
        .iter()
        .map(String::as_str)
        .filter(|name| !SCOPE_COLUMNS.contains(name))
        .collect();
    let base = table::select(assertions, &names)?;
    let mut registry_names = vec!["disease_id"];
    registry_names.extend(SCOPE_COLUMNS);
    let registry = table::select(registry, &registry_names)?;
    let joined = pandas_polars::left_join_text_rows(&base, &registry, &["disease_id"], true)?;
    let mut missing = Vec::new();
    let mut seen = HashSet::new();
    for &(left, right) in &joined {
        if right.is_none_or(|row| table::text(&registry, "priva_scope", row).is_empty()) {
            let disease = table::text(&base, "disease_id", left);
            if seen.insert(disease) && missing.len() < 5 {
                missing.push(disease);
            }
        }
    }
    if !missing.is_empty() {
        return Err(invalid(format!(
            "Disease scope missing for HPO assertions: {}",
            table::repr_list(&missing)
        )));
    }
    let mut scoped = base.take_rows(&joined.iter().map(|(row, _)| *row).collect::<Vec<_>>());
    for name in SCOPE_COLUMNS {
        scoped.set_column(
            name,
            Column::Text(
                joined
                    .iter()
                    .map(|(_, right)| right.map(|row| table::text(&registry, name, row).to_owned()))
                    .collect(),
            ),
        );
    }
    Ok(scoped)
}

const USAGE: &str = "priva mondo-scope build --mondo-obo FILE --hpo-assertions TSV --registry-output TSV[.gz] [--manual-overrides TSV] [--annotated-assertions-output TSV[.gz]]\npriva mondo-scope annotate --hpo-assertions TSV --registry TSV --output TSV[.gz]\nBuild or attach disease scope using existing MONDO/HPO rules. Success prints JSON coverage metadata. Installer callers are unchanged.";

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args.is_empty() || args.iter().any(|arg| arg == "--help" || arg == "-h") {
        println!("{USAGE}");
        return Ok(());
    }
    let names: &[&str] = match args[0].as_str() {
        "build" => &[
            "mondo-obo",
            "hpo-assertions",
            "registry-output",
            "manual-overrides",
            "annotated-assertions-output",
        ],
        "annotate" => &["hpo-assertions", "registry", "output"],
        _ => return Err(invalid(USAGE)),
    };
    let Parsed::Values(options) =
        python_args::parse_long_options(&args[1..], names).map_err(invalid)?
    else {
        println!("{USAGE}");
        return Ok(());
    };
    let required = |name: &str| {
        options
            .get(name)
            .map(|path| Path::new(path))
            .ok_or_else(|| invalid(format!("missing --{name}")))
    };
    if args[0] == "annotate" {
        let assertions_path = required("hpo-assertions")?;
        let registry_path = required("registry")?;
        let output_path = required("output")?;
        return crate::builder_error::finish((|| {
            let assertions = table::read(assertions_path, true)?;
            let registry = table::read(registry_path, true)?;
            let output = annotate(&assertions, &registry)?;
            table::publish(&output, output_path, Publication::Mondo)?;
            table::print_summary(&table::json_counts([
                ("assertion_rows", output.rows),
                ("registry_rows", registry.rows),
                ("annotated_rows", output.rows),
            ]))
        })());
    }
    let obo = required("mondo-obo")?;
    let assertion_path = required("hpo-assertions")?;
    let output = required("registry-output")?;
    crate::builder_error::finish((|| {
        let (registry, mut summary) = build(
            obo,
            assertion_path,
            options.get("manual-overrides").map(Path::new),
        )?;
        table::publish(&registry, output, Publication::Mondo)?;
        if let Some(output) = options.get("annotated-assertions-output") {
            let scoped = annotate(&table::read(assertion_path, true)?, &registry)?;
            table::publish(&scoped, Path::new(output), Publication::Mondo)?;
            if let Json::Object(fields) = &mut summary {
                fields.push(("annotated_rows".to_owned(), Json::Int(scoped.rows as i128)));
            }
        }
        table::print_summary(&summary)
    })())
}
