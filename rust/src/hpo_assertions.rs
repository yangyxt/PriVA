//! Opt-in equivalent of `scripts/collapse_HPO_anno.py`.
//! Gene links retain the complete many-to-many HPOA join. The release match
//! threshold is evaluated on expanded rows, before NOT-only removal or final
//! deduplication. The reference's unmatched `astype(str)` values stay `nan`.
use crate::builder_table::{self as table, Publication};
use crate::cli::{Result, invalid};
use crate::json_text::Json;
use crate::pandas_polars;
use crate::pandas_table::{Column, Table};
use crate::python_args::{self, Parsed};
use crate::python_text::{parse_float, python_strip};
use std::collections::HashSet;
use std::path::Path;

pub(crate) const COLUMNS: [&str; 6] = [
    "gene_symbol",
    "disease_id",
    "hpo_id",
    "frequency",
    "evidence",
    "reference",
];
const GENE_COLUMNS: [&str; 4] = ["disease_id", "frequency", "gene_symbol", "hpo_id"];
const HPOA_COLUMNS: [&str; 6] = [
    "database_id",
    "evidence",
    "frequency",
    "hpo_id",
    "qualifier",
    "reference",
];

fn stripped(mut input: Table) -> Table {
    for column in &mut input.columns {
        let Column::Text(values) = column else {
            unreachable!("dtype=str input")
        };
        for value in values {
            *value = Some(python_strip(value.as_deref().unwrap_or("nan")).to_owned());
        }
    }
    input
}

pub(crate) fn build(
    genes_path: &Path,
    hpoa_path: &Path,
    minimum_match_rate: f64,
) -> Result<(Table, Json)> {
    let genes = table::read(genes_path, false)?;
    let hpoa = table::read(hpoa_path, true)?;
    table::require(&genes, &GENE_COLUMNS, "genes_to_phenotype")?;
    table::require(&hpoa, &HPOA_COLUMNS, "phenotype.hpoa")?;
    let genes = stripped(table::select(&genes, &GENE_COLUMNS)?);
    let hpoa = stripped(table::select(&hpoa, &HPOA_COLUMNS)?);
    let input_gene_rows = genes.rows;
    let invalid_identifiers = (0..genes.rows)
        .filter(|&row| {
            ["disease_id", "hpo_id"]
                .iter()
                .any(|name| ["", "-"].contains(&table::text(&genes, name, row)))
        })
        .count();
    if invalid_identifiers > 0 {
        return Err(invalid(format!(
            "genes_to_phenotype contains {invalid_identifiers} rows with an empty disease_id or hpo_id"
        )));
    }
    let mapped: Vec<usize> = (0..genes.rows)
        .filter(|&row| !["", "-"].contains(&table::text(&genes, "gene_symbol", row)))
        .collect();
    let unmapped_gene_rows = genes.rows - mapped.len();
    let genes = genes.take_rows(&mapped);
    let negative_pairs: HashSet<(&str, &str)> = (0..hpoa.rows)
        .filter(|&row| table::text(&hpoa, "qualifier", row) == "NOT")
        .map(|row| {
            (
                table::text(&hpoa, "database_id", row),
                table::text(&hpoa, "hpo_id", row),
            )
        })
        .collect();
    let positives = (0..hpoa.rows)
        .filter(|&row| table::text(&hpoa, "qualifier", row).is_empty())
        .map(|row| {
            vec![
                table::text(&hpoa, "database_id", row).to_owned(),
                table::text(&hpoa, "hpo_id", row).to_owned(),
                table::optional(table::text(&hpoa, "frequency", row)),
                table::optional(table::text(&hpoa, "evidence", row)),
                table::optional(table::text(&hpoa, "reference", row)),
            ]
        })
        .collect();
    let positives = table::strings(
        &["disease_id", "hpo_id", "frequency", "evidence", "reference"],
        positives,
    );
    let positives = positives.take_rows(&pandas_polars::drop_duplicates_first(&positives, None)?);
    let joined =
        pandas_polars::left_join_text_rows(&genes, &positives, &["disease_id", "hpo_id"], false)?;
    let mut negative_only = 0;
    let mut unknown = 0;
    let mut matched = 0;
    for &(left, right) in &joined {
        if right.is_some() {
            matched += 1;
        } else if negative_pairs.contains(&(
            table::text(&genes, "disease_id", left),
            table::text(&genes, "hpo_id", left),
        )) {
            negative_only += 1;
        } else {
            unknown += 1;
        }
    }
    let match_rate = 1.0
        - if joined.is_empty() {
            0.0
        } else {
            unknown as f64 / joined.len() as f64
        };
    if match_rate < minimum_match_rate {
        return Err(invalid(format!(
            "Only {:.3}% of gene-disease-HPO rows matched phenotype.hpoa; the source files are probably from different HPO releases",
            match_rate * 100.0
        )));
    }
    let mut rows = Vec::with_capacity(joined.len() - negative_only);
    for &(left, right) in &joined {
        let disease = table::text(&genes, "disease_id", left);
        let hpo = table::text(&genes, "hpo_id", left);
        if right.is_none() && negative_pairs.contains(&(disease, hpo)) {
            continue;
        }
        let (frequency, evidence, reference) = match right {
            Some(right) => (
                table::text(&positives, "frequency", right).to_owned(),
                table::optional(table::text(&positives, "evidence", right)),
                table::optional(table::text(&positives, "reference", right)),
            ),
            None => (
                table::optional(table::text(&genes, "frequency", left)),
                "nan".to_owned(),
                "nan".to_owned(),
            ),
        };
        rows.push(vec![
            table::text(&genes, "gene_symbol", left).to_owned(),
            disease.to_owned(),
            hpo.to_owned(),
            frequency,
            evidence,
            reference,
        ]);
    }
    let assertions = table::strings(&COLUMNS, rows);
    let assertions =
        assertions.take_rows(&pandas_polars::drop_duplicates_first(&assertions, None)?);
    let assertions = assertions.take_rows(&pandas_polars::stable_sort(
        &assertions,
        &COLUMNS,
        &[false; 6],
    )?);
    let mut summary = table::json_counts([
        ("input_gene_rows", input_gene_rows),
        ("input_hpoa_rows", hpoa.rows),
        ("unmapped_gene_rows", unmapped_gene_rows),
        ("positive_hpoa_assertions", positives.rows),
        ("joined_rows", joined.len()),
        ("matched_rows", matched),
        ("negative_only_rows", negative_only),
        ("unknown_rows", unknown),
        ("assertion_rows", assertions.rows),
        (
            "genes",
            pandas_polars::drop_duplicates_first(&assertions, Some(&["gene_symbol"]))?.len(),
        ),
        (
            "diseases",
            pandas_polars::drop_duplicates_first(&assertions, Some(&["disease_id"]))?.len(),
        ),
    ]);
    if let Json::Object(fields) = &mut summary {
        fields.push(("match_rate".to_owned(), Json::Float(match_rate)));
    }
    Ok((assertions, summary))
}

const USAGE: &str = "priva hpo-assertions build --genes-to-phenotype FILE --phenotype-hpoa FILE --output TSV[.gz] [--minimum-match-rate 0.999] [--disease-scope-registry TSV]\nBuild release-matched, disease-specific HPO assertions. Existing installer callers are unchanged. Coverage metadata is printed as JSON after successful publication.";

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args.is_empty() || args.iter().any(|arg| arg == "--help" || arg == "-h") {
        println!("{USAGE}");
        return Ok(());
    }
    if args[0] != "build" {
        return Err(invalid(USAGE));
    }
    let Parsed::Values(options) = python_args::parse_long_options(
        &args[1..],
        &[
            "genes-to-phenotype",
            "phenotype-hpoa",
            "output",
            "minimum-match-rate",
            "disease-scope-registry",
        ],
    )
    .map_err(invalid)?
    else {
        println!("{USAGE}");
        return Ok(());
    };
    let required = |name: &str| {
        options
            .get(name)
            .ok_or_else(|| invalid(format!("missing --{name}")))
    };
    let genes = required("genes-to-phenotype")?;
    let hpoa = required("phenotype-hpoa")?;
    let output = required("output")?;
    let minimum = parse_float(
        options
            .get("minimum-match-rate")
            .map(String::as_str)
            .unwrap_or("0.999"),
    )
    .ok_or_else(|| invalid("--minimum-match-rate requires a float"))?;
    crate::builder_error::finish((|| {
        let (mut assertions, mut summary) = build(Path::new(genes), Path::new(hpoa), minimum)?;
        if let Some(path) = options.get("disease-scope-registry") {
            assertions =
                crate::mondo_scope::annotate(&assertions, &table::read(Path::new(path), false)?)?;
            if let Json::Object(fields) = &mut summary {
                fields.push((
                    "annotated_rows".to_owned(),
                    Json::Int(assertions.rows as i128),
                ));
            }
        }
        table::publish(&assertions, Path::new(output), Publication::Hpo)?;
        table::print_summary(&summary)
    })())
}
