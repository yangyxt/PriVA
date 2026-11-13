//! The original per-gene procedure, including pre-filter symbol selection and
//! ordered mutation of all seven gene/transcript/domain variant sets.
use super::maps::{Alleles, Key, TargetMap};
use crate::cache_convert::clinvar_gene_stat::{PlainGene, PlainSets, PlainStats};
use crate::cli::{Result, invalid};
use crate::combine_annotations::gene_stats::{Scalar, Table};
use crate::python_set_order::StringSet;
use std::collections::{BTreeMap, BTreeSet, HashMap};

const REVIEWS: [&str; 3] = [
    "practice_guideline",
    "reviewed_by_expert_panel",
    "criteria_provided,_multiple_submitters,_no_conflicts",
];
type Sets = [BTreeSet<String>; 7];

fn plain(sets: Sets) -> PlainSets {
    sets.map(|set| set.into_iter().collect())
}
fn blank() -> Sets {
    std::array::from_fn(|_| BTreeSet::new())
}

fn variant(table: &Table, row: usize) -> Result<Option<String>> {
    let mut parts = Vec::new();
    for name in ["chrom", "ref", "alt"] {
        let value = table.get(name, row)?;
        if value.missing() {
            return Ok(None);
        }
        let Scalar::Text(text) = value else {
            return Err(invalid("TypeError: a variant-ID component is not text"));
        };
        parts.push(text);
    }
    Ok(Some(format!(
        "{}_{}_{}_{}",
        parts[0],
        table.get("pos", row)?.text(),
        parts[1],
        parts[2]
    )))
}

fn magnitude(text: &str) -> Result<u64> {
    let text = text.trim_start_matches(['+', '-']);
    if !crate::python_text::python_int_is_valid(text) {
        return Err(invalid("ValueError: splice-event offset is not an integer"));
    }
    let text = crate::python_text::python_strip(text).trim_start_matches(['+', '-']);
    let ascii: String = text
        .chars()
        .filter(|c| *c != '_')
        .map(|c| char::from(b'0' + crate::python_text::decimal_digit_value(c).unwrap() as u8))
        .collect();
    Ok(ascii.parse().unwrap_or(u64::MAX))
}

/// [large, small, frameshift]; no simplification of the source's float sums or
/// threshold ordering. A NaN delta does not satisfy the early `<0.3` return.
fn splice_event(event: &str, delta: f64) -> Result<Option<[f64; 3]>> {
    if [
        "NA", "None", "none", "NONE", "Null", "null", "NULL", "nan", "NaN", "NAN", "N/A",
    ]
    .contains(&event)
        || delta < 0.3
    {
        return Ok(None);
    }
    let mut result = [0.0; 3];
    for event in event.split('&') {
        let parts: Vec<_> = event.split(':').collect();
        if parts.len() != 5 {
            return Err(invalid(
                "ValueError: splice event does not have five colon-separated fields",
            ));
        }
        let freq = parts[3];
        let text = if freq.contains('%') {
            freq.trim_end_matches('%')
        } else if freq.contains('#') {
            freq.trim_start_matches('#')
        } else {
            return Err(invalid("ValueError: unrecognized splice-event frequency"));
        };
        let freq = Scalar::Text(text).number_or_nan()? / 100.0;
        if parts[1] == "ES" {
            result[if parts[4] == "Frameshift" { 2 } else { 0 }] += freq;
        } else if ["CD", "CA"].contains(&parts[1]) {
            let offset = magnitude(parts[2])?;
            if parts[4] == "Frameshift" {
                result[2] += freq;
            } else if offset > 15 {
                result[0] += freq;
            } else if (3..=9).contains(&offset) {
                result[1] += freq;
            }
        }
    }
    Ok(Some(result))
}

fn per_gene(
    table: &Table,
    rows: &[usize],
    variants: &[Option<String>],
    seed: u32,
    cutoff: f64,
) -> Result<Option<PlainGene>> {
    let first = rows[0];
    let selected: Vec<usize> = rows.iter().copied().filter(|row| {
        matches!(table.get("BIOTYPE", *row), Ok(Scalar::Text("protein_coding")))
            && matches!(table.get("CLNREVSTAT", *row), Ok(Scalar::Text(value)) if REVIEWS.contains(&value))
    }).collect();
    if selected.is_empty() {
        return Ok(None);
    }
    let Scalar::Text(gene) = table.get("Gene", first)? else {
        return Err(invalid(
            "a non-text gene identifier has no schema-5 representation",
        ));
    };
    let symbol = table.get("SYMBOL", first)?.symbol()?;
    let mut transcripts = StringSet::new(seed);
    for &row in &selected {
        let Scalar::Text(transcript) = table.get("Feature", row)? else {
            return Err(invalid(
                "a non-text transcript has no schema-5 representation",
            ));
        };
        transcripts.insert(transcript.to_owned());
    }
    // .str is attempted even when dropna has removed every numeric value.
    if table.numeric_column("DOMAINS")? {
        return Err(invalid(
            "AttributeError: the DOMAINS column has no .str accessor",
        ));
    }
    let mut domains = StringSet::new(seed);
    for &row in &selected {
        if let Scalar::Text(text) = table.get("DOMAINS", row)? {
            for domain in text.split('&').filter(|domain| !domain.contains("nan")) {
                domains.insert(format!("{gene}:{domain}"));
            }
        }
    }
    let mut transcripts: Vec<_> = transcripts
        .into_order()
        .into_iter()
        .map(|key| (key, blank()))
        .collect();
    let mut domains: Vec<_> = domains
        .into_order()
        .into_iter()
        .map(|key| (key, blank()))
        .collect();
    let transcript_index: HashMap<_, _> = transcripts
        .iter()
        .enumerate()
        .map(|(at, (key, _))| (key.clone(), at))
        .collect();
    let domain_index: HashMap<_, _> = domains
        .iter()
        .enumerate()
        .map(|(at, (key, _))| (key.clone(), at))
        .collect();
    let mut sets = blank();
    for &row in &selected {
        sets[0].insert(variants[row].clone().unwrap_or_else(|| "nan".into()));
    }
    for row in selected {
        let variant = variants[row].as_deref().unwrap_or("nan");
        let transcript = table.get("Feature", row)?.text();
        let Some(&tranx) = transcript_index.get(&transcript) else {
            continue;
        };
        let domain_text = table.get("DOMAINS", row)?.text_or_nan();
        let domain_rows: Vec<_> = domain_text.split('&').filter(|value| !value.contains("nan"))
            .map(|domain| domain_index.get(&format!("{gene}:{domain}")).copied()
                .ok_or_else(|| invalid("KeyError: a row domain is absent from the gene's .str-derived domain dictionary")))
            .collect::<Result<_>>()?;
        transcripts[tranx].1[0].insert(variant.to_owned());
        for &domain in &domain_rows {
            domains[domain].1[0].insert(variant.to_owned());
        }
        let mut add = |index: usize| {
            sets[index].insert(variant.to_owned());
            transcripts[tranx].1[index].insert(variant.to_owned());
            for &domain in &domain_rows {
                domains[domain].1[index].insert(variant.to_owned());
            }
        };
        add(if table.get("CLNSIG", row)?.text().contains("athogenic") {
            1
        } else {
            2
        });
        let delta = table
            .get("SpliceVault_SpliceAI_delta", row)?
            .number_or_nan()?;
        let splicing = splice_event(
            &table.get("SpliceVault_top_events", row)?.text_or_nan(),
            delta,
        )?;
        let lof = table
            .get("LoF", row)?
            .text_or_nan()
            .split(['&', ',', ';', '|'])
            .map(|value| crate::python_text::python_strip(value).to_uppercase())
            .any(|value| ["HC", "OS"].contains(&value.as_str()));
        let vep_lof = lof && !table.get("NMD", row)?.text().contains("escape");
        let donor = table.get("SpliceAI_pred_DS_DL", row)?.number_or_nan()? > cutoff;
        let nmd = vep_lof || donor;
        if nmd {
            add(3);
        }
        let consequence = table.get("Consequence", row)?.text();
        let inframe =
            consequence.contains("inframe_deletion") || consequence.contains("inframe_insertion");
        let length = if inframe {
            let (Scalar::Text(reference), Scalar::Text(alternative)) =
                (table.get("ref", row)?, table.get("alt", row)?)
            else {
                return Err(invalid(
                    "TypeError: an in-frame variant allele has no length",
                ));
            };
            reference
                .chars()
                .count()
                .abs_diff(alternative.chars().count())
        } else {
            0
        };
        let large_indel = inframe && length >= 15;
        let small_indel = inframe && length <= 9;
        let (large_splice, small_splice) =
            splicing.map_or((false, false), |[large, small, frameshift]| {
                (
                    ((frameshift + large) >= 0.2 && delta >= 0.5)
                        || ((frameshift + large) >= 0.4 && delta >= 0.3),
                    (small >= 0.2 && delta >= 0.3) || (small >= 0.4 && delta >= 0.1),
                )
            });
        if (large_indel || large_splice) && !nmd {
            add(5);
        }
        let missense = consequence.contains("missense_variant");
        if (small_indel || small_splice || missense) && !nmd && !large_indel && !large_splice {
            add(4);
        }
        if consequence.contains("synonymous_variant")
            && !nmd
            && !large_indel
            && !large_splice
            && !small_splice
            && !missense
            && !small_indel
        {
            add(6);
        }
    }
    Ok(Some(PlainGene {
        ensembl_id: gene.to_owned(),
        gene_symbol: symbol,
        sets: plain(sets),
        domains: domains
            .into_iter()
            .map(|(key, sets)| (key, plain(sets)))
            .collect(),
        transcripts: transcripts
            .into_iter()
            .map(|(key, sets)| (key, plain(sets)))
            .collect(),
    }))
}

pub(super) fn build(
    table: &Table,
    targets: &TargetMap,
    alleles: &Alleles,
    seed: u32,
    cutoff: f64,
) -> Result<PlainStats> {
    // Series string concatenation rejects a numeric nullable column even if
    // every cell is missing; it does not turn that whole column into NaN IDs.
    for name in ["chrom", "ref", "alt"] {
        if table.numeric_column(name)? {
            return Err(invalid("TypeError: a variant-ID column is numeric"));
        }
    }
    let variants = (0..table.rows)
        .map(|row| variant(table, row))
        .collect::<Result<Vec<_>>>()?;
    let promote = alleles.numeric
        && (alleles.nullable
            || variants.iter().any(|id| {
                id.as_ref()
                    .is_none_or(|id| !alleles.values.contains_key(id))
            }));
    let mut groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (row, variant) in variants.iter().enumerate() {
        let allele = variant
            .as_ref()
            .and_then(|variant| alleles.values.get(variant))
            .unwrap_or(&Key::Missing);
        let promoted;
        let allele = if promote {
            promoted = allele.promoted();
            &promoted
        } else {
            allele
        };
        let Some(Some(target)) = targets.get(allele, promote) else {
            continue;
        };
        if !matches!(table.get("Feature", row)?, Scalar::Text(feature) if feature == target) {
            continue;
        }
        match table.get("Gene", row)? {
            Scalar::Text(gene) => groups.entry(gene.to_owned()).or_default().push(row),
            value if value.missing() => {}
            _ if !matches!(table.get("BIOTYPE", row)?, Scalar::Text("protein_coding"))
                || !matches!(table.get("CLNREVSTAT", row)?, Scalar::Text(value) if REVIEWS.contains(&value)) =>
            {
                // A wholly filtered numeric-key group returns None before the
                // source concatenates its gene identifier with domain text.
            }
            _ => {
                return Err(invalid(
                    "UFuncTypeError: a surviving numeric gene identifier cannot prefix a domain",
                ));
            }
        }
    }
    let mut result = Vec::new();
    for (gene, rows) in groups {
        if let Some(stat) = per_gene(table, &rows, &variants, seed, cutoff)? {
            result.push((gene, stat));
        }
    }
    Ok(result)
}
