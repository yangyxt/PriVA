//! Whole coldspot procedure over native metadata and the qualified TSV reader.
//! The caller supplies domain iteration order. Returned-set table history is
//! separate; this boundary emits each requested domain's scientific summary.
use crate::alpha_domain::metadata;
use crate::builder_error;
use crate::cache_convert::interpro_mapping::InterproCache;
use crate::cli::{Result, invalid};
use crate::das_coldspot_table::{self, Gene};
use crate::das_domain::DomainLookup;
use crate::das_integer_ratio;
use crate::json_text::{self, Json};
use crate::protein_domain_mapping::DomainNormalizer;
use num_bigint::BigInt;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

struct Interval {
    start: BigInt,
    end: BigInt,
}

fn fraction(a: &BigInt, b: &BigInt) -> Result<f64> {
    das_integer_ratio::divide(a, b).map_err(|error| {
        builder_error::category(error.category(), invalid("integer overlap division failed"))
    })
}

fn summary(
    domain: &str,
    metadata: &HashMap<String, metadata::Entry>,
    coldspots: &HashMap<String, Vec<Interval>>,
    lookup: &DomainLookup,
    threshold: f64,
) -> Result<Json> {
    let entry = lookup.first(domain).map_err(invalid)?;
    let mut a = 0.0;
    let mut b = 0.0;
    let mut removed = false;
    if entry.is_some_and(|value| {
        matches!(
            value.entry.kind.as_str(),
            "Family" | "Homologous_superfamily"
        )
    }) {
        if let Some(meta) = metadata.get(domain) {
            if let Some(spots) = coldspots.get(&meta.symbol) {
                let domain_start = BigInt::from(meta.start);
                let domain_end = BigInt::from(meta.end);
                let domain_length = &domain_end - &domain_start + 1;
                for spot in spots {
                    let start = (&domain_start).max(&spot.start);
                    let end = (&domain_end).min(&spot.end);
                    if start > end {
                        continue;
                    }
                    let overlap = end - start + 1;
                    let coldspot_length = &spot.end - &spot.start + 1;
                    let x = fraction(&overlap, &coldspot_length)?;
                    let y = fraction(&overlap, &domain_length)?;
                    if x > a {
                        a = x;
                    }
                    if y > b {
                        b = y;
                    }
                }
                removed = a > threshold || b > threshold;
            }
        }
    }
    Ok(Json::Object(vec![
        ("domain".into(), Json::Str(domain.into())),
        ("coldspot_max_frac_of_coldspot".into(), Json::Float(a)),
        ("coldspot_max_frac_of_domain".into(), Json::Float(b)),
        ("removed".into(), Json::Bool(removed)),
    ]))
}

pub(crate) fn lines(options: &HashMap<&str, &str>, seed: u32, threshold: f64) -> Result<()> {
    let metadata = metadata::read(Path::new(options["--metadata"]))?
        .entries
        .into_iter()
        .collect();
    let normalizer = DomainNormalizer::from_parts(
        &std::fs::read_to_string(options["--normalization-json"])?,
        HashMap::new(),
    )?;
    let lookup = DomainLookup::new(
        normalizer,
        InterproCache::open(options["--interpro"])?,
        seed,
    );
    // Parse/convert every selected row before looking up any domain, including
    // numeric/Boolean/missing gene labels that cannot equal metadata's text.
    let rows = das_coldspot_table::read(Path::new(options["--spots"]))?;
    let mut coldspots = HashMap::<String, Vec<Interval>>::new();
    for row in rows {
        let interval = Interval {
            start: row
                .start
                .parse()
                .map_err(|_| invalid("invalid parsed coldspot start"))?,
            end: row
                .end
                .parse()
                .map_err(|_| invalid("invalid parsed coldspot end"))?,
        };
        if let Gene::Text(gene) = row.gene {
            coldspots.entry(gene).or_default().push(interval);
        }
    }
    let input = BufReader::new(std::fs::File::open(options["--domains"])?);
    let mut output = Vec::new();
    for line in input.lines() {
        let Json::Str(domain) = json_text::parse(&line?, true).map_err(invalid)? else {
            return Err(invalid("each domain-path line must be a JSON string"));
        };
        output.push(summary(&domain, &metadata, &coldspots, &lookup, threshold)?);
    }
    // The original returns only after the complete procedure succeeds. Avoid
    // publishing a scientific prefix when a later domain raises an exception.
    let mut stream = std::io::BufWriter::new(std::io::stdout().lock());
    for row in output {
        writeln!(
            stream,
            "{}",
            json_text::dumps(&row, false).map_err(invalid)?
        )?;
    }
    stream.flush()?;
    Ok(())
}
