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

/// Coldspot intervals by gene label, from the qualified TSV reader.
pub(crate) struct Coldspots(HashMap<String, Vec<Interval>>);

pub(crate) struct Summary {
    pub(crate) of_coldspot: f64,
    pub(crate) of_domain: f64,
    pub(crate) removed: bool,
}

pub(crate) fn summary(
    domain: &str,
    metadata: &HashMap<String, metadata::Entry>,
    coldspots: &Coldspots,
    lookup: &DomainLookup,
    threshold: f64,
) -> Result<Summary> {
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
            if let Some(spots) = coldspots.0.get(&meta.symbol) {
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
    Ok(Summary {
        of_coldspot: a,
        of_domain: b,
        removed,
    })
}

/// Parse/convert every selected row before looking up any domain, including
/// numeric/Boolean/missing gene labels that cannot equal metadata's text.
pub(crate) fn read_coldspots(path: &Path) -> Result<Coldspots> {
    let rows = das_coldspot_table::read(path)?;
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
    Ok(Coldspots(coldspots))
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
    let coldspots = read_coldspots(Path::new(options["--spots"]))?;
    let input = BufReader::new(std::fs::File::open(options["--domains"])?);
    let mut output = Vec::new();
    for line in input.lines() {
        let Json::Str(domain) = json_text::parse(&line?, true).map_err(invalid)? else {
            return Err(invalid("each domain-path line must be a JSON string"));
        };
        let row = summary(&domain, &metadata, &coldspots, &lookup, threshold)?;
        output.push(Json::Object(vec![
            ("domain".into(), Json::Str(domain)),
            ("coldspot_max_frac_of_coldspot".into(), Json::Float(row.of_coldspot)),
            ("coldspot_max_frac_of_domain".into(), Json::Float(row.of_domain)),
            ("removed".into(), Json::Bool(row.removed)),
        ]));
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
