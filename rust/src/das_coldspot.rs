//! DAS coldspot filtering after the caller's exact gene lookup and TSV parse.
//! Inclusive overlap and strict thresholds follow the original Python loop.
use crate::das_reference::integer_ratio;

#[derive(Clone, Copy)]
pub(crate) struct Interval {
    pub(crate) start: i64,
    pub(crate) end: i64,
}

pub(crate) struct Summary {
    pub(crate) fraction_of_coldspot: f64,
    pub(crate) fraction_of_domain: f64,
    pub(crate) removed: bool,
    pub(crate) evaluated: bool,
    pub(crate) exempt: bool,
}

pub(crate) fn filter(
    kind: Option<&str>,
    domain: Option<Interval>,
    coldspots: Option<&[Interval]>,
    threshold: f64,
) -> Summary {
    let mut result = Summary {
        fraction_of_coldspot: 0.0,
        fraction_of_domain: 0.0,
        removed: false,
        evaluated: false,
        exempt: !matches!(kind, Some("Family" | "Homologous_superfamily")),
    };
    if result.exempt {
        return result;
    }
    let (Some(domain), Some(coldspots)) = (domain, coldspots) else {
        return result;
    };
    result.evaluated = true;
    let domain_length = i128::from(domain.end) - i128::from(domain.start) + 1;
    for coldspot in coldspots {
        let start = domain.start.max(coldspot.start);
        let end = domain.end.min(coldspot.end);
        if start > end {
            continue;
        }
        // A nonempty intersection implies both inclusive lengths are positive,
        // even when supplied reversed/empty intervals were not prevalidated.
        let overlap = i128::from(end) - i128::from(start) + 1;
        let coldspot_length = i128::from(coldspot.end) - i128::from(coldspot.start) + 1;
        let a = integer_ratio(overlap, coldspot_length)
            .expect("nonempty intersection has positive length");
        let b = integer_ratio(overlap, domain_length)
            .expect("nonempty intersection has positive length");
        if a > result.fraction_of_coldspot {
            result.fraction_of_coldspot = a;
        }
        if b > result.fraction_of_domain {
            result.fraction_of_domain = b;
        }
    }
    result.removed =
        result.fraction_of_coldspot > threshold || result.fraction_of_domain > threshold;
    result
}
