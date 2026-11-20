//! Sampling procedures from `am_pick_intolerant_domains.py` for one-dimensional
//! float64 score arrays. Preserve reference-array order and the caller's PCG64
//! state, including draws consumed before a later input fails.
use crate::{dip_stat, exon_af_stats, ks_greater, numpy_rng::Pcg64};

#[derive(Debug, PartialEq)]
pub(crate) struct Error {
    pub(crate) category: &'static str,
    pub(crate) message: &'static str,
}

fn value_error(message: &'static str) -> Error {
    Error {
        category: "ValueError",
        message,
    }
}

/// The Python procedure validates each reference just before sampling it.
/// In particular, an empty reference still fails when its allocated take is 0.
pub(crate) fn fixed_composite(
    references: &[Vec<f64>],
    target: i64,
    rng: &mut Pcg64,
) -> Result<Vec<f64>, Error> {
    if references.is_empty() {
        return Err(value_error("No reference arrays provided for composite."));
    }
    let count =
        i64::try_from(references.len()).map_err(|_| value_error("too many reference arrays"))?;
    let base = target.div_euclid(count);
    let remainder = target.rem_euclid(count);
    let mut output = Vec::new();
    for (index, values) in references.iter().enumerate() {
        let take = base + i64::from((index as i64) < remainder);
        if values.is_empty() {
            return Err(value_error("One of the reference arrays is empty."));
        }
        let take = usize::try_from(take)
            .map_err(|_| value_error("negative dimensions are not allowed"))?;
        if values.len() >= take {
            for index in rng.choice(values.len(), take, false).map_err(value_error)? {
                output.push(values[index]);
            }
        } else {
            let padding = rng
                .choice(values.len(), take - values.len(), true)
                .map_err(value_error)?;
            output.extend_from_slice(values);
            output.extend(padding.into_iter().map(|index| values[index]));
        }
    }
    Ok(output)
}

#[derive(Debug)]
pub(crate) struct KsResult {
    pub(crate) statistic: f64,
    pub(crate) pvalue: f64,
    pub(crate) n_query: usize,
    pub(crate) n_composite: usize,
    pub(crate) median_query: f64,
    pub(crate) median_composite: f64,
}

fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        // NumPy's empty median reaches its float64 0/0 reduction.
        f64::from_bits(0xfff8_0000_0000_0000)
    } else {
        exon_af_stats::calculate(values)[0]
    }
}

pub(crate) fn ks_sample(
    query: &[f64],
    composite: &[f64],
    rng: &mut Pcg64,
    cap: Option<i64>,
) -> Result<KsResult, Error> {
    if query.len() < 2 {
        return Ok(KsResult {
            statistic: f64::NAN,
            pvalue: f64::NAN,
            n_query: query.len(),
            n_composite: 0,
            median_query: f64::NAN,
            median_composite: f64::NAN,
        });
    }
    let mut take = query.len().max(50);
    if let Some(cap) = cap {
        if cap < 0 && composite.is_empty() {
            return Err(value_error(
                "a must be a positive integer unless no samples are taken",
            ));
        }
        let cap =
            usize::try_from(cap).map_err(|_| value_error("negative dimensions are not allowed"))?;
        take = take.min(cap);
    }
    take = take.min(composite.len());
    let indices = rng
        .choice(composite.len(), take, false)
        .map_err(value_error)?;
    let sample: Vec<_> = indices.into_iter().map(|index| composite[index]).collect();
    let ks = ks_greater::greater(query, &sample);
    Ok(KsResult {
        statistic: ks.statistic,
        pvalue: ks.pvalue,
        n_query: query.len(),
        n_composite: take,
        median_query: median(query),
        median_composite: median(&sample),
    })
}

#[derive(Debug)]
pub(crate) struct DipResult {
    pub(crate) query: f64,
    pub(crate) reference_mean: f64,
    pub(crate) reference_maximum: f64,
    pub(crate) pvalue: f64,
    pub(crate) n_query: usize,
}

/// The complete numerical body of `_dip_rescue_worker`: query first, then
/// draws in permutation order, inclusive `>=`, and the original reductions.
pub(crate) fn dip_rescue(
    query: &[f64],
    composite: &[f64],
    permutations: i64,
    rng: &mut Pcg64,
) -> Result<DipResult, Error> {
    let dip_query = dip_stat::dip(query).map_err(|message| Error {
        category: "RuntimeError",
        message,
    })?;
    let permutations = usize::try_from(permutations)
        .map_err(|_| value_error("negative dimensions are not allowed"))?;
    let mut references = Vec::new();
    for _ in 0..permutations {
        if composite.is_empty() {
            return Err(value_error("a cannot be empty unless no samples are taken"));
        }
        let indices = rng
            .choice(composite.len(), query.len(), false)
            .map_err(value_error)?;
        let sample: Vec<_> = indices.into_iter().map(|index| composite[index]).collect();
        references.push(dip_stat::dip(&sample).map_err(|message| Error {
            category: "RuntimeError",
            message,
        })?);
    }
    if references.is_empty() {
        return Err(value_error(
            "zero-size array to reduction operation maximum which has no identity",
        ));
    }
    let at_least_query = references
        .iter()
        .filter(|&&value| value >= dip_query)
        .count();
    let summaries = exon_af_stats::calculate(&references);
    Ok(DipResult {
        query: dip_query,
        reference_mean: summaries[1],
        reference_maximum: summaries[2],
        pvalue: at_least_query as f64 / references.len() as f64,
        n_query: query.len(),
    })
}
