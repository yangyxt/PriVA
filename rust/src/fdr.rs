//! One-dimensional float64 Benjamini-Hochberg adjustment used by the PM1
//! producer's scipy.stats.false_discovery_control(..., method="bh") call.
//!
//! Keep division, multiplication, reverse cumulative minimum and clipping as
//! separate operations in the original order. In particular, p * (m / rank)
//! must not become (p * m) / rank: that can change a threshold decision.
//! No multidimensional, float32 or BY-method behavior is implied by this API.

pub(crate) fn benjamini_hochberg(probabilities: &[f64]) -> Result<Vec<f64>, &'static str> {
    if probabilities
        .iter()
        .any(|p| !p.is_finite() || !(0.0..=1.0).contains(p))
    {
        return Err("`ps` must include only numbers between 0 and 1.");
    }
    if probabilities.len() <= 1 {
        // SciPy returns before clipping, so a singleton negative zero stays
        // negative zero. Longer arrays normalize zero during np.clip(ps,0,1).
        return Ok(probabilities.to_vec());
    }
    let mut ordered: Vec<_> = probabilities.iter().copied().enumerate().collect();
    ordered.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    let length = ordered.len() as f64;
    for (at, (_, value)) in ordered.iter_mut().enumerate() {
        *value *= length / (at + 1) as f64;
    }
    for at in (0..ordered.len() - 1).rev() {
        ordered[at].1 = ordered[at].1.min(ordered[at + 1].1);
    }
    let mut adjusted = vec![0.0; ordered.len()];
    for (original, value) in ordered {
        adjusted[original] = if value == 0.0 { 0.0 } else { value.min(1.0) };
    }
    Ok(adjusted)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fdr_keeps_source_order_and_exact_operation_order() {
        assert_eq!(benjamini_hochberg(&[]).unwrap(), Vec::<f64>::new());
        let actual = benjamini_hochberg(&[0.3, 0.2, 0.1]).unwrap();
        assert!(actual.iter().all(|x| x.to_bits() == 0.3_f64.to_bits()));
        assert_eq!(
            benjamini_hochberg(&[-0.0]).unwrap()[0].to_bits(),
            (-0.0_f64).to_bits()
        );
        assert!(
            benjamini_hochberg(&[-0.0, 0.0, -0.0])
                .unwrap()
                .iter()
                .all(|x| x.to_bits() == 0.0_f64.to_bits())
        );
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.1, 1.1] {
            assert_eq!(
                benjamini_hochberg(&[invalid]),
                Err("`ps` must include only numbers between 0 and 1.")
            );
        }
    }
}
