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

/// The two outputs used by DAS from statsmodels 0.14.4's
/// `multipletests(..., method="fdr_bh")`. Its p/(rank/n) operation order,
/// rejection boundary and nonfinite propagation differ from the SciPy API
/// above. The unused Sidak/Bonferroni alpha outputs are outside this API.
/// Adapted from statsmodels/stats/multitest.py; see STATSMODELS_LICENSE.
#[allow(dead_code)] // Qualified independently before the opt-in DAS caller is added.
pub(crate) fn statsmodels_bh(
    probabilities: &[f64],
    alpha: f64,
) -> Result<(Vec<bool>, Vec<f64>), &'static str> {
    if probabilities.is_empty() {
        return Err("float division by zero");
    }
    let order = statsmodels_order(probabilities);
    let length = probabilities.len() as f64;
    let mut last_rejected = None;
    let mut values = Vec::with_capacity(order.len());
    for (at, &index) in order.iter().enumerate() {
        let fraction = (at + 1) as f64 / length;
        if probabilities[index] <= fraction * alpha {
            last_rejected = Some(at);
        }
        values.push(probabilities[index] / fraction);
    }
    for at in (0..values.len() - 1).rev() {
        // np.minimum.accumulate on the reversed array keeps an encountered
        // NaN, and chooses the new operand on an equal (including zero) tie.
        let accumulated = values[at + 1];
        if accumulated.is_nan() || accumulated < values[at] {
            values[at] = accumulated;
        }
    }
    let mut adjusted = vec![0.0; values.len()];
    let mut rejected = vec![false; values.len()];
    for (at, (&index, value)) in order.iter().zip(values).enumerate() {
        adjusted[index] = if value > 1.0 { 1.0 } else { value };
        rejected[index] = last_rejected.is_some_and(|last| at <= last);
    }
    Ok((rejected, adjusted))
}

// NumPy 1.26.4 scalar float64 argsort traversal, needed to select the same
// NaN payload when several payloads compare equal. Adapted from
// numpy/core/src/npysort/{quicksort.cpp,npysort_heapsort.h}; NUMPY_LICENSE.
fn statsmodels_order(values: &[f64]) -> Vec<usize> {
    let less =
        |a: usize, b: usize| values[a] < values[b] || (!values[a].is_nan() && values[b].is_nan());
    let mut order: Vec<_> = (0..values.len()).collect();
    if order.len() < 2 {
        return order;
    }
    let mut stack = vec![(0, order.len() - 1, 2 * order.len().ilog2() as i32)];
    while let Some((mut low, mut high, mut depth)) = stack.pop() {
        if depth < 0 {
            let slice = &mut order[low..=high];
            let sift = |slice: &mut [usize], mut root: usize, end: usize| {
                let value = slice[root];
                while root * 2 + 1 < end {
                    let mut child = root * 2 + 1;
                    if child + 1 < end && less(slice[child], slice[child + 1]) {
                        child += 1;
                    }
                    if !less(value, slice[child]) {
                        break;
                    }
                    slice[root] = slice[child];
                    root = child;
                }
                slice[root] = value;
            };
            let length = slice.len();
            for root in (0..length / 2).rev() {
                sift(slice, root, length);
            }
            for end in (1..length).rev() {
                slice.swap(0, end);
                sift(slice, 0, end);
            }
            continue;
        }
        while high - low > 15 {
            let middle = low + (high - low) / 2;
            for (left, right) in [(low, middle), (middle, high), (low, middle)] {
                if less(order[right], order[left]) {
                    order.swap(left, right);
                }
            }
            let pivot = order[middle];
            order.swap(middle, high - 1);
            let (mut left, mut right) = (low, high - 1);
            loop {
                left += 1;
                while less(order[left], pivot) {
                    left += 1;
                }
                right -= 1;
                while less(pivot, order[right]) {
                    right -= 1;
                }
                if left >= right {
                    break;
                }
                order.swap(left, right);
            }
            order.swap(left, high - 1);
            depth -= 1;
            if left - low < high - left {
                stack.push((left + 1, high, depth));
                high = left - 1;
            } else {
                stack.push((low, left - 1, depth));
                low = left + 1;
            }
        }
        for index in low + 1..=high {
            let value = order[index];
            let mut at = index;
            while at > low && less(value, order[at - 1]) {
                order[at] = order[at - 1];
                at -= 1;
            }
            order[at] = value;
        }
    }
    order
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
