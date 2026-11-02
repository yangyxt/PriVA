//! Hartigan's dip statistic, implemented from the convex-envelope procedure
//! described in P. M. Hartigan (1985), JRSS C 34:320-325,
//! DOI 10.2307/2347485, and Hartigan & Hartigan (1985), Annals of Statistics
//! 13:70-84, DOI 10.1214/aos/1176346577.
//!
//! This implementation constructs monotone hulls afresh within each shrinking
//! modal interval. It does not incorporate the published Fortran listing or
//! the GPL diptest implementation. The installed diptest 0.11.0 package is an
//! independent development oracle only. Its allow_zero=True convention is
//! the convention required by the DAS caller.
//!
//! The pinned Linux reference chooses intervals using x87 extended precision
//! for the interpolation numerator, with a binary64 denominator. Final dip
//! deviations remain binary64. Use software arithmetic for that mixed profile
//! so host `long double` widths cannot change it. See RUSTC_APFLOAT_LICENSE.
use rustc_apfloat::{
    Float, FloatConvert,
    ieee::{Double, X87DoubleExtended},
};

fn extended(value: f64) -> X87DoubleExtended {
    Double::from_bits(value.to_bits() as u128)
        .convert(&mut false)
        .value
}

fn integer(value: usize) -> X87DoubleExtended {
    X87DoubleExtended::from_u128(value as u128).value
}

fn projection(values: &[f64], i: usize, a: usize, b: usize) -> X87DoubleExtended {
    let numerator = (extended(values[i]) - extended(values[a])).value;
    let numerator = (numerator * integer(b - a)).value;
    (numerator / extended(values[b] - values[a])).value
}

fn envelope(values: &[f64], low: usize, high: usize, convex: bool) -> Vec<usize> {
    let mut vertices: Vec<usize> = Vec::new();
    for next in low..=high {
        while vertices.len() >= 2 {
            let a = vertices[vertices.len() - 2];
            let b = vertices[vertices.len() - 1];
            let left = (values[b] - values[a]) * (next - b) as f64;
            let right = (values[next] - values[b]) * (b - a) as f64;
            if if convex { left > right } else { left < right } {
                break;
            }
            vertices.pop();
        }
        vertices.push(next);
    }
    vertices
}

fn deviation(values: &[f64], vertices: &[usize], convex: bool) -> f64 {
    let mut largest = 0.0_f64;
    for segment in vertices.windows(2) {
        let (a, b) = (segment[0], segment[1]);
        // Each empirical step contributes one count even when a segment is
        // constant, adjacent, or its binary64 slope overflows.
        largest = largest.max(1.0);
        if b - a <= 1 || values[b] == values[a] {
            continue;
        }
        let slope = (b - a) as f64 / (values[b] - values[a]);
        for i in a..=b {
            let distance = if convex {
                (i - a + 1) as f64 - (values[i] - values[a]) * slope
            } else {
                (values[i] - values[a]) * slope - ((i - a) as f64 - 1.0)
            };
            if distance > largest {
                largest = distance;
            }
        }
    }
    largest
}

pub(crate) fn dip(values: &[f64]) -> Result<f64, &'static str> {
    if values.is_empty() {
        return Err("N must be >= 1.");
    }
    let mut sorted = values.to_vec();
    // NumPy sorts NaNs after finite values. Their payloads do not enter the
    // real-valued envelope computation; no sample is silently dropped.
    sorted.sort_unstable_by(|a, b| {
        a.partial_cmp(b)
            .unwrap_or_else(|| a.is_nan().cmp(&b.is_nan()))
    });
    let (mut low, mut high) = (0, sorted.len() - 1);
    if !(sorted[high] > sorted[low]) {
        return Ok(0.0);
    }
    let mut twice_dip = 0.0_f64;
    while low < high {
        let lower = envelope(&sorted, low, high, true);
        let upper = envelope(&sorted, low, high, false);
        let (mut new_low, mut new_high) = (low, high);
        let mut separation = X87DoubleExtended::ZERO;
        // The maximum gap between polygonal envelopes occurs at a vertex.
        // Visit the merged contact sequence from left to right and retain
        // the last equal gap and its bracketing endpoints (AS 217, steps 3–4).
        let (mut g, mut h) = (1, 1);
        while lower[g] != upper[h] {
            let (gap, candidate_low, candidate_high) = if lower[g] < upper[h] {
                let (i, a, b) = (lower[g], upper[h - 1], upper[h]);
                let projected = projection(&sorted, i, a, b);
                g = (g + 1).min(lower.len() - 1);
                ((projected - integer(i - a - 1)).value, i, b)
            } else {
                let (i, a, b) = (upper[h], lower[g - 1], lower[g]);
                let projected = projection(&sorted, i, a, b);
                h = (h + 1).min(upper.len() - 1);
                ((integer(i - a + 1) - projected).value, a, i)
            };
            if gap >= separation {
                separation = gap;
                new_low = candidate_low;
                new_high = candidate_high;
            }
        }
        // Equal bounds may still shrink the interval; unchanged endpoints
        // below also terminate the iteration.
        if separation < extended(twice_dip) {
            break;
        }
        let lower_end = lower.partition_point(|&i| i <= new_low);
        let upper_start = upper.partition_point(|&i| i < new_high);
        twice_dip = twice_dip.max(deviation(&sorted, &lower[..lower_end], true));
        twice_dip = twice_dip.max(deviation(&sorted, &upper[upper_start..], false));
        if (new_low, new_high) == (low, high) {
            break;
        }
        (low, high) = (new_low, new_high);
    }
    Ok(twice_dip / (2.0 * values.len() as f64))
}

#[cfg(test)]
mod tests {
    use super::dip;

    #[test]
    fn tied_interval_uses_reference_extended_precision() {
        let values = [0.1, 0.3, 0.4, 0.5, 0.7, 0.8, 0.9];
        assert_eq!(dip(&values).unwrap().to_bits(), 0x3fbb_6db6_db6d_b6d9);
    }

    #[test]
    fn empirical_step_survives_an_overflowing_slope() {
        assert_eq!(dip(&[0.0, -0.0, f64::from_bits(1)]).unwrap(), 1.0 / 6.0);
    }

    #[test]
    fn empty_and_constant_samples_follow_allow_zero_profile() {
        assert_eq!(dip(&[]), Err("N must be >= 1."));
        assert_eq!(dip(&[1.0; 10]).unwrap().to_bits(), 0);
    }
}
