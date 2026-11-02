//! SciPy 1.14.1's one-dimensional `ks_2samp(alternative="greater")` with
//! the default `method="auto"`, as used by the DAS domain procedure.
//!
//! Adapted from scipy/stats/_stats_py.py and the positive-integer branches
//! of special/{binom.h,cephes/{beta.h,gamma.h}}. Keep the arithmetic order,
//! exact-grid rounding and the original exact-to-asymptotic fallback.
//! This is not a replacement for other KS alternatives or general special
//! functions. See SCIPY_LICENSE for SciPy and Cephes notices.

#[derive(Clone, Copy, Debug)]
pub(crate) struct Outcome {
    pub statistic: f64,
    pub pvalue: f64,
    pub fallback: bool,
}

fn polynomial(x: f64, coefficients: &[f64]) -> f64 {
    coefficients[1..]
        .iter()
        .fold(coefficients[0], |value, &coefficient| {
            value * x + coefficient
        })
}

// Every argument to Gamma in this KS path is a positive integer. The
// rational interpolation on (2,3) and negative-argument reflection are unused.
fn gamma_integer(mut x: f64) -> f64 {
    if x <= 33.0 {
        let mut product = 1.0;
        while x >= 3.0 {
            x -= 1.0;
            product *= x;
        }
        return product;
    }
    let reciprocal = 1.0 / x;
    let correction = 1.0
        + reciprocal
            * polynomial(
                reciprocal,
                &[
                    7.87311395793093628397e-4,
                    -2.29549961613378126380e-4,
                    -2.68132617805781232825e-3,
                    3.47222221605458667310e-3,
                    8.33333333333482257126e-2,
                ],
            );
    let exponential = x.exp();
    let value = if x > 143.01608 {
        let power = x.powf(0.5 * x - 0.25);
        power * (power / exponential)
    } else {
        x.powf(x - 0.5) / exponential
    };
    2.50662827463100050242 * value * correction
}

fn log_gamma_integer(x: f64) -> f64 {
    if x < 13.0 {
        let (mut product, mut shift, mut value) = (1.0_f64, 0.0, x);
        while value >= 3.0 {
            shift -= 1.0;
            value = x + shift;
            product *= value;
        }
        return product.ln();
    }
    let leading = (x - 0.5) * x.ln() - x + 0.91893853320467274178;
    let reciprocal_square = 1.0 / (x * x);
    if x >= 1000.0 {
        let correction = ((7.9365079365079365079365e-4 * reciprocal_square
            - 2.7777777777777777777778e-3)
            * reciprocal_square
            + 0.0833333333333333333333)
            / x;
        leading + correction
    } else {
        leading
            + polynomial(
                reciprocal_square,
                &[
                    8.11614167470508450300e-4,
                    -5.95061904284301438324e-4,
                    7.93650340457716943945e-4,
                    -2.77777777730099687205e-3,
                    8.33333333333331927722e-2,
                ],
            ) / x
    }
}

/// Nonnegative integer binomial coefficients needed by the exact KS path.
/// The automatic exact path bounds n by 20,000; it does not use the general
/// binom asymptotic branches for enormous or noninteger arguments.
pub(crate) fn binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let reduced = k.min(n - k);
    if reduced < 20 {
        let (mut numerator, mut denominator) = (1.0_f64, 1.0_f64);
        for i in 1..=reduced {
            numerator *= i as f64 + n as f64 - reduced as f64;
            denominator *= i as f64;
            if numerator.abs() > 1e50 {
                numerator /= denominator;
                denominator = 1.0;
            }
        }
        return numerator / denominator;
    }
    let (mut a, mut b) = ((1 + n - k) as f64, (1 + k) as f64);
    if a < b {
        std::mem::swap(&mut a, &mut b);
    }
    let mut sum = a + b;
    let beta = if sum > 171.624376956302725 {
        sum = log_gamma_integer(sum);
        sum = log_gamma_integer(b) - sum;
        sum = log_gamma_integer(a) + sum;
        sum.exp()
    } else {
        sum = gamma_integer(sum);
        a = gamma_integer(a);
        b = gamma_integer(b);
        if (a - sum).abs() > (b - sum).abs() {
            (b / sum) * a
        } else {
            (a / sum) * b
        }
    };
    1.0 / (n as f64 + 1.0) / beta
}

fn gcd(mut left: usize, mut right: usize) -> usize {
    while right != 0 {
        (left, right) = (right, left % right);
    }
    left
}

fn outside_paths(mut m: usize, mut n: usize, g: usize, h: usize) -> f64 {
    if m < n {
        std::mem::swap(&mut m, &mut n);
    }
    let (mg, ng) = (m / g, n / g);
    // Python // is floor division, including the negative numerator here.
    let length = (n as i64 + (mg as i64 - h as i64).div_euclid(mg as i64)) as usize;
    if length == 0 {
        return binomial(m + n, n);
    }
    let x: Vec<_> = (0..length).map(|j| (h + mg * j).div_ceil(ng)).collect();
    let mut paths = vec![0.0; length];
    paths[0] = 1.0;
    for j in 1..length {
        let mut value = binomial(x[j] + j, j);
        for i in 0..j {
            value -= binomial(x[j] - x[i] + j - i, j - i) * paths[i];
        }
        if !value.is_finite() {
            return f64::NAN;
        }
        paths[j] = value;
    }
    let mut total = 0.0;
    for j in 0..length {
        total += paths[j] * binomial((m - x[j]) + (n - j), n - j);
        if !total.is_finite() {
            return f64::NAN;
        }
    }
    total
}

fn asymptotic(n1: usize, n2: usize, d: f64) -> f64 {
    let m = n1.max(n2) as f64;
    let n = n1.min(n2) as f64;
    let effective = m * n / (m + n);
    let z = effective.sqrt() * d;
    let exponent = -2.0 * z.powi(2) - 2.0 * z * (m + 2.0 * n) / (m * n * (m + n)).sqrt() / 3.0;
    exponent.exp().clamp(0.0, 1.0)
}

pub(crate) fn greater(first: &[f64], second: &[f64]) -> Outcome {
    if first.is_empty() || second.is_empty() || first.iter().chain(second).any(|x| x.is_nan()) {
        // The public SciPy axis/nan decorator handles these before ks_2samp.
        return Outcome {
            statistic: f64::NAN,
            pvalue: f64::NAN,
            fallback: false,
        };
    }
    let (mut first, mut second) = (first.to_vec(), second.to_vec());
    first.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    second.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    let (n1, n2) = (first.len(), second.len());
    let mut d = f64::NEG_INFINITY;
    for &value in first.iter().chain(&second) {
        let left = first.partition_point(|x| *x <= value) as f64 / n1 as f64;
        let right = second.partition_point(|x| *x <= value) as f64 / n2 as f64;
        d = d.max(left - right);
    }
    if n1.max(n2) > 10_000 {
        return Outcome {
            statistic: d,
            pvalue: asymptotic(n1, n2, d),
            fallback: false,
        };
    }
    let divisor = gcd(n1, n2);
    let lcm = (n1 / divisor) * n2;
    let h = (d * lcm as f64).round_ties_even() as usize;
    d = h as f64 * 1.0 / lcm as f64;
    let probability = if h == 0 {
        1.0
    } else if n1 == n2 {
        let mut product = 1.0;
        for j in 0..h {
            product *= (n1 - j) as f64 / ((n1 + j) as f64 + 1.0);
        }
        product
    } else {
        let paths = outside_paths(n1, n2, divisor, h);
        let all = binomial(n1 + n2, n1);
        if paths > all || all.is_infinite() {
            f64::NAN
        } else {
            paths / all
        }
    };
    if (0.0..=1.0).contains(&probability) {
        Outcome {
            statistic: d,
            pvalue: probability.clamp(0.0, 1.0),
            fallback: false,
        }
    } else {
        Outcome {
            statistic: d,
            pvalue: asymptotic(n1, n2, d),
            fallback: true,
        }
    }
}
