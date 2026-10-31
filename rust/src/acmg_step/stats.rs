//! The four SciPy 1.14.1 functions the ACMG step calls. Thread T3 (`acmg_scoring`) and thread T6
//! (`acmg_pm1_pp2_bp1_missense`) call them by these signatures.
//!
//! The owner's decision (2026-10-30): every result that only feeds a threshold test must give
//! the same decision as SciPy; `beta_ppf` may differ from SciPy in the last printed digit (it
//! becomes the printed `control_common_index`).
//!
//! Callers (current line numbers of the Python scripts):
//! - `binom.sf(k, n, p)`: `acmg_pm2_bs1_ba1_frequency.py` 145 and 643;
//! - `binom.cdf(k, n, p)`: `acmg_pm2_bs1_ba1_frequency.py` 240, `acmg_scoring.py` 431 and 439;
//! - `beta.ppf(q, a, b)`: `acmg_scoring.py` 459;
//! - `chi2_contingency(table)[0:2]` and `fisher_exact(table, alternative="greater")`:
//!   `acmg_pm1_pp2_bp1_missense.py` 206 and 212, on 2x2 tables of counts.
//!
//! How the values are made. SciPy 1.14.1 (the conda-forge build in `priva_acmg`, GCC 13,
//! x86-64 baseline without FMA) computes the binomial tails and the hypergeometric cdf with the
//! Boost.Math it vendors (commit a53b013 of boostorg/math), in double precision (its `StatsPolicy`
//! turns off promotion to long double), and the chi-square tail with its C++ copy of Cephes
//! (`igamc`). This module is a line-by-line port of exactly the code paths those functions take for
//! the arguments the step passes (integer counts), keeping Boost's and Cephes's order of
//! floating-point operations and their constants: Boost's own rational `log1p` and `expm1` (the
//! standalone build does not use the C library's), the SSE2 specialisation of the Lanczos sum that
//! x86-64 compiles, the factorial table, and the C library's `exp`, `log`, `pow` and `sqrt` (the same
//! glibc on the same node). So the values are SciPy's values bit for bit, not only the decisions;
//! the comparison harness measures this on millions of inputs.
//!
//! `beta.ppf` is the exception the owner accepted: SciPy runs Boost's `ibeta_inv` in x87 long
//! double, which after rounding is (up to its own last-bit error) the correctly rounded quantile.
//! Here the quantile is found by safeguarded Newton steps on the ported incomplete beta; for
//! integer shape parameters (allele counts) it is then refined to the correctly rounded double
//! with the incomplete beta summed as a binomial tail in double-double arithmetic.
//!
//! Every function takes one element; the callers map them over columns. NaN and out-of-domain
//! arguments give what SciPy's distribution wrappers give (NaN, 0 or 1).

use super::types::HarnessCall;
use crate::cli::{Result, invalid};
use crate::json_text::Json;
use crate::pandas_table::{Column, Table};
use std::sync::OnceLock;

// ------------------------------------------------------------------------------------------------
// Public functions (SciPy's distribution wrappers around the Boost and Cephes cores)

/// `scipy.stats.binom.sf(k, n, p)`: P(X > k) for X ~ Binomial(n, p).
pub(crate) fn binom_sf(k: f64, n: f64, p: f64) -> f64 {
    if k.is_nan() || !binom_args_valid(n, p) {
        return f64::NAN;
    }
    // rv_discrete.sf: below the support 1, at or above n 0, inside Boost's complement cdf.
    if k < 0.0 {
        return 1.0;
    }
    if !(k < n) {
        return 0.0;
    }
    clip01(boost_binomial_cdf(k.floor(), n, p, true))
}

/// `scipy.stats.binom.cdf(k, n, p)`: P(X <= k) for X ~ Binomial(n, p).
pub(crate) fn binom_cdf(k: f64, n: f64, p: f64) -> f64 {
    if k.is_nan() || !binom_args_valid(n, p) {
        return f64::NAN;
    }
    // rv_discrete.cdf: at or above n 1, below the support 0, inside Boost's cdf.
    if k >= n {
        return 1.0;
    }
    if k < 0.0 {
        return 0.0;
    }
    clip01(boost_binomial_cdf(k.floor(), n, p, false))
}

/// `scipy.stats.beta.ppf(q, a, b)`: the q-quantile of Beta(a, b).
pub(crate) fn beta_ppf(q: f64, a: f64, b: f64) -> f64 {
    // rv_continuous.ppf: shape parameters must be positive, q in [0, 1].
    if q.is_nan() || a.is_nan() || b.is_nan() || !(a > 0.0 && b > 0.0) {
        return f64::NAN;
    }
    if q == 0.0 {
        return 0.0;
    }
    if q == 1.0 {
        return 1.0;
    }
    if !(q > 0.0 && q < 1.0) || !a.is_finite() || !b.is_finite() {
        return f64::NAN;
    }
    let estimate = ibeta_quantile(a, b, q);
    // Integer shape parameters (allele counts): refine to the correctly rounded quantile with
    // the incomplete beta as a binomial tail summed in double-double arithmetic, which is what
    // SciPy's long-double computation comes to after rounding, up to its own last-bit error.
    let integral = |value: f64| value.fract() == 0.0;
    if integral(a) && integral(b) && a + b <= 4e6 && estimate > 0.0 && estimate < 1.0 {
        if let Some(refined) = correctly_rounded_quantile(a as u64, b as u64, q, estimate) {
            return refined;
        }
    }
    estimate
}

/// `scipy.stats.chi2_contingency(table)[0:2]` for a 2x2 table with the default Yates
/// continuity correction: `(statistic, p-value)`. A table with a negative cell or a zero
/// margin (where SciPy raises ValueError) gives `(NaN, NaN)`; the callers never build one.
pub(crate) fn chi2_contingency_2x2(table: [[i64; 2]; 2]) -> (f64, f64) {
    if table.iter().flatten().any(|cell| *cell < 0) {
        return (f64::NAN, f64::NAN);
    }
    // expected_freq: margins of the table as float64, outer product over the total.
    let observed = table.map(|row| row.map(|cell| cell as f64));
    let rows = [
        observed[0][0] + observed[0][1],
        observed[1][0] + observed[1][1],
    ];
    let columns = [
        observed[0][0] + observed[1][0],
        observed[0][1] + observed[1][1],
    ];
    let total = ((observed[0][0] + observed[0][1]) + observed[1][0]) + observed[1][1];
    let mut expected = [[0.0; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            expected[i][j] = rows[i] * columns[j] / total;
            if expected[i][j] == 0.0 || expected[i][j].is_nan() {
                return (f64::NAN, f64::NAN);
            }
        }
    }
    // Yates: move each observed count half a unit (at most) towards its expected count.
    let mut corrected = [[0.0; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            let diff = expected[i][j] - observed[i][j];
            let direction = if diff > 0.0 {
                1.0
            } else if diff < 0.0 {
                -1.0
            } else {
                0.0
            };
            let magnitude = diff.abs().min(0.5);
            corrected[i][j] = observed[i][j] + magnitude * direction;
        }
    }
    // power_divergence (Pearson): numpy sums the four terms in order.
    let mut statistic = 0.0;
    for i in 0..2 {
        for j in 0..2 {
            let delta = corrected[i][j] - expected[i][j];
            statistic += delta * delta / expected[i][j];
        }
    }
    (statistic, cephes_chdtrc_1(statistic))
}

/// `scipy.stats.fisher_exact(table, alternative="greater")` for a 2x2 table:
/// `(odds ratio, p-value)`. A negative cell (where SciPy raises ValueError) gives `(NaN, NaN)`.
pub(crate) fn fisher_exact_greater(table: [[i64; 2]; 2]) -> (f64, f64) {
    let c = table;
    if c.iter().flatten().any(|cell| *cell < 0) {
        return (f64::NAN, f64::NAN);
    }
    let column_sums = [c[0][0] + c[1][0], c[0][1] + c[1][1]];
    let row_sums = [c[0][0] + c[0][1], c[1][0] + c[1][1]];
    if column_sums.contains(&0) || row_sums.contains(&0) {
        return (f64::NAN, 1.0);
    }
    let odds_ratio = if c[1][0] > 0 && c[0][1] > 0 {
        (c[0][0] * c[1][1]) as f64 / (c[1][0] * c[0][1]) as f64
    } else {
        f64::INFINITY
    };
    let n1 = row_sums[0];
    let n2 = row_sums[1];
    // hypergeom.cdf(c[0, 1], n1 + n2, n1, c[0, 1] + c[1, 1])
    let pvalue = hypergeom_cdf(c[0][1], n1 + n2, n1, c[0][1] + c[1][1]);
    (odds_ratio, pvalue.min(1.0))
}

// ------------------------------------------------------------------------------------------------
// SciPy's wrappers

/// `binom._argcheck`: n >= 0 and integral (`n == np.round(n)`, true for infinity), p in [0, 1].
fn binom_args_valid(n: f64, p: f64) -> bool {
    !n.is_nan() && !p.is_nan() && n >= 0.0 && n == n.round_ties_even() && (0.0..=1.0).contains(&p)
}

/// `np.clip(value, 0, 1)` (NaN stays NaN).
fn clip01(value: f64) -> f64 {
    if value.is_nan() {
        value
    } else {
        value.clamp(0.0, 1.0)
    }
}

/// `scipy.stats.hypergeom.cdf(k, M, n, N)` for integers (rv_discrete.cdf around Boost's
/// `cdf(hypergeometric_distribution(r = n, n = N, N = M), k)`).
fn hypergeom_cdf(k: i64, total: i64, good: i64, draws: i64) -> f64 {
    if !(total > 0 && good >= 0 && draws >= 0 && good <= total && draws <= total) {
        return f64::NAN;
    }
    let low = (draws - (total - good)).max(0);
    let high = good.min(draws);
    if k >= high {
        return 1.0;
    }
    if k < low {
        return 0.0;
    }
    clip01(boost_hypergeometric_cdf(
        k as u64,
        good as u64,
        draws as u64,
        total as u64,
    ))
}

// ------------------------------------------------------------------------------------------------
// Boost.Math: constants and elementary pieces (double, 53-bit tag)

const EPSILON: f64 = f64::EPSILON;
const MIN_VALUE: f64 = f64::MIN_POSITIVE;
const MAX_VALUE: f64 = f64::MAX;
const LOG_MAX_VALUE: f64 = 709.0;
const LOG_MIN_VALUE: f64 = -708.0;
const MAX_SERIES_ITERATIONS: u64 = 1_000_000;

/// `lanczos13m53::g()`.
const LANCZOS_G: f64 = 6.024680040776729583740234375;

/// `lanczos13m53::lanczos_sum_expG_scaled` coefficients (the long double literals round to these
/// doubles, checked by exact rational rounding through the 64-bit significand).
#[allow(clippy::excessive_precision)]
const LANCZOS_EXPG_SCALED_NUM: [f64; 13] = [
    56906521.91347156388090791033559122686859,
    103794043.1163445451906271053616070238554,
    86363131.28813859145546927288977868422342,
    43338889.32467613834773723740590533316085,
    14605578.08768506808414169982791359218571,
    3481712.15498064590882071018964774556468,
    601859.6171681098786670226533699352302507,
    75999.29304014542649875303443598909137092,
    6955.999602515376140356310115515198987526,
    449.9445569063168119446858607650988409623,
    19.51992788247617482847860966235652136208,
    0.5098416655656676188125178644804694509993,
    0.006061842346248906525783753964555936883222,
];
const LANCZOS_DENOM: [f64; 13] = [
    0.0,
    39916800.0,
    120543840.0,
    150917976.0,
    105258076.0,
    45995730.0,
    13339535.0,
    2637558.0,
    357423.0,
    32670.0,
    1925.0,
    66.0,
    1.0,
];

/// `lanczos13m53::lanczos_sum_expG_scaled<double>` as compiled on x86-64: Boost's SSE2
/// specialization (`detail/lanczos_sse2.hpp`), which runs the numerator and the denominator in
/// the two lanes of one second-order Horner scheme in x (not in 1/x as the generic rational
/// evaluation does for x > 1), and plain Horner in 1/x only above 4.76886e25.
fn lanczos_sum_expg_scaled(x: f64) -> f64 {
    let (num, den) = (&LANCZOS_EXPG_SCALED_NUM, &LANCZOS_DENOM);
    if x > 4.76886e+25 {
        let z = 1.0 / x;
        let horner = |c: &[f64; 13]| {
            c[1..]
                .iter()
                .fold(c[0], |sum, coefficient| sum * z + coefficient)
        };
        return horner(num) / horner(den);
    }
    let x2 = x * x;
    let lane = |c: &[f64; 13]| {
        let mut even = c[12];
        let mut odd = c[11];
        for (e, o) in [(10, 9), (8, 7), (6, 5), (4, 3), (2, 1)] {
            even = even * x2 + c[e];
            odd = odd * x2 + c[o];
        }
        even = even * x2 + c[0];
        odd *= x;
        even + odd
    };
    lane(num) / lane(den)
}

/// `tools::evaluate_polynomial` for 8 terms (`polynomial_horner3_8.hpp`).
fn evaluate_polynomial_8(a: &[f64; 8], x: f64) -> f64 {
    let x2 = x * x;
    let mut t0 = a[7] * x2 + a[5];
    let mut t1 = a[6] * x2 + a[4];
    t0 *= x2;
    t1 *= x2;
    t0 += a[3];
    t1 += a[2];
    t0 *= x2;
    t1 *= x2;
    t0 += a[1];
    t1 += a[0];
    t0 *= x;
    t0 + t1
}

/// `tools::evaluate_polynomial` for 6 terms (`polynomial_horner3_6.hpp`).
fn evaluate_polynomial_6(a: &[f64; 6], x: f64) -> f64 {
    let x2 = x * x;
    let mut t0 = a[5] * x2 + a[3];
    let mut t1 = a[4] * x2 + a[2];
    t0 *= x2;
    t1 *= x2;
    t0 += a[1];
    t1 += a[0];
    t0 *= x;
    t0 + t1
}

/// `boost::math::log1p` for double (`log1p_imp`, 53-bit tag; the standalone build has no
/// `BOOST_HAS_LOG1P`, so the C library's log1p is not used).
#[allow(clippy::excessive_precision)]
fn boost_log1p(x: f64) -> f64 {
    const P: [f64; 8] = [
        0.15141069795941984e-16,
        0.35495104378055055e-15,
        0.33333333333332835,
        0.99249063543365859,
        1.1143969784156509,
        0.58052937949269651,
        0.13703234928513215,
        0.011294864812099712,
    ];
    const Q: [f64; 8] = [
        1.0,
        3.7274719063011499,
        5.5387948649720334,
        4.159201143419005,
        1.6423855110312755,
        0.31706251443180914,
        0.022665554431410243,
        -0.29252538135177773e-5,
    ];
    if x < -1.0 || x.is_nan() {
        return f64::NAN;
    }
    if x == -1.0 {
        return f64::NEG_INFINITY;
    }
    let a = x.abs();
    if a > 0.5 {
        return (1.0 + x).ln();
    }
    if a < EPSILON {
        return x;
    }
    let result = 1.0 - x / 2.0 + evaluate_polynomial_8(&P, x) / evaluate_polynomial_8(&Q, x);
    result * x
}

/// `boost::math::expm1` for double (`expm1_imp`, 53-bit tag).
#[allow(clippy::excessive_precision)]
fn boost_expm1(x: f64) -> f64 {
    const N: [f64; 6] = [
        -0.28127670288085937e-1,
        0.51278186299064534e0,
        -0.6310029069350198e-1,
        0.11638457975729296e-1,
        -0.52143390687521003e-3,
        0.21491399776965688e-4,
    ];
    const D: [f64; 6] = [
        1.0,
        -0.45442309511354755e0,
        0.90850389570911714e-1,
        -0.10088963629815502e-1,
        0.63003407478692265e-3,
        -0.17976570003654402e-4,
    ];
    let a = x.abs();
    if a > 0.5 {
        if a >= LOG_MAX_VALUE {
            return if x > 0.0 { f64::INFINITY } else { -1.0 };
        }
        return x.exp() - 1.0;
    }
    if a < EPSILON {
        return x;
    }
    let y = 0.10281276702880859e1_f32 as f64;
    x * y + x * evaluate_polynomial_6(&N, x) / evaluate_polynomial_6(&D, x)
}

/// `boost::math::powm1(x, y)` for x >= 0.
fn boost_powm1(x: f64, y: f64) -> f64 {
    if x > 0.0 && ((y * (x - 1.0)).abs() < 0.5 || y.abs() < 0.2) {
        let l = y * x.ln();
        if l < 0.5 {
            return boost_expm1(l);
        }
        if l > LOG_MAX_VALUE {
            return f64::INFINITY;
        }
    }
    x.powf(y) - 1.0
}

/// `itrunc`: truncation towards zero to int.
fn itrunc(value: f64) -> i32 {
    value.trunc() as i32
}

/// `unchecked_factorial<double>(i)`: Boost's table (every entry is the double nearest i!).
fn unchecked_factorial(i: u64) -> f64 {
    FACTORIALS[i as usize]
}

const MAX_FACTORIAL: u64 = 170;

#[rustfmt::skip]
#[allow(clippy::excessive_precision)]
const FACTORIALS: [f64; 171] = [
    1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0, 39916800.0,
    479001600.0, 6227020800.0, 87178291200.0, 1307674368000.0, 20922789888000.0,
    355687428096000.0, 6402373705728000.0, 1.21645100408832e+17, 2.43290200817664e+18,
    5.109094217170944e+19, 1.1240007277776077e+21, 2.585201673888498e+22, 6.204484017332394e+23,
    1.5511210043330986e+25, 4.0329146112660565e+26, 1.0888869450418352e+28,
    3.0488834461171387e+29, 8.841761993739702e+30, 2.6525285981219107e+32, 8.222838654177922e+33,
    2.631308369336935e+35, 8.683317618811886e+36, 2.9523279903960416e+38, 1.0333147966386145e+40,
    3.7199332678990125e+41, 1.3763753091226346e+43, 5.230226174666011e+44, 2.0397882081197444e+46,
    8.159152832478977e+47, 3.345252661316381e+49, 1.40500611775288e+51, 6.041526306337383e+52,
    2.658271574788449e+54, 1.1962222086548019e+56, 5.502622159812089e+57, 2.5862324151116818e+59,
    1.2413915592536073e+61, 6.082818640342675e+62, 3.0414093201713376e+64, 1.5511187532873822e+66,
    8.065817517094388e+67, 4.2748832840600255e+69, 2.308436973392414e+71, 1.2696403353658276e+73,
    7.109985878048635e+74, 4.0526919504877214e+76, 2.3505613312828785e+78, 1.3868311854568984e+80,
    8.32098711274139e+81, 5.075802138772248e+83, 3.146997326038794e+85, 1.98260831540444e+87,
    1.2688693218588417e+89, 8.247650592082472e+90, 5.443449390774431e+92, 3.647111091818868e+94,
    2.4800355424368305e+96, 1.711224524281413e+98, 1.1978571669969892e+100,
    8.504785885678623e+101, 6.1234458376886085e+103, 4.4701154615126844e+105,
    3.307885441519386e+107, 2.48091408113954e+109, 1.8854947016660504e+111,
    1.4518309202828587e+113, 1.1324281178206297e+115, 8.946182130782976e+116,
    7.156945704626381e+118, 5.797126020747368e+120, 4.753643337012842e+122,
    3.945523969720659e+124, 3.314240134565353e+126, 2.81710411438055e+128,
    2.4227095383672734e+130, 2.107757298379528e+132, 1.8548264225739844e+134,
    1.650795516090846e+136, 1.4857159644817615e+138, 1.352001527678403e+140,
    1.2438414054641308e+142, 1.1567725070816416e+144, 1.087366156656743e+146,
    1.032997848823906e+148, 9.916779348709496e+149, 9.619275968248212e+151,
    9.426890448883248e+153, 9.332621544394415e+155, 9.332621544394415e+157, 9.42594775983836e+159,
    9.614466715035127e+161, 9.90290071648618e+163, 1.0299016745145628e+166,
    1.081396758240291e+168, 1.1462805637347084e+170, 1.226520203196138e+172,
    1.324641819451829e+174, 1.4438595832024937e+176, 1.588245541522743e+178,
    1.7629525510902446e+180, 1.974506857221074e+182, 2.2311927486598138e+184,
    2.5435597334721877e+186, 2.925093693493016e+188, 3.393108684451898e+190,
    3.969937160808721e+192, 4.684525849754291e+194, 5.574585761207606e+196,
    6.689502913449127e+198, 8.094298525273444e+200, 9.875044200833601e+202,
    1.214630436702533e+205, 1.506141741511141e+207, 1.882677176888926e+209,
    2.372173242880047e+211, 3.0126600184576594e+213, 3.856204823625804e+215,
    4.974504222477287e+217, 6.466855489220474e+219, 8.47158069087882e+221,
    1.1182486511960043e+224, 1.4872707060906857e+226, 1.9929427461615188e+228,
    2.6904727073180504e+230, 3.659042881952549e+232, 5.012888748274992e+234,
    6.917786472619489e+236, 9.615723196941089e+238, 1.3462012475717526e+241,
    1.898143759076171e+243, 2.695364137888163e+245, 3.854370717180073e+247,
    5.5502938327393044e+249, 8.047926057471992e+251, 1.1749972043909107e+254,
    1.727245890454639e+256, 2.5563239178728654e+258, 3.80892263763057e+260,
    5.713383956445855e+262, 8.62720977423324e+264, 1.3113358856834524e+267,
    2.0063439050956823e+269, 3.0897696138473508e+271, 4.789142901463394e+273,
    7.471062926282894e+275, 1.1729568794264145e+278, 1.853271869493735e+280,
    2.9467022724950384e+282, 4.7147236359920616e+284, 7.590705053947219e+286,
    1.2296942187394494e+289, 2.0044015765453026e+291, 3.287218585534296e+293,
    5.423910666131589e+295, 9.003691705778438e+297, 1.503616514864999e+300,
    2.5260757449731984e+302, 4.269068009004705e+304, 7.257415615307999e+306,
];

// ------------------------------------------------------------------------------------------------
// Boost.Math: the beta function and the incomplete beta (beta.hpp)

/// `beta_imp` with the Lanczos approximation, for a, b > 0.
fn boost_beta(a: f64, b: f64) -> f64 {
    let (mut a, mut b) = (a, b);
    if a <= 0.0 || b <= 0.0 {
        return f64::NAN;
    }
    let c = a + b;
    if c == a && b < EPSILON {
        return 1.0 / b;
    } else if c == b && a < EPSILON {
        return 1.0 / a;
    }
    if b == 1.0 {
        return 1.0 / a;
    } else if a == 1.0 {
        return 1.0 / b;
    } else if c < EPSILON {
        let mut result = c / a;
        result /= b;
        return result;
    }
    if a < b {
        std::mem::swap(&mut a, &mut b);
    }
    let agh = a + LANCZOS_G - 0.5;
    let bgh = b + LANCZOS_G - 0.5;
    let cgh = c + LANCZOS_G - 0.5;
    let mut result =
        lanczos_sum_expg_scaled(a) * (lanczos_sum_expg_scaled(b) / lanczos_sum_expg_scaled(c));
    let ambh = a - 0.5 - b;
    if (b * ambh).abs() < cgh * 100.0 && a > 100.0 {
        result *= (ambh * boost_log1p(-b / cgh)).exp();
    } else {
        result *= (agh / cgh).powf(a - 0.5 - b);
    }
    if cgh > 1e10_f32 as f64 {
        result *= ((agh / cgh) * (bgh / cgh)).powf(b);
    } else {
        result *= ((agh * bgh) / (cgh * cgh)).powf(b);
    }
    result *= (std::f64::consts::E / bgh).sqrt();
    result
}

/// `binomial_coefficient<double>(n, k)`.
fn boost_binomial_coefficient(n: u32, k: u32) -> f64 {
    if k > n {
        return f64::NAN;
    }
    if k == 0 || k == n {
        return 1.0;
    }
    if k == 1 || k == n - 1 {
        return n as f64;
    }
    let mut result;
    if u64::from(n) <= MAX_FACTORIAL {
        result = unchecked_factorial(u64::from(n));
        result /= unchecked_factorial(u64::from(n - k));
        result /= unchecked_factorial(u64::from(k));
    } else {
        if k < n - k {
            result = k as f64 * boost_beta(k as f64, (n - k + 1) as f64);
        } else {
            result = (n - k) as f64 * boost_beta((k + 1) as f64, (n - k) as f64);
        }
        if result == 0.0 {
            // SciPy's overflow handler returns 0 (and sets a Python OverflowError).
            return 0.0;
        }
        result = 1.0 / result;
    }
    (result - 0.5).ceil()
}

/// `ibeta_power_terms` (Lanczos, normalised, prefix 1): x^a y^b / B(a, b).
fn ibeta_power_terms(a: f64, b: f64, x: f64, y: f64) -> f64 {
    let c = a + b;
    let agh = a + LANCZOS_G - 0.5;
    let bgh = b + LANCZOS_G - 0.5;
    let cgh = c + LANCZOS_G - 0.5;
    let mut result = if a < MIN_VALUE || b < MIN_VALUE {
        0.0
    } else {
        lanczos_sum_expg_scaled(c) / (lanczos_sum_expg_scaled(a) * lanczos_sum_expg_scaled(b))
    };
    result *= 1.0;
    result *= (bgh / std::f64::consts::E).sqrt();
    result *= (agh / cgh).sqrt();
    let mut l1 = (x * b - y * agh) / agh;
    let mut l2 = (y * a - x * bgh) / bgh;
    if l1.abs().min(l2.abs()) < 0.2 {
        if l1 * l2 > 0.0 || a.min(b) < 1.0 {
            if l1.abs() < 0.1 {
                result *= (a * boost_log1p(l1)).exp();
            } else {
                result *= ((x * cgh) / agh).powf(a);
            }
            if l2.abs() < 0.1 {
                result *= (b * boost_log1p(l2)).exp();
            } else {
                result *= ((y * cgh) / bgh).powf(b);
            }
        } else if l1.abs().max(l2.abs()) < 0.5 {
            let small_a = a < b;
            let ratio = b / a;
            if (small_a && ratio * l2 < 0.1) || (!small_a && l1 / ratio > 0.1) {
                let mut l3 = boost_expm1(ratio * boost_log1p(l2));
                l3 = l1 + l3 + l3 * l1;
                l3 = a * boost_log1p(l3);
                result *= l3.exp();
            } else {
                let mut l3 = boost_expm1(boost_log1p(l1) / ratio);
                l3 = l2 + l3 + l3 * l2;
                l3 = b * boost_log1p(l3);
                result *= l3.exp();
            }
        } else if l1.abs() < l2.abs() {
            let mut l = a * boost_log1p(l1) + b * ((y * cgh) / bgh).ln();
            if l <= LOG_MIN_VALUE || l >= LOG_MAX_VALUE {
                l += result.ln();
                if l >= LOG_MAX_VALUE {
                    return f64::INFINITY;
                }
                result = l.exp();
            } else {
                result *= l.exp();
            }
        } else {
            let mut l = b * boost_log1p(l2) + a * ((x * cgh) / agh).ln();
            if l <= LOG_MIN_VALUE || l >= LOG_MAX_VALUE {
                l += result.ln();
                if l >= LOG_MAX_VALUE {
                    return f64::INFINITY;
                }
                result = l.exp();
            } else {
                result *= l.exp();
            }
        }
    } else {
        let b1 = (x * cgh) / agh;
        let b2 = (y * cgh) / bgh;
        l1 = a * b1.ln();
        l2 = b * b2.ln();
        if l1 >= LOG_MAX_VALUE || l1 <= LOG_MIN_VALUE || l2 >= LOG_MAX_VALUE || l2 <= LOG_MIN_VALUE
        {
            if a < b {
                let p1 = b2.powf(b / a);
                let l3 = if b1 != 0.0 && p1 != 0.0 {
                    a * (b1.ln() + p1.ln())
                } else {
                    MAX_VALUE
                };
                if l3 < LOG_MAX_VALUE && l3 > LOG_MIN_VALUE {
                    result *= (p1 * b1).powf(a);
                } else {
                    l2 += l1 + result.ln();
                    if l2 >= LOG_MAX_VALUE {
                        return f64::INFINITY;
                    }
                    result = l2.exp();
                }
            } else {
                let p1 = if b1 < 1.0 && b < 1.0 && MAX_VALUE * b < a {
                    0.0
                } else {
                    b1.powf(a / b)
                };
                let l3 = if p1 != 0.0 && b2 != 0.0 {
                    (p1.ln() + b2.ln()) * b
                } else {
                    MAX_VALUE
                };
                if l3 < LOG_MAX_VALUE && l3 > LOG_MIN_VALUE {
                    result *= (p1 * b2).powf(b);
                } else if result != 0.0 {
                    l2 += l1 + result.ln();
                    if l2 >= LOG_MAX_VALUE {
                        return f64::INFINITY;
                    }
                    result = l2.exp();
                }
            }
        } else {
            result *= b1.powf(a) * b2.powf(b);
        }
    }
    // An underflow to zero is returned as zero (the default underflow policy ignores it).
    result
}

/// `ibeta_series` (Lanczos, normalised, no derivative).
fn ibeta_series(a: f64, b: f64, x: f64, s0: f64) -> f64 {
    let c = a + b;
    let agh = a + LANCZOS_G - 0.5;
    let bgh = b + LANCZOS_G - 0.5;
    let cgh = c + LANCZOS_G - 0.5;
    let mut result = if a < MIN_VALUE || b < MIN_VALUE {
        0.0
    } else {
        lanczos_sum_expg_scaled(c) / (lanczos_sum_expg_scaled(a) * lanczos_sum_expg_scaled(b))
    };
    if !result.is_finite() {
        result = 0.0;
    }
    let l1 = (cgh / bgh).ln() * (b - 0.5);
    let l2 = (x * cgh / agh).ln() * a;
    if l1 > LOG_MIN_VALUE && l1 < LOG_MAX_VALUE && l2 > LOG_MIN_VALUE && l2 < LOG_MAX_VALUE {
        if a * b < bgh * 10.0 {
            result *= ((b - 0.5) * boost_log1p(a / bgh)).exp();
        } else {
            result *= (cgh / bgh).powf(b - 0.5);
        }
        result *= (x * cgh / agh).powf(a);
        result *= (agh / std::f64::consts::E).sqrt();
    } else if result != 0.0 {
        result = result.ln() + l1 + l2 + (agh.ln() - 1.0) / 2.0;
        result = result.exp();
    }
    if result < MIN_VALUE {
        return s0;
    }
    // ibeta_series_t summed by tools::sum_series from s0.
    let mut term_result = result;
    let mut apn = a;
    let mut poch = 1.0 - b;
    let mut n = 1i32;
    let mut sum = s0;
    let mut counter = MAX_SERIES_ITERATIONS;
    loop {
        let r = term_result / apn;
        apn += 1.0;
        term_result *= poch * x / f64::from(n);
        n += 1;
        poch += 1.0;
        sum += r;
        counter -= 1;
        if !((EPSILON * sum).abs() < r.abs()) || counter == 0 {
            break;
        }
    }
    sum
}

/// `ibeta_fraction2`: the continued fraction (modified Lentz, `continued_fraction_b`).
fn ibeta_fraction2(a: f64, b: f64, x: f64, y: f64) -> f64 {
    let result = ibeta_power_terms(a, b, x, y);
    if result == 0.0 {
        return result;
    }
    let mut m = 0i32;
    let mut next = || {
        let mf = f64::from(m);
        let mut a_n = (a + mf - 1.0) * (a + b + mf - 1.0) * mf * (b - mf) * x * x;
        let denom = a + f64::from(2 * m) - 1.0;
        a_n /= denom * denom;
        let mut b_n = mf;
        b_n += (mf * (b - mf) * x) / (a + f64::from(2 * m) - 1.0);
        b_n += ((a + mf) * (a * y - b * x + 1.0 + mf * (2.0 - x))) / (a + f64::from(2 * m) + 1.0);
        m += 1;
        (a_n, b_n)
    };
    let tiny = 16.0 * MIN_VALUE;
    let (_, b0) = next();
    let mut f = b0;
    if f == 0.0 {
        f = tiny;
    }
    let mut c = f;
    let mut d = 0.0;
    loop {
        let (a_n, b_n) = next();
        d = b_n + a_n * d;
        if d == 0.0 {
            d = tiny;
        }
        c = b_n + a_n / c;
        if c == 0.0 {
            c = tiny;
        }
        d = 1.0 / d;
        let delta = c * d;
        f *= delta;
        if !((delta - 1.0).abs() > EPSILON) {
            break;
        }
    }
    result / f
}

/// `binomial_ccdf(n, k, x, y)`: the finite sum for integer arguments.
fn binomial_ccdf(n: f64, k: f64, x: f64, y: f64) -> f64 {
    let mut result = x.powf(n);
    if result > MIN_VALUE {
        let mut term = result;
        let mut i = itrunc(n - 1.0) as u32;
        while f64::from(i) > k {
            term *= (f64::from(i + 1) * y) / ((n - f64::from(i)) * x);
            result += term;
            i = i.wrapping_sub(1);
        }
    } else {
        let mut start = itrunc(n * x);
        if f64::from(start) <= k + 1.0 {
            start = itrunc(k + 2.0);
        }
        result = x.powf(f64::from(start))
            * y.powf(n - f64::from(start))
            * boost_binomial_coefficient(itrunc(n) as u32, start as u32);
        if result == 0.0 {
            let mut i = (start - 1) as u32;
            while f64::from(i) > k {
                result += x.powf(f64::from(i))
                    * y.powf(n - f64::from(i))
                    * boost_binomial_coefficient(itrunc(n) as u32, i);
                i = i.wrapping_sub(1);
            }
        } else {
            let mut term = result;
            let start_term = result;
            let mut i = (start - 1) as u32;
            while f64::from(i) > k {
                term *= (f64::from(i + 1) * y) / ((n - f64::from(i)) * x);
                result += term;
                i = i.wrapping_sub(1);
            }
            term = start_term;
            let mut i = (start + 1) as u32;
            while f64::from(i) <= n {
                term *= (n - f64::from(i) + 1.0) * x / (f64::from(i) * y);
                result += term;
                i += 1;
            }
        }
    }
    result
}

/// The regularised incomplete beta by its continued fraction on the convergent side. Used only
/// where Boost's `ibeta_imp` would call `ibeta_a_step` and `beta_small_b_large_a_series` (BGRAT),
/// which happens only for non-integer a or b: never for the binomial tails, only for a `beta_ppf`
/// with non-integer allele counts. Accurate to a few units in the last place, not bit-identical.
fn ibeta_fraction_either_side(a: f64, b: f64, x: f64, y: f64) -> f64 {
    if x < (a + 1.0) / (a + b + 2.0) {
        ibeta_fraction2(a, b, x, y)
    } else {
        1.0 - ibeta_fraction2(b, a, y, x)
    }
}

/// `ibeta_imp(a, b, x, inv, normalised = true)` without the derivative.
fn ibeta_imp(a: f64, b: f64, x: f64, inv: bool) -> f64 {
    let (mut a, mut b, mut x) = (a, b, x);
    let mut invert = inv;
    let mut y = 1.0 - x;
    if !a.is_finite() || !b.is_finite() || !(0.0 <= x && x <= 1.0) || a < 0.0 || b < 0.0 {
        return f64::NAN;
    }
    if a == 0.0 {
        if b == 0.0 {
            return f64::NAN;
        }
        return if inv { 0.0 } else { 1.0 };
    } else if b == 0.0 {
        return if inv { 1.0 } else { 0.0 };
    }
    if x == 0.0 {
        return if invert { 1.0 } else { 0.0 };
    }
    if x == 1.0 {
        return if invert { 0.0 } else { 1.0 };
    }
    if a == 0.5 && b == 0.5 {
        let half_pi = std::f64::consts::FRAC_PI_2;
        return if invert {
            y.sqrt().asin() / half_pi
        } else {
            x.sqrt().asin() / half_pi
        };
    }
    if a == 1.0 {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut x, &mut y);
        invert = !invert;
    }
    if b == 1.0 {
        if a == 1.0 {
            return if invert { y } else { x };
        }
        return if y < 0.5 {
            if invert {
                -boost_expm1(a * boost_log1p(-y))
            } else {
                (a * boost_log1p(-y)).exp()
            }
        } else if invert {
            -boost_powm1(x, a)
        } else {
            x.powf(a)
        };
    }
    let mut fract;
    if a.min(b) <= 1.0 {
        // One or both of a, b below 1 (a, b != 1): only non-integer arguments get here.
        if x > 0.5 {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut x, &mut y);
            invert = !invert;
        }
        let series_first = if a.max(b) <= 1.0 {
            a >= b.min(0.2) || x.powf(a) <= 0.9
        } else {
            b <= 1.0 || (x < 0.1 && (b * x).powf(a) <= 0.7)
        };
        if series_first {
            fract = series_with_invert(a, b, x, &mut invert);
        } else {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut x, &mut y);
            invert = !invert;
            if y >= 0.3 {
                fract = series_with_invert(a, b, x, &mut invert);
            } else {
                fract = ibeta_fraction_either_side(a, b, x, y);
            }
        }
    } else {
        // Both a, b > 1.
        let lambda = if a < b {
            a - (a + b) * x
        } else {
            (a + b) * y - b
        };
        if lambda < 0.0 {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut x, &mut y);
            invert = !invert;
        }
        if b < 40.0 {
            if a.floor() == a && b.floor() == b && a < f64::from(i32::MAX - 100) && y != 1.0 {
                let k = a - 1.0;
                let n = b + k;
                fract = binomial_ccdf(n, k, x, y);
            } else if b * x <= 0.7 {
                fract = series_with_invert(a, b, x, &mut invert);
            } else {
                fract = ibeta_fraction_either_side(a, b, x, y);
            }
        } else {
            fract = ibeta_fraction2(a, b, x, y);
        }
    }
    if invert {
        fract = 1.0 - fract;
    }
    fract
}

/// The `ibeta_series` call of `ibeta_imp` with its handling of `invert` (normalised).
fn series_with_invert(a: f64, b: f64, x: f64, invert: &mut bool) -> f64 {
    if !*invert {
        ibeta_series(a, b, x, 0.0)
    } else {
        *invert = false;
        -ibeta_series(a, b, x, -1.0)
    }
}

/// Boost's `cdf(binomial_distribution(n, p), k)` (`complement` false) or
/// `cdf(complement(...))` (true), with SciPy's StatsPolicy (domain errors give NaN).
fn boost_binomial_cdf(k: f64, n: f64, p: f64, complement: bool) -> f64 {
    if !(0.0..=1.0).contains(&p) || !p.is_finite() || n < 0.0 || !n.is_finite() {
        return f64::NAN;
    }
    if k < 0.0 || !k.is_finite() || k > n {
        return f64::NAN;
    }
    if k == n {
        return if complement { 0.0 } else { 1.0 };
    }
    if p == 0.0 {
        return if complement { 0.0 } else { 1.0 };
    }
    if p == 1.0 {
        return if complement { 1.0 } else { 0.0 };
    }
    // cdf = ibetac(k + 1, n - k, p); complement = ibeta(k + 1, n - k, p).
    ibeta_imp(k + 1.0, n - k, p, !complement)
}

// ------------------------------------------------------------------------------------------------
// Beta quantile (beta.ppf): Newton steps on the ported incomplete beta

fn ibeta_quantile(a: f64, b: f64, q: f64) -> f64 {
    let cdf = |x: f64| ibeta_imp(a, b, x, false);
    // Bracket [low, high] with cdf(low) < q <= cdf(high); start from the mean.
    let mut low = 0.0_f64;
    let mut high = 1.0_f64;
    let mut x = a / (a + b);
    for _ in 0..400 {
        let value = cdf(x);
        if value.is_nan() {
            return f64::NAN;
        }
        if value < q {
            low = x;
        } else {
            high = x;
        }
        // The density x^(a-1) y^(b-1) / B(a, b) from the power terms.
        let y = 1.0 - x;
        let density = ibeta_power_terms(a, b, x, y) / (x * y);
        let mut next = if density > 0.0 && density.is_finite() {
            x - (value - q) / density
        } else {
            f64::NAN
        };
        if !(next > low && next < high) {
            next = low + (high - low) / 2.0;
        }
        if next == x || high - low <= f64::EPSILON * x.abs() * 0.5 {
            x = next;
            break;
        }
        x = next;
    }
    // Polish: move to a neighbouring double while that brings cdf(x) nearer to q.
    let distance = |x: f64| (cdf(x) - q).abs();
    let mut best = distance(x);
    for _ in 0..64 {
        let up = next_up(x);
        let down = next_down(x);
        let (candidate, value) = {
            let up_distance = if up < 1.0 {
                distance(up)
            } else {
                f64::INFINITY
            };
            let down_distance = if down > 0.0 {
                distance(down)
            } else {
                f64::INFINITY
            };
            if up_distance < down_distance {
                (up, up_distance)
            } else {
                (down, down_distance)
            }
        };
        if value < best {
            best = value;
            x = candidate;
        } else {
            break;
        }
    }
    x
}

fn next_up(x: f64) -> f64 {
    if x.is_nan() || x == f64::INFINITY {
        return x;
    }
    if x == 0.0 {
        return f64::from_bits(1);
    }
    let bits = x.to_bits();
    f64::from_bits(if x > 0.0 { bits + 1 } else { bits - 1 })
}

fn next_down(x: f64) -> f64 {
    -next_up(-x)
}

// ------------------------------------------------------------------------------------------------
// Double-double arithmetic with a separate binary exponent, for the beta quantile's refinement

/// A value `(hi + lo) * 2^exp` with `|lo| <= ulp(hi) / 2`: about 106 bits of precision and an
/// exponent range wide enough for binomial terms of any size.
#[derive(Clone, Copy, Debug)]
struct Wide {
    hi: f64,
    lo: f64,
    exp: i64,
}

fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    (s, (a - (s - bb)) + (b - bb))
}

fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    (p, a.mul_add(b, -p))
}

/// `a - b` of two double-double pairs.
fn pair_sub(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    let (s, e) = two_sum(a.0, -b.0);
    two_sum(s, e + (a.1 - b.1))
}

/// `a * b` of a double-double pair and a double.
fn pair_mul(a: (f64, f64), b: f64) -> (f64, f64) {
    let (p, e) = two_prod(a.0, b);
    two_sum(p, e + a.1 * b)
}

impl Wide {
    fn from_f64(value: f64) -> Wide {
        Wide {
            hi: value,
            lo: 0.0,
            exp: 0,
        }
        .normalized()
    }

    fn from_pair(hi: f64, lo: f64) -> Wide {
        let (hi, lo) = two_sum(hi, lo);
        Wide { hi, lo, exp: 0 }.normalized()
    }

    /// Moves the binary exponent of `hi` into `exp` (scaling by a power of two is exact).
    fn normalized(self) -> Wide {
        if self.hi == 0.0 || !self.hi.is_finite() {
            return Wide {
                hi: self.hi,
                lo: 0.0,
                exp: if self.hi == 0.0 { 0 } else { self.exp },
            };
        }
        let shift = ((self.hi.to_bits() >> 52) & 0x7ff) as i64 - 1023;
        if shift.abs() < 64 {
            return self;
        }
        let scale = 2f64.powi(-shift as i32);
        Wide {
            hi: self.hi * scale,
            lo: self.lo * scale,
            exp: self.exp + shift,
        }
    }

    fn mul(self, other: Wide) -> Wide {
        let (p, e) = two_prod(self.hi, other.hi);
        let e = e + (self.hi * other.lo + self.lo * other.hi);
        let (hi, lo) = two_sum(p, e);
        Wide {
            hi,
            lo,
            exp: self.exp + other.exp,
        }
        .normalized()
    }

    fn div(self, other: Wide) -> Wide {
        let x = (self.hi, self.lo);
        let y = (other.hi, other.lo);
        let q1 = x.0 / y.0;
        let r = pair_sub(x, pair_mul(y, q1));
        let q2 = r.0 / y.0;
        let r = pair_sub(r, pair_mul(y, q2));
        let q3 = r.0 / y.0;
        let (s, e) = two_sum(q1, q2);
        let (hi, lo) = two_sum(s, e + q3);
        Wide {
            hi,
            lo,
            exp: self.exp - other.exp,
        }
        .normalized()
    }

    /// `(hi, lo)` of the value expressed with binary exponent `exp`.
    fn scaled_to(self, exp: i64) -> (f64, f64) {
        let shift = self.exp - exp;
        if shift < -1000 {
            return (0.0, 0.0);
        }
        let factor = 2f64.powi(shift.min(1000) as i32);
        (self.hi * factor, self.lo * factor)
    }

    fn add(self, other: Wide) -> Wide {
        if self.hi == 0.0 {
            return other;
        }
        if other.hi == 0.0 {
            return self;
        }
        let exp = self.exp.max(other.exp);
        let (ah, al) = self.scaled_to(exp);
        let (bh, bl) = other.scaled_to(exp);
        let (s, e) = two_sum(ah, bh);
        let (t, f) = two_sum(al, bl);
        let (s, e) = two_sum(s, e + t);
        let (hi, lo) = two_sum(s, e + f);
        Wide { hi, lo, exp }.normalized()
    }

    fn sub(self, other: Wide) -> Wide {
        self.add(Wide {
            hi: -other.hi,
            lo: -other.lo,
            exp: other.exp,
        })
    }

    fn powi(self, mut n: u64) -> Wide {
        let mut result = Wide::from_f64(1.0);
        let mut base = self;
        while n > 0 {
            if n & 1 == 1 {
                result = result.mul(base);
            }
            base = base.mul(base);
            n >>= 1;
        }
        result
    }

    /// `self / other` as a double (zero when it underflows).
    fn ratio_to(self, other: Wide) -> f64 {
        let (hi, lo) = self.div(other).scaled_to(0);
        hi + lo
    }

    fn cmp_value(self, other: Wide) -> std::cmp::Ordering {
        self.sub(other)
            .hi
            .partial_cmp(&0.0)
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}

/// `I_x(a, b)` for integer a, b >= 1 as the binomial tail P(Bin(a + b - 1, x) >= a), summed in
/// double-double arithmetic from the first term (exact to about 1e-28 relative).
fn ibeta_integer_wide(a: u64, b: u64, x: f64, binomial: Wide) -> Wide {
    let n = a + b - 1;
    let (y_hi, y_lo) = two_sum(1.0, -x);
    let xw = Wide::from_f64(x);
    let yw = Wide::from_pair(y_hi, y_lo);
    let mut term = binomial.mul(xw.powi(a)).mul(yw.powi(n - a));
    let ratio = xw.div(yw);
    let mut sum = term;
    let mut j = a;
    while j < n {
        let factor = Wide::from_f64((n - j) as f64).div(Wide::from_f64((j + 1) as f64));
        term = term.mul(factor).mul(ratio);
        sum = sum.add(term);
        j += 1;
        if term.hi == 0.0 || term.ratio_to(sum).abs() < 1e-34 {
            break;
        }
    }
    sum
}

/// C(n, k) in double-double arithmetic.
fn binomial_wide(n: u64, k: u64) -> Wide {
    let k = k.min(n - k);
    let mut numerator = Wide::from_f64(1.0);
    let mut denominator = Wide::from_f64(1.0);
    for i in 1..=k {
        numerator = numerator.mul(Wide::from_f64((n - k + i) as f64));
        denominator = denominator.mul(Wide::from_f64(i as f64));
    }
    numerator.div(denominator)
}

/// The double nearest to the q-quantile of Beta(a, b) for integer a, b, starting from an
/// estimate within a few hundred units in the last place. `None` when the bracket is not found.
fn correctly_rounded_quantile(a: u64, b: u64, q: f64, estimate: f64) -> Option<f64> {
    let binomial = binomial_wide(a + b - 1, a);
    let target = Wide::from_f64(q);
    let cdf = |x: f64| ibeta_integer_wide(a, b, x, binomial);
    let below = |x: f64| cdf(x).cmp_value(target) == std::cmp::Ordering::Less;
    let step = |x: f64, ulps: i64| f64::from_bits((x.to_bits() as i64 + ulps) as u64);
    // Bracket: low with cdf < q, high with cdf >= q, widening by powers of two.
    let (mut low, mut high);
    if below(estimate) {
        low = estimate;
        let mut distance = 1i64;
        loop {
            let candidate = step(estimate, distance);
            if !(candidate < 1.0) {
                return None;
            }
            if !below(candidate) {
                high = candidate;
                break;
            }
            low = candidate;
            distance *= 2;
            if distance > 1 << 20 {
                return None;
            }
        }
    } else {
        high = estimate;
        let mut distance = 1i64;
        loop {
            let candidate = step(estimate, -distance);
            if !(candidate > 0.0) {
                return None;
            }
            if below(candidate) {
                low = candidate;
                break;
            }
            high = candidate;
            distance *= 2;
            if distance > 1 << 20 {
                return None;
            }
        }
    }
    while high.to_bits() - low.to_bits() > 1 {
        let middle = f64::from_bits(low.to_bits() + (high.to_bits() - low.to_bits()) / 2);
        if below(middle) {
            low = middle;
        } else {
            high = middle;
        }
    }
    // The root lies between the adjacent doubles low and high: take the nearer one.
    let at_low = cdf(low);
    let at_high = cdf(high);
    let fraction = target.sub(at_low).ratio_to(at_high.sub(at_low));
    Some(if fraction <= 0.5 { low } else { high })
}

// ------------------------------------------------------------------------------------------------
// Boost.Math: the hypergeometric distribution (detail/hypergeometric_pdf.hpp and _cdf.hpp)

/// `hypergeometric_pdf<double>(x, r, n, N)`, clipped to [0, 1].
fn boost_hypergeometric_pdf(x: u64, r: u64, n: u64, big_n: u64) -> f64 {
    let result = if big_n <= MAX_FACTORIAL {
        hypergeometric_pdf_factorial(x, r, n, big_n)
    } else if big_n <= PRIME_METHOD_LIMIT {
        hypergeometric_pdf_prime(x, r, n, big_n)
    } else {
        hypergeometric_pdf_lanczos(x, r, n, big_n)
    };
    result.clamp(0.0, 1.0)
}

fn hypergeometric_pdf_factorial(x: u64, r: u64, n: u64, big_n: u64) -> f64 {
    let mut result = unchecked_factorial(n);
    let num = [
        unchecked_factorial(r),
        unchecked_factorial(big_n - n),
        unchecked_factorial(big_n - r),
    ];
    let denom = [
        unchecked_factorial(big_n),
        unchecked_factorial(x),
        unchecked_factorial(n - x),
        unchecked_factorial(r - x),
        unchecked_factorial(big_n - n - r + x),
    ];
    let mut i = 0;
    let mut j = 0;
    while i < 3 || j < 5 {
        while j < 5 && (result >= 1.0 || i >= 3) {
            result /= denom[j];
            j += 1;
        }
        while i < 3 && (result <= 1.0 || j >= 5) {
            result *= num[i];
            i += 1;
        }
    }
    result
}

/// Boost's prime table: `prime(0)` = 2 to `prime(9999)` = 104729; the pdf uses prime factorisation
/// up to `prime(max_prime - 1)` with `max_prime` = 9999, that is N <= 104723.
const LARGEST_TABLE_PRIME: u64 = 104_729;
const PRIME_METHOD_LIMIT: u64 = 104_723;

fn primes() -> &'static [u64] {
    static PRIMES: OnceLock<Vec<u64>> = OnceLock::new();
    PRIMES.get_or_init(|| {
        let limit = LARGEST_TABLE_PRIME as usize;
        let mut composite = vec![false; limit + 1];
        let mut list = Vec::with_capacity(10_000);
        for value in 2..=limit {
            if !composite[value] {
                list.push(value as u64);
                let mut multiple = value * value;
                while multiple <= limit {
                    composite[multiple] = true;
                    multiple += value;
                }
            }
        }
        list
    })
}

/// `integer_power(x, ex)` with Boost's `pow<N>` for 4 to 8 and the C library's pow above.
fn integer_power(x: f64, ex: i32) -> f64 {
    if ex < 0 {
        return 1.0 / integer_power(x, -ex);
    }
    match ex {
        0 => 1.0,
        1 => x,
        2 => x * x,
        3 => x * x * x,
        4 => {
            let p2 = x * x;
            p2 * p2
        }
        5 => {
            let p2 = x * x;
            x * p2 * p2
        }
        6 => {
            let p3 = x * x * x;
            p3 * p3
        }
        7 => {
            let p3 = x * x * x;
            x * p3 * p3
        }
        8 => {
            let p2 = x * x;
            let p4 = p2 * p2;
            p4 * p4
        }
        _ => x.powf(f64::from(ex)),
    }
}

fn hypergeometric_pdf_prime(x: u64, r: u64, n: u64, big_n: u64) -> f64 {
    // The recursion of `hypergeometric_pdf_prime_loop_imp` as a list of partial products, the
    // last one current; a factor that would overflow or underflow starts a new partial product.
    let mut partials: Vec<f64> = vec![1.0];
    for &prime in primes() {
        if prime > big_n {
            break;
        }
        let mut base = prime;
        let mut prime_powers: i64 = 0;
        while base <= big_n {
            prime_powers += (n / base) as i64;
            prime_powers += (r / base) as i64;
            prime_powers += ((big_n - n) / base) as i64;
            prime_powers += ((big_n - r) / base) as i64;
            prime_powers -= (big_n / base) as i64;
            prime_powers -= (x / base) as i64;
            prime_powers -= ((n - x) / base) as i64;
            prime_powers -= ((r - x) / base) as i64;
            prime_powers -= ((big_n - n - r + x) / base) as i64;
            base *= prime;
        }
        if prime_powers != 0 {
            let p = integer_power(prime as f64, prime_powers as i32);
            let current = *partials.last().expect("one partial product");
            if (p > 1.0 && MAX_VALUE / p < current) || (p < 1.0 && MIN_VALUE / p > current) {
                partials.push(p);
                continue;
            }
            *partials.last_mut().expect("one partial product") *= p;
        }
    }
    // Multiply the partial products, newest first, alternating between those >= 1 and those < 1.
    let newer_at_least_one = |from: Option<usize>| -> Option<usize> {
        let mut index = from;
        while let Some(at) = index {
            if partials[at] >= 1.0 {
                return Some(at);
            }
            index = at.checked_sub(1);
        }
        None
    };
    let newer_below_one = |from: Option<usize>| -> Option<usize> {
        let mut index = from;
        while let Some(at) = index {
            if partials[at] < 1.0 {
                return Some(at);
            }
            index = at.checked_sub(1);
        }
        None
    };
    let newest = Some(partials.len() - 1);
    let mut i = newer_at_least_one(newest);
    let mut j = newer_below_one(newest);
    let mut product = 1.0;
    while i.is_some() || j.is_some() {
        while let Some(at) = i {
            if !(product <= 1.0 || j.is_none()) {
                break;
            }
            product *= partials[at];
            i = newer_at_least_one(at.checked_sub(1));
        }
        while let Some(at) = j {
            if !(product >= 1.0 || i.is_none()) {
                break;
            }
            product *= partials[at];
            j = newer_below_one(at.checked_sub(1));
        }
    }
    product
}

fn hypergeometric_pdf_lanczos(x: u64, r: u64, n: u64, big_n: u64) -> f64 {
    let g = LANCZOS_G;
    let mut bases = [
        n as f64 + g + 0.5,
        r as f64 + g + 0.5,
        (big_n - n) as f64 + g + 0.5,
        (big_n - r) as f64 + g + 0.5,
        1.0 / (big_n as f64 + g + 0.5),
        1.0 / (x as f64 + g + 0.5),
        1.0 / ((n - x) as f64 + g + 0.5),
        1.0 / ((r - x) as f64 + g + 0.5),
        1.0 / ((big_n - n - r + x) as f64 + g + 0.5),
    ];
    let mut exponents = [
        n as f64 + 0.5,
        r as f64 + 0.5,
        (big_n - n) as f64 + 0.5,
        (big_n - r) as f64 + 0.5,
        big_n as f64 + 0.5,
        x as f64 + 0.5,
        (n - x) as f64 + 0.5,
        (r - x) as f64 + 0.5,
        (big_n - n - r + x) as f64 + 0.5,
    ];
    let mut base_e_factors = [-1i32, -1, -1, -1, 1, 1, 1, 1, 1];
    // std::sort with "larger exponent first": libstdc++ uses insertion sort for 9 elements,
    // which keeps equal exponents in their original order, as a stable sort does.
    let mut sorted: Vec<usize> = (0..9).collect();
    sorted.sort_by(|&i, &j| {
        exponents[j]
            .partial_cmp(&exponents[i])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let bubble_down_one = |sorted: &mut [usize], exponents: &[f64; 9]| {
        let mut first = 0;
        while first + 1 < sorted.len() && !(exponents[sorted[first]] > exponents[sorted[first + 1]])
        {
            sorted.swap(first, first + 1);
            first += 1;
        }
    };
    loop {
        exponents[sorted[0]] -= exponents[sorted[1]];
        bases[sorted[1]] *= bases[sorted[0]];
        if bases[sorted[1]] < MIN_VALUE && exponents[sorted[1]] != 0.0 {
            return 0.0;
        }
        base_e_factors[sorted[1]] += base_e_factors[sorted[0]];
        bubble_down_one(&mut sorted[..], &exponents);
        if !(exponents[sorted[1]] > 1.0) {
            break;
        }
    }
    let mut j = 8usize;
    while exponents[sorted[j]] == 0.0 {
        j -= 1;
    }
    while j > 0 {
        while j > 0 && exponents[sorted[j - 1]] == exponents[sorted[j]] {
            bases[sorted[j - 1]] *= bases[sorted[j]];
            exponents[sorted[j]] = 0.0;
            base_e_factors[sorted[j - 1]] += base_e_factors[sorted[j]];
            bubble_down_one(&mut sorted[j..], &exponents);
            j -= 1;
        }
        if j == 0 {
            // `--j` on an unsigned zero would wrap and end the outer loop.
            break;
        }
        j -= 1;
    }
    let factor = |index: usize| bases[index] * f64::from(base_e_factors[index]).exp();
    let mut result = factor(sorted[0]).powf(exponents[sorted[0]]);
    for &index in sorted.iter().skip(1) {
        if !(exponents[index] > 0.0) {
            break;
        }
        if result < MIN_VALUE {
            return 0.0;
        }
        if exponents[index] == 1.0 {
            result *= factor(index);
        } else if exponents[index] == 0.5 {
            result *= factor(index).sqrt();
        } else {
            result *= factor(index).powf(exponents[index]);
        }
    }
    let l = lanczos_sum_expg_scaled;
    result *= l((n + 1) as f64)
        * l((r + 1) as f64)
        * l((big_n - n + 1) as f64)
        * l((big_n - r + 1) as f64)
        / (l((big_n + 1) as f64)
            * l((x + 1) as f64)
            * l((n - x + 1) as f64)
            * l((r - x + 1) as f64)
            * l((big_n - n - r + x + 1) as f64));
    result
}

/// `hypergeometric_cdf<double>(x, r, n, N, invert = false)`, clipped to [0, 1].
fn boost_hypergeometric_cdf(x: u64, r: u64, n: u64, big_n: u64) -> f64 {
    let mut x = x;
    let mut invert = false;
    let mut result = 0.0;
    let mode = ((r + 1) as f64 * (n + 1) as f64 / (big_n + 2) as f64).floor();
    if (x as f64) < mode {
        result = boost_hypergeometric_pdf(x, r, n, big_n);
        let mut diff = result;
        let lower_limit = ((n + r) as i64 - big_n as i64).max(0) as u64;
        while diff > (if invert { 1.0 } else { result }) * EPSILON {
            diff = x as f64 * ((big_n + x) - n - r) as f64 * diff
                / ((1 + n - x) as f64 * (1 + r - x) as f64);
            result += diff;
            if x == lower_limit {
                break;
            }
            x -= 1;
        }
    } else {
        invert = !invert;
        let upper_limit = r.min(n);
        if x != upper_limit {
            x += 1;
            result = boost_hypergeometric_pdf(x, r, n, big_n);
            let mut diff = result;
            while x <= upper_limit && diff > (if invert { 1.0 } else { result }) * EPSILON {
                diff = (n - x) as f64 * (r - x) as f64 * diff
                    / ((x + 1) as f64 * ((big_n + x + 1) - n - r) as f64);
                result += diff;
                x += 1;
            }
        }
    }
    if invert {
        result = 1.0 - result;
    }
    result.clamp(0.0, 1.0)
}

// ------------------------------------------------------------------------------------------------
// Cephes (SciPy's special/cephes headers): chdtrc(1, x) = igamc(0.5, x / 2)

const MACHEP: f64 = 1.11022302462515654042E-16;
const MAXLOG: f64 = 7.09782712893383996732E2;
const IGAM_MAXITER: i32 = 2000;
const IGAM_BIG: f64 = 4.503599627370496e15;
const IGAM_BIGINV: f64 = 2.22044604925031308085e-16;
/// cephes `lgam(0.5)` (SciPy's `gammaln(0.5)`, 0x1.250d048e7a1bcp-1).
const LGAM_HALF: f64 = f64::from_bits(0x3fe2_50d0_48e7_a1bc);
/// cephes `lgam1p(0.5)` (SciPy's `_lgam1p(0.5)`, -0x1.eeb95b094c296p-4).
const LGAM1P_HALF: f64 = f64::from_bits(0xbfbe_eb95_b094_c296);
/// cephes `lanczos_sum_expg_scaled(0.5)` (SciPy's `_lanczos_sum_expg_scaled(0.5)`,
/// 0x1.c5bf891b4ef6bp+0).
const LANCZOS_EXPG_SCALED_HALF: f64 = f64::from_bits(0x3ffc_5bf8_91b4_ef6b);

/// `chdtrc(1, x)`.
fn cephes_chdtrc_1(x: f64) -> f64 {
    if x < 0.0 {
        return 1.0;
    }
    igamc_half(x / 2.0)
}

/// cephes `igamc(0.5, x)` (no asymptotic regime for a = 0.5).
fn igamc_half(x: f64) -> f64 {
    let a = 0.5;
    if x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return 1.0;
    }
    if x.is_infinite() {
        return 0.0;
    }
    if x > 1.1 {
        if x < a {
            1.0 - igam_series(a, x)
        } else {
            igamc_continued_fraction(a, x)
        }
    } else if x <= 0.5 {
        if -0.4 / x.ln() < a {
            1.0 - igam_series(a, x)
        } else {
            igamc_series(a, x)
        }
    } else if x * 1.1 < a {
        1.0 - igam_series(a, x)
    } else {
        igamc_series(a, x)
    }
}

/// cephes `igam_fac(0.5, x)`: x^a exp(-x) / gamma(a).
fn igam_fac(a: f64, x: f64) -> f64 {
    if (a - x).abs() > 0.4 * a.abs() {
        let ax = a * x.ln() - x - LGAM_HALF;
        if ax < -MAXLOG {
            return 0.0;
        }
        return ax.exp();
    }
    let fac = a + LANCZOS_G - 0.5;
    let mut res = (fac / std::f64::consts::E).sqrt() / LANCZOS_EXPG_SCALED_HALF;
    // a = 0.5 and x within 0.2 of it: always the a < 200, x < 200 branch.
    res *= (a - x).exp() * (x / fac).powf(a);
    res
}

fn igamc_continued_fraction(a: f64, x: f64) -> f64 {
    let ax = igam_fac(a, x);
    if ax == 0.0 {
        return 0.0;
    }
    let mut y = 1.0 - a;
    let mut z = x + y + 1.0;
    let mut c = 0.0;
    let mut pkm2 = 1.0;
    let mut qkm2 = x;
    let mut pkm1 = x + 1.0;
    let mut qkm1 = z * x;
    let mut ans = pkm1 / qkm1;
    for _ in 0..IGAM_MAXITER {
        c += 1.0;
        y += 1.0;
        z += 2.0;
        let yc = y * c;
        let pk = pkm1 * z - pkm2 * yc;
        let qk = qkm1 * z - qkm2 * yc;
        let t;
        if qk != 0.0 {
            let r = pk / qk;
            t = ((ans - r) / r).abs();
            ans = r;
        } else {
            t = 1.0;
        }
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if pk.abs() > IGAM_BIG {
            pkm2 *= IGAM_BIGINV;
            pkm1 *= IGAM_BIGINV;
            qkm2 *= IGAM_BIGINV;
            qkm1 *= IGAM_BIGINV;
        }
        if t <= MACHEP {
            break;
        }
    }
    ans * ax
}

fn igam_series(a: f64, x: f64) -> f64 {
    let ax = igam_fac(a, x);
    if ax == 0.0 {
        return 0.0;
    }
    let mut r = a;
    let mut c = 1.0;
    let mut ans = 1.0;
    for _ in 0..IGAM_MAXITER {
        r += 1.0;
        c *= x / r;
        ans += c;
        if c <= MACHEP * ans {
            break;
        }
    }
    ans * ax / a
}

fn igamc_series(a: f64, x: f64) -> f64 {
    let mut fac = 1.0;
    let mut sum = 0.0;
    for n in 1..IGAM_MAXITER {
        fac *= -x / f64::from(n);
        let term = fac / (a + f64::from(n));
        sum += term;
        if term.abs() <= MACHEP * sum.abs() {
            break;
        }
    }
    let logx = x.ln();
    let term = -cephes_expm1(a * logx - LGAM1P_HALF);
    term - (a * logx - LGAM_HALF).exp() * sum
}

/// cephes `expm1` (unity.h).
#[allow(clippy::excessive_precision)]
fn cephes_expm1(x: f64) -> f64 {
    const EP: [f64; 3] = [
        1.2617719307481059087798E-4,
        3.0299440770744196129956E-2,
        9.9999999999999999991025E-1,
    ];
    const EQ: [f64; 4] = [
        3.0019850513866445504159E-6,
        2.5244834034968410419224E-3,
        2.2726554820815502876593E-1,
        2.0000000000000000000897E0,
    ];
    if !x.is_finite() {
        if x.is_nan() || x > 0.0 {
            return x;
        }
        return -1.0;
    }
    if !(-0.5..=0.5).contains(&x) {
        return x.exp() - 1.0;
    }
    let polevl = |value: f64, coefficients: &[f64]| {
        coefficients[1..]
            .iter()
            .fold(coefficients[0], |answer, coefficient| {
                answer * value + coefficient
            })
    };
    let xx = x * x;
    let mut r = x * polevl(xx, &EP);
    r /= polevl(xx, &EQ) - r;
    r + r
}

// ------------------------------------------------------------------------------------------------
// Comparison harness entry

/// `priva acmg-step stats FUNCTION --input FRAME --output FRAME --options JSON`: bulk evaluation
/// for the comparison with SciPy. The values travel as raw little-endian binary files named in
/// the options (the frame files would be too slow for millions of inputs):
/// - `binom_sf`, `binom_cdf`, `beta_ppf`: `{"input": path, "output": path}`, input rows of three
///   float64 (k, n, p or q, a, b), output one float64 per row;
/// - `chi2_contingency_2x2`, `fisher_exact_greater`: input rows of four int64 (a, b, c, d) or
///   `{"max_total": T}` for every table with a total up to T, in the order total, a, b, c;
///   output two float64 per row (statistic or odds ratio, p-value).
///
/// The returned frame has one column, `rows`, with the number of rows evaluated.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    const FUNCTIONS: [&str; 5] = [
        "binom_sf",
        "binom_cdf",
        "beta_ppf",
        "chi2_contingency_2x2",
        "fisher_exact_greater",
    ];
    if !FUNCTIONS.contains(&function) {
        return Err(invalid(format!("stats has no harness function {function}")));
    }
    let options = &call.options;
    let path = |name: &str| -> Result<String> {
        match options.get(name) {
            Some(Json::Str(text)) => Ok(text.clone()),
            _ => Err(invalid(format!("option {name} (a path) is required"))),
        }
    };
    let output_path = path("output")?;
    let mut output: Vec<u8> = Vec::new();
    let rows;
    match function {
        "binom_sf" | "binom_cdf" | "beta_ppf" => {
            let bytes = std::fs::read(path("input")?)?;
            if bytes.len() % 24 != 0 {
                return Err(invalid("input is not rows of three float64"));
            }
            let f: fn(f64, f64, f64) -> f64 = match function {
                "binom_sf" => binom_sf,
                "binom_cdf" => binom_cdf,
                _ => beta_ppf,
            };
            rows = bytes.len() / 24;
            output.reserve(rows * 8);
            for row in bytes.chunks_exact(24) {
                let value = |at: usize| {
                    f64::from_le_bytes(row[at..at + 8].try_into().expect("eight bytes"))
                };
                output.extend_from_slice(&f(value(0), value(8), value(16)).to_le_bytes());
            }
        }
        "chi2_contingency_2x2" | "fisher_exact_greater" => {
            let f: fn([[i64; 2]; 2]) -> (f64, f64) = if function == "fisher_exact_greater" {
                fisher_exact_greater
            } else {
                chi2_contingency_2x2
            };
            let mut push = |table: [[i64; 2]; 2]| {
                let (first, second) = f(table);
                output.extend_from_slice(&first.to_le_bytes());
                output.extend_from_slice(&second.to_le_bytes());
            };
            if let Some(Json::Int(max_total)) = options.get("max_total") {
                let max_total = i64::try_from(*max_total).map_err(|_| invalid("max_total"))?;
                let mut count = 0;
                for total in 0..=max_total {
                    for a in 0..=total {
                        for b in 0..=total - a {
                            for c in 0..=total - a - b {
                                push([[a, b], [c, total - a - b - c]]);
                                count += 1;
                            }
                        }
                    }
                }
                rows = count;
            } else {
                let bytes = std::fs::read(path("input")?)?;
                if bytes.len() % 32 != 0 {
                    return Err(invalid("input is not rows of four int64"));
                }
                rows = bytes.len() / 32;
                for row in bytes.chunks_exact(32) {
                    let value = |at: usize| {
                        i64::from_le_bytes(row[at..at + 8].try_into().expect("eight bytes"))
                    };
                    push([[value(0), value(8)], [value(16), value(24)]]);
                }
            }
        }
        _ => {
            return Err(invalid(format!("stats has no harness function {function}")));
        }
    }
    std::fs::write(output_path, output)?;
    Ok(Table {
        names: vec!["rows".to_owned()],
        columns: vec![Column::Int(vec![rows as i64])],
        rows: 1,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hex(value: f64) -> String {
        format!("{:016x}", value.to_bits())
    }

    /// Parses Python's `float.hex()` text.
    fn from_python_hex(text: &str) -> f64 {
        let (negative, text) = match text.strip_prefix('-') {
            Some(rest) => (true, rest),
            None => (false, text),
        };
        if text == "inf" {
            return if negative {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            };
        }
        let text = text.strip_prefix("0x").expect("hex float");
        let (mantissa, exponent) = text.split_once('p').expect("exponent");
        let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
        let mut value = u64::from_str_radix(whole, 16).expect("whole") as f64;
        let mut scale = 1.0 / 16.0;
        for digit in fraction.chars() {
            value += f64::from(digit.to_digit(16).expect("digit")) * scale;
            scale /= 16.0;
        }
        let exponent: i32 = exponent.parse().expect("exponent value");
        let value = value * 2f64.powi(exponent);
        if negative { -value } else { value }
    }

    /// SciPy 1.14.1 (priva_acmg environment, HGSC) values, printed with `float.hex()`.
    #[test]
    fn binomial_tails_are_scipys_values() {
        let cases: [(&str, f64, f64, f64, &str); 15] = [
            ("sf", 0.0, 10.0, 0.1, "0x1.4d7a04fde6815p-1"),
            ("sf", 3.0, 100.0, 0.01, "0x1.2d0a4b6804a08p-6"),
            ("sf", 9.0, 1000.0, 0.001, "0x1.cd66b5211b7fap-24"),
            ("sf", 150.0, 20000.0, 0.005, "0x1.35bef204ba427p-20"),
            ("sf", 1.0, 250000.0, 1e-6, "0x1.b228a4eb28faep-6"),
            ("sf", 37.0, 1614000.0, 1e-5, "0x1.5093e8388e0cdp-19"),
            ("sf", 5.0, 6.0, 0.3, "0x1.7e34b945308bap-11"),
            ("sf", 12000.0, 1500000.0, 0.0079, "0x1.530c463eecca4p-4"),
            ("cdf", 0.0, 1.0, 0.01, "0x1.fae147ae147aep-1"),
            ("cdf", 2.0, 50.0, 0.0001, "0x1.ffffff583acbap-1"),
            ("cdf", 10.0, 100000.0, 0.01, "0x0.0p+0"),
            ("cdf", 990.0, 100000.0, 0.01, "0x1.8859cdd570ecap-2"),
            ("cdf", 30.0, 1200.0, 0.01, "0x1.ffff9e86c6b5ap-1"),
            ("cdf", 4.0, 800000.0, 3e-6, "0x1.ceea5a00fd93bp-1"),
            ("cdf", 64000.0, 129000.0, 0.5, "0x1.62a66ead3cf9ep-9"),
        ];
        for (kind, k, n, p, expected) in cases {
            let value = if kind == "sf" {
                binom_sf(k, n, p)
            } else {
                binom_cdf(k, n, p)
            };
            assert_eq!(
                hex(value),
                hex(from_python_hex(expected)),
                "{kind}({k}, {n}, {p})"
            );
        }
    }

    #[test]
    fn binomial_edges_follow_scipys_wrappers() {
        assert_eq!(binom_sf(-1.0, 10.0, 0.5), 1.0);
        assert_eq!(binom_sf(10.0, 10.0, 0.5), 0.0);
        assert_eq!(binom_sf(3.0, 10.0, 0.0), 0.0);
        assert_eq!(binom_sf(3.0, 10.0, 1.0), 1.0);
        assert_eq!(binom_cdf(10.0, 10.0, 0.5), 1.0);
        assert_eq!(binom_cdf(-1.0, 10.0, 0.5), 0.0);
        assert_eq!(binom_cdf(3.0, 10.0, 0.0), 1.0);
        assert_eq!(binom_cdf(3.0, 10.0, 1.0), 0.0);
        assert!(binom_sf(3.0, 10.5, 0.5).is_nan());
        assert!(binom_cdf(3.0, 10.0, 1.5).is_nan());
        assert!(binom_cdf(f64::NAN, 10.0, 0.5).is_nan());
        assert!(binom_sf(3.0, -1.0, 0.5).is_nan());
        assert_eq!(binom_cdf(2.7, 10.0, 0.3), binom_cdf(2.0, 10.0, 0.3));
    }

    #[test]
    fn beta_quantiles_are_near_scipys() {
        let cases: [(f64, f64, f64, &str); 5] = [
            (0.01, 1001.0, 9.0, "0x1.f735dd88b6a3fp-1"),
            (0.01, 61.0, 41.0, "0x1.ef1e632029e43p-2"),
            (0.01, 9001.0, 1001.0, "0x1.c9200ad3f47e6p-1"),
            (0.01, 3.0, 2.0, "0x1.207f298787156p-3"),
            (0.01, 1.0, 101.0, "0x1.a158e1b5421c8p-14"),
        ];
        for (q, a, b, expected) in cases {
            let value = beta_ppf(q, a, b);
            let expected = from_python_hex(expected);
            let ulps = (value.to_bits() as i64 - expected.to_bits() as i64).abs();
            assert!(
                ulps <= 1,
                "beta_ppf({q}, {a}, {b}) = {value:e}, SciPy {expected:e}"
            );
        }
        assert!(beta_ppf(0.01, 0.0, 3.0).is_nan());
        assert!(beta_ppf(0.01, 2.0, -1.0).is_nan());
        assert_eq!(beta_ppf(0.0, 2.0, 3.0), 0.0);
        assert_eq!(beta_ppf(1.0, 2.0, 3.0), 1.0);
    }

    #[test]
    fn two_by_two_tests_are_scipys_values() {
        let chi2: [([i64; 4], &str, &str); 6] = [
            (
                [10, 3, 2, 15],
                "0x1.4ea91f46a91f4p+3",
                "0x1.401a88c560cb6p-10",
            ),
            ([5, 5, 5, 5], "0x0.0p+0", "0x1.0000000000000p+0"),
            (
                [30, 12, 7, 40],
                "0x1.ae78202487a82p+4",
                "0x1.cb169e3ad58fep-23",
            ),
            ([0, 4, 9, 1], "0x1.a2a1907f6e5d4p+2", "0x1.5967e43b95fd1p-7"),
            (
                [2, 9, 40, 3],
                "0x1.838a7a0202a10p+4",
                "0x1.cd0c54f29e806p-21",
            ),
            (
                [100, 7, 220, 31],
                "0x1.0b9a199ea5b73p+1",
                "0x1.2f858d0fc4747p-3",
            ),
        ];
        for (cells, statistic, pvalue) in chi2 {
            let (s, p) = chi2_contingency_2x2([[cells[0], cells[1]], [cells[2], cells[3]]]);
            assert_eq!(
                hex(s),
                hex(from_python_hex(statistic)),
                "chi2 statistic {cells:?}"
            );
            assert_eq!(hex(p), hex(from_python_hex(pvalue)), "chi2 p {cells:?}");
        }
        let fisher: [([i64; 4], &str, &str); 7] = [
            (
                [10, 3, 2, 15],
                "0x1.9000000000000p+4",
                "0x1.e7c70f03f508ep-12",
            ),
            ([5, 5, 5, 5], "0x1.0000000000000p+0", "0x1.57fdea81a3466p-1"),
            (
                [30, 12, 7, 40],
                "0x1.c924924924925p+3",
                "0x1.ada6cc8b41e0cp-25",
            ),
            ([1, 0, 3, 2], "inf", "0x1.5555555555556p-1"),
            ([0, 4, 9, 1], "0x0.0p+0", "0x1.0000000000000p+0"),
            (
                [2, 9, 40, 3],
                "0x1.1111111111111p-6",
                "0x1.ffffff062bafdp-1",
            ),
            (
                [100, 7, 220, 31],
                "0x1.01a98ef606a64p+1",
                "0x1.1fd762d35c2cfp-4",
            ),
        ];
        for (cells, odds, pvalue) in fisher {
            let (o, p) = fisher_exact_greater([[cells[0], cells[1]], [cells[2], cells[3]]]);
            assert_eq!(hex(o), hex(from_python_hex(odds)), "fisher odds {cells:?}");
            assert_eq!(hex(p), hex(from_python_hex(pvalue)), "fisher p {cells:?}");
        }
        let (o, p) = fisher_exact_greater([[0, 0], [2, 3]]);
        assert!(o.is_nan());
        assert_eq!(p, 1.0);
    }

    #[test]
    fn double_double_tail_is_exact_on_simple_cases() {
        let value = |w: Wide| {
            let (hi, lo) = w.scaled_to(0);
            hi + lo
        };
        // I_x(1, 1) = x; I_x(2, 1) = x^2; I_x(1, 2) = 1 - (1 - x)^2.
        assert_eq!(
            value(ibeta_integer_wide(1, 1, 0.3, binomial_wide(1, 1))),
            0.3
        );
        assert_eq!(
            value(ibeta_integer_wide(2, 1, 0.75, binomial_wide(2, 2))),
            0.5625
        );
        assert_eq!(
            value(ibeta_integer_wide(1, 2, 0.5, binomial_wide(2, 1))),
            0.75
        );
        assert_eq!(value(binomial_wide(52, 5)), 2598960.0);
        // Large exponents stay representable.
        let tiny = Wide::from_f64(1e-300).mul(Wide::from_f64(1e-300));
        assert_eq!(tiny.ratio_to(Wide::from_f64(1e-300)), 1e-300);
    }

    #[test]
    fn boost_pieces_match_known_values() {
        // The table's primes and factorials.
        assert_eq!(primes().len(), 10_000);
        assert_eq!(primes()[9_999], LARGEST_TABLE_PRIME);
        assert_eq!(primes()[9_998], PRIME_METHOD_LIMIT);
        assert_eq!(unchecked_factorial(20), 2432902008176640000.0);
        assert_eq!(boost_binomial_coefficient(10, 3), 120.0);
        assert_eq!(boost_binomial_coefficient(200, 2), 19900.0);
        // The cephes constants agree with their definitions to a few units in the last place.
        assert!((LGAM_HALF - std::f64::consts::PI.sqrt().ln()).abs() < 1e-15);
        assert_eq!(cephes_expm1(0.0), 0.0);
        assert!((boost_log1p(0.25) - 1.25f64.ln()).abs() < 1e-16);
        assert!((boost_expm1(0.25) - 0.25f64.exp_m1()).abs() < 1e-16);
    }
}
