//! The four SciPy 1.14.1 functions the ACMG step calls. Thread T5 fills this
//! module first; thread T3 (`acmg_scoring`) and thread T6
//! (`acmg_pm1_pp2_bp1_missense`) call it by these signatures.
//!
//! The owner's decision (2026-10-30): every result that only feeds a
//! threshold test must give the same decision as SciPy on all inputs the
//! tests can enumerate; `beta_ppf` may differ from SciPy in the last printed
//! digit (it becomes the printed `control_common_index`).
//!
//! Callers (current line numbers of the Python scripts):
//! - `binom.sf(k, n, p)`: `acmg_pm2_bs1_ba1_frequency.py` 145 and 643;
//! - `binom.cdf(k, n, p)`: `acmg_pm2_bs1_ba1_frequency.py` 240, `acmg_scoring.py` 451 and 459;
//! - `beta.ppf(q, a, b)`: `acmg_scoring.py` 479;
//! - `chi2_contingency(table)[0:2]` and `fisher_exact(table, alternative="greater")`:
//!   `acmg_pm1_pp2_bp1_missense.py` 206 and 212, on 2x2 tables of counts.
//!
//! SciPy is called on numpy arrays there; these functions take one element,
//! and the caller maps them over the column (NaN in gives NaN out, as SciPy does).

/// `scipy.stats.binom.sf(k, n, p)`: P(X > k) for X ~ Binomial(n, p).
pub(crate) fn binom_sf(k: f64, n: f64, p: f64) -> f64 {
    let _ = (k, n, p);
    unimplemented!("binom_sf is filled in by the frequency thread (T5)")
}

/// `scipy.stats.binom.cdf(k, n, p)`: P(X <= k) for X ~ Binomial(n, p).
pub(crate) fn binom_cdf(k: f64, n: f64, p: f64) -> f64 {
    let _ = (k, n, p);
    unimplemented!("binom_cdf is filled in by the frequency thread (T5)")
}

/// `scipy.stats.beta.ppf(q, a, b)`: the q-quantile of Beta(a, b).
pub(crate) fn beta_ppf(q: f64, a: f64, b: f64) -> f64 {
    let _ = (q, a, b);
    unimplemented!("beta_ppf is filled in by the frequency thread (T5)")
}

/// `scipy.stats.chi2_contingency(table)[0:2]` for a 2x2 table with the
/// default Yates continuity correction: `(statistic, p-value)`.
pub(crate) fn chi2_contingency_2x2(table: [[i64; 2]; 2]) -> (f64, f64) {
    let _ = table;
    unimplemented!("chi2_contingency_2x2 is filled in by the frequency thread (T5)")
}

/// `scipy.stats.fisher_exact(table, alternative="greater")` for a 2x2 table:
/// `(odds ratio, p-value)`.
pub(crate) fn fisher_exact_greater(table: [[i64; 2]; 2]) -> (f64, f64) {
    let _ = table;
    unimplemented!("fisher_exact_greater is filled in by the frequency thread (T5)")
}
