//! Corresponds to acmg_scoring.calculate_posterior_probability, lines 300–387.
//! Qualified domain: validated integer strengths (0..=5), DEFAULT prior=0.1,
//! exp_base=2, odds_pvst=350. No alternative prior/model is claimed.

#[derive(Debug, PartialEq)]
pub enum Probability {
    Float(f64),
    Integer(i64),
}

impl Probability {
    pub fn as_f64(&self) -> f64 {
        match self {
            Self::Float(value) => *value,
            Self::Integer(value) => *value as f64,
        }
    }
}

pub fn calculate_posterior_probability(
    criteria_assignment: &[(&str, i64)],
) -> (Probability, &'static str) {
    let patho_matches: Vec<i64> = criteria_assignment
        .iter()
        .filter(|(name, n)| name.starts_with('P') && *n > 0)
        .map(|(_, n)| *n)
        .collect();
    let ps_matches: Vec<i64> = patho_matches.iter().copied().filter(|n| *n == 3).collect();
    let pm_matches: Vec<i64> = patho_matches.iter().copied().filter(|n| *n == 2).collect();
    let pp_matches: Vec<i64> = patho_matches.iter().copied().filter(|n| *n == 1).collect();
    let pvs_matches: Vec<i64> = patho_matches.iter().copied().filter(|n| *n == 4).collect();
    let benign_matches: Vec<i64> = criteria_assignment
        .iter()
        .filter(|(name, n)| name.starts_with('B') && *n > 0)
        .map(|(_, n)| *n)
        .collect();
    let bs_matches: Vec<i64> = benign_matches.iter().copied().filter(|n| *n == 3).collect();
    let bp_matches: Vec<i64> = benign_matches.iter().copied().filter(|n| *n == 1).collect();
    let ba_matches: Vec<i64> = benign_matches.iter().copied().filter(|n| *n == 5).collect();
    let bm_matches: Vec<i64> = benign_matches.iter().copied().filter(|n| *n == 2).collect();
    let weight = |n: &i64| match n {
        4 => 1.0,
        3 => 0.5,
        2 => 0.25,
        1 => 0.125,
        5 => 2.0,
        _ => unreachable!(),
    };
    let power = |n: f64| 350.0_f64.powf(n);

    // Preserve category concatenation and separate exponentiation. Combining
    // weights into one exponent changes floating-point evaluation order.
    let odds_benign = power(
        -bs_matches
            .iter()
            .chain(&bp_matches)
            .chain(&ba_matches)
            .chain(&bm_matches)
            .map(weight)
            .sum::<f64>(),
    );
    let ps_odds = if ps_matches.len() >= 2 {
        1514.0 * power(ps_matches[..ps_matches.len() - 2].iter().map(weight).sum())
    } else {
        power(ps_matches.iter().map(weight).sum())
    };
    let pv_pm_odds = if !pvs_matches.is_empty() && !pm_matches.is_empty() && pp_matches.is_empty() {
        350.0 * power(pm_matches[..pm_matches.len() - 1].iter().map(weight).sum())
    } else {
        power(pvs_matches.iter().chain(&pm_matches).map(weight).sum())
    };
    let pp_odds = power(pp_matches.iter().map(weight).sum());
    let odds_path = ps_odds * pp_odds * pv_pm_odds;
    let odds_path = odds_path * odds_benign;
    let posterior_probability = (odds_path * 0.1) / ((odds_path - 1.0) * 0.1 + 1.0);
    // round(np.float64, 4) delegates to NumPy scalar.round: scale, ties-to-even,
    // unscale. Rust f64::round instead rounds half values away from zero.
    let mut probability =
        Probability::Float((posterior_probability * 10_000.0).round_ties_even() / 10_000.0);
    if !ba_matches.is_empty() {
        // Source assigns Python int(0) here. Preserve this distinction for
        // whole-column pandas inference during final serialization.
        probability = Probability::Integer(0);
    }
    let posterior_probability = probability.as_f64();
    let acmg_class = if posterior_probability > 0.99 {
        "Pathogenic"
    } else if posterior_probability > 0.899 {
        "Likely Pathogenic"
    } else if posterior_probability < 0.0005 {
        "Benign"
    } else if posterior_probability < 0.05 {
        "Likely Benign"
    } else {
        "Uncertain Significance"
    };
    (probability, acmg_class)
}
