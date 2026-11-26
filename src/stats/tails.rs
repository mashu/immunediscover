//! Tail probabilities in log10, accurate far into the tails.

use std::f64::consts::{LN_10, PI};

/// log10 of the standard normal lower tail Φ(z).
pub fn log10_normal_cdf(z: f64) -> f64 {
    if z < -6.0 {
        // Asymptotic expansion of the Mills ratio.
        let w = 1.0 / (z * z);
        let series = (-w + 3.0 * w * w - 15.0 * w * w * w).ln_1p();
        return (-z * z / 2.0 - (-z).ln() - 0.5 * (2.0 * PI).ln() + series) / LN_10;
    }
    (0.5 * libm::erfc(-z / std::f64::consts::SQRT_2)).log10()
}

/// Natural log of the probability of `k` under a negative binomial with
/// mean `mean` and shape `shape`; a Poisson when the shape is infinite.
fn ln_negative_binomial(k: f64, mean: f64, shape: f64) -> f64 {
    if shape.is_infinite() {
        return k * mean.ln() - mean - libm::lgamma(k + 1.0);
    }
    libm::lgamma(k + shape) - libm::lgamma(shape) - libm::lgamma(k + 1.0)
        + shape * (shape / (shape + mean)).ln()
        + k * (mean / (shape + mean)).ln()
}

/// log10 P(X >= k) for a negative binomial (or Poisson) count.
pub fn log10_negative_binomial_upper_tail(k: u64, mean: f64, shape: f64) -> f64 {
    if k == 0 {
        return 0.0;
    }
    if mean <= 0.0 {
        return f64::NEG_INFINITY;
    }
    // Sum the tail relative to its first term; the terms decrease from the
    // mode onwards and the sum stops when they no longer matter.
    let first = ln_negative_binomial(k as f64, mean, shape);
    let mut total = 1.0;
    for i in (k + 1)..(k + 200_000) {
        let term = (ln_negative_binomial(i as f64, mean, shape) - first).exp();
        total += term;
        if term < 1e-15 * total {
            break;
        }
    }
    ((first + total.ln()) / LN_10).min(0.0)
}

/// log10 P(X <= k) for a Poisson count with mean `mean`.
pub fn log10_poisson_cdf(k: u64, mean: f64) -> f64 {
    if mean <= 0.0 {
        return 0.0;
    }
    let terms: Vec<f64> = (0..=k)
        .map(|i| i as f64 * mean.ln() - mean - libm::lgamma(i as f64 + 1.0))
        .collect();
    let largest = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let sum: f64 = terms.iter().map(|t| (t - largest).exp()).sum();
    ((largest + sum.ln()) / LN_10).min(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tolerance: f64) -> bool {
        (a - b).abs() <= tolerance * b.abs().max(1.0)
    }

    #[test]
    fn normal_tail() {
        assert!(close(log10_normal_cdf(0.0), 0.5f64.log10(), 1e-12));
        // Φ(-1.96) = 0.0249979, Φ(-10) = 7.6199e-24, Φ(3) = 0.9986501
        assert!(close(log10_normal_cdf(-1.96), 0.024997895f64.log10(), 1e-7));
        assert!(close(
            log10_normal_cdf(-10.0),
            7.619853e-24f64.log10(),
            1e-6
        ));
        assert!(close(log10_normal_cdf(3.0), 0.998650102f64.log10(), 1e-8));
        // Continuous across the switch to the asymptotic form.
        assert!((log10_normal_cdf(-6.0) - log10_normal_cdf(-6.000001)).abs() < 1e-4);
    }

    #[test]
    fn poisson_tails_add_up() {
        // P(X <= 2 | 3) = e^-3 (1 + 3 + 4.5)
        let expected = (-3.0f64).exp() * 8.5;
        assert!(close(log10_poisson_cdf(2, 3.0), expected.log10(), 1e-12));
        let upper = 10f64.powf(log10_negative_binomial_upper_tail(3, 3.0, f64::INFINITY));
        assert!(close(upper, 1.0 - expected, 1e-12));
        assert_eq!(log10_negative_binomial_upper_tail(0, 3.0, 2.0), 0.0);
    }

    #[test]
    fn negative_binomial_tail_matches_direct_summation() {
        // Shape 2, mean 4: p = 1/3, pmf(k) = (k + 1) p^2 (1 - p)^k
        let p: f64 = 1.0 / 3.0;
        let below: f64 = (0..5)
            .map(|k| (k as f64 + 1.0) * p * p * (1.0 - p).powi(k))
            .sum();
        let tail = 10f64.powf(log10_negative_binomial_upper_tail(5, 4.0, 2.0));
        assert!(close(tail, 1.0 - below, 1e-10));
        // Heavier tail than a Poisson of the same mean.
        assert!(
            log10_negative_binomial_upper_tail(20, 4.0, 2.0)
                > log10_negative_binomial_upper_tail(20, 4.0, f64::INFINITY)
        );
    }
}
