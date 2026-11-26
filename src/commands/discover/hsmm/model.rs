//! Model components of the D-gene hidden semi-Markov model.
//!
//! The model describes a D gene embedded in its two recombination signals:
//!
//! ```text
//! nonamer(9) spacer(12) heptamer(7)  D gene (variable)  heptamer(7) spacer(12) nonamer(9)
//! ```
//!
//! Signal elements are fixed-length states with position-specific emission
//! probabilities; the gene is a variable-duration state with i.i.d. base
//! emissions and an explicit length distribution.

use crate::error::{Error, Result};
use crate::seq::dna::base_index;

/// Position-specific log-probabilities of a fixed-length motif.
#[derive(Debug, Clone, PartialEq)]
pub struct Pwm<const L: usize> {
    log_p: [[f64; 4]; L],
}

impl<const L: usize> Pwm<L> {
    /// Estimate from equally long example motifs with additive smoothing.
    /// Characters other than A, C, G, T contribute no count.
    pub fn fit<S: AsRef<str>>(motifs: &[S], pseudocount: f64) -> Result<Self> {
        if motifs.is_empty() {
            return Err(Error::invalid("Cannot fit a motif model without examples"));
        }
        let mut counts = [[pseudocount; 4]; L];
        for motif in motifs {
            let motif = motif.as_ref().as_bytes();
            if motif.len() != L {
                return Err(Error::invalid(format!(
                    "Training motif '{}' is {} nt long, expected {L}; \
                     the known D match it came from lies too close to a read end",
                    String::from_utf8_lossy(motif),
                    motif.len()
                )));
            }
            for (position, &base) in motif.iter().enumerate() {
                if let Some(index) = base_index(base) {
                    counts[position][index] += 1.0;
                }
            }
        }
        let mut log_p = [[0.0; 4]; L];
        for (position, row) in counts.iter().enumerate() {
            let inverse = 1.0 / row.iter().sum::<f64>();
            for (base, count) in row.iter().enumerate() {
                log_p[position][base] = (count * inverse).ln();
            }
        }
        Ok(Pwm { log_p })
    }

    /// Log-probability of the window of `observations` starting at `start`
    /// (0-based); negative infinity if it contains an unknown base.
    #[inline]
    pub fn score(&self, observations: &[Option<u8>], start: usize) -> f64 {
        let mut total = 0.0;
        for (position, row) in self.log_p.iter().enumerate() {
            match observations[start + position] {
                Some(base) => total += row[base as usize],
                None => return f64::NEG_INFINITY,
            }
        }
        total
    }

    /// Log-probability of a motif given as text.
    pub fn log_prob(&self, motif: &str) -> f64 {
        debug_assert_eq!(motif.len(), L);
        let observations = encode(motif);
        self.score(&observations, 0)
    }
}

/// Position-independent base log-probabilities.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IidEmission {
    pub log_p: [f64; 4],
}

impl IidEmission {
    /// Estimate base frequencies over all example sequences.
    pub fn fit<S: AsRef<str>>(sequences: &[S], pseudocount: f64) -> Self {
        let mut counts = [pseudocount; 4];
        for sequence in sequences {
            for &base in sequence.as_ref().as_bytes() {
                if let Some(index) = base_index(base) {
                    counts[index] += 1.0;
                }
            }
        }
        let inverse = 1.0 / counts.iter().sum::<f64>();
        IidEmission {
            log_p: counts.map(|c| (c * inverse).ln()),
        }
    }
}

/// Length distribution of a variable-duration state.
pub trait DurationModel {
    fn min_len(&self) -> usize;
    fn max_len(&self) -> usize;
    /// Log-probability of duration `len` (`min_len <= len <= max_len`).
    fn log_pmf(&self, len: usize) -> f64;
}

/// Smoothed empirical length distribution on `min_len..=max_len`.
#[derive(Debug, Clone, PartialEq)]
pub struct EmpiricalDuration {
    min_len: usize,
    max_len: usize,
    log_p: Vec<f64>,
}

impl EmpiricalDuration {
    pub fn fit<S: AsRef<str>>(
        sequences: &[S],
        min_len: usize,
        max_len: usize,
        pseudocount: f64,
    ) -> Result<Self> {
        if min_len > max_len {
            return Err(Error::invalid(format!(
                "Minimum gene length {min_len} exceeds maximum gene length {max_len}"
            )));
        }
        let mut counts = vec![pseudocount; max_len - min_len + 1];
        for sequence in sequences {
            let len = sequence.as_ref().len();
            if (min_len..=max_len).contains(&len) {
                counts[len - min_len] += 1.0;
            }
        }
        let inverse = 1.0 / counts.iter().sum::<f64>();
        Ok(EmpiricalDuration {
            min_len,
            max_len,
            log_p: counts.into_iter().map(|c| (c * inverse).ln()).collect(),
        })
    }
}

impl DurationModel for EmpiricalDuration {
    fn min_len(&self) -> usize {
        self.min_len
    }

    fn max_len(&self) -> usize {
        self.max_len
    }

    #[inline]
    fn log_pmf(&self, len: usize) -> f64 {
        self.log_p[len - self.min_len]
    }
}

/// Normal distribution truncated to `min_len..=max_len` and discretised
/// to integer lengths (upstream's `discretize_truncated_normal`).
///
/// Length `k` receives the truncated normal's mass on `k - 0.5 .. k + 0.5`,
/// floored at machine epsilon, and the masses are renormalised.
#[derive(Debug, Clone, PartialEq)]
pub struct TruncatedNormalDuration {
    min_len: usize,
    max_len: usize,
    log_p: Vec<f64>,
}

impl TruncatedNormalDuration {
    pub fn new(min_len: usize, max_len: usize, mean: f64, sd: f64) -> Result<Self> {
        if !(sd > 0.0 && sd.is_finite()) {
            return Err(Error::invalid(format!(
                "Duration standard deviation must be positive and finite, got {sd}"
            )));
        }
        if min_len > max_len {
            return Err(Error::invalid(format!(
                "Minimum gene length {min_len} exceeds maximum gene length {max_len}"
            )));
        }
        let (lower, upper) = (min_len as f64, max_len as f64);
        let normal_cdf = |x: f64| 0.5 * libm::erfc(-(x - mean) / (sd * std::f64::consts::SQRT_2));
        let lower_cdf = normal_cdf(lower);
        let mass = normal_cdf(upper) - lower_cdf;
        // CDF of the truncated distribution, as Distributions.jl defines it.
        let cdf = |x: f64| {
            if x <= lower {
                0.0
            } else if x >= upper {
                1.0
            } else {
                (normal_cdf(x) - lower_cdf) / mass
            }
        };
        let pmf: Vec<f64> = (min_len..=max_len)
            .map(|k| (cdf(k as f64 + 0.5) - cdf(k as f64 - 0.5)).max(f64::EPSILON))
            .collect();
        let total: f64 = pmf.iter().sum();
        Ok(TruncatedNormalDuration {
            min_len,
            max_len,
            log_p: pmf.into_iter().map(|p| (p / total).ln()).collect(),
        })
    }
}

impl DurationModel for TruncatedNormalDuration {
    fn min_len(&self) -> usize {
        self.min_len
    }

    fn max_len(&self) -> usize {
        self.max_len
    }

    #[inline]
    fn log_pmf(&self, len: usize) -> f64 {
        self.log_p[len - self.min_len]
    }
}

/// Encode a read as base indices; `None` for anything but A, C, G, T.
pub fn encode(sequence: &str) -> Vec<Option<u8>> {
    sequence
        .bytes()
        .map(|b| base_index(b).map(|i| i as u8))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encoding_marks_unknown_bases() {
        assert_eq!(encode("ACGT"), vec![Some(0), Some(1), Some(2), Some(3)]);
        assert_eq!(encode("N"), vec![None]);
    }

    #[test]
    fn pwm_rows_are_normalised_and_score_unknown_as_impossible() {
        let pwm: Pwm<3> = Pwm::fit(&["ACG", "ACT", "AAG"], 0.01).unwrap();
        for row in pwm.log_p {
            let total: f64 = row.iter().map(|l| l.exp()).sum();
            assert!((total - 1.0).abs() < 1e-12);
        }
        assert!(pwm.log_prob("ACG") > pwm.log_prob("TTT"));
        assert_eq!(pwm.log_prob("ANG"), f64::NEG_INFINITY);
        assert!(Pwm::<3>::fit(&["AC"], 0.01).is_err());
    }

    #[test]
    fn iid_emission_with_uniform_counts() {
        let e = IidEmission::fit(&["ACGT"], 0.0);
        for l in e.log_p {
            assert!((l - 0.25f64.ln()).abs() < 1e-15);
        }
    }

    #[test]
    fn empirical_duration_prefers_observed_lengths() {
        let d = EmpiricalDuration::fit(&["AAAA", "AAAA", "AAAAAA"], 3, 8, 0.01).unwrap();
        assert!(d.log_pmf(4) > d.log_pmf(6));
        assert!(d.log_pmf(6) > d.log_pmf(5));
        let total: f64 = (3..=8).map(|l| d.log_pmf(l).exp()).sum();
        assert!((total - 1.0).abs() < 1e-12);
    }

    // Reference: `discretize_truncated_normal` of immunediscover.jl
    // (Distributions.jl 0.25).
    #[test]
    fn truncated_normal_duration_matches_upstream() {
        let d = TruncatedNormalDuration::new(10, 20, 15.0, 3.0).unwrap();
        let expected = [
            -3.861967699798511,
            -2.802411924622478,
            -2.417101410493167,
            -2.1418827388911152,
            -1.976752789258199,
            -1.9217096819679351,
            -1.976752789258199,
            -2.1418827388911152,
            -2.417101410493167,
            -2.8024119246224797,
            -3.8619676997985084,
        ];
        for (len, want) in (10..=20).zip(expected) {
            assert!((d.log_pmf(len) / want - 1.0).abs() < 1e-12, "{len}");
        }
        // Far tail: masses below machine epsilon are floored, as upstream.
        let tail = TruncatedNormalDuration::new(5, 8, 30.0, 1.0).unwrap();
        let expected = [
            -36.04365338911715,
            -34.190705066153086,
            -11.147382952901559,
            -1.4413061033419906e-5,
        ];
        for (len, want) in (5..=8).zip(expected) {
            assert!((tail.log_pmf(len) / want - 1.0).abs() < 1e-9, "{len}");
        }
        assert!(TruncatedNormalDuration::new(5, 8, 6.0, 0.0).is_err());
        assert!(TruncatedNormalDuration::new(9, 8, 6.0, 1.0).is_err());
    }
}
