//! Is this donor's amount of a recurrent artifact out of the ordinary?
//!
//! Some artifacts are not sporadic: the same variant of the same parent
//! turns up in most donors, at a low and fairly steady ratio to that
//! parent. Such a sequence is not judged against the sporadic artifacts of
//! the plate but against itself: the reads expected here are the parent's
//! reads times the usual ratio in the other donors, and the question is
//! whether this donor holds so much more that an ordinary fluctuation does
//! not explain it. A real allele whose sequence is also such an artifact
//! does: its carriers hold a hundred times the usual ratio.
//!
//! How much the ratio varies between donors is measured on the plate (a
//! robust spread of the log ratio, so that carriers of real alleles do not
//! widen it), and sampling noise is added for low counts. The evidence is
//! -log10 of the number of recurrent artifacts on the plate expected to
//! reach this far above their usual ratio.
//!
//! Only sequences whose usual ratio is at most
//! [`MAX_ARTIFACT_RATIO`](super::artifacts::MAX_ARTIFACT_RATIO) of the
//! parent are treated this way. A sequence that most donors hold in a
//! larger amount is not an artifact of its parent; whether it segregates
//! is reported by `ubiquitous`.

use super::artifacts::MAX_ARTIFACT_RATIO;
use super::Observation;
use crate::stats::robust::mad;
use crate::stats::tails::log10_normal_cdf;

/// Reads a recurrent artifact must be expected to give for its ratio to
/// say something about the spread between donors.
const MIN_EXPECTED_READS: f64 = 5.0;
/// Recurrent artifacts needed for the spread to be measured.
const MIN_ROWS: usize = 30;
const DEFAULT_SPREAD: f64 = 0.5;
const MIN_SPREAD: f64 = 0.2;
const MAX_EVIDENCE: f64 = 9.0;

/// The variation of recurrent artifacts around their usual ratio.
#[derive(Debug, Clone)]
pub struct Recurrence {
    /// SD of the log ratio between donors.
    pub spread: f64,
    /// Recurrent artifacts on the plate (sequence-in-donor rows).
    pub rows: usize,
}

/// Reads expected from the usual ratio, for a recurrent artifact.
fn expected_reads(observations: &[Observation], o: &Observation) -> Option<f64> {
    let (parent, _) = o.evidence.parent?;
    let ratio = o
        .evidence
        .recurrent_ratio
        .filter(|r| *r > 0.0 && *r <= MAX_ARTIFACT_RATIO)?;
    Some(observations[parent].count as f64 * ratio)
}

impl Recurrence {
    pub fn fit(observations: &[Observation]) -> Self {
        let expected: Vec<(f64, f64)> = observations
            .iter()
            .filter(|o| !o.is_database_allele())
            .filter_map(|o| Some((o.count as f64, expected_reads(observations, o)?)))
            .collect();
        let log_ratios: Vec<f64> = expected
            .iter()
            .filter(|(_, mean)| *mean >= MIN_EXPECTED_READS)
            .map(|(count, mean)| (count / mean).ln())
            .collect();
        let spread = if log_ratios.len() >= MIN_ROWS {
            mad(&log_ratios).max(MIN_SPREAD)
        } else {
            DEFAULT_SPREAD
        };
        Recurrence {
            spread,
            rows: expected.len(),
        }
    }

    /// Evidence that `count` reads are more than a fluctuation around
    /// `expected` reads.
    pub fn evidence(&self, count: u64, expected: f64) -> f64 {
        if count as f64 <= expected {
            return 0.0;
        }
        // Spread between donors plus sampling noise of the count.
        let sd = (self.spread.powi(2) + 1.0 / expected).sqrt();
        let z = (count as f64 / expected).ln() / sd;
        let log10_expected_rows = (self.rows.max(1) as f64).log10() + log10_normal_cdf(-z);
        (-log10_expected_rows).clamp(0.0, MAX_EVIDENCE) + 0.0
    }

    /// Lower `score` where a recurrent artifact explains the observation.
    pub fn limit_scores(&self, observations: &mut [Observation]) {
        let limits: Vec<Option<f64>> = observations
            .iter()
            .map(|o| Some(self.evidence(o.count, expected_reads(observations, o)?)))
            .collect();
        for (o, limit) in observations.iter_mut().zip(limits) {
            if let (Some(score), Some(limit)) = (o.evidence.score, limit) {
                o.evidence.score = Some(score.min(limit));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::testing::observation;

    /// 40 donors with 1000 reads of a parent and a variant at about 5 % of
    /// it; in the last donor the variant has `last` reads.
    fn plate(last: u64) -> Vec<Observation> {
        let mut observations = Vec::new();
        for donor in 0..40 {
            observations.push(observation(donor, "A", "AAAA", 1000));
            let count = if donor == 39 {
                last
            } else {
                40 + (donor as u64 % 5) * 5
            };
            let mut variant = observation(donor, "A", "AAAT", count);
            variant.db_distance = 1;
            variant.evidence.parent = Some((2 * donor, 1));
            variant.evidence.recurrent_ratio = Some(0.05);
            variant.evidence.score = Some(3.0);
            observations.push(variant);
        }
        observations
    }

    #[test]
    fn a_fluctuation_is_no_evidence_and_a_carrier_is() {
        let mut observations = plate(150);
        let model = Recurrence::fit(&observations);
        assert_eq!(model.rows, 40);
        assert_eq!(model.spread, MIN_SPREAD); // counts of 40 to 60 around 50
                                              // Three times the usual ratio happens among 40 rows.
        assert!(model.evidence(150, 50.0) > 2.0);
        assert_eq!(model.evidence(70, 50.0), 0.0);
        assert_eq!(model.evidence(30, 50.0), 0.0);
        // A carrier holds as much as the parent: twenty times the usual.
        assert_eq!(model.evidence(1000, 50.0), MAX_EVIDENCE);
        model.limit_scores(&mut observations);
        // Ordinary rows lose their score, the parent rows are untouched.
        assert_eq!(observations[1].evidence.score, Some(0.0));
        assert_eq!(observations[0].evidence.score, None);
        assert_eq!(observations[79].evidence.score, Some(3.0));
    }

    #[test]
    fn few_reads_are_noisy() {
        let model = Recurrence {
            spread: 0.3,
            rows: 2000,
        };
        // Three times the usual ratio is nothing among 2000 rows.
        assert_eq!(model.evidence(9, 3.0), 0.0);
        assert_eq!(model.evidence(90, 30.0), 0.0);
        // Ten times is, and more so where the counts are not small.
        let (few, many) = (model.evidence(30, 3.0), model.evidence(300, 30.0));
        assert!(few > 0.0 && few < 1.0, "{few}");
        assert!(many > 5.0, "{many}");
    }

    #[test]
    fn a_sequence_most_donors_hold_in_a_larger_amount_is_not_treated() {
        let mut observations = plate(50);
        for o in observations.iter_mut() {
            if o.evidence.recurrent_ratio.is_some() {
                o.evidence.recurrent_ratio = Some(0.45);
            }
        }
        let model = Recurrence::fit(&observations);
        assert_eq!(model.rows, 0);
        model.limit_scores(&mut observations);
        assert_eq!(observations[1].evidence.score, Some(3.0));
    }
}
