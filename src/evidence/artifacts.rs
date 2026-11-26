//! How many artifacts as abundant as an observation are expected.
//!
//! Two kinds of artifact are told apart, each with a rate learned from the
//! plate:
//!
//! * **Sporadic errors** hit any position. The rate is that of reads of a
//!   sequence turning into ONE specific variant at a given distance;
//!   variants never observed count as zeros, and the spread between
//!   variants is a negative binomial shape. An observation is compared with
//!   the number of such variants expected at its abundance or above, over
//!   all variants at that distance.
//! * **Recurrent artifacts** are the same variant of the same parent in
//!   donor after donor (error hotspots, paralogous products). Such a
//!   variant is compared with its own usual ratio to the parent, learned
//!   from the other donors.
//!
//! A sequence above [`MAX_ARTIFACT_RATIO`] of its parent is never used to
//! learn a rate, so real alleles do not inflate the estimates.

use super::parents::DISTANCE_LIMIT;
use super::Observation;
use crate::stats::tails::log10_negative_binomial_upper_tail;

/// Reads a sequence needs to serve as an error source when fitting.
const MIN_SOURCE_READS: u64 = 50;
/// Above this share of the parent a sequence is not treated as an artifact.
pub const MAX_ARTIFACT_RATIO: f64 = 0.1;
/// Largest distance with a sporadic error fit.
pub const MAX_SPORADIC_DISTANCE: usize = 3;
/// Shapes below this describe a fit to a handful of outliers, not a rate.
const MIN_SHAPE: f64 = 0.05;

/// Number of sequences at `distance` substitutions from one of `length`.
pub fn variants_at(length: usize, distance: usize) -> f64 {
    (0..distance).fold(1.0, |n, i| n * (length - i) as f64 / (i + 1) as f64)
        * 3f64.powi(distance as i32)
}

/// Sporadic error fit at one distance.
#[derive(Debug, Clone, Copy)]
pub struct SporadicFit {
    pub distance: usize,
    /// Reads of a specific variant per read of the source.
    pub rate: f64,
    /// Negative binomial shape; infinite for a Poisson.
    pub shape: f64,
    pub sources: usize,
    pub satellites: usize,
    pub satellite_reads: u64,
}

/// Between-donor spread of recurrent artifacts.
#[derive(Debug, Clone, Copy)]
pub struct RecurrentFit {
    pub shape: f64,
    pub observations: usize,
}

fn recurs(o: &Observation) -> Option<f64> {
    o.evidence.recurrent_ratio.filter(|r| *r > 0.0)
}

/// Is `o` a satellite usable for the sporadic fit at `distance`?
pub fn is_sporadic_satellite(
    observations: &[Observation],
    o: &Observation,
    distance: usize,
) -> bool {
    let Some((parent, d)) = o.evidence.parent else {
        return false;
    };
    let source = observations[parent].count;
    d == distance
        && !o.is_database_allele()
        && source >= MIN_SOURCE_READS
        && o.count as f64 <= MAX_ARTIFACT_RATIO * source as f64
        && recurs(o).is_none()
}

/// Fit the sporadic error rate for 1..=[`MAX_SPORADIC_DISTANCE`] differences.
pub fn fit_sporadic(observations: &[Observation]) -> Vec<SporadicFit> {
    let sources: Vec<usize> = (0..observations.len())
        .filter(|&i| observations[i].count >= MIN_SOURCE_READS)
        .collect();
    (1..=MAX_SPORADIC_DISTANCE)
        .map(|distance| {
            let mut satellites_of: Vec<Vec<u64>> = vec![Vec::new(); observations.len()];
            for o in observations {
                if is_sporadic_satellite(observations, o, distance) {
                    let (parent, _) = o.evidence.parent.expect("satellites have a parent");
                    satellites_of[parent].push(o.count);
                }
            }
            let (mut exposure, mut exposure_squared, mut reads) = (0.0, 0.0, 0u64);
            for &i in &sources {
                let n = observations[i].count as f64;
                let variants = variants_at(observations[i].sequence.len(), distance);
                exposure += n * variants;
                exposure_squared += n * n * variants;
                reads += satellites_of[i].iter().sum::<u64>();
            }
            let rate = if exposure > 0.0 {
                reads as f64 / exposure
            } else {
                0.0
            };
            // Method of moments over every possible variant of every source:
            // Var = mean + mean^2 / shape.
            let mut squares = 0.0;
            for &i in &sources {
                let mean = observations[i].count as f64 * rate;
                let variants = variants_at(observations[i].sequence.len(), distance);
                let seen = &satellites_of[i];
                squares += seen.iter().map(|&k| (k as f64 - mean).powi(2)).sum::<f64>()
                    + (variants - seen.len() as f64) * mean * mean;
            }
            let excess = squares - exposure * rate;
            let shape = if rate > 0.0 && excess > 0.0 {
                (rate * rate * exposure_squared / excess).max(MIN_SHAPE)
            } else {
                f64::INFINITY
            };
            SporadicFit {
                distance,
                rate,
                shape,
                sources: sources.len(),
                satellites: satellites_of.iter().map(Vec::len).sum(),
                satellite_reads: reads,
            }
        })
        .collect()
}

/// How far the reads of a recurring variant stray, between donors, from
/// parent reads x its usual ratio.
pub fn fit_recurrent(observations: &[Observation]) -> RecurrentFit {
    let (mut squares, mut excess, mut n) = (0.0, 0.0, 0usize);
    for o in observations {
        let (Some((parent, _)), Some(ratio)) = (o.evidence.parent, recurs(o)) else {
            continue;
        };
        if ratio > MAX_ARTIFACT_RATIO || o.is_database_allele() {
            continue;
        }
        let mean = observations[parent].count as f64 * ratio;
        squares += mean * mean;
        excess += (o.count as f64 - mean).powi(2) - mean;
        n += 1;
    }
    RecurrentFit {
        shape: if excess > 0.0 {
            (squares / excess).max(0.5)
        } else {
            f64::INFINITY
        },
        observations: n,
    }
}

/// Fill in `log10_artifacts`.
///
/// A variant that recurs across donors is tested against its own usual
/// ratio (no correction for the number of possible variants is needed
/// then). Any other variant is tested against the sporadic rate, times the
/// number of variants at its distance. A sequence without a parent, or far
/// from it and not recurrent, has no artifact expectation.
pub fn set_expectation(
    observations: &mut [Observation],
    sporadic: &[SporadicFit],
    recurrent: &RecurrentFit,
) {
    let expectations: Vec<Option<f64>> = observations
        .iter()
        .map(|o| {
            let (parent, distance) = o.evidence.parent?;
            let source = &observations[parent];
            let n = source.count as f64;
            let fit = sporadic
                .iter()
                .find(|f| f.distance == distance && f.rate > 0.0);
            let sporadic_rate = fit.map_or(0.0, |f| f.rate);
            if let Some(ratio) = recurs(o) {
                let mean = n * (ratio.min(MAX_ARTIFACT_RATIO) + sporadic_rate);
                return Some(log10_negative_binomial_upper_tail(
                    o.count,
                    mean,
                    recurrent.shape,
                ));
            }
            let fit = fit?;
            debug_assert!(distance <= DISTANCE_LIMIT);
            Some(
                variants_at(source.sequence.len(), distance).log10()
                    + log10_negative_binomial_upper_tail(o.count, n * fit.rate, fit.shape),
            )
        })
        .collect();
    for (o, expectation) in observations.iter_mut().zip(expectations) {
        o.evidence.log10_artifacts = expectation;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::testing::observation;
    use crate::evidence::{parents, plate};

    #[test]
    fn variant_counts() {
        assert_eq!(variants_at(10, 1), 30.0);
        assert_eq!(variants_at(10, 2), 45.0 * 9.0);
        assert_eq!(variants_at(300, 3), 4_455_100.0 * 27.0);
    }

    /// Twenty donors, each with a 1000-read sequence of length 20 and
    /// `errors` single-read variants of it at distinct positions.
    fn plate_with_errors(errors: usize) -> Vec<Observation> {
        let base = "A".repeat(20);
        let mut observations = Vec::new();
        for donor in 0..20 {
            let mut source = observation(donor, "A", &base, 1000);
            source.evidence.dosage = Some(1.0);
            observations.push(source);
            for e in 0..errors {
                // different variants in different donors: nothing recurs
                let position = (donor + e) % 20;
                let letter = ["C", "G", "T"][(donor / 20 + e / 20) % 3];
                let mut variant = base.clone();
                variant.replace_range(position..=position, letter);
                let mut o = observation(donor, "A", &variant, 1);
                o.db_distance = 1;
                observations.push(o);
            }
        }
        parents::assign(&mut observations);
        observations
    }

    #[test]
    fn sporadic_rate_is_reads_per_possible_variant() {
        let observations = plate_with_errors(6);
        let fits = fit_sporadic(&observations);
        // 6 error reads per 1000 reads over 60 possible variants.
        assert!((fits[0].rate - 6.0 / 1000.0 / 60.0).abs() < 1e-12);
        assert_eq!((fits[0].sources, fits[0].satellites), (20, 120));
        assert_eq!(fits[1].rate, 0.0);
    }

    #[test]
    fn abundant_variant_is_not_explained_but_a_single_read_is() {
        let mut observations = plate_with_errors(6);
        // In donor 0 one more variant has 300 reads: too many, relative to
        // its parent, to be used for learning the error rate.
        let mut candidate = observation(0, "A", &("A".repeat(19) + "G"), 300);
        candidate.db_distance = 1;
        observations.push(candidate);
        parents::assign(&mut observations);
        plate::set_context(&mut observations);
        let sporadic = fit_sporadic(&observations);
        let recurrent = fit_recurrent(&observations);
        set_expectation(&mut observations, &sporadic, &recurrent);
        let expectation = |o: &Observation| o.evidence.log10_artifacts;
        let candidate = observations.last().unwrap();
        assert!(expectation(candidate).unwrap() < -50.0);
        // A single-read variant: several such are expected per source.
        let single = observations.iter().find(|o| o.count == 1).unwrap();
        assert!(expectation(single).unwrap() > 0.0);
        // The source itself has no parent, hence no artifact expectation.
        assert_eq!(expectation(&observations[0]), None);
    }
}
