//! Provisional combined evidence.
//!
//! In log10 units:
//!
//! * plus how badly the artifact models fail to explain the reads, capped
//!   at [`CAP`] (the cap is also what a sequence gets that has no more
//!   abundant sequence to be an artifact of, or lies far from it and does
//!   not recur);
//! * plus log10 P(this few reads | one copy), which is 0 for a full copy's
//!   worth and strongly negative well below it;
//! * minus [`CHIMERA_PENALTY`] for a chimera of two more abundant sequences
//!   of the donor.
//!
//! The form and the constants are a first attempt to be judged on plates
//! with known truth, not a calibrated probability.

use super::Observation;

pub const CAP: f64 = 30.0;
pub const CHIMERA_PENALTY: f64 = 10.0;

pub fn score(log10_p_copy: f64, log10_artifacts: Option<f64>, chimera: bool) -> f64 {
    let unexplained = log10_artifacts.map_or(CAP, |e| (-e).clamp(0.0, CAP));
    unexplained + log10_p_copy.max(-CAP) - if chimera { CHIMERA_PENALTY } else { 0.0 }
}

/// Fill in `score` wherever the copy tail is known.
pub fn set(observations: &mut [Observation]) {
    for o in observations.iter_mut() {
        o.evidence.score = o
            .evidence
            .log10_p_copy
            .map(|p| score(p, o.evidence.log10_artifacts, o.evidence.chimera));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn full_copy_that_no_artifact_explains_scores_highest() {
        assert_eq!(score(0.0, None, false), 30.0);
        assert_eq!(score(-0.3, Some(-80.0), false), 29.7);
        // Explained by errors: nothing is added.
        assert_eq!(score(-0.3, Some(1.5), false), -0.3);
        // Far too few reads for a copy cancels "unexplained".
        assert_eq!(score(-45.0, None, false), 0.0);
        assert_eq!(score(0.0, None, true), 20.0);
    }
}
