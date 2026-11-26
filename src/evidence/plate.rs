//! How a sequence behaves across the donors of the plate.
//!
//! Being seen in many donors is not, by itself, evidence: a recurrent
//! artifact is seen wherever its parent is. What tells the two apart is
//! the pattern. A real allele is carried by some donors in a copy-like
//! amount and absent from others, also where the parent is present. A
//! recurrent artifact follows its parent everywhere, at a steady low ratio.

use std::collections::HashMap;

use super::Observation;
use crate::stats::median;

/// Dosage from which a donor counts as carrying a sequence.
pub const COPY_LIKE: f64 = 0.25;
/// Typical dosage above which a sequence has more copies than one gene on
/// two chromosomes can hold.
pub const MULTI_COPY: f64 = 2.8;
/// Share of the donors from which a sequence counts as not segregating.
pub const UBIQUITOUS_SHARE: f64 = 0.9;
/// Two ratios to the parent within this factor are the same level.
const SAME_LEVEL: f64 = 3.0;
/// Dosage from which a donor counts as carrying the parent.
const PARENT_CARRIER: f64 = 0.5;
/// Parent reads another donor needs for its ratio to be informative.
const MIN_PARENT_READS: u64 = 20;
/// Other donors needed before a usual ratio is reported.
const MIN_OTHER_DONORS: usize = 3;

/// Whether no other donor of the plate carries the sequence.
pub fn is_private(o: &Observation) -> bool {
    let own = o.evidence.dosage.is_some_and(|d| d >= COPY_LIKE);
    o.evidence.carriers == usize::from(own)
}

/// Fill in `carriers`, `tracking`, `recurrent_ratio`, `steady` and
/// `ubiquitous`. Needs dosage and parents.
pub fn set_context(observations: &mut [Observation], n_donors: usize) {
    // sequence => donor => (reads, dosage)
    let mut seen: HashMap<&str, HashMap<usize, (u64, Option<f64>)>> = HashMap::new();
    for o in observations.iter() {
        seen.entry(o.sequence.as_str())
            .or_default()
            .insert(o.donor, (o.count, o.evidence.dosage));
    }
    let mut context = Vec::with_capacity(observations.len());
    for o in observations.iter() {
        let own = &seen[o.sequence.as_str()];
        let carriers = own
            .values()
            .filter(|(_, dosage)| dosage.is_some_and(|d| d >= COPY_LIKE))
            .count();
        let (mut tracking, mut recurrent, mut steady) = (None, None, None);
        if let Some((parent, _)) = o.evidence.parent {
            let of_parent = &seen[observations[parent].sequence.as_str()];
            let parent_carriers: Vec<usize> = of_parent
                .iter()
                .filter(|(_, (_, dosage))| dosage.is_some_and(|d| d >= PARENT_CARRIER))
                .map(|(donor, _)| *donor)
                .collect();
            if !parent_carriers.is_empty() {
                let with_it = parent_carriers
                    .iter()
                    .filter(|d| own.contains_key(d))
                    .count();
                tracking = Some(with_it as f64 / parent_carriers.len() as f64);
            }
            // Ratio to the parent in every OTHER donor with enough parent reads.
            let ratios: Vec<f64> = of_parent
                .iter()
                .filter(|(donor, (reads, _))| **donor != o.donor && *reads >= MIN_PARENT_READS)
                .map(|(donor, (reads, _))| own.get(donor).map_or(0, |x| x.0) as f64 / *reads as f64)
                .collect();
            if ratios.len() >= MIN_OTHER_DONORS {
                recurrent = Some(median(&ratios));
                let here = o.count as f64 / observations[parent].count as f64;
                let same_level = ratios
                    .iter()
                    .filter(|r| (here / SAME_LEVEL..=here * SAME_LEVEL).contains(*r))
                    .count();
                steady = Some(same_level as f64 / ratios.len() as f64);
            }
        }
        context.push((carriers, tracking, recurrent, steady));
    }
    for (o, (carriers, tracking, recurrent, steady)) in observations.iter_mut().zip(context) {
        o.evidence.carriers = carriers;
        o.evidence.tracking = tracking;
        o.evidence.recurrent_ratio = recurrent;
        o.evidence.steady = steady;
        o.evidence.ubiquitous = carriers as f64 >= UBIQUITOUS_SHARE * n_donors as f64
            || steady.is_some_and(|s| s >= UBIQUITOUS_SHARE);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::parents;
    use crate::evidence::testing::observation;

    /// Five donors all carry the parent (100 reads, dosage 1). `hotspot`
    /// follows it at 5 % everywhere; `allele` is carried by two donors only.
    fn plate() -> Vec<Observation> {
        let mut observations = Vec::new();
        for donor in 0..5 {
            let mut add = |sequence: &str, count: u64| {
                let mut o = observation(donor, "A", sequence, count);
                o.evidence.dosage = Some(count as f64 / 100.0);
                observations.push(o);
            };
            add("AAAAAAAA", 100);
            add("AAAAAAAC", 5);
            if donor < 2 {
                add("AAAAAAAG", 90);
            }
        }
        parents::assign(&mut observations);
        set_context(&mut observations, 5);
        observations
    }

    #[test]
    fn recurrent_artifact_follows_its_parent_everywhere() {
        let observations = plate();
        let hotspot = observations
            .iter()
            .find(|o| o.sequence == "AAAAAAAC")
            .unwrap();
        assert_eq!(hotspot.evidence.tracking, Some(1.0));
        assert_eq!(hotspot.evidence.recurrent_ratio, Some(0.05));
        assert_eq!(hotspot.evidence.carriers, 0);
        // The same 5 % of the parent in every other donor: it does not
        // segregate, although no donor carries a copy of it.
        assert_eq!(hotspot.evidence.steady, Some(1.0));
        assert!(hotspot.evidence.ubiquitous);
        // The parent itself is carried by every donor.
        assert!(observations[0].evidence.ubiquitous);
    }

    #[test]
    fn real_allele_is_carried_by_some_donors_and_absent_from_the_others() {
        let observations = plate();
        let allele = observations
            .iter()
            .find(|o| o.sequence == "AAAAAAAG")
            .unwrap();
        assert_eq!(allele.evidence.carriers, 2);
        assert_eq!(allele.evidence.tracking, Some(0.4));
        // In most of the other donors with the parent it is not there at all.
        assert_eq!(allele.evidence.recurrent_ratio, Some(0.0));
        // One of the four other donors holds it at the same level.
        assert_eq!(allele.evidence.steady, Some(0.25));
        assert!(!allele.evidence.ubiquitous);
        assert!(!is_private(allele));
    }
}
