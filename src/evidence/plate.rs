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
const COPY_LIKE: f64 = 0.25;
/// Dosage from which a donor counts as carrying the parent.
const PARENT_CARRIER: f64 = 0.5;
/// Parent reads another donor needs for its ratio to be informative.
const MIN_PARENT_READS: u64 = 20;
/// Other donors needed before a usual ratio is reported.
const MIN_OTHER_DONORS: usize = 3;

/// Fill in `carriers`, `tracking` and `recurrent_ratio`. Needs dosage and
/// parents.
pub fn set_context(observations: &mut [Observation]) {
    // sequence => donor => (reads, dosage)
    let mut seen: HashMap<&str, HashMap<usize, (u64, Option<f64>)>> = HashMap::new();
    for o in observations.iter() {
        seen.entry(o.sequence.as_str())
            .or_default()
            .insert(o.donor, (o.count, o.evidence.dosage));
    }
    let mut context: Vec<(usize, Option<f64>, Option<f64>)> =
        Vec::with_capacity(observations.len());
    for o in observations.iter() {
        let own = &seen[o.sequence.as_str()];
        let carriers = own
            .values()
            .filter(|(_, dosage)| dosage.is_some_and(|d| d >= COPY_LIKE))
            .count();
        let (mut tracking, mut recurrent) = (None, None);
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
            }
        }
        context.push((carriers, tracking, recurrent));
    }
    for (o, (carriers, tracking, recurrent)) in observations.iter_mut().zip(context) {
        o.evidence.carriers = carriers;
        o.evidence.tracking = tracking;
        o.evidence.recurrent_ratio = recurrent;
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
        set_context(&mut observations);
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
    }
}
