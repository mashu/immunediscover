//! The sequence an observation could be an artifact of.
//!
//! Errors and chimeras are copies of something more abundant in the same
//! reaction. For every observation the more abundant sequences of the same
//! gene in the same donor are examined: the nearest one becomes its
//! *parent*, and a sequence that starts like one of them and ends like
//! another is flagged as a chimera.

use std::collections::HashMap;

use super::Observation;
use crate::align::levenshtein;

/// More abundant sequences examined per observation, most abundant first.
const MAX_PARENTS: usize = 8;
/// Distances above this are not told apart.
pub const DISTANCE_LIMIT: usize = 9;

/// Differences between two sequences: substitutions when the lengths agree,
/// edit distance otherwise; anything above [`DISTANCE_LIMIT`] is reported as
/// `DISTANCE_LIMIT + 1`.
pub fn distance(a: &str, b: &str) -> usize {
    let (a, b) = (a.as_bytes(), b.as_bytes());
    if a.len() == b.len() {
        let mut differences = 0;
        for (x, y) in a.iter().zip(b) {
            differences += usize::from(x != y);
            if differences > DISTANCE_LIMIT {
                break;
            }
        }
        return differences.min(DISTANCE_LIMIT + 1);
    }
    if a.len().abs_diff(b.len()) > DISTANCE_LIMIT {
        return DISTANCE_LIMIT + 1;
    }
    levenshtein(a, b).min(DISTANCE_LIMIT + 1)
}

/// Is `sequence` a prefix of `head` joined to a suffix of `tail`, without
/// being either of them? Only sequences of one length are compared.
pub fn is_chimera_of(sequence: &str, head: &str, tail: &str) -> bool {
    let (s, h, t) = (sequence.as_bytes(), head.as_bytes(), tail.as_bytes());
    if s.len() != h.len() || s.len() != t.len() || s == h || s == t {
        return false;
    }
    let prefix = s.iter().zip(h).take_while(|(x, y)| x == y).count();
    let suffix = s
        .iter()
        .rev()
        .zip(t.iter().rev())
        .take_while(|(x, y)| x == y)
        .count();
    prefix + suffix >= s.len()
}

/// Fill in `parent` and `chimera` for every observation.
pub fn assign(observations: &mut [Observation]) {
    let mut groups: HashMap<(usize, &str), Vec<usize>> = HashMap::new();
    for (i, o) in observations.iter().enumerate() {
        groups
            .entry((o.donor, o.gene.as_str()))
            .or_default()
            .push(i);
    }
    struct Found {
        observation: usize,
        parent: Option<(usize, usize)>,
        chimera: bool,
    }
    let mut found: Vec<Found> = Vec::new();
    for members in groups.values_mut() {
        members.sort_by(|&a, &b| {
            (observations[b].count.cmp(&observations[a].count))
                .then_with(|| observations[a].sequence.cmp(&observations[b].sequence))
        });
        for (rank, &i) in members.iter().enumerate() {
            let o = &observations[i];
            let candidates: Vec<usize> = members[..rank]
                .iter()
                .copied()
                .filter(|&j| observations[j].count > o.count)
                .take(MAX_PARENTS)
                .collect();
            // The nearest candidate; among equally near ones the most abundant.
            let mut parent: Option<(usize, usize)> = None;
            for &j in &candidates {
                let d = distance(&o.sequence, &observations[j].sequence);
                if parent.is_none_or(|(_, best)| d < best) {
                    parent = Some((j, d));
                }
            }
            let chimera = candidates.iter().any(|&j| {
                candidates.iter().any(|&k| {
                    j != k
                        && is_chimera_of(
                            &o.sequence,
                            &observations[j].sequence,
                            &observations[k].sequence,
                        )
                })
            });
            found.push(Found {
                observation: i,
                parent,
                chimera,
            });
        }
    }
    for f in found {
        observations[f.observation].evidence.parent = f.parent;
        observations[f.observation].evidence.chimera = f.chimera;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::testing::observation;

    #[test]
    fn distance_counts_substitutions_and_indels() {
        assert_eq!(distance("ACGTACGT", "ACGTACGT"), 0);
        assert_eq!(distance("ACGTACGT", "ACGAACGA"), 2);
        assert_eq!(distance("ACGTACGT", "ACGACGT"), 1);
        assert_eq!(
            distance(&"A".repeat(30), &"C".repeat(30)),
            DISTANCE_LIMIT + 1
        );
        assert_eq!(
            distance(&"A".repeat(30), &"A".repeat(5)),
            DISTANCE_LIMIT + 1
        );
    }

    #[test]
    fn chimera_needs_two_different_parents() {
        assert!(is_chimera_of("AAAATTTT", "AAAAAAAA", "TTTTTTTT"));
        assert!(!is_chimera_of("AAAATTTT", "AAAATTTT", "TTTTTTTT"));
        assert!(!is_chimera_of("AACATTTT", "AAAAAAAA", "TTTTTTTT"));
        assert!(!is_chimera_of("AAAATTT", "AAAAAAAA", "TTTTTTTT"));
    }

    #[test]
    fn parent_is_the_nearest_more_abundant_sequence_of_the_gene() {
        let mut observations = vec![
            observation(0, "A", "AAAAAAAA", 100),
            observation(0, "A", "TTTTAAAA", 80),
            observation(0, "A", "AAAAAAAT", 5), // one from the first
            observation(0, "A", "AAAATTTT", 3), // nearest to the third; not a chimera
            observation(0, "A", "AAAAAAAA", 0),
            observation(1, "A", "AAAAAAAT", 9), // alone in its donor
            observation(0, "B", "AAAAAAAT", 9), // alone in its gene
        ];
        observations[4].sequence = "TTTTAAAT".to_string(); // head of the second, tail of the third
        observations[4].count = 2;
        assign(&mut observations);
        let parent = |i: usize| observations[i].evidence.parent;
        assert_eq!(parent(0), None);
        assert_eq!(parent(1), Some((0, 4)));
        assert_eq!(parent(2), Some((0, 1)));
        assert_eq!(parent(3), Some((2, 3)));
        assert_eq!(parent(4), Some((1, 1)));
        assert!(observations[4].evidence.chimera);
        assert!(!observations[3].evidence.chimera && !observations[2].evidence.chimera);
        assert_eq!((parent(5), parent(6)), (None, None));
    }
}
