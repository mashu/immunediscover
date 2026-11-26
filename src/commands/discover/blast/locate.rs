//! Finding database alleles in reads, exactly or within a few differences.
//!
//! The context extension needs to know where an allele lies in a read. With
//! `max_distance == 0` that is an exact search. With a tolerance, a read in
//! which no allele occurs exactly is assigned to the allele(s) *nearest* to
//! it by edit distance, provided that distance does not exceed the
//! tolerance. This is what makes extension possible when the database comes
//! from another species or lacks the donor's alleles.
//!
//! Approximate search is seeded: an allele cut into `max_distance + 1`
//! pieces cannot be within `max_distance` edits of a read without one piece
//! occurring exactly, so the pieces of all alleles are searched in one pass
//! and only the places they point to are verified, with a banded alignment.

use std::collections::HashMap;

use crate::align::{place_near_diagonal, Placement};
use crate::error::Result;
use crate::seq::search::PatternSet;

/// A sequence found in a read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Located {
    /// Index of the sequence in the set given to [`Locator::new`].
    pub target: usize,
    /// 0-based position of the first read base of the match.
    pub start: usize,
    /// 0-based position one past the last read base of the match.
    pub end: usize,
    pub distance: usize,
}

/// Searches reads for a fixed set of sequences.
pub struct Locator<'a> {
    targets: &'a [String],
    exact: PatternSet,
    approximate: Option<Seeds>,
}

/// Exact pieces of the targets and what they belong to.
struct Seeds {
    pieces: PatternSet,
    /// For every piece: the targets containing it and its offset in each.
    owners: Vec<Vec<(usize, usize)>>,
    max_distance: usize,
}

impl Seeds {
    fn new(targets: &[String], max_distance: usize) -> Result<Self> {
        let mut index: HashMap<&[u8], usize> = HashMap::new();
        let mut pieces: Vec<&[u8]> = Vec::new();
        let mut owners: Vec<Vec<(usize, usize)>> = Vec::new();
        for (target, sequence) in targets.iter().enumerate() {
            let bytes = sequence.as_bytes();
            let n = max_distance + 1;
            for k in 0..n {
                let (from, to) = (k * bytes.len() / n, (k + 1) * bytes.len() / n);
                if from == to {
                    continue; // target shorter than the number of pieces
                }
                let piece = &bytes[from..to];
                let slot = *index.entry(piece).or_insert_with(|| {
                    pieces.push(piece);
                    owners.push(Vec::new());
                    pieces.len() - 1
                });
                owners[slot].push((target, from));
            }
        }
        Ok(Seeds {
            pieces: PatternSet::new(pieces)?,
            owners,
            max_distance,
        })
    }

    /// The targets nearest to `read`, if within the tolerance.
    fn nearest(&self, targets: &[String], read: &str) -> Vec<Located> {
        // Where a seed says a target would start in the read.
        let mut candidates: Vec<(usize, i64)> = self
            .pieces
            .all_occurrences(read)
            .into_iter()
            .flat_map(|hit| {
                self.owners[hit.pattern]
                    .iter()
                    .map(move |&(target, offset)| (target, hit.start as i64 - 1 - offset as i64))
            })
            .collect();
        candidates.sort_unstable();
        candidates.dedup();

        let band = self.max_distance;
        let mut cutoff = self.max_distance;
        let mut found: Vec<Located> = Vec::new();
        let mut verified: Option<(usize, i64)> = None;
        for (target, diagonal) in candidates {
            // A band around one seed of an alignment contains the whole
            // alignment, so neighbouring diagonals need no second look.
            if verified.is_some_and(|(t, d)| t == target && diagonal - d <= band as i64) {
                continue;
            }
            verified = Some((target, diagonal));
            let placement = place_near_diagonal(
                targets[target].as_bytes(),
                read.as_bytes(),
                diagonal,
                band,
                cutoff,
            );
            let Some(Placement {
                distance,
                start,
                end,
            }) = placement
            else {
                continue;
            };
            if distance < cutoff {
                cutoff = distance;
                found.clear();
            }
            if found.iter().all(|located| located.target != target) {
                found.push(Located {
                    target,
                    start,
                    end,
                    distance,
                });
            }
        }
        found
    }
}

impl<'a> Locator<'a> {
    /// `max_distance` is the largest edit distance at which a read is still
    /// assigned to a target; 0 means exact matches only.
    pub fn new(targets: &'a [String], max_distance: usize) -> Result<Self> {
        Ok(Locator {
            targets,
            exact: PatternSet::new(targets.iter().map(|t| t.as_bytes()))?,
            approximate: if max_distance > 0 {
                Some(Seeds::new(targets, max_distance)?)
            } else {
                None
            },
        })
    }

    /// Every target occurring exactly in `read` (leftmost occurrence each,
    /// in target order). If there is none and a tolerance was given, the
    /// nearest target(s) within the tolerance instead.
    pub fn locate(&self, read: &str) -> Vec<Located> {
        let exact: Vec<Located> = self
            .exact
            .first_occurrences(read)
            .into_iter()
            .map(|hit| Located {
                target: hit.pattern,
                start: hit.start - 1,
                end: hit.end,
                distance: 0,
            })
            .collect();
        match &self.approximate {
            Some(seeds) if exact.is_empty() => seeds.nearest(self.targets, read),
            _ => exact,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn targets() -> Vec<String> {
        vec![
            "GATTACAGATTACACC".to_string(), // A*01
            "GATTACAGATTAGACC".to_string(), // A*02, one base from A*01
            "TTGCATGCATGCAAGG".to_string(), // B*01
        ]
    }

    fn hits(max_distance: usize, read: &str) -> Vec<(usize, usize, usize, usize)> {
        let targets = targets();
        Locator::new(&targets, max_distance)
            .unwrap()
            .locate(read)
            .into_iter()
            .map(|l| (l.target, l.start, l.end, l.distance))
            .collect()
    }

    #[test]
    fn exact_matches_need_no_tolerance_and_take_precedence() {
        let read = "AAAGATTACAGATTAGACCTTT";
        assert_eq!(hits(0, read), [(1, 3, 19, 0)]);
        // A*01 is one base away, but the read has an exact allele.
        assert_eq!(hits(2, read), [(1, 3, 19, 0)]);
    }

    #[test]
    fn without_an_exact_match_the_nearest_allele_is_taken() {
        // A*01 with one substitution (position 5 of the allele): distance 1
        // to A*01, distance 2 to A*02.
        let read = "AAAGATTAGAGATTACACCTTT";
        assert_eq!(hits(0, read), []);
        assert_eq!(hits(1, read), [(0, 3, 19, 1)]);
        assert_eq!(hits(3, read), [(0, 3, 19, 1)]);
    }

    #[test]
    fn equally_near_alleles_are_all_reported() {
        // Differs from both A alleles at their distinguishing base.
        let read = "AAAGATTACAGATTATACCTTT";
        assert_eq!(hits(1, read), [(0, 3, 19, 1), (1, 3, 19, 1)]);
    }

    #[test]
    fn indels_move_the_boundaries_not_the_flanks() {
        // B*01 with one base deleted: the match is 15 bases long and the
        // context starts right after it.
        let read = "CCCCTTGCATGATGCAAGGAAAA";
        assert_eq!(hits(1, read), [(2, 4, 19, 1)]);
        // Too far from everything.
        assert_eq!(hits(2, "ACGTACGTACGTACGTACGTACGT"), []);
    }
}
