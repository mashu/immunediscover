//! Exact multi-pattern search.
//!
//! The original implementation calls `findfirst(pattern, read)` for every
//! (read, allele) pair. [`PatternSet`] returns the same answers — the
//! leftmost occurrence of each pattern — from a single Aho-Corasick pass
//! over the read, which is what makes large databases practical.

use aho_corasick::{AhoCorasick, MatchKind};

use crate::error::{Error, Result};

/// A leftmost occurrence of one pattern in a text.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Occurrence {
    /// Index of the pattern in the set.
    pub pattern: usize,
    /// 1-based position of the first matched character.
    pub start: usize,
    /// 1-based position of the last matched character.
    pub end: usize,
}

/// A fixed set of non-empty patterns searched simultaneously.
#[derive(Debug)]
pub struct PatternSet {
    automaton: AhoCorasick,
    n_patterns: usize,
}

impl PatternSet {
    pub fn new<I, P>(patterns: I) -> Result<Self>
    where
        I: IntoIterator<Item = P>,
        P: AsRef<[u8]>,
    {
        let patterns: Vec<P> = patterns.into_iter().collect();
        if let Some(i) = patterns.iter().position(|p| p.as_ref().is_empty()) {
            return Err(Error::invalid(format!(
                "Sequence number {} is empty and cannot be searched for",
                i + 1
            )));
        }
        let automaton = AhoCorasick::builder()
            .match_kind(MatchKind::Standard)
            .build(&patterns)
            .map_err(|e| Error::invalid(format!("Cannot build search index: {e}")))?;
        Ok(PatternSet {
            automaton,
            n_patterns: patterns.len(),
        })
    }

    pub fn len(&self) -> usize {
        self.n_patterns
    }

    pub fn is_empty(&self) -> bool {
        self.n_patterns == 0
    }

    /// Every occurrence of every pattern in `text`, overlapping ones
    /// included, in order of their end position.
    pub fn all_occurrences(&self, text: &str) -> Vec<Occurrence> {
        self.automaton
            .find_overlapping_iter(text)
            .map(|m| Occurrence {
                pattern: m.pattern().as_usize(),
                start: m.start() + 1,
                end: m.end(),
            })
            .collect()
    }

    /// Call `found(pattern, start, end)` for every occurrence of every
    /// pattern in `text`, overlapping ones included, in order of their end
    /// position; `start..end` counts from 0.
    pub fn each_occurrence(&self, text: &[u8], mut found: impl FnMut(usize, usize, usize)) {
        for m in self.automaton.find_overlapping_iter(text) {
            found(m.pattern().as_usize(), m.start(), m.end());
        }
    }

    /// Leftmost occurrence of every pattern present in `text`, ordered by
    /// pattern index (i.e. database order).
    pub fn first_occurrences(&self, text: &str) -> Vec<Occurrence> {
        let mut first: Vec<Option<(usize, usize)>> = vec![None; self.n_patterns];
        let mut hits = 0usize;
        for m in self.automaton.find_overlapping_iter(text) {
            let slot = &mut first[m.pattern().as_usize()];
            if slot.is_none_or(|(start, _)| m.start() < start) {
                hits += usize::from(slot.is_none());
                *slot = Some((m.start(), m.end()));
            }
        }
        let mut out = Vec::with_capacity(hits);
        for (pattern, slot) in first.into_iter().enumerate() {
            if let Some((start, end)) = slot {
                out.push(Occurrence {
                    pattern,
                    start: start + 1,
                    end,
                });
            }
        }
        out
    }
}

/// 1-based inclusive range of the first occurrence of `needle` in `haystack`.
pub fn find_first(needle: &str, haystack: &str) -> Option<(usize, usize)> {
    memchr::memmem::find(haystack.as_bytes(), needle.as_bytes())
        .map(|start| (start + 1, start + needle.len()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reports_leftmost_occurrence_per_pattern_in_pattern_order() {
        let set = PatternSet::new(["GATT", "ACA", "TTT", "A"]).unwrap();
        let hits = set.first_occurrences("CGATTACAGATTACA");
        assert_eq!(
            hits,
            vec![
                Occurrence {
                    pattern: 0,
                    start: 2,
                    end: 5
                },
                Occurrence {
                    pattern: 1,
                    start: 6,
                    end: 8
                },
                Occurrence {
                    pattern: 3,
                    start: 3,
                    end: 3
                },
            ]
        );
    }

    #[test]
    fn duplicate_and_nested_patterns_are_all_reported() {
        let set = PatternSet::new(["ACGT", "ACGT", "CG"]).unwrap();
        let hits = set.first_occurrences("TTACGTACGT");
        assert_eq!(hits.len(), 3);
        assert!(hits.iter().all(|h| h.start == 3 || h.start == 4));
    }

    #[test]
    fn agrees_with_naive_search() {
        let patterns = ["ACG", "CGT", "TTTT", "GA", "ACGTACGT"];
        let set = PatternSet::new(patterns).unwrap();
        let text = "GGACGTACGTTTGA";
        for hit in set.first_occurrences(text) {
            assert_eq!(
                find_first(patterns[hit.pattern], text),
                Some((hit.start, hit.end))
            );
        }
        assert!(find_first("TTTT", text).is_none());
    }

    #[test]
    fn empty_patterns_are_rejected() {
        assert!(PatternSet::new(["A", ""]).is_err());
    }
}
