//! Banded approximate matching of a pattern inside a text.
//!
//! Used to locate a database allele in a read when a few differences are
//! allowed. The caller knows roughly where the pattern lies (a diagonal,
//! from an exact seed), so only a band of `2 * band + 1` cells per pattern
//! position is filled: `O(pattern length * band)` instead of the full
//! `O(pattern length * text length)` table.

/// Where a pattern lies in a text.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Placement {
    /// Edit (Levenshtein) distance between the pattern and `text[start..end]`.
    pub distance: usize,
    /// 0-based position of the first text byte of the match.
    pub start: usize,
    /// 0-based position one past the last text byte of the match.
    pub end: usize,
}

const UNREACHABLE: usize = usize::MAX / 2;

/// One cell of the table: distance, and where in the text the alignment
/// that achieves it starts.
#[derive(Clone, Copy)]
struct Cell {
    distance: usize,
    start: usize,
}

const OUTSIDE: Cell = Cell {
    distance: UNREACHABLE,
    start: 0,
};

/// Best placement of the whole `pattern` in `text` among alignments that
/// stay within `band` cells of `diagonal` (text position minus pattern
/// position), or `None` if every such alignment costs more than `cutoff`.
///
/// The text may start and end anywhere at no cost; the pattern must be
/// matched completely. Among placements of equal distance the one whose
/// length is closest to the pattern's wins (so a difference at the first or
/// last base is read as a substitution, not as a shorter match), then the
/// leftmost.
pub fn place_near_diagonal(
    pattern: &[u8],
    text: &[u8],
    diagonal: i64,
    band: usize,
    cutoff: usize,
) -> Option<Placement> {
    let width = 2 * band + 1;
    let text_len = text.len() as i64;
    // Text column of band slot `k` in the row of pattern prefix length `i`.
    let column = |i: usize, k: usize| i as i64 + diagonal - band as i64 + k as i64;

    let mut previous: Vec<Cell> = (0..width)
        .map(|k| {
            let j = column(0, k);
            if (0..=text_len).contains(&j) {
                Cell {
                    distance: 0,
                    start: j as usize,
                }
            } else {
                OUTSIDE
            }
        })
        .collect();
    let mut current = vec![OUTSIDE; width];

    for i in 1..=pattern.len() {
        let mut row_minimum = UNREACHABLE;
        for k in 0..width {
            let j = column(i, k);
            current[k] = OUTSIDE;
            if !(0..=text_len).contains(&j) {
                continue;
            }
            let mut best = OUTSIDE;
            // Pattern base against text base (same slot of the row above).
            if j >= 1 && previous[k].distance < UNREACHABLE {
                let mismatch = usize::from(pattern[i - 1] != text[j as usize - 1]);
                best = Cell {
                    distance: previous[k].distance + mismatch,
                    start: previous[k].start,
                };
            }
            // Pattern base missing from the text.
            if k + 1 < width && previous[k + 1].distance + 1 < best.distance {
                best = Cell {
                    distance: previous[k + 1].distance + 1,
                    start: previous[k + 1].start,
                };
            }
            // Extra base in the text.
            if k >= 1 && current[k - 1].distance + 1 < best.distance {
                best = Cell {
                    distance: current[k - 1].distance + 1,
                    start: current[k - 1].start,
                };
            }
            current[k] = best;
            row_minimum = row_minimum.min(best.distance);
        }
        if row_minimum > cutoff {
            return None;
        }
        std::mem::swap(&mut previous, &mut current);
    }

    let length_error = |cell: &Cell, end: usize| (end - cell.start).abs_diff(pattern.len());
    (0..width)
        .filter_map(|k| {
            let cell = previous[k];
            (cell.distance <= cutoff).then(|| (cell, column(pattern.len(), k) as usize))
        })
        .min_by_key(|(cell, end)| (cell.distance, length_error(cell, *end), *end))
        .map(|(cell, end)| Placement {
            distance: cell.distance,
            start: cell.start,
            end,
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::levenshtein;

    fn place(pattern: &str, text: &str, diagonal: i64, band: usize) -> Option<Placement> {
        place_near_diagonal(pattern.as_bytes(), text.as_bytes(), diagonal, band, band)
    }

    #[test]
    fn exact_occurrence_on_the_diagonal() {
        let p = place("GATTACA", "CCGATTACACC", 2, 0).unwrap();
        assert_eq!((p.distance, p.start, p.end), (0, 2, 9));
        assert_eq!(place("GATTACA", "CCGATTACACC", 3, 0), None);
        // Off the diagonal but inside the band.
        let p = place("GATTACA", "CCGATTACACC", 4, 2).unwrap();
        assert_eq!((p.distance, p.start, p.end), (0, 2, 9));
    }

    #[test]
    fn substitution_insertion_and_deletion() {
        let p = place("GATTACA", "CCGATGACACC", 2, 1).unwrap();
        assert_eq!((p.distance, p.start, p.end), (1, 2, 9));
        // One extra base in the text.
        let p = place("GATTACA", "CCGATTTACACC", 2, 1).unwrap();
        assert_eq!((p.distance, p.start, p.end), (1, 2, 10));
        // One pattern base missing from the text.
        let p = place("GATTACA", "CCGATACACC", 2, 1).unwrap();
        assert_eq!((p.distance, p.start, p.end), (1, 2, 8));
        assert_eq!(place("GATTACA", "CCGCTGACACC", 2, 1), None);
    }

    #[test]
    fn a_difference_at_either_end_is_a_substitution() {
        // First and last base differ: the match still spans all seven bases,
        // so the flanks outside it are not shifted.
        let p = place("GATTACA", "CCTATTACGCC", 2, 2).unwrap();
        assert_eq!((p.distance, p.start, p.end), (2, 2, 9));
    }

    #[test]
    fn pattern_hanging_over_the_text_ends() {
        // The read starts inside the allele: two pattern bases are missing.
        let p = place("GATTACA", "TTACACC", -2, 2).unwrap();
        assert_eq!((p.distance, p.start, p.end), (2, 0, 5));
        assert_eq!(place("GATTACA", "TTACACC", -2, 1), None);
    }

    #[test]
    fn distance_is_the_edit_distance_of_the_reported_stretch() {
        // Deterministic pseudo-random cases, cross-checked with the plain
        // Levenshtein implementation.
        let mut state = 12345u64;
        let mut next = move |bound: u64| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) % bound
        };
        for _ in 0..60 {
            let text: Vec<u8> = (0..60).map(|_| b"ACGT"[next(4) as usize]).collect();
            let start = next(30) as usize + 5;
            let mut pattern: Vec<u8> = text[start..start + 20].to_vec();
            let edits = next(4) as usize;
            for _ in 0..edits {
                let at = next(pattern.len() as u64) as usize;
                match next(3) {
                    0 => pattern[at] = b"ACGT"[next(4) as usize],
                    1 => pattern.insert(at, b"ACGT"[next(4) as usize]),
                    _ => {
                        pattern.remove(at);
                    }
                }
            }
            let band = 3;
            let found = place_near_diagonal(&pattern, &text, start as i64, band, band)
                .expect("at most three edits were made");
            assert!(found.distance <= edits);
            assert_eq!(
                found.distance,
                levenshtein(&pattern, &text[found.start..found.end])
            );
            // No stretch of the text is closer than what was reported
            // (unless it lies outside the band, which the edits made here
            // cannot cause to be better than `edits`).
            let closest_anywhere = (0..text.len())
                .flat_map(|s| (s..=text.len()).map(move |e| (s, e)))
                .map(|(s, e)| levenshtein(&pattern, &text[s..e]))
                .min()
                .unwrap();
            assert!(found.distance >= closest_anywhere);
        }
    }
}
