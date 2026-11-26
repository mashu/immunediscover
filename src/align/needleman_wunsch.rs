//! Affine-gap Needleman-Wunsch with position-dependent gap scores.
//!
//! This is a line-by-line port of BioAlignments.jl's `NeedlemanWunsch`
//! (`run!` and `traceback`). Reproducing its recurrences *and* its traceback
//! priorities is what guarantees that, among equally scoring alignments,
//! the same one is chosen as in the original implementation.

const TRACE_NONE: u8 = 0b00000;
const TRACE_MATCH: u8 = 0b00001;
const TRACE_DELETE: u8 = 0b00010;
const TRACE_INSERT: u8 = 0b00100;
const TRACE_EXTDEL: u8 = 0b01000;
const TRACE_EXTINS: u8 = 0b10000;

/// Gap scores (non-positive) for one sequence, by position of the gap.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GapScores {
    pub start_open: i32,
    pub start_extend: i32,
    pub middle_open: i32,
    pub middle_extend: i32,
    pub end_open: i32,
    pub end_extend: i32,
}

impl GapScores {
    /// The same affine score everywhere.
    pub fn uniform(open: i32, extend: i32) -> Self {
        GapScores {
            start_open: open,
            start_extend: extend,
            middle_open: open,
            middle_extend: extend,
            end_open: open,
            end_extend: extend,
        }
    }

    /// Affine score inside the sequence, free gaps before and after it.
    pub fn free_ends(open: i32, extend: i32) -> Self {
        GapScores {
            start_open: 0,
            start_extend: 0,
            middle_open: open,
            middle_extend: extend,
            end_open: 0,
            end_extend: 0,
        }
    }
}

/// One alignment column; `None` is a gap.
pub type AlignedPair = (Option<u8>, Option<u8>);

/// An optimal alignment of `a` against `b`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub score: i32,
    /// Columns from left to right, covering all of `a` and all of `b`.
    pub pairs: Vec<AlignedPair>,
}

struct TraceMatrix {
    cells: Vec<u8>,
    rows: usize,
}

impl TraceMatrix {
    fn new(m: usize, n: usize) -> Self {
        TraceMatrix {
            cells: vec![0xff; (m + 1) * (n + 1)],
            rows: m + 1,
        }
    }

    #[inline]
    fn at(&mut self, i: usize, j: usize) -> &mut u8 {
        &mut self.cells[j * self.rows + i]
    }

    #[inline]
    fn get(&self, i: usize, j: usize) -> u8 {
        self.cells[j * self.rows + i]
    }
}

/// Align `a` and `b` (sequences of substitution-matrix codes).
///
/// `gaps_a` scores gaps inserted into `a` (columns consuming only `b`),
/// `gaps_b` scores gaps inserted into `b`.
pub fn align(
    a: &[u8],
    b: &[u8],
    score: impl Fn(u8, u8) -> i32,
    gaps_a: GapScores,
    gaps_b: GapScores,
) -> Alignment {
    let m = a.len();
    let n = b.len();
    let mut trace = TraceMatrix::new(m, n);
    let score = fill(a, b, &score, gaps_a, gaps_b, &mut trace);
    Alignment {
        score,
        pairs: traceback(a, b, &trace),
    }
}

fn fill(
    a: &[u8],
    b: &[u8],
    score: &impl Fn(u8, u8) -> i32,
    ga: GapScores,
    gb: GapScores,
    trace: &mut TraceMatrix,
) -> i32 {
    let m = a.len();
    let n = b.len();
    let middle_init_a = ga.middle_open + ga.middle_extend;
    let end_init_a = ga.end_open + ga.end_extend;
    let middle_init_b = gb.middle_open + gb.middle_extend;
    let end_init_b = gb.end_open + gb.end_extend;

    // h[i] is H[i+1] of the original (row i of the current column).
    let mut h = vec![0i32; m + 1];
    let mut e = vec![0i32; m];
    *trace.at(0, 0) = TRACE_NONE;

    if m == 0 && n == 0 {
        return 0;
    }
    if m == 0 {
        for j in 1..=n {
            *trace.at(0, j) = TRACE_DELETE;
        }
        return ga.start_open + ga.start_extend * n as i32;
    }
    if n == 0 {
        for i in 1..=m {
            *trace.at(i, 0) = TRACE_INSERT;
        }
        return gb.start_open + gb.start_extend * m as i32;
    }

    for i in 1..=m {
        h[i] = gb.start_open + gb.start_extend * i as i32;
        e[i - 1] = h[i] + if i == m { end_init_a } else { middle_init_a };
        *trace.at(i, 0) = TRACE_INSERT;
        *trace.at(i, 1) = TRACE_NONE;
    }

    for j in 1..=n {
        let last_column = j == n;
        let b_j = b[j - 1];
        let mut h_diag = h[0];
        h[0] = ga.start_open + ga.start_extend * j as i32;
        let mut f = h[0]
            + if last_column {
                end_init_b
            } else {
                middle_init_b
            };
        let mut ft = TRACE_NONE;
        *trace.at(0, j) = TRACE_DELETE;

        for i in 1..=m {
            let last_row = i == m;
            // Cell update.
            let e_i = e[i - 1];
            let g = h_diag + score(a[i - 1], b_j);
            let best = e_i.max(f).max(g);
            h_diag = h[i];
            h[i] = best;
            let mut t = trace.get(i, j) | ft;
            if e_i == best {
                t |= TRACE_DELETE;
            }
            if f == best {
                t |= TRACE_INSERT;
            }
            if g == best {
                t |= TRACE_MATCH;
            }
            *trace.at(i, j) = t;

            // Prepare the gap-in-`a` score for the next column.
            if !last_column {
                let (init, extend) = if last_row {
                    (end_init_a, ga.end_extend)
                } else {
                    (middle_init_a, ga.middle_extend)
                };
                let extended = e_i + extend;
                let next = extended.max(best + init);
                e[i - 1] = next;
                *trace.at(i, j + 1) = if next == extended {
                    TRACE_EXTDEL
                } else {
                    TRACE_NONE
                };
            }
            // Prepare the gap-in-`b` score for the next row.
            if !last_row {
                let (init, extend) = if last_column {
                    (end_init_b, gb.end_extend)
                } else {
                    (middle_init_b, gb.middle_extend)
                };
                let extended = f + extend;
                f = extended.max(best + init);
                ft = if f == extended {
                    TRACE_EXTINS
                } else {
                    TRACE_NONE
                };
            }
        }
    }
    h[m]
}

fn traceback(a: &[u8], b: &[u8], trace: &TraceMatrix) -> Vec<AlignedPair> {
    let mut i = a.len();
    let mut j = b.len();
    let mut pairs = Vec::with_capacity(i + j);
    let mut extending_insert = false;
    let mut extending_delete = false;
    while i >= 1 || j >= 1 {
        let t = trace.get(i, j);
        if extending_insert {
            extending_insert = t & TRACE_EXTINS > 0;
            pairs.push((Some(a[i - 1]), None));
            i -= 1;
        } else if extending_delete {
            extending_delete = t & TRACE_EXTDEL > 0;
            pairs.push((None, Some(b[j - 1])));
            j -= 1;
        } else if t & TRACE_INSERT > 0 {
            extending_insert = t & TRACE_EXTINS > 0;
            pairs.push((Some(a[i - 1]), None));
            i -= 1;
        } else if t & TRACE_DELETE > 0 {
            extending_delete = t & TRACE_EXTDEL > 0;
            pairs.push((None, Some(b[j - 1])));
            j -= 1;
        } else if t & TRACE_MATCH > 0 {
            pairs.push((Some(a[i - 1]), Some(b[j - 1])));
            i -= 1;
            j -= 1;
        } else {
            unreachable!("alignment trace matrix is inconsistent");
        }
    }
    pairs.reverse();
    pairs
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple(a: u8, b: u8) -> i32 {
        if a == b {
            0
        } else {
            -6
        }
    }

    fn render(pairs: &[AlignedPair]) -> (String, String) {
        let side = |pick: fn(&AlignedPair) -> Option<u8>| -> String {
            pairs
                .iter()
                .map(|p| pick(p).map_or('-', char::from))
                .collect()
        };
        (side(|p| p.0), side(|p| p.1))
    }

    #[test]
    fn global_alignment_with_gap() {
        let out = align(
            b"ACGT",
            b"ACT",
            simple,
            GapScores::uniform(-5, -3),
            GapScores::uniform(-5, -3),
        );
        assert_eq!(out.score, -8);
        let (a, b) = render(&out.pairs);
        assert_eq!(a.replace('-', ""), "ACGT");
        assert_eq!(b.replace('-', ""), "ACT");
        assert_eq!(a.len(), 4);
    }

    #[test]
    fn free_end_gaps_place_short_sequence_inside_long_one() {
        let out = align(
            b"ACGT",
            b"AACGTTT",
            simple,
            GapScores::free_ends(-5, -3),
            GapScores::uniform(-5, -3),
        );
        assert_eq!(out.score, 0);
        assert_eq!(
            render(&out.pairs),
            ("-ACGT--".to_string(), "AACGTTT".to_string())
        );
    }

    #[test]
    fn empty_inputs() {
        let g = GapScores::uniform(-5, -3);
        assert!(align(b"", b"", simple, g, g).pairs.is_empty());
        assert_eq!(align(b"", b"AC", simple, g, g).pairs.len(), 2);
        assert_eq!(align(b"AC", b"", simple, g, g).score, -11);
    }
}
