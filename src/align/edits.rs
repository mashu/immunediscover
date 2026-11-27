//! Quick answers about how far a read is from a subject.
//!
//! A full alignment ([`super::gotoh`]) is the expensive step of a search.
//! Most subjects can be judged without it: a subject that lies in the
//! read with at most one mismatch needs no alignment at all, and a subject
//! that cannot be turned into any part of the read with a given number of
//! changes cannot be the best match. The two functions here answer these
//! questions in a small fraction of the time of an alignment.

use super::gotoh::is_same_base;

/// Number of positions at which `subject` and an equally long piece of
/// the read differ, or `None` if there are more than `limit`. Bases are
/// compared as in an alignment: ambiguous bases never match.
pub fn mismatches(subject: &[u8], read: &[u8], limit: usize) -> Option<usize> {
    debug_assert_eq!(subject.len(), read.len());
    let mut total = 0;
    for (a, b) in subject.chunks(32).zip(read.chunks(32)) {
        total += a
            .iter()
            .zip(b)
            .filter(|&(&x, &y)| !is_same_base(x, y))
            .count();
        if total > limit {
            return None;
        }
    }
    Some(total)
}

/// Memory reused between calls of [`differences`].
#[derive(Debug, Default)]
pub struct Reach {
    now: Vec<i32>,
    next: Vec<i32>,
}

/// A diagonal that has not been reached.
const NOWHERE: i32 = i32::MIN / 2;

/// The least number of differences (mismatches, inserted and deleted
/// bases) with which `subject` can be aligned inside `read`, starting on
/// one of the diagonals `first..=last` (read position minus subject
/// position). `None` if more than `limit` are needed.
///
/// The whole subject is aligned, except that up to `hang` of its bases
/// may be left out where it reaches over the start or the end of the
/// read.
///
/// The method follows every diagonal as far as the two sequences agree
/// and counts how many differences it takes to get further (Landau and
/// Vishkin). Its cost grows with the number of differences, not with the
/// length of the sequences. Bases are compared as written, so the result
/// is never larger than the differences of any alignment.
pub fn differences(
    subject: &[u8],
    read: &[u8],
    first: isize,
    last: isize,
    limit: usize,
    hang: usize,
    reach: &mut Reach,
) -> Option<usize> {
    let m = subject.len() as i32;
    let n = read.len() as isize;
    if m == 0 || first > last {
        return None;
    }
    // How far subject and read agree from subject position `i` on
    // diagonal `d`.
    let slide = |mut i: i32, d: isize| -> i32 {
        let mut j = i as isize + d;
        while i < m && j < n && subject[i as usize] == read[j as usize] {
            i += 1;
            j += 1;
        }
        i
    };
    // Done when the subject is used up, or when the read is and no more
    // of the subject is left than may hang over its end.
    let hang = hang.min(subject.len().saturating_sub(1)) as i32;
    let done = |i: i32, d: isize| i == m || (i as isize + d == n && m - i <= hang);
    let limit = limit.min(subject.len() + read.len()) as isize;
    // Diagonals `first - limit ..= last + limit`, with one unreachable
    // neighbour on either side.
    let offset = first - limit - 1;
    let width = (last - first + 2 * limit + 3) as usize;
    reach.now.clear();
    reach.now.resize(width, NOWHERE);
    reach.next.clear();
    reach.next.resize(width, NOWHERE);
    let slot = |d: isize| (d - offset) as usize;

    for d in first.max(-(hang as isize))..=last.min(n) {
        // On a diagonal below 0 the first `-d` subject bases hang over
        // the start of the read.
        let i = slide((-d).max(0) as i32, d);
        if done(i, d) {
            return Some(0);
        }
        reach.now[slot(d)] = i;
    }
    for edits in 1..=limit {
        for d in first - edits..=last + edits {
            let at = slot(d);
            let same = reach.now[at];
            let left = reach.now[at - 1];
            let right = reach.now[at + 1];
            let mut i = NOWHERE;
            // A mismatch: one base of each.
            if same >= 0 && (same as isize + d) < n {
                i = same + 1;
            }
            // A base only in the read: from the diagonal below.
            if left >= 0 && (left as isize + d - 1) < n {
                i = i.max(left);
            }
            // A base only in the subject: from the diagonal above.
            if right >= 0 {
                i = i.max(right + 1);
            }
            if i < 0 || i as isize + d < 0 {
                reach.next[at] = NOWHERE;
                continue;
            }
            let i = slide(i.min(m), d);
            if done(i, d) {
                return Some(edits as usize);
            }
            reach.next[at] = i;
        }
        std::mem::swap(&mut reach.now, &mut reach.next);
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::gotoh::{align, Band, Ends, Scoring, Workspace};

    struct Dice(u64);

    impl Dice {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }

        fn below(&mut self, limit: usize) -> usize {
            (self.next() % limit as u64) as usize
        }

        fn sequence(&mut self, length: usize) -> Vec<u8> {
            (0..length).map(|_| b"ACGT"[self.below(4)]).collect()
        }
    }

    /// Fewest differences by the plain dynamic programme: read free at
    /// both ends, start anywhere, and up to `hang` subject bases left out
    /// over either end of the read.
    fn reference(subject: &[u8], read: &[u8], hang: usize) -> usize {
        let (m, n) = (subject.len(), read.len());
        let mut row: Vec<usize> = vec![0; n + 1];
        let mut least = usize::MAX;
        for (i, &a) in subject.iter().enumerate() {
            // Before the read: `hang` bases are left out for free, the
            // others are deleted.
            let start = (i + 1).saturating_sub(hang);
            let mut next = vec![start; n + 1];
            for (j, &b) in read.iter().enumerate() {
                next[j + 1] = (row[j] + usize::from(a != b))
                    .min(row[j + 1] + 1)
                    .min(next[j] + 1);
            }
            row = next;
            // The read is used up and the rest of the subject may hang.
            if m - (i + 1) <= hang.min(m - 1) && i + 1 < m {
                least = least.min(row[n]);
            }
        }
        least.min(*row.iter().min().unwrap())
    }

    #[test]
    fn mismatches_are_counted_up_to_the_limit() {
        assert_eq!(mismatches(b"GATTACA", b"GATTACA", 0), Some(0));
        assert_eq!(mismatches(b"GATTACA", b"GACTACA", 1), Some(1));
        assert_eq!(mismatches(b"GATTACA", b"GACTACC", 1), None);
        // An ambiguous base is a mismatch even against itself.
        assert_eq!(mismatches(b"GANTACA", b"GANTACA", 3), Some(1));
        let long = vec![b'A'; 100];
        let mut other = long.clone();
        other[99] = b'C';
        assert_eq!(mismatches(&long, &other, 1), Some(1));
    }

    #[test]
    fn differences_agree_with_the_plain_dynamic_programme() {
        let mut dice = Dice(0x0BAD_C0DE_1234_5678);
        let mut reach = Reach::default();
        for round in 0..400 {
            let length = 20 + dice.below(60);
            let subject = dice.sequence(length);
            let mut core = subject.clone();
            for _ in 0..dice.below(7) {
                let at = dice.below(core.len());
                match dice.below(3) {
                    0 => core[at] = b"ACGT"[dice.below(4)],
                    1 => core.insert(at, b"ACGT"[dice.below(4)]),
                    _ => {
                        core.remove(at);
                    }
                }
            }
            let before = dice.below(25);
            let mut read = dice.sequence(before);
            read.extend(core);
            let after = dice.below(25);
            read.extend(dice.sequence(after));
            // Some reads are cut at either end, inside the subject.
            if round % 3 == 0 {
                let cut = before + dice.below(12);
                read.drain(..cut.min(read.len() - 1));
            }
            if round % 4 == 0 {
                let keep = read.len() - dice.below(read.len().min(30));
                read.truncate(keep.max(1));
            }
            // Every start diagonal is allowed: the answer is the true one.
            let all = (-(subject.len() as isize), read.len() as isize);
            for hang in [0, 3, 15] {
                let expected = reference(&subject, &read, hang);
                for limit in [0, 2, 5, 40] {
                    let found = differences(&subject, &read, all.0, all.1, limit, hang, &mut reach);
                    let wanted = (expected <= limit).then_some(expected);
                    assert_eq!(found, wanted, "round {round}, limit {limit}, hang {hang}");
                }
                // With few start diagonals the answer can only be larger.
                let near = before as isize;
                let few = differences(&subject, &read, near - 2, near + 2, 40, hang, &mut reach);
                assert!(few.is_none_or(|found| found >= expected));
            }
        }
    }

    #[test]
    fn no_alignment_has_fewer_differences() {
        let mut dice = Dice(0x7777_1111_3333_9999);
        let mut reach = Reach::default();
        let mut workspace = Workspace::default();
        let scoring = Scoring::default();
        for _ in 0..200 {
            let subject = dice.sequence(60);
            let mut core = subject.clone();
            for _ in 0..dice.below(5) {
                let at = dice.below(core.len());
                core[at] = b"ACGT"[dice.below(4)];
            }
            if dice.below(2) == 0 {
                core.remove(30);
            }
            let mut read = dice.sequence(20);
            read.extend(core);
            read.extend(dice.sequence(20));
            let band = Band {
                diagonal: 20,
                half_width: 6,
            };
            let ends = Ends::SubjectInsideRead;
            let path = align(&subject, &read, &scoring, band, ends, &mut workspace).unwrap();
            let least = differences(&subject, &read, 14, 26, 60, 0, &mut reach).unwrap();
            assert!(least <= path.differences());
        }
    }
}
