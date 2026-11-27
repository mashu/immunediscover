//! Banded alignment of a database sequence inside a read, with affine gaps.
//!
//! The *subject* (an allele, usually with the bases recorded around it) is
//! aligned from its first to its last base. The *read* is longer and only
//! the part that matches takes part: the bases before and after it cost
//! nothing. This is the semi-global alignment of Gotoh's algorithm, the
//! right model for a gene sitting somewhere inside a read.
//!
//! A read may also end in the middle of the subject. With
//! [`Ends::ReadMayEndInside`] the part of the subject that hangs over the
//! end of the read is left out for free, but only there: a subject is never
//! cut short while the read still goes on, so differences near the ends of
//! an allele are reported instead of being clipped away.
//!
//! Only the cells near one diagonal are computed. A read position `j`
//! paired with subject position `i` lies on diagonal `j - i`; the band
//! holds the diagonals within `half_width` of the expected one. A path
//! that reaches the border of the band says the band was too narrow, and
//! [`align_widening`] then repeats the alignment in a wider band.

/// Scores of an alignment. A gap of `n` bases costs
/// `gap_open + n * gap_extend`, as in BLAST.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Scoring {
    /// Score of two equal bases (positive).
    pub matched: i32,
    /// Score of two different bases (negative).
    pub mismatched: i32,
    /// Cost of starting a gap (zero or positive).
    pub gap_open: i32,
    /// Cost of every base of a gap (positive).
    pub gap_extend: i32,
}

impl Default for Scoring {
    /// The scores `blastn` uses by default.
    fn default() -> Self {
        Scoring {
            matched: 2,
            mismatched: -3,
            gap_open: 5,
            gap_extend: 2,
        }
    }
}

impl Scoring {
    /// Whether the values make sense for an alignment.
    pub fn is_valid(&self) -> bool {
        self.matched > 0 && self.mismatched < 0 && self.gap_open >= 0 && self.gap_extend > 0
    }

    #[inline]
    fn pair(&self, subject: u8, read: u8) -> i32 {
        if is_same_base(subject, read) {
            self.matched
        } else {
            self.mismatched
        }
    }
}

/// Two bases count as equal when they are the same unambiguous nucleotide.
/// `N` and the other ambiguity codes never match, not even themselves.
#[inline]
pub fn is_same_base(a: u8, b: u8) -> bool {
    a == b && matches!(a, b'A' | b'C' | b'G' | b'T')
}

/// What happens where the read ends.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ends {
    /// The whole subject must lie inside the read.
    SubjectInsideRead,
    /// The read may start or end inside the subject; the part of the
    /// subject beyond the read is left out.
    ReadMayEndInside,
}

/// The diagonals that are computed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Band {
    /// Expected read position minus subject position.
    pub diagonal: isize,
    /// Number of diagonals computed on either side of it.
    pub half_width: usize,
}

impl Band {
    /// A band that holds every cell of the alignment of a subject of
    /// `subject_length` with a read of `read_length`.
    pub fn whole(subject_length: usize, read_length: usize) -> Band {
        Band {
            diagonal: (read_length as isize - subject_length as isize) / 2,
            half_width: subject_length + read_length,
        }
    }

    /// Half width at which a band around `diagonal` holds every cell.
    fn covering(diagonal: isize, subject_length: usize, read_length: usize) -> usize {
        let below = diagonal + subject_length as isize;
        let above = read_length as isize - diagonal;
        below.max(above).max(0) as usize
    }
}

/// One column of an alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    /// The same base in subject and read.
    Match,
    /// Different bases.
    Mismatch,
    /// A read base without a partner in the subject.
    Insertion,
    /// A subject base without a partner in the read.
    Deletion,
}

impl Op {
    fn letter(self) -> char {
        match self {
            Op::Match => '=',
            Op::Mismatch => 'X',
            Op::Insertion => 'I',
            Op::Deletion => 'D',
        }
    }
}

/// The best alignment found in a band.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Path {
    pub score: i32,
    /// Aligned part of the subject, `subject_start..subject_end`.
    pub subject_start: usize,
    pub subject_end: usize,
    /// Aligned part of the read, `read_start..read_end`.
    pub read_start: usize,
    pub read_end: usize,
    /// Columns from left to right.
    pub ops: Vec<Op>,
    /// The path runs along the border of the band: a wider band might
    /// hold a better alignment.
    pub at_band_edge: bool,
}

impl Path {
    fn count(&self, op: Op) -> usize {
        self.ops.iter().filter(|&&o| o == op).count()
    }

    pub fn matches(&self) -> usize {
        self.count(Op::Match)
    }

    pub fn mismatches(&self) -> usize {
        self.count(Op::Mismatch)
    }

    /// Bases that lie in a gap, in the subject or in the read.
    pub fn gap_bases(&self) -> usize {
        self.count(Op::Insertion) + self.count(Op::Deletion)
    }

    /// Number of separate gaps.
    pub fn gap_openings(&self) -> usize {
        let is_gap = |op: Op| matches!(op, Op::Insertion | Op::Deletion);
        let first = usize::from(self.ops.first().is_some_and(|&op| is_gap(op)));
        first
            + self
                .ops
                .windows(2)
                .filter(|pair| is_gap(pair[1]) && pair[0] != pair[1])
                .count()
    }

    /// Mismatches plus gap bases.
    pub fn differences(&self) -> usize {
        self.mismatches() + self.gap_bases()
    }

    /// Number of alignment columns.
    pub fn columns(&self) -> usize {
        self.ops.len()
    }

    /// The columns in CIGAR notation with `=` for a match, `X` for a
    /// mismatch, `I` for a base only in the read and `D` for a base only
    /// in the subject, for example `120=1X30=2D144=`.
    pub fn cigar(&self) -> String {
        let mut text = String::new();
        let mut ops = self.ops.iter().copied().peekable();
        while let Some(op) = ops.next() {
            let mut run = 1usize;
            while ops.next_if_eq(&op).is_some() {
                run += 1;
            }
            text.push_str(&run.to_string());
            text.push(op.letter());
        }
        text
    }

    /// The aligned read bases, with `-` where the subject has a base the
    /// read lacks.
    pub fn gapped_read(&self, read: &[u8]) -> String {
        let mut position = self.read_start;
        let mut text = String::with_capacity(self.ops.len());
        for op in &self.ops {
            match op {
                Op::Deletion => text.push('-'),
                _ => {
                    text.push(read[position] as char);
                    position += 1;
                }
            }
        }
        text
    }
}

/// Memory reused between alignments.
#[derive(Debug, Default)]
pub struct Workspace {
    best: Vec<i32>,
    deletion: Vec<i32>,
    trace: Vec<u8>,
}

/// Score of a cell that cannot be reached; far from overflow.
const UNREACHABLE: i32 = i32::MIN / 4;

/// Where the best score of a cell comes from (two lowest trace bits).
const FROM_PAIR: u8 = 0;
const FROM_INSERTION: u8 = 1;
const FROM_DELETION: u8 = 2;
const FROM_START: u8 = 3;
const SOURCE: u8 = 0b11;
/// The insertion ending here continues one that was already open.
const INSERTION_EXTENDS: u8 = 0b0100;
/// The deletion ending here continues one that was already open.
const DELETION_EXTENDS: u8 = 0b1000;

/// Align `subject` inside `read` within `band`.
///
/// Both sequences are upper-case nucleotide text. `None` if no alignment
/// with a positive score exists in the band.
pub fn align(
    subject: &[u8],
    read: &[u8],
    scoring: &Scoring,
    band: Band,
    ends: Ends,
    workspace: &mut Workspace,
) -> Option<Path> {
    let m = subject.len();
    let n = read.len();
    if m == 0 || n == 0 {
        return None;
    }
    // Wider than every diagonal is the same as all of them.
    let half_width = band.half_width.min(m + n) as isize;
    let diagonal = band.diagonal.clamp(-(m as isize), n as isize);
    let lowest = diagonal - half_width;
    let width = (2 * half_width + 1) as usize;
    // Read position of band column `c` in row `i` (may be out of range).
    let read_position = |i: usize, c: usize| i as isize + lowest + c as isize;
    let overhang = ends == Ends::ReadMayEndInside;
    let open_and_extend = scoring.gap_open + scoring.gap_extend;

    workspace.best.clear();
    workspace.best.resize(width, UNREACHABLE);
    workspace.deletion.clear();
    workspace.deletion.resize(width, UNREACHABLE);
    workspace.trace.clear();
    workspace.trace.resize((m + 1) * width, FROM_START);
    let best = &mut workspace.best;
    let deletion = &mut workspace.deletion;
    let trace = &mut workspace.trace;

    // Row 0: the alignment may start at any read position.
    for (c, cell) in best.iter_mut().enumerate() {
        let j = read_position(0, c);
        if (0..=n as isize).contains(&j) {
            *cell = 0;
        }
    }

    // Best end where the read stops inside the subject: (score, row, column).
    let mut hanging: Option<(i32, usize, usize)> = None;
    for i in 1..=m {
        let base = subject[i - 1];
        let row = &mut trace[i * width..(i + 1) * width];
        let mut left = UNREACHABLE;
        let mut insertion = UNREACHABLE;
        for c in 0..width {
            let j = read_position(i, c);
            if j < 0 || j > n as isize {
                best[c] = UNREACHABLE;
                deletion[c] = UNREACHABLE;
                left = UNREACHABLE;
                insertion = UNREACHABLE;
                continue;
            }
            if j == 0 {
                // No read base yet: the subject starts before the read.
                let score = if overhang {
                    0
                } else {
                    -(scoring.gap_open + i as i32 * scoring.gap_extend)
                };
                best[c] = score;
                deletion[c] = UNREACHABLE;
                row[c] = FROM_START;
                left = score;
                insertion = UNREACHABLE;
                continue;
            }
            let j = j as usize;
            // `best[c]` still holds the row above: the cell on the diagonal.
            let pair = best[c] + scoring.pair(base, read[j - 1]);
            let (above, above_deletion) = if c + 1 < width {
                (best[c + 1], deletion[c + 1])
            } else {
                (UNREACHABLE, UNREACHABLE)
            };
            let mut code = FROM_PAIR;
            let open_deletion = above - open_and_extend;
            let extend_deletion = above_deletion - scoring.gap_extend;
            let delete = if extend_deletion > open_deletion {
                code |= DELETION_EXTENDS;
                extend_deletion
            } else {
                open_deletion
            };
            let open_insertion = left - open_and_extend;
            let extend_insertion = insertion - scoring.gap_extend;
            insertion = if extend_insertion > open_insertion {
                code |= INSERTION_EXTENDS;
                extend_insertion
            } else {
                open_insertion
            };
            // A pair is preferred to a gap, a deletion to an insertion.
            let mut score = pair;
            if delete > score {
                score = delete;
                code = (code & !SOURCE) | FROM_DELETION;
            }
            if insertion > score {
                score = insertion;
                code = (code & !SOURCE) | FROM_INSERTION;
            }
            let score = score.max(UNREACHABLE);
            best[c] = score;
            deletion[c] = delete.max(UNREACHABLE);
            row[c] = code;
            left = score;
        }
        if overhang && i < m {
            let c = n as isize - i as isize - lowest;
            if (0..width as isize).contains(&c) {
                let c = c as usize;
                // Later rows cover more of the subject: they win ties.
                if hanging.is_none_or(|(score, _, _)| best[c] >= score) {
                    hanging = Some((best[c], i, c));
                }
            }
        }
    }

    // Best end with the whole subject aligned; the leftmost wins ties.
    let mut end: Option<(i32, usize, usize)> = None;
    for (c, &score) in best.iter().enumerate() {
        let j = read_position(m, c);
        if (0..=n as isize).contains(&j) && end.is_none_or(|(top, _, _)| score > top) {
            end = Some((score, m, c));
        }
    }
    if let Some(candidate) = hanging {
        if end.is_none_or(|(top, _, _)| candidate.0 > top) {
            end = Some(candidate);
        }
    }
    let (score, end_row, end_column) = end?;
    if score <= 0 {
        return None;
    }

    // Walk back from the end.
    #[derive(PartialEq)]
    enum State {
        Best,
        Insertion,
        Deletion,
    }
    let mut ops: Vec<Op> = Vec::with_capacity(m + 8);
    let mut at_band_edge = false;
    let (mut i, mut c) = (end_row, end_column);
    let mut state = State::Best;
    let subject_end = end_row;
    let read_end = read_position(end_row, end_column) as usize;
    let (subject_start, read_start) = loop {
        let j = read_position(i, c) as usize;
        if (c == 0 && j > 0) || (c + 1 == width && j < n) {
            at_band_edge = true;
        }
        let code = trace[i * width + c];
        match state {
            State::Best => {
                if i == 0 {
                    break (0, j);
                }
                if j == 0 {
                    if overhang {
                        break (i, 0);
                    }
                    ops.extend(std::iter::repeat_n(Op::Deletion, i));
                    break (0, 0);
                }
                match code & SOURCE {
                    FROM_PAIR => {
                        ops.push(if is_same_base(subject[i - 1], read[j - 1]) {
                            Op::Match
                        } else {
                            Op::Mismatch
                        });
                        i -= 1;
                    }
                    FROM_DELETION => state = State::Deletion,
                    _ => state = State::Insertion,
                }
            }
            State::Insertion => {
                ops.push(Op::Insertion);
                if code & INSERTION_EXTENDS == 0 {
                    state = State::Best;
                }
                c -= 1;
            }
            State::Deletion => {
                ops.push(Op::Deletion);
                if code & DELETION_EXTENDS == 0 {
                    state = State::Best;
                }
                i -= 1;
                c += 1;
            }
        }
    };
    ops.reverse();
    Some(Path {
        score,
        subject_start,
        subject_end,
        read_start,
        read_end,
        ops,
        at_band_edge,
    })
}

/// Align in a band around `diagonal`, doubling its width for as long as
/// the path runs along its border, up to `max_half_width`.
#[allow(clippy::too_many_arguments)]
pub fn align_widening(
    subject: &[u8],
    read: &[u8],
    scoring: &Scoring,
    diagonal: isize,
    half_width: usize,
    max_half_width: usize,
    ends: Ends,
    workspace: &mut Workspace,
) -> Option<Path> {
    let covering = Band::covering(diagonal, subject.len(), read.len());
    let limit = max_half_width.max(half_width).min(covering);
    let mut half_width = half_width.max(1).min(limit.max(1));
    loop {
        let band = Band {
            diagonal,
            half_width,
        };
        let path = align(subject, read, scoring, band, ends, workspace);
        let cramped = path.as_ref().is_none_or(|p| p.at_band_edge);
        if !cramped || half_width >= limit {
            return path;
        }
        half_width = (half_width * 2).min(limit);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::needleman_wunsch::{self, GapScores};

    /// Small deterministic generator, so that failures can be repeated.
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

        fn base(&mut self) -> u8 {
            b"ACGT"[self.below(4)]
        }

        fn sequence(&mut self, length: usize) -> Vec<u8> {
            (0..length).map(|_| self.base()).collect()
        }

        /// `sequence` with substitutions, insertions and deletions.
        fn mutated(&mut self, sequence: &[u8], changes: usize) -> Vec<u8> {
            let mut text = sequence.to_vec();
            for _ in 0..changes {
                if text.is_empty() {
                    break;
                }
                let at = self.below(text.len());
                match self.below(3) {
                    0 => text[at] = self.base(),
                    1 => text.insert(at, self.base()),
                    _ => {
                        text.remove(at);
                    }
                }
            }
            text
        }
    }

    /// The best score according to the unbanded reference implementation.
    fn reference_score(subject: &[u8], read: &[u8], scoring: &Scoring, ends: Ends) -> i32 {
        let gaps = GapScores::uniform(-scoring.gap_open, -scoring.gap_extend);
        let free = GapScores::free_ends(-scoring.gap_open, -scoring.gap_extend);
        let in_read = match ends {
            Ends::SubjectInsideRead => gaps,
            Ends::ReadMayEndInside => free,
        };
        needleman_wunsch::align(subject, read, |a, b| scoring.pair(a, b), free, in_read).score
    }

    /// Score of the columns of `path`, counted independently.
    fn rescored(path: &Path, scoring: &Scoring) -> i32 {
        let pairs =
            path.matches() as i32 * scoring.matched + path.mismatches() as i32 * scoring.mismatched;
        pairs
            - path.gap_openings() as i32 * scoring.gap_open
            - path.gap_bases() as i32 * scoring.gap_extend
    }

    fn check_consistency(path: &Path, subject: &[u8], read: &[u8], scoring: &Scoring) {
        let on_subject = path.matches() + path.mismatches() + path.count(Op::Deletion);
        let on_read = path.matches() + path.mismatches() + path.count(Op::Insertion);
        assert_eq!(path.subject_end - path.subject_start, on_subject);
        assert_eq!(path.read_end - path.read_start, on_read);
        assert_eq!(rescored(path, scoring), path.score);
        // Every column says the truth about its two bases.
        let (mut i, mut j) = (path.subject_start, path.read_start);
        for op in &path.ops {
            match op {
                Op::Match => assert!(is_same_base(subject[i], read[j])),
                Op::Mismatch => assert!(!is_same_base(subject[i], read[j])),
                _ => {}
            }
            i += usize::from(*op != Op::Insertion);
            j += usize::from(*op != Op::Deletion);
        }
    }

    #[test]
    fn an_allele_inside_a_read_is_found_with_its_position() {
        let subject = b"GATTACAGATTACA";
        let read = b"CCCCCGATTACAGATTACATTTTT";
        let path = align(
            subject,
            read,
            &Scoring::default(),
            Band::whole(subject.len(), read.len()),
            Ends::SubjectInsideRead,
            &mut Workspace::default(),
        )
        .unwrap();
        assert_eq!((path.read_start, path.read_end), (5, 19));
        assert_eq!((path.subject_start, path.subject_end), (0, 14));
        assert_eq!(path.cigar(), "14=");
        assert_eq!(path.score, 28);
        assert_eq!(path.differences(), 0);
    }

    #[test]
    fn a_difference_at_the_end_of_the_allele_is_reported_not_clipped() {
        let subject = b"GATTACAGATTACC";
        let read = b"CCCCCGATTACAGATTACATTTTT";
        let path = align(
            subject,
            read,
            &Scoring::default(),
            Band::whole(subject.len(), read.len()),
            Ends::ReadMayEndInside,
            &mut Workspace::default(),
        )
        .unwrap();
        assert_eq!(path.cigar(), "13=1X");
        assert_eq!((path.subject_start, path.subject_end), (0, 14));
        assert_eq!(path.mismatches(), 1);
    }

    #[test]
    fn gaps_are_placed_and_counted() {
        let subject = b"ACGTACGTTTGGCCAAGGTTACGATCGA";
        //                         deletion of GG, insertion of T
        let read = b"TTTTACGTACGTTTCCAAGGTTACTGATCGATTTT";
        let path = align(
            subject,
            read,
            &Scoring::default(),
            Band::whole(subject.len(), read.len()),
            Ends::SubjectInsideRead,
            &mut Workspace::default(),
        )
        .unwrap();
        assert_eq!(path.cigar(), "10=2D10=1I6=");
        assert_eq!(path.gap_openings(), 2);
        assert_eq!(path.gap_bases(), 3);
        assert_eq!(path.gapped_read(read), "ACGTACGTTT--CCAAGGTTACTGATCGA");
    }

    #[test]
    fn a_read_that_ends_inside_the_allele_aligns_only_when_allowed() {
        let subject = b"ACGGTCATTGCAGGCTTAACGTTAGC";
        let read = b"TTTTTACGGTCATTGCAGG";
        let scoring = Scoring::default();
        let mut workspace = Workspace::default();
        let band = Band::whole(subject.len(), read.len());
        let hanging = align(
            subject,
            read,
            &scoring,
            band,
            Ends::ReadMayEndInside,
            &mut workspace,
        )
        .unwrap();
        assert_eq!((hanging.subject_start, hanging.subject_end), (0, 14));
        assert_eq!((hanging.read_start, hanging.read_end), (5, 19));
        assert_eq!(hanging.cigar(), "14=");

        // The other end: the read starts inside the subject.
        let read = b"TTAACGTTAGCAAAAA";
        let hanging = align(
            subject,
            read,
            &scoring,
            band,
            Ends::ReadMayEndInside,
            &mut workspace,
        )
        .unwrap();
        assert_eq!((hanging.subject_start, hanging.subject_end), (15, 26));
        assert_eq!((hanging.read_start, hanging.read_end), (0, 11));
    }

    #[test]
    fn ambiguous_bases_never_match() {
        let path = align(
            b"ACGTNACGT",
            b"ACGTNACGT",
            &Scoring::default(),
            Band::whole(9, 9),
            Ends::SubjectInsideRead,
            &mut Workspace::default(),
        )
        .unwrap();
        assert_eq!(path.cigar(), "4=1X4=");
    }

    #[test]
    fn scores_agree_with_the_unbanded_reference() {
        let mut dice = Dice(0x9E37_79B9_7F4A_7C15);
        let mut workspace = Workspace::default();
        let scorings = [
            Scoring::default(),
            Scoring {
                matched: 1,
                mismatched: -2,
                gap_open: 0,
                gap_extend: 2,
            },
            Scoring {
                matched: 5,
                mismatched: -4,
                gap_open: 10,
                gap_extend: 1,
            },
        ];
        for round in 0..600 {
            let scoring = scorings[round % scorings.len()];
            let subject_length = 1 + dice.below(60);
            let subject = dice.sequence(subject_length);
            // Related reads, unrelated reads and reads cut short.
            let read = match round % 4 {
                0 => {
                    let length = 1 + dice.below(80);
                    dice.sequence(length)
                }
                1 => {
                    let cut = dice.below(subject.len());
                    let mut read = dice.sequence(10);
                    read.extend_from_slice(&subject[..=cut]);
                    read
                }
                _ => {
                    let changes = dice.below(6);
                    let core = dice.mutated(&subject, changes);
                    let (before, after) = (dice.below(30), dice.below(30));
                    let mut read = dice.sequence(before);
                    read.extend_from_slice(&core);
                    read.extend(dice.sequence(after));
                    read
                }
            };
            if read.is_empty() {
                continue;
            }
            for ends in [Ends::SubjectInsideRead, Ends::ReadMayEndInside] {
                let expected = reference_score(&subject, &read, &scoring, ends);
                let band = Band::whole(subject.len(), read.len());
                match align(&subject, &read, &scoring, band, ends, &mut workspace) {
                    Some(path) => {
                        assert_eq!(path.score, expected, "round {round}, {ends:?}");
                        assert!(!path.at_band_edge);
                        check_consistency(&path, &subject, &read, &scoring);
                    }
                    None => assert!(expected <= 0, "round {round}, {ends:?}: {expected}"),
                }
            }
        }
    }

    #[test]
    fn a_narrow_band_on_the_right_diagonal_finds_the_same_alignment() {
        let mut dice = Dice(0x2545_F491_4F6C_DD1D);
        let mut workspace = Workspace::default();
        let scoring = Scoring::default();
        for round in 0..200 {
            let subject = dice.sequence(300);
            let changes = dice.below(8);
            let core = dice.mutated(&subject, changes);
            let before = dice.below(120);
            let mut read = dice.sequence(before);
            read.extend_from_slice(&core);
            read.extend(dice.sequence(100));
            let whole = align(
                &subject,
                &read,
                &scoring,
                Band::whole(subject.len(), read.len()),
                Ends::SubjectInsideRead,
                &mut workspace,
            )
            .unwrap();
            // The diagonal is known only roughly, as from shared words.
            let guess = before as isize + dice.below(9) as isize - 4;
            let banded = align_widening(
                &subject,
                &read,
                &scoring,
                guess,
                8,
                64,
                Ends::SubjectInsideRead,
                &mut workspace,
            )
            .unwrap();
            assert_eq!(banded.score, whole.score, "round {round}");
            check_consistency(&banded, &subject, &read, &scoring);
        }
    }

    #[test]
    fn a_band_on_the_wrong_diagonal_is_widened_until_it_fits() {
        let mut dice = Dice(7);
        let subject = dice.sequence(120);
        let mut read = dice.sequence(70);
        read.extend_from_slice(&subject);
        read.extend(dice.sequence(40));
        let scoring = Scoring::default();
        let mut workspace = Workspace::default();
        let narrow = align(
            &subject,
            &read,
            &scoring,
            Band {
                diagonal: 40,
                half_width: 4,
            },
            Ends::SubjectInsideRead,
            &mut workspace,
        );
        assert!(narrow.is_none_or(|p| p.at_band_edge || p.score < 240));
        let widened = align_widening(
            &subject,
            &read,
            &scoring,
            40,
            4,
            64,
            Ends::SubjectInsideRead,
            &mut workspace,
        )
        .unwrap();
        assert_eq!(widened.score, 240);
        assert_eq!((widened.read_start, widened.read_end), (70, 190));
        assert!(!widened.at_band_edge);
    }

    #[test]
    fn empty_and_unrelated_sequences_give_no_alignment() {
        let mut workspace = Workspace::default();
        let scoring = Scoring::default();
        let band = Band::whole(4, 4);
        let ends = Ends::SubjectInsideRead;
        assert!(align(b"", b"ACGT", &scoring, band, ends, &mut workspace).is_none());
        assert!(align(b"ACGT", b"", &scoring, band, ends, &mut workspace).is_none());
        assert!(align(b"AAAA", b"CCCC", &scoring, band, ends, &mut workspace).is_none());
    }
}
