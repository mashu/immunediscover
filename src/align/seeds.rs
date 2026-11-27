//! Finding the database sequences worth aligning to a read.
//!
//! Aligning a read to every sequence of a database would take far too
//! long. Sequences that resemble the read share short words with it, and
//! shared words also tell where in the read the sequence lies: a word at
//! read position `j` and subject position `i` votes for diagonal `j - i`.
//! The subjects are ranked by their votes, and the diagonal with most
//! votes is where the band of the alignment is placed.

/// Length of the words that are looked up.
pub const WORD: usize = 10;
const WORDS: usize = 1 << (2 * WORD);
const WORD_MASK: u32 = (WORDS - 1) as u32;
/// Diagonals are counted in groups of this many (a power of two).
const GROUP_BITS: u32 = 3;
/// Width of one group of diagonals.
pub const GROUP: usize = 1 << GROUP_BITS;

#[inline]
fn base_code(base: u8) -> Option<u32> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// The words of `sequence` as `(start, code)`; words holding anything but
/// A, C, G and T are left out.
pub fn words(sequence: &[u8]) -> impl Iterator<Item = (usize, u32)> + '_ {
    let mut code = 0u32;
    let mut valid = 0usize;
    sequence.iter().enumerate().filter_map(move |(end, &base)| {
        match base_code(base) {
            Some(c) => {
                code = ((code << 2) | c) & WORD_MASK;
                valid += 1;
            }
            None => valid = 0,
        }
        (valid >= WORD).then(|| (end + 1 - WORD, code))
    })
}

/// One place where a word occurs in the database.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Place {
    subject: u32,
    position: u32,
}

/// Every word of a set of sequences with the places where it occurs.
///
/// Sequences that are nearly the same (the alleles of a gene) are indexed
/// as one group: a word found at the same position in several of them is
/// noted once. A read then votes for the group, with at least as many
/// votes as any single member would get.
#[derive(Debug)]
pub struct WordIndex {
    /// `places[starts[w]..starts[w + 1]]` are the places of word `w`.
    starts: Vec<u32>,
    places: Vec<Place>,
    subjects: usize,
    longest: usize,
}

impl WordIndex {
    /// Index groups of sequences; the sequences of a group have the same
    /// length and are numbered by their group.
    pub fn new(groups: &[Vec<&[u8]>]) -> Self {
        // Per group: the distinct (word, position) pairs of its members.
        let distinct: Vec<Vec<(u32, u32)>> = groups
            .iter()
            .map(|members| {
                let mut found: Vec<(u32, u32)> = members
                    .iter()
                    .flat_map(|sequence| {
                        words(sequence).map(|(position, word)| (word, position as u32))
                    })
                    .collect();
                found.sort_unstable();
                found.dedup();
                found
            })
            .collect();
        let mut starts = vec![0u32; WORDS + 1];
        for &(word, _) in distinct.iter().flatten() {
            starts[word as usize + 1] += 1;
        }
        for word in 0..WORDS {
            starts[word + 1] += starts[word];
        }
        let mut next = starts.clone();
        let mut places = vec![
            Place {
                subject: 0,
                position: 0
            };
            starts[WORDS] as usize
        ];
        for (group, found) in distinct.iter().enumerate() {
            for &(word, position) in found {
                let slot = &mut next[word as usize];
                places[*slot as usize] = Place {
                    subject: group as u32,
                    position,
                };
                *slot += 1;
            }
        }
        let longest = groups
            .iter()
            .flatten()
            .map(|sequence| sequence.len())
            .max()
            .unwrap_or(0);
        WordIndex {
            starts,
            places,
            subjects: groups.len(),
            longest,
        }
    }

    #[inline]
    fn places(&self, word: u32) -> &[Place] {
        let word = word as usize;
        &self.places[self.starts[word] as usize..self.starts[word + 1] as usize]
    }

    /// Whether no word was indexed.
    pub fn is_empty(&self) -> bool {
        self.places.is_empty()
    }
}

/// A group of subjects that shares words with a read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Candidate {
    /// Number of the group in the index.
    pub subject: usize,
    /// Shared words, on any diagonal.
    pub votes: u32,
    /// Read position minus subject position where most words agree.
    pub diagonal: isize,
    /// The lowest and the highest diagonal on which words agree (two
    /// votes or more close together). They differ from `diagonal` when
    /// read and subject are shifted against each other part of the way:
    /// an inserted or deleted stretch.
    pub spread: (isize, isize),
}

/// Counters reused between reads.
#[derive(Debug, Default)]
pub struct Votes {
    /// Votes per subject and group of diagonals.
    counts: Vec<u16>,
    /// Cells of `counts` that are not zero.
    touched: Vec<u32>,
    /// Per subject: all votes, the best two neighbouring groups, and the
    /// lowest and highest group with two votes or more.
    totals: Vec<u32>,
    strongest: Vec<(u32, u32)>,
    spread: Vec<(u32, u32)>,
    voted: Vec<u32>,
}

impl Votes {
    /// Count the words `read` shares with every indexed subject, looking
    /// up every `step`-th word of the read, and add the subjects found to
    /// `candidates`.
    pub fn collect(
        &mut self,
        index: &WordIndex,
        read: &[u8],
        step: usize,
        candidates: &mut Vec<Candidate>,
    ) {
        if index.is_empty() || read.len() < WORD {
            return;
        }
        // Diagonals run from -longest to the read length; shifted to start at 0.
        let shift = index.longest;
        let groups = ((read.len() + shift) >> GROUP_BITS) + 2;
        let cells = index.subjects * groups;
        if self.counts.len() < cells {
            self.counts.resize(cells, 0);
        }
        if self.totals.len() < index.subjects {
            self.totals.resize(index.subjects, 0);
            self.strongest.resize(index.subjects, (0, 0));
            self.spread.resize(index.subjects, (u32::MAX, 0));
        }

        for (position, word) in words(read).filter(|(position, _)| position % step == 0) {
            for place in index.places(word) {
                let diagonal = position + shift - place.position as usize;
                let cell = place.subject as usize * groups + (diagonal >> GROUP_BITS);
                if self.counts[cell] == 0 {
                    self.touched.push(cell as u32);
                }
                self.counts[cell] = self.counts[cell].saturating_add(1);
            }
        }

        // A diagonal near the border of its group shares its votes with
        // the neighbour, so two neighbouring groups are judged together.
        for &cell in &self.touched {
            let cell = cell as usize;
            let subject = cell / groups;
            let group = cell % groups;
            let here = u32::from(self.counts[cell]);
            let next = if group + 1 < groups {
                u32::from(self.counts[cell + 1])
            } else {
                0
            };
            if self.totals[subject] == 0 {
                self.voted.push(subject as u32);
                self.strongest[subject] = (0, 0);
                self.spread[subject] = (u32::MAX, 0);
            }
            if here >= 2 {
                let (low, high) = self.spread[subject];
                self.spread[subject] = (low.min(group as u32), high.max(group as u32));
            }
            self.totals[subject] += here;
            let (best, best_group) = self.strongest[subject];
            let pair = here + next;
            if pair > best || (pair == best && (group as u32) < best_group) {
                self.strongest[subject] = (pair, group as u32);
            }
        }
        for &subject in &self.voted {
            let subject = subject as usize;
            let (pair, group) = self.strongest[subject];
            let group = group as usize;
            let here = u64::from(self.counts[subject * groups + group]);
            let next = u64::from(pair) - here;
            // Vote-weighted middle of the two groups.
            let middle = |g: usize| ((g << GROUP_BITS) + GROUP / 2) as u64;
            let centre = (middle(group) * here + middle(group + 1) * next) / (here + next);
            let diagonal = centre as isize - shift as isize;
            let (low, high) = self.spread[subject];
            let spread = if low == u32::MAX {
                (diagonal, diagonal)
            } else {
                let first = ((low as usize) << GROUP_BITS) as isize - shift as isize;
                let last = (((high as usize + 1) << GROUP_BITS) - 1) as isize - shift as isize;
                (first.min(diagonal), last.max(diagonal))
            };
            candidates.push(Candidate {
                subject,
                votes: self.totals[subject],
                diagonal,
                spread,
            });
            self.totals[subject] = 0;
        }
        for &cell in &self.touched {
            self.counts[cell as usize] = 0;
        }
        self.touched.clear();
        self.voted.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn index(sequences: &[&[u8]]) -> WordIndex {
        let groups: Vec<Vec<&[u8]>> = sequences.iter().map(|s| vec![*s]).collect();
        WordIndex::new(&groups)
    }

    #[test]
    fn words_skip_ambiguous_bases() {
        let found: Vec<usize> = words(b"ACGTACGTACGNACGTACGTACG")
            .map(|(start, _)| start)
            .collect();
        // One word before the N, three after it.
        assert_eq!(found, vec![0, 1, 12, 13]);
        let (_, first) = words(b"AAAAAAAAAC").next().unwrap();
        assert_eq!(first, 1);
    }

    #[test]
    fn the_subject_inside_the_read_gets_the_votes_and_the_diagonal() {
        let gene: &[u8] = b"ACGGTCATTGCAGGCTTAACGTTAGCCATGCAAGTCCGATAGGCTTAA";
        let other: &[u8] = b"TTGACCGTAAGGCTAGCTAGGATCCGATCGGATTCGAGCTAGGCTAAC";
        let index = index(&[other, gene]);
        let mut read = b"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG".to_vec();
        read.extend_from_slice(gene);
        read.extend_from_slice(b"CCCCCCCCCCCC");
        let mut candidates = Vec::new();
        Votes::default().collect(&index, &read, 1, &mut candidates);
        assert_eq!(candidates.len(), 1);
        let found = candidates[0];
        assert_eq!(found.subject, 1);
        assert_eq!(found.votes as usize, gene.len() - WORD + 1);
        assert!((found.diagonal - 31).unsigned_abs() <= GROUP);
    }

    #[test]
    fn a_shifted_part_of_the_read_widens_the_spread() {
        let gene: &[u8] = b"ACGGTCATTGCAGGCTTAACGTTAGCCATGCAAGTCCGATAGGCTTAATTGACCGTAAGGCTAGCTAGGATCCGATCGGATTCGAGCTAGGC";
        let index = index(&[gene]);
        // The read holds 30 bases that the gene lacks, in its middle.
        let mut read = gene[..46].to_vec();
        read.extend_from_slice(b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
        read.extend_from_slice(&gene[46..]);
        let mut candidates = Vec::new();
        Votes::default().collect(&index, &read, 1, &mut candidates);
        let (first, last) = candidates[0].spread;
        assert!(first <= 0 && last >= 30, "{first}..{last}");
        // Without the shift the spread stays within one group of diagonals.
        candidates.clear();
        Votes::default().collect(&index, gene, 1, &mut candidates);
        let (first, last) = candidates[0].spread;
        assert!(first <= 0 && last >= 0 && last - first < 2 * GROUP as isize);
    }

    #[test]
    fn counters_are_clean_for_the_next_read() {
        let gene: &[u8] = b"ACGGTCATTGCAGGCTTAACGTTAGCCATGCAAGTCCGATAGGCTTAA";
        let index = index(&[gene]);
        let mut votes = Votes::default();
        let mut first = Vec::new();
        votes.collect(&index, gene, 2, &mut first);
        let mut second = Vec::new();
        votes.collect(&index, gene, 2, &mut second);
        assert_eq!(first, second);
        assert_eq!(first[0].votes as usize, (gene.len() - WORD + 1).div_ceil(2));
        let mut none = Vec::new();
        votes.collect(&index, b"TTTTTTTTTTTTTTTTTTTTTTTT", 1, &mut none);
        assert!(none.is_empty());
    }

    #[test]
    fn a_group_gets_the_votes_of_its_best_member_or_more() {
        let gene: &[u8] = b"ACGGTCATTGCAGGCTTAACGTTAGCCATGCAAGTCCGATAGGCTTAA";
        let mut allele = gene.to_vec();
        allele[20] = b'T';
        let alone = index(&[gene]);
        let together = WordIndex::new(&[vec![gene, allele.as_slice()]]);
        let mut votes = Votes::default();
        for read in [gene, allele.as_slice()] {
            let (mut single, mut grouped) = (Vec::new(), Vec::new());
            votes.collect(&alone, read, 1, &mut single);
            votes.collect(&together, read, 1, &mut grouped);
            assert_eq!(grouped.len(), 1);
            // Either member, read whole, shares all its words with the group.
            assert_eq!(grouped[0].votes as usize, gene.len() - WORD + 1);
            assert!(grouped[0].votes >= single[0].votes);
        }
    }
}
