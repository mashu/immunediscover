//! Finding where sequences lie in a genome.
//!
//! The question is the reverse of [`super::search`]: there, millions of
//! short reads are compared with a small database; here, a few sequences
//! (candidate alleles) are looked for in a genome of billions of bases.
//! The sequences are indexed, in both orientations, and the genome is read
//! through once. Where a stretch of the genome shares several words with a
//! sequence on one diagonal, the least number of differences between the
//! two is counted ([`super::edits`]). For every sequence the nearest place
//! on the wanted chromosome and the nearest place anywhere else are kept,
//! and aligned at the end ([`super::gotoh`]).
//!
//! Nothing of the genome is kept but the stretches of those places, so
//! the genome can come from a stream and need not fit in memory as a
//! whole: one record (chromosome) at a time is enough.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use rayon::prelude::*;

use super::edits::{differences, Reach};
use super::gotoh::{Ends, Op, Path, Workspace};
use super::search::{
    align_in, Database, Settings, HALF_WIDTH, MAX_COUNTED, MAX_HALF_WIDTH, SHORT_SUBJECT,
};
use super::seeds::{words, GROUP, WORD};
use crate::error::Result;
use crate::seq::search::PatternSet;

/// The genome is searched in pieces of at least this many bases.
const PIECE: usize = 4096;
/// Words on one place agree on diagonals at most this far apart.
const NEIGHBOURS: i32 = 2 * GROUP as i32;
/// Shared words that make a place worth counting, for sequences shorter
/// than `LONG` and for the others. Only words that do not overlap are
/// counted: a stretch of a dozen equal bases holds three words and is
/// still chance, in a genome of billions of bases.
const MIN_WORDS: usize = 2;
const MIN_WORDS_LONG: usize = 3;
const LONG: usize = 100;

/// Where a sequence lies in a genome.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Place {
    /// Number of the genome record, as given to [`GenomeSearch::search`].
    pub record: usize,
    /// First base of the place in the record, counted from 0.
    pub start: usize,
    /// The reverse complement of the sequence lies on the record.
    pub minus_strand: bool,
    /// Fewest mismatches, inserted and deleted bases between the sequence
    /// and the genome at this place.
    pub differences: usize,
    /// The alignment along the record: `=` equal, `X` different, `I` base
    /// only in the sequence, `D` base only in the genome.
    pub cigar: String,
    /// The bases of the genome at this place, written like the sequence.
    pub reference: String,
}

/// The nearest places of one sequence.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Places {
    /// On a record that counts as the wanted chromosome.
    pub on_target: Option<Place>,
    /// On any other record.
    pub elsewhere: Option<Place>,
}

/// A place that was found, with what an alignment of it needs.
#[derive(Debug)]
struct Found {
    sequence: usize,
    minus_strand: bool,
    on_target: bool,
    differences: usize,
    record: usize,
    /// Start of `stretch` in the record.
    start: usize,
    /// Diagonals of `stretch` on which the sequence may start.
    band: (isize, isize),
    stretch: Vec<u8>,
}

impl Found {
    /// Order in which places are preferred.
    fn rank(&self) -> (usize, usize, usize, bool) {
        (self.differences, self.record, self.start, self.minus_strand)
    }
}

/// Buffers of one thread.
#[derive(Default)]
struct Scratch {
    /// Shared words of one piece: (group, diagonal, position in the piece).
    hits: Vec<(u32, i32, u32)>,
    reach: Reach,
}

/// A set of sequences prepared for looking them up in genomes.
pub struct GenomeSearch {
    /// Sequence `i` as written is subject `2 * i`, its reverse complement
    /// subject `2 * i + 1`.
    database: Database,
    /// Subjects too short to be found by words, looked for unchanged.
    short: Option<(PatternSet, Vec<usize>)>,
    longest: usize,
    /// Fewest differences so far, per sequence: on target, elsewhere.
    fewest: [Vec<AtomicUsize>; 2],
    found: Mutex<Vec<Found>>,
    settings: Settings,
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

impl GenomeSearch {
    /// Prepare `sequences`. A place counts for a sequence when the two
    /// differ in at most `1 - min_identity` of the bases of the sequence.
    pub fn new(sequences: &[&str], min_identity: f64) -> Result<Self> {
        let mut oriented: Vec<(String, String)> = Vec::with_capacity(2 * sequences.len());
        for (number, sequence) in sequences.iter().enumerate() {
            let forward = sequence.trim().to_ascii_uppercase();
            let reverse = reverse_complement(forward.as_bytes());
            oriented.push((format!("{number}+"), forward));
            oriented.push((
                format!("{number}-"),
                String::from_utf8_lossy(&reverse).into_owned(),
            ));
        }
        let database = Database::new(oriented)?;
        let subjects = database.subjects();
        let short_subjects: Vec<usize> = (0..subjects.len())
            .filter(|&s| subjects[s].sequence.len() < SHORT_SUBJECT)
            .collect();
        let short = if short_subjects.is_empty() {
            None
        } else {
            let patterns: Vec<Vec<u8>> = short_subjects
                .iter()
                .map(|&s| subjects[s].sequence.clone())
                .collect();
            Some((PatternSet::new(patterns)?, short_subjects))
        };
        let longest = subjects.iter().map(|s| s.sequence.len()).max().unwrap_or(0);
        let accepted = |sequence: &&str| {
            let length = sequence.trim().len();
            AtomicUsize::new((((1.0 - min_identity) * length as f64) as usize).min(MAX_COUNTED))
        };
        let fewest = [
            sequences.iter().map(accepted).collect(),
            sequences.iter().map(accepted).collect(),
        ];
        let settings = Settings {
            min_identity,
            min_coverage: 1.0,
            ends: Ends::SubjectInsideRead,
            ..Settings::default()
        };
        Ok(GenomeSearch {
            database,
            short,
            longest,
            fewest,
            found: Mutex::new(Vec::new()),
            settings,
        })
    }

    /// Bases two neighbouring pieces share, so that every place lies whole
    /// in one of them.
    fn overlap(&self) -> usize {
        self.longest + MAX_COUNTED + 2 * HALF_WIDTH
    }

    /// Look for the sequences in one record of a genome. `bases` are
    /// upper-case; `on_target` tells whether the record counts as the
    /// wanted chromosome. Records can be searched in any order, with any
    /// number of threads: the result is the same.
    pub fn search(&self, record: usize, on_target: bool, bases: &[u8]) {
        let step = PIECE.max(4 * self.overlap());
        let starts: Vec<usize> = (0..bases.len()).step_by(step).collect();
        starts
            .par_iter()
            .for_each_init(Scratch::default, |scratch, &start| {
                let end = (start + step + self.overlap()).min(bases.len());
                self.search_piece(record, on_target, start, &bases[start..end], scratch);
            });
    }

    /// Note a place of subject `subject` with `count` differences, if no
    /// nearer one of its kind is known.
    fn note(
        &self,
        subject: usize,
        (record, on_target): (usize, bool),
        count: usize,
        start: usize,
        band: (isize, isize),
        stretch: &[u8],
    ) {
        let sequence = subject / 2;
        let fewest = &self.fewest[usize::from(!on_target)][sequence];
        if fewest.fetch_min(count, Ordering::Relaxed) < count {
            return;
        }
        let found = Found {
            sequence,
            minus_strand: subject % 2 == 1,
            on_target,
            differences: count,
            record,
            start,
            band,
            stretch: stretch.to_vec(),
        };
        self.found
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .push(found);
    }

    fn search_piece(
        &self,
        record: usize,
        on_target: bool,
        offset: usize,
        piece: &[u8],
        scratch: &mut Scratch,
    ) {
        // Short sequences: only where they lie unchanged.
        if let Some((patterns, subjects)) = &self.short {
            patterns.each_occurrence(piece, |pattern, start, end| {
                let place = (record, on_target);
                let stretch = &piece[start..end];
                self.note(subjects[pattern], place, 0, offset + start, (0, 0), stretch);
            });
        }

        // Words the piece shares with the groups of sequences.
        scratch.hits.clear();
        for (position, word) in words(piece) {
            for place in self.database.index.places(word) {
                let diagonal = position as i32 - place.position as i32;
                scratch
                    .hits
                    .push((place.subject, diagonal, position as u32));
            }
        }
        scratch.hits.sort_unstable();

        // Runs of words of one group on neighbouring diagonals.
        let hits = std::mem::take(&mut scratch.hits);
        let mut first = 0;
        while first < hits.len() {
            let (group, lowest, _) = hits[first];
            let mut highest = lowest;
            let mut next = first + 1;
            while next < hits.len() && hits[next].0 == group && hits[next].1 - highest <= NEIGHBOURS
            {
                highest = hits[next].1;
                next += 1;
            }
            let members = &self.database.groups[group as usize];
            let needed = if members.length >= LONG {
                MIN_WORDS_LONG
            } else {
                MIN_WORDS
            };
            if separate_words(&hits[first..next]) >= needed {
                let half = HALF_WIDTH as isize;
                let band = (lowest as isize - half, highest as isize + half);
                for &subject in &members.members {
                    self.count(subject, (record, on_target), offset, piece, band, scratch);
                }
            }
            first = next;
        }
        scratch.hits = hits;
    }

    /// Count the differences between subject `subject` and the piece on
    /// the diagonals `band`, and note the place if it is near enough.
    fn count(
        &self,
        subject: usize,
        place: (usize, bool),
        offset: usize,
        piece: &[u8],
        band: (isize, isize),
        scratch: &mut Scratch,
    ) {
        let sequence = &self.database.subject(subject).sequence;
        let limit = self.fewest[usize::from(!place.1)][subject / 2].load(Ordering::Relaxed);
        let reach = &mut scratch.reach;
        let Some(count) = differences(sequence, piece, band.0, band.1, limit, 0, reach) else {
            return;
        };
        // What an alignment will need: the piece from the lowest start to
        // the furthest end the sequence can have.
        let from = band.0.max(0) as usize;
        let to = (band.1.max(0) as usize + sequence.len() + count + 1).min(piece.len());
        let shifted = (band.0 - from as isize, band.1 - from as isize);
        self.note(
            subject,
            place,
            count,
            offset + from,
            shifted,
            &piece[from..to],
        );
    }

    /// Align the nearest places and return them, one entry per sequence
    /// in the order they were given.
    pub fn finish(mut self) -> Vec<Places> {
        let sequences = self.fewest[0].len();
        let mut nearest: Vec<[Option<Found>; 2]> = (0..sequences).map(|_| [None, None]).collect();
        let found = std::mem::take(
            self.found
                .get_mut()
                .unwrap_or_else(|poisoned| poisoned.into_inner()),
        );
        for place in found {
            let slot = &mut nearest[place.sequence][usize::from(!place.on_target)];
            if slot.as_ref().is_none_or(|held| place.rank() < held.rank()) {
                *slot = Some(place);
            }
        }
        let mut workspace = Workspace::default();
        nearest
            .into_iter()
            .map(|[on_target, elsewhere]| {
                let mut place = |found: Option<Found>| {
                    found.and_then(|found| self.place(found, &mut workspace))
                };
                Places {
                    on_target: place(on_target),
                    elsewhere: place(elsewhere),
                }
            })
            .collect()
    }

    /// Align one place.
    fn place(&self, found: Found, workspace: &mut Workspace) -> Option<Place> {
        let subject = 2 * found.sequence + usize::from(found.minus_strand);
        let sequence = &self.database.subject(subject).sequence;
        let path = if found.differences == 0 && found.stretch.len() == sequence.len() {
            Path {
                score: 0,
                subject_start: 0,
                subject_end: sequence.len(),
                read_start: 0,
                read_end: sequence.len(),
                ops: vec![Op::Match; sequence.len()],
                at_band_edge: false,
            }
        } else {
            // Wide enough for the gaps the differences may hold.
            let wide = MAX_HALF_WIDTH as isize;
            let band = (
                found.band.0.max(found.band.1 - 2 * wide),
                found.band.1.min(found.band.0 + 2 * wide),
            );
            align_in(band, sequence, &found.stretch, &self.settings, workspace)?
        };
        let on_genome = &found.stretch[path.read_start..path.read_end];
        let reference = if found.minus_strand {
            reverse_complement(on_genome)
        } else {
            on_genome.to_vec()
        };
        Some(Place {
            record: found.record,
            start: found.start + path.read_start,
            minus_strand: found.minus_strand,
            differences: found.differences,
            cigar: genome_cigar(&path),
            reference: String::from_utf8_lossy(&reference).into_owned(),
        })
    }
}

/// Number of words among `hits` (of one group, sorted) that do not
/// overlap: words that follow each other base by base on one diagonal
/// are one stretch of equal bases, and a stretch counts for as many words
/// as fit in it side by side.
fn separate_words(hits: &[(u32, i32, u32)]) -> usize {
    let mut words = 0;
    let mut first = 0;
    while first < hits.len() {
        let mut last = first;
        while last + 1 < hits.len()
            && hits[last + 1].1 == hits[first].1
            && hits[last + 1].2 == hits[last].2 + 1
        {
            last += 1;
        }
        // `last - first + 1` words starting one base apart cover that many
        // bases plus the rest of the last word.
        words += (last - first + WORD) / WORD;
        first = last + 1;
    }
    words
}

/// The alignment as a CIGAR along the genome, with the sequence as the
/// query: `I` is a base only in the sequence, `D` a base only in the
/// genome.
fn genome_cigar(path: &Path) -> String {
    let letter = |op: Op| match op {
        Op::Match => '=',
        Op::Mismatch => 'X',
        // The genome is the "read" of the alignment.
        Op::Insertion => 'D',
        Op::Deletion => 'I',
    };
    let mut text = String::new();
    let mut ops = path.ops.iter().copied().peekable();
    while let Some(op) = ops.next() {
        let mut run = 1usize;
        while ops.next_if_eq(&op).is_some() {
            run += 1;
        }
        text.push_str(&run.to_string());
        text.push(letter(op));
    }
    text
}

#[cfg(test)]
mod tests {
    use super::*;

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

        fn substituted(&mut self, sequence: &[u8], positions: &[usize]) -> Vec<u8> {
            let mut changed = sequence.to_vec();
            for &at in positions {
                changed[at] = match changed[at] {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    _ => b'A',
                };
            }
            changed
        }
    }

    fn text(bases: &[u8]) -> &str {
        std::str::from_utf8(bases).unwrap()
    }

    /// A record of random bases with `inserts` written at their positions.
    fn record(dice: &mut Dice, length: usize, inserts: &[(usize, &[u8])]) -> Vec<u8> {
        let mut bases = dice.sequence(length);
        for (at, insert) in inserts {
            bases[*at..*at + insert.len()].copy_from_slice(insert);
        }
        bases
    }

    #[test]
    fn overlapping_words_count_once() {
        // Three words one base apart on a diagonal: twelve equal bases.
        assert_eq!(separate_words(&[(0, 5, 10), (0, 5, 11), (0, 5, 12)]), 1);
        // Eleven such words: twenty equal bases, two words side by side.
        let run: Vec<(u32, i32, u32)> = (0..11).map(|i| (0, 5, 40 + i)).collect();
        assert_eq!(separate_words(&run), 2);
        // Two stretches on one diagonal, and one on the next.
        assert_eq!(
            separate_words(&[(0, 5, 10), (0, 5, 11), (0, 5, 60), (0, 6, 100)]),
            3
        );
        assert_eq!(separate_words(&[]), 0);
    }

    #[test]
    fn a_sequence_is_found_on_either_strand_with_its_position() {
        let mut dice = Dice(0xC0FF_EE00_1234_5678);
        let gene = dice.sequence(300);
        let turned = reverse_complement(&gene);
        // Far inside the record, and across the border of two pieces.
        let border = PIECE.max(4 * (300 + MAX_COUNTED + 2 * HALF_WIDTH)) - 100;
        let genome = record(&mut dice, 3 * border, &[(1000, &gene), (border, &turned)]);
        let search = GenomeSearch::new(&[text(&gene)], 0.8).unwrap();
        search.search(0, true, &genome[..border / 2]);
        search.search(1, false, &genome);
        let places = search.finish();
        let on = places[0].on_target.as_ref().unwrap();
        assert_eq!((on.record, on.start, on.minus_strand), (0, 1000, false));
        assert_eq!((on.differences, on.cigar.as_str()), (0, "300="));
        assert_eq!(on.reference.as_bytes(), gene);
        // Elsewhere the first of two equally near places is kept.
        let off = places[0].elsewhere.as_ref().unwrap();
        assert_eq!((off.record, off.start, off.minus_strand), (1, 1000, false));

        // Only the turned copy: found on the other strand, and the genome
        // there is written like the sequence.
        let search = GenomeSearch::new(&[text(&gene)], 0.8).unwrap();
        search.search(0, true, &genome[border - 500..]);
        let places = search.finish();
        let on = places[0].on_target.as_ref().unwrap();
        assert_eq!((on.start, on.minus_strand, on.differences), (500, true, 0));
        assert_eq!(on.reference.as_bytes(), gene);
        assert!(places[0].elsewhere.is_none());
    }

    #[test]
    fn the_nearest_place_is_kept_for_the_chromosome_and_for_elsewhere() {
        let mut dice = Dice(0xBEEF_0000_FACE_1111);
        let allele = dice.sequence(296);
        // The genome holds the gene with two other bases, a copy further
        // away on the same chromosome, and a distant relative elsewhere.
        let in_genome = dice.substituted(&allele, &[40, 200]);
        let copy = dice.substituted(&allele, &[10, 50, 90, 130, 170, 210, 250]);
        let mut relative = allele.clone();
        for at in (5..290).step_by(12) {
            relative = dice.substituted(&relative, &[at]);
        }
        let chromosome = record(&mut dice, 30_000, &[(4_000, &copy), (20_000, &in_genome)]);
        let other = record(&mut dice, 20_000, &[(7_000, &relative)]);
        let search = GenomeSearch::new(&[text(&allele)], 0.8).unwrap();
        search.search(0, false, &other);
        search.search(1, true, &chromosome);
        let places = search.finish();
        let on = places[0].on_target.as_ref().unwrap();
        assert_eq!((on.record, on.start, on.differences), (1, 20_000, 2));
        assert_eq!(on.cigar, "40=1X159=1X95=");
        assert_eq!(on.reference.as_bytes(), in_genome);
        let off = places[0].elsewhere.as_ref().unwrap();
        assert_eq!((off.record, off.start, off.differences), (0, 7_000, 24));
    }

    #[test]
    fn gaps_are_counted_and_written_along_the_genome() {
        let mut dice = Dice(0x1357_9BDF_2468_ACE0);
        let allele = dice.sequence(200);
        // The genome lacks three bases of the allele and holds two more.
        let mut in_genome = allele.clone();
        in_genome.splice(150..150, *b"AC");
        in_genome.drain(60..63);
        let genome = record(&mut dice, 9_000, &[(5_000, &in_genome)]);
        let search = GenomeSearch::new(&[text(&allele)], 0.8).unwrap();
        search.search(0, true, &genome);
        let places = search.finish();
        let on = places[0].on_target.as_ref().unwrap();
        assert_eq!((on.start, on.differences), (5_000, 5));
        assert_eq!(on.reference.as_bytes(), in_genome);
        // Three bases only in the sequence, two only in the genome.
        let count = |letter: char| {
            let mut total = 0;
            let mut run = 0;
            for c in on.cigar.chars() {
                match c.to_digit(10) {
                    Some(digit) => run = run * 10 + digit as usize,
                    None => {
                        if c == letter {
                            total += run;
                        }
                        run = 0;
                    }
                }
            }
            total
        };
        assert_eq!((count('I'), count('D'), count('X')), (3, 2, 0));
    }

    #[test]
    fn a_place_too_far_from_the_sequence_does_not_count() {
        let mut dice = Dice(0xAAAA_5555_AAAA_5555);
        let allele = dice.sequence(200);
        // Three changed bases in every sixteen: 36 differences in 200
        // bases, beyond one in ten and within one in five.
        let changed: Vec<usize> = (0..200).filter(|at| (5..8).contains(&(at % 16))).collect();
        let distant = dice.substituted(&allele, &changed);
        let genome = record(&mut dice, 9_000, &[(2_000, &distant)]);
        let strict = GenomeSearch::new(&[text(&allele)], 0.9).unwrap();
        strict.search(0, true, &genome);
        assert!(strict.finish()[0].on_target.is_none());
        let lenient = GenomeSearch::new(&[text(&allele)], 0.8).unwrap();
        lenient.search(0, true, &genome);
        let places = lenient.finish();
        let found = places[0].on_target.as_ref().unwrap();
        assert_eq!(found.start, 2_000);
        assert!(
            (30..=39).contains(&found.differences),
            "{}",
            found.differences
        );
    }

    #[test]
    fn short_sequences_are_found_where_they_lie_unchanged() {
        let mut dice = Dice(0x0F1E_2D3C_4B5A_6978);
        let d_gene = b"GGTATAACTGGAACGAC";
        let turned = reverse_complement(d_gene);
        let changed = dice.substituted(d_gene, &[8]);
        let genome = record(&mut dice, 12_000, &[(3_000, &changed), (8_000, &turned)]);
        let search = GenomeSearch::new(&[text(d_gene)], 0.8).unwrap();
        search.search(0, true, &genome);
        let places = search.finish();
        let on = places[0].on_target.as_ref().unwrap();
        // The changed copy is not looked for; the turned one is found.
        assert_eq!(
            (on.start, on.minus_strand, on.differences),
            (8_000, true, 0)
        );
        assert_eq!(
            (on.cigar.as_str(), on.reference.as_bytes()),
            ("17=", &d_gene[..])
        );
    }

    #[test]
    fn several_sequences_keep_their_order_and_their_own_places() {
        let mut dice = Dice(0x7777_8888_9999_AAAA);
        let first = dice.sequence(120);
        let second = dice.substituted(&first, &[60]);
        let third = dice.sequence(150);
        let genome = record(&mut dice, 20_000, &[(2_000, &second), (9_000, &first)]);
        let sequences = [text(&first), text(&second), text(&third)];
        let search = GenomeSearch::new(&sequences, 0.8).unwrap();
        search.search(0, true, &genome);
        let places = search.finish();
        let start = |i: usize| {
            let place = places[i].on_target.as_ref();
            place.map(|p| (p.start, p.differences))
        };
        assert_eq!(start(0), Some((9_000, 0)));
        assert_eq!(start(1), Some((2_000, 0)));
        assert_eq!(start(2), None);
    }
}
