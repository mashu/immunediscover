//! Finding, for every read, the database sequence it matches best.
//!
//! A read is compared with the database in both orientations. Subjects
//! that share words with the read are aligned in a band around the
//! diagonal the words point to ([`super::seeds`], [`super::gotoh`]);
//! subjects too short to be found by words are aligned in whole. The
//! alignments are ranked like BLAST hits are ranked by `discover blast`:
//! highest identity first, then the larger covered part of the subject,
//! then the higher score, then the order of the database.
//!
//! Not every subject that shares a word is aligned. A difference between
//! read and subject removes at most a known number of shared words, so a
//! subject with too few shared words cannot match the read as well as the
//! best alignment found so far, and is passed over.

use std::collections::HashMap;

use rayon::prelude::*;

use super::gotoh::{align, align_widening, Band, Ends, Op, Path, Scoring, Workspace};
use super::seeds::{Candidate, Votes, WordIndex, GROUP, WORD};
use crate::error::{Error, Result};

/// Subjects shorter than this are aligned to every read in whole: they
/// hold too few words to be found by them.
const SHORT_SUBJECT: usize = 3 * WORD;
/// Half width of the band an alignment starts with, and the widest tried.
const HALF_WIDTH: usize = 12;
const MAX_HALF_WIDTH: usize = 96;
/// Fewer shared words than this are chance.
const MIN_VOTES: u32 = 2;
/// Most subjects found by words that are aligned to one read.
const MAX_ALIGNED: usize = 256;
/// A read whose subjects with the most shared words all fail the
/// settings matches nothing; the search stops after this many.
const MAX_FAILED: usize = 8;
/// Most words of a read that are looked up between two that are skipped.
const MAX_STEP: usize = 4;
/// A subject of this length is covered by about this many looked-up words.
const WORDS_PER_SUBJECT: usize = 40;

/// What counts as a match.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Settings {
    pub scoring: Scoring,
    /// Smallest share of equal bases among the aligned columns (0-1).
    pub min_identity: f64,
    /// Smallest aligned share of the subject (0-1).
    pub min_coverage: f64,
    /// Whether a read may end inside a subject.
    pub ends: Ends,
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            scoring: Scoring::default(),
            min_identity: 0.8,
            min_coverage: 0.5,
            ends: Ends::ReadMayEndInside,
        }
    }
}

/// A database sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Subject {
    pub name: String,
    /// Upper-case nucleotides.
    pub sequence: Vec<u8>,
}

/// The best match of a read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Hit {
    /// Number of the subject in the database.
    pub subject: usize,
    /// The reverse complement of the read matches the subject.
    pub minus_strand: bool,
    pub read_length: usize,
    /// The alignment; read positions refer to the read in the orientation
    /// that matches the subject.
    pub path: Path,
}

impl Hit {
    /// Share of equal bases among the aligned columns.
    pub fn identity(&self) -> f64 {
        self.path.matches() as f64 / self.path.columns() as f64
    }

    /// Aligned part of the read as written in the input,
    /// `start..end` counted from 0.
    pub fn read_range(&self) -> (usize, usize) {
        if self.minus_strand {
            (
                self.read_length - self.path.read_end,
                self.read_length - self.path.read_start,
            )
        } else {
            (self.path.read_start, self.path.read_end)
        }
    }

    /// Number of aligned subject bases.
    pub fn subject_span(&self) -> usize {
        self.path.subject_end - self.path.subject_start
    }

    /// The aligned part of `read` (the read this hit was found for),
    /// written like the subject: upper case, and reverse-complemented when
    /// the read matches the other strand.
    pub fn segment(&self, read: &str) -> String {
        let (start, end) = self.read_range();
        let part = &read.trim().as_bytes()[start..end];
        let bytes: Vec<u8> = if self.minus_strand {
            let mut turned = Vec::new();
            reverse_complement_into(&part.to_ascii_uppercase(), &mut turned);
            turned
        } else {
            part.to_ascii_uppercase()
        };
        String::from_utf8_lossy(&bytes).into_owned()
    }

    /// [`Hit::segment`] with `-` where the subject has a base the read
    /// lacks.
    pub fn gapped_segment(&self, read: &str) -> String {
        let segment = self.segment(read);
        let mut bases = segment.chars();
        self.path
            .ops
            .iter()
            .map(|op| match op {
                Op::Deletion => '-',
                _ => bases.next().unwrap_or('N'),
            })
            .collect()
    }

    /// The aligned read bases as BLAST writes them (`qseq`): in the
    /// orientation of the read, with gaps.
    pub fn gapped_read(&self, read: &str) -> String {
        let gapped = self.gapped_segment(read);
        if !self.minus_strand {
            return gapped;
        }
        gapped
            .bytes()
            .rev()
            .map(|base| match base {
                b'A' => 'T',
                b'C' => 'G',
                b'G' => 'C',
                b'T' => 'A',
                b'-' => '-',
                _ => 'N',
            })
            .collect()
    }
}

/// Counters and buffers of one thread.
#[derive(Debug, Default)]
pub struct Scratch {
    workspace: Workspace,
    votes: Votes,
    candidates: Vec<Candidate>,
    /// Numbers of `candidates`, best first.
    order: Vec<usize>,
    forward: Vec<u8>,
    reverse: Vec<u8>,
}

/// Upper-case reverse complement; anything but A, C, G and T becomes N.
fn reverse_complement_into(sequence: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.extend(sequence.iter().rev().map(|base| match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }));
}

/// Rank of an alignment; the larger, the better.
#[derive(Debug, Clone, Copy)]
struct Rank {
    matches: u64,
    columns: u64,
    span: u64,
    subject_length: u64,
    score: i32,
}

impl Rank {
    fn of(path: &Path, subject_length: usize) -> Self {
        Rank {
            matches: path.matches() as u64,
            columns: path.columns() as u64,
            span: (path.subject_end - path.subject_start) as u64,
            subject_length: subject_length as u64,
            score: path.score,
        }
    }

    /// Whether `self` is strictly better than `other`. Fractions are
    /// compared exactly, so that the answer never depends on rounding.
    fn beats(&self, other: &Rank) -> bool {
        let identity = (self.matches * other.columns).cmp(&(other.matches * self.columns));
        let coverage = (self.span * other.subject_length).cmp(&(other.span * self.subject_length));
        identity
            .then(coverage)
            .then(self.score.cmp(&other.score))
            .is_gt()
    }

    fn passes(&self, settings: &Settings) -> bool {
        self.matches as f64 >= settings.min_identity * self.columns as f64
            && self.span as f64 >= settings.min_coverage * self.subject_length as f64
    }
}

/// A set of sequences prepared for searching.
#[derive(Debug)]
pub struct Database {
    subjects: Vec<Subject>,
    index: WordIndex,
    /// Subjects aligned to every read in whole.
    short: Vec<usize>,
    /// Every `step`-th word of a read is looked up.
    step: usize,
}

impl Database {
    /// Prepare `(name, sequence)` pairs for searching.
    pub fn new<N, S>(sequences: impl IntoIterator<Item = (N, S)>) -> Result<Self>
    where
        N: Into<String>,
        S: AsRef<str>,
    {
        let subjects: Vec<Subject> = sequences
            .into_iter()
            .map(|(name, sequence)| Subject {
                name: name.into(),
                sequence: sequence.as_ref().trim().to_ascii_uppercase().into_bytes(),
            })
            .collect();
        if subjects.is_empty() {
            return Err(Error::invalid("The database holds no sequences"));
        }
        if let Some(empty) = subjects.iter().find(|s| s.sequence.is_empty()) {
            return Err(Error::invalid(format!(
                "Database sequence '{}' is empty",
                empty.name
            )));
        }
        let is_short = |subject: &Subject| subject.sequence.len() < SHORT_SUBJECT;
        let index = WordIndex::new(
            subjects
                .iter()
                .map(|s| (!is_short(s)).then_some(s.sequence.as_slice())),
        );
        let short: Vec<usize> = (0..subjects.len())
            .filter(|&i| is_short(&subjects[i]))
            .collect();
        // Short subjects need more of the read's words to be looked up.
        let shortest = subjects
            .iter()
            .filter(|s| !is_short(s))
            .map(|s| s.sequence.len())
            .min()
            .unwrap_or(SHORT_SUBJECT);
        let step = (shortest / WORDS_PER_SUBJECT).clamp(1, MAX_STEP);
        Ok(Database {
            subjects,
            index,
            short,
            step,
        })
    }

    pub fn subjects(&self) -> &[Subject] {
        &self.subjects
    }

    pub fn subject(&self, number: usize) -> &Subject {
        &self.subjects[number]
    }

    /// The best match of one read, if any passes `settings`.
    pub fn best_hit(&self, read: &str, settings: &Settings, scratch: &mut Scratch) -> Option<Hit> {
        let mut forward = std::mem::take(&mut scratch.forward);
        let mut reverse = std::mem::take(&mut scratch.reverse);
        forward.clear();
        forward.extend(read.trim().bytes().map(|b| b.to_ascii_uppercase()));
        reverse_complement_into(&forward, &mut reverse);
        let best = self.search_strands([&forward, &reverse], settings, scratch);
        scratch.forward = forward;
        scratch.reverse = reverse;
        best.map(|(_, hit)| hit)
    }

    /// Keep `path` in `best` if it passes `settings` and is better than
    /// what is there. Of equal alignments the earlier subject, then the
    /// read as written, stays.
    fn consider(
        &self,
        best: &mut Option<(Rank, Hit)>,
        settings: &Settings,
        subject: usize,
        minus_strand: bool,
        read_length: usize,
        path: Option<Path>,
    ) {
        let Some(path) = path else { return };
        let rank = Rank::of(&path, self.subjects[subject].sequence.len());
        if !rank.passes(settings) {
            return;
        }
        let better = match best {
            None => true,
            Some((top, hit)) => {
                rank.beats(top)
                    || (!top.beats(&rank)
                        && (subject, minus_strand) < (hit.subject, hit.minus_strand))
            }
        };
        if better {
            *best = Some((
                rank,
                Hit {
                    subject,
                    minus_strand,
                    read_length,
                    path,
                },
            ));
        }
    }

    /// Align the subjects worth aligning to the two orientations of a
    /// read (as written, and reverse-complemented) and return the best
    /// alignment.
    fn search_strands(
        &self,
        strands: [&[u8]; 2],
        settings: &Settings,
        scratch: &mut Scratch,
    ) -> Option<(Rank, Hit)> {
        let mut best: Option<(Rank, Hit)> = None;
        let read_length = strands[0].len();

        for &subject in &self.short {
            let sequence = &self.subjects[subject].sequence;
            for (minus_strand, read) in [(false, strands[0]), (true, strands[1])] {
                let band = Band::whole(sequence.len(), read.len());
                let path = align(
                    sequence,
                    read,
                    &settings.scoring,
                    band,
                    settings.ends,
                    &mut scratch.workspace,
                );
                self.consider(
                    &mut best,
                    settings,
                    subject,
                    minus_strand,
                    read_length,
                    path,
                );
            }
        }

        // Subjects sharing words with either orientation, most words first.
        scratch.candidates.clear();
        scratch
            .votes
            .collect(&self.index, strands[0], self.step, &mut scratch.candidates);
        let forward_candidates = scratch.candidates.len();
        scratch
            .votes
            .collect(&self.index, strands[1], self.step, &mut scratch.candidates);
        scratch.order.clear();
        scratch.order.extend(0..scratch.candidates.len());
        let candidates = &scratch.candidates;
        scratch.order.sort_unstable_by(|&a, &b| {
            (candidates[b].votes.cmp(&candidates[a].votes))
                .then(candidates[a].subject.cmp(&candidates[b].subject))
                .then(a.cmp(&b))
        });
        let top_votes = scratch.order.first().map_or(0, |&i| candidates[i].votes);
        // A difference removes at most this many of the looked-up words.
        let lost_per_difference = WORD.div_ceil(self.step);
        let mut top: Option<Rank> = best.as_ref().map(|(rank, _)| *rank);
        let mut aligned = 0;
        for &number in &scratch.order {
            let candidate = candidates[number];
            let hopeless = top.is_none() && aligned >= MAX_FAILED;
            if candidate.votes < MIN_VOTES || aligned >= MAX_ALIGNED || hopeless {
                break;
            }
            let minus_strand = number >= forward_candidates;
            let read = strands[usize::from(minus_strand)];
            let sequence = &self.subjects[candidate.subject].sequence;
            // Part of the read the subject would lie on.
            let start = candidate.diagonal.max(0);
            let end = (candidate.diagonal + sequence.len() as isize).min(read.len() as isize);
            let span = (end - start).max(0) as usize;
            // Words an identical subject would share, counted cautiously,
            // and differences this subject may have to match the best.
            let expected = (span.saturating_sub(WORD + GROUP) / self.step).min(top_votes as usize);
            let allowed = match &top {
                Some(top) => {
                    let differences = (top.columns - top.matches) as usize;
                    (differences * span).div_ceil(top.columns as usize) + 1
                }
                None => ((1.0 - settings.min_identity) * span as f64) as usize + 1,
            };
            if candidate.votes as usize + lost_per_difference * allowed < expected {
                continue;
            }
            aligned += 1;
            let path = align_widening(
                sequence,
                read,
                &settings.scoring,
                candidate.diagonal,
                HALF_WIDTH,
                MAX_HALF_WIDTH,
                settings.ends,
                &mut scratch.workspace,
            );
            self.consider(
                &mut best,
                settings,
                candidate.subject,
                minus_strand,
                read_length,
                path,
            );
            top = best.as_ref().map(|(rank, _)| *rank);
        }
        best
    }

    /// The best match of every read. Equal reads are aligned once; the
    /// work is spread over all processor cores, and the result does not
    /// depend on their number.
    pub fn search(&self, reads: &[&str], settings: &Settings) -> Hits {
        let mut number_of: HashMap<&str, usize> = HashMap::new();
        let mut distinct: Vec<&str> = Vec::new();
        let of_read: Vec<usize> = reads
            .iter()
            .map(|&read| {
                *number_of.entry(read).or_insert_with(|| {
                    distinct.push(read);
                    distinct.len() - 1
                })
            })
            .collect();
        let hits: Vec<Option<Hit>> = distinct
            .par_iter()
            .map_init(Scratch::default, |scratch, read| {
                self.best_hit(read, settings, scratch)
            })
            .collect();
        Hits { hits, of_read }
    }
}

/// The matches of a list of reads.
#[derive(Debug)]
pub struct Hits {
    /// One entry per distinct read sequence.
    hits: Vec<Option<Hit>>,
    of_read: Vec<usize>,
}

impl Hits {
    /// The match of read number `read`.
    pub fn get(&self, read: usize) -> Option<&Hit> {
        self.hits[self.of_read[read]].as_ref()
    }

    /// Number of reads.
    pub fn len(&self) -> usize {
        self.of_read.len()
    }

    pub fn is_empty(&self) -> bool {
        self.of_read.is_empty()
    }

    /// Number of distinct read sequences.
    pub fn distinct(&self) -> usize {
        self.hits.len()
    }

    /// Number of reads with a match.
    pub fn matched(&self) -> usize {
        self.of_read
            .iter()
            .filter(|&&number| self.hits[number].is_some())
            .count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::dna::reverse_complement;

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

        fn sequence(&mut self, length: usize) -> String {
            (0..length)
                .map(|_| b"ACGT"[self.below(4)] as char)
                .collect()
        }

        fn substituted(&mut self, sequence: &str, at: &[usize]) -> String {
            let mut bytes = sequence.as_bytes().to_vec();
            for &position in at {
                bytes[position] = match bytes[position] {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    _ => b'A',
                };
            }
            String::from_utf8(bytes).unwrap()
        }
    }

    /// Five genes of 300 bases with three alleles each, one or two
    /// substitutions apart.
    fn alleles(dice: &mut Dice) -> Vec<(String, String)> {
        alleles_of_length(dice, 300)
    }

    fn alleles_of_length(dice: &mut Dice, length: usize) -> Vec<(String, String)> {
        let mut alleles = Vec::new();
        for gene in 1..=5 {
            let first = dice.sequence(length);
            let second = dice.substituted(&first, &[length / 5]);
            let third = dice.substituted(&first, &[length / 5, 2 * length / 3]);
            for (number, sequence) in [first, second, third].into_iter().enumerate() {
                alleles.push((format!("IGHV{gene}-1*0{}", number + 1), sequence));
            }
        }
        alleles
    }

    fn read_of(dice: &mut Dice, allele: &str) -> String {
        let (before, after) = (20 + dice.below(60), 20 + dice.below(60));
        format!("{}{allele}{}", dice.sequence(before), dice.sequence(after))
    }

    /// Every subject aligned in whole: what the search must agree with.
    fn exhaustive(database: &Database, read: &str, settings: &Settings) -> Option<Hit> {
        let forward = read.as_bytes().to_vec();
        let reverse = reverse_complement(read).unwrap().into_bytes();
        let mut workspace = Workspace::default();
        let mut best: Option<(Rank, Hit)> = None;
        for (number, subject) in database.subjects().iter().enumerate() {
            for (minus_strand, oriented) in [(false, &forward), (true, &reverse)] {
                let band = Band::whole(subject.sequence.len(), oriented.len());
                let Some(path) = align(
                    &subject.sequence,
                    oriented,
                    &settings.scoring,
                    band,
                    settings.ends,
                    &mut workspace,
                ) else {
                    continue;
                };
                let rank = Rank::of(&path, subject.sequence.len());
                if rank.passes(settings) && best.as_ref().is_none_or(|(top, _)| rank.beats(top)) {
                    best = Some((
                        rank,
                        Hit {
                            subject: number,
                            minus_strand,
                            read_length: oriented.len(),
                            path,
                        },
                    ));
                }
            }
        }
        best.map(|(_, hit)| hit)
    }

    #[test]
    fn every_read_finds_its_allele_on_either_strand() {
        let mut dice = Dice(0xA5A5_5A5A_1234_5678);
        let alleles = alleles(&mut dice);
        let database = Database::new(alleles.iter().map(|(n, s)| (n.as_str(), s))).unwrap();
        let mut scratch = Scratch::default();
        let settings = Settings::default();
        for (number, (_, allele)) in alleles.iter().enumerate() {
            let read = read_of(&mut dice, allele);
            let hit = database.best_hit(&read, &settings, &mut scratch).unwrap();
            assert_eq!(hit.subject, number);
            assert!(!hit.minus_strand);
            assert_eq!(hit.path.differences(), 0);
            let (start, end) = hit.read_range();
            assert_eq!(&read[start..end], allele);

            let turned = reverse_complement(&read).unwrap();
            let hit = database.best_hit(&turned, &settings, &mut scratch).unwrap();
            assert_eq!(hit.subject, number);
            assert!(hit.minus_strand);
            let (start, end) = hit.read_range();
            assert_eq!(reverse_complement(&turned[start..end]).unwrap(), *allele);
        }
    }

    #[test]
    fn a_read_between_two_alleles_goes_to_the_nearer_one() {
        let mut dice = Dice(0x0123_4567_89AB_CDEF);
        let alleles = alleles(&mut dice);
        let database = Database::new(alleles.iter().map(|(n, s)| (n.as_str(), s))).unwrap();
        // The third allele of the second gene, with one new substitution.
        let novel = dice.substituted(&alleles[5].1, &[120]);
        let read = read_of(&mut dice, &novel);
        let hit = database
            .best_hit(&read, &Settings::default(), &mut Scratch::default())
            .unwrap();
        assert_eq!(database.subject(hit.subject).name, "IGHV2-1*03");
        assert_eq!(hit.path.mismatches(), 1);
        assert_eq!(hit.path.cigar(), "120=1X179=");
    }

    #[test]
    fn equally_near_alleles_are_decided_by_database_order() {
        let mut dice = Dice(0x1111_2222_3333_4444);
        let alleles = alleles(&mut dice);
        let database = Database::new(alleles.iter().map(|(n, s)| (n.as_str(), s))).unwrap();
        // Half way between the second and the third allele of a gene: the
        // read has neither the base of *02 nor that of *03 at position 200.
        let mut bytes = alleles[1].1.as_bytes().to_vec();
        let (second, third) = (bytes[200], alleles[2].1.as_bytes()[200]);
        bytes[200] = *b"ACGT"
            .iter()
            .find(|&&b| b != second && b != third)
            .unwrap();
        let read = read_of(&mut dice, std::str::from_utf8(&bytes).unwrap());
        let hit = database
            .best_hit(&read, &Settings::default(), &mut Scratch::default())
            .unwrap();
        assert_eq!(database.subject(hit.subject).name, "IGHV1-1*02");
    }

    #[test]
    fn the_search_agrees_with_aligning_every_subject_in_whole() {
        let mut dice = Dice(0xDEAD_BEEF_CAFE_F00D);
        let alleles = alleles_of_length(&mut dice, 150);
        let database = Database::new(alleles.iter().map(|(n, s)| (n.as_str(), s))).unwrap();
        let settings = Settings::default();
        let mut scratch = Scratch::default();
        for round in 0..60 {
            let (_, allele) = &alleles[dice.below(alleles.len())];
            // Reads with substitutions, with a gap, and reads cut short.
            let mut sequence = allele.clone();
            for _ in 0..dice.below(6) {
                let at = dice.below(sequence.len());
                sequence = dice.substituted(&sequence, &[at]);
            }
            if round % 5 == 0 {
                let at = 30 + dice.below(90);
                sequence.replace_range(at..at + 2, "");
            }
            let mut read = read_of(&mut dice, &sequence);
            if round % 7 == 0 {
                read.truncate(read.len() - 70);
            }
            if round % 2 == 0 {
                read = reverse_complement(&read).unwrap();
            }
            let found = database.best_hit(&read, &settings, &mut scratch);
            let expected = exhaustive(&database, &read, &settings);
            let summary = |hit: &Option<Hit>| {
                hit.as_ref().map(|h| {
                    (
                        h.subject,
                        h.minus_strand,
                        h.path.score,
                        h.path.cigar(),
                        h.read_range(),
                    )
                })
            };
            assert_eq!(summary(&found), summary(&expected), "round {round}");
        }
    }

    #[test]
    fn short_subjects_are_found_without_shared_words() {
        // D genes: shorter than three words, one differing in one base.
        let database = Database::new([
            ("IGHD1*01", "GGTATAACTGGAACGAC"),
            ("IGHD2*01", "AGGATATTGTAGTGGTGGT"),
            ("IGHD2*02", "AGGATATTGTACTGGTGGT"),
        ])
        .unwrap();
        let mut dice = Dice(99);
        let read = format!(
            "{}AGGATATTGTACTGGTGGT{}",
            dice.sequence(150),
            dice.sequence(150)
        );
        let settings = Settings::default();
        let hit = database
            .best_hit(&read, &settings, &mut Scratch::default())
            .unwrap();
        assert_eq!(database.subject(hit.subject).name, "IGHD2*02");
        assert_eq!(hit.read_range(), (150, 169));
        let turned = reverse_complement(&read).unwrap();
        let hit = database
            .best_hit(&turned, &settings, &mut Scratch::default())
            .unwrap();
        assert!(hit.minus_strand);
        assert_eq!(hit.read_range(), (150, 169));
    }

    #[test]
    fn segments_are_written_like_the_subject_and_like_the_read() {
        let database =
            Database::new([("gene", "ACGGTCATTGCAGGCTTAACGTTAGCCATGCAAGTCCGATAGGC")]).unwrap();
        // The read lacks the two bases GG of the gene and lies on the
        // other strand.
        let forward = "ttttttACGGTCATTGCACTTAACGTTAGCCATGCAAGTCCGATAGGCaaaaa";
        let read = reverse_complement(forward).unwrap().to_lowercase();
        let hit = database
            .best_hit(&read, &Settings::default(), &mut Scratch::default())
            .unwrap();
        assert!(hit.minus_strand);
        assert_eq!(hit.path.cigar(), "12=2D30=");
        assert_eq!(
            hit.segment(&read),
            "ACGGTCATTGCACTTAACGTTAGCCATGCAAGTCCGATAGGC"
        );
        assert_eq!(
            hit.gapped_segment(&read),
            "ACGGTCATTGCA--CTTAACGTTAGCCATGCAAGTCCGATAGGC"
        );
        assert_eq!(
            hit.gapped_read(&read),
            "GCCTATCGGACTTGCATGGCTAACGTTAAG--TGCAATGACCGT"
        );
        assert_eq!(hit.read_range(), (5, 47));
    }

    #[test]
    fn unrelated_reads_have_no_match() {
        let mut dice = Dice(0x5555_AAAA_5555_AAAA);
        let alleles = alleles(&mut dice);
        let database = Database::new(alleles.iter().map(|(n, s)| (n.as_str(), s))).unwrap();
        let read = dice.sequence(400);
        assert!(database
            .best_hit(&read, &Settings::default(), &mut Scratch::default())
            .is_none());
    }

    #[test]
    fn equal_reads_are_aligned_once_and_all_get_the_hit() {
        let mut dice = Dice(0x0F0F_F0F0_0F0F_F0F0);
        let alleles = alleles(&mut dice);
        let database = Database::new(alleles.iter().map(|(n, s)| (n.as_str(), s))).unwrap();
        let read = read_of(&mut dice, &alleles[7].1);
        let other = dice.sequence(300);
        let reads = [read.as_str(), other.as_str(), read.as_str()];
        let hits = database.search(&reads, &Settings::default());
        assert_eq!((hits.len(), hits.distinct(), hits.matched()), (3, 2, 2));
        assert_eq!(hits.get(0).unwrap().subject, 7);
        assert!(hits.get(1).is_none());
        assert_eq!(hits.get(0), hits.get(2));
    }

    #[test]
    fn an_empty_database_is_refused() {
        assert!(Database::new(Vec::<(String, String)>::new()).is_err());
        assert!(Database::new([("empty", "")]).is_err());
    }
}
