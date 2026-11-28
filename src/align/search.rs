//! Finding, for every read, the database sequence nearest to it.
//!
//! A read is compared with the database in both orientations. The nearest
//! subject is the one that can be turned into a part of the read with the
//! fewest differences (mismatches, inserted and deleted bases) per base;
//! between equally near ones the longer subject wins, then the order of
//! the database. A subject that reaches over an end of the read is
//! compared as far as the read goes; every base beyond counts as a
//! fraction of a difference, so that a subject that fits is preferred.
//!
//! Finding the nearest subject takes no alignment:
//!
//! * Subjects that share words with the read are candidates
//!   ([`super::seeds`]). Subjects that are nearly the same (the alleles of
//!   a gene) form a group and share one set of votes.
//! * A subject that lies in the read with at most one mismatch is seen by
//!   comparing it base by base.
//! * For the others the least number of differences is counted
//!   ([`super::edits`]), which costs little when the number is small, and
//!   only up to the number that the nearest subject so far has.
//!
//! The nearest subject is then aligned once ([`super::gotoh`]) to place
//! its gaps. Subjects too short to be found by words are aligned to every
//! read.

use std::collections::HashMap;

use rayon::prelude::*;

use super::edits::{differences, mismatches, Reach};
use super::gotoh::{align, align_widening, is_same_base, Band, Ends, Op, Path, Scoring, Workspace};
use super::seeds::{Candidate, Votes, WordIndex, GROUP, WORD};
use crate::error::{Error, Result};
use crate::logging::Progress;
use crate::seq::search::PatternSet;

/// Subjects shorter than this are aligned to every read in whole: they
/// hold too few words to be found by them.
pub(super) const SHORT_SUBJECT: usize = 3 * WORD;
/// Half width of the band an alignment starts with, and the widest tried.
/// The diagonal found from shared words is off by less than this.
pub(super) const HALF_WIDTH: usize = 12;
pub(super) const MAX_HALF_WIDTH: usize = 96;
/// Fewer shared words than this are chance.
const MIN_VOTES: u32 = 2;
/// Most words of a read that are looked up between two that are skipped.
const MAX_STEP: usize = 4;
/// A subject of this length is covered by about this many looked-up words.
const WORDS_PER_SUBJECT: usize = 40;
/// Most mismatches between a subject and the first subject of its group.
const GROUP_RADIUS: usize = 8;
/// More differences than this are not counted: a read this far from a
/// subject does not match it.
pub(super) const MAX_COUNTED: usize = 128;

/// What counts as a match.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Settings {
    /// Scores that place the gaps of the reported alignment. They also
    /// set what a base beyond the end of the read counts for: the score
    /// of a match relative to the loss of a mismatch.
    pub scoring: Scoring,
    /// A read matches a subject with at most `1 - min_identity`
    /// differences per subject base (0-1).
    pub min_identity: f64,
    /// Smallest share of the subject that must lie inside the read (0-1).
    pub min_coverage: f64,
    /// Whether a read may end inside a subject.
    pub ends: Ends,
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            scoring: Scoring::default(),
            min_identity: 0.9,
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
    reach: Reach,
    votes: Votes,
    candidates: Vec<Candidate>,
    /// Numbers of `candidates`, best first.
    order: Vec<usize>,
    /// Members of a group still to be judged.
    pending: Vec<usize>,
    /// Diagonals on which members of a group may lie nearly unchanged.
    places: Vec<isize>,
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

/// The alignment of a subject that lies in the read from `start` without
/// a gap.
fn gapless(subject: &[u8], read: &[u8], start: usize, scoring: &Scoring) -> Path {
    let ops: Vec<Op> = subject
        .iter()
        .zip(&read[start..start + subject.len()])
        .map(|(&a, &b)| {
            if is_same_base(a, b) {
                Op::Match
            } else {
                Op::Mismatch
            }
        })
        .collect();
    let mismatched = ops.iter().filter(|&&op| op == Op::Mismatch).count();
    let score =
        scoring.matched * (ops.len() - mismatched) as i32 + scoring.mismatched * mismatched as i32;
    Path {
        score,
        subject_start: 0,
        subject_end: subject.len(),
        read_start: start,
        read_end: start + subject.len(),
        ops,
        at_band_edge: false,
    }
}

/// Align `subject` to `read` on the diagonals `band`, and beyond them
/// for as long as the alignment runs along their border.
pub(super) fn align_in(
    band: (isize, isize),
    subject: &[u8],
    read: &[u8],
    settings: &Settings,
    workspace: &mut Workspace,
) -> Option<Path> {
    let half_width = ((band.1 - band.0 + 1) / 2).max(1) as usize;
    align_widening(
        subject,
        read,
        &settings.scoring,
        (band.0 + band.1) / 2,
        half_width,
        MAX_HALF_WIDTH.max(half_width),
        settings.ends,
        workspace,
    )
}

/// How near a subject is to a read: the smaller `cost / length`, the
/// nearer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Nearness {
    /// Differences and bases beyond the read, weighted (see [`Costs`]).
    cost: i64,
    /// Length of the subject.
    length: i64,
}

impl Nearness {
    /// Whether `self` is nearer than `other`; of equally near subjects
    /// the longer is nearer. Shares are compared exactly, so that the
    /// answer never depends on rounding.
    fn nearer_than(&self, other: &Nearness) -> bool {
        (self.cost * other.length)
            .cmp(&(other.cost * self.length))
            .then(other.length.cmp(&self.length))
            .is_lt()
    }

    /// The largest cost at which a subject of `length` is as near.
    fn budget_for(&self, length: usize) -> i64 {
        self.cost * length as i64 / self.length
    }
}

/// What differences and bases beyond the read count for.
#[derive(Debug, Clone, Copy)]
struct Costs {
    /// A mismatch, an inserted or a deleted base.
    difference: i64,
    /// A subject base beyond the end of the read.
    beyond: i64,
}

impl Costs {
    fn of(scoring: &Scoring) -> Self {
        Costs {
            difference: i64::from(scoring.matched - scoring.mismatched),
            beyond: i64::from(scoring.matched),
        }
    }

    /// The cost of an alignment.
    fn of_path(&self, path: &Path, length: usize) -> i64 {
        let beyond = length - (path.subject_end - path.subject_start);
        self.difference * path.differences() as i64 + self.beyond * beyond as i64
    }
}

/// The nearest subject found so far.
#[derive(Debug)]
struct Nearest {
    nearness: Nearness,
    subject: usize,
    minus_strand: bool,
    /// Where to align: the lowest and the highest diagonal (read position
    /// minus subject position) that were tried.
    band: (isize, isize),
    /// The alignment, where it is already known.
    path: Option<Path>,
}

/// How many differences and bases beyond the read a subject may have.
#[derive(Debug, Clone, Copy)]
struct Limits {
    cost: i64,
    differences: usize,
    beyond: usize,
}

/// Subjects of one length that differ from the first of them in few
/// positions.
#[derive(Debug)]
pub(super) struct Group {
    /// Numbers of the subjects, in the order of the database.
    pub(super) members: Vec<usize>,
    pub(super) length: usize,
    /// Most mismatches between a member and the first member.
    radius: usize,
}

/// A set of sequences prepared for searching.
#[derive(Debug)]
pub struct Database {
    subjects: Vec<Subject>,
    pub(super) groups: Vec<Group>,
    /// Words of the groups.
    pub(super) index: WordIndex,
    /// Subjects aligned to every read in whole.
    short: Vec<usize>,
    /// Every `step`-th word of a read is looked up.
    step: usize,
    /// Every subject and its reverse complement, to find reads that hold
    /// a subject unchanged; pattern `2 * s` is subject `s` as written,
    /// `2 * s + 1` its reverse complement.
    unchanged: PatternSet,
    /// Shortest subject found by words.
    shortest: usize,
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
        let short: Vec<usize> = (0..subjects.len())
            .filter(|&i| is_short(&subjects[i]))
            .collect();

        // Every longer subject joins the first group it is close to.
        let mut groups: Vec<Group> = Vec::new();
        for (number, subject) in subjects.iter().enumerate() {
            if is_short(subject) {
                continue;
            }
            let sequence = &subject.sequence;
            let home = groups.iter_mut().find_map(|group| {
                if group.length != sequence.len() {
                    return None;
                }
                let first = &subjects[group.members[0]].sequence;
                mismatches(first, sequence, GROUP_RADIUS).map(|distance| (group, distance))
            });
            match home {
                Some((group, distance)) => {
                    group.members.push(number);
                    group.radius = group.radius.max(distance);
                }
                None => groups.push(Group {
                    members: vec![number],
                    length: sequence.len(),
                    radius: 0,
                }),
            }
        }
        let grouped: Vec<Vec<&[u8]>> = groups
            .iter()
            .map(|group| {
                group
                    .members
                    .iter()
                    .map(|&number| subjects[number].sequence.as_slice())
                    .collect()
            })
            .collect();
        let index = WordIndex::new(&grouped);
        // Short subjects need more of the read's words to be looked up.
        let shortest = groups
            .iter()
            .map(|g| g.length)
            .min()
            .unwrap_or(SHORT_SUBJECT);
        let step = (shortest / WORDS_PER_SUBJECT).clamp(1, MAX_STEP);
        // A subject with an ambiguous base never lies unchanged in a
        // read; a pattern that cannot occur stands in for it.
        let mut patterns: Vec<Vec<u8>> = Vec::with_capacity(2 * subjects.len());
        for subject in &subjects {
            let plain = subject.sequence.iter().all(|&b| is_same_base(b, b));
            let mut turned = Vec::new();
            reverse_complement_into(&subject.sequence, &mut turned);
            if plain {
                patterns.push(subject.sequence.clone());
                patterns.push(turned);
            } else {
                patterns.push(vec![0]);
                patterns.push(vec![0]);
            }
        }
        let unchanged = PatternSet::new(patterns)?;
        Ok(Database {
            subjects,
            groups,
            index,
            short,
            step,
            unchanged,
            shortest,
        })
    }

    /// Number of groups of nearly equal subjects.
    pub fn groups(&self) -> usize {
        self.groups.len()
    }

    pub fn subjects(&self) -> &[Subject] {
        &self.subjects
    }

    pub fn subject(&self, number: usize) -> &Subject {
        &self.subjects[number]
    }

    /// The nearest subject of one read with its alignment, if any subject
    /// is within `settings`.
    pub fn best_hit(&self, read: &str, settings: &Settings, scratch: &mut Scratch) -> Option<Hit> {
        let mut forward = std::mem::take(&mut scratch.forward);
        let mut reverse = std::mem::take(&mut scratch.reverse);
        forward.clear();
        forward.extend(read.trim().bytes().map(|b| b.to_ascii_uppercase()));
        if let Some(hit) = self.unchanged_subject(&forward, &settings.scoring) {
            scratch.forward = forward;
            scratch.reverse = reverse;
            return Some(hit);
        }
        reverse_complement_into(&forward, &mut reverse);
        let hit = self
            .nearest([&forward, &reverse], settings, scratch)
            .and_then(|nearest| {
                let read = if nearest.minus_strand {
                    &reverse
                } else {
                    &forward
                };
                let sequence = &self.subjects[nearest.subject].sequence;
                let path = nearest.path.or_else(|| {
                    align_in(
                        nearest.band,
                        sequence,
                        read,
                        settings,
                        &mut scratch.workspace,
                    )
                })?;
                Some(Hit {
                    subject: nearest.subject,
                    minus_strand: nearest.minus_strand,
                    read_length: forward.len(),
                    path,
                })
            });
        scratch.forward = forward;
        scratch.reverse = reverse;
        hit
    }

    /// The subject that lies in the read without any difference, in either
    /// orientation: the longest one, then the first of the database, then
    /// the read as written. Found by one pass over the read, without
    /// looking for candidates.
    fn unchanged_subject(&self, read: &[u8], scoring: &Scoring) -> Option<Hit> {
        // (length, subject, reverse-complemented) and where the match ends.
        let mut best: Option<((usize, usize, bool), usize)> = None;
        self.unchanged.each_occurrence(read, |pattern, start, end| {
            let found = (end - start, pattern / 2, pattern % 2 == 1);
            let better = match &best {
                None => true,
                Some((held, _)) => {
                    (found.0, held.1, held.2) > (held.0, found.1, found.2)
                        // The same subject again: on the other strand the
                        // last occurrence is the first of the turned read.
                        || (found == *held && found.2)
                }
            };
            if better {
                best = Some((found, end));
            }
        });
        let ((length, subject, minus_strand), end) = best?;
        // Position in the read as it matches the subject.
        let start = if minus_strand {
            read.len() - end
        } else {
            end - length
        };
        Some(Hit {
            subject,
            minus_strand,
            read_length: read.len(),
            path: Path {
                score: scoring.matched * length as i32,
                subject_start: 0,
                subject_end: length,
                read_start: start,
                read_end: start + length,
                ops: vec![Op::Match; length],
                at_band_edge: false,
            },
        })
    }

    /// Keep a subject in `nearest` if it is nearer than what is there. Of
    /// equally near subjects the earlier one, then the read as written,
    /// stays.
    fn offer(nearest: &mut Option<Nearest>, offered: Nearest) {
        let wins = match nearest {
            None => true,
            Some(held) => {
                offered.nearness.nearer_than(&held.nearness)
                    || (!held.nearness.nearer_than(&offered.nearness)
                        && (offered.subject, offered.minus_strand)
                            < (held.subject, held.minus_strand))
            }
        };
        if wins {
            *nearest = Some(offered);
        }
    }

    /// What a subject of `length` may cost to be accepted by `settings`
    /// and to be as near as `nearest`.
    fn limits(&self, nearest: &Option<Nearest>, settings: &Settings, length: usize) -> Limits {
        let costs = Costs::of(&settings.scoring);
        let accepted = ((1.0 - settings.min_identity) * length as f64) as usize;
        let accepted = accepted.min(MAX_COUNTED);
        let outside = if settings.ends == Ends::ReadMayEndInside {
            ((1.0 - settings.min_coverage) * length as f64) as usize
        } else {
            0
        };
        let cost = match nearest {
            Some(held) => held.nearness.budget_for(length),
            None => i64::MAX / 4,
        };
        Limits {
            cost,
            differences: accepted.min((cost / costs.difference) as usize),
            beyond: outside.min((cost / costs.beyond) as usize),
        }
    }

    /// The nearest subject of the two orientations of a read (as written,
    /// and reverse-complemented).
    fn nearest(
        &self,
        strands: [&[u8]; 2],
        settings: &Settings,
        scratch: &mut Scratch,
    ) -> Option<Nearest> {
        let costs = Costs::of(&settings.scoring);
        let mut nearest: Option<Nearest> = None;

        // Short subjects: aligned in whole to both orientations.
        for &subject in &self.short {
            let sequence = &self.subjects[subject].sequence;
            for (minus_strand, read) in [(false, strands[0]), (true, strands[1])] {
                let band = Band::whole(sequence.len(), read.len());
                let scoring = &settings.scoring;
                let workspace = &mut scratch.workspace;
                let aligned = align(sequence, read, scoring, band, settings.ends, workspace);
                let Some(path) = aligned else { continue };
                let limits = self.limits(&None, settings, sequence.len());
                let beyond = sequence.len() - (path.subject_end - path.subject_start);
                if path.differences() > limits.differences || beyond > limits.beyond {
                    continue;
                }
                let nearness = Nearness {
                    cost: costs.of_path(&path, sequence.len()),
                    length: sequence.len() as i64,
                };
                let offered = Nearest {
                    nearness,
                    subject,
                    minus_strand,
                    band: (0, 0),
                    path: Some(path),
                };
                Self::offer(&mut nearest, offered);
            }
        }

        // Groups sharing words with either orientation. A read that
        // shares few words with every group is far from all of them:
        // every word of it is then looked up, not every `step`-th.
        let mut step = self.step;
        let mut forward_candidates;
        loop {
            scratch.candidates.clear();
            scratch
                .votes
                .collect(&self.index, strands[0], step, &mut scratch.candidates);
            forward_candidates = scratch.candidates.len();
            scratch
                .votes
                .collect(&self.index, strands[1], step, &mut scratch.candidates);
            let most = scratch
                .candidates
                .iter()
                .map(|c| c.votes)
                .max()
                .unwrap_or(0);
            let unchanged = self.shortest.saturating_sub(WORD + GROUP) / step;
            if step == 1 || 2 * most as usize >= unchanged {
                break;
            }
            step = 1;
        }
        // Most words first.
        let mut order = std::mem::take(&mut scratch.order);
        let candidates = std::mem::take(&mut scratch.candidates);
        order.clear();
        order.extend(0..candidates.len());
        order.sort_unstable_by(|&a, &b| {
            (candidates[b].votes.cmp(&candidates[a].votes))
                .then(candidates[a].subject.cmp(&candidates[b].subject))
                .then(a.cmp(&b))
        });
        let top_votes = order.first().map_or(0, |&i| candidates[i].votes);
        for &number in &order {
            let candidate = candidates[number];
            if candidate.votes < MIN_VOTES {
                break;
            }
            let minus_strand = number >= forward_candidates;
            let read = strands[usize::from(minus_strand)];
            self.judge_group(
                candidate,
                (top_votes, step),
                read,
                minus_strand,
                settings,
                scratch,
                &mut nearest,
            );
        }
        scratch.order = order;
        scratch.candidates = candidates;
        nearest
    }

    /// Judge the members of one group against one orientation of a read
    /// and keep the nearest in `nearest`.
    #[allow(clippy::too_many_arguments)]
    fn judge_group(
        &self,
        candidate: Candidate,
        (top_votes, step): (u32, usize),
        read: &[u8],
        minus_strand: bool,
        settings: &Settings,
        scratch: &mut Scratch,
        nearest: &mut Option<Nearest>,
    ) {
        let costs = Costs::of(&settings.scoring);
        let group = &self.groups[candidate.subject];
        let length = group.length;
        let (m, n) = (length as isize, read.len() as isize);
        let diagonal = candidate.diagonal;
        let limits = self.limits(nearest, settings, length);

        // Shared words: a difference removes at most `lost` of the
        // looked-up words, so few shared words mean many differences.
        let span = ((diagonal + m).min(n) - diagonal.max(0)).max(0) as usize;
        let ideal = span.saturating_sub(WORD + GROUP) / step;
        let expected = ideal.min(top_votes as usize);
        let lost = WORD.div_ceil(step);
        let votes = candidate.votes as usize;
        if expected.saturating_sub(votes).div_ceil(lost) > limits.differences {
            return;
        }
        // Diagonals that are tried, and those among them on which the
        // whole subject lies inside the read.
        // Around the diagonal with most words, and as far as words agree
        // on other diagonals (a read with an inserted or deleted stretch).
        let half = HALF_WIDTH as isize;
        let widest = MAX_HALF_WIDTH as isize;
        let (first_voted, last_voted) = candidate.spread;
        let band = (
            (diagonal - half)
                .min(first_voted - half / 2)
                .max(diagonal - widest),
            (diagonal + half)
                .max(last_voted + half / 2)
                .min(diagonal + widest),
        );
        let whole = (band.0.max(0), band.1.min(n - m));
        // Whether, on a tried diagonal, between one and `bases` bases of
        // the subject lie beyond an end of the read.
        let reaches_beyond = |bases: usize| {
            let bases = bases.min(length) as isize;
            bases > 0
                && ((band.0 < 0 && band.1 >= -bases) || (band.1 > n - m && band.0 <= n - m + bases))
        };

        // The first member tells for the whole group: no member is more
        // than `radius` nearer to the read.
        if group.members.len() > 1 {
            let leader = &self.subjects[group.members[0]].sequence;
            let within = limits.differences + group.radius;
            let beyond = if reaches_beyond(limits.beyond) {
                limits.beyond
            } else {
                0
            };
            let reach = &mut scratch.reach;
            if differences(leader, read, band.0, band.1, within, beyond, reach).is_none() {
                return;
            }
        }

        let offer = |nearest: &mut Option<Nearest>, subject: usize, cost: i64, path| {
            let nearness = Nearness {
                cost,
                length: length as i64,
            };
            let offered = Nearest {
                nearness,
                subject,
                minus_strand,
                band,
                path,
            };
            Self::offer(nearest, offered);
        };

        // Members that lie in the read with at most one mismatch. They can
        // do so only where the first member has few mismatches itself.
        let leader = &self.subjects[group.members[0]].sequence;
        let mut places = std::mem::take(&mut scratch.places);
        places.clear();
        places.extend((whole.0..=whole.1).filter(|&diagonal| {
            let start = diagonal as usize;
            mismatches(leader, &read[start..start + length], group.radius + 1).is_some()
        }));
        scratch.pending.clear();
        for &subject in &group.members {
            let sequence = &self.subjects[subject].sequence;
            match self.nearly_unchanged(sequence, read, &places) {
                Some((0, start)) => {
                    let path = gapless(sequence, read, start, &settings.scoring);
                    offer(nearest, subject, 0, Some(path));
                }
                // Leaving out a base or two beyond the read may be nearer
                // than the mismatch: the count below decides.
                Some((_, start))
                    if self.limits(nearest, settings, length).differences >= 1
                        && !reaches_beyond((costs.difference / costs.beyond) as usize) =>
                {
                    let path = gapless(sequence, read, start, &settings.scoring);
                    offer(nearest, subject, costs.difference, Some(path));
                }
                _ => scratch.pending.push(subject),
            }
        }

        scratch.places = places;

        // The others: the least differences, counted only as far as the
        // nearest subject so far allows.
        let pending = std::mem::take(&mut scratch.pending);
        for &subject in &pending {
            let sequence = &self.subjects[subject].sequence;
            let limits = self.limits(nearest, settings, length);
            let reach = &mut scratch.reach;
            // With the whole subject inside the read.
            let inside = (limits.differences > 0)
                .then(|| differences(sequence, read, band.0, band.1, limits.differences, 0, reach))
                .flatten();
            let inside_cost = inside.map(|count| count as i64 * costs.difference);
            // With some of it beyond the read: nearer only with fewer
            // differences, and with few enough bases beyond.
            let cost_below = inside_cost.unwrap_or(limits.cost.saturating_add(1));
            let beyond = limits
                .beyond
                .min(((cost_below - 1) / costs.beyond).max(0) as usize);
            let fewer = match inside {
                Some(0) => None,
                Some(count) => Some(count - 1),
                None => Some(limits.differences),
            };
            let outside = match fewer {
                Some(within) if reaches_beyond(beyond) => {
                    differences(sequence, read, band.0, band.1, within, beyond, reach)
                }
                _ => None,
            };
            // At least this many bases lie beyond the read on every tried
            // diagonal that is not whole.
            let surely_beyond = (-band.1).max(band.0 - (n - m)).max(1) as i64;
            let outside = outside.filter(|&count| {
                let least = count as i64 * costs.difference + surely_beyond * costs.beyond;
                let held_earlier = nearest.as_ref().is_some_and(|held| {
                    held.nearness.length == length as i64
                        && (held.subject, held.minus_strand) < (subject, minus_strand)
                });
                least < limits.cost || (least == limits.cost && !held_earlier)
            });
            if outside.is_some() {
                // How many bases lie beyond is known from an alignment.
                let path = align_in(band, sequence, read, settings, &mut scratch.workspace);
                if let Some(path) = path {
                    let left_out = length - (path.subject_end - path.subject_start);
                    let cost = costs.of_path(&path, length);
                    if path.differences() <= limits.differences
                        && left_out <= limits.beyond
                        && cost <= limits.cost
                        && inside_cost.is_none_or(|inside| cost < inside)
                    {
                        offer(nearest, subject, cost, Some(path));
                        continue;
                    }
                }
            }
            if let Some(cost) = inside_cost {
                offer(nearest, subject, cost, None);
            }
        }
        scratch.pending = pending;
    }

    /// Where `subject` lies whole in the read with at most one mismatch:
    /// the number of mismatches and the start in the read, on the lowest
    /// of the diagonals `places` with the fewest mismatches.
    fn nearly_unchanged(
        &self,
        subject: &[u8],
        read: &[u8],
        places: &[isize],
    ) -> Option<(usize, usize)> {
        let m = subject.len();
        let mut found: Option<(usize, usize)> = None;
        for &diagonal in places {
            let limit = match found {
                Some((0, _)) => break,
                Some((fewest, _)) => fewest - 1,
                None => 1,
            };
            let start = diagonal as usize;
            if let Some(count) = mismatches(subject, &read[start..start + m], limit) {
                found = Some((count, start));
            }
        }
        found
    }

    /// The best match of every read. Equal reads are aligned once; the
    /// work is spread over all processor cores, and the result does not
    /// depend on their number. Progress is shown on stderr.
    pub fn search(&self, reads: &[&str], settings: &Settings) -> Hits {
        let mut number_of: HashMap<&str, usize> = HashMap::new();
        let mut distinct: Vec<&str> = Vec::new();
        // How many reads hold each distinct sequence.
        let mut copies: Vec<usize> = Vec::new();
        let of_read: Vec<usize> = reads
            .iter()
            .map(|&read| {
                let number = *number_of.entry(read).or_insert_with(|| {
                    distinct.push(read);
                    copies.push(0);
                    distinct.len() - 1
                });
                copies[number] += 1;
                number
            })
            .collect();
        let progress = Progress::new("Aligning", "Aligned", "reads", reads.len());
        let hits: Vec<Option<Hit>> = distinct
            .par_iter()
            .zip(&copies)
            .map_init(Scratch::default, |scratch, (read, &copies)| {
                let hit = self.best_hit(read, settings, scratch);
                progress.add(copies);
                hit
            })
            .collect();
        progress.finish();
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

    /// The nearest subject by the plain dynamic programme over every
    /// subject and both orientations: (subject, minus strand, differences,
    /// bases beyond the read).
    fn reference(
        database: &Database,
        read: &str,
        settings: &Settings,
    ) -> (usize, bool, usize, usize) {
        let costs = Costs::of(&settings.scoring);
        let forward = read.as_bytes().to_vec();
        let reverse = reverse_complement(read).unwrap().into_bytes();
        let mut nearest: Option<(Nearness, (usize, bool, usize, usize))> = None;
        for (number, subject) in database.subjects().iter().enumerate() {
            for (minus_strand, oriented) in [(false, &forward), (true, &reverse)] {
                let (cost, beyond) = least_cost(&subject.sequence, oriented, &costs);
                let nearness = Nearness {
                    cost,
                    length: subject.sequence.len() as i64,
                };
                let differences =
                    ((cost - beyond as i64 * costs.beyond) / costs.difference) as usize;
                if nearest
                    .as_ref()
                    .is_none_or(|(held, _)| nearness.nearer_than(held))
                {
                    nearest = Some((nearness, (number, minus_strand, differences, beyond)));
                }
            }
        }
        nearest.unwrap().1
    }

    /// Least cost of laying `subject` on `read`, and the subject bases
    /// left beyond the read at that cost: every difference costs
    /// `costs.difference`, every base beyond an end of the read
    /// `costs.beyond`.
    fn least_cost(subject: &[u8], read: &[u8], costs: &Costs) -> (i64, usize) {
        let (m, n) = (subject.len(), read.len());
        // (cost, bases beyond the start of the read)
        let mut row: Vec<(i64, usize)> = vec![(0, 0); n + 1];
        let mut least = (i64::MAX, 0);
        for i in 1..=m {
            let mut next = vec![(costs.beyond * i as i64, i); n + 1];
            for j in 1..=n {
                let pair = i64::from(!is_same_base(subject[i - 1], read[j - 1])) * costs.difference;
                next[j] = [
                    (row[j - 1].0 + pair, row[j - 1].1),
                    (row[j].0 + costs.difference, row[j].1),
                    (next[j - 1].0 + costs.difference, next[j - 1].1),
                    next[j],
                ]
                .into_iter()
                .min()
                .unwrap();
            }
            row = next;
            // The rest of the subject beyond the end of the read.
            let rest = m - i;
            least = least.min((row[n].0 + costs.beyond * rest as i64, row[n].1 + rest));
        }
        least.min(*row.iter().min().unwrap())
    }

    /// Check that the columns of a hit say the truth about the two
    /// sequences.
    fn assert_consistent(hit: &Hit, database: &Database, read: &str) {
        let subject = &database.subject(hit.subject).sequence;
        let oriented = if hit.minus_strand {
            reverse_complement(read).unwrap()
        } else {
            read.to_string()
        };
        let oriented = oriented.as_bytes();
        let (mut i, mut j) = (hit.path.subject_start, hit.path.read_start);
        for op in &hit.path.ops {
            match op {
                Op::Match => assert!(is_same_base(subject[i], oriented[j])),
                Op::Mismatch => assert!(!is_same_base(subject[i], oriented[j])),
                _ => {}
            }
            i += usize::from(*op != Op::Insertion);
            j += usize::from(*op != Op::Deletion);
        }
        assert_eq!((i, j), (hit.path.subject_end, hit.path.read_end));
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

    /// Families of genes a tenth apart, each gene with alleles a few
    /// substitutions apart and one allele that lacks two bases: what
    /// makes real databases hard.
    fn families(dice: &mut Dice) -> Vec<(String, String)> {
        let mut alleles = Vec::new();
        for family in 1..=3 {
            let ancestor = dice.sequence(150);
            for gene in 1..=3 {
                let mut first = ancestor.clone();
                for _ in 0..15 {
                    let at = dice.below(first.len());
                    first = dice.substituted(&first, &[at]);
                }
                let mut sequences = vec![first.clone()];
                for _ in 0..5 {
                    let mut allele = first.clone();
                    for _ in 0..1 + dice.below(3) {
                        let at = dice.below(allele.len());
                        allele = dice.substituted(&allele, &[at]);
                    }
                    sequences.push(allele);
                }
                let mut shorter = first.clone();
                let at = 20 + dice.below(100);
                shorter.replace_range(at..at + 2, "");
                sequences.push(shorter);
                for (number, sequence) in sequences.into_iter().enumerate() {
                    alleles.push((format!("IGHV{family}-{gene}*0{}", number + 1), sequence));
                }
            }
        }
        alleles
    }

    #[test]
    fn the_nearest_allele_is_the_one_the_plain_programme_finds() {
        let mut dice = Dice(0x600D_F00D_0BAD_CAFE);
        let alleles = families(&mut dice);
        let database = Database::new(alleles.iter().map(|(n, s)| (n.as_str(), s))).unwrap();
        // The alleles of a gene share a group; the shorter one stands alone.
        assert_eq!(database.groups(), 18);
        let rounds: usize = std::env::var("IMMUNEDISCOVER_STRESS")
            .ok()
            .and_then(|text| text.parse().ok())
            .unwrap_or(40);
        let settings = Settings {
            min_identity: 0.8,
            min_coverage: 0.5,
            ..Settings::default()
        };
        let mut scratch = Scratch::default();
        let mut checked = 0;
        for round in 0..rounds {
            let (_, allele) = &alleles[dice.below(alleles.len())];
            let mut sequence = allele.clone();
            // Substitutions: none, a few, or so many that nothing matches.
            let changes = [0, 0, 1, 1, 2, 3, 5, 12, 30][dice.below(9)];
            for _ in 0..changes {
                let at = dice.below(sequence.len());
                sequence = dice.substituted(&sequence, &[at]);
            }
            // A base too many or too few, also at the very ends.
            match dice.below(8) {
                0 => {
                    let at = dice.below(sequence.len() + 1);
                    sequence.insert(at, b"ACGT"[dice.below(4)] as char);
                }
                1 => {
                    let at = dice.below(sequence.len() - 2);
                    let length = 1 + dice.below(2);
                    sequence.replace_range(at..at + length, "");
                }
                2 => {
                    let at = dice.below(sequence.len());
                    sequence.replace_range(at..at + 1, "N");
                }
                _ => {}
            }
            // The gene lies at the end of the read, near it, or far inside;
            // some reads stop inside the gene.
            let margins = [0, 1, 2, 5, 13, 40];
            let (left, right) = (dice.below(margins.len()), dice.below(margins.len()));
            let before = dice.sequence(margins[left]);
            let after = dice.sequence(margins[right]);
            let mut read = format!("{before}{sequence}{after}");
            match dice.below(10) {
                0 => read.truncate(read.len() - after.len() - 1 - dice.below(60)),
                1 => {
                    read.drain(..before.len() + 1 + dice.below(60));
                }
                _ => {}
            }
            if round % 2 == 0 {
                read = reverse_complement(&read).unwrap();
            }
            let found = database.best_hit(&read, &settings, &mut scratch);
            if let Some(hit) = &found {
                assert_consistent(hit, &database, &read);
            }
            // Far from every allele (more than one difference in ten
            // bases of what the read holds of it) few words are shared and
            // the nearest allele may be
            // missed; at the limits of `settings` rounding decides.
            // Well inside both, the search must name the nearest allele.
            let (subject, minus_strand, differences, beyond) =
                reference(&database, &read, &settings);
            let length = database.subject(subject).sequence.len();
            if 10 * differences > length - beyond || 5 * beyond > 2 * length {
                continue;
            }
            checked += 1;
            let hit = found.unwrap_or_else(|| panic!("round {round}: no match for {read}"));
            let name = |number: usize| database.subject(number).name.as_str();
            let context = format!(
                "round {round}: {differences} differences, {beyond} beyond, read {read}, found {} {}",
                name(hit.subject),
                hit.path.cigar()
            );
            if beyond == 0 {
                // The allele lies whole inside the read: the nearest one,
                // exactly.
                assert_eq!(
                    (name(hit.subject), hit.minus_strand),
                    (name(subject), minus_strand),
                    "{context}"
                );
            } else {
                // The read ends inside the allele. How much of the allele
                // lies beyond is taken from an alignment, which may leave
                // out a base more or less than the plain programme: the
                // allele found is at most one difference further.
                let costs = Costs::of(&settings.scoring);
                let oriented = if hit.minus_strand {
                    reverse_complement(&read).unwrap()
                } else {
                    read.clone()
                };
                let found = &database.subject(hit.subject).sequence;
                let (cost, _) = least_cost(found, oriented.as_bytes(), &costs);
                let nearest =
                    (differences as i64 + 1) * costs.difference + beyond as i64 * costs.beyond;
                assert!(
                    cost * length as i64 <= nearest * found.len() as i64,
                    "{context}"
                );
            }
        }
        assert!(checked * 2 > rounds);
    }

    #[test]
    fn an_allele_that_lies_unchanged_in_the_read_is_named_at_once() {
        let mut dice = Dice(0x0DDB_A11A_D5EE_D5ED);
        let gene = dice.sequence(120);
        let part = gene[10..90].to_string();
        let other = dice.sequence(120);
        let with_n = format!("{}N{}", &other[..60], &other[61..]);
        let database = Database::new([
            ("part", part.as_str()),
            ("gene*01", gene.as_str()),
            ("gene*02", gene.as_str()),
            ("other", with_n.as_str()),
        ])
        .unwrap();
        let scoring = Scoring::default();
        let settings = Settings::default();
        let mut scratch = Scratch::default();

        // The longest allele that fits wins, then the first of the database.
        let read = format!("{}{gene}{}", dice.sequence(33), dice.sequence(21));
        let hit = database
            .unchanged_subject(read.as_bytes(), &scoring)
            .unwrap();
        assert_eq!(database.subject(hit.subject).name, "gene*01");
        assert_eq!((hit.read_range(), hit.minus_strand), ((33, 153), false));
        assert_eq!(hit.path.cigar(), "120=");
        assert_eq!(database.best_hit(&read, &settings, &mut scratch), Some(hit));

        // On the other strand, with the position in the read as written.
        let turned = reverse_complement(&read).unwrap();
        let hit = database
            .unchanged_subject(turned.as_bytes(), &scoring)
            .unwrap();
        assert_eq!(database.subject(hit.subject).name, "gene*01");
        assert_eq!((hit.read_range(), hit.minus_strand), ((21, 141), true));
        assert_eq!(hit.segment(&turned), gene);

        // Only the shorter allele lies in this read unchanged.
        let read = format!("{}{part}{}", dice.sequence(15), dice.sequence(15));
        let hit = database
            .unchanged_subject(read.as_bytes(), &scoring)
            .unwrap();
        assert_eq!(database.subject(hit.subject).name, "part");

        // A read with a difference is not named here, and not lost: the
        // search that follows finds its allele.
        let changed = dice.substituted(&gene, &[50]);
        let read = format!("{}{changed}{}", dice.sequence(33), dice.sequence(21));
        assert!(database
            .unchanged_subject(read.as_bytes(), &scoring)
            .is_none());
        let hit = database.best_hit(&read, &settings, &mut scratch).unwrap();
        assert_eq!(database.subject(hit.subject).name, "gene*01");
        assert_eq!(hit.path.cigar(), "50=1X69=");

        // An allele with an ambiguous base never counts as unchanged.
        let read = format!("{}{with_n}{}", dice.sequence(20), dice.sequence(20));
        assert!(database
            .unchanged_subject(read.as_bytes(), &scoring)
            .is_none());
        let hit = database.best_hit(&read, &settings, &mut scratch).unwrap();
        assert_eq!(database.subject(hit.subject).name, "other");
        assert_eq!(hit.path.mismatches(), 1);
    }

    #[test]
    fn a_read_with_a_long_inserted_stretch_still_finds_its_allele() {
        let mut dice = Dice(0x05EE_D0F1_DEA5);
        let alleles = alleles(&mut dice);
        let database = Database::new(alleles.iter().map(|(n, s)| (n.as_str(), s))).unwrap();
        // 20 bases that the allele lacks, in its middle: more than the
        // band around one diagonal holds.
        let allele = &alleles[4].1;
        let inserted = format!("{}{}{}", &allele[..150], dice.sequence(20), &allele[150..]);
        let read = read_of(&mut dice, &inserted);
        let settings = Settings::default();
        let hit = database
            .best_hit(&read, &settings, &mut Scratch::default())
            .unwrap();
        assert_eq!(hit.subject, 4);
        assert_eq!(hit.path.cigar(), "150=20I150=");
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
