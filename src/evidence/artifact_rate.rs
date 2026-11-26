//! Among the sequences like this one, which share is expected to be artifacts?
//!
//! A real allele gives about one or two copies' worth of reads; an artifact
//! usually gives far less. The question for a candidate is how often
//! artifacts of its [`Kind`] reach its dosage on this plate, and that is
//! measured on the plate itself:
//!
//! 1. **Where artifacts are alone.** Below [`ARTIFACT_ONLY`] of a copy the
//!    sequences of a plate are artifacts, a few exceptions aside. Those
//!    private to one donor are counted there in two windows of half a
//!    decade of dosage each.
//! 2. **How fast they thin out.** The decline from the lower window to the
//!    upper one is continued upwards as a power-law tail, the usual model
//!    for the far end of a distribution. It gives the number of artifacts
//!    expected in any range of dosage above the windows.
//! 3. **Against what was seen.** The artifacts expected within half a
//!    decade around a dosage, divided by the candidates observed there, is
//!    the expected share of artifacts among them (a local false discovery
//!    rate). The evidence is its -log10: 1 means one in ten, 2 one in a
//!    hundred. More copies are never weaker evidence: a sequence takes
//!    the best rate found at or below its dosage.
//!
//! There is no dosage below which a sequence is declared an artifact: a
//! real allele that amplifies poorly, at half a copy, is judged by how many
//! artifacts the plate has at half a copy.
//!
//! Sporadic artifacts arise in one donor at a time, so the tail is learned
//! from sequences that no other donor carries and that most other donors
//! do not show at all, and a candidate private to its donor is compared
//! with the private candidates only. A sequence other donors carry is
//! compared with all candidates, which bounds its rate from above.
//! Sequences at the same level in every donor are systematic, not
//! sporadic: they are left out of the counts and are told apart by
//! `ubiquitous`, not by this evidence.
//!
//! The extrapolation is kept from running away: artifacts are taken to
//! thin out by at least half and at most twenty-fold per half decade of
//! dosage, half an artifact is added to every count, and a kind with too
//! few artifacts to measure its own decline uses the decline of all kinds
//! together.

use super::kinds::Kind;
use super::plate::is_private;
use super::Observation;

/// Reads one copy must be expected to give for an observation to be
/// assessed. Below it, a few stray reads are a copy-like dosage.
pub const MIN_COPY_READS: f64 = 20.0;
/// Dosage below which the sequences of a plate are taken to be artifacts
/// when their tail is measured.
pub const ARTIFACT_ONLY: f64 = 0.2;
/// Artifacts a kind needs in the two windows for its own decline.
const MIN_TAIL_ROWS: usize = 20;
/// Half a decade of dosage: the width of a window.
const STEP: f64 = 3.162_277_660_168_379_5;
/// Bounds on the decline per step.
const DECLINE: (f64, f64) = (0.05, 0.5);
const PSEUDO_COUNT: f64 = 0.5;
/// Evidence is not reported beyond this.
const MAX_EVIDENCE: f64 = 9.0;

/// Decline per window, from the lower window to the upper one.
fn decline(windows: [usize; 2]) -> f64 {
    let [lower, upper] = windows.map(|n| n as f64 + PSEUDO_COUNT);
    (upper / lower).clamp(DECLINE.0, DECLINE.1)
}

/// The artifacts of one kind below [`ARTIFACT_ONLY`], and their tail above.
#[derive(Debug, Clone)]
pub struct ArtifactTail {
    /// Artifacts in the two windows below [`ARTIFACT_ONLY`], lower first.
    pub windows: [usize; 2],
    /// Factor by which artifacts thin out per half decade of dosage.
    pub decline: f64,
    /// Whether the decline was measured on this kind alone.
    pub own_decline: bool,
}

impl ArtifactTail {
    /// Artifacts of this kind expected on the plate with at least `dosage`.
    pub fn expected(&self, dosage: f64) -> f64 {
        let exponent = -self.decline.ln() / STEP.ln();
        (self.windows[1] as f64 + PSEUDO_COUNT) * self.decline / (1.0 - self.decline)
            * (dosage / ARTIFACT_ONLY).max(1.0).powf(-exponent)
    }

    /// Artifacts expected between two dosages.
    fn expected_within(&self, low: f64, high: f64) -> f64 {
        self.expected(low) - self.expected(high)
    }
}

/// Half a decade around a dosage, not reaching into the counting windows.
fn neighbourhood(dosage: f64) -> (f64, f64) {
    let half = STEP.sqrt();
    ((dosage / half).max(ARTIFACT_ONLY), dosage * half)
}

/// The candidates of one kind above [`ARTIFACT_ONLY`], with the expected
/// share of artifacts at every one of them.
#[derive(Debug, Clone, Default)]
struct Candidates {
    /// Ascending.
    dosages: Vec<f64>,
    /// Lowest share of artifacts found around this or any weaker candidate.
    rate: Vec<f64>,
}

impl Candidates {
    fn new(mut dosages: Vec<f64>, tail: &ArtifactTail) -> Self {
        dosages.sort_by(f64::total_cmp);
        let mut candidates = Candidates {
            dosages,
            rate: Vec::new(),
        };
        let mut lowest = 1.0f64;
        for i in 0..candidates.dosages.len() {
            lowest = lowest.min(candidates.local_rate(candidates.dosages[i], tail));
            candidates.rate.push(lowest);
        }
        candidates
    }

    /// Artifacts expected around a dosage over candidates seen there. A
    /// sequence not among the candidates counts as one.
    fn local_rate(&self, dosage: f64, tail: &ArtifactTail) -> f64 {
        let (low, high) = neighbourhood(dosage);
        let seen = self.dosages.partition_point(|d| *d <= high)
            - self.dosages.partition_point(|d| *d < low);
        (tail.expected_within(low, high) / seen.max(1) as f64).min(1.0)
    }

    fn rate(&self, dosage: f64, tail: &ArtifactTail) -> f64 {
        let here = self.local_rate(dosage, tail);
        match self.dosages.partition_point(|d| *d < dosage) {
            0 => here,
            weaker => here.min(self.rate[weaker - 1]),
        }
    }
}

/// The fitted tails and candidates of every kind.
#[derive(Debug, Clone)]
pub struct ArtifactRate {
    tails: Vec<ArtifactTail>,
    private: Vec<Candidates>,
    all: Vec<Candidates>,
}

/// Seen in most of the other donors holding its parent: whatever it is,
/// it did not arise by chance in this donor.
fn recurs(o: &Observation) -> bool {
    o.evidence.recurrent_ratio.is_some_and(|r| r > 0.0)
}

/// Dosage of an observation, if the donor is deep enough for the gene.
fn assessed_dosage(o: &Observation) -> Option<f64> {
    o.evidence
        .per_copy_reads
        .filter(|reads| *reads >= MIN_COPY_READS)
        .and(o.evidence.dosage)
}

impl ArtifactRate {
    pub fn fit(observations: &[Observation]) -> Self {
        // Sequences that could be sporadic artifacts or unknown alleles.
        let open: Vec<(&Observation, Kind, bool, f64)> = observations
            .iter()
            .filter(|o| !o.is_database_allele() && !o.evidence.ubiquitous)
            .filter_map(|o| Some((o, Kind::of(o), is_private(o), assessed_dosage(o)?)))
            .collect();
        let count = |kind: Option<Kind>, low: f64, high: f64| {
            open.iter()
                .filter(|(o, k, private, dosage)| {
                    *private
                        && !recurs(o)
                        && kind.is_none_or(|kind| kind == *k)
                        && (low..high).contains(dosage)
                })
                .count()
        };
        let windows = |kind| {
            [2, 1].map(|i| {
                count(
                    kind,
                    ARTIFACT_ONLY / STEP.powi(i),
                    ARTIFACT_ONLY / STEP.powi(i - 1),
                )
            })
        };
        let pooled = decline(windows(None));
        let tails: Vec<ArtifactTail> = Kind::ALL
            .iter()
            .map(|kind| {
                let windows = windows(Some(*kind));
                let own_decline = windows.iter().sum::<usize>() >= MIN_TAIL_ROWS;
                ArtifactTail {
                    windows,
                    decline: if own_decline {
                        decline(windows)
                    } else {
                        pooled
                    },
                    own_decline,
                }
            })
            .collect();
        let candidates = |only_private: bool| -> Vec<Candidates> {
            Kind::ALL
                .iter()
                .zip(&tails)
                .map(|(kind, tail)| {
                    let dosages = open
                        .iter()
                        .filter(|(_, k, private, dosage)| {
                            k == kind && *dosage >= ARTIFACT_ONLY && (*private || !only_private)
                        })
                        .map(|(_, _, _, dosage)| *dosage)
                        .collect();
                    Candidates::new(dosages, tail)
                })
                .collect()
        };
        ArtifactRate {
            private: candidates(true),
            all: candidates(false),
            tails,
        }
    }

    pub fn tail(&self, kind: Kind) -> &ArtifactTail {
        &self.tails[kind.index()]
    }

    /// Candidates above [`ARTIFACT_ONLY`]: private to their donor, and all.
    pub fn candidates(&self, kind: Kind) -> (usize, usize) {
        (
            self.private[kind.index()].dosages.len(),
            self.all[kind.index()].dosages.len(),
        )
    }

    /// Expected share of artifacts among the sequences of this kind (and,
    /// for a private one, private like it) around this dosage.
    pub fn rate(&self, kind: Kind, private: bool, dosage: f64) -> f64 {
        if dosage < ARTIFACT_ONLY {
            return 1.0;
        }
        let candidates = if private { &self.private } else { &self.all };
        candidates[kind.index()].rate(dosage, self.tail(kind))
    }

    /// The evidence: -log10 of [`ArtifactRate::rate`].
    pub fn evidence(&self, kind: Kind, private: bool, dosage: f64) -> f64 {
        // Adding zero turns the -0.0 of a rate of one into 0.0.
        (-self.rate(kind, private, dosage).log10()).clamp(0.0, MAX_EVIDENCE) + 0.0
    }

    /// Fill in `score` for every observation deep enough to be assessed.
    pub fn set_scores(&self, observations: &mut [Observation]) {
        for o in observations.iter_mut() {
            let (kind, private) = (Kind::of(o), is_private(o));
            o.evidence.score = assessed_dosage(o).map(|d| self.evidence(kind, private, d));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::testing::observation;

    /// A sequence one difference from a more abundant one.
    fn variant(dosage: f64, carriers_elsewhere: usize) -> Observation {
        let mut o = observation(0, "A", "variant", (dosage * 100.0).round() as u64);
        o.db_distance = 1;
        o.evidence.per_copy_reads = Some(100.0);
        o.evidence.dosage = Some(dosage);
        o.evidence.parent = Some((0, 1));
        o.evidence.carriers = carriers_elsewhere + usize::from(dosage >= 0.25);
        o
    }

    /// 1000 and 100 artifacts in the two windows below a fifth of a copy;
    /// three private artifacts just above it; five private candidates and
    /// fifty shared ones at one copy.
    fn plate() -> Vec<Observation> {
        let mut observations = Vec::new();
        observations.extend((0..1000).map(|_| variant(0.03, 0)));
        observations.extend((0..100).map(|_| variant(0.1, 0)));
        observations.extend((0..3).map(|_| variant(0.22, 0)));
        observations.extend((0..5).map(|_| variant(1.0, 0)));
        observations.extend((0..50).map(|_| variant(1.0, 3)));
        observations
    }

    #[test]
    fn tail_continues_the_decline_between_the_windows() {
        let model = ArtifactRate::fit(&plate());
        let tail = model.tail(Kind::Near(1));
        assert_eq!(tail.windows, [1000, 100]);
        assert!(tail.own_decline);
        let decline = 100.5 / 1000.5;
        assert!((tail.decline - decline).abs() < 1e-12);
        // Above the windows: the next half decade holds `decline` times the
        // upper one, the one after it `decline` squared, and so on.
        let beyond = 100.5 * decline / (1.0 - decline);
        assert!((tail.expected(ARTIFACT_ONLY) - beyond).abs() < 1e-9);
        assert!((tail.expected(ARTIFACT_ONLY * STEP) - beyond * decline).abs() < 1e-9);
        let next = tail.expected_within(ARTIFACT_ONLY, ARTIFACT_ONLY * STEP);
        assert!((next - 100.5 * decline).abs() < 1e-9);
        assert_eq!(model.candidates(Kind::Near(1)), (8, 58));
    }

    #[test]
    fn rate_is_expected_artifacts_over_candidates_seen_around_a_dosage() {
        let mut observations = plate();
        let model = ArtifactRate::fit(&observations);
        let near = Kind::Near(1);
        let tail = model.tail(near);
        // Just above the windows ten artifacts are expected and three
        // sequences seen: nothing can be said for them.
        assert_eq!(model.rate(near, true, 0.22), 1.0);
        assert_eq!(model.evidence(near, true, 0.22), 0.0);
        // Around one copy: private candidates are compared with each other,
        // shared ones with all candidates.
        let (low, high) = neighbourhood(1.0);
        let around = tail.expected_within(low, high);
        assert!((1.2..1.4).contains(&around), "{around}");
        assert!((model.rate(near, true, 1.0) - around / 5.0).abs() < 1e-12);
        assert!((model.rate(near, false, 1.0) - around / 55.0).abs() < 1e-12);
        // Inside the windows everything counts as an artifact.
        assert_eq!(model.rate(near, false, 0.19), 1.0);
        // More copies are never weaker evidence, also where no candidate is.
        assert!(model.rate(near, true, 5.0) < model.rate(near, true, 1.0));
        assert!(model.rate(near, true, 0.5) > model.rate(near, true, 1.0));
        model.set_scores(&mut observations);
        let shared = observations.last().unwrap().evidence.score.unwrap();
        assert!((shared - -(around / 55.0).log10()).abs() < 1e-12);
        assert_eq!(observations[0].evidence.score, Some(0.0));
    }

    #[test]
    fn half_a_copy_is_judged_by_the_artifacts_at_half_a_copy() {
        // No artifact anywhere near: a lone sequence at half a copy stands.
        let mut observations: Vec<Observation> = (0..300).map(|_| variant(0.03, 0)).collect();
        observations.push(variant(0.5, 0));
        let clean = ArtifactRate::fit(&observations);
        assert!(clean.evidence(Kind::Near(1), true, 0.5) > 1.5);
        // Many artifacts right below it: it does not.
        observations.extend((0..150).map(|_| variant(0.15, 0)));
        let noisy = ArtifactRate::fit(&observations);
        assert_eq!(noisy.evidence(Kind::Near(1), true, 0.5), 0.0);
    }

    #[test]
    fn rare_kinds_borrow_the_decline_and_shallow_cells_are_skipped() {
        let mut observations = plate();
        // A lone sequence far from everything, at one copy in one donor.
        let mut far = variant(1.0, 0);
        far.evidence.parent = Some((0, 6));
        observations.push(far);
        // Two reads where a copy gives four.
        let mut shallow = variant(0.5, 0);
        shallow.evidence.per_copy_reads = Some(4.0);
        observations.push(shallow);
        let model = ArtifactRate::fit(&observations);
        let tail = model.tail(Kind::Far);
        assert_eq!((tail.windows, tail.own_decline), ([0; 2], false));
        assert_eq!(tail.decline, model.tail(Kind::Near(1)).decline);
        // Half a pseudo-artifact in the upper window, thinning out as elsewhere.
        let rate = model.rate(Kind::Far, true, 1.0);
        assert!((0.001..0.02).contains(&rate), "{rate}");
        model.set_scores(&mut observations);
        assert_eq!(observations.last().unwrap().evidence.score, None);
    }

    #[test]
    fn decline_stays_within_bounds() {
        assert_eq!(decline([10, 10]), 0.5);
        assert_eq!(decline([100_000, 3]), 0.05);
        assert!((decline([1000, 100]) - 100.5 / 1000.5).abs() < 1e-12);
    }

    #[test]
    fn systematic_and_recurrent_sequences_are_left_out_of_the_tail() {
        let mut observations = plate();
        let mut systematic: Vec<Observation> = (0..500).map(|_| variant(0.1, 40)).collect();
        for o in &mut systematic {
            o.evidence.ubiquitous = true;
        }
        observations.extend(systematic);
        // An unusually strong copy of a variant most donors show faintly.
        let mut recurrent = variant(0.1, 0);
        recurrent.evidence.recurrent_ratio = Some(0.02);
        observations.push(recurrent);
        let model = ArtifactRate::fit(&observations);
        assert_eq!(model.tail(Kind::Near(1)).windows, [1000, 100]);
    }
}
