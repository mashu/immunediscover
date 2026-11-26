//! Among the sequences like this one, which share is expected to be artifacts?
//!
//! A real allele gives one or two copies' worth of reads; an artifact
//! usually gives far less. The question for a candidate is how often
//! artifacts of its [`Kind`] reach its dosage on this plate, and that is
//! measured on the plate itself:
//!
//! 1. **Where artifacts are alone.** One copy varies around a dosage of 1
//!    by the allele-level spread of the plate; 99 % of single copies lie
//!    above the *edge* (see [`copy_edge`]). A sequence below the edge is
//!    not a copy, whatever else it is.
//! 2. **How fast they thin out.** The artifacts of a kind below the edge
//!    are counted in three windows of half a decade of dosage each. Their
//!    decline from window to window gives a power-law tail, the usual
//!    model for the far end of a distribution, and it extrapolates the
//!    number of artifacts expected above any dosage beyond the edge. Two
//!    declines are measured, over the last window and on average over the
//!    last two, and the slower one is used: abundant faint errors in the
//!    lowest window must not hide a rarer kind of artifact with a longer
//!    tail.
//! 3. **Against what was seen.** The expected number of artifacts with at
//!    least a dosage, divided by the number of candidates observed with at
//!    least that dosage, is the expected share of artifacts among them (a
//!    false discovery rate). The evidence is its -log10: 1 means one in
//!    ten, 2 one in a hundred.
//!
//! Sporadic artifacts arise in one donor at a time, so the tail is learned
//! from sequences that no other donor carries and that most other donors
//! do not show at all, and a candidate private to its donor is compared
//! with the private candidates only. A sequence other
//! donors carry is compared with all candidates, which bounds its rate
//! from above. Sequences at the same level in every donor are systematic,
//! not sporadic: they are left out of the counts and are told apart by
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
/// Artifacts a kind needs below the edge for its own decline.
const MIN_TAIL_ROWS: usize = 20;
/// Half a decade of dosage: the width of the three counting windows.
const STEP: f64 = 3.162_277_660_168_379_5;
/// Bounds on the decline per step.
const DECLINE: (f64, f64) = (0.05, 0.5);
const PSEUDO_COUNT: f64 = 0.5;
/// Artifacts' worth of weight given to the average decline when the
/// decline over the last window is measured.
const SHRINKAGE: f64 = 5.0;
/// Standard normal deviate of the lowest 1 %.
const LOWEST_COPIES: f64 = 2.326;
/// Bounds on the edge, in copies.
const EDGE: (f64, f64) = (0.3, 0.75);
/// Evidence is not reported beyond this.
const MAX_EVIDENCE: f64 = 9.0;

/// The dosage below which a sequence is not a copy: the first percentile of
/// one copy, given the allele-level spread (SD of the log dosage).
pub fn copy_edge(allele_noise: f64) -> f64 {
    (-LOWEST_COPIES * allele_noise).exp().clamp(EDGE.0, EDGE.1)
}

/// Decline per window: the slower of the average decline over two windows
/// and the decline over the last one. With few artifacts the latter is
/// mostly noise, and is drawn towards the former.
fn decline(windows: [usize; 3]) -> f64 {
    let [lowest, middle, highest] = windows.map(|n| n as f64);
    let average = ((highest + PSEUDO_COUNT) / (lowest + PSEUDO_COUNT)).sqrt();
    let last = (highest + SHRINKAGE * average) / (middle + SHRINKAGE);
    average.max(last).clamp(DECLINE.0, DECLINE.1)
}

/// The artifacts of one kind below the edge, and their tail beyond it.
#[derive(Debug, Clone)]
pub struct ArtifactTail {
    /// Artifacts in the three windows below the edge, lowest dosage first.
    pub windows: [usize; 3],
    /// Factor by which artifacts thin out per half decade of dosage.
    pub decline: f64,
    /// Whether the decline was measured on this kind alone.
    pub own_decline: bool,
    edge: f64,
}

impl ArtifactTail {
    /// Artifacts of this kind expected on the plate with at least `dosage`.
    pub fn expected(&self, dosage: f64) -> f64 {
        let exponent = -self.decline.ln() / STEP.ln();
        (self.windows[2] as f64 + PSEUDO_COUNT) * self.decline / (1.0 - self.decline)
            * (dosage / self.edge).max(1.0).powf(-exponent)
    }
}

/// The candidates of one kind above the edge, with the expected share of
/// artifacts at every one of them.
#[derive(Debug, Clone, Default)]
struct Candidates {
    /// Ascending.
    dosages: Vec<f64>,
    /// Lowest share of artifacts among the candidates holding at least the
    /// dosage of this or any weaker candidate.
    rate: Vec<f64>,
    /// The same for all candidates above the edge.
    rate_at_edge: f64,
}

impl Candidates {
    fn new(mut dosages: Vec<f64>, tail: &ArtifactTail) -> Self {
        dosages.sort_by(f64::total_cmp);
        let n = dosages.len();
        let rate_at_edge = (tail.expected(tail.edge) / n.max(1) as f64).min(1.0);
        let mut lowest = rate_at_edge;
        let rate = dosages
            .iter()
            .enumerate()
            .map(|(i, dosage)| {
                lowest = lowest.min(tail.expected(*dosage) / (n - i) as f64);
                lowest
            })
            .collect();
        Candidates {
            dosages,
            rate,
            rate_at_edge,
        }
    }

    /// For any dosage above the edge, also of a sequence not among the
    /// candidates (it then counts as one more).
    fn rate(&self, dosage: f64, tail: &ArtifactTail) -> f64 {
        let weaker = self.dosages.partition_point(|d| *d < dosage);
        let at_least = (self.dosages.len() - weaker).max(1);
        let before = match weaker {
            0 => self.rate_at_edge,
            i => self.rate[i - 1],
        };
        before.min(tail.expected(dosage) / at_least as f64)
    }
}

/// The fitted tails and candidates of every kind.
#[derive(Debug, Clone)]
pub struct ArtifactRate {
    /// Dosage below which a sequence is not a copy.
    pub edge: f64,
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
    pub fn fit(observations: &[Observation], allele_noise: f64) -> Self {
        let edge = copy_edge(allele_noise);
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
        let windows =
            |kind| [3, 2, 1].map(|i| count(kind, edge / STEP.powi(i), edge / STEP.powi(i - 1)));
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
                    edge,
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
                            k == kind && *dosage >= edge && (*private || !only_private)
                        })
                        .map(|(_, _, _, dosage)| *dosage)
                        .collect();
                    Candidates::new(dosages, tail)
                })
                .collect()
        };
        ArtifactRate {
            edge,
            private: candidates(true),
            all: candidates(false),
            tails,
        }
    }

    pub fn tail(&self, kind: Kind) -> &ArtifactTail {
        &self.tails[kind.index()]
    }

    /// Candidates above the edge: private to their donor, and all.
    pub fn candidates(&self, kind: Kind) -> (usize, usize) {
        (
            self.private[kind.index()].dosages.len(),
            self.all[kind.index()].dosages.len(),
        )
    }

    /// Expected share of artifacts among the sequences of this kind (and,
    /// for a private one, private like it) with at least this dosage.
    pub fn rate(&self, kind: Kind, private: bool, dosage: f64) -> f64 {
        if dosage < self.edge {
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

    /// Allele-level spread 0.2: the edge is exp(-2.326 * 0.2) = 0.628.
    /// 1000, 100 and 20 artifacts in the three windows below it; five
    /// private candidates and fifty shared ones at one copy.
    fn plate() -> Vec<Observation> {
        let mut observations = Vec::new();
        observations.extend((0..1000).map(|_| variant(0.03, 0)));
        observations.extend((0..100).map(|_| variant(0.1, 0)));
        observations.extend((0..20).map(|_| variant(0.4, 0)));
        observations.extend((0..5).map(|_| variant(1.0, 0)));
        observations.extend((0..50).map(|_| variant(1.0, 3)));
        observations
    }

    #[test]
    fn edge_follows_the_spread_of_one_copy() {
        assert!((copy_edge(0.2) - 0.628).abs() < 1e-3);
        assert_eq!(copy_edge(0.05), 0.75);
        assert_eq!(copy_edge(1.0), 0.3);
    }

    #[test]
    fn tail_is_extrapolated_from_the_windows_below_the_edge() {
        let model = ArtifactRate::fit(&plate(), 0.2);
        let tail = model.tail(Kind::Near(1));
        assert_eq!(tail.windows, [1000, 100, 20]);
        assert!(tail.own_decline);
        // Over two windows the decline averages sqrt(20.5 / 1000.5) = 0.143;
        // over the last one it is slower, and that one counts.
        let average = (20.5f64 / 1000.5).sqrt();
        let decline = (20.0 + 5.0 * average) / 105.0;
        assert!(decline > average);
        assert!((tail.decline - decline).abs() < 1e-12);
        // Beyond the edge: the next half decade holds `decline` times the
        // last one, the one after it `decline` squared, and so on.
        let beyond = 20.5 * decline / (1.0 - decline);
        assert!((tail.expected(model.edge) - beyond).abs() < 1e-9);
        assert!((tail.expected(model.edge * STEP) - beyond * decline).abs() < 1e-9);
        assert_eq!(model.candidates(Kind::Near(1)), (5, 55));
    }

    #[test]
    fn rate_is_expected_artifacts_over_candidates_seen() {
        let mut observations = plate();
        let model = ArtifactRate::fit(&observations, 0.2);
        let near = Kind::Near(1);
        let expected = model.tail(near).expected(1.0);
        assert!((expected - 2.615).abs() < 0.01, "{expected}");
        // Private candidates are compared with each other, shared ones
        // with all candidates.
        assert!((model.rate(near, true, 1.0) - expected / 5.0).abs() < 1e-12);
        assert!((model.rate(near, false, 1.0) - expected / 55.0).abs() < 1e-12);
        // Below the edge nothing is a copy.
        assert_eq!(model.rate(near, false, 0.5), 1.0);
        assert_eq!(model.evidence(near, false, 0.5), 0.0);
        // More copies are never weaker evidence, also where no candidate is.
        assert!(model.rate(near, true, 2.0) <= model.rate(near, true, 1.0));
        assert!(model.rate(near, true, 5.0) < model.rate(near, true, 1.0));
        assert!(model.rate(near, true, 0.7) >= model.rate(near, true, 1.0));
        model.set_scores(&mut observations);
        let shared = observations.last().unwrap().evidence.score.unwrap();
        assert!((shared - -(expected / 55.0).log10()).abs() < 1e-12);
        assert_eq!(observations[0].evidence.score, Some(0.0));
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
        let model = ArtifactRate::fit(&observations, 0.2);
        let tail = model.tail(Kind::Far);
        assert_eq!((tail.windows, tail.own_decline), ([0; 3], false));
        assert_eq!(tail.decline, model.tail(Kind::Near(1)).decline);
        // Half a pseudo-artifact below the edge, thinning out as elsewhere.
        let rate = model.rate(Kind::Far, true, 1.0);
        assert!((0.03..0.08).contains(&rate), "{rate}");
        model.set_scores(&mut observations);
        assert_eq!(observations.last().unwrap().evidence.score, None);
    }

    #[test]
    fn decline_is_the_slower_of_two() {
        // Steep on average, nothing in the last windows: the average.
        assert!((decline([400, 0, 0]) - (0.5f64 / 400.5).sqrt().max(0.05)).abs() < 1e-12);
        // Faint errors below, a long tail above them: the last window.
        let long_tail = decline([1000, 80, 20]);
        assert!((long_tail - (20.0 + 5.0 * (20.5f64 / 1000.5).sqrt()) / 85.0).abs() < 1e-12);
        // Never slower than by half, never faster than twenty-fold.
        assert_eq!(decline([10, 10, 10]), 0.5);
        assert_eq!(decline([100_000, 3, 0]), 0.05);
    }

    #[test]
    fn systematic_and_recurrent_sequences_are_left_out_of_the_tail() {
        let mut observations = plate();
        let mut systematic: Vec<Observation> = (0..500).map(|_| variant(0.4, 40)).collect();
        for o in &mut systematic {
            o.evidence.ubiquitous = true;
        }
        observations.extend(systematic);
        // An unusually strong copy of a variant most donors show faintly.
        let mut recurrent = variant(0.4, 0);
        recurrent.evidence.recurrent_ratio = Some(0.02);
        observations.push(recurrent);
        let model = ArtifactRate::fit(&observations, 0.2);
        assert_eq!(model.tail(Kind::Near(1)).windows, [1000, 100, 20]);
    }
}
