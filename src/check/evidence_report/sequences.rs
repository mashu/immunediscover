//! Every sequence judged once, on the donor where its evidence is strongest.

use std::cmp::Ordering;
use std::collections::{BTreeMap, HashSet};
use std::fmt::Write as _;

use super::{number, optional, Class, Plate};
use crate::evidence::plate::{COPY_LIKE, MULTI_COPY};
use crate::evidence::{Kind, Observation};
use crate::stats::median;
use crate::stats::robust::auc;

/// Evidence from which a donor counts as supporting a sequence.
const SUPPORT: f64 = 1.0;
const CUTS: [f64; 6] = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0];
const MAX_COST_LINES: usize = 10;

pub(super) struct Sequence<'a> {
    pub class: Class,
    /// The donor where the evidence is strongest.
    pub best: &'a Observation,
    pub max_count: u64,
    pub max_ratio: f64,
    /// Donors showing the sequence at all.
    pub donors: usize,
    /// Donors with evidence of at least [`SUPPORT`].
    pub supporting: usize,
    /// Strongest evidence among the donors passing the count thresholds;
    /// `None` when no donor passes them. A donor too shallow to be assessed
    /// counts as infinite: nothing speaks against the sequence there.
    pub passing_evidence: Option<f64>,
    /// Median dosage over the donors carrying the sequence.
    pub typical_dosage: f64,
}

impl Sequence<'_> {
    fn evidence(&self) -> Option<f64> {
        self.best.evidence.score
    }

    fn at_least(&self, cut: f64) -> bool {
        self.evidence().is_some_and(|e| e >= cut)
    }

    fn passes(&self) -> bool {
        self.passing_evidence.is_some()
    }

    fn ubiquitous(&self) -> bool {
        self.best.evidence.ubiquitous
    }

    fn multi_copy(&self) -> bool {
        self.typical_dosage > MULTI_COPY
    }

    fn note(&self) -> &'static str {
        match (self.ubiquitous(), self.multi_copy()) {
            (true, true) => "does not segregate, several copies",
            (true, false) => "does not segregate",
            (false, true) => "several copies",
            (false, false) => "",
        }
    }
}

fn evidence_or_lowest(o: &Observation) -> f64 {
    o.evidence.score.unwrap_or(f64::NEG_INFINITY)
}

/// Evidence first; dosage and count order the observations not assessed.
fn by_strength(a: &Observation, b: &Observation) -> Ordering {
    evidence_or_lowest(a)
        .total_cmp(&evidence_or_lowest(b))
        .then(
            a.evidence
                .dosage
                .unwrap_or(0.0)
                .total_cmp(&b.evidence.dosage.unwrap_or(0.0)),
        )
        .then(a.count.cmp(&b.count))
}

pub(super) fn summarise<'a>(plate: &Plate<'a>) -> Vec<Sequence<'a>> {
    let mut by_sequence: BTreeMap<&str, Vec<(&Observation, Class)>> = BTreeMap::new();
    for (o, class) in plate.rows() {
        by_sequence
            .entry(o.sequence.as_str())
            .or_default()
            .push((o, class));
    }
    by_sequence
        .into_values()
        .map(|rows| {
            let observations = || rows.iter().map(|(o, _)| *o);
            let best = observations()
                .max_by(|a, b| by_strength(a, b))
                .expect("a row");
            let carried: Vec<f64> = observations()
                .filter_map(|o| o.evidence.dosage)
                .filter(|d| *d >= COPY_LIKE)
                .collect();
            let supporting: HashSet<usize> = observations()
                .filter(|o| o.evidence.score.is_some_and(|s| s >= SUPPORT))
                .map(|o| o.donor)
                .collect();
            Sequence {
                class: rows.iter().map(|(_, c)| *c).min().expect("a row"),
                best,
                max_count: observations().map(|o| o.count).max().unwrap_or(0),
                max_ratio: observations().map(|o| o.ratio).fold(0.0, f64::max),
                donors: observations()
                    .map(|o| o.donor)
                    .collect::<HashSet<_>>()
                    .len(),
                supporting: supporting.len(),
                passing_evidence: observations()
                    .filter(|o| plate.thresholds.passes(o))
                    .map(|o| o.evidence.score.unwrap_or(f64::INFINITY))
                    .reduce(f64::max),
                typical_dosage: if carried.is_empty() {
                    best.evidence.dosage.unwrap_or(f64::NAN)
                } else {
                    median(&carried)
                },
            }
        })
        .collect()
}

pub(super) fn write_overview(out: &mut String, plate: &Plate, sequences: &[Sequence]) {
    let _ = writeln!(
        out,
        "\n== Sequences, each judged on the donor where its evidence is strongest"
    );
    let _ = writeln!(
        out,
        "  class         sequences  pass {}  evidence: none    <1   1-2   >=2  not assessed   does not segregate  >{MULTI_COPY} copies",
        plate.thresholds.label()
    );
    let pad = plate.thresholds.label().len() + 5;
    for class in Class::ALL {
        let items: Vec<&Sequence> = sequences.iter().filter(|s| s.class == class).collect();
        if items.is_empty() {
            continue;
        }
        let count = |f: &dyn Fn(&Sequence) -> bool| items.iter().filter(|s| f(s)).count();
        let within = |low: f64, high: f64| {
            count(&|s| s.evidence().is_some_and(|e| (low..high).contains(&e)))
        };
        let _ = writeln!(
            out,
            "  {:<12} {:10}  {:>pad$}           {:5} {:5} {:5} {:5}  {:12}   {:18}  {:10}",
            class.label(),
            items.len(),
            count(&|s| s.passes()),
            count(&|s| s.evidence() == Some(0.0)),
            within(f64::MIN_POSITIVE, 1.0),
            within(1.0, 2.0),
            within(2.0, f64::INFINITY),
            count(&|s| s.evidence().is_none()),
            count(&|s| s.ubiquitous()),
            count(&|s| s.multi_copy()),
        );
    }
}

pub(super) fn write_separation(out: &mut String, plate: &Plate, sequences: &[Sequence]) {
    let _ = writeln!(
        out,
        "\n== Targets against other non-database sequences, per sequence (AUC; 0.5 = none, 1 = perfect)"
    );
    let _ = writeln!(
        out,
        "  sequences considered                 targets  others  max count  max ratio  dosage  evidence"
    );
    type Keep<'k> = Box<dyn Fn(&Sequence) -> bool + 'k>;
    let subsets: [(String, Keep); 3] = [
        ("assessed".to_string(), Box::new(|_| true)),
        (
            "assessed, segregating".to_string(),
            Box::new(|s| !s.ubiquitous()),
        ),
        (
            format!("passing {} in a donor", plate.thresholds.label()),
            Box::new(|s| s.passes()),
        ),
    ];
    for (label, keep) in &subsets {
        let pick = |class: Class| -> Vec<&Sequence> {
            sequences
                .iter()
                .filter(|s| s.class == class && s.evidence().is_some() && keep(s))
                .collect()
        };
        let (targets, others) = (pick(Class::Target), pick(Class::Other));
        let separation = |feature: &dyn Fn(&Sequence) -> f64| {
            let t: Vec<f64> = targets.iter().map(|s| feature(s)).collect();
            let n: Vec<f64> = others.iter().map(|s| feature(s)).collect();
            number(auc(&t, &n))
        };
        let _ = writeln!(
            out,
            "  {label:<36} {:7} {:7}     {}     {}  {}    {}",
            targets.len(),
            others.len(),
            separation(&|s| s.max_count as f64),
            separation(&|s| s.max_ratio),
            separation(&|s| s.best.evidence.dosage.unwrap_or(f64::NAN)),
            separation(&|s| s.evidence().unwrap_or(f64::NAN)),
        );
    }
}

/// Sequences kept by the count thresholds and by cuts on the evidence.
pub(super) fn write_rules(out: &mut String, plate: &Plate, sequences: &[Sequence]) {
    let of = |class: Class| -> Vec<&Sequence> {
        sequences.iter().filter(|s| s.class == class).collect()
    };
    let (targets, other_cohort, others) =
        (of(Class::Target), of(Class::OtherTarget), of(Class::Other));
    let _ = writeln!(
        out,
        "\n== Count thresholds against evidence: sequences kept (targets on the plate: {}, other-cohort: {}, other: {})",
        targets.len(),
        other_cohort.len(),
        others.len()
    );
    let _ = writeln!(
        out,
        "  rule, met in at least one donor               targets  other-cohort  other, segregating  other, not segregating"
    );
    let mut line = |label: &str, keep: &dyn Fn(&Sequence) -> bool| {
        let kept = |items: &[&Sequence]| items.iter().filter(|s| keep(s)).count();
        let _ = writeln!(
            out,
            "  {label:<44} {:7}  {:12}  {:18}  {:18}",
            kept(&targets),
            kept(&other_cohort),
            others.iter().filter(|s| keep(s) && !s.ubiquitous()).count(),
            others.iter().filter(|s| keep(s) && s.ubiquitous()).count(),
        );
    };
    let thresholds = plate.thresholds.label();
    line(&thresholds, &|s| s.passes());
    for cut in [1.0, 2.0] {
        line(&format!("{thresholds} and evidence >= {cut}"), &|s| {
            s.passing_evidence.is_some_and(|e| e >= cut)
        });
    }
    line("fewer artifacts than candidates (evidence > 0)", &|s| {
        s.evidence().is_some_and(|e| e > 0.0)
    });
    for cut in CUTS {
        line(&format!("evidence >= {cut}"), &|s| s.at_least(cut));
    }
}

/// How many other sequences a cut on the evidence removes from today's
/// calls, by the number of targets it costs.
pub(super) fn write_cost(out: &mut String, plate: &Plate, sequences: &[Sequence]) {
    let passing = |class: Class| -> Vec<&Sequence> {
        sequences
            .iter()
            .filter(|s| s.class == class && s.passes())
            .collect()
    };
    let (targets, others) = (passing(Class::Target), passing(Class::Other));
    let mut of_targets: Vec<f64> = targets.iter().filter_map(|s| s.passing_evidence).collect();
    of_targets.sort_by(f64::total_cmp);
    let _ = writeln!(
        out,
        "\n== {} and evidence >= cut: what the cut removes ({} targets, {} other sequences pass today)",
        plate.thresholds.label(),
        targets.len(),
        others.len()
    );
    let _ = writeln!(
        out,
        "  targets lost    cut   other, segregating  other, not segregating   of them 4+ differences from the database"
    );
    // Every cut that would cost one more target, weakest target first.
    let mut cuts = of_targets.clone();
    cuts.dedup();
    for cut in cuts.into_iter().take(MAX_COST_LINES) {
        let kept: Vec<&&Sequence> = others
            .iter()
            .filter(|s| s.passing_evidence.is_some_and(|e| e >= cut))
            .collect();
        let _ = writeln!(
            out,
            "  {:12}  {}   {:18}  {:22}   {:7}",
            of_targets.iter().filter(|e| **e < cut).count(),
            number(cut),
            kept.iter().filter(|s| !s.ubiquitous()).count(),
            kept.iter().filter(|s| s.ubiquitous()).count(),
            kept.iter().filter(|s| s.best.db_distance >= 4).count(),
        );
    }
}

/// Differences from the nearest database allele, for today's calls.
pub(super) fn write_distances(out: &mut String, plate: &Plate, sequences: &[Sequence]) {
    let _ = writeln!(
        out,
        "\n== Differences from the nearest database allele, sequences passing {}",
        plate.thresholds.label()
    );
    let _ = writeln!(
        out,
        "  class         sequences      1      2      3    4-5   6-10    >10"
    );
    for class in [Class::Target, Class::OtherTarget, Class::Other] {
        let items: Vec<&Sequence> = sequences
            .iter()
            .filter(|s| s.class == class && s.passes())
            .collect();
        if items.is_empty() {
            continue;
        }
        let within = |low: u64, high: u64| {
            items
                .iter()
                .filter(|s| (low..=high).contains(&s.best.db_distance))
                .count()
        };
        let _ = writeln!(
            out,
            "  {:<12} {:10} {:6} {:6} {:6} {:6} {:6} {:6}",
            class.label(),
            items.len(),
            within(0, 1),
            within(2, 2),
            within(3, 3),
            within(4, 5),
            within(6, 10),
            within(11, u64::MAX),
        );
    }
}

fn write_sequence_lines(out: &mut String, items: &[&Sequence]) {
    let _ = writeln!(
        out,
        "  gene            nearest allele        d(db)  kind            count  dosage  evidence  reads per copy  donors  with evidence >= {SUPPORT}  carriers  typical copies  usual % of parent  passes  note"
    );
    for s in items {
        let o = s.best;
        let _ = writeln!(
            out,
            "  {:<15} {:<22} {:4}  {:<14} {:6}  {}    {}          {}  {:6}  {:20}  {:8}          {}             {}  {:<6}  {}",
            o.gene,
            o.nearest,
            o.db_distance,
            Kind::of(o).label(),
            o.count,
            optional(o.evidence.dosage),
            optional(s.evidence()),
            optional(o.evidence.per_copy_reads),
            s.donors,
            s.supporting,
            o.evidence.carriers,
            number(s.typical_dosage),
            optional(o.evidence.recurrent_ratio.map(|r| 100.0 * r)),
            if s.passes() { "yes" } else { "no" },
            s.note()
        );
    }
}

pub(super) fn write_strongest_others(out: &mut String, sequences: &[Sequence], list: usize) {
    if list == 0 {
        return;
    }
    let _ = writeln!(
        out,
        "\n== Other sequences with the strongest evidence (unknown alleles, sequences of other loci, or what the evidence gets wrong)"
    );
    let mut others: Vec<&Sequence> = sequences
        .iter()
        .filter(|s| s.class == Class::Other && s.evidence().is_some_and(|e| e > 0.0))
        .collect();
    if others.is_empty() {
        let _ = writeln!(out, "  none with evidence above 0");
        return;
    }
    others.sort_by(|a, b| by_strength(b.best, a.best));
    let total = others.len();
    others.truncate(list);
    write_sequence_lines(out, &others);
    if total > others.len() {
        let _ = writeln!(out, "  ... and {} more", total - others.len());
    }
}

pub(super) fn write_weakest_targets(out: &mut String, sequences: &[Sequence], list: usize) {
    if list == 0 {
        return;
    }
    let _ = writeln!(
        out,
        "\n== Targets with evidence below {SUPPORT} in every donor"
    );
    let mut weak: Vec<&Sequence> = sequences
        .iter()
        .filter(|s| {
            // The donors the sequence is called in today count, where there
            // are any; one too shallow to assess does not speak against it.
            s.class == Class::Target
                && s.passing_evidence
                    .or(s.evidence())
                    .is_some_and(|e| e < SUPPORT)
        })
        .collect();
    weak.sort_by(|a, b| by_strength(a.best, b.best));
    let total = weak.len();
    weak.truncate(list);
    if weak.is_empty() {
        let _ = writeln!(out, "  none");
        return;
    }
    write_sequence_lines(out, &weak);
    if total > weak.len() {
        let _ = writeln!(out, "  ... and {} more", total - weak.len());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::check::evidence_report::Thresholds;
    use crate::evidence::assess;
    use crate::evidence::testing::observation;

    #[test]
    fn a_sequence_is_judged_on_its_strongest_donor() {
        let scored = |donor: usize, sequence: &str, count: u64, score: Option<f64>| {
            let mut o = observation(donor, "A", sequence, count);
            o.ratio = 1.0;
            o.best_full_count = count;
            o.best_full_ratio = 1.0;
            o.evidence.score = score;
            o.evidence.dosage = Some(count as f64 / 100.0);
            o
        };
        let observations = vec![
            scored(0, "AAAA", 100, Some(2.5)),
            scored(1, "AAAA", 3, Some(-2.0)),
            scored(2, "AAAA", 200, None),
            scored(0, "CCCC", 4, Some(-1.0)),
        ];
        let classes = [Class::Other, Class::Target, Class::Other, Class::Other];
        let models = assess(&mut Vec::new(), 3);
        let plate = Plate {
            observations: &observations,
            classes: &classes,
            n_donors: 3,
            models: &models,
            thresholds: Thresholds {
                min_count: 50,
                min_ratio: 0.1,
            },
        };
        let sequences = summarise(&plate);
        assert_eq!(sequences.len(), 2);
        let a = &sequences[0];
        assert_eq!(a.class, Class::Target);
        assert_eq!(a.best.donor, 0);
        assert_eq!((a.donors, a.supporting, a.max_count), (3, 1, 200));
        // Donors 0 and 2 pass the thresholds; donor 2 is too shallow to be
        // assessed, which does not speak against the sequence.
        assert_eq!(a.passing_evidence, Some(f64::INFINITY));
        // Carried at 1 and 2 copies (the 3-read donor does not carry it).
        assert_eq!(a.typical_dosage, 1.5);
        assert!(!sequences[1].passes());
    }
}
