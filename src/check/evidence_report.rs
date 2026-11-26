//! Does the evidence of [`crate::evidence`] separate held-out alleles from
//! everything else on a plate with known truth?
//!
//! The report is aggregate only: no donor identifiers and no sequences.

use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use std::path::Path;

use crate::error::Result;
use crate::evidence::artifacts::{is_sporadic_satellite, variants_at};
use crate::evidence::{Models, Observation};
use crate::io::fasta::read_fasta;
use crate::seq::gene_of;
use crate::stats::median;
use crate::stats::robust::{auc, mad, quantile};
use crate::stats::tails::log10_negative_binomial_upper_tail;

/// What a sequence is, according to the truth set.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Class {
    /// A database allele.
    Known,
    /// A part of a database allele of the same gene.
    Partial,
    /// A held-out allele that should be rediscovered.
    Target,
    /// Anything else: artifacts, and real alleles nobody knows yet.
    Other,
}

impl Class {
    const ALL: [Class; 4] = [Class::Known, Class::Partial, Class::Target, Class::Other];

    fn label(self) -> &'static str {
        match self {
            Class::Known => "known",
            Class::Partial => "partial",
            Class::Target => "target",
            Class::Other => "other",
        }
    }
}

/// Classify every observation by exact sequence match against the database
/// FASTA and the FASTA of held-out target alleles.
pub fn classify(observations: &[Observation], db: &Path, targets: &Path) -> Result<Vec<Class>> {
    let mut known: HashSet<String> = HashSet::new();
    let mut by_gene: HashMap<String, Vec<String>> = HashMap::new();
    for record in read_fasta(db)? {
        let sequence = record.sequence.to_ascii_uppercase();
        by_gene
            .entry(gene_of(record.identifier()).to_string())
            .or_default()
            .push(sequence.clone());
        known.insert(sequence);
    }
    let targets: HashSet<String> = read_fasta(targets)?
        .into_iter()
        .map(|r| r.sequence.to_ascii_uppercase())
        .collect();
    let mut cache: HashMap<(&str, &str), Class> = HashMap::new();
    Ok(observations
        .iter()
        .map(|o| {
            *cache
                .entry((o.gene.as_str(), o.sequence.as_str()))
                .or_insert_with(|| {
                    if known.contains(&o.sequence) {
                        Class::Known
                    } else if targets.contains(&o.sequence) {
                        Class::Target
                    } else if by_gene
                        .get(&o.gene)
                        .is_some_and(|alleles| alleles.iter().any(|a| a.contains(&o.sequence)))
                    {
                        Class::Partial
                    } else {
                        Class::Other
                    }
                })
        })
        .collect())
}

/// The count and ratio thresholds of the run being compared against.
#[derive(Debug, Clone, Copy)]
pub struct Thresholds {
    pub min_count: u64,
    pub min_ratio: f64,
}

impl Thresholds {
    /// Both levels, as `discover blast` applies them.
    fn passes(&self, o: &Observation) -> bool {
        self.passes_allele_level(o)
            && o.best_full_count >= self.min_count
            && o.best_full_ratio >= self.min_ratio
    }

    fn passes_allele_level(&self, o: &Observation) -> bool {
        o.count >= self.min_count && o.ratio >= self.min_ratio
    }
}

fn number(x: f64) -> String {
    if x.is_nan() {
        "   n/a".to_string()
    } else if x != 0.0 && (x.abs() < 1e-3 || x.abs() >= 1e5) {
        format!("{x:6.1e}")
    } else {
        let digits = (2 - x.abs().max(1e-9).log10().floor() as i32).clamp(0, 6) as usize;
        format!("{x:6.digits$}")
    }
}

fn unexplained(o: &Observation) -> f64 {
    o.evidence.log10_artifacts.map_or(30.0, |e| -e)
}

fn recurs(o: &Observation) -> bool {
    o.evidence
        .recurrent_ratio
        .is_some_and(|r| r > 0.0 && r <= 0.1)
}

/// The evidence report for one plate.
pub fn format_evidence_report(
    observations: &[Observation],
    classes: &[Class],
    n_donors: usize,
    models: &Models,
    thresholds: Thresholds,
    list: usize,
) -> String {
    let mut out = String::new();
    let log2 = std::f64::consts::LN_2;
    let (min_count, min_ratio) = (thresholds.min_count, thresholds.min_ratio);
    let rows: Vec<(&Observation, Class)> =
        observations.iter().zip(classes.iter().copied()).collect();
    let of_class = |class: Class| {
        rows.iter()
            .filter(move |(_, c)| *c == class)
            .map(|(o, _)| *o)
    };
    let distinct = |items: &[&Observation]| {
        items
            .iter()
            .map(|o| o.sequence.as_str())
            .collect::<HashSet<_>>()
            .len()
    };

    let _ = writeln!(out, "== Input");
    let _ = writeln!(
        out,
        "donors (well + case): {n_donors}   genes: {}   sequence-in-donor rows: {}",
        models.amplification.genes.len(),
        observations.len()
    );
    for class in Class::ALL {
        let items: Vec<&Observation> = of_class(class).collect();
        let _ = writeln!(
            out,
            "  {:<8} rows {:7}   distinct sequences {:6}   rows with >= {min_count} reads {:6}",
            class.label(),
            items.len(),
            distinct(&items),
            items.iter().filter(|o| o.count >= min_count).count()
        );
    }

    let _ = writeln!(
        out,
        "\n== Amplification model: log(reads of gene in donor) = donor + gene"
    );
    let residuals = &models.amplification.residuals;
    let _ = writeln!(
        out,
        "cells fitted {}   residual spread (robust SD, log2) {:.2}   5%..95% residual (log2) {:.2} .. {:.2}",
        residuals.len(),
        mad(residuals) / log2,
        quantile(residuals, 0.05) / log2,
        quantile(residuals, 0.95) / log2
    );
    let genes: Vec<f64> = models
        .amplification
        .gene
        .iter()
        .copied()
        .filter(|g| !g.is_nan())
        .collect();
    let donors: Vec<f64> = models
        .amplification
        .donor
        .iter()
        .copied()
        .filter(|d| !d.is_nan())
        .collect();
    let _ = writeln!(
        out,
        "primer efficiency between genes (log2, 5%..95%): {:.2} .. {:.2}   donor depth (log2, 5%..95%): {:.2} .. {:.2}",
        quantile(&genes, 0.05) / log2,
        quantile(&genes, 0.95) / log2,
        quantile(&donors, 0.05) / log2,
        quantile(&donors, 0.95) / log2
    );
    let _ = writeln!(
        out,
        "allele-level spread around 1 or 2 copies (SD of log, database alleles): {:.2}",
        models.allele_noise
    );
    let _ = writeln!(
        out,
        "estimated copies (dosage) by class, rows with >= {min_count} reads:"
    );
    let _ = writeln!(
        out,
        "  class      n      5%    25%    50%    75%    95%   share in 0.5-2.8"
    );
    for class in Class::ALL {
        let dosage: Vec<f64> = of_class(class)
            .filter(|o| o.count >= min_count)
            .filter_map(|o| o.evidence.dosage)
            .collect();
        if dosage.is_empty() {
            continue;
        }
        let _ = writeln!(
            out,
            "  {:<8} {:5}  {} {} {} {} {}   {:.2}",
            class.label(),
            dosage.len(),
            number(quantile(&dosage, 0.05)),
            number(quantile(&dosage, 0.25)),
            number(quantile(&dosage, 0.5)),
            number(quantile(&dosage, 0.75)),
            number(quantile(&dosage, 0.95)),
            dosage.iter().filter(|d| (0.5..=2.8).contains(*d)).count() as f64 / dosage.len() as f64
        );
    }

    let _ = writeln!(
        out,
        "\n== Error model: reads of a sequence becoming one specific variant"
    );
    for fit in &models.sporadic {
        let shape = if fit.shape.is_infinite() {
            "Poisson".to_string()
        } else {
            format!("{:.3}", fit.shape)
        };
        let _ = writeln!(
            out,
            "distance {}: rate per read {:.3e}   shape {shape}   error sources {}   satellites {} ({} reads)",
            fit.distance, fit.rate, fit.sources, fit.satellites, fit.satellite_reads
        );
        if fit.rate <= 0.0 {
            continue;
        }
        let _ = write!(out, "   satellites with >= k reads, observed / expected:");
        for k in [1u64, 2, 3, 5, 10, 20] {
            let observed = observations
                .iter()
                .filter(|o| o.count >= k && is_sporadic_satellite(observations, o, fit.distance))
                .count();
            let expected: f64 = observations
                .iter()
                .filter(|o| o.count >= 50)
                .map(|o| {
                    variants_at(o.sequence.len(), fit.distance)
                        * 10f64.powf(log10_negative_binomial_upper_tail(
                            k,
                            o.count as f64 * fit.rate,
                            fit.shape,
                        ))
                })
                .sum();
            let _ = write!(out, "  k={k}: {observed} / {expected:.1}");
        }
        let _ = writeln!(out);
    }
    let recurring: Vec<&Observation> = observations
        .iter()
        .filter(|o| o.evidence.parent.is_some() && recurs(o) && !o.is_database_allele())
        .collect();
    let ratios: Vec<f64> = recurring
        .iter()
        .filter_map(|o| o.evidence.recurrent_ratio)
        .collect();
    let _ = writeln!(
        out,
        "recurrent artifacts (same variant of the same parent in at least half of the other donors): {} rows, {} distinct variants, usual ratio to parent 25%..75%: {} .. {}, between-donor shape {}",
        recurring.len(),
        distinct(&recurring),
        number(quantile(&ratios, 0.25)),
        number(quantile(&ratios, 0.75)),
        if models.recurrent.shape.is_infinite() { "Poisson".to_string() } else { format!("{:.3}", models.recurrent.shape) }
    );

    let scored: Vec<(&Observation, Class)> = rows
        .iter()
        .copied()
        .filter(|(o, _)| o.evidence.score.is_some())
        .collect();
    let score = |o: &Observation| o.evidence.score.unwrap_or(f64::NAN);
    let _ = writeln!(
        out,
        "\n== Separation of held-out targets from other non-database sequences (AUC; 0.5 = none, 1 = perfect)"
    );
    let _ = writeln!(
        out,
        "  rows considered                targets  others   count   ratio  dosage  -log10E   score"
    );
    type Keep = Box<dyn Fn(&Observation) -> bool>;
    let subsets: [(String, Keep); 3] = [
        (
            "all with >= 2 reads".to_string(),
            Box::new(|o| o.count >= 2),
        ),
        (
            format!("all with >= {min_count} reads"),
            Box::new(move |o| o.count >= min_count),
        ),
        (
            format!("passing -c {min_count} -f {min_ratio}"),
            Box::new(move |o| thresholds.passes(o)),
        ),
    ];
    for (label, keep) in &subsets {
        let pick = |class: Class| -> Vec<&Observation> {
            scored
                .iter()
                .filter(|(o, c)| *c == class && keep(o))
                .map(|(o, _)| *o)
                .collect()
        };
        let (targets, others) = (pick(Class::Target), pick(Class::Other));
        let separation = |feature: &dyn Fn(&Observation) -> f64| {
            let t: Vec<f64> = targets.iter().map(|o| feature(o)).collect();
            let n: Vec<f64> = others.iter().map(|o| feature(o)).collect();
            number(auc(&t, &n))
        };
        let _ = writeln!(
            out,
            "  {label:<30} {:7} {:7}  {}  {}  {}  {}  {}",
            targets.len(),
            others.len(),
            separation(&|o| o.count as f64),
            separation(&|o| o.ratio),
            separation(&|o| o.evidence.dosage.unwrap_or(f64::NAN)),
            separation(&|o| o.evidence.log10_artifacts.map_or(1e3, |e| -e)),
            separation(&score),
        );
    }

    let _ = writeln!(out, "\n== Thresholds versus score, on rows with >= 2 reads");
    let pool: Vec<(&Observation, Class)> = scored
        .iter()
        .copied()
        .filter(|(o, c)| o.count >= 2 && matches!(c, Class::Target | Class::Other))
        .collect();
    let mut describe = |label: &str, kept: &[(&Observation, Class)]| {
        let of = |class: Class| -> Vec<&Observation> {
            kept.iter()
                .filter(|(_, c)| *c == class)
                .map(|(o, _)| *o)
                .collect()
        };
        let (targets, others) = (of(Class::Target), of(Class::Other));
        let _ = writeln!(
            out,
            "  {label:<44} target rows {:5} (sequences {:4})   other rows {:5} (sequences {:4})",
            targets.len(),
            distinct(&targets),
            others.len(),
            distinct(&others)
        );
    };
    let current: Vec<_> = pool
        .iter()
        .copied()
        .filter(|(o, _)| thresholds.passes(o))
        .collect();
    describe(
        &format!("current: -c {min_count} -f {min_ratio} (both levels)"),
        &current,
    );
    let allele_level: Vec<_> = pool
        .iter()
        .copied()
        .filter(|(o, _)| thresholds.passes_allele_level(o))
        .collect();
    describe("count and ratio only (allele level)", &allele_level);
    let mut by_score = pool.clone();
    by_score.sort_by(|a, b| score(b.0).total_cmp(&score(a.0)));
    let targets_now = current.iter().filter(|(_, c)| *c == Class::Target).count();
    let others_now = current.len() - targets_now;
    let (mut targets_so_far, mut others_so_far) = (0, 0);
    let (mut same_targets, mut same_others) = (None, None);
    for (i, (_, class)) in by_score.iter().enumerate() {
        match class {
            Class::Target => targets_so_far += 1,
            _ => others_so_far += 1,
        }
        if same_targets.is_none() && targets_now > 0 && targets_so_far >= targets_now {
            same_targets = Some(i + 1);
        }
        if others_so_far <= others_now {
            same_others = Some(i + 1);
        }
    }
    if let Some(n) = same_targets {
        describe("score, same number of target rows", &by_score[..n]);
    }
    if let Some(n) = same_others {
        describe("score, same number of other rows", &by_score[..n]);
    }
    for cut in [25.0, 20.0, 15.0, 10.0, 5.0] {
        let kept: Vec<_> = pool
            .iter()
            .copied()
            .filter(|(o, _)| score(o) >= cut)
            .collect();
        describe(&format!("score >= {cut}"), &kept);
    }

    let _ = writeln!(
        out,
        "\n== Evidence by class (medians), rows with >= {min_count} reads"
    );
    let _ = writeln!(
        out,
        "  class      n   count  dosage  -log10E  log10 p(copy)  carriers  tracking  chimera  top of gene  recurs"
    );
    for class in Class::ALL {
        let items: Vec<&Observation> = scored
            .iter()
            .filter(|(o, c)| *c == class && o.count >= min_count)
            .map(|(o, _)| *o)
            .collect();
        if items.is_empty() {
            continue;
        }
        let med = |f: &dyn Fn(&Observation) -> Option<f64>| -> f64 {
            median(&items.iter().filter_map(|o| f(o)).collect::<Vec<f64>>())
        };
        let share = |f: &dyn Fn(&Observation) -> bool| {
            items.iter().filter(|o| f(o)).count() as f64 / items.len() as f64
        };
        let _ = writeln!(
            out,
            "  {:<8} {:5}  {:6.0}  {}   {}        {}    {:6.0}    {}    {:5.2}       {:5.2}    {:5.2}",
            class.label(),
            items.len(),
            med(&|o| Some(o.count as f64)),
            number(med(&|o| o.evidence.dosage)),
            number(med(&|o| Some(unexplained(o)))),
            number(med(&|o| o.evidence.log10_p_copy)),
            med(&|o| Some(o.evidence.carriers as f64)),
            number(med(&|o| o.evidence.tracking)),
            share(&|o| o.evidence.chimera),
            share(&|o| o.evidence.parent.is_none()),
            share(&|o| recurs(o)),
        );
    }

    if list > 0 {
        let _ = writeln!(
            out,
            "\n== Highest-scoring sequences that are neither database nor target (unknown real alleles, or what the score gets wrong)"
        );
        let _ = writeln!(
            out,
            "  gene            nearest allele        d(db) d(parent)  count  dosage  -log10E  carriers  tracking  recurs  chimera"
        );
        let mut others: Vec<&Observation> = scored
            .iter()
            .filter(|(o, c)| *c == Class::Other && o.count >= 2)
            .map(|(o, _)| *o)
            .collect();
        others.sort_by(|a, b| score(b).total_cmp(&score(a)));
        let mut seen: HashSet<&str> = HashSet::new();
        for o in others {
            if seen.len() >= list {
                break;
            }
            if !seen.insert(o.sequence.as_str()) {
                continue;
            }
            let _ = writeln!(
                out,
                "  {:<15} {:<22} {:4} {:>8} {:6}  {}   {}    {:5}    {}  {}  {}",
                o.gene,
                o.nearest,
                o.db_distance,
                o.evidence
                    .parent
                    .map_or("none".to_string(), |(_, d)| d.to_string()),
                o.count,
                number(o.evidence.dosage.unwrap_or(f64::NAN)),
                number(unexplained(o)),
                o.evidence.carriers,
                number(o.evidence.tracking.unwrap_or(f64::NAN)),
                number(o.evidence.recurrent_ratio.unwrap_or(f64::NAN)),
                if o.evidence.chimera { "yes" } else { "no" }
            );
        }
    }
    out
}
