//! What was read, what the plate-level models learned, and the resulting
//! evidence scale.

use std::collections::{BTreeMap, HashSet};
use std::fmt::Write as _;

use super::{number, Class, Plate};
use crate::evidence::artifact_rate::MIN_COPY_READS;
use crate::evidence::artifacts::{is_sporadic_satellite, variants_at};
use crate::evidence::plate::COPY_LIKE;
use crate::evidence::{Kind, Observation};
use crate::stats::robust::{mad, quantile};
use crate::stats::tails::log10_negative_binomial_upper_tail;

/// Share of a gene's database-allele rows below the edge that makes its
/// unit suspect.
const SUSPECT_SHARE: f64 = 0.2;
const MAX_SUSPECT_GENES: usize = 15;

fn distinct(items: &[&Observation]) -> usize {
    items
        .iter()
        .map(|o| o.sequence.as_str())
        .collect::<HashSet<_>>()
        .len()
}

fn shape(shape: f64) -> String {
    if shape.is_infinite() {
        "Poisson".to_string()
    } else {
        format!("{shape:.3}")
    }
}

pub(super) fn write_input(out: &mut String, plate: &Plate) {
    let min_count = plate.thresholds.min_count;
    let _ = writeln!(out, "== Input");
    let _ = writeln!(
        out,
        "donors (well + case): {}   genes: {}   sequence-in-donor rows: {}",
        plate.n_donors,
        plate.models.amplification.genes.len(),
        plate.observations.len()
    );
    for class in Class::ALL {
        let items: Vec<&Observation> = plate.of_class(class).collect();
        let _ = writeln!(
            out,
            "  {:<12} rows {:7}   distinct sequences {:6}   rows with >= {min_count} reads {:6}",
            class.label(),
            items.len(),
            distinct(&items),
            items.iter().filter(|o| o.count >= min_count).count()
        );
    }
    let shallow = plate
        .observations
        .iter()
        .filter(|o| o.evidence.score.is_none())
        .count();
    let _ = writeln!(
        out,
        "rows not assessed (gene too shallow in the donor: fewer than {MIN_COPY_READS} reads expected per copy): {shallow}"
    );
}

pub(super) fn write_amplification(out: &mut String, plate: &Plate) {
    let log2 = std::f64::consts::LN_2;
    let min_count = plate.thresholds.min_count;
    let amplification = &plate.models.amplification;
    let _ = writeln!(
        out,
        "\n== Amplification model: log(reads of gene in donor) = donor + gene"
    );
    let residuals = &amplification.residuals;
    let _ = writeln!(
        out,
        "cells fitted {}   residual spread (robust SD, log2) {:.2}   5%..95% residual (log2) {:.2} .. {:.2}",
        residuals.len(),
        mad(residuals) / log2,
        quantile(residuals, 0.05) / log2,
        quantile(residuals, 0.95) / log2
    );
    let finite =
        |values: &[f64]| -> Vec<f64> { values.iter().copied().filter(|v| !v.is_nan()).collect() };
    let (genes, donors) = (finite(&amplification.gene), finite(&amplification.donor));
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
        plate.models.allele_noise
    );
    let _ = writeln!(
        out,
        "estimated copies (dosage) by class, rows with >= {min_count} reads:"
    );
    let _ = writeln!(
        out,
        "  class            n      5%    25%    50%    75%    95%   share in 0.5-2.8"
    );
    for class in Class::ALL {
        let dosage: Vec<f64> = plate
            .of_class(class)
            .filter(|o| o.count >= min_count)
            .filter_map(|o| o.evidence.dosage)
            .collect();
        if dosage.is_empty() {
            continue;
        }
        let in_range = dosage.iter().filter(|d| (0.5..=2.8).contains(*d)).count();
        let _ = writeln!(
            out,
            "  {:<12} {:5}  {} {} {} {} {}   {:.2}",
            class.label(),
            dosage.len(),
            number(quantile(&dosage, 0.05)),
            number(quantile(&dosage, 0.25)),
            number(quantile(&dosage, 0.5)),
            number(quantile(&dosage, 0.75)),
            number(quantile(&dosage, 0.95)),
            in_range as f64 / dosage.len() as f64
        );
    }
}

/// Genes whose database alleles are often below a copy-like amount: the
/// gene has more than two copies, or another sequence inflates its total.
/// A real allele of such a gene can fall below the edge.
pub(super) fn write_gene_units(out: &mut String, plate: &Plate) {
    let edge = plate.models.artifact_rate.edge;
    // gene => dosages of the database alleles present beyond stray reads
    let mut by_gene: BTreeMap<&str, Vec<f64>> = BTreeMap::new();
    for o in plate.of_class(Class::Known) {
        if let Some(dosage) = o.evidence.dosage.filter(|d| *d >= COPY_LIKE) {
            by_gene.entry(o.gene.as_str()).or_default().push(dosage);
        }
    }
    let below = |dosages: &[f64]| {
        dosages.iter().filter(|d| **d < edge).count() as f64 / dosages.len() as f64
    };
    let mut suspect: Vec<(&str, &Vec<f64>)> = by_gene
        .iter()
        .filter(|(_, dosages)| dosages.len() >= 5 && below(dosages) >= SUSPECT_SHARE)
        .map(|(gene, dosages)| (*gene, dosages))
        .collect();
    suspect.sort_by(|a, b| below(b.1).total_cmp(&below(a.1)));
    let _ = writeln!(
        out,
        "\n== Genes where one copy may not be half the gene: at least {:.0}% of database-allele rows between {COPY_LIKE} and {edge:.2} copies ({} of {} genes)",
        100.0 * SUSPECT_SHARE,
        suspect.len(),
        by_gene.len()
    );
    if suspect.is_empty() {
        return;
    }
    let _ = writeln!(
        out,
        "  gene             rows  share below {edge:.2}    25%    50%    75% of dosage"
    );
    for (gene, dosages) in suspect.iter().take(MAX_SUSPECT_GENES) {
        let _ = writeln!(
            out,
            "  {gene:<15} {:5}              {:.2} {} {} {}",
            dosages.len(),
            below(dosages),
            number(quantile(dosages, 0.25)),
            number(quantile(dosages, 0.5)),
            number(quantile(dosages, 0.75)),
        );
    }
    if suspect.len() > MAX_SUSPECT_GENES {
        let _ = writeln!(out, "  ... and {} more", suspect.len() - MAX_SUSPECT_GENES);
    }
}

/// The error models are not part of the evidence; they describe the plate.
pub(super) fn write_errors(out: &mut String, plate: &Plate) {
    let observations = plate.observations;
    let _ = writeln!(
        out,
        "\n== Errors on this plate (description only, not part of the evidence)"
    );
    for fit in plate.models.sporadic.iter().filter(|fit| fit.distance == 1) {
        let _ = writeln!(
            out,
            "one specific single-base variant: rate per read {:.3e}   shape {}   error sources {}   satellites {} ({} reads)",
            fit.rate,
            shape(fit.shape),
            fit.sources,
            fit.satellites,
            fit.satellite_reads
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
                    let tail =
                        log10_negative_binomial_upper_tail(k, o.count as f64 * fit.rate, fit.shape);
                    variants_at(o.sequence.len(), fit.distance) * 10f64.powf(tail)
                })
                .sum();
            let _ = write!(out, "  k={k}: {observed} / {expected:.1}");
        }
        let _ = writeln!(out);
    }
    let recurring: Vec<&Observation> = observations
        .iter()
        .filter(|o| o.evidence.parent.is_some() && !o.is_database_allele())
        .filter(|o| {
            o.evidence
                .recurrent_ratio
                .is_some_and(|r| r > 0.0 && r <= 0.1)
        })
        .collect();
    let ratios: Vec<f64> = recurring
        .iter()
        .filter_map(|o| o.evidence.recurrent_ratio)
        .collect();
    let _ = writeln!(
        out,
        "recurrent variants (same variant of the same parent in at least half of the other donors, at most 10% of it): {} rows, {} distinct, usual ratio to parent 25%..75%: {} .. {}, between-donor shape {}",
        recurring.len(),
        distinct(&recurring),
        number(quantile(&ratios, 0.25)),
        number(quantile(&ratios, 0.75)),
        shape(plate.models.recurrent.shape)
    );
}

pub(super) fn write_evidence_scale(out: &mut String, plate: &Plate) {
    let model = &plate.models.artifact_rate;
    let _ = writeln!(
        out,
        "\n== Evidence: -log10 of the share of artifacts expected among sequences of the same kind with at least the same dosage (1 = one in ten, 2 = one in a hundred)"
    );
    let _ = writeln!(
        out,
        "not a copy below a dosage of {:.2} (lowest 1% of single copies); artifacts counted in three windows of half a decade from {:.3} up to it, among sequences private to one donor and absent from most others",
        model.edge,
        model.edge / 10f64.powf(1.5)
    );
    let _ = writeln!(
        out,
        "  kind (vs. more abundant sequence)  artifacts: lowest  middle  highest  thinning per half decade   expected above: edge  1 copy  2 copies   candidates: private    all   evidence at 1 copy: private  shared"
    );
    for kind in Kind::ALL {
        let tail = model.tail(kind);
        let (private, all) = model.candidates(kind);
        let _ = writeln!(
            out,
            "  {:<33} {:17} {:7} {:8}  {:>8} {:<15}                {}  {}    {}   {:19} {:6}                       {:5.2}   {:5.2}",
            kind.label(),
            tail.windows[0],
            tail.windows[1],
            tail.windows[2],
            format!("{:.3}", tail.decline),
            if tail.own_decline { "" } else { "(all kinds)" },
            number(tail.expected(model.edge)),
            number(tail.expected(1.0)),
            number(tail.expected(2.0)),
            private,
            all,
            model.evidence(kind, true, 1.0),
            model.evidence(kind, false, 1.0),
        );
    }
}
