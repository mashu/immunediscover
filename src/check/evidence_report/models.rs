//! What was read, what the plate-level models learned, and the resulting
//! evidence scale.

use std::collections::HashSet;
use std::fmt::Write as _;

use super::{number, Class, Plate};
use crate::evidence::artifact_rate::{ARTIFACT_ONLY, MIN_COPY_READS};
use crate::evidence::artifacts::{is_sporadic_satellite, variants_at};
use crate::evidence::{Kind, Observation};
use crate::stats::robust::{mad, quantile};
use crate::stats::tails::log10_negative_binomial_upper_tail;

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
    let own_unit = plate
        .observations
        .iter()
        .filter(|o| o.evidence.score.is_some() && !plate.models.amplification.is_fitted(o))
        .count();
    let _ = writeln!(
        out,
        "rows assessed with the unit of their own donor (gene in too few donors): {own_unit}"
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

/// Genes whose unit is not the default of two copies fitted on database
/// alleles.
pub(super) fn write_gene_units(out: &mut String, plate: &Plate) {
    let amplification = &plate.models.amplification;
    let named = |keep: &dyn Fn(usize) -> bool| -> String {
        let names: Vec<String> = (0..amplification.genes.len())
            .filter(|g| keep(*g))
            .map(|g| {
                if amplification.copies[g] > 2.0 {
                    format!("{} ({})", amplification.genes[g], amplification.copies[g])
                } else {
                    amplification.genes[g].clone()
                }
            })
            .collect();
        if names.is_empty() {
            "none".to_string()
        } else {
            names.join(", ")
        }
    };
    let fitted = |g: usize| !amplification.gene[g].is_nan();
    let _ = writeln!(out, "\n== The unit of one copy, by gene");
    let _ = writeln!(
        out,
        "more than two copies in most donors (copies): {}",
        named(&|g| amplification.copies[g] > 2.0)
    );
    let _ = writeln!(
        out,
        "no database allele in enough donors, fitted on all main sequences: {}",
        named(&|g| fitted(g) && !amplification.from_database[g])
    );
    let _ = writeln!(
        out,
        "in too few donors to fit, unit taken from the donor itself: {}",
        named(&|g| !fitted(g))
    );
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
    let _ = writeln!(
        out,
        "recurrent artifacts (usual ratio at most 10% of the parent): {} rows, ratio varies between donors by a factor of {:.2} (1 SD); a row is judged against its own usual ratio",
        plate.models.recurrence.rows,
        plate.models.recurrence.spread.exp()
    );
}

const SCALE_DOSAGES: [f64; 4] = [0.3, 0.5, 1.0, 2.0];

pub(super) fn write_evidence_scale(out: &mut String, plate: &Plate) {
    let model = &plate.models.artifact_rate;
    let _ = writeln!(
        out,
        "\n== Evidence: -log10 of the share of artifacts expected among sequences of the same kind around the same dosage (1 = one in ten, 2 = one in a hundred)"
    );
    let _ = writeln!(
        out,
        "artifacts counted in two windows of half a decade below {ARTIFACT_ONLY} copies, among sequences private to one donor and absent from most others; their decline is continued upwards"
    );
    let _ = writeln!(
        out,
        "  kind (vs. more abundant sequence)  artifacts: lower  upper  thinning per half decade   expected above: {ARTIFACT_ONLY}     0.6       2   candidates: private    all   evidence, private at: 0.3   0.5     1     2   shared at: 0.3   0.5     1     2"
    );
    for kind in Kind::ALL {
        let tail = model.tail(kind);
        let (private, all) = model.candidates(kind);
        let _ = write!(
            out,
            "  {:<33} {:16} {:6}  {:>8} {:<15}                {}  {}  {}   {:19} {:6}                      ",
            kind.label(),
            tail.windows[0],
            tail.windows[1],
            format!("{:.3}", tail.decline),
            if tail.own_decline { "" } else { "(all kinds)" },
            number(tail.expected(ARTIFACT_ONLY)),
            number(tail.expected(0.6)),
            number(tail.expected(2.0)),
            private,
            all,
        );
        for is_private in [true, false] {
            for dosage in SCALE_DOSAGES {
                let _ = write!(out, " {:5.2}", model.evidence(kind, is_private, dosage));
            }
            let _ = write!(out, "            ");
        }
        let _ = writeln!(out);
    }
}
