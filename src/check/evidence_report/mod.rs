//! Does the evidence of [`crate::evidence`] separate held-out alleles from
//! everything else on a plate with known truth?
//!
//! The truth is known per sequence, not per donor: a held-out allele is
//! real in the donors that carry it and an error copy in others. The
//! comparison is therefore made per sequence, on the donor where its
//! evidence is strongest, which is also the unit of the allele recovery
//! check.
//!
//! The report is aggregate only: no donor identifiers and no sequences.

mod classify;
mod models;
mod rows;
mod sequences;

use crate::evidence::{Models, Observation};

pub use classify::{classify, Class, Truth};

/// The count and ratio thresholds of the run being compared against.
#[derive(Debug, Clone, Copy)]
pub struct Thresholds {
    pub min_count: u64,
    pub min_ratio: f64,
}

impl Thresholds {
    /// Both levels, as `discover blast` applies them.
    pub fn passes(&self, o: &Observation) -> bool {
        o.count >= self.min_count
            && o.ratio >= self.min_ratio
            && o.best_full_count >= self.min_count
            && o.best_full_ratio >= self.min_ratio
    }

    fn label(&self) -> String {
        format!("-c {} -f {}", self.min_count, self.min_ratio)
    }
}

/// An assessed plate with the truth about its sequences.
pub struct Plate<'a> {
    pub observations: &'a [Observation],
    pub classes: &'a [Class],
    pub n_donors: usize,
    pub models: &'a Models,
    pub thresholds: Thresholds,
}

impl<'a> Plate<'a> {
    fn rows(&self) -> impl Iterator<Item = (&'a Observation, Class)> + '_ {
        self.observations.iter().zip(self.classes.iter().copied())
    }

    fn of_class(&self, class: Class) -> impl Iterator<Item = &'a Observation> + '_ {
        self.rows()
            .filter(move |(_, c)| *c == class)
            .map(|(o, _)| o)
    }
}

/// Six characters wide, three significant digits.
fn number(x: f64) -> String {
    if x.is_nan() {
        "   n/a".to_string()
    } else if x == 0.0 {
        "     0".to_string()
    } else if x != 0.0 && (x.abs() < 1e-3 || x.abs() >= 1e5) {
        format!("{x:6.1e}")
    } else {
        let digits = (2 - x.abs().max(1e-9).log10().floor() as i32).clamp(0, 6) as usize;
        format!("{x:6.digits$}")
    }
}

fn optional(x: Option<f64>) -> String {
    number(x.unwrap_or(f64::NAN))
}

/// The evidence report for one plate; `list` bounds the two sequence lists.
pub fn format_evidence_report(plate: &Plate, list: usize) -> String {
    let sequences = sequences::summarise(plate);
    let mut out = String::new();
    models::write_input(&mut out, plate);
    models::write_amplification(&mut out, plate);
    models::write_errors(&mut out, plate);
    models::write_evidence_scale(&mut out, plate);
    models::write_gene_units(&mut out, plate);
    sequences::write_overview(&mut out, plate, &sequences);
    sequences::write_separation(&mut out, plate, &sequences);
    sequences::write_rules(&mut out, plate, &sequences);
    sequences::write_cost(&mut out, plate, &sequences);
    sequences::write_distances(&mut out, plate, &sequences);
    rows::write_passing_rows(&mut out, plate);
    sequences::write_strongest_others(&mut out, &sequences, list);
    sequences::write_weakest_targets(&mut out, &sequences, list);
    out
}
