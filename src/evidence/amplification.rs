//! How many reads one copy of a gene gives in a donor.
//!
//! `ln(reads of gene g in donor d) = donor[d] + gene[g]`, fitted by median
//! polish, which a deleted or duplicated gene in a few donors does not
//! disturb. `exp(donor + gene) / 2` is then the expectation for ONE copy.
//!
//! The total of a gene is taken over its main sequences only (at least
//! [`MAIN_RATIO`] of the best one in the donor), so that reads carrying
//! errors, a sizeable share of a 300 nt amplicon, do not enter the unit.

use std::collections::HashMap;

use super::Observation;
use crate::stats::median;
use crate::stats::robust::median_polish;
use crate::stats::tails::{log10_normal_cdf, log10_poisson_cdf};

/// Share of the best sequence a sequence needs to count towards its gene.
pub const MAIN_RATIO: f64 = 0.2;
/// Reads a gene needs in a donor for the cell to enter the fit.
const MIN_GENE_READS: f64 = 20.0;
/// Donors a gene needs for its efficiency to be estimated.
const MIN_DONORS_PER_GENE: usize = 3;

/// The fitted donor and gene factors (natural log scale).
#[derive(Debug, Clone)]
pub struct Amplification {
    pub donor: Vec<f64>,
    pub gene: Vec<f64>,
    pub genes: Vec<String>,
    gene_index: HashMap<String, usize>,
    /// Residuals of the fitted cells.
    pub residuals: Vec<f64>,
}

impl Amplification {
    pub fn fit(observations: &[Observation], n_donors: usize) -> Self {
        let mut genes: Vec<String> = observations.iter().map(|o| o.gene.clone()).collect();
        genes.sort();
        genes.dedup();
        let gene_index: HashMap<String, usize> = genes
            .iter()
            .cloned()
            .enumerate()
            .map(|(i, g)| (g, i))
            .collect();

        let mut best = vec![vec![0u64; genes.len()]; n_donors];
        for o in observations {
            let cell = &mut best[o.donor][gene_index[&o.gene]];
            *cell = (*cell).max(o.count);
        }
        let mut totals = vec![vec![0.0f64; genes.len()]; n_donors];
        for o in observations {
            let g = gene_index[&o.gene];
            if o.count as f64 >= MAIN_RATIO * best[o.donor][g] as f64 {
                totals[o.donor][g] += o.count as f64;
            }
        }
        let cells: Vec<Vec<Option<f64>>> = totals
            .iter()
            .map(|row| {
                row.iter()
                    .map(|&t| (t >= MIN_GENE_READS).then(|| t.ln()))
                    .collect()
            })
            .collect();
        let mut fit = median_polish(&cells, genes.len(), 30);
        for (g, value) in fit.column.iter_mut().enumerate() {
            if cells.iter().filter(|row| row[g].is_some()).count() < MIN_DONORS_PER_GENE {
                *value = f64::NAN;
            }
        }
        let residuals = cells
            .iter()
            .zip(&fit.row)
            .flat_map(|(row, d)| {
                row.iter()
                    .zip(&fit.column)
                    .filter_map(move |(cell, g)| cell.map(|x| x - d - g))
            })
            .filter(|r| !r.is_nan())
            .collect();
        Amplification {
            donor: fit.row,
            gene: fit.column,
            genes,
            gene_index,
            residuals,
        }
    }

    /// Expected reads of one copy of the observation's gene in its donor.
    pub fn per_copy_reads(&self, o: &Observation) -> Option<f64> {
        let reads = (self.donor[o.donor] + self.gene[*self.gene_index.get(&o.gene)?]).exp() / 2.0;
        (reads.is_finite() && reads > 0.0).then_some(reads)
    }

    pub fn set_dosage(&self, observations: &mut [Observation]) {
        for o in observations.iter_mut() {
            o.evidence.dosage = self.per_copy_reads(o).map(|reads| o.count as f64 / reads);
        }
    }

    /// Fill in `log10_p_copy` given the allele-level noise.
    pub fn set_copy_tail(&self, observations: &mut [Observation], sigma: f64) {
        for o in observations.iter_mut() {
            o.evidence.log10_p_copy = self
                .per_copy_reads(o)
                .map(|reads| log10_copy_tail(o.count, reads, sigma));
        }
    }
}

/// Spread of the log dosage of database alleles around one or two copies.
pub fn allele_noise(observations: &[Observation]) -> f64 {
    let residuals: Vec<f64> = observations
        .iter()
        .filter(|o| o.is_database_allele() && o.count >= 10)
        .filter_map(|o| o.evidence.dosage)
        .filter(|d| (0.3..=3.0).contains(d))
        .map(|d| d.ln().abs().min((d / 2.0).ln().abs()))
        .collect();
    if residuals.is_empty() {
        return 0.3;
    }
    // The median of |x| for a normal x is 0.6745 standard deviations.
    (median(&residuals) / 0.6745).max(0.1)
}

/// log10 P(at most `k` reads | one copy): Poisson sampling around a
/// per-copy expectation that itself varies log-normally between alleles.
pub fn log10_copy_tail(k: u64, per_copy: f64, sigma: f64) -> f64 {
    if k > 300 {
        // Sampling noise is negligible next to sigma.
        let z = ((k as f64 + 0.5) / per_copy).ln() / (sigma * sigma + 1.0 / per_copy).sqrt();
        return log10_normal_cdf(z);
    }
    let step = 0.25;
    let total: f64 = (-24..=24)
        .map(|i| {
            let z = i as f64 * step;
            let weight = (-z * z / 2.0).exp() / (2.0 * std::f64::consts::PI).sqrt() * step;
            weight * 10f64.powf(log10_poisson_cdf(k, per_copy * (sigma * z).exp()))
        })
        .sum();
    if total > 0.0 {
        total.log10().min(0.0)
    } else {
        -300.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::testing::observation;

    #[test]
    fn dosage_is_in_copies_whatever_the_depth_and_the_primer() {
        // Three donors (depth x1, x2, x4), three genes (efficiency x1, x3, x0.5);
        // one copy gives 100 * depth * efficiency reads.
        let (depths, efficiencies) = ([1.0, 2.0, 4.0], [1.0, 3.0, 0.5]);
        let mut observations = Vec::new();
        for (d, depth) in depths.iter().enumerate() {
            for (g, efficiency) in efficiencies.iter().enumerate() {
                let one_copy = (100.0 * depth * efficiency) as u64;
                let gene = format!("G{g}");
                if (d, g) == (0, 0) {
                    // heterozygous, plus an artifact at 5 % of a copy
                    observations.push(observation(d, &gene, "AAAA", one_copy));
                    observations.push(observation(d, &gene, "AAAT", one_copy));
                    observations.push(observation(d, &gene, "AATT", one_copy / 20));
                } else {
                    observations.push(observation(d, &gene, "AAAA", 2 * one_copy));
                }
            }
        }
        let model = Amplification::fit(&observations, 3);
        model.set_dosage(&mut observations);
        let dosage: Vec<f64> = observations
            .iter()
            .map(|o| o.evidence.dosage.unwrap())
            .collect();
        assert!((dosage[0] - 1.0).abs() < 1e-9 && (dosage[1] - 1.0).abs() < 1e-9);
        assert!((dosage[2] - 0.05).abs() < 1e-9);
        assert!(dosage[3..].iter().all(|d| (d - 2.0).abs() < 1e-9));
        assert!(model.residuals.iter().all(|r| r.abs() < 1e-9));
    }

    #[test]
    fn copy_tail_is_small_only_well_below_one_copy() {
        // 100 reads expected per copy, 20 % spread between alleles.
        assert!(log10_copy_tail(100, 100.0, 0.2) > -0.5);
        assert!(log10_copy_tail(250, 100.0, 0.2) > -0.01);
        assert!(log10_copy_tail(50, 100.0, 0.2) < -2.5);
        assert!(log10_copy_tail(5, 100.0, 0.2) < -12.0);
        // Monotone in the count, and continuous across the approximation switch.
        assert!(log10_copy_tail(10, 100.0, 0.2) < log10_copy_tail(20, 100.0, 0.2));
        let (below, above) = (
            log10_copy_tail(300, 400.0, 0.2),
            log10_copy_tail(301, 400.0, 0.2),
        );
        assert!((below - above).abs() < 0.05, "{below} {above}");
    }
}
