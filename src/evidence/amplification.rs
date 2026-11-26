//! How many reads one copy of a gene gives in a donor.
//!
//! `ln(reads of gene g in donor d) = donor[d] + gene[g]`, fitted by median
//! polish, which a deleted or duplicated gene in a few donors does not
//! disturb. `exp(donor + gene) / copies[g]` is then the expectation for
//! ONE copy.
//!
//! The total of a gene is taken over its database alleles, and over the
//! main ones only (at least [`MAIN_RATIO`] of the best one in the donor).
//! Reads carrying errors, a sizeable share of a 300 nt amplicon, do not
//! enter the unit, and neither does a sequence of another locus that the
//! search assigned to the gene, however abundant.
//!
//! `copies[g]` is two unless the donors show otherwise: where a donor
//! holds two or more database alleles of a gene, the smallest of them is
//! one copy, and its share of the total says how many copies there are.
//! A gene that is duplicated in most donors has four, and one copy of it
//! is a quarter of its reads, not half.
//!
//! A gene without database alleles in enough donors is fitted on all its
//! main sequences. Where even that fails, the unit is half of what the
//! donor itself shows for the gene.

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
/// Donors with two or more database alleles a gene needs for its number
/// of copies to be estimated.
const MIN_RESOLVED_DONORS: usize = 5;
const MAX_COPIES: f64 = 8.0;

/// The fitted donor and gene factors (natural log scale).
#[derive(Debug, Clone)]
pub struct Amplification {
    pub donor: Vec<f64>,
    pub gene: Vec<f64>,
    pub genes: Vec<String>,
    /// Copies of every gene a donor usually holds.
    pub copies: Vec<f64>,
    /// Whether the gene was fitted on its database alleles.
    pub from_database: Vec<bool>,
    gene_index: HashMap<String, usize>,
    /// Reads per copy from the donor's own reads, per donor and gene, for
    /// the genes without a fitted factor.
    own_unit: Vec<Vec<f64>>,
    /// Residuals of the fitted cells.
    pub residuals: Vec<f64>,
}

/// Reads of the main sequences of every gene in every donor.
struct Mains {
    /// `[donor][gene]`: the sequences counted, most abundant first.
    counts: Vec<Vec<Vec<f64>>>,
}

impl Mains {
    fn of<'a>(
        observations: impl Iterator<Item = &'a Observation> + Clone,
        gene_index: &HashMap<String, usize>,
        n_donors: usize,
    ) -> Self {
        let mut best = vec![vec![0u64; gene_index.len()]; n_donors];
        for o in observations.clone() {
            let cell = &mut best[o.donor][gene_index[&o.gene]];
            *cell = (*cell).max(o.count);
        }
        let mut counts = vec![vec![Vec::new(); gene_index.len()]; n_donors];
        for o in observations {
            let g = gene_index[&o.gene];
            if o.count as f64 >= MAIN_RATIO * best[o.donor][g] as f64 {
                counts[o.donor][g].push(o.count as f64);
            }
        }
        Mains { counts }
    }

    fn total(&self, donor: usize, gene: usize) -> f64 {
        self.counts[donor][gene].iter().sum()
    }

    /// Donors in which the gene has enough reads to enter the fit.
    fn donors_with(&self, gene: usize) -> usize {
        (0..self.counts.len())
            .filter(|d| self.total(*d, gene) >= MIN_GENE_READS)
            .count()
    }

    /// Copies of a gene, from the donors holding two or more of its
    /// sequences: the smallest is one copy.
    fn copies(&self, gene: usize) -> f64 {
        let resolved: Vec<f64> = self
            .counts
            .iter()
            .map(|row| &row[gene])
            .filter(|cell| cell.len() >= 2 && cell.iter().sum::<f64>() >= MIN_GENE_READS)
            .map(|cell| cell.iter().sum::<f64>() / cell.iter().copied().fold(f64::MAX, f64::min))
            .collect();
        if resolved.len() < MIN_RESOLVED_DONORS {
            return 2.0;
        }
        median(&resolved).round().clamp(2.0, MAX_COPIES)
    }
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

        let known = Mains::of(
            observations.iter().filter(|o| o.is_database_allele()),
            &gene_index,
            n_donors,
        );
        let all = Mains::of(observations.iter(), &gene_index, n_donors);
        // Database alleles where enough donors have them.
        let by_known: Vec<bool> = (0..genes.len())
            .map(|g| known.donors_with(g) >= MIN_DONORS_PER_GENE)
            .collect();
        let total = |d: usize, g: usize| {
            if by_known[g] {
                known.total(d, g)
            } else {
                all.total(d, g)
            }
        };
        let cells: Vec<Vec<Option<f64>>> = (0..n_donors)
            .map(|d| {
                (0..genes.len())
                    .map(|g| {
                        Some(total(d, g))
                            .filter(|t| *t >= MIN_GENE_READS)
                            .map(f64::ln)
                    })
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
        let copies = (0..genes.len())
            .map(|g| if by_known[g] { known.copies(g) } else { 2.0 })
            .collect();
        let own_unit = (0..n_donors)
            .map(|d| (0..genes.len()).map(|g| all.total(d, g) / 2.0).collect())
            .collect();
        Amplification {
            donor: fit.row,
            gene: fit.column,
            genes,
            copies,
            from_database: by_known,
            gene_index,
            own_unit,
            residuals,
        }
    }

    /// Expected reads of one copy of the observation's gene in its donor.
    pub fn per_copy_reads(&self, o: &Observation) -> Option<f64> {
        let g = *self.gene_index.get(&o.gene)?;
        let fitted = (self.donor[o.donor] + self.gene[g]).exp() / self.copies[g];
        let reads = if fitted.is_finite() {
            fitted
        } else {
            self.own_unit[o.donor][g]
        };
        (reads.is_finite() && reads > 0.0).then_some(reads)
    }

    /// Whether the unit of the observation's gene comes from the plate.
    pub fn is_fitted(&self, o: &Observation) -> bool {
        self.gene_index
            .get(&o.gene)
            .is_some_and(|g| (self.donor[o.donor] + self.gene[*g]).is_finite())
    }

    pub fn set_dosage(&self, observations: &mut [Observation]) {
        for o in observations.iter_mut() {
            o.evidence.per_copy_reads = self.per_copy_reads(o);
            o.evidence.dosage = o
                .evidence
                .per_copy_reads
                .map(|reads| o.count as f64 / reads);
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

    /// Ten donors, one gene, 100 reads per copy.
    fn dosages(cell: impl Fn(usize) -> Vec<Observation>) -> Vec<Observation> {
        let mut observations: Vec<Observation> = (0..10).flat_map(cell).collect();
        // A second, ordinary gene keeps the donor factors apart from the gene's.
        observations.extend((0..10).map(|d| observation(d, "H", "CCCC", 200)));
        let model = Amplification::fit(&observations, 10);
        model.set_dosage(&mut observations);
        observations
    }

    fn unknown(donor: usize, gene: &str, sequence: &str, count: u64) -> Observation {
        let mut o = observation(donor, gene, sequence, count);
        o.db_distance = 5;
        o
    }

    #[test]
    fn a_foreign_sequence_does_not_set_the_unit() {
        // Every donor: the database allele twice, and ten times as many
        // reads of a sequence of another locus assigned to the gene.
        let observations = dosages(|d| {
            vec![
                observation(d, "G", "AAAA", 200),
                unknown(d, "G", "TTTT", 2000),
            ]
        });
        assert!((observations[0].evidence.dosage.unwrap() - 2.0).abs() < 1e-9);
        assert!((observations[1].evidence.dosage.unwrap() - 20.0).abs() < 1e-9);
    }

    #[test]
    fn one_copy_of_a_duplicated_gene_is_a_quarter_of_it() {
        // Four copies per donor: one allele three times and another once in
        // most donors, two and two in the others. Donor 0 also holds an
        // unknown allele on one more copy.
        let observations = dosages(|d| {
            let mut cell = if d % 3 == 0 {
                vec![
                    observation(d, "G", "AAAA", 200),
                    observation(d, "G", "AAAT", 200),
                ]
            } else {
                vec![
                    observation(d, "G", "AAAA", 300),
                    observation(d, "G", "AAAT", 100),
                ]
            };
            if d == 0 {
                cell.push(unknown(d, "G", "AATT", 100));
            }
            cell
        });
        let novel = observations.iter().find(|o| o.sequence == "AATT").unwrap();
        assert!((novel.evidence.dosage.unwrap() - 1.0).abs() < 1e-9);
        let dosage_of = |donor: usize, sequence: &str| {
            observations
                .iter()
                .find(|o| o.donor == donor && o.sequence == sequence)
                .and_then(|o| o.evidence.dosage)
                .unwrap()
        };
        assert!((dosage_of(1, "AAAA") - 3.0).abs() < 1e-9);
        assert!((dosage_of(1, "AAAT") - 1.0).abs() < 1e-9);
    }

    #[test]
    fn a_gene_seen_in_one_donor_takes_its_unit_from_that_donor() {
        let mut observations = dosages(|d| vec![observation(d, "G", "AAAA", 200)]);
        observations.push(unknown(4, "R", "GGGG", 180));
        observations.push(unknown(4, "R", "GGGT", 9));
        let model = Amplification::fit(&observations, 10);
        model.set_dosage(&mut observations);
        let rare = &observations[observations.len() - 2];
        assert!(!model.is_fitted(rare));
        assert_eq!(rare.evidence.dosage, Some(2.0));
        assert!(model.is_fitted(&observations[0]));
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
