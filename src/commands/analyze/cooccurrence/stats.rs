//! Pairwise co-occurrence statistics of alleles across donors.

use std::collections::{BTreeSet, HashMap};

use crate::error::Result;
use crate::frame::Frame;
use crate::stats::hypergeometric::upper_tail;

/// Presence of alleles in donors plus pairwise association matrices.
#[derive(Debug, Clone)]
pub struct CooccurrenceStats {
    /// Alleles in sorted order; index into the matrices.
    pub alleles: Vec<String>,
    /// Donors carrying each allele (sorted).
    pub donors_of: HashMap<String, BTreeSet<String>>,
    /// Number of donors.
    pub n_donors: usize,
    /// Phi coefficient (row-major, `n` x `n`).
    pub rho: Vec<f64>,
    pub jaccard: Vec<f64>,
    /// Number of donors carrying both alleles.
    pub support: Vec<usize>,
    /// Enrichment p-value.
    pub p_value: Vec<f64>,
}

impl CooccurrenceStats {
    pub fn n(&self) -> usize {
        self.alleles.len()
    }

    pub fn at(&self, i: usize, j: usize) -> usize {
        i * self.n() + j
    }

    /// `max(rho, 0)` as a similarity matrix.
    pub fn similarity(&self) -> Vec<f64> {
        self.rho.iter().map(|r| r.max(0.0)).collect()
    }
}

/// Phi (Matthews) coefficient of a 2x2 presence table.
pub fn phi_coefficient(n11: i64, n10: i64, n01: i64, n00: i64) -> f64 {
    let denominator = (n11 + n10) * (n11 + n01) * (n10 + n00) * (n01 + n00);
    if denominator <= 0 {
        return 0.0;
    }
    (n11 * n00 - n10 * n01) as f64 / (denominator as f64).sqrt()
}

/// Enrichment p-value for `n11` shared donors: the probability that two
/// alleles carried by `n_a` and `n_b` of `n_donors` donors share at least
/// `n11` donors by chance (upper hypergeometric tail, `n_a` successes,
/// `n_donors - n_a` failures, `n_b` draws).
pub fn enrichment_p_value(n_donors: usize, n_a: usize, n_b: usize, n11: usize) -> f64 {
    if n_donors == 0 || n_a == 0 || n_b == 0 || n11 == 0 {
        return 1.0;
    }
    if n11 > n_a.min(n_b) {
        return 0.0;
    }
    let failures = n_donors.saturating_sub(n_a);
    upper_tail(n_a as u64, failures as u64, n_b as u64, n11 as i64)
}

/// Build presence sets and association matrices.
///
/// Alleles carried by fewer than `min_donors` distinct donors are dropped
/// first.
pub fn compute_full_stats(
    table: &Frame,
    case_column: &str,
    allele_column: &str,
    min_donors: usize,
) -> Result<CooccurrenceStats> {
    let groups = table.group_by(&[allele_column])?;
    let all_cases = table.column(case_column)?;
    let mut keep = vec![false; table.nrows()];
    for rows in groups.iter() {
        let carriers: BTreeSet<String> = rows
            .iter()
            .map(|&row| all_cases.get(row).julia_string())
            .collect();
        if carriers.len() >= min_donors {
            for &row in rows {
                keep[row] = true;
            }
        }
    }
    let filtered = table.filter(&keep);

    let distinct_sorted = |column: &str| -> Result<Vec<String>> {
        let unique = filtered.select(&[column])?.unique().sorted_asc(&[column])?;
        Ok(unique
            .column(column)?
            .iter()
            .map(|v| v.julia_string())
            .collect())
    };
    let donors = distinct_sorted(case_column)?;
    let alleles = distinct_sorted(allele_column)?;

    let mut donors_of: HashMap<String, BTreeSet<String>> = alleles
        .iter()
        .map(|a| (a.clone(), BTreeSet::new()))
        .collect();
    let case_cells = filtered.column(case_column)?;
    let allele_cells = filtered.column(allele_column)?;
    for row in 0..filtered.nrows() {
        donors_of
            .get_mut(&allele_cells.get(row).julia_string())
            .expect("every allele was registered")
            .insert(case_cells.get(row).julia_string());
    }

    let n = alleles.len();
    let n_donors = donors.len();
    let mut rho = vec![0.0; n * n];
    let mut jaccard = vec![0.0; n * n];
    let mut support = vec![0usize; n * n];
    let mut p_value = vec![1.0; n * n];
    for i in 0..n {
        rho[i * n + i] = 1.0;
        jaccard[i * n + i] = 1.0;
        p_value[i * n + i] = 0.0;
    }
    for i in 0..n {
        let a = &donors_of[&alleles[i]];
        for j in (i + 1)..n {
            let b = &donors_of[&alleles[j]];
            let n11 = a.intersection(b).count();
            let n10 = a.len() - n11;
            let n01 = b.len() - n11;
            let n00 = n_donors as i64 - (n11 + n10 + n01) as i64;
            let r = phi_coefficient(n11 as i64, n10 as i64, n01 as i64, n00);
            let union = n11 + n10 + n01;
            let jac = if union == 0 {
                0.0
            } else {
                n11 as f64 / union as f64
            };
            let p = enrichment_p_value(n_donors, a.len(), b.len(), n11);
            for (x, y) in [(i, j), (j, i)] {
                rho[x * n + y] = r;
                jaccard[x * n + y] = jac;
                support[x * n + y] = n11;
                p_value[x * n + y] = p;
            }
        }
    }
    Ok(CooccurrenceStats {
        alleles,
        donors_of,
        n_donors,
        rho,
        jaccard,
        support,
        p_value,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frame::Column;

    fn toy() -> Frame {
        Frame::from_columns([
            (
                "case",
                Column::from_strings([
                    "D1", "D1", "D1", "D2", "D2", "D2", "D3", "D3", "D3", "D4", "D4", "D4",
                ]),
            ),
            (
                "db_name",
                Column::from_strings([
                    "A*01", "B*01", "C*01", "A*01", "B*01", "C*01", "A*01", "B*01", "D*01", "A*01",
                    "C*01", "D*01",
                ]),
            ),
        ])
        .unwrap()
    }

    #[test]
    fn presence_and_matrices() {
        let stats = compute_full_stats(&toy(), "case", "db_name", 1).unwrap();
        assert_eq!(stats.alleles, ["A*01", "B*01", "C*01", "D*01"]);
        assert_eq!(stats.n_donors, 4);
        // B and C share D1, D2; B also D3, C also D4.
        let bc = stats.at(1, 2);
        assert_eq!(stats.support[bc], 2);
        assert_eq!(stats.jaccard[bc], 0.5);
        assert!((stats.rho[bc] - (-1.0 / 3.0)).abs() < 1e-12);
        assert_eq!(stats.rho[stats.at(2, 1)], stats.rho[bc]);
        // A is present everywhere: no variance, phi defined as 0.
        assert_eq!(stats.rho[stats.at(0, 1)], 0.0);
    }

    #[test]
    fn phi_of_perfect_association() {
        assert_eq!(phi_coefficient(2, 0, 0, 2), 1.0);
        assert_eq!(phi_coefficient(0, 2, 2, 0), -1.0);
        assert_eq!(phi_coefficient(2, 0, 0, 0), 0.0);
    }

    #[test]
    fn enrichment_is_the_overlap_test() {
        // All five carriers of one allele among the five of the other, 20 donors:
        // exactly one of C(20, 5) = 15504 placements.
        let p = enrichment_p_value(20, 5, 5, 5);
        assert!((p - 1.0 / 15504.0).abs() < 1e-18, "{p}");
        // Reference values: ccdf(Hypergeometric(Ka, N - Ka, Kb), k - 1) in Distributions.jl.
        assert!((enrichment_p_value(50, 10, 10, 8) / 3.4560006468360775e-6 - 1.0).abs() < 1e-11);
        assert!((enrichment_p_value(50, 45, 45, 45) / 4.719741735732219e-7 - 1.0).abs() < 1e-11);
        // Symmetric in the two alleles; everything shared by necessity is not enrichment.
        assert!(
            (enrichment_p_value(30, 7, 12, 4) - enrichment_p_value(30, 12, 7, 4)).abs() < 1e-15
        );
        assert_eq!(enrichment_p_value(10, 10, 4, 4), 1.0);
    }

    #[test]
    fn min_donors_counts_donors_not_rows() {
        // E*01 occurs three times (three flank variants) in a single donor.
        let table = Frame::from_columns([
            (
                "case",
                Column::from_strings(["D1", "D2", "D3", "D1", "D1", "D1"]),
            ),
            (
                "db_name",
                Column::from_strings(["A*01", "A*01", "A*01", "E*01", "E*01", "E*01"]),
            ),
        ])
        .unwrap();
        let stats = compute_full_stats(&table, "case", "db_name", 3).unwrap();
        assert_eq!(stats.alleles, ["A*01"]);
    }
}
