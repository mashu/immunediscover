//! Calling a diploid genotype from allele read counts of one gene.

use crate::frame::round_digits;

/// Genotype class of a gene in one donor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Zygosity {
    Homozygous,
    Heterozygous,
    /// More than two alleles reach the ratio threshold.
    Duplication,
}

impl Zygosity {
    pub fn as_str(self) -> &'static str {
        match self {
            Zygosity::Homozygous => "homozygous",
            Zygosity::Heterozygous => "heterozygous",
            Zygosity::Duplication => "duplication",
        }
    }
}

/// Genotype call for one donor and gene.
#[derive(Debug, Clone, PartialEq)]
pub struct GenotypeCall {
    pub zygosity: Zygosity,
    pub allele_1: String,
    /// Empty for homozygous calls.
    pub allele_2: String,
    pub count_1: i64,
    /// Count of the second most frequent allele, also when it is not called.
    pub count_2: i64,
    pub freq_1: f64,
    pub freq_2: f64,
    /// `count_2 / count_1`.
    pub ratio: f64,
    pub total_count: i64,
    /// Comma-separated further alleles.
    pub other_alleles: String,
    /// Positions (in the input order) of the two most frequent alleles.
    pub top: [Option<usize>; 2],
}

/// Call a genotype from `(label, count)` pairs of one gene.
///
/// Alleles are ranked by count (ties keep input order). An allele is
/// "effective" when its count is at least `min_ratio` of the top count.
pub fn call_genotype(alleles: &[(String, i64)], min_ratio: f64) -> GenotypeCall {
    let mut order: Vec<usize> = (0..alleles.len()).collect();
    order.sort_by(|&a, &b| alleles[b].1.cmp(&alleles[a].1));
    let counts: Vec<i64> = order.iter().map(|&i| alleles[i].1).collect();
    let labels: Vec<&str> = order.iter().map(|&i| alleles[i].0.as_str()).collect();
    let total: i64 = counts.iter().sum();

    let top = counts.first().copied().unwrap_or(0);
    let ratios: Vec<f64> = counts
        .iter()
        .map(|&c| if top > 0 { c as f64 / top as f64 } else { 0.0 })
        .collect();
    let effective = ratios.iter().filter(|r| **r >= min_ratio).count();
    let zygosity = match effective {
        0 | 1 => Zygosity::Homozygous,
        2 => Zygosity::Heterozygous,
        _ => Zygosity::Duplication,
    };

    let label = |i: usize| labels.get(i).map_or(String::new(), |l| l.to_string());
    let count_1 = counts.first().copied().unwrap_or(0);
    let count_2 = counts.get(1).copied().unwrap_or(0);
    let others: Vec<&str> = match zygosity {
        Zygosity::Duplication => (2..labels.len())
            .filter(|&i| ratios[i] >= min_ratio)
            .map(|i| labels[i])
            .collect(),
        _ => labels.iter().skip(2).copied().collect(),
    };
    let share = |count: i64| {
        if total > 0 {
            count as f64 / total as f64
        } else {
            0.0
        }
    };

    GenotypeCall {
        zygosity,
        allele_1: label(0),
        allele_2: if zygosity == Zygosity::Homozygous {
            String::new()
        } else {
            label(1)
        },
        count_1,
        count_2,
        freq_1: round_digits(share(count_1), 3),
        freq_2: round_digits(share(count_2), 3),
        ratio: round_digits(
            if count_1 > 0 {
                count_2 as f64 / count_1 as f64
            } else {
                0.0
            },
            3,
        ),
        total_count: total,
        other_alleles: others.join(","),
        top: [order.first().copied(), order.get(1).copied()],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alleles(counts: &[i64]) -> Vec<(String, i64)> {
        counts
            .iter()
            .enumerate()
            .map(|(i, c)| (format!("a{i}"), *c))
            .collect()
    }

    #[test]
    fn homozygous_keeps_second_count_but_not_second_allele() {
        let call = call_genotype(&alleles(&[100, 5]), 0.1);
        assert_eq!(call.zygosity, Zygosity::Homozygous);
        assert_eq!((call.allele_1.as_str(), call.allele_2.as_str()), ("a0", ""));
        assert_eq!((call.count_1, call.count_2), (100, 5));
        assert_eq!(call.ratio, 0.05);
        assert_eq!(call.freq_1, 0.952);
    }

    #[test]
    fn heterozygous_lists_minor_alleles_as_other() {
        let call = call_genotype(&alleles(&[8, 10, 1]), 0.2);
        assert_eq!(call.zygosity, Zygosity::Heterozygous);
        assert_eq!(
            (call.allele_1.as_str(), call.allele_2.as_str()),
            ("a1", "a0")
        );
        assert_eq!(call.other_alleles, "a2");
        assert_eq!(call.top, [Some(1), Some(0)]);
    }

    #[test]
    fn duplication_lists_only_effective_other_alleles() {
        let call = call_genotype(&alleles(&[10, 9, 8, 1]), 0.5);
        assert_eq!(call.zygosity, Zygosity::Duplication);
        assert_eq!(call.other_alleles, "a2");
        assert_eq!(call.total_count, 28);
    }

    #[test]
    fn ties_keep_input_order() {
        let call = call_genotype(&alleles(&[7, 7]), 0.1);
        assert_eq!(
            (call.allele_1.as_str(), call.allele_2.as_str()),
            ("a0", "a1")
        );
    }
}
