//! What every sequence of the plate is, according to the truth set.

use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::Path;

use crate::check::allele_recovery::normalize_allele;
use crate::error::Result;
use crate::evidence::Observation;
use crate::io::fasta::read_fasta;
use crate::seq::gene_of;

/// What a sequence is. Where a sequence has several classes (under two
/// genes), the first in this order wins.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Class {
    /// A database allele.
    Known,
    /// A held-out allele reported for this cohort: it should be found.
    Target,
    /// A held-out allele reported for other cohorts only: real, but maybe
    /// not on this plate.
    OtherTarget,
    /// A part of a database allele of the same gene.
    Partial,
    /// Anything else: artifacts, sequences of other loci, and real alleles
    /// nobody knows yet.
    Other,
}

impl Class {
    pub const ALL: [Class; 5] = [
        Class::Known,
        Class::Target,
        Class::OtherTarget,
        Class::Partial,
        Class::Other,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Class::Known => "database",
            Class::Target => "target",
            Class::OtherTarget => "other-cohort",
            Class::Partial => "partial",
            Class::Other => "other",
        }
    }
}

/// The database and the held-out alleles, by sequence.
#[derive(Debug, Default)]
pub struct Truth {
    known: HashSet<String>,
    known_by_gene: HashMap<String, Vec<String>>,
    /// Held-out sequence -> reported for this cohort.
    targets: HashMap<String, bool>,
}

impl Truth {
    pub fn add_known(&mut self, name: &str, sequence: &str) {
        let sequence = sequence.to_ascii_uppercase();
        self.known_by_gene
            .entry(gene_of(name).to_string())
            .or_default()
            .push(sequence.clone());
        self.known.insert(sequence);
    }

    pub fn add_target(&mut self, sequence: &str, in_cohort: bool) {
        *self
            .targets
            .entry(sequence.to_ascii_uppercase())
            .or_default() |= in_cohort;
    }

    /// `cohort_targets` are the normalised names of the held-out alleles
    /// reported for the cohort of the plate.
    pub fn read(db: &Path, targets: &Path, cohort_targets: &BTreeSet<String>) -> Result<Self> {
        let mut truth = Truth::default();
        for record in read_fasta(db)? {
            truth.add_known(record.identifier(), &record.sequence);
        }
        for record in read_fasta(targets)? {
            let name = normalize_allele(record.identifier());
            truth.add_target(&record.sequence, cohort_targets.contains(&name));
        }
        Ok(truth)
    }

    pub fn class_of(&self, gene: &str, sequence: &str) -> Class {
        if self.known.contains(sequence) {
            return Class::Known;
        }
        match self.targets.get(sequence) {
            Some(true) => Class::Target,
            Some(false) => Class::OtherTarget,
            None if self
                .known_by_gene
                .get(gene)
                .is_some_and(|alleles| alleles.iter().any(|a| a.contains(sequence))) =>
            {
                Class::Partial
            }
            None => Class::Other,
        }
    }
}

/// Classify every observation by exact sequence match.
pub fn classify(observations: &[Observation], truth: &Truth) -> Vec<Class> {
    let mut cache: HashMap<(&str, &str), Class> = HashMap::new();
    observations
        .iter()
        .map(|o| {
            *cache
                .entry((o.gene.as_str(), o.sequence.as_str()))
                .or_insert_with(|| truth.class_of(&o.gene, &o.sequence))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::testing::observation;

    #[test]
    fn classes_by_exact_sequence() {
        let mut truth = Truth::default();
        truth.add_known("IGHV1-2*01", "acgtacgtac");
        truth.add_target("ACGTACGTAA", false);
        truth.add_target("ACGTACGTAA", true); // also listed for this cohort
        truth.add_target("ACGTACGTTT", false);
        let observations: Vec<_> = [
            ("IGHV1-2", "ACGTACGTAC"),
            ("IGHV1-2", "ACGTACGTAA"),
            ("IGHV1-2", "ACGTACGTTT"),
            ("IGHV1-2", "CGTACG"),
            ("IGHV3-7", "CGTACG"),
            ("IGHV1-2", "TTTTTTTTTT"),
        ]
        .iter()
        .map(|(gene, sequence)| observation(0, gene, sequence, 10))
        .collect();
        assert_eq!(
            classify(&observations, &truth),
            [
                Class::Known,
                Class::Target,
                Class::OtherTarget,
                Class::Partial,
                Class::Other,
                Class::Other
            ]
        );
    }
}
