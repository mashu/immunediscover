//! Evidence that a sequence seen in a donor is a real germline allele.
//!
//! A read count alone says little: it has to be compared with what a real
//! allele would give in this donor and with what the known kinds of
//! artifact would give. This module turns the clusters of `discover blast`
//! into one [`Observation`] per sequence and donor and attaches that
//! comparison to it:
//!
//! * [`amplification`]: genomic DNA carries every gene in the same number
//!   of copies, so the reads of a gene in a donor are (amount of DNA) x
//!   (primer efficiency) x copies. Fitting the first two over the plate
//!   gives every sequence a **dosage** in copies, and the probability of
//!   seeing so few reads from one real copy.
//! * [`parents`], [`artifacts`]: the nearest more abundant sequence of the
//!   same gene in the same donor is what the sequence could be a copy of
//!   with errors. Two rates are learned from the plate: sporadic errors
//!   (any variant) and recurrent artifacts (this variant of this parent,
//!   at its usual ratio in the other donors). They give the expected
//!   number of artifacts at least this abundant.
//! * [`plate`]: how the sequence behaves across the donors of the plate.
//!
//! Results stay per donor: the rest of the plate is only used to learn
//! rates, as it is for the context extension.
//!
//! The combined [`score`] is provisional. It exists to find out, on plates
//! with known truth, whether these quantities separate real alleles from
//! artifacts better than count thresholds do.

pub mod amplification;
pub mod artifacts;
pub mod observations;
pub mod parents;
pub mod plate;
pub mod score;

pub use amplification::Amplification;
pub use artifacts::{RecurrentFit, SporadicFit};
pub use observations::{from_clusters, read_clusters};

/// One sequence in one donor, summed over its extension variants.
#[derive(Debug, Clone)]
pub struct Observation {
    /// Index of the donor (well + case) on the plate.
    pub donor: usize,
    pub gene: String,
    /// Database allele of the best BLAST hit.
    pub nearest: String,
    pub sequence: String,
    /// Reads with this sequence in any extension context.
    pub count: u64,
    /// Reads of its best-supported extension variant.
    pub best_full_count: u64,
    pub best_full_ratio: f64,
    /// `count` relative to the best sequence of the gene in the donor.
    pub ratio: f64,
    /// Differences from `nearest`; 0 for a database allele.
    pub db_distance: u64,
    pub evidence: Evidence,
}

impl Observation {
    pub fn is_database_allele(&self) -> bool {
        self.db_distance == 0
    }
}

/// What the models say about an [`Observation`]. `None` where a quantity
/// could not be computed.
#[derive(Debug, Clone, Default)]
pub struct Evidence {
    /// Estimated copies: 1 is one chromosome's worth of reads.
    pub dosage: Option<f64>,
    /// log10 P(at most this many reads | one real copy).
    pub log10_p_copy: Option<f64>,
    /// The more abundant sequence it could be an artifact of (index into
    /// the observations), and the differences from it.
    pub parent: Option<(usize, usize)>,
    /// log10 of the expected number of artifacts at least this abundant.
    pub log10_artifacts: Option<f64>,
    /// Starts like one more abundant sequence of the donor, ends like another.
    pub chimera: bool,
    /// Donors on the plate holding a copy-like amount of the sequence.
    pub carriers: usize,
    /// Share of the parent's carriers in which the sequence is seen at all.
    pub tracking: Option<f64>,
    /// Its usual ratio to the parent in the other donors (median).
    pub recurrent_ratio: Option<f64>,
    /// Provisional combined evidence, log10 units.
    pub score: Option<f64>,
}

/// The models fitted on a plate.
#[derive(Debug, Clone)]
pub struct Models {
    pub amplification: Amplification,
    /// Spread (SD of log) of database alleles around one or two copies.
    pub allele_noise: f64,
    /// Sporadic error fits for 1, 2 and 3 differences.
    pub sporadic: Vec<SporadicFit>,
    pub recurrent: RecurrentFit,
}

/// Fit the plate models and fill in the evidence of every observation.
pub fn assess(observations: &mut [Observation], n_donors: usize) -> Models {
    let amplification = Amplification::fit(observations, n_donors);
    amplification.set_dosage(observations);
    let allele_noise = amplification::allele_noise(observations);
    amplification.set_copy_tail(observations, allele_noise);
    parents::assign(observations);
    plate::set_context(observations);
    let sporadic = artifacts::fit_sporadic(observations);
    let recurrent = artifacts::fit_recurrent(observations);
    artifacts::set_expectation(observations, &sporadic, &recurrent);
    score::set(observations);
    Models {
        amplification,
        allele_noise,
        sporadic,
        recurrent,
    }
}

#[cfg(test)]
pub(crate) mod testing {
    use super::{Evidence, Observation};

    /// A database allele observed `count` times.
    pub fn observation(donor: usize, gene: &str, sequence: &str, count: u64) -> Observation {
        Observation {
            donor,
            gene: gene.to_string(),
            nearest: format!("{gene}*01"),
            sequence: sequence.to_string(),
            count,
            best_full_count: count,
            best_full_ratio: 1.0,
            ratio: 1.0,
            db_distance: 0,
            evidence: Evidence::default(),
        }
    }
}
