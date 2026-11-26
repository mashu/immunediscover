//! Synthetic plate libraries for tests and demonstrations.
//!
//! Produces a random reference allele, mutated "novel" variants of it, and
//! barcoded reads for two donors: the same construction as the original
//! `Simulate` module (the random streams themselves necessarily differ).

use std::path::Path;

use rand::seq::SliceRandom;
use rand::Rng;

use crate::error::{Error, Result};
use crate::frame::{Column, Frame};
use crate::io::fasta::{write_fasta, DEFAULT_LINE_WIDTH};
use crate::io::tsv::write_tsv;

const NUCLEOTIDES: [char; 4] = ['A', 'C', 'G', 'T'];
const RSS_SIGNAL: &str = "CACAGTG";
const LEADER_SIGNAL: &str = "GTTTTTGT";
const BARCODE_LENGTH: usize = 10;
const MIN_BARCODE_DISTANCE: usize = 4;
const EXTRA_VARIANTS_PER_DONOR: usize = 17;

/// Kind of mutation introduced into the reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mutation {
    Insertion,
    Deletion,
    Substitution,
}

impl Mutation {
    pub const ALL: [Mutation; 3] = [
        Mutation::Insertion,
        Mutation::Deletion,
        Mutation::Substitution,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Mutation::Insertion => "insertion",
            Mutation::Deletion => "deletion",
            Mutation::Substitution => "substitution",
        }
    }

    /// Apply the mutation at 1-based position `pos`.
    fn apply(self, rng: &mut impl Rng, seq: &str, pos: usize, length: usize) -> String {
        match self {
            Mutation::Insertion => {
                format!(
                    "{}{}{}",
                    &seq[..pos],
                    random_sequence(rng, length, length),
                    &seq[pos..]
                )
            }
            Mutation::Deletion => format!("{}{}", &seq[..pos - 1], &seq[pos - 1 + length..]),
            Mutation::Substitution => {
                let replaced: String = seq[pos - 1..pos - 1 + length]
                    .chars()
                    .map(|original| {
                        let choices: Vec<char> = NUCLEOTIDES
                            .iter()
                            .copied()
                            .filter(|n| *n != original)
                            .collect();
                        *choices.choose(rng).expect("three alternatives remain")
                    })
                    .collect();
                format!(
                    "{}{}{}",
                    &seq[..pos - 1],
                    replaced,
                    &seq[pos - 1 + length..]
                )
            }
        }
    }
}

fn random_sequence(rng: &mut impl Rng, min_length: usize, max_length: usize) -> String {
    let length = rng.gen_range(min_length..=max_length);
    (0..length)
        .map(|_| *NUCLEOTIDES.choose(rng).expect("alphabet is non-empty"))
        .collect()
}

fn hamming(a: &str, b: &str) -> usize {
    a.bytes().zip(b.bytes()).filter(|(x, y)| x != y).count()
}

fn distant_barcode(rng: &mut impl Rng, existing: &[String]) -> Result<String> {
    (0..1000)
        .map(|_| random_sequence(rng, BARCODE_LENGTH, BARCODE_LENGTH))
        .find(|candidate| {
            existing
                .iter()
                .all(|e| hamming(candidate, e) >= MIN_BARCODE_DISTANCE)
        })
        .ok_or_else(|| {
            Error::invalid(
                "Could not generate barcode with required minimum distance after 1000 attempts",
            )
        })
}

fn apply_random_mutation(
    rng: &mut impl Rng,
    seq: &str,
    mutation: Mutation,
    length: usize,
) -> String {
    let len = seq.len();
    let mut start = 50usize;
    let mut end = len.saturating_sub(50).min(len - length);
    if start >= end {
        start = (len / 4).max(1);
        end = (len - length).min(3 * (len / 4));
    }
    let pos = rng.gen_range(start..=end);
    mutation.apply(rng, seq, pos, length)
}

/// A read: barcodes around random padding, leader, allele and RSS signal.
fn barcoded_read(rng: &mut impl Rng, donor: &DonorIndex, core: &str) -> String {
    let prefix = random_sequence(rng, 20, 30);
    let suffix = random_sequence(rng, 20, 30);
    format!(
        "{}{}{}{}{}{}{}",
        donor.forward, prefix, LEADER_SIGNAL, core, RSS_SIGNAL, suffix, donor.reverse
    )
}

/// Donor barcodes written to the index sheet.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DonorIndex {
    pub forward: String,
    pub reverse: String,
    pub case: String,
}

/// Everything a simulation run produced.
#[derive(Debug, Clone)]
pub struct Simulation {
    /// `(read name, read sequence)`.
    pub reads: Vec<(String, String)>,
    pub indices: Vec<DonorIndex>,
    pub reference: String,
    /// `(allele name, sequence)` of every simulated variant.
    pub novel: Vec<(String, String)>,
}

/// Simulate a two-donor library and write reads, index sheet, reference and
/// variant FASTA files.
pub fn generate_library(
    rng: &mut impl Rng,
    reads_fasta: &Path,
    indices_tsv: &Path,
    reference_fasta: &Path,
    novel_fasta: &Path,
    n_reads: usize,
    base_length: usize,
) -> Result<Simulation> {
    let mut barcodes: Vec<String> = Vec::new();
    for _ in 0..4 {
        let barcode = distant_barcode(rng, &barcodes)?;
        barcodes.push(barcode);
    }
    let indices: Vec<DonorIndex> = (0..2)
        .map(|i| DonorIndex {
            forward: barcodes[i].clone(),
            reverse: barcodes[i + 2].clone(),
            case: format!("Donor{}", i + 1),
        })
        .collect();
    let sheet = Frame::from_columns([
        (
            "forward_index",
            Column::from_strings(indices.iter().map(|d| d.forward.as_str())),
        ),
        (
            "reverse_index",
            Column::from_strings(indices.iter().map(|d| d.reverse.as_str())),
        ),
        (
            "case",
            Column::from_strings(indices.iter().map(|d| d.case.as_str())),
        ),
    ])?;
    write_tsv(indices_tsv, &sheet)?;

    let reference = random_sequence(rng, base_length, base_length);
    write_fasta(
        reference_fasta,
        [("REF*01", reference.as_str())],
        Some(DEFAULT_LINE_WIDTH),
    )?;

    let lengths = [1usize, 3, 5];
    let mut reads: Vec<(String, String)> = Vec::new();
    let mut novel: Vec<(String, String)> = Vec::new();
    for (donor_number, donor) in indices.iter().enumerate() {
        for mutation in Mutation::ALL {
            for length in lengths {
                let mutated = apply_random_mutation(rng, &reference, mutation, length);
                novel.push((
                    format!("REF*{}_{}_{}", donor_number + 1, mutation.name(), length),
                    mutated.clone(),
                ));
                for read in 1..=n_reads {
                    let name =
                        format!("{}_{}_{}_read{}", donor.case, mutation.name(), length, read);
                    reads.push((name, barcoded_read(rng, donor, &mutated)));
                }
            }
        }
        for read in 1..=n_reads {
            let name = format!("{}_reference_read{}", donor.case, read);
            reads.push((name, barcoded_read(rng, donor, &reference)));
        }
    }
    for donor_number in 1..=indices.len() {
        for i in 1..=EXTRA_VARIANTS_PER_DONOR {
            let mutation = *Mutation::ALL
                .choose(rng)
                .expect("mutation list is non-empty");
            let length = *lengths.choose(rng).expect("length list is non-empty");
            let variant = apply_random_mutation(rng, &reference, mutation, length);
            novel.push((format!("REF*{donor_number}_extra_{i}"), variant));
        }
    }

    write_fasta(
        reads_fasta,
        reads.iter().map(|(n, s)| (n.as_str(), s.as_str())),
        Some(DEFAULT_LINE_WIDTH),
    )?;
    write_fasta(
        novel_fasta,
        novel.iter().map(|(n, s)| (n.as_str(), s.as_str())),
        Some(DEFAULT_LINE_WIDTH),
    )?;
    Ok(Simulation {
        reads,
        indices,
        reference,
        novel,
    })
}
