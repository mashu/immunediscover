//! Database extension: the read context ("affixes") around known alleles.
//!
//! For every database allele the bases immediately upstream and downstream
//! of its exact occurrences in the reads are collected and reduced to a
//! consensus. The allele extended by this consensus context is what reads
//! are BLASTed against.
//!
//! The consensus is a majority vote per position, walking outwards from
//! the allele: a side is extended for as long as more than half of all
//! occurrences of the allele carry the same base at that distance. The
//! result does not depend on the order of the reads, and a few reads with
//! a sequencing error or a truncated context do not shorten it.
//! Occurrences on the reverse strand are counted in allele orientation.
//!
//! Which reads count as occurrences of an allele is decided by
//! [`super::locate`]: exact matches, or, with a tolerance, the nearest
//! allele of a read that matches none exactly.
//!
//! The context is collected per allele. An allele that occurs in no read
//! can take the consensus over all occurrences of the other alleles of its
//! gene instead (`gene_fallback`), so that it is searched with context like
//! its siblings rather than bare.

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use rayon::prelude::*;

use super::locate::Locator;
use crate::error::{Error, Result};
use crate::frame::{Column, Frame};
use crate::info;
use crate::io::compress;
use crate::io::tsv::{read_tsv, write_tsv};
use crate::seq::dna::{base_index, reverse_complement};
use crate::seq::{gene_of, Allele};

const CHUNK: usize = 50_000;

/// How the read context of the alleles is collected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExtensionSettings {
    /// Bases of context taken upstream of an allele.
    pub forward: usize,
    /// Bases of context taken downstream of an allele.
    pub reverse: usize,
    /// Largest edit distance at which a read without any exact allele is
    /// still assigned to its nearest allele(s); 0: exact occurrences only.
    pub max_distance: usize,
    /// Give an allele that occurs in no read the context of its gene.
    pub gene_fallback: bool,
}

/// An allele with its consensus read context.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExtendedAllele {
    pub name: String,
    /// `prefix + allele + suffix`.
    pub sequence: String,
    pub prefix: String,
    pub suffix: String,
}

/// Base counts per distance from the allele on one side.
#[derive(Debug, Clone, Default)]
struct SideCounts {
    /// `counts[k]`: A, C, G, T, other at distance `k + 1` from the allele.
    counts: Vec<[u32; 5]>,
}

impl SideCounts {
    /// Add one occurrence; `bases` runs outwards from the allele.
    fn add(&mut self, bases: impl Iterator<Item = u8>) {
        for (distance, base) in bases.enumerate() {
            if self.counts.len() <= distance {
                self.counts.push([0; 5]);
            }
            self.counts[distance][base_index(base).unwrap_or(4)] += 1;
        }
    }

    fn merge(&mut self, other: &SideCounts) {
        if self.counts.len() < other.counts.len() {
            self.counts.resize(other.counts.len(), [0; 5]);
        }
        for (mine, theirs) in self.counts.iter_mut().zip(&other.counts) {
            for (a, b) in mine.iter_mut().zip(theirs) {
                *a += b;
            }
        }
    }

    /// Majority bases, outwards from the allele, among `occurrences`.
    fn consensus(&self, occurrences: u32) -> Vec<u8> {
        let mut bases = Vec::new();
        for counts in &self.counts {
            let (best, support) = (0..4).fold((0, 0), |(best, support), base| {
                if counts[base] > support {
                    (base, counts[base])
                } else {
                    (best, support)
                }
            });
            if u64::from(support) * 2 <= u64::from(occurrences) {
                break;
            }
            bases.push(b"ACGT"[best]);
        }
        bases
    }
}

/// Running context statistics of one allele.
#[derive(Debug, Clone, Default)]
struct AffixConsensus {
    occurrences: u32,
    upstream: SideCounts,
    downstream: SideCounts,
}

impl AffixConsensus {
    /// `prefix` and `suffix` are the context in allele orientation.
    fn add(&mut self, prefix: &str, suffix: &str) {
        self.occurrences += 1;
        self.upstream.add(prefix.bytes().rev());
        self.downstream.add(suffix.bytes());
    }

    /// Pool the occurrences of another allele into this one.
    fn merge(&mut self, other: &AffixConsensus) {
        self.occurrences += other.occurrences;
        self.upstream.merge(&other.upstream);
        self.downstream.merge(&other.downstream);
    }

    fn prefix(&self) -> String {
        let mut bases = self.upstream.consensus(self.occurrences);
        bases.reverse();
        String::from_utf8(bases).expect("consensus bases are ASCII")
    }

    fn suffix(&self) -> String {
        String::from_utf8(self.downstream.consensus(self.occurrences))
            .expect("consensus bases are ASCII")
    }
}

/// Compute the consensus context of every allele over all reads.
pub fn accumulate_affixes(
    db: &[Allele],
    reads: &[String],
    settings: &ExtensionSettings,
) -> Result<Vec<ExtendedAllele>> {
    let (forward_extension, reverse_extension) = (settings.forward, settings.reverse);
    let max_distance = settings.max_distance;
    // Targets `0..n` are the alleles, `n..2n` their reverse complements
    // (an allele that cannot be complemented is only searched forward).
    let n = db.len();
    let targets: Vec<String> = db
        .iter()
        .map(|a| a.sequence.clone())
        .chain(
            db.iter()
                .map(|a| reverse_complement(&a.sequence).unwrap_or_else(|| a.sequence.clone())),
        )
        .collect();
    let locator = Locator::new(&targets, max_distance)?;
    let mut consensus: Vec<AffixConsensus> = vec![AffixConsensus::default(); n];
    let mut approximate = 0usize;

    for chunk in reads.chunks(CHUNK) {
        let per_read: Vec<Vec<(usize, String, String, bool)>> = chunk
            .par_iter()
            .map(|read| {
                locator
                    .locate(read)
                    .into_iter()
                    .filter_map(|hit| {
                        let before =
                            |length: usize| &read[hit.start.saturating_sub(length)..hit.start];
                        let after =
                            |length: usize| &read[hit.end..(hit.end + length).min(read.len())];
                        let inexact = hit.distance > 0;
                        if hit.target < n {
                            return Some((
                                hit.target,
                                before(forward_extension).to_string(),
                                after(reverse_extension).to_string(),
                                inexact,
                            ));
                        }
                        let allele = hit.target - n;
                        if targets[hit.target] == targets[allele] {
                            return None; // already counted as a forward match
                        }
                        // What follows the match in the read precedes the allele.
                        Some((
                            allele,
                            reverse_complement(after(forward_extension))?,
                            reverse_complement(before(reverse_extension))?,
                            inexact,
                        ))
                    })
                    .collect()
            })
            .collect();
        for hits in per_read {
            approximate += usize::from(hits.iter().any(|hit| hit.3));
            for (allele, prefix, suffix, _) in hits {
                consensus[allele].add(&prefix, &suffix);
            }
        }
    }
    if max_distance > 0 {
        info!(
            "{approximate} reads without an exact allele were assigned to their nearest allele (distance <= {max_distance})"
        );
    }
    let unseen = consensus.iter().filter(|c| c.occurrences == 0).count();
    info!(
        "Context found for {} of {n} alleles; {unseen} were not found in any read",
        n - unseen
    );
    if settings.gene_fallback && unseen > 0 {
        let borrowed = fill_from_gene(db, &mut consensus);
        info!(
            "{borrowed} of those take the context of the other alleles of their gene; {} stay unextended",
            unseen - borrowed
        );
    }

    Ok(db
        .iter()
        .zip(consensus)
        .map(|(allele, c)| {
            let (prefix, suffix) = (c.prefix(), c.suffix());
            ExtendedAllele {
                name: allele.name.clone(),
                sequence: format!("{prefix}{}{suffix}", allele.sequence),
                prefix,
                suffix,
            }
        })
        .collect())
}

/// Give every allele without occurrences the pooled statistics of the
/// alleles of its gene that have some. Returns how many were filled.
fn fill_from_gene(db: &[Allele], consensus: &mut [AffixConsensus]) -> usize {
    let mut by_gene: HashMap<&str, AffixConsensus> = HashMap::new();
    for (allele, seen) in db.iter().zip(consensus.iter()) {
        if seen.occurrences > 0 {
            by_gene
                .entry(gene_of(&allele.name))
                .or_default()
                .merge(seen);
        }
    }
    let mut filled = 0;
    for (allele, own) in db.iter().zip(consensus.iter_mut()) {
        if own.occurrences == 0 {
            if let Some(pooled) = by_gene.get(gene_of(&allele.name)) {
                *own = pooled.clone();
                filled += 1;
            }
        }
    }
    filled
}

/// Write the extended alleles as FASTA and their affixes as a TSV.
pub fn save_extended(extended: &[ExtendedAllele], fasta: &Path, affixes: &Path) -> Result<()> {
    let mut out = compress::create(fasta, false)?;
    for allele in extended {
        writeln!(out, ">{}\n{}", allele.name, allele.sequence).map_err(|e| Error::io(fasta, e))?;
    }
    out.flush().map_err(|e| Error::io(fasta, e))?;
    write_affixes(
        affixes,
        extended
            .iter()
            .map(|a| (a.name.as_str(), a.prefix.as_str(), a.suffix.as_str())),
    )
}

/// Write an affix table (`name`, `prefix`, `suffix`).
pub fn write_affixes<'a>(
    path: &Path,
    rows: impl Iterator<Item = (&'a str, &'a str, &'a str)>,
) -> Result<()> {
    let rows: Vec<(&str, &str, &str)> = rows.collect();
    let table = Frame::from_columns([
        ("name", Column::from_strings(rows.iter().map(|r| r.0))),
        ("prefix", Column::from_strings(rows.iter().map(|r| r.1))),
        ("suffix", Column::from_strings(rows.iter().map(|r| r.2))),
    ])?;
    write_tsv(path, &table)
}

/// Read an affix table into `name -> (prefix, suffix)`.
pub fn load_affixes(path: &Path) -> Result<HashMap<String, (String, String)>> {
    let table = read_tsv(path)?;
    let text = |column: &str, row: usize| -> Result<String> {
        let value = table.value(row, column)?;
        Ok(if value.is_missing() {
            String::new()
        } else {
            value.julia_string()
        })
    };
    let mut affixes = HashMap::new();
    for row in 0..table.nrows() {
        let name = text("name", row)?.trim().to_string();
        if name.is_empty() {
            continue;
        }
        affixes.insert(name, (text("prefix", row)?, text("suffix", row)?));
    }
    Ok(affixes)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn collect(reads: &[String], max_distance: usize) -> Vec<ExtendedAllele> {
        let settings = ExtensionSettings {
            forward: 5,
            reverse: 5,
            max_distance,
            gene_fallback: false,
        };
        accumulate_affixes(&db(), reads, &settings).unwrap()
    }

    fn db() -> Vec<Allele> {
        vec![
            Allele::new("A*01", "GATTACAGATTACA"),
            Allele::new("B*01", "CCCCGGGGTTTT"),
        ]
    }

    #[test]
    fn consensus_context_is_shared_part_adjacent_to_the_allele() {
        let reads = vec![
            "TTAAGCGATTACAGATTACATTGCAA".to_string(),
            "GGAAGCGATTACAGATTACATTGGGG".to_string(),
        ];
        let out = collect(&reads, 0);
        assert_eq!(out[0].prefix, "AAGC");
        assert_eq!(out[0].suffix, "TTG");
        assert_eq!(out[0].sequence, "AAGCGATTACAGATTACATTG");
        // Allele never observed: unchanged, empty affixes.
        assert_eq!(out[1].sequence, "CCCCGGGGTTTT");
        assert_eq!((out[1].prefix.as_str(), out[1].suffix.as_str()), ("", ""));
    }

    #[test]
    fn consensus_is_a_majority_vote_independent_of_read_order() {
        let mut reads = vec![
            "ACGTTGATTACAGATTACATGCAT".to_string(),
            "ACGTTGATTACAGATTACATGCAT".to_string(),
            "ACGTTGATTACAGATTACATGCAT".to_string(),
            // One read with an error next to the allele on both sides.
            "ACGTAGATTACAGATTACAAGCAT".to_string(),
            // One read in which the allele starts at the first base.
            "GATTACAGATTACATGCAT".to_string(),
        ];
        let forward = collect(&reads, 0);
        assert_eq!(forward[0].prefix, "ACGTT");
        assert_eq!(forward[0].suffix, "TGCAT");
        reads.reverse();
        let backward = collect(&reads, 0);
        assert_eq!(forward, backward);
    }

    #[test]
    fn consensus_stops_where_no_base_has_a_majority() {
        let reads = vec![
            "AAAAAGATTACAGATTACACCCCC".to_string(),
            "TTTTAGATTACAGATTACACGGGG".to_string(),
        ];
        let out = collect(&reads, 0);
        assert_eq!(out[0].prefix, "A");
        assert_eq!(out[0].suffix, "C");
    }

    #[test]
    fn reverse_strand_occurrences_count_in_allele_orientation() {
        // Reverse complement of ACGTT-GATTACAGATTACA-TGCAT.
        let reads = vec!["ATGCATGTAATCTGTAATCAACGT".to_string(); 3];
        let out = collect(&reads, 0);
        assert_eq!(out[0].prefix, "ACGTT");
        assert_eq!(out[0].suffix, "TGCAT");
    }

    #[test]
    fn tolerance_lets_reads_of_a_related_allele_supply_the_context() {
        // The donor's allele differs from A*01 by one base (C -> G at the
        // sixth position): no exact match anywhere.
        let reads = vec!["ACGTTGATTAGAGATTACATGCAT".to_string(); 3];
        let exact = collect(&reads, 0);
        assert_eq!(
            (exact[0].prefix.as_str(), exact[0].suffix.as_str()),
            ("", "")
        );
        let tolerant = collect(&reads, 1);
        assert_eq!(tolerant[0].prefix, "ACGTT");
        assert_eq!(tolerant[0].suffix, "TGCAT");
        // The database sequence itself is what gets extended.
        assert_eq!(tolerant[0].sequence, "ACGTTGATTACAGATTACATGCAT");
        assert_eq!(tolerant[1].sequence, "CCCCGGGGTTTT");
    }

    #[test]
    fn unseen_allele_takes_the_context_of_its_gene() {
        let db = vec![
            Allele::new("A*01", "GATTACAGATTACA"),
            Allele::new("A*02", "GATTACAGATTAGA"), // in no read
            Allele::new("A*03", "GATTACAGATTATA"),
            Allele::new("B*01", "CCCCGGGGTTTT"), // in no read, no sibling
        ];
        let reads = vec![
            "ACGTTGATTACAGATTACATGCAT".to_string(),
            "ACGTTGATTACAGATTACATGCAT".to_string(),
            "ACGTTGATTACAGATTATATGCAA".to_string(),
        ];
        let mut settings = ExtensionSettings {
            forward: 5,
            reverse: 5,
            max_distance: 0,
            gene_fallback: false,
        };
        let own = accumulate_affixes(&db, &reads, &settings).unwrap();
        assert_eq!(own[1].sequence, "GATTACAGATTAGA");
        settings.gene_fallback = true;
        let filled = accumulate_affixes(&db, &reads, &settings).unwrap();
        // Alleles with reads keep their own context.
        assert_eq!(filled[0], own[0]);
        assert_eq!(filled[2], own[2]);
        assert_eq!(filled[2].suffix, "TGCAA");
        // A*02 gets the majority over all three reads of gene A: the last
        // downstream base is T in two of them.
        assert_eq!(filled[1].prefix, "ACGTT");
        assert_eq!(filled[1].suffix, "TGCAT");
        assert_eq!(filled[1].sequence, "ACGTTGATTACAGATTAGATGCAT");
        // B has no allele with reads.
        assert_eq!(filled[3].sequence, "CCCCGGGGTTTT");
    }

    #[test]
    fn affix_table_round_trips_empty_values() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.affixes");
        write_affixes(&path, [("A*01 ", "ACG", ""), ("", "", "")].into_iter()).unwrap();
        let loaded = load_affixes(&path).unwrap();
        assert_eq!(loaded.len(), 1);
        assert_eq!(loaded["A*01"], ("ACG".to_string(), String::new()));
    }
}
