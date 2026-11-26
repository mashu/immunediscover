//! `search bwa`: keep candidate sequences whose best genome alignment lies
//! on the expected chromosome.

pub mod aligner;
pub mod cigar;

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use self::aligner::{align_all, GenomeAlignment};
use self::cigar::Cigar;
use crate::align::{hamming, levenshtein};
use crate::commands::preprocess::pattern::compile;
use crate::error::Result;
use crate::frame::{Column, Frame};
use crate::info;
use crate::io::fasta::read_fasta;
use crate::io::path::always_gz;
use crate::io::tsv::{read_tsv, write_tsv, write_tsv_gz};
use crate::seq::dna::reverse_complement;

#[derive(Debug, Clone)]
pub struct BwaParams {
    pub tsv: PathBuf,
    pub output: PathBuf,
    /// FASTA files that have a BWA index next to them.
    pub genomes: Vec<PathBuf>,
    /// Text the reference description must contain.
    pub chromosome: String,
    /// Column with allele names.
    pub name_column: String,
    /// Columns concatenated into the sequence to align.
    pub sequence_columns: Vec<String>,
    /// Regex the reference description must match.
    pub tag: String,
}

/// Reference sequences of one genome FASTA.
struct Genome {
    path: PathBuf,
    sequences: HashMap<String, String>,
    /// Header text after the sequence name.
    descriptions: HashMap<String, String>,
}

impl Genome {
    fn load(path: &Path) -> Result<Self> {
        info!("Loading reference sequences from {}", path.display());
        let mut sequences = HashMap::new();
        let mut descriptions = HashMap::new();
        for record in read_fasta(path)? {
            let name = record.identifier().to_string();
            descriptions.insert(name.clone(), record.comment().to_string());
            sequences.insert(name, record.sequence);
        }
        info!("Loaded {} reference sequences", sequences.len());
        Ok(Genome {
            path: path.to_path_buf(),
            sequences,
            descriptions,
        })
    }

    fn description(&self, reference: &str) -> &str {
        self.descriptions.get(reference).map_or("", String::as_str)
    }

    /// Up to `length` bases starting at 0-based `position`; the part of the
    /// window that falls outside the sequence is left out.
    fn window(&self, reference: &str, position: i64, length: usize) -> &str {
        let Some(sequence) = self.sequences.get(reference) else {
            return "";
        };
        let end = (position + length as i64).clamp(0, sequence.len() as i64) as usize;
        let start = position.clamp(0, sequence.len() as i64) as usize;
        &sequence[start.min(end)..end]
    }
}

/// Per-sequence outcome of the genome check.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct GenomeCheck {
    pub on_target: bool,
    pub position: String,
    pub edit_distance: i64,
    pub ref_sequence: String,
    pub orientation: String,
    pub cigar: String,
}

/// A candidate whose best alignment lies outside the target.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Discarded {
    pub name: String,
    pub genome: PathBuf,
    /// Description of the reference sequence it aligned to instead.
    pub description: String,
}

/// Fill `check` from the accepted alignment.
///
/// `ref_sequence` is the stretch of the genome the *whole* candidate lies
/// on: the aligned span widened by the soft-clipped ends, in the
/// candidate's orientation and in upper case. `edit_distance` compares the candidate with
/// that stretch (Hamming when the lengths agree, Levenshtein otherwise), so
/// clipped ends count as differences instead of being ignored.
fn describe_hit(genome: &Genome, sequence: &str, best: &GenomeAlignment, check: &mut GenomeCheck) {
    check.on_target = true;
    check.position = format!("{}:{}", best.reference, best.position);
    check.orientation = if best.reverse { "-" } else { "+" }.to_string();
    check.cigar = best.cigar.clone();
    let cigar = Cigar::parse(&best.cigar);
    let (start, length) = if cigar.is_empty() {
        (best.position, sequence.len())
    } else {
        (
            best.position - cigar.leading_clip() as i64,
            cigar.leading_clip() + cigar.reference_span() + cigar.trailing_clip(),
        )
    };
    // Soft-masked genomes mark repeats in lower case; case carries no
    // sequence information and must not count as a difference.
    let window = genome
        .window(&best.reference, start, length)
        .to_ascii_uppercase();
    check.ref_sequence = if best.reverse {
        // Sequences with non-nucleotide characters are passed through as is.
        reverse_complement(&window).unwrap_or(window)
    } else {
        window
    };
    let candidate = sequence.to_ascii_uppercase();
    let (candidate, reference) = (candidate.as_bytes(), check.ref_sequence.as_bytes());
    check.edit_distance =
        hamming(candidate, reference).unwrap_or_else(|| levenshtein(candidate, reference)) as i64;
}

/// Align every sequence to every genome and decide whether it belongs to
/// the target chromosome.
///
/// The alignments of all genomes are pooled: a sequence is on target when
/// one of its top-scoring alignments, over every genome given, lies on a
/// reference whose description contains `chromosome` and matches `tag`.
/// The reported columns describe the first such alignment.
pub fn check_sequences(
    genomes: &[PathBuf],
    names: &[String],
    sequences: &[String],
    chromosome: &str,
    tag: &str,
) -> Result<(Vec<GenomeCheck>, Vec<Discarded>)> {
    let tag_regex = compile(tag)?;
    let loaded: Vec<Genome> = genomes
        .iter()
        .map(|p| Genome::load(p))
        .collect::<Result<_>>()?;
    let alignments: Vec<Vec<Vec<GenomeAlignment>>> = genomes
        .iter()
        .map(|genome| align_all(genome, sequences))
        .collect::<Result<_>>()?;

    let mut checks = vec![
        GenomeCheck {
            edit_distance: -1,
            ..GenomeCheck::default()
        };
        sequences.len()
    ];
    let mut discarded = Vec::new();

    for (i, (name, sequence)) in names.iter().zip(sequences).enumerate() {
        let pooled: Vec<(&Genome, &GenomeAlignment)> = loaded
            .iter()
            .zip(&alignments)
            .flat_map(|(genome, per_query)| per_query[i].iter().map(move |a| (genome, a)))
            .collect();
        let Some(best_score) = pooled.iter().map(|(_, a)| a.score).max() else {
            info!("{name} does **NOT** align to any genome (skipping)");
            continue;
        };
        let top: Vec<&(&Genome, &GenomeAlignment)> = pooled
            .iter()
            .filter(|(_, a)| a.score == best_score)
            .collect();
        let is_target = |(genome, a): &(&Genome, &GenomeAlignment)| -> bool {
            let description = genome.description(&a.reference);
            description.contains(chromosome) && tag_regex.is_match(description).unwrap_or(false)
        };
        match top.iter().find(|hit| is_target(hit)) {
            Some((genome, best)) => describe_hit(genome, sequence, best, &mut checks[i]),
            None => {
                let (genome, other) = top[0];
                let description = genome.description(&other.reference).to_string();
                info!(
                    "Discarded {name} matching {description} in {}",
                    genome.path.display()
                );
                discarded.push(Discarded {
                    name: name.clone(),
                    genome: genome.path.clone(),
                    description,
                });
            }
        }
    }
    Ok((checks, discarded))
}

/// `<output>.discarded.tsv` next to the result table.
fn discarded_path(output: &Path) -> PathBuf {
    let text = output.to_string_lossy();
    let base = text.strip_suffix(".gz").unwrap_or(&text);
    PathBuf::from(format!("{base}.discarded.tsv"))
}

/// Run the command. The result is gzip-compressed and written to `output`
/// with a `.gz` suffix; the candidates that were dropped are listed in
/// `<output>.discarded.tsv`.
pub fn run(params: &BwaParams) -> Result<()> {
    info!("BWA search to filter candidates if they match correct chromosome");
    let mut table = read_tsv(&params.tsv)?;
    let name_column = table.column(&params.name_column)?;
    let names: Vec<String> = name_column.iter().map(|v| v.julia_string()).collect();
    let sequence_columns: Vec<&Column> = params
        .sequence_columns
        .iter()
        .map(|c| table.column(c))
        .collect::<Result<_>>()?;
    let sequences: Vec<String> = (0..table.nrows())
        .map(|row| {
            sequence_columns
                .iter()
                // A missing cell contributes nothing to the sequence.
                .map(|c| c.get(row))
                .filter(|cell| !cell.is_missing())
                .map(|cell| cell.julia_string())
                .collect()
        })
        .collect();

    let (checks, discarded) = check_sequences(
        &params.genomes,
        &names,
        &sequences,
        &params.chromosome,
        &params.tag,
    )?;
    table.set_column(
        "position",
        Column::from_strings(checks.iter().map(|c| c.position.as_str())),
    )?;
    table.set_column(
        "edit_distance",
        Column::from_ints(checks.iter().map(|c| c.edit_distance)),
    )?;
    table.set_column(
        "ref_sequence",
        Column::from_strings(checks.iter().map(|c| c.ref_sequence.as_str())),
    )?;
    table.set_column(
        "orientation",
        Column::from_strings(checks.iter().map(|c| c.orientation.as_str())),
    )?;
    table.set_column(
        "cigar",
        Column::from_strings(checks.iter().map(|c| c.cigar.as_str())),
    )?;

    let mask: Vec<bool> = checks.iter().map(|c| c.on_target).collect();
    info!(
        "{} sequences processed, {} matched chromosome {}",
        table.nrows(),
        mask.iter().filter(|m| **m).count(),
        params.chromosome
    );
    let output = always_gz(&params.output);
    write_tsv_gz(&output, &table.filter(&mask))?;
    info!("Filtered result saved in {}", output.display());

    let discarded_file = discarded_path(&params.output);
    let genomes: Vec<String> = discarded
        .iter()
        .map(|d| d.genome.to_string_lossy().into_owned())
        .collect();
    let report = Frame::from_columns([
        (
            "name",
            Column::from_strings(discarded.iter().map(|d| d.name.as_str())),
        ),
        (
            "description",
            Column::from_strings(discarded.iter().map(|d| d.description.as_str())),
        ),
        (
            "genome",
            Column::from_strings(genomes.iter().map(String::as_str)),
        ),
    ])?;
    write_tsv(&discarded_file, &report)?;
    info!(
        "{} discarded sequence names written to {}",
        discarded.len(),
        discarded_file.display()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    //                 0         1         2         3
    //                 0123456789012345678901234567890123456789
    const REFERENCE: &str = "TTGACCGTAAGCTAGGCATCGATTACGGATCCAGTTACGA";

    fn genome() -> Genome {
        Genome {
            path: PathBuf::from("genome.fasta"),
            sequences: [("chr14".to_string(), REFERENCE.to_string())].into(),
            descriptions: [("chr14".to_string(), "chromosome 14".to_string())].into(),
        }
    }

    fn hit(position: i64, reverse: bool, cigar: &str) -> GenomeAlignment {
        GenomeAlignment {
            reference: "chr14".into(),
            position,
            reverse,
            cigar: cigar.into(),
            score: 0,
            edit_distance: 0,
        }
    }

    fn check(sequence: &str, alignment: &GenomeAlignment) -> GenomeCheck {
        let mut check = GenomeCheck::default();
        describe_hit(&genome(), sequence, alignment, &mut check);
        check
    }

    #[test]
    fn unclipped_hit_reports_the_aligned_reference() {
        let c = check(&REFERENCE[10..30], &hit(10, false, "20M"));
        assert_eq!(c.ref_sequence, &REFERENCE[10..30]);
        assert_eq!((c.edit_distance, c.cigar.as_str()), (0, "20M"));
        assert_eq!(c.position, "chr14:10");
    }

    #[test]
    fn soft_clipped_ends_are_compared_with_the_genome_not_skipped() {
        // Candidate = reference 10..30 with its first three bases replaced.
        let candidate = format!("AAA{}", &REFERENCE[13..30]);
        let c = check(&candidate, &hit(13, false, "3S17M"));
        // The window is widened by the clip, so it covers the whole candidate.
        assert_eq!(c.ref_sequence, &REFERENCE[10..30]);
        assert_eq!(c.edit_distance, 3);
        assert_eq!(c.cigar, "3S17M");
        // Trailing clip.
        let candidate = format!("{}CC", &REFERENCE[10..28]);
        let c = check(&candidate, &hit(10, false, "18M2S"));
        assert_eq!(c.ref_sequence, &REFERENCE[10..30]);
        assert_eq!(c.edit_distance, 2);
    }

    #[test]
    fn reverse_hits_are_reported_in_candidate_orientation() {
        let candidate = reverse_complement(&REFERENCE[10..30]).unwrap();
        let c = check(&candidate, &hit(10, true, "20M"));
        assert_eq!(c.ref_sequence, candidate);
        assert_eq!((c.edit_distance, c.orientation.as_str()), (0, "-"));
    }

    #[test]
    fn letter_case_is_not_a_difference() {
        let mut genome = genome();
        genome
            .sequences
            .insert("chr14".to_string(), REFERENCE.to_ascii_lowercase());
        for (reverse, candidate) in [
            (false, REFERENCE[10..30].to_string()),
            (true, reverse_complement(&REFERENCE[10..30]).unwrap()),
        ] {
            let mut check = GenomeCheck::default();
            describe_hit(&genome, &candidate, &hit(10, reverse, "20M"), &mut check);
            assert_eq!(check.edit_distance, 0);
            assert_eq!(check.ref_sequence, candidate);
        }
    }

    #[test]
    fn indels_use_edit_distance() {
        // One reference base (index 20) missing from the candidate.
        let candidate = format!("{}{}", &REFERENCE[10..20], &REFERENCE[21..30]);
        let c = check(&candidate, &hit(10, false, "10M1D9M"));
        assert_eq!(c.ref_sequence, &REFERENCE[10..30]);
        assert_eq!(c.edit_distance, 1);
    }

    #[test]
    fn window_is_cut_at_the_reference_end() {
        let candidate = format!("GGGG{}", &REFERENCE[..10]);
        let c = check(&candidate, &hit(0, false, "4S10M"));
        assert_eq!(c.ref_sequence, &REFERENCE[..10]);
        assert_eq!(c.edit_distance, 4);
    }

    #[test]
    fn discarded_list_sits_next_to_the_output() {
        assert_eq!(
            discarded_path(Path::new("out/v_bwa.tsv.gz")),
            PathBuf::from("out/v_bwa.tsv.discarded.tsv")
        );
        assert_eq!(
            discarded_path(Path::new("v_bwa.tsv")),
            PathBuf::from("v_bwa.tsv.discarded.tsv")
        );
    }
}
