//! `discover blast`: BLAST-based discovery of candidate alleles.
//!
//! Reads are assigned to (context-extended) database alleles with BLAST,
//! identical aligned segments are clustered per donor, the extension is
//! trimmed off again, and clusters passing the count, ratio, length and
//! distance filters are named: known alleles keep their name, anything
//! else receives a content-derived `_Sxxxx` name.

pub mod affixes;
pub mod assign;
pub mod evidence_columns;
pub mod flank_split;
pub mod locate;
pub mod runner;
pub mod trim;
pub mod workdir;

use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};

use rayon::prelude::*;

use self::affixes::ExtensionSettings;
pub use self::assign::{Aligner, AssignSettings, DEFAULT_BUILTIN_IDENTITY};
use self::trim::{trim_and_align, AlignmentStats};
use self::workdir::{content_digest, file_identity, Fingerprint, RunLock, WorkLayout};
use crate::commands::variants::{
    allele_counts, rank_variants, ratio_to_gene_max, sort_variants, VariantColumns,
};
use crate::error::{Error, Result};
use crate::filters::{GermlineFilter, MaxThreshold, MinStringLength, MinThreshold, NonNegative};
use crate::frame::{format_float, Column, Frame};
use crate::io::compress;
use crate::io::path::always_gz;
use crate::io::tsv::{read_tsv, write_tsv_gz};
use crate::seq::{load_database, unique_name, Allele};
use crate::{info, warn};

/// Parameters of the `discover blast` command (after preset resolution).
#[derive(Debug, Clone)]
pub struct BlastParams {
    pub input: PathBuf,
    pub fasta: PathBuf,
    pub output: PathBuf,
    /// FASTA of pseudogenes; hits to them are discarded.
    pub pseudo: Option<PathBuf>,
    pub min_full_count: i64,
    pub min_full_ratio: f64,
    pub subject_coverage: f64,
    pub max_distance: i64,
    /// Minimum length of the BLAST-aligned read segment (allele plus
    /// extension, before trimming).
    pub min_length: usize,
    pub edge: i64,
    /// Which program compares the reads with the database (with the
    /// arguments for `blastn`, if it is BLAST).
    pub aligner: Aligner,
    pub overwrite: bool,
    pub forward_extension: usize,
    pub reverse_extension: usize,
    /// Largest edit distance at which a read without an exact allele still
    /// supplies context to its nearest allele (0: exact matches only).
    pub extension_max_distance: usize,
    /// Extend an allele that occurs in no read with the context of the
    /// other alleles of its gene.
    pub gene_context: bool,
    pub min_quality: f64,
    /// Name non-identical sequences that are substrings of a known allele
    /// after that allele instead of as novel.
    pub name_substrings_as_known: bool,
    /// Keep clusters whose extension could not be trimmed off; they are
    /// exempt from the core coverage and distance filters.
    pub keep_failed: bool,
    pub min_core_coverage: f64,
    /// Maximum number of extension variants reported per allele.
    pub top: usize,
    pub verbose: bool,
    /// Add the evidence columns (copies, evidence, donors, ...) to the
    /// output.
    pub evidence: bool,
    pub work_dir: String,
}

/// Map a BLAST subject id to its database key; novel alleles map to their
/// base name (`TRGV2*01_S2223` -> `TRGV2*01`) when only that is known.
fn subject_to_db_key<'a>(subject: &'a str, keys: &HashSet<&str>) -> &'a str {
    if keys.contains(subject) {
        return subject;
    }
    if let Some((base, hash)) = subject.rsplit_once("_S") {
        let is_hash = !hash.is_empty() && hash.bytes().all(|b| b.is_ascii_digit());
        if is_hash && !base.is_empty() && keys.contains(base) {
            return base;
        }
    }
    subject
}

fn warn_if_short(label: &str, extension: usize) {
    if extension > 0 && extension < 7 {
        warn!("{label} extension {extension} is short and may lead to false positives");
    }
}

fn warn_if_strict(affix_length: usize, min_quality: f64) {
    if affix_length > 0 && affix_length <= 20 && min_quality > 0.5 {
        warn!(
            "Quality threshold {}% might be too strict for short affixes ({affix_length} nt). Consider lowering --minquality.",
            format_float(crate::frame::round_digits(min_quality * 100.0, 1))
        );
    }
}

/// Load a FASTA for BLAST: every allele is named by the first word of its
/// header, which is what BLAST reports as subject id.
fn load_subjects(path: &Path) -> Result<Vec<Allele>> {
    Ok(load_database(path)?
        .into_iter()
        .map(|a| {
            let identifier = a.name.split_whitespace().next().unwrap_or_default();
            Allele::new(identifier, a.sequence)
        })
        .collect())
}

/// Pseudogenes under their BLAST names (prefixed `P`), then the alleles.
fn blast_subjects(pseudo: &[Allele], db: &[Allele]) -> Vec<Allele> {
    pseudo
        .iter()
        .map(|a| Allele::new(format!("P{}", a.name), a.sequence.clone()))
        .chain(db.iter().cloned())
        .collect()
}

/// Write the BLAST database FASTA.
fn write_combined_fasta(path: &Path, subjects: &[Allele]) -> Result<()> {
    let mut out = compress::create(path, false)?;
    let records = subjects
        .iter()
        .map(|a| (a.name.as_str(), a.sequence.as_str()));
    for (name, sequence) in records {
        // The trailing blank after the name is part of the original format.
        writeln!(out, ">{name} \n{sequence}").map_err(|e| Error::io(path, e))?;
    }
    out.flush().map_err(|e| Error::io(path, e))
}

/// Trim the database extension off every cluster and add `aln_qseq`,
/// `aln_mismatch`, `db_length` and `corecov`.
fn trim_clusters(
    clusters: &mut Frame,
    db: &[Allele],
    layout: &WorkLayout,
    params: &BlastParams,
) -> Result<()> {
    info!("Loading affixes from {}", layout.affixes.display());
    let affixes = affixes::load_affixes(&layout.affixes)?;
    warn_if_strict(params.forward_extension, params.min_quality);
    warn_if_strict(params.reverse_extension, params.min_quality);

    let by_name: HashMap<&str, &str> = db
        .iter()
        .map(|a| (a.name.as_str(), a.sequence.as_str()))
        .collect();
    let keys: HashSet<&str> = by_name.keys().copied().collect();
    let subjects: Vec<String> = clusters
        .column("sseqid")?
        .iter()
        .map(|v| v.julia_string().trim().to_string())
        .collect();
    let segments: Vec<String> = clusters
        .column("qseq")?
        .iter()
        .map(|v| v.julia_string())
        .collect();
    let no_affix = (String::new(), String::new());

    let results: Vec<(trim::Trimmed, AlignmentStats, usize)> = subjects
        .par_iter()
        .zip(&segments)
        .map(|(subject, segment)| {
            let key = subject_to_db_key(subject, &keys);
            let reference = by_name.get(key).ok_or_else(|| {
                Error::invalid(format!(
                    "BLAST subject '{subject}' is not a name in the database FASTA"
                ))
            })?;
            let (prefix, suffix) = affixes.get(key).unwrap_or(&no_affix);
            let mut stats = AlignmentStats::default();
            let trimmed = trim_and_align(
                segment,
                prefix,
                suffix,
                reference,
                params.min_quality,
                subject,
                &mut stats,
            )?;
            Ok((trimmed, stats, reference.len()))
        })
        .collect::<Result<_>>()?;

    let mut stats = AlignmentStats::default();
    let mut sequences = Vec::with_capacity(results.len());
    let mut distances = Vec::with_capacity(results.len());
    let mut db_lengths = Vec::with_capacity(results.len());
    let mut coverage = Vec::with_capacity(results.len());
    for (trimmed, local, db_length) in results {
        stats.merge(local);
        coverage.push(if db_length > 0 {
            trimmed.sequence.len() as f64 / db_length as f64
        } else {
            0.0
        });
        db_lengths.push(db_length as i64);
        distances.push(trimmed.distance);
        sequences.push(trimmed.sequence);
    }
    clusters.set_column("aln_qseq", Column::from_strings(sequences))?;
    clusters.set_column("aln_mismatch", Column::from_ints(distances))?;
    clusters.set_column("db_length", Column::from_ints(db_lengths))?;
    clusters.set_column("corecov", Column::from_floats(coverage.iter().copied()))?;

    info!(
        "Alignment stats: {} attempts, {} prefix, {} suffix failures",
        stats.total_attempts, stats.prefix_failures, stats.suffix_failures
    );
    if params.verbose {
        for (subject, message) in &stats.failures {
            info!("  Failed [{subject}]: {message}");
        }
    }
    Ok(())
}

/// Drop clusters whose trimmed core covers too little of the database
/// allele.
///
/// A cluster whose trimming failed has no core; with `keep_failed` it is
/// passed on instead of being removed by the coverage threshold.
fn core_coverage_filter(clusters: &mut Frame, params: &BlastParams) -> Result<()> {
    let before = clusters.nrows();
    let coverage = clusters.floats("corecov")?;
    let failed = clusters.ints("aln_mismatch")?;
    let mask: Vec<bool> = coverage
        .iter()
        .zip(&failed)
        .map(|(c, distance)| {
            c.is_some_and(|c| c >= params.min_core_coverage)
                || (params.keep_failed && distance.is_some_and(|d| d < 0))
        })
        .collect();
    *clusters = clusters.filter(&mask);
    info!(
        "Core coverage filter (>= {}): kept {} / {before}",
        format_float(params.min_core_coverage),
        clusters.nrows()
    );
    Ok(())
}

/// Temporary column holding the sequence that identifies a cluster's
/// allele: the trimmed sequence, or the aligned segment where trimming
/// failed (such clusters stay on their own).
const ALLELE_KEY: &str = "__allele";

/// How `discover blast` tables name donor, gene and allele.
const COLUMNS: VariantColumns<'static> = VariantColumns {
    well: "well",
    case: "case",
    gene: "gene",
    name: "sseqid",
    sequence: ALLELE_KEY,
};

/// Add `count` and `ratio`: reads per trimmed sequence, summed over the
/// clusters that differ only in their extension, and that count relative
/// to the best allele of the gene. Same meaning as in `search exact`.
fn add_allele_counts(clusters: &mut Frame) -> Result<()> {
    let key: Vec<String> = clusters
        .strs("aln_qseq")?
        .iter()
        .zip(clusters.strs("qseq")?)
        .map(|(core, segment)| match core.as_deref() {
            Some(core) if !core.is_empty() => core.to_string(),
            _ => segment.clone().unwrap_or_default(),
        })
        .collect();
    clusters.set_column(ALLELE_KEY, Column::from_strings(key))?;
    let full_count: Vec<i64> = clusters
        .ints("full_count")?
        .iter()
        .map(|c| c.unwrap_or(0))
        .collect();
    let count = allele_counts(clusters, &COLUMNS, &full_count)?;
    let ratio = ratio_to_gene_max(clusters, &COLUMNS, &count)?;
    clusters.set_column("count", Column::from_ints(count))?;
    clusters.set_column("ratio", Column::from_floats(ratio))?;
    Ok(())
}

/// Order clusters by donor, gene and allele, number the extension variants
/// of every allele by support (`flank_index`) and keep the best `top`.
fn rank_clusters(clusters: &Frame, top: usize) -> Result<Frame> {
    let ranked = rank_variants(&sort_variants(clusters, &COLUMNS)?, &COLUMNS, top)?;
    let columns: Vec<String> = ranked
        .names()
        .iter()
        .filter(|name| name.as_str() != ALLELE_KEY)
        .cloned()
        .collect();
    ranked.select(&columns)
}

/// Name every remaining cluster.
fn allele_names(
    clusters: &Frame,
    db: &[Allele],
    name_substrings_as_known: bool,
) -> Result<Vec<String>> {
    let distance = clusters.ints("aln_mismatch")?;
    let subjects: Vec<String> = clusters
        .column("sseqid")?
        .iter()
        .map(|v| v.julia_string())
        .collect();
    let sequences: Vec<String> = clusters
        .column("aln_qseq")?
        .iter()
        .map(|v| v.julia_string())
        .collect();
    Ok((0..clusters.nrows())
        .map(|row| {
            if distance[row] == Some(0) {
                return subjects[row].clone();
            }
            if name_substrings_as_known {
                if let Some(known) = db
                    .iter()
                    .find(|a| a.sequence.contains(sequences[row].as_str()))
                {
                    return known.name.clone();
                }
            }
            unique_name(&subjects[row], &sequences[row])
        })
        .collect())
}

/// Run the command.
pub fn run(params: &BlastParams) -> Result<()> {
    info!(
        "Discovery with read assignment by {}",
        params.aligner.label()
    );
    let work_dir = workdir::resolve_work_dir(&params.work_dir)?;
    std::fs::create_dir_all(&work_dir).map_err(|e| Error::io(&work_dir, e))?;
    info!(
        "Work directory (caches, intermediates): {}",
        work_dir.display()
    );
    let layout = WorkLayout::new(
        &work_dir,
        &params.input,
        &params.fasta,
        params.pseudo.as_deref(),
    )?;
    layout.create_directories()?;
    // Held until the run ends: a second run on the same reads and
    // database would write the same intermediate files.
    let _lock = RunLock::acquire(&layout.run_dir)?;
    let db = load_subjects(&params.fasta)?;
    let all_reads = read_tsv(&params.input)?;
    warn_if_short("Forward", params.forward_extension);
    warn_if_short("Reverse", params.reverse_extension);

    let pseudo = match &params.pseudo {
        Some(path) => load_subjects(path)?,
        None => Vec::new(),
    };
    let subjects = blast_subjects(&pseudo, &db);
    write_combined_fasta(&layout.combined_fasta, &subjects)?;

    let extending = params.forward_extension != 0 || params.reverse_extension != 0;
    let blast_fasta = if extending {
        // The extended database depends on the reads, on every subject
        // sequence and on the extension lengths. Pseudogenes are extended
        // like the alleles, so that they compete for reads on equal terms.
        let fingerprint = Fingerprint::of(&[
            file_identity(&params.input)?,
            content_digest(&layout.combined_fasta)?,
            params.forward_extension.to_string(),
            params.reverse_extension.to_string(),
            params.extension_max_distance.to_string(),
            params.gene_context.to_string(),
        ]);
        let current = fingerprint.matches(&layout.extended_fasta) && layout.affixes.is_file();
        if !current || params.overwrite {
            info!(
                "Extending gene sequences by {} forward and {} reverse nucleotides",
                params.forward_extension, params.reverse_extension
            );
            let sequences: Vec<String> = all_reads
                .column("genomic_sequence")?
                .iter()
                .map(|v| v.julia_string())
                .collect();
            let extended = affixes::accumulate_affixes(
                &subjects,
                &sequences,
                &ExtensionSettings {
                    forward: params.forward_extension,
                    reverse: params.reverse_extension,
                    max_distance: params.extension_max_distance,
                    gene_fallback: params.gene_context,
                },
            )?;
            affixes::save_extended(&extended, &layout.extended_fasta, &layout.affixes)?;
            fingerprint.record(&layout.extended_fasta)?;
            info!("Saved affixes in {}", layout.affixes.display());
        } else {
            info!(
                "Using existing extended sequences from {}",
                layout.extended_fasta.display()
            );
        }
        &layout.extended_fasta
    } else {
        info!("No sequence extension requested, using original sequences");
        affixes::write_affixes(&layout.affixes, [("", "", "")].into_iter())?;
        &layout.combined_fasta
    };

    let mut clusters = assign::blast_discover(
        &params.input,
        &all_reads,
        blast_fasta,
        &layout,
        &AssignSettings {
            aligner: params.aligner.clone(),
            max_distance: params.max_distance,
            min_edge: params.edge,
            min_subject_coverage: params.subject_coverage,
            verbose: params.verbose,
            overwrite: params.overwrite,
        },
    )?;

    if extending {
        trim_clusters(&mut clusters, &db, &layout, params)?;
    } else {
        clusters.set_column("aln_qseq", clusters.column("qseq")?.clone())?;
        clusters.set_column("aln_mismatch", clusters.column("mismatch")?.clone())?;
    }
    add_allele_counts(&mut clusters)?;
    if extending {
        core_coverage_filter(&mut clusters, params)?;
    }

    let output = always_gz(&params.output);
    // Alleles that clear -c/-f on count/ratio but lose every extension to
    // full_count/full_ratio (support split across flanks).
    let flank_split = flank_split::find_flank_split_alleles(
        &clusters,
        params.min_full_count,
        params.min_full_ratio,
    )?;
    flank_split::report_flank_split_losses(
        &flank_split,
        params.min_full_count,
        params.min_full_ratio,
        &output,
    )?;

    if params.verbose {
        let unfiltered = flank_split::unfiltered_sidecar(&output);
        let export = clusters_without_allele_key(&clusters)?;
        write_tsv_gz(&unfiltered, &export)?;
        info!(
            "Unfiltered clusters (before abundance cuts) written to {}",
            unfiltered.display()
        );
    }

    let mut filter = GermlineFilter::new()
        .with(MinThreshold::new(
            "full_count",
            params.min_full_count as f64,
            "Min cluster size",
        ))
        .with(MinThreshold::new(
            "full_ratio",
            params.min_full_ratio,
            "Min allelic ratio",
        ))
        // The same thresholds at allele level, as in `search exact`.
        .with(MinThreshold::new(
            "count",
            params.min_full_count as f64,
            "Min allele count",
        ))
        .with(MinThreshold::new(
            "ratio",
            params.min_full_ratio,
            "Min allele ratio",
        ))
        .with(MinStringLength::new(
            "qseq",
            params.min_length,
            "Min read length",
        ));
    if !params.keep_failed {
        filter.push(NonNegative::new("aln_mismatch", "Exclude failed trimming"));
    }
    filter.push(MaxThreshold::new(
        "aln_mismatch",
        params.max_distance as f64,
        "Max distance",
    ));
    // The evidence is learned from everything the run saw, cut or not.
    let all_clusters = params.evidence.then(|| clusters.clone());
    filter.apply(&mut clusters)?;
    let mut clusters = rank_clusters(&clusters, params.top)?;

    let names = allele_names(&clusters, &db, params.name_substrings_as_known)?;
    clusters.set_column("allele_name", Column::from_strings(names))?;
    if let Some(all_clusters) = &all_clusters {
        info!(
            "Adding evidence columns: {}",
            evidence_columns::COLUMNS.join(", ")
        );
        evidence_columns::add_evidence_columns(&mut clusters, all_clusters)?;
    }

    info!(
        "Saving discovered gene sequences in compressed {} file",
        output.display()
    );
    write_tsv_gz(&output, &clusters)
}

/// Drop the temporary `__allele` column used for counting.
fn clusters_without_allele_key(clusters: &Frame) -> Result<Frame> {
    let columns: Vec<String> = clusters
        .names()
        .iter()
        .filter(|name| name.as_str() != ALLELE_KEY)
        .cloned()
        .collect();
    clusters.select(&columns)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn novel_subject_ids_map_to_their_base_allele() {
        let keys: HashSet<&str> = ["TRGV2*01", "IGHV1-2*02_S1234"].into_iter().collect();
        assert_eq!(subject_to_db_key("TRGV2*01", &keys), "TRGV2*01");
        assert_eq!(subject_to_db_key("TRGV2*01_S2223", &keys), "TRGV2*01");
        assert_eq!(
            subject_to_db_key("IGHV1-2*02_S1234", &keys),
            "IGHV1-2*02_S1234"
        );
        assert_eq!(subject_to_db_key("X*01_S12", &keys), "X*01_S12");
        assert_eq!(subject_to_db_key("TRGV2*01_Sab", &keys), "TRGV2*01_Sab");
    }

    /// Donor D1, subject A*01: the same trimmed sequence in three extension
    /// contexts (6, 3 and 2 reads), a second sequence (4 reads), and two
    /// clusters whose trimming failed.
    fn clusters() -> Frame {
        Frame::from_columns([
            ("well", Column::from_ints([1; 6])),
            ("case", Column::from_strings(["D1"; 6])),
            ("sseqid", Column::from_strings(["A*01"; 6])),
            (
                "qseq",
                Column::from_strings([
                    "ttACGTcc", "ttACGTca", "tgACGTcc", "ttACGAcc", "xxxx", "yyyy",
                ]),
            ),
            ("full_count", Column::from_ints([3, 6, 2, 4, 5, 1])),
            ("gene", Column::from_strings(["A"; 6])),
            (
                "aln_qseq",
                Column::from_strings(["ACGT", "ACGT", "ACGT", "ACGA", "", ""]),
            ),
            ("aln_mismatch", Column::from_ints([0, 0, 0, 1, -1, -1])),
        ])
        .unwrap()
    }

    #[test]
    fn count_is_the_allele_total_over_extension_variants() {
        let mut frame = clusters();
        add_allele_counts(&mut frame).unwrap();
        // Failed clusters have no trimmed sequence and are not pooled.
        assert_eq!(
            frame.ints("count").unwrap(),
            [11, 11, 11, 4, 5, 1].map(Some)
        );
        let ratio: Vec<f64> = frame
            .floats("ratio")
            .unwrap()
            .iter()
            .map(|r| r.unwrap())
            .collect();
        assert_eq!(ratio[0], 1.0);
        assert_eq!(ratio[3], 4.0 / 11.0);
    }

    #[test]
    fn best_supported_variant_represents_the_allele() {
        let mut frame = clusters();
        add_allele_counts(&mut frame).unwrap();
        let best = rank_clusters(&frame, 1).unwrap();
        assert!(!best.has_column(ALLELE_KEY));
        let segments = best.strs("qseq").unwrap();
        let segments: Vec<&str> = segments.iter().map(|s| s.as_deref().unwrap()).collect();
        // One row per sequence: ACGA, ACGT (its 6-read variant), and the
        // two failed clusters, each on its own.
        assert_eq!(segments, ["ttACGAcc", "ttACGTca", "xxxx", "yyyy"]);
        assert_eq!(best.ints("flank_index").unwrap(), [Some(1); 4]);
        let all = rank_clusters(&frame, 10).unwrap();
        assert_eq!(all.nrows(), 6);
        assert_eq!(
            all.ints("flank_index").unwrap(),
            [1, 1, 2, 3, 1, 1].map(Some)
        );
    }
}
