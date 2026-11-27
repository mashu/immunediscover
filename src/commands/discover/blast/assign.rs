//! BLAST assignment of reads to database alleles and clustering of the
//! aligned read segments.

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;
use std::time::Instant;

use super::runner::{blastn, ensure_blast_db, verify_blastn_version, COLUMNS};
use super::workdir::{content_digest, file_identity, replace_extension, Fingerprint, WorkLayout};
use crate::align::search::{Hit, Settings};
use crate::commands::align::{bit_score, e_value, load_database};
use crate::error::{Error, Result};
use crate::frame::{Column, Frame, Value};
use crate::io::compress;
use crate::io::tsv::{read_table, write_table, ReadOptions, WriteOptions};
use crate::seq::dna::reverse_complement;
use crate::seq::gene_of;
use crate::seq::search::find_first;
use crate::{info, warn};

/// The program that compares reads with the database.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Aligner {
    /// NCBI `blastn`, which must be installed, with extra arguments for it.
    Blast { args: String },
    /// The aligner that is part of this program ([`crate::align::search`]).
    Builtin,
}

impl Aligner {
    /// Name for messages.
    pub fn label(&self) -> &'static str {
        match self {
            Aligner::Blast { .. } => "BLAST",
            Aligner::Builtin => "the built-in aligner",
        }
    }
}

/// Thresholds of the assignment step.
#[derive(Debug, Clone)]
pub struct AssignSettings {
    pub aligner: Aligner,
    /// Maximum BLAST mismatches of a cluster.
    pub max_distance: i64,
    /// Minimum distance of the aligned segment from both read ends.
    pub min_edge: i64,
    /// Minimum alignment length relative to the subject length.
    pub min_subject_coverage: f64,
    /// Also write intermediate tables into the run directory.
    pub verbose: bool,
    /// Ignore a cached BLAST result.
    pub overwrite: bool,
}

fn text_cells(frame: &Frame, column: &str) -> Result<Vec<String>> {
    Ok(frame
        .column(column)?
        .iter()
        .map(|v| v.julia_string())
        .collect())
}

fn numbers(frame: &Frame, column: &str) -> Result<Vec<f64>> {
    frame
        .floats(column)?
        .into_iter()
        .map(|v| {
            v.ok_or_else(|| Error::invalid(format!("BLAST column {column} has an empty value")))
        })
        .collect()
}

fn write_query_fasta(path: &Path, reads: &Frame) -> Result<()> {
    let unique = reads
        .select(&["well", "case", "name", "genomic_sequence"])?
        .unique();
    let mut out = compress::create(path, false)?;
    for row in 0..unique.nrows() {
        let cell =
            |column: &str| -> Result<String> { Ok(unique.value(row, column)?.julia_string()) };
        writeln!(
            out,
            ">{} {} {}\n{}",
            cell("name")?,
            cell("well")?,
            cell("case")?,
            cell("genomic_sequence")?
        )
        .map_err(|e| Error::io(path, e))?;
    }
    out.flush().map_err(|e| Error::io(path, e))
}

/// Row index of the best hit of every query, queries in group order.
/// "Best" is the highest (pident, qcovhsp, qcovs, bitscore), first hit wins ties.
fn best_hits(hits: &Frame) -> Result<Vec<usize>> {
    let keys = [
        numbers(hits, "pident")?,
        numbers(hits, "qcovhsp")?,
        numbers(hits, "qcovs")?,
        numbers(hits, "bitscore")?,
    ];
    let better = |a: usize, b: usize| -> bool {
        keys.iter()
            .map(|k| k[a].total_cmp(&k[b]))
            .find(|o| o.is_ne())
            .is_some_and(|o| o.is_gt())
    };
    Ok(hits
        .group_by(&["qseqid"])?
        .iter()
        .map(|rows| {
            rows.iter()
                .copied()
                .reduce(|best, row| if better(row, best) { row } else { best })
                .expect("groups are non-empty")
        })
        .collect())
}

/// Number of read bases before and after the aligned segment.
fn edges(segment: &str, read: &str) -> Option<(i64, i64)> {
    find_first(segment, read).map(|(start, end)| (start as i64 - 1, (read.len() - end) as i64))
}

/// The aligned read segment in the orientation of the subject allele, and
/// the read bases on the allele's 5' and 3' side
/// (`five_prime_edge`, `three_prime_edge`).
///
/// BLAST reports `qseq` in read orientation; for a hit on the subject's
/// minus strand the segment is reverse-complemented and the two edges
/// swap. `None` if the segment is not part of the read.
fn oriented_segment(segment: &str, read: &str, minus_strand: bool) -> Option<(String, i64, i64)> {
    let (before, after) = edges(segment, read)?;
    if minus_strand {
        if let Some(complemented) = reverse_complement(segment) {
            return Some((complemented, after, before));
        }
    }
    Some((segment.to_string(), before, after))
}

/// Smallest share of equal bases at which the built-in aligner reports a
/// match. Low on purpose: `--maxdist` decides later what is too far.
const BUILTIN_MIN_IDENTITY: f64 = 0.7;

/// Best BLAST hit of every read: run `blastn` (or reuse its cached result)
/// and keep the best row per read. Columns as in [`COLUMNS`].
fn blast_best_hits(
    input_tsv: &Path,
    reads: &Frame,
    database_fasta: &Path,
    layout: &WorkLayout,
    blast_args: &str,
    settings: &AssignSettings,
) -> Result<Frame> {
    verify_blastn_version()?;
    // The cached hits depend on the reads, the database and the arguments.
    let fingerprint = Fingerprint::of(&[
        file_identity(input_tsv)?,
        content_digest(database_fasta)?,
        blast_args.to_string(),
    ]);
    let cached = layout.hits.is_file();
    let current = fingerprint.matches(&layout.hits);
    info!(
        "BLAST cache check: file='{}', exists={cached}, current={current}, overwrite={}",
        layout.hits.display(),
        settings.overwrite
    );
    if cached && !current && !settings.overwrite {
        info!(
            "Cached BLASTn results were made from other inputs or arguments; running BLASTn again."
        );
    }
    if current && !settings.overwrite {
        info!(
            "BLASTn results already exist {}. Skipping BLASTn.",
            layout.hits.display()
        );
    } else {
        info!("Running BLASTn. This may take a while.");
        let database = ensure_blast_db(database_fasta)?;
        let legacy = replace_extension(input_tsv, "blast");
        if legacy.is_file() {
            std::fs::remove_file(&legacy).map_err(|e| Error::io(&legacy, e))?;
            info!(
                "Removed legacy uncompressed BLAST file beside input: {}",
                legacy.display()
            );
        }
        // The reads are written out only when BLAST is going to read them,
        // and the result takes its place only once it is complete.
        write_query_fasta(&layout.query_fasta, reads)?;
        let partial = layout.run_dir.join("hits.partial.blast.gz");
        let searched = blastn(&layout.query_fasta, &database, &partial, blast_args);
        let _ = std::fs::remove_file(&layout.query_fasta);
        searched?;
        std::fs::rename(&partial, &layout.hits).map_err(|e| Error::io(&layout.hits, e))?;
        fingerprint.record(&layout.hits)?;
    }

    let hits = read_table(&layout.hits, &ReadOptions::tsv().header(&COLUMNS))?;
    info!(
        "BLASTn results read from {}: {} rows",
        layout.hits.display(),
        hits.nrows()
    );
    if hits.nrows() == 0 {
        return Err(Error::invalid(
            "No BLASTn results found (wrong BLAST parameters?). Cannot proceed.",
        ));
    }
    let best = hits.take(&best_hits(&hits)?);
    info!("BLASTn results after best hits: {} rows", best.nrows());
    Ok(best)
}

/// Best match of every read according to the built-in aligner, in the
/// columns BLAST would give ([`COLUMNS`]).
///
/// The whole database sequence (allele plus recorded surroundings) is
/// aligned inside the read, so `qseq` reaches from its first to its last
/// base unless the read ends before.
fn builtin_best_hits(
    reads: &Frame,
    database_fasta: &Path,
    settings: &AssignSettings,
) -> Result<Frame> {
    let database = load_database(database_fasta)?;
    let names = text_cells(reads, "name")?;
    let sequences = text_cells(reads, "genomic_sequence")?;
    let search = Settings {
        min_identity: BUILTIN_MIN_IDENTITY,
        min_coverage: settings.min_subject_coverage,
        ..Settings::default()
    };
    let started = Instant::now();
    let texts: Vec<&str> = sequences.iter().map(String::as_str).collect();
    let hits = database.search(&texts, &search);
    info!(
        "Built-in aligner: {} of {} reads matched ({} distinct read sequences, {:.1} seconds)",
        hits.matched(),
        hits.len(),
        hits.distinct(),
        started.elapsed().as_secs_f64()
    );
    if hits.matched() == 0 {
        return Err(Error::invalid(
            "No read matches a database sequence. Cannot proceed.",
        ));
    }
    let database_length: usize = database.subjects().iter().map(|s| s.sequence.len()).sum();
    let matched: Vec<(usize, &Hit)> = (0..hits.len())
        .filter_map(|read| hits.get(read).map(|hit| (read, hit)))
        .collect();
    let subject = |hit: &Hit| database.subject(hit.subject);
    let ints = |value: &dyn Fn(usize, &Hit) -> usize| {
        Column::from_ints(matched.iter().map(|&(read, hit)| value(read, hit) as i64))
    };
    let floats = |value: &dyn Fn(usize, &Hit) -> f64| {
        Column::from_floats(matched.iter().map(|&(read, hit)| value(read, hit)))
    };
    let read_coverage = |_: usize, hit: &Hit| {
        let (start, end) = hit.read_range();
        (100.0 * (end - start) as f64 / hit.read_length as f64).round()
    };
    Frame::from_columns([
        (
            "qseqid",
            Column::from_strings(matched.iter().map(|&(read, _)| names[read].as_str())),
        ),
        (
            "sseqid",
            Column::from_strings(matched.iter().map(|&(_, hit)| subject(hit).name.as_str())),
        ),
        ("pident", floats(&|_, hit| 100.0 * hit.identity())),
        ("nident", ints(&|_, hit| hit.path.matches())),
        ("length", ints(&|_, hit| hit.path.columns())),
        ("mismatch", ints(&|_, hit| hit.path.mismatches())),
        ("gapopen", ints(&|_, hit| hit.path.gap_openings())),
        ("qcovs", floats(&read_coverage)),
        ("qcovhsp", floats(&read_coverage)),
        ("qstart", ints(&|_, hit| hit.read_range().0 + 1)),
        ("qend", ints(&|_, hit| hit.read_range().1)),
        // On the other strand BLAST counts the subject positions downwards.
        (
            "sstart",
            ints(&|_, hit| {
                if hit.minus_strand {
                    hit.path.subject_end
                } else {
                    hit.path.subject_start + 1
                }
            }),
        ),
        (
            "send",
            ints(&|_, hit| {
                if hit.minus_strand {
                    hit.path.subject_start + 1
                } else {
                    hit.path.subject_end
                }
            }),
        ),
        ("qlen", ints(&|_, hit| hit.read_length)),
        ("slen", ints(&|_, hit| subject(hit).sequence.len())),
        (
            "evalue",
            floats(&|_, hit| e_value(hit.path.score, hit.read_length, database_length)),
        ),
        ("bitscore", floats(&|_, hit| bit_score(hit.path.score))),
        (
            "sstrand",
            Column::from_strings(matched.iter().map(|&(_, hit)| {
                if hit.minus_strand {
                    "minus"
                } else {
                    "plus"
                }
            })),
        ),
        (
            "qseq",
            Column::from_strings(
                matched
                    .iter()
                    .map(|&(read, hit)| hit.gapped_read(&sequences[read])),
            ),
        ),
    ])
}

/// Assign every read to a database sequence (with BLAST or the built-in
/// aligner) and build clusters of identical aligned segments per well,
/// donor and subject allele.
///
/// Returns columns `well`, `case`, `sseqid`, `qseq`, `mismatch`,
/// `full_count`, `gene`, `full_ratio`, sorted by well, case and subject.
pub fn blast_discover(
    input_tsv: &Path,
    all_reads: &Frame,
    database_fasta: &Path,
    layout: &WorkLayout,
    settings: &AssignSettings,
) -> Result<Frame> {
    layout.create_directories()?;

    info!(
        "Read {} rows from {} before assignment",
        all_reads.nrows(),
        input_tsv.display()
    );
    let reads = all_reads.unique_by(&["name"])?;
    info!("Unique reads: {}", reads.nrows());
    if reads.nrows() != all_reads.nrows() {
        warn!("Duplicated names found in input. Using unique ones but this is likely user error!");
    }

    let best = match &settings.aligner {
        Aligner::Blast { args } => {
            blast_best_hits(input_tsv, &reads, database_fasta, layout, args, settings)?
        }
        Aligner::Builtin => builtin_best_hits(&reads, database_fasta, settings)?,
    };

    // Attach read information.
    let read_row: HashMap<String, usize> = text_cells(&reads, "name")?
        .into_iter()
        .enumerate()
        .map(|(row, name)| (name, row))
        .collect();
    let rows: Vec<usize> = text_cells(&best, "qseqid")?
        .iter()
        .map(|id| {
            read_row.get(id).copied().ok_or_else(|| {
                Error::invalid(format!(
                    "Hit for unknown read '{id}'; stale cache? Re-run with --overwrite"
                ))
            })
        })
        .collect::<Result<_>>()?;
    let sequences = text_cells(&reads, "genomic_sequence")?;
    let segments = text_cells(&best, "qseq")?;
    let ungapped: Vec<String> = segments.iter().map(|s| s.replace('-', "")).collect();

    // Orient every segment like its subject allele; edge filter.
    let strands = text_cells(&best, "sstrand")?;
    let oriented: Vec<(String, i64, i64)> = ungapped
        .iter()
        .zip(&rows)
        .zip(&strands)
        .map(|((segment, &row), strand)| {
            oriented_segment(segment, &sequences[row], strand == "minus").ok_or_else(|| {
                Error::invalid(
                    "An aligned read segment was not found in its read; stale cache? Re-run with --overwrite",
                )
            })
        })
        .collect::<Result<_>>()?;
    let mut keep: Vec<bool> = oriented
        .iter()
        .map(|(_, five, three)| *five >= settings.min_edge && *three >= settings.min_edge)
        .collect();
    let count = |keep: &[bool]| keep.iter().filter(|k| **k).count();
    info!(
        "After filtering edge >= {}: {} rows",
        settings.min_edge,
        count(&keep)
    );

    // Subject coverage filter.
    let length = numbers(&best, "length")?;
    let subject_length = numbers(&best, "slen")?;
    for (row, flag) in keep.iter_mut().enumerate() {
        *flag &= length[row] / subject_length[row] >= settings.min_subject_coverage;
    }
    info!(
        "After filtering scov >= {}: {} rows",
        crate::frame::format_float(settings.min_subject_coverage),
        count(&keep)
    );

    // Pseudogene hits (database names prefixed with `P`).
    let subjects = text_cells(&best, "sseqid")?;
    for (row, flag) in keep.iter_mut().enumerate() {
        *flag &= !subjects[row].starts_with('P');
    }
    info!("After filtering pseudo genes: {} rows", count(&keep));

    let kept_rows: Vec<usize> = (0..best.nrows()).filter(|&r| keep[r]).collect();
    let read_rows: Vec<usize> = kept_rows.iter().map(|&r| rows[r]).collect();
    let mut assigned = Frame::from_columns([
        ("well", reads.column("well")?.take(&read_rows)),
        ("case", reads.column("case")?.take(&read_rows)),
        ("sseqid", best.column("sseqid")?.take(&kept_rows)),
        (
            "qseq",
            Column::from_strings(kept_rows.iter().map(|&r| oriented[r].0.as_str())),
        ),
        ("mismatch", best.column("mismatch")?.take(&kept_rows)),
    ])?;
    if settings.verbose {
        // Diagnostic table: best hits joined with their reads.
        let mut detail = best.take(&kept_rows);
        detail.set_column("qseq", assigned.column("qseq")?.clone())?;
        for name in reads.names().iter().filter(|n| n.as_str() != "name") {
            detail.set_column(name, reads.column(name)?.take(&read_rows))?;
        }
        detail.set_column(
            "five_prime_edge",
            Column::from_ints(kept_rows.iter().map(|&r| oriented[r].1)),
        )?;
        detail.set_column(
            "three_prime_edge",
            Column::from_ints(kept_rows.iter().map(|&r| oriented[r].2)),
        )?;
        detail.set_column(
            "scov",
            Column::from_floats(kept_rows.iter().map(|&r| length[r] / subject_length[r])),
        )?;
        write_table(
            &layout.run_dir.join("pseudo.tsv"),
            &detail,
            WriteOptions::tsv(),
        )?;
    }

    // Cluster identical segments.
    let by = ["well", "case", "sseqid", "qseq", "mismatch"];
    let groups = assigned.group_by(&by)?;
    let sizes: Vec<i64> = groups.iter().map(|rows| rows.len() as i64).collect();
    assigned = assigned.take(&groups.first_rows());
    assigned.push_column("full_count", Column::from_ints(sizes))?;
    let mut clusters = assigned;
    info!("Clusters after grouping: {} rows", clusters.nrows());
    if settings.verbose {
        write_table(
            &layout.run_dir.join("clusters.tsv"),
            &clusters,
            WriteOptions::tsv(),
        )?;
    }

    let within_distance: Vec<bool> = clusters
        .column("mismatch")?
        .iter()
        .map(|v| {
            v.as_f64()
                .is_some_and(|m| m <= settings.max_distance as f64)
        })
        .collect();
    clusters = clusters.filter(&within_distance);
    info!(
        "After filtering mismatches <= {}: {} rows",
        settings.max_distance,
        clusters.nrows()
    );
    if settings.verbose {
        write_table(
            &layout.run_dir.join("clusters-mismatch.tsv"),
            &clusters,
            WriteOptions::tsv(),
        )?;
    }

    let genes: Vec<String> = clusters
        .column("sseqid")?
        .iter()
        .map(|v| match v {
            Value::Str(s) => gene_of(s).to_string(),
            other => other.julia_string(),
        })
        .collect();
    clusters.push_column("gene", Column::from_strings(genes))?;
    let full_count: Vec<i64> = clusters
        .ints("full_count")?
        .iter()
        .map(|c| c.unwrap_or(0))
        .collect();
    let full_ratio: Vec<f64> = clusters
        .group_by(&["well", "case", "gene"])?
        .map_rows(|rows| {
            let max = rows.iter().map(|&r| full_count[r]).max().unwrap_or(1) as f64;
            rows.iter().map(|&r| full_count[r] as f64 / max).collect()
        });
    clusters.push_column("full_ratio", Column::from_floats(full_ratio))?;
    clusters.sorted_asc(&["well", "case", "sseqid"])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn edges_count_the_bases_around_the_segment() {
        assert_eq!(edges("GATTACA", "GATTACACC"), Some((0, 2)));
        assert_eq!(edges("GATTACA", "TTGATTACA"), Some((2, 0)));
        assert_eq!(edges("GATTACA", "TTTT"), None);
    }

    #[test]
    fn minus_strand_segments_are_turned_to_allele_orientation() {
        // The read carries the reverse complement of the allele GATTACA.
        let read = "CCCTGTAATCA";
        assert_eq!(
            oriented_segment("TGTAATC", read, true),
            Some(("GATTACA".to_string(), 1, 3))
        );
        assert_eq!(
            oriented_segment("TGTAATC", read, false),
            Some(("TGTAATC".to_string(), 3, 1))
        );
    }
}
