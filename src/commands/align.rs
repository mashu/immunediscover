//! `align`: the best matching database sequence of every read.
//!
//! Reads come from a FASTA file or from the table written by
//! `preprocess demultiplex`; the database is a FASTA file of alleles.
//! Every read is compared with the database in both orientations
//! ([`crate::align::search`]) and its best match is written as one row.

use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

use crate::align::gotoh::Scoring;
use crate::align::search::{Database, Hit, Hits, Settings};
use crate::error::{Error, Result};
use crate::info;
use crate::io::compress;
use crate::io::fasta::read_fasta;
use crate::io::tsv::read_tsv;
use crate::logging::with_thousands;

/// Layout of the output.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    /// Named columns, one row per matched read.
    Table,
    /// The twelve columns of `blastn -outfmt 6`, without a header.
    Blast,
}

/// Parameters of the `align` command.
#[derive(Debug, Clone)]
pub struct AlignParams {
    /// Reads: FASTA, or a table with the columns `name` and
    /// `genomic_sequence`.
    pub reads: PathBuf,
    /// Database sequences (FASTA).
    pub database: PathBuf,
    pub output: PathBuf,
    pub format: OutputFormat,
    pub settings: Settings,
}

/// Columns of the table format, after `well` and `case` when the reads
/// come from a demultiplexed table.
pub const COLUMNS: [&str; 18] = [
    "read",
    "allele",
    "strand",
    "identity",
    "differences",
    "mismatches",
    "gap_bases",
    "gaps",
    "score",
    "read_start",
    "read_end",
    "read_length",
    "allele_start",
    "allele_end",
    "allele_length",
    "coverage",
    "cigar",
    "sequence",
];

/// Statistical parameters of gapped `blastn` alignments with the default
/// scores (match 2, mismatch -3, gap 5 + 2 per base), as printed by BLAST.
const LAMBDA: f64 = 0.625;
const K: f64 = 0.41;

/// BLAST's bit score of a raw score.
pub fn bit_score(score: i32) -> f64 {
    (LAMBDA * f64::from(score) - K.ln()) / std::f64::consts::LN_2
}

/// Expected number of chance alignments with at least this score between
/// a read and a database of the given sizes. BLAST shortens both lengths
/// a little before multiplying; this value is therefore slightly larger.
pub fn e_value(score: i32, read_length: usize, database_length: usize) -> f64 {
    K * read_length as f64 * database_length as f64 * (-LAMBDA * f64::from(score)).exp()
}

/// The reads of an input file.
#[derive(Debug, Default)]
pub struct Reads {
    pub names: Vec<String>,
    pub sequences: Vec<String>,
    /// Donor of every read (`well`, `case`), for a demultiplexed table.
    pub donors: Option<Vec<(String, String)>>,
}

/// Read a FASTA file or a demultiplexed table; which one is decided by
/// the first character of the content.
pub fn read_reads(path: &Path) -> Result<Reads> {
    let mut reader = compress::open_maybe_gzip(path)?;
    let is_fasta = reader
        .fill_buf()
        .map_err(|e| Error::io(path, e))?
        .iter()
        .find(|b| !b.is_ascii_whitespace())
        .is_some_and(|&b| b == b'>');
    drop(reader);
    if is_fasta {
        let records = read_fasta(path)?;
        return Ok(Reads {
            names: records.iter().map(|r| r.identifier().to_string()).collect(),
            sequences: records.into_iter().map(|r| r.sequence).collect(),
            donors: None,
        });
    }
    let table = read_tsv(path)?;
    for column in ["name", "genomic_sequence"] {
        if !table.has_column(column) {
            return Err(Error::invalid(format!(
                "{} is neither a FASTA file nor a table of reads: column '{column}' is missing",
                path.display()
            )));
        }
    }
    let text = |column: &str| -> Result<Vec<String>> {
        Ok(table
            .column(column)?
            .iter()
            .map(|v| v.julia_string())
            .collect())
    };
    let donors = if table.has_column("well") && table.has_column("case") {
        Some(text("well")?.into_iter().zip(text("case")?).collect())
    } else {
        None
    };
    Ok(Reads {
        names: text("name")?,
        sequences: text("genomic_sequence")?,
        donors,
    })
}

/// Load a FASTA file as a search database; sequences are named by the
/// first word of their header.
pub fn load_database(path: &Path) -> Result<Database> {
    let records = read_fasta(path)?;
    Database::new(
        records
            .iter()
            .map(|r| (r.identifier().to_string(), r.sequence.as_str())),
    )
    .map_err(|e| Error::invalid(format!("{}: {e}", path.display())))
}

/// The values of one row of the table format, in the order of [`COLUMNS`].
fn table_row(name: &str, read: &str, hit: &Hit, database: &Database) -> Vec<String> {
    let subject = database.subject(hit.subject);
    let (start, end) = hit.read_range();
    let path = &hit.path;
    vec![
        name.to_string(),
        subject.name.clone(),
        if hit.minus_strand { "-" } else { "+" }.to_string(),
        format!("{:.4}", hit.identity()),
        path.differences().to_string(),
        path.mismatches().to_string(),
        path.gap_bases().to_string(),
        path.gap_openings().to_string(),
        path.score.to_string(),
        (start + 1).to_string(),
        end.to_string(),
        hit.read_length.to_string(),
        (path.subject_start + 1).to_string(),
        path.subject_end.to_string(),
        subject.sequence.len().to_string(),
        format!(
            "{:.4}",
            hit.subject_span() as f64 / subject.sequence.len() as f64
        ),
        path.cigar(),
        hit.segment(read),
    ]
}

/// The values of one row of the BLAST format.
fn blast_row(name: &str, hit: &Hit, database: &Database, database_length: usize) -> Vec<String> {
    let subject = database.subject(hit.subject);
    let (start, end) = hit.read_range();
    let path = &hit.path;
    // On the other strand BLAST counts the subject positions downwards.
    let (subject_start, subject_end) = if hit.minus_strand {
        (path.subject_end, path.subject_start + 1)
    } else {
        (path.subject_start + 1, path.subject_end)
    };
    vec![
        name.to_string(),
        subject.name.clone(),
        format!("{:.3}", 100.0 * hit.identity()),
        path.columns().to_string(),
        path.mismatches().to_string(),
        path.gap_openings().to_string(),
        (start + 1).to_string(),
        end.to_string(),
        subject_start.to_string(),
        subject_end.to_string(),
        format!(
            "{:.2e}",
            e_value(path.score, hit.read_length, database_length)
        ),
        format!("{:.1}", bit_score(path.score)),
    ]
}

/// Write the matches of `reads`.
fn write_hits(
    path: &Path,
    reads: &Reads,
    hits: &Hits,
    database: &Database,
    format: OutputFormat,
) -> Result<()> {
    let gzip = path.to_string_lossy().ends_with(".gz");
    let mut out = compress::create(path, gzip)?;
    let failed = |e| Error::io(path, e);
    if format == OutputFormat::Table {
        let donor_columns = reads.donors.as_ref().map(|_| "well\tcase\t");
        writeln!(
            out,
            "{}{}",
            donor_columns.unwrap_or_default(),
            COLUMNS.join("\t")
        )
        .map_err(failed)?;
    }
    let database_length: usize = database.subjects().iter().map(|s| s.sequence.len()).sum();
    for (number, name) in reads.names.iter().enumerate() {
        let Some(hit) = hits.get(number) else {
            continue;
        };
        let row = match format {
            OutputFormat::Table => {
                if let Some(donors) = &reads.donors {
                    let (well, case) = &donors[number];
                    write!(out, "{well}\t{case}\t").map_err(failed)?;
                }
                table_row(name, &reads.sequences[number], hit, database)
            }
            OutputFormat::Blast => blast_row(name, hit, database, database_length),
        };
        writeln!(out, "{}", row.join("\t")).map_err(failed)?;
    }
    out.flush().map_err(failed)
}

/// Run the command.
pub fn run(params: &AlignParams) -> Result<()> {
    let settings = &params.settings;
    if !settings.scoring.is_valid() {
        return Err(Error::invalid(
            "Scores must be: match above 0, mismatch above 0, gap opening 0 or more, gap extension above 0",
        ));
    }
    if params.format == OutputFormat::Blast && settings.scoring != Scoring::default() {
        return Err(Error::invalid(
            "The BLAST format holds bit scores and E-values, which are only defined for the default scores. Use the table format with other scores.",
        ));
    }
    let database = load_database(&params.database)?;
    let reads = read_reads(&params.reads)?;
    info!(
        "Aligning {} reads from {} to {} sequences from {} ({} groups of nearly equal sequences)",
        with_thousands(reads.sequences.len()),
        params.reads.display(),
        with_thousands(database.subjects().len()),
        params.database.display(),
        with_thousands(database.groups())
    );
    let sequences: Vec<&str> = reads.sequences.iter().map(String::as_str).collect();
    let hits = database.search(&sequences, settings);
    info!(
        "{} of {} reads matched an allele ({} distinct read sequences)",
        with_thousands(hits.matched()),
        with_thousands(hits.len()),
        with_thousands(hits.distinct())
    );
    write_hits(&params.output, &reads, &hits, &database, params.format)?;
    info!("Matches written to {}", params.output.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::dna::reverse_complement;

    const GENE: &str = "ACGGTCATTGCAGGCTTAACGTTAGCCATGCAAGTCCGATAGGCTTAACCGGATTACAGGT";
    const OTHER: &str = "TTGACCGTAAGGCTAGCTAGGATCCGATCGGATTCGAGCTAGGCTAACGATCGATTTGCA";

    fn write(path: &Path, text: &str) {
        std::fs::write(path, text).unwrap();
    }

    fn rows(path: &Path) -> Vec<Vec<String>> {
        compress::read_to_string(path)
            .unwrap()
            .lines()
            .map(|line| line.split('\t').map(str::to_string).collect())
            .collect()
    }

    fn params(dir: &Path, reads: &str, format: OutputFormat) -> AlignParams {
        AlignParams {
            reads: dir.join(reads),
            database: dir.join("alleles.fasta"),
            output: dir.join("hits.tsv"),
            format,
            settings: Settings::default(),
        }
    }

    #[test]
    fn fasta_reads_give_one_row_per_matched_read() {
        let dir = tempfile::tempdir().unwrap();
        write(
            &dir.path().join("alleles.fasta"),
            &format!(">IGHV1-1*01 first\n{GENE}\n>IGHV2-1*01\n{OTHER}\n"),
        );
        // One read of the second gene with a mismatch, turned around; one
        // read of nothing.
        let mut changed = OTHER.to_string();
        changed.replace_range(20..21, "T");
        let read = reverse_complement(&format!("GGGGGGG{changed}CCCC")).unwrap();
        write(
            &dir.path().join("reads.fasta"),
            &format!(
                ">r1 comment\n{read}\n>r2\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
            ),
        );
        run(&params(dir.path(), "reads.fasta", OutputFormat::Table)).unwrap();
        let rows = rows(&dir.path().join("hits.tsv"));
        assert_eq!(rows[0], COLUMNS);
        assert_eq!(rows.len(), 2);
        let value = |column: &str| {
            let at = COLUMNS.iter().position(|c| *c == column).unwrap();
            rows[1][at].as_str()
        };
        assert_eq!(value("read"), "r1");
        assert_eq!(value("allele"), "IGHV2-1*01");
        assert_eq!(value("strand"), "-");
        assert_eq!(value("differences"), "1");
        assert_eq!(value("identity"), "0.9833");
        assert_eq!((value("read_start"), value("read_end")), ("5", "64"));
        assert_eq!((value("allele_start"), value("allele_end")), ("1", "60"));
        assert_eq!(value("coverage"), "1.0000");
        assert_eq!(value("cigar"), "20=1X39=");
        assert_eq!(value("sequence"), changed);
    }

    #[test]
    fn a_demultiplexed_table_keeps_well_and_case() {
        let dir = tempfile::tempdir().unwrap();
        write(
            &dir.path().join("alleles.fasta"),
            &format!(">IGHV1-1*01\n{GENE}\n"),
        );
        write(
            &dir.path().join("reads.tsv"),
            &format!("well\tcase\tname\tgenomic_sequence\n3\tD7\tr1\tTTTT{GENE}AAAA\n"),
        );
        run(&params(dir.path(), "reads.tsv", OutputFormat::Table)).unwrap();
        let rows = rows(&dir.path().join("hits.tsv"));
        assert_eq!(&rows[0][..3], ["well", "case", "read"]);
        assert_eq!(&rows[1][..5], ["3", "D7", "r1", "IGHV1-1*01", "+"]);
    }

    #[test]
    fn the_blast_format_has_twelve_columns_and_counts_down_on_the_other_strand() {
        let dir = tempfile::tempdir().unwrap();
        write(
            &dir.path().join("alleles.fasta"),
            &format!(">IGHV1-1*01\n{GENE}\n"),
        );
        let read = reverse_complement(&format!("TTTT{GENE}AAAA")).unwrap();
        write(&dir.path().join("reads.fasta"), &format!(">r1\n{read}\n"));
        run(&params(dir.path(), "reads.fasta", OutputFormat::Blast)).unwrap();
        let rows = rows(&dir.path().join("hits.tsv"));
        assert_eq!(rows.len(), 1);
        assert_eq!(
            &rows[0][..10],
            [
                "r1",
                "IGHV1-1*01",
                "100.000",
                "61",
                "0",
                "0",
                "5",
                "65",
                "61",
                "1"
            ]
        );
        // Raw score 122 is 111.3 bits with BLAST's parameters.
        assert_eq!(rows[0][11], "111.3");
    }

    #[test]
    fn other_scores_cannot_be_written_as_blast() {
        let dir = tempfile::tempdir().unwrap();
        let mut params = params(dir.path(), "reads.fasta", OutputFormat::Blast);
        params.settings.scoring.matched = 1;
        assert!(run(&params).is_err());
    }

    #[test]
    fn a_table_without_reads_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        write(&dir.path().join("reads.tsv"), "a\tb\n1\t2\n");
        assert!(read_reads(&dir.path().join("reads.tsv")).is_err());
    }
}
