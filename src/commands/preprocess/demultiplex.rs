//! `preprocess demultiplex`: assign reads to plate wells and donors by
//! their forward and reverse index sequences.

use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

use super::pattern::{compile, Anchored};
use crate::commands::report::print_counts_table;
use crate::error::{Error, Result};
use crate::frame::{round_digits, Column, Frame, SortKey};
use crate::io::fastq::{for_each_record, FastqRecord};
use crate::io::path::always_gz;
use crate::io::tsv::{read_table, write_tsv, write_tsv_gz, ReadOptions};
use crate::io::{compress, fastq};
use crate::stats::{mean, sample_std};
use crate::{info, warn};

const INDEX_COLUMNS: [&str; 3] = ["forward_index", "reverse_index", "case"];

/// Parameters of the demultiplex command.
#[derive(Debug, Clone)]
pub struct DemultiplexParams {
    pub fastq: PathBuf,
    pub indices: PathBuf,
    pub output: PathBuf,
    /// Reads shorter than this are dropped.
    pub min_length: usize,
    /// Also write one FASTQ per donor next to the input.
    pub split: bool,
    /// Optional array index expected after the first two bases.
    pub forward_array_index: String,
    /// Keep only donors whose identifier starts with a match of this regex.
    pub case_filter_regex: Option<String>,
}

/// One row of the index sheet.
struct IndexPattern {
    forward: Anchored,
    reverse: Anchored,
    case: String,
}

fn load_patterns(path: &Path, array_index: &str) -> Result<Vec<IndexPattern>> {
    let options = ReadOptions::default().text_column("case");
    let sheet = read_table(path, &options)?;
    if sheet.names() != INDEX_COLUMNS {
        return Err(Error::invalid(
            "File must contain following columns forward_index, reverse_index, case",
        ));
    }
    info!("Loaded {} indices", sheet.nrows());
    let text =
        |name: &str, row: usize| -> Result<String> { Ok(sheet.value(row, name)?.julia_string()) };
    (0..sheet.nrows())
        .map(|row| {
            let forward = text("forward_index", row)?;
            let forward = if array_index.is_empty() {
                Anchored::prefix(&forward)?
            } else {
                Anchored::prefix(&format!("^.{{2}}{array_index}.{{24}}{forward}"))?
            };
            Ok(IndexPattern {
                forward,
                reverse: Anchored::suffix(&text("reverse_index", row)?)?,
                case: text("case", row)?,
            })
        })
        .collect()
}

/// Lazily opened per-donor FASTQ writers.
struct SplitWriters {
    directory: PathBuf,
    writers: HashMap<String, Box<dyn Write + Send>>,
}

impl SplitWriters {
    fn new(fastq: &Path) -> Self {
        let mut directory = fastq.as_os_str().to_os_string();
        directory.push("_split");
        SplitWriters {
            directory: PathBuf::from(directory),
            writers: HashMap::new(),
        }
    }

    fn write(&mut self, case: &str, record: &FastqRecord) -> Result<()> {
        if !self.writers.contains_key(case) {
            std::fs::create_dir_all(&self.directory).map_err(|e| Error::io(&self.directory, e))?;
            let path = self.directory.join(format!("{case}.fastq"));
            self.writers
                .insert(case.to_string(), compress::create(&path, false)?);
        }
        let writer = self
            .writers
            .get_mut(case)
            .expect("writer was just inserted");
        record.write_to(writer).map_err(Error::from)
    }

    fn finish(self) -> Result<()> {
        for mut writer in self.writers.into_values() {
            writer.flush()?;
        }
        Ok(())
    }
}

/// Demultiplexed reads plus per-well statistics.
pub struct Demultiplexed {
    /// Columns `well`, `case`, `name`, `genomic_sequence`.
    pub reads: Frame,
    /// Per-well read count and length statistics.
    pub stats: Frame,
}

/// Assign every read of `fastq` to the first index pair it carries.
pub fn demultiplex(
    fastq: &Path,
    indices: &Path,
    array_index: &str,
    min_length: usize,
    split: bool,
) -> Result<Demultiplexed> {
    let patterns = load_patterns(indices, array_index)?;
    let mut wells: Vec<i64> = Vec::new();
    let mut cases: Vec<String> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    let mut sequences: Vec<String> = Vec::new();
    let mut writers = split.then(|| SplitWriters::new(fastq));

    for_each_record(fastq, |record| {
        let sequence = record.sequence.as_str();
        for (i, pattern) in patterns.iter().enumerate() {
            if pattern.forward.matches_start(sequence) && pattern.reverse.matches_end(sequence) {
                if sequence.len() < min_length {
                    // A short read may still be claimed by a later index pair.
                    continue;
                }
                wells.push(i as i64 + 1);
                cases.push(pattern.case.clone());
                names.push(record.identifier().to_string());
                sequences.push(sequence.to_string());
                if let Some(writers) = writers.as_mut() {
                    writers.write(&pattern.case, record)?;
                }
                break;
            }
        }
        Ok(())
    })?;
    if let Some(writers) = writers {
        writers.finish()?;
    }

    if wells.is_empty() {
        warn!("No reads were demultiplexed. Check your indices and FASTQ file.");
        let empty_ints = || Column::Int(Vec::new());
        let empty_text = || Column::Str(Vec::new());
        return Ok(Demultiplexed {
            reads: Frame::from_columns([
                ("well", empty_ints()),
                ("case", empty_text()),
                ("name", empty_text()),
                ("genomic_sequence", empty_text()),
            ])?,
            stats: Frame::from_columns([
                ("well", empty_ints()),
                ("total", empty_ints()),
                ("short", empty_ints()),
            ])?,
        });
    }

    let reads = Frame::from_columns([
        ("well", Column::from_ints(wells)),
        ("case", Column::from_strings(cases)),
        ("name", Column::from_strings(names)),
        (
            "genomic_sequence",
            Column::from_strings(sequences.iter().map(String::as_str)),
        ),
    ])?;
    let stats = well_statistics(&reads, &sequences)?;
    Ok(Demultiplexed { reads, stats })
}

fn well_statistics(reads: &Frame, sequences: &[String]) -> Result<Frame> {
    let groups = reads.group_by(&["well"])?;
    let wells = reads.ints("well")?;
    let mut well = Vec::new();
    let mut count = Vec::new();
    let mut mean_length = Vec::new();
    let mut std_length = Vec::new();
    for rows in groups.iter() {
        let lengths: Vec<f64> = rows.iter().map(|&r| sequences[r].len() as f64).collect();
        well.push(wells[rows[0]].expect("well is never missing"));
        count.push(rows.len() as i64);
        mean_length.push(round_digits(mean(&lengths), 1));
        std_length.push(round_digits(sample_std(&lengths), 1));
    }
    Frame::from_columns([
        ("well", Column::from_ints(well)),
        ("count", Column::from_ints(count)),
        ("mean_length", Column::from_floats(mean_length)),
        ("std_length", Column::from_floats(std_length)),
    ])?
    .sorted(&[SortKey::asc("well")])
}

fn count_distinct(frame: &Frame, column: &str) -> Result<usize> {
    Ok(frame.unique_rows(&[column])?.len())
}

/// Run the command: demultiplex, report, and write the read table.
pub fn run(params: &DemultiplexParams) -> Result<()> {
    info!("Demultiplexing");
    let Demultiplexed { mut reads, stats } = demultiplex(
        &params.fastq,
        &params.indices,
        &params.forward_array_index,
        params.min_length,
        params.split,
    )?;

    if let Some(pattern) = &params.case_filter_regex {
        info!("Filtering cases with regex: {pattern}");
        let before = reads.nrows();
        let regex = compile(&format!("^(?:{pattern})"))?;
        let cases = reads.strs("case")?.to_vec();
        reads.retain_rows(|row| {
            cases[row]
                .as_deref()
                .is_some_and(|case| regex.is_match(case).unwrap_or(false))
        });
        info!("Filtered from {before} to {} rows", reads.nrows());
    }

    info!("Demultiplexing summary:");
    info!("  - Total sequences: {}", reads.nrows());
    info!("  - Unique cases: {}", count_distinct(&reads, "case")?);
    info!("  - Unique wells: {}", count_distinct(&reads, "well")?);

    let mut log_path = params.output.as_os_str().to_os_string();
    log_path.push(".log");
    write_tsv(Path::new(&log_path), &stats)?;

    let groups = reads.group_by(&["case"])?;
    let cases = reads.strs("case")?;
    let mut per_case: Vec<(String, usize)> = groups
        .iter()
        .map(|rows| (cases[rows[0]].clone().unwrap_or_default(), rows.len()))
        .collect();
    per_case.sort_by_key(|(_, count)| std::cmp::Reverse(*count));
    let (labels, counts): (Vec<String>, Vec<usize>) = per_case.into_iter().unzip();
    print_counts_table(&labels, &counts);

    let output = always_gz(&params.output);
    write_tsv_gz(&output, &reads)?;
    info!(
        "Demultiplexing data saved in compressed {} file",
        output.display()
    );
    Ok(())
}

/// Convert FASTA-style records into a FASTQ file with constant quality.
/// Used to build test inputs.
pub fn write_constant_quality_fastq(path: &Path, records: &[(String, String)]) -> Result<()> {
    let records: Vec<FastqRecord> = records
        .iter()
        .map(|(name, sequence)| FastqRecord {
            description: name.clone(),
            sequence: sequence.clone(),
            quality: "I".repeat(sequence.len()),
        })
        .collect();
    fastq::write_fastq(path, &records)
}
