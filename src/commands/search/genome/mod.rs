//! `search genome`: keep candidate sequences whose nearest place in a
//! genome lies on the expected chromosome.
//!
//! Every candidate is looked for in the whole genome, in both
//! orientations ([`crate::align::genome`]). Its nearest place on the
//! wanted chromosome is compared with its nearest place anywhere else: a
//! candidate is kept when no place elsewhere is nearer. The genome is read
//! as a stream, from a file or from the internet, and is not stored.

pub mod chromosome;
pub mod source;

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use self::chromosome::{chromosome_of, wanted_chromosome};
use self::source::{Origin, Source};
use crate::align::genome::{GenomeSearch, Place, Places};
use crate::commands::preprocess::pattern::compile;
use crate::error::{Error, Result};
use crate::frame::{Column, Frame};
use crate::io::fasta::{identifier_of, FastaStream};
use crate::io::path::always_gz;
use crate::io::tsv::{read_tsv, write_tsv, write_tsv_gz};
use crate::logging::{with_thousands, Progress};
use crate::{info, warn};

#[derive(Debug, Clone)]
pub struct GenomeParams {
    pub tsv: PathBuf,
    pub output: PathBuf,
    /// Genomes: FASTA files, `-` for the standard input, names of known
    /// assemblies or addresses.
    pub genomes: Vec<String>,
    /// The chromosome the candidates should lie on (`14`, `chr14`,
    /// `chromosome 14`).
    pub chromosome: String,
    /// Column with allele names.
    pub name_column: String,
    /// Columns concatenated into the sequence to look for.
    pub sequence_columns: Vec<String>,
    /// Regex that the header of a record must also match to count as the
    /// wanted chromosome.
    pub tag: Option<String>,
    /// A place counts when it differs from the candidate in at most
    /// `1 - min_identity` of the candidate's bases.
    pub min_identity: f64,
}

/// A record of a genome that was searched.
#[derive(Debug, Clone)]
struct Record {
    /// Full header line.
    header: String,
    genome: String,
}

impl Record {
    fn name(&self) -> &str {
        identifier_of(&self.header)
    }
}

/// What the search found for one candidate.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct GenomeCheck {
    /// No place off the chromosome is nearer than the one on it.
    pub on_target: bool,
    /// Nearest place on the wanted chromosome, as `record:start` (from 0).
    pub position: Option<String>,
    pub edit_distance: Option<i64>,
    pub ref_sequence: Option<String>,
    pub orientation: Option<String>,
    pub cigar: Option<String>,
    /// Nearest place on any other record.
    pub elsewhere: Option<String>,
    pub elsewhere_distance: Option<i64>,
    /// Header and genome of the record of `elsewhere`.
    elsewhere_record: Option<(String, String)>,
}

fn position(place: &Place, records: &[Record]) -> String {
    format!("{}:{}", records[place.record].name(), place.start)
}

impl GenomeCheck {
    fn of(places: &Places, records: &[Record]) -> Self {
        let on = places.on_target.as_ref();
        let off = places.elsewhere.as_ref();
        GenomeCheck {
            // Equally near on and off the chromosome counts as on it.
            on_target: on.is_some_and(|on| off.is_none_or(|off| on.differences <= off.differences)),
            position: on.map(|place| position(place, records)),
            edit_distance: on.map(|place| place.differences as i64),
            ref_sequence: on.map(|place| place.reference.clone()),
            orientation: on.map(|place| if place.minus_strand { "-" } else { "+" }.to_string()),
            cigar: on.map(|place| place.cigar.clone()),
            elsewhere: off.map(|place| position(place, records)),
            elsewhere_distance: off.map(|place| place.differences as i64),
            elsewhere_record: off.map(|place| {
                let record = &records[place.record];
                (record.header.clone(), record.genome.clone())
            }),
        }
    }
}

/// Search every genome for `sequences` and tell for each whether it
/// belongs to `chromosome`.
///
/// The genomes are pooled: the nearest place on the chromosome and the
/// nearest place elsewhere are the nearest over all of them.
pub fn check_sequences(
    genomes: &[String],
    sequences: &[String],
    chromosome: &str,
    tag: Option<&str>,
    min_identity: f64,
) -> Result<Vec<GenomeCheck>> {
    let tag = tag.map(compile).transpose()?;
    let wanted = wanted_chromosome(chromosome);
    // Equal sequences are looked for once.
    let mut number_of: HashMap<String, usize> = HashMap::new();
    let mut distinct: Vec<String> = Vec::new();
    let of_row: Vec<usize> = sequences
        .iter()
        .map(|sequence| {
            let key = sequence.trim().to_ascii_uppercase();
            *number_of.entry(key.clone()).or_insert_with(|| {
                distinct.push(key);
                distinct.len() - 1
            })
        })
        .collect();
    let texts: Vec<&str> = distinct.iter().map(String::as_str).collect();
    let search = GenomeSearch::new(&texts, min_identity)?;

    let mut records: Vec<Record> = Vec::new();
    let mut on_chromosome = 0usize;
    let mut bases: Vec<u8> = Vec::new();
    for genome in genomes {
        let origin = Origin::of(genome)?;
        let label = origin.label();
        match &origin {
            Origin::Internet { url, .. } => {
                info!("Reading {label} from {url}. It is searched while it arrives and not stored.")
            }
            _ => info!("Reading the genome from {label}"),
        }
        let mut source = Source::open(&origin)?;
        let progress = Progress::new("Searching", "Searched", "bases", 0);
        let mut seen = 0usize;
        let read = (|| -> Result<()> {
            let mut stream = FastaStream::new(&mut source.reader, label.clone());
            while let Some(header) = stream.next_record(&mut bases)? {
                let on_target = chromosome_of(&header).is_some_and(|found| found == wanted)
                    && tag
                        .as_ref()
                        .is_none_or(|tag| tag.is_match(&header).unwrap_or(false));
                on_chromosome += usize::from(on_target);
                seen += 1;
                search.search(records.len(), on_target, &bases);
                progress.add(bases.len());
                records.push(Record {
                    header,
                    genome: label.clone(),
                });
            }
            Ok(())
        })();
        if let Err(error) = read {
            // A download that broke off shows as a damaged stream; the
            // reason curl gives is the better message.
            return Err(source.abandon().unwrap_or(error));
        }
        source.finish()?;
        progress.finish();
        if seen == 0 {
            return Err(Error::invalid(format!(
                "The genome from {label} holds no sequences"
            )));
        }
    }
    if on_chromosome == 0 {
        let examples: Vec<&str> = records.iter().take(3).map(|r| r.header.as_str()).collect();
        return Err(Error::invalid(format!(
            "No record of the genome belongs to chromosome {wanted}, so nothing can be said about the candidates. Its records start with: {}. Check --chromosome{}",
            examples.join(" | "),
            if tag.is_some() { " and --tag" } else { "" }
        )));
    }
    info!(
        "{} genome records searched, {} of them of chromosome {wanted}",
        with_thousands(records.len()),
        with_thousands(on_chromosome)
    );
    let places = search.finish();
    Ok(of_row
        .iter()
        .map(|&number| GenomeCheck::of(&places[number], &records))
        .collect())
}

/// `<output>.discarded.tsv` next to the result table.
fn discarded_path(output: &Path) -> PathBuf {
    let text = output.to_string_lossy();
    let base = text.strip_suffix(".gz").unwrap_or(&text);
    PathBuf::from(format!("{base}.discarded.tsv"))
}

fn text_column<'a>(values: impl Iterator<Item = Option<&'a String>>) -> Column {
    Column::Str(values.map(|value| value.cloned()).collect())
}

/// Run the command. The result is gzip-compressed and written to `output`
/// with a `.gz` suffix; the candidates that were dropped are listed in
/// `<output>.discarded.tsv`.
pub fn run(params: &GenomeParams) -> Result<()> {
    info!("Genome search: do the candidates lie on the expected chromosome?");
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
    if let Some(row) = sequences.iter().position(|s| s.trim().is_empty()) {
        return Err(Error::invalid(format!(
            "Candidate '{}' has no sequence in column(s) {}",
            names[row],
            params.sequence_columns.join(", ")
        )));
    }

    let checks = check_sequences(
        &params.genomes,
        &sequences,
        &params.chromosome,
        params.tag.as_deref(),
        params.min_identity,
    )?;
    let wanted = wanted_chromosome(&params.chromosome);

    let ints =
        |value: fn(&GenomeCheck) -> Option<i64>| Column::Int(checks.iter().map(value).collect());
    table.set_column(
        "position",
        text_column(checks.iter().map(|c| c.position.as_ref())),
    )?;
    table.set_column("edit_distance", ints(|c| c.edit_distance))?;
    table.set_column(
        "ref_sequence",
        text_column(checks.iter().map(|c| c.ref_sequence.as_ref())),
    )?;
    table.set_column(
        "orientation",
        text_column(checks.iter().map(|c| c.orientation.as_ref())),
    )?;
    table.set_column(
        "cigar",
        text_column(checks.iter().map(|c| c.cigar.as_ref())),
    )?;
    table.set_column(
        "elsewhere",
        text_column(checks.iter().map(|c| c.elsewhere.as_ref())),
    )?;
    table.set_column("elsewhere_distance", ints(|c| c.elsewhere_distance))?;

    let mask: Vec<bool> = checks.iter().map(|c| c.on_target).collect();
    let kept = mask.iter().filter(|kept| **kept).count();
    info!(
        "{} sequences processed, {kept} lie on chromosome {wanted}",
        table.nrows()
    );
    // Kept, but with a place elsewhere that is just as near.
    let ambiguous = checks
        .iter()
        .filter(|c| c.on_target && c.elsewhere_distance == c.edit_distance)
        .count();
    if ambiguous > 0 {
        warn!(
            "{ambiguous} of them have an equally near place outside chromosome {wanted} (columns elsewhere, elsewhere_distance)"
        );
    }
    let output = always_gz(&params.output);
    write_tsv_gz(&output, &table.filter(&mask))?;
    info!("Filtered result saved in {}", output.display());

    // The candidates that were dropped, with where they lie instead.
    let dropped: Vec<usize> = (0..checks.len()).filter(|&row| !mask[row]).collect();
    for &row in &dropped {
        match &checks[row].elsewhere_record {
            Some((header, genome)) => {
                info!("Discarded {} matching {header} in {genome}", names[row])
            }
            None => info!(
                "Discarded {}: no place in the genome is near enough",
                names[row]
            ),
        }
    }
    let no_place = "no place near enough".to_string();
    let report = Frame::from_columns([
        (
            "name",
            Column::from_strings(dropped.iter().map(|&row| names[row].as_str())),
        ),
        (
            "description",
            Column::from_strings(dropped.iter().map(|&row| {
                let record = checks[row].elsewhere_record.as_ref();
                record.map_or(no_place.as_str(), |(header, _)| header.as_str())
            })),
        ),
        (
            "genome",
            text_column(dropped.iter().map(|&row| {
                checks[row]
                    .elsewhere_record
                    .as_ref()
                    .map(|(_, genome)| genome)
            })),
        ),
        (
            "position",
            text_column(dropped.iter().map(|&row| checks[row].elsewhere.as_ref())),
        ),
        (
            "edit_distance",
            Column::Int(
                dropped
                    .iter()
                    .map(|&row| checks[row].elsewhere_distance)
                    .collect(),
            ),
        ),
        (
            "chromosome_position",
            text_column(dropped.iter().map(|&row| checks[row].position.as_ref())),
        ),
        (
            "chromosome_distance",
            Column::Int(
                dropped
                    .iter()
                    .map(|&row| checks[row].edit_distance)
                    .collect(),
            ),
        ),
    ])?;
    let discarded_file = discarded_path(&params.output);
    write_tsv(&discarded_file, &report)?;
    info!(
        "{} discarded sequences listed in {}",
        dropped.len(),
        discarded_file.display()
    );
    Ok(())
}
