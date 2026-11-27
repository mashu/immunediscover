//! `align` command.

use std::path::PathBuf;

use clap::{Args, ValueEnum};

use super::validators::unit_interval;
use crate::align::gotoh::{Ends, Scoring};
use crate::align::search::Settings;
use crate::commands::align::{self, AlignParams, OutputFormat};
use crate::error::Result;

/// Layout of the output file.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum Format {
    /// One row per read with named columns: allele, strand, identity,
    /// differences, positions, CIGAR and the matching part of the read
    Table,
    /// The twelve columns of `blastn -outfmt 6`, without a header
    Blast,
}

#[derive(Debug, Args)]
pub struct AlignArgs {
    /// Reads: a FASTA file, or the table written by `preprocess
    /// demultiplex` (columns name and genomic_sequence). May be gzipped
    reads: PathBuf,
    /// Sequences to look for in the reads, for example known alleles
    /// (FASTA)
    database: PathBuf,
    /// Where to write the matches (tab-separated; gzipped if the name ends
    /// in .gz)
    output: PathBuf,
    /// Layout of the output
    #[arg(short = 'f', long = "format", value_enum, default_value_t = Format::Table)]
    format: Format,
    /// Smallest share of equal bases between a read and its allele (0-1).
    /// Reads that match no sequence this well are left out
    #[arg(short = 'i', long = "min-identity", default_value_t = 0.8, value_parser = unit_interval)]
    min_identity: f64,
    /// Smallest share of the allele that must lie inside the read (0-1).
    /// 1 keeps only reads that hold the allele from end to end
    #[arg(short = 'c', long = "min-coverage", default_value_t = 0.5, value_parser = unit_interval)]
    min_coverage: f64,
    /// Score of two equal bases
    #[arg(long = "match", default_value_t = 2, value_parser = clap::value_parser!(i32).range(1..))]
    matched: i32,
    /// Cost of two different bases
    #[arg(long = "mismatch", default_value_t = 3, value_parser = clap::value_parser!(i32).range(1..))]
    mismatch: i32,
    /// Cost of starting a gap (an inserted or deleted stretch)
    #[arg(long = "gap-open", default_value_t = 5, value_parser = clap::value_parser!(i32).range(0..))]
    gap_open: i32,
    /// Cost of every base of a gap
    #[arg(long = "gap-extend", default_value_t = 2, value_parser = clap::value_parser!(i32).range(1..))]
    gap_extend: i32,
}

impl AlignArgs {
    pub fn run(self) -> Result<()> {
        align::run(&AlignParams {
            reads: self.reads,
            database: self.database,
            output: self.output,
            format: match self.format {
                Format::Table => OutputFormat::Table,
                Format::Blast => OutputFormat::Blast,
            },
            settings: Settings {
                scoring: Scoring {
                    matched: self.matched,
                    mismatched: -self.mismatch,
                    gap_open: self.gap_open,
                    gap_extend: self.gap_extend,
                },
                min_identity: self.min_identity,
                min_coverage: self.min_coverage,
                ends: Ends::ReadMayEndInside,
            },
        })
    }
}
