//! `preprocess` command group.

use std::path::PathBuf;

use clap::{Args, Subcommand};

use crate::commands::preprocess::demultiplex::{self, DemultiplexParams};
use crate::error::Result;

#[derive(Debug, Subcommand)]
pub enum PreprocessCommand {
    /// Demultiplex indexed plate libraries into a TSV with per-read metadata
    Demultiplex(DemultiplexArgs),
}

impl PreprocessCommand {
    pub fn run(self) -> Result<()> {
        match self {
            PreprocessCommand::Demultiplex(args) => args.run(),
        }
    }
}

#[derive(Debug, Args)]
pub struct DemultiplexArgs {
    /// Input FASTQ file with reads (single-end)
    fastq: PathBuf,
    /// Index sheet (TSV) with the columns forward_index, reverse_index and
    /// case; one row per well
    indices: PathBuf,
    /// Output TSV (gz auto-enabled) with demultiplexed reads and metadata
    output: PathBuf,
    /// Reads shorter than this are dropped
    #[arg(short = 'l', long = "length", default_value_t = 200)]
    length: usize,
    /// Also write one FASTQ file per donor
    #[arg(short = 's', long = "split")]
    split: bool,
    /// Name of the forward array index to use for demultiplexing (if present)
    #[arg(short = 'f', long = "forwardarrayindex", default_value = "")]
    forwardarrayindex: String,
    /// Regex to keep only cases matching pattern (e.g., '[ACDERF]')
    #[arg(long = "case-filter-regex")]
    case_filter_regex: Option<String>,
}

impl DemultiplexArgs {
    fn run(self) -> Result<()> {
        demultiplex::run(&DemultiplexParams {
            fastq: self.fastq,
            indices: self.indices,
            output: self.output,
            min_length: self.length,
            split: self.split,
            forward_array_index: self.forwardarrayindex,
            case_filter_regex: self.case_filter_regex,
        })
    }
}
