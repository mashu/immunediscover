//! `fasta` command group.

use std::path::PathBuf;

use clap::{Args, Subcommand};

use crate::commands::fasta::diff;
use crate::commands::fasta::merge::{self, MergeParams};
use crate::error::Result;

#[derive(Debug, Subcommand)]
pub enum FastaCommand {
    /// Merge multiple FASTA files, keeping only unique sequences
    Merge(MergeArgs),
    /// Diff two FASTA files based on sequence identity but keep associated names
    Diff(DiffArgs),
    /// Add hash based _S suffix to all allele names in the FASTA file
    Hash(HashArgs),
}

impl FastaCommand {
    pub fn run(self) -> Result<()> {
        match self {
            FastaCommand::Merge(args) => merge::run(&MergeParams {
                output: args.output,
                inputs: args.inputs,
                sort_by_name: !args.no_sort,
                cleanup: args.cleanup,
                prefer_first: !args.prefer_last,
                add_source_prefix: args.add_source_prefix,
            })
            .map(|_| ()),
            FastaCommand::Diff(args) => diff::diff(&args.fasta),
            FastaCommand::Hash(args) => diff::hash(&args.fastain),
        }
    }
}

#[derive(Debug, Args)]
pub struct MergeArgs {
    /// Output merged FASTA file
    output: PathBuf,
    /// Input FASTA files to merge (2 or more files)
    #[arg(required = true, num_args = 1..)]
    inputs: Vec<PathBuf>,
    /// Optional regex pattern to remove from sequence names (e.g., ' Novel')
    #[arg(short = 'c', long = "cleanup", allow_hyphen_values = true)]
    cleanup: Option<String>,
    /// Disable sorting sequences by name (default: sort enabled)
    #[arg(long = "no-sort")]
    no_sort: bool,
    /// When duplicate sequences have different names, prefer the last encountered (default: prefer first)
    #[arg(long = "prefer-last")]
    prefer_last: bool,
    /// Add source filename as prefix to sequence names
    #[arg(long = "add-source-prefix")]
    add_source_prefix: bool,
}

#[derive(Debug, Args)]
pub struct DiffArgs {
    /// FASTA files with sequences
    #[arg(required = true, num_args = 1..)]
    fasta: Vec<PathBuf>,
}

#[derive(Debug, Args)]
pub struct HashArgs {
    /// Input FASTA file path
    fastain: PathBuf,
}
