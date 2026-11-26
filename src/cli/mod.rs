//! Command-line interface: argument definitions and dispatch.
//!
//! The command tree, option names and defaults mirror `immunediscover.jl`.
//! Each argument struct only translates the command line into the plain
//! parameter struct of the corresponding [`crate::commands`] module.

mod analyze;
mod discover;
mod fasta;
mod invocation_log;
mod preprocess;
pub mod presets;
mod search;
mod table;
mod validators;

use clap::{Parser, Subcommand};

use crate::error::Result;

/// Version reported by `--version` and recorded in the invocation log.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Debug, Parser)]
#[command(
    name = "immunediscover",
    version,
    about = "Tool for processing immune NGS data",
    after_help = "GKHLab",
    subcommand_required = true,
    arg_required_else_help = true
)]
pub struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// De novo allele discovery (blast, hsmm)
    #[command(subcommand)]
    Discover(discover::DiscoverCommand),
    /// Search against known references (exact, heptamer, bwa)
    #[command(subcommand)]
    Search(search::SearchCommand),
    /// Downstream analysis (cooccurrence, haplotype)
    #[command(subcommand)]
    Analyze(analyze::AnalyzeCommand),
    /// Data preparation (demultiplex)
    #[command(subcommand)]
    Preprocess(preprocess::PreprocessCommand),
    /// Table utilities (outerjoin, leftjoin, transform, aggregate, unique, sort, filter, select, fasta, collect, exclude)
    #[command(subcommand)]
    Table(table::TableCommand),
    /// FASTA utilities (merge, diff, hash)
    #[command(subcommand)]
    Fasta(fasta::FastaCommand),
}

impl Cli {
    /// Execute the selected command.
    pub fn run(self) -> Result<()> {
        match self.command {
            Command::Discover(command) => command.run(),
            Command::Search(command) => command.run(),
            Command::Analyze(command) => command.run(),
            Command::Preprocess(command) => command.run(),
            Command::Table(command) => command.run(),
            Command::Fasta(command) => command.run(),
        }
    }
}

/// Parse the process arguments, record them, and run the command.
pub fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    invocation_log::record(VERSION, &args);
    Cli::parse().run()
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn command_tree_is_consistent() {
        Cli::command().debug_assert();
    }

    // The same invocations the upstream precompile workload and tests parse.
    #[test]
    fn upstream_invocations_parse() {
        let cases: &[&[&str]] = &[
            &["discover", "blast", "i.tsv", "d.fa", "o.tsv", "-g", "V"],
            &[
                "discover",
                "blast",
                "i.tsv",
                "d.fa",
                "o.tsv",
                "-a",
                "-task blastn -word_size 7",
                "-i",
            ],
            &["discover", "hsmm", "i.tsv", "d.fa", "o.tsv.gz"],
            &["search", "exact", "i.tsv.gz", "d.fa", "o.tsv.gz"],
            &[
                "search", "exact", "i.tsv.gz", "d.fa", "o.tsv.gz", "-r", "IGHV3-23", "IGHV1-2",
                "-g", "J",
            ],
            &[
                "search", "heptamer", "i.tsv.gz", "d.fa", "o.tsv.gz", "s.tsv",
            ],
            &["search", "bwa", "i.tsv", "o.tsv", "g.fa"],
            &[
                "search", "bwa", "i.tsv", "o.tsv", "g1.fa", "g2.fa", "-s", "seq",
            ],
            &["analyze", "cooccurrence", "i.tsv"],
            &["analyze", "haplotype", "i.tsv", "o.tsv"],
            &["preprocess", "demultiplex", "i.fq", "idx.tsv", "o.tsv"],
            &["table", "outerjoin", "l.tsv", "r.tsv", "o.tsv", "-k", "key"],
            &["table", "leftjoin", "l.tsv", "r.tsv", "o.tsv", "-k", "key"],
            &[
                "table",
                "transform",
                "i.tsv",
                "o.tsv",
                "-c",
                "col",
                "-p",
                "(.*)",
                "-r",
                "\\1",
            ],
            &["table", "aggregate", "i.tsv", "o.tsv", "-g", "col"],
            &["table", "unique", "i.tsv", "o.tsv", "-c", "col"],
            &["table", "sort", "i.tsv", "o.tsv", "-c", "col"],
            &[
                "table",
                "filter",
                "i.tsv",
                "o.tsv",
                "-c",
                "col",
                "--pattern",
                "x",
            ],
            &[
                "table",
                "filter",
                "i.tsv",
                "o.tsv",
                "-c",
                "col",
                "--operator",
                ">=",
                "--threshold",
                "-2.5",
            ],
            &["table", "select", "i.tsv", "o.tsv", "-c", "col"],
            &["table", "fasta", "i.tsv", "o.fa"],
            &["table", "collect", "*.tsv", "o.tsv"],
            &["table", "exclude", "i.tsv", "o.tsv", "r.fa"],
            &["fasta", "merge", "o.fa", "a.fa", "b.fa"],
            &["fasta", "diff", "a.fa", "b.fa"],
            &["fasta", "hash", "i.fa"],
        ];
        for args in cases {
            let full: Vec<&str> = std::iter::once("immunediscover")
                .chain(args.iter().copied())
                .collect();
            Cli::try_parse_from(&full).unwrap_or_else(|e| panic!("{args:?} failed to parse: {e}"));
        }
    }

    #[test]
    fn out_of_range_values_are_rejected() {
        let bad: &[&[&str]] = &[
            &["search", "exact", "i", "d", "o", "--minratio", "1.5"],
            &["search", "exact", "i", "d", "o", "--mincount", "0"],
            &["search", "exact", "i", "d", "o", "-g", "X"],
            &["search", "heptamer", "i", "d", "o", "s", "-c", "TRAV"],
            &["discover", "blast", "i", "d", "o", "-g", "Q"],
            &["analyze", "cooccurrence", "i", "--cluster-method", "ward"],
            &["table", "filter", "i", "o", "-c", "x", "--operator", "=="],
        ];
        for args in bad {
            let full: Vec<&str> = std::iter::once("immunediscover")
                .chain(args.iter().copied())
                .collect();
            assert!(
                Cli::try_parse_from(&full).is_err(),
                "{args:?} should be rejected"
            );
        }
    }
}
