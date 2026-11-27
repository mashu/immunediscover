//! `analyze` command group.

use std::path::PathBuf;

use clap::{Args, Subcommand};

use super::validators::{at_least_one, positive_size, unit_interval};
use crate::commands::analyze::cooccurrence::{self, ClusterMethod, CooccurrenceParams};
use crate::commands::analyze::haplotype::{self, HaplotypeParams};
use crate::error::Result;

#[derive(Debug, Subcommand)]
pub enum AnalyzeCommand {
    /// Analyze cross-donor co-occurrence among alleles across cases
    Cooccurrence(CooccurrenceArgs),
    /// Infer approximate haplotypes per case using diploid assumptions
    Haplotype(HaplotypeArgs),
}

impl AnalyzeCommand {
    pub fn run(self) -> Result<()> {
        match self {
            AnalyzeCommand::Cooccurrence(args) => args.run(),
            AnalyzeCommand::Haplotype(args) => args.run(),
        }
    }
}

#[derive(Debug, Args)]
pub struct CooccurrenceArgs {
    /// TSV/TSV.GZ with columns: case and db_name (or specify columns)
    input: PathBuf,
    /// Name of column with donor/case id
    #[arg(short = 'C', long = "case-col", default_value = "case")]
    case_col: String,
    /// Name of column with allele name
    #[arg(short = 'A', long = "allele-col", default_value = "db_name")]
    allele_col: String,
    /// Alleles found in fewer donors are left out
    #[arg(short = 'm', long = "min-donors", default_value_t = 2, value_parser = positive_size)]
    min_donors: usize,
    /// How alleles are grouped: components, complete, average or single
    #[arg(long = "cluster-method", default_value = "components")]
    cluster_method: ClusterMethod,
    /// How strongly alleles must go together (0-1) to be grouped
    #[arg(long = "cluster-threshold", default_value_t = 0.5, value_parser = unit_interval)]
    cluster_threshold: f64,
    /// Print rho-based triangle diagnostics at cluster-threshold
    #[arg(long = "debug-triangles")]
    debug_triangles: bool,
    /// Smallest group of alleles reported
    #[arg(long = "min-cluster-size", default_value_t = 3, value_parser = positive_size)]
    min_cluster_size: usize,
    /// Optional path to save clusters (TSV)
    #[arg(long = "clusters")]
    clusters: Option<PathBuf>,
}

impl CooccurrenceArgs {
    fn run(self) -> Result<()> {
        cooccurrence::run(&CooccurrenceParams {
            input: self.input,
            case_column: self.case_col,
            allele_column: self.allele_col,
            min_donors: self.min_donors,
            cluster_method: self.cluster_method,
            cluster_threshold: self.cluster_threshold,
            debug_triangles: self.debug_triangles,
            min_cluster_size: self.min_cluster_size,
            clusters: self.clusters,
        })
    }
}

#[derive(Debug, Args)]
pub struct HaplotypeArgs {
    /// TSV/TSV.GZ with columns: case and allele (or specify columns)
    input: PathBuf,
    /// TSV file to save haplotype inference results
    output: PathBuf,
    /// Name of column with donor/case id
    #[arg(short = 'C', long = "case-col", default_value = "case")]
    case_col: String,
    /// Name of column with allele name
    #[arg(short = 'A', long = "allele-col", default_value = "db_name")]
    allele_col: String,
    /// Name of column with gene name (for grouping alleles by gene)
    #[arg(short = 'G', long = "gene-col", default_value = "gene")]
    gene_col: String,
    /// Minimum count for an allele to be considered
    #[arg(short = 'c', long = "mincount", default_value_t = 5, value_parser = at_least_one)]
    mincount: i64,
    /// Minimum reads of the second allele relative to the first for a donor
    /// to be called heterozygous (0-1)
    #[arg(short = 'r', long = "min-ratio", default_value_t = 0.1, value_parser = unit_interval)]
    min_ratio: f64,
    /// Optional FASTA file with novel alleles to mark in results
    #[arg(short = 'f', long = "novel-fasta")]
    novel_fasta: Option<PathBuf>,
}

impl HaplotypeArgs {
    fn run(self) -> Result<()> {
        haplotype::run(&HaplotypeParams {
            input: self.input,
            output: self.output,
            case_column: self.case_col,
            allele_column: self.allele_col,
            gene_column: self.gene_col,
            mincount: self.mincount,
            min_ratio: self.min_ratio,
            novel_fasta: self.novel_fasta,
        })
    }
}
