//! `search` command group.

use std::path::PathBuf;

use clap::{Args, Subcommand};

use super::validators::{
    at_least_one, non_negative, positive_size, positive_unit_interval, unit_interval,
};
use crate::commands::search::bwa::{self, BwaParams};
use crate::commands::search::exact::{
    self, ExactParams, ExactSearchSettings, FlankMode, RssElement,
};
use crate::commands::search::heptamer::{self, HeptamerParams, DEFAULT_HEPTAMER_FILE};
use crate::error::Result;
use crate::seq::GeneType;

#[derive(Debug, Subcommand)]
pub enum SearchCommand {
    /// Exact match search of reads to database alleles with robust filters
    Exact(Box<ExactArgs>),
    /// Identify heptamer RSS positions and extend/trim V reads accordingly
    Heptamer(HeptamerArgs),
    /// Genome mapping QC with BWA to retain sequences mapped to the target chromosome
    Bwa(BwaArgs),
}

impl SearchCommand {
    pub fn run(self) -> Result<()> {
        match self {
            SearchCommand::Exact(args) => args.run(),
            SearchCommand::Heptamer(args) => args.run(),
            SearchCommand::Bwa(args) => args.run(),
        }
    }
}

#[derive(Debug, Args)]
pub struct ExactArgs {
    /// TSV file with demultiplexed data
    tsv: PathBuf,
    /// FASTA file with query alleles
    fasta: PathBuf,
    /// TSV file to save ouput
    output: PathBuf,
    /// Minimum allelic ratio applied within each gene group
    #[arg(short = 'f', long = "minratio", default_value_t = 0.1, value_parser = unit_interval)]
    minratio: f64,
    /// Minimum cluster size
    #[arg(short = 'c', long = "mincount", default_value_t = 5, value_parser = at_least_one)]
    mincount: i64,
    /// Minimum allelic ratio applied within each gene group for the allele against median
    #[arg(long = "min-allele-mratio", default_value_t = 0.05, value_parser = non_negative)]
    min_allele_mratio: f64,
    /// Minimum allelic ratio applied within each gene group for the gene against median
    #[arg(long = "min-gene-mratio", default_value_t = 0.05, value_parser = non_negative)]
    min_gene_mratio: f64,
    /// Comma-separated list of rss fragments: heptamer, spacer, nonamer
    #[arg(long = "rss", default_value = "heptamer")]
    rss: String,
    /// Length of extension on RSS side instead of RSS elements
    #[arg(long = "extension")]
    extension: Option<usize>,
    /// Number of nucleotides from both read ends forming a border; reject if extension overlaps
    #[arg(long = "border", default_value_t = 0)]
    border: usize,
    /// Auto-reduce extension per gene and side to avoid crossing the border
    #[arg(long = "adjust-per-gene-extension")]
    adjust_per_gene_extension: bool,
    /// Target fraction (0-1] of reads per gene that should be safe (no border overlap) when calibrating per-gene extension
    #[arg(long = "adjust-percent", default_value_t = 1.0, value_parser = positive_unit_interval)]
    adjust_percent: f64,
    /// Unfiltered exact search results for diagnostics
    #[arg(long = "raw")]
    raw: Option<PathBuf>,
    /// Disable unicode gene plot
    #[arg(short = 'n', long = "noplot")]
    noplot: bool,
    /// gene; must be one of V, D or J
    #[arg(short = 'g', long = "gene", default_value = "V")]
    gene: GeneType,
    /// Number of bases to extract from the non-RSS side of the sequence
    #[arg(short = 'a', long = "affix", default_value_t = 13, value_parser = positive_size)]
    affix: usize,
    /// Saves at most N records of flank and sequence.
    #[arg(short = 't', long = "top", default_value_t = 1, value_parser = positive_size)]
    top: usize,
    /// Space separated reference genes to use for computing ratio
    #[arg(short = 'r', long = "refgene", num_args = 0..)]
    refgene: Vec<String>,
    /// Limit to this number of sequences, zero means no limit
    #[arg(short = 'l', long = "limit", default_value_t = 0)]
    limit: usize,
    /// TSV file containing gene names and their corresponding allele_freq threshold, with two columns: name and ratio
    #[arg(short = 'e', long = "expect")]
    expect: Option<PathBuf>,
    /// TSV file containing gene names and their corresponding gene_case_freq threshold, with two columns: name and ratio
    #[arg(short = 'd', long = "deletion")]
    deletion: Option<PathBuf>,
    /// Locus to filter genes to start with this string (e.g. IGHV) excluding other genes from the analysis (i.e control genes)
    #[arg(long = "locus", default_value = "IG")]
    locus: String,
    /// Optional reference FASTA file to check if sequences are in database (adds isin_db column)
    #[arg(long = "ref-fasta")]
    ref_fasta: Option<PathBuf>,
}

impl ExactArgs {
    fn run(self) -> Result<()> {
        let mode = match self.extension {
            Some(length) => FlankMode::Extension(length),
            None => FlankMode::Rss(RssElement::parse_list(&self.rss)?),
        };
        let mut search = ExactSearchSettings::new(self.gene);
        search.mincount = self.mincount;
        search.minratio = self.minratio;
        search.affix = self.affix;
        search.mode = mode;
        search.top = self.top;
        search.border = self.border;
        search.adjust_per_gene_extension = self.adjust_per_gene_extension;
        search.adjust_percent = self.adjust_percent;
        search.raw = self.raw;
        exact::run(&ExactParams {
            tsv: self.tsv,
            fasta: self.fasta,
            output: self.output,
            search,
            min_allele_median_ratio: self.min_allele_mratio,
            min_gene_median_ratio: self.min_gene_mratio,
            noplot: self.noplot,
            refgenes: self.refgene,
            limit: self.limit,
            expect: self.expect,
            deletion: self.deletion,
            locus: self.locus,
            ref_fasta: self.ref_fasta,
        })
    }
}

#[derive(Debug, Args)]
pub struct HeptamerArgs {
    /// TSV file with demultiplexed reads
    tsv: PathBuf,
    /// FASTA file with query alleles
    fasta: PathBuf,
    /// TSV file to save data with identified heptamers
    output: PathBuf,
    /// TSV file to save summary collapsed alleles with statistics
    summary: PathBuf,
    /// JSON file with dictionary containing haptamers
    #[arg(short = 'j', long = "json", default_value = DEFAULT_HEPTAMER_FILE)]
    json: PathBuf,
    /// chain; must be one of IGKV, IGLV or IGHV
    #[arg(short = 'c', long = "chain", default_value = "IGHV", value_parser = ["IGKV", "IGLV", "IGHV"])]
    chain: String,
    /// A positive integer indicating maximum Hamming distance from any of heptamers in JSON file
    #[arg(short = 'd', long = "maxdist", default_value_t = 1)]
    maxdist: usize,
    /// How much to trim from the 5' beginning of the query sequence
    #[arg(short = 'b', long = "begin", default_value_t = 0)]
    begin: usize,
    /// How much to trim from the 3' end of the query sequence
    #[arg(short = 'e', long = "end", default_value_t = 8)]
    end: usize,
    /// Minimum count allowed in summary
    #[arg(short = 'm', long = "mincount", default_value_t = 1, value_parser = at_least_one)]
    mincount: i64,
    /// Lowest allowed ratio between counts of full allele sequence and trimmed allele sequence
    #[arg(short = 'r', long = "ratio", default_value_t = 0.25, value_parser = unit_interval)]
    ratio: f64,
}

impl HeptamerArgs {
    fn run(self) -> Result<()> {
        heptamer::run(&HeptamerParams {
            tsv: self.tsv,
            fasta: self.fasta,
            output: self.output,
            summary: self.summary,
            json: self.json,
            chain: self.chain,
            max_distance: self.maxdist,
            trim_begin: self.begin,
            trim_end: self.end,
            mincount: self.mincount,
            ratio: self.ratio,
        })
    }
}

#[derive(Debug, Args)]
pub struct BwaArgs {
    /// TSV file with candidate alleles (e.g. `discover blast` output)
    tsv: PathBuf,
    /// TSV file to save filtered input (gzip-compressed, `.gz` is appended)
    output: PathBuf,
    /// FASTA file(s) with a BWA index; alignments to all of them are pooled
    #[arg(required = true, num_args = 1..)]
    genome: Vec<PathBuf>,
    /// Chromosome string to filter by
    #[arg(short = 'c', long = "chromosome", default_value = "chromosome 14")]
    chromosome: String,
    /// Name of the column with allele names
    #[arg(short = 'n', long = "colname", default_value = "allele_name")]
    colname: String,
    /// Columns concatenated into the sequence to align
    #[arg(short = 's', long = "colseq", num_args = 0.., default_values = ["aln_qseq"])]
    colseq: Vec<String>,
    /// Regex to filter valid descriptions of chromosomes
    #[arg(
        short = 't',
        long = "tag",
        default_value = "(.*Primary Assembly.*)|(.*alternate locus.*)"
    )]
    tag: String,
}

impl BwaArgs {
    fn run(self) -> Result<()> {
        bwa::run(&BwaParams {
            tsv: self.tsv,
            output: self.output,
            genomes: self.genome,
            chromosome: self.chromosome,
            name_column: self.colname,
            sequence_columns: self.colseq,
            tag: self.tag,
        })
    }
}
