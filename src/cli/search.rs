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
    /// Reads of the plate: the table written by `preprocess demultiplex`
    tsv: PathBuf,
    /// Known alleles (FASTA, names like IGHV1-69*01)
    fasta: PathBuf,
    /// Where to write the result table (.gz is added)
    output: PathBuf,
    /// Minimum reads of an allele relative to the best allele of the same
    /// gene in the donor (0-1)
    #[arg(short = 'f', long = "minratio", default_value_t = 0.1, value_parser = unit_interval)]
    minratio: f64,
    /// Minimum reads of an allele in a donor
    #[arg(short = 'c', long = "mincount", default_value_t = 5, value_parser = at_least_one)]
    mincount: i64,
    /// Minimum reads of an allele in a donor relative to what that allele
    /// usually has across donors (its median). Removes traces carried over
    /// from other wells
    #[arg(long = "min-allele-mratio", default_value_t = 0.05, value_parser = non_negative)]
    min_allele_mratio: f64,
    /// Minimum reads of a gene in a donor relative to what that gene usually
    /// has across donors (its median)
    #[arg(long = "min-gene-mratio", default_value_t = 0.05, value_parser = non_negative)]
    min_gene_mratio: f64,
    /// Parts of the recombination signal to report next to the allele,
    /// comma-separated: heptamer, spacer, nonamer
    #[arg(long = "rss", default_value = "heptamer")]
    rss: String,
    /// Report this many bases on the recombination-signal side instead of
    /// the signal parts
    #[arg(long = "extension")]
    extension: Option<usize>,
    /// Number of bases at each end of a read that the bases reported next to
    /// the allele must not reach into; reads where they do are rejected
    #[arg(long = "border", default_value_t = 0)]
    border: usize,
    /// Shorten the extension per gene automatically so that it stays clear
    /// of the border
    #[arg(long = "adjust-per-gene-extension")]
    adjust_per_gene_extension: bool,
    /// With --adjust-per-gene-extension: share (0-1] of a gene's reads that
    /// must be clear of the border
    #[arg(long = "adjust-percent", default_value_t = 1.0, value_parser = positive_unit_interval)]
    adjust_percent: f64,
    /// Also write the result before filtering to this file, for troubleshooting
    #[arg(long = "raw")]
    raw: Option<PathBuf>,
    /// Do not draw the gene plot on the screen
    #[arg(short = 'n', long = "noplot")]
    noplot: bool,
    /// Gene type: V, D or J. Decides on which side of the allele the
    /// recombination signal is expected
    #[arg(short = 'g', long = "gene", default_value = "V")]
    gene: GeneType,
    /// Number of bases to report on the side of the allele away from the
    /// recombination signal
    #[arg(short = 'a', long = "affix", default_value_t = 13, value_parser = positive_size)]
    affix: usize,
    /// How many versions of the bases next to an allele to report per allele
    /// and donor
    #[arg(short = 't', long = "top", default_value_t = 1, value_parser = positive_size)]
    top: usize,
    /// Reference gene(s), space-separated; adds columns with counts relative
    /// to them
    #[arg(short = 'r', long = "refgene", num_args = 0..)]
    refgene: Vec<String>,
    /// Use only the first N reads (0 = all). For quick tests
    #[arg(short = 'l', long = "limit", default_value_t = 0)]
    limit: usize,
    /// Table with the columns name and ratio: a minimum allele_freq for each
    /// name listed
    #[arg(short = 'e', long = "expect")]
    expect: Option<PathBuf>,
    /// Table with the columns name and ratio: a gene_case_freq threshold for
    /// each gene listed
    #[arg(short = 'd', long = "deletion")]
    deletion: Option<PathBuf>,
    /// Analyse only genes whose name starts with this text (for example
    /// IGHV); other genes, such as control genes, are left out
    #[arg(long = "locus", default_value = "IG")]
    locus: String,
    /// A second FASTA; adds the column isin_db, saying whether each sequence
    /// is in it
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
    /// Known alleles (FASTA, names like IGHV1-69*01)
    fasta: PathBuf,
    /// TSV file to save data with identified heptamers
    output: PathBuf,
    /// TSV file to save summary collapsed alleles with statistics
    summary: PathBuf,
    /// JSON file listing the accepted heptamers per chain
    #[arg(short = 'j', long = "json", default_value = DEFAULT_HEPTAMER_FILE)]
    json: PathBuf,
    /// Chain: IGHV, IGKV or IGLV. Selects the heptamer list
    #[arg(short = 'c', long = "chain", default_value = "IGHV", value_parser = ["IGKV", "IGLV", "IGHV"])]
    chain: String,
    /// Number of differences allowed between a heptamer in a read and the
    /// listed ones
    #[arg(short = 'd', long = "maxdist", default_value_t = 1)]
    maxdist: usize,
    /// Number of bases to ignore at the start of each database allele
    #[arg(short = 'b', long = "begin", default_value_t = 0)]
    begin: usize,
    /// Number of bases to ignore at the end of each database allele; the end
    /// is rebuilt from the read up to the heptamer
    #[arg(short = 'e', long = "end", default_value_t = 8)]
    end: usize,
    /// Minimum reads for a sequence to appear in the summary
    #[arg(short = 'm', long = "mincount", default_value_t = 1, value_parser = at_least_one)]
    mincount: i64,
    /// Minimum reads of the rebuilt sequence relative to reads of the
    /// shortened one (0-1)
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
    /// Text that the description of a genome record must contain for hits
    /// on it to be kept
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
