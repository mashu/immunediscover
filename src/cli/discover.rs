//! `discover` command group.

use std::path::PathBuf;

use clap::{ArgAction, Args, Subcommand};

use crate::error::Error;

use super::presets::{describe_presets, resolve, Preset, PresetValues};
use super::validators::{at_least_one, positive_size, unit_interval};
use crate::commands::discover::blast::workdir::DEFAULT_WORK_DIR;
use crate::commands::discover::blast::{self, BlastParams};
use crate::commands::discover::hsmm::{self, HsmmParams};
use crate::error::Result;
use crate::info;

#[derive(Debug, Subcommand)]
pub enum DiscoverCommand {
    /// BLAST-based candidate discovery with trimming, filtering, and identity clustering
    Blast(Box<BlastArgs>),
    /// Detect D genes using an HSMM trained on RSS flanks (V/J masked)
    Hsmm(HsmmArgs),
}

impl DiscoverCommand {
    pub fn run(self) -> Result<()> {
        match self {
            DiscoverCommand::Blast(args) => args.run(),
            DiscoverCommand::Hsmm(args) => args.run(),
        }
    }
}

#[derive(Debug, Args)]
pub struct BlastArgs {
    /// TSV file with demultiplex data
    #[arg(required_unless_present = "show_presets")]
    input: Option<PathBuf>,
    /// FASTA file with database sequences
    #[arg(required_unless_present = "show_presets")]
    fasta: Option<PathBuf>,
    /// TSV file to save discovery results
    #[arg(required_unless_present = "show_presets")]
    output: Option<PathBuf>,
    /// FASTA file with pseudo-genes
    #[arg(short = 'p', long = "pseudo", default_value = "")]
    pseudo: String,
    /// Minimum number of reads, required both of the reported variant
    /// (`full_count`) and of the allele (`count`) [default: 5]
    #[arg(short = 'c', long = "minfullcount")]
    minfullcount: Option<i64>,
    /// Minimum ratio to the best of the gene, required both of
    /// `full_ratio` and of `ratio` [default: 0.1]
    #[arg(short = 'f', long = "minfullratio", value_parser = unit_interval)]
    minfullratio: Option<f64>,
    /// Minimum subject (database) coverage [default: 0.1]
    #[arg(short = 's', long = "subjectcov", value_parser = unit_interval)]
    subjectcov: Option<f64>,
    /// Maximum distance allowed for alleles [default: 20]
    #[arg(short = 'd', long = "maxdist", value_parser = clap::value_parser!(i64).range(0..))]
    maxdist: Option<i64>,
    /// Minimum length of the read segment aligned by BLAST (allele plus
    /// extension) [default: 290]
    #[arg(short = 'l', long = "length", value_parser = positive_size)]
    length: Option<usize>,
    /// Minimum number of nucleotides required between target gene and end of the read [default: 0]
    #[arg(short = 'e', long = "edge", value_parser = clap::value_parser!(i64).range(0..))]
    edge: Option<i64>,
    /// Additional arguments to pass to blastn
    /// [default: "-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 50"]
    #[arg(short = 'a', long = "args", allow_hyphen_values = true)]
    args: Option<String>,
    /// Overwrite existing files (i.e BLAST cache)
    #[arg(short = 'o', long = "overwrite")]
    overwrite: bool,
    /// Use gene preset parameters for V, D, or J analysis. Any option given
    /// explicitly takes precedence over the preset.
    #[arg(short = 'g', long = "gene")]
    gene: Option<Preset>,
    /// Show preset parameters for V, D, and J analysis and exit
    #[arg(short = 'G', long = "show-presets")]
    show_presets: bool,
    /// Forward extension length [default: 20]
    #[arg(long = "forward")]
    forward: Option<usize>,
    /// Reverse extension length [default: 20]
    #[arg(long = "reverse")]
    reverse: Option<usize>,
    /// How far a read may be from a database allele and still supply the
    /// context (flanks) used to extend it. A read containing no allele
    /// exactly is assigned to the allele(s) at the smallest edit distance,
    /// if that distance is at most this value. 0 uses exact matches only;
    /// raise it when the database is from another species or lacks the
    /// donors' alleles.
    #[arg(long = "extension-maxdist", default_value_t = 0)]
    extension_maxdist: usize,
    /// Leave an allele that is found in no read without extension. By
    /// default such an allele takes the context of the other alleles of
    /// its gene, so that all alleles of a gene are searched alike.
    #[arg(long = "no-gene-context")]
    no_gene_context: bool,
    /// Minimum fraction (0-1) of affix positions that match the read in the semi-global affix-read alignment
    /// used for 5'/3' trimming; prefix and suffix each must meet this or the candidate is dropped.
    /// [default: 0.75]
    #[arg(short = 'q', long = "minquality", value_parser = unit_interval)]
    minquality: Option<f64>,
    /// Report a sequence that differs from its best hit but is contained in
    /// a known allele as novel. Without this flag it takes the name of the
    /// known allele that contains it.
    #[arg(short = 'i', long = "isin", action = ArgAction::SetFalse)]
    isin: bool,
    /// Keep rows where trimming failed (aln_qseq empty, aln_mismatch -1). By
    /// default such rows are dropped. Kept rows are exempt from --min-corecov
    /// and --maxdist, which cannot be evaluated for them.
    #[arg(long = "keep-failed")]
    keep_failed: bool,
    /// Minimum ratio length(aln_qseq)/length(db_seq) after trimming [default: 0.6]
    #[arg(long = "min-corecov", value_parser = unit_interval)]
    min_corecov: Option<f64>,
    /// Maximum number of extension variants (rows with the same trimmed
    /// sequence but different read context) reported per allele
    #[arg(short = 't', long = "top", default_value_t = 1, value_parser = positive_size)]
    top: usize,
    /// Print verbose output and save intermediate files, including an
    /// unfiltered cluster table (`*.unfiltered.tsv.gz`) beside the output.
    /// Alleles lost only to full_count/full_ratio (flank-split) are always
    /// summarised in `*.flank-split.tsv.gz` when any exist.
    #[arg(short = 'v', long = "verbose")]
    verbose: bool,
    /// Leave out the evidence columns (db_distance, copies, evidence,
    /// relation, donors, carriers, segregates, usual_parent_ratio). They
    /// describe how much of every reported sequence a donor holds and what
    /// the other donors of the plate say about it; see docs/evidence.md.
    #[arg(long = "no-evidence")]
    no_evidence: bool,
    /// Directory for BLAST cache, temporary query FASTA, combined/extended DB, and affix files
    /// (relative paths use the current directory). Nothing is written beside the input TSV.
    #[arg(long = "work-dir", default_value = DEFAULT_WORK_DIR)]
    work_dir: String,
}

impl BlastArgs {
    fn run(self) -> Result<()> {
        if self.show_presets {
            print!("{}", describe_presets());
            return Ok(());
        }
        let positional = |value: Option<PathBuf>, name: &str| {
            value.ok_or_else(|| Error::invalid(format!("Missing argument <{name}>")))
        };
        let resolved = resolve(
            self.gene,
            PresetValues {
                forward: self.forward,
                reverse: self.reverse,
                min_full_ratio: self.minfullratio,
                length: self.length,
                max_distance: self.maxdist,
                min_full_count: self.minfullcount,
                edge: self.edge,
                subject_coverage: self.subjectcov,
                min_quality: self.minquality,
                min_core_coverage: self.min_corecov,
                blast_args: self.args,
            },
        );
        let params = BlastParams {
            input: positional(self.input, "INPUT")?,
            fasta: positional(self.fasta, "FASTA")?,
            output: positional(self.output, "OUTPUT")?,
            pseudo: (!self.pseudo.is_empty()).then(|| PathBuf::from(&self.pseudo)),
            min_full_count: resolved.min_full_count,
            min_full_ratio: resolved.min_full_ratio,
            subject_coverage: resolved.subject_coverage,
            max_distance: resolved.max_distance,
            min_length: resolved.length,
            edge: resolved.edge,
            blast_args: resolved.blast_args,
            overwrite: self.overwrite,
            forward_extension: resolved.forward,
            reverse_extension: resolved.reverse,
            extension_max_distance: self.extension_maxdist,
            gene_context: !self.no_gene_context,
            min_quality: resolved.min_quality,
            name_substrings_as_known: self.isin,
            keep_failed: self.keep_failed,
            min_core_coverage: resolved.min_core_coverage,
            top: self.top,
            verbose: self.verbose,
            evidence: !self.no_evidence,
            work_dir: self.work_dir,
        };
        info!("Running with the following parameters:");
        info!("  {params:?}");
        blast::run(&params)
    }
}

#[derive(Debug, Args)]
pub struct HsmmArgs {
    /// TSV/TSV.GZ demultiplex file with columns well, case, name, genomic_sequence
    tsv: PathBuf,
    /// FASTA file with D alleles (known reference)
    fasta: PathBuf,
    /// TSV.GZ file to save detected D alleles (novel and/or known) with flanks
    output: PathBuf,
    /// Allelic ratio threshold for known D selection per donor and gene
    #[arg(short = 'r', long = "ratio", default_value_t = 0.2, value_parser = unit_interval)]
    ratio: f64,
    /// Minimum count for a known D allele to be considered in training
    #[arg(short = 'c', long = "mincount", default_value_t = 10, value_parser = at_least_one)]
    mincount: i64,
    /// Minimum posterior probability for accepting an HSMM detection
    #[arg(long = "min-posterior", default_value_t = 0.7, value_parser = unit_interval)]
    min_posterior: f64,
    /// Minimum D gene length for HSMM duration model (auto if 0)
    #[arg(long = "min-gene-len", default_value_t = 10)]
    min_gene_len: usize,
    /// Maximum D gene length for HSMM duration model (auto if 0)
    #[arg(long = "max-gene-len", default_value_t = 70)]
    max_gene_len: usize,
    /// Limit number of demultiplexed reads to process (0 means no limit)
    #[arg(short = 'l', long = "limit", default_value_t = 0)]
    limit: usize,
    /// Minimum count for an extracted D (after HSMM) to keep in output
    #[arg(long = "out-mincount", default_value_t = 10, value_parser = at_least_one)]
    out_mincount: i64,
    /// Minimum allelic ratio within gene (per donor) for output filtering
    #[arg(long = "out-minratio", default_value_t = 0.2, value_parser = unit_interval)]
    out_minratio: f64,
    /// Minimum probability under pre-heptamer PWM to keep detection (0 disables)
    #[arg(long = "min-heptamer-prob-pre", default_value_t = 0.05, value_parser = unit_interval)]
    min_heptamer_prob_pre: f64,
    /// Minimum probability under post-heptamer PWM to keep detection (0 disables)
    #[arg(long = "min-heptamer-prob-post", default_value_t = 0.05, value_parser = unit_interval)]
    min_heptamer_prob_post: f64,
}

impl HsmmArgs {
    fn run(self) -> Result<()> {
        hsmm::run(&HsmmParams {
            tsv: self.tsv,
            fasta: self.fasta,
            output: self.output,
            ratio: self.ratio,
            mincount: self.mincount,
            min_posterior: self.min_posterior,
            min_gene_len: self.min_gene_len,
            max_gene_len: self.max_gene_len,
            limit: self.limit,
            out_mincount: self.out_mincount,
            out_minratio: self.out_minratio,
            min_heptamer_prob_pre: self.min_heptamer_prob_pre,
            min_heptamer_prob_post: self.min_heptamer_prob_post,
        })
        .map(|_| ())
    }
}
