//! `discover` command group.

use std::path::PathBuf;

use clap::{ArgAction, Args, Subcommand, ValueEnum};

use crate::error::Error;

use super::presets::{describe_presets, resolve, Preset, PresetValues};
use super::validators::{at_least_one, positive_size, unit_interval};
use crate::commands::discover::blast::workdir::DEFAULT_WORK_DIR;
use crate::commands::discover::blast::{self, Aligner, BlastParams};
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

/// The program that compares the reads with the known alleles.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum AlignerChoice {
    /// NCBI blastn, which must be installed
    Blast,
    /// The aligner built into this program. Needs no other software and
    /// compares the whole allele, so differences at its ends are counted
    Builtin,
}

#[derive(Debug, Args)]
pub struct BlastArgs {
    /// Reads of the plate, one row per read: the table written by
    /// `preprocess demultiplex` (columns well, case, name, genomic_sequence)
    #[arg(required_unless_present = "show_presets")]
    input: Option<PathBuf>,
    /// Known alleles (FASTA, names like IGHV1-69*01)
    #[arg(required_unless_present = "show_presets")]
    fasta: Option<PathBuf>,
    /// Where to write the result table (.gz is added)
    #[arg(required_unless_present = "show_presets")]
    output: Option<PathBuf>,
    /// Pseudogenes (FASTA). Reads that resemble one of them more than any
    /// allele are discarded
    #[arg(short = 'p', long = "pseudo", default_value = "")]
    pseudo: String,
    /// Minimum number of reads for a sequence to be reported in a donor.
    /// Both `full_count` and `count` must reach it [default: 5]
    #[arg(short = 'c', long = "minfullcount")]
    minfullcount: Option<i64>,
    /// Minimum number of reads relative to the best allele of the same gene
    /// in the donor (0-1). Both `full_ratio` and `ratio` must reach it.
    /// Sequencing errors stay far below it; a second allele of the gene is
    /// near 1 [default: 0.1]
    #[arg(short = 'f', long = "minfullratio", value_parser = unit_interval)]
    minfullratio: Option<f64>,
    /// Minimum share (0-1) of the database allele that a read must cover to
    /// be assigned to it [default: 0.1]
    #[arg(short = 's', long = "subjectcov", value_parser = unit_interval)]
    subjectcov: Option<f64>,
    /// Maximum number of differences from the nearest database allele.
    /// Sequences further away are not reported [default: 20]
    #[arg(short = 'd', long = "maxdist", value_parser = clap::value_parser!(i64).range(0..))]
    maxdist: Option<i64>,
    /// Minimum length of the matching part of a read (the allele plus the
    /// bases recorded around it). Shorter matches are fragments
    /// [default: 290]
    #[arg(short = 'l', long = "length", value_parser = positive_size)]
    length: Option<usize>,
    /// Minimum number of bases between the allele and the end of the read.
    /// Removes alleles cut off by the read end [default: 0]
    #[arg(short = 'e', long = "edge", value_parser = clap::value_parser!(i64).range(0..))]
    edge: Option<i64>,
    /// Settings passed on to blastn. Change only if you know BLAST
    /// [default: "-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 50"]
    #[arg(short = 'a', long = "args", allow_hyphen_values = true)]
    args: Option<String>,
    /// Which program compares the reads with the known alleles. `builtin`
    /// needs no BLAST installation; --args does not apply to it
    #[arg(long = "aligner", value_enum, default_value_t = AlignerChoice::Blast)]
    aligner: AlignerChoice,
    /// Search again instead of reusing the BLAST results of an earlier run
    #[arg(short = 'o', long = "overwrite")]
    overwrite: bool,
    /// Use the settings tuned for V, D or J genes (see --show-presets). Any
    /// option given explicitly takes precedence over the preset.
    #[arg(short = 'g', long = "gene")]
    gene: Option<Preset>,
    /// Print the V, D and J presets and exit
    #[arg(short = 'G', long = "show-presets")]
    show_presets: bool,
    /// Number of bases recorded before each allele. Alleles are searched
    /// together with the bases around them, which anchors their ends
    /// [default: 20]
    #[arg(long = "forward")]
    forward: Option<usize>,
    /// Number of bases recorded after each allele [default: 20]
    #[arg(long = "reverse")]
    reverse: Option<usize>,
    /// How many differences a read may have from a database allele and still
    /// be used to learn the bases around that allele. 0 uses only reads
    /// that contain the allele exactly. Raise it (for example to 3) when the
    /// database comes from another species or lacks most alleles of the
    /// donors.
    #[arg(long = "extension-maxdist", default_value_t = 0)]
    extension_maxdist: usize,
    /// Do not lend surrounding bases to alleles found in no read. By default
    /// such an allele takes the surroundings of the other alleles of its
    /// gene, so that all alleles of a gene are searched alike.
    #[arg(long = "no-gene-context")]
    no_gene_context: bool,
    /// How well the bases recorded around an allele must match the read
    /// when they are cut off again (0-1). A sequence whose surroundings
    /// match worse on either side is dropped [default: 0.75]
    #[arg(short = 'q', long = "minquality", value_parser = unit_interval)]
    minquality: Option<f64>,
    /// Report a sequence that is a part of a known allele as new. Without
    /// this flag it takes the name of the known allele that contains it.
    #[arg(short = 'i', long = "isin", action = ArgAction::SetFalse)]
    isin: bool,
    /// Keep sequences whose surrounding bases could not be cut off (their
    /// aln_qseq is empty and aln_mismatch is -1). For troubleshooting;
    /// --min-corecov and --maxdist cannot be applied to them.
    #[arg(long = "keep-failed")]
    keep_failed: bool,
    /// Minimum length of the found sequence relative to the database allele
    /// (0-1). Removes truncated sequences [default: 0.6]
    #[arg(long = "min-corecov", value_parser = unit_interval)]
    min_corecov: Option<f64>,
    /// How many versions of the surrounding bases to report per allele (rows
    /// with the same sequence in different surroundings)
    #[arg(short = 't', long = "top", default_value_t = 1, value_parser = positive_size)]
    top: usize,
    /// Also write every sequence before the thresholds to
    /// `<output>.unfiltered.tsv.gz`. Alleles whose reads are spread over
    /// several surroundings, none passing on its own, are always listed in
    /// `<output>.flank-split.tsv.gz`.
    #[arg(short = 'v', long = "verbose")]
    verbose: bool,
    /// Leave out the evidence columns (db_distance, copies, evidence,
    /// relation, donors, carriers, segregates, usual_parent_ratio). They
    /// describe how much of every reported sequence a donor holds and what
    /// the other donors of the plate say about it. Explained in the user
    /// guide: https://mashu.github.io/immunediscover/evidence/
    #[arg(long = "no-evidence")]
    no_evidence: bool,
    /// Folder for intermediate files. BLAST results are kept there and
    /// reused by later runs; it is safe to delete.
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
            aligner: match self.aligner {
                AlignerChoice::Blast => Aligner::Blast,
                AlignerChoice::Builtin => Aligner::Builtin,
            },
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
    /// Minimum reads of a known D allele relative to the best allele of its
    /// gene in the donor (0-1), for the allele to be used for learning
    #[arg(short = 'r', long = "ratio", default_value_t = 0.2, value_parser = unit_interval)]
    ratio: f64,
    /// Minimum reads of a known D allele for it to be used for learning
    #[arg(short = 'c', long = "mincount", default_value_t = 10, value_parser = at_least_one)]
    mincount: i64,
    /// Minimum confidence (0-1) of a detected D gene
    #[arg(long = "min-posterior", default_value_t = 0.7, value_parser = unit_interval)]
    min_posterior: f64,
    /// Shortest D gene considered (0 = set automatically)
    #[arg(long = "min-gene-len", default_value_t = 10)]
    min_gene_len: usize,
    /// Longest D gene considered (0 = set automatically)
    #[arg(long = "max-gene-len", default_value_t = 70)]
    max_gene_len: usize,
    /// Use only the first N reads (0 = all). For quick tests
    #[arg(short = 'l', long = "limit", default_value_t = 0)]
    limit: usize,
    /// Minimum reads for a found D sequence to be reported
    #[arg(long = "out-mincount", default_value_t = 10, value_parser = at_least_one)]
    out_mincount: i64,
    /// Minimum reads of a found D sequence relative to the best of its gene
    /// in the donor (0-1), to be reported
    #[arg(long = "out-minratio", default_value_t = 0.2, value_parser = unit_interval)]
    out_minratio: f64,
    /// Minimum quality (0-1) of the heptamer before the gene; 0 switches the
    /// check off
    #[arg(long = "min-heptamer-prob-pre", default_value_t = 0.05, value_parser = unit_interval)]
    min_heptamer_prob_pre: f64,
    /// Minimum quality (0-1) of the heptamer after the gene; 0 switches the
    /// check off
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
