//! `table` command group.

use std::path::PathBuf;

use clap::{Args, Subcommand};

use super::validators::positive_size;
use crate::commands::table::common::split_list;
use crate::commands::table::filter::{Comparison, RowCondition};
use crate::commands::table::join::JoinParams;
use crate::commands::table::to_fasta::ToFastaParams;
use crate::commands::table::{collect, exclude, filter, join, simple, to_fasta, transform};
use crate::error::Result;
use crate::frame::JoinKind;
use crate::info;

#[derive(Debug, Subcommand)]
pub enum TableCommand {
    /// Outer-join two TSVs on specified key columns
    Outerjoin(JoinArgs),
    /// Left-join two TSVs on specified key columns
    Leftjoin(JoinArgs),
    /// Regex-based column transform with capture-group replacement and optional new column
    Transform(TransformArgs),
    /// Group by selected columns and count unique groups (keep optional columns)
    Aggregate(AggregateArgs),
    /// Select distinct rows by specified columns
    Unique(ColumnsArgs),
    /// Sort TSV by one or more columns (asc/desc)
    Sort(SortArgs),
    /// Filter rows by string regex or numeric threshold operation
    Filter(FilterArgs),
    /// Project a subset of columns from a TSV
    Select(ColumnsArgs),
    /// Export sequences from TSV to FASTA with optional filtering/cleanup
    Fasta(Box<ToFastaArgs>),
    /// Concatenate many TSVs into a single table (schema must match)
    Collect(CollectArgs),
    /// Exclude sequences whose names or sequences overlap with a FASTA reference
    Exclude(ExcludeArgs),
}

impl TableCommand {
    pub fn run(self) -> Result<()> {
        match self {
            TableCommand::Outerjoin(args) => args.run(JoinKind::Outer),
            TableCommand::Leftjoin(args) => args.run(JoinKind::Left),
            TableCommand::Transform(args) => args.run(),
            TableCommand::Aggregate(args) => args.run(),
            TableCommand::Unique(args) => {
                simple::unique(&args.input, &args.output, &split_list(&args.columns)).map(|_| ())
            }
            TableCommand::Sort(args) => simple::sort(
                &args.input,
                &args.output,
                &split_list(&args.columns),
                args.reverse,
            )
            .map(|_| ()),
            TableCommand::Filter(args) => args.run(),
            TableCommand::Select(args) => {
                simple::select(&args.input, &args.output, &split_list(&args.columns)).map(|_| ())
            }
            TableCommand::Fasta(args) => args.run(),
            TableCommand::Collect(args) => collect::run(&args.pattern, &args.output),
            TableCommand::Exclude(args) => exclude::run(
                &args.input,
                &args.output,
                &args.fasta,
                &args.colname,
                &args.colseq,
            ),
        }
    }
}

#[derive(Debug, Args)]
pub struct JoinArgs {
    /// Left TSV file path
    left: PathBuf,
    /// Right TSV file path
    right: PathBuf,
    /// Output TSV (gz auto-enabled)
    output: PathBuf,
    /// Comma-separated column names to join on
    #[arg(short = 'k', long = "keys", required_unless_present = "left_keys")]
    keys: Option<String>,
    /// Comma-separated join keys from left file (defaults to --keys)
    #[arg(long = "left-keys")]
    left_keys: Option<String>,
    /// Comma-separated join keys from right file (defaults to --keys)
    #[arg(long = "right-keys")]
    right_keys: Option<String>,
    /// Optional prefix for left non-key columns
    #[arg(long = "left-prefix")]
    left_prefix: Option<String>,
    /// Optional prefix for right non-key columns
    #[arg(long = "right-prefix")]
    right_prefix: Option<String>,
    /// Comma-separated subset of columns to keep from left file
    #[arg(long = "left-select")]
    left_select: Option<String>,
    /// Comma-separated subset of columns to keep from right file
    #[arg(long = "right-select")]
    right_select: Option<String>,
}

impl JoinArgs {
    fn run(self, kind: JoinKind) -> Result<()> {
        // Left keys: `--left-keys`, else `--keys`. Right keys: `--right-keys`,
        // else `--keys`, else the left keys.
        let left_keys = self
            .left_keys
            .or_else(|| self.keys.clone())
            .expect("clap requires --keys or --left-keys");
        join::run(&JoinParams {
            left: self.left,
            right: self.right,
            output: self.output,
            kind,
            left_keys: split_list(&left_keys),
            right_keys: self.right_keys.or(self.keys).as_deref().map(split_list),
            left_prefix: self.left_prefix,
            right_prefix: self.right_prefix,
            left_select: self.left_select.as_deref().map(split_list),
            right_select: self.right_select.as_deref().map(split_list),
        })
        .map(|_| ())
    }
}

#[derive(Debug, Args)]
pub struct TransformArgs {
    /// Input TSV file path
    input: PathBuf,
    /// Output TSV (gz auto-enabled)
    output: PathBuf,
    /// Target column to transform
    #[arg(short = 'c', long = "column")]
    column: String,
    /// Regex with capture groups (e.g., 'ID_(\d+)_(\w+)')
    #[arg(short = 'p', long = "pattern", allow_hyphen_values = true)]
    pattern: String,
    /// Replacement string using capture groups (e.g., '\1-\2')
    #[arg(short = 'r', long = "replacement", allow_hyphen_values = true)]
    replacement: String,
    /// Optional name for a new column storing captured groups
    #[arg(long = "new-column")]
    new_column: Option<String>,
}

impl TransformArgs {
    fn run(self) -> Result<()> {
        info!("Transforming TSV file");
        transform::run(
            &self.input,
            &self.output,
            &self.column,
            &self.pattern,
            &self.replacement,
            self.new_column.as_deref(),
        )
        .map(|_| ())
    }
}

#[derive(Debug, Args)]
pub struct AggregateArgs {
    /// Input TSV file path
    input: PathBuf,
    /// Output TSV (gz auto-enabled)
    output: PathBuf,
    /// Comma-separated column names to group by
    #[arg(short = 'g', long = "group-by")]
    group_by: String,
    /// Comma-separated additional columns to keep (defaults to all non-group columns)
    #[arg(short = 'k', long = "keep-columns")]
    keep_columns: Option<String>,
    /// Name of the count column
    #[arg(short = 'c', long = "count-column", default_value = "count")]
    count_column: String,
}

impl AggregateArgs {
    fn run(self) -> Result<()> {
        info!("Aggregating TSV file");
        let keep = self.keep_columns.as_deref().map(split_list);
        simple::aggregate(
            &self.input,
            &self.output,
            &split_list(&self.group_by),
            keep.as_deref(),
            &self.count_column,
        )
        .map(|_| ())
    }
}

#[derive(Debug, Args)]
pub struct ColumnsArgs {
    /// Input TSV file path
    input: PathBuf,
    /// Output TSV (gz auto-enabled)
    output: PathBuf,
    /// Comma-separated list of columns
    #[arg(short = 'c', long = "columns")]
    columns: String,
}

#[derive(Debug, Args)]
pub struct SortArgs {
    /// Input TSV file path
    input: PathBuf,
    /// Output TSV (gz auto-enabled)
    output: PathBuf,
    /// Comma-separated columns to sort by (priority order)
    #[arg(short = 'c', long = "columns")]
    columns: String,
    /// Sort in descending order (default ascending)
    #[arg(short = 'r', long = "reverse")]
    reverse: bool,
}

#[derive(Debug, Args)]
pub struct FilterArgs {
    /// Input TSV file path
    input: PathBuf,
    /// Output TSV (gz auto-enabled)
    output: PathBuf,
    /// Target column to filter on
    #[arg(short = 'c', long = "column")]
    column: String,
    /// Regex pattern for string filtering (for text columns)
    #[arg(long = "pattern", allow_hyphen_values = true)]
    pattern: Option<String>,
    /// Numeric operator: <, <=, >=, > (requires --threshold)
    #[arg(long = "operator")]
    operator: Option<Comparison>,
    /// Numeric threshold value (used with --operator for numeric columns)
    #[arg(long = "threshold", allow_negative_numbers = true)]
    threshold: Option<f64>,
}

impl FilterArgs {
    fn run(self) -> Result<()> {
        let condition = RowCondition::from_options(self.pattern, self.operator, self.threshold)?;
        filter::run(&self.input, &self.output, &self.column, &condition).map(|_| ())
    }
}

#[derive(Debug, Args)]
pub struct ToFastaArgs {
    /// Input TSV file path
    input: PathBuf,
    /// Output FASTA file path
    output: PathBuf,
    /// Column name with sequence names/IDs
    #[arg(short = 'n', long = "colname", default_value = "allele_name")]
    colname: String,
    /// Column name(s) with nucleotide sequences (comma-separated for concatenation)
    #[arg(short = 's', long = "colseq", default_value = "seq")]
    colseq: String,
    /// Optional column with descriptions appended to FASTA headers
    #[arg(short = 'd', long = "coldesc")]
    coldesc: Option<String>,
    /// Regex to filter sequence names (e.g., 'Novel')
    #[arg(short = 'f', long = "filter", allow_hyphen_values = true)]
    filter: Option<String>,
    /// Regex to remove from sequence names (e.g., ' Novel')
    #[arg(short = 'c', long = "cleanup", allow_hyphen_values = true)]
    cleanup: Option<String>,
    /// Regex to filter description column; capture group 1 (if present) is appended
    #[arg(long = "desc-filter", allow_hyphen_values = true)]
    desc_filter: Option<String>,
    /// Do not sort records by sequence name
    #[arg(long = "no-sort")]
    no_sort: bool,
    /// Minimum number of cases that must include the allele to export
    #[arg(long = "mincase", default_value_t = 1, value_parser = positive_size)]
    mincase: usize,
    /// Column name with donor/case identifiers
    #[arg(long = "case-col", default_value = "case")]
    case_col: String,
    /// Keep only unique sequences (ignore sequence names, use first encountered name)
    #[arg(long = "unique-sequences")]
    unique_sequences: bool,
}

impl ToFastaArgs {
    fn run(self) -> Result<()> {
        info!("Exporting TSV to FASTA");
        to_fasta::run(&ToFastaParams {
            input: self.input,
            output: self.output,
            name_column: self.colname,
            sequence_columns: self.colseq,
            description_column: self.coldesc,
            filter: self.filter,
            cleanup: self.cleanup,
            description_filter: self.desc_filter,
            sort_by_name: !self.no_sort,
            min_cases: self.mincase,
            case_column: self.case_col,
            unique_sequences: self.unique_sequences,
        })
    }
}

#[derive(Debug, Args)]
pub struct CollectArgs {
    /// Glob pattern of TSV files to concatenate (columns must match)
    pattern: String,
    /// Output TSV file path for concatenated table
    output: PathBuf,
}

#[derive(Debug, Args)]
pub struct ExcludeArgs {
    /// Input TSV with allele_name and seq columns
    input: PathBuf,
    /// Output TSV with excluded rows removed
    output: PathBuf,
    /// Reference FASTA with names/sequences to exclude against
    fasta: PathBuf,
    /// Column name containing allele names
    #[arg(short = 'n', long = "colname", default_value = "allele_name")]
    colname: String,
    /// Column name containing nucleotide sequences
    #[arg(short = 's', long = "colseq", default_value = "seq")]
    colseq: String,
}
