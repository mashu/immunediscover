//! `search exact`: count exact occurrences of database alleles in reads.

pub mod border;
pub mod flanks;
pub mod frequencies;
pub mod refgene;
pub mod scan;
pub mod summarize;

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

pub use flanks::{FlankMode, RssElement};
pub use frequencies::Thresholds;
pub use scan::BorderRow;

use crate::commands::report::print_counts_table;
use crate::error::{Error, Result};
use crate::frame::{round_digits, Column, Frame, SortKey};
use crate::io::fasta::read_fasta;
use crate::io::path::always_gz;
use crate::io::tsv::{read_tsv, write_tsv, write_tsv_gz};
use crate::seq::{load_database, load_demultiplex, Allele, DemuxTable, GeneType};
use crate::{info, warn};

/// Settings of the core search (shared with `discover hsmm`).
#[derive(Debug, Clone)]
pub struct ExactSearchSettings {
    pub gene: GeneType,
    /// Minimum number of reads for an allele (and for a flank variant).
    pub mincount: i64,
    /// Minimum count relative to the best allele of the gene.
    pub minratio: f64,
    /// Bases reported on the side without a recombination signal.
    pub affix: usize,
    pub mode: FlankMode,
    /// Maximum number of flank variants reported per allele.
    pub top: usize,
    /// Alleles/genes exempt from ratio thresholds.
    pub expected: Thresholds,
    /// Sequences of a reference set, to label unknown ones as `Novel`.
    pub known_sequences: Option<HashSet<String>>,
    pub border: usize,
    pub adjust_per_gene_extension: bool,
    pub adjust_percent: f64,
    /// Write every unfiltered match to this path (gzip, `.gz` appended).
    pub raw: Option<PathBuf>,
}

impl ExactSearchSettings {
    /// The defaults of the original library function `exact_search`.
    pub fn new(gene: GeneType) -> Self {
        ExactSearchSettings {
            gene,
            mincount: 10,
            minratio: 0.01,
            affix: 13,
            mode: FlankMode::all_rss(),
            top: 10,
            expected: HashMap::new(),
            known_sequences: None,
            border: 0,
            adjust_per_gene_extension: false,
            adjust_percent: 1.0,
            raw: None,
        }
    }
}

/// Result of the core search.
pub struct ExactSearch {
    /// One row per allele (and flank variant) passing the thresholds.
    pub counts: Frame,
    /// Border-filter statistics per donor and gene (extension mode only).
    pub border: Vec<BorderRow>,
}

fn log_border_statistics(border: &[BorderRow]) {
    if border.is_empty() {
        return;
    }
    let matched: usize = border.iter().map(|r| r.matched_total).sum();
    let rejected: usize = border.iter().map(BorderRow::rejected).sum();
    let percent = if matched > 0 {
        round_digits(100.0 * rejected as f64 / matched as f64, 2)
    } else {
        0.0
    };
    info!(
        "Border filter rejected {rejected} of {matched} potential matches ({}%)",
        crate::frame::format_float(percent)
    );
}

/// Border-filter statistics per donor and gene as a table.
pub fn border_table(border: &[BorderRow]) -> Result<Frame> {
    let ratio = |row: &BorderRow| match row.matched_total {
        0 => 0.0,
        matched => row.rejected() as f64 / matched as f64,
    };
    Frame::from_columns([
        (
            "case",
            Column::from_strings(border.iter().map(|r| r.case.as_str())),
        ),
        (
            "gene",
            Column::from_strings(border.iter().map(|r| r.gene.as_str())),
        ),
        (
            "matched_total",
            Column::from_ints(border.iter().map(|r| r.matched_total as i64)),
        ),
        (
            "accepted_total",
            Column::from_ints(border.iter().map(|r| r.accepted_total as i64)),
        ),
        (
            "rejected_border",
            Column::from_ints(border.iter().map(|r| r.rejected() as i64)),
        ),
        (
            "rejected_ratio",
            Column::from_floats(border.iter().map(ratio)),
        ),
    ])
}

/// `<output>.border.tsv` next to the result table.
fn border_path(output: &Path) -> PathBuf {
    let text = output.to_string_lossy();
    let base = text.strip_suffix(".gz").unwrap_or(&text);
    PathBuf::from(format!("{base}.border.tsv"))
}

/// Search every read for every allele and summarise per donor.
pub fn exact_search(
    table: &DemuxTable,
    db: &[Allele],
    settings: &ExactSearchSettings,
) -> Result<ExactSearch> {
    info!(
        "Using mincount: {}, minratio: {} for genes: {}",
        settings.mincount,
        crate::frame::format_float(settings.minratio),
        settings.gene
    );
    let scanned = scan::scan(
        table,
        db,
        &scan::ScanSettings {
            gene: settings.gene,
            affix: settings.affix,
            mode: &settings.mode,
            border: settings.border,
            adjust_per_gene_extension: settings.adjust_per_gene_extension,
            adjust_percent: settings.adjust_percent,
            keep_raw: settings.raw.is_some(),
        },
    )?;
    log_border_statistics(&scanned.border);
    if scanned.collapsed.nrows() == 0 {
        return Err(Error::invalid(
            "No exact matches found: none of the database sequences occurs in any read",
        ));
    }
    if let (Some(path), Some(raw)) = (&settings.raw, &scanned.raw) {
        let mut gz = path.as_os_str().to_os_string();
        gz.push(".gz");
        write_tsv_gz(Path::new(&gz), raw)?;
    }

    let with_ratios =
        summarize::add_counts_and_ratios(&scanned.collapsed, settings.known_sequences.as_ref())?;
    let counts = summarize::filter_and_rank(
        &with_ratios,
        &summarize::CountFilters {
            mincount: settings.mincount,
            minratio: settings.minratio,
            expected: &settings.expected,
            top: settings.top,
        },
    )?;
    Ok(ExactSearch {
        counts,
        border: scanned.border,
    })
}

/// Parameters of the `search exact` command.
#[derive(Debug, Clone)]
pub struct ExactParams {
    pub tsv: PathBuf,
    pub fasta: PathBuf,
    pub output: PathBuf,
    pub search: ExactSearchSettings,
    pub min_allele_median_ratio: f64,
    pub min_gene_median_ratio: f64,
    pub noplot: bool,
    pub refgenes: Vec<String>,
    /// Use only the first `limit` reads (0 = all).
    pub limit: usize,
    /// TSV (`name`, `ratio`) of per-allele/gene minimum `allele_freq`.
    pub expect: Option<PathBuf>,
    /// TSV (`name`, `ratio`) of per-allele/gene minimum `gene_case_freq`.
    pub deletion: Option<PathBuf>,
    /// Only alleles starting with this prefix enter frequency totals.
    pub locus: String,
    pub ref_fasta: Option<PathBuf>,
}

/// Read a two-column `name`/`ratio` threshold table.
pub fn load_thresholds(path: &Path, what: &str) -> Result<Thresholds> {
    let table = read_tsv(path)?;
    let names = table.column("name")?;
    let ratios = table.floats("ratio")?;
    if table.nrows() > 0 {
        info!("Using {what} file with {} entries", table.nrows());
    }
    names
        .iter()
        .zip(ratios)
        .map(|(name, ratio)| {
            let ratio = ratio.ok_or_else(|| Error::parse(path, "missing ratio value"))?;
            Ok((name.julia_string(), ratio))
        })
        .collect()
}

/// Sequences of a reference FASTA, for the `isin_db` column.
pub fn load_known_sequences(path: &Path) -> Result<HashSet<String>> {
    info!(
        "Building sequence lookup from reference FASTA: {}",
        path.display()
    );
    let sequences: HashSet<String> = read_fasta(path)?.into_iter().map(|r| r.sequence).collect();
    info!("Loaded {} sequences from reference FASTA", sequences.len());
    Ok(sequences)
}

fn plot_genes(counts: &Frame) -> Result<()> {
    let groups = counts.group_by(&["gene"])?;
    let genes = counts.strs("gene")?;
    let mut rows: Vec<(String, usize)> = groups
        .iter()
        .map(|rows| (genes[rows[0]].clone().unwrap_or_default(), rows.len()))
        .collect();
    rows.sort_by_key(|(_, count)| std::cmp::Reverse(*count));
    let (labels, values): (Vec<String>, Vec<usize>) = rows.into_iter().unzip();
    print_counts_table(&labels, &values);
    Ok(())
}

/// Run the command.
pub fn run(params: &ExactParams) -> Result<()> {
    info!("Exact search");
    if !params.refgenes.is_empty() {
        info!("Using reference genes {:?}", params.refgenes);
    }
    let mut table = load_demultiplex(&params.tsv, None)?;
    if params.limit > 0 {
        info!("Limiting reads to {}", params.limit);
        table = table.head(params.limit)?;
    }
    let db = load_database(&params.fasta)?;
    if params.search.mincount < 5 {
        warn!("Decreasing mincount below 5 may lead to false positives");
    }
    match &params.search.mode {
        FlankMode::Extension(length) => info!("Using extension mode with length {length}"),
        FlankMode::Rss(elements) => {
            let names: Vec<&str> = elements.iter().map(|e| e.column()).collect();
            info!("Extract RSS: {}", names.join(","));
        }
    }
    if params.search.top != 1 {
        info!(
            "Uncollapsed mode; at most {} full records returned.",
            params.search.top
        );
    }

    let mut settings = params.search.clone();
    if let Some(path) = &params.expect {
        settings.expected = load_thresholds(path, "expect")?;
    }
    let deletion = match &params.deletion {
        Some(path) => load_thresholds(path, "deletion")?,
        None => Thresholds::new(),
    };
    if let Some(path) = &params.ref_fasta {
        settings.known_sequences = Some(load_known_sequences(path)?);
    }

    let found = exact_search(&table, &db, &settings)?;
    let mut counts = found
        .counts
        .sorted(&[SortKey::asc("case"), SortKey::asc("db_name")])?;

    if !params.noplot {
        if counts.nrows() > 0 {
            plot_genes(&counts)?;
        } else {
            warn!("No exact matches to plot");
        }
    }

    frequencies::add_frequencies(&mut counts, &params.locus)?;
    frequencies::apply_frequency_filters(
        &mut counts,
        &params.locus,
        &settings.expected,
        &deletion,
        params.min_allele_median_ratio,
        params.min_gene_median_ratio,
    )?;

    for refgene in &params.refgenes {
        counts = refgene::add_reference_ratios(&counts, refgene)?;
    }

    let output = always_gz(&params.output);
    if !found.border.is_empty() {
        let path = border_path(&params.output);
        write_tsv(&path, &border_table(&found.border)?)?;
        info!("Border filter statistics saved in {}", path.display());
    }
    write_tsv_gz(&output, &counts)?;
    info!(
        "Exact search data saved in compressed {} file",
        output.display()
    );
    Ok(())
}
