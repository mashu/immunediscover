//! `discover hsmm`: detect D genes with a hidden semi-Markov model of the
//! gene and its recombination signals.
//!
//! The model is trained on the donor-supported known D alleles found by
//! exact search, then every read is scanned for the most probable
//! signal-gene-signal placement. Detected genes are collapsed per well and
//! donor and named after the nearest database allele.

pub mod model;
pub mod scan;

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use rayon::prelude::*;

use self::scan::{DGeneRssHsmm, TrainingExample};
use crate::align::levenshtein;
use crate::commands::search::exact::{exact_search, ExactSearchSettings};
use crate::error::{Error, Result};
use crate::filters::{GermlineFilter, MinThreshold};
use crate::frame::{Column, Frame, SortKey};
use crate::io::path::always_gz;
use crate::io::tsv::write_tsv_gz;
use crate::seq::slice::range;
use crate::seq::{
    gene_of, load_database, load_demultiplex, unique_name, Allele, DemuxTable, GeneType,
};
use crate::{info, warn};

const FLANK_COLUMNS: [&str; 6] = [
    "pre_nonamer",
    "pre_spacer",
    "pre_heptamer",
    "post_heptamer",
    "post_spacer",
    "post_nonamer",
];

/// Lengths of the signal elements, in `FLANK_COLUMNS` order.
const FLANK_LENGTHS: [usize; 6] = [9, 12, 7, 7, 12, 9];

#[derive(Debug, Clone)]
pub struct HsmmParams {
    pub tsv: PathBuf,
    pub fasta: PathBuf,
    pub output: PathBuf,
    /// Minimum count relative to the donor's best allele of the gene for a
    /// known allele to be used in training.
    pub ratio: f64,
    /// Minimum read count of a known allele to be used in training.
    pub mincount: i64,
    pub min_posterior: f64,
    /// Shortest gene considered (0 = shortest training gene).
    pub min_gene_len: usize,
    /// Longest gene considered (0 = longest training gene).
    pub max_gene_len: usize,
    /// Use only the first `limit` reads (0 = all).
    pub limit: usize,
    pub out_mincount: i64,
    pub out_minratio: f64,
    /// Minimum probability of the upstream heptamer (0 disables).
    pub min_heptamer_prob_pre: f64,
    /// Minimum probability of the downstream heptamer (0 disables).
    pub min_heptamer_prob_post: f64,
}

fn text(frame: &Frame, column: &str) -> Result<Vec<String>> {
    Ok(frame
        .column(column)?
        .iter()
        .map(|v| v.julia_string())
        .collect())
}

/// Select training examples from the exact matches of known D alleles.
///
/// Per donor and gene, alleles reaching `ratio` of the best one are
/// candidates; rows of the two best-supported sequences are kept and paired
/// with every distinct flank combination observed for that allele.
///
/// An allele's support in a donor is its read count (`count` is the allele
/// total and repeats on each flank-variant row, so it is taken once per
/// well). Flank combinations cut short by a read end cannot train the
/// fixed-length signal models and are skipped.
fn training_examples(known: &Frame, ratio: f64) -> Result<Vec<TrainingExample>> {
    let by = ["case", "gene", "sequence", "db_name"];
    let groups = known.group_by(&by)?;
    let count: Vec<i64> = known
        .ints("count")?
        .iter()
        .map(|c| c.unwrap_or(0))
        .collect();
    let wells = text(known, "well")?;
    let case_count: Vec<i64> = groups
        .iter()
        .map(|rows| {
            let mut seen = HashSet::new();
            rows.iter()
                .filter(|&&r| seen.insert(wells[r].as_str()))
                .map(|&r| count[r])
                .sum()
        })
        .collect();
    let mut aggregated = known.take(&groups.first_rows()).select(&by)?;
    aggregated.push_column("case_count", Column::from_ints(case_count.iter().copied()))?;

    let case_ratio: Vec<f64> = aggregated.group_by(&["case", "gene"])?.map_rows(|rows| {
        let max = rows.iter().map(|&r| case_count[r]).max().unwrap_or(1) as f64;
        rows.iter().map(|&r| case_count[r] as f64 / max).collect()
    });
    let supported: Vec<bool> = case_ratio.iter().map(|r| *r >= ratio).collect();
    let aggregated = aggregated.filter(&supported);

    // Keep rows of the two best-supported sequences per donor and gene.
    let sequences = text(&aggregated, "sequence")?;
    let names = text(&aggregated, "db_name")?;
    let mut selected: Vec<(String, String)> = Vec::new();
    for rows in aggregated.group_by(&["case", "gene"])?.iter() {
        let sorted = aggregated
            .take(rows)
            .sort_permutation(&[SortKey::desc("case_count")])?;
        let ordered: Vec<usize> = sorted.into_iter().map(|i| rows[i]).collect();
        let mut top: Vec<&str> = Vec::new();
        for &row in &ordered {
            if top.len() < 2 && !top.contains(&sequences[row].as_str()) {
                top.push(&sequences[row]);
            }
        }
        for &row in &ordered {
            if top.contains(&sequences[row].as_str()) {
                selected.push((names[row].clone(), sequences[row].clone()));
            }
        }
    }

    // Distinct flank combinations per (allele name, sequence).
    let mut flank_columns: Vec<&str> = vec!["db_name", "sequence"];
    flank_columns.extend(FLANK_COLUMNS);
    let flanks = known.select(&flank_columns)?.unique();
    let column = |name: &str| text(&flanks, name);
    let (names, genes) = (column("db_name")?, column("sequence")?);
    let flank_text: Vec<Vec<String>> = FLANK_COLUMNS
        .iter()
        .map(|c| column(c))
        .collect::<Result<_>>()?;
    let mut by_allele: HashMap<(&str, &str), Vec<usize>> = HashMap::new();
    for (row, key) in names.iter().zip(&genes).enumerate() {
        by_allele
            .entry((key.0.as_str(), key.1.as_str()))
            .or_default()
            .push(row);
    }

    let complete = |row: usize| {
        FLANK_LENGTHS
            .iter()
            .zip(&flank_text)
            .all(|(length, column)| column[row].len() == *length)
    };
    let mut examples = Vec::new();
    let mut skipped = 0usize;
    for (name, sequence) in &selected {
        for &row in by_allele
            .get(&(name.as_str(), sequence.as_str()))
            .into_iter()
            .flatten()
        {
            if !complete(row) {
                skipped += 1;
                continue;
            }
            let flank = |i: usize| flank_text[i][row].clone();
            examples.push(TrainingExample {
                pre_nonamer: flank(0),
                pre_spacer: flank(1),
                pre_heptamer: flank(2),
                gene: sequence.clone(),
                post_heptamer: flank(3),
                post_spacer: flank(4),
                post_nonamer: flank(5),
            });
        }
    }
    if skipped > 0 {
        info!("Skipped {skipped} training examples with flanks cut short by a read end");
    }
    Ok(examples)
}

/// One detected D segment.
struct ScanRow {
    read: usize,
    sequence: String,
    flanks: [String; 6],
    heptamer_logp_pre: f64,
    heptamer_logp_post: f64,
    log_path_prob: f64,
    log_total_prob: f64,
    posterior_prob: f64,
    isin_db: bool,
    db_name: String,
    nearest_db: String,
    nearest_db_dist: i64,
}

fn scan_read(
    read: usize,
    sequence: &str,
    model: &DGeneRssHsmm,
    db: &[Allele],
    by_sequence: &HashMap<&str, &str>,
) -> Option<ScanRow> {
    let detection = model.scan(sequence)?;
    let part = |first: usize, len: usize| {
        range(sequence, first as i64, (first + len - 1) as i64).to_string()
    };
    let (ps, ss) = (detection.prefix_start, detection.suffix_start);
    let flanks = [
        part(ps, 9),
        part(ps + 9, 12),
        part(ps + 21, 7),
        part(ss, 7),
        part(ss + 7, 12),
        part(ss + 19, 9),
    ];
    let gene = range(
        sequence,
        detection.gene_start as i64,
        detection.gene_end as i64,
    )
    .to_string();

    let known_name = by_sequence.get(gene.as_str());
    let mut nearest = "";
    let mut nearest_dist = i64::MAX;
    for allele in db {
        let distance = levenshtein(gene.as_bytes(), allele.sequence.as_bytes()) as i64;
        if distance < nearest_dist {
            nearest_dist = distance;
            nearest = &allele.name;
            if distance == 0 {
                break;
            }
        }
    }
    Some(ScanRow {
        read,
        heptamer_logp_pre: model.pre_heptamer.log_prob(&flanks[2]),
        heptamer_logp_post: model.post_heptamer.log_prob(&flanks[3]),
        sequence: gene,
        flanks,
        log_path_prob: detection.log_path_prob,
        log_total_prob: detection.log_total_prob,
        posterior_prob: detection.posterior_prob,
        isin_db: known_name.is_some(),
        db_name: known_name.map_or("Novel", |n| n).to_string(),
        nearest_db: nearest.to_string(),
        nearest_db_dist: nearest_dist,
    })
}

fn detections_to_frame(table: &DemuxTable, rows: &[ScanRow]) -> Result<Frame> {
    let mut frame = Frame::from_columns([
        (
            "well",
            Column::from_strings(rows.iter().map(|r| table.well[r.read].as_str())),
        ),
        (
            "case",
            Column::from_strings(rows.iter().map(|r| table.case[r.read].as_str())),
        ),
        (
            "sequence",
            Column::from_strings(rows.iter().map(|r| r.sequence.as_str())),
        ),
    ])?;
    for (i, name) in FLANK_COLUMNS.iter().enumerate() {
        frame.push_column(
            *name,
            Column::from_strings(rows.iter().map(|r| r.flanks[i].as_str())),
        )?;
    }
    frame.push_column(
        "heptamer_logp_pre",
        Column::from_floats(rows.iter().map(|r| r.heptamer_logp_pre)),
    )?;
    frame.push_column(
        "heptamer_logp_post",
        Column::from_floats(rows.iter().map(|r| r.heptamer_logp_post)),
    )?;
    frame.push_column(
        "log_path_prob",
        Column::from_floats(rows.iter().map(|r| r.log_path_prob)),
    )?;
    frame.push_column(
        "log_total_prob",
        Column::from_floats(rows.iter().map(|r| r.log_total_prob)),
    )?;
    frame.push_column(
        "posterior_prob",
        Column::from_floats(rows.iter().map(|r| r.posterior_prob)),
    )?;
    frame.push_column(
        "isin_db",
        Column::from_bools(rows.iter().map(|r| r.isin_db)),
    )?;
    frame.push_column(
        "db_name",
        Column::from_strings(rows.iter().map(|r| r.db_name.as_str())),
    )?;
    frame.push_column(
        "nearest_db",
        Column::from_strings(rows.iter().map(|r| r.nearest_db.as_str())),
    )?;
    frame.push_column(
        "nearest_db_dist",
        Column::from_ints(rows.iter().map(|r| r.nearest_db_dist)),
    )?;
    Ok(frame)
}

/// Collapse detections per well, donor and gene sequence. Per-read
/// columns take the value of the first read of the group.
fn collapse(detections: &Frame) -> Result<Frame> {
    let groups = detections.group_by(&["well", "case", "sequence"])?;
    let isin: Vec<bool> = match detections.column("isin_db")? {
        Column::Bool(values) => values.iter().map(|v| v.unwrap_or(false)).collect(),
        _ => return Err(Error::invalid("isin_db must be boolean")),
    };
    let distance: Vec<i64> = detections
        .ints("nearest_db_dist")?
        .iter()
        .map(|d| d.unwrap_or(i64::MAX))
        .collect();

    let first = detections.take(&groups.first_rows());
    let mut out = first.select(&["well", "case", "sequence"])?;
    out.push_column(
        "count",
        Column::from_ints(groups.iter().map(|rows| rows.len() as i64)),
    )?;
    for name in FLANK_COLUMNS.iter().chain(&[
        "heptamer_logp_pre",
        "heptamer_logp_post",
        "log_path_prob",
        "log_total_prob",
        "posterior_prob",
    ]) {
        out.push_column(*name, first.column(name)?.clone())?;
    }
    out.push_column(
        "isin_db",
        Column::from_bools(groups.iter().map(|rows| rows.iter().any(|&r| isin[r]))),
    )?;
    out.push_column("db_name", first.column("db_name")?.clone())?;
    out.push_column("nearest_db", first.column("nearest_db")?.clone())?;
    out.push_column(
        "nearest_db_dist",
        Column::from_ints(
            groups
                .iter()
                .map(|rows| rows.iter().map(|&r| distance[r]).min().unwrap_or(i64::MAX)),
        ),
    )?;
    Ok(out)
}

/// Add names, heptamer probabilities and gene; apply the output filters.
fn annotate_and_filter(collapsed: &mut Frame, params: &HsmmParams) -> Result<()> {
    let nearest = text(collapsed, "nearest_db")?;
    let db_name = text(collapsed, "db_name")?;
    let sequence = text(collapsed, "sequence")?;
    let distance = collapsed.ints("nearest_db_dist")?.to_vec();
    let isin: Vec<bool> = match collapsed.column("isin_db")? {
        Column::Bool(values) => values.iter().map(|v| v.unwrap_or(false)).collect(),
        _ => return Err(Error::invalid("isin_db must be boolean")),
    };
    let names: Vec<String> = (0..collapsed.nrows())
        .map(|row| {
            let base = if nearest[row].is_empty() {
                &db_name[row]
            } else {
                &nearest[row]
            };
            if distance[row] == Some(0) || isin[row] {
                base.clone()
            } else {
                unique_name(base, &sequence[row])
            }
        })
        .collect();
    collapsed.push_column("allele_name", Column::from_strings(names))?;

    for (source, target) in [
        ("heptamer_logp_pre", "heptamer_prob_pre"),
        ("heptamer_logp_post", "heptamer_prob_post"),
    ] {
        let probabilities: Vec<f64> = collapsed
            .floats(source)?
            .into_iter()
            .map(|l| l.filter(|l| l.is_finite()).map_or(0.0, f64::exp))
            .collect();
        collapsed.push_column(target, Column::from_floats(probabilities))?;
    }

    let genes: Vec<String> = nearest
        .iter()
        .map(|n| {
            if n.is_empty() {
                String::new()
            } else {
                gene_of(n).to_string()
            }
        })
        .collect();
    let any_gene = genes.iter().any(|g| !g.is_empty());
    collapsed.push_column("gene", Column::from_strings(genes))?;
    if any_gene {
        let count: Vec<i64> = collapsed
            .ints("count")?
            .iter()
            .map(|c| c.unwrap_or(0))
            .collect();
        let ratio: Vec<f64> = collapsed
            .group_by(&["well", "case", "gene"])?
            .map_rows(|rows| {
                let max = rows.iter().map(|&r| count[r]).max().unwrap_or(1) as f64;
                rows.iter().map(|&r| count[r] as f64 / max).collect()
            });
        collapsed.push_column("ratio", Column::from_floats(ratio))?;
        GermlineFilter::new()
            .with(MinThreshold::new(
                "count",
                params.out_mincount as f64,
                "Min output count",
            ))
            .with(MinThreshold::new(
                "ratio",
                params.out_minratio,
                "Min output ratio",
            ))
            .apply(collapsed)?;
    }

    let mut heptamer_filter = GermlineFilter::new();
    if params.min_heptamer_prob_pre > 0.0 {
        heptamer_filter.push(MinThreshold::new(
            "heptamer_prob_pre",
            params.min_heptamer_prob_pre,
            "Min pre-heptamer prob",
        ));
    }
    if params.min_heptamer_prob_post > 0.0 {
        heptamer_filter.push(MinThreshold::new(
            "heptamer_prob_post",
            params.min_heptamer_prob_post,
            "Min post-heptamer prob",
        ));
    }
    heptamer_filter.apply(collapsed)
}

/// Run the command. Returns the output table, or `None` when no D segment
/// could be detected (nothing is written in that case).
pub fn run(params: &HsmmParams) -> Result<Option<Frame>> {
    info!("HSMM D detection");
    info!("Loading demultiplex: {}", params.tsv.display());
    let table = load_demultiplex(&params.tsv, (params.limit > 0).then_some(params.limit))?;
    info!("Loaded {} rows", table.len());
    info!("Loading D allele FASTA: {}", params.fasta.display());
    let db = load_database(&params.fasta)?;
    // For identical sequences the last database name wins.
    let by_sequence: HashMap<&str, &str> = db
        .iter()
        .map(|a| (a.sequence.as_str(), a.name.as_str()))
        .collect();

    info!("Searching known D alleles");
    let mut settings = ExactSearchSettings::new(GeneType::D);
    settings.mincount = params.mincount;
    settings.minratio = 0.0;
    settings.top = 1000;
    let known = exact_search(&table, &db, &settings)?.counts;
    if known.nrows() == 0 {
        warn!("No exact D matches");
        return Ok(None);
    }

    let examples = training_examples(&known, params.ratio)?;
    if examples.is_empty() {
        return Err(Error::invalid(
            "No known D allele qualifies as a training example",
        ));
    }
    let lengths: HashSet<usize> = examples.iter().map(|e| e.gene.len()).collect();
    let min_len = if params.min_gene_len > 0 {
        params.min_gene_len
    } else {
        *lengths.iter().min().expect("examples are non-empty")
    };
    let max_len = if params.max_gene_len > 0 {
        params.max_gene_len
    } else {
        *lengths.iter().max().expect("examples are non-empty")
    };
    info!(
        "Training HSMM [{min_len}, {max_len}] from {} examples",
        examples.len()
    );
    let model = DGeneRssHsmm::fit(&examples, min_len, max_len)?;

    info!("Scanning reads with HSMM");
    let detected: Vec<ScanRow> = table
        .sequence
        .par_iter()
        .enumerate()
        .filter_map(|(read, sequence)| scan_read(read, sequence, &model, &db, &by_sequence))
        .collect();
    if detected.is_empty() {
        warn!("No D segments detected");
        return Ok(None);
    }

    let mut detections = detections_to_frame(&table, &detected)?;
    let before = detections.nrows();
    let confident: Vec<bool> = detected
        .iter()
        .map(|r| r.posterior_prob >= params.min_posterior)
        .collect();
    detections = detections.filter(&confident);
    info!("Posterior filter: {}/{before}", detections.nrows());

    let mut collapsed = collapse(&detections)?;
    annotate_and_filter(&mut collapsed, params)?;

    let output = always_gz(&params.output);
    write_tsv_gz(&output, &collapsed)?;
    info!("Saved to {} ({} rows)", output.display(), collapsed.nrows());
    Ok(Some(collapsed))
}
