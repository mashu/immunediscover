//! Scanning reads for exact allele matches and collapsing identical hits.

use std::collections::HashMap;

use rayon::prelude::*;

use super::border::{
    allowed_lengths, extension_overlaps_border, rejects_with_calibration, safe_quantile, MatchSpan,
    PerGeneExtension,
};
use super::flanks::{extension_flanks, rss_flanks, select_rss_fields, FlankMode, Flanks};
use crate::error::Result;
use crate::frame::{Column, Frame};
use crate::info;
use crate::seq::search::PatternSet;
use crate::seq::{gene_of, Allele, DemuxTable, GeneType};

/// Reads processed per parallel batch; bounds peak memory.
const CHUNK: usize = 50_000;

/// How matches are located and what is extracted around them.
pub struct ScanSettings<'a> {
    pub gene: GeneType,
    pub affix: usize,
    pub mode: &'a FlankMode,
    pub border: usize,
    pub adjust_per_gene_extension: bool,
    pub adjust_percent: f64,
    /// Keep every individual match (for the `--raw` diagnostic table).
    pub keep_raw: bool,
}

impl ScanSettings<'_> {
    /// Border filtering only applies in extension mode with a border.
    fn border_extension(&self) -> Option<usize> {
        self.mode.extension().filter(|_| self.border > 0)
    }
}

/// Border-filter bookkeeping for one donor and gene.
#[derive(Debug, Clone, PartialEq)]
pub struct BorderRow {
    pub case: String,
    pub gene: String,
    pub matched_total: usize,
    pub accepted_total: usize,
}

impl BorderRow {
    pub fn rejected(&self) -> usize {
        self.matched_total - self.accepted_total
    }
}

/// Distinct match rows with their multiplicities.
pub struct ScanResult {
    /// One row per distinct `(well, case, db_name, flanks...)` combination,
    /// in order of first appearance, with a `full_count` column.
    pub collapsed: Frame,
    /// Every match in read order (only when requested).
    pub raw: Option<Frame>,
    pub border: Vec<BorderRow>,
}

/// `(allele, longest safe upstream extension, longest safe downstream extension)`.
type AllowedLengths = (usize, Option<usize>, Option<usize>);

/// One match found in a read.
struct Hit {
    allele: usize,
    flanks: Flanks,
    /// `Some(accepted)` when the border filter looked at this match.
    border_verdict: Option<bool>,
}

struct Collapser {
    columns: Vec<&'static str>,
    index: HashMap<(usize, usize, Vec<String>), usize>,
    /// `(read, allele, flank values, count)` per distinct row.
    rows: Vec<(usize, usize, Vec<String>, i64)>,
    read_keys: HashMap<(String, String), usize>,
    name_keys: HashMap<String, usize>,
}

impl Collapser {
    fn new() -> Self {
        Collapser {
            columns: Vec::new(),
            index: HashMap::new(),
            rows: Vec::new(),
            read_keys: HashMap::new(),
            name_keys: HashMap::new(),
        }
    }

    fn add(&mut self, table: &DemuxTable, db: &[Allele], read: usize, hit: Hit) {
        if self.columns.is_empty() {
            self.columns = hit.flanks.iter().map(|(name, _)| *name).collect();
        }
        // Rows are identical when well, case, allele *name* and flanks agree.
        let next = self.read_keys.len();
        let read_key = *self
            .read_keys
            .entry((table.well[read].clone(), table.case[read].clone()))
            .or_insert(next);
        let next = self.name_keys.len();
        let name_key = *self
            .name_keys
            .entry(db[hit.allele].name.clone())
            .or_insert(next);
        let values: Vec<String> = hit.flanks.into_iter().map(|(_, value)| value).collect();
        match self.index.get(&(read_key, name_key, values.clone())) {
            Some(&row) => self.rows[row].3 += 1,
            None => {
                self.index
                    .insert((read_key, name_key, values.clone()), self.rows.len());
                self.rows.push((read, hit.allele, values, 1));
            }
        }
    }
}

fn rows_to_frame(
    table: &DemuxTable,
    db: &[Allele],
    columns: &[&'static str],
    rows: &[(usize, usize, Vec<String>, i64)],
    with_lengths: bool,
    with_counts: bool,
) -> Result<Frame> {
    let mut frame = Frame::from_columns([
        (
            "well",
            Column::from_strings(rows.iter().map(|r| table.well[r.0].as_str())),
        ),
        (
            "case",
            Column::from_strings(rows.iter().map(|r| table.case[r.0].as_str())),
        ),
        (
            "db_name",
            Column::from_strings(rows.iter().map(|r| db[r.1].name.as_str())),
        ),
    ])?;
    for (i, name) in columns.iter().enumerate() {
        frame.push_column(
            *name,
            Column::from_strings(rows.iter().map(|r| r.2[i].as_str())),
        )?;
    }
    if with_lengths {
        for (source, target) in [("prefix", "prefix_len"), ("suffix", "suffix_len")] {
            let i = columns
                .iter()
                .position(|c| *c == source)
                .expect("extension mode has prefix and suffix");
            frame.push_column(
                target,
                Column::from_ints(rows.iter().map(|r| r.2[i].len() as i64)),
            )?;
        }
    }
    if with_counts {
        frame.push_column("full_count", Column::from_ints(rows.iter().map(|r| r.3)))?;
    }
    Ok(frame)
}

/// Choose per-gene extension lengths so that the requested fraction of
/// matches stays clear of the read borders.
fn calibrate(
    table: &DemuxTable,
    db: &[Allele],
    patterns: &PatternSet,
    settings: &ScanSettings<'_>,
    extension: usize,
) -> PerGeneExtension {
    info!(
        "Calibrating per-gene extension targeting ≥ {}% safe reads",
        (settings.adjust_percent * 100.0).round() as i64
    );
    let per_read: Vec<Vec<AllowedLengths>> = table
        .sequence
        .par_iter()
        .map(|read| {
            patterns
                .first_occurrences(read)
                .into_iter()
                .map(|hit| {
                    let span = MatchSpan::new(hit.start, hit.end, read.len());
                    let (pre, suf) =
                        allowed_lengths(span, settings.gene, extension, settings.border);
                    (hit.pattern, pre, suf)
                })
                .collect()
        })
        .collect();

    let mut prefix: HashMap<String, Vec<usize>> = HashMap::new();
    let mut suffix: HashMap<String, Vec<usize>> = HashMap::new();
    for (allele, pre, suf) in per_read.into_iter().flatten() {
        let gene = gene_of(&db[allele].name);
        if let Some(pre) = pre {
            prefix.entry(gene.to_string()).or_default().push(pre);
        }
        if let Some(suf) = suf {
            suffix.entry(gene.to_string()).or_default().push(suf);
        }
    }
    let quantiles = |lengths: HashMap<String, Vec<usize>>| -> HashMap<String, usize> {
        lengths
            .into_iter()
            .map(|(gene, values)| {
                (
                    gene,
                    safe_quantile(&values, settings.adjust_percent, extension),
                )
            })
            .collect()
    };
    let calibrated = PerGeneExtension {
        prefix: quantiles(prefix),
        suffix: quantiles(suffix),
    };
    info!(
        "Per-gene extension calibrated for {} genes",
        calibrated.n_genes()
    );
    calibrated
}

fn hits_in_read(
    read: &str,
    db: &[Allele],
    patterns: &PatternSet,
    settings: &ScanSettings<'_>,
    calibrated: Option<&PerGeneExtension>,
) -> Vec<Hit> {
    let mut hits = Vec::new();
    for occurrence in patterns.first_occurrences(read) {
        let gene_name = gene_of(&db[occurrence.pattern].name);
        let span = MatchSpan::new(occurrence.start, occurrence.end, read.len());
        let mut border_verdict = None;
        if let Some(extension) = settings.border_extension() {
            let rejected = match calibrated {
                Some(calibrated) => rejects_with_calibration(
                    span,
                    settings.gene,
                    calibrated,
                    gene_name,
                    extension,
                    settings.border,
                ),
                None => extension_overlaps_border(span, settings.gene, extension, settings.border),
            };
            border_verdict = Some(!rejected);
        }
        let flanks = match settings.mode {
            FlankMode::Rss(selected) => {
                let all = rss_flanks(
                    read,
                    occurrence.start,
                    occurrence.end,
                    settings.gene,
                    settings.affix,
                );
                select_rss_fields(all, selected, settings.gene)
            }
            FlankMode::Extension(extension) => {
                let (pre, suf) = match calibrated {
                    Some(c) => (
                        c.prefix_for(gene_name, *extension),
                        c.suffix_for(gene_name, *extension),
                    ),
                    None => (*extension, *extension),
                };
                extension_flanks(
                    read,
                    occurrence.start,
                    occurrence.end,
                    settings.gene,
                    settings.affix,
                    pre,
                    suf,
                )
            }
        };
        hits.push(Hit {
            allele: occurrence.pattern,
            flanks,
            border_verdict,
        });
    }
    hits
}

/// Find the first occurrence of every database allele in every read.
pub fn scan(table: &DemuxTable, db: &[Allele], settings: &ScanSettings<'_>) -> Result<ScanResult> {
    let patterns = PatternSet::new(db.iter().map(|a| a.sequence.as_bytes()))?;
    let calibrated = match settings.border_extension() {
        Some(extension) if settings.adjust_per_gene_extension => {
            Some(calibrate(table, db, &patterns, settings, extension))
        }
        _ => None,
    };

    let mut collapser = Collapser::new();
    let mut raw_rows: Vec<(usize, usize, Vec<String>, i64)> = Vec::new();
    let mut border_index: HashMap<(String, String), usize> = HashMap::new();
    let mut border: Vec<BorderRow> = Vec::new();

    for (chunk_number, chunk) in table.sequence.chunks(CHUNK).enumerate() {
        let offset = chunk_number * CHUNK;
        let per_read: Vec<Vec<Hit>> = chunk
            .par_iter()
            .map(|read| hits_in_read(read, db, &patterns, settings, calibrated.as_ref()))
            .collect();
        for (i, hits) in per_read.into_iter().enumerate() {
            let read = offset + i;
            for hit in hits {
                if let Some(accepted) = hit.border_verdict {
                    let key = (
                        table.case[read].clone(),
                        gene_of(&db[hit.allele].name).to_string(),
                    );
                    let slot = *border_index.entry(key.clone()).or_insert_with(|| {
                        border.push(BorderRow {
                            case: key.0,
                            gene: key.1,
                            matched_total: 0,
                            accepted_total: 0,
                        });
                        border.len() - 1
                    });
                    border[slot].matched_total += 1;
                    if !accepted {
                        continue;
                    }
                    border[slot].accepted_total += 1;
                }
                if settings.keep_raw {
                    let values = hit.flanks.iter().map(|(_, v)| v.clone()).collect();
                    raw_rows.push((read, hit.allele, values, 1));
                }
                collapser.add(table, db, read, hit);
            }
        }
    }

    let with_lengths = settings.mode.extension().is_some();
    let collapsed = rows_to_frame(
        table,
        db,
        &collapser.columns,
        &collapser.rows,
        with_lengths,
        true,
    )?;
    let raw = settings
        .keep_raw
        .then(|| {
            rows_to_frame(
                table,
                db,
                &collapser.columns,
                &raw_rows,
                with_lengths,
                false,
            )
        })
        .transpose()?;
    Ok(ScanResult {
        collapsed,
        raw,
        border,
    })
}
