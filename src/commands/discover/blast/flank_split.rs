//! Detect alleles that pass allele-level abundance but lose every
//! extension variant to `full_count` / `full_ratio`.
//!
//! That pattern appears when a germline sequence has enough total reads
//! (`count`) yet no single extension context reaches `-c` / `-f` — the
//! support is split across flanks (KI: `IGHV1-46*01_S9407`).

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use crate::error::Result;
use crate::frame::{Column, Frame};
use crate::io::path::always_gz;
use crate::io::tsv::write_tsv_gz;
use crate::{info, warn};

/// One allele whose total support clears `-c`/`-f` but every extension
/// fails the same cut on `full_count` / `full_ratio`.
#[derive(Debug, Clone)]
pub struct FlankSplitAllele {
    pub well: String,
    pub case: String,
    pub gene: String,
    pub sseqid: String,
    pub count: i64,
    pub ratio: f64,
    pub best_full_count: i64,
    pub best_full_ratio: f64,
    pub n_variants: usize,
    /// Why the best extension still fails: `full_count`, `full_ratio`, or both.
    pub reason: String,
}

/// Alleles that would be kept by allele-level cuts but dropped by the
/// extension-level cuts (same numeric thresholds).
pub fn find_flank_split_alleles(
    clusters: &Frame,
    min_count: i64,
    min_ratio: f64,
) -> Result<Vec<FlankSplitAllele>> {
    let wells = string_column(clusters, "well")?;
    let cases = string_column(clusters, "case")?;
    let genes = string_column(clusters, "gene")?;
    let subjects = string_column(clusters, "sseqid")?;
    let full_count = clusters
        .ints("full_count")?
        .into_iter()
        .map(|v| v.unwrap_or(0))
        .collect::<Vec<_>>();
    let full_ratio = clusters
        .floats("full_ratio")?
        .into_iter()
        .map(|v| v.unwrap_or(0.0))
        .collect::<Vec<_>>();
    let count = clusters
        .ints("count")?
        .into_iter()
        .map(|v| v.unwrap_or(0))
        .collect::<Vec<_>>();
    let ratio = clusters
        .floats("ratio")?
        .into_iter()
        .map(|v| v.unwrap_or(0.0))
        .collect::<Vec<_>>();

    // (well, case, gene, allele_key) -> row indices
    let groups = clusters.group_by(&["well", "case", "gene", "__allele"])?;
    let mut lost = Vec::new();
    for rows in groups.iter() {
        let first = rows[0];
        let allele_count = count[first];
        let allele_ratio = ratio[first];
        if allele_count < min_count || allele_ratio < min_ratio {
            continue;
        }
        let best_full_count = rows.iter().map(|&r| full_count[r]).max().unwrap_or(0);
        let best_full_ratio = rows
            .iter()
            .map(|&r| full_ratio[r])
            .fold(f64::NEG_INFINITY, f64::max);
        let fail_count = best_full_count < min_count;
        let fail_ratio = best_full_ratio < min_ratio;
        if !fail_count && !fail_ratio {
            continue;
        }
        let reason = match (fail_count, fail_ratio) {
            (true, true) => "full_count+full_ratio",
            (true, false) => "full_count",
            (false, true) => "full_ratio",
            (false, false) => unreachable!(),
        };
        // Prefer the subject id of the best-supported extension.
        let best_row = *rows
            .iter()
            .max_by(|&&a, &&b| {
                full_count[a]
                    .cmp(&full_count[b])
                    .then(full_ratio[a].total_cmp(&full_ratio[b]))
            })
            .unwrap_or(&first);
        lost.push(FlankSplitAllele {
            well: wells[first].clone(),
            case: cases[first].clone(),
            gene: genes[first].clone(),
            sseqid: subjects[best_row].clone(),
            count: allele_count,
            ratio: allele_ratio,
            best_full_count,
            best_full_ratio,
            n_variants: rows.len(),
            reason: reason.to_string(),
        });
    }
    lost.sort_by(|a, b| {
        a.gene
            .cmp(&b.gene)
            .then(a.case.cmp(&b.case))
            .then(b.count.cmp(&a.count))
    });
    Ok(lost)
}

/// Log a gene-wise summary and write a sidecar table next to `output`.
///
/// Returns the path written, if any losses were found.
pub fn report_flank_split_losses(
    lost: &[FlankSplitAllele],
    min_count: i64,
    min_ratio: f64,
    output: &Path,
) -> Result<Option<PathBuf>> {
    if lost.is_empty() {
        return Ok(None);
    }

    let mut by_gene: BTreeMap<&str, Vec<&FlankSplitAllele>> = BTreeMap::new();
    for allele in lost {
        by_gene.entry(allele.gene.as_str()).or_default().push(allele);
    }

    warn!(
        "{} allele(s) meet -c/-f on count/ratio but every extension fails full_count/full_ratio \
         (support split across flanks; threshold count>={min_count}, ratio>={})",
        lost.len(),
        format_ratio(min_ratio)
    );
    for (gene, alleles) in &by_gene {
        let total_count: i64 = alleles.iter().map(|a| a.count).sum();
        warn!(
            "  {gene}: {} allele(s), summed count={total_count}",
            alleles.len()
        );
        for allele in alleles.iter().take(8) {
            warn!(
                "    {} {}: count={} ratio={} best_full_count={} best_full_ratio={} variants={} ({})",
                allele.case,
                allele.sseqid,
                allele.count,
                format_ratio(allele.ratio),
                allele.best_full_count,
                format_ratio(allele.best_full_ratio),
                allele.n_variants,
                allele.reason
            );
        }
        if alleles.len() > 8 {
            warn!("    … and {} more", alleles.len() - 8);
        }
    }

    let path = flank_split_sidecar(output);
    write_tsv_gz(&path, &flank_split_frame(lost)?)?;
    info!(
        "Flank-split allele report ({} rows) written to {}",
        lost.len(),
        path.display()
    );
    Ok(Some(path))
}

/// Path of the flank-split sidecar for a discovery output.
pub fn flank_split_sidecar(output: &Path) -> PathBuf {
    sidecar_gz(output, "flank-split")
}

/// Path of the unfiltered (pre-abundance-filter) sidecar.
pub fn unfiltered_sidecar(output: &Path) -> PathBuf {
    sidecar_gz(output, "unfiltered")
}

fn sidecar_gz(output: &Path, tag: &str) -> PathBuf {
    let output = always_gz(output);
    let text = output.to_string_lossy();
    let path = if let Some(stem) = text.strip_suffix(".tsv.gz") {
        format!("{stem}.{tag}.tsv.gz")
    } else if let Some(stem) = text.strip_suffix(".gz") {
        format!("{stem}.{tag}.gz")
    } else {
        format!("{text}.{tag}.tsv.gz")
    };
    PathBuf::from(path)
}

fn flank_split_frame(lost: &[FlankSplitAllele]) -> Result<Frame> {
    Frame::from_columns([
        (
            "well",
            Column::from_strings(lost.iter().map(|a| a.well.as_str())),
        ),
        (
            "case",
            Column::from_strings(lost.iter().map(|a| a.case.as_str())),
        ),
        (
            "gene",
            Column::from_strings(lost.iter().map(|a| a.gene.as_str())),
        ),
        (
            "sseqid",
            Column::from_strings(lost.iter().map(|a| a.sseqid.as_str())),
        ),
        (
            "count",
            Column::from_ints(lost.iter().map(|a| a.count)),
        ),
        (
            "ratio",
            Column::from_floats(lost.iter().map(|a| a.ratio)),
        ),
        (
            "best_full_count",
            Column::from_ints(lost.iter().map(|a| a.best_full_count)),
        ),
        (
            "best_full_ratio",
            Column::from_floats(lost.iter().map(|a| a.best_full_ratio)),
        ),
        (
            "n_variants",
            Column::from_ints(lost.iter().map(|a| a.n_variants as i64)),
        ),
        (
            "reason",
            Column::from_strings(lost.iter().map(|a| a.reason.as_str())),
        ),
    ])
}

fn string_column(frame: &Frame, name: &str) -> Result<Vec<String>> {
    Ok(frame
        .column(name)?
        .iter()
        .map(|v| v.julia_string())
        .collect())
}

fn format_ratio(x: f64) -> String {
    crate::frame::format_float(crate::frame::round_digits(x, 4))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> Frame {
        // Allele A: count 24, best full_count 15 (flank-split at -c 18).
        // Allele B: count 30, best full_count 26 (kept).
        // Allele C: count 10 (fails allele-level -c 18).
        Frame::from_columns([
            ("well", Column::from_ints([3, 3, 3, 3, 3])),
            (
                "case",
                Column::from_strings(["D1", "D1", "D1", "D1", "D1"]),
            ),
            (
                "gene",
                Column::from_strings(["IGHV1-46", "IGHV1-46", "IGHV1-46", "IGHV1-2", "IGHV1-2"]),
            ),
            (
                "sseqid",
                Column::from_strings(["A*01", "A*01", "A*01", "B*01", "C*01"]),
            ),
            (
                "__allele",
                Column::from_strings(["seqA", "seqA", "seqA", "seqB", "seqC"]),
            ),
            ("full_count", Column::from_ints([15, 5, 4, 26, 10])),
            (
                "full_ratio",
                Column::from_floats([0.5, 0.17, 0.13, 1.0, 0.3]),
            ),
            ("count", Column::from_ints([24, 24, 24, 30, 10])),
            (
                "ratio",
                Column::from_floats([0.8, 0.8, 0.8, 1.0, 0.3]),
            ),
        ])
        .unwrap()
    }

    #[test]
    fn detects_allele_with_count_but_weak_extensions() {
        let lost = find_flank_split_alleles(&sample(), 18, 0.083).unwrap();
        assert_eq!(lost.len(), 1);
        assert_eq!(lost[0].gene, "IGHV1-46");
        assert_eq!(lost[0].count, 24);
        assert_eq!(lost[0].best_full_count, 15);
        assert_eq!(lost[0].n_variants, 3);
        assert_eq!(lost[0].reason, "full_count");
    }

    #[test]
    fn no_loss_when_threshold_at_best_full_count() {
        let lost = find_flank_split_alleles(&sample(), 15, 0.083).unwrap();
        assert!(lost.is_empty());
    }

    #[test]
    fn sidecar_names() {
        assert_eq!(
            flank_split_sidecar(Path::new("out/KI-blast-V.tsv.gz")),
            PathBuf::from("out/KI-blast-V.flank-split.tsv.gz")
        );
        assert_eq!(
            unfiltered_sidecar(Path::new("out/KI-blast-V.tsv")),
            PathBuf::from("out/KI-blast-V.unfiltered.tsv.gz")
        );
    }
}
