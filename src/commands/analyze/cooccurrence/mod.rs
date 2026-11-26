//! `analyze cooccurrence`: which alleles tend to occur in the same donors.

pub mod cluster;
pub mod stats;

use std::path::{Path, PathBuf};

pub use self::cluster::ClusterMethod;
pub use self::stats::{compute_full_stats, CooccurrenceStats};
use crate::error::Result;
use crate::frame::{Column, Frame};
use crate::info;
use crate::io::tsv::{read_tsv, write_tsv};
use crate::stats::benjamini_hochberg;

#[derive(Debug, Clone)]
pub struct CooccurrenceParams {
    pub input: PathBuf,
    pub case_column: String,
    pub allele_column: String,
    pub min_donors: usize,
    pub cluster_method: ClusterMethod,
    pub cluster_threshold: f64,
    pub debug_triangles: bool,
    pub min_cluster_size: usize,
    /// Where to write clusters; no clustering when absent.
    pub clusters: Option<PathBuf>,
}

/// Edge table of every allele pair sharing at least one donor.
pub fn edges_frame(stats: &CooccurrenceStats) -> Result<Frame> {
    let n = stats.n();
    let pairs: Vec<(usize, usize)> = (0..n)
        .flat_map(|i| ((i + 1)..n).map(move |j| (i, j)))
        .filter(|&(i, j)| stats.support[stats.at(i, j)] > 0)
        .collect();
    let index: Vec<usize> = pairs.iter().map(|&(i, j)| stats.at(i, j)).collect();
    Frame::from_columns([
        (
            "allele_a",
            Column::from_strings(pairs.iter().map(|p| stats.alleles[p.0].as_str())),
        ),
        (
            "allele_b",
            Column::from_strings(pairs.iter().map(|p| stats.alleles[p.1].as_str())),
        ),
        (
            "rho",
            Column::from_floats(index.iter().map(|&k| stats.rho[k])),
        ),
        (
            "jaccard",
            Column::from_floats(index.iter().map(|&k| stats.jaccard[k])),
        ),
        (
            "support",
            Column::from_ints(index.iter().map(|&k| stats.support[k] as i64)),
        ),
        (
            "p_value",
            Column::from_floats(index.iter().map(|&k| stats.p_value[k])),
        ),
    ])
}

/// Edge table restricted by support and Jaccard index, with
/// Benjamini-Hochberg adjusted p-values (`q_enrich`).
pub fn filtered_edges(
    stats: &CooccurrenceStats,
    min_support: usize,
    min_jaccard: f64,
) -> Result<Frame> {
    let n = stats.n();
    let pairs: Vec<(usize, usize)> = (0..n)
        .flat_map(|i| ((i + 1)..n).map(move |j| (i, j)))
        .filter(|&(i, j)| {
            let k = stats.at(i, j);
            stats.support[k] >= min_support && stats.jaccard[k] >= min_jaccard
        })
        .collect();
    let index: Vec<usize> = pairs.iter().map(|&(i, j)| stats.at(i, j)).collect();
    let p: Vec<f64> = index.iter().map(|&k| stats.p_value[k]).collect();
    let donors = |allele: usize| stats.donors_of[&stats.alleles[allele]].len() as i64;
    Frame::from_columns([
        (
            "allele_a",
            Column::from_strings(pairs.iter().map(|p| stats.alleles[p.0].as_str())),
        ),
        (
            "allele_b",
            Column::from_strings(pairs.iter().map(|p| stats.alleles[p.1].as_str())),
        ),
        ("n_a", Column::from_ints(pairs.iter().map(|p| donors(p.0)))),
        ("n_b", Column::from_ints(pairs.iter().map(|p| donors(p.1)))),
        (
            "n_shared",
            Column::from_ints(index.iter().map(|&k| stats.support[k] as i64)),
        ),
        (
            "jaccard",
            Column::from_floats(index.iter().map(|&k| stats.jaccard[k])),
        ),
        (
            "rho",
            Column::from_floats(index.iter().map(|&k| stats.rho[k])),
        ),
        ("p_enrich", Column::from_floats(p.iter().copied())),
        ("q_enrich", Column::from_floats(benjamini_hochberg(&p))),
    ])
}

/// One row per allele of every cluster, with the donors carrying it.
pub fn clusters_frame(stats: &CooccurrenceStats, groups: &[Vec<String>]) -> Result<Frame> {
    let rows: Vec<(i64, &String)> = groups
        .iter()
        .enumerate()
        .flat_map(|(id, group)| group.iter().map(move |allele| (id as i64 + 1, allele)))
        .collect();
    let donors = |allele: &String| -> Vec<&str> {
        stats
            .donors_of
            .get(allele)
            .map(|set| set.iter().map(String::as_str).collect())
            .unwrap_or_default()
    };
    Frame::from_columns([
        ("group_id", Column::from_ints(rows.iter().map(|r| r.0))),
        (
            "allele",
            Column::from_strings(rows.iter().map(|r| r.1.as_str())),
        ),
        (
            "donors",
            Column::from_strings(rows.iter().map(|r| donors(r.1).join(","))),
        ),
        (
            "n_donors",
            Column::from_ints(rows.iter().map(|r| donors(r.1).len() as i64)),
        ),
    ])
}

/// `x.tsv` / `x.tsv.gz` -> `x_edges.tsv`; anything else gets the suffix appended.
fn edges_path(input: &Path) -> PathBuf {
    let text = input.to_string_lossy();
    for extension in [".tsv.gz", ".tsv"] {
        if let Some(stem) = text.strip_suffix(extension) {
            return PathBuf::from(format!("{stem}_edges.tsv"));
        }
    }
    PathBuf::from(format!("{text}_edges.tsv"))
}

fn log_pairs_above_threshold(stats: &CooccurrenceStats, similarity: &[f64], threshold: f64) {
    const MAX_PRINTED: usize = 20;
    let n = stats.n();
    let mut count = 0;
    for i in 0..n {
        for j in (i + 1)..n {
            let s = similarity[i * n + j];
            if s >= threshold {
                count += 1;
                if count <= MAX_PRINTED {
                    info!("  {} — {}: {s:.3}", stats.alleles[i], stats.alleles[j]);
                }
            }
        }
    }
    info!(
        "Total edges above threshold {}: {count}",
        crate::frame::format_float(threshold)
    );
}

/// Run the command.
pub fn run(params: &CooccurrenceParams) -> Result<()> {
    info!("Co-occurrence analysis (Jaccard + support)");
    let table = read_tsv(&params.input)?;
    info!("Loaded {} rows from input file", table.nrows());
    let stats = compute_full_stats(
        &table,
        &params.case_column,
        &params.allele_column,
        params.min_donors,
    )?;

    let edges = edges_frame(&stats)?;
    let edges_output = edges_path(&params.input);
    write_tsv(&edges_output, &edges)?;
    info!(
        "Edges saved to {} ({} edges from {} alleles)",
        edges_output.display(),
        edges.nrows(),
        stats.n()
    );

    let similarity = stats.similarity();
    if params.debug_triangles {
        log_pairs_above_threshold(&stats, &similarity, params.cluster_threshold);
    }

    if let Some(path) = &params.clusters {
        let groups = match params.cluster_method {
            ClusterMethod::Components => cluster::components(
                &similarity,
                &stats.alleles,
                params.cluster_threshold,
                params.min_cluster_size,
            ),
            ClusterMethod::Hierarchical(linkage) => cluster::hierarchical(
                &similarity,
                &stats.alleles,
                params.cluster_threshold,
                params.min_cluster_size,
                linkage,
            )?,
        };
        let frame = clusters_frame(&stats, &groups)?;
        write_tsv(path, &frame)?;
        info!(
            "Clusters saved in {} ({} groups, {} alleles)",
            path.display(),
            groups.len(),
            frame.nrows()
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn edge_output_path() {
        assert_eq!(
            edges_path(Path::new("d/x.tsv")),
            PathBuf::from("d/x_edges.tsv")
        );
        assert_eq!(
            edges_path(Path::new("x.tsv.gz")),
            PathBuf::from("x_edges.tsv")
        );
        assert_eq!(
            edges_path(Path::new("x.txt")),
            PathBuf::from("x.txt_edges.tsv")
        );
    }

    // Same toy data and expectation as the upstream `cooccurrence.jl` test.
    #[test]
    fn shared_donor_counts() {
        let table = Frame::from_columns([
            (
                "case",
                Column::from_strings([
                    "D1", "D1", "D2", "D2", "D3", "D3", "D4", "D4", "D5", "D6", "D7", "D8",
                ]),
            ),
            (
                "db_name",
                Column::from_strings([
                    "A*01", "B*01", "A*01", "B*01", "A*01", "C*01", "B*01", "C*01", "B*01", "C*01",
                    "X*01", "Y*01",
                ]),
            ),
        ])
        .unwrap();
        let stats = compute_full_stats(&table, "case", "db_name", 1).unwrap();
        let edges = filtered_edges(&stats, 1, 0.0).unwrap();
        assert!(edges.nrows() > 0);
        let a = edges.strs("allele_a").unwrap();
        let b = edges.strs("allele_b").unwrap();
        let row = (0..edges.nrows())
            .find(|&r| a[r].as_deref() == Some("A*01") && b[r].as_deref() == Some("B*01"))
            .unwrap();
        assert_eq!(edges.ints("n_shared").unwrap()[row], Some(2));
        assert_eq!(edges.names().last().map(String::as_str), Some("q_enrich"));
    }
}
