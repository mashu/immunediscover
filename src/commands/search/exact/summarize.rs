//! Turning collapsed match rows into per-allele counts and ratios.

use std::collections::{HashMap, HashSet};

use crate::commands::variants::{
    allele_counts, rank_variants, ratio_to_gene_max, sort_variants, VariantColumns,
};
use crate::error::Result;
use crate::frame::{Column, Frame};
use crate::info;
use crate::seq::gene_of;

/// Thresholds applied to the per-allele counts.
pub struct CountFilters<'a> {
    pub mincount: i64,
    pub minratio: f64,
    /// Alleles or genes exempt from the ratio thresholds.
    pub expected: &'a HashMap<String, f64>,
    /// Maximum number of flank variants reported per allele sequence.
    pub top: usize,
}

fn text(frame: &Frame, name: &str) -> Result<Vec<String>> {
    Ok(frame
        .strs(name)?
        .iter()
        .map(|c| c.clone().unwrap_or_default())
        .collect())
}

fn ints(frame: &Frame, name: &str) -> Result<Vec<i64>> {
    Ok(frame.ints(name)?.iter().map(|c| c.unwrap_or(0)).collect())
}

/// How `search exact` tables name donor, gene and allele.
pub const COLUMNS: VariantColumns<'static> = VariantColumns {
    well: "well",
    case: "case",
    gene: "gene",
    name: "db_name",
    sequence: "sequence",
};

/// Add `count`, `gene`, optional `isin_db`, `full_ratio` and `ratio` to the
/// collapsed rows (which already carry `full_count`).
pub fn add_counts_and_ratios(
    collapsed: &Frame,
    known_sequences: Option<&HashSet<String>>,
) -> Result<Frame> {
    let mut frame = collapsed.clone();
    let full_count = ints(&frame, "full_count")?;

    // `count`: reads per allele sequence, regardless of flanks.
    let count = allele_counts(&frame, &COLUMNS, &full_count)?;
    frame.push_column("count", Column::from_ints(count.iter().copied()))?;

    let gene: Vec<String> = text(&frame, "db_name")?
        .iter()
        .map(|n| gene_of(n).to_string())
        .collect();
    frame.push_column("gene", Column::from_strings(gene))?;

    if let Some(known) = known_sequences {
        info!("Adding isin_db column based on reference FASTA");
        let labels = text(&frame, "sequence")?.into_iter().map(|s| {
            if known.contains(&s) {
                ""
            } else {
                "Novel"
            }
        });
        frame.push_column("isin_db", Column::from_strings(labels))?;
    }

    let full_ratio = ratio_to_gene_max(&frame, &COLUMNS, &full_count)?;
    let ratio = ratio_to_gene_max(&frame, &COLUMNS, &count)?;
    frame.push_column("full_ratio", Column::from_floats(full_ratio))?;
    frame.push_column("ratio", Column::from_floats(ratio))?;
    Ok(frame)
}

/// Ratio threshold for one row: alleles or genes listed in the
/// expectation table are exempt (threshold 0).
fn ratio_threshold(filters: &CountFilters<'_>, db_name: &str, gene: &str, case: &str) -> f64 {
    if filters.expected.contains_key(db_name) {
        info!("Skipping allelic ratio filters for {db_name} in case {case}");
        0.0
    } else if filters.expected.contains_key(gene) {
        info!("Skipping gene ratio filters for {db_name} in case {case}");
        0.0
    } else {
        filters.minratio
    }
}

/// Sort, apply count/ratio thresholds, and number flank variants.
pub fn filter_and_rank(frame: &Frame, filters: &CountFilters<'_>) -> Result<Frame> {
    let mut kept = sort_variants(frame, &COLUMNS)?;

    // First on full (flank-aware) counts, then on per-sequence counts.
    for (count_column, ratio_column) in [("full_count", "full_ratio"), ("count", "ratio")] {
        let counts = ints(&kept, count_column)?;
        let ratios = kept.floats(ratio_column)?;
        let db_name = text(&kept, "db_name")?;
        let gene = text(&kept, "gene")?;
        let case = text(&kept, "case")?;
        let mask: Vec<bool> = (0..kept.nrows())
            .map(|row| {
                counts[row] >= filters.mincount
                    && ratios[row].unwrap_or(f64::NAN)
                        >= ratio_threshold(filters, &db_name[row], &gene[row], &case[row])
            })
            .collect();
        kept = kept.filter(&mask);
    }

    let ranked = rank_variants(&kept, &COLUMNS, filters.top)?;
    main_columns_first(&ranked)
}

/// Identification and count columns first, flanks after them.
fn main_columns_first(frame: &Frame) -> Result<Frame> {
    const MAIN: [&str; 8] = [
        "well",
        "case",
        "gene",
        "db_name",
        "count",
        "full_count",
        "ratio",
        "full_ratio",
    ];
    let mut order: Vec<String> = MAIN.iter().map(|c| c.to_string()).collect();
    order.extend(
        frame
            .names()
            .iter()
            .filter(|name| !MAIN.contains(&name.as_str()))
            .cloned(),
    );
    frame.select(&order)
}
