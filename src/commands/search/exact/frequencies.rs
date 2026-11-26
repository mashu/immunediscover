//! Cross-donor frequency columns and the filters built on them.

use std::collections::HashMap;

use super::summarize::COLUMNS;
use crate::commands::variants::first_variant_rows;
use crate::error::{Error, Result};
use crate::filters::{CustomFilter, GermlineFilter};
use crate::frame::{Column, Frame};
use crate::info;
use crate::stats::median;

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

/// Per-group sum over the rows whose allele belongs to the locus.
fn locus_sum(frame: &Frame, by: &[&str], in_locus: &[bool], values: &[i64]) -> Result<Vec<i64>> {
    Ok(frame.group_by(by)?.broadcast(|rows| {
        rows.iter()
            .filter(|&&r| in_locus[r])
            .map(|&r| values[r])
            .sum::<i64>()
    }))
}

/// Per-group median over the rows whose allele belongs to the locus;
/// `None` for groups without such rows.
fn locus_median(
    frame: &Frame,
    by: &[&str],
    in_locus: &[bool],
    values: &[i64],
) -> Result<Vec<Option<f64>>> {
    Ok(frame.group_by(by)?.broadcast(|rows| {
        let selected: Vec<f64> = rows
            .iter()
            .filter(|&&r| in_locus[r])
            .map(|&r| values[r] as f64)
            .collect();
        (!selected.is_empty()).then(|| median(&selected))
    }))
}

/// A median column: floating point, except that it degenerates to integer
/// zeros when no group had any locus row (as the original's column does).
fn median_column(medians: &[Option<f64>]) -> Column {
    if medians.iter().all(Option::is_none) {
        Column::from_ints(medians.iter().map(|_| 0))
    } else {
        Column::from_floats(medians.iter().map(|m| m.unwrap_or(0.0)))
    }
}

/// Which rows belong to the locus (allele name starts with `locus`).
fn locus_rows(counts: &Frame, locus: &str) -> Result<Vec<bool>> {
    Ok(text(counts, "db_name")?
        .iter()
        .map(|n| n.starts_with(locus))
        .collect())
}

/// Add gene/donor totals, cross-donor medians and the derived frequencies.
///
/// Rows whose allele name does not start with `locus` (control genes) are
/// excluded from every total and median but still receive the columns.
pub fn add_frequencies(counts: &mut Frame, locus: &str) -> Result<()> {
    info!("Excluding genes not starting with {locus} for frequency calculation");
    // `count` repeats on every flank variant of an allele: totals and
    // medians take each allele once.
    let once = first_variant_rows(counts, &COLUMNS)?;
    let in_locus: Vec<bool> = locus_rows(counts, locus)?
        .iter()
        .zip(&once)
        .map(|(inside, first)| *inside && *first)
        .collect();
    let count = ints(counts, "count")?;

    let gene_count = locus_sum(counts, &["well", "case", "gene"], &in_locus, &count)?;
    let case_count = locus_sum(counts, &["well", "case"], &in_locus, &count)?;
    counts.push_column("gene_count", Column::from_ints(gene_count.iter().copied()))?;
    counts.push_column("case_count", Column::from_ints(case_count.iter().copied()))?;

    let median_count = locus_median(counts, &["gene"], &in_locus, &count)?;
    let median_gene_count = locus_median(counts, &["gene"], &in_locus, &gene_count)?;
    let median_allele_count = locus_median(counts, &["db_name"], &in_locus, &count)?;
    counts.push_column("cross_case_median_count", median_column(&median_count))?;
    counts.push_column(
        "cross_case_median_gene_count",
        median_column(&median_gene_count),
    )?;
    counts.push_column(
        "cross_case_median_allele_count",
        median_column(&median_allele_count),
    )?;

    let divide = |numerator: &[i64], denominator: &[f64]| -> Column {
        Column::from_floats(
            numerator
                .iter()
                .zip(denominator)
                .map(|(n, d)| *n as f64 / d),
        )
    };
    let as_f64 = |values: &[i64]| -> Vec<f64> { values.iter().map(|v| *v as f64).collect() };
    let or_zero =
        |values: &[Option<f64>]| -> Vec<f64> { values.iter().map(|v| v.unwrap_or(0.0)).collect() };

    counts.push_column("allele_case_freq", divide(&count, &as_f64(&case_count)))?;
    counts.push_column("gene_case_freq", divide(&gene_count, &as_f64(&case_count)))?;
    counts.push_column(
        "allele_to_cross_case_median_ratio",
        divide(&count, &or_zero(&median_allele_count)),
    )?;
    counts.push_column(
        "gene_to_cross_case_median_ratio",
        divide(&gene_count, &or_zero(&median_gene_count)),
    )?;

    let allele_freq: Vec<f64> = counts
        .group_by(&["well", "case", "gene"])?
        .map_rows(|rows| {
            let total: i64 = rows.iter().filter(|&&r| once[r]).map(|&r| count[r]).sum();
            rows.iter()
                .map(|&r| count[r] as f64 / total as f64)
                .collect()
        });
    counts.push_column("allele_freq", Column::from_floats(allele_freq))?;
    Ok(())
}

/// Per-allele or per-gene thresholds read from an expectation table.
pub type Thresholds = HashMap<String, f64>;

/// Threshold for a row: allele entry first, then gene entry, else 0.
fn threshold_for(thresholds: &Thresholds, frame: &Frame, row: usize, kind: &str) -> Result<f64> {
    let db_name = frame.value(row, "db_name")?.julia_string();
    let gene = frame.value(row, "gene")?.julia_string();
    match thresholds.get(&db_name).or_else(|| thresholds.get(&gene)) {
        Some(value) => {
            info!("Applying {kind} >= {value} for {db_name}");
            Ok(*value)
        }
        None => Ok(0.0),
    }
}

fn float_at(frame: &Frame, row: usize, column: &str) -> Result<f64> {
    frame
        .value(row, column)?
        .as_f64()
        .ok_or_else(|| Error::invalid(format!("Column '{column}' is not numeric")))
}

/// Apply the expectation, deletion and cross-donor median filters.
///
/// The two cross-donor median filters apply to locus genes only: rows
/// outside the locus (control genes) take no part in the medians, so the
/// ratios are undefined for them and they are kept.
pub fn apply_frequency_filters(
    counts: &mut Frame,
    locus: &str,
    expected: &Thresholds,
    deletion: &Thresholds,
    min_allele_median_ratio: f64,
    min_gene_median_ratio: f64,
) -> Result<()> {
    fn outside(frame: &Frame, row: usize, locus: &str) -> Result<bool> {
        Ok(!frame
            .value(row, "db_name")?
            .julia_string()
            .starts_with(locus))
    }
    GermlineFilter::new()
        .with(CustomFilter::new(
            "Allele frequency threshold",
            |frame: &Frame, row| {
                Ok(float_at(frame, row, "allele_freq")?
                    >= threshold_for(expected, frame, row, "allele_freq")?)
            },
        ))
        .with(CustomFilter::new(
            "Gene case frequency threshold",
            |frame: &Frame, row| {
                Ok(float_at(frame, row, "gene_case_freq")?
                    >= threshold_for(deletion, frame, row, "gene_case_freq")?)
            },
        ))
        .with(CustomFilter::new(
            "Min allele median ratio",
            |frame: &Frame, row| {
                Ok(outside(frame, row, locus)?
                    || float_at(frame, row, "allele_to_cross_case_median_ratio")?
                        >= min_allele_median_ratio)
            },
        ))
        .with(CustomFilter::new(
            "Min gene median ratio",
            |frame: &Frame, row| {
                Ok(outside(frame, row, locus)?
                    || float_at(frame, row, "gene_to_cross_case_median_ratio")?
                        >= min_gene_median_ratio)
            },
        ))
        .apply(counts)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Donor D1: A*01 with two flank variants (count 30 on both rows),
    /// A*02 (10), B*01 (20) and the control CTRL*01 (50). Donor D2: A*01 (8).
    fn counts() -> Frame {
        Frame::from_columns([
            ("well", Column::from_ints([1, 1, 1, 1, 1, 2])),
            (
                "case",
                Column::from_strings(["D1", "D1", "D1", "D1", "D1", "D2"]),
            ),
            (
                "db_name",
                Column::from_strings(["A*01", "A*01", "A*02", "B*01", "CTRL*01", "A*01"]),
            ),
            (
                "sequence",
                Column::from_strings(["AAAA", "AAAA", "AAAT", "CCCC", "GGGG", "AAAA"]),
            ),
            (
                "gene",
                Column::from_strings(["A", "A", "A", "B", "CTRL", "A"]),
            ),
            ("count", Column::from_ints([30, 30, 10, 20, 50, 8])),
        ])
        .unwrap()
    }

    #[test]
    fn each_allele_is_counted_once_whatever_its_number_of_flank_variants() {
        let mut frame = counts();
        add_frequencies(&mut frame, "").unwrap();
        assert_eq!(
            frame.ints("gene_count").unwrap(),
            [40, 40, 40, 20, 50, 8].map(Some)
        );
        assert_eq!(
            frame.ints("case_count").unwrap(),
            [110, 110, 110, 110, 110, 8].map(Some)
        );
        let allele_freq: Vec<f64> = frame
            .floats("allele_freq")
            .unwrap()
            .iter()
            .map(|v| v.unwrap())
            .collect();
        assert_eq!(allele_freq, [0.75, 0.75, 0.25, 1.0, 1.0, 1.0]);
    }

    #[test]
    fn control_genes_stay_out_of_the_totals_but_in_the_table() {
        let mut frame = counts();
        add_frequencies(&mut frame, "A").unwrap();
        // Only locus alleles enter the donor total, each once: 30 + 10.
        assert_eq!(frame.ints("case_count").unwrap()[0], Some(40));
        let none = Thresholds::new();
        apply_frequency_filters(&mut frame, "A", &none, &none, 0.1, 0.1).unwrap();
        let kept = frame.strs("db_name").unwrap();
        let kept: Vec<&str> = kept.iter().map(|n| n.as_deref().unwrap()).collect();
        // B and CTRL are outside the locus "A": the median filters do not
        // apply to them, so the control gene is still there for --refgene.
        assert_eq!(kept, ["A*01", "A*01", "A*02", "B*01", "CTRL*01", "A*01"]);
    }
}
