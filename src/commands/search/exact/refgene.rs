//! Ratios of allele counts relative to a reference gene (`--refgene`).
//!
//! A reference is named either as a gene (`CTRL9`) or as one allele
//! (`CTRL9*01`). Names are compared exactly: the gene `IGHV3-2` does not
//! match alleles of `IGHV3-23`.
//!
//! For every donor two ratios are added:
//!
//! * `count_<gene>_ratio`: allele count over the reference count (the
//!   reference allele's count, or the reference gene's total when a gene
//!   was named);
//! * `ref_gene_count_<gene>_ratio`: gene total over the reference gene's
//!   total.

use std::collections::HashSet;

use crate::error::Result;
use crate::frame::{Column, Frame};
use crate::seq::gene_of;
use crate::{info, warn};

fn ints(frame: &Frame, name: &str) -> Result<Vec<i64>> {
    Ok(frame.ints(name)?.iter().map(|c| c.unwrap_or(0)).collect())
}

fn text(frame: &Frame, name: &str) -> Result<Vec<String>> {
    Ok(frame
        .strs(name)?
        .iter()
        .map(|c| c.clone().unwrap_or_default())
        .collect())
}

/// Does `db_name` belong to the reference named `reference`?
pub fn is_reference(reference: &str, db_name: &str) -> bool {
    if reference.contains('*') {
        db_name == reference
    } else {
        gene_of(db_name) == reference
    }
}

/// Sum of `count` over the distinct alleles among `rows` that satisfy
/// `select`. `count` repeats on every flank variant of an allele, so each
/// allele name is counted once.
fn allele_total(
    rows: &[usize],
    db_name: &[String],
    count: &[i64],
    select: impl Fn(&str) -> bool,
) -> i64 {
    let mut seen = HashSet::new();
    rows.iter()
        .filter(|&&r| select(&db_name[r]) && seen.insert(db_name[r].as_str()))
        .map(|&r| count[r])
        .sum()
}

/// Add `ref_gene_count` and both reference ratios for `refgene`.
/// Rows come back grouped by (well, case) in order of first appearance.
///
/// In a donor without the reference the ratios are undefined: they are
/// left empty (missing) and a warning is logged.
pub fn add_reference_ratios(counts: &Frame, refgene: &str) -> Result<Frame> {
    let reference_gene = gene_of(refgene).to_string();
    let db_name = text(counts, "db_name")?;
    let gene = text(counts, "gene")?;
    let count = ints(counts, "count")?;

    let mut order = Vec::with_capacity(counts.nrows());
    let mut gene_totals = Vec::with_capacity(counts.nrows());
    let mut allele_ratio = Vec::with_capacity(counts.nrows());
    let mut gene_ratio = Vec::with_capacity(counts.nrows());
    for rows in counts.group_by(&["well", "case"])?.iter() {
        let well = counts.value(rows[0], "well")?.julia_string();
        let case = counts.value(rows[0], "case")?.julia_string();
        let reference_count =
            allele_total(rows, &db_name, &count, |name| is_reference(refgene, name));
        let reference_gene_count = allele_total(rows, &db_name, &count, |name| {
            gene_of(name) == reference_gene
        });
        if reference_count == 0 {
            warn!("Reference name {refgene} not found in well {well} and case {case}");
        } else {
            info!("Applying name {refgene} to well {well} and case {case}");
        }
        let relative = |value: i64, reference: i64| -> Option<f64> {
            (reference > 0).then(|| value as f64 / reference as f64)
        };
        for &row in rows {
            let gene_total =
                allele_total(rows, &db_name, &count, |name| gene_of(name) == gene[row]);
            order.push(row);
            gene_totals.push(gene_total);
            allele_ratio.push(relative(count[row], reference_count));
            gene_ratio.push(relative(gene_total, reference_gene_count));
        }
    }

    let mut out = counts.take(&order);
    out.set_column("ref_gene_count", Column::from_ints(gene_totals))?;
    // Column order: allele ratio, gene total, gene ratio.
    let allele_column = format!("count_{reference_gene}_ratio");
    let gene_column = format!("ref_gene_count_{reference_gene}_ratio");
    out.set_column(&allele_column, Column::Float(allele_ratio))?;
    out.set_column(&gene_column, Column::Float(gene_ratio))?;
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counts() -> Frame {
        // Donor D1: CTRL9*01 (two flank variants), CTRL9*02, IGHV3-2*01, IGHV3-23*01.
        Frame::from_columns([
            ("well", Column::from_ints([1, 1, 1, 1, 1, 2])),
            (
                "case",
                Column::from_strings(["D1", "D1", "D1", "D1", "D1", "D2"]),
            ),
            (
                "db_name",
                Column::from_strings([
                    "CTRL9*01",
                    "CTRL9*01",
                    "CTRL9*02",
                    "IGHV3-2*01",
                    "IGHV3-23*01",
                    "IGHV3-23*01",
                ]),
            ),
            (
                "gene",
                Column::from_strings([
                    "CTRL9", "CTRL9", "CTRL9", "IGHV3-2", "IGHV3-23", "IGHV3-23",
                ]),
            ),
            ("count", Column::from_ints([40, 40, 10, 20, 100, 7])),
        ])
        .unwrap()
    }

    /// Ratios with `NAN` standing for a missing value.
    fn floats(frame: &Frame, column: &str) -> Vec<f64> {
        frame
            .floats(column)
            .unwrap()
            .iter()
            .map(|v| v.unwrap_or(f64::NAN))
            .collect()
    }

    #[test]
    fn names_are_matched_exactly() {
        assert!(is_reference("IGHV3-2", "IGHV3-2*01"));
        assert!(!is_reference("IGHV3-2", "IGHV3-23*01"));
        assert!(is_reference("CTRL9*01", "CTRL9*01"));
        assert!(!is_reference("CTRL9*01", "CTRL9*011"));
        assert!(!is_reference("CTRL9*01", "CTRL9*02"));
    }

    #[test]
    fn gene_reference_uses_the_gene_total_once_per_allele() {
        let out = add_reference_ratios(&counts(), "CTRL9").unwrap();
        // CTRL9 total = 40 + 10; the second flank row of *01 is not added again.
        assert_eq!(
            floats(&out, "count_CTRL9_ratio")[..5],
            [0.8, 0.8, 0.2, 0.4, 2.0]
        );
        assert_eq!(
            out.ints("ref_gene_count").unwrap(),
            [50, 50, 50, 20, 100, 7].map(Some)
        );
        assert_eq!(
            floats(&out, "ref_gene_count_CTRL9_ratio")[..5],
            [1.0, 1.0, 1.0, 0.4, 2.0]
        );
    }

    #[test]
    fn allele_reference_divides_alleles_by_the_allele_and_genes_by_the_gene() {
        let out = add_reference_ratios(&counts(), "CTRL9*01").unwrap();
        assert_eq!(
            floats(&out, "count_CTRL9_ratio")[..5],
            [1.0, 1.0, 0.25, 0.5, 2.5]
        );
        assert_eq!(
            floats(&out, "ref_gene_count_CTRL9_ratio")[..5],
            [1.0, 1.0, 1.0, 0.4, 2.0]
        );
    }

    #[test]
    fn donor_without_the_reference_has_no_ratio() {
        // D2 (last row) carries no CTRL9 allele.
        let out = add_reference_ratios(&counts(), "CTRL9").unwrap();
        assert_eq!(out.floats("count_CTRL9_ratio").unwrap()[5], None);
        assert_eq!(out.floats("ref_gene_count_CTRL9_ratio").unwrap()[5], None);
        // The gene total itself is still reported.
        assert_eq!(out.ints("ref_gene_count").unwrap()[5], Some(7));
    }
}
