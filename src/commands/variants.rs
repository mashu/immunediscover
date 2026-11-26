//! Counting alleles and their flank variants: the conventions shared by
//! every command that reports read counts (`search exact`, `discover
//! blast`).
//!
//! A row is one *variant*: an allele sequence together with the read
//! context it was seen in (recombination signal, extension). Two counts
//! describe it:
//!
//! * `full_count`: reads showing exactly this variant;
//! * `count`: reads showing the allele, whatever the context, i.e. the sum
//!   of `full_count` over the variants of the allele. It repeats on every
//!   variant row of the allele.
//!
//! `full_ratio` and `ratio` are those counts relative to the largest value
//! among the rows of the same gene and donor. Thresholds are applied to
//! both levels, variants are ranked by `full_count` within their allele
//! (`flank_index`) and the best `top` are reported.

use crate::error::Result;
use crate::frame::{Column, Frame, SortKey};

/// The columns that identify donor, gene and allele in a table of variants.
#[derive(Debug, Clone, Copy)]
pub struct VariantColumns<'a> {
    pub well: &'a str,
    pub case: &'a str,
    pub gene: &'a str,
    /// Allele (or nearest database allele) name.
    pub name: &'a str,
    /// Allele sequence without read context.
    pub sequence: &'a str,
}

impl<'a> VariantColumns<'a> {
    fn donor_gene(&self) -> [&'a str; 3] {
        [self.well, self.case, self.gene]
    }

    fn allele(&self) -> [&'a str; 4] {
        [self.well, self.case, self.name, self.sequence]
    }

    fn sort_order(&self) -> [&'a str; 5] {
        [self.well, self.case, self.gene, self.name, self.sequence]
    }
}

/// `count` for every row: `full_count` summed over the allele's variants.
pub fn allele_counts(
    frame: &Frame,
    columns: &VariantColumns<'_>,
    full_count: &[i64],
) -> Result<Vec<i64>> {
    Ok(frame
        .group_by(&columns.allele())?
        .broadcast(|rows| rows.iter().map(|&r| full_count[r]).sum::<i64>()))
}

/// `values ./ maximum(values)` within each donor and gene.
pub fn ratio_to_gene_max(
    frame: &Frame,
    columns: &VariantColumns<'_>,
    values: &[i64],
) -> Result<Vec<f64>> {
    Ok(frame.group_by(&columns.donor_gene())?.map_rows(|rows| {
        let max = rows.iter().map(|&r| values[r]).max().unwrap_or(1) as f64;
        rows.iter().map(|&r| values[r] as f64 / max).collect()
    }))
}

/// Order rows by donor, gene, allele; within an allele the best-supported
/// variant comes first.
pub fn sort_variants(frame: &Frame, columns: &VariantColumns<'_>) -> Result<Frame> {
    frame
        .sorted(&[SortKey::desc("full_count"), SortKey::desc("count")])?
        .sorted_asc(&columns.sort_order())
}

/// Number the variants of every allele in row order (`flank_index`, from 1)
/// and keep the first `top`. Rows must already be in [`sort_variants`]
/// order for the numbering to mean "by support".
pub fn rank_variants(frame: &Frame, columns: &VariantColumns<'_>, top: usize) -> Result<Frame> {
    let flank_index: Vec<i64> = frame
        .group_by(&columns.sort_order())?
        .map_rows(|rows| (1..=rows.len() as i64).collect());
    let mask: Vec<bool> = flank_index.iter().map(|&i| i as usize <= top).collect();
    let mut ranked = frame.clone();
    ranked.set_column("flank_index", Column::from_ints(flank_index))?;
    Ok(ranked.filter(&mask))
}

/// One row per allele of a donor: `true` on the first of its variant rows.
///
/// Totals and medians of `count` must use each allele once, because
/// `count` repeats on every variant row.
pub fn first_variant_rows(frame: &Frame, columns: &VariantColumns<'_>) -> Result<Vec<bool>> {
    let mut first = vec![false; frame.nrows()];
    for row in frame.group_by(&columns.allele())?.first_rows() {
        first[row] = true;
    }
    Ok(first)
}

#[cfg(test)]
mod tests {
    use super::*;

    const COLUMNS: VariantColumns<'static> = VariantColumns {
        well: "well",
        case: "case",
        gene: "gene",
        name: "name",
        sequence: "sequence",
    };

    /// Donor D1: A*01 in three contexts (5, 9 and 2 reads), A*02 in one (4).
    fn variants() -> Frame {
        Frame::from_columns([
            ("well", Column::from_ints([1, 1, 1, 1])),
            ("case", Column::from_strings(["D1", "D1", "D1", "D1"])),
            ("gene", Column::from_strings(["A", "A", "A", "A"])),
            (
                "name",
                Column::from_strings(["A*01", "A*02", "A*01", "A*01"]),
            ),
            (
                "sequence",
                Column::from_strings(["ACGT", "ACGA", "ACGT", "ACGT"]),
            ),
            ("context", Column::from_strings(["x", "x", "y", "z"])),
            ("full_count", Column::from_ints([5, 4, 9, 2])),
        ])
        .unwrap()
    }

    #[test]
    fn count_sums_the_variants_of_an_allele() {
        let frame = variants();
        let count = allele_counts(&frame, &COLUMNS, &[5, 4, 9, 2]).unwrap();
        assert_eq!(count, [16, 4, 16, 16]);
        assert_eq!(
            ratio_to_gene_max(&frame, &COLUMNS, &count).unwrap(),
            [1.0, 0.25, 1.0, 1.0]
        );
        assert_eq!(
            first_variant_rows(&frame, &COLUMNS).unwrap(),
            [true, true, false, false]
        );
    }

    #[test]
    fn variants_are_ranked_by_support_within_their_allele() {
        let mut frame = variants();
        frame
            .push_column("count", Column::from_ints([16, 4, 16, 16]))
            .unwrap();
        let sorted = sort_variants(&frame, &COLUMNS).unwrap();
        let ranked = rank_variants(&sorted, &COLUMNS, 2).unwrap();
        let context = ranked.strs("context").unwrap();
        let context: Vec<&str> = context.iter().map(|c| c.as_deref().unwrap()).collect();
        // A*01: y (9 reads) then x (5); z (2) is beyond the top two. Then A*02.
        assert_eq!(context, ["y", "x", "x"]);
        assert_eq!(
            ranked.ints("flank_index").unwrap(),
            [Some(1), Some(2), Some(1)]
        );
        assert_eq!(rank_variants(&sorted, &COLUMNS, 1).unwrap().nrows(), 2);
    }
}
