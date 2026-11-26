//! Read-border handling for extension mode.
//!
//! The outermost `border` bases of a read (typically primer or index
//! sequence) must not end up inside an extracted extension. Matches whose
//! extension would reach into a border are rejected, optionally after the
//! extension length has been reduced per gene (`--adjust-per-gene-extension`).

use std::collections::HashMap;

use crate::seq::GeneType;

/// Per-gene extension lengths chosen by calibration.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PerGeneExtension {
    pub prefix: HashMap<String, usize>,
    pub suffix: HashMap<String, usize>,
}

impl PerGeneExtension {
    pub fn prefix_for(&self, gene: &str, default: usize) -> usize {
        self.prefix.get(gene).copied().unwrap_or(default)
    }

    pub fn suffix_for(&self, gene: &str, default: usize) -> usize {
        self.suffix.get(gene).copied().unwrap_or(default)
    }

    pub fn n_genes(&self) -> usize {
        let mut genes: Vec<&String> = self.prefix.keys().chain(self.suffix.keys()).collect();
        genes.sort();
        genes.dedup();
        genes.len()
    }
}

/// Geometry of one match inside its read (1-based, inclusive).
#[derive(Debug, Clone, Copy)]
pub struct MatchSpan {
    pub start: i64,
    pub end: i64,
    pub read_length: i64,
}

impl MatchSpan {
    pub fn new(start: usize, end: usize, read_length: usize) -> Self {
        MatchSpan {
            start: start as i64,
            end: end as i64,
            read_length: read_length as i64,
        }
    }

    /// Last position of the left border.
    fn left_border_end(&self, border: i64) -> i64 {
        border.min(self.read_length)
    }

    /// First position of the right border.
    fn right_border_start(&self, border: i64) -> i64 {
        (self.read_length - border + 1).max(1)
    }

    fn upstream_reaches_border(&self, extension: i64, border: i64) -> bool {
        (self.start - extension).max(1) <= self.left_border_end(border)
    }

    fn downstream_reaches_border(&self, extension: i64, border: i64) -> bool {
        (self.end + extension).min(self.read_length) >= self.right_border_start(border)
    }
}

/// Whether a fixed-length extension overlaps a read border.
///
/// Unlike [`rejects_with_calibration`], an extension that is empty because
/// the match touches the read end does not count as an overlap.
pub fn extension_overlaps_border(
    span: MatchSpan,
    gene: GeneType,
    extension: usize,
    border: usize,
) -> bool {
    let (ext, border) = (extension as i64, border as i64);
    if border <= 0 || ext <= 0 {
        return false;
    }
    let upstream = || {
        let first = (span.start - ext).max(1);
        first < span.start && first <= span.left_border_end(border)
    };
    let downstream = || {
        let last = (span.end + ext).min(span.read_length);
        span.end < last && last >= span.right_border_start(border)
    };
    match gene {
        GeneType::V => downstream(),
        GeneType::J => upstream(),
        GeneType::D => upstream() || downstream(),
    }
}

/// Border test used with calibrated per-gene extension lengths.
pub fn rejects_with_calibration(
    span: MatchSpan,
    gene: GeneType,
    calibrated: &PerGeneExtension,
    gene_name: &str,
    extension: usize,
    border: usize,
) -> bool {
    let border = border as i64;
    let prefix = calibrated.prefix_for(gene_name, extension) as i64;
    let suffix = calibrated.suffix_for(gene_name, extension) as i64;
    match gene {
        GeneType::V => span.downstream_reaches_border(suffix, border),
        GeneType::J => span.upstream_reaches_border(prefix, border),
        GeneType::D => {
            span.upstream_reaches_border(prefix, border)
                || span.downstream_reaches_border(suffix, border)
        }
    }
}

/// Longest border-safe extension lengths for one match: `(upstream, downstream)`.
/// A side that is irrelevant for the gene type is `None`.
pub fn allowed_lengths(
    span: MatchSpan,
    gene: GeneType,
    extension: usize,
    border: usize,
) -> (Option<usize>, Option<usize>) {
    let (ext, border) = (extension as i64, border as i64);
    let upstream = || ext.min((span.start - 1 - span.left_border_end(border)).max(0)) as usize;
    let downstream = || ext.min((span.right_border_start(border) - 1 - span.end).max(0)) as usize;
    match gene {
        GeneType::V => (None, Some(downstream())),
        GeneType::J => (Some(upstream()), None),
        GeneType::D => (Some(upstream()), Some(downstream())),
    }
}

/// The extension length that keeps at least a fraction `safe_fraction` of
/// the observed matches clear of the border.
pub fn safe_quantile(values: &[usize], safe_fraction: f64, fallback: usize) -> usize {
    if values.is_empty() {
        return fallback;
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    let rank = ((1.0 - safe_fraction) * sorted.len() as f64).ceil() as usize;
    sorted[rank.clamp(1, sorted.len()) - 1]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn v_extension_into_right_border_is_rejected() {
        // Read of 100, border 10 -> right border starts at 91.
        let near = MatchSpan::new(40, 85, 100);
        assert!(extension_overlaps_border(near, GeneType::V, 10, 10));
        assert!(!extension_overlaps_border(near, GeneType::V, 5, 10));
        assert!(!extension_overlaps_border(near, GeneType::V, 10, 0));
        assert!(!extension_overlaps_border(near, GeneType::J, 10, 10));
    }

    #[test]
    fn j_and_d_check_the_upstream_side() {
        let near_start = MatchSpan::new(15, 60, 100);
        assert!(extension_overlaps_border(near_start, GeneType::J, 6, 10));
        assert!(!extension_overlaps_border(near_start, GeneType::J, 4, 10));
        assert!(extension_overlaps_border(near_start, GeneType::D, 6, 10));
    }

    #[test]
    fn allowed_lengths_are_capped_by_the_border() {
        let span = MatchSpan::new(15, 85, 100);
        assert_eq!(allowed_lengths(span, GeneType::V, 20, 10), (None, Some(5)));
        assert_eq!(allowed_lengths(span, GeneType::J, 20, 10), (Some(4), None));
        assert_eq!(
            allowed_lengths(span, GeneType::D, 3, 10),
            (Some(3), Some(3))
        );
    }

    #[test]
    fn quantile_targets_safe_fraction() {
        let values = [10, 2, 8, 6, 4, 10, 10, 10, 10, 10];
        assert_eq!(safe_quantile(&values, 1.0, 99), 2);
        assert_eq!(safe_quantile(&values, 0.9, 99), 2);
        assert_eq!(safe_quantile(&values, 0.8, 99), 4);
        assert_eq!(safe_quantile(&[], 0.9, 99), 99);
    }
}
