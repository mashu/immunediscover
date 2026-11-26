//! Robust summaries: quantiles, spread, two-way median fit, rank statistics.

use super::descriptive::median;

/// The `q` quantile (0..=1) as an observed value (`NaN` when empty).
pub fn quantile(values: &[f64], q: f64) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let rank = (q * sorted.len() as f64).ceil() as usize;
    sorted[rank.clamp(1, sorted.len()) - 1]
}

/// Median absolute deviation, scaled to estimate a normal standard deviation.
pub fn mad(values: &[f64]) -> f64 {
    let centre = median(values);
    let deviations: Vec<f64> = values.iter().map(|v| (v - centre).abs()).collect();
    1.4826 * median(&deviations)
}

/// Additive fit `x[r][c] ≈ row[r] + column[c]` by median polish.
#[derive(Debug, Clone)]
pub struct TwoWayFit {
    pub row: Vec<f64>,
    pub column: Vec<f64>,
}

/// Median polish of a table with missing cells. Rows or columns without any
/// cell get `NaN`.
pub fn median_polish(cells: &[Vec<Option<f64>>], n_columns: usize, iterations: usize) -> TwoWayFit {
    let column_of = |c: usize, row: &[f64]| -> f64 {
        let values: Vec<f64> = cells
            .iter()
            .zip(row)
            .filter_map(|(cells, r)| cells[c].map(|x| x - r))
            .filter(|v| !v.is_nan())
            .collect();
        median(&values)
    };
    let mut row = vec![0.0; cells.len()];
    let mut column: Vec<f64> = (0..n_columns).map(|c| column_of(c, &row)).collect();
    for _ in 0..iterations {
        for (r, cells) in cells.iter().enumerate() {
            let values: Vec<f64> = cells
                .iter()
                .zip(&column)
                .filter_map(|(cell, c)| cell.map(|x| x - c))
                .filter(|v| !v.is_nan())
                .collect();
            row[r] = median(&values);
        }
        for (c, value) in column.iter_mut().enumerate() {
            *value = column_of(c, &row);
        }
    }
    TwoWayFit { row, column }
}

/// Area under the ROC curve: the probability that a positive scores above a
/// negative, ties counting half (`NaN` if either group is empty).
pub fn auc(positive: &[f64], negative: &[f64]) -> f64 {
    if positive.is_empty() || negative.is_empty() {
        return f64::NAN;
    }
    let mut all: Vec<(f64, bool)> = positive
        .iter()
        .map(|&v| (v, true))
        .chain(negative.iter().map(|&v| (v, false)))
        .collect();
    all.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut rank_sum = 0.0;
    let mut i = 0;
    while i < all.len() {
        let mut j = i;
        while j + 1 < all.len() && all[j + 1].0 == all[i].0 {
            j += 1;
        }
        let mean_rank = (i + j) as f64 / 2.0 + 1.0;
        rank_sum += mean_rank * all[i..=j].iter().filter(|x| x.1).count() as f64;
        i = j + 1;
    }
    let (np, nn) = (positive.len() as f64, negative.len() as f64);
    (rank_sum - np * (np + 1.0) / 2.0) / (np * nn)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantile_and_mad() {
        let v = [5.0, 1.0, 3.0, 2.0, 4.0];
        assert_eq!(quantile(&v, 0.5), 3.0);
        assert_eq!(quantile(&v, 0.05), 1.0);
        assert_eq!(quantile(&v, 1.0), 5.0);
        assert!((mad(&v) - 1.4826).abs() < 1e-12);
        assert!(quantile(&[], 0.5).is_nan());
    }

    #[test]
    fn median_polish_recovers_an_additive_table() {
        let (rows, columns) = ([0.0, 1.0, -2.0, 0.5], [10.0, 11.0, 13.0]);
        let mut cells: Vec<Vec<Option<f64>>> = rows
            .iter()
            .map(|r| columns.iter().map(|c| Some(r + c)).collect())
            .collect();
        cells[1][2] = None; // a missing cell
        cells[2][0] = Some(50.0); // an outlier
        let fit = median_polish(&cells, 3, 30);
        for (r, row) in rows.iter().enumerate() {
            for (c, column) in columns.iter().enumerate() {
                if (r, c) != (2, 0) {
                    assert!((fit.row[r] + fit.column[c] - (row + column)).abs() < 1e-9);
                }
            }
        }
    }

    #[test]
    fn auc_counts_ties_as_half() {
        assert_eq!(auc(&[3.0, 4.0], &[1.0, 2.0]), 1.0);
        assert_eq!(auc(&[1.0, 2.0], &[3.0, 4.0]), 0.0);
        assert_eq!(auc(&[1.0], &[1.0]), 0.5);
        assert_eq!(auc(&[2.0, 1.0], &[1.0]), 0.75);
        assert!(auc(&[], &[1.0]).is_nan());
    }
}
