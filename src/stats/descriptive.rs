//! Summary statistics with Julia `Statistics` conventions.

/// Arithmetic mean (`NaN` for an empty slice).
pub fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

/// Sample standard deviation with Bessel's correction (`NaN` for n < 2).
pub fn sample_std(values: &[f64]) -> f64 {
    let n = values.len();
    if n < 2 {
        return f64::NAN;
    }
    let m = mean(values);
    let sum_sq: f64 = values.iter().map(|v| (v - m) * (v - m)).sum();
    (sum_sq / (n - 1) as f64).sqrt()
}

/// Median; the mean of the two middle values for an even count.
pub fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 1 {
        sorted[mid]
    } else {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basics() {
        assert_eq!(mean(&[1.0, 2.0, 3.0]), 2.0);
        assert_eq!(median(&[1.0, 3.0, 5.0]), 3.0);
        assert_eq!(median(&[1.0, 2.0]), 1.5);
        assert!(
            (sample_std(&[2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]) - 2.138089935299395).abs()
                < 1e-15
        );
        assert!(sample_std(&[1.0]).is_nan());
    }
}
