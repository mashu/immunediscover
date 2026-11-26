//! Multiple-testing corrections.

/// Benjamini-Hochberg adjusted p-values (q-values).
pub fn benjamini_hochberg(pvalues: &[f64]) -> Vec<f64> {
    let m = pvalues.len();
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&a, &b| pvalues[a].total_cmp(&pvalues[b]));
    let mut adjusted = vec![0.0; m];
    let mut running_min = 1.0f64;
    for rank in (1..=m).rev() {
        let index = order[rank - 1];
        let q = ((pvalues[index] * m as f64) / rank as f64).min(1.0);
        running_min = running_min.min(q);
        adjusted[index] = running_min;
    }
    adjusted
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adjusts_and_keeps_input_order() {
        let q = benjamini_hochberg(&[0.01, 0.04, 0.03, 0.5]);
        assert_eq!(q, vec![0.04, 0.05333333333333334, 0.05333333333333334, 0.5]);
        assert!(benjamini_hochberg(&[]).is_empty());
    }
}
