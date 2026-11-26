//! Position frequency profiles of equally long motifs.

use crate::seq::dna::base_index;

/// Smoothing constant added to every count.
const PSEUDOCOUNT: f64 = 0.0001;

/// A 4 x L matrix of per-position base values (rows A, C, G, T).
pub type Profile = Vec<[f64; 4]>;

/// Base counts per position; gaps and unknown characters are not counted.
pub fn counts<S: AsRef<str>>(motifs: &[S]) -> Vec<[u64; 4]> {
    let length = motifs.first().map_or(0, |m| m.as_ref().len());
    let mut table = vec![[0u64; 4]; length];
    for motif in motifs {
        for (position, &base) in motif.as_ref().as_bytes().iter().enumerate().take(length) {
            if let Some(index) = base_index(base) {
                table[position][index] += 1;
            }
        }
    }
    table
}

/// Smoothed base probabilities per position.
pub fn motif_profile<S: AsRef<str>>(motifs: &[S]) -> Profile {
    counts(motifs)
        .into_iter()
        .map(|column| {
            let smoothed = column.map(|c| c as f64 + PSEUDOCOUNT);
            let total: f64 = smoothed.iter().sum();
            smoothed.map(|c| c / total)
        })
        .collect()
}

/// Probability of `motif` under `profile`; a gap makes it impossible.
pub fn motif_prob(motif: &str, profile: &Profile) -> f64 {
    let mut log2_sum = 0.0;
    for (position, &base) in motif.as_bytes().iter().enumerate() {
        if base == b'-' {
            return 0.0;
        }
        if let Some(index) = base_index(base) {
            log2_sum += profile[position][index].log2();
        }
    }
    2f64.powf(log2_sum)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn profile_columns_sum_to_one() {
        let motifs = ["CACAGTG", "CACAATG", "CAC-GTG"];
        let table = counts(&motifs);
        assert_eq!(table[0], [0, 3, 0, 0]);
        assert_eq!(table[3], [2, 0, 0, 0]);
        let profile = motif_profile(&motifs);
        for column in &profile {
            assert!((column.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        }
        assert!(motif_prob("CACAGTG", &profile) > motif_prob("CACAATG", &profile));
        assert_eq!(motif_prob("CAC-GTG", &profile), 0.0);
    }
}
