//! Content-derived allele names (`_Sxxxx` suffixes).

use md5::{Digest, Md5};

/// Number of decimal digits in a sequence hash suffix.
pub const HASH_DIGITS: usize = 4;

/// `S` followed by a zero-padded decimal number derived from the last
/// `digits` hex characters of the sequence's MD5.
pub fn sequence_hash(sequence: &str, digits: usize) -> String {
    let digest = Md5::digest(sequence.as_bytes());
    let hex: String = digest.iter().map(|b| format!("{b:02x}")).collect();
    let tail = &hex[hex.len() - digits..];
    let value = u64::from_str_radix(tail, 16).expect("md5 hex digits are valid hex");
    let reduced = value % 10u64.pow(digits as u32);
    format!("S{reduced:0digits$}")
}

/// Name for a novel allele: the base name up to its first underscore,
/// followed by `_` and the sequence hash.
pub fn unique_name(name: &str, sequence: &str) -> String {
    let base = name.split('_').next().unwrap_or(name);
    format!("{base}_{}", sequence_hash(sequence, HASH_DIGITS))
}

#[cfg(test)]
mod tests {
    use super::*;

    // Reference values produced by immunediscover.jl `Data.sequence_hash`.
    #[test]
    fn hash_matches_julia() {
        assert_eq!(sequence_hash("ACGT", 4), "S1086");
        assert_eq!(sequence_hash("", 4), "S7022");
        assert_eq!(sequence_hash("GATTACA", 4), "S5291");
    }

    #[test]
    fn unique_name_replaces_existing_suffix() {
        assert_eq!(unique_name("IGHV1-2*02_S1234", "ACGT"), "IGHV1-2*02_S1086");
        assert_eq!(unique_name("A_B_C", "ACGT"), "A_S1086");
    }
}
