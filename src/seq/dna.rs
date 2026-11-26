//! Nucleotide helpers: IUPAC validation and reverse complement.

/// Complement of an IUPAC nucleotide code, preserving case.
fn complement(base: u8) -> Option<u8> {
    let upper = match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'M' => b'K',
        b'R' => b'Y',
        b'W' => b'W',
        b'S' => b'S',
        b'Y' => b'R',
        b'K' => b'M',
        b'V' => b'B',
        b'H' => b'D',
        b'D' => b'H',
        b'B' => b'V',
        b'N' => b'N',
        _ => return None,
    };
    Some(upper)
}

/// Whether `sequence` consists only of IUPAC nucleotide codes (any case).
pub fn is_iupac_dna(sequence: &str) -> bool {
    sequence.bytes().all(|b| complement(b).is_some())
}

/// Upper-case reverse complement, or `None` if a non-IUPAC character occurs.
pub fn reverse_complement(sequence: &str) -> Option<String> {
    sequence
        .bytes()
        .rev()
        .map(|b| complement(b).map(char::from))
        .collect()
}

/// Index of an unambiguous base (A, C, G, T -> 0..=3), case-sensitive.
#[inline]
pub fn base_index(base: u8) -> Option<usize> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_complement_handles_ambiguity_codes() {
        assert_eq!(reverse_complement("ACGTN").as_deref(), Some("NACGT"));
        assert_eq!(reverse_complement("acgR").as_deref(), Some("YCGT"));
        assert_eq!(reverse_complement("AC-T"), None);
    }
}
