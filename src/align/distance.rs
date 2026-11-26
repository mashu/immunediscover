//! Edit and Hamming distances on byte strings.

/// Levenshtein distance (unit cost substitutions, insertions, deletions).
pub fn levenshtein(a: &[u8], b: &[u8]) -> usize {
    let (short, long) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    let mut row: Vec<usize> = (0..=short.len()).collect();
    for (j, &lc) in long.iter().enumerate() {
        let mut diagonal = row[0];
        row[0] = j + 1;
        for (i, &sc) in short.iter().enumerate() {
            let substitution = diagonal + usize::from(sc != lc);
            diagonal = row[i + 1];
            row[i + 1] = substitution.min(row[i] + 1).min(diagonal + 1);
        }
    }
    row[short.len()]
}

/// Number of mismatching positions; `None` when the lengths differ.
pub fn hamming(a: &[u8], b: &[u8]) -> Option<usize> {
    (a.len() == b.len()).then(|| a.iter().zip(b).filter(|(x, y)| x != y).count())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn levenshtein_basics() {
        assert_eq!(levenshtein(b"kitten", b"sitting"), 3);
        assert_eq!(levenshtein(b"", b"abc"), 3);
        assert_eq!(levenshtein(b"ACGT", b"ACGT"), 0);
        assert_eq!(levenshtein(b"ACGT", b"AGT"), 1);
    }

    #[test]
    fn hamming_requires_equal_length() {
        assert_eq!(hamming(b"CACAGTG", b"CACAATG"), Some(1));
        assert_eq!(hamming(b"AC", b"A"), None);
    }
}
