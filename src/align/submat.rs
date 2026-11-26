//! Nucleotide substitution scores.

/// IUPAC nucleotide codes in the row/column order of the NUC.4.4 matrix.
const SYMBOLS: &[u8; 15] = b"ATGCSWRYKMBVHDN";

/// EDNAFULL / NUC.4.4 substitution matrix.
#[rustfmt::skip]
const EDNAFULL: [[i32; 15]; 15] = [
    //  A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N
    [   5, -4, -4, -4, -4,  1,  1, -4, -4,  1, -4, -1, -1, -1, -2], // A
    [  -4,  5, -4, -4, -4,  1, -4,  1,  1, -4, -1, -4, -1, -1, -2], // T
    [  -4, -4,  5, -4,  1, -4,  1, -4,  1, -4, -1, -1, -4, -1, -2], // G
    [  -4, -4, -4,  5,  1, -4, -4,  1, -4,  1, -1, -1, -1, -4, -2], // C
    [  -4, -4,  1,  1, -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1], // S
    [   1,  1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1], // W
    [   1, -4,  1, -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1], // R
    [  -4,  1, -4,  1, -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1], // Y
    [  -4,  1,  1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1], // K
    [   1, -4, -4,  1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1], // M
    [  -4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1], // B
    [  -1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1], // V
    [  -1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1], // H
    [  -1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1], // D
    [  -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1], // N
];

/// A nucleotide sequence encoded as indices into the substitution matrix.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EncodedDna(Vec<u8>);

impl EncodedDna {
    /// Encode IUPAC text (case-insensitive); `None` on any other character.
    pub fn new(sequence: &str) -> Option<Self> {
        sequence
            .bytes()
            .map(|b| {
                let upper = b.to_ascii_uppercase();
                SYMBOLS.iter().position(|&s| s == upper).map(|i| i as u8)
            })
            .collect::<Option<Vec<u8>>>()
            .map(EncodedDna)
    }

    /// Wrap codes previously obtained from an [`EncodedDna`].
    pub fn from_codes(codes: Vec<u8>) -> Self {
        debug_assert!(codes.iter().all(|&c| (c as usize) < SYMBOLS.len()));
        EncodedDna(codes)
    }

    pub fn codes(&self) -> &[u8] {
        &self.0
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Upper-case text of a slice of codes.
    pub fn decode(codes: &[u8]) -> String {
        codes.iter().map(|&c| SYMBOLS[c as usize] as char).collect()
    }
}

/// EDNAFULL score of two encoded nucleotides.
#[inline]
pub fn ednafull(a: u8, b: u8) -> i32 {
    EDNAFULL[a as usize][b as usize]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matrix_is_symmetric_with_expected_diagonal() {
        for (i, row) in EDNAFULL.iter().enumerate() {
            for (j, score) in row.iter().enumerate() {
                assert_eq!(*score, EDNAFULL[j][i]);
            }
        }
        let acgt = EncodedDna::new("acgt").unwrap();
        for &c in acgt.codes() {
            assert_eq!(ednafull(c, c), 5);
        }
        assert_eq!(EncodedDna::decode(acgt.codes()), "ACGT");
        assert!(EncodedDna::new("AC-T").is_none());
    }
}
