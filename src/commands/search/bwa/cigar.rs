//! SAM CIGAR handling for the BWA genome check.

/// A parsed SAM CIGAR string.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Cigar {
    operations: Vec<(usize, char)>,
}

impl Cigar {
    /// Parse SAM text (`*` and malformed input give an empty CIGAR).
    pub fn parse(text: &str) -> Self {
        let mut operations = Vec::new();
        let mut length = String::new();
        for ch in text.chars() {
            if ch.is_ascii_digit() {
                length.push(ch);
            } else {
                if let Ok(n) = length.parse() {
                    operations.push((n, ch));
                }
                length.clear();
            }
        }
        Cigar { operations }
    }

    pub fn is_empty(&self) -> bool {
        self.operations.is_empty()
    }

    fn clip(operation: Option<&(usize, char)>) -> usize {
        match operation {
            Some(&(length, 'S' | 'H')) => length,
            _ => 0,
        }
    }

    /// Query bases clipped before the first aligned base (reference order).
    pub fn leading_clip(&self) -> usize {
        Self::clip(self.operations.first())
    }

    /// Query bases clipped after the last aligned base (reference order).
    pub fn trailing_clip(&self) -> usize {
        if self.operations.len() < 2 {
            return 0;
        }
        Self::clip(self.operations.last())
    }

    /// Reference bases covered by the alignment (`M`, `D`, `N`, `=`, `X`).
    /// Clips and insertions consume no reference.
    pub fn reference_span(&self) -> usize {
        self.operations
            .iter()
            .filter(|(_, op)| matches!(op, 'M' | 'D' | 'N' | 'X' | '='))
            .map(|(length, _)| length)
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clips_consume_no_reference() {
        let cigar = Cigar::parse("5S90M2D3I5M7S");
        assert_eq!(cigar.leading_clip(), 5);
        assert_eq!(cigar.trailing_clip(), 7);
        assert_eq!(cigar.reference_span(), 97);
    }

    #[test]
    fn unclipped_and_empty() {
        let cigar = Cigar::parse("100M");
        assert_eq!((cigar.leading_clip(), cigar.trailing_clip()), (0, 0));
        assert_eq!(cigar.reference_span(), 100);
        assert!(Cigar::parse("*").is_empty());
        assert_eq!(Cigar::parse("*").reference_span(), 0);
        // A lone clip is not both leading and trailing.
        assert_eq!(Cigar::parse("10S").trailing_clip(), 0);
    }
}
