//! Anchored index (barcode) patterns.
//!
//! Index sequences are regular expressions upstream. Almost all real
//! indices are plain nucleotide strings, so literals take a fast path and
//! only patterns with metacharacters go through the regex engine.

use fancy_regex::Regex;

use crate::error::{Error, Result};

const METACHARACTERS: &[char] = &[
    '\\', '.', '^', '$', '|', '?', '*', '+', '(', ')', '[', ']', '{', '}',
];

/// Compile a user-supplied pattern, reporting failures with the pattern text.
pub fn compile(pattern: &str) -> Result<Regex> {
    Regex::new(pattern).map_err(|e| Error::Regex {
        pattern: pattern.to_string(),
        message: e.to_string(),
    })
}

/// A pattern that must match at the very start or very end of a read.
pub enum Anchored {
    Literal(String),
    Regex(Regex),
}

impl Anchored {
    /// Pattern that has to match a prefix of the text.
    pub fn prefix(pattern: &str) -> Result<Self> {
        if pattern.contains(METACHARACTERS) {
            Ok(Anchored::Regex(compile(&format!("^(?:{pattern})"))?))
        } else {
            Ok(Anchored::Literal(pattern.to_string()))
        }
    }

    /// Pattern that has to match a suffix of the text.
    pub fn suffix(pattern: &str) -> Result<Self> {
        if pattern.contains(METACHARACTERS) {
            Ok(Anchored::Regex(compile(&format!("(?:{pattern})\\z"))?))
        } else {
            Ok(Anchored::Literal(pattern.to_string()))
        }
    }

    pub fn matches_start(&self, text: &str) -> bool {
        match self {
            Anchored::Literal(literal) => text.starts_with(literal.as_str()),
            Anchored::Regex(regex) => regex.is_match(text).unwrap_or(false),
        }
    }

    pub fn matches_end(&self, text: &str) -> bool {
        match self {
            Anchored::Literal(literal) => text.ends_with(literal.as_str()),
            Anchored::Regex(regex) => regex.is_match(text).unwrap_or(false),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn literal_and_regex_anchoring() {
        let literal = Anchored::prefix("ACG").unwrap();
        assert!(literal.matches_start("ACGTT") && !literal.matches_start("TACG"));
        let regex = Anchored::prefix("^.{2}AC.{2}GG").unwrap();
        assert!(regex.matches_start("TTACTTGGA") && !regex.matches_start("ATTACTTGG"));
        let tail = Anchored::suffix("G[AT]").unwrap();
        assert!(tail.matches_end("CCGA") && !tail.matches_end("GAC"));
    }
}
