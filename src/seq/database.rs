//! Allele databases loaded from FASTA.

use std::collections::HashSet;
use std::path::Path;
use std::sync::OnceLock;

use fancy_regex::Regex;

use crate::error::{Error, Result};
use crate::io::fasta::read_fasta;

/// One database entry. `name` is the complete FASTA header line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Allele {
    pub name: String,
    pub sequence: String,
}

impl Allele {
    pub fn new(name: impl Into<String>, sequence: impl Into<String>) -> Self {
        Allele {
            name: name.into(),
            sequence: sequence.into(),
        }
    }
}

/// Load a database without consistency checks.
pub fn load_database(path: &Path) -> Result<Vec<Allele>> {
    Ok(read_fasta(path)?
        .into_iter()
        .map(|r| Allele::new(r.description, r.sequence))
        .collect())
}

/// Load a database and require well-formed, unique names and sequences.
pub fn load_validated_database(path: &Path) -> Result<Vec<Allele>> {
    let alleles = load_database(path)?;
    let mut names = HashSet::new();
    let mut sequences = HashSet::new();
    for allele in &alleles {
        if allele.sequence.is_empty() {
            return Err(Error::invalid(format!(
                "Empty sequence found: {}",
                allele.name
            )));
        }
        if !is_valid_identifier(&allele.name) {
            return Err(Error::invalid(format!(
                "Invalid identifier found: {}",
                allele.name
            )));
        }
        if !names.insert(allele.name.as_str()) {
            return Err(Error::invalid(format!(
                "Duplicate identifier found: {}",
                allele.name
            )));
        }
        if !sequences.insert(allele.sequence.as_str()) {
            return Err(Error::invalid(format!(
                "Duplicate sequence found: {}",
                allele.sequence
            )));
        }
    }
    Ok(alleles)
}

/// Allele naming rule: `GENE*ALLELE`, optionally with one `_S<hash>` suffix.
pub fn is_valid_identifier(id: &str) -> bool {
    static PATTERN: OnceLock<Regex> = OnceLock::new();
    let pattern = PATTERN.get_or_init(|| {
        Regex::new(r"[A-Z].*\d.*\*(.*S.*\d.*|.*\d+$)").expect("identifier pattern is valid")
    });
    let underscores = id.matches('_').count();
    let stars = id.matches('*').count();
    pattern.is_match(id).unwrap_or(false) && underscores <= 1 && stars <= 1
}

/// Whether a sequence contains only unambiguous bases.
pub fn is_unambiguous_dna(sequence: &str) -> bool {
    sequence
        .bytes()
        .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identifier_rule() {
        assert!(is_valid_identifier("IGHV1-2*02"));
        assert!(is_valid_identifier("IGHV1-2*02_S1234"));
        assert!(!is_valid_identifier("IGHV1-2"));
        assert!(!is_valid_identifier("IGHV1-2*02_S1_S2"));
        assert!(!is_valid_identifier("ighv1-2*02x"));
    }
}
