//! Gene segment types and allele-name conventions.

use std::fmt;
use std::str::FromStr;

use crate::error::Error;

/// Immunoglobulin / T-cell receptor gene segment type.
///
/// The segment type decides where recombination signal sequences sit
/// relative to the coding sequence, which every flank-related routine
/// dispatches on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GeneType {
    V,
    D,
    J,
}

impl GeneType {
    pub const ALL: [GeneType; 3] = [GeneType::V, GeneType::D, GeneType::J];

    pub fn as_str(self) -> &'static str {
        match self {
            GeneType::V => "V",
            GeneType::D => "D",
            GeneType::J => "J",
        }
    }

    /// Infer the segment type from an allele or gene name.
    ///
    /// The first of `V`, `D`, `J` (in that priority) occurring anywhere in
    /// the name wins, e.g. `IGHV1-2*02` is V and `IGHD3-10` is D.
    pub fn from_name(name: &str) -> Option<GeneType> {
        [('V', GeneType::V), ('D', GeneType::D), ('J', GeneType::J)]
            .into_iter()
            .find_map(|(letter, gene_type)| name.contains(letter).then_some(gene_type))
    }
}

impl fmt::Display for GeneType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl FromStr for GeneType {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "V" => Ok(GeneType::V),
            "D" => Ok(GeneType::D),
            "J" => Ok(GeneType::J),
            other => Err(Error::invalid(format!("Invalid gene type: {other}"))),
        }
    }
}

/// Gene part of an allele name: everything before the first `*`.
pub fn gene_of(allele: &str) -> &str {
    allele.split('*').next().unwrap_or(allele)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gene_is_prefix_before_star() {
        assert_eq!(gene_of("IGHV1-2*02"), "IGHV1-2");
        assert_eq!(gene_of("IGHV1-2*02_S1234"), "IGHV1-2");
        assert_eq!(gene_of("noallele"), "noallele");
    }

    #[test]
    fn type_from_name_prefers_v_then_d_then_j() {
        assert_eq!(GeneType::from_name("IGHV1-2"), Some(GeneType::V));
        assert_eq!(GeneType::from_name("IGHD3-10"), Some(GeneType::D));
        assert_eq!(GeneType::from_name("IGHJ4"), Some(GeneType::J));
        assert_eq!(GeneType::from_name("TRAC"), None);
    }
}
