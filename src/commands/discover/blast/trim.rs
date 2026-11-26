//! Removing database-extension affixes from BLAST query alignments.
//!
//! Database alleles are extended with the read context they are usually
//! found in, so BLAST alignments run past the allele ends. Each affix is
//! semi-globally aligned to the aligned read segment; if enough of the
//! affix matches, everything up to (prefix) or from (suffix) the affix is
//! cut off, leaving the candidate allele sequence.

use crate::align::{levenshtein, semiglobal, AlignedPair, EncodedDna};
use crate::error::{Error, Result};
use crate::frame::{format_float, round_digits};

const GAP_OPEN: i32 = -10;
const GAP_EXTEND: i32 = -1;
/// Shortest sequence accepted after trimming.
const MIN_CORE_LENGTH: usize = 10;

/// Which end of the read segment an affix belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AffixSide {
    Prefix,
    Suffix,
}

impl AffixSide {
    fn label(self) -> &'static str {
        match self {
            AffixSide::Prefix => "prefix",
            AffixSide::Suffix => "suffix",
        }
    }

    fn title(self) -> &'static str {
        match self {
            AffixSide::Prefix => "Prefix",
            AffixSide::Suffix => "Suffix",
        }
    }
}

/// Counters and messages describing trimming outcomes.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct AlignmentStats {
    pub total_attempts: usize,
    pub prefix_failures: usize,
    pub suffix_failures: usize,
    pub successful: usize,
    /// `(subject id, reason)` for every failed trim.
    pub failures: Vec<(String, String)>,
}

impl AlignmentStats {
    fn fail(&mut self, subject: &str, side: Option<AffixSide>, message: String) {
        self.failures.push((subject.to_string(), message));
        match side {
            Some(AffixSide::Prefix) => self.prefix_failures += 1,
            Some(AffixSide::Suffix) => self.suffix_failures += 1,
            None => {}
        }
    }

    /// Fold another worker's statistics into this one.
    pub fn merge(&mut self, other: AlignmentStats) {
        self.total_attempts += other.total_attempts;
        self.prefix_failures += other.prefix_failures;
        self.suffix_failures += other.suffix_failures;
        self.successful += other.successful;
        self.failures.extend(other.failures);
    }
}

/// The part of the query on the far side of the aligned affix.
fn core_columns(
    pairs: &[AlignedPair],
    side: AffixSide,
) -> std::result::Result<&[AlignedPair], &'static str> {
    match side {
        AffixSide::Prefix => {
            let boundary = pairs
                .iter()
                .rposition(|p| p.0.is_some())
                .ok_or("No match found")?;
            if boundary + 1 >= pairs.len() {
                return Err("Query too short after trimming");
            }
            Ok(&pairs[boundary + 1..])
        }
        AffixSide::Suffix => {
            let boundary = pairs
                .iter()
                .position(|p| p.0.is_some())
                .ok_or("No match found")?;
            if boundary == 0 {
                return Err("Starts too early");
            }
            Ok(&pairs[..boundary])
        }
    }
}

/// Align `affix` to `query` and cut it off. Returns the remaining query
/// codes, or `None` (with the reason recorded) when trimming fails.
fn trim_one_affix(
    affix: &EncodedDna,
    query: EncodedDna,
    side: AffixSide,
    min_quality: f64,
    subject: &str,
    stats: &mut AlignmentStats,
) -> Option<EncodedDna> {
    if affix.is_empty() {
        return Some(query);
    }
    let label = side.label();
    let alignment = semiglobal(affix, &query, GAP_OPEN, GAP_EXTEND);
    let affix_columns: Vec<&AlignedPair> =
        alignment.pairs.iter().filter(|p| p.0.is_some()).collect();
    if affix_columns.is_empty() {
        stats.fail(
            subject,
            Some(side),
            format!("No {label} content in alignment"),
        );
        return None;
    }
    let matches = affix_columns.iter().filter(|p| p.0 == p.1).count();
    let quality = matches as f64 / affix_columns.len() as f64;
    if quality < min_quality {
        let percent = format_float(round_digits(quality * 100.0, 1));
        stats.fail(
            subject,
            Some(side),
            format!("Poor {label} alignment quality ({percent}% match)"),
        );
        return None;
    }
    let core = match core_columns(&alignment.pairs, side) {
        Ok(core) => core,
        Err(reason) => {
            stats.fail(subject, Some(side), format!("{}: {reason}", side.title()));
            return None;
        }
    };
    let remaining: Vec<u8> = core.iter().filter_map(|p| p.1).collect();
    if remaining.is_empty() {
        stats.fail(
            subject,
            Some(side),
            format!("Empty core after {label} trimming"),
        );
        return None;
    }
    Some(EncodedDna::from_codes(remaining))
}

/// Outcome of trimming one clustered read segment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Trimmed {
    /// Trimmed sequence; empty when trimming failed.
    pub sequence: String,
    /// Edit distance to the reference allele; -1 when trimming failed.
    pub distance: i64,
}

impl Trimmed {
    fn failed() -> Self {
        Trimmed {
            sequence: String::new(),
            distance: -1,
        }
    }
}

/// Trim `prefix` and `suffix` off `query` and measure the edit distance of
/// what remains to `reference`.
pub fn trim_and_align(
    query: &str,
    prefix: &str,
    suffix: &str,
    reference: &str,
    min_quality: f64,
    subject: &str,
    stats: &mut AlignmentStats,
) -> Result<Trimmed> {
    let encode = |what: &str, text: &str| {
        EncodedDna::new(text).ok_or_else(|| {
            Error::invalid(format!(
                "{what} for {subject} contains a non-nucleotide character: {text}"
            ))
        })
    };
    let query = encode("Query sequence", query)?;
    let prefix = encode("Prefix", prefix)?;
    let suffix = encode("Suffix", suffix)?;

    if query.is_empty() {
        stats.fail(subject, None, "Empty query sequence".to_string());
        return Ok(Trimmed::failed());
    }
    stats.total_attempts += 1;

    let Some(without_prefix) = trim_one_affix(
        &prefix,
        query,
        AffixSide::Prefix,
        min_quality,
        subject,
        stats,
    ) else {
        return Ok(Trimmed::failed());
    };
    let Some(core) = trim_one_affix(
        &suffix,
        without_prefix,
        AffixSide::Suffix,
        min_quality,
        subject,
        stats,
    ) else {
        return Ok(Trimmed::failed());
    };
    if core.len() < MIN_CORE_LENGTH {
        // Counted as a prefix failure, as in the original.
        stats.fail(
            subject,
            Some(AffixSide::Prefix),
            format!(
                "Core sequence too short ({} < {MIN_CORE_LENGTH} nt)",
                core.len()
            ),
        );
        return Ok(Trimmed::failed());
    }
    stats.successful += 1;
    let sequence = EncodedDna::decode(core.codes());
    let distance = levenshtein(sequence.as_bytes(), reference.as_bytes()) as i64;
    Ok(Trimmed { sequence, distance })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn trim(query: &str, prefix: &str, suffix: &str, reference: &str) -> (Trimmed, AlignmentStats) {
        let mut stats = AlignmentStats::default();
        let out =
            trim_and_align(query, prefix, suffix, reference, 0.75, "S*01", &mut stats).unwrap();
        (out, stats)
    }

    #[test]
    fn exact_affixes_are_removed() {
        let core = "ACGTTGCAAGGCTTAACCGG";
        let query = format!("TTTTTTGGGG{core}CCCCAAAAAT");
        let (out, stats) = trim(&query, "TTTTTTGGGG", "CCCCAAAAAT", core);
        assert_eq!(out.sequence, core);
        assert_eq!(out.distance, 0);
        assert_eq!((stats.total_attempts, stats.successful), (1, 1));
    }

    #[test]
    fn empty_affixes_leave_the_query_untouched() {
        let (out, _) = trim("ACGTACGTACGT", "", "", "ACGTACGTACGA");
        assert_eq!(out.sequence, "ACGTACGTACGT");
        assert_eq!(out.distance, 1);
    }

    #[test]
    fn poor_affix_alignment_fails() {
        let (out, stats) = trim("ACGTACGTACGTACGTACGT", "GGGGGGGGGG", "", "ACGT");
        assert_eq!(out, Trimmed::failed());
        assert_eq!(stats.prefix_failures, 1);
        assert!(stats.failures[0]
            .1
            .starts_with("Poor prefix alignment quality"));
    }

    #[test]
    fn short_core_and_bad_characters() {
        let (out, stats) = trim(
            "TTTTTTGGGGACGTACCCCAAAAAT",
            "TTTTTTGGGG",
            "CCCCAAAAAT",
            "ACGTA",
        );
        assert_eq!(out.distance, -1);
        assert_eq!(stats.prefix_failures, 1);
        let mut s = AlignmentStats::default();
        assert!(trim_and_align("AC!T", "", "", "ACGT", 0.75, "x", &mut s).is_err());
    }
}
