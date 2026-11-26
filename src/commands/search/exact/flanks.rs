//! Flanking-sequence extraction around an exact allele match.
//!
//! Which flanks exist depends on the gene segment type:
//!
//! * V: recombination signal (heptamer, 23-bp spacer, nonamer) downstream;
//! * J: recombination signal (nonamer, 23-bp spacer, heptamer) upstream;
//! * D: signals with 12-bp spacers on both sides.
//!
//! All coordinates are 1-based and inclusive, as in the original code.
//! An element that lies partly outside the read is truncated; one that lies
//! entirely outside it is empty.

use std::str::FromStr;

use crate::error::{Error, Result};
use crate::seq::slice::range;
use crate::seq::GeneType;

/// A recombination signal element selectable with `--rss`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RssElement {
    Heptamer,
    Spacer,
    Nonamer,
}

impl RssElement {
    pub fn column(self) -> &'static str {
        match self {
            RssElement::Heptamer => "heptamer",
            RssElement::Spacer => "spacer",
            RssElement::Nonamer => "nonamer",
        }
    }

    /// Parse a comma-separated list such as `heptamer,spacer`.
    pub fn parse_list(list: &str) -> Result<Vec<RssElement>> {
        let elements = list
            .split(',')
            .map(RssElement::from_str)
            .collect::<Result<Vec<_>>>()?;
        if elements.is_empty() {
            return Err(Error::invalid("At least one type must be specified."));
        }
        Ok(elements)
    }
}

impl FromStr for RssElement {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        match s {
            "heptamer" => Ok(RssElement::Heptamer),
            "spacer" => Ok(RssElement::Spacer),
            "nonamer" => Ok(RssElement::Nonamer),
            other => Err(Error::invalid(format!(
                "Invalid type: {other}. Allowed types are: heptamer, spacer, nonamer."
            ))),
        }
    }
}

/// What to extract next to a match.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FlankMode {
    /// Recombination signal elements (the listed ones are reported for V/J).
    Rss(Vec<RssElement>),
    /// A fixed-length extension on the signal side(s).
    Extension(usize),
}

impl FlankMode {
    pub fn all_rss() -> Self {
        FlankMode::Rss(vec![
            RssElement::Heptamer,
            RssElement::Spacer,
            RssElement::Nonamer,
        ])
    }

    pub fn extension(&self) -> Option<usize> {
        match self {
            FlankMode::Extension(length) => Some(*length),
            FlankMode::Rss(_) => None,
        }
    }
}

/// Named flank strings in output-column order (includes `sequence`).
pub type Flanks = Vec<(&'static str, String)>;

/// Column order of the V/J recombination-signal fields (upstream's
/// `full_tuple_keys_as_strings`).
const RSS_FIELD_ORDER: [&str; 6] = [
    "prefix", "sequence", "heptamer", "spacer", "nonamer", "suffix",
];

fn own(s: &str) -> String {
    s.to_string()
}

/// The `affix` bases before the match (fewer if the read starts earlier).
fn leading_affix(read: &str, start: i64, affix: i64) -> &str {
    if start > affix {
        range(read, start - affix, start - 1)
    } else {
        range(read, 1, start - 1)
    }
}

/// Recombination-signal flanks of the match `read[start..=end]`.
pub fn rss_flanks(read: &str, start: usize, end: usize, gene: GeneType, affix: usize) -> Flanks {
    let (sp, ep, n, len) = (start as i64, end as i64, affix as i64, read.len() as i64);
    let sequence = own(range(read, sp, ep));
    match gene {
        GeneType::V => vec![
            ("prefix", own(leading_affix(read, sp, n))),
            ("sequence", sequence),
            ("heptamer", own(range(read, ep + 1, (ep + 7).min(len)))),
            ("spacer", own(range(read, ep + 8, (ep + 30).min(len)))),
            ("nonamer", own(range(read, ep + 31, (ep + 39).min(len)))),
        ],
        GeneType::J => vec![
            ("nonamer", own(range(read, (sp - 39).max(1), sp - 31))),
            ("spacer", own(range(read, (sp - 30).max(1), sp - 8))),
            ("heptamer", own(range(read, (sp - 7).max(1), sp - 1))),
            ("suffix", own(range(read, ep + 1, (ep + n).min(len)))),
            ("sequence", sequence),
        ],
        GeneType::D => vec![
            ("pre_nonamer", own(range(read, (sp - 28).max(1), sp - 20))),
            ("pre_spacer", own(range(read, (sp - 19).max(1), sp - 8))),
            ("pre_heptamer", own(range(read, (sp - 7).max(1), sp - 1))),
            ("sequence", sequence),
            ("post_heptamer", own(range(read, ep + 1, (ep + 7).min(len)))),
            ("post_spacer", own(range(read, ep + 8, (ep + 19).min(len)))),
            (
                "post_nonamer",
                own(range(read, ep + 20, (ep + 28).min(len))),
            ),
        ],
    }
}

/// Keep the fields selected with `--rss` (V and J only; D keeps everything).
pub fn select_rss_fields(flanks: Flanks, selected: &[RssElement], gene: GeneType) -> Flanks {
    let mandatory = match gene {
        GeneType::V => "prefix",
        GeneType::J => "suffix",
        GeneType::D => return flanks,
    };
    let keep = |name: &str| {
        name == "sequence" || name == mandatory || selected.iter().any(|e| e.column() == name)
    };
    let mut kept: Flanks = flanks.into_iter().filter(|(name, _)| keep(name)).collect();
    kept.sort_by_key(|(name, _)| RSS_FIELD_ORDER.iter().position(|field| field == name));
    kept
}

/// Extension-mode flanks with separate lengths for the two sides.
///
/// On the side without a recombination signal (upstream of V, downstream
/// of J) the fixed `affix` length is used instead.
pub fn extension_flanks(
    read: &str,
    start: usize,
    end: usize,
    gene: GeneType,
    affix: usize,
    prefix_extension: usize,
    suffix_extension: usize,
) -> Flanks {
    let (sp, ep, n, len) = (start as i64, end as i64, affix as i64, read.len() as i64);
    let (pre, suf) = (prefix_extension as i64, suffix_extension as i64);
    let upstream = |length: i64| own(range(read, (sp - length).max(1), sp - 1));
    let downstream = |length: i64| own(range(read, ep + 1, (ep + length).min(len)));
    let (prefix, suffix) = match gene {
        GeneType::V => (own(leading_affix(read, sp, n)), downstream(suf)),
        GeneType::J => (upstream(pre), downstream(n)),
        GeneType::D => (upstream(pre), downstream(suf)),
    };
    vec![
        ("prefix", prefix),
        ("sequence", own(range(read, sp, ep))),
        ("suffix", suffix),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read() -> String {
        "ACGTACGTAC".repeat(12)
    }

    fn names(flanks: &Flanks) -> Vec<&'static str> {
        flanks.iter().map(|(n, _)| *n).collect()
    }

    fn value<'a>(flanks: &'a Flanks, name: &str) -> &'a str {
        &flanks.iter().find(|(n, _)| *n == name).unwrap().1
    }

    // Expected values produced by immunediscover.jl `Exact.extract_flanking`.
    #[test]
    fn v_flanks_near_read_end() {
        let f = rss_flanks(&read(), 115, 118, GeneType::V, 13);
        assert_eq!(value(&f, "prefix"), "CGTACGTACACGT");
        assert_eq!(value(&f, "sequence"), "ACGT");
        assert_eq!(value(&f, "heptamer"), "AC");
        assert_eq!(value(&f, "spacer"), "");
        assert_eq!(value(&f, "nonamer"), "");
    }

    #[test]
    fn v_flanks_near_read_start() {
        let f = rss_flanks(&read(), 3, 6, GeneType::V, 13);
        assert_eq!(value(&f, "prefix"), "AC");
        assert_eq!(value(&f, "heptamer"), "GTACACG");
        assert_eq!(value(&f, "spacer"), "TACGTACACGTACGTACACGTAC");
        assert_eq!(value(&f, "nonamer"), "GTACACGTA");
    }

    #[test]
    fn elements_before_the_read_start_are_empty() {
        let j = rss_flanks(&read(), 3, 6, GeneType::J, 13);
        assert_eq!(value(&j, "nonamer"), "");
        assert_eq!(value(&j, "spacer"), "");
        assert_eq!(value(&j, "heptamer"), "AC");
        assert_eq!(value(&j, "suffix"), "GTACACGTACGTA");
        let d = rss_flanks(&read(), 3, 6, GeneType::D, 13);
        assert_eq!(value(&d, "pre_nonamer"), "");
        assert_eq!(value(&d, "pre_spacer"), "");
        assert_eq!(value(&d, "pre_heptamer"), "AC");
        // Partly inside the read: truncated, not dropped.
        let j = rss_flanks(&read(), 35, 40, GeneType::J, 13);
        assert_eq!(value(&j, "nonamer"), "ACGT");
        assert_eq!(value(&j, "spacer"), "ACGTACACGTACGTACACGTACG");
        assert_eq!(value(&d, "post_heptamer"), "GTACACG");
        assert_eq!(value(&d, "post_spacer"), "TACGTACACGTA");
        assert_eq!(value(&d, "post_nonamer"), "CGTACACGT");
    }

    #[test]
    fn rss_field_selection_and_order() {
        use RssElement::*;
        let pick = |gene, selected: &[RssElement]| {
            names(&select_rss_fields(
                rss_flanks(&read(), 50, 60, gene, 13),
                selected,
                gene,
            ))
        };
        assert_eq!(
            pick(GeneType::V, &[Heptamer]),
            ["prefix", "sequence", "heptamer"]
        );
        assert_eq!(
            pick(GeneType::V, &[Spacer]),
            ["prefix", "sequence", "spacer"]
        );
        assert_eq!(
            pick(GeneType::V, &[Heptamer, Spacer, Nonamer]),
            ["prefix", "sequence", "heptamer", "spacer", "nonamer"]
        );
        assert_eq!(
            pick(GeneType::J, &[Heptamer]),
            ["sequence", "heptamer", "suffix"]
        );
        assert_eq!(
            pick(GeneType::J, &[Nonamer, Heptamer, Spacer]),
            ["sequence", "heptamer", "spacer", "nonamer", "suffix"]
        );
        assert_eq!(pick(GeneType::D, &[Heptamer]).len(), 7);
    }

    #[test]
    fn extension_flanks_use_affix_on_the_signal_free_side() {
        let v = extension_flanks(&read(), 50, 60, GeneType::V, 5, 9, 7);
        assert_eq!(
            (value(&v, "prefix").len(), value(&v, "suffix").len()),
            (5, 7)
        );
        let j = extension_flanks(&read(), 50, 60, GeneType::J, 5, 9, 7);
        assert_eq!(
            (value(&j, "prefix").len(), value(&j, "suffix").len()),
            (9, 5)
        );
        let d = extension_flanks(&read(), 50, 60, GeneType::D, 5, 9, 7);
        assert_eq!(
            (value(&d, "prefix").len(), value(&d, "suffix").len()),
            (9, 7)
        );
    }

    #[test]
    fn rss_list_parsing() {
        assert_eq!(
            RssElement::parse_list("heptamer,nonamer").unwrap(),
            vec![RssElement::Heptamer, RssElement::Nonamer]
        );
        assert!(RssElement::parse_list("heptamer,foo").is_err());
    }
}
