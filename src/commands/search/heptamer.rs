//! `search heptamer`: locate the recombination-signal heptamer downstream
//! of V alleles and report the allele extended up to it.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use rayon::prelude::*;

use crate::align::hamming;
use crate::error::{Error, Result};
use crate::frame::{Column, Frame, SortKey};
use crate::info;
use crate::io::compress::read_to_string;
use crate::io::path::always_gz;
use crate::io::tsv::{write_tsv, write_tsv_gz};
use crate::seq::search::{find_first, PatternSet};
use crate::seq::slice::range;
use crate::seq::{load_database, load_demultiplex, Allele, DemuxTable};

const HEPTAMER_LENGTH: usize = 7;

/// Heptamer sets shipped with the original repository (`heptamers.json`),
/// used when no file is given and none exists in the working directory.
const BUILTIN_HEPTAMERS: &str = r#"{"IGHV":["CACAGTG","CACAATG","CACAGAG","CACGGTG","CACAGCG"],"IGLV":["CACAGTG","CACGGTG","CATGGTG","CACGCTG","CACAGCG","CACAGTA","CATAGTG","CACAATG"],"IGKV":["CACAGTG","CACACTG","CACTGTG","CACGGTG","CACAATG","CACATTG"]}"#;

/// Default heptamer file name, looked up in the working directory.
pub const DEFAULT_HEPTAMER_FILE: &str = "heptamers.json";

#[derive(Debug, Clone)]
pub struct HeptamerParams {
    pub tsv: PathBuf,
    pub fasta: PathBuf,
    pub output: PathBuf,
    pub summary: PathBuf,
    pub json: PathBuf,
    pub chain: String,
    pub max_distance: usize,
    /// Bases trimmed from the 5' end of each allele before searching.
    pub trim_begin: usize,
    /// Bases trimmed from the 3' end of each allele before searching.
    pub trim_end: usize,
    pub mincount: i64,
    pub ratio: f64,
}

/// Load `chain -> heptamers` from JSON.
pub fn load_heptamers(path: &Path) -> Result<HashMap<String, Vec<String>>> {
    let text = if path.exists() || path != Path::new(DEFAULT_HEPTAMER_FILE) {
        read_to_string(path)?
    } else {
        info!(
            "{DEFAULT_HEPTAMER_FILE} not found in the working directory; using built-in heptamers"
        );
        BUILTIN_HEPTAMERS.to_string()
    };
    serde_json::from_str(&text).map_err(|e| Error::parse(path, e.to_string()))
}

/// 1-based position of the first window of `suffix` within `max_distance`
/// mismatches of any heptamer, preferring smaller distances. Every window
/// is examined, including the one that ends at the last base.
pub fn find_heptamer(suffix: &str, heptamers: &[String], max_distance: usize) -> Option<usize> {
    let bytes = suffix.as_bytes();
    let n_windows = (bytes.len() + 1).saturating_sub(HEPTAMER_LENGTH);
    for distance in 0..=max_distance {
        for k in 0..n_windows {
            let window = &bytes[k..k + HEPTAMER_LENGTH];
            if heptamers
                .iter()
                .any(|h| hamming(h.as_bytes(), window).is_some_and(|d| d <= distance))
            {
                return Some(k + 1);
            }
        }
    }
    None
}

struct HeptamerHit {
    full_match: String,
    trimmed_match: String,
    sequence: String,
    heptamer: String,
}

fn julia_range(first: usize, last: usize) -> String {
    format!("{first}:{last}")
}

/// For every allele and every read containing its trimmed sequence, find
/// the heptamer downstream and extend the allele up to it.
pub fn extract_heptamers(
    table: &DemuxTable,
    db: &[Allele],
    heptamers: &[String],
    max_distance: usize,
    trim_begin: usize,
    trim_end: usize,
) -> Result<Frame> {
    if let Some(bad) = heptamers.iter().find(|h| h.len() != HEPTAMER_LENGTH) {
        return Err(Error::invalid(format!(
            "Heptamer '{bad}' is not {HEPTAMER_LENGTH} nt long"
        )));
    }
    let trimmed: Vec<&str> = db
        .iter()
        .map(|allele| {
            let len = allele.sequence.len() as i64;
            let query = range(
                &allele.sequence,
                trim_begin as i64 + 1,
                len - trim_end as i64,
            );
            if query.is_empty() {
                Err(Error::invalid(format!(
                    "Allele {} is too short to trim {trim_begin} + {trim_end} bases",
                    allele.name
                )))
            } else {
                Ok(query)
            }
        })
        .collect::<Result<_>>()?;
    let patterns = PatternSet::new(trimmed.iter().map(|q| q.as_bytes()))?;

    let per_read: Vec<Vec<(usize, HeptamerHit)>> = table
        .sequence
        .par_iter()
        .map(|read| {
            patterns
                .first_occurrences(read)
                .into_iter()
                .map(|hit| {
                    let suffix = &read[hit.end..];
                    let (extra, heptamer) = match find_heptamer(suffix, heptamers, max_distance) {
                        Some(loc) => (
                            &suffix[..loc - 1],
                            &suffix[loc - 1..loc - 1 + HEPTAMER_LENGTH],
                        ),
                        None => ("", ""),
                    };
                    let full_match = find_first(&db[hit.pattern].sequence, read)
                        .map(|(first, last)| julia_range(first, last))
                        .unwrap_or_default();
                    (
                        hit.pattern,
                        HeptamerHit {
                            full_match,
                            trimmed_match: julia_range(hit.start, hit.end),
                            sequence: format!("{}{}", trimmed[hit.pattern], extra),
                            heptamer: heptamer.to_string(),
                        },
                    )
                })
                .collect()
        })
        .collect();

    // Output is ordered by allele (database order), then by read.
    let mut per_allele: Vec<Vec<HeptamerHit>> = (0..db.len()).map(|_| Vec::new()).collect();
    for (allele, hit) in per_read.into_iter().flatten() {
        per_allele[allele].push(hit);
    }
    let rows: Vec<(usize, &HeptamerHit)> = per_allele
        .iter()
        .enumerate()
        .flat_map(|(allele, hits)| hits.iter().map(move |hit| (allele, hit)))
        .collect();

    Frame::from_columns([
        (
            "db_name",
            Column::from_strings(rows.iter().map(|(a, _)| db[*a].name.as_str())),
        ),
        (
            "full_match",
            Column::from_strings(rows.iter().map(|(_, h)| h.full_match.as_str())),
        ),
        (
            "trimmed_match",
            Column::from_strings(rows.iter().map(|(_, h)| h.trimmed_match.as_str())),
        ),
        (
            "db_length",
            Column::from_ints(rows.iter().map(|(a, _)| db[*a].sequence.len() as i64)),
        ),
        (
            "full_length",
            Column::from_ints(rows.iter().map(|(_, h)| h.sequence.len() as i64)),
        ),
        (
            "sequence",
            Column::from_strings(rows.iter().map(|(_, h)| h.sequence.as_str())),
        ),
        (
            "heptamer",
            Column::from_strings(rows.iter().map(|(_, h)| h.heptamer.as_str())),
        ),
    ])
}

/// Collapse identical (allele, extended sequence, heptamer) rows and keep
/// those passing the count and within-allele ratio thresholds.
pub fn summarize(hits: &Frame, min_ratio: f64, min_count: i64) -> Result<Frame> {
    if hits.nrows() == 0 {
        return Ok(Frame::new());
    }
    let by = ["db_name", "sequence", "heptamer"];
    let groups = hits.group_by(&by)?;
    let mut summary = hits.take(&groups.first_rows()).select(&by)?;
    let counts: Vec<i64> = groups.iter().map(|rows| rows.len() as i64).collect();
    summary.push_column("count", Column::from_ints(counts.iter().copied()))?;
    let ratios: Vec<f64> = summary.group_by(&["db_name"])?.map_rows(|rows| {
        let max = rows.iter().map(|&r| counts[r]).max().unwrap_or(1) as f64;
        rows.iter().map(|&r| counts[r] as f64 / max).collect()
    });
    summary.push_column("ratio", Column::from_floats(ratios.iter().copied()))?;
    let mask: Vec<bool> = (0..summary.nrows())
        .map(|row| counts[row] >= min_count && ratios[row] >= min_ratio)
        .collect();
    summary
        .filter(&mask)
        .sorted(&[SortKey::asc("db_name"), SortKey::desc("count")])
}

/// Run the command.
pub fn run(params: &HeptamerParams) -> Result<()> {
    info!("Heptamer search");
    let table = load_demultiplex(&params.tsv, None)?;
    let db = load_database(&params.fasta)?;
    let all = load_heptamers(&params.json)?;
    let heptamers = all.get(&params.chain).ok_or_else(|| {
        Error::invalid(format!(
            "No heptamers for chain {} in {}",
            params.chain,
            params.json.display()
        ))
    })?;

    let hits = extract_heptamers(
        &table,
        &db,
        heptamers,
        params.max_distance,
        params.trim_begin,
        params.trim_end,
    )?;
    let output = always_gz(&params.output);
    write_tsv_gz(&output, &hits)?;
    info!("Heptamer search results saved to {}", output.display());

    let summary = summarize(&hits, params.ratio, params.mincount)?;
    write_tsv(&params.summary, &summary)?;
    info!("Heptamer summary saved to {}", params.summary.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn heptamers() -> Vec<String> {
        vec!["CACAGTG".to_string(), "CACAATG".to_string()]
    }

    #[test]
    fn exact_heptamer_is_preferred_over_closer_inexact_one() {
        // Position 2 is one mismatch away, position 10 is exact.
        let suffix = "TCACTGTGGCACAGTGAAAA";
        assert_eq!(find_heptamer(suffix, &heptamers(), 1), Some(10));
        assert_eq!(find_heptamer("TCACTGTGGAAAAAAAA", &heptamers(), 1), Some(2));
        assert_eq!(find_heptamer("TCACTGTGGAAAAAAAA", &heptamers(), 0), None);
    }

    #[test]
    fn window_ending_at_the_last_base_is_examined() {
        assert_eq!(find_heptamer("AACACAGTG", &heptamers(), 0), Some(3));
        assert_eq!(find_heptamer("CACAGT", &heptamers(), 1), None);
        assert_eq!(find_heptamer("AACACAGTGA", &heptamers(), 0), Some(3));
    }
}
