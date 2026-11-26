//! Running `bwa mem` and reading its alignments.
//!
//! The original links libbwa in-process. This port runs the `bwa`
//! executable instead and reads the SAM stream; `-a` reports every
//! alignment (the original inspects all of them), `-Y` keeps clipping soft
//! and `-T 0` disables the score cutoff that only the command-line tool
//! applies.

use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};

use crate::error::{Error, Result};

/// One alignment of a query sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomeAlignment {
    pub reference: String,
    /// 0-based leftmost reference position.
    pub position: i64,
    pub reverse: bool,
    pub cigar: String,
    /// Alignment score (`AS`).
    pub score: i64,
    /// Edit distance (`NM`).
    pub edit_distance: i64,
}

fn tag_value(fields: &[&str], tag: &str) -> Option<i64> {
    fields
        .iter()
        .skip(11)
        .find_map(|field| field.strip_prefix(tag))
        .and_then(|value| value.parse().ok())
}

/// Parse SAM text into per-query alignment lists (queries are named by
/// their index). Unmapped records contribute no alignment.
pub fn parse_sam(sam: &str, n_queries: usize) -> Result<Vec<Vec<GenomeAlignment>>> {
    let malformed = |line: &str| Error::tool("bwa", format!("unexpected SAM line: {line}"));
    let mut alignments: Vec<Vec<GenomeAlignment>> = vec![Vec::new(); n_queries];
    for line in sam.lines().filter(|l| !l.is_empty() && !l.starts_with('@')) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            return Err(malformed(line));
        }
        let query: usize = fields[0].parse().map_err(|_| malformed(line))?;
        let flag: u32 = fields[1].parse().map_err(|_| malformed(line))?;
        if flag & 0x4 != 0 {
            continue;
        }
        let position: i64 = fields[3].parse().map_err(|_| malformed(line))?;
        alignments
            .get_mut(query)
            .ok_or_else(|| malformed(line))?
            .push(GenomeAlignment {
                reference: fields[2].to_string(),
                position: position - 1,
                reverse: flag & 0x10 != 0,
                cigar: fields[5].to_string(),
                score: tag_value(&fields, "AS:i:").unwrap_or(0),
                edit_distance: tag_value(&fields, "NM:i:").unwrap_or(0),
            });
    }
    Ok(alignments)
}

/// Align `sequences` against the BWA index at `genome`.
/// Returns, for each sequence, its alignments in BWA's order.
pub fn align_all(genome: &Path, sequences: &[String]) -> Result<Vec<Vec<GenomeAlignment>>> {
    let directory = tempdir()?;
    let query_path = directory.join("queries.fasta");
    {
        let mut file = std::fs::File::create(&query_path).map_err(|e| Error::io(&query_path, e))?;
        for (i, sequence) in sequences.iter().enumerate() {
            // bwa cannot read empty records; a single N never aligns.
            let sequence = if sequence.is_empty() { "N" } else { sequence };
            writeln!(file, ">{i}\n{sequence}").map_err(|e| Error::io(&query_path, e))?;
        }
    }
    let threads = std::thread::available_parallelism().map_or(1, |n| n.get());
    let output = Command::new("bwa")
        .args(["mem", "-a", "-Y", "-T", "0", "-v", "1", "-t"])
        .arg(threads.to_string())
        .arg(genome)
        .arg(&query_path)
        .stdin(Stdio::null())
        .output()
        .map_err(|e| {
            Error::tool(
                "bwa",
                format!("could not be started ({e}); is bwa installed and on PATH?"),
            )
        })?;
    let _ = std::fs::remove_dir_all(&directory);
    if !output.status.success() {
        return Err(Error::tool(
            "bwa",
            format!(
                "exited with {}: {}",
                output.status,
                String::from_utf8_lossy(&output.stderr).trim()
            ),
        ));
    }
    parse_sam(&String::from_utf8_lossy(&output.stdout), sequences.len())
}

fn tempdir() -> Result<std::path::PathBuf> {
    let unique = format!(
        "immunediscover-bwa-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos())
    );
    let path = std::env::temp_dir().join(unique);
    std::fs::create_dir_all(&path).map_err(|e| Error::io(&path, e))?;
    Ok(path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sam_records_are_grouped_by_query() {
        let sam = "@SQ\tSN:chr14\tLN:1000\n\
                   0\t0\tchr14\t101\t60\t50M\t*\t0\t0\tACGT\t*\tNM:i:1\tMD:Z:50\tAS:i:45\n\
                   0\t272\tchr2\t7\t0\t5S45M\t*\t0\t0\t*\t*\tNM:i:0\tAS:i:45\n\
                   1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\tAS:i:0\n";
        let parsed = parse_sam(sam, 2).unwrap();
        assert_eq!(parsed[0].len(), 2);
        assert_eq!(parsed[0][0].position, 100);
        assert!(!parsed[0][0].reverse);
        assert_eq!((parsed[0][0].score, parsed[0][0].edit_distance), (45, 1));
        assert!(parsed[0][1].reverse);
        assert_eq!(parsed[0][1].cigar, "5S45M");
        assert!(parsed[1].is_empty());
    }
}
