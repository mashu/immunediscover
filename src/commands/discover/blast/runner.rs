//! Running NCBI BLAST+ (`blastn`, `makeblastdb`).

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::Instant;

use crate::error::{Error, Result};
use crate::info;
use crate::io::compress;
use crate::io::path::stem_before_first_dot;

use super::workdir::{content_digest, Fingerprint};

/// Columns requested from `blastn -outfmt 6`, in order.
pub const COLUMNS: [&str; 19] = [
    "qseqid", "sseqid", "pident", "nident", "length", "mismatch", "gapopen", "qcovs", "qcovhsp",
    "qstart", "qend", "sstart", "send", "qlen", "slen", "evalue", "bitscore", "sstrand", "qseq",
];

/// Supported `blastn` versions (inclusive).
pub const MIN_VERSION: (u32, u32, u32) = (2, 15, 0);
pub const MAX_VERSION: (u32, u32, u32) = (2, 17, 0);

fn parse_version(output: &str) -> Option<(u32, u32, u32)> {
    let rest = output.split("blastn:").nth(1)?.trim_start();
    let token: String = rest
        .chars()
        .take_while(|c| c.is_ascii_digit() || *c == '.')
        .collect();
    let mut parts = token.split('.').map(|p| p.parse::<u32>());
    Some((
        parts.next()?.ok()?,
        parts.next()?.ok()?,
        parts.next()?.ok()?,
    ))
}

/// Check that `blastn` is installed and within the supported version range.
pub fn verify_blastn_version() -> Result<()> {
    let version = |v: (u32, u32, u32)| format!("{}.{}.{}", v.0, v.1, v.2);
    let install_hint = format!(
        "Please install BLAST version {} - {}",
        version(MIN_VERSION),
        version(MAX_VERSION)
    );
    let output = Command::new("blastn")
        .arg("-version")
        .output()
        .map_err(|_| {
            Error::tool(
                "blastn",
                format!("blastn command not found in PATH. {install_hint}"),
            )
        })?;
    let current = parse_version(&String::from_utf8_lossy(&output.stdout)).ok_or_else(|| {
        Error::tool(
            "blastn",
            format!("Could not parse blastn version. {install_hint}"),
        )
    })?;
    if current < MIN_VERSION || current > MAX_VERSION {
        return Err(Error::tool(
            "blastn",
            format!(
                "blastn version {} is outside acceptable range ({} - {}). {install_hint}",
                version(current),
                version(MIN_VERSION),
                version(MAX_VERSION)
            ),
        ));
    }
    info!(
        "blastn version {} is within acceptable range",
        version(current)
    );
    Ok(())
}

/// Path (without extension) of the BLAST database for `fasta`, building it
/// with `makeblastdb` when missing or older than the FASTA.
pub fn ensure_blast_db(fasta: &Path) -> Result<PathBuf> {
    let db = fasta
        .parent()
        .unwrap_or_else(|| Path::new(""))
        .join(stem_before_first_dot(fasta));
    let mut index = db.as_os_str().to_os_string();
    index.push(".nin");
    // Built from exactly this FASTA? File times are not trusted: they are
    // too coarse on some file systems to tell two quick runs apart.
    let index = PathBuf::from(index);
    let fingerprint = Fingerprint::of(&[content_digest(fasta)?]);
    if fingerprint.matches(&index) {
        info!("Using existing BLAST database {}", db.display());
        return Ok(db);
    }
    info!(
        "Building BLAST database from {} (enables multithreaded blastn)",
        fasta.display()
    );
    let status = Command::new("makeblastdb")
        .arg("-in")
        .arg(fasta)
        .args(["-dbtype", "nucl", "-parse_seqids", "-out"])
        .arg(&db)
        .status()
        .map_err(|e| Error::tool("makeblastdb", format!("could not be started: {e}")))?;
    if !status.success() {
        return Err(Error::tool("makeblastdb", format!("exited with {status}")));
    }
    fingerprint.record(&index)?;
    Ok(db)
}

/// Thread count for `blastn`: `BLAST_NUM_THREADS` if set, else all CPUs.
pub fn blastn_threads() -> Result<usize> {
    match std::env::var("BLAST_NUM_THREADS") {
        Ok(value) if !value.trim().is_empty() => match value.trim().parse::<i64>() {
            Ok(n) if n >= 1 => Ok(n as usize),
            Ok(n) => Err(Error::invalid(format!(
                "BLAST_NUM_THREADS must be >= 1, got {n}"
            ))),
            Err(_) => Err(Error::invalid(format!(
                "BLAST_NUM_THREADS must be a positive integer, got {value:?}"
            ))),
        },
        _ => Ok(std::thread::available_parallelism().map_or(1, |n| n.get())),
    }
}

/// NCBI tools take single-dash long options; accept GNU-style `--task` too.
fn normalise_option(token: &str) -> String {
    match token.strip_prefix("--") {
        Some(rest) if !rest.is_empty() && !rest.starts_with('-') => format!("-{rest}"),
        _ => token.to_string(),
    }
}

/// Run `blastn`, streaming the tabular result gzip-compressed to `output_gz`.
pub fn blastn(query: &Path, database: &Path, output_gz: &Path, extra_args: &str) -> Result<()> {
    if !output_gz.to_string_lossy().ends_with(".gz") {
        return Err(Error::invalid(format!(
            "BLAST cache path must end with .gz, got {:?}",
            output_gz.display().to_string()
        )));
    }
    let threads = blastn_threads()?;
    let started = Instant::now();
    let mut child = Command::new("blastn")
        .arg("-num_threads")
        .arg(threads.to_string())
        .arg("-query")
        .arg(query)
        .arg("-db")
        .arg(database)
        .args(["-out", "-", "-outfmt"])
        .arg(format!("6 {}", COLUMNS.join(" ")))
        .args(extra_args.split_whitespace().map(normalise_option))
        .stdout(Stdio::piped())
        .spawn()
        .map_err(|e| Error::tool("blastn", format!("could not be started: {e}")))?;
    let mut stdout = child.stdout.take().expect("stdout was requested as a pipe");
    let copied = compress::create(output_gz, true).and_then(|mut sink| {
        std::io::copy(&mut stdout, &mut sink).map_err(|e| Error::io(output_gz, e))
    });
    if let Err(error) = copied {
        // The result cannot be written (disk full?): do not leave blastn
        // running, nor half a result behind.
        let _ = child.kill();
        let _ = child.wait();
        let _ = std::fs::remove_file(output_gz);
        return Err(error);
    }
    let status = child
        .wait()
        .map_err(|e| Error::tool("blastn", e.to_string()))?;
    if !status.success() {
        let _ = std::fs::remove_file(output_gz);
        return Err(Error::tool("blastn", format!("exited with {status}")));
    }
    info!(
        "BLASTn completed in {:.2} seconds (blastn num_threads={threads}, streamed to gzip)",
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_parsing_and_option_normalisation() {
        assert_eq!(
            parse_version("blastn: 2.17.0+\n Package: blast 2.17.0"),
            Some((2, 17, 0))
        );
        assert_eq!(parse_version("garbage"), None);
        assert_eq!(normalise_option("--task"), "-task");
        assert_eq!(normalise_option("-task"), "-task");
        assert_eq!(normalise_option("---x"), "---x");
        assert_eq!(normalise_option("megablast"), "megablast");
    }
}
