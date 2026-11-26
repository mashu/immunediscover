//! `fasta merge`: merge FASTA files, keeping each distinct sequence once.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use crate::commands::table::common::{is_match, regex, remove_matches};
use crate::error::Result;
use crate::io::fasta::{read_fasta, write_fasta, DEFAULT_LINE_WIDTH};
use crate::{info, warn};

#[derive(Debug, Clone)]
pub struct MergeParams {
    pub output: PathBuf,
    pub inputs: Vec<PathBuf>,
    pub sort_by_name: bool,
    /// Regex removed from sequence names.
    pub cleanup: Option<String>,
    /// For duplicate sequences keep the first name (else the last).
    pub prefer_first: bool,
    /// Prefix names with the source file name.
    pub add_source_prefix: bool,
}

/// One distinct sequence and where its name came from.
struct Entry {
    sequence: String,
    name: String,
    /// Every `(name, file)` this sequence was seen under, if more than one.
    duplicates: Vec<(String, String)>,
    source: String,
}

fn file_stem(path: &Path) -> String {
    path.file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default()
}

/// Name given to a record: its identifier, followed by the rest of the
/// header line when there is one.
fn record_name(identifier: &str, comment: &str) -> String {
    if comment.is_empty() {
        identifier.to_string()
    } else {
        format!("{identifier} {comment}")
    }
}

/// Merge the inputs and write the result; returns the number of distinct
/// sequences. Output keeps first-seen order unless sorted by name.
pub fn run(params: &MergeParams) -> Result<usize> {
    info!(
        "Merging {} FASTA files into {}",
        params.inputs.len(),
        params.output.display()
    );
    let cleanup = params.cleanup.as_deref().map(regex).transpose()?;
    let mut entries: Vec<Entry> = Vec::new();
    let mut index: HashMap<String, usize> = HashMap::new();
    let mut total = 0usize;

    for (i, input) in params.inputs.iter().enumerate() {
        info!(
            "Processing file {}/{}: {}",
            i + 1,
            params.inputs.len(),
            input.display()
        );
        if !input.is_file() {
            warn!("File not found: {} - skipping", input.display());
            continue;
        }
        let source = input.display().to_string();
        let records = read_fasta(input)?;
        let count = records.len();
        for record in records {
            let mut name = record_name(record.identifier(), record.comment());
            if let Some(regex) = cleanup.as_ref().filter(|r| is_match(r, &name)) {
                name = remove_matches(regex, &name).trim().to_string();
            }
            if params.add_source_prefix {
                name = format!("{}_{name}", file_stem(input));
            }
            match index.get(&record.sequence) {
                Some(&slot) => {
                    let entry = &mut entries[slot];
                    if entry.duplicates.is_empty() {
                        entry
                            .duplicates
                            .push((entry.name.clone(), entry.source.clone()));
                    }
                    entry.duplicates.push((name.clone(), source.clone()));
                    if !params.prefer_first {
                        entry.name = name;
                        entry.source = source.clone();
                    }
                }
                None => {
                    index.insert(record.sequence.clone(), entries.len());
                    entries.push(Entry {
                        sequence: record.sequence,
                        name,
                        duplicates: Vec::new(),
                        source: source.clone(),
                    });
                }
            }
            total += 1;
        }
        info!("  Loaded {count} sequences from {}", input.display());
    }
    info!("Total sequences processed: {total}");
    info!("Unique sequences found: {}", entries.len());

    let mut order: Vec<usize> = (0..entries.len()).collect();
    if params.sort_by_name {
        info!("Sorting sequences by name");
        order.sort_by(|&a, &b| entries[a].name.cmp(&entries[b].name));
    }
    info!(
        "Writing {} unique sequences to {}",
        entries.len(),
        params.output.display()
    );
    write_fasta(
        &params.output,
        order
            .iter()
            .map(|&i| (entries[i].name.as_str(), entries[i].sequence.as_str())),
        Some(DEFAULT_LINE_WIDTH),
    )?;
    info!(
        "Successfully merged FASTA files: {} duplicates removed",
        total - entries.len()
    );

    if entries.iter().any(|e| !e.duplicates.is_empty()) {
        info!("Duplicate sequences found:");
        for entry in entries.iter().filter(|e| !e.duplicates.is_empty()) {
            info!("  Sequence with {} duplicates:", entry.duplicates.len());
            for (name, source) in &entry.duplicates {
                let file = Path::new(source)
                    .file_name()
                    .map(|f| f.to_string_lossy().into_owned())
                    .unwrap_or_default();
                info!("    - {name} (from {file})");
            }
        }
    }
    Ok(entries.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_text_is_kept_once() {
        assert_eq!(record_name("A*01", ""), "A*01");
        assert_eq!(record_name("A*01", "functional"), "A*01 functional");
    }
}
