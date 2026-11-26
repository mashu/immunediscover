//! `table collect`: concatenate tables matched by a glob pattern.

use std::path::Path;

use crate::error::{Error, Result};
use crate::frame::{Column, Frame};
use crate::io::path::always_gz;
use crate::io::tsv::{read_tsv, write_tsv_gz};
use crate::{info, warn};

/// Concatenate all TSV files matching `pattern`, adding a `file` column.
/// The result is gzip-compressed and written to `output` with a `.gz`
/// suffix.
pub fn run(pattern: &str, output: &Path) -> Result<()> {
    info!("Collecting TSV files");
    let files: Vec<String> = glob::glob(pattern)
        .map_err(|e| Error::invalid(format!("Invalid glob pattern '{pattern}': {e}")))?
        .filter_map(|entry| entry.ok())
        .map(|path| path.to_string_lossy().into_owned())
        .collect();
    info!("Found {} files matching pattern {pattern}", files.len());
    if files.is_empty() {
        warn!("No files found matching pattern {pattern}");
        return Ok(());
    }
    let mut tables: Vec<Frame> = Vec::with_capacity(files.len());
    for file in &files {
        let mut table = read_tsv(Path::new(file))?;
        table.set_column(
            "file",
            Column::from_strings(std::iter::repeat_n(file.as_str(), table.nrows())),
        )?;
        if let Some(first) = tables.first() {
            if first.names() != table.names() {
                return Err(Error::invalid(format!(
                    "Column names in {file} do not match first file"
                )));
            }
        }
        tables.push(table);
    }
    let collected = Frame::vcat(tables)?;
    let output = always_gz(output);
    write_tsv_gz(&output, &collected)?;
    info!(
        "Collected {} rows into {}",
        collected.nrows(),
        output.display()
    );
    Ok(())
}
