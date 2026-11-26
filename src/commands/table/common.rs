//! Helpers shared by the table utilities.

use std::path::{Path, PathBuf};

use fancy_regex::Regex;

use crate::error::{Error, Result};
use crate::frame::Frame;
use crate::info;
use crate::io::path::always_gz;
use crate::io::tsv::{read_tsv, write_tsv_gz};

/// Split a comma-separated option value.
pub fn split_list(list: &str) -> Vec<String> {
    list.split(',').map(str::to_string).collect()
}

/// Load a TSV, logging its shape.
pub fn load(path: &Path) -> Result<Frame> {
    info!("Loading input file: {}", path.display());
    let table = read_tsv(path)?;
    info!(
        "Input file: {} rows, {} columns",
        table.nrows(),
        table.ncols()
    );
    Ok(table)
}

/// Fail early with the list of available columns.
pub fn require_columns(table: &Frame, columns: &[String]) -> Result<()> {
    columns
        .iter()
        .try_for_each(|c| table.index_of(c).map(|_| ()))
}

/// Write gzip-compressed, adding `.gz` to the name when missing.
pub fn save_gz(path: &Path, table: &Frame) -> Result<PathBuf> {
    let output = always_gz(path);
    write_tsv_gz(&output, table)?;
    Ok(output)
}

/// Compile a user-supplied regular expression.
pub fn regex(pattern: &str) -> Result<Regex> {
    Regex::new(pattern).map_err(|e| Error::Regex {
        pattern: pattern.to_string(),
        message: e.to_string(),
    })
}

/// `occursin(regex, text)`.
pub fn is_match(regex: &Regex, text: &str) -> bool {
    regex.is_match(text).unwrap_or(false)
}

/// `replace(text, regex => "")`: delete every match.
pub fn remove_matches(regex: &Regex, text: &str) -> String {
    regex.replace_all(text, "").into_owned()
}
