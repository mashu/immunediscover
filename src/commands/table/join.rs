//! `table outerjoin` / `table leftjoin`.

use std::path::{Path, PathBuf};

use super::common::save_gz;
use crate::error::Result;
use crate::frame::{join, Frame, JoinKind};
use crate::info;
use crate::io::tsv::read_tsv;

#[derive(Debug, Clone)]
pub struct JoinParams {
    pub left: PathBuf,
    pub right: PathBuf,
    pub output: PathBuf,
    pub kind: JoinKind,
    pub left_keys: Vec<String>,
    /// Defaults to `left_keys`.
    pub right_keys: Option<Vec<String>>,
    /// Prefix (joined with `_`) for non-key columns of the left table.
    pub left_prefix: Option<String>,
    pub right_prefix: Option<String>,
    /// Columns to keep from the left table before joining.
    pub left_select: Option<Vec<String>>,
    pub right_select: Option<Vec<String>>,
}

fn prepare(path: &Path, side: &str, select: &Option<Vec<String>>) -> Result<Frame> {
    info!("Loading {side} file: {}", path.display());
    let table = read_tsv(path)?;
    match select {
        Some(columns) => {
            info!("Selecting columns from {side} file: {columns:?}");
            table.select(columns)
        }
        None => Ok(table),
    }
}

fn add_prefix(table: &mut Frame, keys: &[String], prefix: &Option<String>) -> Result<()> {
    let Some(prefix) = prefix.as_deref().filter(|p| !p.is_empty()) else {
        return Ok(());
    };
    for name in table.names().to_vec() {
        if !keys.contains(&name) {
            table.rename(&name, format!("{prefix}_{name}"))?;
        }
    }
    Ok(())
}

/// Join two TSV files and write the result.
pub fn run(params: &JoinParams) -> Result<Frame> {
    info!("Performing {} join", params.kind.label());
    let mut left = prepare(&params.left, "left", &params.left_select)?;
    let mut right = prepare(&params.right, "right", &params.right_select)?;
    let right_keys = params
        .right_keys
        .clone()
        .unwrap_or_else(|| params.left_keys.clone());
    info!("Left file: {} rows, {} columns", left.nrows(), left.ncols());
    info!(
        "Right file: {} rows, {} columns",
        right.nrows(),
        right.ncols()
    );
    info!(
        "Join keys - Left: {:?}, Right: {right_keys:?}",
        params.left_keys
    );

    add_prefix(&mut left, &params.left_keys, &params.left_prefix)?;
    add_prefix(&mut right, &right_keys, &params.right_prefix)?;
    let result = join(&left, &right, &params.left_keys, &right_keys, params.kind)?;
    info!(
        "Result: {} rows, {} columns",
        result.nrows(),
        result.ncols()
    );
    let output = save_gz(&params.output, &result)?;
    info!("Output saved to: {}", output.display());
    Ok(result)
}
