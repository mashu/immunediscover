//! `table unique`, `table sort`, `table select`, `table aggregate`.

use std::path::Path;

use super::common::{load, require_columns, save_gz};
use crate::error::Result;
use crate::frame::{Column, Frame, SortKey};
use crate::info;

/// Distinct rows over `columns` (only those columns are kept).
pub fn unique(input: &Path, output: &Path, columns: &[String]) -> Result<Frame> {
    let table = load(input)?;
    require_columns(&table, columns)?;
    let distinct = table.select(columns)?.unique();
    info!(
        "Unique rows: {} (removed {} duplicates)",
        distinct.nrows(),
        table.nrows() - distinct.nrows()
    );
    save_gz(output, &distinct)?;
    Ok(distinct)
}

/// Stable sort by `columns`, all ascending or all descending.
pub fn sort(input: &Path, output: &Path, columns: &[String], descending: bool) -> Result<Frame> {
    let table = load(input)?;
    require_columns(&table, columns)?;
    let keys: Vec<SortKey> = columns
        .iter()
        .map(|c| SortKey {
            column: c.clone(),
            descending,
        })
        .collect();
    let sorted = table.sorted(&keys)?;
    info!(
        "Sorted {} rows by {columns:?} ({})",
        sorted.nrows(),
        if descending { "desc" } else { "asc" }
    );
    save_gz(output, &sorted)?;
    Ok(sorted)
}

/// Projection onto `columns`.
pub fn select(input: &Path, output: &Path, columns: &[String]) -> Result<Frame> {
    let table = load(input)?;
    require_columns(&table, columns)?;
    let selected = table.select(columns)?;
    info!(
        "Selected {} columns from {}",
        selected.ncols(),
        table.ncols()
    );
    save_gz(output, &selected)?;
    Ok(selected)
}

/// Count rows per group; additional columns carry the group's first value.
///
/// `keep_columns = None` keeps every non-group column.
pub fn aggregate(
    input: &Path,
    output: &Path,
    group_by: &[String],
    keep_columns: Option<&[String]>,
    count_column: &str,
) -> Result<Frame> {
    let table = load(input)?;
    require_columns(&table, group_by)?;
    let keep: Vec<String> = match keep_columns {
        Some(columns) => {
            require_columns(&table, columns)?;
            columns.to_vec()
        }
        None => table
            .names()
            .iter()
            .filter(|n| !group_by.contains(n))
            .cloned()
            .collect(),
    };

    let groups = table.group_by(group_by)?;
    let first_rows = groups.first_rows();
    let mut grouped = table.take(&first_rows).select(group_by)?;
    grouped.push_column(
        count_column,
        Column::from_ints(groups.iter().map(|rows| rows.len() as i64)),
    )?;

    let carried: Vec<&String> = keep.iter().filter(|c| !grouped.has_column(c)).collect();
    if !carried.is_empty() {
        // Each kept column takes the value of the group's first row; groups
        // with missing keys are groups like any other.
        let firsts = table.take(&first_rows);
        for column in carried {
            if !grouped.has_column(column) {
                grouped.push_column(column.as_str(), firsts.column(column)?.clone())?;
            }
        }
    }
    info!(
        "Aggregated from {} rows to {} unique groups",
        table.nrows(),
        grouped.nrows()
    );
    save_gz(output, &grouped)?;
    Ok(grouped)
}
