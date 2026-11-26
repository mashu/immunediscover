//! `Frame`: a small column-oriented table.
//!
//! This is deliberately not a general data-frame library. It implements the
//! handful of relational operations the pipelines need, with the row-order
//! and missing-value semantics of Julia's DataFrames.jl so that outputs line
//! up with the original implementation.

use std::collections::HashSet;
use std::sync::Arc;

use super::column::Column;
use super::value::{Key, Value};
use crate::error::{Error, Result};

/// Distinct values of a pooled text column, in order of first appearance.
///
/// CSV.jl stores low-cardinality text columns as pooled arrays, and
/// DataFrames orders groups by pool position when every grouping column is
/// pooled or integer-like. The pool is therefore part of a column's
/// observable behaviour: it survives row selection (including levels whose
/// rows were filtered out) and is dropped when a column is rebuilt.
pub type Pool = Arc<Vec<Option<String>>>;

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Frame {
    names: Vec<String>,
    columns: Vec<Column>,
    pools: Vec<Option<Pool>>,
}

impl Frame {
    pub fn new() -> Self {
        Self::default()
    }

    /// Build from `(name, column)` pairs; all columns must be equally long.
    pub fn from_columns<I, S>(columns: I) -> Result<Self>
    where
        I: IntoIterator<Item = (S, Column)>,
        S: Into<String>,
    {
        let mut frame = Frame::new();
        for (name, column) in columns {
            frame.push_column(name, column)?;
        }
        Ok(frame)
    }

    pub fn nrows(&self) -> usize {
        self.columns.first().map_or(0, Column::len)
    }

    pub fn ncols(&self) -> usize {
        self.columns.len()
    }

    pub fn names(&self) -> &[String] {
        &self.names
    }

    pub fn has_column(&self, name: &str) -> bool {
        self.names.iter().any(|n| n == name)
    }

    pub fn index_of(&self, name: &str) -> Result<usize> {
        self.names
            .iter()
            .position(|n| n == name)
            .ok_or_else(|| Error::MissingColumn {
                column: name.to_string(),
                available: self.names.join(", "),
            })
    }

    pub fn column(&self, name: &str) -> Result<&Column> {
        Ok(&self.columns[self.index_of(name)?])
    }

    pub fn column_at(&self, index: usize) -> &Column {
        &self.columns[index]
    }

    pub fn columns(&self) -> &[Column] {
        &self.columns
    }

    /// Text column, with a typed error otherwise.
    pub fn strs(&self, name: &str) -> Result<&[Option<String>]> {
        let column = self.column(name)?;
        column.as_strs().ok_or_else(|| Error::ColumnType {
            column: name.to_string(),
            expected: "String",
            found: column.type_name(),
        })
    }

    /// Integer column, with a typed error otherwise. A column without any
    /// value (for instance in a table with no rows) reads as all-missing.
    pub fn ints(&self, name: &str) -> Result<Vec<Option<i64>>> {
        let column = self.column(name)?;
        if let Some(values) = column.as_ints() {
            return Ok(values.to_vec());
        }
        if column.iter().all(|v| v.is_missing()) {
            return Ok(vec![None; column.len()]);
        }
        Err(Error::ColumnType {
            column: name.to_string(),
            expected: "Int",
            found: column.type_name(),
        })
    }

    /// Numeric column widened to `f64` (all-missing columns are accepted).
    pub fn floats(&self, name: &str) -> Result<Vec<Option<f64>>> {
        let column = self.column(name)?;
        if let Some(values) = column.to_f64() {
            return Ok(values);
        }
        if column.iter().all(|v| v.is_missing()) {
            return Ok(vec![None; column.len()]);
        }
        Err(Error::ColumnType {
            column: name.to_string(),
            expected: "numeric",
            found: column.type_name(),
        })
    }

    pub fn value(&self, row: usize, name: &str) -> Result<Value<'_>> {
        Ok(self.column(name)?.get(row))
    }

    /// Append a column; fails on duplicate names or length mismatch.
    pub fn push_column(&mut self, name: impl Into<String>, column: Column) -> Result<()> {
        let name = name.into();
        if self.has_column(&name) {
            return Err(Error::invalid(format!("Duplicate column name '{name}'")));
        }
        if !self.columns.is_empty() && column.len() != self.nrows() {
            return Err(Error::invalid(format!(
                "Column '{name}' has {} rows, table has {}",
                column.len(),
                self.nrows()
            )));
        }
        self.names.push(name);
        self.columns.push(column);
        self.pools.push(None);
        Ok(())
    }

    /// Pool of a text column read from a file, if it qualified for pooling.
    pub fn pool_at(&self, index: usize) -> Option<&Pool> {
        self.pools[index].as_ref()
    }

    /// Mark column `name` as pooled with the given levels.
    pub fn set_pool(&mut self, name: &str, pool: Pool) -> Result<()> {
        let index = self.index_of(name)?;
        self.pools[index] = Some(pool);
        Ok(())
    }

    /// Replace an existing column or append a new one.
    pub fn set_column(&mut self, name: &str, column: Column) -> Result<()> {
        match self.names.iter().position(|n| n == name) {
            Some(index) => {
                if column.len() != self.nrows() {
                    return Err(Error::invalid(format!(
                        "Column '{name}' has {} rows, table has {}",
                        column.len(),
                        self.nrows()
                    )));
                }
                self.columns[index] = column;
                self.pools[index] = None;
                Ok(())
            }
            None => self.push_column(name, column),
        }
    }

    pub fn rename(&mut self, from: &str, to: impl Into<String>) -> Result<()> {
        let index = self.index_of(from)?;
        self.names[index] = to.into();
        Ok(())
    }

    /// Projection, in the requested order.
    pub fn select<S: AsRef<str>>(&self, names: &[S]) -> Result<Frame> {
        let mut out = Frame::new();
        for name in names {
            let name = name.as_ref();
            let index = self.index_of(name)?;
            out.push_column(name, self.columns[index].clone())?;
            *out.pools.last_mut().expect("a column was just pushed") = self.pools[index].clone();
        }
        Ok(out)
    }

    /// Rows `rows`, in that order.
    pub fn take(&self, rows: &[usize]) -> Frame {
        Frame {
            names: self.names.clone(),
            columns: self.columns.iter().map(|c| c.take(rows)).collect(),
            pools: self.pools.clone(),
        }
    }

    /// Rows whose mask entry is `true`.
    pub fn filter(&self, mask: &[bool]) -> Frame {
        Frame {
            names: self.names.clone(),
            columns: self.columns.iter().map(|c| c.filter(mask)).collect(),
            pools: self.pools.clone(),
        }
    }

    /// In-place row filter driven by a predicate over row indices.
    pub fn retain_rows(&mut self, mut keep: impl FnMut(usize) -> bool) {
        let mask: Vec<bool> = (0..self.nrows()).map(&mut keep).collect();
        *self = self.filter(&mask);
    }

    pub fn head(&self, n: usize) -> Frame {
        let rows: Vec<usize> = (0..self.nrows().min(n)).collect();
        self.take(&rows)
    }

    /// Hashable identity of `row` over the given column indices.
    pub fn row_key(&self, row: usize, column_indices: &[usize]) -> Vec<Key<'_>> {
        column_indices
            .iter()
            .map(|&c| self.columns[c].get(row).key())
            .collect()
    }

    /// Indices of the first occurrence of every distinct row (over `names`).
    pub fn unique_rows<S: AsRef<str>>(&self, names: &[S]) -> Result<Vec<usize>> {
        let indices = self.indices_of(names)?;
        let mut seen: HashSet<Vec<Key<'_>>> = HashSet::with_capacity(self.nrows());
        Ok((0..self.nrows())
            .filter(|&row| seen.insert(self.row_key(row, &indices)))
            .collect())
    }

    /// `unique(df, cols)`: keep the first row of every distinct key.
    pub fn unique_by<S: AsRef<str>>(&self, names: &[S]) -> Result<Frame> {
        Ok(self.take(&self.unique_rows(names)?))
    }

    /// `unique(df)`: distinct over all columns.
    pub fn unique(&self) -> Frame {
        let names = self.names.clone();
        self.unique_by(&names)
            .expect("a frame's own column names always resolve")
    }

    pub fn indices_of<S: AsRef<str>>(&self, names: &[S]) -> Result<Vec<usize>> {
        names.iter().map(|n| self.index_of(n.as_ref())).collect()
    }

    /// Vertical concatenation; column names must match exactly.
    pub fn vcat(frames: Vec<Frame>) -> Result<Frame> {
        let mut iter = frames.into_iter();
        let Some(first) = iter.next() else {
            return Ok(Frame::new());
        };
        let names = first.names;
        let mut columns = first.columns;
        for frame in iter {
            if frame.names != names {
                return Err(Error::invalid(
                    "Cannot concatenate tables with different column names",
                ));
            }
            columns = columns
                .into_iter()
                .zip(frame.columns)
                .map(|(a, b)| a.append(b))
                .collect();
        }
        let pools = vec![None; names.len()];
        Ok(Frame {
            names,
            columns,
            pools,
        })
    }

    /// Consume into `(name, column)` pairs.
    pub fn into_columns(self) -> impl Iterator<Item = (String, Column)> {
        self.names.into_iter().zip(self.columns)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> Frame {
        Frame::from_columns([
            ("a", Column::from_strings(["x", "x", "y"])),
            ("b", Column::from_ints([1, 1, 2])),
        ])
        .unwrap()
    }

    #[test]
    fn unique_keeps_first_occurrence() {
        let out = sample().unique();
        assert_eq!(out.nrows(), 2);
        assert_eq!(out.strs("a").unwrap()[1].as_deref(), Some("y"));
    }

    #[test]
    fn vcat_promotes_int_and_float() {
        let a = Frame::from_columns([("v", Column::from_ints([1]))]).unwrap();
        let b = Frame::from_columns([("v", Column::from_floats([2.5]))]).unwrap();
        let out = Frame::vcat(vec![a, b]).unwrap();
        assert_eq!(out.column("v").unwrap(), &Column::from_floats([1.0, 2.5]));
    }

    #[test]
    fn missing_column_error_lists_available() {
        let err = sample().column("zzz").unwrap_err().to_string();
        assert!(err.contains("zzz") && err.contains("a, b"));
    }
}
