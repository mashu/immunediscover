//! Row grouping with DataFrames.jl's group ordering.
//!
//! DataFrames orders groups by first appearance, except when *every*
//! grouping column is integer-like with a compact value range; then a
//! reference-pool fast path is taken and groups come out sorted by value.
//! Output row order of several commands depends on this, so both paths are
//! reproduced here (rules taken from DataFrames 1.8 `refpool_and_array` and
//! `row_group_slots!`).

use std::cmp::Ordering;
use std::collections::HashMap;

use super::column::Column;
use super::table::Frame;
use super::value::{Key, Value};
use crate::error::Result;

/// Partition of a frame's rows into groups.
#[derive(Debug, Clone, PartialEq)]
pub struct Groups {
    rows: Vec<Vec<usize>>,
    nrows: usize,
}

impl Groups {
    pub fn len(&self) -> usize {
        self.rows.len()
    }

    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// Row indices of each group, groups in DataFrames order, rows ascending.
    pub fn iter(&self) -> impl Iterator<Item = &[usize]> + '_ {
        self.rows.iter().map(Vec::as_slice)
    }

    /// Group index of every row.
    pub fn group_of_row(&self) -> Vec<usize> {
        let mut out = vec![0; self.nrows];
        for (group, rows) in self.rows.iter().enumerate() {
            for &row in rows {
                out[row] = group;
            }
        }
        out
    }

    /// Compute one value per group and spread it back over the rows
    /// (DataFrames `transform(groupby(...), ...)` with a scalar result).
    pub fn broadcast<T: Clone>(&self, mut per_group: impl FnMut(&[usize]) -> T) -> Vec<T> {
        let mut out: Vec<Option<T>> = vec![None; self.nrows];
        for rows in &self.rows {
            let value = per_group(rows);
            for &row in rows {
                out[row] = Some(value.clone());
            }
        }
        out.into_iter()
            .map(|v| v.expect("every row belongs to exactly one group"))
            .collect()
    }

    /// Compute one value per row, group by group, written back in row order
    /// (DataFrames `transform` with a vector result).
    pub fn map_rows<T: Clone + Default>(
        &self,
        mut per_group: impl FnMut(&[usize]) -> Vec<T>,
    ) -> Vec<T> {
        let mut out = vec![T::default(); self.nrows];
        for rows in &self.rows {
            let values = per_group(rows);
            debug_assert_eq!(values.len(), rows.len());
            for (&row, value) in rows.iter().zip(values) {
                out[row] = value;
            }
        }
        out
    }

    /// First row of every group (the group's representative).
    pub fn first_rows(&self) -> Vec<usize> {
        self.rows.iter().map(|rows| rows[0]).collect()
    }
}

impl Frame {
    /// `groupby(df, cols)` with the default `sort=nothing`.
    pub fn group_by<S: AsRef<str>>(&self, names: &[S]) -> Result<Groups> {
        let indices = self.indices_of(names)?;
        let nrows = self.nrows();

        // Assign dense ids in order of first appearance, one column at a time.
        let mut ids = vec![0u32; nrows];
        let mut ngroups = usize::from(nrows > 0);
        for &c in &indices {
            let column = self.column_at(c);
            let mut seen: HashMap<(u32, Key<'_>), u32> = HashMap::new();
            for (row, id) in ids.iter_mut().enumerate() {
                let next = seen.len() as u32;
                *id = *seen.entry((*id, column.get(row).key())).or_insert(next);
            }
            ngroups = seen.len();
        }

        let mut rows: Vec<Vec<usize>> = vec![Vec::new(); ngroups];
        for (row, &id) in ids.iter().enumerate() {
            rows[id as usize].push(row);
        }

        if let Some(ranks) = self.refpool_ranks(&indices) {
            rows.sort_by(|a, b| {
                ranks
                    .iter()
                    .map(|rank| rank[a[0]].cmp(&rank[b[0]]))
                    .find(|o| *o != Ordering::Equal)
                    .unwrap_or(Ordering::Equal)
            });
        }
        Ok(Groups { rows, nrows })
    }

    /// Per-column reference-pool rank of every row, when DataFrames would
    /// take the reference-pool path for these grouping columns; `None`
    /// when it falls back to hashing (order of first appearance).
    fn refpool_ranks(&self, indices: &[usize]) -> Option<Vec<Vec<i128>>> {
        if indices.is_empty() {
            return None;
        }
        let nrows = self.nrows() as i128;
        let mut combinations: i128 = 1;
        let mut ranks = Vec::with_capacity(indices.len());
        for &c in indices {
            let (pool_len, rank) = match self.pool_at(c) {
                Some(pool) => pooled_ranks(self.column_at(c), pool)?,
                None => integer_ranks(self.column_at(c))?,
            };
            combinations = combinations.saturating_mul(pool_len);
            ranks.push(rank);
        }
        (combinations <= 2 * nrows).then_some(ranks)
    }
}

/// Pool length and per-row pool position of a pooled text column.
fn pooled_ranks(column: &Column, pool: &[Option<String>]) -> Option<(i128, Vec<i128>)> {
    let position: HashMap<Option<&str>, i128> = pool
        .iter()
        .enumerate()
        .map(|(i, level)| (level.as_deref(), i as i128))
        .collect();
    let cells = column.as_strs()?;
    let rank = cells
        .iter()
        .map(|cell| position.get(&cell.as_deref()).copied())
        .collect::<Option<Vec<i128>>>()?;
    Some((pool.len() as i128, rank))
}

/// Pool length and per-row rank of an integer-like column, or `None` when
/// DataFrames would not build an integer reference pool for it.
fn integer_ranks(column: &Column) -> Option<(i128, Vec<i128>)> {
    let len = column.len() as i128;
    if len == 0 {
        return None;
    }
    let mut values: Vec<Option<i128>> = Vec::with_capacity(column.len());
    for value in column.iter() {
        values.push(match value {
            Value::Missing => None,
            Value::Int(v) => Some(v as i128),
            Value::Bool(v) => Some(v as i128),
            Value::Float(v) => {
                let integral = v.is_finite() && v.fract() == 0.0;
                let negative_zero = v == 0.0 && v.is_sign_negative();
                if !integral || negative_zero || v.abs() >= 9.0e18 {
                    return None;
                }
                Some(v as i128)
            }
            Value::Str(_) => return None,
        });
    }
    let min = values.iter().flatten().min().copied()?;
    let max = values.iter().flatten().max().copied()?;
    let ngroups = max - min + 1;
    if ngroups >= 2 * len {
        return None;
    }
    // Missing values occupy the slot after the largest value.
    let rank = values
        .iter()
        .map(|v| v.map_or(ngroups, |v| v - min))
        .collect();
    Some((ngroups + i128::from(column.has_missing()), rank))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn order(frame: &Frame, cols: &[&str]) -> Vec<usize> {
        frame.group_by(cols).unwrap().first_rows()
    }

    // Expected orders below were observed with DataFrames 1.8.3.
    #[test]
    fn strings_group_in_order_of_appearance() {
        let f =
            Frame::from_columns([("s", Column::from_strings(["b", "a", "b", "c", "a"]))]).unwrap();
        assert_eq!(order(&f, &["s"]), vec![0, 1, 3]);
    }

    #[test]
    fn compact_integers_group_sorted() {
        let f = Frame::from_columns([("n", Column::from_ints([3, 1, 3, 2, 1]))]).unwrap();
        assert_eq!(order(&f, &["n"]), vec![1, 3, 0]);
    }

    #[test]
    fn sparse_integers_group_in_order_of_appearance() {
        let f = Frame::from_columns([("n", Column::from_ints([100, 1, 50, 1]))]).unwrap();
        assert_eq!(order(&f, &["n"]), vec![0, 1, 2]);
    }

    #[test]
    fn mixed_string_and_integer_uses_appearance_order() {
        let f = Frame::from_columns([
            ("s", Column::from_strings(["b", "a", "b", "c", "a"])),
            ("n", Column::from_ints([3, 1, 3, 2, 1])),
        ])
        .unwrap();
        assert_eq!(order(&f, &["s", "n"]), vec![0, 1, 3]);
    }

    #[test]
    fn integral_floats_and_bools_group_sorted() {
        let f = Frame::from_columns([("f", Column::from_floats([2.0, 1.0, 2.0, 3.0]))]).unwrap();
        assert_eq!(order(&f, &["f"]), vec![1, 0, 3]);
        let g = Frame::from_columns([("f", Column::from_floats([2.0, 1.0, 2.0, 1.5]))]).unwrap();
        assert_eq!(order(&g, &["f"]), vec![0, 1, 3]);
        let b = Frame::from_columns([("b", Column::from_bools([true, false, true]))]).unwrap();
        assert_eq!(order(&b, &["b"]), vec![1, 0]);
    }

    #[test]
    fn missing_integers_sort_last() {
        let f = Frame::from_columns([("n", Column::Int(vec![Some(3), None, Some(1), Some(3)]))])
            .unwrap();
        assert_eq!(order(&f, &["n"]), vec![2, 0, 1]);
    }

    #[test]
    fn pooled_text_columns_group_in_pool_order() {
        use std::sync::Arc;
        let mut f = Frame::from_columns([
            ("s", Column::from_strings(["b", "a", "b", "c", "a", "c"])),
            ("t", Column::from_strings(["y", "x", "x", "y", "y", "x"])),
        ])
        .unwrap();
        // Unpooled: order of first appearance.
        assert_eq!(order(&f, &["s", "t"]), vec![0, 1, 2, 3, 4, 5]);
        let levels = |values: &[&str]| {
            Arc::new(
                values
                    .iter()
                    .map(|v| Some(v.to_string()))
                    .collect::<Vec<_>>(),
            )
        };
        f.set_pool("s", levels(&["b", "a", "c"])).unwrap();
        f.set_pool("t", levels(&["y", "x"])).unwrap();
        // Pooled: (b,y) (b,x) (a,y) (a,x) (c,y) (c,x), as DataFrames 1.8.3 does.
        assert_eq!(order(&f, &["s", "t"]), vec![0, 2, 4, 1, 3, 5]);
        // Too many empty combinations: back to order of appearance.
        f.set_pool("s", levels(&["b", "a", "c", "d", "e", "f", "g"]))
            .unwrap();
        assert_eq!(order(&f, &["s", "t"]), vec![0, 1, 2, 3, 4, 5]);
    }

    #[test]
    fn broadcast_spreads_group_value() {
        let f = Frame::from_columns([("s", Column::from_strings(["b", "a", "b"]))]).unwrap();
        let sizes = f.group_by(&["s"]).unwrap().broadcast(|rows| rows.len());
        assert_eq!(sizes, vec![2, 1, 2]);
    }
}
