//! Stable multi-column sorting (`sort(df, cols; rev)`).

use std::cmp::Ordering;

use super::table::Frame;
use crate::error::Result;

/// One sort criterion.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SortKey {
    pub column: String,
    pub descending: bool,
}

impl SortKey {
    pub fn asc(column: impl Into<String>) -> Self {
        SortKey {
            column: column.into(),
            descending: false,
        }
    }

    pub fn desc(column: impl Into<String>) -> Self {
        SortKey {
            column: column.into(),
            descending: true,
        }
    }
}

impl Frame {
    /// Permutation that sorts the rows; ties keep their original order.
    pub fn sort_permutation(&self, keys: &[SortKey]) -> Result<Vec<usize>> {
        let resolved: Vec<(usize, bool)> = keys
            .iter()
            .map(|k| Ok((self.index_of(&k.column)?, k.descending)))
            .collect::<Result<_>>()?;
        let mut rows: Vec<usize> = (0..self.nrows()).collect();
        rows.sort_by(|&a, &b| {
            for &(c, descending) in &resolved {
                let column = self.column_at(c);
                let ordering = column.get(a).total_cmp(&column.get(b));
                let ordering = if descending {
                    ordering.reverse()
                } else {
                    ordering
                };
                if ordering != Ordering::Equal {
                    return ordering;
                }
            }
            Ordering::Equal
        });
        Ok(rows)
    }

    /// Sorted copy.
    pub fn sorted(&self, keys: &[SortKey]) -> Result<Frame> {
        Ok(self.take(&self.sort_permutation(keys)?))
    }

    /// Ascending sort on the given columns.
    pub fn sorted_asc<S: AsRef<str>>(&self, columns: &[S]) -> Result<Frame> {
        let keys: Vec<SortKey> = columns.iter().map(|c| SortKey::asc(c.as_ref())).collect();
        self.sorted(&keys)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frame::Column;

    #[test]
    fn missing_sorts_last_ascending_and_first_descending() {
        let f = Frame::from_columns([
            ("a", Column::Int(vec![Some(2), None, Some(1)])),
            ("b", Column::from_strings(["x", "y", "z"])),
        ])
        .unwrap();
        let asc = f.sorted(&[SortKey::asc("a")]).unwrap();
        assert_eq!(
            asc.column("a").unwrap(),
            &Column::Int(vec![Some(1), Some(2), None])
        );
        let desc = f.sorted(&[SortKey::desc("a")]).unwrap();
        assert_eq!(
            desc.column("a").unwrap(),
            &Column::Int(vec![None, Some(2), Some(1)])
        );
    }

    #[test]
    fn sort_is_stable() {
        let f = Frame::from_columns([
            ("k", Column::from_ints([1, 0, 1, 0])),
            ("v", Column::from_strings(["a", "b", "c", "d"])),
        ])
        .unwrap();
        let out = f.sorted(&[SortKey::desc("k")]).unwrap();
        let v: Vec<_> = out.strs("v").unwrap().iter().flatten().cloned().collect();
        assert_eq!(v, ["a", "c", "b", "d"]);
    }
}
