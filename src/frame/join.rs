//! Equi-joins with DataFrames.jl's result layout and row order.
//!
//! Result rows are: matched rows, then left-only rows, then right-only rows.
//! The order *within* the matched block follows DataFrames' algorithm choice
//! (`find_inner_rows`): a merge when both key columns are already sorted,
//! otherwise the longer table is scanned against a hash of the shorter one.

use std::cmp::Ordering;
use std::collections::HashMap;

use super::table::Frame;
use super::value::Key;
use crate::error::{Error, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JoinKind {
    Inner,
    Left,
    Outer,
}

impl JoinKind {
    pub fn label(self) -> &'static str {
        match self {
            JoinKind::Inner => "inner",
            JoinKind::Left => "left",
            JoinKind::Outer => "outer",
        }
    }
}

/// Join `left` and `right` on pairwise-matched key columns.
///
/// Key columns keep the left names; clashing non-key column names from the
/// right table get a `_1`, `_2`, ... suffix (`makeunique=true`).
pub fn join(
    left: &Frame,
    right: &Frame,
    left_on: &[String],
    right_on: &[String],
    kind: JoinKind,
) -> Result<Frame> {
    if left_on.len() != right_on.len() || left_on.is_empty() {
        return Err(Error::invalid(
            "Join needs the same, non-zero number of key columns on both sides",
        ));
    }
    let left_keys = left.indices_of(left_on)?;
    let right_keys = right.indices_of(right_on)?;
    check_keys(left, &left_keys, left_on)?;
    check_keys(right, &right_keys, right_on)?;
    for (&l, &r) in left_keys.iter().zip(&right_keys) {
        if left.column_at(l).is_numeric() != right.column_at(r).is_numeric() {
            return Err(Error::invalid(format!(
                "Join key '{}' is {} on the left but {} on the right",
                left.names()[l],
                left.column_at(l).type_name(),
                right.column_at(r).type_name()
            )));
        }
    }

    let (left_matched, right_matched) = inner_rows(left, right, &left_keys, &right_keys);

    let left_only: Vec<usize> = match kind {
        JoinKind::Inner => Vec::new(),
        _ => unmatched(&left_matched, left.nrows()),
    };
    let right_only: Vec<usize> = match kind {
        JoinKind::Outer => unmatched(&right_matched, right.nrows()),
        _ => Vec::new(),
    };

    let left_rows: Vec<Option<usize>> = left_matched
        .iter()
        .chain(&left_only)
        .map(|&r| Some(r))
        .chain(right_only.iter().map(|_| None))
        .collect();
    let right_rows: Vec<Option<usize>> = right_matched
        .iter()
        .map(|&r| Some(r))
        .chain(left_only.iter().map(|_| None))
        .chain(right_only.iter().map(|&r| Some(r)))
        .collect();

    let mut out = Frame::new();
    for (c, name) in left.names().iter().enumerate() {
        let column = match left_keys.iter().position(|&k| k == c) {
            // Key columns are filled from the right table for right-only rows.
            Some(k) => {
                let from_left: Vec<usize> = left_rows.iter().flatten().copied().collect();
                left.column_at(c)
                    .take(&from_left)
                    .append(right.column_at(right_keys[k]).take(&right_only))
            }
            None => left.column_at(c).take_optional(&left_rows),
        };
        out.push_column(name.clone(), column)?;
    }
    for (c, name) in right.names().iter().enumerate() {
        if right_keys.contains(&c) {
            continue;
        }
        let unique = make_unique(&out, name);
        out.push_column(unique, right.column_at(c).take_optional(&right_rows))?;
    }
    Ok(out)
}

fn check_keys(frame: &Frame, keys: &[usize], names: &[String]) -> Result<()> {
    for (&k, name) in keys.iter().zip(names) {
        if frame.column_at(k).has_missing() {
            return Err(Error::invalid(format!(
                "Missing values in key column '{name}' are not allowed in joins"
            )));
        }
    }
    Ok(())
}

fn make_unique(frame: &Frame, name: &str) -> String {
    if !frame.has_column(name) {
        return name.to_string();
    }
    (1..)
        .map(|i| format!("{name}_{i}"))
        .find(|candidate| !frame.has_column(candidate))
        .expect("an unused suffix always exists")
}

fn unmatched(matched: &[usize], nrows: usize) -> Vec<usize> {
    let mut seen = vec![false; nrows];
    for &row in matched {
        seen[row] = true;
    }
    (0..nrows).filter(|&r| !seen[r]).collect()
}

fn cmp_keys(
    a: &Frame,
    a_keys: &[usize],
    a_row: usize,
    b: &Frame,
    b_keys: &[usize],
    b_row: usize,
) -> Ordering {
    a_keys
        .iter()
        .zip(b_keys)
        .map(|(&ac, &bc)| {
            a.column_at(ac)
                .get(a_row)
                .total_cmp(&b.column_at(bc).get(b_row))
        })
        .find(|o| *o != Ordering::Equal)
        .unwrap_or(Ordering::Equal)
}

fn is_sorted(frame: &Frame, keys: &[usize]) -> bool {
    (1..frame.nrows()).all(|r| cmp_keys(frame, keys, r - 1, frame, keys, r) != Ordering::Greater)
}

/// Matching `(left_row, right_row)` pairs in DataFrames order.
fn inner_rows(
    left: &Frame,
    right: &Frame,
    left_keys: &[usize],
    right_keys: &[usize],
) -> (Vec<usize>, Vec<usize>) {
    if left.nrows() == 0 || right.nrows() == 0 {
        return (Vec::new(), Vec::new());
    }
    if is_sorted(left, left_keys) && is_sorted(right, right_keys) {
        return merge_sorted(left, right, left_keys, right_keys);
    }
    if right.nrows() < left.nrows() {
        hash_probe(left, left_keys, right, right_keys)
    } else {
        let (right_rows, left_rows) = hash_probe(right, right_keys, left, left_keys);
        (left_rows, right_rows)
    }
}

/// Scan `probe` in row order against a hash of `build`; for each probe row
/// the matching build rows are emitted in ascending order.
fn hash_probe(
    probe: &Frame,
    probe_keys: &[usize],
    build: &Frame,
    build_keys: &[usize],
) -> (Vec<usize>, Vec<usize>) {
    let mut index: HashMap<Vec<Key<'_>>, Vec<usize>> = HashMap::with_capacity(build.nrows());
    for row in 0..build.nrows() {
        index
            .entry(build.row_key(row, build_keys))
            .or_default()
            .push(row);
    }
    let mut probe_rows = Vec::new();
    let mut build_rows = Vec::new();
    for row in 0..probe.nrows() {
        if let Some(matches) = index.get(&probe.row_key(row, probe_keys)) {
            for &m in matches {
                probe_rows.push(row);
                build_rows.push(m);
            }
        }
    }
    (probe_rows, build_rows)
}

fn run_end(frame: &Frame, keys: &[usize], start: usize) -> usize {
    let mut end = start + 1;
    while end < frame.nrows() && cmp_keys(frame, keys, start, frame, keys, end) == Ordering::Equal {
        end += 1;
    }
    end
}

fn merge_sorted(
    left: &Frame,
    right: &Frame,
    left_keys: &[usize],
    right_keys: &[usize],
) -> (Vec<usize>, Vec<usize>) {
    let (mut l, mut r) = (0, 0);
    let mut left_rows = Vec::new();
    let mut right_rows = Vec::new();
    while l < left.nrows() && r < right.nrows() {
        let l_end = run_end(left, left_keys, l);
        let r_end = run_end(right, right_keys, r);
        match cmp_keys(left, left_keys, l, right, right_keys, r) {
            Ordering::Equal => {
                for right_row in r..r_end {
                    for left_row in l..l_end {
                        left_rows.push(left_row);
                        right_rows.push(right_row);
                    }
                }
                l = l_end;
                r = r_end;
            }
            Ordering::Less => l = l_end,
            Ordering::Greater => r = r_end,
        }
    }
    (left_rows, right_rows)
}

/// Rows of `left` whose key occurs in `right`, in `left` order (semi-join).
pub fn semi_join_rows(left: &Frame, right: &Frame, on: &[String]) -> Result<Vec<usize>> {
    let left_keys = left.indices_of(on)?;
    let right_keys = right.indices_of(on)?;
    let present: std::collections::HashSet<Vec<Key<'_>>> = (0..right.nrows())
        .map(|r| right.row_key(r, &right_keys))
        .collect();
    Ok((0..left.nrows())
        .filter(|&r| present.contains(&left.row_key(r, &left_keys)))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frame::Column;

    fn keys(frame: &Frame, column: &str) -> Vec<String> {
        frame
            .column(column)
            .unwrap()
            .iter()
            .map(|v| v.julia_string())
            .collect()
    }

    fn on() -> Vec<String> {
        vec!["key".to_string()]
    }

    // Expected layouts below were observed with DataFrames 1.8.3.
    #[test]
    fn left_join_puts_unmatched_rows_last() {
        let l = Frame::from_columns([
            ("key", Column::from_strings(["a", "b", "c", "b"])),
            ("val", Column::from_ints([1, 2, 3, 4])),
        ])
        .unwrap();
        let r = Frame::from_columns([
            ("key", Column::from_strings(["c", "b", "d"])),
            ("score", Column::from_ints([10, 20, 30])),
        ])
        .unwrap();
        let out = join(&l, &r, &on(), &on(), JoinKind::Left).unwrap();
        assert_eq!(keys(&out, "key"), ["b", "c", "b", "a"]);
        assert_eq!(keys(&out, "score"), ["20", "10", "20", "missing"]);
        let outer = join(&l, &r, &on(), &on(), JoinKind::Outer).unwrap();
        assert_eq!(keys(&outer, "key"), ["b", "c", "b", "a", "d"]);
        assert_eq!(keys(&outer, "val"), ["2", "3", "4", "1", "missing"]);
    }

    #[test]
    fn shorter_left_table_yields_right_order_and_unique_names() {
        let l = Frame::from_columns([
            ("key", Column::from_strings(["a", "b"])),
            ("val", Column::from_ints([1, 2])),
        ])
        .unwrap();
        let r = Frame::from_columns([
            ("key", Column::from_strings(["c", "b", "d", "b", "a"])),
            ("val", Column::from_ints([10, 20, 30, 40, 50])),
        ])
        .unwrap();
        let out = join(&l, &r, &on(), &on(), JoinKind::Outer).unwrap();
        assert_eq!(out.names(), ["key", "val", "val_1"]);
        assert_eq!(keys(&out, "key"), ["b", "b", "a", "c", "d"]);
        assert_eq!(keys(&out, "val_1"), ["20", "40", "50", "10", "30"]);
    }

    #[test]
    fn missing_keys_are_rejected() {
        let l = Frame::from_columns([("key", Column::Str(vec![None]))]).unwrap();
        let r = Frame::from_columns([("key", Column::from_strings(["a"]))]).unwrap();
        assert!(join(&l, &r, &on(), &on(), JoinKind::Left).is_err());
    }
}
