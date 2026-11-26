//! Agglomerative hierarchical clustering.
//!
//! A port of the parts of Clustering.jl (`hclust`, `cutree`) that the
//! co-occurrence analysis uses: single linkage via the nearest-neighbour
//! array algorithm, complete and average linkage via the nearest-neighbour
//! chain with Lance-Williams updates, R-compatible merge ordering and
//! height-based tree cutting. Tie-breaking follows the original code so
//! that flat clusters agree when distances repeat, which is the norm for
//! correlation-derived distances on few donors.
//!
//! Tree and leaf identifiers are 1-based as in the original: a negative id
//! `-k` is leaf `k`, a positive id `k` is the tree created by merge `k`.

use crate::error::{Error, Result};

/// Cluster linkage criterion.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Linkage {
    Single,
    Average,
    Complete,
}

/// A dendrogram as a sequence of merges ordered by height.
#[derive(Debug, Clone, PartialEq)]
pub struct Dendrogram {
    n_leaves: usize,
    left: Vec<i64>,
    right: Vec<i64>,
    heights: Vec<f64>,
}

/// Dense square matrix addressed with 1-based indices.
#[derive(Clone)]
struct Matrix {
    n: usize,
    data: Vec<f64>,
}

impl Matrix {
    fn get(&self, i: usize, j: usize) -> f64 {
        self.data[(i - 1) * self.n + (j - 1)]
    }

    fn set(&mut self, i: usize, j: usize, value: f64) {
        self.data[(i - 1) * self.n + (j - 1)] = value;
    }
}

struct Merges {
    left: Vec<i64>,
    right: Vec<i64>,
    heights: Vec<f64>,
}

impl Merges {
    fn new() -> Self {
        Merges {
            left: Vec::new(),
            right: Vec::new(),
            heights: Vec::new(),
        }
    }

    fn push(&mut self, left: i64, right: i64, height: f64) -> i64 {
        self.left.push(left);
        self.right.push(right);
        self.heights.push(height);
        self.heights.len() as i64
    }
}

/// How the distance to a freshly merged cluster is derived.
trait LanceWilliams {
    /// New distance between cluster `k` and the union of `i` and `j`.
    fn merged(&self, d_ki: f64, d_kj: f64, ni: usize, nj: usize) -> f64;
}

struct Minimum;
struct Maximum;
struct Average;

impl LanceWilliams for Minimum {
    fn merged(&self, d_ki: f64, d_kj: f64, _: usize, _: usize) -> f64 {
        if d_ki > d_kj {
            d_kj
        } else {
            d_ki
        }
    }
}

impl LanceWilliams for Maximum {
    fn merged(&self, d_ki: f64, d_kj: f64, _: usize, _: usize) -> f64 {
        if d_ki < d_kj {
            d_kj
        } else {
            d_ki
        }
    }
}

impl LanceWilliams for Average {
    fn merged(&self, d_ki: f64, d_kj: f64, ni: usize, nj: usize) -> f64 {
        let combined = (ni as f64 * d_ki + nj as f64 * d_kj) / (ni + nj) as f64;
        // Guard reducibility against floating-point error.
        combined.max(d_ki.min(d_kj))
    }
}

/// Update the upper triangle of `d` when cluster `j` is merged into `i`
/// (`i < j`) and the last active cluster `n_active` moves into slot `j`.
fn update_distances(
    d: &mut Matrix,
    metric: &dyn LanceWilliams,
    size: &dyn Fn(usize) -> usize,
    i: usize,
    j: usize,
    n_active: usize,
) {
    let ni = size(i);
    let nj = size(j);
    for k in 1..=i {
        let value = metric.merged(d.get(k, i), d.get(k, j), ni, nj);
        d.set(k, i, value);
    }
    for k in (i + 1)..j {
        let value = metric.merged(d.get(i, k), d.get(k, j), ni, nj);
        d.set(i, k, value);
    }
    for k in (j + 1)..=n_active {
        let value = metric.merged(d.get(i, k), d.get(j, k), ni, nj);
        d.set(i, k, value);
    }
    if j < n_active {
        d.set(j, j, d.get(n_active, n_active));
        for k in 1..j {
            d.set(k, j, d.get(k, n_active));
        }
        for k in (j + 1)..n_active {
            d.set(j, k, d.get(k, n_active));
        }
    }
}

/// Nearest active neighbour of `i` among clusters `1..=n_active`.
fn nearest_neighbor(d: &Matrix, i: usize, n_active: usize) -> (usize, f64) {
    if n_active <= 1 {
        return (0, f64::NAN);
    }
    let (mut best, mut best_dist) = if i > 1 {
        (1, d.get(1, i))
    } else {
        (2, d.get(i, 2))
    };
    for j in (best + 1)..i {
        if best_dist > d.get(j, i) {
            best = j;
            best_dist = d.get(j, i);
        }
    }
    for j in (i + 1)..=n_active {
        if best_dist > d.get(i, j) {
            best = j;
            best_dist = d.get(i, j);
        }
    }
    (best, best_dist)
}

fn single_linkage(mut d: Matrix) -> Merges {
    let n = d.n;
    let mut merges = Merges::new();
    // nn[k-1] is the nearest neighbour of active cluster k.
    let mut nn: Vec<usize> = (1..=n).map(|i| nearest_neighbor(&d, i, n).0).collect();
    let mut trees: Vec<i64> = (1..=n as i64).map(|i| -i).collect();
    let upper = |d: &Matrix, a: usize, b: usize| if a < b { d.get(a, b) } else { d.get(b, a) };

    while trees.len() > 1 {
        let mut i = 1;
        let mut min_dist = upper(&d, 1, nn[0]);
        for k in 2..=trees.len() {
            let dist = upper(&d, k, nn[k - 1]);
            if dist < min_dist {
                min_dist = dist;
                i = k;
            }
        }
        let mut j = nn[i - 1];
        if i > j {
            std::mem::swap(&mut i, &mut j);
        }
        let last = trees.len();
        update_distances(&mut d, &Minimum, &|_| 0, i, j, last);
        trees[i - 1] = merges.push(trees[i - 1], trees[j - 1], min_dist);
        trees[j - 1] = trees[last - 1];
        nn[j - 1] = nn[last - 1];
        nn.pop();
        trees.pop();
        for neighbor in nn.iter_mut() {
            if *neighbor == j {
                *neighbor = i;
            } else if *neighbor == last {
                *neighbor = j;
            }
        }
        // Recompute the neighbour of the merged cluster.
        let mut best = 0;
        let mut best_dist = f64::MAX;
        for k in 1..i {
            if best == 0 || d.get(k, i) < best_dist {
                best_dist = d.get(k, i);
                best = k;
            }
        }
        for k in (i + 1)..=trees.len() {
            if best == 0 || d.get(i, k) < best_dist {
                best_dist = d.get(i, k);
                best = k;
            }
        }
        nn[i - 1] = best;
    }
    merges
}

fn nn_chain(mut d: Matrix, metric: &dyn LanceWilliams) -> Merges {
    let n = d.n;
    let mut merges = Merges::new();
    let mut ids: Vec<i64> = (1..=n as i64).map(|i| -i).collect();
    // Number of leaves in each created tree, indexed by merge id - 1.
    let mut tree_sizes: Vec<usize> = Vec::new();
    let mut chain: Vec<usize> = vec![1];

    while ids.len() > 1 {
        let mut min_dist;
        loop {
            let top = *chain.last().expect("chain is never empty here");
            let (next, dist) = nearest_neighbor(&d, top, ids.len());
            min_dist = dist;
            if chain.len() > 1 && next == chain[chain.len() - 2] {
                break;
            }
            chain.push(next);
        }
        let mut lo = chain.pop().expect("chain holds a mutual pair");
        let mut hi = chain.pop().expect("chain holds a mutual pair");
        if lo > hi {
            std::mem::swap(&mut lo, &mut hi);
        }
        let last = ids.len();
        {
            let size = |slot: usize| -> usize {
                let id = ids[slot - 1];
                if id > 0 {
                    tree_sizes[(id - 1) as usize]
                } else {
                    1
                }
            };
            update_distances(&mut d, metric, &size, lo, hi, last);
        }
        let size_of = |id: i64, sizes: &[usize]| if id > 0 { sizes[(id - 1) as usize] } else { 1 };
        let merged_size = size_of(ids[lo - 1], &tree_sizes) + size_of(ids[hi - 1], &tree_sizes);
        let new_id = merges.push(ids[lo - 1], ids[hi - 1], min_dist);
        tree_sizes.push(merged_size);
        ids[lo - 1] = new_id;
        ids[hi - 1] = ids[last - 1];
        ids.pop();

        let mut restart = false;
        for slot in chain.iter_mut() {
            if *slot == lo || *slot == hi {
                // With duplicate distances a merged cluster can still sit in
                // the chain; restart the search in that case.
                restart = true;
                break;
            } else if *slot == last {
                *slot = hi;
            }
        }
        if restart {
            chain.clear();
        }
        if chain.is_empty() {
            chain.push(1);
        }
    }
    merges
}

/// R's ordering of the two subtrees of a merge.
fn is_r_ordered(i: i64, j: i64) -> bool {
    (i < 0 && j < 0 && i > j) || (i > 0 && j > 0 && i < j) || (i < 0 && j > 0)
}

/// Sort merges by height (stably) and renumber subtree references.
fn order_like_r(merges: &mut Merges) {
    let count = merges.heights.len();
    let mut order: Vec<usize> = (0..count).collect();
    order.sort_by(|&a, &b| merges.heights[a].total_cmp(&merges.heights[b]));
    let mut position = vec![0i64; count];
    for (new, &old) in order.iter().enumerate() {
        position[old] = new as i64 + 1;
    }
    for i in 0..count {
        if merges.left[i] > 0 {
            merges.left[i] = position[(merges.left[i] - 1) as usize];
        }
        if merges.right[i] > 0 {
            merges.right[i] = position[(merges.right[i] - 1) as usize];
        }
        if !is_r_ordered(merges.left[i], merges.right[i]) {
            std::mem::swap(&mut merges.left[i], &mut merges.right[i]);
        }
    }
    merges.left = order.iter().map(|&o| merges.left[o]).collect();
    merges.right = order.iter().map(|&o| merges.right[o]).collect();
    merges.heights = order.iter().map(|&o| merges.heights[o]).collect();
}

/// Cluster the points of a symmetric distance matrix (row-major, `n` x `n`).
pub fn hclust(distances: &[f64], n: usize, linkage: Linkage) -> Result<Dendrogram> {
    if distances.len() != n * n {
        return Err(Error::invalid("Distance matrix should be square."));
    }
    let matrix = Matrix {
        n,
        data: distances.to_vec(),
    };
    for i in 1..=n {
        for j in (i + 1)..=n {
            if matrix.get(i, j) != matrix.get(j, i) {
                return Err(Error::invalid("Distance matrix should be symmetric."));
            }
        }
    }
    let mut merges = match linkage {
        Linkage::Single => single_linkage(matrix),
        Linkage::Complete => nn_chain(matrix, &Maximum),
        Linkage::Average => nn_chain(matrix, &Average),
    };
    order_like_r(&mut merges);
    Ok(Dendrogram {
        n_leaves: n,
        left: merges.left,
        right: merges.right,
        heights: merges.heights,
    })
}

impl Dendrogram {
    pub fn heights(&self) -> &[f64] {
        &self.heights
    }

    /// Cut the tree at `height`; returns a 1-based cluster label per leaf.
    /// Labels are numbered by each cluster's smallest leaf index.
    pub fn cut_at_height(&self, height: f64) -> Vec<usize> {
        let n = self.n_leaves;
        let m = self.heights.len();
        let mut by_height: Vec<usize> = (0..m).collect();
        by_height.sort_by(|&a, &b| self.heights[a].total_cmp(&self.heights[b]));

        let mut n_merges = m;
        match by_height.iter().rposition(|&i| self.heights[i] <= height) {
            Some(last) if last + 1 < n_merges => n_merges = last + 1,
            Some(_) => {}
            None => {
                if m >= 1 {
                    n_merges = 0;
                }
            }
        }

        let mut position = vec![0usize; m];
        for (pos, &merge) in by_height.iter().enumerate() {
            position[merge] = pos;
        }
        let mut clusters: Vec<Vec<usize>> = Vec::with_capacity(n_merges);
        let mut unmerged = vec![true; n];
        for &merge in &by_height[..n_merges] {
            let mut members = Vec::new();
            for side in [self.left[merge], self.right[merge]] {
                if side < 0 {
                    let leaf = (-side) as usize;
                    unmerged[leaf - 1] = false;
                    members.push(leaf);
                } else {
                    let child = position[(side - 1) as usize];
                    members.append(&mut clusters[child]);
                }
            }
            clusters.push(members);
        }

        let mut clusters: Vec<Vec<usize>> =
            clusters.into_iter().filter(|c| !c.is_empty()).collect();
        clusters.sort_by_key(|c| *c.iter().min().expect("cluster is non-empty"));

        let mut labels = vec![0usize; n];
        let mut singles = (1..=n).filter(|&leaf| unmerged[leaf - 1]).peekable();
        let mut next = 1;
        for cluster in &clusters {
            let smallest = *cluster.iter().min().expect("cluster is non-empty");
            while let Some(&leaf) = singles.peek() {
                if leaf >= smallest {
                    break;
                }
                labels[leaf - 1] = next;
                next += 1;
                singles.next();
            }
            for &leaf in cluster {
                labels[leaf - 1] = next;
            }
            next += 1;
        }
        for leaf in singles {
            labels[leaf - 1] = next;
            next += 1;
        }
        labels
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square(rows: &[&[f64]]) -> Vec<f64> {
        rows.iter().flat_map(|r| r.iter().copied()).collect()
    }

    #[test]
    fn two_tight_pairs() {
        let d = square(&[
            &[0.0, 0.1, 0.9, 0.8],
            &[0.1, 0.0, 0.85, 0.9],
            &[0.9, 0.85, 0.0, 0.2],
            &[0.8, 0.9, 0.2, 0.0],
        ]);
        for linkage in [Linkage::Single, Linkage::Average, Linkage::Complete] {
            let tree = hclust(&d, 4, linkage).unwrap();
            assert_eq!(tree.cut_at_height(0.5), vec![1, 1, 2, 2], "{linkage:?}");
            assert_eq!(tree.cut_at_height(0.05), vec![1, 2, 3, 4], "{linkage:?}");
            assert_eq!(tree.cut_at_height(1.0), vec![1, 1, 1, 1], "{linkage:?}");
        }
    }

    #[test]
    fn asymmetric_input_is_rejected() {
        assert!(hclust(&[0.0, 1.0, 2.0, 0.0], 2, Linkage::Single).is_err());
    }

    #[test]
    fn single_point() {
        let tree = hclust(&[0.0], 1, Linkage::Complete).unwrap();
        assert_eq!(tree.cut_at_height(0.5), vec![1]);
    }
}
