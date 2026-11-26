//! Grouping alleles by co-occurrence similarity.

use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};
use std::str::FromStr;

use crate::error::{Error, Result};
use crate::stats::{hclust, Linkage};

/// How alleles are grouped from the similarity matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClusterMethod {
    /// Connected components of the graph of pairs above the threshold.
    Components,
    /// Hierarchical clustering cut at the threshold.
    Hierarchical(Linkage),
}

impl ClusterMethod {
    pub const NAMES: [&'static str; 4] = ["components", "complete", "average", "single"];
}

impl FromStr for ClusterMethod {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        match s {
            "components" => Ok(ClusterMethod::Components),
            "complete" => Ok(ClusterMethod::Hierarchical(Linkage::Complete)),
            "average" => Ok(ClusterMethod::Hierarchical(Linkage::Average)),
            "single" => Ok(ClusterMethod::Hierarchical(Linkage::Single)),
            other => Err(Error::invalid(format!(
                "Unknown cluster method '{other}', expected one of: {}",
                Self::NAMES.join(", ")
            ))),
        }
    }
}

/// Connected components among items linked by `similarity >= threshold`.
/// Components are ordered by their first item; members are sorted.
pub fn components(
    similarity: &[f64],
    items: &[String],
    threshold: f64,
    min_size: usize,
) -> Vec<Vec<String>> {
    let n = items.len();
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        for j in (i + 1)..n {
            if similarity[i * n + j] >= threshold {
                neighbors[i].push(j);
                neighbors[j].push(i);
            }
        }
    }
    let mut visited = vec![false; n];
    let mut out = Vec::new();
    for start in 0..n {
        if visited[start] || neighbors[start].is_empty() {
            continue;
        }
        let mut members = Vec::new();
        let mut queue = VecDeque::from([start]);
        visited[start] = true;
        while let Some(u) = queue.pop_front() {
            members.push(items[u].clone());
            for &v in &neighbors[u] {
                if !visited[v] {
                    visited[v] = true;
                    queue.push_back(v);
                }
            }
        }
        if members.len() >= min_size {
            members.sort();
            out.push(members);
        }
    }
    out
}

/// Hierarchical clusters of `1 - similarity`, cut at `1 - threshold`.
///
/// Groups are returned in order of their first member. (Upstream numbers
/// them in the iteration order of a hash table; membership is identical.)
pub fn hierarchical(
    similarity: &[f64],
    items: &[String],
    threshold: f64,
    min_size: usize,
    linkage: Linkage,
) -> Result<Vec<Vec<String>>> {
    let n = items.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let distances: Vec<f64> = similarity.iter().map(|s| 1.0 - s).collect();
    let labels = hclust(&distances, n, linkage)?.cut_at_height(1.0 - threshold);
    let mut groups: BTreeMap<usize, Vec<String>> = BTreeMap::new();
    for (item, label) in items.iter().zip(labels) {
        groups.entry(label).or_default().push(item.clone());
    }
    Ok(groups
        .into_values()
        .filter(|g| g.len() >= min_size)
        .map(|mut g| {
            g.sort();
            g
        })
        .collect())
}

/// Connected components of an explicit edge list.
pub fn groups_from_edges(edges: &[(String, String)], min_size: usize) -> Vec<Vec<String>> {
    let mut adjacency: HashMap<&str, BTreeSet<&str>> = HashMap::new();
    let mut nodes: BTreeSet<&str> = BTreeSet::new();
    for (a, b) in edges {
        nodes.insert(a);
        nodes.insert(b);
        if a != b {
            adjacency.entry(a).or_default().insert(b);
            adjacency.entry(b).or_default().insert(a);
        }
    }
    let mut visited: BTreeSet<&str> = BTreeSet::new();
    let mut out = Vec::new();
    for &start in &nodes {
        if !visited.insert(start) {
            continue;
        }
        let mut members = vec![start.to_string()];
        let mut queue = VecDeque::from([start]);
        while let Some(u) = queue.pop_front() {
            for &v in adjacency.get(u).into_iter().flatten() {
                if visited.insert(v) {
                    members.push(v.to_string());
                    queue.push_back(v);
                }
            }
        }
        if members.len() >= min_size {
            members.sort();
            out.push(members);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn items(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("a{i}")).collect()
    }

    #[test]
    fn components_follow_threshold() {
        // 0-1 and 1-2 linked, 3 isolated, 4-5 linked.
        let n = 6;
        let mut s = vec![0.0; n * n];
        for (i, j) in [(0, 1), (1, 2), (4, 5)] {
            s[i * n + j] = 0.8;
            s[j * n + i] = 0.8;
        }
        let groups = components(&s, &items(n), 0.5, 2);
        assert_eq!(groups, vec![vec!["a0", "a1", "a2"], vec!["a4", "a5"]]);
        assert_eq!(components(&s, &items(n), 0.5, 3).len(), 1);
        assert!(components(&s, &items(n), 0.9, 2).is_empty());
    }

    #[test]
    fn hierarchical_groups_respect_linkage() {
        let n = 3;
        // 0-1 strongly similar, 1-2 similar, 0-2 dissimilar.
        let mut s = vec![0.0; n * n];
        for i in 0..n {
            s[i * n + i] = 1.0;
        }
        for (i, j, v) in [(0, 1, 0.9), (1, 2, 0.7), (0, 2, 0.1)] {
            s[i * n + j] = v;
            s[j * n + i] = v;
        }
        let single = hierarchical(&s, &items(n), 0.5, 2, Linkage::Single).unwrap();
        assert_eq!(single, vec![vec!["a0", "a1", "a2"]]);
        let complete = hierarchical(&s, &items(n), 0.5, 2, Linkage::Complete).unwrap();
        assert_eq!(complete, vec![vec!["a0", "a1"]]);
    }

    #[test]
    fn method_names_parse() {
        assert_eq!(
            "components".parse::<ClusterMethod>().unwrap(),
            ClusterMethod::Components
        );
        assert_eq!(
            "average".parse::<ClusterMethod>().unwrap(),
            ClusterMethod::Hierarchical(Linkage::Average)
        );
        assert!("ward".parse::<ClusterMethod>().is_err());
    }
}
