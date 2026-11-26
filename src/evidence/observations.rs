//! Building observations from the cluster table of `discover blast`
//! (the `*.unfiltered.tsv.gz` written with `-v`).

use std::collections::HashMap;
use std::path::Path;

use super::{Evidence, Observation};
use crate::error::Result;
use crate::frame::Frame;
use crate::io::tsv::{read_table, ReadOptions};
use crate::seq::gene_of;

/// Read a cluster table and collapse it to one observation per sequence
/// and donor. Returns the observations and the number of donors.
pub fn read_clusters(path: &Path) -> Result<(Vec<Observation>, usize)> {
    let options = ["well", "case", "sseqid", "aln_qseq"]
        .iter()
        .fold(ReadOptions::tsv(), |options, column| {
            options.text_column(column)
        });
    from_clusters(&read_table(path, &options)?)
}

/// Collapse clusters (one row per extension variant) to one observation
/// per trimmed sequence and donor. Clusters whose trimming failed carry no
/// sequence and are left out. A sequence hit under several database
/// alleles takes the best supported one as `nearest`.
pub fn from_clusters(clusters: &Frame) -> Result<(Vec<Observation>, usize)> {
    let text = |name: &str| -> Result<Vec<String>> {
        Ok(clusters
            .column(name)?
            .iter()
            .map(|v| {
                if v.is_missing() {
                    String::new()
                } else {
                    v.julia_string()
                }
            })
            .collect())
    };
    let (wells, cases, hits, sequences) = (
        text("well")?,
        text("case")?,
        text("sseqid")?,
        text("aln_qseq")?,
    );
    let full_counts = clusters.ints("full_count")?;
    let distances = clusters.ints("aln_mismatch")?;
    let full_ratios = clusters.floats("full_ratio")?;
    let ratios = clusters.floats("ratio")?;

    let mut donors: HashMap<(&str, &str), usize> = HashMap::new();
    let mut index: HashMap<(usize, &str), usize> = HashMap::new();
    let mut observations: Vec<Observation> = Vec::new();
    let mut reads_by_hit: Vec<HashMap<&str, u64>> = Vec::new();

    for row in 0..clusters.nrows() {
        let sequence = sequences[row].as_str();
        if sequence.is_empty() {
            continue;
        }
        let next = donors.len();
        let donor = *donors
            .entry((wells[row].as_str(), cases[row].as_str()))
            .or_insert(next);
        let full = full_counts[row].unwrap_or(0).max(0) as u64;
        let distance = distances[row].unwrap_or(0).max(0) as u64;
        let i = *index.entry((donor, sequence)).or_insert_with(|| {
            observations.push(Observation {
                donor,
                gene: String::new(),
                nearest: String::new(),
                sequence: sequence.to_string(),
                count: 0,
                best_full_count: 0,
                best_full_ratio: 0.0,
                ratio: 0.0,
                db_distance: distance,
                evidence: Evidence::default(),
            });
            reads_by_hit.push(HashMap::new());
            observations.len() - 1
        });
        let o = &mut observations[i];
        o.count += full;
        o.best_full_count = o.best_full_count.max(full);
        o.best_full_ratio = o.best_full_ratio.max(full_ratios[row].unwrap_or(0.0));
        o.ratio = o.ratio.max(ratios[row].unwrap_or(0.0));
        o.db_distance = o.db_distance.min(distance);
        *reads_by_hit[i].entry(hits[row].as_str()).or_insert(0) += full;
    }
    // A sequence belongs to one allele, and so to one gene, on the whole
    // plate: the one most of its reads were assigned to. Assigned donor by
    // donor, a sequence two genes share would move between them and put
    // its reads into a different gene's total in every donor.
    let mut reads_on_plate: HashMap<&str, HashMap<&str, u64>> = HashMap::new();
    for (o, by_hit) in observations.iter().zip(&reads_by_hit) {
        let of_sequence = reads_on_plate.entry(o.sequence.as_str()).or_default();
        for (hit, reads) in by_hit {
            *of_sequence.entry(hit).or_insert(0) += reads;
        }
    }
    let nearest: Vec<String> = observations
        .iter()
        .map(|o| {
            let (hit, _) = reads_on_plate[o.sequence.as_str()]
                .iter()
                .max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0)))
                .expect("every observation has at least one cluster");
            hit.to_string()
        })
        .collect();
    for (o, nearest) in observations.iter_mut().zip(nearest) {
        o.gene = gene_of(&nearest).to_string();
        o.nearest = nearest;
    }
    Ok((observations, donors.len()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frame::Column;

    #[test]
    fn variants_collapse_to_one_observation_per_sequence_and_donor() {
        let clusters = Frame::from_columns([
            ("well", Column::from_ints([1, 1, 1, 2, 1])),
            ("case", Column::from_strings(["D1", "D1", "D1", "D2", "D1"])),
            (
                "sseqid",
                Column::from_strings(["A*01", "A*01", "A*02", "A*01", "A*01"]),
            ),
            ("full_count", Column::from_ints([7, 3, 5, 4, 9])),
            ("full_ratio", Column::from_floats([1.0, 0.4, 0.7, 1.0, 1.0])),
            ("ratio", Column::from_floats([1.0, 1.0, 0.5, 1.0, 0.9])),
            (
                "aln_qseq",
                Column::from_strings(["ACGT", "ACGT", "ACGT", "ACGT", ""]),
            ),
            ("aln_mismatch", Column::from_ints([0, 0, 1, 0, -1])),
        ])
        .unwrap();
        let (observations, donors) = from_clusters(&clusters).unwrap();
        assert_eq!(donors, 2);
        assert_eq!(observations.len(), 2); // the failed cluster is left out
        let first = &observations[0];
        assert_eq!((first.count, first.best_full_count), (15, 7));
        assert_eq!(first.nearest, "A*01"); // 10 reads under A*01, 5 under A*02
        assert_eq!(first.gene, "A");
        assert!(first.is_database_allele());
        assert_eq!(observations[1].count, 4);
    }

    #[test]
    fn a_sequence_keeps_one_gene_on_the_whole_plate() {
        // Donor 2 has the sequence under gene B, donors 1 and 3 under gene
        // A, with more reads altogether.
        let clusters = Frame::from_columns([
            ("well", Column::from_ints([1, 2, 3])),
            ("case", Column::from_strings(["D1", "D2", "D3"])),
            ("sseqid", Column::from_strings(["A*01", "B*01", "A*01"])),
            ("full_count", Column::from_ints([6, 8, 5])),
            ("full_ratio", Column::from_floats([1.0, 1.0, 1.0])),
            ("ratio", Column::from_floats([1.0, 1.0, 1.0])),
            ("aln_qseq", Column::from_strings(["ACGT", "ACGT", "ACGT"])),
            ("aln_mismatch", Column::from_ints([1, 1, 1])),
        ])
        .unwrap();
        let (observations, _) = from_clusters(&clusters).unwrap();
        assert!(observations
            .iter()
            .all(|o| o.gene == "A" && o.nearest == "A*01"));
    }
}
