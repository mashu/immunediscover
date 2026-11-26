//! `fasta diff` and `fasta hash`.

use std::path::{Path, PathBuf};

use super::keyed_set::KeyedSet;
use crate::error::Result;
use crate::info;
use crate::seq::{load_database, unique_name};

/// Load a FASTA as a sequence-keyed set of names.
pub fn sequence_set(path: &Path) -> Result<KeyedSet> {
    Ok(KeyedSet::from_pairs(
        load_database(path)?
            .into_iter()
            .map(|a| (a.sequence, a.name)),
    ))
}

/// Compare every pair of FASTA files by sequence content and print the
/// sizes of union, intersection and both differences (with allele names).
pub fn diff(paths: &[PathBuf]) -> Result<()> {
    let sets: Vec<KeyedSet> = paths
        .iter()
        .map(|p| sequence_set(p))
        .collect::<Result<_>>()?;
    for i in 0..sets.len() {
        for j in (i + 1)..sets.len() {
            println!("Comparing {} vs {}", paths[i].display(), paths[j].display());
            println!("Union: {}", sets[i].union(&sets[j]).len());
            println!("Intersection: {}", sets[i].intersection(&sets[j]).len());
            for (a, b) in [(i, j), (j, i)] {
                let only = sets[a].difference(&sets[b]);
                let names: Vec<&str> = only.iter().map(|(_, name)| name).collect();
                println!("Difference ({}): \n{}", names.len(), names.join("\n"));
            }
        }
    }
    Ok(())
}

/// Print the FASTA with every name replaced by its content-derived name.
pub fn hash(path: &Path) -> Result<()> {
    info!("Hashing alleles");
    for allele in load_database(path)? {
        println!(
            ">{}\n{}",
            unique_name(&allele.name, &allele.sequence),
            allele.sequence
        );
    }
    Ok(())
}
