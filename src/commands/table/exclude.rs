//! `table exclude`: drop rows whose sequence contains, or is contained in,
//! a sequence of a reference FASTA.

use std::collections::HashSet;
use std::path::Path;

use memchr::memmem;

use super::common::save_gz;
use crate::error::{Error, Result};
use crate::info;
use crate::io::tsv::read_tsv;
use crate::seq::load_database;

pub fn run(
    input: &Path,
    output: &Path,
    fasta: &Path,
    name_column: &str,
    sequence_column: &str,
) -> Result<()> {
    info!("Exclude");
    let db = load_database(fasta)?;
    let table = read_tsv(input)?;
    let names: Vec<String> = table
        .column(name_column)?
        .iter()
        .map(|v| v.julia_string())
        .collect();
    let sequences: Vec<String> = table
        .column(sequence_column)?
        .iter()
        .map(|v| {
            if v.is_missing() {
                Err(Error::invalid(format!(
                    "Column '{sequence_column}' has an empty value"
                )))
            } else {
                Ok(v.julia_string())
            }
        })
        .collect::<Result<_>>()?;

    let mut discard: HashSet<&str> = HashSet::new();
    for (name, sequence) in names.iter().zip(&sequences) {
        for allele in &db {
            if memmem::find(allele.sequence.as_bytes(), sequence.as_bytes()).is_some() {
                info!("Allele {name} is a substring of {}", allele.name);
                discard.insert(sequence);
            }
            if memmem::find(sequence.as_bytes(), allele.sequence.as_bytes()).is_some() {
                info!("Allele {} is a substring of {name}", allele.name);
                discard.insert(sequence);
            }
        }
    }
    for sequence in &discard {
        info!("Discarding sequences matching: {sequence}");
    }
    let keep: Vec<bool> = sequences
        .iter()
        .map(|s| !discard.contains(s.as_str()))
        .collect();
    save_gz(output, &table.filter(&keep))?;
    Ok(())
}
