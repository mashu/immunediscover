//! `analyze haplotype`: approximate per-gene genotypes for every donor
//! under a diploid assumption.

pub mod flanks;
pub mod genotype;

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use self::flanks::{extended_sequence, FlankColumns};
use self::genotype::{call_genotype, GenotypeCall, Zygosity};
use crate::error::{Error, Result};
use crate::frame::{format_float, round_digits, Column, Frame};
use crate::info;
use crate::io::fasta::read_fasta;
use crate::io::tsv::{read_tsv, write_tsv};

#[derive(Debug, Clone)]
pub struct HaplotypeParams {
    pub input: PathBuf,
    pub output: PathBuf,
    pub case_column: String,
    pub allele_column: String,
    pub gene_column: String,
    /// Minimum count of an allele observation.
    pub mincount: i64,
    /// Minimum count relative to the major allele to be called.
    pub min_ratio: f64,
    /// FASTA whose identifiers are flagged as novel in the output.
    pub novel_fasta: Option<PathBuf>,
}

fn load_novel_names(path: &Path) -> Result<HashSet<String>> {
    if !path.is_file() {
        return Ok(HashSet::new());
    }
    Ok(read_fasta(path)?
        .iter()
        .map(|r| r.identifier().to_string())
        .collect())
}

fn counts_of(table: &Frame, column: &str) -> Result<Vec<i64>> {
    table
        .ints(column)?
        .iter()
        .map(|c| c.ok_or_else(|| Error::invalid(format!("Column '{column}' has an empty value"))))
        .collect()
}

/// Sum of `counts` per distinct (allele, extended sequence), in order of
/// first appearance among `rows`.
fn variant_totals<'a>(
    rows: &[usize],
    alleles: &'a [String],
    extended: &'a [String],
    counts: &[i64],
) -> Vec<((&'a str, &'a str), i64)> {
    let mut index: HashMap<(&str, &str), usize> = HashMap::new();
    let mut totals: Vec<((&str, &str), i64)> = Vec::new();
    for &row in rows {
        let key = (alleles[row].as_str(), extended[row].as_str());
        match index.get(&key) {
            Some(&i) => totals[i].1 += counts[row],
            None => {
                index.insert(key, totals.len());
                totals.push((key, counts[row]));
            }
        }
    }
    totals
}

/// Stable display labels `<allele>_<rank>` for every (allele, extended
/// sequence) variant, ranked by total count across all donors (ties broken
/// by extended sequence).
fn variant_labels<'a>(
    alleles: &'a [String],
    extended: &'a [String],
    counts: &[i64],
) -> HashMap<(&'a str, &'a str), String> {
    let all_rows: Vec<usize> = (0..alleles.len()).collect();
    let mut by_allele: HashMap<&str, Vec<(&str, i64)>> = HashMap::new();
    for ((allele, ext), total) in variant_totals(&all_rows, alleles, extended, counts) {
        by_allele.entry(allele).or_default().push((ext, total));
    }
    let mut labels = HashMap::new();
    for (allele, mut variants) in by_allele {
        variants.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(b.0)));
        for (rank, (ext, _)) in variants.into_iter().enumerate() {
            labels.insert((allele, ext), format!("{allele}_{}", rank + 1));
        }
    }
    labels
}

/// Infer genotypes for every donor and gene of a search result table.
pub fn infer_haplotypes(params: &HaplotypeParams) -> Result<Frame> {
    let table = read_tsv(&params.input)?;
    info!("Loaded {} rows from input file", table.nrows());

    let novel = match &params.novel_fasta {
        Some(path) => {
            let names = load_novel_names(path)?;
            info!("Loaded {} novel alleles from FASTA", names.len());
            Some(names)
        }
        None => None,
    };

    if !table.has_column("count") {
        return Err(Error::invalid(
            "Input file must have columns: count. This command expects output from exact search or similar commands.",
        ));
    }
    info!("Using existing count and ratio columns from exact search output");
    let count_column = if table.has_column("full_count") {
        info!("Using full_count as allele count source");
        "full_count"
    } else {
        info!("Using count as allele count source");
        "count"
    };

    let all_counts = counts_of(&table, count_column)?;
    let enough: Vec<bool> = all_counts.iter().map(|c| *c >= params.mincount).collect();
    let filtered = table.filter(&enough);
    info!(
        "Filtered to {} allele observations meeting mincount >= {}",
        filtered.nrows(),
        params.mincount
    );

    let present = FlankColumns::detect(&filtered);
    let extended: Vec<String> = (0..filtered.nrows())
        .map(|row| extended_sequence(&filtered, row, &params.gene_column, &present))
        .collect::<Result<_>>()?;
    let alleles: Vec<String> = filtered
        .column(&params.allele_column)?
        .iter()
        .map(|v| v.julia_string())
        .collect();
    let counts = counts_of(&filtered, count_column)?;
    let labels = variant_labels(&alleles, &extended, &counts);

    let groups = filtered.group_by(&[params.case_column.as_str(), params.gene_column.as_str()])?;
    let mut cases = Vec::new();
    let mut genes = Vec::new();
    let mut calls: Vec<GenotypeCall> = Vec::new();
    let mut novel_flags: Vec<(bool, bool)> = Vec::new();
    for rows in groups.iter() {
        cases.push(filtered.value(rows[0], &params.case_column)?.julia_string());
        genes.push(filtered.value(rows[0], &params.gene_column)?.julia_string());
        let totals = variant_totals(rows, &alleles, &extended, &counts);
        let labelled: Vec<(String, i64)> = totals
            .iter()
            .map(|(key, total)| {
                (
                    labels
                        .get(key)
                        .cloned()
                        .unwrap_or_else(|| key.0.to_string()),
                    *total,
                )
            })
            .collect();
        let call = call_genotype(&labelled, params.min_ratio);
        let is_novel = |position: Option<usize>| -> bool {
            match (&novel, position) {
                (Some(names), Some(i)) => names.contains(totals[i].0 .0),
                _ => false,
            }
        };
        novel_flags.push((is_novel(call.top[0]), is_novel(call.top[1])));
        calls.push(call);
    }

    let results = Frame::from_columns([
        ("case", Column::from_strings(cases)),
        ("gene", Column::from_strings(genes)),
        (
            "genotype",
            Column::from_strings(calls.iter().map(|c| c.zygosity.as_str())),
        ),
        (
            "allele_1",
            Column::from_strings(calls.iter().map(|c| c.allele_1.as_str())),
        ),
        (
            "allele_2",
            Column::from_strings(calls.iter().map(|c| c.allele_2.as_str())),
        ),
        (
            "count_1",
            Column::from_ints(calls.iter().map(|c| c.count_1)),
        ),
        (
            "count_2",
            Column::from_ints(calls.iter().map(|c| c.count_2)),
        ),
        (
            "freq_1",
            Column::from_floats(calls.iter().map(|c| c.freq_1)),
        ),
        (
            "freq_2",
            Column::from_floats(calls.iter().map(|c| c.freq_2)),
        ),
        ("ratio", Column::from_floats(calls.iter().map(|c| c.ratio))),
        (
            "total_count",
            Column::from_ints(calls.iter().map(|c| c.total_count)),
        ),
        (
            "other_alleles",
            Column::from_strings(calls.iter().map(|c| c.other_alleles.as_str())),
        ),
        (
            "novel_1",
            Column::from_bools(novel_flags.iter().map(|f| f.0)),
        ),
        (
            "novel_2",
            Column::from_bools(novel_flags.iter().map(|f| f.1)),
        ),
    ])?;
    write_tsv(&params.output, &results)?;

    let total = calls.len();
    info!("Haplotype inference completed:");
    info!("  Total case-gene combinations: {total}");
    if total > 0 {
        let share = |zygosity: Option<Zygosity>| -> String {
            let n = calls
                .iter()
                .filter(|c| Some(c.zygosity) == zygosity)
                .count();
            format!(
                "{n} ({}%)",
                format_float(round_digits(100.0 * n as f64 / total as f64, 1))
            )
        };
        info!("  Homozygous: {}", share(Some(Zygosity::Homozygous)));
        info!("  Heterozygous: {}", share(Some(Zygosity::Heterozygous)));
        info!("  Duplication: {}", share(Some(Zygosity::Duplication)));
        info!("  Uncertain: {}", share(None));
    }
    info!("  Results saved to: {}", params.output.display());
    Ok(results)
}

/// Run the command.
pub fn run(params: &HaplotypeParams) -> Result<()> {
    info!("Haplotype inference");
    infer_haplotypes(params).map(|_| ())
}
