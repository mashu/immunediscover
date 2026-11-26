//! Match `discover blast` allele names against KI `_S` target alleles.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::path::Path;

use crate::error::{Error, Result};
use crate::frame::{Frame, Value};
use crate::io::fasta::{read_fasta, write_fasta};
use crate::io::tsv::{read_table, ReadOptions};

/// Gene segment used to select KI targets and BLAST tables.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GeneKind {
    V,
    D,
    J,
}

impl GeneKind {
    pub fn as_str(self) -> &'static str {
        match self {
            GeneKind::V => "V",
            GeneKind::D => "D",
            GeneKind::J => "J",
        }
    }

    pub fn prefix(self) -> &'static str {
        match self {
            GeneKind::V => "IGHV",
            GeneKind::D => "IGHD",
            GeneKind::J => "IGHJ",
        }
    }

    pub fn parse(s: &str) -> Result<Self> {
        match s {
            "V" | "v" => Ok(GeneKind::V),
            "D" | "d" => Ok(GeneKind::D),
            "J" | "j" => Ok(GeneKind::J),
            other => Err(Error::invalid(format!(
                "gene kind '{other}' is not one of V, D, J"
            ))),
        }
    }
}

/// Per-gene and overall recovery / false-positive summary.
#[derive(Debug, Clone)]
pub struct RecoveryReport {
    pub gene: GeneKind,
    pub target_alleles: BTreeSet<String>,
    pub found_targets: BTreeSet<String>,
    pub missing_targets: BTreeSet<String>,
    pub all_s_alleles: BTreeSet<String>,
    /// Extra `_S` calls that did not match any KI target (by signature).
    pub false_positives: BTreeSet<String>,
    /// Target count and found count per gene name (slash-alleles may appear
    /// under more than one gene, matching the Julia checker).
    pub per_gene: BTreeMap<String, GeneStats>,
}

#[derive(Debug, Clone)]
pub struct GeneStats {
    pub total_target: usize,
    pub found: usize,
    pub found_alleles: BTreeSet<String>,
    pub missing_alleles: BTreeSet<String>,
}

impl RecoveryReport {
    pub fn total_target(&self) -> usize {
        self.per_gene.values().map(|g| g.total_target).sum()
    }

    pub fn total_found(&self) -> usize {
        self.per_gene.values().map(|g| g.found).sum()
    }

    pub fn unique_target_count(&self) -> usize {
        self.target_alleles.len()
    }

    pub fn unique_found_count(&self) -> usize {
        self.found_targets.len()
    }

    pub fn recovery_pct(&self) -> f64 {
        let total = self.total_target();
        if total == 0 {
            0.0
        } else {
            100.0 * self.total_found() as f64 / total as f64
        }
    }

    pub fn false_positive_rate(&self) -> f64 {
        let n = self.all_s_alleles.len();
        if n == 0 {
            0.0
        } else {
            100.0 * self.false_positives.len() as f64 / n as f64
        }
    }

    pub fn all_recovered(&self) -> bool {
        self.missing_targets.is_empty()
    }
}

/// Strip trailing `" Novel"` and surrounding whitespace.
pub fn normalize_allele(name: &str) -> String {
    name.replace(" Novel", "").trim().to_string()
}

/// Gene names before each `*` (handles `/`-joined alleles).
pub fn extract_genes(allele_name: &str) -> BTreeSet<String> {
    let allele_name = normalize_allele(allele_name);
    let mut genes = BTreeSet::new();
    for part in allele_name.split('/') {
        let part = part.trim();
        if let Some(gene) = part.split('*').next() {
            if !gene.is_empty() {
                genes.insert(gene.to_string());
            }
        }
    }
    genes
}

/// Drop a trailing `D` on the gene (e.g. `IGHV1-69D` → `IGHV1-69`).
pub fn normalize_gene_name(gene: &str) -> String {
    gene.strip_suffix('D').unwrap_or(gene).to_string()
}

/// `(normalized_gene, _S…)` ignoring the allele number.
pub fn extract_gene_s_signature(allele_name: &str) -> Option<(String, String)> {
    let allele_name = normalize_allele(allele_name);
    let suffix = allele_name
        .rfind("_S")
        .map(|i| &allele_name[i..])
        .filter(|s| !s.contains('/'))?;
    let gene = allele_name.split('*').next()?;
    Some((normalize_gene_name(gene), suffix.to_string()))
}

/// Expand `A*01/B*02_S1234` into `{A*01_S1234, B*02_S1234}`.
pub fn expand_allele_combinations(allele_name: &str) -> BTreeSet<String> {
    let allele_name = normalize_allele(allele_name);
    let mut expanded = BTreeSet::new();
    if !allele_name.contains('/') {
        expanded.insert(allele_name);
        return expanded;
    }
    if !allele_name.contains("_S") {
        expanded.insert(allele_name);
        return expanded;
    }
    let Some(suffix_start) = allele_name.rfind("_S") else {
        expanded.insert(allele_name);
        return expanded;
    };
    let suffix = &allele_name[suffix_start..];
    if suffix.contains('/') {
        expanded.insert(allele_name);
        return expanded;
    }
    let parts: Vec<&str> = allele_name.split('/').collect();
    for (i, part) in parts.iter().enumerate() {
        let part = part.trim();
        if i + 1 == parts.len() {
            expanded.insert(part.to_string());
        } else {
            let base = part.split("_S").next().unwrap_or(part);
            if !base.is_empty() {
                expanded.insert(format!("{base}{suffix}"));
            }
        }
    }
    expanded
}

/// KI `_S` (not `_F`) alleles of this gene whose cohort contains `KI`.
pub fn load_ki_target_alleles(alleles_tsv: &Path, gene: GeneKind) -> Result<HashSet<String>> {
    let frame = read_table(
        alleles_tsv,
        &ReadOptions::tsv()
            .text_column("db_name")
            .text_column("cohort"),
    )?;
    let db_names = column_as_strings(&frame, "db_name")?;
    let cohorts = column_as_strings(&frame, "cohort")?;
    let prefix = gene.prefix();
    let mut out = HashSet::new();
    for (db_name, cohort) in db_names.iter().zip(cohorts.iter()) {
        let (Some(db_name), Some(cohort)) = (db_name.as_deref(), cohort.as_deref()) else {
            continue;
        };
        if db_name.starts_with(prefix)
            && db_name.contains("_S")
            && !db_name.contains("_F")
            && cohort.contains("KI")
        {
            out.insert(normalize_allele(db_name));
        }
    }
    Ok(out)
}

/// KI targets present in both `alleles.tsv` and the target FASTA.
pub fn load_target_alleles(
    fasta_path: &Path,
    alleles_tsv: &Path,
    gene: GeneKind,
) -> Result<(BTreeSet<String>, BTreeMap<String, BTreeSet<String>>)> {
    let ki = load_ki_target_alleles(alleles_tsv, gene)?;
    let mut target_alleles = BTreeSet::new();
    let mut target_genes: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for record in read_fasta(fasta_path)? {
        let allele_name = normalize_allele(record.identifier());
        if !ki.contains(&allele_name) {
            continue;
        }
        let genes = extract_genes(&allele_name);
        if genes.is_empty() {
            continue;
        }
        target_alleles.insert(allele_name.clone());
        for g in genes {
            target_genes.entry(g).or_default().insert(allele_name.clone());
        }
    }
    Ok((target_alleles, target_genes))
}

/// Compare a BLAST result table against KI targets.
pub fn analyze_blast_results(
    blast_path: &Path,
    target_fasta: &Path,
    alleles_tsv: &Path,
    gene: GeneKind,
) -> Result<RecoveryReport> {
    let (target_alleles, target_genes) = load_target_alleles(target_fasta, alleles_tsv, gene)?;

    let mut target_signatures: HashMap<(String, String), HashSet<String>> = HashMap::new();
    let mut target_suffixes: HashMap<String, HashSet<String>> = HashMap::new();
    for allele in &target_alleles {
        if let Some((g, suffix)) = extract_gene_s_signature(allele) {
            target_signatures
                .entry((g, suffix.clone()))
                .or_default()
                .insert(allele.clone());
            target_suffixes
                .entry(suffix)
                .or_default()
                .insert(allele.clone());
        }
    }

    let frame = read_table(
        blast_path,
        &ReadOptions {
            delimiter: None,
            text_columns: vec!["allele_name".into()],
            ..ReadOptions::default()
        },
    )?;
    if !frame.has_column("allele_name") {
        return Err(Error::invalid(format!(
            "{}: column 'allele_name' not found",
            blast_path.display()
        )));
    }
    let allele_names = column_as_strings(&frame, "allele_name")?;

    let mut found_targets = BTreeSet::new();
    let mut all_s_alleles = BTreeSet::new();

    for cell in &allele_names {
        let Some(raw) = cell.as_deref() else {
            continue;
        };
        let allele = normalize_allele(raw);
        for expanded in expand_allele_combinations(&allele) {
            if expanded.contains("_S") {
                all_s_alleles.insert(expanded.clone());
            }
            let Some((g, suffix)) = extract_gene_s_signature(&expanded) else {
                continue;
            };
            if let Some(matched) = target_signatures.get(&(g, suffix.clone())) {
                found_targets.extend(matched.iter().cloned());
            } else if let Some(matched) = target_suffixes.get(&suffix) {
                found_targets.extend(matched.iter().cloned());
            }
        }
    }

    let false_positives: BTreeSet<String> = all_s_alleles
        .iter()
        .filter(|a| {
            let Some((g, suffix)) = extract_gene_s_signature(a) else {
                return true;
            };
            !target_signatures.contains_key(&(g, suffix.clone()))
                && !target_suffixes.contains_key(&suffix)
        })
        .cloned()
        .collect();

    let mut per_gene = BTreeMap::new();
    for (gene_name, alleles) in &target_genes {
        let found_alleles: BTreeSet<String> = alleles
            .iter()
            .filter(|a| found_targets.contains(*a))
            .cloned()
            .collect();
        let missing_alleles: BTreeSet<String> = alleles
            .iter()
            .filter(|a| !found_targets.contains(*a))
            .cloned()
            .collect();
        per_gene.insert(
            gene_name.clone(),
            GeneStats {
                total_target: alleles.len(),
                found: found_alleles.len(),
                found_alleles,
                missing_alleles,
            },
        );
    }

    let missing_targets: BTreeSet<String> = target_alleles
        .difference(&found_targets)
        .cloned()
        .collect();

    Ok(RecoveryReport {
        gene,
        target_alleles,
        found_targets,
        missing_targets,
        all_s_alleles,
        false_positives,
        per_gene,
    })
}

/// Human-readable report matching the Julia checker's layout.
pub fn format_report(report: &RecoveryReport) -> String {
    use std::fmt::Write;
    let mut out = String::new();
    let _ = writeln!(out, "{}", "=".repeat(80));
    let _ = writeln!(out, "Target Allele Recovery Report (per gene)");
    let _ = writeln!(out, "{}", "=".repeat(80));
    let _ = writeln!(
        out,
        "{:<20} {:>12} {:>15} {:>10}",
        "Gene", "Total Target", "Found in BLAST", "Recovery %"
    );
    let _ = writeln!(out, "{}", "-".repeat(80));

    let mut total_target = 0usize;
    let mut total_found = 0usize;
    for (gene, stats) in &report.per_gene {
        let recovery = if stats.total_target == 0 {
            0.0
        } else {
            100.0 * stats.found as f64 / stats.total_target as f64
        };
        total_target += stats.total_target;
        total_found += stats.found;
        let _ = writeln!(
            out,
            "{:<20} {:>12} {:>15} {:>9.1}%",
            gene, stats.total_target, stats.found, recovery
        );
    }
    let overall = if total_target == 0 {
        0.0
    } else {
        100.0 * total_found as f64 / total_target as f64
    };
    let _ = writeln!(out, "{}", "-".repeat(80));
    let _ = writeln!(
        out,
        "{:<20} {:>12} {:>15} {:>9.1}%",
        "TOTAL", total_target, total_found, overall
    );
    let _ = writeln!(out, "{}", "=".repeat(80));

    let fp_pct = report.false_positive_rate();
    let _ = writeln!(out);
    let _ = writeln!(out, "{}", "=".repeat(80));
    let _ = writeln!(out, "FALSE POSITIVE ANALYSIS");
    let _ = writeln!(out, "{}", "=".repeat(80));
    let _ = writeln!(
        out,
        "Total _S alleles found in BLAST:     {}",
        report.all_s_alleles.len()
    );
    let _ = writeln!(
        out,
        "Target _S alleles found:              {}",
        report.found_targets.len()
    );
    let _ = writeln!(
        out,
        "Additional _S alleles (non-target):  {}",
        report.false_positives.len()
    );
    let _ = writeln!(out, "False positive rate:                 {fp_pct:.1}%");
    if fp_pct > 50.0 {
        let _ = writeln!(out);
        let _ = writeln!(
            out,
            "⚠ WARNING: High false positive rate ({fp_pct:.1}%)"
        );
        let _ = writeln!(
            out,
            "  This suggests BLAST parameters may be too relaxed, allowing many non-target alleles."
        );
    } else if fp_pct > 20.0 {
        let _ = writeln!(out);
        let _ = writeln!(
            out,
            "⚠ NOTE: Moderate false positive rate ({fp_pct:.1}%)"
        );
        let _ = writeln!(
            out,
            "  Some non-target alleles are being detected. This may be acceptable for maximum recovery."
        );
    } else {
        let _ = writeln!(out);
        let _ = writeln!(out, "✓ Low false positive rate - parameters appear well-tuned.");
    }
    let _ = writeln!(out, "{}", "=".repeat(80));

    let _ = writeln!(out);
    let _ = writeln!(out, "{}", "=".repeat(80));
    let _ = writeln!(out, "Missing Alleles (genes with < 100% recovery)");
    let _ = writeln!(out, "{}", "=".repeat(80));
    let mut has_missing = false;
    for (gene, stats) in &report.per_gene {
        if stats.found < stats.total_target {
            has_missing = true;
            let _ = writeln!(
                out,
                "\n{gene} ({}/{} found, {} missing):",
                stats.found,
                stats.total_target,
                stats.missing_alleles.len()
            );
            for allele in &stats.missing_alleles {
                let _ = writeln!(out, "  - {allele}");
            }
        }
    }
    if !has_missing {
        let _ = writeln!(out);
        let _ = writeln!(out, "All target alleles were found in BLAST results!");
    }
    let _ = writeln!(out, "{}", "=".repeat(80));

    if !report.false_positives.is_empty() {
        let _ = writeln!(out);
        let _ = writeln!(out, "{}", "=".repeat(80));
        let _ = writeln!(
            out,
            "Extra _S alleles (not KI targets), {} total:",
            report.false_positives.len()
        );
        let _ = writeln!(out, "{}", "=".repeat(80));
        for allele in &report.false_positives {
            let _ = writeln!(out, "  - {allele}");
        }
        let _ = writeln!(out, "{}", "=".repeat(80));
    }
    out
}

fn column_as_strings(frame: &Frame, name: &str) -> Result<Vec<Option<String>>> {
    let column = frame.column(name)?;
    if let Some(strs) = column.as_strs() {
        return Ok(strs.to_vec());
    }
    Ok(column
        .iter()
        .map(|v| match v {
            Value::Missing => None,
            other => Some(other.julia_string()),
        })
        .collect())
}

/// Split a combined FASTA into origin (no `_S`) and target (has `_S`).
pub fn split_fasta_by_s_suffix(
    input: &Path,
    origin_out: &Path,
    target_out: &Path,
) -> Result<(usize, usize)> {
    let records = read_fasta(input)?;
    let mut origin: Vec<(&str, &str)> = Vec::new();
    let mut target: Vec<(&str, &str)> = Vec::new();
    for record in &records {
        if record.identifier().contains("_S") {
            target.push((&record.description, &record.sequence));
        } else {
            origin.push((&record.description, &record.sequence));
        }
    }
    let origin_n = origin.len();
    let target_n = target.len();
    write_fasta(origin_out, origin, Some(crate::io::fasta::DEFAULT_LINE_WIDTH))?;
    write_fasta(target_out, target, Some(crate::io::fasta::DEFAULT_LINE_WIDTH))?;
    Ok((origin_n, target_n))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalize_strips_novel() {
        assert_eq!(
            normalize_allele("IGHV3-23*04_S7240 Novel"),
            "IGHV3-23*04_S7240"
        );
    }

    #[test]
    fn expand_slash_alleles() {
        let expanded =
            expand_allele_combinations("IGHV3-30*02/IGHV3-30-5*02_S9143");
        assert!(expanded.contains("IGHV3-30*02_S9143"));
        assert!(expanded.contains("IGHV3-30-5*02_S9143"));
    }

    #[test]
    fn signature_drops_d_and_allele_number() {
        assert_eq!(
            extract_gene_s_signature("IGHV1-69D*01_S9407"),
            Some(("IGHV1-69".into(), "_S9407".into()))
        );
        assert_eq!(
            extract_gene_s_signature("IGHV1-46*03_S9407"),
            Some(("IGHV1-46".into(), "_S9407".into()))
        );
    }

    #[test]
    fn extract_genes_handles_slash() {
        let genes = extract_genes("IGHD4-17*01/IGHD4-4*01_S0251");
        assert!(genes.contains("IGHD4-17"));
        assert!(genes.contains("IGHD4-4"));
    }
}
