//! KI target-allele recovery self-check.
//!
//! Ports the matching logic of
//! `ImmuneDiscover-IGH-2026-Check/check_target_alleles.jl`: load KI `_S`
//! targets from `alleles.tsv` ∩ target FASTA, match them against
//! `discover blast` `allele_name` values (flexible gene / `_S` / slash
//! expansion), and report recovery plus extra `_S` calls (false
//! positives).

pub mod allele_recovery;
pub mod evidence_report;

pub use allele_recovery::{
    analyze_blast_results, format_report, load_ki_target_alleles, load_target_alleles,
    split_fasta_by_s_suffix, GeneKind, RecoveryReport,
};
