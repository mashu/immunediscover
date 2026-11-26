//! KI target-allele recovery self-check.
//!
//! Runs `discover blast` on the KI demultiplex tables (or reuses existing
//! BLAST outputs) and reports whether every KI `_S` target allele is
//! recovered, plus how many extra `_S` calls appear (false positives).
//!
//! Matching logic mirrors `ImmuneDiscover-IGH-2026-Check/check_target_alleles.jl`.
//!
//! ```text
//! cargo run --release --bin ki_check -- \
//!   --check-root /path/to/ImmuneDiscover-IGH-2026-Check
//!
//! # Analyze existing BLAST tables only:
//! cargo run --release --bin ki_check -- \
//!   --check-root /path/to/ImmuneDiscover-IGH-2026-Check --skip-discover
//! ```

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;

use immunediscover::check::allele_recovery::{
    analyze_blast_results, format_report, split_fasta_by_s_suffix, GeneKind, RecoveryReport,
};
use immunediscover::check::evidence_report::{classify, format_evidence_report, Thresholds};
use immunediscover::cli::presets::{resolve, Preset, PresetValues, Resolved};
use immunediscover::commands::discover::blast::flank_split::unfiltered_sidecar;
use immunediscover::commands::discover::blast::{self, BlastParams};
use immunediscover::error::{Error, Result};
use immunediscover::evidence;

#[derive(Debug, Parser)]
#[command(
    name = "ki-check",
    about = "Self-check: recover KI _S target alleles with discover blast and report misses / false positives"
)]
struct Args {
    /// Root of ImmuneDiscover-IGH-2026-Check (contains data/KI and data/db)
    #[arg(long = "check-root")]
    check_root: PathBuf,

    /// Directory for BLAST outputs and analysis reports
    #[arg(long = "out-dir")]
    out_dir: Option<PathBuf>,

    /// Skip discover blast; only analyze existing KI-blast-{V,D,J}.tsv.gz
    #[arg(long = "skip-discover")]
    skip_discover: bool,

    /// Overwrite BLAST cache / outputs when running discovery
    #[arg(short = 'o', long = "overwrite")]
    overwrite: bool,

    /// Gene segments to check (comma-separated: V,D,J)
    #[arg(long = "genes", default_value = "V,D,J")]
    genes: String,

    /// Fail if any target allele is missing (default: true)
    #[arg(long = "allow-misses", default_value_t = false)]
    allow_misses: bool,

    /// Fail when false-positive rate (%) exceeds this (omit to never fail on FP)
    #[arg(long = "max-fp-rate")]
    max_fp_rate: Option<f64>,

    /// Work directory for discover blast intermediates
    #[arg(long = "work-dir", default_value = "blast")]
    work_dir: String,

    /// Also evaluate the evidence statistics (dosage, artifact expectation,
    /// provisional score) against the KI targets. Needs the unfiltered
    /// cluster table, which discovery writes when this flag is set; with
    /// --skip-discover it must already exist. The report is aggregate only.
    #[arg(long = "evidence")]
    evidence: bool,

    /// Top unexplained candidates listed in the evidence report (gene and
    /// nearest allele names only; 0 leaves the list out)
    #[arg(long = "evidence-list", default_value_t = 20)]
    evidence_list: usize,
}

struct GenePaths {
    kind: GeneKind,
    reads: PathBuf,
    origin_fasta: PathBuf,
    target_fasta: PathBuf,
    blast_out: PathBuf,
}

fn main() -> ExitCode {
    match run() {
        Ok(ok) => {
            if ok {
                ExitCode::SUCCESS
            } else {
                ExitCode::FAILURE
            }
        }
        Err(error) => {
            eprintln!("ERROR: {error}");
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<bool> {
    let args = Args::parse();
    let root = args.check_root.canonicalize().map_err(|e| {
        Error::io(
            &args.check_root,
            std::io::Error::new(e.kind(), format!("cannot resolve --check-root: {e}")),
        )
    })?;
    let out_dir = args
        .out_dir
        .clone()
        .unwrap_or_else(|| root.join("analysis_results_rs"));
    std::fs::create_dir_all(&out_dir).map_err(|e| Error::io(&out_dir, e))?;

    let db_dir = root.join("data/db");
    let alleles_tsv = db_dir.join("alleles.tsv");
    ensure_file(&alleles_tsv)?;
    ensure_origin_target(&db_dir)?;

    let genes = parse_genes(&args.genes)?;
    let mut reports: Vec<RecoveryReport> = Vec::new();

    for kind in genes {
        let paths = gene_paths(&root, &out_dir, kind)?;
        if !args.skip_discover {
            println!(
                "==> discover blast -g {} → {}",
                kind.as_str(),
                paths.blast_out.display()
            );
            run_discover(&paths, &args)?;
        } else {
            ensure_file(&paths.blast_out)?;
            println!("==> analyzing existing {}", paths.blast_out.display());
        }

        let report =
            analyze_blast_results(&paths.blast_out, &paths.target_fasta, &alleles_tsv, kind)?;
        let text = format_report(&report);
        print!("{text}");
        let analysis_path = out_dir.join(format!("analysis-{}.txt", kind.as_str()));
        std::fs::write(&analysis_path, &text).map_err(|e| Error::io(&analysis_path, e))?;
        reports.push(report);

        if args.evidence {
            let text = evidence_report(&paths, &args)?;
            print!("{text}");
            let path = out_dir.join(format!("evidence-{}.txt", kind.as_str()));
            std::fs::write(&path, &text).map_err(|e| Error::io(&path, e))?;
            println!("Evidence report saved to: {}", path.display());
        }
    }

    let summary = format_overall_summary(&reports);
    print!("{summary}");
    let summary_path = out_dir.join("final_report.txt");
    std::fs::write(&summary_path, &summary).map_err(|e| Error::io(&summary_path, e))?;
    println!("Full report saved to: {}", summary_path.display());

    let mut ok = true;
    if !args.allow_misses {
        for report in &reports {
            if !report.all_recovered() {
                eprintln!(
                    "FAIL: {} missed {} / {} unique KI targets",
                    report.gene.as_str(),
                    report.missing_targets.len(),
                    report.unique_target_count()
                );
                ok = false;
            }
        }
    }
    if let Some(max_fp) = args.max_fp_rate {
        for report in &reports {
            let fp = report.false_positive_rate();
            if fp > max_fp {
                eprintln!(
                    "FAIL: {} false-positive rate {fp:.1}% exceeds --max-fp-rate {max_fp}",
                    report.gene.as_str()
                );
                ok = false;
            }
        }
    }
    if ok {
        println!("PASS: all checked KI targets recovered within FP limits.");
    }
    Ok(ok)
}

/// Evaluate the evidence statistics on the unfiltered clusters of a gene
/// segment against its database and target FASTAs.
fn evidence_report(paths: &GenePaths, args: &Args) -> Result<String> {
    let clusters = unfiltered_sidecar(&paths.blast_out);
    if !clusters.is_file() {
        return Err(Error::invalid(format!(
            "missing {}: run without --skip-discover, or run discover blast with -v",
            clusters.display()
        )));
    }
    let (mut observations, n_donors) = evidence::read_clusters(&clusters)?;
    let models = evidence::assess(&mut observations, n_donors);
    let classes = classify(&observations, &paths.origin_fasta, &paths.target_fasta)?;
    let resolved = ki_resolved(paths.kind);
    let thresholds = Thresholds {
        min_count: resolved.min_full_count.max(0) as u64,
        min_ratio: resolved.min_full_ratio,
    };
    Ok(format!(
        "\n======== Evidence, {} ========\n{}",
        paths.kind.as_str(),
        format_evidence_report(
            &observations,
            &classes,
            n_donors,
            &models,
            thresholds,
            args.evidence_list
        )
    ))
}

fn parse_genes(s: &str) -> Result<Vec<GeneKind>> {
    let mut out = Vec::new();
    for part in s.split(',') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        out.push(GeneKind::parse(part)?);
    }
    if out.is_empty() {
        return Err(Error::invalid("--genes is empty"));
    }
    Ok(out)
}

fn gene_paths(root: &Path, out_dir: &Path, kind: GeneKind) -> Result<GenePaths> {
    let (reads_name, fasta_name) = match kind {
        GeneKind::V => (
            "KI_IMD_IGHV_old38-42.demultiplexed.combined.tsv.gz",
            "KI+1KGP-IGHV-SHORT.fasta",
        ),
        GeneKind::D => (
            "KI_IMD_IGHDJ_old38-42.demultiplexed.combined.tsv.gz",
            "KI+1KGP-IGHD-SHORT.fasta",
        ),
        GeneKind::J => (
            "KI_IMD_IGHDJ_old38-42.demultiplexed.combined.tsv.gz",
            "KI+1KGP-IGHJ-SHORT.fasta",
        ),
    };
    let reads = root.join("data/KI").join(reads_name);
    let origin_fasta = root.join("data/db/origin").join(fasta_name);
    let target_fasta = root.join("data/db/target").join(fasta_name);
    ensure_file(&reads)?;
    ensure_file(&origin_fasta)?;
    ensure_file(&target_fasta)?;
    Ok(GenePaths {
        kind,
        reads,
        origin_fasta,
        target_fasta,
        blast_out: out_dir.join(format!("KI-blast-{}.tsv.gz", kind.as_str())),
    })
}

fn ensure_file(path: &Path) -> Result<()> {
    if path.is_file() {
        Ok(())
    } else {
        Err(Error::invalid(format!("missing file: {}", path.display())))
    }
}

fn ensure_origin_target(db_dir: &Path) -> Result<()> {
    let origin = db_dir.join("origin");
    let target = db_dir.join("target");
    let need_split = !origin.is_dir()
        || !target.is_dir()
        || std::fs::read_dir(&origin)
            .map(|mut d| d.next().is_none())
            .unwrap_or(true)
        || std::fs::read_dir(&target)
            .map(|mut d| d.next().is_none())
            .unwrap_or(true);
    if !need_split {
        // Still ensure the three SHORT FASTAs exist; if not, rebuild from combined.
        for name in [
            "KI+1KGP-IGHV-SHORT.fasta",
            "KI+1KGP-IGHD-SHORT.fasta",
            "KI+1KGP-IGHJ-SHORT.fasta",
        ] {
            if !origin.join(name).is_file() || !target.join(name).is_file() {
                split_one(db_dir, name)?;
            }
        }
        return Ok(());
    }
    std::fs::create_dir_all(&origin).map_err(|e| Error::io(&origin, e))?;
    std::fs::create_dir_all(&target).map_err(|e| Error::io(&target, e))?;
    for name in [
        "KI+1KGP-IGHV-SHORT.fasta",
        "KI+1KGP-IGHD-SHORT.fasta",
        "KI+1KGP-IGHJ-SHORT.fasta",
    ] {
        split_one(db_dir, name)?;
    }
    Ok(())
}

fn split_one(db_dir: &Path, name: &str) -> Result<()> {
    let input = db_dir.join(name);
    ensure_file(&input)?;
    let (o, t) = split_fasta_by_s_suffix(
        &input,
        &db_dir.join("origin").join(name),
        &db_dir.join("target").join(name),
    )?;
    println!("  split {name}: origin={o}, target={t}");
    Ok(())
}

/// Parameters tuned for full KI target recovery under Rust affixes.
///
/// Julia used `-c 18` for V; Rust's 20 nt affix consensus splits some
/// alleles across extension variants (e.g. IGHV1-46*01_S9407 best
/// `full_count` 15), so V uses 15 while keeping the `full_count` filter.
fn ki_resolved(kind: GeneKind) -> Resolved {
    match kind {
        GeneKind::V => resolve(
            Some(Preset::V),
            PresetValues::<String> {
                min_full_ratio: Some(0.083102),
                min_full_count: Some(15),
                length: Some(299),
                min_core_coverage: Some(0.98),
                max_distance: Some(10),
                min_quality: Some(0.5),
                ..PresetValues::default()
            },
        ),
        GeneKind::D | GeneKind::J => {
            let preset = match kind {
                GeneKind::D => Preset::D,
                GeneKind::J => Preset::J,
                GeneKind::V => unreachable!(),
            };
            resolve(
                Some(preset),
                PresetValues::<String> {
                    forward: Some(20),
                    reverse: Some(20),
                    min_full_ratio: Some(0.01),
                    min_full_count: Some(10),
                    ..PresetValues::default()
                },
            )
        }
    }
}

fn run_discover(paths: &GenePaths, args: &Args) -> Result<()> {
    let resolved = ki_resolved(paths.kind);
    let params = BlastParams {
        input: paths.reads.clone(),
        fasta: paths.origin_fasta.clone(),
        output: paths.blast_out.clone(),
        pseudo: None,
        min_full_count: resolved.min_full_count,
        min_full_ratio: resolved.min_full_ratio,
        subject_coverage: resolved.subject_coverage,
        max_distance: resolved.max_distance,
        min_length: resolved.length,
        edge: resolved.edge,
        blast_args: resolved.blast_args,
        overwrite: args.overwrite,
        forward_extension: resolved.forward,
        reverse_extension: resolved.reverse,
        extension_max_distance: 0,
        gene_context: true,
        min_quality: resolved.min_quality,
        name_substrings_as_known: true,
        keep_failed: false,
        min_core_coverage: resolved.min_core_coverage,
        top: 1,
        // The evidence report reads the unfiltered clusters, written with -v.
        verbose: args.evidence,
        work_dir: args.work_dir.clone(),
    };
    println!("  params: {params:?}");
    blast::run(&params)
}

fn format_overall_summary(reports: &[RecoveryReport]) -> String {
    use std::fmt::Write;
    let mut out = String::new();
    let _ = writeln!(out, "{}", "=".repeat(80));
    let _ = writeln!(
        out,
        "TARGET ALLELE RECOVERY ANALYSIS - FINAL SUMMARY REPORT"
    );
    let _ = writeln!(out, "{}", "=".repeat(80));
    let _ = writeln!(out, "Generated: {}", chrono::Local::now());
    let _ = writeln!(out);
    let _ = writeln!(
        out,
        "Gene Type    Total Target    Found in BLAST    Recovery %    Unique targets    Unique found    Extra _S    FP %"
    );
    let _ = writeln!(out, "{}", "-".repeat(100));

    let mut sum_total = 0usize;
    let mut sum_found = 0usize;
    let mut sum_unique_t = 0usize;
    let mut sum_unique_f = 0usize;
    let mut sum_all_s = 0usize;
    let mut sum_fp = 0usize;

    for report in reports {
        let total = report.total_target();
        let found = report.total_found();
        let pct = report.recovery_pct();
        let fp = report.false_positive_rate();
        sum_total += total;
        sum_found += found;
        sum_unique_t += report.unique_target_count();
        sum_unique_f += report.unique_found_count();
        sum_all_s += report.all_s_alleles.len();
        sum_fp += report.false_positives.len();
        let _ = writeln!(
            out,
            "{:<8} {:>14} {:>16} {:>11.1}% {:>16} {:>14} {:>10} {:>6.1}%",
            format!("{} Genes", report.gene.as_str()),
            total,
            found,
            pct,
            report.unique_target_count(),
            report.unique_found_count(),
            report.false_positives.len(),
            fp
        );
    }
    let overall = if sum_total == 0 {
        0.0
    } else {
        100.0 * sum_found as f64 / sum_total as f64
    };
    let fp_all = if sum_all_s == 0 {
        0.0
    } else {
        100.0 * sum_fp as f64 / sum_all_s as f64
    };
    let _ = writeln!(out, "{}", "-".repeat(100));
    let _ = writeln!(
        out,
        "{:<8} {:>14} {:>16} {:>11.1}% {:>16} {:>14} {:>10} {:>6.1}%",
        "TOTAL", sum_total, sum_found, overall, sum_unique_t, sum_unique_f, sum_fp, fp_all
    );
    let _ = writeln!(out, "{}", "=".repeat(80));
    out
}
