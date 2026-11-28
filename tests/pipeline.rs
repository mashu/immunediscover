//! End-to-end runs of the command-line program on a small synthetic plate.
//!
//! Nothing here comes from real samples: genes, alleles, donors and reads
//! are generated from a fixed seed. The plate is small enough for the whole
//! file to run in a few seconds.
//!
//! `discover blast` with BLAST needs BLAST+. Without `blastn` on the `PATH`
//! those parts are skipped, unless `IMMUNEDISCOVER_REQUIRE_BLAST` is set (as
//! in CI), in which case a missing `blastn` fails the test. The built-in
//! aligner is tested everywhere.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use immunediscover::frame::Frame;
use immunediscover::io::tsv::read_tsv;

const DONORS: usize = 12;
const GENES: usize = 5;
const READS_PER_COPY: usize = 30;
/// The donors carrying the allele that is missing from the database.
const NOVEL_CARRIERS: [usize; 3] = [0, 3, 6];
const NOVEL_GENE: usize = 1;
const HEPTAMER: &str = "CACAGTG";

/// Small deterministic generator (xorshift64*), so that the plate is the
/// same on every machine without depending on a random number crate.
struct Generator(u64);

impl Generator {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }

    fn sequence(&mut self, length: usize) -> String {
        (0..length)
            .map(|_| b"ACGT"[self.below(4)] as char)
            .collect()
    }
}

fn substitute(sequence: &str, position: usize) -> String {
    let mut bytes = sequence.as_bytes().to_vec();
    bytes[position] = match bytes[position] {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    };
    String::from_utf8(bytes).expect("nucleotides")
}

struct Gene {
    name: String,
    /// Bases before the gene in every read of it.
    upstream: String,
    /// Bases after the heptamer.
    downstream: String,
    alleles: [String; 2],
}

/// The synthetic plate and what is true about it.
struct Plate {
    genes: Vec<Gene>,
    novel: String,
    /// `(forward index, reverse index, case)` per well.
    indices: Vec<(String, String, String)>,
    /// `[donor][gene]`: the sequences on the two chromosomes.
    genotypes: Vec<Vec<[String; 2]>>,
    reads: Vec<String>,
}

impl Plate {
    fn generate() -> Plate {
        let mut rng = Generator(0x9E37_79B9_7F4A_7C15);
        let genes: Vec<Gene> = (0..GENES)
            .map(|g| {
                let first = rng.sequence(300);
                let second = [60, 150, 240]
                    .iter()
                    .fold(first.clone(), |s, p| substitute(&s, *p));
                Gene {
                    name: format!("IGHV{}-1", g + 1),
                    upstream: rng.sequence(30),
                    downstream: rng.sequence(40),
                    alleles: [first, second],
                }
            })
            .collect();
        // One difference from the first allele of one gene.
        let novel = substitute(&genes[NOVEL_GENE].alleles[0], 100);
        let indices: Vec<(String, String, String)> = (0..DONORS)
            .map(|d| (rng.sequence(12), rng.sequence(12), format!("D{:02}", d + 1)))
            .collect();
        let genotypes: Vec<Vec<[String; 2]>> = (0..DONORS)
            .map(|d| {
                (0..GENES)
                    .map(|g| {
                        let [first, second] = genes[g].alleles.clone();
                        if g == NOVEL_GENE && NOVEL_CARRIERS.contains(&d) {
                            [first, novel.clone()]
                        } else {
                            match (d + g) % 3 {
                                0 => [first, second],
                                1 => [first.clone(), first],
                                _ => [second.clone(), second],
                            }
                        }
                    })
                    .collect()
            })
            .collect();
        let mut reads = Vec::new();
        for (d, (forward, reverse, _)) in indices.iter().enumerate() {
            for (g, gene) in genes.iter().enumerate() {
                let read = |allele: &str| {
                    format!(
                        "{forward}{}{allele}{HEPTAMER}{}{reverse}",
                        gene.upstream, gene.downstream
                    )
                };
                for copy in &genotypes[d][g] {
                    reads.extend((0..READS_PER_COPY).map(|_| read(copy)));
                    // Two reads with a sequencing error each.
                    for _ in 0..2 {
                        reads.push(read(&substitute(copy, 5 + rng.below(290))));
                    }
                }
            }
        }
        Plate {
            genes,
            novel,
            indices,
            genotypes,
            reads,
        }
    }

    fn write(&self, dir: &Path) {
        let mut fastq = String::new();
        for (i, read) in self.reads.iter().enumerate() {
            let _ = writeln!(fastq, "@read{i}\n{read}\n+\n{}", "I".repeat(read.len()));
        }
        std::fs::write(dir.join("plate.fastq"), fastq).unwrap();
        let mut sheet = String::from("forward_index\treverse_index\tcase\n");
        for (forward, reverse, case) in &self.indices {
            let _ = writeln!(sheet, "{forward}\t{reverse}\t{case}");
        }
        std::fs::write(dir.join("indices.tsv"), sheet).unwrap();
        let mut database = String::new();
        for gene in &self.genes {
            for (i, allele) in gene.alleles.iter().enumerate() {
                let _ = writeln!(database, ">{}*0{}\n{allele}", gene.name, i + 1);
            }
        }
        std::fs::write(dir.join("db.fasta"), database).unwrap();
    }

    fn allele_name(&self, gene: usize, sequence: &str) -> Option<String> {
        let allele = self.genes[gene]
            .alleles
            .iter()
            .position(|a| a == sequence)?;
        Some(format!("{}*0{}", self.genes[gene].name, allele + 1))
    }
}

fn run(dir: &Path, arguments: &[&str]) -> Output {
    let output = Command::new(env!("CARGO_BIN_EXE_immunediscover"))
        .args(arguments)
        .current_dir(dir)
        .output()
        .expect("the program starts");
    assert!(
        output.status.success(),
        "immunediscover {} failed:\n{}{}",
        arguments.join(" "),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn printed(output: &Output) -> String {
    format!(
        "{}{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    )
}

fn column(frame: &Frame, name: &str) -> Vec<String> {
    frame
        .column(name)
        .unwrap_or_else(|_| panic!("column {name} is missing"))
        .iter()
        .map(|v| {
            if v.is_missing() {
                String::new()
            } else {
                v.julia_string()
            }
        })
        .collect()
}

fn table(path: PathBuf) -> Frame {
    read_tsv(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// The plate, demultiplexed in a fresh directory.
fn demultiplexed() -> (tempfile::TempDir, Plate) {
    let dir = tempfile::tempdir().unwrap();
    let plate = Plate::generate();
    plate.write(dir.path());
    run(
        dir.path(),
        &[
            "preprocess",
            "demultiplex",
            "plate.fastq",
            "indices.tsv",
            "demux.tsv",
        ],
    );
    (dir, plate)
}

fn blastn_available() -> bool {
    let found = Command::new("blastn").arg("-version").output().is_ok();
    if !found && std::env::var_os("IMMUNEDISCOVER_REQUIRE_BLAST").is_some() {
        panic!("blastn is required for this run (IMMUNEDISCOVER_REQUIRE_BLAST) but was not found");
    }
    found
}

#[test]
fn every_read_goes_to_its_donor() {
    let (dir, plate) = demultiplexed();
    let demux = table(dir.path().join("demux.tsv.gz"));
    assert_eq!(demux.nrows(), plate.reads.len());
    let mut per_case: BTreeMap<String, usize> = BTreeMap::new();
    for case in column(&demux, "case") {
        *per_case.entry(case).or_default() += 1;
    }
    assert_eq!(per_case.len(), DONORS);
    let per_donor = plate.reads.len() / DONORS;
    assert!(per_case.values().all(|n| *n == per_donor), "{per_case:?}");
    // Wells are numbered by the rows of the index sheet.
    let wells: BTreeSet<String> = column(&demux, "well").into_iter().collect();
    assert_eq!(wells.len(), DONORS);
}

#[test]
fn known_alleles_are_found_with_their_reads_and_genotypes_follow() {
    let (dir, plate) = demultiplexed();
    run(
        dir.path(),
        &[
            "search",
            "exact",
            "demux.tsv.gz",
            "db.fasta",
            "exact.tsv",
            "-g",
            "V",
            "-n",
        ],
    );
    let exact = table(dir.path().join("exact.tsv.gz"));
    let found: BTreeMap<(String, String), i64> = column(&exact, "case")
        .into_iter()
        .zip(column(&exact, "db_name"))
        .zip(column(&exact, "count"))
        .map(|((case, allele), count)| ((case, allele), count.parse().unwrap()))
        .collect();
    // What the genotypes say: reads per database allele and donor. Reads
    // with an error do not contain the allele and are not counted.
    let mut expected: BTreeMap<(String, String), i64> = BTreeMap::new();
    for (d, donor) in plate.genotypes.iter().enumerate() {
        for (g, copies) in donor.iter().enumerate() {
            for copy in copies {
                if let Some(name) = plate.allele_name(g, copy) {
                    *expected
                        .entry((plate.indices[d].2.clone(), name))
                        .or_default() += READS_PER_COPY as i64;
                }
            }
        }
    }
    assert_eq!(found, expected);
    let heptamers: BTreeSet<String> = column(&exact, "heptamer").into_iter().collect();
    assert_eq!(heptamers, BTreeSet::from([HEPTAMER.to_string()]));

    run(
        dir.path(),
        &["analyze", "haplotype", "exact.tsv.gz", "genotypes.tsv"],
    );
    let path = ["genotypes.tsv", "genotypes.tsv.gz"]
        .iter()
        .map(|name| dir.path().join(name))
        .find(|path| path.is_file())
        .expect("the genotype table was written");
    let genotypes = table(path);
    let called: BTreeMap<(String, String), String> = column(&genotypes, "case")
        .into_iter()
        .zip(column(&genotypes, "gene"))
        .zip(column(&genotypes, "genotype"))
        .collect();
    for (d, donor) in plate.genotypes.iter().enumerate() {
        for (g, copies) in donor.iter().enumerate() {
            // Only database alleles are visible to the exact search.
            let known: BTreeSet<String> = copies
                .iter()
                .filter_map(|c| plate.allele_name(g, c))
                .collect();
            let expected = if known.len() == 2 {
                "heterozygous"
            } else {
                "homozygous"
            };
            let key = (plate.indices[d].2.clone(), plate.genes[g].name.clone());
            assert_eq!(
                called.get(&key).map(String::as_str),
                Some(expected),
                "{key:?}"
            );
        }
    }
}

/// What `discover blast` must report for the plate, whichever aligner
/// compared the reads with the database.
fn check_discovery(plate: &Plate, result: &Frame) {
    let names = column(result, "allele_name");
    let cases = column(result, "case");
    let sequences = column(result, "aln_qseq");
    let new: Vec<usize> = (0..result.nrows())
        .filter(|r| names[*r].contains("_S"))
        .collect();
    // The new allele in its three carriers, and nothing else that is new.
    let carriers: BTreeSet<&str> = new.iter().map(|r| cases[*r].as_str()).collect();
    let expected: BTreeSet<&str> = NOVEL_CARRIERS
        .iter()
        .map(|d| plate.indices[*d].2.as_str())
        .collect();
    assert_eq!(carriers, expected);
    assert_eq!(new.len(), NOVEL_CARRIERS.len());
    let nearest = format!("{}*01_S", plate.genes[NOVEL_GENE].name);
    for row in &new {
        assert_eq!(sequences[*row], plate.novel);
        assert!(names[*row].starts_with(&nearest), "{}", names[*row]);
    }
    // Every database allele of every donor is reported under its own name.
    let reported: BTreeSet<(&str, &str)> = (0..result.nrows())
        .map(|r| (cases[r].as_str(), names[r].as_str()))
        .collect();
    for (d, donor) in plate.genotypes.iter().enumerate() {
        for (g, copies) in donor.iter().enumerate() {
            for name in copies.iter().filter_map(|c| plate.allele_name(g, c)) {
                let key = (plate.indices[d].2.as_str(), name.as_str());
                assert!(reported.contains(&key), "{key:?} is missing");
            }
        }
    }

    // The evidence columns of the new allele.
    let number = |name: &str, row: usize| -> f64 {
        column(result, name)[row]
            .parse()
            .unwrap_or_else(|_| panic!("{name} is empty in row {row}"))
    };
    for row in &new {
        assert_eq!(number("db_distance", *row), 1.0);
        assert_eq!(number("carriers", *row), NOVEL_CARRIERS.len() as f64);
        let copies = number("copies", *row);
        assert!((0.7..1.4).contains(&copies), "copies {copies}");
        assert!(number("evidence", *row) > 0.0);
        assert_eq!(column(result, "segregates")[*row], "true");
    }
}

#[test]
fn the_allele_missing_from_the_database_is_discovered() {
    if !blastn_available() {
        eprintln!("skipped: blastn is not on the PATH");
        return;
    }
    let (dir, plate) = demultiplexed();
    let discover = [
        "discover",
        "blast",
        "demux.tsv.gz",
        "db.fasta",
        "novel.tsv",
        "-g",
        "V",
    ];
    let first = run(dir.path(), &discover);
    assert!(!printed(&first).contains("Skipping BLASTn"));
    let result_path = dir.path().join("novel.tsv.gz");
    check_discovery(&plate, &table(result_path.clone()));

    // A second run reuses the search and gives the same table.
    let before = std::fs::read(&result_path).unwrap();
    let second = run(dir.path(), &discover);
    assert!(printed(&second).contains("Skipping BLASTn"));
    assert_eq!(std::fs::read(&result_path).unwrap(), before);

    // Another database on the same reads has its own intermediate files
    // and leaves the first run's search in place.
    let database = std::fs::read_to_string(dir.path().join("db.fasta")).unwrap();
    let half: String = database
        .lines()
        .take(2 * GENES)
        .map(|l| format!("{l}\n"))
        .collect();
    std::fs::write(dir.path().join("half.fasta"), half).unwrap();
    run(
        dir.path(),
        &[
            "discover",
            "blast",
            "demux.tsv.gz",
            "half.fasta",
            "half.tsv",
            "-g",
            "V",
        ],
    );
    let run_directories = std::fs::read_dir(dir.path().join(".immunediscover/blast"))
        .unwrap()
        .count();
    assert_eq!(run_directories, 2);
    let third = run(dir.path(), &discover);
    assert!(printed(&third).contains("Skipping BLASTn"));

    // Export of the new allele: one record, found in three donors.
    run(
        dir.path(),
        &[
            "table",
            "fasta",
            "novel.tsv.gz",
            "novel.fasta",
            "--colname",
            "allele_name",
            "--colseq",
            "aln_qseq",
            "--filter",
            "_S",
            "--mincase",
            "2",
        ],
    );
    let exported = std::fs::read_to_string(dir.path().join("novel.fasta")).unwrap();
    assert_eq!(exported.matches('>').count(), 1);
    assert!(exported.replace('\n', "").ends_with(&plate.novel));
}

/// Rows of a discovery table as (case, sequence, name, reads), sorted.
fn discovered(result: &Frame) -> Vec<(String, String, String, String)> {
    let mut rows: Vec<_> = column(result, "case")
        .into_iter()
        .zip(column(result, "aln_qseq"))
        .zip(column(result, "allele_name"))
        .zip(column(result, "count"))
        .map(|(((case, sequence), name), count)| (case, sequence, name, count))
        .collect();
    rows.sort();
    rows
}

#[test]
fn the_builtin_aligner_discovers_the_same_alleles_without_blast() {
    let (dir, plate) = demultiplexed();
    let output = run(
        dir.path(),
        &[
            "discover",
            "blast",
            "demux.tsv.gz",
            "db.fasta",
            "builtin.tsv",
            "-g",
            "V",
            "--aligner",
            "builtin",
        ],
    );
    assert!(printed(&output).contains("Built-in aligner"));
    // The time the alignment took and its speed are reported.
    assert!(printed(&output).contains("reads per second"));
    // The presets are applied, but nothing of BLAST is used or reported.
    assert!(printed(&output).contains("Applying V preset length"));
    for word in ["blastn", "megablast", "preset args"] {
        assert!(!printed(&output).contains(word), "{word} in the log");
    }
    let refused = Command::new(env!("CARGO_BIN_EXE_immunediscover"))
        .args([
            "discover",
            "blast",
            "demux.tsv.gz",
            "db.fasta",
            "refused.tsv",
            "--aligner",
            "builtin",
            "--args",
            "-task blastn",
        ])
        .current_dir(dir.path())
        .output()
        .unwrap();
    assert!(!refused.status.success());
    assert!(printed(&refused).contains("--args sets options of blastn"));
    let builtin = table(dir.path().join("builtin.tsv.gz"));
    check_discovery(&plate, &builtin);

    // A stricter aligner gives the same alleles here, where every read is
    // near one; its option belongs to the built-in aligner alone.
    run(
        dir.path(),
        &[
            "discover",
            "blast",
            "demux.tsv.gz",
            "db.fasta",
            "strict.tsv",
            "-g",
            "V",
            "--aligner",
            "builtin",
            "--min-identity",
            "0.9",
        ],
    );
    let strict = table(dir.path().join("strict.tsv.gz"));
    assert_eq!(discovered(&strict), discovered(&builtin));
    let refused = Command::new(env!("CARGO_BIN_EXE_immunediscover"))
        .args([
            "discover",
            "blast",
            "demux.tsv.gz",
            "db.fasta",
            "refused.tsv",
            "--min-identity",
            "0.9",
        ])
        .current_dir(dir.path())
        .output()
        .unwrap();
    assert!(!refused.status.success());
    assert!(printed(&refused).contains("--min-identity is an option of the built-in aligner"));

    if !blastn_available() {
        eprintln!("comparison with BLAST skipped: blastn is not on the PATH");
        return;
    }
    run(
        dir.path(),
        &[
            "discover",
            "blast",
            "demux.tsv.gz",
            "db.fasta",
            "blast.tsv",
            "-g",
            "V",
        ],
    );
    let blast = table(dir.path().join("blast.tsv.gz"));
    assert_eq!(discovered(&builtin), discovered(&blast));
}

#[test]
fn a_difference_at_the_first_base_of_an_allele_is_not_clipped() {
    let (dir, plate) = demultiplexed();
    // A database whose third gene starts with another base than the
    // alleles of the plate: every allele of that gene is new, and differs
    // from the database exactly where a local alignment would stop early.
    let gene = 2;
    let mut database = String::new();
    for (g, known) in plate.genes.iter().enumerate() {
        let allele = &known.alleles[0];
        let sequence = if g == gene {
            substitute(allele, 0)
        } else {
            allele.clone()
        };
        let _ = writeln!(database, ">{}*01\n{sequence}", known.name);
    }
    std::fs::write(dir.path().join("changed.fasta"), database).unwrap();
    // No surrounding bases: the first base of the allele is the first
    // base that is compared.
    run(
        dir.path(),
        &[
            "discover",
            "blast",
            "demux.tsv.gz",
            "changed.fasta",
            "ends.tsv",
            "-g",
            "V",
            "--forward",
            "0",
            "--reverse",
            "0",
            "--aligner",
            "builtin",
        ],
    );
    let result = table(dir.path().join("ends.tsv.gz"));
    let rows = discovered(&result);
    let prefix = format!("{}*01_S", plate.genes[gene].name);
    for (d, donor) in plate.genotypes.iter().enumerate() {
        let case = &plate.indices[d].2;
        let copies: BTreeSet<&String> = donor[gene].iter().collect();
        for copy in copies {
            let found = rows
                .iter()
                .find(|(c, sequence, _, _)| c == case && sequence == copy);
            let (_, _, name, _) = found.unwrap_or_else(|| panic!("{case}: allele not found whole"));
            assert!(name.starts_with(&prefix), "{name}");
        }
    }
}

#[test]
fn align_names_the_allele_of_every_read() {
    let (dir, plate) = demultiplexed();
    run(
        dir.path(),
        &["align", "demux.tsv.gz", "db.fasta", "hits.tsv.gz"],
    );
    let hits = table(dir.path().join("hits.tsv.gz"));
    assert_eq!(hits.nrows(), plate.reads.len());
    let cases = column(&hits, "case");
    let alleles = column(&hits, "allele");
    let differences = column(&hits, "differences");
    let sequences = column(&hits, "sequence");
    let coverage = column(&hits, "coverage");
    assert!(coverage.iter().all(|c| c == "1.0"), "every allele is whole");

    // Reads per donor and allele that match it exactly.
    let mut exact: BTreeMap<(&str, &str), usize> = BTreeMap::new();
    for row in 0..hits.nrows() {
        if differences[row] == "0" {
            *exact
                .entry((cases[row].as_str(), alleles[row].as_str()))
                .or_default() += 1;
        }
    }
    let mut expected: BTreeMap<(&str, String), usize> = BTreeMap::new();
    for (d, donor) in plate.genotypes.iter().enumerate() {
        for (g, copies) in donor.iter().enumerate() {
            for name in copies.iter().filter_map(|c| plate.allele_name(g, c)) {
                *expected
                    .entry((plate.indices[d].2.as_str(), name))
                    .or_default() += READS_PER_COPY;
            }
        }
    }
    for ((case, name), reads) in &expected {
        // A read with a sequencing error may happen to equal the other
        // allele of the gene; never fewer than the error-free reads.
        let found = exact.get(&(*case, name.as_str())).copied().unwrap_or(0);
        assert!(found >= *reads, "{case} {name}: {found} of {reads}");
    }

    // Reads of the new allele: one difference from its nearest allele,
    // and the sequence that is reported is the new allele itself.
    let nearest = format!("{}*01", plate.genes[NOVEL_GENE].name);
    let carriers: BTreeSet<&str> = NOVEL_CARRIERS
        .iter()
        .map(|d| plate.indices[*d].2.as_str())
        .collect();
    let of_novel = (0..hits.nrows())
        .filter(|&row| sequences[row] == plate.novel)
        .inspect(|&row| {
            assert_eq!(alleles[row], nearest);
            assert_eq!(differences[row], "1");
            assert!(carriers.contains(cases[row].as_str()));
        })
        .count();
    assert!(of_novel >= NOVEL_CARRIERS.len() * READS_PER_COPY);

    // The same matches in the layout of BLAST.
    run(
        dir.path(),
        &[
            "align",
            "demux.tsv.gz",
            "db.fasta",
            "hits.blast",
            "--format",
            "blast",
        ],
    );
    let text = std::fs::read_to_string(dir.path().join("hits.blast")).unwrap();
    assert_eq!(text.lines().count(), plate.reads.len());
    assert!(text.lines().all(|line| line.split('\t').count() == 12));
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'A',
        })
        .collect()
}

/// `text` with `insert` written over it from `at`.
fn written_over(text: &str, at: usize, insert: &str) -> String {
    format!("{}{insert}{}", &text[..at], &text[at + insert.len()..])
}

/// A small genome and candidates of every kind, written to `dir`.
///
/// Chromosome 14 holds one gene as written and one turned around;
/// chromosome 15 holds a distant relative of a third gene. `headers` are
/// the header lines of chromosomes 1, 14 and 15.
fn genome_and_candidates(dir: &Path, headers: [&str; 3]) -> [String; 3] {
    let mut generator = Generator(0x6E0E_5EA2_C4A7_11D3);
    let genes: [String; 3] = std::array::from_fn(|_| generator.sequence(300));
    let mut relative = genes[2].clone();
    for at in (7..300).step_by(15) {
        relative = substitute(&relative, at);
    }
    let chromosome_1 = generator.sequence(30_000);
    let chromosome_14 = written_over(&generator.sequence(40_000), 10_000, &genes[0]);
    let chromosome_14 = written_over(&chromosome_14, 25_000, &reverse_complement(&genes[1]));
    let chromosome_15 = written_over(&generator.sequence(30_000), 5_000, &relative);
    // Lower case marks repeats in some genome files and means nothing.
    let chromosome_14 = written_over(
        &chromosome_14,
        10_050,
        &chromosome_14[10_050..10_150].to_lowercase(),
    );
    let mut genome = String::new();
    for (header, bases) in headers
        .iter()
        .zip([&chromosome_1, &chromosome_14, &chromosome_15])
    {
        let _ = writeln!(genome, ">{header}");
        for line in bases.as_bytes().chunks(60) {
            let _ = writeln!(genome, "{}", std::str::from_utf8(line).unwrap());
        }
    }
    std::fs::write(dir.join("genome.fa"), genome).unwrap();

    let new_allele = substitute(&substitute(&genes[0], 120), 200);
    let candidates = [
        ("known", genes[0].clone()),
        ("new", new_allele),
        ("turned", genes[1].clone()),
        ("elsewhere", genes[2].clone()),
        ("nowhere", generator.sequence(300)),
    ];
    let mut text = String::from("allele_name\taln_qseq\n");
    for (name, sequence) in &candidates {
        let _ = writeln!(text, "{name}\t{sequence}");
    }
    std::fs::write(dir.join("candidates.tsv"), text).unwrap();
    genes
}

/// What `search genome` must report for the genome above, whose
/// chromosome 14 record is called `name`.
fn check_genome_result(dir: &Path, output: &str, name: &str, genes: &[String; 3]) {
    let kept = table(dir.join(format!("{output}.gz")));
    assert_eq!(column(&kept, "allele_name"), ["known", "new", "turned"]);
    assert_eq!(
        column(&kept, "position"),
        [
            format!("{name}:10000"),
            format!("{name}:10000"),
            format!("{name}:25000")
        ]
    );
    assert_eq!(column(&kept, "edit_distance"), ["0", "2", "0"]);
    assert_eq!(column(&kept, "orientation"), ["+", "+", "-"]);
    assert_eq!(column(&kept, "cigar"), ["300=", "120=1X79=1X99=", "300="]);
    // The genome at the place, written like the candidate and in upper case.
    assert_eq!(
        column(&kept, "ref_sequence"),
        [genes[0].clone(), genes[0].clone(), genes[1].clone()]
    );

    let discarded = table(dir.join(format!("{output}.discarded.tsv")));
    assert_eq!(column(&discarded, "name"), ["elsewhere", "nowhere"]);
    assert_eq!(column(&discarded, "edit_distance")[0], "20");
    assert!(column(&discarded, "position")[0].ends_with(":5000"));
    assert_eq!(column(&discarded, "description")[1], "no place near enough");
}

#[test]
fn candidates_are_kept_when_they_lie_on_the_right_chromosome() {
    let dir = tempfile::tempdir().unwrap();
    // Names as UCSC writes them.
    let genes = genome_and_candidates(dir.path(), ["chr1", "chr14", "chr15"]);
    let output = run(
        dir.path(),
        &[
            "search",
            "genome",
            "candidates.tsv",
            "kept.tsv",
            "genome.fa",
        ],
    );
    assert!(printed(&output).contains("5 sequences processed, 3 lie on chromosome 14"));
    assert!(printed(&output).contains("bases per second"));
    check_genome_result(dir.path(), "kept.tsv", "chr14", &genes);

    // The same genome gzipped, under the name of the command this one
    // replaces, and asked for chromosome 15: only its candidate is kept.
    let text = std::fs::read(dir.path().join("genome.fa")).unwrap();
    let mut encoder = flate2::write::GzEncoder::new(
        std::fs::File::create(dir.path().join("genome.fa.gz")).unwrap(),
        flate2::Compression::fast(),
    );
    std::io::Write::write_all(&mut encoder, &text).unwrap();
    encoder.finish().unwrap();
    run(
        dir.path(),
        &[
            "search",
            "bwa",
            "candidates.tsv",
            "fifteen.tsv",
            "genome.fa.gz",
            "-c",
            "chr15",
        ],
    );
    let kept = table(dir.path().join("fifteen.tsv.gz"));
    assert_eq!(column(&kept, "allele_name"), ["elsewhere"]);
    assert_eq!(column(&kept, "edit_distance"), ["20"]);
    // Where the other candidates lie is told in the list of the dropped.
    let discarded = table(dir.path().join("fifteen.tsv.discarded.tsv"));
    assert_eq!(
        column(&discarded, "name"),
        ["known", "new", "turned", "nowhere"]
    );
    assert_eq!(column(&discarded, "position")[0], "chr14:10000");

    // A chromosome the genome does not have is an error, not an empty table.
    let refused = Command::new(env!("CARGO_BIN_EXE_immunediscover"))
        .args([
            "search",
            "genome",
            "candidates.tsv",
            "none.tsv",
            "genome.fa",
        ])
        .args(["-c", "22"])
        .current_dir(dir.path())
        .output()
        .unwrap();
    assert!(!refused.status.success());
    assert!(printed(&refused).contains("No record of the genome belongs to chromosome 22"));
}

#[test]
fn genomes_with_descriptions_and_from_a_stream_give_the_same_answer() {
    let dir = tempfile::tempdir().unwrap();
    // Headers as NCBI writes them; chromosome 1 must not pass for 14.
    let genes = genome_and_candidates(
        dir.path(),
        [
            "NC_000001.11 Homo sapiens chromosome 1, GRCh38.p14 Primary Assembly",
            "NC_000014.9 Homo sapiens chromosome 14, GRCh38.p14 Primary Assembly",
            "NC_000015.10 Homo sapiens chromosome 15, GRCh38.p14 Primary Assembly",
        ],
    );
    run(
        dir.path(),
        &[
            "search",
            "genome",
            "candidates.tsv",
            "kept.tsv",
            "genome.fa",
            "-c",
            "chromosome 14",
            "--tag",
            "Primary Assembly|alternate locus",
        ],
    );
    check_genome_result(dir.path(), "kept.tsv", "NC_000014.9", &genes);

    // From the standard input.
    let genome = std::fs::File::open(dir.path().join("genome.fa")).unwrap();
    let piped = Command::new(env!("CARGO_BIN_EXE_immunediscover"))
        .args(["search", "genome", "candidates.tsv", "piped.tsv", "-"])
        .current_dir(dir.path())
        .stdin(genome)
        .output()
        .unwrap();
    assert!(piped.status.success(), "{}", printed(&piped));
    check_genome_result(dir.path(), "piped.tsv", "NC_000014.9", &genes);

    // From an address, as a genome on a server is read. `curl` fetches
    // it; a file address needs no network.
    let has_curl = Command::new("curl")
        .arg("--version")
        .output()
        .is_ok_and(|out| out.status.success());
    if !has_curl {
        eprintln!("reading from an address skipped: curl is not on the PATH");
        return;
    }
    let address = format!("file://{}", dir.path().join("genome.fa").display());
    let output = run(
        dir.path(),
        &[
            "search",
            "genome",
            "candidates.tsv",
            "fetched.tsv",
            &address,
        ],
    );
    assert!(printed(&output).contains("not stored"));
    check_genome_result(dir.path(), "fetched.tsv", "NC_000014.9", &genes);
    // A download that fails is said to have failed.
    let missing = format!("file://{}", dir.path().join("missing.fa").display());
    let failed = Command::new(env!("CARGO_BIN_EXE_immunediscover"))
        .args(["search", "genome", "candidates.tsv", "failed.tsv", &missing])
        .current_dir(dir.path())
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert!(printed(&failed).contains("curl"), "{}", printed(&failed));
}
