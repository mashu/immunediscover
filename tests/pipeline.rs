//! End-to-end runs of the command-line program on a small synthetic plate.
//!
//! Nothing here comes from real samples: genes, alleles, donors and reads
//! are generated from a fixed seed. The plate is small enough for the whole
//! file to run in a few seconds.
//!
//! `discover blast` needs BLAST+. Without `blastn` on the `PATH` its test is
//! skipped, unless `IMMUNEDISCOVER_REQUIRE_BLAST` is set (as in CI), in
//! which case a missing `blastn` fails the test.

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
    let result = table(result_path.clone());

    let names = column(&result, "allele_name");
    let cases = column(&result, "case");
    let sequences = column(&result, "aln_qseq");
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
        column(&result, name)[row]
            .parse()
            .unwrap_or_else(|_| panic!("{name} is empty in row {row}"))
    };
    for row in &new {
        assert_eq!(number("db_distance", *row), 1.0);
        assert_eq!(number("carriers", *row), NOVEL_CARRIERS.len() as f64);
        let copies = number("copies", *row);
        assert!((0.7..1.4).contains(&copies), "copies {copies}");
        assert!(number("evidence", *row) > 0.0);
        assert_eq!(column(&result, "segregates")[*row], "true");
    }

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
