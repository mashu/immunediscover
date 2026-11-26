# immunediscover (Rust)

A Rust port of [immunediscover.jl](https://github.com/mashu/immunediscover.jl)
(commit `b89d5cc`): germline immune receptor allele discovery from plate
based NGS libraries.

The port keeps the command line and the file formats of the Julia tool.

- **0.1.0** reproduced the Julia tool exactly, bugs included (148 of 163
  compared outputs byte-identical, the rest explained).
- **0.2.0** fixed the upstream bugs found during the port.
- **0.3.0** made read counts mean the same thing in every command:
  `full_count`, `count`, `full_ratio`, `ratio`, `flank_index` and `--top`
  behave identically in `search exact` and `discover blast`; added
  `--extension-maxdist`.
- **0.4.0** (this version): alleles without reads take the context of their
  gene; `search bwa` ignores letter case; `--refgene` ratios are empty
  where the reference is absent.

[BUGS.md](BUGS.md) lists every fix, the choices made where a fix was not
obvious, and what changes in the output. No known issue is open.

## Build

```
cargo build --release          # binary: target/release/immunediscover
cargo build --release --bin ki-check
cargo test                     # unit + integration tests
```

Rust 1.82 or newer. External programs, looked up on `PATH` when the
corresponding command runs:

| Command | Needs |
|---|---|
| `discover blast` | `blastn`, `makeblastdb` (BLAST+ 2.15–2.17) |
| `search bwa` | `bwa` (0.7.x) and a `bwa index` of the genome |
| `ki-check` | same as `discover blast` |

## KI allele recovery self-check

`ki-check` runs `discover blast` on the KI demultiplex tables and checks
that every KI `_S` target allele is recovered (matching logic of
`ImmuneDiscover-IGH-2026-Check/check_target_alleles.jl`), while counting
extra `_S` calls as false positives.

```
cargo build --release --bin ki-check

# Full run (discovery + recovery report):
./target/release/ki-check \
  --check-root /path/to/ImmuneDiscover-IGH-2026-Check

# Re-score existing BLAST tables only:
./target/release/ki-check \
  --check-root /path/to/ImmuneDiscover-IGH-2026-Check \
  --out-dir /path/to/ImmuneDiscover-IGH-2026-Check/analysis_results \
  --skip-discover
```

Defaults write under `{check-root}/analysis_results_rs/`. Exit status is
non-zero if any target is missed; pass `--max-fp-rate N` to also fail when
the false-positive rate exceeds N percent.

## Commands

Same tree, option names and defaults as upstream (`--help` at every level):

```
immunediscover preprocess demultiplex <fastq> <indices.tsv> <out.tsv>
immunediscover search     exact | heptamer | bwa
immunediscover discover   blast | hsmm
immunediscover analyze    cooccurrence | haplotype
immunediscover table      outerjoin | leftjoin | transform | aggregate | unique
                          | sort | filter | select | fasta | collect | exclude
immunediscover fasta      merge | diff | hash
```

Typical run:

```
immunediscover preprocess demultiplex lib.fastq.gz indices.tsv demux.tsv
immunediscover search exact  demux.tsv.gz IGHV.fasta exact.tsv -g V
immunediscover discover blast demux.tsv.gz IGHV.fasta novel.tsv -g V
immunediscover analyze haplotype exact.tsv.gz haplotype.tsv
```

Every invocation is appended to `immunediscover.log` in the working
directory, as upstream.

`discover blast` warns when alleles clear `-c`/`-f` on allele-level
`count`/`ratio` but every extension fails the same cut on `full_count` /
`full_ratio` (support split across flanks). A small
`*.flank-split.tsv.gz` sidecar lists them. With `-v`, a full
`*.unfiltered.tsv.gz` is also written (clusters before abundance cuts).
The final filtered table alone cannot show these losses.

## Layout

The crate is a library plus a thin binary. Dependencies point downwards
only: `cli → commands → {frame, io, seq, align, stats, filters}`.

```
src/
  main.rs, lib.rs, error.rs, logging.rs
  cli/        clap definitions only: one file per command group, BLAST
              presets, validators, invocation log. Converts arguments into
              the commands' *Params structs; contains no analysis code.
  commands/   one module per sub-command, each exposing `Params` + `run`
    preprocess/  demultiplex, pattern (barcode matching), simulate
    search/      exact/{scan, flanks, border, summarize, frequencies, refgene}
                 heptamer, bwa/{aligner, cigar}
    discover/    blast/{runner, workdir, locate, affixes, trim, assign}
                 hsmm/{model, scan}, profile
    variants.rs  the one definition of full_count / count / ratios /
                 flank_index / top, used by search exact and discover blast
    analyze/     cooccurrence/{stats, cluster}, haplotype/{flanks, genotype}
    table/       join, filter, transform, simple, to_fasta, collect, exclude
    fasta/       keyed_set, merge, diff
  frame/      small typed data frame with the DataFrames.jl semantics the
              tool depends on: column type inference, group order, stable
              sorts with missing values, join row order, float printing
  io/         TSV/CSV (CSV.jl 0.10 compatible), FASTA, FASTQ, gzip sniffing
  seq/        gene types, allele naming (`_S####` hashes), DNA helpers,
              multi-pattern search, database and demultiplex-table loading
  align/      Needleman-Wunsch (port of BioAlignments.jl, same traceback
              priority), EDNAFULL, semi-global mode, edit distance, banded
              approximate matching
  stats/      hypergeometric tail, Benjamini-Hochberg, hierarchical
              clustering (port of Clustering.jl), descriptive statistics
  filters/    `RowFilter` trait and the threshold filters combined by
              `GermlineFilter` (reports rows removed per criterion)
tests/
  reference_values.rs   replays values generated by the Julia packages
  pipeline.rs           end-to-end runs on a simulated library
  fixtures/             Julia-generated cases + the script that made them
parity/                 Julia-vs-Rust comparison harness (see below)
```

Design rules followed throughout:

- **Commands are plain functions of a parameter struct.** Nothing below
  `cli/` knows about clap, so every command can be called from tests or from
  other code.
- **One responsibility per file**; the largest source file is about 500
  lines, most are under 200.
- **Behaviour that must match Julia is isolated and named.** DataFrames/CSV
  semantics live in `frame/` and `io/`; every bug fix is one function or
  one small module, named in `BUGS.md`.
- **Traits for open-ended variation** (`RowFilter`); enums for closed sets
  (gene type, linkage, join kind); plain structs elsewhere.
- **Errors are values**: failures on user input are reported through
  `error::Error` and printed as one `ERROR:` line.

## Context extension in `discover blast`

Reads are BLASTed against database alleles extended by the read context
around them (`--forward`, `--reverse`), which is trimmed off again
afterwards.

- **The context is per allele.** Each database allele gets the consensus of
  the reads in which *it* is found.
- **An allele found in no read takes the context of its gene**: the
  consensus over all occurrences of the other alleles of that gene, so that
  all alleles of a gene are searched alike. `--no-gene-context` turns this
  off and leaves such alleles bare. The log says how many alleles had reads
  of their own, how many borrowed, and how many stayed unextended (genes of
  which no allele was seen).
- **Consensus**: per position, outwards from the allele, while more than
  half of the allele's occurrences agree.
- **`--extension-maxdist N`** (default 0) decides which reads count as
  occurrences. With 0, only reads containing the allele exactly. With
  N > 0, a read that contains *no* allele exactly is assigned to the
  allele(s) at the smallest edit distance, if that distance is at most N.
  Reads with an exact allele are never reassigned. Use it when the database
  is from another species or lacks the donors' alleles.

The tolerant search is seeded (an allele cut into N + 1 pieces must share
one piece exactly with any read within N edits) and verified with a banded
alignment, so it costs little: on the 2800-read test library, with a
database in which every allele carries four substitutions, no allele has an
exact match and none is extended at tolerance 0; at tolerance 4, 21 of 23
alleles are extended and discovery returns exactly the 53 donor–sequence
pairs obtained with the library's own database. Run time is unchanged
(0.8 s).

## Comparison with immunediscover.jl

`parity/runall.sh` generates synthetic V, D and J plate libraries and runs
73 scenarios covering every sub-command through both implementations.
Reference environment: Julia 1.13.4, CSV.jl **0.10.15** (upstream HEAD does
not run with CSV 1.x), DataFrames 1.8.3, BLAST 2.17.0, bwa 0.7.19.

Since 0.2.0 fixes bugs, its output is *meant* to differ from Julia's. The
check is therefore run per command: `runall.sh isolated` gives every Rust
command the **Julia** intermediates as input (a fresh copy per scenario), so
a difference can only come from that command itself.

Result: **119 of 163 outputs byte-identical**; the other 44 are all
accounted for by a fix:

| Outputs | Difference | Fix |
|---|---|---|
| 12 × `search exact` tables, 1 raw table | column order; rows for the control gene; values of locus rows identical, except the frequency columns of the one `--top 5` scenario | D1, A5, counts |
| 1 × `search exact --refgene` wide table | only the `ref_gene_count_*_ratio` columns | A5 |
| 9 × `discover blast` tables, affixes, extended FASTA | affix consensus and what follows from it; `count`, `ratio`, `flank_index`; one row per allele (without extension: same rows, three more columns) | A10, A3, A4, edge, coverage, counts |
| 6 × `discover blast` stdout | `makeblastdb` timestamps and paths | — |
| 4 × co-occurrence edge tables | `p_value` only | A1 |
| 2 × hierarchical cluster tables | same clusters, numbering (hash order upstream) | D3 |
| 2 × `discover hsmm` | Julia crashes on integer wells; log-probabilities differ by ≤ 4.3e-16 relative | B2 |
| `table transform`, `table fasta` (2), `fasta merge`, `fasta diff`, `search bwa` | `missing` text; names (and one export on which Julia crashes on empty cells); headers; listing order; soft clip | A12, A6, A7, D3, A2 |
| `table collect` | none: the glob also picks up a file a later scenario wrote into the Julia directory | — |

Against the planted truth of the test libraries, `discover blast -g V` went
from 5 rows (2 of 16 planted novel alleles) to 53 rows, one per
donor–allele pair (9 of 16), with no wrong sequence; `-g D` went from 27 of
35 known pairs, 9 of 13 novel pairs and 12 calls that are not in the truth
to 33, 13 and 0. BUGS.md has the full table.

`tests/reference_values.rs` still checks 4000 affix trims, 4500
dendrograms, 3000 flank extractions, 6000 border decisions and the
hypergeometric tail against values produced by the Julia packages; the one
place where a fix changes a reference value (flanks before the read start)
is corrected inside the test and asserted to be exactly that case.

To reproduce:

```
cargo build --release
export JULIA="julia --project=/path/to/immunediscover.jl"
parity/runall.sh jl          # Julia reference (slow)
parity/runall.sh isolated    # every Rust command on the Julia intermediates
parity/runall.sh rs          # the Rust pipeline end to end, then a full diff
```

See the header of `parity/runall.sh` for the environment variables.

## Speed

281 000 reads, 2 cores: demultiplex 3.9 s, exact search 1.1–1.4 s,
heptamer 2.4 s, BLAST discovery 57 s (48 s of it inside `blastn`).
