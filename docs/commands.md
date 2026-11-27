# Commands, options and output columns

Every command follows the same pattern:

```
immunediscover <group> <command> [options] <input files> <output file>
```

`align` stands on its own, without a group. `--help` after any command
lists its options. Tables are tab-separated; `.gz` output is compressed.

Contents: [preprocess demultiplex](#preprocess-demultiplex) ·
[search exact](#search-exact) · [search heptamer](#search-heptamer) ·
[search bwa](#search-bwa) · [align](#align) ·
[discover blast](#discover-blast) ·
[discover hsmm](#discover-hsmm) ·
[analyze haplotype](#analyze-haplotype) ·
[analyze cooccurrence](#analyze-cooccurrence) · [table](#table) ·
[fasta](#fasta)

Columns that mean the same everywhere:

| Column | Meaning |
|---|---|
| `well` | Well of the plate (row number of the index sheet). |
| `case` | Donor identifier. |
| `db_name` | Name of a database allele (`GENE*ALLELE`). |
| `gene` | The part of the allele name before `*`. |
| `count` | Reads supporting the sequence in this donor. |
| `ratio` | `count` relative to the best allele of the same gene in this donor. |

---

## preprocess demultiplex

Sorts the reads of a plate into donors by the index sequences at the two
ends of every read.

```
immunediscover preprocess demultiplex plate.fastq.gz indices.tsv demux.tsv
```

| Option | Default | What it is for |
|---|---|---|
| `-l`, `--length` | 200 | Reads shorter than this are dropped. Set it a little below your amplicon length. |
| `-s`, `--split` | off | Also write one FASTQ file per donor. |
| `-f`, `--forwardarrayindex` | none | Name of an additional forward array index, if your library design has one. |
| `--case-filter-regex` | none | Keep only donors whose identifier matches this pattern. |

**Output** `demux.tsv.gz`: one row per read.

| Column | Meaning |
|---|---|
| `well`, `case` | Where the read belongs. |
| `name` | Read name from the FASTQ. |
| `genomic_sequence` | The sequence of the read. |

**Output** `demux.tsv.log`: one row per well with `count` (reads),
`mean_length` and `std_length` (read length). Check it for empty or weak
wells.

---

## search exact

Finds the database alleles that occur **exactly** in the reads and reports,
per donor, how many reads support each, together with the bases next to
the allele (the recombination signal). Use it to genotype donors with a
database you trust.

```
immunediscover search exact demux.tsv.gz IGHV.fasta exact.tsv -g V
```

### Options

| Option | Default | What it is for |
|---|---|---|
| `-g`, `--gene` | V | Gene type: `V`, `D` or `J`. Decides on which side of the allele the recombination signal is expected. |
| `-c`, `--mincount` | 5 | Minimum reads for an allele in a donor. |
| `-f`, `--minratio` | 0.1 | Minimum reads relative to the best allele of the same gene in the donor. Removes weak second alleles that are errors. |
| `--min-allele-mratio` | 0.05 | Minimum for the allele's reads relative to what this allele usually has across donors. Removes alleles present in a trace amount compared with other donors (cross-contamination). |
| `--min-gene-mratio` | 0.05 | The same for the whole gene. |
| `--rss` | heptamer | Which parts of the recombination signal to report next to the allele: any of `heptamer`, `spacer`, `nonamer`, comma-separated. |
| `--extension` | none | Report a fixed number of bases on the signal side instead of the signal parts. |
| `-a`, `--affix` | 13 | Bases to report on the other side of the allele. |
| `-t`, `--top` | 1 | How many versions of the neighbouring bases to report per allele. |
| `--border` | 0 | Bases at each end of the read that the reported neighbourhood must not reach into; reads where it does are rejected. |
| `--adjust-per-gene-extension` | off | Shorten the extension per gene automatically so that it stays clear of the border. |
| `--adjust-percent` | 1 | With the option above: share of a gene's reads that must be clear of the border. |
| `-r`, `--refgene` | none | Reference genes; adds columns with counts relative to them (see below). |
| `-e`, `--expect` | none | Table of allele names with a minimum `allele_freq` for each ([format](inputs.md)). |
| `-d`, `--deletion` | none | Table of gene names with a `gene_case_freq` limit for each. |
| `--locus` | IG | Only genes whose name starts with this text are analysed; others (control genes) are left out. |
| `--ref-fasta` | none | A second FASTA; adds `isin_db`, saying whether each sequence is in it. |
| `--raw` | none | Also write the unfiltered result to this file, for troubleshooting. |
| `-l`, `--limit` | 0 | Use only the first N reads (0 = all). For quick tests. |
| `-n`, `--noplot` | off | Do not draw the gene plot on the screen. |

### Output

One row per allele found in a donor.

| Column | Meaning |
|---|---|
| `well`, `case`, `db_name`, `gene` | Donor and allele. |
| `sequence` | The allele sequence. |
| `count` | Reads containing the allele, whatever is next to it. |
| `ratio` | `count` relative to the best allele of the gene in this donor. |
| `full_count` | Reads containing the allele with exactly the neighbouring bases shown in this row. |
| `full_ratio` | `full_count` relative to the best of the gene. |
| `flank_index` | Rank of this version of the neighbouring bases (1 = most reads). |
| `allele_freq` | Share of the gene's reads that belong to this allele (about 0.5 for a heterozygous donor, 1 for a homozygous one). |
| `gene_count` | Reads of the gene in this donor. |
| `case_count` | Reads of all genes in this donor. |
| `allele_case_freq` | `count` as a share of all the donor's reads. |
| `gene_case_freq` | `gene_count` as a share of all the donor's reads. A gene far below its usual share may be deleted on one chromosome. |
| `cross_case_median_count`, `cross_case_median_allele_count`, `cross_case_median_gene_count` | What is usual across donors (medians) for this allele and gene. |
| `allele_to_cross_case_median_ratio` | This donor's reads of the allele relative to the usual. Far below 1 suggests contamination from another well. |
| `gene_to_cross_case_median_ratio` | The same for the gene. |

Neighbouring bases, depending on `-g` and `--rss`:

| Gene type | Columns |
|---|---|
| V | `prefix` (bases before the allele), then `heptamer`, `spacer`, `nonamer` after it |
| J | `nonamer`, `spacer`, `heptamer` before the allele, then `suffix` (bases after it) |
| D | `pre_nonamer`, `pre_spacer`, `pre_heptamer` before; `post_heptamer`, `post_spacer`, `post_nonamer` after |
| any, with `--extension` | `prefix`, `suffix`, `prefix_len`, `suffix_len` |

Optional columns: `isin_db` (with `--ref-fasta`);
`count_<gene>_ratio`, `ref_gene_count`, `ref_gene_count_<gene>_ratio`
(with `--refgene`): counts relative to the reference gene.

---

## search heptamer

For V genes: finds where each allele ends by looking for the heptamer
that follows it, and reports alleles whose real end differs from the
database sequence.

```
immunediscover search heptamer demux.tsv.gz IGHV.fasta hept.tsv summary.tsv
```

| Option | Default | What it is for |
|---|---|---|
| `-c`, `--chain` | IGHV | `IGHV`, `IGKV` or `IGLV`: which heptamer list to use. |
| `-j`, `--json` | heptamers.json | File with the accepted heptamers ([format](inputs.md)). |
| `-d`, `--maxdist` | 1 | Differences allowed between a heptamer in the read and the listed ones. |
| `-b`, `--begin` | 0 | Bases to ignore at the start of each database allele when searching. |
| `-e`, `--end` | 8 | Bases to ignore at the end of each database allele; the end is then rebuilt from the read up to the heptamer. |
| `-m`, `--mincount` | 1 | Minimum reads for a row of the summary. |
| `-r`, `--ratio` | 0.25 | Minimum reads of the rebuilt sequence relative to reads of the shortened one. |

**Output** `hept.tsv.gz`, one row per read with a match: `db_name`,
`full_match` and `trimmed_match` (positions of the whole and of the
shortened allele in the read), `db_length`, `full_length` (length up to
the heptamer), `sequence` (the allele as found, up to the heptamer),
`heptamer`.

**Output** `summary.tsv`, one row per distinct sequence: `db_name`,
`sequence`, `heptamer`, `count`, `ratio`.

---

## search bwa

Maps candidate sequences to a genome and keeps those that land on the
expected chromosome. A check that new alleles come from the right locus.

```
immunediscover search bwa novel.tsv.gz novel_chr14.tsv genome.fasta
```

| Option | Default | What it is for |
|---|---|---|
| `-c`, `--chromosome` | chromosome 14 | Text that the description of the genome record must contain. Use "chromosome 2" for kappa, "chromosome 22" for lambda. |
| `-n`, `--colname` | allele_name | Column with the names. |
| `-s`, `--colseq` | aln_qseq | Column(s) with the sequence; several are joined. |
| `-t`, `--tag` | primary assembly or alternate locus | Pattern that the record description must match. |

**Output**: the input rows that map to the chromosome, with added columns
`position` (where), `edit_distance` (differences from the genome),
`ref_sequence` (the genome sequence there), `orientation` (`+` or `-`),
`cigar` (alignment).

---

## align

Writes, for every read, the allele of a database it resembles most, with
the number of differences, the position of the match in the read and the
matching part of the read. No other software is needed. Options, columns
and method are described in [The built-in aligner](align.md).

```
immunediscover align demux.tsv.gz IGHV.fasta hits.tsv.gz
```

---

## discover blast

Finds sequences that resemble database alleles without being identical:
candidates for new alleles. Described in full in
[discover-blast.md](discover-blast.md); the evidence columns in
[evidence.md](evidence.md).

```
immunediscover discover blast demux.tsv.gz IGHV.fasta novel.tsv -g V
```

Reads are compared with the alleles by BLAST+, or with
`--aligner builtin` by immunediscover itself.

---

## discover hsmm

Finds D genes without relying on their sequence: it looks for the pattern
"recombination signal, short gene, recombination signal" in the reads.
The pattern is learned from the known D alleles present in your donors.
Use it when D alleles may differ too much from the database for
`discover blast`.

```
immunediscover discover hsmm demux.tsv.gz IGHD.fasta d_found.tsv.gz
```

| Option | Default | What it is for |
|---|---|---|
| `-c`, `--mincount` | 10 | Minimum reads for a known D allele to be used for learning. |
| `-r`, `--ratio` | 0.2 | Minimum ratio to the gene's best allele for a known D allele to be used for learning. |
| `--min-posterior` | 0.7 | Minimum confidence (0 to 1) of a detection. |
| `--min-gene-len`, `--max-gene-len` | 10, 70 | Shortest and longest D gene considered. |
| `--min-heptamer-prob-pre`, `--min-heptamer-prob-post` | 0.05 | Minimum quality of the heptamer before and after the gene (0 switches the check off). |
| `--out-mincount` | 10 | Minimum reads for a found D to be reported. |
| `--out-minratio` | 0.2 | Minimum ratio to the gene's best allele for a found D to be reported. |
| `-l`, `--limit` | 0 | Use only the first N reads. |

**Output**, one row per D sequence per donor:

| Column | Meaning |
|---|---|
| `well`, `case` | Donor. |
| `sequence` | The D gene found. |
| `count`, `ratio` | Reads, and reads relative to the best allele of the gene. |
| `allele_name` | Database name, or nearest allele plus `_S<digits>` for a new one. |
| `isin_db` | Whether the sequence is in the database. |
| `db_name` | Database name if it is. |
| `nearest_db`, `nearest_db_dist` | Nearest database allele and the differences from it. |
| `gene` | Gene of the nearest allele. |
| `pre_nonamer`, `pre_spacer`, `pre_heptamer`, `post_heptamer`, `post_spacer`, `post_nonamer` | The recombination signals found on both sides. |
| `posterior_prob` | Confidence of the detection, 0 to 1. |
| `heptamer_prob_pre`, `heptamer_prob_post` | How typical the two heptamers are, 0 to 1. |
| `heptamer_logp_pre`, `heptamer_logp_post`, `log_path_prob`, `log_total_prob` | The same quantities on a log scale, for diagnostics. |

---

## analyze haplotype

Decides per donor and gene whether the donor is homozygous or
heterozygous, from a table of alleles with counts (normally the output of
`search exact`).

```
immunediscover analyze haplotype exact.tsv.gz genotypes.tsv
```

| Option | Default | What it is for |
|---|---|---|
| `-C`, `--case-col` | case | Column with the donor. |
| `-A`, `--allele-col` | db_name | Column with the allele name. |
| `-G`, `--gene-col` | gene | Column with the gene. |
| `-c`, `--mincount` | 5 | Alleles with fewer reads are ignored. |
| `-r`, `--min-ratio` | 0.1 | Minimum reads of the second allele relative to the first for the donor to be called heterozygous. |
| `-f`, `--novel-fasta` | none | FASTA of new alleles; marks them in the result. |

**Output**, one row per donor and gene:

| Column | Meaning |
|---|---|
| `case`, `gene` | Donor and gene. |
| `genotype` | `homozygous`, `heterozygous` or `uncertain`. |
| `allele_1`, `allele_2` | The two alleles, most reads first (`allele_2` empty if homozygous). |
| `count_1`, `count_2` | Their reads. |
| `freq_1`, `freq_2` | Their shares of the gene's reads. |
| `ratio` | `count_2` / `count_1`. |
| `total_count` | Reads of the gene in the donor. |
| `other_alleles` | Further alleles seen, if more than two. A sign of a duplicated gene or of contamination. |
| `novel_1`, `novel_2` | Whether the alleles are in `--novel-fasta`. |

---

## analyze cooccurrence

Finds alleles that tend to be carried by the same donors (for example
alleles of neighbouring genes inherited together).

```
immunediscover analyze cooccurrence exact.tsv.gz --clusters groups.tsv
```

| Option | Default | What it is for |
|---|---|---|
| `-C`, `--case-col` | case | Column with the donor. |
| `-A`, `--allele-col` | db_name | Column with the allele name. |
| `-m`, `--min-donors` | 2 | Alleles in fewer donors are left out. |
| `--cluster-method` | components | How alleles are grouped: `components`, `complete`, `average` or `single`. |
| `--cluster-threshold` | 0.5 | How strongly alleles must go together (0 to 1) to be grouped. |
| `--min-cluster-size` | 3 | Smallest group reported. |
| `--clusters` | none | File for the groups. |
| `--debug-triangles` | off | Print diagnostics. |

**Output** `<input>_edges.tsv`, one row per pair of alleles: `allele_a`,
`allele_b`, `rho` (correlation between the two across donors, -1 to 1),
`jaccard` (donors with both / donors with either), `support` (donors with
both), `p_value`.

**Output** `--clusters` file, one row per allele: `group_id`, `allele`,
`donors` (list), `n_donors`.

---

## table

Small tools for tables, so that results can be combined without a
spreadsheet. All read and write TSV (`.gz` allowed).

| Command | What it does | Main options |
|---|---|---|
| `table filter` | Keep rows by text pattern or number. | `-c` column; `--pattern` for text; `--operator` (`<`, `<=`, `>=`, `>`) with `--threshold` for numbers |
| `table select` | Keep some columns. | `-c` columns, comma-separated |
| `table sort` | Sort rows. | `-c` columns; `-r` descending |
| `table unique` | Keep one row per distinct combination. | `-c` columns |
| `table aggregate` | Count rows per group. | `-g` grouping columns; `-k` columns to keep; `-c` name of the count column |
| `table transform` | Rewrite the text of a column with a pattern. | `-c` column; `-p` pattern; `-r` replacement; `--new-column` |
| `table leftjoin`, `table outerjoin` | Combine two tables by shared columns. `leftjoin` keeps all rows of the first table, `outerjoin` all rows of both. | `-k` key columns; `--left-keys`, `--right-keys`; `--left-prefix`, `--right-prefix`; `--left-select`, `--right-select` |
| `table collect` | Stack many tables with the same columns (for example one per plate). | a file pattern such as `"plates/*.tsv.gz"` |
| `table fasta` | Write sequences of a table to FASTA. | `-n` name column; `-s` sequence column(s); `-f` keep names matching a pattern; `-c` remove a pattern from names; `--mincase` minimum donors; `--unique-sequences`; `-d` description column |
| `table exclude` | Remove rows whose name or sequence is in a FASTA. | `-n` name column; `-s` sequence column |

Examples:

```
# calls with at least one copy
immunediscover table filter novel.tsv.gz strong.tsv -c copies --operator ">=" --threshold 0.8

# only sequences not in the database
immunediscover table filter novel.tsv.gz new.tsv -c allele_name --pattern "_S"

# in how many donors was each sequence found?
immunediscover table aggregate new.tsv.gz per_allele.tsv -g allele_name -k gene -c n_donors

# several plates into one table
immunediscover table collect "plates/*_novel.tsv.gz" all_novel.tsv
```

---

## fasta

| Command | What it does |
|---|---|
| `fasta merge <out> <in1> <in2> ...` | Joins FASTA files, keeping every distinct sequence once. `--prefer-last` decides which name wins when two files name the same sequence differently; `-c` removes a pattern from names; `--add-source-prefix` puts the file name in front of each name. |
| `fasta diff <a> <b>` | Compares two FASTA files by sequence and prints what they share and what is only in one. |
| `fasta hash <in>` | Prints the FASTA with `_S<digits>` added to every name. The digits are derived from the sequence, so the same sequence always gets the same suffix. |
