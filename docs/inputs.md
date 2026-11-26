# Input files

All tables are tab-separated text (TSV) with a header line. Files ending
in `.gz` are read and written compressed; most commands add `.gz` to the
output name themselves.

## Reads: FASTQ

Used by: `preprocess demultiplex`.

Single-end reads of amplicons from **genomic DNA**, one library per plate.
Every read starts with the index of its well (forward) and ends with the
second index (reverse). Plain or gzipped FASTQ.

## Index sheet: which indices belong to which donor

Used by: `preprocess demultiplex`.

Exactly three columns, in this order:

| Column | Content |
|---|---|
| `forward_index` | Index sequence at the start of the read. |
| `reverse_index` | Index sequence at the end of the read. |
| `case` | Donor identifier. Free text; it is carried into every later table. |

```
forward_index	reverse_index	case
TTGGAAATAG	CAACGGGCAC	D1
ACCTGGTTCA	GGTACATTGC	D2
```

The row number becomes the `well` (the first donor row is well 1). Give
every well its own row even if two wells hold the same donor: results are
kept per well and donor.

## Demultiplexed reads

Written by `preprocess demultiplex`; used by `search exact`,
`search heptamer`, `discover blast`, `discover hsmm`.

| Column | Content |
|---|---|
| `well` | Well number (row of the index sheet). |
| `case` | Donor identifier. |
| `name` | Read name from the FASTQ. |
| `genomic_sequence` | The sequence of the read. |

You can also build this table yourself from already separated samples, as
long as it has these four columns.

## Allele database: FASTA

Used by: `search exact`, `search heptamer`, `discover blast`,
`discover hsmm`.

One record per known allele. The name must have the form `GENE*ALLELE`,
for example `IGHV1-69*01`: everything before the `*` is taken as the gene,
and alleles of one gene are compared with each other.

```
>IGHV1-69*01
CAGGTGCAGCTGGTGCAGTCTGGGGCT...
>IGHV1-69*02
CAGGTCCAGCTGGTGCAATCTGGGGCT...
```

- Use one database per gene type (V, D, J) and run the commands once for
  each.
- Sequences should cover the same region as your amplicon (the coding
  region without the recombination signal).
- A name containing `_S` and digits marks an allele that was itself found
  by this program (see `fasta hash`).

## Pseudogenes: FASTA (optional)

Used by: `discover blast --pseudo`.

Same format as the database. Reads that resemble one of these more than
any allele are discarded.

## Per-allele and per-gene limits (optional)

Used by: `search exact --expect` and `--deletion`.

Two columns, `name` and `ratio`:

```
name	ratio
IGHV1-2*02	0.6
```

- `--expect`: allele names with the minimum `allele_freq` each must reach.
- `--deletion`: gene names with a `gene_case_freq` limit for each.

## Heptamer list: JSON

Used by: `search heptamer` (default file name `heptamers.json`).

The accepted heptamer sequences per chain:

```
{"IGHV": ["CACAGTG", "CACAATG", "CACAGAG"],
 "IGKV": ["CACAGTG", "CACACTG"],
 "IGLV": ["CACAGTG", "CACGGTG"]}
```

## Genome: FASTA with a BWA index

Used by: `search bwa`.

One or more genome FASTA files, each indexed with `bwa index`. The record
descriptions must name the chromosome (for example "chromosome 14"),
because candidates are kept by the chromosome they map to.

## Tables for `analyze`

- `analyze cooccurrence`: any table with a donor column and an allele
  column (defaults `case` and `db_name`; for `discover blast` output use
  `-A allele_name`).
- `analyze haplotype`: the output of `search exact`, or any table with
  donor, allele, gene and `count` columns.
