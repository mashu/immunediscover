# The built-in aligner

immunediscover can compare reads with alleles itself, without BLAST. The
aligner is made for one job: finding a known gene, with the bases around
it, somewhere inside a read, and counting how the read differs from it.

It is used in two places:

- `immunediscover align` writes the best matching allele of every read to
  a table.
- `discover blast --aligner builtin` uses it instead of BLAST to assign
  reads to alleles (see [discover blast](discover-blast.md#search)).

## `align`: the best allele of every read

```
immunediscover align demux.tsv.gz IGHV.fasta hits.tsv.gz
```

- `demux.tsv.gz`: the reads. Either the table written by
  `preprocess demultiplex`, or a FASTA file.
- `IGHV.fasta`: the sequences to look for, for example the known alleles.
- `hits.tsv.gz`: the result, one row per read that matches an allele.
  Leave out `.gz` for an uncompressed file.

Reads are compared in both orientations. A read that matches no allele
well enough is left out of the result; the log says how many reads
matched.

### Options

| Option | Default | What it is for |
|---|---|---|
| `-i`, `--min-identity` | 0.8 | Smallest share of equal bases between the read and the allele. 0.8 means at most one difference in five positions. |
| `-c`, `--min-coverage` | 0.5 | Smallest share of the allele that must lie inside the read. A read can end in the middle of a gene; with 0.5 at least half of the allele must be there. 1 keeps only reads that hold the allele from its first to its last base. |
| `-f`, `--format` | `table` | `table` for the columns below, `blast` for the twelve columns of `blastn -outfmt 6`. |
| `--match` | 2 | Score of two equal bases. |
| `--mismatch` | 3 | Cost of two different bases. |
| `--gap-open` | 5 | Cost of starting a gap (an inserted or deleted stretch). |
| `--gap-extend` | 2 | Cost of every base of a gap. |

The scores are those `blastn` uses. They decide where gaps are placed in
an alignment, and rarely need changing.

### Output columns

| Column | Meaning |
|---|---|
| `well`, `case` | Well and donor of the read. Only when the reads come from a demultiplexed table. |
| `read` | Name of the read. |
| `allele` | Name of the best matching sequence of the database. |
| `strand` | `+` if the read matches as written, `-` if its reverse complement does. |
| `identity` | Share of equal bases among the aligned positions (0 to 1). |
| `differences` | Mismatches plus bases in gaps: how far the read is from the allele. |
| `mismatches` | Positions where read and allele hold different bases. |
| `gap_bases` | Bases that the read has and the allele lacks, or the other way round. |
| `gaps` | Number of separate gaps. |
| `score` | Alignment score. |
| `read_start`, `read_end` | Where the match lies in the read as written (first base is 1). |
| `read_length` | Length of the read. |
| `allele_start`, `allele_end` | Which part of the allele is matched. 1 and the allele length when the read holds the whole allele. |
| `allele_length` | Length of the allele. |
| `coverage` | Matched share of the allele (0 to 1). |
| `cigar` | The alignment, position by position: `=` equal, `X` different, `I` base only in the read, `D` base only in the allele. `120=1X175=` is 120 equal bases, one mismatch, 175 equal bases. |
| `sequence` | The matching part of the read, written in the orientation of the allele. With `differences` 0 this is the allele itself; otherwise it is what the donor's read holds in its place. |

In the `blast` format the columns are `qseqid`, `sseqid`, `pident`,
`length`, `mismatch`, `gapopen`, `qstart`, `qend`, `sstart`, `send`,
`evalue` and `bitscore`, as BLAST writes them, without a header. Bit
scores are computed as BLAST computes them. E-values are slightly larger
than BLAST's, because BLAST shortens the sequence lengths a little before
it computes them. This format is available with the default scores only.

### Examples

```
# which allele does every read of a plate resemble most?
immunediscover align demux.tsv.gz IGHV.fasta hits.tsv.gz

# only reads that hold a whole allele, nearly unchanged
immunediscover align demux.tsv.gz IGHV.fasta hits.tsv.gz -c 1 -i 0.98

# reads in a FASTA file, result in the layout of BLAST
immunediscover align reads.fasta IGHD.fasta hits.blast --format blast
```

## How it differs from BLAST

BLAST looks for the stretch where two sequences agree best and stops
where agreement ends. A difference in the first or last bases of an
allele is therefore often cut off: the match is reported a few bases
shorter, and the difference is not counted.

The built-in aligner always compares the allele from its first to its
last base. A difference at the end of an allele is counted like one in
the middle, and the sequence that is reported has the full length of the
allele. The only place where an allele may be left incomplete is the end
of the read: if the read stops inside the gene, the rest of the allele is
not there to compare, and `coverage` tells how much was.

| | BLAST | Built-in aligner |
|---|---|---|
| Compared part of the allele | The stretch that agrees best | The whole allele, unless the read ends inside it |
| Difference at the first or last base | Often cut off | Counted |
| Short genes (D) | Need special settings | Aligned without special settings |
| Other software needed | BLAST+ 2.15 to 2.17 | None |
| Same result on every computer | Depends on the BLAST version | Yes |

## How it works

1. **Candidates.** Every stretch of 10 bases of every allele is noted in
   an index. A read is cut into such stretches too, and the alleles that
   share stretches with it are its candidates. The shared stretches also
   show where in the read the allele lies. Alleles shorter than 30 bases
   hold too few such stretches and are compared with every read.
2. **Alignment.** Each candidate is aligned to the read around the
   position found in step 1. The alignment allows mismatches and gaps
   (the method of Gotoh, with separate costs for starting and extending a
   gap). Only positions near the expected place are tried; if the
   alignment runs to the border of that region, the region is widened and
   the alignment repeated.
3. **Choice.** The allele with the highest `identity` wins. Between equal
   ones the allele with the larger matched part wins, then the higher
   score, then the allele that comes first in the database file. The
   result is therefore the same in every run and on every computer.

Not every candidate has to be aligned. One difference between a read and
an allele can remove only a limited number of shared stretches, so an
allele that shares too few of them cannot match the read as well as the
best alignment found so far, and is skipped. The alleles that could match
equally well are all aligned, up to 256 per read.

Equal reads are aligned once. The work is spread over all processor
cores; set the environment variable `RAYON_NUM_THREADS` to use fewer.

## Limits

- Only A, C, G and T count as equal. An `N` in a read or an allele is a
  mismatch.
- One row per read: the best allele. Reads that span two genes are
  reported for the better matching one.
- A read is expected to hold the gene in one piece. Alignments in which
  read and allele drift apart by more than about 100 bases of gaps are not
  found.
- An allele of 30 bases or more is found through the stretches of 10
  bases it shares with the read. A read that differs from it so often
  that no such stretch is left is not matched to it.
