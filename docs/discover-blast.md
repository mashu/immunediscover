# `discover blast`: finding new alleles

`discover blast` looks in every donor for sequences that resemble the
alleles of a database without being identical to them, and reports them
together with the known alleles it finds. It is the command to use when
you want to know which alleles a cohort carries that your database does
not have yet.

```
immunediscover discover blast demux.tsv.gz IGHV.fasta novel.tsv -g V
```

- `demux.tsv.gz`: the reads, one row per read, from
  `preprocess demultiplex` (see [inputs](inputs.md)).
- `IGHV.fasta`: the database of known alleles.
- `novel.tsv`: where to write the result (saved as `novel.tsv.gz`).
- `-g V`: use the settings tuned for V genes (`-g D`, `-g J` likewise).

## What it does, step by step

1. **Learns the surroundings of every allele.** For each database allele
   the reads that contain it are collected, and the bases just before and
   after the allele (20 on each side by default) are recorded. Reads are
   then compared with "allele plus surroundings". This anchors both ends
   of the allele, so that a read is not mistaken for a shorter or longer
   allele.
2. **Assigns every read to its most similar allele** with BLAST. Reads
   that resemble a pseudogene more (if you give `--pseudo`) are dropped.
3. **Groups identical reads** within each donor. The number of reads in a
   group is its `full_count`.
4. **Cuts the surroundings off again**, leaving the allele-sized sequence
   (`aln_qseq`), and counts its differences from the database allele
   (`aln_mismatch`).
5. **Counts per allele.** The same allele sequence can sit in slightly
   different surroundings; `count` adds them all up.
6. **Applies the thresholds** (below).
7. **Names the sequences.** An exact match keeps the database name. Any
   other sequence is named after its nearest allele plus `_S` and four
   digits derived from the sequence, for example `IGHV1-69*01_S1234`. The
   same sequence always gets the same name, in every donor and every run.
8. **Adds the evidence columns** (see [evidence](evidence.md)).

Results are per donor: one row per sequence found in a donor. The rest of
the plate is used only to learn the surroundings (step 1) and the rates
behind the evidence (step 8).

## The thresholds that decide what is reported

| Option | Default | What it is for |
|---|---|---|
| `-c`, `--minfullcount` | 5 | Minimum reads. Applied twice: to the reads of the allele in its best surroundings (`full_count`) and to all reads of the allele (`count`). Raise it to remove rare errors; lower it for shallow donors. |
| `-f`, `--minfullratio` | 0.1 | Minimum reads relative to the best allele of the same gene in the donor (`full_ratio` and `ratio`). 0.1 keeps sequences with at least a tenth of the reads of the gene's top allele. Errors usually stay far below; an allele on one chromosome next to an allele on the other is near 1. |
| `-d`, `--maxdist` | 20 | Maximum differences from the database allele. Lower it to keep only close relatives of known alleles. |
| `-l`, `--length` | 290 | Minimum length of the part of the read that matches (allele plus surroundings). Removes fragments. Must be much smaller for D and J genes (see presets). |
| `--min-corecov` | 0.6 | Minimum length of the found sequence relative to the database allele. Removes truncated sequences. |
| `-t`, `--top` | 1 | How many versions of the surroundings to report per allele. 1 reports the best supported one. |

!!! note
    The evidence columns do **not** remove anything. They are there so
    that you can set the count thresholds generously and judge the calls
    afterwards.

## Presets for V, D and J

`-g V`, `-g D` and `-g J` replace the defaults with values that suit the
gene type. Any option you give yourself takes precedence. Show them with
`immunediscover discover blast -G`:

| Option | default | `-g V` | `-g D` | `-g J` |
|---|---|---|---|---|
| `--forward` | 20 | 20 | 40 | 12 |
| `--reverse` | 20 | 20 | 40 | 12 |
| `--minfullratio` | 0.1 | 0.035 | 0.2 | 0.1 |
| `--length` | 290 | 283 | 5 | 10 |
| `--maxdist` | 20 | 14 | 20 | 10 |
| `--minfullcount` | 5 | 5 | 5 | 10 |
| `--edge` | 0 | 0 | 10 | 0 |
| `--subjectcov` | 0.1 | 0.1 | 0.25 | 0.1 |
| `--minquality` | 0.75 | 0.62 | 0.5 | 0.75 |
| `--min-corecov` | 0.6 | 0.5 | 0.6 | 0.6 |

D genes are short, so they need longer surroundings and a more sensitive
BLAST search; the D preset also changes the BLAST settings (`--args`).

## All options

### Files

| Option | What it is for |
|---|---|
| `-p`, `--pseudo <FASTA>` | Pseudogenes. Reads that match one of them better than any allele are discarded, so that pseudogenes are not reported as new alleles. |
| `--work-dir <DIR>` | Where intermediate files are kept (default `.immunediscover`). The BLAST results are cached there, so a second run with other thresholds is fast. |
| `-o`, `--overwrite` | Ignore the cache and search again. Use it after changing the reads or the database in place. |
| `-v`, `--verbose` | Also write `<output>.unfiltered.tsv.gz`: every sequence before the thresholds. Useful to see what a threshold removed. |
| `--no-evidence` | Leave out the evidence columns. |

### Surroundings (context) of the alleles

| Option | Default | What it is for |
|---|---|---|
| `--forward` | 20 | Bases recorded before the allele. |
| `--reverse` | 20 | Bases recorded after the allele. |
| `--extension-maxdist` | 0 | How different a read may be from an allele and still be used to learn its surroundings. 0 uses only reads containing the allele exactly. Raise it (for example to 3) when the database comes from another species or lacks most alleles of your donors, so that no read contains a database allele exactly. |
| `--no-gene-context` | off | By default an allele that occurs in no read borrows the surroundings of the other alleles of its gene, so that all alleles of a gene are searched alike. This switches the borrowing off. |
| `-q`, `--minquality` | 0.75 | How well the recorded surroundings must match the read when they are cut off again (step 4). A sequence whose surroundings match worse is dropped. Lower it if many reads are lost in "failed trimming". |
| `--keep-failed` | off | Keep sequences whose surroundings could not be cut off (their `aln_qseq` is empty and `aln_mismatch` is -1). For troubleshooting. |

### Search

| Option | Default | What it is for |
|---|---|---|
| `-s`, `--subjectcov` | 0.1 | Minimum share of the database allele that the read must cover. |
| `-e`, `--edge` | 0 | Minimum bases between the allele and the end of the read. Removes alleles cut off by the read end. |
| `-a`, `--args` | megablast settings | Extra settings passed to `blastn`. Change only if you know BLAST. |
| `-i`, `--isin` | off | A found sequence that is a piece of a known allele normally takes that allele's name. With `-i` it is reported as new. |
| `-g`, `--gene` | none | Preset for `V`, `D` or `J`. |
| `-G`, `--show-presets` | | Print the presets and stop. |

## The output table

One row per sequence found in a donor.

| Column | Meaning |
|---|---|
| `well` | Well of the plate (row number in the index sheet). |
| `case` | Donor identifier. |
| `allele_name` | Name of the sequence: a database name, or `<nearest allele>_S<digits>` for a sequence not in the database. **This is the column to look at first.** |
| `gene` | Gene of the nearest database allele. |
| `sseqid` | The nearest database allele. |
| `aln_qseq` | The sequence found, allele-sized. This is what to export as a new allele. |
| `aln_mismatch` | Differences between `aln_qseq` and the nearest database allele. 0 = known allele. |
| `count` | Reads of this sequence in this donor, in any surroundings. |
| `ratio` | `count` relative to the best allele of the gene in this donor (1 = it is the best). |
| `full_count` | Reads of this sequence in the surroundings shown in `qseq`. |
| `full_ratio` | `full_count` relative to the best of the gene in this donor. |
| `flank_index` | Rank of this version of the surroundings (1 = best supported). Only 1 appears unless `--top` is raised. |
| `qseq` | The sequence with its surroundings, as found in the reads. |
| `mismatch` | Differences between `qseq` and the database allele with its surroundings. |
| `db_length` | Length of the database allele. |
| `corecov` | Length of `aln_qseq` relative to `db_length`. |
| `db_distance`, `copies`, `evidence`, `relation`, `donors`, `carriers`, `segregates`, `usual_parent_ratio` | The evidence columns, explained in [evidence](evidence.md). |

### Other files written

| File | Content |
|---|---|
| `<output>.flank-split.tsv.gz` | Written only when needed: alleles with enough reads in total (`count`) whose reads are spread over several surroundings, none of which passes on its own. They are missing from the main table; this file tells you so. Lowering `-c`/`-f` or checking the surroundings recovers them. |
| `<output>.unfiltered.tsv.gz` | With `-v`: all sequences before the thresholds. |
| `immunediscover.log` | Every command you ran, with date and settings, appended in the current folder. |

## From the table to a list of new alleles

```
# keep the new ones (names with _S), one FASTA record per sequence
immunediscover table fasta novel.tsv.gz novel.fasta \
    --colname allele_name --colseq aln_qseq --filter "_S" --mincase 2
```

`--mincase 2` keeps sequences found in at least two donors. Before
exporting, sort the table by `allele_name` and look at `copies`,
`evidence`, `db_distance`, `donors` and `segregates` for each candidate.

A check against the genome is available with
[`search bwa`](commands.md#search-bwa): it keeps the candidates that map
to the expected chromosome.

## Troubleshooting

**A known allele of my donor is missing.** Look in
`<output>.flank-split.tsv.gz`, then run with `-v` and search for it in
the unfiltered table: its `count`, `ratio` and `aln_mismatch` show which
threshold removed it.

**Too many new alleles.** Raise `-c` and `-f`, give `--pseudo`, and use
the evidence columns: sequences with evidence 0, with `segregates` false,
or many differences from the database are the first to review.

**Nothing is found with a database from another species.** No read
contains a database allele exactly, so no surroundings are learned. Set
`--extension-maxdist` to a few differences.

**The second run is not faster / uses old results.** The cache is in
`--work-dir`. It is reused only when reads, database and surroundings
settings are unchanged; `-o` forces a new search.
