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
2. **Assigns every read to its most similar allele**, with BLAST or with
   the aligner built into the program (`--aligner`, see
   [Search](#search)). Reads that resemble a pseudogene more (if you give
   `--pseudo`) are dropped.
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

D genes are short, so they need longer surroundings. These values apply
with BLAST and with the built-in aligner alike.

With BLAST, each preset also chooses the settings of `blastn` itself
(`--args`); for D genes that is a more sensitive search for short
matches. `-G` lists them. The built-in aligner has no such settings: it
treats short and long genes the same way, and `--args` cannot be combined
with `--aligner builtin`.

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
| `--aligner` | `blast` | Which program compares the reads with the alleles: `blast` (NCBI BLAST+, which must be installed) or `builtin` (part of immunediscover, nothing else needed). |
| `-a`, `--args` | megablast settings | Extra settings passed to `blastn`. Change only if you know BLAST. Not allowed with `--aligner builtin`. |
| `--min-identity` | 0.7 | Only with `--aligner builtin`: how near a read must be to an allele (with its surroundings) to be given to it. 0.7 allows three differences per ten bases, so no read that resembles an allele is left out. Reads further from every allele get no allele. Raising it (for example to 0.9, one difference per ten bases) saves time on plates with many reads that belong to no gene of the database. |
| `-i`, `--isin` | off | A found sequence that is a piece of a known allele normally takes that allele's name. With `-i` it is reported as new. |
| `-g`, `--gene` | none | Preset for `V`, `D` or `J`. |
| `-G`, `--show-presets` | | Print the presets and stop. |

#### BLAST or the built-in aligner

```
immunediscover discover blast demux.tsv.gz IGHV.fasta novel.tsv -g V --aligner builtin
```

Both programs answer the same question, which allele a read resembles
most, and everything after that step is the same. They differ in how they
treat the ends of an allele. BLAST reports the stretch where read and
allele agree best and may stop a few bases before the end when the last
bases differ. The built-in aligner always compares the allele, with its
recorded surroundings, from the first to the last base, so a difference
at the very end of an allele is counted and the sequence is reported at
full length. This matters most when few or no surrounding bases are
recorded (`--forward 0 --reverse 0`) and for short genes.

The built-in aligner needs no installation and gives the same result on
every computer. It shows its progress, the reads handled per second and
the time left while it runs, and the total time at the end. Its results
are not cached: every run aligns again. [The built-in aligner](align.md)
describes how it works.

It gives every read to the nearest allele (fewest differences per base).
A read that holds an allele with its surroundings unchanged is recognised
at once; every other read is compared with the alleles it resembles. No
read is set aside at this step because it is far from the database: a
read is left without an allele only if it differs from every allele in
more than three of ten bases. What is reported is decided afterwards, by
`--maxdist` and the other thresholds.

Reads far from every allele take the most time, because each is compared
with every allele it shares anything with. If a plate holds many of them
and the run is slow, `--min-identity` sets how far a read is followed:

```
immunediscover discover blast demux.tsv.gz IGHV.fasta novel.tsv -g V --aligner builtin --min-identity 0.9
```

With 0.9 a read more than one difference in ten bases away from every
allele (counted over the allele and its recorded surroundings) gets no
allele. Keep it well below what `--maxdist` asks for: with the V preset
(14 differences in an allele of about 300 bases, plus differences in the
40 bases around it) 0.9 leaves room; 0.97 would not.

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
settings are unchanged; `-o` forces a new search. Only BLAST results are
cached; with `--aligner builtin` every run aligns again.

**"blastn command not found".** Install BLAST+ (see
[Installing](install.md#what-else-is-needed)), or run with
`--aligner builtin`, which does not need it.

**"Another discover blast run is working on the same reads and
database".** Every combination of reads, database and pseudogenes has its
own folder of intermediate files, and one run uses it at a time. Runs on
other reads or with another database (V, D and J on the same plate, for
example) can go on side by side. For two runs on the same reads and
database, wait for the first or give the second another `--work-dir`. If
the message appears although no run is active, delete the `lock` file it
names.
