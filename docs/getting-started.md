# Getting started

## What you need

- The program: `immunediscover` (see [Installing](install.md)).
- For `discover blast`: NCBI BLAST+ (`blastn` and `makeblastdb` on your
  `PATH`), unless you run it with `--aligner builtin`.
- For `search bwa`: `bwa`.
- Your data: a FASTQ file per plate, an index sheet, and a FASTA database
  of known alleles for each gene type. See [inputs](inputs.md).

## Getting help

Every command and sub-command explains itself with `--help`:

```
immunediscover --help
immunediscover discover blast --help
```

## A complete run for V genes

```
# 1. Split the plate into donors
immunediscover preprocess demultiplex plate.fastq.gz indices.tsv demux.tsv

# 2. Which known alleles does each donor carry?
immunediscover search exact demux.tsv.gz IGHV.fasta exact.tsv -g V

# 3. Which alleles are not in the database yet?
immunediscover discover blast demux.tsv.gz IGHV.fasta novel.tsv -g V

# 4. Export the new alleles found in at least two donors
immunediscover table fasta novel.tsv.gz novel.fasta \
    --colname allele_name --colseq aln_qseq --filter "_S" --mincase 2

# 5. Add them to the database and repeat step 2 with the larger database
immunediscover fasta merge IGHV_plus.fasta IGHV.fasta novel.fasta
immunediscover search exact demux.tsv.gz IGHV_plus.fasta exact_plus.tsv -g V

# 6. Genotype per donor
immunediscover analyze haplotype exact_plus.tsv.gz genotypes.tsv
```

For D and J genes, repeat steps 2 to 6 with the D or J database and
`-g D` or `-g J`.

### What to look at after each step

1. **Demultiplex.** `demux.tsv.log` lists the number of reads per well.
   Wells with very few reads will give weak results later.
2. **Exact search.** The table has one row per known allele per donor.
   `count` and `ratio` say how strongly the allele is supported; a gene
   plot is printed to the screen (switch it off with `-n`).
3. **Discovery.** Sort by `allele_name`. Names with `_S` are sequences
   not in the database. For each, read `db_distance`, `copies`,
   `evidence`, `donors`, `carriers` and `segregates`
   ([what they mean](evidence.md)).
4. **Export.** Open the FASTA and check the number of records against
   your expectations before adding them to the database.

## Where things are written

- Results go to the file you name; `.gz` is added where the table is
  compressed.
- `discover blast` keeps its intermediate files in `.immunediscover/`
  (change with `--work-dir`). Deleting the folder is always safe; the next
  run is just slower.
- `immunediscover.log` in the current folder records every command with
  its settings. Keep it with your results.

## Next

- [Input files](inputs.md)
- [All commands, options and output columns](commands.md)
- [`discover blast` in depth](discover-blast.md)
- [The evidence columns](evidence.md)
