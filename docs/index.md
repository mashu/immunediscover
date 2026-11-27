# immunediscover: user guide

immunediscover finds the germline alleles of immune receptor genes (V, D
and J) in amplicon libraries made from genomic DNA, one plate of donors at
a time. It tells you which known alleles every donor carries, and which
sequences in your donors are not in the database yet.

```
FASTQ of a plate ──demultiplex──▶ reads per donor ──┬─ search exact ──▶ known alleles per donor ──▶ genotypes
                                                    └─ discover blast ─▶ new alleles, with evidence
```

!!! tip "New here?"
    Start with [Getting started](getting-started.md): it runs the whole
    analysis once, from the FASTQ of a plate to a list of new alleles.

| Page | Read it when |
|---|---|
| [Installing](install.md) | you do not have the program yet |
| [Getting started](getting-started.md) | you want to run the whole analysis once, from FASTQ to new alleles |
| [Input files](inputs.md) | you are preparing the files, or a command complains about a column |
| [Commands](commands.md) | you need an option or the meaning of an output column |
| [`discover blast` in depth](discover-blast.md) | you are looking for new alleles and need to choose thresholds |
| [Evidence](evidence.md) | you have a list of candidates and must decide which to believe |
| [For developers](developers.md) | you want to change the code or check it against the Julia original |

## The commands at a glance

| Command | What it does |
|---|---|
| `preprocess demultiplex` | Sorts the reads of a plate into donors. |
| `search exact` | Reports the known alleles each donor carries and how many reads support them. |
| `discover blast` | Finds sequences that resemble known alleles without being in the database. |
| `discover hsmm` | Finds D genes by their recombination signals alone. |
| `align` | Reports the most similar allele of every read, with its differences. |
| `search heptamer` | Finds where each V allele ends from the heptamer that follows it. |
| `search bwa` | Checks that candidates map to the expected chromosome of the genome. |
| `analyze haplotype` | Calls each donor homozygous or heterozygous for each gene. |
| `analyze cooccurrence` | Finds alleles that tend to occur in the same donors. |
| `table ...` | Joins, filters, sorts and exports tables. |
| `fasta ...` | Merges, compares and names FASTA files. |

## Words used in this guide

| Word | Meaning |
|---|---|
| donor, case | One individual. `case` is the column name. |
| well | One position of the plate; normally one donor. |
| allele | One sequence variant of a gene, named `GENE*NUMBER`. |
| database | Your FASTA file of known alleles. |
| new (novel) allele | A sequence close to a database allele but not identical to any; named `<nearest allele>_S<digits>`. |
| count | Number of reads supporting a sequence in a donor. |
| ratio | A count relative to the best allele of the same gene in the same donor. |
| flank, surroundings, context | The bases immediately before and after the allele in the read. |
| RSS | Recombination signal sequence next to a gene: heptamer, spacer, nonamer. |
| copies | Estimated number of copies of a sequence in a donor's genome (1 = one chromosome). |
| evidence | How unlikely it is that a sequence is an artifact; see [evidence](evidence.md). |
