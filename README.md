# immunediscover

[![Tests](https://github.com/mashu/immunediscover/actions/workflows/ci.yml/badge.svg)](https://github.com/mashu/immunediscover/actions/workflows/ci.yml)
[![Documentation](https://github.com/mashu/immunediscover/actions/workflows/docs.yml/badge.svg)](https://mashu.github.io/immunediscover/)
[![Release](https://img.shields.io/github/v/release/mashu/immunediscover)](https://github.com/mashu/immunediscover/releases/latest)

Finds the germline alleles of immune receptor genes (V, D and J) in
amplicon libraries made from genomic DNA, one plate of donors at a time:
which known alleles every donor carries, and which sequences are not in
the database yet.

Install on Linux or macOS:

```
curl -fsSL https://github.com/mashu/immunediscover/releases/latest/download/install.sh | sh
```

```
immunediscover preprocess demultiplex plate.fastq.gz indices.tsv demux.tsv
immunediscover search exact   demux.tsv.gz IGHV.fasta exact.tsv -g V
immunediscover discover blast demux.tsv.gz IGHV.fasta novel.tsv -g V
```

`discover blast` needs BLAST+ (`blastn`, `makeblastdb`) on the `PATH`, or
`--aligner builtin`. Nothing else is needed.
