# immunediscover

[![Tests](https://github.com/mashu/immunediscover/actions/workflows/ci.yml/badge.svg)](https://github.com/mashu/immunediscover/actions/workflows/ci.yml)
[![Documentation](https://github.com/mashu/immunediscover/actions/workflows/docs.yml/badge.svg)](https://mashu.github.io/immunediscover/)

Finds the germline alleles of immune receptor genes (V, D and J) in
amplicon libraries made from genomic DNA, one plate of donors at a time:
which known alleles every donor carries, and which sequences are not in
the database yet.

```
cargo build --release
./target/release/immunediscover --help
```

```
immunediscover preprocess demultiplex plate.fastq.gz indices.tsv demux.tsv
immunediscover search exact   demux.tsv.gz IGHV.fasta exact.tsv -g V
immunediscover discover blast demux.tsv.gz IGHV.fasta novel.tsv -g V
```

`discover blast` needs BLAST+ (`blastn`, `makeblastdb`) on the `PATH`;
`search bwa` needs `bwa`.
