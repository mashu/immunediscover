# Evidence columns

`discover blast` reports every sequence that passes the read-count
thresholds. A count alone is hard to judge: 40 reads are a lot for a gene
that amplifies poorly in a donor with little DNA, and very little for a
gene that amplifies well in a deep donor. The evidence columns put every
reported sequence on a common scale and say what the rest of the plate
knows about it.

## What they measure

The `evidence` column measures one thing: whether a sequence is present in
a donor's DNA in a real amount, or is an artifact (a PCR or sequencing
error, a chimera) of something else in the well.

It does not establish which gene a real sequence belongs to. A sequence
from another locus or a pseudogene is real DNA and gets strong evidence.
The other columns are there for that decision: the distance to the
database, the number of donors, the number of copies and whether the
sequence segregates.

## The columns

| Column | What it tells you |
|---|---|
| `db_distance` | Differences (substitutions, insertions, deletions) from the nearest database allele. 0 is a known allele. Most new alleles are 1 to 3 differences away; a sequence 5 or 10 away is more often another gene. |
| `copies` | How many copies of the sequence the donor holds, estimated from the reads. About 1 for an allele on one chromosome, about 2 for an allele on both. Far below 1 for artifacts. Above 2 when the sequence sits in more than two places in the genome. |
| `evidence` | How unlikely it is that the sequence is an artifact. Read it as "1 in 10 to the power of evidence": **1** means about one in ten sequences like this one is an artifact, **2** one in a hundred, **3** one in a thousand. **0** means artifacts are as common as real sequences at this amount. Empty when the donor has too few reads of the gene to judge. |
| `relation` | What the sequence could be an artifact of, in this donor: `1 difference`, `2 differences`, `3 differences` or `4+ differences` from a more abundant sequence of the same gene, a `chimera` of two of them, or `top of gene` when nothing in the gene is more abundant. |
| `donors` | Donors on the plate in which the sequence was seen at all, even as a single read. |
| `carriers` | Donors on the plate holding at least a quarter of a copy of it. |
| `segregates` | `true` when some donors carry the sequence and others do not, as alleles do. `false` when nearly every donor (90 % or more) holds it at the same level. Such a sequence is not a chance artifact, and is not shown to be an allele of this gene either: think of another locus, a duplicated gene or a systematic artifact. |
| `usual_parent_ratio` | For a sequence that most donors show: its usual amount relative to the more abundant sequence it resembles (0.01 = 1 %). Empty when most donors do not show it. |

### Reading a row

| allele_name | count | db_distance | copies | evidence | relation | donors | carriers | segregates |
|---|---|---|---|---|---|---|---|---|
| IGHV1-69*01_S1234 | 163 | 1 | 1.05 | 2.4 | 1 difference | 6 | 5 | true |
| IGHV3-15*01_S5678 | 648 | 11 | 9.8 | 6.4 | top of gene | 95 | 95 | false |
| IGHV2-26*02_S9012 | 19 | 1 | 0.33 | 0 | 1 difference | 2 | 1 | true |

Example values for three typical cases.

- **First row.** One copy, one difference from the database, evidence 2.4
  (fewer than one in a hundred such sequences is an artifact), carried by
  5 donors and absent from the others. This is what a new allele looks
  like.
- **Second row.** Very strong evidence that the molecule is real, but it
  is 11 differences from the database, every donor has it, and it comes
  in about ten copies. It is real DNA and almost certainly not an allele
  of IGHV3-15.
- **Third row.** A third of a copy and evidence 0: at this amount the
  plate has as many artifacts as real sequences. The call may still be
  right (a gene present in several copies gives less than a full "copy"
  per allele), but the reads cannot support it. Look for it in other
  plates.

### Strength of evidence

| evidence | meaning |
|---|---|
| 2 or more | strong: artifact unlikely |
| 1 to 2 | moderate |
| above 0, below 1 | weak |
| 0 | none: indistinguishable from the artifacts of the plate |
| empty | not judged: fewer than 20 reads expected per copy of this gene in this donor |

!!! warning "Evidence is about artifacts, not about genes"
    Evidence says nothing about which gene a sequence belongs to. Always
    read it together with `db_distance`, `copies` and `segregates`.

## How it is computed

The evidence comes from three models, each fitted to the plate being
analysed. Nothing is trained in advance and no reference data are used.
All sequences of the run enter the fit, including those the count
thresholds remove: they show what the artifacts of the plate look like.

### 1. From reads to copies

The library is made from genomic DNA, so every gene is present in the
same number of copies in every cell, whatever the cell does. The number
of reads of a gene in a donor then depends on only three things:

```
reads = amount of DNA of the donor  x  efficiency of the primers for the gene  x  copies
```

The first factor is one number per donor, the second one number per gene.
Both are estimated from the whole plate (by median polish of the log read
counts, a robust fit that a deleted or duplicated gene in a few donors
does not disturb). Dividing a sequence's reads by "DNA amount x primer
efficiency" leaves its number of copies.

Details that matter:

- The reads of a gene are taken from its **database alleles**. An
  abundant sequence of another locus that the search assigned to the gene
  therefore cannot distort the gene's unit.
- A gene has **two copies** unless the donors show otherwise. Where a
  donor holds two or more database alleles of a gene, the smallest is
  taken as one copy; if it is typically a third or a quarter of the gene,
  the gene has three or four copies in most donors (IGHV1-69, IGHV3-30
  and the like), and one copy is a third or a quarter of its reads. This
  estimate cannot tell an extra copy from an allele that amplifies
  poorly; it errs towards giving sequences more copies, never fewer.
- A gene seen in fewer than three donors takes its unit from the donor
  itself (half of what the donor shows for the gene).
- Where one copy would give fewer than 20 reads, `evidence` is left empty:
  a few stray reads would already look like a copy.

### 2. How often do artifacts reach this amount?

Artifacts arise by chance, in one donor at a time, and nearly all of them
stay far below one copy. The question for a sequence at, say, 0.8 copies
is how many artifacts of its kind the plate has at about 0.8 copies.

1. Sequences are grouped by `relation`, because each kind of artifact has
   its own rate: single-base errors are common, sequences four
   differences from everything else are rare.
2. Below a fifth of a copy the sequences of a plate are artifacts, a few
   exceptions aside. Those private to one donor are counted there in two
   ranges (0.02 to 0.063 and 0.063 to 0.2 copies).
3. Artifacts thin out quickly with amount. The drop from the lower range
   to the upper one is continued upwards (a power-law tail, the usual way
   to extend the far end of a distribution). This gives the number of
   artifacts expected around any amount.
4. That expected number, divided by the number of sequences actually seen
   around the same amount, is the expected share of artifacts among them:
   a local false discovery rate. `evidence` is minus its log10.

```
evidence = -log10( artifacts expected around this amount / sequences seen around this amount )
```

Safeguards:

- A sequence private to one donor is compared with other private
  sequences only. One that other donors also carry is compared with all,
  which is the generous side: chance artifacts do not repeat across
  donors, so carried-by-several is itself evidence.
- More copies never give weaker evidence than fewer.
- The drop per range is kept between 2-fold and 20-fold, half an artifact
  is added to every count, and a kind with fewer than 20 artifacts uses
  the drop of all kinds together. With no artifacts to learn from, the
  evidence stays modest instead of becoming certain.
- Sequences that do not segregate are left out of these counts.

### 3. Artifacts that recur

Some errors are not chance events: the same variant of the same parent
sequence appears in most donors, at a low and steady ratio to the parent
(typically a fraction of a percent). Its amount in one
donor is then judged against **its own** usual ratio in the other donors.
A donor that really carries the sequence holds a hundred times the usual
ratio; a fluctuation holds two or three times.

- The usual ratio is the median ratio to the parent in the other donors.
- How much the ratio varies between donors is measured on the plate, and
  sampling noise is added for small counts.
- Which sequences count as recurrent artifacts is also measured: up to
  three spreads above the middle of the usual ratios on the plate, and
  never above 10 % of the parent. A sequence that most donors show at a
  higher ratio keeps the evidence of step 2, and its ratio is reported in
  `usual_parent_ratio`.

This step can only lower the evidence.

## What it can and cannot do

- **It needs a plate.** The rates are learned from the donors analysed
  together. With a handful of donors there is little to learn from, units
  fall back to the donor itself, and the evidence stays modest.
- **It needs genomic DNA.** The "same copies in every cell" premise does
  not hold for expressed repertoires.
- **It judges amounts.** A real allele that amplifies at a fifth of the
  usual efficiency looks like an artifact. A real molecule from another
  locus gets strong evidence. Neither is a mistake of the method; both
  are reasons to read the other columns.
- **Below a fifth of a copy it has nothing to say** (evidence 0).

## Switching it off

`--no-evidence` leaves the eight columns out. The other columns of the
output are identical with and without them.

## In the source code

`src/evidence/`: `amplification.rs` (copies), `parents.rs` and `kinds.rs`
(relation), `plate.rs` (donors, carriers, segregation, usual ratio),
`artifact_rate.rs` (step 2), `recurrence.rs` (step 3). The evidence can
be checked against alleles held out of a database with
`ki-check --evidence`, described under [For developers](developers.md).
