# Specificity annotation: how it works

`immvoke annotate` and `annotateSpecificity()` match a query AIRR table
against a reference AIRR table and report, for every query row, which
reference rows it matches. This page explains the flow and the rules behind
each step; [README.md](README.md) covers how to run it.

## Workflow

Each query row goes through these steps in order. Steps marked *automatic*
always run; the others are set with the option shown (the
`annotateSpecificity()` argument has the same name, with `_` for `-`).

| # | Step | Option | Choices | Default |
| --- | --- | --- | --- | --- |
| 1 | Check the reference is AIRR | *automatic* | -- | -- |
| 2 | Paired or single chain? | `--paired` | on, off | off |
| 2a | Paired: one heavy and one light chain per query cell | `--count-field` | an AIRR count column | `duplicate_count` |
| 3 | Keep heavy chain rows | *automatic* | `IGH` only | -- |
| 4 | Sequence to compare | `--sequence-region` | `CDR3H`, `VDJ`, `CDRs` | `CDR3H` |
| 5 | Gene gate | `--partition-by` | none, `V`, `VJ`, `VJL` | none |
| 6 | How close is a hit | `--identity-threshold` | a fraction in (0, 1] | exact |
| 7 | Distance | `--distance-metric` | `hamming`, `levenshtein` | `hamming` |
| 8 | Search | *automatic* | symscan; rows too large for it are scored directly instead, and counted in `annotation_metadata.yml` | -- |
| 9 | Output | `--outdir` (`outdir` in Python), `--full-dist-matrix` | `annotation.tsv` (each query's best hits, also returned as a DataFrame) and `annotation_metadata.yml`, plus `annotation_full_matrix.txt` (every hit within the threshold) with `--full-dist-matrix`, each prefixed with the reference's database and table (`iedb_bcr_annotation.tsv`) | -- |

With `--paired` on, steps 3-7 are replaced by one fixed comparison: an
antibody (`cell_id`) matches a reference antibody only if both its heavy and
its light chain have the same V gene, the same J gene and an exact CDR3 --
classified `both`, `heavy_only`, `light_only`, `conflicting` or `neither`.
Without it, every query heavy chain row that matched something gets one
result row -- an unmatched row is left out, counted in the metadata
instead: its best hits -- the highest identity, all of them on a tie --
with `n_hits`, `identity`, `hit_ids` and each hit's `antigen`, `epitope`
and `antigen_source_species` from the reference, `;` between hits.
`annotation_full_matrix.txt` (`--full-dist-matrix`) lists every hit within
the threshold, one (query, hit, identity) per row.

The output files are named after the reference, `<db>_<table>_` -- the
database read off the `download_metadata.yml` `immvoke specificity` writes
beside the reference table, the table off its file name -- so runs against
different databases stay apart, even in one `--outdir`. A reference without
that record is named by its file name alone.

`annotation_metadata.yml` records the run: the options, the query and
reference (with the reference's source database, read from the
`download_metadata.yml` beside it), the counts -- rows scored directly
included -- and the immvoke version, time and command.

## Input contract

The reference must always be an **AIRR rearrangement table**, whichever
database it comes from: every field the AIRR schema requires (`sequence_id`,
`sequence`, `rev_comp`, `productive`, `v_call`, `d_call`, `j_call`,
`sequence_alignment`, `germline_alignment`, `junction`, `junction_aa`,
`v_cigar`, `d_cigar`, `j_cigar`) present as a column, though values may be
blank. Any table immvoke writes as AIRR already is one; anything else is
refused with an error naming the missing fields.

On top of that, each step reads its own AIRR fields: `cdr3_aa`,
`sequence_aa` or `cdr1_aa`/`cdr2_aa` for the sequence region, `v_call` and
`j_call` for partitions and paired matching, `locus` for the chain filter,
and `cell_id` to link the two chains of a paired antibody. The query only
needs the fields the chosen options read. A table without a `locus` column
is compared whole; a table with one but no `IGH` row left is an error, never
an empty result.

`cdr3_aa` must follow the AIRR definition in both tables: the CDR3 without
its two flanking conserved residues, i.e. `junction_aa` minus its first and
last residue (`ARGLDS`, not `CARGLDSW`). IEDB's does; annotate does not check
or convert a query's, so a query whose `cdr3_aa` holds junctions will find
no exact match.

A blank `locus` is kept as the reference has it -- nothing is inferred
from the gene calls -- so such a row is in neither the heavy nor the light
chain pool: single-chain matching never compares it, and paired matching
cannot pair it. (Of IEDB's human chains, 277 have a blank locus, nearly all
of them light chains.)

## Matching rules

**Sequence region.** `CDR3H` compares the heavy chain CDR3 (`cdr3_aa`), `VDJ` compares
`sequence_aa`, and `CDRs` pools `cdr1_aa`, `cdr2_aa` and `cdr3_aa` under one
shared tolerance. Pooling only works with Hamming distance, because edit
distance does not add up across columns.

**Blanks.** A row with a blank value in any compared column (say, no
`cdr1_aa` under `CDRs`) has nothing to compare and never matches -- not
even another blank. A blank query row gets no hit; a blank reference row is
never hit.

**Partition.** `V` and `VJ` only compare row pairs whose V-gene (and
J-gene) agree; `VJL` also requires the same CDR3 length, and so needs
`cdr3_aa` in both tables whatever the sequence region is.

**Gene names.** `v_call`/`j_call` cells are not always one clean name: IEDB
mixes allele calls (`IGHV1-2*02`), gene calls (`IGHV1-2`), several calls in
one cell (`IGHV1-2*02, IGHV1-2*04`) and blanks. They are parsed with
Change-O's gene pattern (`Genes.py`) and compared by exact gene name,
whatever the allele: `IGHV1-2*02`, `IGHV1-2*04` and `IGHV1-2` all agree,
while `IGHV1-2` and `IGHV1-69` do not.

A cell with several calls agrees if any one of them does. A blank or
unrecognized call **agrees with nothing**, not even another blank; the
number of unrecognized values is logged. So under `VJ`, a row with no J
call can never match -- use `V`, or no partition, to include such rows.

**Threshold and metric.** An unset (or 1.0) `identity_threshold` is an
exact match. Below that, a query may differ from a hit by up to
floor((1 - threshold) x query length) positions: substitutions only under
`hamming`, or substitutions, insertions and deletions under `levenshtein`,
which also lets the two differ in length. Hamming reports which positions
differ; Levenshtein cannot, since an indel shifts every position after it.

**X residues.** An `X` is an unresolved residue. Under exact and Hamming
matching it matches any residue, so `CARXYW` equals `CARDYW`; under
Levenshtein it is an ordinary letter.

**Paired.** An antibody is its heavy (`IGH`) and light (`IGK`/`IGL`) chain
rows sharing one `cell_id`. It matches a reference antibody only if that
*same* antibody agrees on V-gene and J-gene (same gene, any allele) and
exact CDR3 on both chains; a heavy hit on one reference antibody and a light
hit on another is `conflicting`, not `both`. Region, partition, threshold
and metric do not apply. A row with a blank `cell_id` has no partner chain,
so it is left out, on either side, with a warning giving the count.

A query cell can carry more than one heavy chain, or more than one light
chain (a doublet, ambient RNA, or a second light chain). Before matching,
each query cell keeps one heavy (`IGH`) and one light (`IGK` or `IGL`)
chain: the one with the highest `--count-field` value (`duplicate_count` by
default, or e.g. `umi_count`). On a tie, a blank count, or a file without
that column, the chain that comes first in the file is kept. The reference
is not reduced this way: a reference antibody with several chains of one
kind (4 heavy and 26 light cases among IEDB's human antibodies, which carry
no count column) keeps them all.

## Engine

- **symscan finds the hits.** It is a required dependency, a
  symmetric-deletion index that only looks at likely matches instead of
  scoring every query x reference pair. There is no install check and no
  pure-Python matcher.
- **Every gated pair through numpy.** symscan never scores a pair beyond
  the threshold, so the pairs file of a Hamming or Levenshtein run is
  scored separately (`Distances.py`): references are indexed by gene and
  gate key so each query only meets the pairs its gates allow, then a
  block of pairs is compared position by position (Hamming, `X` never
  counting) or through an edit distance table filled for the whole block
  at once (Levenshtein). Queries go in blocks, so the file is written as it
  goes, and its `within_threshold` pairs are exactly symscan's hits.
- **`X` through symscan.** symscan has no wildcard and compares `X` as a
  letter, so the `X`s are taken out of the comparison instead: reference
  rows are grouped by where their `X`s sit, and those positions are cut out
  of them and of every same-length query before symscan sees them. A
  query's own `X`s (rarer) are covered by slack -- one extra difference
  each -- and every pair involving an `X` is re-scored with the wildcard
  rule, keeping only those within the real threshold.
- **Rows beyond symscan's limits are scored directly.** symscan aborts the
  process above 14 differences, or once its index -- every deletion variant
  of every sequence, which grows steeply with length and tolerance --
  outgrows memory. Both the queries and the references are chunked to a
  fixed budget, so table size alone never matters; what cannot fit at all
  (a threshold allowing over 14 differences, `X` slack included, or one
  sequence too long to index at that threshold) is kept away from symscan,
  and those query rows' pairs are scored directly instead (`Distances.py`),
  so no row misses a hit; their count is in `annotation_metadata.yml`. On
  today's IEDB table that is a handful of outliers, such as reference
  "CDR3s" of 60+ residues at a loose threshold. A VDJ comparison (~110-130
  residues) fits symscan up to 4 differences -- exact matching and
  thresholds of about 0.97 and up; below that most VDJ rows are scored
  directly, which is correct but slower, so `--partition-by` matters there.
  (Under `--paired`, which is always exact, only a query with over 14 `X`s
  goes beyond symscan; such a row is skipped, with a warning.)
- **Order of steps.** symscan measures distance first; the exact-column and
  gene checks are applied to the pairs it returns.

## Current limits

- Only IEDB's `bcr` table comes in AIRR so far (IEDB publishes it that
  way); CoV-AbDab and ASD are still downloaded in their own columns, so they
  cannot be used as a reference until they are converted.
- Antibodies only: TCR tables are not supported yet (no `IGH` rows, so they
  are refused).
- No substitution-matrix (BLOSUM) scoring.
- Fuzzy VDJ matching below a threshold of about 0.97 is out of symscan's
  reach; those rows are scored directly, which is slow without
  `--partition-by`.
