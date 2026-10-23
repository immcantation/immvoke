# Specificity annotations

Annotate your own receptors with antigen specificity from public databases:
download a reference table in AIRR format, then match your receptors against
it. How the matching works, step by step, is in [idea.md](idea.md).

## IEDB: download

[IEDB](https://www.iedb.org/) (`Iedb.py`) publishes its BCR (antibody)
receptors in AIRR format, as `bcr_full_airr.tsv` inside its bulk
`airr_full.zip` export, refreshed weekly. `immvoke specificity bcr` downloads
that table, keeping **human receptors only**: a receptor is kept when it has
a human chain (`locus_species` `NCBITAXON:9606`) and no chain of another
species -- a chain with no species recorded rides along with its human
partner, so pairs are never split. The rest is kept as is -- one
rearrangement record per chain, linked by `cell_id`,
built from IEDB's calculated (IMGT standardized) gene calls and CDRs, with
the epitope and assay context (`iedb_epitope_name`, `iedb_assay_ids`, ...)
alongside -- and writes it with a `.validation.txt` report next to it.

```bash
immvoke specificity bcr --outdir databases
```

Every download lands in its database's own folder under `--outdir`
(`databases` by default for `specificity bcr`), so more databases can share
one outdir without mixing their files:

```
databases/
└── iedb/
    ├── bcr.tsv                  # the reference table for annotate
    ├── bcr.tsv.validation.txt   # AIRR validation report
    ├── download_metadata.yml    # source, license, citation, when, command
    └── raw/                     # IEDB's airr_full.zip and the file taken from it
```

The written `databases/iedb/bcr.tsv` keeps every column IEDB provides and
always carries every field the AIRR schema requires; one IEDB lacks is
written blank.

## Annotate

`immvoke annotate` matches a query AIRR table -- your own receptors --
against a reference AIRR table, such as the `databases/iedb/bcr.tsv`
downloaded above, and reports which reference rows each query row matches.
The reference must be an AIRR rearrangement table; see
[idea.md](idea.md#input-contract) for exactly what that requires. From
Python, `immvoke.Sources.Specificity.Annotate.annotateSpecificity()` takes
the same options, as file paths or DataFrames, returns the result table,
and, given `outdir`, writes the same files the command does.

### Options

| `immvoke annotate` | `annotateSpecificity()` | Values | Meaning |
| --- | --- | --- | --- |
| `--query`, `--reference` | `query`, `reference` | AIRR TSV (or DataFrame in Python) | the tables to match |
| `--sequence-region` | `sequence_region` | `CDR3H` (default), `VDJ`, `CDRs` | compare the heavy chain CDR3 (`cdr3_aa`); `sequence_aa`; or CDR1-3 pooled |
| `--partition-by` | `partition_by` | unset (default), `V`, `VJ`, `VJL` | only compare pairs with the same V-gene; V- and J-gene; or V-gene, J-gene and CDR3 length (same gene, any allele) |
| `--identity-threshold` | `identity_threshold` | unset (exact), or a fraction in (0, 1] | e.g. `0.9` allows up to 10% of the query's positions to differ |
| `--distance-metric` | `distance_metric` | `hamming` (default), `levenshtein` | substitutions only, or substitutions, insertions and deletions |
| `--paired` | `paired` | flag | match heavy and light chain of the same reference antibody |
| `--count-field` | `count_field` | an AIRR count column; `duplicate_count` (default) | `--paired`: which of a query cell's several heavy (or light) chains to keep -- the highest count |
| `--threads` | `threads` | `1` (default); `-1` for every core | symscan's threads |
| `--outdir` | `outdir` | folder, created if needed | where to write this run's files (see Output) |

`--partition-by` is off by default, but an unpartitioned CDR3 comparison at
a loose threshold is prone to coincidental hits, especially on short CDR3s,
so `VJ` is usually what you want there.

### Examples

```bash
# exact CDR3 match
immvoke annotate --query mine.tsv --reference databases/iedb/bcr.tsv \
    --outdir annotation/exact

# exact full VDJ match within the same V- and J-gene
immvoke annotate --query mine.tsv --reference databases/iedb/bcr.tsv \
    --sequence-region VDJ --partition-by VJ --outdir annotation/vdj_exact

# same V+J gene, CDR3 at 85% amino acid identity (Hamming distance)
immvoke annotate --query mine.tsv --reference databases/iedb/bcr.tsv \
    --partition-by VJ --identity-threshold 0.85 --outdir annotation/hamming_0.85

# CDR1-3 pooled at 90% identity, same V-gene, J-gene and CDR3 length
immvoke annotate --query mine.tsv --reference databases/iedb/bcr.tsv \
    --sequence-region CDRs --partition-by VJL --identity-threshold 0.9 \
    --outdir annotation/cdrs_0.9

# CDR3 at 90% identity allowing insertions/deletions
immvoke annotate --query mine.tsv --reference databases/iedb/bcr.tsv \
    --partition-by VJ --identity-threshold 0.9 --distance-metric levenshtein \
    --outdir annotation/levenshtein_0.9

# paired heavy+light match
immvoke annotate --query mine.tsv --reference databases/iedb/bcr.tsv \
    --paired --outdir annotation/paired
```

### Output

Each run writes three files into its `--outdir`, kept apart from the
downloaded databases:

```
annotation/exact/
├── annotation.tsv            # one row per query: its hits
├── annotation_pairs.tsv      # one row per (query, reference) pair, scored
└── annotation_metadata.yml   # the options, query, reference and its source
                              # database, counts, version, time and command
```

`annotation_metadata.yml` names the reference's database (`source: iedb`)
by reading the `download_metadata.yml` beside the reference table, so a run
stays traceable to the database it was annotated against.

`annotation.tsv`, without `--paired`, has one row per query heavy chain
row, including rows with no hit:

| Column | Meaning |
| --- | --- |
| `query_row` | the query row's position among the heavy chain rows |
| `n_hits_total` | how many reference rows matched |
| `min_mismatches` | the closest hit's distance (edit distance under `levenshtein`) |
| `is_unique_hit` | 0 or 1 hit |
| `hit_ids` | the matching reference `sequence_id`s, closest first, `;`-joined |
| `mismatch_positions` | Hamming only: the differing positions per hit |

With `--paired`, one row per query antibody (`cell_id`): `classification`
(`both`, `heavy_only`, `light_only`, `conflicting` or `neither`),
`heavy_hit_ids`, `light_hit_ids` and `both_hit_ids`.

`annotation_pairs.tsv` lists the same hits one (query, reference) pair per
row, closest first, each with its own score; a query with no hit has no
row. Without `--paired`:

| Column | Meaning |
| --- | --- |
| `query_row` | as above, to join back to the result table |
| `hit_id` | one matching reference `sequence_id` |
| `mismatches` | this pair's distance (edit distance under `levenshtein`) |
| `identity` | `1 - mismatches / query length` (an `X` counts as identical under `hamming`) |
| `mismatch_positions` | Hamming only: this pair's differing positions |

With `--paired`: `cell_id`, `hit_id` (a reference `cell_id`) and
`matched_chains` (`both`, `heavy` or `light`).

IEDB lists a receptor once per epitope it was tested against: each record
gets its own `cell_id`, but the chains keep the receptor's `sequence_id`, so
one id in `hit_ids` can stand for several reference rows. Join `hit_ids`
back to the reference on `sequence_id` to see every epitope behind a hit.

Reference command: `immvoke specificity bcr --help`,
`immvoke annotate --help`.
