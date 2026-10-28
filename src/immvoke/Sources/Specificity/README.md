# Specificity annotations

Annotate your own receptors with antigen specificity from public databases:
download a reference table in AIRR format, then match your receptors against
it. How the matching works, step by step, is in [idea.md](idea.md).

## IEDB: download

[IEDB](https://www.iedb.org/) (`Iedb.py`) publishes its BCR (antibody)
receptors in AIRR format, as `bcr_full_airr.tsv` inside its bulk
`airr_full.zip` export, refreshed weekly. `immvoke specificity iedb download
bcr` downloads that table, keeping **every receptor species** unless
`--air_species` names one. With `--air_species "Homo sapiens"`, a receptor
is kept when it has a human chain (`locus_species` `NCBITAXON:9606`) and no
chain of another species -- a chain with no species recorded rides along
with its partner, so pairs are never split. A species includes its strains:
`"Mus musculus"` also keeps BALB/c, C57BL/6, Xenomouse and the rest, while
`"Mus musculus BALB/c"` keeps that strain alone. The names are IEDB's;
`immvoke specificity iedb download bcr --help` lists them all. Likewise
every assay outcome is kept unless `--outcome` names one (`Positive` or
`Negative`, matched against `iedb_assay_outcome`), and every antigen unless
`--antigen_species` names its organism: any name IEDB uses, ignoring case,
with or without its bracketed common name (`"Homo sapiens"` for
`Homo sapiens (human)`). An organism includes those under it, so
`"Influenza A virus"` keeps every flu strain. The AIRR export does not say
which organism an epitope comes from, so this asks IEDB's API about the
table's epitopes at download time (about 15 seconds for bcr); an unknown name fails,
naming the closest IEDB names. The rest is
kept as is -- one
rearrangement record per chain, linked by `cell_id`,
built from IEDB's calculated (IMGT standardized) gene calls and CDRs, with
the epitope and assay context (`iedb_epitope_name`, `iedb_assay_ids`, ...)
alongside -- and writes it with a `bcr.airr_validation.log` report next to it.

```bash
immvoke specificity iedb download bcr --air_species "Homo sapiens" \
    --antigen_species "Influenza A virus" --outcome Positive --outdir databases
```

Every download lands in its database's own folder under `--outdir`, so
more databases can share one outdir without mixing their files:

```
databases/
└── iedb/
    ├── bcr.tsv                  # the reference table for annotate
    ├── bcr.airr_validation.log  # AIRR validation report
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
downloaded databases, each named after the reference's database and table,
`<db>_<table>_`:

```
annotation/exact/
├── iedb_bcr_annotation.tsv            # one row per query: its hits
├── iedb_bcr_annotation_pairs.tsv      # one row per (query, reference) pair, scored
└── iedb_bcr_annotation_metadata.yml   # the options, query, reference and its
                                       # database, counts, version, time and command
```

The database is read off the `download_metadata.yml` that `immvoke
specificity` writes beside the reference table, and the table off its file
name, so runs against different databases can share one `--outdir` without
overwriting each other. A reference without that record is named by its
file name alone (`myref.tsv` gives `myref_annotation.tsv`), and a DataFrame
passed from Python by `reference_`.

`<db>_<table>_annotation.tsv`, without `--paired`, has one row per query
heavy chain row, including rows with no hit:

| Column | Meaning |
| --- | --- |
| `sequence_id` | the query's own `sequence_id`, so the table joins straight back to the query file (`query_row`, its position among the heavy chain rows, for a query without that column) |
| `n_hits_total` | how many reference rows matched |
| `min_mismatches` | the closest hit's distance (edit distance under `levenshtein`) |
| `is_unique_hit` | 0 or 1 hit |
| `hit_ids` | the matching reference `sequence_id`s, closest first, `;`-joined |
| `mismatch_positions` | Hamming only: the differing positions per hit |

With `--paired`, one row per query antibody (`cell_id`): `classification`
(`both`, `heavy_only`, `light_only`, `conflicting` or `neither`),
`heavy_hit_ids`, `light_hit_ids` and `both_hit_ids`.

`<db>_<table>_annotation_pairs.tsv` has one (query, reference) pair per
row, each with its own score, by query and closest first:

- an exact match (no `--identity-threshold`) lists its hits;
- a Hamming or Levenshtein run lists **every pair through the filters you
  chose** -- `--partition-by`'s gene gate (and CDR3 length under `VJL`), no
  blank sequences, and under Hamming same-length pairs only -- whatever its
  distance, with `within_threshold` marking the hits. Without
  `--partition-by` that is every query against every reference: for one
  repertoire against IEDB, millions of rows under Hamming and some 10^8
  under Levenshtein, so `--partition-by VJ` keeps it to a manageable size.

| Column | Meaning |
| --- | --- |
| `sequence_id` | the query, as above |
| `reference_id` | the reference `sequence_id` |
| `mismatches` | this pair's distance (edit distance under `levenshtein`) |
| `identity` | `1 - mismatches / query length` (an `X` counts as identical under `hamming`) |
| `within_threshold` | Hamming/Levenshtein runs: whether the pair is a hit |
| `mismatch_positions` | Hamming only: this pair's differing positions |

With `--paired`, the hits: `cell_id`, `hit_id` (a reference `cell_id`) and
`matched_chains` (`both`, `heavy` or `light`).

IEDB lists a receptor once per epitope it was tested against: each record
gets its own `cell_id`, but the chains keep the receptor's `sequence_id`, so
one id in `hit_ids` can stand for several reference rows. Join `hit_ids`
back to the reference on `sequence_id` to see every epitope behind a hit.

Reference command: `immvoke specificity iedb download bcr --help`,
`immvoke annotate --help`.
