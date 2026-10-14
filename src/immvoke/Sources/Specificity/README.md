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
immvoke specificity bcr --outdir tmp
```

The written `tmp/specificity/bcr.tsv` keeps every column IEDB provides and
always carries every field the AIRR schema requires; one IEDB lacks is
written blank.

## Annotate

`immvoke annotate` matches a query AIRR table -- your own receptors --
against a reference AIRR table, such as the `tmp/specificity/bcr.tsv`
written above, and reports which reference rows each query row matches. The
reference must be an AIRR rearrangement table; see
[idea.md](idea.md#input-contract) for exactly what that requires. From
Python, `immvoke.Sources.Specificity.Annotate.annotateSpecificity()` takes
the same options, as file paths or DataFrames, and returns the result table.

### Options

| `immvoke annotate` | `annotateSpecificity()` | Values | Meaning |
| --- | --- | --- | --- |
| `--query`, `--reference` | `query`, `reference` | AIRR TSV (or DataFrame in Python) | the tables to match |
| `--sequence-region` | `sequence_region` | `CDR3` (default), `VDJ`, `CDRs` | compare `cdr3_aa`; `sequence_aa`; or CDR1-3 pooled |
| `--partition-by` | `partition_by` | unset (default), `V`, `VJ`, `VJL` | only compare pairs with the same V-gene; V- and J-gene; or V-gene, J-gene and CDR3 length (same gene, any allele) |
| `--identity-threshold` | `identity_threshold` | unset (exact), or a fraction in (0, 1] | e.g. `0.9` allows up to 10% of the query's positions to differ |
| `--distance-metric` | `distance_metric` | `hamming` (default), `levenshtein` | substitutions only, or substitutions, insertions and deletions |
| `--paired` | `paired` | flag | match heavy and light chain of the same reference antibody |
| `--threads` | `threads` | `1` (default); `-1` for every core | symscan's threads |
| `-o`, `--out` | `output` | path | where to write the result TSV |

`--partition-by` is off by default, but an unpartitioned CDR3 comparison at
a loose threshold is prone to coincidental hits, especially on short CDR3s,
so `VJ` is usually what you want there.

### Examples

```bash
# exact CDR3 match
immvoke annotate --query mine.tsv --reference tmp/specificity/bcr.tsv \
    -o annotated.tsv

# exact full VDJ match within the same V- and J-gene
immvoke annotate --query mine.tsv --reference tmp/specificity/bcr.tsv \
    --sequence-region VDJ --partition-by VJ -o annotated.tsv

# same V+J gene, CDR3 at 85% amino acid identity (Hamming distance)
immvoke annotate --query mine.tsv --reference tmp/specificity/bcr.tsv \
    --partition-by VJ --identity-threshold 0.85 -o annotated.tsv

# CDR1-3 pooled at 90% identity, same V-gene, J-gene and CDR3 length
immvoke annotate --query mine.tsv --reference tmp/specificity/bcr.tsv \
    --sequence-region CDRs --partition-by VJL --identity-threshold 0.9 \
    -o annotated.tsv

# CDR3 at 90% identity allowing insertions/deletions
immvoke annotate --query mine.tsv --reference tmp/specificity/bcr.tsv \
    --partition-by VJ --identity-threshold 0.9 --distance-metric levenshtein \
    -o annotated.tsv

# paired heavy+light match
immvoke annotate --query mine.tsv --reference tmp/specificity/bcr.tsv \
    --paired -o annotated.tsv
```

### Output

Without `--paired`, one row per query heavy chain row, including rows with
no hit:

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

IEDB lists a receptor once per epitope it was tested against: each record
gets its own `cell_id`, but the chains keep the receptor's `sequence_id`, so
one id in `hit_ids` can stand for several reference rows. Join `hit_ids`
back to the reference on `sequence_id` to see every epitope behind a hit.

Reference command: `immvoke specificity bcr --help`,
`immvoke annotate --help`.
