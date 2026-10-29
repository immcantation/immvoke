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
chain of another species, and only its human chains are written: a chain
with no species recorded is left out (`search` would show it as
`(not recorded)`), so a heavy chain whose light chain IEDB gave no species
is kept alone. A species includes its strains:
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
naming the closest IEDB names. Every epitope structure is kept unless
`--antigen_structure_type` names one: `Linear peptide`, `Discontinuous
peptide` (which includes `Discontinuous peptide on multi chain`) or
`Non-peptidic` (haptens, glycans, lipids...) -- IEDB's `structure_type`,
the same choice as its own search's Epitope Structure, looked up in the
same API call. `--assay_type` keeps records with an assay of that type
(`neutralization`, `3D structure`, `qualitative binding`, ... -- the first
part of IEDB's assay name), and `--disease` those with an assay studying
that disease or one under it (`"viral infectious disease"` takes in
COVID-19, influenza and the rest), both any name IEDB uses, ignoring case;
given both, one assay must satisfy both. They ask IEDB's API about the
table's assays, adding about 40 seconds, and only when given. IEDB records
no disease for most assays, so for a pathogen `--antigen_species` finds
more (SARS-CoV2 neutralizing human antibodies: 63 receptor groups by
`--antigen_species SARS-CoV2`, 39 by `--disease COVID-19`). `--paired_only`
keeps only receptors with both chains of a pair -- heavy with kappa or
lambda, beta with alpha, delta with gamma -- naming a chain IEDB gave no
locus by its V or J gene. The rest is
kept as is -- one
rearrangement record per chain, linked by `cell_id`,
built from IEDB's calculated (IMGT standardized) gene calls and CDRs, with
the epitope and assay context (`iedb_epitope_name`, `iedb_assay_ids`, ...)
alongside -- and writes it as `iedb_bcr.tsv` with an `iedb_bcr.airr_validation.log` report next to it.

```bash
immvoke specificity iedb download bcr --air_species "Homo sapiens" \
    --antigen_species "Influenza A virus" --antigen_structure_type "Discontinuous peptide" \
    --outcome Positive --outdir databases

# neutralizing, paired human antibodies
immvoke specificity iedb download bcr --air_species "Homo sapiens" \
    --assay_type neutralization --paired_only --outdir databases
```

`search` takes the same filters and counts what a download would give,
writing nothing (unless `--out` asks for the counts as a TSV): rows,
receptor groups, epitopes, and the most common loci, receptor species,
antigen species, antigen structure types and outcomes (`--top` of each, 10
by default).

```bash
immvoke specificity iedb search bcr --air_species "Homo sapiens"
```

```
iedb bcr: BCR (antibody) receptor sequences, from the AIRR bulk export
filters: air_species = Homo sapiens

  rows             10,861 of 18,650
  receptor groups  3,417
  epitopes         3,128

  locus
    5,917  IGH
    2,920  IGK
    1,747  IGL
      277  (not recorded)
...
```

IEDB's bulk ZIP is kept in `~/.cache/immvoke/iedb/` (`$XDG_CACHE_HOME` if
set) for a day, so a search, a download, and `bcr` then `tcr`, share one
copy; after a day the current export is fetched again.

Every download lands in its database's own folder under `--outdir`, so
more databases can share one outdir without mixing their files:

```
databases/
└── iedb/
    ├── iedb_bcr.tsv                  # the reference table for annotate
    ├── iedb_bcr.airr_validation.log  # AIRR validation report
    └── download_metadata.yml         # source, license, citation, when, command
```

Each record also gets its epitope's antigen, that antigen's source
species and the epitope's structure type, `iedb_antigen`,
`iedb_antigen_source_species` and `iedb_antigen_structure_type` (after
`iedb_epitope_name`), which IEDB's AIRR export lacks: the download looks the
table's epitopes up on IEDB's API (`curated_source_antigens`, several
` | ` joined, and `structure_type`), adding about 15 seconds.

The written `databases/iedb/iedb_bcr.tsv` keeps every column IEDB provides and
always carries every field the AIRR schema requires; one IEDB lacks is
written blank.

## Annotate

`immvoke annotate` matches a query AIRR table -- your own receptors --
against a reference AIRR table, such as the `databases/iedb/iedb_bcr.tsv`
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
| `--full-dist-matrix` | `full_dist_matrix` | flag | also write every hit within the threshold, not only the best ones (see Output) |

`--partition-by` is off by default, but an unpartitioned CDR3 comparison at
a loose threshold is prone to coincidental hits, especially on short CDR3s,
so `VJ` is usually what you want there.

### Examples

```bash
# exact CDR3 match
immvoke annotate --query mine.tsv --reference databases/iedb/iedb_bcr.tsv \
    --outdir annotation/exact

# exact full VDJ match within the same V- and J-gene
immvoke annotate --query mine.tsv --reference databases/iedb/iedb_bcr.tsv \
    --sequence-region VDJ --partition-by VJ --outdir annotation/vdj_exact

# same V+J gene, CDR3 at 85% amino acid identity (Hamming distance)
immvoke annotate --query mine.tsv --reference databases/iedb/iedb_bcr.tsv \
    --partition-by VJ --identity-threshold 0.85 --outdir annotation/hamming_0.85

# CDR1-3 pooled at 90% identity, same V-gene, J-gene and CDR3 length
immvoke annotate --query mine.tsv --reference databases/iedb/iedb_bcr.tsv \
    --sequence-region CDRs --partition-by VJL --identity-threshold 0.9 \
    --outdir annotation/cdrs_0.9

# CDR3 at 90% identity allowing insertions/deletions
immvoke annotate --query mine.tsv --reference databases/iedb/iedb_bcr.tsv \
    --partition-by VJ --identity-threshold 0.9 --distance-metric levenshtein \
    --outdir annotation/levenshtein_0.9

# paired heavy+light match
immvoke annotate --query mine.tsv --reference databases/iedb/iedb_bcr.tsv \
    --paired --outdir annotation/paired
```

### Output

Each run writes two files into its `--outdir`, kept apart from the
downloaded databases, each named after the reference's database and table,
`<db>_<table>_`, and a third with `--full-dist-matrix`:

```
annotation/exact/
├── iedb_bcr_annotation.tsv               # one row per query: its best hits
├── iedb_bcr_annotation_metadata.yml      # the options, query, reference and its
│                                         # database, counts, version, time and command
└── iedb_bcr_annotation_full_matrix.txt   # --full-dist-matrix: every hit, one per row
```

The database is read off the `download_metadata.yml` that `immvoke
specificity` writes beside the reference table, and the table off its file
name, so runs against different databases can share one `--outdir` without
overwriting each other. A reference without that record is named by its
file name alone (`myref.tsv` gives `myref_annotation.tsv`), and a DataFrame
passed from Python by `reference_`.

`<db>_<table>_annotation.tsv`, without `--paired`, has one row per query
heavy chain row, including rows with no hit. Of a query's hits it keeps the
**best**: the ones with the highest identity, all of them when several tie.

| Column | Meaning |
| --- | --- |
| `sequence_id` | the query's own `sequence_id`, so the table joins straight back to the query file (`query_row`, its position among the heavy chain rows, for a query without that column) |
| `n_hits` | how many best hits (0 without a hit) |
| `identity` | their identity, `1 - mismatches / query length` (an `X` counts as identical under `hamming`) |
| `hit_ids` | their reference `sequence_id`s, `;`-joined |
| `antigen` | each hit's antigen (`iedb_antigen`), in `hit_ids` order |
| `epitope` | each hit's epitope (`iedb_epitope_name`) |
| `antigen_source_species` | the organism each hit's antigen comes from (`iedb_antigen_source_species`) |

Hits are separated by `;` rather than commas, since IEDB names a
discontinuous epitope by its residues (`G58, L59, S61`). IEDB lists a
receptor once per epitope it was tested against, all under its
`sequence_id`, so one hit can stand for several epitopes; their values are
joined with ` | ` within that hit's entry. A reference without IEDB's
columns leaves those three blank.

With `--paired`, one row per query antibody (`cell_id`): `classification`
(`both`, `heavy_only`, `light_only`, `conflicting` or `neither`),
`heavy_hit_ids`, `light_hit_ids` and `both_hit_ids`.

`<db>_<table>_annotation_full_matrix.txt`, with `--full-dist-matrix`, lists
every hit within `--identity-threshold`, not only the best ones: one row
per (query, hit), by query and closest first -- `sequence_id`, `hit_id`
and `identity`. With `--paired`, each antibody's hits: `cell_id`, `hit_id`
(a reference `cell_id`) and `matched_chains` (`both`, `heavy` or `light`).

Reference command: `immvoke specificity iedb download bcr --help`,
`immvoke annotate --help`.
