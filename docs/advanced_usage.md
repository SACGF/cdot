# Advanced usage

Two things you'll likely want once you're past the [basic examples](../README.md#examples):

1. [Fixing messy HGVS input](#fixing-messy-hgvs-input) before handing it to biocommons HGVS.
2. [Read-ahead batch retrieval](#read-ahead-batch-retrieval) to make bulk processing fast over the
   REST data provider.

## Fixing messy HGVS input

Real-world HGVS strings — from spreadsheets, lab reports, literature, decades-old archives — are
frequently *almost* valid: stray whitespace, a `p.` suffix glued on, an uppercase `DEL`, a missing
underscore (`NM000059`), or just a gene symbol instead of a transcript (`BRCA2:c.36del`). biocommons
HGVS will reject these. cdot's `cdot.hgvs` module fixes the common cases and tells you what it changed.

### `fix_hgvs` — the one-call entry point

`fix_hgvs()` is the recommended function. It chains string cleaning and (optionally) gene→transcript
resolution, returning the cleaned string plus a list of `HGVSFix` records describing every change:

```python
from cdot.hgvs import fix_hgvs
from cdot.hgvs.dataproviders import JSONDataProvider

hdp = JSONDataProvider(["./cdot-0.2.32.refseq.GRCh38.json.gz"])

# Gene-only input → resolved to the MANE/canonical transcript
result, fixes = fix_hgvs("BRCA2:c.36DEL", hdp, "GRCh38")
# result = "NM_000059.4:c.36del"
# fixes  = [HGVSFix(WARNING, LOWERCASED_MUTATION_TYPE, ...),
#           HGVSFix(WARNING, RESOLVED_GENE_TO_TRANSCRIPT, ...)]

for fix in fixes:
    print(fix.severity.value, fix.code.value, fix.message)
```

Each [`HGVSFix`](../cdot/hgvs/clean.py) has a `severity` (`WARNING` = fixed something unambiguously,
`ERROR` = could not fix), a `code` (`HGVSFixCode`), a human-readable `message`, and the
`original`/`fixed` values. Inspect the list to decide whether to trust the result, log a warning,
or reject the variant.

**Cleaning only (no data provider).** If your input already names a transcript, omit the provider and
build — only string normalisation runs, so no data load is needed:

```python
result, fixes = fix_hgvs("NM_000059.4 c.316+5G>A")
# result = "NM_000059.4:c.316+5G>A"
```

Useful options:

| Argument | Default | Effect |
|----------|---------|--------|
| `raise_on_errors` | `False` | Raise `HGVSInputError` on the first ERROR-level fix instead of returning it. |
| `prefer_consortium` | `Consortium.REFSEQ` | Hard-filter resolved transcripts to RefSeq (default) or `ENSEMBL`; pass `None` to allow either. |
| `fallback_to_longest` | `False` | If no tagged/canonical transcript exists for a gene, fall back to the longest. |
| `tag_priority` | MANE_Select > MANE_Plus_Clinical > RefSeq_Select > Ensembl_canonical | Order used when resolving a gene to a transcript. |
| `ops` | `None` (all) | Restrict which cleaning operations run (see below). |
| `version_fallback` | `None` (off) | Opt in to adjacent transcript-version fallback (see below). |

### Transcript version fallback

Real-world HGVS often names a transcript version your data release doesn't have (e.g. a report cites
`NM_000059.2` but you only ship `.3`/`.4`). This is **off by default** — pass a `VersionStrategy` to
opt in. When set, `fix_hgvs` substitutes the best available version and reports a WARNING:

```python
from cdot.hgvs import fix_hgvs
from cdot.hgvs.clean import VersionStrategy

result, fixes = fix_hgvs("NM_000059.2:c.36del", hdp,
                         version_fallback=VersionStrategy.UP_THEN_DOWN)
# result = "NM_000059.4:c.36del"  (if .2 absent but .4 present)
# fixes[-1].code == HGVSFixCode.USED_ADJACENT_VERSION
```

`version_fallback` only needs a `data_provider` (not a `genome_build`). The strategies are:

| `VersionStrategy` | Behaviour |
|-------------------|-----------|
| `UP_THEN_DOWN` | Prefer higher versions first, then lower (default for the helper). |
| `CLOSEST` | Nearest version by absolute distance (ties prefer the higher one). |
| `LATEST` | Always the highest available, ignoring the requested version. |

If the requested version exists it's left untouched (non-destructive); if the accession has no
versions at all an ERROR fix is returned. Call `resolve_transcript_version()` directly for the same
behaviour outside `fix_hgvs`. Both rely on the data provider's `get_tx_versions(accession)`
(implemented by `JSONDataProvider` and `RESTDataProvider`; the REST version uses the versionless
`/transcript/<ac>` lookup and warms the cache with every returned version).

### `clean_hgvs` — string cleaning with op selection

`fix_hgvs` calls `clean_hgvs()` for the string-fixing half. Call it directly when you want only the
pure-string normalisation (no provider, no gene resolution) and fine control over which operations
run via the `ops` set of `HGVSCleanOp`:

```python
from cdot.hgvs.clean import clean_hgvs, HGVSCleanOp, ALL_CLEAN_OPS

# Allowlist: only strip whitespace
cleaned, fixes = clean_hgvs(messy, ops={HGVSCleanOp.STRIP_WHITESPACE})

# Blocklist via set algebra: everything except stripping the protein suffix
cleaned, fixes = clean_hgvs(messy, ops=ALL_CLEAN_OPS - {HGVSCleanOp.STRIP_PROTEIN_SUFFIX})
```

Ops always run in the canonical pipeline order regardless of set order. Pass `validate=False` to skip
the ERROR-level validation pass when you only want normalisation.

## Read-ahead batch retrieval

The biocommons HGVS `Interface` is **one transcript at a time**: each `c_to_g` (etc.) triggers a
single `/transcript/<ac>` lookup that's then reused for `get_tx_info` / `get_tx_exons` / .... For a
local `JSONDataProvider` that's all in memory, so it's already fast. But for the `RESTDataProvider`,
processing N variants means N sequential HTTP round-trips — the latency dominates.

When you know the transcripts up front (bulk processing a file of variants), warm the cache first
with `prefetch()`. Every later lookup is then a cache hit — **no change to how you call biocommons HGVS**:

```python
from cdot.hgvs.dataproviders import RESTDataProvider

hdp = RESTDataProvider()  # cdotlib.org

# Warm the cache in one round-trip
hdp.prefetch(["NM_001637.3", "NM_000059.4", "NM_007294.4"])

# ... now run your AssemblyMapper / c_to_g loop as usual; all cache hits
```

### From a list of HGVS strings

If you have HGVS strings rather than bare accessions, `prefetch_from_hgvs()` extracts the accessions
for you (running `clean_hgvs` first so messy input still yields a usable accession). Gene-only HGVS
and genomic `NC_` references are skipped:

```python
variants = ["NM_001637.3:c.1582G>A", "NM_000059.4:c.316+5G>A", "BRCA2:c.36del"]
hdp.prefetch_from_hgvs(variants)   # prefetches NM_001637.3 and NM_000059.4; skips the gene-only one
```

### How it fetches (batch vs concurrent)

`prefetch()` has two strategies:

- **Batch (default).** POSTs the whole accession list to the cdot_rest batch `POST /transcripts`
  endpoint in a **single round-trip**. Versionless accessions (e.g. `NM_000059`) are expanded
  server-side to *every* available version — handy for warming the version-bump path.
- **Concurrent fallback.** Servers without the batch endpoint (older cdot_rest) automatically fall
  back to a thread-pool of single `/transcript/<ac>` requests, turning N sequential round-trips into
  `ceil(N / max_workers)` waves. Force this with `batch=False`; tune the pool with `max_workers`
  (default 10). Note: the fallback does **not** expand versionless accessions — pass exact versions.

```python
hdp.prefetch(accessions, batch=False, max_workers=20)
```

Already-cached accessions are skipped, missing/404 transcripts are cached as `None` (matching normal
lookup behaviour, so a prefetch miss is never fatal), and the call returns the number of transcripts
added to the cache.

> Batch prefetch needs a cdot_rest server with the `POST /transcripts` endpoint
> ([SACGF/cdot_rest#9](https://github.com/SACGF/cdot_rest/pull/9)). Older servers transparently use
> the concurrent fallback.

## See also

- [Biocommons HGVS examples](examples_biocommons.md)
- [JSON data format reference](json_data_format.md)
- [Coordinates & exon alignments](coordinates_and_exons.md)
