# cdot — Deep Research Notes

## Project Identity

- **Package name:** `cdot` ("Complete Dict Of Transcripts")
- **Version:** `0.2.26` (code); data schema version `0.2.33`
- **GitHub:** https://github.com/SACGF/cdot
- **Author:** Dave Lawrence (SACGF)
- **Purpose:** Convert RefSeq/Ensembl GTF/GFF3 annotation files into compact JSON, then serve that data to the two major Python HGVS variant-annotation libraries: [biocommons/hgvs](https://github.com/biocommons/hgvs) and [pyhgvs](https://github.com/counsyl/hgvs).
- **Performance:** Local JSON.gz gives ~500–1000 transcripts/s vs UTA public DB at ~1–1.5 s/transcript.

---

## Repository Layout

```
cdot/                          # Installed Python package
  __init__.py                  # __version__, get_data_schema_int()
  data_release.py              # GitHub API helpers for downloading data releases
  hgvs/dataproviders/
    __init__.py                # Re-exports all providers
    json_data_provider.py      # AbstractJSONDataProvider, LocalDataProvider,
                               #   JSONDataProvider, RESTDataProvider
    ensembl_tark_data_provider.py  # EnsemblTarkDataProvider + EnsemblTarkSeqFetcher
    seqfetcher.py              # PrefixSeqFetcher, ChainedSeqFetcher,
                               #   VerifyMultipleSeqFetcher, AlwaysFailSeqFetcher
    fasta_seqfetcher.py        # GenomeFastaSeqFetcher,
                               #   ExonsFromGenomeFastaSeqFetcher, FastaSeqFetcher
  pyhgvs/
    pyhgvs_transcript.py       # AbstractPyHGVSTranscriptFactory + concrete classes

generate_transcript_data/      # NOT part of the pip package
  cdot_json.py                 # CLI: gtf_to_json, gff3_to_json, uta_to_json,
                               #   merge_historical, combine_builds, release_notes
  gff_parser.py                # GFFParser (base), GTFParser, GFF3Parser
  cdot_gene_info.py            # NCBI Entrez gene-summary fetcher
  json_schema_version.py       # JSON_SCHEMA_VERSION constant + changelog comments
  json_encoders.py             # SortedSetEncoder (sets → sorted lists)
  cdot_transcripts.yaml        # All known GTF/GFF3 source URLs + build mappings
  ensembl_grch37_canonical_transcripts.py  # Canonical transcript fetcher

tests/
  test_gff_parsers.py
  test_json_data_provider_refseq.py
  test_json_data_provider_ensembl.py
  test_pyhgvs.py
  test_uta_conversion.py
  benchmark_hgvs.py
  mock_seqfetcher.py
  mock_ensembl_tark.py
  genome.py                    # MockGenomeTestFile
  test_data/                   # Small JSON, GTF, GFF fixtures
```

---

## JSON Data Format (Schema v0.2.33)

The on-disk format is gzipped JSON with this top-level shape:

```json
{
  "cdot_version": "0.2.33",
  "genome_builds": ["GRCh37", "GRCh38"],
  "transcripts": { "<id>": { ... } },
  "genes": { "<gene_version_id>": { ... } },
  "refseq_gene_summary_api_retrieval_date": "...",
  "metadata": {
    "method": "...",
    "input_files": [],
    "sys.argv": "...",
    "url_counts": {},
    "cdot_version": "...",
    "genome_builds": []
  }
}
```

### Transcript object

Fields shared across all genome builds (top-level in the transcript dict):
- `id` — accession (e.g. `"NM_001637.3"`)
- `gene_name` — HGNC symbol
- `gene_version` — GeneID (RefSeq) or Ensembl gene version
- `protein` — protein accession, or absent if non-coding
- `start_codon`, `stop_codon` — 0-based cDNA positions of codon start/end
- `biotype` — set → sorted list, e.g. `["protein_coding"]`
- `source` — GTF column 2 value(s)

Build-specific block under `genome_builds["GRCh37"]` (or GRCh38 / T2T-CHM13v2.0):
- `contig` — NCBI accession (e.g. `"NC_000007.13"`)
- `strand` — `"+"` or `"-"`
- `cds_start`, `cds_end` — genomic coordinates of CDS (on contig)
- `exons` — array of 6-tuples (see below)
- `url` — source URL
- `source` — GFF column 2 text
- `tag` — e.g. `"MANE_Select"`, `"Ensembl_canonical"` (optional)
- `ccds` — CCDS identifier (optional)
- `transcript_support_level` — Ensembl TSL (optional)
- `note` — RefSeq note field (optional)
- `other_chroms` — for transcripts spanning multiple chromosomes (rare)

### Exon tuple

```
[alt_start_i, alt_end_i, exon_id, cds_start_i, cds_end_i, gap_or_null]
```
- Positions 0–1: genomic coordinates, 0-based half-open
- Position 2: exon order number
- Positions 3–4: transcript/cDNA coordinates (1-based, inclusive end)
- Position 5: GFF3-style gap string (e.g. `"M196 I1 M61"`) or `null`

### Gene object

```json
{
  "gene_symbol": "AOAH",
  "description": "acyloxyacyl hydrolase",
  "hgnc": "14472",
  "aliases": "...",
  "summary": "..."
}
```

### Schema versioning

`cdot/__init__.py`:
```python
def get_data_schema_int(version: str) -> int:
    major, minor, *_ = version.split(".")
    return 1000 * int(major) + int(minor)
```
Only major.minor compatibility is enforced; patch mismatches are allowed.
`JSONDataProvider.__init__` raises on incompatibility.

---

## Class Hierarchy: biocommons HGVS Providers

```
hgvs.dataproviders.interface.Interface
└── AbstractJSONDataProvider          (json_data_provider.py)
    ├── LocalDataProvider             (abstract, adds interval-tree queries)
    │   ├── JSONDataProvider          (in-memory JSON.gz)
    │   └── (no RedisDataProvider yet, but the abstraction exists for it)
    └── RESTDataProvider              (queries cdot.cc REST API)

EnsemblTarkDataProvider               (ensembl_tark_data_provider.py,
                                       independent hierarchy)
```

### AbstractJSONDataProvider — key internals

- **`NCBI_ALN_METHOD = "splign"`** — only alignment method; raises `HGVSError` if called with anything else.
- `assemblies` arg (default `["GRCh37", "GRCh38"]`) — built into `assembly_by_contig` dict via bioutils.
- **`get_tx_exons(tx_ac, alt_ac, alt_aln_method)`** — converts raw 6-tuple exons to the biocommons exon-dict format with `cigar` key; calls `_convert_gap_to_cigar()`.
- **`get_tx_info(tx_ac, alt_ac, alt_aln_method)`** — returns `{hgnc, cds_start_i, cds_end_i, tx_ac, alt_ac, alt_aln_method}`.
- **`get_tx_identity_info(tx_ac)`** — returns `lengths` array: `[cds_end + 1 - cds_start for each exon]` plus `{tx_ac, origin, hgnc}`.
- **`get_gene_info(gene)`** — returns `{hgnc, maploc, descr, summary, aliases, added}`; aliases formatted UTA-style: `{ALIAS1,ALIAS2}`.
- **`_convert_gap_to_cigar(gap)`** — GFF3 gap ops are inverted for HGVS CIGAR (I↔D, M→=).

### LocalDataProvider — region/gene queries

- **`get_tx_for_gene(gene)`** — sorts transcripts by total exon length descending.
- **`get_tx_for_region(alt_ac, alt_aln_method, start_i, end_i)`** — queries an `IntervalTree`; stores only transcript start/end (not full exon list) to save memory.
- **`_get_tx_by_gene_and_intervals()`** — static helper returning `(tx_by_gene: dict[str,set], tx_intervals: dict[str,IntervalTree])`.

### JSONDataProvider — file loading

- Accepts single path or list of paths; handles `.gz` transparently.
- Validates schema compatibility on init.
- `_tx_by_gene_and_intervals` is a `lazy` (deferred) property — built on first access.
- Stores `cdot_data_version` as a tuple.

### RESTDataProvider

- Default base URL: `https://cdot.cc` (or `http://cdot.cc` if `secure=False`).
- Caches transcripts and genes locally in dicts.
- `_get_from_url()` — checks `Content-Type` header, returns `None` on 404 (cached as missing).
- Endpoints: `/transcripts/gene/{gene}`, `/transcripts/region/{alt_ac}/{method}/{start}/{end}`, `/transcript/{tx_ac}`.

### EnsemblTarkDataProvider

- Base URL: `https://tark.ensembl.org/api`
- Pagination: follows `"next"` links via `_get_all_paginated_transcript_results()`.
- De-duplication: `_filter_dupes_take_most_recent()` — picks transcript with most recent `transcript_release_set` date.
- No alignment gaps in Tark data — CIGAR always `"="`.
- `_get_cds_start_end()` — calculates from `five_prime_utr_seq` / `three_prime_utr_seq` lengths (not stored in the main JSON format).
- `get_gene_info()` — **not implemented** (TARK lacks summary data).
- Has its own `EnsemblTarkSeqFetcher` built from a `PrefixSeqFetcher`:
  - `"NC_"` contigs → genome FASTA
  - `"ENST"` → TARK API
  - `"NM_"`, `"NR_"` → RefSeq, cross-validated with `VerifyMultipleSeqFetcher`

---

## Sequence Fetcher Hierarchy

```
seqfetcher.py
  PrefixSeqFetcher             — routes by accession prefix
  MultiSeqFetcher (abstract)
    ChainedSeqFetcher          — first fetcher to succeed wins
    VerifyMultipleSeqFetcher   — all fetchers must agree (used for validation)
  AlwaysFailSeqFetcher         — raises with a custom message
  AbstractTranscriptSeqFetcher — caches; subclasses implement _get_transcript_seq()

fasta_seqfetcher.py
  GenomeFastaSeqFetcher        — direct pysam FASTA access, contig→FastaFile map
  ExonsFromGenomeFastaSeqFetcher — reconstructs transcript from exons + CIGAR
    CIGAR 'D' → fills with "N" (deletion in ref)
    CIGAR 'I' → skips (insertion in ref)
    CIGAR '='/'X' → includes sequence
    Handles strand reversal
    Validates reconstructed length against expected
  FastaSeqFetcher              — PrefixSeqFetcher wrapper:
    default → ExonsFromGenomeFastaSeqFetcher
    "NC_" → GenomeFastaSeqFetcher
```

---

## PyHGVS Integration

`cdot/pyhgvs/pyhgvs_transcript.py`

```
AbstractPyHGVSTranscriptFactory
  PyHGVSTranscriptFactory       — in-memory transcripts dict
  JSONPyHGVSTranscriptFactory   — loads from JSON.gz or plain JSON
  RESTPyHGVSTranscriptFactory   — REST API, caches locally
                                   endpoint: /transcript/{tx_ac}
```

- `get_transcript(tx_ac, genome_build, sacgf_pyhgvs_fork=False)` — wraps `get_pyhgvs_data()`.
- `is_sacgf_pyhgvs_fork()` — checks installed pyhgvs version >= 0.12.0.

### PyHGVS data dict format

Standard PyHGVS:
```python
{
  "id": tx_ac,
  "chrom": contig,
  "start": first_exon_start,
  "end": last_exon_end,
  "strand": "+"/"-",
  "cds_start": cds_start,   # = tx end for non-coding
  "cds_end": cds_end,       # = tx end for non-coding
  "gene_name": gene_name,
  "exons": [[start, end], ...]
}
```

SACGF fork extras (pyhgvs >= 0.12.0):
- `cdna_match` — exons without 3rd element
- `start_codon_transcript_pos`
- `stop_codon_transcript_pos`
- `other_chroms` (optional)

---

## Data Generation Pipeline

### GFF/GTF Parsing (`gff_parser.py`)

**`GFFParser` (abstract base)**

Constructor params:
- `filename`, `genome_build`, `url`
- `discard_contigs_with_underscores=True` — filters alt-loci
- `no_contig_conversion=False` — skip bioutils name mapping
- `skip_missing_parents=False`

Key internal dicts:
```python
gene_data_by_accession: dict[str, dict]    # gene features
transcript_data_by_accession: dict[str, dict]  # transcript features
```

Processing stages:
1. `_parse()` — HTSeq.GFF_Reader iteration → `handle_feature()` (implemented by subclass)
2. `_add_transcript_data()` — stores exon/CDS; detects `cDNA_match` for gapped alignments
3. `_finish_process_features()` — converts raw features into final exon 6-tuples
4. `_finish()` — moves build-specific fields into `genome_builds` dict

Gap handling via `get_cdna_match_offset(gap_string, offset)`:
- Parses `"M185 I3 M250"` style strings
- M=match, I=ref gap (skip in ref), D=transcript gap
- Raises `ValueError` if `offset` falls inside a deletion

**`GTFParser(GFFParser)` — Ensembl GTF**
- FEATURE_ALLOW_LIST: `gene`, `transcript`, `exon`, `CDS`, `start_codon`, `stop_codon`
- Extracts gene/transcript version suffixes
- Gets protein version from CDS `protein_id` attribute

**`GFF3Parser(GFFParser)` — RefSeq GFF3**
- GFF3_GENES: `gene`, `pseudogene`, `ncRNA_gene`
- GFF3_TRANSCRIPTS_DATA: `exon`, `CDS`, `cDNA_match`, `five_prime_UTR`, `three_prime_UTR`
- Parses `Dbxref` for HGNC and CCDS IDs
- **Mitochondrial special case:** GRCh37 has only CDS (no mRNA), so fake transcript IDs are generated with prefix `"fake-"`.

### `cdot_json.py` — CLI subcommands

1. `gtf_to_json` / `gff3_to_json` — parse → enrich → write JSON.gz
2. `uta_to_json` — reads UTA CSV dumps; `_convert_uta_exons()` maps 0-based UTA to 1-based cdot; `_cigar_to_gap_and_length()` converts CIGAR to GFF3 gap (op flip: `=`→M, `D`→I, `I`→D, `X`→M).
3. `merge_historical` — merges multiple JSON files, latest transcript version wins; handles fake UTA gene accessions (prefix `"_"`) by mapping to real accessions from newer versions.
4. `combine_builds` — merges GRCh37 + GRCh38 + T2T-CHM13v2.0; checks codon consistency, removes stale older versions on mismatch; uses `ijson` for memory efficiency.
5. `release_notes` — prints release metadata.

Enrichment functions:
- `add_gene_info(gene_info_filename, genes)` — merges Entrez gene summaries (description, aliases, summary).
- `add_gencode_hgnc(gencode_hgnc_filename, genes, transcripts)` — adds HGNC IDs from GENCODE TSV.
- `add_canonical_transcripts(csv, genome_build, transcripts)` — tags Ensembl canonical transcripts.

### `cdot_gene_info.py` — Entrez fetcher

- Batches gene IDs (1000/batch, NCBI limit 10k); retries up to 3 times.
- Output JSON.gz:
  ```json
  {
    "cdot_version": "...",
    "api_retrieval_date": "...",
    "gene_info": {
      "GeneID": { "gene_symbol", "map_location", "description", "aliases", "summary" }
    }
  }
  ```

### `json_schema_version.py`

`JSON_SCHEMA_VERSION = "0.2.33"`

Changelog in comments:
- 0.2.29: Ensembl HGNC from outside GTFs
- 0.2.30: Ensembl GRCh37 canonical transcripts
- 0.2.31: `metadata` block with method/urls
- 0.2.32: `source` field (GTF column 2)
- 0.2.33: `ccds`, `transcript_support_level`

---

## Data Release Management (`data_release.py`)

- Queries GitHub API: `https://api.github.com/repos/SACGF/cdot/releases`
- `get_latest_data_version_and_release()` — finds release compatible with current schema (major.minor match).
- `get_latest_combo_file_urls(annotation_consortia, genome_builds)` — regex-parses release asset names:
  `r"cdot-(\d+\.\d+\.\d+)\.(refseq|ensembl)\.(.+)\.json\.gz"`
- Case-insensitive matching for consortium and build names.
- Supports: GRCh37, GRCh38, T2T-CHM13v2.0.

---

## Data Sources (`cdot_transcripts.yaml`)

| Source | Build | Coverage |
|--------|-------|----------|
| Ensembl GTF | GRCh37 | releases 82+ |
| Ensembl GTF | GRCh38 | releases 81–115 |
| Ensembl GTF | T2T-CHM13v2.0 | 2022-06 snapshot |
| RefSeq GFF3 | GRCh37 | multiple release dates |
| RefSeq GFF3 | GRCh38 | multiple release dates up to RS_2025_08 |
| RefSeq GFF3 | T2T-CHM13v2.0 | multiple releases |
| UTA | GRCh37 | uta_20241220 |
| UTA | GRCh38 | uta_20241220 |

Within a merge, newer source entries override older ones.

---

## Test Suite

### Test files

| File | What it tests |
|------|--------------|
| `test_gff_parsers.py` | GTF + GFF3 parsing, MT genes, contig conversion, tags |
| `test_json_data_provider_refseq.py` | RefSeq JSON provider, HGVS c_to_g, gene/region/protein queries |
| `test_json_data_provider_ensembl.py` | Ensembl JSON + Tark provider, both via same abstract test class |
| `test_pyhgvs.py` | PyHGVS transcript factory, coding + non-coding transcripts |
| `test_uta_conversion.py` | UTA exon/CIGAR conversion |
| `benchmark_hgvs.py` | Performance benchmarking (not a pytest test) |

### Key mock classes

- `MockSeqFetcher` — loads sequences from `test_data/transcript_sequences.json`
- `MockEnsemblTarkDataProvider` — overrides `_get_from_url()` to serve files from `test_data/ensembl_tark/`
- `MockGenomeTestFile` — pysam-compatible mock genome for PyHGVS tests

### Test data files

```
test_data/
  cdot.refseq.grch37.json          # RefSeq test transcripts
  cdot.ensembl.grch38.json         # Ensembl test transcripts
  transcript_sequences.json        # Nucleotide sequences for MockSeqFetcher
  ensembl_test.GRCh38.104.gtf      # Ensembl GTF fixture
  ensembl_test.GRCh38.111.gtf
  refseq_test.GRCh38.p13_genomic.109.20210514.gff
  refseq_test.GRCh38.p14_genomic.RS_2023_03.gff
  refseq_grch37_mt.gff             # MT-specific GRCh37
  refseq_grch38.p14_mt.gff         # MT-specific GRCh38
  hg19_chrY_300kb_genes.gtf        # UCSC format
  ensembl_tark/                    # Tark API JSON responses
```

---

## Alignment Gap Encoding

### GFF3 gap format (stored in JSON)

`"M196 I1 M61 I1 M181"`
- M = match (both transcript and genome advance)
- I = insertion in reference (transcript has extra bases; genome does not advance)
- D = deletion from reference (genome has extra bases; transcript does not advance)

### HGVS CIGAR format (computed on read)

Ops are **inverted** relative to GFF3 when returned to the HGVS library (because "to" and "from" sequences are swapped):
- GFF3 M → CIGAR `=`
- GFF3 I → CIGAR `D`
- GFF3 D → CIGAR `I`

### UTA CIGAR → GFF3 gap (in cdot_json.py)

During `uta_to_json`, CIGAR ops are translated:
- `=` → M
- `D` → I
- `I` → D
- `X` → M (mismatch treated as match)

Adjacent matches are merged before storage.

---

## Edge Cases and Special Handling

### Mitochondrial genes

- GRCh37 RefSeq: Only CDS features present (no mRNA records). GFF3Parser creates fake transcript IDs with prefix `"fake-"`.
- GRCh38 RefSeq: Both mRNA and CDS records exist; handled normally.
- MT contig: `NC_012920.1`

### Non-coding transcripts

- Detected by absence of coding features (CDS, start_codon, stop_codon).
- PyHGVS: `cds_start` and `cds_end` set equal to transcript end (zero-length CDS).

### Contig name normalization

- Input chromosome names ("1", "MT", "X") are mapped to NCBI accessions ("NC_000001.11", "NC_012920.1") via bioutils.
- `--no-contig-conversion` flag for non-standard or custom genomes.
- `discard_contigs_with_underscores=True` silently drops alt-loci (e.g. `"NT_"`, `"NW_"` contig names in input).

### HGNC ID sources

- **RefSeq GFF3:** Parsed from `Dbxref` attribute (may have `"HGNC:"` prefix that is stripped).
- **Ensembl GTF:** Parsed from `description` field via regex, or supplied via GENCODE external file.
- Stored without the `"HGNC:"` prefix.

### Gene accession edge cases

- RefSeq: numeric GeneID as `gene_version`.
- Ensembl: gene version number as `gene_version`.
- UTA: fake gene accessions prefixed with `"_"` — `merge_historical` attempts to remap these to real accessions from a newer source.

### `combine_builds` codon consistency check

When merging GRCh37 and GRCh38, if `start_codon`/`stop_codon` disagree between versions of the same transcript, the older (lower version) is dropped and a warning is logged.

### Ensembl TARK limitations

- No alignment gaps; CIGAR is always `"="`.
- `get_gene_info()` is not implemented (TARK API has no gene summary endpoint).
- Deduplication of duplicate transcript entries by most recent `transcript_release_set`.

---

## Known Unimplemented / Deferred Items

- `get_acs_for_protein_seq()` — always returns `[]`; comment in code says "TODO: drop"
- `EnsemblTarkDataProvider.get_gene_info()` — not implemented
- `EnsemblTarkDataProvider.get_tx_for_region()` — partially implemented in excerpt
- `RedisDataProvider` — abstraction exists (`LocalDataProvider`) but no Redis implementation is present

---

## Known Bugs / Correctness Gaps (from issues)

### CDS phase and ribosomal frameshifting (#76)

Some transcripts use −1 ribosomal frameshifting where one base is read twice, expressed in GFF3 as overlapping CDS exon coordinates, e.g. `NM_015068.3` (PEG10): `join(480..1436,1436..2605)`. cdot currently ignores this — it stores exon boundaries without the overlap, so the encoded protein cannot be correctly reconstructed from cdot data alone.

The GFF3 spec's column 8 `phase` field is the correct mechanism to track this. Dave's proposed fix: add an optional `phasing` array to the JSON exon data (omit if all zeros to avoid bloating unaffected transcripts). Blocked in practice until biocommons HGVS adds phase support. Other affected transcripts include `NM_001172437.2`, `NM_001184961.1`, `NM_004152.3`, `NM_001301020.1`, `NM_002537.3`, `NM_001301302.1`.

Note: GRCh38 GFF3 shows phase=0 for `NM_015068.3` even though the frameshift is well-established — the phasing may only appear in older GFF3 versions.

### FastaSeqFetcher sequence mismatch vs RefSeq (#55)

`NM_000399.3` (GRCh37 and GRCh38) has a CIGAR gap `M2190 D1 M282`. The RefSeq-stated sequence has a "T" at position 790, but the FASTA-reconstructed sequence has "A" at position 2697. Both sequences are the same total length. The RefSeq GFF3 flags the transcript with `gap_count=1` and a `Note` about a non-frameshifting indel vs the genomic sequence. Root cause unresolved — may be a RefSeq/genome divergence rather than a cdot bug.

### Ensembl GFF3 protein version parsing is broken (#101)

`GFF3Parser` does not correctly extract protein accession versions for Ensembl GFF3 files. Only `GTFParser` handles Ensembl protein IDs reliably. There are no explicit tests for the Ensembl GFF3 path, and the error message when it fails is not informative.

### `gene_name` is `None` for symbol-less transcripts (#89)

Some Ensembl transcripts (particularly those associated with BAC-derived or novel loci) have no gene symbol in the source GTF/GFF3. These are stored in cdot with `gene_name: None`. This is correct behaviour but can surprise users expecting a gene symbol. Transcripts like this exist because Ensembl stores `description=novel transcript` with no `Name` attribute on the gene.

---

## Transcript Source Priority

Within a merged JSON file, sources are processed in order from `cdot_transcripts.yaml`. The general rule is: **UTA first (baseline), then official annotation releases override in chronological order, with later releases taking precedence**. This means:

1. UTA (provides many historical transcript alignments)
2. Older RefSeq/Ensembl GFF/GTF releases (overrides UTA)
3. Newer RefSeq/Ensembl GFF/GTF releases (overrides older)

A consequence: the same transcript accession (e.g. `NM_001017995.2`) may have different `cds_start`/`cds_end` values depending on which source "won". The UTA 20241220 vs older RefSeq GFFs differ by 1 base for some transcripts — this is a genuine alignment difference between sources, not a cdot bug (#95). Needs documentation in the data files and/or a summary statistics file from the pipeline.

---

## Proposed JSON Schema Changes (#78)

If a breaking schema change is ever made (requiring a major version bump), the desired changes are:

1. **Exons as dicts instead of fixed-length arrays** — current 6-tuple format forces positional parsing; dicts would allow adding new fields (e.g. `phase`) without breaking existing consumers.
2. **Phasing support** — add per-exon phase as a new field (only possible cleanly once exons are dicts).
3. A suggestion to rename `cds_start`/`cds_end` to `tx_start`/`tx_end` was raised externally but rejected — those fields *are* coding coordinates; transcript extents are recoverable from `exons[0][0]` / `exons[-1][1]`.

---

## VEP-Compatible Releases (#63)

VEP (Variant Effect Predictor) uses specific annotation versions tied to Ensembl releases. It would be valuable to build and release VEP-compatible cdot JSON files from CI. The main constraint is memory; `orjson` was suggested as a more memory-efficient JSON library to support this. The VEP cache table at https://asia.ensembl.org/info/docs/tools/vep/script/vep_cache.html maps VEP versions to annotation sources and builds.

---

## Dependencies

### Installed with package (`setup.cfg`)

| Package | Role |
|---------|------|
| `requests` | REST API calls |
| `intervaltree` | Efficient region overlap queries |
| `more_itertools` | Utility iteration helpers |
| `bioutils>=0.5.8` | Assembly/contig name mapping |
| `lazy` | Deferred property evaluation |

### Required for data generation or optional features

| Package | Role |
|---------|------|
| `HTSeq` | GTF/GFF3 file parsing |
| `pysam` | FASTA file access |
| `hgvs` (biocommons) | HGVS library integration |
| `pyhgvs` | PyHGVS integration |
| `Bio` (BioPython) | Entrez API access |
| `ijson` | Streaming JSON for large combine_builds files |

---

## Project Context and Motivation (from docs/)

### Why cdot exists

- HGVS nomenclature (e.g. `NM_000492.3(CFTR):c.1438G>T`) is the de-facto standard for variant descriptions in clinical reports and publications.
- High-throughput sequencing works in genomic coordinates, so HGVS↔genomic conversion is a fundamental requirement.
- Transcripts are versioned (e.g. `NM_153253.3` = version 3); exon coordinates can change between versions, so historical HGVS resolution requires many transcript versions.

### Limitations of prior art that cdot addresses

| Tool | Problem |
|------|---------|
| **UTA (biocommons)** | Local install requires PostgreSQL; remote access connects to PostgreSQL (firewalls block this); no Ensembl transcript version support ([issue #233, June 2021](https://github.com/biocommons/uta/issues/233)); REST API planned since 2014 ([issue #164](https://github.com/biocommons/uta/issues/164)) but not delivered |
| **Counsyl pyHGVS** | Coordinate conversion bugs; no alignment gap support |

### Transcript count comparison

- UTA `uta_20210129`: 314,227 total transcripts; **141,262 with versioned accessions** (those with a `.` in the accession).
- VariantGrid (a related clinical platform): 788,252 distinct transcript versions.
- cdot target: **1.3 million+** transcript versions (exact figures TBD for paper).

### REST service URL

- Currently hosted at `http://cdot.cc/`
- New URL `cdotlib.org` has been obtained — migration planned.
- REST server code: https://github.com/SACGF/cdot_rest

### Related projects

- **pyreference** (https://github.com/SACGF/pyreference) — uses cdot JSON files for general bioinformatics analysis (not HGVS-specific). Demonstrates the broader utility of the JSON format.
- **gtf-to-genes** — discontinued predecessor that used a binary file format. cdot's gzipped JSON is faster and human-readable.
- **VariantGrid** (https://github.com/SACGF/variantgrid) — clinical genomics platform that contains HGVS utility code (in `genes/hgvs/`) that could potentially be extracted into a standalone library.

### Benchmark idea (from notes)

Resolve all HGVS entries in ClinVar and validate against the genomic coordinates in the VCF — a meaningful real-world accuracy benchmark.

### Key references

- HGVS nomenclature standard: *HGVS Recommendations for the Description of Sequence Variants: 2016 Update* — https://onlinelibrary.wiley.com/doi/pdf/10.1002/humu.22981
- biocommons hgvs library: Wang M et al. (2018), *Hum Mutat* 39(12):1803–1813, PMID 30129167

---

## Planned / In-Progress Work (from docs/todo.txt)

- **LRG transcripts** — currently handled in VariantGrid; not yet in cdot.
- **MANE transcript lookup by gene name** — e.g. resolving `GENENAME:c.XXX` by finding the MANE Select transcript. Could be implemented in the HGVS data provider by accepting a MANE file, or by searching transcript tags already stored in cdot JSON. Open question: does the genome build need to be passed in for this?
- **Snakemake run summary statistics** — the data generation pipeline (Snakemake) should emit a summary file (transcript counts, etc.) for use in the paper.
- **UTA in Snakemake** — UTA incorporation may already be done.
- **Combined build release files** — GitHub releases should include `grch37_grch38` and `grch37_grch38_t2t` combined files (used by the cdot REST service). Question open on whether the GRCh37+GRCh38 combined files should also include T2T.

---

## Key Algorithmic Notes

### Interval tree storage (memory optimization)

`LocalDataProvider` stores only `[tx_start, tx_end]` per transcript in the interval tree, not the full exon list. The full exon data lives in the main transcript dict and is fetched on demand.

### Sorted JSON output

`SortedSetEncoder` converts Python `set` objects to sorted lists before serialization. This ensures stable, diffable JSON output across runs.

### Lazy transcript-interval map

`JSONDataProvider._tx_by_gene_and_intervals` uses the `lazy` decorator so the interval tree is built only when a region or gene query first arrives, not at load time.

### Exon sort order

Exons are stored in **genomic order** (ascending `alt_start_i`), not transcript order. The `exon_id` field (position 2 in the tuple) preserves the original exon number for strand-aware reconstruction.
