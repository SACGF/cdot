# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**cdot** ("Complete Dict Of Transcripts") converts RefSeq/Ensembl GTF/GFF3 annotation files into JSON format, providing high-performance transcript data for two Python HGVS variant annotation libraries: **biocommons HGVS** and **PyHGVS**.

Performance comparison: UTA public DB ~1s/transcript vs cdot JSON.gz ~500-1000 transcripts/second.

## Important Constraints

**Never run code against full datasets.** Processing full GTF/GFF3 files or complete JSON.gz data releases takes a very long time. Always use test data (in `tests/`) when verifying changes.

**PyHGVS is abandoned — prefer biocommons HGVS.** The `cdot/pyhgvs/` integration exists for legacy compatibility but PyHGVS is no longer maintained. Do not write new features that require significant PyHGVS-specific work. Focus new development on the biocommons HGVS path (`cdot/hgvs/dataproviders/`). If a feature is straightforward to support in both libraries, fine; if it requires real effort for PyHGVS, skip it and biocommons-only is acceptable.

## Commands

```bash
# Install for development
pip install -e .

# Run all tests
python -m pytest tests/

# Run a single test file
python -m pytest tests/test_json_data_provider_refseq.py

# Run a single test method
python -m pytest tests/test_json_data_provider_refseq.py::TestJsonDataProvider::test_method_name

# Convert GTF/GFF3 to JSON (data generation - not part of installed package)
python generate_transcript_data/cdot_json.py gtf_or_gff [options]
```

## Architecture

### Package Layout

- **`cdot/hgvs/dataproviders/`** — Biocommons HGVS integration. Core class hierarchy: `AbstractJSONDataProvider` (implements `hgvs.dataproviders.interface.Interface`) → `LocalDataProvider` → `JSONDataProvider`. Also `EnsemblTarkDataProvider` for the Ensembl TARK REST API.
- **`cdot/pyhgvs/`** — PyHGVS integration. `AbstractPyHGVSTranscriptFactory` → `JSONPyHGVSTranscriptFactory` / `RESTPyHGVSTranscriptFactory`. Handles both standard PyHGVS and the SACGF fork.
- **`generate_transcript_data/`** — Data generation pipeline (not distributed as part of the pip package). Converts GTF/GFF3 files using HTSeq, enriches with gene info, and serializes to JSON.gz. Entry point: `cdot_json.py`.

### JSON Data Format

Each transcript contains `id`, `gene_name`, `protein`, `start_codon`, `stop_codon`, and a `genome_builds` dict mapping build names (e.g., `"GRCh38"`) to coordinates:
- `contig`: chromosome accession (e.g., `"NC_000007.14"`)
- `strand`: `"+"` or `"-"`
- `exons`: array of `[alt_start, alt_end, exon_id, cds_start, cds_end, gap_or_null]`

Schema versioning (`cdot/__init__.py`) uses major.minor; clients validate compatibility on load. Current schema: v0.2.

### Key Data Flows

**GTF/GFF3 → JSON** (`generate_transcript_data/`):
1. Parse with HTSeq (`gff_parser.py`: `GFF3Parser` or `GTFParser`)
2. Extract transcript/gene/exon features; compute CDS coords from start/stop codon features
3. Normalize contig names via bioutils assemblies
4. Convert CIGAR strings for alignment gaps
5. Optionally merge sources (`merge_historical`) or combine genome builds (`combine_builds`)
6. Serialize to JSON.gz using `SortedSetEncoder`

**JSON → HGVS queries** (`cdot/hgvs/dataproviders/json_data_provider.py`):
1. Load JSON.gz into memory; build interval trees for region queries and gene/contig→transcript mappings
2. Respond to biocommons HGVS interface methods: `get_tx_exons`, `get_tx_info`, `get_tx_for_region`, etc.
3. Sequence fetching: `ChainedSeqFetcher` tries multiple sources; `FastaSeqFetcher` for local FASTA files

### Data Sources

`generate_transcript_data/cdot_transcripts.yaml` defines all transcript sources: Ensembl GTF (releases 81–115, GRCh37/GRCh38), RefSeq GFF3, UTA, and T2T-CHM13v2.0. When sources are merged, newer entries override older ones.

### Release Management

`cdot/data_release.py` queries the GitHub API for the latest compatible data release (filtered by schema version), returning download URLs for RefSeq/Ensembl JSON.gz files across GRCh37, GRCh38, and T2T-CHM13v2.0.

## GitHub Comments
When writing any comment on a GitHub issue or pull request, always preface it with 🤖 Written by Claude.
Do NOT close GitHub issues. Issues must go through a testing lifecycle before being closed by the user.
