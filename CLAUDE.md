# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**cdot** converts RefSeq/Ensembl GTF/GFF3 annotation files into JSON format, providing high-performance transcript data for two Python HGVS variant annotation libraries: **biocommons HGVS** and **PyHGVS**.

Performance comparison: UTA public DB ~1s/transcript vs cdot JSON.gz ~500-1000 transcripts/second.

## Important Constraints

**Never run code against full datasets.** Processing full GTF/GFF3 files or complete JSON.gz data releases takes a very long time. Always use test data (in `tests/`) when verifying changes.

**PyHGVS is abandoned — prefer biocommons HGVS.** The `cdot/pyhgvs/` integration exists for legacy compatibility but PyHGVS is no longer maintained. Do not write new features that require significant PyHGVS-specific work. Focus new development on the biocommons HGVS path (`cdot/hgvs/dataproviders/`). If a feature is straightforward to support in both libraries, fine; if it requires real effort for PyHGVS, skip it and biocommons-only is acceptable.

**Keep the changelog up to date.** `CHANGELOG.md` is for **users of the client code** — only add entries for changes such a user would care about: client API/behaviour changes, dependency/compatibility changes (eg dropping a Python version), and changes to the **data content** they consume (note when a change only affects data and not client code). Add the entry to the `[unreleased]` section under `### Added`, `### Changed`, or `### Fixed`, and reference the GitHub issue it is done against (eg `#27`).

Do **NOT** add changelog entries for things a client-code user never sees: documentation, the paper (`paper/`), benchmark/analysis tooling (`analysis/`, `tests/benchmark_*.py`), CI/dev infrastructure, or the data-generation/build pipeline (`generate_transcript_data/`, Snakemake) when it doesn't change the released data. When in doubt, ask "would someone who only does `pip install cdot` and calls the library (or downloads a data release) notice this?" — if no, leave it out.

**The `cdot_private` repo holds real-world data — never copy examples from it into this repo.** It is checked out as a sibling directory (`../cdot_private`, github.com/SACGF/cdot_private) and contains a corpus of real search-bar HGVS strings (`combined_search_hgvs.csv`, built by its `process_search_hgvs.py`) plus `analyze_cleaning.py`, which runs that corpus through `cdot.hgvs.clean.clean_hgvs` to report rescue rates and failure patterns (issue #112). Use it to find what cleaning should handle, but **no example string sourced from `cdot_private` may ever appear in this repo's tests, comments, docstrings, or changelog** — synthesise equivalents from the standard public examples already used in the tests (eg `NM_000059.4` BRCA2, `NM_001754.5` RUNX1).

**Benchmarking.** `paper/scripts/benchmark_resolution.py` resolves real ClinVar (g.HGVS, c.HGVS) pairs (in `tests/test_data/clinvar_hgvs/`) through a pluggable provider (REST/JSON/UTA) to measure resolution accuracy, recovery (cleaning + version-bump), and speed; supports `--prefetch` and a local `--fasta`. This is how the README performance numbers are produced.

**Avoid AI writing tells in documentation and paper edits.** When writing or editing prose (docs, `README`, `paper/`, comments, commit/PR text), do not use em-dashes. Prefer plain punctuation: commas, periods, parentheses, or a rephrase. Also avoid other common LLM tells: "It's not just X, it's Y" constructions, hedging filler ("it's worth noting", "importantly", "delve"), over-use of bold for emphasis, and formulaic tricolon lists. Write plainly and directly, matching the voice of the surrounding text.

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

### REST API

The same team works on the REST server - http://github.com/SACGF/cdot_rest the REST API docs are here: https://cdotlib.org/static/api-docs.html 

## Git Commits
Commit only as the configured git user. Do NOT add `Co-Authored-By: Claude` trailers or any `Claude-Session` / Claude attribution lines to commit messages.

## GitHub Comments
When writing any comment on a GitHub issue or pull request, always preface it with 🤖 Written by Claude.
Do NOT close GitHub issues. Issues must go through a testing lifecycle before being closed by the user.
