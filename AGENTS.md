# AGENTS.md

This file provides guidance to AI coding agents (Claude Code, Codex, etc.) when working with code in this repository.

## Project Overview

**cdot** converts RefSeq/Ensembl GTF/GFF3 annotation files into JSON format, providing high-performance transcript data for two Python HGVS variant annotation libraries: **biocommons HGVS** and **PyHGVS**.

Performance comparison: UTA public DB ~1s/transcript vs cdot JSON.gz ~500-1000 transcripts/second.

## Important Constraints

**Never run code against full datasets.** Processing full GTF/GFF3 files or complete JSON.gz data releases takes a very long time. Always use test data (in `tests/`) when verifying changes. The test suite itself is not covered by this rule: it runs in a few seconds, so run it as often as you like (see [Test Suite](#test-suite)).

**PyHGVS is abandoned — prefer biocommons HGVS.** The `cdot/pyhgvs/` integration exists for legacy compatibility but PyHGVS is no longer maintained. Do not write new features that require significant PyHGVS-specific work. Focus new development on the biocommons HGVS path (`cdot/hgvs/dataproviders/`). If a feature is straightforward to support in both libraries, fine; if it requires real effort for PyHGVS, skip it and biocommons-only is acceptable.

**Keep the changelogs up to date.** There are two, because code and data are released separately, on independent version numbers:

- `CHANGELOG.md` is for users of the client code (`pip install cdot`), versioned by `cdot.__version__` and released to PyPI. Client API/behaviour changes, dependency/compatibility changes (eg dropping a Python version).
- `CHANGELOG-data.md` is for users of the published transcript data, versioned by `generate_transcript_data/json_schema_version.py` and released as a `data_v<version>` GitHub release. Changes to the data content they consume.

Put the entry in whichever one the change actually ships in; a change affecting both gets an entry in each. Add it to that file's `[unreleased]` section under `### Added`, `### Changed`, or `### Fixed`, and reference the GitHub issue it is done against (eg `#27`).

Do **NOT** add changelog entries for things a client-code user never sees: documentation, the paper (`paper/`), benchmark/analysis tooling (`paper/scripts/`, `tests/benchmark_*.py`), CI/dev infrastructure, or the data-generation/build pipeline (`generate_transcript_data/`, Snakemake) when it doesn't change the released data. When in doubt, ask "would someone who only does `pip install cdot` and calls the library (or downloads a data release) notice this?" — if no, leave it out.

**The `cdot_private` repo holds real-world data — never copy examples from it into this repo.** It is checked out as a sibling directory (`../cdot_private`, github.com/SACGF/cdot_private) and contains a corpus of real search-bar HGVS strings (`combined_search_hgvs.csv`, built by its `process_search_hgvs.py`) plus `analyze_cleaning.py`, which runs that corpus through `cdot.hgvs.clean.clean_hgvs` to report rescue rates and failure patterns (issue #112). Use it to find what cleaning should handle, but **no example string sourced from `cdot_private` may ever appear in this repo's tests, comments, docstrings, or changelog** — synthesise equivalents from the standard public examples already used in the tests (eg `NM_000059.4` BRCA2, `NM_001754.5` RUNX1).

**Benchmarking.** `paper/scripts/benchmark_resolution.py` resolves real ClinVar (g.HGVS, c.HGVS) pairs (in `tests/test_data/clinvar_hgvs/`) through a pluggable provider (REST/JSON/UTA) to measure resolution accuracy, recovery (cleaning + version-bump), and speed; supports `--prefetch` and a local `--fasta`. This is how the README performance numbers are produced.

**Sequence fetching: use `FastaSeqFetcher` (local genome FASTA) whenever Ensembl transcripts are involved.** Resolving a `c.`/`n.` HGVS needs the transcript sequence (biocommons `c_to_g` fetches it), served by SeqRepo, cdot's `FastaSeqFetcher` (splices the transcript's exon ranges out of a local genome FASTA), or `ChainedSeqFetcher` (tries several in order). SeqRepo is RefSeq-centric and does not hold every Ensembl (`ENST`) transcript sequence, so a SeqRepo-only run fails on those with `Failed to fetch ENST... from SeqRepo`, which surfaces as `no_data` and looks like a cdot coverage gap when it is really a sequence-layer gap. `FastaSeqFetcher` avoids this: it reproduces an Ensembl transcript exactly (Ensembl transcript sequence == genome), and it is also faster and offline (no remote NCBI fetch). Use it (or `ChainedSeqFetcher(SeqRepo, FastaSeqFetcher)`) for any Ensembl work or benchmark, e.g. `benchmark_resolution.py ... --fasta GENOME.fna.gz`. Caveats: (1) the FASTA must use the same contig names as cdot's data (NC_ accessions), so use the NCBI GRCh38 assembly (`GCF_000001405.39`, `NC_000001.11`...), not an Ensembl-named FASTA (`1`,`2`,...); (2) `FastaSeqFetcher` is not guaranteed to match a curated RefSeq transcript exactly (RefSeq sequence can differ from the genome at a few positions), so for an exact RefSeq sequence prefer SeqRepo or chain SeqRepo first. On this dev machine the GRCh38 FASTA is `/data/annotation/fasta/GCF_000001405.39_GRCh38.p13_genomic.fna.gz` (indexed).

**Avoid AI writing tells in documentation and paper edits.** When writing or editing prose (docs, `README`, `paper/`, comments, commit/PR text), do not use em-dashes. Prefer plain punctuation: commas, periods, parentheses, or a rephrase. Also avoid other common LLM tells: "It's not just X, it's Y" constructions, hedging filler ("it's worth noting", "importantly", "delve"), over-use of bold for emphasis, and formulaic tricolon lists. Write plainly and directly, matching the voice of the surrounding text.

**Docstrings, schema descriptions and comments describe the current state only.** Do not record how something used to behave, which version fixed it, or which issue changed it (eg "data before 0.2.37 had the gene biotype here", "previously this returned X"). That history belongs in the changelogs and the issue tracker. Notes about obsolete behaviour accumulate and quickly stop matching the code. The same goes when you change existing code: remove notes that are no longer true rather than adding a correction next to them.

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

# Regenerate docs/json_data_format.md and docs/cdot-json-schema.json after editing cdot/models.py
python generate_transcript_data/generate_json_docs.py
```

## Test Suite

`python -m pytest tests/` runs about 1300 tests in about 3 seconds, so run it freely and after every change. It needs the `test` extra (`pip install -e '.[test]'`: pytest, msgspec, pysam) and `pyhgvs` (in the `build` extra) for `tests/test_pyhgvs.py`. The handful of skips reported are the abstract base class `AbstractEnsemblTestCase` in `tests/test_json_data_provider_ensembl.py` skipping itself; its subclasses run the same tests against the JSON provider and the mocked TARK provider.

Almost everything in the suite runs offline:

- REST provider tests build `RESTDataProvider(url="https://example.org")` and either monkeypatch its `_get_from_url` / `_post_to_url` seams or hand it a `mock.Mock(spec=requests.Session)`; `tests/test_data_release.py` patches `requests.get`.
- Ensembl TARK tests run the real `EnsemblTarkDataProvider` code against recorded responses on disk (see `mock_ensembl_tark.py` below).
- Sequence fetching uses the mocks below, or a tiny FASTA written to a temp dir (`tests/test_fasta_seqfetcher.py`).

The exception is the `test_transcript` methods in `tests/test_json_data_provider_refseq.py` and `tests/test_json_data_provider_ensembl.py` (and any test built the same way). Their seqfetcher is `ChainedSeqFetcher(MockSeqFetcher, SeqFetcher())` and `c_to_g` with `replace_reference=True` fetches one genomic base, which the mock does not hold, so it falls through to biocommons `SeqFetcher`: a local SeqRepo when `HGVS_SEQREPO_DIR` is set, otherwise NCBI over the network. Everything else passes with `SeqFetcher.fetch_seq` patched to raise.

Helper modules in `tests/` (no `test_` prefix, not collected by pytest):

- `tests/genome.py`: `MockGenome` / `MockGenomeTestFile`, copied from counsyl pyhgvs. A pygr-style genome that answers `genome[chrom][start:end]` from a tab-separated lookup table of `(chrom, start, end, seq)` rows, `tests/test_data/grch37.genome`. Used by `test_pyhgvs.py`. Constructing it with `create_data=True` and a real FASTA `db_filename` records every query to the table, which is how the table is extended.
- `tests/mock_seqfetcher.py`: `MockSeqFetcher`, a biocommons-style seqfetcher that serves transcript sequences from a JSON `{accession: sequence}` file (`tests/test_data/transcript_sequences.json`, `translation_transcript_sequences.json`) and raises `HGVSDataNotAvailableError` for anything else. Used by the JSON provider and translation tests.
- `tests/mock_ensembl_tark.py`: `MockEnsemblTarkDataProvider`, a subclass of `EnsemblTarkDataProvider` whose `_get_from_url` maps the request URL to `tests/test_data/ensembl_tark/<path>/<query string>.json` and returns that recorded response. Add a new JSON file there to cover a new TARK query.
- `tests/benchmark_hgvs.py`: not a test. A CLI that resolves a file of (g.HGVS, c.HGVS) pairs through a UTA, REST, JSON or TARK provider and reports accuracy and timing, with notes at the end of the file on how the ClinVar input files were made. `paper/scripts/benchmark_resolution.py` is the current benchmark; this is the older one.

Test data lives in `tests/test_data/`: small GTF/GFF3 slices (single genes, the mitochondrial genome, a chr21 slice), small cdot JSON files, `clinvar_hgvs/` TSVs of (g., c.) pairs for the benchmarks, `ensembl_tark/` recorded responses, `version_safety/` real ClinVar records, and `gff_parser_golden/`, the expected JSON output for every GTF/GFF3 fixture. `tests/test_gff_parser_golden.py` compares parser output to the goldens byte for byte; after an intentional change to the output, regenerate them with `CDOT_UPDATE_GOLDEN=1 python -m pytest tests/test_gff_parser_golden.py` and review the diff.

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

The typed models are in `cdot/models.py`. `docs/json_data_format.md` and `docs/cdot-json-schema.json` are generated from them by `generate_transcript_data/generate_json_docs.py` (field descriptions come from the docstring under each field), so never edit those two files by hand. After touching `models.py`, run `python generate_transcript_data/generate_json_docs.py` and commit the result; a test checks that the generated files match the models.

### Key Data Flows

**GTF/GFF3 → JSON** (`generate_transcript_data/`):
1. Parse with HTSeq (`gff_parser.py`: `GFF3Parser` or `GTFParser` work out the gene/transcript/exon hierarchy for their file format)
2. Consortium-specific conventions (`annotation_consortium.py`: `RefSeq` / `Ensembl`, detected from the file or `--annotation-consortium`) supply HGNC, CCDS, protein versions, translation exceptions, RefSeq mito fakes
3. `transcript_builder.py` collects the features and computes exons, CDS coords and codon positions (coordinate maths in `transcript_coordinates.py`). Output must stay byte-identical to `tests/test_data/gff_parser_golden/` (`tests/test_gff_parser_golden.py`)
3. Normalize contig names via bioutils assemblies
4. Convert CIGAR strings for alignment gaps
5. Optionally merge sources (`merge_historical`) or combine genome builds (`combine_builds`)
6. Serialize to JSON.gz using `SortedSetEncoder`

**JSON → HGVS queries** (`cdot/hgvs/dataproviders/json_data_provider.py`):
1. Load JSON.gz into memory; build interval trees for region queries and gene/contig→transcript mappings
2. Respond to biocommons HGVS interface methods: `get_tx_exons`, `get_tx_info`, `get_tx_for_region`, etc.
3. Sequence fetching: `ChainedSeqFetcher` tries multiple sources; `FastaSeqFetcher` for local FASTA files

### Data Sources

`generate_transcript_data/cdot_transcripts.yaml` defines all transcript sources: Ensembl GTF (releases 81–116, GRCh37/GRCh38), RefSeq GFF3, UTA, and T2T-CHM13v2.0. When sources are merged, newer entries override older ones.

Adding a source touches more than the yaml. Checklist:

1. `generate_transcript_data/cdot_transcripts.yaml` (and `ensembl_gencode.csv` for a new Ensembl release).
2. The per-source release assets: the `files` array in the "Assemble release assets and notes" step of `.github/workflows/data_release.yml`, and the same list in `generate_transcript_data/github_release_upload.sh`. The two must match.
3. The file tables in `docs/release_files.md`, which list every asset a data release ships.
4. The example release number in `docs/create_data_from_scratch.md`.
5. `CHANGELOG-data.md`.

`docs/checking_for_new_sources.md` has the detail for steps 1, 2, 4 and 5 (which files are shipped individually, the GENCODE mapping, the parse check); `release_files.md` is the one it does not mention.

### Release Management

**Releases are done manually by the maintainer.** Do not tag, publish to PyPI, publish a GitHub release (including the data-release draft), or start the data-release workflow unless the user explicitly asks for that specific step in that message. Preparing a release (version bumps, changelog headings, building to check) is fine when asked; answering a question about release tooling is not permission to release.

To check for new upstream annotation releases (Ensembl, RefSeq, UTA, GENCODE) and add them, follow `docs/checking_for_new_sources.md`.

`cdot/data_release.py` queries the GitHub API for the latest compatible data release (filtered by schema version), returning download URLs for RefSeq/Ensembl JSON.gz files across GRCh37, GRCh38, and T2T-CHM13v2.0.

### REST API

The same team works on the REST server - http://github.com/SACGF/cdot_rest the REST API docs are here: https://cdotlib.org/static/api-docs.html 

## Branches, Commits and Pull Requests

Work on a branch off `origin/main` named `<topic>-<issue number>`, eg `ensembl-transcript-biotype-153`. Commit subjects are `#NNN - description`, with the issue number first (`#150, #152 - clean_hgvs: ...` when a change covers several). PR titles use the same `#NNN - description` form. The PR body starts with `🤖 By Claude <model>` on its own line, then `Fixes #NNN`, then what changed and how it was verified. The maintainer merges PRs (with a merge commit) and closes issues; do not merge, and do not close issues unless told to.

Commit only as the configured git user. Do NOT add `Co-Authored-By: Claude` trailers or any `Claude-Session` / Claude attribution lines to commit messages.

## GitHub Comments
When writing any comment on a GitHub issue or pull request, always preface it with 🤖 By Claude followed by the model that wrote it, eg `🤖 By Claude Opus 5.5`.
Do NOT close GitHub issues unless instructed by the user.
