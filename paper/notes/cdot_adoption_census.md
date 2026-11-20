# Census of public projects that use cdot

Backs `empirical_results/adoption.csv`. Snapshot of public GitHub, 2026-11-20.

## Rule

A project counts if its own code (not only a dependency file) imports the `cdot` Python
package or loads cdot JSON, and it is not by SACGF (the cdot authors). Each repository
below was opened and the cited file read. A project's data pipeline repository counts
with the project it feeds.

This only sees public GitHub repositories, so private, institutional and non-GitHub use is
missing. The count is a lower bound.

## Method

- GitHub code search (`gh api search/code`): `cdot.hgvs.dataproviders language:Python`,
  `"from cdot" language:Python`, `"import cdot" language:Python`, `cdot language:Rust`, and
  `cdot` in `requirements.txt`, `pyproject.toml`, `setup.py`, `environment.yml` and
  `Cargo.toml`. "cdot" is a common substring, so most hits are unrelated and were
  discarded on reading.
- GitHub's dependents page for SACGF/cdot (repositories and packages).
- deps.dev reverse dependencies for the PyPI package (none recorded).
- The survey of biocommons/hgvs data-provider users in biocommons/hgvs#958, each entry
  re-checked here.

## Counted (12)

Stars and last push as of 2026-11-20.

### Python, using cdot's data providers (9)

| Project | Stars | Evidence |
|---|---|---|
| VariantEffect/mavedb-api (MaveDB) | 19 | `src/mavedb/data_providers/services.py`: `StrictRESTDataProvider` subclasses cdot `RESTDataProvider`, with a `ChainedSeqFetcher` of FASTA fetchers |
| phiweger/faltwerk | 9 | `scripts/map_clinvar_to_proteins.py` |
| phiweger/primer4 | 4 | `primer4/stream.py` |
| KarchinLab/open-cravat-extras | 1 | `hgvs-api/hgvs_api/service/coordinates.py` |
| bihealth/dotty | 0 | `dotty/core.py`: one `JSONDataProvider` per assembly |
| VCCRI/heartvar-open | 0 | `backend/clients/cdot_sqlite.py`: subclasses `LocalDataProvider` over cdot JSON loaded into SQLite |
| MarcoMernberger/mutlibrary | 0 | `src/mutlibrary/models/annotator.py` |
| biocommons/uta-clients | 0 | `src/uta_clients/cdot.py` |
| varfish-org/mehari-annotation-comparison | 0 | `hgvs-comparison/workflow/scripts/check_hgvs.py` |

### Rust, reading cdot JSON directly (3)

| Project | Stars | Evidence |
|---|---|---|
| varfish-org/mehari | 33 | `mehari/src/db/transcripts/create/cdot.rs` builds its transcript database from cdot JSON; its data pipeline varfish-org/mehari-data-tx (`workflow/rules/cdot.smk`) counts with it |
| fulcrumgenomics/ferro-hgvs | 32 | `src/data/mod.rs`, `src/data/mapping.rs` |
| varfish-org/hgvs-rs | 18 | `src/data/cdot/json.rs`, a cdot JSON data provider for the Rust port of biocommons/hgvs |

## Not counted

- SACGF's own projects: SACGF/variantgrid, SACGF/cdot_rest, SACGF/hgvs_shim,
  SACGF/too-many-transcripts, SACGF/clinvar-hidden-structures.
- Declared but not used in the project's code: VariantEffect/mavetools (`pyproject.toml`),
  ave-dcd/dcd_mapping (`notebooks/analysis/requirements.txt` only), cancervariants/fusor
  (transitive, listed on the dependents page).
- yxuil/hgvs-ext: declares `cdot>=0.2.26` but imports a module path that does not exist
  in cdot (`cdot.hgvs.dataproviders.uta`).
- populationgenomics/hgvs-weaver: implements the biocommons/hgvs data-provider interface
  over its own data; cdot appears only in a benchmark container script that prepares
  ferro-hgvs inputs.
