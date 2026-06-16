# cdot documentation

Reference and how-to docs for [cdot](../README.md). These live in the repo, versioned with the code.

## Getting started

- **[Using local downloaded JSON.gz files](local_json_files.md)** — load release files into the HGVS libraries.
- **[Biocommons HGVS examples](examples_biocommons.md)** — c→g / g→c, plus a T2T-CHM13v2.0 example.
- **[PyHGVS examples](examples_pyhgvs.md)** — legacy PyHGVS integration (prefer biocommons).
- **[FastaSeqFetcher](fasta_seqfetcher.md)** — local FASTA sequence fetching (SeqRepo replacement).

## Advanced usage

- **[Advanced usage](advanced_usage.md)** — fixing messy HGVS input (`fix_hgvs` / `clean_hgvs`) and
  read-ahead batch retrieval for bulk processing (`RESTDataProvider.prefetch`).

## Data format reference

- **[JSON data format](json_data_format.md)** — every field in a cdot JSON(.gz) file, auto-generated
  from the typed models in [`cdot/models.py`](../cdot/models.py). Machine-readable
  [JSON Schema](cdot-json-schema.json) alongside it.
- **[Coordinates & exon alignments](coordinates_and_exons.md)** — how exon coordinates, exon IDs and
  the alignment `gap` (CIGAR-like) strings work, with worked examples.

## Data files & generation

- **[GitHub release file details](release_files.md)** — what each released `.json.gz` file contains.
- **[Create data from scratch](create_data_from_scratch.md)** — build the JSON files yourself from GTF/GFF3.

## Background

- **[cdot vs UTA](cdot_vs_uta.md)** — how cdot compares to the Universal Transcript Archive.
- **[Design notes & project direction](design_notes.md)** — why JSON, known issues, project goals.
