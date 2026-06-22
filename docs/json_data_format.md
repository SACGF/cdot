# cdot JSON data format

> Auto-generated from the typed models in [`cdot/models.py`](../cdot/models.py) by `generate_transcript_data/generate_json_docs.py`. Do not edit by hand.

Generated from cdot **0.2.28**. A machine-readable [JSON Schema](cdot-json-schema.json) is generated alongside this file.

See [Coordinates & exon alignments](coordinates_and_exons.md) for a conceptual walk-through of exon coordinates, exon ordering and the alignment gap strings.

```python
from cdot import models

data = models.load("cdot-0.2.32.refseq.GRCh38.json.gz")
tx = data.transcripts["NM_001637.3"]
print(tx.gene_name, tx.protein)
for build_name, build in tx.genome_builds.items():
    for exon in build.exons:
        print(exon.alt_start, exon.alt_end, exon.cds_start, exon.gap)
```

## CdotData

The top-level contents of a cdot JSON(.gz) file.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `cdot_version` | string | yes |  |
| `genome_builds` | array of string | yes |  |
| `transcripts` | object ([Transcript](#transcript) values) | no | Keyed by transcript accession (e.g. `NM_001637.3`). |
| `genes` | object ([Gene](#gene) values) | no | Keyed by **gene ID** (e.g. RefSeq `80167`), *not* the symbol - the symbol is in each gene's `gene_symbol` field. For UTA-sourced data the gene ID is unknown, so a placeholder symbol is used as the key instead. |
| `metadata` | object (any values) or null | no | Release provenance (input_files, method, sys.argv, url_counts, ...); freeform, present in merged release files. |

## Transcript

A single transcript and its per-build coordinates.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `id` | string | yes | Transcript accession, e.g. `'NM_001637.3'`. |
| `genome_builds` | object ([GenomeBuild](#genomebuild) values) | yes | Build name (`'GRCh37'` / `'GRCh38'` / `'T2T-CHM13v2.0'`) -> coordinates. |
| `gene_name` | string or null | no |  |
| `gene_version` | string or null | no |  |
| `biotype` | array of string or null | no |  |
| `protein` | string or null | no | Protein accession (coding transcripts only). |
| `start_codon` | integer or null | no | 1-based transcript (cDNA) coordinate of the CDS start (coding only). |
| `stop_codon` | integer or null | no | 1-based transcript (cDNA) coordinate of the CDS end (coding only). |
| `hgnc` | string or null | no |  |
| `cdot` | string or null | no | cdot version that generated/last touched this transcript record. |
| `source` | array of string or null | no | Annotation source(s) this transcript came from (e.g. `['NCBI']`). |
| `partial` | integer or null | no | Non-zero if the transcript is annotated as partial/incomplete. |

## GenomeBuild

A transcript's coordinates on one genome build (e.g. `GRCh38`).

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `contig` | string | yes | Chromosome/contig accession, e.g. `'NC_000007.14'`. |
| `strand` | string | yes | `'+'` or `'-'`. |
| `exons` | array of [Exon](#exon) | yes |  |
| `url` | string or null | no | Source annotation file the transcript was extracted from. |
| `cds_start` | integer or null | no | 0-based genomic CDS start on the contig (coding transcripts only). |
| `cds_end` | integer or null | no | 0-based genomic CDS end on the contig (coding transcripts only). |
| `start` | integer or null | no | 0-based genomic start of the transcript on the contig. |
| `stop` | integer or null | no | 0-based genomic end of the transcript on the contig. |
| `tag` | string or null | no | Comma-separated tags (e.g. `'MANE_Select,Ensembl_canonical'`), verbatim from the source GTF/GFF. Spelling differs by consortium: RefSeq uses spaces (`'MANE Select'`, `'RefSeq Select'`) while Ensembl uses underscores (`'MANE_Select'`, `'Ensembl_canonical'`), so a raw-JSON consumer must handle both. To rank/compare tags across sources use `cdot.hgvs.gene_hgvs`, which normalises spelling before comparison. |
| `note` | string or null | no |  |
| `other_chroms` | array of string or null | no | Other contigs this transcript also aligns to (e.g. PAR/alt loci). |
| `source` | string or array of string or null | no | Annotation source (GTF/GFF column 2, e.g. `'BestRefSeq'`); data schema >= 0.2.32. A single string in early 0.2.32 data, a list (e.g. `['BestRefSeq']`) from 0.2.33 on. |
| `ccds` | string or null | no | CCDS id, when present; data schema >= 0.2.33. |
| `transcript_support_level` | string or null | no | Ensembl transcript support level (TSL); data schema >= 0.2.33. |

## Exon

A single exon alignment.

Serialised in JSON as a positional array, e.g.:

    [36552548, 36552986, 20, 2001, 2440, "M196 I1 M61 I1 M181"]

Coordinates: `alt_start`/`alt_end` are 0-based genomic coordinates on the
contig; `cds_start`/`cds_end` are 1-based transcript (cDNA) coordinates.
`exon_id` is the ordinal of the exon in stranded (transcript) order.

The original positional access (`exon[0]`, tuple unpacking) still works, so
this is a drop-in for the plain list it replaces.

Encoded as a JSON array (positional):

| # | Field | Type | Description |
|---|-------|------|-------------|
| 0 | `alt_start` | integer | 0-based genomic start of the exon on the contig. |
| 1 | `alt_end` | integer | 0-based genomic end of the exon on the contig (exclusive). |
| 2 | `exon_id` | integer | Exon ordinal in stranded (transcript) order, starting at 0. |
| 3 | `cds_start` | integer | 1-based transcript (cDNA) start coordinate of the exon. |
| 4 | `cds_end` | integer | 1-based transcript (cDNA) end coordinate of the exon. |
| 5 | `gap` | string or null | Alignment gap as a cdot 'gap' string (e.g. `'M196 I1 M61'`) or `None` if the exon aligns cleanly. |

## Gene

Gene-level metadata (present when the source provided gene info).

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `gene_symbol` | string or null | no |  |
| `aliases` | string or null | no |  |
| `biotype` | string or array of string or null | no |  |
| `description` | string or null | no |  |
| `hgnc` | string or null | no |  |
| `map_location` | string or null | no |  |
| `summary` | string or null | no |  |
| `url` | string or null | no |  |
| `source` | array of string or null | no |  |
| `transcripts` | array of string or null | no | Transcript accessions belonging to this gene (when provided). |

## Source-specific notes

* **Gene map keys.** The top-level `genes` map is keyed by gene ID (e.g. RefSeq `80167`),
  not by gene symbol. The human-readable symbol is stored in each record's `gene_symbol` field.
* **UTA-derived data.** When data is built from UTA (a `url` like
  `postgresql://uta.biocommons.org/uta_20210129`), the gene ID is not available - only the symbol.
  In that case a placeholder symbol is used as the `genes` map key.
* **Coordinate systems.** Genomic coordinates (`alt_start`/`alt_end`, build `cds_start`/`cds_end`,
  `start`/`stop`) are 0-based. Transcript (cDNA) coordinates inside each exon
  (`cds_start`/`cds_end`) are 1-based.
* **Tags are verbatim.** The build `tag` field is passed through unchanged from the source
  GTF/GFF, so MANE/canonical spelling differs by consortium: RefSeq uses spaces (`MANE Select`,
  `RefSeq Select`) while Ensembl uses underscores (`MANE_Select`, `Ensembl_canonical`). A consumer
  reading the raw JSON must handle both spellings. To rank or compare tags across sources, use
  `cdot.hgvs.gene_hgvs`, which normalises spelling (spaces to underscores) before comparison.
