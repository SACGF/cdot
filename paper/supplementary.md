# Supplementary Material

*Reference in abstract: "Supplementary data are available at Bioinformatics online."*
*Submit as a single PDF labelled "Supplementary Data".*

---

## Supplementary Tables

### Table S1 — RefSeq GFF3 annotation releases

| Genome build | Annotation release | Transcripts added |
|-------------|-------------------|------------------|
| GRCh37 | ... | ... |
| GRCh38 | ... | ... |
| T2T-CHM13v2.0 | ... | ... |

*Generated from Snakemake pipeline summary stats.*

### Table S2 — Ensembl GTF releases

| Genome build | Ensembl release | Transcript count |
|-------------|-----------------|-----------------|
| GRCh37 | 75–81 | ... |
| GRCh38 | 76–115 | ... |
| T2T-CHM13v2.0 | ... | ... |

### Table S3 — JSON schema fields (v0.2.33)

| Field | Level | Description |
|-------|-------|-------------|
| `schema_version` | root | Schema compatibility version |
| `transcripts` | root | Dict keyed by transcript accession |
| `gene_name` | transcript | HGNC gene symbol |
| `biotype` | transcript | e.g. `"protein_coding"` |
| `genome_builds` | transcript | Dict keyed by build name |
| `contig` | build | Chromosome/contig accession |
| `strand` | build | `+1` or `-1` |
| `exons` | build | List of `[alt_start, alt_end, exon_id, cds_start, cds_end, gap]` |
| `cds_start` / `cds_end` | build | CDS coordinates |
| `mane_select` | build | MANE Select accession if applicable |
| `mane_plus_clinical` | build | MANE Plus Clinical accession if applicable |
| `refseq_select` | build | RefSeq Select flag |
| `ensembl_canonical` | build | Ensembl canonical flag |

### Table S4 — ClinVar benchmark details

Breakdown of [N] ClinVar HGVS variants tested: counts by resolution source (RefSeq GRCh38, Ensembl GRCh38, GRCh37, unresolved), and failure reason (unknown accession, unknown version, parse error).

---

## Supplementary Figures

**Figure S1** — Speed benchmark: throughput (transcripts/second) for cdot local, cdot REST, UTA remote, UTA local. Log-scale bar chart.

**Figure S2** — Cumulative unique NM_ accession versions by RefSeq annotation release year. Motivates ingestion of multiple historical releases.

**Figure S3** — T2T-CHM13v2.0 unique transcripts: genes present in T2T but absent from GRCh37/GRCh38, annotated by chromosomal location.
