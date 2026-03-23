# Supplementary Material

*Supplementary data should be referenced in the abstract.*
*Submit as a single PDF labelled "Supplementary Data" at the journal.*

---

## Supplementary Tables

### Table S1 — RefSeq GFF3 annotation releases ingested by cdot

| Genome build | Annotation release | Date | Transcripts |
|-------------|-------------------|------|-------------|
| GRCh37 | [earliest release] | [date] | [N] |
| GRCh37 | ... | ... | ... |
| GRCh38 | [earliest release] | [date] | [N] |
| GRCh38 | ... | ... | ... |
| T2T-CHM13v2.0 | [release] | [date] | [N] |

*Source: Snakemake pipeline summary statistics (to be generated)*

### Table S2 — Ensembl GTF releases ingested by cdot

| Genome build | Ensembl release | Transcript count |
|-------------|-----------------|-----------------|
| GRCh37 | 75 | [N] |
| GRCh37 | ... | ... |
| GRCh38 | 76 | [N] |
| GRCh38 | ... | ... |
| T2T-CHM13v2.0 | [release] | [N] |

### Table S3 — cdot JSON schema fields (v0.2.33)

| Field | Level | Type | Description |
|-------|-------|------|-------------|
| `schema_version` | root | string | JSON schema version |
| `transcripts` | root | dict | keyed by transcript accession |
| `gene_name` | transcript | string | HGNC gene symbol |
| `gene_version` | transcript | int | |
| `hgnc` | transcript | string | HGNC ID |
| `biotype` | transcript | string | e.g. "protein_coding" |
| `genome_builds` | transcript | dict | keyed by build name (e.g. "GRCh38") |
| `contig` | build | string | chromosome/contig name |
| `strand` | build | int | +1 or -1 |
| `exons` | build | list | list of 6-tuples |
| `cds_start` | build | int | CDS start in genomic coords |
| `cds_end` | build | int | CDS end in genomic coords |
| `mane_select` | build | string | MANE Select accession (if applicable) |
| `mane_plus_clinical` | build | string | MANE Plus Clinical accession (if applicable) |
| `refseq_select` | build | bool | RefSeq Select flag |
| `ensembl_canonical` | build | bool | Ensembl canonical flag |

*Exon 6-tuple format: `[alt_start, alt_end, exon_id, cds_start, cds_end, gap]` where gap is a GFF3 gap string.*

### Table S4 — ClinVar HGVS benchmark details

*[Full breakdown of [N] ClinVar variants tested: counts by transcript source (RefSeq/Ensembl), genome build, failure reason (unknown transcript accession, unknown version, parse error). To be generated from benchmark script.]*

---

## Supplementary Figures

### Figure S1 — Transcript biotype distribution

Bar chart showing transcript counts by biotype (protein_coding, lincRNA, miRNA, pseudogene, etc.) across cdot sources and genome builds.

### Figure S2 — Historical RefSeq coverage by annotation release year

Cumulative count of unique NM_ accession versions added per RefSeq GFF3 annotation release. Demonstrates the value of ingesting multiple historical releases for resolving older HGVS strings.

### Figure S3 — T2T-CHM13v2.0 unique genes

Venn diagram or table showing genes with transcripts in T2T-CHM13v2.0 that have no equivalent (by gene name) in GRCh37 or GRCh38. Includes gene names, biotypes, and chromosomal locations (particularly acrocentric chromosome short arms).

### Figure S4 — cdot REST API response times

Distribution of response times for transcript lookups via cdotlib.org REST API, measured from [location] under [load conditions].

---

*Notes:*
- *Tables S1 and S2 require Snakemake pipeline to output summary statistics*
- *Table S3 should be verified against the actual JSON files before submission*
- *Table S4 and Figure S2 require the ClinVar benchmark script (#5)*
- *Reference supplementary material in the abstract ("Supplementary data are available at Bioinformatics online")*
