# Abstract

*Application Note abstract format: Summary / Availability and Implementation / Contact / Supplementary Information*
*Target: ~150 words. One paragraph per heading — no sub-bullets.*

---

**Summary:**

HGVS transcript notation is the clinical standard for variant description, underpinning >{{ literature.clinvar_variants | commas }} ClinVar submissions and ACMG/AMP classification guidelines. Resolving HGVS to genomic coordinates requires versioned transcript data, yet the dominant resource (UTA) requires PostgreSQL infrastructure, lacks Ensembl support, and covers only ~{{ literature.uta_count | commas }} transcript alignments. cdot (Complete Dict Of Transcripts) provides {{ coverage.total_count | commas }} versioned transcript/genome alignments from RefSeq and Ensembl across GRCh37, GRCh38, and T2T-CHM13v2.0 as compact JSON files loadable at {{ benchmark.cdot_local_min_tps | commas }}–{{ benchmark.cdot_local_max_tps | commas }} transcripts/second — no database required. Python integrations are provided for both major HGVS libraries (biocommons/hgvs and PyHGVS). cdot stores MANE Select and Ensembl canonical tags enabling programmatic canonical transcript selection, and is the first HGVS resource to support the T2T-CHM13v2.0 assembly.

**Availability and Implementation:**

https://github.com/SACGF/cdot; `pip install cdot`; MIT licence. Data files (JSON.gz) at cdotlib.org and [Zenodo DOI].

**Contact:** [email]

**Supplementary information:** Supplementary data are available at *Bioinformatics* online.

---

*Fill in: Zenodo DOI, contact email. Run analysis scripts to populate coverage.csv and benchmark.csv.*
