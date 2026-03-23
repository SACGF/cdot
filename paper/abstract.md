# Abstract

*Application Note abstract format: Summary / Availability and Implementation / Contact / Supplementary Information*
*Target: ~150 words. One paragraph per heading — no sub-bullets.*

---

**Summary:**

HGVS transcript notation is the clinical standard for variant description, underpinning >3 million ClinVar submissions and ACMG/AMP classification guidelines. Resolving HGVS to genomic coordinates requires versioned transcript data, yet the dominant resource (UTA) requires PostgreSQL infrastructure, lacks Ensembl support, and covers only ~141,000 transcript alignments. cdot (Complete Dict Of Transcripts) provides [X] versioned transcript/genome alignments from RefSeq and Ensembl across GRCh37, GRCh38, and T2T-CHM13v2.0 as compact JSON files loadable at [500–1000] transcripts/second — no database required. Python integrations are provided for both major HGVS libraries (biocommons/hgvs and PyHGVS). cdot stores MANE Select and Ensembl canonical tags enabling programmatic canonical transcript selection, and is the first HGVS resource to support the T2T-CHM13v2.0 assembly.

**Availability and Implementation:**

https://github.com/SACGF/cdot; `pip install cdot`; MIT licence. Data files (JSON.gz) at cdotlib.org and [Zenodo DOI].

**Contact:** [email]

**Supplementary information:** Supplementary data are available at *Bioinformatics* online.

---

*Fill in: [X] transcript count, speed numbers, Zenodo DOI, contact email.*
