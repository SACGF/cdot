# Abstract

*Bioinformatics structured abstract — Motivation / Results / Availability / Contact*
*Target length: ~150–200 words*

---

**Motivation:**

HGVS transcript nomenclature is the clinical and research standard for describing sequence variants, underpinning ACMG/AMP variant classification guidelines and the >3 million variants archived in ClinVar. Converting HGVS descriptions to genomic coordinates requires comprehensive, versioned transcript data, yet existing resources require database infrastructure (PostgreSQL), lack Ensembl support, and provide limited historical transcript coverage — leaving tools unable to resolve a substantial fraction of clinical HGVS strings.

**Results:**

cdot (Complete Dict Of Transcripts) provides a compact JSON format for transcript coordinate data, generated from all available RefSeq GFF3 and Ensembl GTF annotation releases across GRCh37, GRCh38, and T2T-CHM13v2.0. cdot covers [X] versioned transcript/genome alignments — [9]× more than the Universal Transcript Archive — with no database infrastructure required. Data can be loaded locally ([500–1000] transcripts/second) or accessed via a REST API. Python integrations are provided for both major HGVS libraries (biocommons/hgvs and PyHGVS). cdot stores MANE Select, MANE Plus Clinical, and Ensembl canonical tags, enabling programmatic canonical transcript selection. On a benchmark of [N] ClinVar HGVS variants, cdot resolves [X]% versus [Y]% for UTA.

**Availability and Implementation:**

https://github.com/SACGF/cdot; `pip install cdot`; MIT licence. Data files available at [cdotlib.org / Zenodo DOI].

**Contact:** [dave.lawrence@sa.gov.au or similar]

**Supplementary information:** Supplementary data are available at *Bioinformatics* online.

---

*Notes for final draft:*
- *Fill in bracketed numbers once benchmarks are run (see results.md)*
- *Add Zenodo DOI once data files are archived*
- *Confirm corresponding author contact*
