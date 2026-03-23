# Figures

*Aim for ≤5 main figures. Additional figures go to Supplementary.*
*Each legend must be self-contained — reader should understand the figure without the main text.*

---

## Figure 1 — Architecture overview

**Type**: Schematic diagram

**Content**: Show the cdot ecosystem:
- Left: data generation pipeline (RefSeq GFF3 + Ensembl GTF + UTA CSV → Snakemake pipeline → JSON.gz files)
- Centre: JSON.gz files and REST API (cdotlib.org)
- Right: two client paths: (a) biocommons/hgvs dataprovider → HGVS variant normalisation; (b) PyHGVS factory → HGVS conversion

Show SeqRepo alongside as the sequence complement.

**Legend**: `cdot architecture. Transcript coordinate data from RefSeq GFF3 and Ensembl GTF annotation releases is processed into compact JSON.gz files. These can be loaded locally via the Python client or queried on demand via the cdotlib.org REST API. Python integrations are provided for biocommons/hgvs and PyHGVS. Sequence data is supplied by SeqRepo or cdot's FastaSeqFetcher.`

---

## Figure 2 — Transcript coverage comparison

**Type**: Bar chart (grouped by genome build)

**Content**: Side-by-side bars for:
- cdot RefSeq (GRCh37, GRCh38, T2T)
- cdot Ensembl (GRCh37, GRCh38, T2T)
- UTA (GRCh37, GRCh38 only)
- TARK (GRCh38 only, Ensembl)

Y-axis: number of versioned transcript alignments (log scale or absolute)

**Legend**: `Transcript coverage by source and genome build. cdot covers [X] versioned transcript/genome alignments across RefSeq and Ensembl for GRCh37, GRCh38, and T2T-CHM13v2.0, compared to approximately 141,000 in UTA. The Universal Transcript Archive (UTA) does not include Ensembl transcripts or T2T assembly data. TARK (Ensembl Transcript Archive) provides Ensembl GRCh38 data only, as an online REST service.`

---

## Figure 3 — Speed benchmark

**Type**: Bar chart (log scale)

**Content**: Throughput (transcripts/second) for:
- cdot local JSON (GRCh38)
- cdot REST (cdotlib.org)
- UTA local (PostgreSQL, same machine)
- UTA remote (public server, typical conditions)

**Legend**: `Transcript retrieval throughput. cdot local JSON access achieves [500–1000] transcripts per second, compared to approximately 1 transcript per second via UTA remote access. All measurements performed on [hardware description]; UTA remote measurements represent median over [N] queries to the public UTA server. Error bars show [SD/IQR].`

---

## Figure 4 — ClinVar resolution benchmark

**Type**: Stacked bar chart or pie chart

**Content**: For [N] ClinVar HGVS variant descriptions:
- Resolved by cdot (GRCh38)
- Resolved by cdot (GRCh37, not GRCh38)
- Resolved by cdot (Ensembl only)
- Resolved by UTA only (not cdot)
- Unresolved by either

Compare cdot vs UTA side-by-side.

**Legend**: `Resolution of ClinVar HGVS variants. [N] HGVS variant descriptions from ClinVar [Landrum 2025] were tested against cdot and UTA. cdot resolves [X]% of variants, compared to [Y]% for UTA, driven by greater historical RefSeq version coverage and inclusion of Ensembl transcripts.`

---

## Figure 5 (optional) — MANE canonical transcript workflow

**Type**: Schematic / workflow diagram

**Content**: Show the gene-name HGVS lookup flow:
`Gene symbol ("BRCA2")` → `CanonicalTranscriptSelector` → `MANE Select transcript` → `biocommons/hgvs` → `HGVS c. notation`

Annotate with the MANE coverage statistic (>97% of protein-coding genes).

**Legend**: `Canonical transcript selection via MANE. cdot's CanonicalTranscriptSelector maps gene symbols to MANE Select (or MANE Plus Clinical) transcript accessions, enabling gene-name HGVS lookup. MANE Select covers >97% of protein-coding genes [Morales 2022].`

---

## Supplementary Figures

**Supp Fig S1** — Transcript biotype breakdown (protein-coding vs non-coding vs pseudogene) by source and build
**Supp Fig S2** — Cumulative transcript version coverage by RefSeq annotation release year (shows why historical releases matter)
**Supp Fig S3** — T2T unique transcripts: genes present in T2T but absent from GRCh37/38

---

*Notes:*
- *Figures 2, 3, 4 require benchmark data — see results.md and pre-submission checklist*
- *Figure 1 can be drafted now (architecture diagram)*
- *Figure 5 requires #36 (canonical transcript API) to be implemented*
- *All figure files should be submitted as separate PDFs or high-res TIFFs at revision*
