# Introduction

*Target: ~800–1000 words. Sets up the problem, existing solutions and their gaps, and cdot's approach.*

---

## 1.1 HGVS nomenclature as the clinical standard

The Human Genome Variation Society (HGVS) nomenclature [den Dunnen 2016] is the established international standard for describing sequence variants in clinical reports, research publications, and variant databases. HGVS strings (e.g., `NM_000492.3(CFTR):c.1521_1523delCTT`) express variants in terms of a named, versioned transcript and a position within that transcript's coordinate system. This representation is human-interpretable and clinically precise, but creates a fundamental computational challenge: converting between HGVS transcript coordinates and the genomic coordinates used by sequencing pipelines requires authoritative, versioned transcript data.

The scale of demand is substantial. ClinVar, the primary archive of clinically interpreted variants, now holds >3 million variants submitted by >2,800 organisations [Landrum 2025], virtually all described in HGVS transcript notation. The ACMG/AMP variant interpretation guidelines [Richards 2015], now the de-facto standard for clinical variant classification worldwide, assume reliable HGVS coordinate conversion as a prerequisite. Mutalyzer, an HGVS syntax checker, processed over 133 million HGVS descriptions in five years, finding approximately 50% contained errors [Lefter 2021] — a rate consistent with the original Mutalyzer study [Wildeman 2008]. Many of these errors trace to missing or misidentified transcript data.

## 1.2 The annotation disagreement problem

Transcript choice has measurable consequences for annotation outcomes. McCarthy et al. [2014] found that only 44% of putative loss-of-function variants agreed between RefSeq and Ensembl annotation sets when processed by ANNOVAR, and only 65% agreed between ANNOVAR and VEP. Park et al. [2022] found similar ~85% agreement between ANNOVAR and SnpEff on clinical variants. Most dramatically, Münz et al. [2015] showed that most variant annotation tools correctly annotated only 32% of BRCA2 variants — the remainder failing due to improper handling of alignment gaps between transcript and reference genome sequences, or incorrect transcript strand orientation.

These numbers establish that transcript data quality and coverage are not minor implementation details but have major consequences for clinical and research conclusions.

## 1.3 Canonical transcript selection

A longstanding problem is the lack of a single "standard" transcript per gene. RefSeq and Ensembl have historically annotated different representative transcripts, causing inconsistency in clinical reports. The MANE (Matched Annotation from NCBI and EMBL-EBI) initiative [Morales 2022] addresses this by jointly selecting one representative transcript per protein-coding gene (MANE Select) and additional clinically relevant transcripts (MANE Plus Clinical), now covering >97% of protein-coding genes. Wright et al. [2023] argue that clinical laboratories should standardise on MANE transcripts; tools must therefore support programmatic canonical transcript lookup.

## 1.4 Existing infrastructure and its limitations

The two dominant Python HGVS libraries — biocommons/hgvs [Hart 2015; Wang 2018] and PyHGVS (Counsyl) — both rely on external data providers for transcript information. The Universal Transcript Archive (UTA), the standard backend for biocommons/hgvs, stores transcript alignments in PostgreSQL. This creates three practical barriers: (i) UTA requires either a local PostgreSQL installation or a remote connection, which is frequently blocked by firewalls in clinical environments; (ii) UTA provides approximately [141,000] versioned transcript alignments, with no Ensembl support; and (iii) remote UTA access is slow (~1 transcript/second), insufficient for production-scale HGVS processing. PyHGVS additionally lacks alignment gap support, causing incorrect coordinate conversion for a meaningful fraction of transcripts.

Reference genome evolution adds further complexity. The T2T-CHM13v2.0 assembly [Nurk 2022] adds ~200 million base pairs not present in GRCh38, including regions containing real genes. The Human Pangenome Reference [Liao 2023] extends this further with 47 phased diploid assemblies. Existing HGVS tools are limited to GRCh37 and GRCh38, leaving clinical and research variants against newer assemblies unresolvable.

## 1.5 cdot

Here we describe cdot (Complete Dict Of Transcripts), a Python package and data resource that addresses these limitations. cdot converts RefSeq GFF3 and Ensembl GTF annotation files — spanning multiple historical releases — into compact JSON files covering [X] versioned transcript alignments across GRCh37, GRCh38, and T2T-CHM13v2.0. These files can be loaded locally at [500–1000] transcripts/second or queried through a REST API (cdotlib.org), with no database infrastructure required. cdot implements the biocommons/hgvs data provider interface and provides a PyHGVS transcript factory, integrating with both major Python HGVS libraries. It stores MANE Select, MANE Plus Clinical, RefSeq Select, and Ensembl canonical tags, enabling programmatic canonical transcript selection.

cdot plays the same role in the biocommons ecosystem for transcript coordinates as SeqRepo [Hart 2020] plays for biological sequences: a local, high-performance, offline-capable alternative to remote database access, achieving [500–1000]× speedup over UTA remote access (compared to SeqRepo's 1300× speedup over remote sequence retrieval).

---

*Notes for final draft:*
- *Fill in bracketed numbers once benchmarks are run*
- *Tighten to ~800 words; trim wherever possible*
- *Ensure all citations match the references.md list*
