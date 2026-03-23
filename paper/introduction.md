# Introduction

*Target: ~400 words. Three tight paragraphs: (1) problem scale, (2) existing tools and their gaps, (3) cdot's approach.*

---

HGVS nomenclature [den Dunnen 2016] is the international standard for describing sequence variants in clinical reports, databases, and publications. Converting HGVS descriptions to genomic coordinates requires authoritative, versioned transcript data — a transcript accession like `NM_000492.3` encodes not just a gene but a specific sequence version, and coordinate conversion is only possible with the matching alignment. The clinical scale of this requirement is large: ClinVar holds >3 million variants described in HGVS notation [Landrum 2025], the ACMG/AMP guidelines assume reliable conversion [Richards 2015], and Mutalyzer found ~50% error rates in submitted HGVS descriptions over five years [Lefter 2021], many attributable to missing transcript data. Transcript choice also has downstream consequences: only 44% of putative loss-of-function variants agreed between RefSeq and Ensembl annotation sets in ANNOVAR [McCarthy 2014], and tools without proper alignment gap support correctly annotate only 32% of BRCA2 variants [Münz 2015].

The two dominant Python HGVS libraries — biocommons/hgvs [Hart 2015; Wang 2018] and PyHGVS — both rely on external data providers. The Universal Transcript Archive (UTA), the standard backend for biocommons/hgvs, stores transcript alignments in PostgreSQL. This creates three barriers: a database installation requirement (or dependence on a remote server frequently blocked by clinical firewalls); limited coverage (~141,000 versioned alignments, RefSeq only); and slow remote access (~1 transcript/second). Ensembl transcripts, widely used in research genomics, are absent from UTA entirely. PyHGVS additionally lacks alignment gap support, causing incorrect coordinate conversion for transcripts with indels relative to the reference genome.

cdot (Complete Dict Of Transcripts) addresses these gaps. We generate compact JSON files from all available RefSeq GFF3 and Ensembl GTF annotation releases, covering [X] versioned transcript alignments across GRCh37, GRCh38, and — for the first time in an HGVS resource — T2T-CHM13v2.0 [Nurk 2022]. Files load locally in seconds and serve queries at [500–1000] transcripts/second; a REST API (cdotlib.org) is available for on-demand access. cdot plays the same role in the biocommons ecosystem for transcript coordinates as SeqRepo [Hart 2020] plays for biological sequences: local, high-performance, offline-capable access without database infrastructure.

---

*Keep this to 3 paragraphs. If it runs long, cut the McCarthy/Münz sentence — move that detail to implementation.*
