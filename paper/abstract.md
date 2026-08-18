# cdot: fast, versioned RefSeq and Ensembl transcript data for HGVS resolution

# Abstract

**Motivation:**

Resolving a transcript-level HGVS variant description to genomic coordinates requires
accurate, versioned transcript-to-genome alignments. The standard source for this,
UTA, requires a PostgreSQL database, covers only
~{{ literature.uta_count | commas }} alignments, retains limited transcript history,
and omits Ensembl, while the descriptions real pipelines receive are often malformed or
cite retired transcript versions.

**Results:**

cdot provides {{ coverage.total_count | commas }} versioned transcript/genome
alignments from RefSeq and Ensembl across GRCh37, GRCh38, and, for the first time for
the Python HGVS libraries, T2T-CHM13v2.0, as a single gzipped JSON file or a REST API
(cdotlib.org), with no database required. On ClinVar's submitted variant descriptions,
{{ clinvar_submitted.version_not_current_pct | dp(1) }}% of which cite a transcript
version that is no longer current, cdot resolved
{{ clinvar_submitted.cdot_resolved_pct | dp(1) }}% versus
{{ clinvar_submitted.uta_resolved_pct | dp(1) }}% for UTA. Loaded in memory, cdot
resolves ~{{ benchmark.cdot_local_tps | commas }} HGVS/second, nearly four orders of
magnitude above the public UTA server. A parser-independent cleaning step
(`clean_hgvs()`) repairs common formatting errors before resolution, and an opt-in
fallback substitutes a retired transcript version only when a coordinate-safety check
confirms the substitution does not move the variant.

**Availability and Implementation:**

https://github.com/SACGF/cdot; `pip install cdot`; MIT licence. Data files (JSON.gz) at
cdotlib.org and [Zenodo DOI].

**Contact:** [email]

**Supplementary information:** Supplementary data are available at *Bioinformatics*
online.
