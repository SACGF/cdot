# cdot: fast, versioned RefSeq and Ensembl transcript data for HGVS resolution

# Abstract

**Summary:**

HGVS nomenclature is the international standard for describing sequence variants in
clinical reports, public databases, and the research literature. A routine but
error-prone task is turning these descriptions back into genomic coordinates, which
depends on accurate, versioned transcript data. cdot supplies that data. It does not
resolve HGVS itself; it provides the versioned transcript-to-genome alignments
that the established Python HGVS libraries already use, so existing pipelines gain
coverage and speed without changing how they work.

The standard resource for this, UTA, requires a PostgreSQL database, has no Ensembl
support, and covers only ~{{ literature.uta_count | commas }} transcript alignments. cdot
provides {{ coverage.total_count | commas }} versioned transcript/genome alignments from
RefSeq and Ensembl across GRCh37, GRCh38, and T2T-CHM13v2.0, as a single gzipped JSON
file or served on demand over HTTP via a REST API (cdotlib.org), with no database
required.

cdot also helps with the malformed strings that clinical search boxes and importers
collect. A parser-independent cleaning step (`clean_hgvs()`) repairs common formatting
errors before resolution. When a cited transcript version has been retired, an opt-in
fallback substitutes the nearest available version, but only when a build-independent
check confirms the substitution does not move the variant. cdot integrates with both
major Python HGVS libraries (biocommons/hgvs and PyHGVS), stores MANE Select and Ensembl
canonical tags for gene-symbol lookup, and is the first transcript data source to bring
T2T-CHM13v2.0 to these libraries.

**Availability and Implementation:**

https://github.com/SACGF/cdot; `pip install cdot`; MIT licence. Data files (JSON.gz) at
cdotlib.org and [Zenodo DOI].

**Contact:** [email]

**Supplementary information:** Supplementary data are available at *Bioinformatics*
online.
