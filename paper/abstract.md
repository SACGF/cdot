# cdot: fast, versioned RefSeq and Ensembl transcript data for HGVS resolution

# Abstract

**Motivation:**

Clinical and research pipelines must resolve transcript-level HGVS variant descriptions
to genomic coordinates, and the practical goal is to resolve as many real-world
descriptions as possible, including the malformed strings and long-retired transcript
versions that accumulate in variant databases and clinical records. The standard
transcript source for the Python HGVS libraries, UTA, forces a tradeoff between a
locally installed PostgreSQL database and a slow public server that many clinical
networks firewall; it covers only ~{{ literature.uta_count | commas }} RefSeq
alignments, retains limited transcript history, and omits Ensembl.

**Results:**

cdot provides {{ coverage.total_count | commas }} versioned transcript/genome
alignments from RefSeq and Ensembl across GRCh37, GRCh38, and, for the first time for
the Python HGVS libraries, T2T-CHM13v2.0, as a single gzipped JSON file or a REST API
(cdotlib.org), with no database required. On ClinVar's submitted variant descriptions,
{{ clinvar_submitted.version_not_current_pct | dp(1) }}% of which cite a transcript
version that is no longer current, cdot resolved
{{ clinvar_submitted.cdot_resolved_pct | dp(1) }}% versus
{{ clinvar_submitted.uta_resolved_pct | dp(1) }}% for UTA. Loaded in memory, cdot
resolves ~{{ benchmark.cdot_local_tps | commas }} HGVS/second, about four times a
locally installed UTA, and compared remote-to-remote its REST API is roughly two orders
of magnitude faster than the public UTA server. A parser-independent cleaning step
(`clean_hgvs()`) repairs common formatting errors before resolution, and an opt-in
version-substitution step supplies a retired transcript version only when a
coordinate-safety check confirms it does not change the variant's coordinate.

**Availability and Implementation:**

https://github.com/SACGF/cdot; `pip install cdot`; MIT licence. Data files (JSON.gz) are
published on the GitHub releases page and archived at [Zenodo DOI]; a REST API is served
at cdotlib.org.

**Contact:** [email]

**Supplementary information:** Supplementary data are available at *Bioinformatics*
online.
