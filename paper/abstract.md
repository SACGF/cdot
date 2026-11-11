# cdot: fast, versioned RefSeq and Ensembl transcript data for HGVS resolution

# Abstract

**Motivation:**

cdot was built for the Australian Genomics Shariant project [@Tudini2022], which pools
variant classifications from clinical laboratories across the country, each recorded
against whatever transcript version the submitting lab used at the time. Making that
shared history usable meant resolving as many real-world transcript-level HGVS
descriptions as possible, including malformed strings and long-retired transcript
versions. The standard transcript source for the Python HGVS libraries, UTA, forces a
tradeoff between a locally installed PostgreSQL database and a slow public server that
many clinical networks firewall; it covered ~{{ literature.uta_count | commas }}
RefSeq alignments in its 2021 release, retains limited transcript history, and omits
Ensembl.

**Results:**

cdot provides {{ coverage.total_count | commas }} versioned transcript/genome
alignments from RefSeq and Ensembl across GRCh37, GRCh38, and, for the first time for
the Python HGVS libraries, T2T-CHM13v2.0, as a single gzipped JSON file or a REST API
(cdotlib.org), with no database required. {# lint-ok: 'no database required' describes the distribution format #}
On ClinVar's submitted variant descriptions,
{{ clinvar_submitted.version_not_current_pct | pct(1) }} of which cite a transcript
{% assert clinvar_submitted.cdot_resolved_pct > clinvar_submitted.uta_resolved_pct, "abstract: cdot X% 'versus' UTA Y% on submitted descriptions" %}
version that is no longer current, cdot reproduced ClinVar's genomic coordinate for
{{ clinvar_submitted.cdot_resolved_pct | pct(1) }} versus
{{ clinvar_submitted.uta_resolved_pct | pct(1) }} for UTA. Loaded in memory, cdot
resolves ~{{ benchmark.cdot_local_tps | commas }} HGVS/second,
{{ benchmark.cdot_local_vs_uta_local_ratio | approx }} times a locally installed UTA,
and its REST API is
{{ benchmark.cdot_rest_vs_uta_remote_ratio | orders_of_magnitude }} faster than the
public UTA server. A parser-independent cleaning step (`clean_hgvs()`) repairs
common formatting errors before resolution, and an opt-in version-substitution step
supplies a retired transcript version only when a coordinate-safety check indicates
the variant's coordinate is unchanged.

**Availability and Implementation:**

https://github.com/SACGF/cdot; `pip install cdot`; MIT licence. Data files (JSON.gz) are
published on the GitHub releases page and archived at [Zenodo DOI]; a REST API is served
at cdotlib.org.

**Contact:** [email]

**Supplementary information:** Supplementary data are available at *Bioinformatics*
online.
