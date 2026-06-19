# Introduction

*Target: ~400 words. Three tight paragraphs: (1) problem scale, (2) existing tools and
their gaps, (3) cdot's approach.*

---

HGVS nomenclature [@DenDunnen2016] is the international standard for describing sequence
variants in clinical reports, databases, and publications. Converting HGVS descriptions
to genomic coordinates requires authoritative, versioned transcript data. A transcript
accession such as `NM_000492.3` names a reference sequence at a specific version, and
that transcript is aligned to a particular genome build; the alignment determines where
the transcript's exons fall in genomic coordinates, so conversion is only possible with
the matching transcript-version alignment for the build in question. The scale is large:
ClinVar alone provides HGVS descriptions for >{{ literature.clinvar_variants | commas }}
variants [@Landrum2025]. Getting HGVS descriptions to resolve reliably is a known
problem, especially for the human-entered strings that reach search boxes and importers:
Mutalyzer, whose users type descriptions into a free-text box, found
~{{ literature.hgvs_error_rate_pct | dp(0) }}% error rates in
submitted HGVS descriptions over five years [@Lefter2021], many attributable to missing
transcript data. Transcript choice also has downstream consequences: only
{{ literature.lof_agreement_pct | dp(0) }}% of putative loss-of-function variants agreed
between RefSeq and Ensembl annotation sets in ANNOVAR [@McCarthy2014], underscoring the
value of covering both annotation sources.

The two dominant Python HGVS libraries are biocommons/hgvs [@Hart2015; @Wang2018] and
PyHGVS, and they obtain transcript data differently. biocommons/hgvs reads from a
pluggable data provider behind a common interface (a de facto protocol that additional
tools such as hgvs-weaver [@HgvsWeaver] now implement too, so a transcript backend can
target it once and serve multiple clients), whereas PyHGVS loads transcripts from its own
data structures, with example code for obtaining them from a UCSC RefGene export. The
standard backend for biocommons/hgvs, the Universal Transcript Archive (UTA), stores
transcript alignments in PostgreSQL, which is awkward to access either way it is reached.
Locally it requires installing and loading a PostgreSQL database. Remotely it depends on a
server reached over the PostgreSQL wire protocol, which many hospital and corporate
networks block while permitting only HTTP/HTTPS, and which is slow over the network
(~{{ literature.uta_remote_tps }} transcript/second). UTA is also limited in coverage
(~{{ literature.uta_count | commas }} versioned alignments, RefSeq only). Ensembl
transcripts,
widely used in research genomics and the predominant annotation source in Europe, are
absent from UTA entirely.

cdot was developed for the Australian Genomics Shariant project [@Tudini2022], a national
platform that pools historical variant classifications from clinical genetic-testing
laboratories. Those classifications span many years and are recorded in HGVS against
whichever transcript version each lab used at the time, so making them usable meant cdot
had to resolve as many of these real-world strings as possible, including ones that
reference transcript versions long retired from the current annotation.

To do this, cdot generates compact JSON files from all available RefSeq GFF3
and Ensembl GTF annotation releases, covering {{ coverage.total_count | commas }}
versioned transcript alignments across GRCh37 and GRCh38, and, for the first time for the
Python HGVS libraries, T2T-CHM13v2.0 [@Nurk2022]. JSON was chosen deliberately: every major
language parses it at high speed with built-in libraries, and it serialises cleanly over
HTTP, so the same files drive both the in-memory local provider and the REST API without
a separate data format. Files load locally in seconds and serve queries at
{{ benchmark.cdot_local_min_tps | commas }}–{{ benchmark.cdot_local_max_tps | commas }}
transcripts/second; a REST API (cdotlib.org) is available for on-demand access.

---

*Now 4 short paragraphs: (1) problem, (2) existing tools/gaps, (3) Shariant motivation +
focus thesis, (4) cdot's approach. If it runs long, the Shariant and approach paragraphs
can be merged. Münz/BRCA2 gap framing removed (gap correctness is a downstream library
concern, not a cdot result, per feedback).*
