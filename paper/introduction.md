# Introduction

HGVS nomenclature [@DenDunnen2016] is the international standard for describing sequence
variants in clinical reports, databases, and publications. Converting a transcript-level
HGVS description (a `c.` or `n.` description) to genomic coordinates requires
authoritative, versioned transcript data. A transcript
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
{{ literature.lof_agreement_pct | dp(0) }}% of putative loss-of-function variants were
classified as loss-of-function by both RefSeq and Ensembl annotation sets in ANNOVAR
[@McCarthy2014], underscoring the
value of covering both annotation sources.

The two dominant Python HGVS libraries are biocommons/hgvs [@Hart2015; @Wang2018] and
PyHGVS, and they obtain transcript data differently. biocommons/hgvs reads from a
pluggable data provider behind a common interface, whereas PyHGVS loads transcripts from
its own data structures, with example code for obtaining them from a UCSC RefGene export.
That data-provider interface has become a de facto standard: other tools such as
hgvs-weaver [@HgvsWeaver] implement it too, so a single transcript backend written to it
serves multiple clients at once, and any improvement to the backend benefits all of them. The
standard backend for biocommons/hgvs is the Universal Transcript Archive (UTA), which
stores transcript alignments in PostgreSQL. Accessing it means either installing and
loading that database locally or reaching a remote server over the PostgreSQL wire
protocol, which many hospital and corporate networks block and which is slow over the
network (~{{ literature.uta_remote_tps }} transcript/second). UTA is also limited in
coverage (~{{ literature.uta_count | commas }} versioned alignments, RefSeq only): it
retains only a limited set of historical transcript versions, and Ensembl transcripts,
widely used in research genomics and the predominant annotation source in Europe, are
absent entirely.

cdot was developed for the Australian Genomics Shariant project [@Tudini2022], a national
platform that pools variant classifications from clinical genetic-testing laboratories,
accumulated over many years and continually updated as laboratories submit new ones.
Those classifications span many years and are recorded in HGVS against
whichever transcript version each lab used at the time, so making them usable meant cdot
had to resolve as many of these real-world strings as possible, including ones that
reference transcript versions long retired from the current annotation.

To do this, cdot generates compact JSON files from all available RefSeq GFF3
and Ensembl GTF annotation releases, covering {{ coverage.total_count | commas }}
versioned transcript alignments across GRCh37, GRCh38, and, for the first time for the
Python HGVS libraries, T2T-CHM13v2.0 [@Nurk2022]. The files load locally in seconds and
can also be served on demand through a REST API (cdotlib.org), so the same data backs
both a fast in-memory provider and lightweight remote access.
