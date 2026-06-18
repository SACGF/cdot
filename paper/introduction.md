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
the matching transcript-version alignment for the build in question. The clinical scale
of this requirement is large: ClinVar holds
>{{ literature.clinvar_variants | commas }} variants described in HGVS notation
[@Landrum2025], the ACMG/AMP guidelines assume reliable conversion [@Richards2015], and
Mutalyzer found ~{{ literature.hgvs_error_rate_pct | dp(0) }}% error rates in submitted
HGVS descriptions over five years [@Lefter2021], many attributable to missing transcript
data. Transcript choice also has downstream consequences: only
{{ literature.lof_agreement_pct | dp(0) }}% of putative loss-of-function variants agreed
between RefSeq and Ensembl annotation sets in ANNOVAR [@McCarthy2014], underscoring the
value of covering both annotation sources.

The two dominant Python HGVS libraries, biocommons/hgvs [@Hart2015; @Wang2018] and
PyHGVS, both rely on external data providers through a common data-provider interface;
biocommons/hgvs defines this interface, and additional tools such as hgvs-weaver
[@HgvsWeaver] now implement it too, making it a de facto protocol that a transcript
backend can target once and serve to multiple clients. The Universal Transcript Archive
(UTA), the standard backend for biocommons/hgvs, stores transcript alignments in
PostgreSQL. UTA imposes three practical barriers: a database installation requirement (or
dependence on a remote server frequently blocked by clinical firewalls); limited coverage
(~{{ literature.uta_count | commas }} versioned alignments, RefSeq only); and slow
remote access (~{{ literature.uta_remote_tps }} transcript/second). Ensembl transcripts,
widely used in research genomics and the predominant annotation source in Europe, are
absent from UTA entirely.

cdot was developed for the Australian Genomics Shariant project [@Tudini2022], a national
platform that pools historical variant classifications from clinical genetic-testing
laboratories. Those classifications span many years and are recorded in HGVS against
whichever transcript version each lab used at the time, so making them usable demanded a
single overriding goal that has shaped cdot ever since: **resolve as many real-world HGVS
strings as possible**, including those naming transcript versions long retired from the
current annotation.

cdot addresses these gaps. We generate compact JSON files from all available RefSeq GFF3
and Ensembl GTF annotation releases, covering {{ coverage.total_count | commas }}
versioned transcript alignments across GRCh37 and GRCh38, and, for the first time in an
HGVS resource, T2T-CHM13v2.0 [@Nurk2022]. JSON was chosen deliberately: every major
language parses it at high speed with built-in libraries, and it serialises cleanly over
HTTP, so the same files drive both the in-memory local provider and the REST API without
a separate data format. Files load locally in seconds and serve queries at
{{ benchmark.cdot_local_min_tps | commas }}–{{ benchmark.cdot_local_max_tps | commas }}
transcripts/second; a REST API (cdotlib.org) is available for on-demand access. cdot
provides the same local, offline access to transcript coordinates that SeqRepo
[@Hart2020] provides for biological sequences, without database infrastructure.

---

*Now 4 short paragraphs: (1) problem, (2) existing tools/gaps, (3) Shariant motivation +
focus thesis, (4) cdot's approach. If it runs long, the Shariant and approach paragraphs
can be merged. Münz/BRCA2 gap framing removed (gap correctness is a downstream library
concern, not a cdot result — feedback).*
