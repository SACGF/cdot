# Introduction

HGVS nomenclature [@DenDunnen2016] is the international standard for describing
sequence variants in clinical reports, databases, and publications. Converting a
transcript-level HGVS description (a `c.` or `n.` description) to genomic coordinates
requires authoritative, versioned transcript data: an accession such as `NM_000492.3`
names a reference sequence at a specific version, and its alignment to a genome build
determines where its exons fall, so conversion needs the matching transcript-version
alignment for the build in question. ClinVar alone provides HGVS descriptions for
>{{ literature.clinvar_variants | commas }} variants [@Landrum2025]. Human-entered
strings are often broken: of the ~26 million unique descriptions submitted to Mutalyzer
over five years of production logs, only
~{{ literature.mutalyzer_correct_pct | dp(0) }}% were correct as written,
~{{ literature.mutalyzer_error_pct | dp(0) }}% contained a syntactic or semantic error,
and Mutalyzer could automatically correct only
~{{ literature.mutalyzer_autocorrect_pct | dp(0) }}% (all rates per unique description)
[@Lefter2021].

The two dominant Python HGVS libraries are biocommons/hgvs [@Wang2018] and
PyHGVS. biocommons/hgvs reads from a pluggable data provider behind a common interface,
which other tools such as hgvs-weaver [@HgvsWeaver] also implement, so a single
transcript backend serves multiple clients; PyHGVS loads transcripts from its own data
structures. The standard backend for biocommons/hgvs is the Universal Transcript
Archive (UTA), which stores transcript alignments in PostgreSQL. Accessing it means
either installing that database locally or reaching a remote server over the PostgreSQL
wire protocol, which many hospital and corporate networks block and which is slow over
the network (~{{ literature.uta_remote_tps }} transcript/second). UTA is also limited
in coverage (~{{ literature.uta_count | commas }} versioned alignments, RefSeq only):
it retains only a limited set of historical transcript versions, and Ensembl
transcripts, widely used in research genomics and the predominant annotation source in
Europe, are absent entirely.

cdot was developed for the Australian Genomics Shariant project [@Tudini2022], a
national platform that pools variant classifications from clinical genetic-testing
laboratories. Classifications accumulate over many years, each recorded in HGVS against
whichever transcript version the lab used at the time, so making them usable meant
resolving as many of these real-world strings as possible, including ones that
reference transcript versions long retired from the current annotation.

To resolve as many of those descriptions as possible, cdot combines three things: broad
versioned transcript coverage, so the cited version is present in the first place;
parser-independent string cleaning, so a malformed description still parses; and an
opt-in, coordinate-safe version substitution for the cases where the exact cited version
is genuinely gone. For coverage, cdot generates compact JSON files from all available
RefSeq GFF3 and Ensembl GTF annotation releases, covering
{{ coverage.total_count | commas }} versioned transcript alignments across GRCh37,
GRCh38, and T2T-CHM13v2.0 [@Nurk2022]. The files load locally in seconds and can also be
served on demand through a REST API (cdotlib.org), so the same data backs both a fast
in-memory provider and lightweight remote access.
