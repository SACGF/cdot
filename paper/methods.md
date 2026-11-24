# Methods

## Data sources and generation

cdot merges three transcript-annotation sources into a single dataset per build
(Figure 1A): the complete run of historical RefSeq GFF3 [@OLeary2016] and Ensembl GTF
[@Martin2023] releases (Supplementary Tables S1 and S2), plus UTA alignments that fill
gaps predating the oldest releases and are overridden by official annotation where
available. Releases are merged chronologically, newer entries overwriting older ones for
the same transcript version; this also backfills fields missing from early files. A
Snakemake pipeline parses each source and serialises to gzip-compressed JSON
(Supplementary Methods).

## JSON format

Each transcript entry stores shared metadata plus a `genome_builds` dict keyed by assembly
name, so one entry carries a transcript's GRCh37, GRCh38, and T2T-CHM13v2.0 alignments and
a single-build file is a drop-in GTF/GFF replacement. Per-build fields include a per-exon
`gap` string recording transcript-vs-genome indels and canonical tags (MANE Select and Plus
Clinical, RefSeq Select, Ensembl canonical); a root-level `schema_version` lets clients
reject incompatible files on load (Supplementary Table S3).

## String cleaning

`clean_hgvs()` is a pure string operation: it takes an HGVS string and returns the
repaired string with structured `HGVSFix` records describing every change. It needs no
genome build, sequence, parser, or data provider, so it can run as a pre-pass anywhere.
Cleaning is an ordered pipeline of single-purpose operations (Supplementary Table S7);
what it cannot fix is flagged, and it never raises. One repair needs transcript data: a
bare-number accession whose RefSeq prefix was dropped is restored only when exactly one
candidate exists in the loaded data.

## Version substitution

A separate, opt-in helper, `get_best_transcript_version()`, addresses transcript-version
drift. Substituting a different version is a heuristic that can be wrong, so it is never
applied automatically: when the requested version is absent, the helper returns the best
available adjacent version as a reported `HGVSFix`, and the caller decides whether to
accept it. A substitution is offered as safe only when the two versions' intrinsic CDS
structure (CDS length plus coding-exon segment lengths in transcript coordinates), which
is the same in every build, matches (Supplementary Methods). The check is validated
against ClinVar in Results R5.

## Access and client libraries

cdot implements the full `biocommons.hgvs.dataproviders.interface.Interface`, so it is a
drop-in replacement for UTA, through three providers (Figure 1B): `JSONDataProvider` (a
JSON.gz file loaded into memory), `RESTDataProvider` (the same data served at
cdotlib.org) and `EnsemblTarkDataProvider` (the Ensembl Transcript Archive).
`resolve_gene_hgvs()` maps a gene symbol to its MANE transcript
[@Morales2022; @Wright2023]. Sequence data comes from SeqRepo [@Hart2020] or cdot's
`FastaSeqFetcher`, which splices transcript sequence out of a local genome FASTA,
enabling fully offline operation (Supplementary Methods).

## Evaluation

Resolution accuracy, throughput, cleaning and version-substitution safety are measured
with scripts committed to the repository (`paper/scripts/`); protocols are in
Supplementary Methods. Resolution is scored against ClinVar's VCF coordinate, on
current-version ClinVar pairs and on the strings submitted to ClinVar. Throughput resolves
one committed ClinVar pair set through every backend with the engine and sequence layer
held constant. Cleaning is evaluated on a production query corpus and on clean ClinVar
strings with injected errors, run through `clean_hgvs()`, the LOVD HGVS syntax checker
[@LovdHgvsChecker], VariantValidator [@Freeman2018] and Mutalyzer [@Lefter2021]. Results
marked **[private data]** are aggregate statistics from private production corpora,
reported as frozen constants; every other result is reproducible from public data.
