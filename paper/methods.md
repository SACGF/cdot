# Methods

## Data sources and generation

cdot merges three transcript-annotation sources into a single dataset per build
(Figure 1A), downloading the complete run of historical releases from the RefSeq and
Ensembl FTP sites because an NM_ version retired years ago is still cited in patient
records and ClinVar submissions:

1. **RefSeq GFF3**: {{ sources.refseq_grch37_releases | int }} GRCh37,
   {{ sources.refseq_grch38_releases | int }} GRCh38, and
   {{ sources.refseq_t2t_releases | int }} T2T-CHM13v2.0 NCBI annotation releases
   through the dated RS_2025_08 release (Supplementary Table S1), plus NCBI's
   historical transcript-alignment file (RS_2024_08), covering GRCh38 versions
   replaced or suppressed before any archived release.

2. **Ensembl GTF**: {{ sources.ensembl_grch37_releases | int }} GRCh37,
   {{ sources.ensembl_grch38_releases | int }} GRCh38, and
   {{ sources.ensembl_t2t_releases | int }} T2T-CHM13v2.0 releases, spanning Ensembl
   releases 81 to 116 (Supplementary Table S2); GTF is used because Ensembl's GFF3
   files omit transcript protein versions.

3. **UTA**: the Universal Transcript Archive computes its own alignments, so it holds
   some accessions never published in any annotation file; these fill gaps predating
   the oldest releases and are overridden by official annotation where available.

Releases are merged chronologically, newer entries overwriting older ones for the same
transcript version; this also backfills fields missing from early files, such as the
alignment-gap CIGAR strings absent from early RefSeq GFF3. A Snakemake pipeline parses
each source with HTSeq, normalises contig names via bioutils [@bioutils], extracts
coding-sequence (CDS) boundaries from start/stop codon features, attaches gene metadata
(symbols, HGNC IDs, biotypes) from HGNC and NCBI gene-info files, and serialises to
gzip-compressed JSON.

## JSON format

JSON parses quickly in every major language and serialises cleanly over HTTP, so the
same files drive both the in-memory provider and the REST API. Each transcript entry stores shared metadata (`gene_name`, `hgnc`, `biotype`,
accession and version) plus a `genome_builds` dict keyed by assembly name. The build is
always a dict key, so a single-build file is a drop-in GTF/GFF replacement and the REST
API can return a transcript's GRCh37, GRCh38, and T2T-CHM13v2.0 alignments in one
response. Build-specific fields (`contig`, `strand`, CDS bounds, exons; Supplementary
Table S3) include a per-exon `gap` string (e.g. `"M196 I1 M61"`) recording
transcript-vs-genome indels, converted to HGVS CIGAR format at query time. A root-level
`schema_version` lets clients reject incompatible files on load.

Canonical transcript tags (`mane_select`, `mane_plus_clinical`, `refseq_select`, and
`ensembl_canonical`) are stored per build where present. MANE Select covers
>{{ literature.mane_coverage_pct | dp(0) }}% of protein-coding genes [@Morales2022]
and is available for GRCh38; Ensembl canonical tags cover GRCh37 and GRCh38.

## String cleaning (`clean_hgvs()`)

Many of the descriptions that reach a curation platform do not parse: copy/paste
whitespace, lost casing, transposed punctuation, a swapped gene symbol and accession, a
trailing protein annotation. The intended variant is usually unmistakable, but a strict
grammar rejects the string; cdot repairs it instead.

`clean_hgvs()` is a pure string operation: it takes an HGVS string and returns the
repaired string with structured `HGVSFix` records describing every change. It needs no
genome build, sequence, parser, or data provider, so it can run as a pre-pass anywhere,
whichever downstream library consumes the result. Cleaning is an ordered pipeline of
single-purpose operations; callers may restrict it to a subset, which filters but never
reorders. The operations strip extraneous characters, repair structural punctuation,
normalise casing and prefixes, and rebuild the canonical `transcript(GENE):c.` shape,
including a transposed gene symbol and accession; the per-operation catalogue, with
examples, is Supplementary Table S7. Each `HGVSFix` carries a severity, a stable code,
a message, and the before/after values, so a caller can audit exactly what changed. A
final validation pass flags what cleaning cannot fix (no colon at all, a bare variant
body, an insertion given as a length instead of a sequence) as `ERROR`-level fixes; by
default cleaning never raises and returns its best attempt.

One repair needs transcript data: a bare-number accession whose RefSeq prefix was
dropped entirely (`000059.4:c.68del`). `resolve_missing_accession_prefix()`, applied by
`fix_hgvs()` when a data provider is supplied, generates the candidate accessions the
kind letter allows and restores the prefix only when exactly one exists in the loaded
data, reported as an `HGVSFix` like every other repair.

## Version fallback

A separate, opt-in helper, `get_best_transcript_version()`, addresses transcript-version
drift. Substituting a different version is a heuristic that can be wrong, so it is
never applied automatically: when the requested version is absent from a caller's data,
the helper returns the best available adjacent version (up-then-down, closest, or
latest strategies) as a reported `HGVSFix`, and the caller decides whether to accept
it. Before a substitution is offered as safe, cdot applies a build-independent check: a
transcript version's intrinsic CDS structure (CDS length plus coding-exon segment
lengths in transcript coordinates) is the same in every build, and a substitution
passes only when the two versions' structures match, their CDS alignment gaps match,
and, for a UTR variant, the relevant UTR length is preserved
(`intrinsic_cds_structure()`, `is_version_substitution_safe()`; Supplementary Methods).
The check is validated against ClinVar in Results R5.

## Access and client libraries

cdot's client stack (Figure 1B) offers three data providers behind the same interface.

**Local JSON**: `JSONDataProvider` loads a JSON.gz file into memory
(~{{ benchmark.grch38_load_time_s | dp(0) }} seconds for GRCh38 RefSeq) with lazy
interval trees and lookup dictionaries; retrieval is O(1) and end-to-end resolution
throughput ~{{ benchmark.cdot_local_tps | commas }} HGVS/second (Results, Table 1).

**REST API**: `cdot_rest` (https://github.com/SACGF/cdot_rest) serves the same JSON
data at cdotlib.org. `RESTDataProvider` fetches one transcript version per request,
suited to occasional lookups without downloading the full file, and can warm its cache
with a single batched `prefetch()` request (Results R3).

**Ensembl TARK**: `EnsemblTarkDataProvider` exposes the Ensembl Transcript Archive
(TARK) REST service through the same interface (Discussion).

**biocommons/hgvs integration**: cdot implements the full
`biocommons.hgvs.dataproviders.interface.Interface`, so it is a drop-in replacement for
UTA in any biocommons/hgvs pipeline:

```python
from cdot.hgvs.dataproviders import JSONDataProvider
hdp = JSONDataProvider(["cdot.0.2.34.refseq.grch38.json.gz"])
```

Sequence data comes from SeqRepo [@Hart2020] or cdot's `FastaSeqFetcher`, which
reconstructs a transcript sequence by splicing its exon ranges out of a local genome
FASTA, enabling fully offline operation. That reproduces an Ensembl transcript exactly
but is not guaranteed to match a curated RefSeq transcript, whose sequence can differ
from the genome at a few positions (Supplementary Methods). `ChainedSeqFetcher` tries
several sequence sources in a caller-defined order, so a pipeline can prefer SeqRepo
and fall back to the FASTA, or the reverse.

**PyHGVS integration**: `JSONPyHGVSTranscriptFactory` exposes the same data to the
Counsyl PyHGVS library. PyHGVS is no longer maintained, so new development targets the
biocommons/hgvs path.

## Canonical transcript selection

`resolve_gene_hgvs()` maps a gene symbol to a MANE Select (or MANE Plus Clinical)
transcript accession [@Morales2022; @Wright2023] via the data provider's
`get_tx_ac_tags_for_gene()`, which ranks a gene's transcripts by tag priority, so a
pipeline given `BRCA1:c.68_69del` can resolve it through `NM_007294.4`. The selection
is a sensible default for search and gene-name lookup, not an authoritative clinical
choice.

## Benchmarking

Resolution accuracy and throughput are measured with scripts committed to the
repository (`paper/scripts/`); protocol detail beyond what follows is in Supplementary
Methods. `benchmark_resolution.py` resolves real (g.HGVS, c.HGVS) pairs through a
pluggable provider (local JSON, REST, or UTA) and reports resolution rate, recovery
from cleaning and version fallback, and speed; the ClinVar pair set is built by
`build_clinvar_pairs.py`. Cleaning is evaluated on a production query corpus and, as a
reproducible control, with `inject_and_clean.py`, which injects each fix category into
clean ClinVar strings. `lovd_head_to_head.py` runs the same injected cases through both
`clean_hgvs()` and the LOVD HGVS syntax checker [@LovdHgvsChecker]
({{ lovd_comparison.lovd_version }}, run locally as a PHP CLI), scoring a case as
recovered when the tool's output (for LOVD, its top-ranked correction) exactly matches
the canonical target; the same false-correction check runs over the uncorrupted
originals.

`vv_mutalyzer_head_to_head.py` extends the same cases and scoring rule to the two
sequence-aware validation services, VariantValidator [@Freeman2018] and Mutalyzer
[@Lefter2021], over their public REST APIs
({{ vv_mutalyzer_comparison.service_date }}; remote services cannot be version-pinned,
so the facts record the service date and the version metadata each API reports).
VariantValidator requests are routed by accession family (its Ensembl transcripts live
on a separate endpoint) and a case counts as recovered when the validated transcript
description matches the target; for Mutalyzer a match by either its
`corrected_description` or its `normalized_description` counts, since the normalizer
may legitimately re-shift a representation. On the uncorrupted originals, a valid input
the service *alters* is a false correction, while one it *rejects* (for example
Mutalyzer's `EINTRONIC` for intronic positions on a transcript reference, or a
transcript absent from VariantValidator's database) is counted separately as a
validity or coverage position, matching how LOVD's flagged-invalid originals are
reported. Version-fallback safety is measured by `compute_version_stability.py` on
GRCh38, over a seeded {{ version_stability.sample_n | commas }}-accession sample of
accessions cdot holds at two or more versions; the same run bins preserved coding bases
by relative CDS position (Supplementary Figure S1).

The submitted-string corpus (Results R2) is built by
`build_clinvar_submitted_pairs.py` from the per-submission HGVS attributes of a ClinVar
VCV XML release (ClinVarVCVRelease_2026-06): each SCV's first transcript c./n.
expression is joined verbatim to the variant's VCF coordinate via its AlleleID and
collapsed to unique (AlleleID, string) pairs (Supplementary Methods). Samples are drawn
with a fixed seed (42): 3,000 pairs for the cdot-versus-UTA comparison and a committed
500-pair sample. Scoring uses the VCF coordinate rather than the g.HGVS string, since a
submitted string may legitimately spell an indel differently from ClinVar's normalised
form. Transcript version age is computed by `compute_submitted_version_age.py` against
the released cdot RefSeq JSON, whose per-transcript source URL identifies whether an
entry survives only through cdot's historical releases.

Throughput (Table 1) is measured by `compute_benchmark.py`: every configuration
resolves the identical committed set of {{ benchmark.n_pairs | commas }} ClinVar pairs
through the same biocommons/hgvs engine with the sequence layer held constant (one
shared local SeqRepo, file-descriptor cache enabled), so only the transcript-data layer
varies between rows. Each configuration is timed over
{{ benchmark.n_repeats | int }} passes and reported as median (IQR), steady-state
resolution only (exclusions in Supplementary Methods). The REST rows run against the
production server at cdotlib.org over the public internet; the public remote UTA row
uses the first {{ benchmark.uta_remote_n | int }} pairs, a full-size pass being
impractical at ~{{ benchmark.uta_remote_tps | dp(1) }} HGVS/s.
