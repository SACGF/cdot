# Results

> **Provenance flags.** Results are reproducible from public data committed to this repo
> unless marked **[Tier 2]**, which denotes aggregate statistics from a private production
> corpus, published as frozen constants and not reproducible by a referee (see Methods /
> data availability).

---

## R1: Transcript coverage

The unit of coverage is the *transcript-version alignment*: a particular
transcript version aligned to a particular genome build. The same transcript version is
counted separately per build, because each alignment is what a resolution against that
build actually needs. cdot covers {{ coverage.total_count | commas }} such alignments
across all builds and sources, compared with ~{{ literature.uta_count | commas }} in UTA,
an {{ coverage.improvement_fold | fmt('.1f') }}× increase. Two things drive the
gain. The first is historical depth: UTA does hold several versions per transcript,
but cdot ingests the complete run of RefSeq and Ensembl releases from the FTP
archives, so it retains many more historical versions. An NM_ version cited in an older
clinical report or ClinVar submission still resolves even after NCBI has retired it from
the current annotation. The second is Ensembl, which cdot covers in
full: {{ coverage.ensembl_unique_count | commas }} Ensembl transcript accessions are
present in cdot but absent from UTA entirely. T2T-CHM13v2.0 adds a further
{{ coverage.t2t_unique_count | commas }} alignments.

Those T2T-CHM13v2.0 alignments make cdot the first transcript data source to bring that
assembly to the Python HGVS libraries (biocommons/hgvs and PyHGVS); they are absent from
UTA and from every other backend these libraries can use. (Ensembl VEP can already
generate HGVS against T2T, but it is a standalone annotation tool, not a transcript
backend for these libraries.) The JSON format also stores per-exon alignment-gap
information (indels of the transcript relative to the genome) so that downstream libraries
can apply it during coordinate conversion; this is a property of the stored data rather
than a separate result.

## R2: ClinVar and clinical resolution accuracy

To measure practical impact rather than raw counts, a seeded sample of
{{ clinvar.n_variants | commas }} ClinVar [@Landrum2025] variant descriptions,
spanning both RefSeq (NM_) and Ensembl (ENST) transcripts, was resolved against
cdot and a locally loaded UTA (release `uta_20241220`) through the identical
biocommons/hgvs code path (genomic→transcript projection, with sequences served from a
shared local SeqRepo so only the transcript-data layer differs). A sample is used here
because this comparison is gated by UTA's throughput; the full ClinVar set is resolved
through cdot alone in R3. cdot resolved
{{ clinvar.cdot_resolution_pct | dp(1) }}% versus
{{ clinvar.uta_resolution_pct | dp(1) }}% for UTA. To show where that gap comes from, we
resolved the committed ClinVar test pairs ({{ clinvar.n_refseq | commas }} RefSeq and
{{ clinvar.n_ensembl | commas }} Ensembl, reproducible from this repo) through each
backend, with cdot serving transcript sequence from a local genome FASTA via its
`FastaSeqFetcher` so no pair is dropped for a sequence SeqRepo happens not to hold. RefSeq
is at parity, with both backends holding every cited current version
(cdot {{ clinvar.cdot_refseq_pct | dp(1) }}%, UTA {{ clinvar.uta_refseq_pct | dp(1) }}%),
whereas on Ensembl cdot resolves {{ clinvar.cdot_ensembl_pct | dp(1) }}% and UTA
{{ clinvar.uta_ensembl_pct | dp(1) }}%, because UTA stores no Ensembl alignments at all
(Supplementary Table S4). The benchmark is reproducible: the ClinVar build script, random
seed, and resolution harness are committed.

At full scale, resolving every RefSeq and Ensembl c.HGVS in ClinVar through cdot alone
({{ clinvar_vcf.n_pairs | commas }} (g.,c.) pairs) reaches
{{ clinvar_vcf.resolved_pct | dp(1) }}% resolution, and of the variants resolved
{{ clinvar_vcf.matched_of_resolved_pct | dp(1) }}% reproduce ClinVar's own VCF
coordinate exactly. We score the projection as a VCF coordinate (CHROM/POS/REF/ALT)
rather than as a g.HGVS string, so equivalent representations (3'-shift, del/dup
spellings) and the ClinVar variants written with tandem-repeat or identity notation are
not miscounted as disagreements. The residual {{ clinvar_vcf.incorrect_pct | dp(2) }}% is
dominated by paralog and copy-number transcripts that map to more than one genomic locus
and by indel-representation differences, not by coordinate errors; the per-source split
and the residual breakdown are in Supplementary Table S4.

**[Tier 2].** The same gap holds on the
historical clinical data that motivated cdot. We resolved the complete set of
{{ historical.n_lines | commas }} unique HGVS descriptions imported into the Australian
Genomics Shariant variant-sharing platform [@Tudini2022]: real classifications submitted
by clinical genetic-testing laboratories over many years, each written against whichever
transcript version was current when the variant was first classified
({{ historical.n_multi_version | commas }} of the {{ historical.n_unique_tx | commas }}
distinct transcripts are cited at more than one version). Run through the identical
biocommons engine with only the transcript-data layer swapped (pure coordinate
projection, `replace_reference=False`, so the sequence layer never differs), cdot
produced a genomic coordinate for {{ historical.cdot_resolved_pct | dp(1) }}% of them
versus {{ historical.uta_resolved_pct | dp(1) }}% for the same locally loaded UTA release
(`uta_20241220`). Of the strings cdot resolved but UTA could not
({{ historical.cdot_only_pct | dp(1) }}% of the corpus),
{{ historical.cdot_only_historical_pct | dp(0) }}% were RefSeq transcript versions for
which UTA holds no GRCh38 alignment and {{ historical.cdot_only_ensembl_pct | dp(0) }}%
were Ensembl transcripts (absent from UTA entirely). The ClinVar comparison above shows
the GTF→JSON→biocommons pipeline works at scale on current transcript versions; this
corpus shows the historical depth matters in practice: it is the
older-version traffic a working clinical lab generates. There is no ground-truth genomic
coordinate for the private corpus, so the metric is resolution rate rather than
correctness (Methods, data availability).

## R3: Throughput

Backends were compared with the sequence layer held constant, so the only thing that
varies across rows of Table 1 is the transcript-data layer.

**Table 1. End-to-end HGVS resolution throughput by transcript backend**, sequence layer
held constant (shared local SeqRepo), identical biocommons/hgvs engine.

| Configuration | Throughput (HGVS/s) |
|---|---|
| UTA: public remote database | ~{{ benchmark.uta_remote_tps | dp(1) }} |
| UTA: local PostgreSQL | ~{{ benchmark.uta_local_tps | int }} |
| cdot REST (one request per transcript) | ~{{ benchmark.cdot_rest_tps | int }} |
| cdot REST (after one batch `prefetch()`) | ~{{ benchmark.cdot_rest_prefetch_tps | int }} |
| cdot local JSON | {{ benchmark.cdot_local_min_tps | int }}–{{ benchmark.cdot_local_max_tps | int }} |

A GRCh38 RefSeq JSON file loads in ~{{ benchmark.grch38_load_time_s | dp(0) }} s and then
resolves at {{ benchmark.cdot_local_min_tps | int }}–{{ benchmark.cdot_local_max_tps |
int }} HGVS/s through the biocommons engine. The REST provider serves
~{{ benchmark.cdot_rest_tps | int }} HGVS/s when each transcript version is fetched in its
own request; batching those lookups into one `prefetch()` request amortises the
per-request network and processing overhead across the whole set and warms the transcript
cache. Later lookups are then in-memory hits, so REST throughput rises to match local. A
locally loaded UTA reached only ~{{ benchmark.uta_local_tps | int }} HGVS/s, and the
public remote UTA database only ~{{ benchmark.uta_remote_tps | dp(1) }} HGVS/s. cdot's
local data layer is thus roughly 30× faster than a local UTA on the identical engine.
Once the data layer is local, the remaining bottleneck is sequence fetching rather than
transcript lookup.

At scale, a single local-JSON process resolved the entire set of 3,660,452
unique ClinVar (g.HGVS, c.HGVS) pairs in ~92 minutes (665 HGVS/s; 99.3% produced a
genomic coordinate, 98.8% matched the ClinVar genomic HGVS exactly). The REST provider
matched this over the network: after one batch cache-warming pass (21,277 distinct
transcripts warmed in ~6 s) it resolved the same set in ~83 minutes (731 HGVS/s),
effectively matching local JSON: once the cache is warm both providers resolve from an
in-memory dict with no further network I/O, so throughput is bounded by the shared
sequence layer. The marginal difference is within run-to-run variance, plausibly because
the warmed REST cache holds only the few thousand transcripts the set actually touches
whereas local JSON holds the entire dataset in memory. The same exhaustive pass is
impractical against the public remote UTA database (extrapolated at hundreds of days from
its ~0.1 HGVS/s). (`paper/scripts/build_clinvar_pairs.py` builds the pair set by joining
ClinVar's variant_summary with the ClinVar VCF.)

## R4: String cleaning recovers malformed real-world HGVS

The main test of cleaning is its effect on a real production query stream; a
reproducible injection benchmark provides a supporting safety check.

**[Tier 2].** The corpus is
N = 32,752 real queries typed into the HGVS search box of production clinical and
research variant-curation platforms based on VariantGrid [@VariantGrid]. Users treat that box as a shortcut to jump straight
to a variant or its classification, so the strings are whatever a clinician or curator
happened to paste or type, not curated HGVS. They arrive carrying the damage of their
route to the box: stray whitespace and non-printable characters from copying out of Word
documents and report PDFs and pasting between systems, lost casing, transposed
punctuation, and trailing protein annotations. Run over this corpus, `clean_hgvs()` raised
the fraction parseable by biocommons/hgvs from 91.5% as-submitted to 96.6% after
cleaning: a +5.1% absolute gain (1,678 additional strings rescued) with zero
regressions (no already-valid string was broken). This measures what cleaning recovers
from messy, human-entered input rather than from synthetic errors. The rescues break down by fix type as shown in
Table 2; they are dominated by whitespace removal and base re-casing, followed by
protein-suffix stripping and gene/transcript-wrapper repair.

**Table 2. Fixes applied across the production corpus (N = 32,752).** Each row is a
`clean_hgvs()` fix category, with the number of rescued queries in which it fired and that
number as a share of the 1,678 rescued queries. Categories overlap (a single query may
need several fixes), so the counts sum to more than the total. *(Tier 2; counts are
frozen constants from a deterministic run of `clean_hgvs()` over the production corpus.)*

| Fix category | Example (→ repaired) | Rescued queries | % of rescued |
|---|---|---|---|
| Whitespace / non-printable removal | `NM_000059.4: c.1A>G` → `NM_000059.4:c.1A>G` | 940 | 56.0% |
| Base re-casing | `NM_000059.4:c.1delg` → `…delG` | 553 | 33.0% |
| Structural-punctuation repair | `NM_000059..4:c.1A>G` → `NM_000059.4:c.1A>G` | 310 | 18.5% |
| Gene/transcript-wrapper repair | `BRCA2(NM_000059.4):c.1A>G` → `NM_000059.4(BRCA2):c.1A>G` | 209 | 12.5% |
| Protein-suffix stripping | `NM_000059.4:c.1A>G p.(Met1?)` → `NM_000059.4:c.1A>G` | 130 | 7.7% |
| Genomic-ref-in-parens removal | `NM_000059.4(NC_000013.11):c.68del` → `NM_000059.4:c.68del` | 57 | 3.4% |
| Other (del/dup count, mutation-type case, …) | `NM_000059.4:c.1_2del2` → `…del` | 21 | 1.3% |
| Prefix / kind restoration | `NM_000059.4:1A>G` → `NM_000059.4:c.1A>G` | 11 | 0.7% |
| **Total unique queries rescued** | | **1,678** | **100%** |

As a reproducible control,
`paper/scripts/inject_and_clean.py` injects each `clean_hgvs()` fix category into a seeded
sample of clean, parseable ClinVar c.HGVS strings committed to this repo and confirms the
cleaner recovers the canonical target with **{{ cleaning.inject_regressions | int }}
regressions**: no already-valid string is ever broken (Supplementary Table S5). Because
this benchmark injects the very errors it then repairs, its recovery rate is not an
independent measure of real-world performance and is reported in the supplement only; its
purpose here is to demonstrate the no-regression guarantee on which the production result
depends.

### Residual errors: the ceiling of cleaning *(Table S6)*

**[Tier 2].** The 3.4% of the production corpus (1,118 queries; 860
unique strings) that still fail to parse after cleaning define the ceiling of what
pure string repair can achieve. Each residual string was assigned to one single-label
error class under a fixed decision-tree taxonomy. Of the eight classes, the seven
repair-relevant ones are shown below with synthesised examples; the eighth was non-HGVS
input (81 queries, 7.2%: pasted URLs, report templates, or prose) and is excluded from
the table, as there is nothing in it for cleaning to repair.

**Residual error classes after cleaning** (counts and % of the 1,118 residual
queries; examples synthesised from public BRCA2 `NM_000059.4`). *(Tier 2; frozen
constants from a deterministic run over the production corpus.)*

| Class | Queries | What it is (*example*) |
|---|---|---|
| Truncated | 284 (25.4%) | cut off before a complete variant: `NM_000059.4:c.68_69` (range, no edit) |
| No reference | 277 (24.8%) | a bare variant body, no transcript/gene/accession: `c.68_69delAG` |
| Bad accession | 167 (14.9%) | missing prefix, or misplaced/truncated version: `000059.4:c.68del` (`NM_` prefix dropped) |
| Edit syntax | 113 (10.1%) | malformed or non-standard edit operation: `NM_000059.4:c.68AG>T` (multi-base reference in a substitution) |
| Trailing / concatenated | 85 (7.6%) | extra characters after a complete variant, or several run together: `NM_000059.4:c.68delAG;c.70A>G` |
| Grammar gap | 81 (7.2%) | legitimate HGVS the biocommons grammar rejects: `NM_000059.4:c.(67+1_68-1)_(70+1_71-1)del` (uncertain-range deletion) |
| Insertion (length only) | 30 (2.7%) | an insertion given as a base count, not a sequence: `NM_000059.4:c.68_69ins5` (position and length recoverable; inserted bases not) |

The residual falls into three groups. Just over half (~50%: Truncated + No reference) is
incomplete or reference-less user input: information the user never supplied, which no
string-level repair can invent; the integer-length insertions (2.7%) belong with these,
since the inserted bases were never given (only the position and length are recoverable, so
the variant is not properly resolvable). About 33% (Bad accession + Edit syntax + Trailing /
concatenated) is in principle fixable and marks the frontier for future cleaning rules.
The remaining ~14% splits into two equal classes of 81 queries (7.2%) each: a grammar gap
(valid HGVS the biocommons grammar rejects rather than the input) and the non-HGVS input
excluded above (strings that should not be parsed at all). Most of what remains is
therefore either out of scope for any repair or a downstream grammar limitation.

*Method and limitation:* classification was performed by a large language model (Claude
Opus 4, Anthropic; 2026-06-17) applying the shared decision tree to each unique string,
single-label and single-rater; no second-rater adjudication was done, so no inter-rater
agreement (κ) is reported. The taxonomy is version `v1`. Synthesised examples (from public
NM_000059.4 / NM_001754.5) illustrate each class; no corpus string is reproduced.

## R5: Transcript version fallback and safe substitution

The opt-in adjacent-version fallback (Methods) lets a variant cited against a retired
transcript version still resolve. An end-to-end ablation
(`paper/scripts/benchmark_resolution.py`) removes the requested version from each test
variant and resolves through the fallback; it recovered the correct genomic coordinate
with no false rescues (a false rescue being a substitution that resolves to a different
coordinate). This is a client-layer feature rather than a backend one: biocommons/hgvs
has no adjacent-version fallback regardless of its data provider, so a UTA-backed pipeline
gains nothing here even though UTA stores several versions per transcript. cdot's
multi-release depth gives the fallback more versions to choose from, and it is never
applied automatically, preserving exact-version semantics by default.

The fallback is only safe if substituting an adjacent version
does not move the variant: a coding `c.` position must still map to the same genomic
coordinate. We measured this directly on the released cdot data (Methods, Benchmarking), exploiting
the fact that
cdot stores the genome *alignment* rather than the sequence, so a version bump can only
move a coordinate if it changes that alignment. Across consecutive RefSeq version bumps
({{ version_stability.refseq_pairs | commas }} pairs),
{{ version_stability.refseq_preserving_pct | dp(1) }}% preserved every coding coordinate;
for Ensembl ({{ version_stability.ensembl_pairs | commas }} pairs)
{{ version_stability.ensembl_preserving_pct | dp(1) }}%. Weighting by coding base (the
chance a *random* variant is unaffected), safety is higher still,
{{ (version_stability.refseq_pervariant_safety * 100) | dp(1) }}% (RefSeq) and
{{ (version_stability.ensembl_pervariant_safety * 100) | dp(1) }}% (Ensembl); and the most
dangerous case, a *partial* bump that mis-places some variants but not others, is rare
({{ version_stability.refseq_partial_drift_pct | dp(1) }}% of RefSeq bumps,
{{ version_stability.ensembl_partial_drift_pct | dp(1) }}% Ensembl); when a coordinate
does move it is almost always the whole CDS, driven by a re-annotation of the coding
region. Drift is also highly concentrated: only
{{ version_stability.refseq_accessions_drift_pct | dp(1) }}% of RefSeq accessions ever
drift, and within that set the moves cluster in a small minority of accessions.
This concentration has a practical consequence: the risky transcripts form a short,
identifiable list that can be flagged or blocklisted, rather than every version bump being
treated as unsafe.

Whether a bump is coordinate-safe is predictable from a single version, with no
flanking pair or genomic context. A transcript version's intrinsic CDS structure (its
CDS length plus the lengths of its coding-exon segments in transcript coordinates)
is build-independent, identical across GRCh37/GRCh38/T2T for
{{ version_stability.refseq_struct_portable_pct | dp(1) }}% of RefSeq and
{{ version_stability.ensembl_struct_portable_pct | dp(1) }}% of Ensembl versions. A
change in that structure is almost equivalent to a genomic-coordinate drift:
every drifting RefSeq bump ({{ version_stability.refseq_drift_struct_flagged_pct | dp(0) }}%)
and {{ version_stability.ensembl_drift_struct_flagged_pct | dp(1) }}% of drifting Ensembl
bumps carry an intrinsic-structure change, while conversely
{{ version_stability.refseq_struct_unchanged_preserved_pct | dp(0) }}% (RefSeq) /
{{ version_stability.ensembl_struct_unchanged_preserved_pct | dp(1) }}% (Ensembl) of
structure-unchanged bumps are genomically preserved. This makes the safety decision
work even when the requested version was never aligned to the target build. To judge
substituting an available version for a requested one in, say, GRCh38, cdot reads the
requested version's intrinsic CDS structure from any build that carries it (the structure
is build-independent, so a GRCh37 or T2T record serves), reads the substitute version's
structure from GRCh38, and substitutes only when the two match. Because the comparison is
on the structure rather than on flanking genomic coordinates, it also catches
transient-revert versions (a coordinate that goes A→B→A) that a genomic-bracket check
between neighbours cannot see.

Identical transcript structure does not on its own guarantee safety, so the gate
adds two refinements. An alignment gap (a transcript-vs-genome indel) present in one
version and not the other shifts every coding base downstream of it, so the gate also
requires the two versions' CDS alignment gaps to match. And because the structure is
build-independent it cannot see the variant's position relative to the UTRs: a 5' or 3'
UTR variant can move when only the UTR length changes (UTR annotations change between
versions far more often than the CDS), so for a cited UTR position the matching UTR length
must be preserved too, while a coding variant on a UTR-only change stays safe. We
validated the gate empirically against ClinVar: of the
{{ version_safety_validation.n_safe_substitutions | commas }} version substitutions the
structure, gap, and UTR checks accept across the full RefSeq corpus, only
{{ version_safety_validation.n_coordinate_changes }} still move the genomic coordinate,
both a single transcript that re-aligns the same CDS structure to a different genomic
locus. The build-independent structure cannot see such a re-placement, so cdot
additionally compares the two versions' genomic CDS maps when both are loaded, and ships a
small precomputed blocklist of these re-placements for the case where the requested
version is absent from the loaded data; with these, no accepted substitution moves a
coordinate.

cdot ships this as the safety gate behind the version fallback
(`intrinsic_cds_structure()` and the data-provider method
`is_version_substitution_safe()`): when the requested version is absent, the substitute
is applied only if it is coordinate-safe by these tests, reported as a coordinate-safe
`HGVSFix`; a structural mismatch (or a version absent from every build, where only a
probabilistic genomic-bracket check remains) is refused by default rather than silently
substituted, preserving exact-variant semantics.
