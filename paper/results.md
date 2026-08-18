# Results

> **Provenance flags.** Results are reproducible from public data committed to this repo
> unless marked **[Tier 2]**, which denotes aggregate statistics from a private production
> corpus, published as frozen constants and not reproducible by a referee (see Methods /
> data availability).

---

## R1: Transcript coverage

The unit of coverage is the *transcript-version alignment*: a particular transcript
version aligned to a particular genome build, counted separately per build because each
alignment is what a resolution against that build needs. cdot covers
{{ coverage.total_count | commas }} such alignments across all builds and sources,
compared with ~{{ literature.uta_count | commas }} in UTA, an increase of
{{ coverage.improvement_fold | fmt('.1f') }}×. The gain has two sources: historical
depth (cdot ingests the complete run of archived RefSeq and Ensembl releases, so an NM_
version cited in an older report or ClinVar submission still resolves after NCBI has
retired it from the current annotation) and Ensembl itself, absent from UTA entirely
({{ coverage.ensembl_unique_count | commas }} accessions present in cdot).
T2T-CHM13v2.0 adds a further {{ coverage.t2t_unique_count | commas }} alignments,
making cdot the first transcript data source to bring that assembly to the Python HGVS
libraries (biocommons/hgvs and PyHGVS); Ensembl VEP can already generate HGVS against
T2T, but it is a standalone annotation tool, not a transcript backend for these
libraries. The JSON format also stores per-exon alignment-gap information (indels of
the transcript relative to the genome) that downstream libraries apply during
coordinate conversion.

## R2: ClinVar and clinical resolution accuracy

To measure practical impact, a seeded sample of
{{ clinvar.n_variants | commas }} ClinVar [@Landrum2025] variant descriptions, spanning
RefSeq (NM_) and Ensembl (ENST) transcripts, was resolved against cdot and a locally
loaded UTA (release `uta_20241220`) through the identical biocommons/hgvs code path,
with sequences from a shared local SeqRepo so only the transcript-data layer differs
(the comparison is gated to a sample by UTA's throughput; the full ClinVar set is
resolved through cdot alone below). cdot resolved
{{ clinvar.cdot_resolution_pct | dp(1) }}% versus
{{ clinvar.uta_resolution_pct | dp(1) }}% for UTA. The gap is Ensembl: on the committed
ClinVar test pairs ({{ clinvar.n_refseq | commas }} RefSeq and
{{ clinvar.n_ensembl | commas }} Ensembl), with cdot serving transcript sequence from a
local genome FASTA so no pair is dropped for a sequence SeqRepo lacks, RefSeq is at
parity (cdot {{ clinvar.cdot_refseq_pct | dp(1) }}%, UTA
{{ clinvar.uta_refseq_pct | dp(1) }}%), while on Ensembl cdot resolves
{{ clinvar.cdot_ensembl_pct | dp(1) }}% and UTA
{{ clinvar.uta_ensembl_pct | dp(1) }}%, UTA storing no Ensembl alignments at all
(Supplementary Table S4).

At full scale, resolving every RefSeq and Ensembl c.HGVS in ClinVar through cdot alone
({{ clinvar_vcf.n_pairs | commas }} (g., c.) pairs) reaches
{{ clinvar_vcf.resolved_pct | dp(1) }}% resolution, and
{{ clinvar_vcf.matched_of_resolved_pct | dp(1) }}% of resolved variants reproduce
ClinVar's own VCF coordinate exactly, scored as a VCF coordinate rather than a g.HGVS
string (Methods). The residual {{ clinvar_vcf.incorrect_pct | dp(2) }}% is dominated by
paralog and copy-number transcripts that map to more than one genomic locus and by
indel-representation differences; the per-source split and residual breakdown are in
Supplementary Table S4.

Both corpora above take their c.HGVS from the `Name` column of ClinVar's
variant_summary, ClinVar's own recomputed preferred-transcript name, always at the
current transcript version. They therefore measure a current-version ceiling and cannot
exercise historical transcript depth.

### What laboratories actually submit

To measure what laboratories actually write, we built a second public corpus from the
per-submission (SCV) HGVS attributes of the ClinVar VCV XML: each submitted string kept
verbatim and joined to the variant's VCF coordinate as ground truth via its AlleleID
({{ clinvar_submitted.n_unique_pairs | commas }} unique submitted-string/variant pairs
from {{ clinvar_submitted.n_scv_tx_strings | commas }} SCV transcript expressions;
Methods). The corpus is entirely RefSeq (not one submitted string cites an Ensembl
transcript), and its version profile confirms that submitted traffic is historical:
{{ clinvar_submitted.version_not_current_pct | dp(1) }}% of submitted strings cite a
transcript version that is no longer the version in the current RefSeq annotation
release ({{ clinvar_submitted.scv_weighted_not_current_pct | dp(1) }}% weighted by
submission count), and {{ clinvar_submitted.base_retired_pct | dp(1) }}% cite a
transcript no longer annotated at any version. cdot's merged release history holds
{{ clinvar_submitted.not_current_in_cdot_pct | dp(1) }}% of those superseded versions;
only {{ clinvar_submitted.absent_cdot_pct | dp(1) }}% of cited versions are absent from
its GRCh38 data.

On a fixed-seed sample of {{ clinvar_submitted.n_sample | commas }} submitted pairs,
cdot resolved
{{ clinvar_submitted.cdot_resolved_pct | dp(1) }}% and reproduced ClinVar's VCF
coordinate for {{ clinvar_submitted.cdot_matched_pct | dp(1) }}%, versus
{{ clinvar_submitted.uta_resolved_pct | dp(1) }}% for the same locally loaded UTA. On
submitted strings, the RefSeq gap invisible at the current-version ceiling (both
backends {{ clinvar.cdot_refseq_pct | dp(1) }}% above) opens to
{{ clinvar_submitted.cdot_only_pct | dp(1) }} points: UTA holds no GRCh38 alignment for
{{ clinvar_submitted.uta_no_data_pct | dp(1) }}% of the cited versions, and
{{ clinvar_submitted.cdot_only | commas }} of the
{{ clinvar_submitted.n_sample | commas }} pairs resolve through cdot alone (one through
UTA alone). The submitted strings are largely well-formed, so cleaning has little to
rescue (`fix_hgvs()` repaired {{ clinvar_submitted.rescued_by_fix | int }} string, with
{{ clinvar_submitted.regressions | int }} regressions). The residual
{{ clinvar_submitted.residual_pct | dp(1) }}% after cleaning and version fallback
({{ clinvar_submitted.residual_n | int }} of
{{ clinvar_submitted.n_sample | commas }}) is dominated by version effects, not
formatting: {{ clinvar_submitted_residual.version_refused | int }} cite a version
absent from the data, where the fallback declines to substitute because coordinate
safety cannot be verified, and
{{ clinvar_submitted_residual.coordinate_drift | int }} resolve through the cited
historical version to a coordinate that differs from ClinVar's current interpretation;
the full breakdown is in Supplementary Table S8.

**[Tier 2].** The same gap holds on the historical clinical data that motivated cdot:
the complete set of {{ historical.n_lines | commas }} unique HGVS descriptions imported
into the Australian Genomics Shariant variant-sharing platform [@Tudini2022],
classifications submitted by clinical laboratories over many years, each written
against whichever transcript version was current at the time
({{ historical.n_multi_version | commas }} of the
{{ historical.n_unique_tx | commas }} distinct transcripts are cited at more than one
version). Through the identical biocommons engine with only the transcript-data layer
swapped (Supplementary Methods), cdot produced a genomic coordinate for
{{ historical.cdot_resolved_pct | dp(1) }}% versus
{{ historical.uta_resolved_pct | dp(1) }}% for the same locally loaded UTA release. Of
the strings cdot resolved but UTA could not
({{ historical.cdot_only_pct | dp(1) }}% of the corpus),
{{ historical.cdot_only_historical_pct | dp(0) }}% were RefSeq transcript versions for
which UTA holds no GRCh38 alignment and
{{ historical.cdot_only_ensembl_pct | dp(0) }}% were Ensembl transcripts. There is no
ground-truth genomic coordinate for the private corpus, so the metric is resolution
rate rather than correctness (Methods, data availability).

## R3: Throughput

Table 1 compares backends over the identical committed ClinVar pair set with the
engine and sequence layer held constant (Methods).

**Table 1. End-to-end HGVS resolution throughput by transcript backend**: median (IQR)
HGVS/s over {{ benchmark.n_repeats | int }} timed passes of the identical
{{ benchmark.n_pairs | commas }}-pair ClinVar set, sequence layer held constant,
identical biocommons/hgvs engine; steady-state resolution only (Methods). The public
remote UTA row uses the first {{ benchmark.uta_remote_n | int }} pairs of the same set.

| Configuration | Throughput (HGVS/s), median (IQR) |
|---|---|
| UTA: public remote database | {{ benchmark.uta_remote_tps | dp(2) }} ({{ benchmark.uta_remote_tps_q1 | dp(2) }}–{{ benchmark.uta_remote_tps_q3 | dp(2) }}) |
| UTA: local PostgreSQL | {{ benchmark.uta_local_tps | int }} ({{ benchmark.uta_local_tps_q1 | int }}–{{ benchmark.uta_local_tps_q3 | int }}) |
| cdot REST (one request per transcript) | {{ benchmark.cdot_rest_tps | int }} ({{ benchmark.cdot_rest_tps_q1 | int }}–{{ benchmark.cdot_rest_tps_q3 | int }}) |
| cdot REST (after one batch `prefetch()`) | {{ benchmark.cdot_rest_prefetch_tps | int }} ({{ benchmark.cdot_rest_prefetch_tps_q1 | int }}–{{ benchmark.cdot_rest_prefetch_tps_q3 | int }}) |
| cdot local JSON | {{ benchmark.cdot_local_tps | int }} ({{ benchmark.cdot_local_tps_q1 | int }}–{{ benchmark.cdot_local_tps_q3 | int }}) |

A GRCh38 RefSeq JSON file loads in ~{{ benchmark.grch38_load_time_s | dp(0) }} s and
resolves at {{ benchmark.cdot_local_tps | int }} HGVS/s (median). Batching the
per-transcript REST lookups into one `prefetch()` request (all transcripts for the set,
under a second, untimed) makes REST throughput equivalent to local JSON, the two
differing by under 1% across repeats: with the transcript data in process memory, both
are bounded by the shared engine and sequence layer, not by the transcript backend. A
locally loaded UTA reached {{ benchmark.uta_local_tps | int }} HGVS/s, about a quarter
of local-JSON throughput, each lookup being a set of SQL queries instead of a dict hit.
The public remote UTA database, at {{ benchmark.uta_remote_tps | dp(2) }} HGVS/s, is
nearly four orders of magnitude slower than any local configuration, paying wide-area
round trips to a shared server on every lookup.

At scale, a single local-JSON process resolved the entire set of 3,660,452 unique
ClinVar (g.HGVS, c.HGVS) pairs in ~92 minutes (665 HGVS/s; 99.3% produced a genomic
coordinate, 98.8% matched the ClinVar genomic HGVS exactly), and the REST provider,
after one batch cache warm (21,277 distinct transcripts in ~6 s), in ~83 minutes
(731 HGVS/s). The two runs' sequence-layer cache conditions differed, so their small
difference is not evidence of REST outrunning local JSON; the controlled comparison is
Table 1. The same exhaustive pass against the public remote UTA database extrapolates
to close to a year.

## R4: String cleaning recovers malformed real-world HGVS

**[Tier 2].** The main test of cleaning is a production query stream: N = 32,752 real
queries typed into the HGVS search box of clinical and research variant-curation
platforms based on VariantGrid [@VariantGrid]. The strings are whatever a clinician or
curator pasted or typed, carrying the damage of their route to the box: whitespace and
non-printable characters from Word documents and report PDFs, lost casing, transposed
punctuation, trailing protein annotations. The cleaning pipeline (`clean_hgvs()` plus
the provider-verified accession-prefix restoration, Methods) raised the fraction
parseable by biocommons/hgvs from 91.5% as-submitted to 96.7%, a +5.3% absolute gain
(1,721 strings rescued, about 62% of the 8.5 percentage points that failed
as-submitted) with zero regressions (no already-valid string was broken). Table 2
breaks the rescues down by fix type: whitespace removal and structural-punctuation
repair dominate, followed by gene/transcript-wrapper repair and structure
reconstruction. The rarest category, restoring a fully dropped accession prefix
against the loaded transcript data, is impossible for a purely string-level checker;
having the transcript data locally is what makes it safe.

The public submitted-string corpus (R2) is the complement, and the two corpora locate
the two failure axes: formal database submissions are largely
well-formed and fail by transcript-version age, whereas interactive human-typed input
carries the formatting damage cleaning repairs. cdot addresses the first with
historical transcript depth and the second with `clean_hgvs()`.

**Table 2. Fixes applied across the production corpus (N = 32,752).** Each row is a
cleaning fix category, with the number of rescued queries in which it fired and its
share of the 1,721 rescued queries. Categories overlap (a single query may need several
fixes), so the counts sum to more than the total. *(Tier 2; counts are frozen constants
from a deterministic run of the cleaning pipeline, `clean_hgvs()` plus the
provider-verified accession-prefix restoration, over the production corpus.)*

| Fix category | Example (→ repaired) | Rescued queries | % of rescued |
|---|---|---|---|
| Whitespace / non-printable removal | `NM_000059.4: c.1A>G` → `NM_000059.4:c.1A>G` | 664 | 38.6% |
| Structural-punctuation repair | `NM_000059.4;c.1A>G` → `NM_000059.4:c.1A>G` | 640 | 37.2% |
| Gene/transcript-wrapper repair | `BRCA2(NM_000059.4):c.1A>G` → `NM_000059.4(BRCA2):c.1A>G` | 260 | 15.1% |
| Structure reconstruction | `NM_000059.4c.1A>G` → `NM_000059.4:c.1A>G` | 248 | 14.4% |
| Protein-suffix stripping | `NM_000059.4:c.1A>G p.(Met1?)` → `NM_000059.4:c.1A>G` | 130 | 7.6% |
| Genomic-ref-in-parens removal | `NM_000059.4(NC_000013.11):c.68del` → `NM_000059.4:c.68del` | 57 | 3.3% |
| Base re-casing | `NM_000059.4:c.1delg` → `…delG` | 42 | 2.4% |
| Other (del/dup count, mutation-type case, …) | `NM_000059.4:c.1_2del2` → `…del` | 11 | 0.6% |
| Prefix / kind restoration | `NM_000059.4:1A>G` → `NM_000059.4:c.1A>G` | 4 | 0.2% |
| Accession prefix restoration (provider-verified) | `000059.4:c.68del` → `NM_000059.4:c.68del` | 1 | 0.1% |
| **Total unique queries rescued** | | **1,721** | **100%** |

As a control, `paper/scripts/inject_and_clean.py` injects each
`clean_hgvs()` fix category into a seeded sample of clean, parseable ClinVar c.HGVS
strings and confirms the cleaner recovers the canonical target with
**{{ cleaning.inject_regressions | int }} regressions**: no already-valid string is
ever broken (Supplementary Table S5). Because this benchmark injects the very errors it
then repairs, its recovery rate is reported in the supplement only; its purpose is the
no-regression guarantee on which the production result depends.

### Residual errors: the ceiling of cleaning *(Table S6)*

**[Tier 2].** The 3.3% of the production corpus (1,075 queries; 826 unique strings)
that still fail to parse after cleaning define the ceiling of pure string repair.
Classified under a fixed decision-tree taxonomy (Supplementary Table S6, with
synthesised examples and the classification method and its limitations), just over half
is incomplete or reference-less input that no string-level repair can invent, about 30%
is in principle fixable and marks the frontier for future cleaning rules, and the
remainder splits evenly between valid HGVS the biocommons grammar rejects and non-HGVS
input (pasted URLs, report templates, prose) that should not be parsed at all.

## R5: Transcript version fallback and safe substitution

When a cited transcript version is absent from the loaded data, the opt-in fallback
(Methods) substitutes an adjacent version only if the coordinate-safety check passes;
a substitution that cannot be verified safe is refused by default, preserving
exact-version semantics. The fallback is a client-layer feature: biocommons/hgvs has no
adjacent-version fallback with any data provider. In an end-to-end ablation
(`paper/scripts/benchmark_resolution.py`) that removes the requested version from each
test variant, the fallback recovered the correct genomic coordinate with no false
rescues (a false rescue being a substitution that resolves to a different coordinate).

Caution is warranted because a version bump can move a variant. Across consecutive
RefSeq version bumps ({{ version_stability.refseq_pairs | commas }} pairs),
{{ version_stability.refseq_preserving_pct | dp(1) }}% preserved every coding
coordinate; for Ensembl ({{ version_stability.ensembl_pairs | commas }} pairs)
{{ version_stability.ensembl_preserving_pct | dp(1) }}%. Weighted by coding base (the
chance a *random* variant is unaffected), safety is
{{ (version_stability.refseq_pervariant_safety * 100) | dp(1) }}% (RefSeq) and
{{ (version_stability.ensembl_pervariant_safety * 100) | dp(1) }}% (Ensembl). When a
coordinate does move it is almost always the whole CDS, driven by a re-annotation of
the coding region; the most dangerous case, a *partial* bump that mis-places some
variants but not others, is rare
({{ version_stability.refseq_partial_drift_pct | dp(1) }}% of RefSeq bumps,
{{ version_stability.ensembl_partial_drift_pct | dp(1) }}% Ensembl). Within that tail
the risk is positional: a partial drift keeps a 5' prefix intact up to its first
alignment change and moves what lies downstream, so preservation falls from
{{ positional_drift.refseq_partial_decile1_pct | dp(0) }}% of coding bases in the
5'-most decile to {{ positional_drift.refseq_partial_decile10_pct | dp(0) }}% in the
3'-most for RefSeq ({{ positional_drift.ensembl_partial_decile1_pct | dp(0) }}% to
{{ positional_drift.ensembl_partial_decile10_pct | dp(0) }}% for Ensembl; Supplementary
Figure S1). Drift is also concentrated: only
{{ version_stability.refseq_accessions_drift_pct | dp(1) }}% of RefSeq accessions ever
drift, so the risky transcripts form a short, identifiable list that can be flagged or
blocklisted.

The intrinsic-structure check (Methods) works even when the requested version was never
aligned to the target build, because the structure is build-independent, identical
across GRCh37/GRCh38/T2T for
{{ version_stability.refseq_struct_portable_pct | dp(1) }}% of RefSeq and
{{ version_stability.ensembl_struct_portable_pct | dp(1) }}% of Ensembl versions
(predictor detail in Supplementary Methods). Validated against ClinVar: of the
{{ version_safety_validation.n_safe_substitutions | commas }} version substitutions the
structure, gap, and UTR checks accept across the full RefSeq corpus, only
{{ version_safety_validation.n_coordinate_changes }} still move the genomic coordinate,
both a single transcript that re-aligns the same CDS structure to a different genomic
locus. A genomic CDS-map comparison and a small precomputed blocklist catch these
(Supplementary Methods); with them, no accepted substitution moves a coordinate.
