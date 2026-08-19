# Results

> **Provenance flags.** Results are reproducible from public data committed to this repo
> unless marked **[private data]**, which denotes aggregate statistics from a private
> production corpus, published as frozen constants and not reproducible by a referee (see
> Methods / data availability).

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
libraries (biocommons/hgvs and PyHGVS). The JSON format also stores per-exon
alignment-gap information (indels of the transcript relative to the genome) that
downstream libraries apply during coordinate conversion; Ensembl VEP, which can itself
generate HGVS, does not apply these transcript-to-genome gaps when converting a c.
description to genomic coordinates [@VepHgvsGaps], so it can misplace variants that fall
downstream of an indel in a gapped RefSeq alignment.

## R2: ClinVar and clinical resolution accuracy

Two comparisons measure resolution rate against a locally loaded UTA (release
`uta_20241220`): a controlled current-version check here, and the historical
submitted-string corpus below that is the real test of transcript depth.

For the current-version check, a seeded ClinVar [@Landrum2025] sample of
preferred-transcript descriptions, deliberately including both RefSeq (NM_) and Ensembl
(ENST) transcripts, was resolved against cdot and UTA through the identical
biocommons/hgvs code path, with sequences from a shared local SeqRepo so only the
transcript-data layer differs. On RefSeq the two backends are at parity (cdot
{{ clinvar.cdot_refseq_pct | dp(1) }}%, UTA {{ clinvar.uta_refseq_pct | dp(1) }}%; cdot
serving transcript sequence from a local genome FASTA so no pair is dropped for a
sequence SeqRepo lacks). The remaining difference is capability, not accuracy: UTA stores
no Ensembl alignments at all, so it resolves {{ clinvar.uta_ensembl_pct | dp(1) }}% of
the Ensembl pairs against cdot's {{ clinvar.cdot_ensembl_pct | dp(1) }}% (Supplementary
Table S4). Because this sample is deliberately Ensembl-enriched to exercise that
capability, its aggregate rate is not a like-for-like accuracy score, so we do not read a
single headline resolution figure from it; the fair head-to-head is RefSeq parity here
and the historical submitted corpus below.

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

### ClinVar submissions as a historical record of transcripts used

To read ClinVar's submission history as a record of the transcript versions labs
actually used over the years, we built a second public
corpus from the per-submission (SCV) HGVS attributes of the ClinVar VCV XML: each
submitted string kept verbatim and joined to the variant's VCF coordinate as ground
truth via its AlleleID ({{ clinvar_submitted.n_unique_pairs | commas }} unique
submitted-string/variant pairs from {{ clinvar_submitted.n_scv_tx_strings | commas }} SCV
transcript expressions; Methods). Unlike the deliberately Ensembl-enriched sample above,
this corpus reflects the natural submission mix, which is essentially all RefSeq (not one
submitted string cites an Ensembl transcript; Ensembl accessions enter ClinVar through
its recomputed names, not submitter descriptions). Its version profile confirms that
submitted traffic is historical:
{{ clinvar_submitted.version_not_current_pct | dp(1) }}% of submitted strings cite a
transcript version that is no longer the version in the current RefSeq annotation
release ({{ clinvar_submitted.scv_weighted_not_current_pct | dp(1) }}% weighted by
submission count), and {{ clinvar_submitted.base_retired_pct | dp(1) }}% cite a
transcript no longer annotated at any version. cdot's merged release history holds
{{ clinvar_submitted.not_current_in_cdot_pct | dp(1) }}% of those superseded versions;
only {{ clinvar_submitted.absent_cdot_pct | dp(1) }}% of cited versions are absent from
its GRCh38 data.

Because ClinVar grows over time, a flat random draw over-represents recent submissions,
so we scored two {{ clinvar_submitted.n_sample | commas }}-pair samples (seed 42): one
balanced evenly across submission-year eras (2008-2026), which represents the historical
record fairly, and one drawn at random, which reflects the current file's recency skew
(Methods).

On the fair era-balanced sample, cdot resolved
{{ clinvar_submitted.cdot_resolved_pct | dp(1) }}% and reproduced ClinVar's VCF
coordinate for {{ clinvar_submitted.after_fix_matched_pct | dp(1) }}%, versus
{{ clinvar_submitted.uta_resolved_pct | dp(1) }}% for the same locally loaded UTA. The
RefSeq gap invisible at the current-version ceiling (both backends
{{ clinvar.cdot_refseq_pct | dp(1) }}% above) opens to
{{ clinvar_submitted.cdot_only_pct | dp(1) }} points:
{{ clinvar_submitted.cdot_only | commas }} of the
{{ clinvar_submitted.n_sample | commas }} pairs resolve through cdot alone (one through
UTA alone), because UTA holds no GRCh38 alignment for
{{ clinvar_submitted.uta_no_data_pct | dp(1) }}% of the cited versions. The random draw
is easier on both backends ({{ clinvar_submitted.rnd_after_fix_matched_pct | dp(1) }}%
cdot versus {{ clinvar_submitted.rnd_uta_resolved_pct | dp(1) }}% UTA, a
{{ clinvar_submitted.rnd_cdot_only_pct | dp(1) }}-point gap): it is dominated by recent
submissions, and the fair sample is harder precisely because older submissions cite the
superseded versions this corpus exists to exercise.

The submitted strings are largely well-formed, so cleaning has little to rescue
(`fix_hgvs()` repaired {{ clinvar_submitted.rescued_by_fix | int }} string, with
{{ clinvar_submitted.regressions | int }} regressions). The residual
{{ clinvar_submitted.residual_pct | dp(1) }}% on the fair sample after cleaning and
version substitution ({{ clinvar_submitted.residual_n | int }} of
{{ clinvar_submitted.n_sample | commas }}) is dominated by historical-version effects,
not formatting: {{ clinvar_submitted_residual.coordinate_drift | int }} resolve through
the cited version to a coordinate that differs from ClinVar's current interpretation,
{{ clinvar_submitted_residual.reference_mismatch | int }} cite a reference base that
does not exist on the cited version, and
{{ clinvar_submitted_residual.version_refused | int }} cite a version absent from the
data where substitution declines to act because coordinate safety cannot be verified;
the full breakdown is in Supplementary Table S8.

**[private data].** The same gap holds on the historical clinical data that motivated cdot:
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
are bounded by the shared engine and sequence layer, not by the transcript backend.
Comparing like with like: locally, a loaded UTA reached
{{ benchmark.uta_local_tps | int }} HGVS/s, about a quarter of local-JSON throughput,
each lookup resolving a set of SQL queries where cdot needs a single JSON object; and
remote to remote, cdot's REST API at {{ benchmark.cdot_rest_tps | int }} HGVS/s (one
request per transcript, no prefetch) is more than two orders of magnitude faster than the
public UTA server's {{ benchmark.uta_remote_tps | dp(2) }} HGVS/s, which pays wide-area
round trips to a shared database on every lookup.

At scale, a single local-JSON process resolved the entire set of 3,660,452 unique
ClinVar (g.HGVS, c.HGVS) pairs in ~92 minutes (665 HGVS/s; 99.3% produced a genomic
coordinate, 98.8% matched the ClinVar genomic HGVS exactly), and the REST provider,
after one batch cache warm (21,277 distinct transcripts in ~6 s), in ~83 minutes
(731 HGVS/s). The two runs' sequence-layer cache conditions differed, so their small
difference is not evidence of REST outrunning local JSON; the controlled comparison is
Table 1. The same exhaustive pass against the public remote UTA database extrapolates
to close to a year.

## R4: String cleaning recovers malformed real-world HGVS

**[private data].** The main test of cleaning is a production query stream:
N = {{ cleaning_corpus.corpus_n | commas }}
real search-box queries from clinical and research variant-curation platforms based on
VariantGrid [@VariantGrid], restricted to the subset that matched a broad HGVS regex
(loosely HGVS-shaped strings a cleaner could plausibly repair, not arbitrary free-text
search terms), with a small residue of non-HGVS input that slipped the regex (pasted
URLs, report templates, prose) excluded as a collection artifact (see Residual errors,
below). The strings are whatever a clinician or curator pasted or typed, carrying
the damage of their route to the box: whitespace and
non-printable characters from Word documents and report PDFs, lost casing, transposed
punctuation, trailing protein annotations. The cleaning pipeline (`clean_hgvs()` plus
the provider-verified accession-prefix restoration, Methods) raised the fraction
parseable by biocommons/hgvs from {{ cleaning_corpus.as_submitted_pct | dp(1) }}%
as-submitted to {{ cleaning_corpus.after_pct | dp(1) }}%, a
+{{ cleaning_corpus.gain_pct | dp(1) }}% absolute gain
({{ cleaning_corpus.rescued | commas }} strings rescued, about
{{ cleaning_corpus.rescued_share_pct | int }}% of the
{{ cleaning_corpus.failed_pp | dp(1) }} percentage points that failed
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

**Table 2. Fixes applied across the production corpus (N = {{ cleaning_corpus.corpus_n | commas }}).** Each row is a
cleaning fix category, with the number of rescued queries in which it fired and its
share of the 1,721 rescued queries. Categories overlap (a single query may need several
fixes), so the counts sum to more than the total. *(Private data; counts are frozen constants
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
| **Total unique queries rescued** | | **{{ cleaning_corpus.rescued | commas }}** | **100%** |

As a control, `paper/scripts/inject_and_clean.py` injects each
`clean_hgvs()` fix category into a seeded sample of clean, parseable ClinVar c.HGVS
strings and confirms the cleaner recovers the canonical target with
**{{ cleaning.inject_regressions | int }} regressions**: no already-valid string is
ever broken (Supplementary Table S5). Because this benchmark injects the very errors it
then repairs, its recovery rate is reported in the supplement only; its purpose is the
no-regression guarantee on which the production result depends.

### Residual errors: the ceiling of cleaning *(Table S6)*

**[private data].** A small residue of non-HGVS input ({{ cleaning_corpus.nonhgvs_n | int }}
queries: pasted URLs, report templates, prose) that slipped the corpus regex is a
data-collection artifact, not something cleaning could ever repair, so it is excluded
from the corpus above. The remaining {{ cleaning_corpus.residual_pct | dp(1) }}%
({{ cleaning_corpus.residual_n | commas }} queries) of genuine HGVS-shaped input that
still fails to parse after cleaning defines the ceiling of pure string repair. Classified under a fixed
decision-tree taxonomy (Supplementary Table S6, with synthesised examples and the
classification method and its limitations), well over half is incomplete or
reference-less input that no string-level repair can invent, about a third is in
principle fixable and marks the frontier for future cleaning rules, and the remaining
~8% is valid HGVS the biocommons grammar rejects.

## R5: Transcript version substitution and coordinate safety

When a cited transcript version is absent from the loaded data, the opt-in substitution
step (Methods) supplies an adjacent version only if the coordinate-safety check passes;
a substitution that cannot be verified safe is refused by default, preserving
exact-version semantics. This is a client-layer feature: biocommons/hgvs has no
adjacent-version substitution with any data provider. In an end-to-end ablation
(`paper/scripts/benchmark_resolution.py`) that removes the requested version from each
test variant, the substitution recovered the correct genomic coordinate with no false
rescues (a false rescue being a substitution that resolves to a different coordinate).

Caution is warranted because substituting one version for another can change the
coordinate a variant projects to. Across consecutive
RefSeq version pairs ({{ version_stability.refseq_pairs | commas }} pairs),
{{ version_stability.refseq_preserving_pct | dp(1) }}% preserved every coding
coordinate; for Ensembl ({{ version_stability.ensembl_pairs | commas }} pairs)
{{ version_stability.ensembl_preserving_pct | dp(1) }}%. Weighted by coding base (the
chance a *random* variant is unaffected), safety is
{{ (version_stability.refseq_pervariant_safety * 100) | dp(1) }}% (RefSeq) and
{{ (version_stability.ensembl_pervariant_safety * 100) | dp(1) }}% (Ensembl). When a
coordinate does move it is almost always the whole CDS, driven by a re-annotation of
the coding region; the most dangerous case, a *partial* substitution that mis-places some
variants but not others, is rare
({{ version_stability.refseq_partial_drift_pct | dp(1) }}% of RefSeq pairs,
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
