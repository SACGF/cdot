# Results

## R1: Transcript coverage

The unit of coverage is the *transcript-version alignment*: a particular transcript
version aligned to a particular genome build, counted separately per build. cdot covers
{{ coverage.total_count | commas }} such alignments across all builds and sources,
compared with ~{{ literature.uta_count | commas }} in UTA's 2021 release [@Hart2020],
{{ coverage.improvement_fold | fold(1) | article }} increase. The gain has two sources:
{% assert clinvar.uta_ensembl_pct == 0 and coverage.ensembl_unique_count > 0, "'Ensembl itself, absent from UTA entirely'" %}
historical depth (cdot ingests the complete run of archived RefSeq and Ensembl releases,
so an NM_ version cited in an older report or ClinVar submission still resolves after
NCBI has retired it from the current annotation) and Ensembl itself, absent from UTA
entirely ({{ coverage.ensembl_unique_count | commas }} versioned Ensembl transcripts on
GRCh38 alone).
T2T-CHM13v2.0 adds a further {{ coverage.t2t_unique_count | commas }} alignments,
making cdot the first transcript data source to bring that assembly to the Python HGVS
libraries. The JSON format also stores per-exon alignment-gap information (indels of
the transcript relative to the genome) that downstream libraries apply during
coordinate conversion; Ensembl VEP, which can itself generate HGVS, does not apply these
gaps when converting a c. description to genomic coordinates [@VepHgvsGaps], so it can
misplace variants downstream of an indel in a gapped RefSeq alignment.

## R2: ClinVar and clinical resolution accuracy

Resolution rate is compared against a locally loaded UTA (release `{{ clinvar_submitted._prov.tools.uta }}`),
first on current-version descriptions and then on the historical submitted-string
corpus that is the real test of transcript depth.

For the current-version check, the committed ClinVar [@Landrum2025] pair sets
({{ clinvar.n_refseq | commas }} RefSeq NM_ and {{ clinvar.n_ensembl | commas }} Ensembl
ENST preferred-transcript descriptions) were resolved through cdot, with transcript
sequence spliced from a local genome FASTA so that no Ensembl pair is dropped for a
sequence SeqRepo lacks. {# lint-ok: 'no Ensembl pair is dropped' states why FastaSeqFetcher is used #}
UTA's rate is membership of the exact cited accession.version
in its GRCh38 alignment set, which is what resolution through UTA requires; the live
{% assert clinvar.cdot_refseq_pct == clinvar.uta_refseq_pct, "'On RefSeq the two backends are at parity'" %}
database gives the same rates. On RefSeq the two backends are at parity (cdot
{{ clinvar.cdot_refseq_pct | pct(1) }}, UTA {{ clinvar.uta_refseq_pct | pct(1) }}). The
{% assert clinvar.uta_ensembl_pct == 0, "'UTA stores no Ensembl alignments at all'" %}
remaining difference is capability, not accuracy: UTA stores no Ensembl alignments at
all, so it resolves {{ clinvar.uta_ensembl_pct | pct(1) }} of the Ensembl pairs against
cdot's {{ clinvar.cdot_ensembl_pct | pct(1) }}. Because the pair sets are
Ensembl-enriched, their aggregate rate is not a like-for-like accuracy score and we do
not report one.

At full scale, resolving every c.HGVS in ClinVar's variant_summary through cdot alone
({{ clinvar_vcf.n_pairs | commas }} (g., c.) pairs, effectively all RefSeq) reaches {# lint-ok: 'every c.HGVS' describes the corpus; clinvar_source_breakdown.csv has no Ensembl rows #}
{{ clinvar_vcf.resolved_pct | pct(1) }} resolution, and
{{ clinvar_vcf.matched_of_resolved_pct | pct(2) }} of resolved variants reproduce
{% assert clinvar_vcf_residual.indel_mismatch * 2 > clinvar_vcf_residual.n_incorrect and clinvar_vcf_residual.snv_diff_pos > clinvar_vcf_residual.unaligned_between_exons, "'dominated by indel-representation differences and by paralog and copy-number transcripts ..., plus N transcripts with bases that align to neither'" %}
ClinVar's own VCF coordinate exactly (Methods). The residual
{{ clinvar_vcf.incorrect_pct | pct(2) }} is dominated by indel-representation
differences and by paralog and copy-number transcripts that map to more than one
genomic locus, plus
{{ clinvar_vcf_residual.unaligned_between_exons_transcripts | number_word }} transcripts
with bases that align to neither of two adjacent exons, where cdot mostly gives the same
coordinate as UTA and ClinVar places the unaligned bases differently (Supplementary Table S4).

### ClinVar submissions as a historical record of transcripts used

Both corpora above take their c.HGVS from the `Name` column of ClinVar's
variant_summary, ClinVar's own recomputed name, always at the current transcript
version, so they measure a current-version ceiling. To read ClinVar's submission
history as a record of the transcript versions labs actually used, we built a second
public corpus from the per-submission (SCV) HGVS attributes of the ClinVar VCV XML:
each submitted string kept verbatim and joined to the variant's VCF coordinate as
ground truth via its AlleleID
({{ clinvar_submitted.n_unique_pairs | commas }} unique submitted-string/variant pairs
from {{ clinvar_submitted.n_scv_tx_strings | commas }} SCV transcript expressions;
{% assert clinvar_submitted.ensembl_pct == 0, "'not one submitted string cites an Ensembl transcript'" %}
Methods). This corpus reflects the natural submission mix, which is essentially all
RefSeq (not one submitted string cites an Ensembl transcript; Ensembl accessions enter
ClinVar through its recomputed names, not submitter descriptions). Its version profile
confirms that submitted traffic is historical:
{% assert clinvar_submitted.version_not_current_pct > 50 and clinvar_submitted.absent_cdot_pct < 1, "'submitted traffic is historical' (most cite a superseded version); 'only X% of cited versions are absent'" %}
{{ clinvar_submitted.version_not_current_pct | pct(1) }} of submitted strings cite a
transcript version that is no longer the version in the current RefSeq annotation
release ({{ clinvar_submitted.scv_weighted_not_current_pct | pct(1) }} weighted by
submission count), and {{ clinvar_submitted.base_retired_pct | pct(1) }} cite a
transcript no longer annotated at any version. cdot's merged release history holds
{{ clinvar_submitted.not_current_in_cdot_pct | pct(1) }} of those superseded versions;
only {{ clinvar_submitted.absent_cdot_pct | pct(1) }} of cited versions are absent from
its GRCh38 data.

{% assert submitter_attribution.n_submitters_with_absent * 10 < submitter_attribution.n_submitters_total, "'only N submitters ever cite a version cdot lacks' (under a tenth)" %}
Because cdot ingests every RefSeq and Ensembl annotation release, a version absent from
its data should be one that never appeared in an official release, and the absent
citations bear that out: of {{ submitter_attribution.n_submitters_total | commas }}
submitters only {{ submitter_attribution.n_submitters_with_absent | int }} ever cite a
version cdot lacks, a single large laboratory accounts for
{{ submitter_attribution.top1_absent_share_pct | pct(1) }} of all absent-version
citations, and {{ submitter_attribution.single_submitter_version_pct | pct(1) }} of the
{{ submitter_attribution.n_absent_versions | int }} distinct absent versions are cited by
one submitter alone.
These are the signature of transcripts a laboratory aligned to the
genome itself rather than drawing from a published release, so they are addressed by
version substitution (R5), not by deeper ingest.

cdot resolved the full set of {{ clinvar_submitted.full_n | commas }} submitted pairs,
reproducing ClinVar's VCF coordinate for
{% assert clinvar_submitted.cdot_resolved_pct > clinvar_submitted.uta_resolved_pct, "cdot X% 'versus' UTA Y% on the submitted pairs" %}
{{ clinvar_submitted.cdot_resolved_pct | pct(1) }}, versus
{% assert clinvar.cdot_refseq_pct == clinvar.uta_refseq_pct, "'The RefSeq gap invisible at the current-version ceiling (both backends ...)'" %}
{{ clinvar_submitted.uta_resolved_pct | pct(1) }} for the same locally loaded UTA. The
RefSeq gap invisible at the current-version ceiling (both backends
{{ clinvar.cdot_refseq_pct | pct(1) }} above) opens to
{{ clinvar_submitted.cdot_only_pct | dp(1) }} points:
{{ clinvar_submitted.cdot_only | commas }} pairs resolve through cdot alone and
{{ clinvar_submitted.uta_only | commas }} through UTA alone, because UTA holds no GRCh38
alignment for {{ clinvar_submitted.uta_no_data_pct | pct(1) }} of the cited versions.
{% assert clinvar_submitted.era_old_uta_pct < clinvar_submitted.era_recent_uta_pct < clinvar_submitted.era_mid_uta_pct, "UTA by era: 'lowest on the oldest submissions ... and lower on the newest than on 2016-2020'" %}
Splitting by submission era locates the gap: cdot reproduces the coordinate for
{{ clinvar_submitted.era_old_cdot_pct | pct(1) }} of the oldest (2008-2015)
submissions, {{ clinvar_submitted.era_mid_cdot_pct | pct(1) }} of 2016-2020 and
{{ clinvar_submitted.era_recent_cdot_pct | pct(1) }} of 2021-2026, while UTA reaches
{{ clinvar_submitted.era_old_uta_pct | pct(1) }},
{{ clinvar_submitted.era_mid_uta_pct | pct(1) }} and
{{ clinvar_submitted.era_recent_uta_pct | pct(1) }}: lowest on the oldest submissions,
which cite the most superseded versions, and lower on the newest than on 2016-2020.

The submitted strings are largely well-formed, so cleaning has little to rescue
(`fix_hgvs()` repaired {{ clinvar_submitted.rescued_by_fix | commas }} strings, with
{% assert clinvar_submitted_residual.version_refused + clinvar_submitted_residual.coordinate_drift > clinvar_submitted_residual.n_residual / 2, "'The residual ... is dominated by historical-version effects'" %}
{{ clinvar_submitted.regressions | int }} regressions). The residual
{{ clinvar_submitted.residual_pct | pct(1) }} after cleaning and version substitution
({{ clinvar_submitted.residual_n | commas }} of {{ clinvar_submitted.full_n | commas }})
is dominated by historical-version effects:
{{ clinvar_submitted_residual.version_refused | commas }} cite a version absent from the
data where substitution declines to act because coordinate safety cannot be verified,
and {{ clinvar_submitted_residual.coordinate_drift | commas }} resolve through the cited
version to a coordinate that differs from ClinVar's current interpretation; the remainder
cite a position or base that does not exist on the cited version, notation the
biocommons grammar rejects, or an accession absent entirely (Supplementary Table S8).

**[private data].** The same gap holds on the historical clinical data that motivated cdot:
the complete set of {{ historical.n_lines | commas }} unique HGVS descriptions imported
into the Shariant platform [@Tudini2022], classifications submitted by clinical
laboratories over many years ({{ historical.n_multi_version | commas }} of the
{{ historical.n_unique_tx | commas }} distinct transcripts are cited at more than one
version). With only the transcript-data layer swapped (Supplementary Methods), cdot
{% assert historical.cdot_resolved_pct > historical.uta_resolved_pct, "historical corpus: cdot X% 'versus' UTA Y%" %}
produced a genomic coordinate for {{ historical.cdot_resolved_pct | pct(1) }} versus
{{ historical.uta_resolved_pct | pct(1) }} for UTA. Of the strings cdot resolved but
UTA could not ({{ historical.cdot_only_pct | pct(1) }} of the corpus),
{{ historical.cdot_only_historical_pct | pct(0) }} were RefSeq transcript versions for
which UTA holds no GRCh38 alignment and
{{ historical.cdot_only_ensembl_pct | pct(0) }} were Ensembl transcripts. There is no
ground-truth genomic coordinate for the private corpus, so the metric is resolution
rate rather than correctness. {# lint-ok(over, more than, only, no): method wording: 'over many years', 'more than one version', 'only the transcript-data layer swapped', 'holds no GRCh38 alignment', 'no ground truth' #}

## R3: Throughput

**Table 1. End-to-end HGVS resolution throughput by transcript backend**: median (IQR)
HGVS/s over {{ benchmark.n_repeats | int }} timed passes of the identical
{{ benchmark.n_pairs | commas }}-pair ClinVar set, engine and sequence layer held
constant (Methods). The public remote UTA row uses the first
{{ benchmark.uta_remote_n | int }} pairs of the same set. {# lint-ok(over): 'over' is a preposition #}

| Configuration | Throughput (HGVS/s), median (IQR) |
|---|---|
| UTA: public remote database | {{ benchmark.uta_remote_tps | dp(2) }} ({{ benchmark.uta_remote_tps_q1 | dp(2) }}–{{ benchmark.uta_remote_tps_q3 | dp(2) }}) |
| UTA: local PostgreSQL | {{ benchmark.uta_local_tps | int }} ({{ benchmark.uta_local_tps_q1 | int }}–{{ benchmark.uta_local_tps_q3 | int }}) |
| cdot REST (one request per transcript) | {{ benchmark.cdot_rest_tps | int }} ({{ benchmark.cdot_rest_tps_q1 | int }}–{{ benchmark.cdot_rest_tps_q3 | int }}) |
| cdot REST (after one batch `prefetch()`) | {{ benchmark.cdot_rest_prefetch_tps | int }} ({{ benchmark.cdot_rest_prefetch_tps_q1 | int }}–{{ benchmark.cdot_rest_prefetch_tps_q3 | int }}) |
| cdot local JSON | {{ benchmark.cdot_local_tps | int }} ({{ benchmark.cdot_local_tps_q1 | int }}–{{ benchmark.cdot_local_tps_q3 | int }}) |

A GRCh38 RefSeq JSON file loads in ~{{ benchmark.grch38_load_time_s | dp(0) }} s.
Batching the per-transcript REST lookups into one `prefetch()` request (all transcripts
for the set, under a second, untimed) makes REST throughput equivalent to local JSON,
the two {{ benchmark.rest_prefetch_vs_local_diff_pct | within_pct }} of each other:
with the transcript data in process memory, both are bounded by the shared engine and sequence layer, not by the transcript backend.
{% assert 3.6 <= benchmark.cdot_local_vs_uta_local_ratio <= 4.4, "local UTA 'about a quarter of local-JSON throughput'" %}
Comparing like with like: locally, a loaded UTA reached
{{ benchmark.uta_local_tps | int }} HGVS/s, about a quarter of local-JSON throughput,
each lookup resolving a set of SQL queries where cdot needs a single JSON object; and
remote to remote, cdot's REST API (one request per transcript, no prefetch) is
{{ benchmark.cdot_rest_vs_uta_remote_ratio | orders_of_magnitude }} faster than the
public UTA server, which pays wide-area
round trips to a shared database on every lookup.

At scale, a single local-JSON process resolved the full set of
{{ benchmark_fullscale.n_pairs | commas }} unique ClinVar (g.HGVS, c.HGVS) pairs in
~{{ benchmark_fullscale.wall_min | int }} minutes,
{{ benchmark_fullscale.resolved_pct | pct(1) }} producing a genomic coordinate. Over
{{ benchmark_fullscale.n_runs | number_word }} hot-cache repeats it averaged
{{ benchmark_fullscale.hot_tps | int }} HGVS/s
({{ benchmark_fullscale.hot_tps_min | int }}–{{ benchmark_fullscale.hot_tps_max | int }}),
{{ benchmark_fullscale.hot_vs_cold_diff_pct | within_pct }} of the cold pass ({{ benchmark_fullscale.cold_tps | int }} HGVS/s), so the state of the
sequence cache is not the bound either. The same exhaustive pass against the public
remote UTA database extrapolates to
{{ benchmark_fullscale.uta_remote_days | approx(per=365, noun='year') }}.

## R4: String cleaning recovers malformed real-world HGVS

**[private data].** The main test of cleaning is a production query stream:
N = {{ cleaning_corpus.corpus_n | commas }}
real search-box queries from clinical and research variant-curation platforms based on
VariantGrid [@VariantGrid], restricted to those matching a broad HGVS regex (loosely
HGVS-shaped strings, not arbitrary free-text search terms). The strings are whatever a
clinician or curator pasted or typed, carrying
the damage of their route to the box: whitespace and non-printable characters from
Word documents and report PDFs, lost casing, transposed punctuation, trailing protein
annotations. The cleaning pipeline (`clean_hgvs()` plus the provider-verified
accession-prefix restoration, Methods) raised the fraction parseable by
biocommons/hgvs from {{ cleaning_corpus.as_submitted_pct | pct(1) }} as-submitted to
{{ cleaning_corpus.after_pct | pct(1) }} ({{ cleaning_corpus.rescued | commas }}
strings rescued, about {{ cleaning_corpus.rescued_share_pct | pct(0) }} of those that
{% assert cleaning_corpus.regressions == 0, "'with zero regressions'" %}
failed as-submitted) with {{ cleaning_corpus.regressions | number_word }} regressions.
Table 2 breaks the rescues down by fix
type: structural-punctuation repair and whitespace removal dominate, followed by
gene/transcript-wrapper repair and structure reconstruction. The rarest category,
restoring a fully dropped accession prefix against the loaded transcript data, is
impossible for a purely string-level checker; having the transcript data locally is
what makes it safe.

The public submitted-string corpus (R2) is the complement, and the two corpora locate
the two failure axes: formal database submissions are largely well-formed and fail by
transcript-version age, whereas interactive human-typed input carries the formatting
damage cleaning repairs.

**Table 2. Fixes applied across the production corpus (N = {{ cleaning_corpus.corpus_n | commas }}).** Each row is a
cleaning fix category, with the number of rescued queries in which it fired and its
share of the {{ cleaning_corpus.rescued | commas }} rescued queries. Categories overlap
(a single query may need several fixes), so the counts sum to more than the total.
*(Private data; frozen constants.)*

<!-- include-csv: empirical_results/cleaning_fix_categories.csv
  {align: left, na_rep: ""}
-->

As a control, `paper/scripts/inject_and_clean.py` injects each `clean_hgvs()` fix
category into a seeded sample of clean, parseable ClinVar c.HGVS strings and confirms
the cleaner recovers the canonical target with
{{ cleaning.inject_regressions | int }} regressions (Supplementary Table S5).
Because this benchmark injects the very errors it then repairs, its recovery rate is
reported in the supplement only; its purpose is the no-regression guarantee on which
the production result depends.

### Residual errors: the ceiling of cleaning *(Table S6)*

**[private data].** A small residue of non-HGVS input ({{ cleaning_corpus.nonhgvs_n | int }}
queries: pasted URLs, report templates, prose) that slipped the corpus regex is a
data-collection artifact and is excluded from the corpus above. The remaining
{{ cleaning_corpus.residual_pct | pct(1) }} ({{ cleaning_corpus.residual_n | commas }}
queries) of genuine HGVS-shaped input that still fails to parse after cleaning defines
the ceiling of pure string repair. Classified under a fixed decision-tree taxonomy
(Supplementary Table S6), well over half is incomplete or reference-less input that no
{% assert cleaning_corpus.residual_incomplete_pct > 55 and 20 < cleaning_corpus.residual_fixable_pct < 30 and (cleaning_corpus.residual_incomplete_pct + cleaning_corpus.residual_fixable_pct + cleaning_corpus.residual_grammar_gap_pct - 100) | abs < 0.5, "Table S6 residual: 'well over half' incomplete, 'about a quarter' fixable, 'the remaining' grammar gap" %}
string-level repair can invent, about a quarter is in principle fixable and marks the
frontier for future cleaning rules, and the remaining
~{{ cleaning_corpus.residual_grammar_gap_pct | pct(0) }} is valid HGVS the biocommons
grammar rejects.

## R5: Transcript version substitution and coordinate safety

Adjacent-version substitution is a client-layer feature: biocommons/hgvs has none with
any data provider.

Caution is warranted because substituting one version for another can change the
coordinate a variant projects to. Across consecutive RefSeq version pairs
({{ version_stability.refseq_pairs | commas }} pairs),
{{ version_stability.refseq_preserving_pct | pct(1) }} preserved every coding {# lint-ok: 'preserved every coding coordinate' defines the preserving category #}
coordinate; for Ensembl ({{ version_stability.ensembl_pairs | commas }} pairs)
{{ version_stability.ensembl_preserving_pct | pct(1) }}. Weighted by coding base (the
chance a *random* variant is unaffected), safety is
{{ version_stability.refseq_pervariant_safety | pct(1) }} (RefSeq) and
{{ version_stability.ensembl_pervariant_safety | pct(1) }} (Ensembl).
{% assert version_stability.refseq_full_drift_pct > 2 * version_stability.refseq_partial_drift_pct and version_stability.ensembl_full_drift_pct > 2 * version_stability.ensembl_partial_drift_pct, "'When a coordinate does move it is usually the whole CDS'" %}
When a coordinate does move it is usually the whole CDS
({{ version_stability.refseq_full_drift_pct | pct(1) }} of RefSeq pairs,
{{ version_stability.ensembl_full_drift_pct | pct(1) }} Ensembl), driven by a
re-annotation of the coding region; the most dangerous case, a *partial* substitution
that mis-places some variants but not others, is rare
({{ version_stability.refseq_partial_drift_pct | pct(1) }} of RefSeq pairs,
{{ version_stability.ensembl_partial_drift_pct | pct(1) }} Ensembl). Within that tail
the risk is positional: a partial drift keeps a 5' prefix intact up to its first
{% assert positional_drift.refseq_partial_decile1_pct > positional_drift.refseq_partial_decile10_pct and positional_drift.ensembl_partial_decile1_pct > positional_drift.ensembl_partial_decile10_pct, "'preservation falls from the 5'-most decile to the 3'-most'" %}
alignment change and moves what lies downstream, so preservation falls from
{{ positional_drift.refseq_partial_decile1_pct | pct(0) }} of coding bases in the
5'-most decile to {{ positional_drift.refseq_partial_decile10_pct | pct(0) }} in the
3'-most for RefSeq ({{ positional_drift.ensembl_partial_decile1_pct | pct(0) }} to
{{ positional_drift.ensembl_partial_decile10_pct | pct(0) }} for Ensembl; Supplementary
Figure S1). Drift is also concentrated: only
{% assert version_stability.refseq_accessions_drift_pct < 10, "'only X% of RefSeq accessions ever drift ... a short, identifiable list'" %}
{{ version_stability.refseq_accessions_drift_pct | pct(1) }} of RefSeq accessions ever
drift, so the risky transcripts form a short, identifiable list that can be flagged or
blocklisted.

The intrinsic-structure check (Methods) works even when the requested version was never {# lint-ok: 'never aligned' describes the case the check handles #}
aligned to the target build, because the structure is build-independent, identical
across GRCh37/GRCh38/T2T for
{{ version_stability.refseq_struct_portable_pct | pct(1) }} of RefSeq and
{{ version_stability.ensembl_struct_portable_pct | pct(1) }} of Ensembl versions
{% assert version_safety_validation.n_coordinate_changes == 0, "'every one that resolves under the substitute gives the same genomic coordinate'" %}
(predictor detail in Supplementary Methods). Validated against ClinVar: of the
{{ version_safety_validation.n_safe_substitutions | commas }} version substitutions the
safety checks accept across the full RefSeq corpus, every one that resolves under the
substitute ({{ version_safety_validation.n_resolved_under_w | commas }}) gives the same
genomic coordinate as the cited version
({{ version_safety_validation.n_coordinate_changes }} coordinate changes).
