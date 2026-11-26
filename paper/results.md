# Results

## R1: Transcript coverage

The following counts are of a transcript version aligned to a particular genome build.
cdot has {{ coverage.total_count | commas }} alignments compared to
{{ coverage.uta_count | commas }} in the 2024 UTA release [@Hart2020],
{{ coverage.improvement_fold | dp(1) }} times as many. This gain comes from ingesting every
archived RefSeq release, and also from bringing in Ensembl transcripts which UTA lacks.
Ensembl GRCh38 alone accounts for {{ coverage.ensembl_unique_count | commas }} transcript
versions. Adding T2T brings in {{ coverage.t2t_unique_count | commas }}.

## R2: ClinVar and clinical resolution accuracy

At full scale, resolving every c.HGVS in ClinVar's variant_summary through cdot alone
({{ clinvar_vcf.n_pairs | commas }} (g., c.) pairs, effectively all RefSeq) reaches
{{ clinvar_vcf.resolved_pct | pct(1) }} resolution, and
{{ clinvar_vcf.matched_of_resolved_pct | pct(2) }} of resolved variants reproduce
ClinVar's own VCF coordinate exactly (Methods). The residual
{{ clinvar_vcf.incorrect_pct | pct(2) }} is dominated by indel-representation
differences and by paralog and copy-number transcripts that map to more than one
genomic locus (Supplementary Table S4). UTA stores no Ensembl alignments at all, so it
resolves {{ clinvar.uta_ensembl_pct | pct(1) }} of the
{{ clinvar.n_ensembl | commas }} Ensembl ClinVar [@Landrum2025] pairs against cdot's
{{ clinvar.cdot_ensembl_pct | pct(1) }} (Supplementary Results).

The corpus above takes its c.HGVS from the `Name` column of ClinVar's
variant_summary, ClinVar's own recomputed name, always at the current transcript
version, so it measures a current-version ceiling. To read ClinVar's submission
history as a record of the transcript versions labs actually used, we built a second
public corpus from the per-submission (SCV) HGVS attributes of the ClinVar VCV XML:
each submitted string kept verbatim and joined to the variant's VCF coordinate as
ground truth via its AlleleID
({{ clinvar_submitted.n_unique_pairs | commas }} unique submitted-string/variant pairs
from {{ clinvar_submitted.n_scv_tx_strings | commas }} SCV transcript expressions;
Methods). Its version profile confirms that submitted traffic is historical:
{{ clinvar_submitted.version_not_current_pct | pct(1) }} of submitted strings cite a
transcript version that is no longer the version in the current RefSeq annotation
release. cdot's merged release history holds
{{ clinvar_submitted.not_current_in_cdot_pct | pct(1) }} of those superseded versions.
cdot resolved the full set of {{ clinvar_submitted.full_n | commas }} submitted pairs,
reproducing ClinVar's VCF coordinate for
{{ clinvar_submitted.cdot_resolved_pct | pct(1) }}, versus
{{ clinvar_submitted.uta_resolved_pct | pct(1) }} for the same locally loaded UTA
(release `{{ clinvar_submitted._prov.tools.uta }}`). **[private data].** The same gap
holds on the historical clinical data that motivated cdot: the complete set of
{{ historical.n_lines | commas }} unique HGVS descriptions imported into the Shariant
platform [@Tudini2022], classifications submitted by clinical laboratories over many
years ({{ historical.n_multi_version | commas }} of the
{{ historical.n_unique_tx | commas }} distinct transcripts are cited at more than one
version). With only the transcript-data layer swapped (Supplementary Methods), cdot
produced a genomic coordinate for {{ historical.cdot_resolved_pct | pct(1) }} versus
{{ historical.uta_resolved_pct | pct(1) }} for UTA. There is no ground-truth genomic
coordinate for the private corpus, so the metric is resolution rate rather than
correctness.

## R3: Throughput

**Table 1. End-to-end HGVS resolution throughput by transcript backend**: median (IQR)
HGVS/s over {{ benchmark.n_repeats | int }} timed passes of the identical
{{ benchmark.n_pairs | commas }}-pair ClinVar set, engine and sequence layer held
constant (Methods). The public remote UTA row uses the first
{{ benchmark.uta_remote_n | int }} pairs of the same set.

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
{{ benchmark_fullscale.resolved_pct | pct(1) }} producing a genomic coordinate. The same
exhaustive pass against the public remote UTA database extrapolates to
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
failed as-submitted) with {{ cleaning_corpus.regressions | number_word }} regressions.
Table 2 breaks the rescues down by fix type. The rarest category,
restoring a fully dropped accession prefix against the loaded transcript data, is
impossible for a purely string-level checker; having the transcript data locally is
what makes it safe.

**Table 2. Fixes applied across the production corpus (N = {{ cleaning_corpus.corpus_n | commas }}).** Each row is a
cleaning fix category, with the number of rescued queries in which it fired and its
share of the {{ cleaning_corpus.rescued | commas }} rescued queries. Categories overlap
(a single query may need several fixes), so the counts sum to more than the total.
*(Private data; frozen constants.)*

<!-- include-csv: empirical_results/cleaning_fix_categories.csv
  {align: left, na_rep: ""}
-->

The LOVD checker is the closest comparator (offline, no reference sequence), so we ran
it head-to-head with `clean_hgvs()` on the injection corpus
({{ lovd_comparison.n_cases | commas }} corrupted strings, checker
{{ lovd_comparison.lovd_version }}; scoring in Methods, per-category breakdown in
Supplementary Table S5). The corpus injects the error classes `clean_hgvs()` targets,
so `clean_hgvs()` recovers every case by construction and the comparison measures how
much of that territory a general syntax checker also covers. Weighted by the production
error mix, LOVD's top-ranked correction restored
{{ lovd_comparison.lovd_top1_weighted_pct | pct(0) }} of cases
({{ lovd_comparison.lovd_top1_pct | pct(0) }} unweighted), matching `clean_hgvs()` on
common single-token errors but not on free-text damage or cross-accession-family
repairs; neither tool altered any valid input.
The same corpus puts the sequence-aware services in context. Run over their public REST
APIs ({{ vv_mutalyzer_comparison.service_date }}; Methods), VariantValidator's validated
description recovered {{ vv_mutalyzer_comparison.vv_weighted_pct | pct(0) }} of cases
weighted by the production error mix, and Mutalyzer
{{ vv_mutalyzer_comparison.mut_weighted_pct | pct(0) }} (Supplementary Table S5). Each
repairs the damage it can interpret and rejects the rest as invalid; neither falsely
corrected a valid input.

## R5: Transcript version substitution and coordinate safety

Adjacent-version substitution is a client-layer feature: biocommons/hgvs has none with
any data provider. Caution is warranted because substituting one version for another
can change the coordinate a variant projects to. Across consecutive RefSeq version pairs
({{ version_stability.refseq_pairs | commas }} pairs),
{{ version_stability.refseq_preserving_pct | pct(1) }} preserved every coding
coordinate; for Ensembl ({{ version_stability.ensembl_pairs | commas }} pairs)
{{ version_stability.ensembl_preserving_pct | pct(1) }} (Supplementary Results,
Supplementary Figure S1).
The intrinsic-structure check (Methods) works even when the requested version was never
aligned to the target build, because the structure is build-independent, identical
across GRCh37/GRCh38/T2T for
{{ version_stability.refseq_struct_portable_pct | pct(1) }} of RefSeq and
{{ version_stability.ensembl_struct_portable_pct | pct(1) }} of Ensembl versions
(predictor detail in Supplementary Methods). Validated against ClinVar: of the
{{ version_safety_validation.n_safe_substitutions | commas }} version substitutions the
safety checks accept across the full RefSeq corpus, every one that resolves under the
substitute ({{ version_safety_validation.n_resolved_under_w | commas }}) gives the same
genomic coordinate as the cited version
({{ version_safety_validation.n_coordinate_changes }} coordinate changes).
