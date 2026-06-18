# Results

*Full Original Paper (Bioinformatics): Results section. New for the full-paper switch.*
*Numbers come from the two-tier fact model (plan §3): Tier 1 = reproducible templated
facts; Tier 2 =
frozen aggregate constants from the private corpus, flagged non-reproducible.
Synthesised examples
only; never corpus strings.*

> **Provenance flags.** Every reported number is tagged **[Tier 1]** (reproducible from
public data
> committed to this repo) or **[Tier 2]** (aggregate statistics from a private
production corpus,
> published as frozen constants, not reproducible by a referee; see Methods /
data-availability).

---

## R1: Transcript coverage and ClinVar resolution  *(Figure 1)*

**[Tier 1]** cdot covers {{ coverage.total_count | commas }} versioned transcript
alignments across all builds and sources, compared with
~{{ literature.uta_count | commas }} in UTA, a
{{ coverage.improvement_fold | fmt('.1f') }}× increase (Figure 1). The gain has two
distinct origins. First, multiple historical RefSeq releases are retained, so an NM_
version cited in an older clinical report or ClinVar submission still resolves even
after NCBI has retired it from the current annotation. Second, Ensembl is covered in
full: {{ coverage.ensembl_unique_count | commas }} Ensembl transcript accessions are
present in cdot but absent from UTA entirely. T2T-CHM13v2.0 contributes a further
{{ coverage.t2t_unique_count | commas }} transcripts with no GRCh37/38 equivalent, in
genomic regions that were inaccessible to short-read sequencing.

**[Tier 1]** To measure practical impact rather than raw counts, a seeded sample of
{{ clinvar.n_variants | commas }} ClinVar [@Landrum2025] variant descriptions,
spanning both RefSeq (NM_) and Ensembl (ENST) transcripts, was resolved against
cdot and a locally loaded UTA (release `uta_20241220`) through the identical
biocommons/hgvs code path (genomic→transcript projection, with sequences served from a
shared local SeqRepo so only the transcript data layer differs). cdot resolved
{{ clinvar.cdot_resolution_pct | dp(1) }}% versus
{{ clinvar.uta_resolution_pct | dp(1) }}% for UTA. On the RefSeq subset the two are at
parity (~99% each); the entire difference is Ensembl, which UTA cannot resolve at all
(0% of ENST inputs, as it stores no Ensembl alignments) whereas cdot resolves it
natively (Figure 1; Supplementary Table S4). The benchmark is reproducible: the ClinVar
build script, random seed, and resolution harness are committed.

## R2: String cleaning recovers malformed real-world HGVS  *(Figure 2)*

The cleaning contribution is evaluated at two levels: a reproducible injection benchmark
that isolates each fix category (Tier 1), and the aggregate behaviour of the same code
on a real production query stream (Tier 2).

**[Tier 1, headline] Injection benchmark.** `analysis/inject_and_clean.py` takes
{{ cleaning.inject_sample_size | commas }} clean, parseable ClinVar c.HGVS strings (a
seeded sample committed to this repo), and for each of the
{{ cleaning.inject_n_categories | int }} `clean_hgvs()` fix categories injects exactly
that error into each clean string, runs `clean_hgvs()`, and tests whether the result is
recovered. Because the biocommons grammar is lenient (it tolerates, for example,
lowercase bases), recovery is scored by the stronger criterion that the cleaned string
*exactly equals* the canonical target a correct cleaner should produce, with that target
independently confirmed to parse. Across {{ cleaning.inject_total_attempted | commas }}
injected cases spanning all {{ cleaning.inject_n_categories | int }} categories,
`clean_hgvs()` recovered {{ cleaning.inject_overall_pct | dp(1) }}%
({{ cleaning.inject_weighted_pct | dp(1) }}% when each category is weighted by its
real-world rescue frequency), with **{{ cleaning.inject_regressions | int }}
regressions**: no string that was already valid was broken by cleaning (Figure 2;
Supplementary Table S5). The per-category injection weights are constants citing the
aggregate production rescue-op distribution (Tier 2); no corpus string is read or copied
(plan §3). For the formatting-error classes that cleaning targets, the repair is
therefore deterministic and lossless by construction.

**[Tier 2: production validation, not reproducible].** Run over a corpus of **N =
32,752** real search-bar HGVS queries pooled from production clinical and research
variant-curation platforms, the same `clean_hgvs()` raised the fraction parseable by
biocommons/hgvs from **91.5%** as-submitted to **96.4%** after cleaning, a **+4.9%**
absolute gain (≈1,600 additional strings rescued) with **0 regressions**. The rescues
are dominated by whitespace removal and base re-casing, followed by protein-suffix
stripping and gene/transcript-wrapper repair, consistent with the failure modes cleaning
was designed for and with the injection weights above. These figures are aggregate
production statistics: the underlying queries contain real patient-derived data and are
not shareable, so, unlike the injection benchmark, they cannot be regenerated by a
referee. They are reported as frozen constants (source: `cdot_private/output/`, issue
#112).

## R3: The ceiling of cleaning, a taxonomy of residual errors  *(Figure 2 / Table S6)*

**[Tier 2: not reproducible].** The 3.6% of the production corpus (1,175 queries; 913
unique strings) that still fail to parse after cleaning define the ceiling of what
pure string repair can achieve. Each unique residual string was assigned to one of eight
single-label error classes under a fixed decision-tree taxonomy (Figure 2; Supplementary
Table S6):

| Class | Rows | What it is |
|---|---|---|
| TRUNCATED | 24.2% | cut off before a complete variant (ends mid-token or with an edit lacking sequence) |
| MISSING_REFERENCE | 23.6% | a bare variant body (`c.`/`n.`/`g.`/…) with no transcript, gene, or accession |
| MALFORMED_ACCESSION | 14.2% | broken reference: missing prefix, misplaced or truncated version |
| EDIT_SYNTAX_ERROR | 12.2% | malformed or non-standard edit operation (typo'd `del`/`dup`, multi-base substitution, count shorthand) |
| TRAILING_OR_CONCATENATED | 7.2% | a complete variant followed by extra characters, or several variants run together |
| NON_HGVS | 6.9% | not a variant attempt at all (URL, pasted report template, prose) |
| UNSUPPORTED_VALID_SYNTAX | 6.9% | legitimate HGVS the biocommons grammar rejects (uncertain ranges, `?` endpoints, allele/compound notation) |
| GENOMIC_REF_IN_PARENS | 4.9% | a genomic `NC_` accession placed in the parenthetical slot meant for a gene symbol |

The residual falls into four regimes. Roughly half (~48%: TRUNCATED
+ MISSING_REFERENCE) is incomplete or reference-less user input:
information the user never supplied, which no string-level repair can invent. About
38% (MALFORMED_ACCESSION + EDIT_SYNTAX_ERROR + TRAILING_OR_CONCATENATED +
GENOMIC_REF_IN_PARENS) is in principle fixable and marks the frontier for
future cleaning rules. The remaining residual is a 6.9% grammar gap, valid HGVS
that the parser rejects rather than the input (UNSUPPORTED_VALID_SYNTAX), plus 6.9%
non-HGVS junk that should not be parsed at all. Most of what remains is therefore
either out of scope for any repair or a downstream grammar limitation rather than a
cdot deficiency.

*Method and limitation (plan §6):* classification was performed by an LLM applying the
shared decision tree to each unique string, single-label and single-rater; no
second-rater adjudication was done, so no inter-rater agreement (κ) is reported. The
taxonomy is version `v1`. Synthesised examples (from public NM_000059.4 / NM_001754.5)
illustrate each class in the supplement; no corpus string is reproduced.

## R4: Transcript version fallback

**[Tier 1]** cdot exposes an opt-in adjacent-version fallback
(`get_best_transcript_version()`, surfaced through `fix_hgvs(...,
version_fallback=...)`): when a requested transcript version is absent from the loaded
data, the nearest available version is substituted under a configurable up-then-down,
closest, or latest policy, and the substitution is always reported as an `HGVSFix` so
the caller decides whether to accept it. In the end-to-end ablation
(`analysis/benchmark_resolution.py`), removing the requested version from each test
variant and resolving through the fallback recovered the correct genomic coordinate with
no false rescues (a false rescue being a substitution that resolves to a different
coordinate). The capability depends on cdot's multi-release version depth, since UTA's
single annotation snapshot has no adjacent versions to fall back to, and it is never
applied automatically, preserving exact-version semantics by default.

## R5: Throughput

**[Tier 1]** A GRCh38 RefSeq JSON file loads in
~{{ benchmark.grch38_load_time_s | dp(0) }} s and then resolves at
{{ benchmark.cdot_local_min_tps | int }}–{{ benchmark.cdot_local_max_tps | int }} HGVS/s
through the biocommons engine. The REST provider serves
~{{ benchmark.cdot_rest_tps | int }} HGVS/s on demand, but a single batch `prefetch()`
round-trip warms the transcript cache and lifts REST to match local throughput. To
compare the data layers on equal terms we held the sequence layer constant (a shared
local SeqRepo) and ran the identical engine over a locally loaded UTA: it reached only
~{{ benchmark.uta_local_tps | int }} HGVS/s, and over the public remote UTA database
~{{ benchmark.uta_remote_tps | dp(1) }} HGVS/s. cdot's local data layer is thus roughly
30× faster than a local UTA and orders of magnitude faster than the remote service most
users reach; once the data layer is local, the remaining bottleneck is sequence fetching
rather than transcript lookup.

**[Tier 1]** At scale, a single local-JSON process resolved the entire set of 3,660,452
unique ClinVar (g.HGVS, c.HGVS) pairs in ~92 minutes (665 HGVS/s; 99.3% produced a
genomic coordinate, 98.8% matched the ClinVar genomic HGVS exactly). The REST provider
matched this over the network: after one batch cache-warming pass (21,277 distinct
transcripts warmed in ~6 s) it resolved the same set in ~83 minutes (731 HGVS/s),
marginally faster than local JSON because warmed lookups hit a flat in-memory cache
instead of the interval-tree index. The same exhaustive pass is impractical against the
public remote UTA database (extrapolated at hundreds of days from its ~0.1 HGVS/s).
(`analysis/build_clinvar_pairs.py` builds the pair set by joining ClinVar's
variant_summary with the ClinVar VCF; the residual ~1% non-exact/error rate is dominated
by genomic-HGVS normalisation differences between the ClinVar source string and the
biocommons output, not resolution failures.)

## R6: Alignment-gap correctness and T2T uniqueness

**[Tier 1]** cdot encodes per-exon alignment gaps (indels of the transcript relative to
the genome) and applies them during coordinate conversion. This corrects the class of
error responsible for the {{ literature.brca2_accuracy_pct | dp(0) }}% BRCA2 accuracy
reported by @Munz2015 in tools lacking gap support: a correctness issue rather than a
coverage one, and the reason the PyHGVS integration adds gap handling that the upstream
library lacks. Separately, cdot is the first HGVS resource to cover T2T-CHM13v2.0,
contributing {{ coverage.t2t_unique_count | commas }} transcript alignments with no
GRCh37/38 equivalent and resolving HGVS in genomic regions that were inaccessible to
prior assemblies and absent from UTA.

---

*Notes: keep Tier-1/Tier-2 provenance explicit in every reported number. Audit
"incorrect" buckets
before reporting (alignment diffs #95). Normalise before string-compare.*
