# Results

*Full Original Paper (Bioinformatics) — Results section. New for the full-paper switch.*
*Numbers come from the two-tier fact model (plan §3): Tier 1 = reproducible `{{ }}` facts; Tier 2 =
frozen aggregate constants from the private corpus, flagged non-reproducible. Synthesised examples
only — never corpus strings.*

> **STATUS: STUB.** Subsection skeleton below; prose written in the implementation step. Some prose
> migrates from the old `implementation.md` "Coverage and ClinVar benchmark" subsection.

---

## R1 — Transcript coverage and ClinVar resolution  *(Figure 1)*

*[Tier 1]* Coverage counts (cdot RefSeq + Ensembl vs UTA, by build) and ClinVar c→g resolution
(% correct / no_data, cdot vs UTA). Migrate + expand the old "Coverage and ClinVar benchmark" prose.
Facts: `coverage.*`, `clinvar.*`. Scale beyond the easy 100-set (plan §5.4).

## R2 — String cleaning recovers malformed real-world HGVS  *(Figure 2)*

*[Tier 1 — headline]* Injection benchmark (`analysis/inject_and_clean.py`): perturb public ClinVar
c.HGVS with each `clean.py` fix category at weights informed by the aggregate real distribution;
report recovery % per class and **0 regressions**. Fact: `cleaning.*`.

*[Tier 2 — production validation, non-reproducible]* Real corpus: N = 32,752 production queries,
parseable **91.5% → 96.4% (+4.9%), 0 regressions**; rescue-op distribution. Constants from
`cdot_private/output/` — flag as aggregate, not reproducible.

## R3 — The ceiling of cleaning: a taxonomy of residual errors  *(Figure 2 / Table S6)*

*[Tier 2]* 8-class taxonomy of what remains unparseable after cleaning (TRUNCATED 24.2%,
MISSING_REFERENCE 23.6%, MALFORMED_ACCESSION 14.2%, …). Frame the ceiling: ~48% irreducible
incomplete/reference-less user input; ~31% further-fixable; 6.9% a biocommons grammar gap; 6.9%
non-HGVS junk. Disclose LLM single-rater method (plan §6).

## R4 — Transcript version fallback  *(Table S?)*

*[Tier 1]* `get_best_transcript_version()`: drop the requested version → rescue % + **false-rescue
rate**. Unique-to-cdot capability (needs version depth). `benchmark_plan.md` Table 3.2.

## R5 — Throughput  *(Supplementary Figure S1)*

*[Tier 1]* Load time, throughput (HGVS/s, tx/s), prefetch. Local JSON vs REST vs UTA. Facts already
in `output/facts/resolution_*.json`. REST+prefetch 451 HGVS/s; sequence-fetch is the bottleneck.

## R6 — Alignment-gap correctness & T2T uniqueness

*[Tier 1]* Münz BRCA2-class gap experiment (32%→~100% with gap encoding); count of transcripts/HGVS
resolvable only on CHM13v2.0. `benchmark_plan.md` Tables 4–5.

---

*Notes: keep Tier-1/Tier-2 provenance explicit in every reported number. Audit "incorrect" buckets
before reporting (alignment diffs #95). Normalise before string-compare.*
