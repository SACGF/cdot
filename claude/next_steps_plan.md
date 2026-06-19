# cdot — Next Steps Plan (Bioinformatics paper target)

Outstanding work only. Items that have shipped (code + paper Results + committed facts
CSVs) have been removed; see "Completed" at the bottom for the IDs, so nothing is lost.

## Conventions for analysis work
- **LLM classification/assignment tasks: emit SHORT class names** (e.g. `Truncated`,
  `No reference`, `Bad accession`, `Grammar gap`, `Genomic ref in parens`), not long
  ALL_CAPS taxonomy keys. Do it at the source so tables need no post-processing.
- **Never put private-corpus strings in this repo** (tests/comments/docstrings/paper).
  Synthesise from public examples (BRCA2 `NM_000059.4`, RUNX1 `NM_001754.5`).

---

## STRENGTHENS THE PAPER — non-blocking but important

### R7. Versionless transcript resolution in `get_best_transcript_version` (new issue)

**Why:** A client / code feature (changelog-worthy), not paper analysis, but motivated by
the same safety machinery. Today `resolve_transcript_version` no-ops on a versionless
accession (`_parse_versioned_transcript` returns `None`). Search-bar users routinely omit
the version (`NM_000059:c.…`). When there is only one available version, or all available
versions are structurally equivalent (same intrinsic CDS structure ⟹ coordinate-safe),
we can pick one with a WARNING.

**What to do:**
- Extend the cleaning/resolution path so a versionless transcript resolves to a concrete
  version when it is unambiguous (single version, or all versions structurally equal via
  `intrinsic_cds_structure`).
- **Ambiguity policy (decided): REFUSE.** When available versions are NOT structurally
  equivalent, return an ERROR fix and leave the string unchanged (mirrors the `REFUSE`
  default of `on_unsafe_version`). No silent pick.
- File a GitHub issue; add a `CHANGELOG.md` `[unreleased]` entry referencing it. Lands in
  `cdot/hgvs/clean.py` / `gene_hgvs.py`. Document under the planned `clean-hgvs.md` docs page (#117).

---

## Optional plumbing

- **Facts-driven R2/R3 tables.** R2 Table 1 and the R3 taxonomy table are hardcoded
  literal markdown (Tier-2 frozen, hand-transcribed — repo convention), NOT vibepaper
  facts. To make future refreshes a one-line CSV edit, promote them to a `residual.csv` /
  extended `cleaning.csv` + Snakefile rule. Optional.

---

## Summary Table

| ID | Task | Effort | Blocks submission |
|----|------|--------|-------------------|
| R7 | Versionless transcript resolution (REFUSE on ambiguity) | Medium | No — client feature + changelog |
| S2 | BRCA2 gap verification (demoted, post-paper) | Medium | No — dropped from paper framing |
| F3 | Automated data-release regeneration | Medium | No — backs a Discussion sentence |

---

## Completed

- **A1. Single ClinVar parse/resolution pass.** `paper/scripts/resolve_clinvar_pass.py`
  resolves the (g.HGVS, c.HGVS) pairs exactly once and writes the per-variant results
  table `g_hgvs, c_hgvs, tx, version, bucket, converted_g, fix_codes` (CSV by default,
  parquet if pyarrow is installed) to a `--out` path under gitignored `output/`. `bucket`
  is the baseline raw resolution (no cleaning/version bump) scored against the ClinVar
  CLNHGVS, which is what R5b and R6 build on; `fix_codes` is opt-in (`--with-fixes`, off
  by default since it doubles the work and fires nothing on clean ClinVar). Pluggable
  provider shared with `benchmark_resolution.py` (REST default, `--json` for the offline
  full run). R5b/R6 consume this table instead of re-resolving. Verified on
  `tests/test_data/clinvar_hgvs/` (10/10 correct via REST).

- **R5b. Empirical validation of safe version substitution.**
  `paper/scripts/validate_safe_versions.py` reads the A1 table, and for every resolved
  variant whose transcript has another version W that `is_version_substitution_safe`
  calls safe, resolves the variant under W and compares to the stored under-V coordinate
  (`converted_g`, so V is not re-resolved). Emits `n_safe_substitutions` and
  `n_coordinate_changes` (the headline K, expect 0) to
  `output/facts/version_safety_validation.csv`. Verified end-to-end on the 10-variant
  table (9 variants with alt versions, 11 safe substitutions, K=0); the change-detector
  was confirmed to fire (K=1) on a deliberately-wrong stored coordinate, so the zero is
  real, not trivial.

- **R6. Categorize the ClinVar genomic-match failures.**
  `paper/scripts/categorize_genomic_mismatches.py` reads the A1 table's `incorrect` rows
  and buckets each with deterministic rules (no LLM): `normalization_equivalent` /
  `ref_base` / `build` / `gap_related` (positions differ after normalisation and the
  cited transcript has an alignment gap) / `genuine_disagreement` / `unparseable`. Emits
  the breakdown plus the corrected genuine-error rate to
  `output/facts/genomic_mismatch.csv`. Replaces the stale PyHGVS/GRCh37
  `investigate_fails.py` stub (still wired into the Snakefile `clinvar` rule — repoint
  that when the facts are promoted into the paper). All six categories verified on
  crafted public-coordinate cases plus the gapped test transcript `NM_001637.3`.

- **S4. Per-source / failure-reason breakdown (Supplementary Table S4).**
  `paper/scripts/summarize_clinvar_pass.py` aggregates the A1 table into the S4 breakdown:
  per source (RefSeq vs Ensembl) the resolved/matched/no_data/error counts, plus an
  optional `--split-no-data` enrichment (needs a provider) that splits unresolved into
  unknown-accession vs unknown-version. Pure aggregation, instant, no provider for the
  base table. Writes `output/facts/clinvar_source_breakdown.csv` (+ a `_failure_reasons`
  file when split). This is the cdot-only full-ClinVar view (the R1 cdot-vs-UTA numbers
  stay sampled because UTA gates them). **Caveat baked into the script and to be stated in
  the paper:** ClinVar is dominated by a few large (largely US) labs citing current RefSeq
  versions, so S4 is a scale sanity check, not an unbiased transcript sample; the Shariant
  historical corpus (R1 Tier 2) is the unbiased complement. Fills the S4 stub in
  `supplementary.md` once the full run produces real numbers (paper-template wiring is the
  remaining step). Verified on the 10-variant table and a synthetic mixed-source table
  (no_data split correctly separates unknown-accession from unknown-version).

Notes from the "what else needs a full ClinVar run" review (2026-06-19): **R4** (version-
fallback false-rescue) is already done. its ablation is `benchmark_resolution.py
--recovery` (`absent_version` degraded mode) and the R4 Results section is written;
running it at full ClinVar scale would only restate R5b. The injection-cleaning benchmark
(R2 Tier-1) could be scaled from its seeded sample to full ClinVar to harden the "0
regressions" claim, but it is supplement-only and low marginal value.

