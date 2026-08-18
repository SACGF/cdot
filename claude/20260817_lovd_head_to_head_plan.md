# LOVD syntax checker head-to-head plan

Follow-up to item 3 of the (now deleted) 20260808 paper feedback plan. The related-work
paragraph in `paper/discussion.md` now cites the LOVD HGVS syntax checker
(github.com/LOVDnl/HGVS-syntax-checker) as a repair tool. A Bioinformatics reviewer is
likely to ask how `clean_hgvs()` compares, so run the comparison pre-emptively: score both
tools on the same corrupted-HGVS corpus and report recovery head-to-head.

## Feasibility facts (assessed 2026-08-17)

- The checker is a self-contained PHP library, PHP >= 7.4 with json and zlib extensions,
  runnable fully offline from the CLI: `php -f HGVS.php "NM_004006.3:c.157C>T"` emits JSON
  with `valid` and `corrected_values` (corrections ranked by likelihood, each with a
  confidence score).
- Installable via composer as `lovd/hgvs-syntax-checker`, or plain clone of the repo.
- Public REST fallback: api.lovd.nl checkHGVS v2 (single or batch), rate-limited to about
  5 variants/s or 1 batch/s, no registration. Local CLI is minutes for the whole corpus,
  REST would be ~12 minutes; use the local CLI.
- Estimated effort: about half a day including the PHP install and a subprocess wrapper.

## Corpus

`paper/scripts/inject_and_clean.py` generates ~3,400 corrupted strings by injecting known
error classes into committed ClinVar c.HGVS (all public data, safe to use anywhere). Each
corrupted string has a known canonical target, so scoring is exact-match recovery. Reuse
its injection machinery and categories; do not invent a new corpus.

## Method

1. Install PHP + the checker locally (pin the checker version and record it in Methods for
   reproducibility; a composer.lock or pinned git SHA in the script docstring is enough).
2. New script `paper/scripts/lovd_head_to_head.py`: subprocess wrapper around the PHP CLI,
   feeding each corrupted string, parsing the JSON, and scoring. Batch via the CLI's own
   batching if it has one, otherwise one call per string is fine at local speed.
3. Scoring must be decided and documented up front, then applied identically to both tools:
   - Primary metric: top-ranked LOVD correction equals the canonical target, versus
     `clean_hgvs()` output equals the canonical target.
   - Secondary (report, do not headline): target appears anywhere in LOVD's ranked
     `corrected_values` list.
4. Fairness stratification, the key design point: LOVD's checker is deliberately
   sequence-agnostic, so injection categories that require transcript data (anything a
   data-provider-aware fix handles, eg accession-prefix restoration, version work) are
   structurally out of its scope. Split results into two strata: string-repairable
   (both tools eligible) and data-requiring (cdot only). Report the head-to-head on the
   first stratum and the second as cdot-only capability. Do not present a blended number
   that makes LOVD look artificially weak.
5. Also record LOVD's behaviour on the uncorrupted originals (false-correction rate), since
   the paper claims `clean_hgvs()` guarantees no regressions; the comparison is only
   meaningful with the same check applied to LOVD.

## Deliverables

- `paper/scripts/lovd_head_to_head.py` (plus whatever tiny install notes it needs in its
  docstring; the PHP install itself stays out of the repo).
- Facts CSV under `paper/empirical_results/` (eg `lovd_comparison.csv`) following the
  existing pattern, wired into `paper/Snakefile` as a frozen-constants rule with a
  provenance docstring, matching how the other measured facts are handled.
- One short paragraph in `paper/discussion.md` extending the existing related-work
  paragraph with the measured comparison (or a sentence there plus a supplementary table
  if the strata need a table). Methods gets the scoring and stratification rules.
- No CHANGELOG entry (paper/analysis only). No em-dashes in prose. Nothing from
  `../cdot_private` (the injection corpus is public ClinVar, so this is naturally safe).

## Caveats

- If the checker cannot be installed offline on this machine (no PHP, composer blocked),
  fall back to the REST API with the rate limit respected, and note the service date in
  Methods instead of a version pin.
- LOVD may legitimately return multiple plausible corrections where the injection is
  ambiguous; that is what the secondary metric captures. Do not count ambiguity as failure
  in prose without saying the top-1 rule caused it.
