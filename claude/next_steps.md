# cdot paper — next steps (deferred from feedback triage, 2026-06-18)

These are the items from `claude/paper_feedback_2026-06-18.md` that are **not** pure paper
edits — they need new analysis runs or code changes (some touching the private corpus or
the frozen Tier-2 constants). The paper-edit items from that feedback are already applied
and pushed (commit `4fafe87`).

## Conventions for future analysis work
- **LLM classification/assignment tasks: emit SHORT class names** (e.g. `Truncated`,
  `No reference`, `Bad accession`, `Grammar gap`, `Genomic ref in parens`), not long
  ALL_CAPS taxonomy keys. The R3 residual table had to be hand-shortened after the fact;
  do it at the source next time so tables need no post-processing.
- **Never put private-corpus strings in this repo** (tests/comments/docstrings/paper).
  Synthesise from public examples (BRCA2 `NM_000059.4`, RUNX1 `NM_001754.5`).

## B3. Fix GENOMIC_REF_IN_PARENS in clean_hgvs (code + small paper number update)
Feedback: "we could probably fix this." A genomic `NC_`/`NG_`/`NW_` accession sitting in
the gene-symbol parenthetical slot (`NM_…(NC_000013.11):c.…`) is detectable and
strippable so the string parses.

**Scope (kept deliberately small per Dave):**
1. Add a cleaning op to `cdot/hgvs/clean.py` (detect a genomic accession in the gene-paren
   slot → drop/relocate it), emitting an `HGVSFix`. Unit-test on **synthesised public**
   examples only.
2. `CHANGELOG.md` `[unreleased]` → `### Fixed`/`### Added`: new cleaning capability,
   client-visible, reference `#112`.
3. Paper update is minimal — **just remove the category, no reclassification.** The
   GENOMIC_REF_IN_PARENS residual class (57 rows / 52 unique) now parses, so:
   - drop the `Genomic ref in parens` row from the R3 table;
   - the other classes' raw counts are unchanged (do **not** re-run the LLM classifier);
   - bump the R2 rescue total (96.4% → slightly higher; +~57 rows rescued) and the
     residual total (1,175 → ~1,118 rows; 913 → ~861 unique);
   - recompute the remaining classes' % against the new residual total (arithmetic, not
     reclassification), or keep counts and note the removed class — pick whichever reads
     cleaner; keep Tier-2 provenance flags.
4. Optional exact-number check (cheap, deterministic, NOT an LLM run): re-run
   `clean_hgvs()` over the corpus via `../cdot_private/analyze_cleaning.py` to get the
   precise new rescue/residual counts. `cdot_private` is the sibling dir `../cdot_private`.

## C1. R4 version-fallback — concrete numbers (LOW expected value; deprioritised)
R4 reports no numbers. **Caveat (Dave):** ClinVar cites very recent transcript versions and
the as-is success rate is already very high, so the version-bump path rarely fires on the
ClinVar set — this measurement probably shows little. If done at all, run the existing
`benchmark_resolution.py` "absent_version" degraded mode over `tests/test_data/clinvar_hgvs/`
and report how often a bump (a) fires, (b) still resolves, (c) lands on a *different* g.
coordinate. Cheap, Tier-1 — but expect a near-null result. The interesting version question
is C2, not this.

## C2. Transcript-version drift study (the real version question — explore next)
How often is a version bump actually *safe*? Walk adjacent transcript-version pairs (same
accession, consecutive versions) across cdot's historical data and measure how many c.
positions map to a **different** g. coordinate (exon-boundary / CDS-offset drift). Output:
distribution of "fraction of c. bases whose g. coordinate is unchanged across a version
bump", per transcript/gene; how often a bump is coordinate-preserving vs not.
- Genuinely novel ("not sure anyone has measured this"); could become its own figure/section.
- **For now: explore only.** Write a standalone analysis script, run it on available local
  data, and **report findings back** — do **NOT** wire it into Snakemake or the paper yet.
- Local data (see memory `benchmark-data-locations`): cdot JSON at `/data/cdot_data/`
  (refseq/ensembl/all-builds), genome FASTA, local UTA, SeqRepo. Use test-scale subsets;
  per CLAUDE.md never run the full GTF/GFF generation pipeline.

## Done / no-action
- **B1** framed in-paper (Methods contrasts unambiguous cleaning vs opt-in, never-automatic
  version/canonical heuristics reported as `HGVSFix`; Discussion has a safety/limitation
  note). No `is_unambiguous()` API unless we later want it.
- **B2** bibliography is author-date (CSL) → alphabetical by surname; "biocommons #2" isn't
  hand-assignable and Dave accepted the journal ordering. No change.
- **Facts plumbing note:** R2 Table 1 and the R3 taxonomy table are hardcoded literal
  markdown (Tier-2 frozen, hand-transcribed — repo convention), NOT vibepaper facts. Only
  the new historical benchmark and single cleaning values are templated. If we want these
  tables facts-driven (so B3-type refreshes are a one-line CSV edit), promote them to a
  `residual.csv`/extended `cleaning.csv` + Snakefile rule — optional, bundle with B3 if desired.
