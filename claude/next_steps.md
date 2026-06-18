# cdot paper — next steps (deferred from feedback triage, 2026-06-18)

These are the items from `claude/paper_feedback_2026-06-18.md` that are **not** pure paper
edits — they need new analysis runs or code changes (some touching the private corpus or
the frozen Tier-2 constants). The paper-edit items from that feedback are already applied.

## C1. R4 version-fallback — concrete numbers (cheap, recommended)
R4 currently states the fallback works but reports **no numbers** — a real weakness the
feedback flagged. Add, over the committed ClinVar (g,c) pair set:
- how often a requested transcript version is **absent** and the fallback fires;
- of those, how often the substituted version still resolves, and how often it lands on a
  **different genomic coordinate** (i.e. the bump silently changed the answer — the risk case).
Run via `paper/scripts/benchmark_resolution.py` (it already has the degraded "absent_version"
mode) over `tests/test_data/clinvar_hgvs/`. Tier 1, reproducible. ~minutes.

## C2. Transcript-version drift study (ambitious, novel — its own analysis)
The deeper question the feedback raised: *how often is a version bump actually safe?*
Walk adjacent transcript-version pairs (same accession, consecutive versions) across cdot's
historical data and measure how many c. positions map to a **different** g. coordinate
(exon-boundary / CDS-offset drift). Output: distribution of "fraction of bases whose g.
coordinate is unchanged across a version bump", by gene/transcript.
- Genuinely novel ("not sure anyone has measured this") — could be its own figure or a
  short standalone section, possibly its own note.
- Cost: hours+ (needs a coordinate-walk harness over many version pairs; design first).
- **Decision needed before doing.** Keep OUT of the paper until run and reviewed.

## B3. Fix GENOMIC_REF_IN_PARENS in clean_hgvs (code + refresh frozen R3 numbers)
Feedback: "we could probably fix this." A genomic `NC_` accession sitting in the
gene-symbol parenthetical slot (`NM_…(NC_000013.11):c.…`) is detectable and strippable.
This is **code**, not a paper edit, and it changes the frozen Tier-2 numbers:
1. add a cleaning op (detect `NC_`/`NG_`/`NW_` in the gene-paren slot → drop or relocate);
2. unit test on synthesised public examples (no corpus strings);
3. re-run `cdot_private/analyze_cleaning.py` (corpus is at `../cdot_private`) to refresh the
   rescue rate (currently 96.4%) and the residual taxonomy — GENOMIC_REF_IN_PARENS (57
   queries, 4.9%) should shrink toward zero;
4. transcribe the new constants into `paper/results.md` R2/R3 + the Snakefile rules.
Greenlit by Dave (cdot_private reachable). Do as the immediate follow-up code task.
Note: changes R2 headline (rescue %) and R3 table — keep Tier-2 provenance flags.

## B1 (done in paper) / B2 (no action)
- B1: framed in-paper — Methods now contrasts unambiguous cleaning vs the opt-in,
  never-automatic version/canonical heuristics (reported as `HGVSFix`); Discussion has a
  matching safety/limitation note. No new `is_unambiguous()` API (redundant with `HGVSFix`
  severity) unless we later decide we want it.
- B2: bibliography is author-date (CSL) → sorted alphabetically by surname; "biocommons #2"
  isn't hand-assignable and Dave accepted the journal's ordering. No change.
