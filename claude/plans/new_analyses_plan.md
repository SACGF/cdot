# New analyses to do for the paper (deferred)

From Dave's feedback (`claude/paper_feedbac_2026-08-19.md`, 2026-08-19). These need
re-running scripts / new data work, so they are parked here. The prose/framing edits
that do **not** need new data have already been applied to the paper.

Constraint reminder: never run against full datasets (CLAUDE.md). Use samples / local
UTA / committed test data. No `cdot_private` example strings in-repo.

---

## 1. Fair time-bucketed sampling of the ClinVar submitted (SCV) corpus  [R2]

**Why:** the current submitted-corpus sample is a single fixed-seed flat draw over all
unique (AlleleID, string) pairs, which has recency bias (ClinVar grows over time, so a
flat draw over-represents recent submissions).

**Do:**
- In `build_clinvar_submitted_pairs.py`, add time-bucketed sampling: bucket SCV records
  by submission date (or release/version era) and draw evenly across buckets, so
  historical submissions are represented fairly.
- Keep a **second** run that is a whole-file random sample, reported alongside and
  explicitly labelled as recency-biased (reflects the live distribution).
- Report both. The time-bucketed one is the fair "what labs submitted over the years"
  number; the random one is "what the current file looks like".

**Touches:** R2 numbers in `results.md`, `clinvar_submitted.csv`, Methods description of
the sampling.

## 2. Run local UTA over (all of) ClinVar  [R2 / R3]

**Why:** we currently only sample the cdot-vs-UTA comparison "because UTA is slow", but
that caveat is for the *remote* server. The **local** `uta_20241220` may be able to do
the whole ClinVar set in a few hours.

**Do:** run the full ClinVar resolution through local UTA (background job), so the
head-to-head is full-scale, not sampled. Removes the "gated to a sample by UTA
throughput" caveat for the local case.

**Touches:** R2 (drop/reduce the sampling caveat), possibly R3.

## 3. Why is version substitution needed at all? (submitter attribution)  [new analysis]

**Why:** we ingest every RefSeq/Ensembl release off the FTP sites, so in principle we
should have every published version. Yet labs cite versions we don't hold. Dave's
hypothesis: those transcripts come from a small number of labs that align transcripts
to the genome themselves (so the version never appeared in an official annotation
release).

**Do (exploratory first, time-boxed):**
- Take the cited transcript versions that are absent from cdot's data.
- Group by submitter (SCV submitter / lab) in the ClinVar VCV XML.
- Test whether the absent-version citations concentrate in a small number of submitters,
  or whether specific labs dominate.
- Characterise: are these self-aligned transcripts, pre-release versions, or genuine
  gaps in our ingest?

**Scope risk:** this is a genuinely new result and could balloon. Decide whether it
earns a place in the paper (probably a short paragraph in R2 or Discussion) before
investing heavily.

## 4. Warm-cache benchmark rerun  [R3]

**Why:** the full-scale R3 paragraph currently says the local-JSON and REST runs' "cache
conditions differed, so their small difference is not evidence of REST outrunning local
JSON." That is a hand-wave. Replace it with a hot-cache protocol.

**Do:** for the full-scale runs, do one throwaway pass to warm the sequence-layer cache,
discard it, then time the next pass. Report the hot-cache time and note it is only
marginally slower than cold. Removes the "cache conditions differed" sentence entirely.

**Touches:** R3 full-scale paragraph, `benchmark.csv` (if regenerated), Methods.

## 5. Drop non-HGVS input from the cleaning residual corpus  [R4]

**Why:** the R4 residual currently includes "non-HGVS input (pasted URLs, report
templates, prose)". Dave's point: that is a data-collection artifact, not a measure of
tool correctness, and should be removed from the denominator before counting the
residual, since it does not reflect what the cleaner could ever be expected to fix.

**Do:** filter non-HGVS input out of the production cleaning corpus before computing the
residual rate, so the residual reflects genuine HGVS-repair ceiling only. Recompute the
R4 residual number and the Table S6 taxonomy denominators.

**Touches:** R4 residual paragraph + Table 2 context, `cleaning.csv`, Supplementary
Table S6.

---

## Also consider (from feedback, smaller)

- ~~**VEP + RefSeq alignment gaps citation:**~~ DONE. R1 now cites Ensembl/ensembl-vep
  issue #1053 (`@VepHgvsGaps` in references.bib) for the claim that VEP does not apply
  transcript-to-genome alignment gaps when converting c. to genomic coordinates.
