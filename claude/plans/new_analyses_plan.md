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

**Code DONE 2026-08-19 (run pending):** `build_clinvar_submitted_pairs.py` now captures
each pair's earliest SCV `SubmissionDate` from the VCV XML and writes it as a trailing
`submit_date` corpus column (downstream `iter_vcf_pairs` reads by name, so it is
ignored there). Two draws: `--sample-out` (the existing whole-file random draw,
recency-biased) and `--time-bucketed-out` (new, even across submission-year eras,
`--bucket-years` width, short eras drawn fully with the deficit redistributed). Both can
be emitted from one build/`--from-pairs`. Validated on synthetic skewed corpora + a tiny
synthetic VCV XML (earliest-date tracking, SCV collapse, even allocation, undated-drop,
and the missing-column error path). STILL TO RUN (needs full VCV XML + VCF, a big job):
rebuild the corpus with dates, emit both samples, resolve each via
`resolve_clinvar_pass.py`, then update the R2 numbers + Methods prose to report both.
Note: the committed `clinvar_submitted_500.tsv` predates the date column and should be
regenerated when the corpus is rebuilt.

**RUN DONE 2026-08-19 (results pending prose):** rebuilt from `ClinVarVCVRelease_2026-06`.
New corpus = 3,198,528 unique (AlleleID, string) pairs, 100% RefSeq / 0 Ensembl, all with
a submission date, 19 eras 2008-2026. (Corpus is larger than the committed 2,933,667: the
XML path scans to the first *transcript* HGVS per assertion, capturing SCVs whose first
HGVS attribute is genomic/protein, which the CSV-extraction path dropped. More complete
and self-contained.) Raw distribution is heavily recency-skewed: 2024 = 1.14M, 2025 =
1.21M pairs vs a few thousand per year in 2010-2013.

cdot resolution (0.2.34 refseq GRCh38, FastaSeqFetcher, --with-fixes) on 3,000-pair
seed-42 samples, VCF-coordinate scored:

| | correct (after fix) | residual | 2008-2015 | 2016-2020 | 2021-2026 |
|---|---|---|---|---|---|
| RANDOM (recency-biased) | 2969 (99.0%) | 31 (1.0%) | 96.9% (n=32) | 99.2% | 99.0% |
| TIME-BUCKETED (fair) | 2939 (98.0%) | 61 (2.0%) | 95.6% (n=1072) | 99.2% | 99.3% |

The fair draw is 2x harder (2.0% vs 1.0% residual) because it up-weights old submissions:
the random draw holds only 32 pre-2016 pairs (1%), the fair draw 1072 (36%), and pre-2016
strings resolve at ~95.6% vs ~99% for recent ones. Residual composition also differs:
random residual is mostly `no_data` (recent strings citing versions cdot lacks), fair
residual is mostly `incorrect`/`error` (old versions projecting to drifted/invalid coords).
Pass CSVs: scratchpad `random_pass_cdot.csv`, `bucketed_pass_cdot.csv`. Corpus + samples
under `/data/clinvar/`.

Dave's decisions 2026-08-19: ADOPT the 3.20M corpus; run UTA head-to-head on the samples
(done); HOLD item 1.

Full measured numbers (new 3.20M corpus, cdot 0.2.34 refseq GRCh38, local uta_20241220):

cdot-vs-UTA head-to-head (VCF-scored):
| draw | cdot correct | UTA correct | UTA no_data | cdot-only | uta-only |
|---|---|---|---|---|---|
| RANDOM (recency-biased) | 2969 (99.0%) | 2466 (82.2%) | 530 (17.7%) | 504 (16.8pts) | 1 |
| TIME-BUCKETED (fair)    | 2939 (98.0%) | 2396 (79.9%) | 586 (19.5%) | 544 (18.1pts) | 1 |

The fair draw WIDENS the cdot advantage (16.8 -> 18.1 pts): historical submissions cite
more superseded versions UTA holds no GRCh38 alignment for.

Version-age on the new corpus (`version_age_new.csv`): version_not_current 75.1% (was
81.8% on the old 2.93M CSV-path corpus - the drop is the corpus redefinition),
scv_weighted_not_current 69.8%, base_retired 0.7%, not_current_in_cdot 99.3%,
absent_cdot_grch38 0.6%.

STILL TODO (R2 prose/CSV rewrite): report BOTH samples (bucketed = fair headline, random
= recency-biased); restructure `clinvar_submitted.csv` to carry both draws; pick which
sample the abstract reports; recompute the residual taxonomy
(`clinvar_submitted_residual.csv`) for the headline sample; regenerate committed
`clinvar_submitted_500.tsv` from the new corpus. Pass CSVs in scratchpad:
{random,bucketed}_pass_{cdot,uta}.csv.

## PR #128 review comments (Dave, 2026-08-19) - addressed in parallel

- DONE: renamed R2 heading "What laboratories submitted" -> "ClinVar submissions as a
  historical record of transcripts used" (+ reworded lead-in).
- DONE: turned the hardcoded R4 cleaning-corpus constants into a vibepaper fact
  (`cleaning_corpus.csv`, wired into Snakefile FROZEN_FACTS); templated N, rescued,
  rates, residual across results.md Table 2 + supplementary S6. Render verified.
- DONE (draft, for Dave review): abstract motivation now leads with the Shariant origin
  story. Intro already had the Shariant paragraph (lines 33-38).

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

**Code DONE 2026-08-19 (Dave: "just do everything"), run pending.** `submitter_attribution.py`
streams the VCV XML, attributes each SCV's first RefSeq transcript citation to its
submitter (`ClinVarAccession@SubmitterName`, else `ClinVarSubmissionID@submitter`), and
flags absent-from-cdot versions using the same version-set test as
`compute_submitted_version_age.py`. Reports concentration (top-N submitter share of
absent citations) and the self-alignment signature (absent versions cited by a single
submitter). Validated on a synthetic XML. Full run deferred until item 1's timed pass
finishes (heavy parse would confound item 1 throughput); will run alongside item 4.

## Item 4 RESULT (full-scale, 2026-08-21)

Ran the full submitted corpus (3,198,528 pairs) through cdot (--fasta --with-fixes, ~20h
single-core) and local UTA (SeqRepo, ~7.5h). Full head-to-head:

| | cdot | UTA |
|---|---|---|
| matched VCF coordinate | 99.0% (3,166,066) | 81.9% (2,618,768) |
| resolves through this backend alone | 548,524 (17.2 pts) | 1,226 |
| no_data | 19,124 | 573,994 (18.0%) |

Per-era: 2008-2015 cdot 97.3% / UTA 76.2%; 2016-2020 99.3% / 88.8%; 2021-2026 99.0% /
81.0%. Residual 1.0% (32,462), 0 regressions. The full-scale numbers land near the
recency-weighted random draw; the per-era split is the fair view without sampling.
R2 rewritten to full-scale + per-era, dropping the two-sample framing (superseded).
Note: cdot --fasta is the bottleneck (~43 HGVS/s, historical versions force genome
reconstruction); UTA-full was actually the fast part.

## Status 2026-08-19 (PR #128 MERGED; items 1/4/5 on new branch paper-benchmarks-2026-08-19)

- Items 2, 3 + all PR #128 feedback: merged to main.
- Item 1 (warm-cache) + Item 4 (full UTA): RUNNING (sequential background runner PID
  357490, monitor b4g9mry9g). Item 1 first (~4h), then item 4 cdot-full + UTA-full (~15h).
- Item 5: code ready, run queued for after item 1's timed pass.
- When results land: update R3 (hot-cache throughput), R2 (full-scale head-to-head), add
  item 5 paragraph; open new PR from paper-benchmarks-2026-08-19.

## 4. Warm-cache benchmark rerun  [R3]

**Why:** the full-scale R3 paragraph currently says the local-JSON and REST runs' "cache
conditions differed, so their small difference is not evidence of REST outrunning local
JSON." That is a hand-wave. Replace it with a hot-cache protocol.

**Do:** for the full-scale runs, do one throwaway pass to warm the sequence-layer cache,
discard it, then time the next pass. Report the hot-cache time and note it is only
marginally slower than cold. Removes the "cache conditions differed" sentence entirely.

**Touches:** R3 full-scale paragraph, `benchmark.csv` (if regenerated), Methods.

## 5. Drop non-HGVS input from the cleaning residual corpus  [R4]  — DONE 2026-08-19

**Why:** the R4 residual currently includes "non-HGVS input (pasted URLs, report
templates, prose)". Dave's point: that is a data-collection artifact, not a measure of
tool correctness, and should be removed from the denominator before counting the
residual, since it does not reflect what the cleaner could ever be expected to fix.

**Do:** filter non-HGVS input out of the production cleaning corpus before computing the
residual rate, so the residual reflects genuine HGVS-repair ceiling only. Recompute the
R4 residual number and the Table S6 taxonomy denominators.

**Touches:** R4 residual paragraph + Table 2 context, `cleaning.csv`, Supplementary
Table S6.

**Done:** the R4/S6 residual and cleaning-rate numbers are literal frozen constants in
the prose (not a facts CSV; `cleaning.csv` is the unrelated injection benchmark and was
not touched). Kept the frozen LLM per-class taxonomy counts and moved the 81 non-HGVS
queries out of the corpus. Reran `cdot_private/analyze_cleaning.py` to confirm the
frozen baseline reproduces (as-is 29,956, after 31,676, residual 1,076 vs frozen 1,075,
±1 code drift). Recomputed against N = 32,752 − 81 = 32,671:

| quantity | was | now |
|---|---|---|
| corpus N | 32,752 | 32,671 |
| parseable as-submitted | 91.5% | 91.7% (29,956/32,671) |
| parseable after cleaning | 96.7% | 97.0% (31,677/32,671) |
| absolute gain | +5.3% | +5.3% (1,721/32,671 = 5.27%) |
| share of as-submitted failures rescued | ~62% of 8.5 pp | ~63% of 8.3 pp (1,721/2,715) |
| residual | 3.3% (1,075 q) | 3.0% (994 q) |
| S6 class %s | of 1,075 | of 994 (counts unchanged) |

The "826 unique strings" figure was dropped (the unique count *within* the 81 non-HGVS
needs the lost per-string LLM labels; residual is now reported in queries only).

Edited: `results.md` (main R4 para, Table 2 header, residual para), `supplementary.md`
(S6 header + table %s), `paper/README.md` (Tier-2 rate note). Paper re-rendered OK.

---

## Also consider (from feedback, smaller)

- ~~**VEP + RefSeq alignment gaps citation:**~~ DONE. R1 now cites Ensembl/ensembl-vep
  issue #1053 (`@VepHgvsGaps` in references.bib) for the claim that VEP does not apply
  transcript-to-genome alignment gaps when converting c. to genomic coordinates.
