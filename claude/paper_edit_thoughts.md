# cdot paper: editorial thoughts and proposed edits

Read target: `paper/build/.last_render.md` (rendered 2026-07-01, the most recent build), which
is the fully assembled paper with facts injected. Notes below are organised from the highest
level (what the paper is and how it's ordered) down to individual sentences, ending with a
quick-wins checklist.

Overall: this is a strong, honest, well-evidenced applications note. The science and the
provenance discipline (Tier 1 / Tier 2) are its best features. The main weaknesses are (1) a
couple of blocking data/placeholder bugs, (2) it is too long and too evenly weighted for a
*Bioinformatics* Applications Note, so the headline results are buried, and (3) sentence-level
density: many sentences carry three or four clauses and have to be re-read.

---

## 0. Blocking issues (fix before anything else)

- **R1 per-source split is broken.** The text reads "on RefSeq the two are close (cdot
  **0.0%**, UTA **0.0%**), whereas on Ensembl cdot resolves **0.0%** and UTA **0.0%**". The
  fact file `output/facts/clinvar.csv` has `n_refseq=0, n_ensembl=0, cdot_refseq_pct=0.0, ...`,
  so the injection produced zeros. This is the single most important sentence in R1 (it
  "locates the gap") and it currently says the opposite of the intended claim. Regenerate the
  fact or fix the pipeline before submission; right now the paragraph is self-contradicting
  (cdot resolved 99.0% overall but "0.0%" per source).
- **Unfilled placeholders**: `[Zenodo DOI]`, `[email]`, `**Contact:** [email]` in the
  abstract; Table S1 / S2 bodies are literal `...`; Table S3 lists schema `v0.2.33` while the
  code example uses `cdot.0.2.33` (consistent, good) but the supplement tables have no data.
  These must be populated.
- **Consistency check on the two headline throughput numbers.** R6 body says local JSON
  resolves the full ClinVar set at "665 HGVS/s" and REST at "731 HGVS/s", but Table 2 lists
  local JSON as "540-665" and REST-prefetch as "~731". The 665 is the top of the range used as
  a point estimate at scale; fine, but state once that the 665 is the large-run figure so a
  reader doesn't think the table and text disagree.

---

## 1. Structure / framing

### Title
"cdot: resolving as many real-world HGVS strings as possible" is informal and sells only the
cleaning/robustness angle. The paper's biggest quantitative wins are coverage (11.5x more
alignments, Ensembl, T2T) and speed (up to ~6000x vs remote UTA, ~30x vs local UTA). A more
conventional and more complete title would foreground the data resource, e.g.
"cdot: fast, versioned RefSeq and Ensembl transcript data for HGVS resolution". Keep the
"real-world" framing as the sub-theme it already is in the intro.

### Results ordering
Seven result sections (R1-R7) is a lot, and the order buries the strongest, most reproducible
selling points. Suggested regrouping (7 -> 5):

1. **Coverage** = current R1 coverage paragraph **+ R7 (T2T)**. T2T is a coverage story; it
   should not be a standalone final section that reads as an afterthought. Fold the 186,644
   T2T alignments into the coverage count discussion.
2. **Resolution accuracy** = R1 ClinVar sample + full-scale + Shariant. (Keep as is, but see
   §2: it is three results under one heading and should be signposted.)
3. **Throughput** (current R6). This is a headline result and sits sixth. Move it up to third,
   right after accuracy, because "as accurate as UTA and 30x-6000x faster" is the cleanest
   one-line pitch and it is fully Tier 1.
4. **String cleaning** = R2 + R3. R3 (residual taxonomy) is a natural "how far can this go"
   coda to R2 and does not need its own top-level heading; make it the closing subsection of
   the cleaning result.
5. **Version fallback and safety** = R4 + R5. These are one feature (the fallback) and its
   safety proof; R4 is only 11 lines and reads as a preamble to R5. Merge.

This turns 7 sections into 5 and puts the two most reproducible, most impressive results
(coverage, throughput) first.

### Tier 1 / Tier 2 provenance flags
The discipline is excellent and should stay, but the *presentation* is heavy: a block-quote
disclaimer opens Results and then almost every paragraph and table caption is prefixed with a
bold `**[Tier 1]**` / `**[Tier 2: production validation, not reproducible]**`. Because Tier 1
is the default, tagging it everywhere adds noise. Suggestion: keep the opening definition,
then tag **only** Tier 2 items inline (a short "(production corpus; not referee-reproducible)"
is enough) and drop the repeated `[Tier 1]` bold prefixes. This alone de-clutters the whole
Results section.

### Add one figure
The paper is all prose and tables. A *Bioinformatics* note benefits from a single schematic:
either (a) the data flow (RefSeq/Ensembl/UTA releases -> merge -> JSON.gz -> {local provider,
REST, TARK} -> biocommons/PyHGVS), or (b) a two-panel bar chart of coverage (cdot vs UTA, by
source/build) and throughput (Table 2 as a log-scale bar). Option (b) makes the two headline
numbers visual and is trivial to generate from the existing fact CSVs.

---

## 2. Section-by-section

### Abstract
- Paragraph 3 is overloaded: cleaning + version fallback + both libraries + MANE + T2T in one
  block. Split or trim. The MANE/gene-symbol lookup and "both libraries" points are
  secondary; consider demoting them to a single trailing sentence so the paragraph's spine is
  cleaning + safe fallback.
- "resolving as many ... as possible" echoes the title; once is enough.

### Introduction
- Strong. The four-paragraph arc (what HGVS resolution needs -> the two libraries + UTA's
  limits -> Shariant motivation -> what cdot does) is the right shape.
- Paragraph 1 is a single long build-up. The ANNOVAR "44% of LoF" stat is a slightly indirect
  argument for covering both annotation sources; it works but consider tightening the framing
  so the reader sees immediately why it is cited (transcript-source choice changes the
  answer, therefore cover both).
- The UTA "~1 transcript/second" and coverage "~141,000" numbers recur in abstract, intro,
  R1, R6. Some overlap is expected, but the intro can state them once and let Results own the
  measured versions.

### Methods
- Longest section (~230 lines) and carries content that is really Results or Reference
  material. Specifically, the **String cleaning** subsection's four bullet groups (Stripping /
  Structural punctuation / Casing / Reconstruction) read as documentation. Keep a short prose
  description of the pipeline (pure string op, ordered, records `HGVSFix`, never raises by
  default) in Methods and move the exhaustive operation-by-operation catalogue to a
  supplementary table. That table is genuinely useful, just not in the main Methods flow.
- Methods contains forward-looking result numbers ("540-665 transcripts/second (Results,
  Table 2)", "~11 seconds"). Fine to mention scale once, but avoid stating the measured
  throughput in Methods and again in R6; cite it once.
- The `get_best_transcript_version()` paragraph and the R5 safety gate describe the same
  machinery from two directions; make sure Methods sets up the mechanism and Results reports
  the measurement, without each re-deriving the "never applied automatically / opt-in /
  preserves exact-variant semantics" point (it appears in Methods twice, R4, R5, and
  Discussion, ~5 times total).

### Results R1
- Three distinct results under one heading (Tier 1 ClinVar sample vs UTA; Tier 1 full-scale
  ClinVar-only; Tier 2 Shariant). Add sub-labels or split. A reader currently has to track
  which population and which comparator each percentage refers to.
- Fix the 0.0% bug (see §0).
- The full-scale ClinVar paragraph and Supplementary Table S4 restate the same
  VCF-coordinate-scoring rationale ("we score as VCF coordinate ... so equivalent
  representations are not miscounted") almost verbatim. State it once in the main text, refer
  to S4 for the breakdown.

### Results R2 / R3 (cleaning)
- R2 is clear and the Table 1 breakdown is good. The "+5.1% absolute (~1,700 strings)" is a
  clean headline; keep it prominent.
- R3's residual taxonomy is thorough and honest, but it is heavily Tier 2 and rests on a
  single-rater LLM classification (disclosed). For a short note this is a lot of real estate
  on the *limits* of a secondary feature. Recommend compressing R3 to: the ceiling number
  (3.4% residual), the three-way grouping (~50% unfixable/reference-less, ~33% fixable
  frontier, ~14% grammar-gap + non-HGVS), and push the full 7-class table to the supplement.
  The LLM-method caveat can move to the supplementary table legend.

### Results R4 / R5 (fallback + safety)
- Merge (see §1). R4 is a short preamble; R5 is the substance.
- R5 is the most technically dense passage in the paper (intrinsic CDS structure,
  build-independence, gap and UTR refinements, the 2-in-4,024,794 residual, the blocklist).
  It is a genuine and clever contribution but currently ~4 dense paragraphs. Tighten to: the
  claim (a version bump is safe iff intrinsic CDS structure + gaps + UTR length are
  preserved), the empirical validation (2 residual movers out of ~4M accepted substitutions,
  both one re-placed transcript, handled by a blocklist), and one sentence on why
  build-independence matters (lets the check work when the requested version was never aligned
  to the target build). The mechanism-vs-measurement split with Methods will help here.

### Results R6 (throughput)
- Strong and fully Tier 1; move earlier (see §1).
- The "sequence layer held constant / shared local SeqRepo so only the transcript layer
  varies" caveat appears in Methods (Benchmarking), the R6 intro paragraph, the Table 2
  caption, and the R6 body. Four times. State once (table caption is the natural home) and
  delete the rest.

### Discussion
- Tight and well-structured ("two gaps" framing is effective).
- It restates UTA's limitations (no DB, no Ensembl, limited history) that the intro already
  made. A discussion can recap, but this one recaps at nearly the original length; trim the
  first paragraph to a sentence and spend the space on the forward-looking automation point,
  which is currently one line at the very end.

### Supplement
- Populate S1/S2 (currently `...`). S3 schema table is good. Consider a supplementary table
  for the full `clean_hgvs()` operation catalogue moved out of Methods.

---

## 3. Repetition to cut (recurring phrases)

These claims are each stated 3+ times; pick one home for each and delete the rest.

| Claim | Appears in | Keep in |
|---|---|---|
| "sequence layer held constant / same local SeqRepo, only transcript layer varies" | Methods, R6 intro, Table 2 caption, R6 body | Table 2 caption |
| "never applied automatically / opt-in / preserves exact-variant semantics" | Methods x2, R4, R5, Discussion | Methods (define) + Discussion (recap) |
| "first to bring T2T-CHM13v2.0 to the Python HGVS libraries" | Abstract, Intro, R7 | Intro + R7 (merged into coverage) |
| "identical biocommons/hgvs code path / engine, only data layer swapped" | Methods, R1 (x2), R6 | R1 (define once) |
| "1,621,073 alignments" / "~141,000 UTA" | Abstract, Intro, R1 | Abstract (claim) + R1 (measured) |
| VCF-coordinate scoring rationale | R1 body, Table S4 | R1 body |
| UTA limitations (no DB, no Ensembl, limited history) | Intro, Discussion | Intro (full) + Discussion (one line) |

---

## 4. Sentence-level patterns

The dominant issue is **multi-clause sentences chained with "so", "because", ";", and "which"**
that pack three ideas where a reader can hold one. A few representative offenders and rewrites:

- Abstract, "cdot also helps ... T2T-CHM13v2.0 assembly to the Python HGVS libraries." — one
  paragraph, five ideas. Split after "before resolution." and again after "does not move the
  variant."

- Methods JSON, lines ~151-156: "The build is always a dict key ... with no reshaping when
  builds are combined." This is one ~90-word sentence with an em-of clauses. Break into:
  (1) the build is always a dict key, even for a single build; (2) so a single-build file is a
  drop-in GTF/GFF replacement and a multi-build file needs no reshaping; (3) the REST API can
  therefore return all builds for a transcript in one response.

- R5, "A transcript version's intrinsic CDS structure ... 99.3% of drifting Ensembl bumps
  carry an intrinsic-structure change ..." — a single sentence spanning the definition, the
  build-independence stat, and the two-directional validation. Split into three: definition;
  build-independence; predictive validity.

- R6, "batching those lookups into one prefetch() request amortises the per-request network
  and processing overhead across the whole set and warms the transcript cache, so later
  lookups are in-memory hits and REST throughput rises to match local." — split at "and warms
  the transcript cache."

General guidance for a pass:
- Target one main clause + at most one subordinate clause per sentence in Results.
- Prefer a period over a semicolon where the two halves are independent claims (most are).
- Cut "so that", "which means", "in practice" connectors where the next sentence can just
  state the consequence.

### Style compliance (CLAUDE.md)
- No em-dashes found in the rendered text (good; the style rule is being followed).
- Watch a couple of near-tricolons and mild hedges ("in practice", "roughly"). Not egregious.
- "underscoring the value of covering both annotation sources" (intro) is a light editorial
  flourish; a plain "so covering both sources matters" is more in-voice.

---

## 5. Small factual / clarity nits

- R1: "a 11.5x increase" -> "an 11.5x increase".
- R1 counts: 1,621,073 total, of which 574,584 Ensembl absent from UTA and 186,644 T2T.
  Make explicit whether the 11.5x is total-vs-total (1.62M / 141k = 11.5x, yes) so the reader
  can reconstruct it; currently the drivers (history + Ensembl + T2T) are listed but the
  arithmetic isn't tied back to the 11.5x.
- R2: "+5.1% absolute gain (~1,700 additional strings rescued)" but Table 1 total is 1,678 and
  the text elsewhere says 1,678. Use 1,678 consistently rather than "~1,700" in one place and
  1,678 in the table.
- R6: "roughly 30x faster than a local UTA" (665 vs 24 = ~28x) and "~6000x" vs remote (665 vs
  0.1) is implied but never stated as a single number; the abstract could use the clean
  "up to ~6000x faster than the public UTA server, ~30x faster than a local UTA" line.
- Table 2 uses "~0.1", "~24", "~39", "~731", "540-665" — mixed tilde/range style. Pick one
  (e.g. all point estimates with a footnote that they are medians of N runs).
- R3: "the eighth was non-HGVS input (81 queries, 7.2%)" then "a further ~7% (excluded above)
  was non-HGVS" — the 7.2% is stated twice; and note the Grammar-gap class is also 81 (7.2%),
  so "7% grammar gap + 7% non-HGVS" collides numerically with two different 81s. Disambiguate.

---

## 6. Quick-wins checklist (in priority order)

1. Fix the R1 `0.0%` per-source numbers (broken fact injection). **Blocking.**
2. Fill placeholders: Zenodo DOI, contact email, Table S1/S2 bodies. **Blocking.**
3. Merge R4+R5 (fallback+safety) and fold R7 (T2T) into R1 (coverage): 7 sections -> 5.
4. Move throughput (R6) up to third; it is the strongest reproducible pitch.
5. Drop the repeated `[Tier 1]` bold prefixes; tag only Tier 2 inline. Keep the opening def.
6. Cut the four repeated "sequence-layer-held-constant" statements down to one.
7. Move the `clean_hgvs()` operation catalogue and the R3 7-class table to the supplement.
8. Add one figure (data-flow schematic or coverage+throughput bar chart).
9. Sentence-splitting pass on the JSON-format, R5, and abstract long sentences.
10. Retitle to foreground the transcript-data resource, keeping "real-world" as sub-theme.

Nothing here changes the paper's claims or evidence; it is about making the strongest results
land first and reducing the reread rate. The underlying work is solid and unusually honest
about its own limits, which is worth preserving through the trimming.
