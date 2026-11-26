# Runbook: rewriting the paper in the author's own words

The main text of the paper was AI-generated. The author (Dave) is rewriting it so that it is
his own prose and reads as human-written. This file is for an agent picking up the work.
Read it, then `outline.md` in this directory, before doing anything.

## Files

- `paper/rewrite/outline.md`: the reverse outline. One entry per target paragraph, in reading
  order, each with Job, Opens with, Lands on, Beats, Must keep (facts, citations, checked
  claims), Flags, then the author's dictation and the accepted text. The structural
  changes (splits, merges, moves) are listed at the top and have already been decided.
- `paper/*.md`: the paper sections. Accepted text is copied into these.
- `paper/*.claims.toml`: claims sidecars. Rewording breaks their quotes, and they are fixed
  after each section goes in.

## The loop, per one or two paragraphs

1. **Present the next paragraph(s)** from `outline.md`. Always start with **"Where you
   are"**: the section, the paragraph number out of how many, what the previous paragraphs
   established and what comes next. Then give the opening, the landing point and the beats
   as short notes, with real numbers filled in (look them up in
   `paper/empirical_results/*.csv`). The author asked for this context. Without it the
   beats are hard to write from.
2. **The author dictates.** The text arrives as speech-to-text, with slips ("HTVS",
   "Ensemble", "CDOT"), false starts and filler.
3. **Record the dictation verbatim** in `outline.md` under "Your version (as dictated)".
4. **Make a draft by cutting and joining only.** The author's rule: *you may chop, cut,
   reorder and join his sentences, but never add words.* Allowed:
   - deleting filler, false starts and repeats ("basically", "So", "And", "a lot of, well,")
   - joining or splitting sentences, and reordering his words
   - fixing speech-to-text slips and names (HGVS, Ensembl, GRCh37, cdot, PostgreSQL)
   - filling an "X" or a rough number he gave with the real value as a `{{ fact }}`, and
     adding `[@citations]`, since that is markup, not prose
   - inserting exact words or phrases he dictates as a correction ("put opt-in and flagged in
     brackets after...")
   Not allowed: rephrasing, tightening, "light tidies", or adding linking words. A Claude
   "light tidy" of his dictation scored 100% AI on Pangram. His own words, cut only, scored
   100% human every time (I1, I2, I3).
5. **Give feedback in this format.** Quote each sentence that needs work, say what's wrong
   with it, and let him re-dictate. For missing content, say where it goes and what it must
   say, with the numbers. Don't write replacement prose. Check:
   - structure: does it open on the topic and land on the point?
   - facts and citations dropped from the beats
   - factual errors. He welcomes these being caught. Examples so far: "p." can't be
     converted to genomic coordinates; UTA isn't the *only* provider, because cdot is one;
     the tool comparison runs on the public injected corpus, not the private search corpus.
   - overlap with paragraphs already written
   - filler that cutting can fix, and tense clashes that can't be fixed without a word from
     him
6. Once he's happy, mark the entry **"Accepted"** in `outline.md`.
7. **When a section is complete**, copy the accepted paragraphs into the section `.md`,
   rewrite its claims sidecar to match the new sentences, and make the checks pass:
   ```
   .venv/bin/snakemake -s paper/Snakefile --cores 1        # refresh output/facts
   .venv/bin/vibepaper lint --config paper/paper.toml --numbers --unused
   .venv/bin/vibepaper check --config paper/paper.toml
   .venv/bin/vibepaper wordcount --config paper/paper.toml
   ```
   A fact that is no longer cited either goes under `[lint] unused_ok` in `paper.toml` or is
   removed. Keep facts that are used elsewhere.

Pangram is not run on every paragraph any more. The method has been shown to work, so it
runs once on the whole main text at the end. He pastes the text in himself. Give it with
the numbers rendered and without `{{ }}` or `[@...]` markup.

## Each Results paragraph needs one central result about cdot

His rule (2026-11-26): a paragraph earns its place only if it can be summarised in one line
as a result, and that result is a claim about cdot that we care about. "We store gaps and
VEP doesn't apply them" fails: biocommons HGVS applies the gaps, not cdot, and VEP's
behaviour isn't our result (R1b was cut for this). A paragraph that only sets up another,
describes a table, explains a residual, or repeats another paragraph's result on a different
corpus gets merged, moved to the supplement, or deleted. Write the one-line result in the
outline entry before presenting the paragraph.

## When he asks for a change to data or facts

Do it properly, in the facts pipeline, not in the prose. For example, the UTA count is now
measured by `paper/scripts/compute_coverage.py` from the uta_20241220 dumps, as
`coverage.uta_count`. Follow `paper/AGENTS.md` and the repository's `AGENTS.md` (branches,
`#112 - ...` commit subjects, no Claude attribution lines in commits).

## Word budget

Bioinformatics allows about 5,000 words of main text. After the Introduction and Methods
rewrites the main text was 4,898 words (Methods grew from about 500 to 966, the Introduction
from about 440 to about 780). After the Results restructure (2026-11-26) it is 3,811, with
Results at about 1,090 and the Discussion at about 650, so the budget is no longer tight.
Keep an eye on `vibepaper wordcount` as Results and the Discussion are rewritten. The detail
of Methods was moved to Supplementary Methods for the same reason.

## Current state

- **Introduction: done.** All 11 paragraphs are accepted and in `introduction.md`, with
  claims updated.
- **Methods: done.** All 6 paragraphs (M1 to M6) are accepted and in `methods.md`. M3 now
  leads with `fix_hgvs()` (the recommended entry point) and ends on `clean_hgvs()` (the pure
  string function); the MANE gene lookup moved from M5 into M3. `methods.claims.toml` has one
  claim (the remote UTA subset). The old Methods detail is in `supplementary.md` under
  Supplementary Methods, as "... (detail)" subsections, which stay AI text for now, to be
  lightly edited later. Some of it overlaps the older supplementary subsections, which needs
  tidying.
  **Results in progress.** Restructured on 2026-11-26 to seven paragraphs (R1, R2a, R2b, R3,
  R4a, R4b, R5, matching the `## R1` to `## R5` headings), each with one central result about
  cdot. The restructure is applied in the files: details not about cdot moved to a new
  Supplementary Results section in `supplementary.md`, the LOVD/VariantValidator/Mutalyzer
  numbers moved from the Discussion into R4b, the two-corpora point moved to the start of the
  Discussion (for D1), and the claims sidecars moved with their sentences. R1 is his accepted
  text. R2a to R5 still hold the old AI wording, cut, as placeholders.
  **Next:** he dictates R2a to R5 from `results_briefs.md` (one brief per paragraph, numbers
  filled in), in one go if he likes. Record each dictation in `outline.md`, cut and join it,
  then give feedback on the set and swap the accepted text into `results.md`.
- **Informality pass: Introduction and Methods done.** The written text keeps spoken register ("you",
  contractions, "So"/"And"/"But" openers, "a bit of", "many, many"). Claude gives three
  informal sentences at a time with the reason, he re-dictates them, and they are swapped in
  by cutting only. He has approved small mechanical fixes on top (deleting a word, "their" →
  "its", expanding a contraction, cutting a repeated sentence, cutting ", it" to fix a
  dangling modifier). Dictations are recorded in `outline.md` under each entry.
  - **Done:** the whole Introduction and the whole of Methods (M3 to M6; M1 and M2 were
    already formal). He keeps the conversational I1b opener ("One of the surprising things
    we've seen").
  - **Left for his paper edit (noted, not fixed):** I1a "convert them" after the singular
    "The c. form"; I3 "went back many years" then "many years old"; I3 list not parallel
    ("cleaning … collecting … and conversion"); I4 "server" twice.
  - Do another pass over Results and the Discussion once they are written. After the passes
    he prints the paper and edits it by hand.
- **Discussion, Abstract: not started.** Their outline entries describe the target
  structure. The LOVD/VariantValidator/Mutalyzer move into Results is done; the Discussion's
  first paragraph is the moved two-corpora point, a placeholder for D1. Do the Abstract last.
- **Supplementary Results** is AI text moved from Results, joined lightly. Tidy it with the
  rest of the supplement later.
- **Figure caption, table captions, supplement:** optional, later.

## Tooling notes

- Don't run `vibepaper wrap` on a section: it splits long lines without rejoining short
  ones, which leaves ragged one-word lines. Wrap edits by hand.

## Open questions for the author

- I3 says "cdot was split out of the Australian Genomics Shariant project". It was strictly
  extracted from VariantGrid, the software Shariant runs on. He hasn't decided whether to
  change it.
- I1b ("One of the surprising things we've seen...") and the I3 search-box paragraph make
  overlapping observations. He hasn't decided whether to trim one.
- `todo.md` lists the remaining before-submission items.
