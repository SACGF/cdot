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

## When he asks for a change to data or facts

Do it properly, in the facts pipeline, not in the prose. For example, the UTA count is now
measured by `paper/scripts/compute_coverage.py` from the uta_20241220 dumps, as
`coverage.uta_count`. Follow `paper/AGENTS.md` and the repository's `AGENTS.md` (branches,
`#112 - ...` commit subjects, no Claude attribution lines in commits).

## Word budget

Bioinformatics allows about 5,000 words of main text. After the Introduction and Methods
rewrites, the main text is 4,898 words (Methods grew from about 500 to 966). The Introduction grew from about 440 words to
about 780. Results and the Discussion must come down as they are rewritten: keep an eye on
`vibepaper wordcount`.
The detail of Methods was moved to Supplementary Methods for the same reason.

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
  **Next: present R1a and R1b** (transcript coverage, alignment gaps).
- **Results, Discussion, Abstract: not started.** Their outline entries describe the target
  structure, including moves not yet made in the files. The biggest is that the
  LOVD/VariantValidator/Mutalyzer numbers move from the Discussion into Results R4.6. Make
  each move when that section's accepted text goes in. Do the Abstract last.
- **Figure caption, table captions, supplement:** optional, later.

## Open questions for the author

- I3 says "cdot was split out of the Australian Genomics Shariant project". It was strictly
  extracted from VariantGrid, the software Shariant runs on. He hasn't decided whether to
  change it.
- I1b ("One of the surprising things we've seen...") and the I3 search-box paragraph make
  overlapping observations. He hasn't decided whether to trim one.
- `todo.md` lists the remaining before-submission items.
