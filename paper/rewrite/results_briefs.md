# Results briefs: dictate all six in one go

Six paragraphs are left in Results. Each brief gives where it sits, the one-line result it
has to land, the numbers (already looked up, so say them roughly or as "X" and they get
filled in as facts), and what not to repeat. Dictate them in any order, all at once or a
few at a time. Your words are only cut and joined, never added to.

`results.md` already has the new shape. Each of these paragraphs currently holds the old
AI wording, cut down, as a placeholder until your version replaces it. The details that are
not about cdot have moved to Supplementary Results, so a pointer ("Supplementary Results")
is all you need for them.

Already done: **R1** (coverage). You said: 1,954,742 alignments against 301,271 in UTA, 6.5
times as many, from every archived RefSeq release plus Ensembl, which UTA lacks; Ensembl
GRCh38 alone is 854,940, and T2T adds 186,644.

Target: three to five sentences each. Main text is 3,811 of 5,000 words, so there is room,
but short is better.

---

## R2a. cdot reproduces ClinVar's coordinates

**Where you are:** first ClinVar paragraph. R1 said how much cdot covers; this says it is
correct.

**One-line result:** across every c.HGVS in ClinVar, cdot reproduces ClinVar's own genomic
coordinate for 99.95% of the variants it resolves.

**Numbers:**
- every c.HGVS in ClinVar's variant_summary: 4,583,798 pairs, effectively all RefSeq
- 99.5% resolve
- 99.95% of those resolved match ClinVar's VCF coordinate exactly
- the 0.05% left over is mostly indel representation and paralog transcripts that map to
  more than one place (Supplementary Table S4)
- Ensembl: 818 Ensembl pairs, cdot 99.8%, UTA 0%

**Say, in order:** the corpus, the two rates, the residual as one clause with the table
pointer, then Ensembl.

**Don't repeat:** Methods already says how scoring works (against ClinVar's VCF coordinate),
and R1 already says UTA has no Ensembl. Give the Ensembl numbers, not the reason.

**Your dictation:**

---

## R2b. Labs cite old transcript versions, which cdot has and UTA doesn't

**Where you are:** second ClinVar paragraph. R2a used ClinVar's own recomputed names, which
are always at the current version. Here are the strings labs actually submitted.

**One-line result:** most transcript versions labs submitted to ClinVar are superseded, so
cdot resolves 99.3% of submitted strings against 81.9% for UTA. The same holds on Shariant.

**Numbers:**
- the R2a corpus is always the current version, so it is a best case
- the submitted strings (from the ClinVar VCV XML, Methods): 3,245,277 pairs
- 74.7% cite a transcript version that is no longer current
- cdot holds 99.3% of those old versions
- cdot resolves 99.3% vs UTA 81.9% (local UTA, release uta_20241220). UTA has no alignment
  for 17.9% of the cited versions, if you want the "because"
- Shariant [private data]: 36,010 strings from clinical labs, cdot 99.6% vs UTA 93.7%. No
  ground truth for Shariant, so this is a resolution rate, not correctness

**Say, in order:** why the R2a corpus is a best case, the submitted corpus, how old the
versions are, cdot vs UTA, then Shariant in a sentence or two.

**Watch:** 99.3% comes up twice (old versions cdot holds, and strings cdot resolves). Say
them so they can't be confused.

**Your dictation:**

---

## R3. cdot is faster than UTA

**Where you are:** throughput. Table 1 sits just above this paragraph with every number in
it, so the text only needs the comparisons.

**One-line result:** cdot is about 4 times faster than UTA locally and about 830 times
faster remotely, and all of ClinVar takes about 2 hours, where the public UTA would take
over a year.

**Numbers:**
- same 500 ClinVar pairs, same HGVS engine and sequence source for every backend
- local: cdot JSON 741 HGVS/s, local UTA 187 (about a quarter)
- remote: cdot REST 108 HGVS/s, public UTA 0.13 (about 830 times)
- after one `prefetch()` call, REST is within 0.7% of local JSON (746 vs 741)
- a GRCh38 JSON file loads in about 5 seconds
- all of ClinVar (4,583,798 pairs) in about 121 minutes in one process
- the same run against public UTA would take about 408 days

**Say, in order:** local comparison, remote comparison, prefetch, then the all-of-ClinVar
point to finish.

**Don't repeat:** the Introduction already gives public UTA at 0.13 HGVS/s. Use the ratio.

**Your dictation:**

---

## R4a. Cleaning makes most failing real-world queries usable [private data]

**Where you are:** string cleaning. Methods described `fix_hgvs()`; this is how it does on
what people really type. Table 2 (fixes by type) sits just after this paragraph.

**One-line result:** on real search-box queries, cleaning makes two thirds of the failing
strings parseable (91.7% to 97.3%), with no regressions.

**Numbers:**
- 32,671 real search queries from VariantGrid-based platforms that look like HGVS
- the damage: whitespace and odd characters from Word and PDF reports, wrong case,
  punctuation, a protein change pasted on the end
- parseable rises from 91.7% as typed to 97.3% after cleaning
- 1,828 strings rescued, 67% of those that failed
- 0 regressions (nothing that parsed before was broken)
- Table 2: punctuation and whitespace fixes do most of the work
- the rarest fix, restoring a dropped accession prefix (`000059.4` to `NM_000059.4`), needs
  the transcript data, which a string-only checker doesn't have
- optional: what still fails is mostly incomplete input that no repair can invent
  (Supplementary Results)

**Say, in order:** the corpus, the damage, the numbers, then Table 2 and the prefix point.

**Don't repeat:** Methods already lists example fixes and the dropped-prefix example. Here
just say it is the rarest and why it needs data.

**Your dictation:**

---

## R4b. Other tools repair only part of the same errors

**Where you are:** after the cleaning result. These numbers used to be in the Discussion;
the Discussion will keep only what they mean.

**One-line result:** on the same common errors, LOVD repairs 69%, VariantValidator 54% and
Mutalyzer 34%, and none of them changes a valid input.

**Numbers:**
- the injected corpus: 998 clean ClinVar strings, each of the 18 error types `clean_hgvs()`
  fixes injected into them, 3,419 broken strings in total (Methods)
- run through LOVD's checker (v1.2.2), and VariantValidator and Mutalyzer over their public
  APIs (August 2026)
- caveat: these are exactly the errors `clean_hgvs()` was built to fix, so it gets 100% by
  construction. The question is how much of the same ground the other tools cover
- weighted by how often each error type turns up in the real search-box data: LOVD 69%
  (34% unweighted), VariantValidator 54%, Mutalyzer 34% (Supplementary Table S5)
- none of the three changed a valid input

**Say, in order:** the corpus and tools, the caveat, the three numbers, then no false
corrections to finish.

**Your dictation:**

---

## R5. Version substitution: risky without a check, safe with cdot's

**Where you are:** last Results paragraph. Methods M4 described the opt-in substitution and
its CDS-structure check. This shows the risk is real and the check catches it.

**One-line result:** swapping to the next version moves coordinates for about 4% of RefSeq
version pairs, and cdot's check accepted 8.4 million ClinVar substitutions with no
coordinate changes.

**Numbers:**
- optional opener: biocommons HGVS has no version substitution with any data provider
- 19,393 consecutive RefSeq version pairs: 96.0% keep every coding coordinate
- 12,377 Ensembl pairs: 98.0%
- how the drift happens is in Supplementary Results and Figure S1
- the structure check works across builds: the structure is identical on GRCh37, GRCh38
  and T2T for 99.4% of RefSeq and 100% of Ensembl versions
- on ClinVar: 8,411,749 substitutions the check accepted, 8,411,000 of them resolve, 0 give
  a different coordinate

**Say, in order:** the risk, the cross-build point, then the ClinVar validation to finish.

**Don't repeat:** M4 already says what the check compares (CDS length plus coding-exon
lengths) and why it works across builds. Here just give the 99.4%/100%.

**Your dictation:**
