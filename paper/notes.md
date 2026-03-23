# Working Notes

---

## Article type decision

**Application Note** (≤4 pages / ~2,600 words, or ~2,000 words + 1 figure).

Rationale: cdot is a tool + data resource — the classic Application Note scope. The Hart 2015 biocommons/hgvs paper (the tool cdot integrates with) was itself a Bioinformatics Application Note. A tight 4-page paper with strong numbers is more impactful than a padded 7-page paper.

---

## Numbers to fill in before submission

- [ ] Total transcript count (cdot combined, all builds) — from Snakemake summary stats
- [ ] RefSeq GRCh37, GRCh38, T2T counts separately
- [ ] Ensembl GRCh37, GRCh38, T2T counts separately
- [ ] UTA count (run `SELECT COUNT(*) FROM tx_info` against uta_20210129)
- [ ] Exact multiplier: cdot / UTA (currently stated as ~9×)
- [ ] Speed benchmark: cdot local JSON, cdot REST, UTA remote, UTA local (transcripts/s)
- [ ] ClinVar benchmark: N variants tested, cdot resolution %, UTA resolution %
- [ ] T2T unique transcript count (genes not in GRCh37/38)
- [ ] JSON file sizes (compressed and uncompressed, GRCh38 RefSeq)
- [ ] Load time for GRCh38 RefSeq JSON.gz

---

## Word count tracking

| Section | Target | Status |
|---------|--------|--------|
| Abstract | ~150 | draft |
| Introduction | ~400 | draft |
| Implementation | ~1,200 | draft |
| Discussion | ~300 | draft |
| **Total** | **~2,050** | **draft** |
| Page budget (4 pages) | ~2,600 words | — |
| Remaining budget | ~550 words | — |

*~550 words of slack. Keep for figure legends, captions, and any result detail that needs expanding once real numbers are in.*

---

## Open questions

- [ ] Does VariantValidator support T2T? Confirm before claiming "first HGVS resource"
- [ ] Does TARK support T2T? Same check
- [ ] #36 (CanonicalTranscriptSelector) — must be implemented before submission; it's in the abstract
- [ ] #100 — cdotlib.org migration from cdot.cc — URL in paper must be stable at submission
- [ ] Zenodo DOI for data files — register before submission
- [ ] Corresponding author and co-author list
- [ ] Cover letter: mention production use in VariantGrid/Shariant

---

## Decisions made

- Application Note, not Original Paper
- cdot_rest = same paper; pyreference = 1 paragraph in Discussion (cut for word count)
- SeqRepo analogy: use explicitly
- Speed benchmark goes to Supplementary Figure S1 to save main figure space for Figure 1 (coverage + ClinVar)
- pyreference mention cut from discussion draft (word count); can add back if under budget

---

## Pre-submission code checklist

- [ ] #36 CanonicalTranscriptSelector implemented
- [ ] #100 cdotlib.org domain live
- [ ] #55 FastaSeqFetcher mismatch bug fixed (affects gap correctness claim)
- [ ] CI running with test suite (#62)
- [ ] Python version bump + build deps (#104)
