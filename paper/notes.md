# Working Notes

*Decisions made, things to resolve, numbers to fill in, open questions.*

---

## Numbers to fill in before submission

- [ ] Total transcript count (cdot combined, all builds) — from Snakemake summary stats
- [ ] RefSeq GRCh37, GRCh38, T2T counts separately
- [ ] Ensembl GRCh37, GRCh38, T2T counts separately
- [ ] UTA transcript count (run query against uta_20210129)
- [ ] Speed benchmark numbers (local JSON, REST, UTA remote, UTA local)
- [ ] ClinVar benchmark: total variants tested, resolution rates for cdot vs UTA
- [ ] T2T unique transcript count (not in GRCh37/38)
- [ ] JSON file sizes (compressed and uncompressed, by build)
- [ ] MANE Select coverage in cdot GRCh38 (% of protein-coding genes)

---

## Decisions made

**Journal**: Bioinformatics (Oxford), Original Paper (7 pages). Application Note (4 pages) as fallback if content is tight.

**Scope**: cdot + cdot_rest in same paper (same data, different delivery). pyreference mentioned in Discussion as non-HGVS use. biocommons/hgvs client code stays in cdot repo.

**SeqRepo analogy**: Use explicitly. "cdot is to transcript coordinates what SeqRepo is to biological sequences."

**T2T novelty claim**: Check that no other HGVS tool supports T2T before asserting "first". Tools to check: VariantValidator, VEP, TARK, TransVar.

**pyreference**: One paragraph in Discussion, not a separate paper section.

---

## Open questions

- [ ] Is the gap correctness benchmark (BRCA2, with/without gap support) feasible before submission?
- [ ] Should TransVar be included in the comparison table? Need to verify T2T and gap support.
- [ ] Does VariantValidator support T2T? Check their docs/GitHub.
- [ ] Will #36 (CanonicalTranscriptSelector) be implemented before submission? This is the most paper-ready feature.
- [ ] Confirm corresponding author and co-author list with SACGF team.
- [ ] Should we invite a biocommons maintainer as co-author?
- [ ] Zenodo DOI for data files — needs to be registered before submission.
- [ ] cdotlib.org domain — confirm migration from cdot.cc (#100) is complete before final draft.

---

## Word count tracking

*Update as sections are drafted*

| Section | Target words | Current estimate |
|---------|-------------|-----------------|
| Abstract | 200 | draft done |
| Introduction | 900 | draft done |
| Methods | 1200 | draft done |
| Results | 1000 | outline only |
| Discussion | 900 | draft done |
| **Total** | **4200** | **~3200 words drafted** |
| Page budget (7 pages) | ~5000 words | — |
| Remaining budget | ~1800 words | — |

*Page limit is 7 pages including figures. Each figure roughly costs 0.5–1 page of page budget.*
*With 5 figures, text should be ~4,000–4,500 words.*

---

## Revision checklist (post-submission)

- [ ] Respond to all reviewer comments
- [ ] Update transcript counts if data releases have advanced
- [ ] Ensure all URLs in the paper are stable (cdotlib.org domain)
- [ ] Add formatting per Bioinformatics style guide at revision stage

---

## Misc notes

- The `FastaSeqFetcher` mismatch bug (issue #55) should be fixed before submission — it affects the accuracy claim for gap handling
- Cover letter: mention that cdot has been used in production clinical pipelines (VariantGrid/Shariant) — cite the system if published
- The 50% HGVS error rate in Mutalyzer [Lefter 2021] is a striking intro statistic — use it in the very first paragraph
