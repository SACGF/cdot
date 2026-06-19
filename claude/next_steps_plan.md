# cdot — Next Steps Plan (Bioinformatics paper target)

Outstanding work only. Items that have shipped (code + paper Results + committed facts
CSVs) have been removed; see "Completed" at the bottom for the IDs, so nothing is lost.

## Conventions for analysis work
- **LLM classification/assignment tasks: emit SHORT class names** (e.g. `Truncated`,
  `No reference`, `Bad accession`, `Grammar gap`, `Genomic ref in parens`), not long
  ALL_CAPS taxonomy keys. Do it at the source so tables need no post-processing.
- **Never put private-corpus strings in this repo** (tests/comments/docstrings/paper).
  Synthesise from public examples (BRCA2 `NM_000059.4`, RUNX1 `NM_001754.5`).

---

## BLOCKING — must exist before submission

### B7. Documentation site (MkDocs + mkdocs-material)

**Why:** The paper needs a citable docs URL. The GitHub wiki is not versioned with code
and isn't citable. Developers integrating cdot expect discoverable, versioned docs. The
`docs/` directory already holds the source material (schema, examples, fasta notes,
version-safety write-up) but there is no `mkdocs.yml` and nothing is published.

**What to do:**
- Add `mkdocs.yml`; configure the mkdocs-material theme over the existing `docs/`.
- Host on GitHub Pages (or Read the Docs) at a stable URL, ideally `docs.cdotlib.org` or
  `cdotlib.readthedocs.io`.
- Migrate useful wiki content into the structure.
- Structure:
  - `index.md` — what cdot is, one-minute install + example
  - `quickstart.md` — biocommons/hgvs integration, PyHGVS integration, full working code
  - `data-files.md` — how to get JSON.gz files (GitHub releases, `data_release.py`)
  - `clean-hgvs.md` — `clean_hgvs()` and `get_best_transcript_version()` with examples
  - `canonical.md` — gene-symbol HGVS resolution (`resolve_gene_hgvs` / `fix_hgvs`),
    `get_tx_ac_tags_for_gene` tag ranking, MANE, resolving gene-only HGVS (#36)
  - `rest-api.md` — cdotlib.org endpoints, `RESTDataProvider`
  - `data-format.md` — JSON schema, exon 6-tuple, gap encoding (for tool authors)
  - `data-generation.md` — building your own JSON.gz from RefSeq/Ensembl GTF/GFF3;
    Snakemake pipeline; source priority; adding a genome build or species (#95)

**Note:** `data-generation.md` matters for users with non-human organisms or custom
assemblies who need to run the pipeline themselves.

---

## STRENGTHENS THE PAPER — non-blocking but important

### P0. Paper/plan accuracy fixes before sending out (do first)

**Why:** Found while stripping author-notes; these are wrong-as-written, not future work.

**What to do:**
- **`CanonicalTranscriptSelector` does not exist.** `paper/methods.md` ("Canonical
  transcript selection") and `paper/abstract.md` name a class that is not in the code.
  #36 *is* implemented, but as `cdot.hgvs.gene_hgvs.resolve_gene_hgvs` / `fix_hgvs` and
  the data-provider method `get_tx_ac_tags_for_gene`. Rename in the paper to the real
  API (or describe the capability without inventing a class name).
- **Orphaned `Munz2015` citation.** The gap-accuracy (Münz 2015) framing was removed from
  the paper body per feedback (see S2), but `references.bib` / `references.md` still carry
  the entry, cited nowhere. Remove it, or re-cite it somewhere legitimate.
- `paper/notes.md` still has a stale unchecked `#36 ... must be implemented before
  submission` line; #36 is done (see CHANGELOG). Update or drop that scratch note.

### S4. FastaSeqFetcher mismatch detection and documentation (#84, #55)

**Why:** FastaSeqFetcher is the offline sequence path; if it silently returns wrong
sequences for gap-flagged transcripts it undermines the correctness story.

**What to do:**
- Use the stored `note` field to detect RefSeq-flagged indel transcripts.
- Emit a warning (not an error) when `ExonsFromGenomeFastaSeqFetcher` is used for such a
  transcript.
- Add a docs page: when to use FastaSeqFetcher, what its limitations are.
- Investigate `NM_000399.3` discrepancy (#55) — cdot bug or genuine RefSeq/genome
  divergence?

### A1. Single ClinVar parse/resolution pass (shared infra for R5b + R6)

**Why:** Several paper scripts each walk all of ClinVar (`build_clinvar_pairs.py`,
`benchmark_resolution.py`, `compute_coverage.py`, `compute_benchmark.py`). R5b and R6
below both need another full resolution pass. Adding two more is wasteful and risks the
numbers drifting out of sync between sections. One pass that emits a per-variant results
table, consumed by every downstream analysis, is faster and guarantees consistency.

**What to do:**
- Build one resolution pass over the ClinVar (g.HGVS, c.HGVS) pairs that writes a
  per-variant results table (parquet/CSV): `g_hgvs, c_hgvs, tx, version, bucket
  (correct/incorrect/no_data/error), converted_g, fix_codes`.
- Point R5b, R6, coverage, and benchmark at that artifact rather than re-resolving.
- Scope the refactor to what R5b/R6 touch; do NOT big-bang rewrite all four scripts.
- Constraint (CLAUDE.md): never run against full datasets in dev — use
  `tests/test_data/clinvar_hgvs/` for verification; the full pass is a dedicated run.

### R5b. Empirical ClinVar validation of safe version substitution

**Why:** R5 (`compute_version_stability.py` + `cdot/hgvs/version_safety.py`) proves
*structurally* that a version bump is coordinate-preserving iff the intrinsic CDS
structure is unchanged. R5b is the empirical complement: demonstrate end-to-end on real
variants that every substitution our check calls "safe" actually preserves the genomic
coordinate. This turns the Tier-1 structural claim into an observed zero-error rate and
is the strongest reviewer defense of the safe-fallback feature. (The benchmark's degraded
`absent_version` mode forces `.99` fallbacks but does NOT isolate the safe-identified
subset, so it cannot make this claim.)

**What to do:**
- From the A1 results table, for each ClinVar c.HGVS where another version exists and
  `is_version_substitution_safe` calls the V→W substitution safe, resolve under W and
  confirm the g.HGVS is byte-identical to resolving under V.
- Headline fact: `N safe substitutions, K coordinate changes` (expect K = 0). Any K > 0
  is a real bug to chase. Emit to `output/facts/`.

### R6. Categorize the ClinVar genomic-match failures (the 1.2%)

**Why:** 98.8% of c→g conversions reproduce ClinVar's CLNHGVS exactly. A reviewer will
ask what the residual 1.2% are. Strong prior: most "incorrect" cases are *representation*
differences (3'-shift vs ClinVar normalization, equivalent del/dup/inv spellings,
ref-base disagreements), not cdot coordinate errors. Categorizing likely shows true
accuracy is well above 98.8%, which strengthens the paper rather than weakening it.

**What to do:**
- From the A1 results table, take the `incorrect` bucket and re-normalize both the cdot
  conversion and the ClinVar g.HGVS, then bucket each: normalization-equivalent /
  ref-base / gap-related / build / genuine-disagreement.
- Use deterministic rules for the bucketing (structured equivalences), not an LLM.
- Note: `investigate_fails.py` is a stale PyHGVS/GRCh37 stub — replace it, don't extend.
- Report the breakdown and the corrected "genuine error" rate as a fact.

### R7. Versionless transcript resolution in `get_best_transcript_version` (new issue)

**Why:** A client / code feature (changelog-worthy), not paper analysis, but motivated by
the same safety machinery. Today `resolve_transcript_version` no-ops on a versionless
accession (`_parse_versioned_transcript` returns `None`). Search-bar users routinely omit
the version (`NM_000059:c.…`). When there is only one available version, or all available
versions are structurally equivalent (same intrinsic CDS structure ⟹ coordinate-safe),
we can pick one with a WARNING.

**What to do:**
- Extend the cleaning/resolution path so a versionless transcript resolves to a concrete
  version when it is unambiguous (single version, or all versions structurally equal via
  `intrinsic_cds_structure`).
- **Ambiguity policy (decided): REFUSE.** When available versions are NOT structurally
  equivalent, return an ERROR fix and leave the string unchanged (mirrors the `REFUSE`
  default of `on_unsafe_version`). No silent pick.
- File a GitHub issue; add a `CHANGELOG.md` `[unreleased]` entry referencing it. Lands in
  `cdot/hgvs/clean.py` / `gene_hgvs.py`. Document under the planned `clean-hgvs.md` (B7).

---

## POST-PAPER / FUTURE WORK

### S2 (demoted). BRCA2 gap support verification — no longer a paper claim

**Status change:** The paper deliberately removed the gap-correctness framing (per
feedback): cdot *stores* per-exon gap data, but applying it is a downstream-library
concern, not promoted as a cdot result (see `results.md` R7, and the now-removed
intro/methods notes). So reproducing the Münz 32% is no longer a paper-strengthening
item — it would re-introduce a claim the paper chose not to make. The correctness-of-gap
concern that *does* still matter (does `FastaSeqFetcher` return wrong sequence for
gap-flagged transcripts) is covered by **S4**. Keep this only as optional post-paper
validation; do not add a benchmark-figure panel for it.

### F3. Automated data-release regeneration (stated as planned in Discussion)

**Why:** `discussion.md` now states that automating regeneration, so each new RefSeq and
Ensembl release is ingested and published as a data release without manual intervention,
is a planned improvement. The plan should track what the paper promises.

**What to do:**
- Wire the Snakemake pipeline + release publishing into a scheduled/CI job that detects a
  new RefSeq/Ensembl release, regenerates the JSON.gz, and publishes a GitHub data release.
- Not blocking for submission; it backs a forward-looking sentence, not a result.

### F1. GFF/GTF parser refactor: split by source not format (#101)

The GTFParser/GFF3Parser split follows file format but logic differs by source (RefSeq vs
Ensembl). This breaks Ensembl GFF3 protein version extraction. A proper fix splits into
`RefSeqParser` and `EnsemblParser`. High effort; needs CI first. Not needed for the paper.

### F2. CDS phase / ribosomal frameshift support (#76)

Genes with −1 frameshifting (e.g. `NM_015068.3` / PEG10) have incorrect protein
reconstruction from cdot data. Fix requires GFF3 phase column parsing and either a 7th
exon tuple element or migration to exon dicts (#78). Biocommons HGVS does not yet consume
phase, so impact is limited to non-HGVS consumers (pyreference). Worth a mention in the
paper as a known limitation.

---

## Optional plumbing

- **Facts-driven R2/R3 tables.** R2 Table 1 and the R3 taxonomy table are hardcoded
  literal markdown (Tier-2 frozen, hand-transcribed — repo convention), NOT vibepaper
  facts. To make future refreshes a one-line CSV edit, promote them to a `residual.csv` /
  extended `cleaning.csv` + Snakefile rule. Optional.

---

## Summary Table

| ID | Task | Effort | Blocks submission |
|----|------|--------|-------------------|
| P0 | Paper accuracy fixes (wrong class name, orphaned cite) | Low | Yes — before sending out |
| B7 | Documentation site (MkDocs) | Medium | Yes — citable URL for paper |
| S4 | FastaSeqFetcher mismatch detection | Medium | No |
| A1 | Single ClinVar parse/resolution pass | Medium | No — shared infra for R5b/R6 |
| R5b | Empirical validation of safe version bumps | Medium | No — but strong reviewer defense |
| R6 | Categorize the 1.2% genomic-match failures | Low | No — likely raises true accuracy |
| R7 | Versionless transcript resolution (REFUSE on ambiguity) | Medium | No — client feature + changelog |
| S2 | BRCA2 gap verification (demoted, post-paper) | Medium | No — dropped from paper framing |
| F3 | Automated data-release regeneration | Medium | No — backs a Discussion sentence |
| F1 | Parser refactor | High | No |
| F2 | CDS phase support | Medium | No |

