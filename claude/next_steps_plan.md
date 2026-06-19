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
  - `canonical.md` — `CanonicalTranscriptSelector`, MANE, resolving gene-only HGVS
  - `rest-api.md` — cdotlib.org endpoints, `RESTDataProvider`
  - `data-format.md` — JSON schema, exon 6-tuple, gap encoding (for tool authors)
  - `data-generation.md` — building your own JSON.gz from RefSeq/Ensembl GTF/GFF3;
    Snakemake pipeline; source priority; adding a genome build or species (#95)

**Note:** `data-generation.md` matters for users with non-human organisms or custom
assemblies who need to run the pipeline themselves.

---

## STRENGTHENS THE PAPER — non-blocking but important

### S2. BRCA2 gap support verification

**Why:** The 32% accuracy claim (Münz 2015) motivates gap support. It is currently cited
as a Tier-2 literature number (`literature.csv`), not reproduced. Reproducing it would
turn it into a Tier-1 result.

**What to do:**
- Run the BRCA2 variant set from Münz et al. through biocommons/hgvs with cdot (gap
  support via `FastaSeqFetcher`) vs without.
- If reproducible, include as a panel in the benchmarks figure.

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

---

## POST-PAPER / FUTURE WORK

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
| B7 | Documentation site (MkDocs) | Medium | Yes — citable URL for paper |
| S2 | BRCA2 gap support verification | Medium | No — but supports key claim |
| S4 | FastaSeqFetcher mismatch detection | Medium | No |
| F1 | Parser refactor | High | No |
| F2 | CDS phase support | Medium | No |

