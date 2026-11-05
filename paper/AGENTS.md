# paper/AGENTS.md

Notes for AI coding agents working on the manuscript. The repo-wide rules in the top-level
`AGENTS.md` still apply. `README.md` in this directory covers the build, the fact CSVs and the
public/private fact model.

## Conventions

- **No private-corpus strings.** The cleaning rescue rates, the per-fix distribution and the
  injection weights in `scripts/inject_and_clean.py` come from the private `cdot_private`
  corpus. Only aggregate numbers are cited, as frozen constants. Synthesise any example from
  the public ones (BRCA2 `NM_000059.4`, RUNX1 `NM_001754.5`).
- **LLM classification tasks emit short class names** (eg `Truncated`, `No reference`,
  `Bad accession`), not long ALL_CAPS taxonomy keys. Do it at the source so the tables need no
  post-processing.
- **UTA coverage is measured by membership** of the exact accession.version in UTA's
  alignment dump (`ac` column of `Homo_sapiens_GRCh38_UTA_20241220.csv`). The cdot JSON
  converted from UTA (`cdot-*_uta.*.json.gz`) under-counts badly and is not a fair proxy.
- Analysis write-ups that back a fact or a PROVENANCE entry go in `notes/`.

## After the next client re-run (0.2.37, #235)

Every resolution-side fact is re-measured on the released client (VCF pass, safe-version
validation, submitted pass, residual categorisation and positions, Table 1, full-scale
throughput). The cleaning, LOVD, VV/Mutalyzer, coverage and version-stability facts do not
move. Prose written against the old numbers to re-check once the CSVs are refreshed:

- `results.md` R2, full-scale paragraph: the residual composition sentence ("dominated by
  indel-representation differences and by paralog and copy-number transcripts"). #235 adds a
  third class: about 230 coordinate differences on seven transcripts (SHANK3 `NM_001372044.2`
  and others) with transcript bases unaligned between two exons, most of them identical to
  UTA's answer. Say so here or in Table S4: cdot and UTA agree, ClinVar distributes the
  unaligned bases differently, so they are alignment-choice differences, not cdot errors.
- `results.md` R2, submitted corpus: headline rates, cdot-only count, the era sentence (cdot
  per era, UTA lowest on 2008-2015 and lower on 2021-2026 than 2016-2020), and the residual
  sentence (about 896 submitted rows on the #235 transcripts move out of `no_data`/`error`).
- `results.md` R3: "within 2%" for REST-prefetch vs local JSON, "nearly three orders of
  magnitude" (also in `abstract.md`), "about four times a locally installed UTA" (abstract),
  the full-scale wall time and hot/cold throughput, and "extrapolates to over a year".
- `results.md` R5, last paragraph: safe-substitution counts and "0 coordinate changes".
- `abstract.md`: the submitted-corpus cdot vs UTA percentages and the throughput figures.
- `discussion.md`, limitations: consider adding biocommons/hgvs's handling of intronic
  offsets anchored inside an alignment gap (14 ClinVar rows raise `IndexError` on 0.2.37,
  UTA takes the same path; davmlaw/dave-admin#253) as an engine limitation.
- `methods.md` JSON format: once released, the per-exon gap string also covers transcript
  bases that align to neither flanking exon (#235); the example release filename in the
  `JSONDataProvider` snippet if a data release ships with the client.
- `supplementary.md` Table S4: the version line ("cdot 0.2.36 with data release 0.2.36"),
  Table S3 heading, the residual breakdown and positional paragraphs; Table S6 "Counts are
  for cdot 0.2.36" (cleaning did not change, so counts stay); Table S8 counts.
- Expected Table S4 movement: `no_data` down about 2,100, matched up about 1,800, incorrect up
  about 230, 14 new `error` rows (the `IndexError` cases).

## Before submission

- Author list and corresponding author. Core: Dave Lawrence and the Shariant team who helped
  debug HGVS issues. Consider inviting a biocommons/hgvs maintainer and a clinical genomics
  collaborator. Acknowledge the NCBI RefSeq, Ensembl and biocommons teams.
- Fill the placeholders: `[Zenodo DOI]` and `[email]` in `abstract.md`, and the `...` rows of
  Tables S1 and S2 in `supplementary.md`.
- The UTA alignment count cited in the abstract, introduction and R1 is Hart 2020's figure
  for the 2021 release; the comparison runs use `uta_20241220`, whose GRCh38 dump alone holds
  148,275 accession versions. Decide whether to re-measure `literature.uta_count` (and so
  `coverage.improvement_fold`) against that release.
- The version-removal ablation of `get_best_transcript_version()` is not reported; if it
  should be, run `benchmark_resolution.py`'s degraded mode and add the counts as a fact.
- HTSeq, Snakemake and PyHGVS have no citation; UTA, RefSeq and Ensembl are cited.
- `bioinformatics.csl` is author-date; the journal wants numbered citations in order of
  appearance (formatting is only required at revision). `rsvg-convert` must be installed for
  the docx to embed the SVG figures.
