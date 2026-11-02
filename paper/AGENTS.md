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

## Before submission

- Author list and corresponding author. Core: Dave Lawrence and the Shariant team who helped
  debug HGVS issues. Consider inviting a biocommons/hgvs maintainer and a clinical genomics
  collaborator. Acknowledge the NCBI RefSeq, Ensembl and biocommons teams.
- Fill the placeholders: `[Zenodo DOI]` and `[email]` in `abstract.md`, and the `...` rows of
  Table S1 in `supplementary.md`.
