# To do

Open work on the paper. Delete an item when it is done. How to work on the paper is in
`AGENTS.md`.

## Review

- [ ] **Read the rewritten R5 sentence** ("only X% of RefSeq accessions ever drift ...").
  It now also cites the Ensembl share and the share of drifting pairs held by the top 5%
  of accessions (`version_stability.ensembl_accessions_drift_pct`,
  `*_top5_drift_share_pct`), the facts "a short, identifiable list" rests on (PR #245).

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
