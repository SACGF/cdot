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
- **Wording that depends on a fact is generated or checked**: a word filter (`approx`,
  `within_pct`, `orders_of_magnitude`, ...), an `{% assert %}` on the line above the
  sentence, or `{# lint-ok: reason #}` for a true hedge. A derived value a filter needs
  (ratio, percent difference) is emitted by the script that owns the fact. A judgement
  word ("essentially flat", "nearly uniform") gets one tolerance in `paper.toml`
  `[claims]`, used as `{% assert x is essentially_flat(y) %}`. Keep
  `vibepaper lint --config paper/paper.toml --numbers` and `vibepaper check` clean (CI
  runs both), and after a re-run read `vibepaper diff` before the paper.
  `lint --unused` lists facts the paper never cites; list ones kept on purpose under
  `[lint] unused_ok`.
- **vibepaper 0.10** (`vibepaper~=0.10.0` in the `paper` extra and in CI). To move to a new
  release, follow vibepaper's `UPGRADE_PROMPT.md`.
- **Facts are typed.** Every facts CSV is vertical (`field,value,type`). Scripts write
  them with `write_typed_facts()` from `scripts/facts_io.py`, which types each fact from
  its name (`*_pct` percent, `n_*` count, `*_ratio` ratio, ...), so name a new fact to
  match. In templates use `{{ x_pct | pct(1) }}`, not `{{ x_pct | dp(1) }}%`; a value out of
  its type's range fails the load.
- **Tool and data versions behind a fact come from its provenance sidecar**
  (`empirical_results/<fact>.provenance.json`, `{{ clinvar_vcf._prov.tools.cdot }}`), not
  typed into the prose. A script that writes a facts CSV whose versions the prose cites
  should call `vibepaper.provenance.write_provenance()`; for a frozen fact, rewrite the
  sidecar with that function when the CSV is re-measured. Lint flags a paragraph that mixes
  versions of one tool. Write any other CSV with LF line endings
  (`csv.writer(f, lineterminator="\n")`).
- **Re-run PRs**: paste
  `vibepaper diff --config paper/paper.toml --facts-dir paper/empirical_results --against main --format md`
  into the PR description (`output/facts/` is not committed, so compare the snapshot).

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
