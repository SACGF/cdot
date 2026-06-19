# cdot Paper — Working Directory

Target journal: **Bioinformatics (Oxford)** — DECIDED 2026-06-17
Article type: **Original Paper** (full article; not an Application Note). See
`claude/paper_plan.md`.
Structure: Abstract · Introduction · Methods · Results · Discussion · Figures.

## Files

| File | Purpose |
|------|---------|
| `bioinformatics_submission_guide.md` | Journal requirements |
| `abstract.md` | Abstract (Summary / Availability / Contact) |
| `introduction.md` | Background and problem statement (~400 words) |
| `methods.md` | Methods: data generation, JSON format, `clean_hgvs()`, access/clients, canonical selection |
| `results.md` | Results: coverage + ClinVar resolution, cleaning recovery + taxonomy, version fallback, speed, gap/T2T |
| `discussion.md` | Positioning, limitations, future (~300 words) |
| `figures.md` | Figure descriptions (target: 2 main figures) |
| `references.md` | Formatted reference list |
| `supplementary.md` | Supplementary tables and figures |

## Building the paper

Two build modes. Both render the Markdown sections to `output/cdot_<date>.docx`
(+ `.md`, + a supplementary doc) with [vibepaper](https://pypi.org/project/vibepaper/);
install the tooling with `pip install -e '.[paper]'` (adds `vibepaper`, `snakemake`,
and `weasyprint` for the optional PDF output).

Numbers in the text are not hard-coded — they are `{{ fact }}` templates filled from
one-row CSVs ("facts", one per column). The render always reads facts from the
**facts dir** (`output/facts/`, set by `facts_dir` in `paper.toml`). The two modes
differ only in **how that dir gets populated**:

- **full** regenerates the facts there from production data / benchmark runs;
- **quick** copies the committed *stored results* from `paper/empirical_results/` into it.

`paper/empirical_results/` is the checked-in snapshot of measured results (it is *not*
the facts dir — "facts" are the generated values). The full build refreshes it after a
successful run.

### Quick build (default) — seconds, no data needed

Copies `paper/empirical_results/*.csv` → `output/facts/`, then renders. No production
JSON, no private corpus, no network, no UTA.

```bash
# from the repo root, with the venv set up; default target is the quick build
snakemake -s paper/Snakefile --cores 1
# equivalent by hand:
cp paper/empirical_results/*.csv output/facts/ && \
    .venv/bin/vibepaper build --config paper/paper.toml
```

Use this for editing prose, checking layout, or any rebuild where the numbers
haven't changed.

### PDF output

By default the build produces `.docx` (+ `.md`) only. To also render
`output/cdot_<date>.pdf`, use the dedicated `pdf` rule (a quick build plus
vibepaper's `--pdf`):

```bash
snakemake -s paper/Snakefile pdf --cores 1
```

Either form adds `--pdf`, which needs the **weasyprint** toolchain (pulled in by
`pip install -e '.[paper]'`). Probe whether it's usable first with:

```bash
.venv/bin/vibepaper build --config paper/paper.toml --is-pdf-available
```

To attach a PDF to a **full** build (or the quick build by hand), pass
`--config pdf=true`:

```bash
snakemake -s paper/Snakefile full --config pdf=true --cores 1 \
    --config data_dir=/path/to/cdot_data ...
```

### Full build — regenerates every fact, then renders

Recomputes the reproducible (Tier-1) facts from scratch, refreshes the
`paper/empirical_results/` snapshot, then renders. Slow, and needs the inputs each fact requires.

```bash
snakemake -s paper/Snakefile full --cores 1 \
    --config data_dir=/path/to/cdot_data \
             clinvar_hgvs_csv=/path/to/clinvar_pairs.tsv \
             uta_uri=postgresql://uta:uta@localhost:5432/uta/uta_20241220
```

What each fact comes from (see the rules in `paper/Snakefile`):

| Fact CSV | Source | Needs |
|---|---|---|
| `coverage.csv` | `paper/scripts/compute_coverage.py` over the release JSON.gz files | `data_dir` |
| `sources.csv` | `paper/scripts/compute_sources.py` over `cdot_transcripts.yaml` | committed (none) |
| `cleaning.csv` | `paper/scripts/inject_and_clean.py` (Tier-1 injection benchmark) | committed test data (none) |
| `benchmark.csv` | measured throughput (`paper/scripts/benchmark_resolution.py`, `paper/scripts/compare_providers.py`) | recorded as frozen measurements in the rule |
| `clinvar.csv` | `paper/scripts/benchmark_resolution.py` (cdot vs UTA resolution) | `data_dir`, `uta_uri`, SeqRepo |
| `literature.csv` | static constants from published papers | none |

The full-scale ClinVar throughput runs take ~1.5 h each — see `claude/benchmark_plan.md`.

### Two-tier facts and the private corpus

- **Tier 1 (reproducible)** lives in the fact CSVs above and regenerates from public
  data committed here.
- **Tier 2 (production validation, not reproducible)** — the cleaning rescue rate
  (91.5% → 96.4%), the per-fix rescue distribution (Results Table 1), and the residual
  error taxonomy — comes from the private `cdot_private` corpus and is written into
  `results.md` as **literal frozen constants**, not regenerable facts. No corpus string
  ever enters this repo. When the corpus is re-analysed (`cdot_private/analyze_cleaning.py`),
  transcribe the new aggregate numbers into `results.md` by hand. These are flagged
  `[Tier 2]` in the text.

To refresh the committed snapshot without a full data run (e.g. after editing one
analysis script), regenerate that one CSV into `output/facts/` and copy it into
`paper/empirical_results/`.

## Key resources

- `claude/paper_thoughts.md` — framing, scope decisions, pre-submission checklist
- `literature/literature_review.md` — literature synthesis with all references
- `claude/research.md` — deep technical notes on the codebase
