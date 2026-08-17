# Talk slides

Source for a 20-minute conference talk on the cdot paper, first built August 2026. The deck
walks the same story as `paper/results.md`: where cdot sits in the Python HGVS stack, how the
data is built, then the five results (coverage, ClinVar resolution, throughput, `clean_hgvs()`
recovery, and safe version fallback), plus two backup slides holding the bands of Figure 1.

This is a one-off. It is deliberately not wired into the paper's Snakemake pipeline, and
nothing else in the repo depends on it. It lives here so it is version controlled rather than
sitting in the ignored `output/` directory.

## Files

| File | Purpose |
|------|---------|
| `make_deck.py` | Builds the whole deck with python-pptx. Slide content, layout helpers, and speaker notes all live here |
| `make_figs.py` | Renders the six matplotlib chart figures (coverage, ClinVar, throughput, cleaning, residual, version safety) |
| `make_figure1.py` | Rasterises `paper/figures/figure1.svg` to slide-ready PNGs: the whole figure and each band separately |

## Building

```bash
pip install python-pptx matplotlib cairosvg pillow

python paper/slides/make_figure1.py    # writes output/talk_figs/figure1*.png
python paper/slides/make_figs.py       # writes the rest of output/talk_figs/
python paper/slides/make_deck.py       # writes output/cdot_talk_20min.pptx
```

Only the scripts are tracked. The figures and the deck they build are artefacts and go to the
ignored `output/` directory, the same way the paper tracks `figure1.svg` but not its
rasterisations. Set `DECK_OUT` to write the deck somewhere other than the default.

## Notes for editing

Numbers in the deck are typed in as literals, copied from `paper/results.md` at the time of
writing. They are not `{{ fact }}` templates like the paper, so if the facts CSVs are
refreshed, the deck does not follow automatically and has to be checked by hand.

`DURATIONS` near the top of `make_deck.py` holds the per-slide speaking budget in seconds, in
slide order, for the content slides only (the backup slides at the end are excluded). Each
`notes()` call pops the next entry, so the running clock in the speaker notes is computed
rather than typed. Adding or removing a slide means adding or removing a `DURATIONS` entry in
the matching position. The script prints the total talk length when it finishes.

Every HGVS string on the slides is synthesised from the public examples used in the paper and
the tests (BRCA1 `NM_007294.4`, BRCA2 `NM_000059.4`). None of them come from the real search
query corpus in `cdot_private`, and none should.
