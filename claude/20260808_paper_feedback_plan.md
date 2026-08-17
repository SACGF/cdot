# Paper feedback plan

Response to the presentation feedback in `claude/feedback.md`. Each item below records what
the paper currently says, what to do about it, and what else it touches. Ordered by item
number; a suggested working order is at the end.

## 1. The Mutalyzer "50%" claim (introduction)

`paper/introduction.md` says Mutalyzer "found ~50% error rates in submitted HGVS
descriptions over five years, many attributable to missing transcript data". Lefter 2021
reports the opposite split: ~50% of submitted descriptions were *correct*, ~41% had a
syntactic or semantic error, and only ~7% could be automatically corrected by Mutalyzer.
The "many attributable to missing transcript data" clause appears unsupported.

To do:

- Rewrite the sentence with the real breakdown (50% correct, 41% error, 7% auto-corrected).
  Change `literature.csv` and `paper/Snakefile` (the `literature` rule) from the single
  `hgvs_error_rate_pct: 50` fact to the three separate facts.
- Expand it into the motivation for R4. An independent, large-scale production stream says
  roughly 4 in 10 human-entered HGVS strings are broken, and the leading checker repairs
  only ~7% of submissions. That is the gap `clean_hgvs()` fills: we rescue 5.1 points of
  the 8.5 that fail, about 60% of the failures. Much stronger framing than the current line.
- Take the numbers and denominators from the paper body (PMC8479679), not the abstract, and
  check whether Mutalyzer's rate is per submission or per unique string, since we report both.

## 2. Figure 1: show UTA in the data providers

UTA appears once in Figure 1, as an input in panel A (alignments absent from any annotation
file). Panel B shows only `JSONDataProvider`, `RESTDataProvider` and `EnsemblTarkDataProvider`.
That undersells the drop-in-replacement claim and leaves the R2/R3 comparator invisible.

To do:

- Add a greyed `UTADataProvider` to `PostgreSQL` node in panel B alongside cdot's providers,
  with both plugging into the same `hgvs.dataproviders.interface`. Grey already means
  external in the existing legend, so no legend change is needed.
- Update the figure legend in `paper/figures.md` to name UTA as the alternative backend.

This makes the "swap only the transcript-data layer" experimental design visible at a glance.

## 3. Literature search on existing HGVS cleaning

The discussion currently describes VariantValidator and Mutalyzer as validators only. Both
in fact perform repair. VariantValidator auto-corrects where it can, including re-mapping
incorrectly reported intron/exon boundary coordinates. Mutalyzer 2 returns a corrected
description for ~7% of submissions. A third tool is missing from the paper entirely: the
LOVD HGVS syntax checker (LOVDnl/HGVS-syntax-checker), which offers ranked,
confidence-scored corrections of invalid descriptions.

To do:

- Search properly: VariantValidator, Mutalyzer 2 and 3, the LOVD syntax checker,
  hgvs-weaver, the ClinGen Allele Registry, and the leniency of the biocommons parser
  itself. Add a short related-work paragraph to `paper/discussion.md`.
- State the differentiator honestly. Those tools are web services and validators that repair
  as a side effect of validation. cdot's contribution is an offline Python function that
  repairs before parsing, reports every change as an inspectable `HGVSFix`, and guarantees
  no regressions. Position on repair scope, auditability and offline use, not on "they do
  not do this".
- Optional and stronger: run the injection corpus (`paper/scripts/inject_and_clean.py`)
  through the LOVD syntax checker for a head-to-head recovery comparison. This is the kind
  of comparison a Bioinformatics reviewer is likely to ask for, so it may be worth doing
  pre-emptively.

Reading:

- Mutalyzer 2 (Lefter 2021): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8479679/
- VariantValidator (Freeman 2018): https://onlinelibrary.wiley.com/doi/full/10.1002/humu.23348
- LOVD HGVS syntax checker: https://github.com/LOVDnl/HGVS-syntax-checker

## 4. ClinVar resolution is too high (highest value item)

Confirmed from `paper/scripts/build_clinvar_pairs.py`: the c.HGVS comes from the `Name`
column of `variant_summary.txt.gz`, which is ClinVar's own recomputed preferred-transcript
name and is always at the current version. So the 99.4% resolved / 98.8% matched figures
measure the easy case, and cannot exercise the historical-depth claim at all. That is why
the interesting version story currently sits in the Tier 2 (private, non-reproducible)
Shariant corpus.

The XML parser exists, but not in `too-many-transcripts`. It is
`../clinvar-hidden-structures/scripts/extract_xml_to_csv.py`, which already walks the VCV
XML, pulls `ClinicalAssertion .//AttributeSet/Attribute[@Type="HGVS"]` values, and has
`_extract_coding_transcript()` to grab the versioned accession the submitting lab actually
used.

To do:

- Add `paper/scripts/build_clinvar_submitted_pairs.py`: SCV-submitted HGVS (as written)
  joined to the VCV VCF coordinate as ground truth, via VariationID / AlleleID.
- Report the version-age distribution. How many submitted strings cite a transcript version
  that is no longer current, and how many cite a version retired from the current RefSeq
  release. This is the number that makes the historical-depth argument concrete.
- Re-run R2 as cdot vs UTA on this corpus. Expect the gap to widen and the absolute numbers
  to drop, which is the point.
- Biggest win in the whole list: submitted SCV strings are also messy, so this yields a
  public, reproducible (Tier 1) counterpart to the Tier 2 cleaning result in R4, plus a
  public residual-error taxonomy. That would demote the private corpus to a confirmation
  rather than the primary evidence for two separate results.
- Keep the existing variant_summary pass as an explicit "current-version ceiling" row, so
  the comparison is stated rather than silently swapped.

Caveats to flag in Methods: the full VCV XML is large, and submitted HGVS is not
deduplicated across SCVs, so the dedup and sampling rules need documenting.

## 5. REST faster than local is a red flag, not a result

`empirical_results/benchmark.csv` has local JSON at 540 to 665 HGVS/s and REST after
prefetch at 731. R3 currently explains this away as run-to-run variance. A reviewer is
unlikely to accept a benchmark whose headline configuration loses to the network.

To do:

- Replace the single-shot numbers with N=5 repeats per configuration over an identical fixed
  HGVS set, reporting median and IQR (or a boxplot as a supplementary figure). This means
  updating `paper/scripts/compute_benchmark.py`, which currently does one run of
  `N_BENCHMARK=500`, and carrying dispersion into the Table 1 cells. Note also that the
  540 to 665 range and the 731 REST figure came from differently sized runs, which is an
  apples-to-oranges problem on its own.
- Investigate two concrete mechanisms before settling on variance as the explanation:
  - `AbstractJSONDataProvider.get_tx_exons` (`cdot/hgvs/dataproviders/json_data_provider.py`)
    rebuilds the full exon dict list on every call, with no memoisation. Both providers pay
    this, and it is the per-variant hot path. An `lru_cache` here is a real, client-visible
    speedup and a genuine changelog entry.
  - GC pressure. Local JSON holds hundreds of thousands of transcript dicts resident, while
    a warmed REST cache holds a few thousand. `gc.freeze()` after load (or `gc.disable()`
    for the duration of the benchmark) is a cheap test. If it explains the inversion, that
    is both a better sentence for the paper and a real fix.
- If REST genuinely matches local after repeats, say so plainly with the dispersion shown,
  and keep the existing explanation (both resolve from an in-memory dict, so the shared
  sequence layer bounds both).

## 6. Fix the dropped `NM_` prefix in `clean_hgvs()`

Agreed, this is the most tractable slice of the residual. The "bad accession" class is 167
queries (14.9% of the residual) and dropped prefix is its headline example. `clean.py`
already has `_op_add_transcript_underscore` (`NM000059` to `NM_000059`), so the missing
sibling is the bare-number case.

To do:

- Measure first. Run the sub-breakdown of the 167 bad-accession residuals through
  `../cdot_private/analyze_cleaning.py` to see how much is genuinely dropped prefix versus
  misplaced or truncated version. Implement rules for whatever dominates.
- The dropped-prefix rule can be disambiguated by the kind letter: a bare
  `\d{6,9}\.\d+:c.` implies `NM_` or `XM_`, and `:n.` implies `NR_` or `XR_`. The NM/XM
  ambiguity is real, so make the rule data-provider aware: try the candidate accessions
  against the loaded provider and apply only when exactly one exists. This fits the existing
  design (reported as an `HGVSFix`, never silent) and demonstrates something no pure
  string-level checker can do, namely that having the transcript data locally enables better
  cleaning.
- Add a new `HGVSFixCode`, plus tests using synthesised public examples only (BRCA2
  `NM_000059.4`, RUNX1 `NM_001754.5`). Nothing sourced from `cdot_private`.
- Add a `CHANGELOG.md` entry under `[unreleased]` / `### Added` referencing #112.
- Re-run the production corpus and regenerate Table 2 and Table S6, since both the rescue
  count and the residual taxonomy shift.

## 7. Quantify positional drift along the transcript

Good hypothesis, and directly computable from data we already have.
`paper/scripts/compute_version_stability.py::_preserved_fraction` already walks the CDS
breakpoint by breakpoint, so binning preserved and total bases by relative CDS position
(deciles) is a small change to an existing loop.

To do:

- Emit a per-decile preservation curve across consecutive version bumps, for RefSeq and
  Ensembl. Prediction: monotonic decline toward the 3' end.
- State an important nuance up front, because it partly contradicts the intuition. The
  current facts show drift is overwhelmingly whole-CDS (RefSeq 3.2% full versus 0.7%
  partial), which is a relocation and therefore position independent. The positional effect
  can only live inside the partial-drift bucket, so the analysis must be conditioned on
  partial drift. Reported honestly this is a better result: most version risk is
  all-or-nothing and detectable, and the position-dependent part is the small partial-drift
  tail.
- Give the 1,766 incorrect ClinVar projections (`empirical_results/clinvar_vcf_residual.csv`)
  the same treatment: bin by relative position in the transcript to see whether errors
  concentrate at the 3' end.
- If the effect is clear, it is a supplementary figure plus one sentence in R5. It also gives
  the version fallback a usable rule of thumb: a 5' coding variant is safer to substitute
  than a 3' one.

## Suggested order

1. Item 4 (submitted-HGVS ClinVar corpus). Largest gain, and it feeds items 6 and 7 with a
   public corpus.
2. Item 5 (benchmark repeats plus the two mechanism checks). Reviewer critical, and may
   produce a genuine performance fix.
3. Items 1 and 3 together (fix the Mutalyzer number, add the repair-tool related work). Both
   cheap, both touch the same motivation paragraph, and item 1 is currently a factual error.
4. Item 6 (`NM_` prefix rule), then regenerate the R4 tables.
5. Item 7 (positional drift), then item 2 (figure). Both small.

Items 1, 2, 3 and 7 are paper only. Items 5 and 6 touch shipped code and need `CHANGELOG.md`
entries. Item 4 is analysis tooling only and does not.
