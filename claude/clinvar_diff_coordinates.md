# ClinVar coordinate differences (VCF-comparison residual)

Analysis of where cdot's resolved genomic coordinate disagrees with ClinVar, when the
ClinVar c.HGVS is resolved through cdot and compared as a **VCF coordinate**
(CHROM/POS/REF/ALT) rather than as a g.HGVS string.

## Why VCF comparison

The original benchmark compared cdot's `c_to_g` output as a g.HGVS *string* against
ClinVar's CLNHGVS. That mis-scores representation differences (3'-shift, equivalent
del/dup spellings) as errors, and cannot even parse the ClinVar variants whose CLNHGVS
uses tandem-repeat `ref[N]` or identity `=` notation (about 17.5k of the original
mismatch bucket). Comparing VCF tuples instead (cdot `c_to_g` -> biocommons Babelfish
`hgvs_to_vcf`, parsimony-normalised, vs the ClinVar VCF's own CHROM/POS/REF/ALT columns)
is representation-robust and defined for every variant.

Run: `paper/scripts/resolve_clinvar_pass.py` on the VCF-format pairs
(`build_clinvar_pairs.py` now emits the VCF ground truth), GRCh38, cdot 0.2.33 RefSeq +
genome FASTA. 4,382,308 pairs.

## Result

| bucket | count | % |
|--------|-------|---|
| correct (matched) | 4,354,691 | 99.37 |
| incorrect | 1,766 | 0.04 |
| no_data | 3,198 | 0.07 |
| error (input unresolvable) | 22,653 | 0.52 |

Of variants cdot resolved, 99.96% match ClinVar's VCF coordinate exactly. Versus the
g.HGVS-string comparison the incorrect bucket collapses from 0.51% to 0.04%, because the
representation-difference and unparseable-notation buckets are no longer counted as
errors.

## Where the residual lives

The non-matching rows (and `error`/`no_data`) are written to
`output/clinvar_pass/vcf_errors.csv` (columns: bucket, category, c_hgvs, tx, version,
gt_chrom, gt_pos, gt_ref, gt_alt, clnhgvs, cdot_vcf). The 1,766 `incorrect`:

| category | count | what it is |
|----------|-------|------------|
| `indel_mismatch` | 1,116 | indel representation / left-alignment differences in repeat regions (471 transcripts; NM_000271 alone 100) |
| `snv_diff_pos` | 406 | SNV whose position differs; concentrated in 42 transcripts (NM_001164277=177, NM_000854.3, NM_015691, ...) |
| `babelfish_megaallele` | 178 | Babelfish left-shuffles an insertion through a repeat into a huge ref/alt; a normaliser artifact, not a cdot coordinate |
| `identity_or_symbolic` | 39 | ClinVar identity (`=`) / symbolic alleles with no comparable VCF ALT |
| `ambiguity_code_allele` | 27 | IUPAC ambiguity ALTs (G>H, G>Y, G>R); position matches, only the degenerate base differs |

Key observation: `snv_diff_pos` is **not** random per-variant error. It clusters by
transcript with a constant per-transcript offset (e.g. every NM_000854.3 variant ~24 kb
off, every NM_001001548.3 variant exactly 3,655 bp off), which means cdot aligned the
*whole transcript* to a different genomic locus than ClinVar. These are
multi-mapping / paralog / copy-number genes (NM_000854 is GSTT-family, a known
copy-number region) where a transcript legitimately maps to more than one place, so cdot
and ClinVar can pick different copies. This is the same class of "same transcript,
different genomic placement" as the version re-placement blocklist (see
`docs/transcript_version_safety.md`), but here within a single version.

So genuine cdot coordinate errors are a small fraction of even the 1,766: most are
paralog/copy-number placement choices, indel/insertion representation, ambiguity codes,
or release-skew ground truth (the variant_summary and VCF are slightly different ClinVar
releases, so a few ALLELEIDs join to a different variant).

## UTA isolation: is it cdot's data?

To separate a cdot-specific problem from a shared-engine or ClinVar-ground-truth one, the
1,766 `incorrect` c.HGVS were re-resolved through **UTA** (a different transcript backend,
the identical biocommons engine + Babelfish, sequence from SeqRepo) and each compared to
both cdot's call and ClinVar's VCF
(`paper/scripts/exploratory/isolate_uta_vs_cdot.py`, needs the local UTA + SeqRepo;
per-variant verdicts in `output/clinvar_pass/uta_vs_cdot.csv`).

| UTA verdict | count | meaning |
|-------------|-------|---------|
| `uta_matches_cdot` | 1,674 (94.8%) | UTA == cdot, both differ from ClinVar -> not cdot-specific |
| `uta_no_data` | 76 | UTA lacks the transcript/version |
| `uta_matches_clinvar` | 10 | UTA == ClinVar, cdot differs -> genuinely cdot-specific |
| `uta_differs_both` | 6 | three-way (all identity/symbolic) |

By category, `snv_diff_pos` is **0%** cdot-specific (334 match cdot, 72 UTA no_data, 0 match
ClinVar): UTA places the paralog/copy-number clusters exactly where cdot does, so the
offset is ClinVar mapping a multi-copy gene to a different locus, not a cdot error. The 10
cdot-specific cases are all indels and are not wrong exon mappings either:

- **Sequence-fetcher N-padding** (NM_002111.8, NM_152503.8): cdot's call contains `N`s
  because this benchmark fed cdot the `FastaSeqFetcher`, which rebuilds transcript sequence
  from the genome and pads gapped regions with `N`, while UTA used SeqRepo. A setup
  artifact; would resolve if cdot used a transcript SeqRepo.
- **Repeat-region left-alignment** (NM_030582.4 x5 in a CCCCCC/GGCCCC tract ~15 bp off;
  NM_001164277.2, NM_001417890.1, NM_032790.3 1-3 bp off): cdot normalises the indel to a
  slightly different position than UTA within a low-complexity tract. Representation, not a
  wrong coordinate.

Conclusion: on this measure cdot's transcript coordinates are as accurate as UTA's. None of
the 1,766 differences is a genuine cdot exon-mapping error; they are shared-engine
representation, ClinVar paralog placement, or (for 10) seqfetcher / left-alignment
artifacts.

## Caveats

- `error` (22,653) is mostly input the c.HGVS engine cannot resolve (c.HGVS that itself
  uses tandem-repeat notation, ambiguity codes that fail conversion, malformed input),
  not a coordinate disagreement.
- The pairs were built from `variant_summary.txt.gz` (newer) joined to
  `clinvar_latest.vcf.gz` 2026-03-21 by ALLELEID; the slight release skew contributes a
  few wrong-ground-truth rows in `incorrect`.

## Files

- Full per-variant table (not in source control): `output/clinvar_pass/refseq_full_vcf.csv`
- Residual / errors (not in source control): `output/clinvar_pass/vcf_errors.csv`
- Headline facts (committed): `paper/empirical_results/clinvar_vcf.csv`
