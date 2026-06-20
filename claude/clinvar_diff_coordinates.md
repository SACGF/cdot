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
