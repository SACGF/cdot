# Transcript-version "bump safety" experiment (C2)

**Status:** exploration only. Not wired into Snakemake, the paper, or the library.
**Script:** `analysis/transcript_version_drift.py`
**Data:** local cdot 0.2.33 single-build JSON (`/data/cdot_data/{refseq,ensembl}/…GRCh38.json.gz`),
which keeps *all* historical versions of each accession as separate keys.
**Runs quoted below:** `--sample 12000 --seed 1`, GRCh38. RefSeq and Ensembl.
Full console output saved at `/tmp/c2_refseq.GRCh38.txt`, `/tmp/c2_ensembl.GRCh38.txt`.

## The question

cdot exposes an opt-in adjacent-version fallback (`get_best_transcript_version`):
when the requested transcript version is absent from the loaded data, substitute
the nearest available version. **Is that safe?** A version bump is "safe" for a
variant if the variant's coding position `c.X` still maps to the *same* genomic
coordinate after the substitution. If it doesn't, the substitution silently
changes what the variant means.

Sub-questions Dave raised:
1. Can you tell a bump will drift *without* holding the other version?
2. Interpolation: if you hold v2 and v4 and they agree, is v3 safe?
3. Does drift oscillate, or only go one way?
4. Do the rules *combine* (e.g. no-gap **and** bracketed)?
5. Could we precompute odds, look them up at query time, and is that
   clinically acceptable?

## Method (why this is exact and cheap)

cdot JSON stores the **alignment**, not the sequence: each exon is
`[g_start, g_end, exon_id, tx_start, tx_end, gap]`, a piecewise-linear map
between spliced-transcript coordinates and genomic coordinates. A version bump
happens for *sequence* reasons, but only the part of that change that alters the
**alignment** can move a coordinate — which is exactly what this measures.

`c.X` anchors at the start codon (`n_pos(c.X) = start_codon + X`), so **UTR-only
edits are automatically counted as preserving** — drift is reported only when the
*coding* alignment moves. Each version's CDS map is summarised as an exact
signature (genomic coordinate at every coding-exon breakpoint). Two maps are
unit-slope piecewise-linear, so they agree on a sub-interval iff they agree at
its left breakpoint → exact comparison and exact preserved-fraction in O(exons),
not O(bases). Gapped alignments (the linear assumption fails) are excluded from
the precise numbers and reported separately; gaps are rare (0.2% RefSeq, 0%
Ensembl).

Validated: self-vs-self = 1.0, cross-accession control = 0.0, and `c.1` of
NM_000059.4 maps to the build's stated `cds_start`.

## Results

### (1) How often is a bump coordinate-preserving?

| | RefSeq (17,162 pairs) | Ensembl (12,812 pairs) |
|---|---|---|
| identical map | 97.8% | 79.2% |
| extended-but-compatible | 0.0% | 5.7% |
| **partial drift** (some c. move, some don't) | **0.1%** | **2.9%** |
| **full drift** (whole CDS relocated) | **2.1%** | **12.1%** |
| relocated contig/strand | 0.0% | 0.0% (a literal handful) |
| **coordinate-preserving overall** | **97.8%** | **85.0%** |
| per-bump mean preserved fraction | 0.979 | 0.868 |
| **per-variant (CDS-base-weighted) safety** | **0.981** | **0.907** |

By accession class: RefSeq `NM_` (curated) 98.2% preserving, `XM_` (predicted)
97.6%; Ensembl `ENST` 85.0%.

Two things to notice:
- When a coordinate *does* move it's almost always **full drift** (the whole CDS
  shifts; `first_diff = c.1`), driven by a CDS re-annotation (`CDS length
  changed` on ~96–99% of drifting pairs). The genuinely insidious case —
  **partial** drift, where a substitution mis-places *some* variants but not
  others — is rare (0.1% RefSeq, 2.9% Ensembl).
- Per-variant safety (a random coding base) is *higher* than per-bump,
  especially for Ensembl (0.907 vs 0.868), i.e. drifting bumps are concentrated
  in shorter transcripts.

### (2b) Single-version predictor: alignment gap

A gapped cDNA→genome alignment (indels vs the reference) is the strongest signal
you can read off the *one* version you hold:

| held version | next bump moves the map |
|---|---|
| RefSeq, has gap | **63.6%** |
| RefSeq, no gap | 2.2% |

≈29× risk multiplier — but gaps are rare (57/29,206 RefSeq versions, **0** in
Ensembl, which carries no cigar gaps in this data). High precision, low coverage.

### (3) Oscillation vs one-way

Walking each accession's full version chain:

| | RefSeq | Ensembl |
|---|---|---|
| chains with an identical map across **all** versions | 97.1% | 78.6% |
| chains with **≥1 revert** (a map returning after changing) | **0.2%** | **0.2%** |
| mean map-changes per chain | 0.032 | 0.222 |

**Drift is essentially one-way / monotone.** Once coordinates move they don't
come back. This is *why* interpolation works: if a mapping never reverts, two
agreeing versions must bracket an interval where nothing moved.

### (4) Bracketing / interpolation — "if v2 and v4 agree, is v3 safe?"

Treat each interior version as the absent "requested" one; the bracket is its two
held neighbours.

| | RefSeq | Ensembl |
|---|---|---|
| base rate P(safe) | 0.976 | 0.722 |
| **P(safe \| bracket agrees)** | **0.995** | **0.961** |
| P(safe \| bracket disagrees) | 0.460 | 0.184 |
| strong test: agreeing (lo,hi) pairs whose **every** intermediate also agrees | 99.4% | 96.0% |

Bracketing is both a strong **positive** signal (agree → ~96–99.5% safe) and a
strong **negative** one (disagree → coin-flip or worse — an excellent risk flag).

### (5) Do the rules stack? (joint conditioning) — answers Q4

P(mid safe | conditions), RefSeq:

| condition | P(safe) | n |
|---|---|---|
| unconditional | 0.975 | 5,206 |
| no gap | 0.976 | 5,191 |
| has gap | 0.800 | 15 |
| bracket agree | 0.995 | 5,013 |
| bracket disagree | 0.456 | 193 |
| **no gap + bracket agree** | **0.995** | 5,002 |
| no gap + bracket disagree | 0.460 | 189 |
| has gap + bracket agree | 1.000 | 11 |
| has gap + bracket disagree | 0.250 | 4 |
| NM_ + bracket agree | 0.994 | 877 |
| XM_ + bracket agree | 0.996 | 4,136 |

Ensembl: bracket agree 0.961, disagree 0.184 (no gapped versions exist, so the
gap cells are empty).

**The headline finding: the features are *not* independent — bracketing dominates
and largely absorbs the others.** "no gap + bracket agree" ≈ "bracket agree"
(both ~0.995), because once you have an agreeing bracket the local region is
demonstrably stable. The gap signal matters mainly when you have **no** bracket
(it's a single-version prior); a gapped version's risk is largely *mediated*
through making the bracket disagree (note `has gap + bracket disagree` = 0.25).
Conditioning on an agreeing bracket, even consortium barely matters (NM 0.994,
XM 0.996, ENST 0.961). So the practical decision variable is **bracketing**, with
consortium/gap as fallbacks when you can't bracket.

## Q5 — A stored-odds model, and is it clinically acceptable?

### Feasible? Yes.

Every feature is observable at query time *from cdot's own data*: which versions
are loaded, whether they bracket the requested one and agree, the consortium, and
the gap flag. So you can precompute a calibrated table `P(safe | feature-bucket)`
(this experiment *is* that table) and return a number per query. Nothing external
is needed.

Two refinements make it much stronger than a single global odds:

1. **Prefer the exact version — usually no bump is needed.** cdot holds the
   historical versions. The fallback only fires when the *exact* requested
   version is genuinely absent. In that case you often still hold its
   *neighbours* (e.g. have v2,v4; want v3) → you can bracket.

2. **Make the check position-specific, not transcript-global.** The real
   question is "does *this variant's* `c.X` map identically across the versions I
   hold", not "is the whole bump preserving". Because drift is monotone (§3), if
   you hold flanking versions and `c.X` maps to the same genomic coordinate in
   both, that variant is safe **regardless of global drift** — even a partial-drift
   bump is fine for positions outside the moved segment. This converts a ~99.5%
   probabilistic statement into a near-deterministic per-variant check. (The
   script currently measures transcript-global preservation; a per-variant
   bracketing check is a small, recommended extension.)

### Clinically acceptable?

Short version: **a stored probability is useful for triage and transparency, but
"≥X% likely safe" should not by itself silently substitute a version in a
clinical context.** Reasons:

- HGVS is a *precise* assertion. Substituting a version changes the identifier;
  a silent ~0.5–15% chance of moving the coordinate is a wrong-variant risk, and
  the failure mode is invisible to the curator.
- The dangerous case (partial drift) is exactly the one a *global* odds smears
  over. A position-specific bracketing check addresses this directly; a global
  probability does not.
- Calibration caveats: these odds are the empirical distribution of *past* bumps
  on GRCh38 in one cdot release; future bumps, GRCh37/T2T, and per-gene
  heterogeneity may differ. Odds are a prior, not a guarantee.

A defensible clinical policy:
1. **Always try the exact requested version first** (cdot usually has it).
2. If absent, **auto-accept only on a passing position-specific bracketing
   check** (flanking held versions both map *this variant's* coordinate
   identically) — deterministic, not probabilistic.
3. Otherwise **refuse / route to human review**; never substitute silently.
   cdot already surfaces every substitution as an `HGVSFix` — keep that, and
   attach the computed odds + bracket status as evidence for the reviewer.
4. Use the stored odds for triage/ranking and transparency, **not** as the sole
   auto-commit gate.

In short: the odds are real and computable, bracketing makes them strong, and a
*position-specific* bracket check is strong enough to auto-accept; but the
appropriate clinical default remains "exact version, else human-visible warning",
with odds as supporting evidence rather than an automatic decision.

## Limitations / next steps

- GRCh38 only; one cdot release. Cross-check GRCh37 and T2T for generality.
- Measures **alignment** drift, not sequence identity (correct for the coordinate
  question, but a bump can change sequence without moving coordinates — those are
  the "identical map" cases and are safe *for coordinates*, which is what matters
  here).
- Gapped-alignment subset is tiny (RefSeq) or absent (Ensembl); the 63.6%
  gap→drift figure rests on 44 RefSeq pairs.
- The script measures transcript-global preservation; the **position-specific
  bracketing check** (the clinically interesting one) is a recommended extension.
- Single deterministic measurement (no model, no sequence fetch). If this becomes
  a paper figure, the per-variant check + GRCh37 replication are the obvious adds.

## Reproduce

```bash
python analysis/transcript_version_drift.py \
    --json /data/cdot_data/refseq/cdot-0.2.33.refseq.GRCh38.json.gz \
    --build GRCh38 --sample 12000 --seed 1
# and the ensembl/…ensembl.GRCh38.json.gz equivalent
```
