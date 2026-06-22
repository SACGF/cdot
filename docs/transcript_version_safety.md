# Transcript-version safety: the safe version fallback

When an HGVS string cites a transcript version you don't have loaded (e.g. a report says
`NM_000059.2:c.36del` but your data only has `NM_000059.4`), cdot can fall back to the
nearest available version. That is only safe if the substitution is *coordinate-preserving*
— i.e. the coding `c.` position still maps to the **same** genomic coordinate under the
substitute. If it doesn't, the substitution silently changes what the variant means.

This document explains how cdot decides whether a version substitution is safe, the public
API, and the empirical study that the check is built on.

- API entry points: [`fix_hgvs` / `resolve_transcript_version`](advanced_usage.md) and the
  data-provider method `is_version_substitution_safe`.
- Reproducible analysis: `analysis/transcript_version_drift.py`,
  `analysis/transcript_version_crossbuild.py`, and the paper facts script
  `paper/scripts/compute_version_stability.py` (paper Results R5).

## The check, in one sentence

> A version bump is coordinate-preserving **iff** the version's *intrinsic CDS structure* —
> its CDS length plus the lengths of its coding-exon segments in transcript coordinates — is
> unchanged.

That structure is **build-independent** (only a version's *genomic placement* changes
between GRCh37/GRCh38/T2T), so it can be read for the requested version from *any* build
that holds it, even when that version is absent from the build you're mapping in. Comparing
the requested version's structure to the substitute's is therefore near-deterministic, needs
no flanking pair, and — crucially — catches the transient-revert versions (a coordinate that
goes A→B→A) that a genomic-neighbour check is blind to.

## Using it

The fallback is **opt-in** and never silent. Through the one-call entry point:

```python
from cdot.hgvs import fix_hgvs, VersionStrategy, UnsafeVersionPolicy

result, fixes = fix_hgvs(
    "NM_000059.2:c.36del",
    data_provider=dp,            # JSONDataProvider or RESTDataProvider
    genome_build="GRCh38",
    version_fallback=VersionStrategy.UP_THEN_DOWN,
)
```

What happens when `.2` is absent:

| Situation | `result` | `HGVSFix` code | Severity |
|---|---|---|---|
| Substitute is structure-verified safe | rewritten to the substitute | `USED_ADJACENT_VERSION_COORD_SAFE` | WARNING |
| Not safe, **`on_unsafe_version=REFUSE`** (default) | **unchanged** | `REFUSED_UNSAFE_VERSION` | ERROR |
| Not safe, `on_unsafe_version=SUBSTITUTE` | rewritten (coordinate may have moved) | `USED_ADJACENT_VERSION_COORD_UNVERIFIED` | WARNING |
| Provider can't assess safety | rewritten | `USED_ADJACENT_VERSION` | WARNING |

The default refuses a substitution that isn't verified coordinate-safe, preserving
exact-variant semantics; pass `on_unsafe_version=UnsafeVersionPolicy.SUBSTITUTE` to apply it
anyway with a "coordinate may have moved" warning. With `raise_on_errors=True`, a refusal is
raised as `HGVSInputError`.

### Lower-level pieces

```python
# Build-independent structure of a transcript record (dict or models.Transcript).
from cdot.hgvs import intrinsic_cds_structure
intrinsic_cds_structure(record)                 # -> ("+", 10257, (228, 4932, ...)) or None

# Provider-level verdict for a specific substitution.
safe, reason = dp.is_version_substitution_safe("NM_000059", 2, 4, genome_build="GRCh38")
```

`is_version_substitution_safe` reads the requested version's structure from any build that
holds it (cdot all-builds JSON files, or the REST `/transcript/<ac>` view, return every
version across builds). Only when the requested version exists in **no** build does it fall
back to a probabilistic genomic-bracket check (do the held versions either side of the gap
agree?). `get_tx_versions(accession, genome_build=...)` restricts the substitution candidate
pool to versions placeable in the target build.

## The evidence

Measured over a 12,000-accession sample on GRCh38 (cdot 0.2.33). cdot stores the genome
*alignment*, not the sequence, so a version bump can only move a coordinate if it changes
that alignment — which is exactly what is measured. Anchoring at the start codon means
UTR-only edits count as preserving; drift is reported only when the *coding* alignment moves.

**How often a bump is coordinate-preserving**

| | RefSeq | Ensembl |
|---|---|---|
| coordinate-preserving (every coding base) | **97.8%** | **85.0%** |
| per-variant safety (a random coding base) | 0.981 | 0.907 |
| partial drift (some variants move, some don't) | 0.1% | 2.9% |
| whole-CDS drift | 2.1% | 12.1% |

By class: RefSeq `NM_` (curated) 98.2% preserving, `XM_` (predicted) 97.6%. When a
coordinate does move it is almost always the **whole CDS** (a re-annotation of the coding
region), never a contig/strand change; the genuinely insidious *partial* drift is rare.
Drift is also highly concentrated: only 2.9% of RefSeq accessions ever drift, and the
worst-affected 5% hold 100% of it.

**The pivotal result — structure ⇔ drift equivalence**

A version's intrinsic CDS structure is identical across builds for **99.5%** of RefSeq /
99.8% of Ensembl versions, and a change in it is almost exactly equivalent to a genomic
drift:

| | RefSeq | Ensembl |
|---|---|---|
| drifting bumps flagged by an intrinsic-structure change | **100%** | 99.3% |
| structure-unchanged bumps that are genomically preserved | **100%** | 99.9% |

Because the structure is build-portable and cdot returns every version in every build, the
requested version's structure is readable even when it's absent from the target build — so
the safety question becomes essentially deterministic.

## Why structural, and not the alternatives

Several other signals were measured and rejected as the *primary* check:

- **Genomic bracketing** ("hold v2 and v4; if they agree, v3 is safe"): a strong signal
  (P(safe | bracket agrees) ≈ 0.995 RefSeq), but it needs a flanking pair and is capped by
  **transient-revert** versions (A→B→A): the two stable flanks agree with each other while
  the requested middle version differs, and a bracket can't see a version it doesn't hold.
  The structural check reads that middle version's changed structure from another build and
  catches it.
- **Single-version alignment-gap flag**: a gapped cDNA→genome alignment is ~29× more likely
  to drift on the next bump — high precision, but gaps are rare (57/29,206 RefSeq versions,
  none in Ensembl), so coverage is tiny.
- **Known-unstable blocklist**: because drift is so concentrated (≤5% of RefSeq accessions
  hold all of it), a precomputed blocklist is a cheap defence-in-depth — but it's a
  denylist, not a positive safety proof.
- **Stored-odds model** (`P(safe | features)`): useful for triage and transparency, but a
  *global* probability smears over the dangerous partial-drift case and is never
  deterministic.

The structural check subsumes the useful part of all of these: it is build-independent,
needs no flanking pair, is near-deterministic, and catches the transient reverts bracketing
misses. The genomic bracket survives only as the fallback for the residual case where the
requested version is in no build at all.

## Limitations

- Measured on GRCh38 from a single cdot release; the concentration/percentages are an
  empirical prior, not a guarantee. GRCh37/T2T replication is future work.
- Measures **alignment** drift, not sequence identity — correct for the coordinate question
  (a bump can change sequence without moving a coordinate; those are the safe "identical map"
  cases), but it is not a statement about sequence equality.
- A bump that only *extends* the CDS (same shared region, longer) counts as a structure
  change and is conservatively refused even though a variant in the shared region is fine —
  an acceptable false-negative for a safety gate.
- Structure is build-portable ~99.5%, not 100%: rare build-specific realignments exist.
- Identical CDS structure is necessary but not sufficient. An alignment gap (a
  transcript-vs-genome indel) present in one version and not the other shifts every coding
  base downstream of it, so the check also requires the two versions' CDS alignment gaps to
  match (`cds_alignment_gaps`). Gaps are build-specific, so this is exact when both versions
  are placed in the target build and conservative across builds.
- The check is position-aware. Coding (`c.`) and CDS-intronic positions are covered by the
  CDS structure plus alignment gaps. A 5'UTR (`c.-N`) or 3'UTR (`c.*N`) position additionally
  requires the matching UTR *length* to be unchanged (`utr_lengths`), since UTR annotations
  change between versions far more often than the CDS. A residual edge case remains: two
  versions with the same total UTR length but a different UTR exon *split* could still move a
  UTR-intronic position; covering that needs full non-coding exon structure, not just length.
- The intrinsic structure is build-independent and so is blind to genomic *placement*: a rare
  version re-aligns the same CDS structure to a different locus (eg NM_018263, a ~9.75 kb
  jump). When both versions are loaded this is caught by comparing their genomic CDS maps; for
  the case where the requested version is absent from all loaded data (nothing to compare) a
  small precomputed blocklist of known re-placements is shipped
  (`generate_transcript_data/generate_version_replacement_blocklist.py`). Re-placement is a
  genomic property and so build-specific (a pair can re-place in one build but not the other),
  and the blocklist is keyed without a build and applied in every build, so it is generated as
  the union of re-placements found across GRCh37 and GRCh38. On the 0.2.33 RefSeq release this
  is 74 substitutions across 22 transcripts.

## Reproduce

```bash
# Single-build drift, concentration, per-variant safety (Q1-Q7 in the analysis script)
python analysis/transcript_version_drift.py \
    --json /path/to/cdot-0.2.33.refseq.GRCh38.json.gz --build GRCh38 --sample 12000 --seed 1

# Cross-build structure portability + structure<=>drift equivalence
python analysis/transcript_version_crossbuild.py \
    --json /path/to/cdot-0.2.33.all-builds-refseq-...json.gz \
    --target GRCh38 --other GRCh37 --sample 12000 --seed 1

# The paper's one-row facts CSV (reuses cdot.hgvs.version_safety)
python paper/scripts/compute_version_stability.py \
    --refseq-grch38 ... --ensembl-grch38 ... \
    --refseq-allbuilds ... --ensembl-allbuilds ... --sample 12000 --seed 1
```
