#!/usr/bin/env python3
"""Compute transcript-version stability facts from cdot release JSON.gz files.

Populates ``output/facts/version_stability.csv`` for the paper's R5 section
("Transcript-version stability and safe version fallback"). This is the
reproducible (Tier-1) write-up of the C2 experiment: how often a transcript
version bump is coordinate-preserving, and the pivotal result that a build-
independent intrinsic-CDS-structure change is (almost exactly) equivalent to a
genomic-coordinate drift — which is what makes cdot's safe-version-fallback check
near-deterministic (see cdot.hgvs.version_safety / is_version_substitution_safe).

The drift / concentration / per-variant numbers come from the single-build GRCh38
files (which keep every historical version as a separate key); the cross-build
structure-portability and structure⇔drift-equivalence numbers come from the
all-builds files (each version record spans GRCh37/GRCh38/T2T). All measurement
reuses the packaged cdot.hgvs.version_safety helpers, so the facts and the shipped
safety check are computed by the same code.

Per CLAUDE.md, run this only against production release files on a dedicated run,
never the GTF/GFF generation pipeline. Sampling (--sample) keeps it cheap and
deterministic (--seed).

Usage:
    python paper/scripts/compute_version_stability.py \
        --refseq-grch38   cdot-0.2.33.refseq.GRCh38.json.gz \
        --ensembl-grch38  cdot-0.2.33.ensembl.GRCh38.json.gz \
        --refseq-allbuilds  cdot-0.2.33.all-builds-refseq-...json.gz \
        --ensembl-allbuilds cdot-0.2.33.all-builds-ensembl-...json.gz \
        --build GRCh38 --other GRCh37 --sample 12000 --seed 1
"""

import argparse
import gzip
import json
import random
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from cdot.hgvs.version_safety import (
    cds_genomic_map,
    intrinsic_cds_structure,
    _genomic_at,
)

ACC_RE = re.compile(r"^(?P<acc>[A-Z]+_?\d+|ENST\d+)\.(?P<ver>\d+)$")


def _load(path):
    with gzip.open(path, "rt") as fh:
        return json.load(fh)["transcripts"]


def _by_acc(transcripts):
    """Group coding (NM_/XM_/ENST) version records by versionless accession."""
    out = defaultdict(dict)
    for key, t in transcripts.items():
        m = ACC_RE.match(key)
        if not m:
            continue
        acc, ver = m.group("acc"), int(m.group("ver"))
        if not acc.startswith(("NM_", "XM_", "ENST")):
            continue
        if t.get("start_codon") is None or t.get("stop_codon") is None:
            continue
        out[acc][ver] = t
    return out


def _sample(accs, n, seed):
    accs = sorted(accs)
    if n and n < len(accs):
        accs = random.Random(seed).sample(accs, n)
    return accs


def _preserved_fraction(a, b):
    """Exact fraction of shared CDS offsets where two genomic maps agree (or None).

    Uses the packaged cds_genomic_map signatures: unit-slope linear, so equality is
    constant between breakpoints. 0.0 if contig/strand differ (whole-CDS relocation).
    """
    if a is None or b is None or a["has_gap"] or b["has_gap"]:
        return None
    if a["contig"] != b["contig"] or a["strand"] != b["strand"]:
        return 0.0
    L = min(a["cds_len"], b["cds_len"])
    if L <= 0:
        return None
    bps = sorted({0, L}
                 | {o for o, _g in a["segs"] if 0 < o < L}
                 | {o for o, _g in b["segs"] if 0 < o < L})
    same = total = 0
    for x0, x1 in zip(bps, bps[1:]):
        length = x1 - x0
        if length <= 0:
            continue
        total += length
        if _genomic_at(a, x0) == _genomic_at(b, x0):
            same += length
    return same / total if total else None


def drift_stats(transcripts, build, sample, seed):
    """Single-build drift, per-variant safety and concentration for one consortium."""
    by_acc = _by_acc(transcripts)
    multi = {a: v for a, v in by_acc.items() if len(v) > 1}
    accs = _sample(multi, sample, seed)

    pairs = 0
    preserving = full_drift = partial_drift = 0
    bases_pres = bases_total = 0
    acc_drifts = Counter()
    acc_seen = set()
    for acc in accs:
        vs = sorted(multi[acc])
        maps = {v: cds_genomic_map(multi[acc][v], build) for v in vs}
        maps = {v: m for v, m in maps.items() if m is not None}
        vs = sorted(maps)
        for lo, hi in zip(vs, vs[1:]):
            a, b = maps[lo], maps[hi]
            if a["has_gap"] or b["has_gap"]:
                continue
            frac = _preserved_fraction(a, b)
            if frac is None:
                continue
            pairs += 1
            acc_seen.add(acc)
            overlap = min(a["cds_len"], b["cds_len"])
            bases_total += overlap
            bases_pres += round(frac * overlap)
            if frac == 1.0:
                preserving += 1
            elif frac == 0.0:
                full_drift += 1
                acc_drifts[acc] += 1
            else:
                partial_drift += 1
                acc_drifts[acc] += 1

    n_acc = len(acc_seen)
    n_drift_acc = len(acc_drifts)
    total_drift = sum(acc_drifts.values())
    ranked = sorted(acc_drifts.values(), reverse=True)
    top5 = max(1, int(0.05 * n_acc))
    top5_share = (sum(ranked[:top5]) / total_drift * 100) if total_drift else 0.0

    def pct(x):
        return round(100 * x / pairs, 1) if pairs else 0.0

    return {
        "pairs": pairs,
        "preserving_pct": pct(preserving),
        "full_drift_pct": pct(full_drift),
        "partial_drift_pct": pct(partial_drift),
        "pervariant_safety": round(bases_pres / bases_total, 3) if bases_total else 0.0,
        "accessions_drift_pct": round(100 * n_drift_acc / n_acc, 1) if n_acc else 0.0,
        "top5_drift_share_pct": round(top5_share, 1),
    }


def crossbuild_stats(transcripts, build, other, sample, seed):
    """Structure portability (X2) and structure⇔drift equivalence (X4)."""
    by_acc = _by_acc(transcripts)
    accs = _sample([a for a, v in by_acc.items() if len(v) > 1], sample, seed)

    x2_same = x2_total = 0
    drift = drift_flagged = 0
    same_struct = same_struct_preserved = 0
    for acc in accs:
        vmap = by_acc[acc]
        for v, t in vmap.items():
            if build in t["genome_builds"] and other in t["genome_builds"]:
                sb = intrinsic_cds_structure(t, build)
                so = intrinsic_cds_structure(t, other)
                if sb is not None and so is not None:
                    x2_total += 1
                    if sb == so:
                        x2_same += 1
        vs = sorted(vmap)
        for lo, hi in zip(vs, vs[1:]):
            a = cds_genomic_map(vmap[lo], build)
            b = cds_genomic_map(vmap[hi], build)
            if a is None or b is None or a["has_gap"] or b["has_gap"]:
                continue
            frac = _preserved_fraction(a, b)
            if frac is None:
                continue
            sa = intrinsic_cds_structure(vmap[lo], build)
            sb = intrinsic_cds_structure(vmap[hi], build)
            if frac < 1.0:
                drift += 1
                if sa != sb:
                    drift_flagged += 1
            if sa == sb:
                same_struct += 1
                if frac == 1.0:
                    same_struct_preserved += 1

    return {
        "struct_portable_pct": round(100 * x2_same / x2_total, 1) if x2_total else 0.0,
        "drift_struct_flagged_pct": round(100 * drift_flagged / drift, 1) if drift else 0.0,
        "struct_unchanged_preserved_pct": (
            round(100 * same_struct_preserved / same_struct, 1) if same_struct else 0.0),
    }


def main():
    ap = argparse.ArgumentParser(description="cdot transcript-version stability facts.")
    ap.add_argument("--refseq-grch38", default="")
    ap.add_argument("--ensembl-grch38", default="")
    ap.add_argument("--refseq-allbuilds", default="")
    ap.add_argument("--ensembl-allbuilds", default="")
    ap.add_argument("--build", default="GRCh38")
    ap.add_argument("--other", default="GRCh37")
    ap.add_argument("--sample", type=int, default=12000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    facts = {"sample_n": args.sample, "build": args.build}

    for label, path in (("refseq", args.refseq_grch38), ("ensembl", args.ensembl_grch38)):
        if path and Path(path).exists():
            print(f"drift: loading {path} ...", flush=True)
            d = drift_stats(_load(path), args.build, args.sample, args.seed)
            for k, v in d.items():
                facts[f"{label}_{k}"] = v

    for label, path in (("refseq", args.refseq_allbuilds), ("ensembl", args.ensembl_allbuilds)):
        if path and Path(path).exists():
            print(f"crossbuild: loading {path} ...", flush=True)
            d = crossbuild_stats(_load(path), args.build, args.other, args.sample, args.seed)
            for k, v in d.items():
                facts[f"{label}_{k}"] = v

    out = Path("output/facts/version_stability.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([facts]).to_csv(out, index=False)
    print(f"Written: {out}")
    for k, v in facts.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
