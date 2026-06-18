#!/usr/bin/env python3
"""
Cross-build transcript-version analysis (cdot paper, item C2 — EXPLORE ONLY).

Follow-up to transcript_version_drift.py.  Single-build bracketing has two blind
spots: (a) the requested version may be absent from the target build, and (b)
transient-revert versions (A->B->A) are invisible to a bracket that doesn't hold
the middle version.  Other genome builds can help with both, because:

  * cdot's data (and the REST API) carries each accession's versions across
    GRCh37 / GRCh38 / T2T, and the version *sets differ by build* — a version
    missing from GRCh38 is often present in GRCh37.
  * a version's *intrinsic* structure (CDS length + the coding-exon lengths in
    transcript coordinates) is build-INDEPENDENT.  Only the genomic placement
    changes between builds.  So another build's copy of a version tells you that
    version's c. structure even when you can't place it in the target build.

This script uses the single all-builds file (each version record has a
genome_builds dict spanning builds) and measures:

  X1  availability: how many versions does pooling builds add, and how often is a
      version absent from GRCh38 present in another build (so REST would surface
      it)?
  X2  same-version cross-build intrinsic-structure identity (is structure really
      build-portable?).
  X3  cross-build concordance of the drift call for a consecutive bump.
  X4  transient detection: are drifting / transient versions flagged by their
      build-independent intrinsic structure alone (CDS-length / exon-length
      change), i.e. can REST's all-build view catch what single-build bracketing
      misses?

Standalone — NOT wired into Snakemake, the paper, or the library.

Usage:
  python analysis/transcript_version_crossbuild.py \
      --json /data/cdot_data/refseq/cdot-0.2.33.all-builds-refseq-grch37_grch38_t2t-chm13v2.0.json.gz \
      --target GRCh38 --other GRCh37 --sample 12000 [--seed 0]
"""
import argparse
import gzip
import json
import random
import re
from collections import Counter, defaultdict

from transcript_version_drift import cds_signature, preserved_fraction

ACC_RE = re.compile(r"^(?P<acc>[A-Z]+_?\d+|ENST\d+)\.(?P<ver>\d+)$")


def load_transcripts(path):
    with gzip.open(path, "rt") as fh:
        return json.load(fh)["transcripts"]


def intrinsic(sig):
    """Build-independent CDS shape: strand, CDS length, coding-exon-segment lengths.

    Genomic placement (intron sizes / locus) is excluded — only what the spliced
    coding sequence looks like, which is identical for a version in any build.
    """
    bounds = list(sig["bounds"]) + [sig["cds_len"]]
    seglens = tuple(b - a for a, b in zip(bounds, bounds[1:]))
    return (sig["strand"], sig["cds_len"], seglens)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True)
    ap.add_argument("--target", default="GRCh38")
    ap.add_argument("--other", default="GRCh37")
    ap.add_argument("--sample", type=int, default=12000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    print(f"Loading {args.json} ...", flush=True)
    tx = load_transcripts(args.json)
    print(f"  {len(tx):,} version records", flush=True)

    # accession -> version -> transcript record (coding only)
    by_acc = defaultdict(dict)
    for key, t in tx.items():
        m = ACC_RE.match(key)
        if not m:
            continue
        acc, ver = m.group("acc"), int(m.group("ver"))
        if not acc.startswith(("NM_", "XM_", "ENST")):
            continue
        if t.get("start_codon") is None or t.get("stop_codon") is None:
            continue
        by_acc[acc][ver] = t

    accs = sorted(a for a, v in by_acc.items() if len(v) > 1)
    if args.sample and args.sample < len(accs):
        accs = random.Random(args.seed).sample(accs, args.sample)
    print(f"  analysing {len(accs):,} multi-version accessions", flush=True)

    T, O = args.target, args.other

    # ---- X1: availability across builds ---------------------------------
    av = Counter()                  # which-builds-set -> n (acc,ver)
    holes_total = holes_filled = 0  # internal target-build gaps; filled by other?
    add_versions = 0                # versions absent in target but present elsewhere
    n_acc_with_target = 0
    for acc in accs:
        vmap = by_acc[acc]
        vt = sorted(v for v, t in vmap.items() if T in t["genome_builds"])
        vo = {v for v, t in vmap.items() if O in t["genome_builds"]}
        vall = set(vmap)
        for v, t in vmap.items():
            bs = tuple(b for b in (O, T) if b in t["genome_builds"])
            av[bs] += 1
        if vt:
            n_acc_with_target += 1
            add_versions += len(vall - set(vt))
            lo, hi = vt[0], vt[-1]
            for v in range(lo + 1, hi):
                if v not in vt:                       # internal gap in target build
                    holes_total += 1
                    if v in vo:                       # filled by the other build
                        holes_filled += 1

    # ---- X2/X3/X4: signatures per (acc, ver, build) ---------------------
    # Precompute target+other signatures.
    sig_t = {}  # (acc, ver) -> sig in target
    sig_o = {}  # (acc, ver) -> sig in other
    for acc in accs:
        for v, t in by_acc[acc].items():
            if T in t["genome_builds"]:
                s = cds_signature(t, T)
                if s and not s["has_gap"]:
                    sig_t[(acc, v)] = s
            if O in t["genome_builds"]:
                s = cds_signature(t, O)
                if s and not s["has_gap"]:
                    sig_o[(acc, v)] = s

    # X2: same-version intrinsic identity across builds
    x2_same = x2_total = 0
    for key in sig_t:
        if key in sig_o:
            x2_total += 1
            if intrinsic(sig_t[key]) == intrinsic(sig_o[key]):
                x2_same += 1

    # X3: cross-build concordance of the consecutive-bump drift call
    x3_concordant = x3_total = x3_both_preserve = x3_both_drift = x3_discordant = 0
    # X4: structural detectability of drift (intrinsic change), per build (target)
    x4_drift = x4_drift_struct = x4_drift_cdslen = 0
    # converse: among bumps whose intrinsic structure is UNCHANGED, are they
    # genomically preserved? (is structure-unchanged sufficient for safety?)
    x4_same_struct = x4_same_struct_preserved = 0
    for acc in accs:
        vs = sorted(by_acc[acc])
        for lo, hi in zip(vs, vs[1:]):
            kt_lo, kt_hi = (acc, lo), (acc, hi)
            # X4 (target build): is a drifting bump flagged by intrinsic change?
            if kt_lo in sig_t and kt_hi in sig_t:
                a, b = sig_t[kt_lo], sig_t[kt_hi]
                preserved = preserved_fraction(a, b) == 1.0
                if not preserved:
                    x4_drift += 1
                    if intrinsic(a) != intrinsic(b):
                        x4_drift_struct += 1
                    if a["cds_len"] != b["cds_len"]:
                        x4_drift_cdslen += 1
                if intrinsic(a) == intrinsic(b):
                    x4_same_struct += 1
                    if preserved:
                        x4_same_struct_preserved += 1
            # X3: both versions present (gap-free) in BOTH builds
            ko_lo, ko_hi = (acc, lo), (acc, hi)
            if (kt_lo in sig_t and kt_hi in sig_t
                    and ko_lo in sig_o and ko_hi in sig_o):
                x3_total += 1
                pt = preserved_fraction(sig_t[kt_lo], sig_t[kt_hi]) == 1.0
                po = preserved_fraction(sig_o[ko_lo], sig_o[ko_hi]) == 1.0
                if pt == po:
                    x3_concordant += 1
                    if pt:
                        x3_both_preserve += 1
                    else:
                        x3_both_drift += 1
                else:
                    x3_discordant += 1

    # ---- report ----------------------------------------------------------
    def pct(x, d):
        return f"{100*x/d:.1f}%" if d else "n/a"

    print("\n" + "=" * 72)
    print(f"CROSS-BUILD VERSION ANALYSIS  ({O} <-> {T})")
    print("=" * 72)

    print("\n(X1) Availability across builds")
    for bs, n in sorted(av.items(), key=lambda kv: -kv[1]):
        label = "+".join(bs) if bs else "(neither)"
        print(f"    versions in {label:18s}: {n:,}")
    print(f"    accessions with >=1 {T} version : {n_acc_with_target:,}")
    print(f"    extra versions pooling builds adds (absent in {T}): {add_versions:,}")
    print(f"    internal version gaps in {T}     : {holes_total:,}")
    print(f"      of those present in {O} (REST could supply): "
          f"{holes_filled:,} ({pct(holes_filled, holes_total)})")

    print("\n(X2) Same version, both builds: intrinsic structure identical?")
    print(f"    versions present (gap-free) in both: {x2_total:,}")
    print(f"    intrinsic CDS structure identical  : {x2_same:,} ({pct(x2_same, x2_total)})")

    print("\n(X3) Cross-build concordance of the bump drift call")
    print(f"    bumps with both versions in both builds: {x3_total:,}")
    print(f"    concordant call                        : {x3_concordant:,} ({pct(x3_concordant, x3_total)})")
    print(f"      both preserve {x3_both_preserve:,} | both drift {x3_both_drift:,}")
    print(f"    DISCORDANT (build-specific realignment) : {x3_discordant:,} ({pct(x3_discordant, x3_total)})")

    print("\n(X4) Are drifting versions flagged by build-independent structure?")
    print(f"    drifting bumps in {T}              : {x4_drift:,}")
    print(f"    intrinsic structure changed (detectable in ANY build): "
          f"{x4_drift_struct:,} ({pct(x4_drift_struct, x4_drift)})")
    print(f"    CDS length changed                 : "
          f"{x4_drift_cdslen:,} ({pct(x4_drift_cdslen, x4_drift)})")
    print(f"    converse - bumps with UNCHANGED intrinsic structure: {x4_same_struct:,}")
    print(f"      of those, genomically preserved  : "
          f"{x4_same_struct_preserved:,} ({pct(x4_same_struct_preserved, x4_same_struct)})")
    print("=" * 72)


if __name__ == "__main__":
    main()
