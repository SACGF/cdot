#!/usr/bin/env python3
"""
Transcript-version drift study (cdot paper, item C2 — EXPLORE ONLY).

Core question: when a transcript accession gets a new version (e.g. NM_000059.3
-> NM_000059.4), how often does a given coding (c.) position still map to the
*same* genomic coordinate?  If a version bump is coordinate-preserving,
substituting an adjacent version (cdot's opt-in `get_best_transcript_version`) is
safe; if it moves c. positions to different g. coordinates, the substitution
silently changes what the variant means.

This script answers four things:

  (1) Pairwise drift over consecutive versions: what fraction of bumps preserve
      every coding coordinate, and when they drift, by what mechanism.
  (2) Bump categorisation: identical map / extended-but-compatible / partial
      drift / full drift — separating "sequence bump that didn't touch the
      alignment" from "bump that moved coordinates".
  (3) Oscillation vs one-way: walking each accession's whole version chain, does
      a coordinate mapping ever *revert* to an earlier one, or does drift only
      ever move forward?  (If it never reverts, two agreeing versions guarantee
      everything between them agrees.)
  (4) Bracketing / interpolation: if you hold versions on *both sides* of an
      absent one and they agree, is the missing one safe?  And can a single held
      version's features (consortium, alignment gaps) predict bump safety?

Method (pure cdot JSON, no sequence, no biocommons needed):
  * cdot's merged single-build JSON keeps *all* historical versions of an
    accession as separate keys ("NM_000059.3", "NM_000059.4", ...).
  * Each transcript's genome_builds[build]["exons"] is a piecewise-linear map
    between spliced-transcript coordinates (1-based n.) and genomic coordinates:
    each exon is [g_start, g_end, exon_id, tx_start, tx_end, gap].
  * c.X anchors at the start codon (transcript-level "start_codon", a 0-based
    count of pre-CDS bases), so n_pos(c.X) = start_codon + X.  Anchoring at the
    start codon means UTR-only edits are automatically "preserving" — drift is
    reported only when the *coding* alignment moves in genomic space.
  * We summarise each version's CDS c.->g. map as an exact signature: the
    genomic coordinate at every coding-exon boundary.  Two maps are piecewise
    linear with unit slope, so they agree on a sub-interval iff they agree at its
    left breakpoint — making exact comparison and exact preserved-fraction cheap
    (O(exons), not O(bases)).

NOTE: cdot stores the genome *alignment*, not the transcript sequence.  A version
bump happens for sequence reasons, but only the part of that change that alters
the alignment can move a coordinate — which is exactly what this measures.

Standalone — NOT wired into Snakemake or the paper.

Usage:
  python analysis/transcript_version_drift.py \
      --json /data/cdot_data/refseq/cdot-0.2.33.refseq.GRCh38.json.gz \
      --build GRCh38 --sample 12000 [--seed 0]
"""
import argparse
import gzip
import json
import random
import re
from collections import Counter, defaultdict

ACC_RE = re.compile(r"^(?P<acc>[A-Z]+_?\d+|ENST\d+)\.(?P<ver>\d+)$")


def load_transcripts(path):
    with gzip.open(path, "rt") as fh:
        return json.load(fh)["transcripts"]


def cds_signature(t, build):
    """Exact signature of a version's CDS coding-coordinate map, or None.

    Returns a dict:
      contig, strand, cds_len,
      segs:  sorted list of (cds_offset_start, genomic_at_that_offset)  — one per
             coding exon segment; within a segment genomic advances by the strand
             sign, so genomic(off) = g0 + sign*(off - off0).
      bounds: tuple of the segment-start offsets (the breakpoints).
      has_gap: True if any exon carries an alignment gap (linear map inexact).
      canon: a hashable canonical form for exact-equality / revert detection.
    """
    gb = t["genome_builds"].get(build)
    if not gb:
        return None
    sc = t.get("start_codon")
    ec = t.get("stop_codon")
    if sc is None or ec is None:
        return None
    cds_len = ec - sc
    if cds_len <= 0:
        return None

    strand = gb["strand"]
    contig = gb["contig"]
    sign = 1 if strand == "+" else -1
    cds_lo = sc + 1          # first CDS tx position (1-based)
    cds_hi = sc + cds_len    # last CDS tx position

    segs = []
    has_gap = False
    for g_start, g_end, _eid, tx_start, tx_end, gap in gb["exons"]:
        if gap:
            has_gap = True
        a = max(tx_start, cds_lo)
        b = min(tx_end, cds_hi)
        if a > b:
            continue
        # genomic at the lower CDS tx position within this exon
        off_in_exon = a - tx_start
        g_at_a = (g_start + off_in_exon) if sign == 1 else (g_end - 1 - off_in_exon)
        cds_off = a - cds_lo
        segs.append((cds_off, g_at_a))

    if not segs:
        return None
    segs.sort()
    bounds = tuple(s[0] for s in segs)
    canon = (contig, strand, cds_len, tuple(segs))
    return {
        "contig": contig, "strand": strand, "cds_len": cds_len,
        "sign": sign, "segs": segs, "bounds": bounds,
        "has_gap": has_gap, "canon": canon,
    }


def genomic_at(sig, off):
    """Genomic coordinate of CDS offset `off` (0-based) under signature `sig`."""
    segs = sig["segs"]
    # last segment whose start <= off
    lo, hi = 0, len(segs) - 1
    idx = 0
    while lo <= hi:
        mid = (lo + hi) // 2
        if segs[mid][0] <= off:
            idx = mid
            lo = mid + 1
        else:
            hi = mid - 1
    off0, g0 = segs[idx]
    return g0 + sig["sign"] * (off - off0)


def preserved_fraction(a, b):
    """Exact fraction of shared CDS offsets where a and b map to the same genome.

    Returns a float in [0,1], or None if not comparable.  0.0 if the contig or
    strand differs (a wholesale relocation — nothing is preserved).
    """
    if a is None or b is None:
        return None
    if a["contig"] != b["contig"] or a["strand"] != b["strand"]:
        return 0.0
    L = min(a["cds_len"], b["cds_len"])
    if L <= 0:
        return None
    bps = sorted({0, L} | {o for o in a["bounds"] if 0 < o < L}
                            | {o for o in b["bounds"] if 0 < o < L})
    same = total = 0
    for x0, x1 in zip(bps, bps[1:]):
        length = x1 - x0
        if length <= 0:
            continue
        total += length
        if genomic_at(a, x0) == genomic_at(b, x0):
            same += length
    return same / total if total else None


def classify_bump(a, b):
    """Categorise a consecutive-version bump from its two signatures."""
    if a["contig"] != b["contig"]:
        return "relocated_contig"
    if a["strand"] != b["strand"]:
        return "relocated_strand"
    frac = preserved_fraction(a, b)
    if frac == 1.0:
        # whole shared range identical; did the CDS length change at the ends?
        return "identical_map" if a["cds_len"] == b["cds_len"] else "extended_compatible"
    if frac == 0.0:
        return "full_drift"
    return "partial_drift"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True)
    ap.add_argument("--build", default="GRCh38")
    ap.add_argument("--sample", type=int, default=12000,
                    help="number of multi-version accessions to sample (0 = all)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    print(f"Loading {args.json} ...", flush=True)
    tx = load_transcripts(args.json)
    print(f"  {len(tx):,} transcript versions", flush=True)

    # Group coding accessions (need a CDS / c. coordinates) by accession.
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

    multi = {a: v for a, v in by_acc.items() if len(v) > 1}
    print(f"  {len(by_acc):,} coding accessions; {len(multi):,} with >1 version", flush=True)

    accs = sorted(multi)
    if args.sample and args.sample < len(accs):
        accs = random.Random(args.seed).sample(accs, args.sample)
    print(f"  analysing {len(accs):,} accessions", flush=True)

    # Precompute per-version signatures.
    sigs = {}  # acc -> {ver: sig}
    for acc in accs:
        d = {}
        for ver, t in multi[acc].items():
            s = cds_signature(t, args.build)
            if s is not None:
                d[ver] = s
        if len(d) > 1:
            sigs[acc] = d

    # ---- (1)+(2) consecutive-pair drift and bump categorisation ----------
    pair_class = Counter()
    fracs = []
    drift_mech = Counter()
    by_prefix = defaultdict(lambda: Counter())
    gap_drift = Counter()       # (held_has_gap, drifted) -> n
    worst = []
    consec_pairs = 0
    for acc in sigs:
        vs = sorted(sigs[acc])
        prefix = acc[:4] if acc.startswith("ENST") else acc[:3]
        for v_lo, v_hi in zip(vs, vs[1:]):
            a, b = sigs[acc][v_lo], sigs[acc][v_hi]
            if a["has_gap"] or b["has_gap"]:
                continue  # linear map inexact for gapped alignments
            consec_pairs += 1
            cls = classify_bump(a, b)
            pair_class[cls] += 1
            frac = preserved_fraction(a, b)
            fracs.append(frac)
            drifted = frac < 1.0
            by_prefix[prefix][cls] += 1
            gap_drift[(a["has_gap"], drifted)] += 1
            if drifted:
                if a["contig"] != b["contig"]:
                    drift_mech["contig changed"] += 1
                if a["strand"] != b["strand"]:
                    drift_mech["strand changed"] += 1
                if a["cds_len"] != b["cds_len"]:
                    drift_mech["CDS length changed"] += 1
                worst.append((frac, acc, v_lo, v_hi, cls))

    # ---- single-version predictor: does an alignment gap flag instability? --
    # A gapped alignment means the transcript sequence doesn't co-linearly match
    # the genome (indels vs the reference) — a plausible a-priori risk signal you
    # can read off the *one* version you hold.  "map changed" here uses gap-aware
    # canonical equality, so it stays valid for gapped pairs (where the precise
    # preserved-fraction would be unreliable).
    gap_pred = Counter()   # (held_has_gap, map_changed) -> n
    n_versions = n_gap_versions = 0
    for acc in sigs:
        vs = sorted(sigs[acc])
        for v in vs:
            n_versions += 1
            if sigs[acc][v]["has_gap"]:
                n_gap_versions += 1
        for v_lo, v_hi in zip(vs, vs[1:]):
            a, b = sigs[acc][v_lo], sigs[acc][v_hi]
            changed = (a["canon"], a["has_gap"]) != (b["canon"], b["has_gap"])
            gap_pred[(a["has_gap"], changed)] += 1

    # ---- (3) oscillation / one-way over the whole version chain ----------
    chains_total = 0
    chains_stable = 0           # mapping identical across ALL versions
    chains_with_revert = 0      # a canonical map reappears after changing away
    transitions = []
    for acc in sigs:
        vs = sorted(sigs[acc])
        if any(sigs[acc][v]["has_gap"] for v in vs):
            continue
        canons = [sigs[acc][v]["canon"] for v in vs]
        chains_total += 1
        ntrans = sum(1 for x, y in zip(canons, canons[1:]) if x != y)
        transitions.append(ntrans)
        if ntrans == 0:
            chains_stable += 1
        # revert: a canon equals one seen before but not the immediately prior one
        seen = set()
        prev = None
        reverted = False
        for c in canons:
            if c != prev and c in seen:
                reverted = True
                break
            seen.add(c)
            prev = c
        if reverted:
            chains_with_revert += 1

    # ---- (4) bracketing / interpolation ---------------------------------
    # For every interior version v_i with both neighbours present, treat v_i as
    # the absent "requested" version.  Does the bracket (v_{i-1}, v_{i+1})
    # agreeing predict that v_i agrees with its lower neighbour (the substitute)?
    br_total = 0
    br_bracket_agree = 0
    br_mid_safe_given_agree = 0
    br_mid_safe_given_disagree = 0
    br_disagree = 0
    # Stronger: across ALL agreeing (lo,hi) pairs, are ALL intermediates consistent?
    strong_brackets = 0
    strong_brackets_consistent = 0
    for acc in sigs:
        vs = sorted(sigs[acc])
        if any(sigs[acc][v]["has_gap"] for v in vs):
            continue
        S = sigs[acc]
        for i in range(1, len(vs) - 1):
            lo, mid, hi = vs[i - 1], vs[i], vs[i + 1]
            br_total += 1
            agree = preserved_fraction(S[lo], S[hi]) == 1.0
            mid_safe = preserved_fraction(S[lo], S[mid]) == 1.0
            if agree:
                br_bracket_agree += 1
                if mid_safe:
                    br_mid_safe_given_agree += 1
            else:
                br_disagree += 1
                if mid_safe:
                    br_mid_safe_given_disagree += 1
        # strong bracket test over all pairs
        for x in range(len(vs)):
            for y in range(x + 2, len(vs)):
                if preserved_fraction(S[vs[x]], S[vs[y]]) == 1.0:
                    strong_brackets += 1
                    inner_ok = all(
                        preserved_fraction(S[vs[x]], S[vs[z]]) == 1.0
                        for z in range(x + 1, y)
                    )
                    if inner_ok:
                        strong_brackets_consistent += 1

    # ---- report ----------------------------------------------------------
    def pct(x, d):
        return f"{100*x/d:.1f}%" if d else "n/a"

    n = consec_pairs
    print("\n" + "=" * 72)
    print("TRANSCRIPT-VERSION DRIFT  —  " + args.json.split("/")[-1])
    print("=" * 72)
    print(f"build {args.build} | accessions {len(sigs):,} | "
          f"consecutive clean pairs {n:,} (gapped pairs excluded)")

    print("\n(1) Consecutive-bump categories")
    for cls in ["identical_map", "extended_compatible", "partial_drift",
                "full_drift", "relocated_contig", "relocated_strand"]:
        if pair_class.get(cls):
            print(f"    {cls:20s} {pair_class[cls]:>7,}  ({pct(pair_class[cls], n)})")
    preserving = pair_class["identical_map"] + pair_class["extended_compatible"]
    drifting = n - preserving
    print(f"    -> coordinate-preserving (any c. base safe): {preserving:,} ({pct(preserving, n)})")
    print(f"    -> drifting                                : {drifting:,} ({pct(drifting, n)})")
    if fracs:
        print(f"    mean preserved fraction {sum(fracs)/len(fracs):.4f} | "
              f"median {sorted(fracs)[len(fracs)//2]:.4f}")

    print("\n    by consortium / accession class:")
    for prefix, c in sorted(by_prefix.items()):
        tot = sum(c.values())
        pres = c["identical_map"] + c["extended_compatible"]
        print(f"      {prefix:5s} pairs={tot:>7,}  "
              f"preserving={pct(pres, tot)}  "
              f"partial={pct(c['partial_drift'], tot)}  "
              f"full={pct(c['full_drift'], tot)}  "
              f"relocated={pct(c['relocated_contig']+c['relocated_strand'], tot)}")

    print("\n    drift mechanism (among drifting pairs):")
    for k, v in drift_mech.most_common():
        print(f"      {k:22s} {v:,}")

    print("\n(2b) Single-version predictor: alignment gap")
    print(f"    versions carrying an alignment gap: {n_gap_versions:,}/{n_versions:,} "
          f"({pct(n_gap_versions, n_versions)})")
    for flag in (True, False):
        chg = gap_pred[(flag, True)]
        tot = chg + gap_pred[(flag, False)]
        print(f"    held version has_gap={str(flag):5s}: next bump changes map "
              f"{chg:,}/{tot:,} ({pct(chg, tot)})")

    print("\n(3) Oscillation vs one-way (full version chains, gap-free)")
    print(f"    chains analysed                : {chains_total:,}")
    print(f"    stable across ALL versions     : {chains_stable:,} ({pct(chains_stable, chains_total)})")
    print(f"    chains with >=1 revert (oscillation): "
          f"{chains_with_revert:,} ({pct(chains_with_revert, chains_total)})")
    if transitions:
        print(f"    mean map-changes per chain     : {sum(transitions)/len(transitions):.3f}")

    print("\n(4) Bracketing / interpolation")
    print(f"    interior versions tested       : {br_total:,}")
    if br_total:
        base = (br_mid_safe_given_agree + br_mid_safe_given_disagree) / br_total
        print(f"    base rate P(mid safe)          : {base:.4f}")
        print(f"    bracket (lo,hi) agree          : {br_bracket_agree:,} ({pct(br_bracket_agree, br_total)})")
        print(f"    P(mid safe | bracket agree)    : "
              f"{br_mid_safe_given_agree/br_bracket_agree:.4f}" if br_bracket_agree else "    n/a")
        print(f"    P(mid safe | bracket disagree) : "
              f"{br_mid_safe_given_disagree/br_disagree:.4f}" if br_disagree else "    n/a")
    print(f"    strong test: agreeing (lo,hi) pairs = {strong_brackets:,}; "
          f"all intermediates also agree = {pct(strong_brackets_consistent, strong_brackets)}")

    print("\n    most-drifted consecutive pairs (public accessions):")
    for frac, acc, v_lo, v_hi, cls in sorted(worst)[:8]:
        print(f"      {acc}.{v_lo}->{v_hi}  preserved={frac:.3f}  {cls}")
    print("=" * 72)


if __name__ == "__main__":
    main()
