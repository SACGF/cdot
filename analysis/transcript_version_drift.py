#!/usr/bin/env python3
"""
Transcript-version drift study (cdot paper, item C2 — EXPLORE ONLY).

Question: when a transcript accession gets a new version (e.g. NM_000059.3 ->
NM_000059.4), how often does a given coding (c.) position still map to the *same*
genomic coordinate?  If a version bump is coordinate-preserving, substituting an
adjacent version (cdot's opt-in `get_best_transcript_version`) is safe; if it
moves c. positions to different g. coordinates, the substitution silently changes
what the variant means.

Method (pure cdot JSON, no sequence, no biocommons needed):
  * cdot's merged single-build JSON keeps *all* historical versions of an
    accession as separate keys ("NM_000059.3", "NM_000059.4", ...).
  * Each transcript's genome_builds[build]["exons"] is a piecewise-linear map
    between spliced-transcript coordinates (1-based n.) and genomic coordinates:
    each exon is [g_start, g_end, exon_id, tx_start, tx_end, gap].
  * c.X anchors at the start codon (transcript-level "start_codon", a 0-based
    count of pre-CDS bases), so n_pos(c.X) = start_codon + X.  Anchoring at the
    start codon means UTR-only edits are automatically "preserving" — drift is
    reported only when the *coding* exon structure moves in genomic space.
  * For each consecutive (v_n, v_{n+1}) pair we map every CDS position c.1 ..
    c.(min CDS length) through both versions and count how many land on a
    different genomic coordinate.

Outputs summary statistics over the sampled accessions.  Standalone — NOT wired
into Snakemake or the paper.

Usage:
  python analysis/transcript_version_drift.py \
      --json /data/cdot_data/refseq/cdot-0.2.33.refseq.GRCh38.json.gz \
      --build GRCh38 --sample 3000 [--seed 0] [--max-cds 30000]
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


def build_segments(gb):
    """Return (segments, has_gap, strand, contig) for one genome build.

    segments: list of (tx_start, tx_end, g_start, g_end) in transcript order.
    Linear within each exon; gaps (alignment indels vs the genome) make the
    linear assumption wrong, so we flag the transcript instead of trusting it.
    """
    strand = gb["strand"]
    contig = gb["contig"]
    segs = []
    has_gap = False
    for exon in gb["exons"]:
        g_start, g_end, _exon_id, tx_start, tx_end, gap = exon
        if gap:
            has_gap = True
        segs.append((tx_start, tx_end, g_start, g_end))
    return segs, has_gap, strand, contig


def tx_to_genomic(segs, strand, p):
    """Map a 1-based spliced-transcript position to a 0-based genomic coordinate."""
    for tx_start, tx_end, g_start, g_end in segs:
        if tx_start <= p <= tx_end:
            offset = p - tx_start
            if strand == "+":
                return g_start + offset
            return g_end - 1 - offset
    return None  # position outside the exon span (shouldn't happen for valid c.)


def cds_length(t):
    sc = t.get("start_codon")
    ec = t.get("stop_codon")
    if sc is None or ec is None:
        return None
    return ec - sc  # includes the stop codon; fine for a relative comparison


def compare_pair(t_lo, t_hi, build, max_cds):
    """Compare two versions of one accession over their shared CDS c. range.

    Returns a dict of metrics, or None if not comparable (missing build, no CDS).
    """
    gb_lo = t_lo["genome_builds"].get(build)
    gb_hi = t_hi["genome_builds"].get(build)
    if not gb_lo or not gb_hi:
        return None

    len_lo = cds_length(t_lo)
    len_hi = cds_length(t_hi)
    if not len_lo or not len_hi or len_lo <= 0 or len_hi <= 0:
        return None

    segs_lo, gap_lo, strand_lo, contig_lo = build_segments(gb_lo)
    segs_hi, gap_hi, strand_hi, contig_hi = build_segments(gb_hi)

    sc_lo = t_lo["start_codon"]
    sc_hi = t_hi["start_codon"]

    overlap = min(len_lo, len_hi)
    if overlap > max_cds:
        overlap = max_cds  # cap pathological lengths for speed

    diff = 0
    compared = 0
    first_diff_c = None
    for x in range(1, overlap + 1):
        g_lo = tx_to_genomic(segs_lo, strand_lo, sc_lo + x)
        g_hi = tx_to_genomic(segs_hi, strand_hi, sc_hi + x)
        if g_lo is None or g_hi is None:
            continue
        compared += 1
        same = (g_lo == g_hi) and (contig_lo == contig_hi) and (strand_lo == strand_hi)
        if not same:
            diff += 1
            if first_diff_c is None:
                first_diff_c = x

    if compared == 0:
        return None

    return {
        "compared": compared,
        "diff": diff,
        "preserved_frac": (compared - diff) / compared,
        "cds_len_changed": len_lo != len_hi,
        "contig_changed": contig_lo != contig_hi,
        "strand_changed": strand_lo != strand_hi,
        "start_codon_changed": sc_lo != sc_hi,
        "has_gap": gap_lo or gap_hi,
        "first_diff_c": first_diff_c,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True)
    ap.add_argument("--build", default="GRCh38")
    ap.add_argument("--sample", type=int, default=3000,
                    help="number of multi-version accessions to sample (0 = all)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-cds", type=int, default=30000)
    args = ap.parse_args()

    print(f"Loading {args.json} ...", flush=True)
    tx = load_transcripts(args.json)
    print(f"  {len(tx):,} transcript versions", flush=True)

    # Group versions by accession (only RefSeq/Ensembl mRNA-style coding accessions
    # with a CDS are interesting; non-coding (NR_/XR_) have no c. coordinates).
    by_acc = defaultdict(dict)
    for key, t in tx.items():
        m = ACC_RE.match(key)
        if not m:
            continue
        acc, ver = m.group("acc"), int(m.group("ver"))
        if not acc.startswith(("NM_", "XM_", "ENST")):
            continue  # need a CDS / c. coordinates
        if t.get("start_codon") is None or t.get("stop_codon") is None:
            continue
        by_acc[acc][ver] = t

    multi = {a: v for a, v in by_acc.items() if len(v) > 1}
    print(f"  {len(by_acc):,} coding accessions; {len(multi):,} with >1 version", flush=True)

    accs = sorted(multi)
    if args.sample and args.sample < len(accs):
        rng = random.Random(args.seed)
        accs = rng.sample(accs, args.sample)
    print(f"  analysing {len(accs):,} accessions", flush=True)

    results = []
    skipped_gap = 0
    for acc in accs:
        versions = sorted(multi[acc])
        for v_lo, v_hi in zip(versions, versions[1:]):
            r = compare_pair(multi[acc][v_lo], multi[acc][v_hi], args.build, args.max_cds)
            if r is None:
                continue
            r["acc"] = acc
            r["v_lo"] = v_lo
            r["v_hi"] = v_hi
            results.append(r)

    n = len(results)
    if n == 0:
        print("No comparable pairs found.")
        return

    # ---- summary ----------------------------------------------------------
    with_gap = [r for r in results if r["has_gap"]]
    clean = [r for r in results if not r["has_gap"]]  # linear map is exact

    fully_preserved = [r for r in clean if r["diff"] == 0]
    fully_diverged = [r for r in clean if r["preserved_frac"] == 0.0]
    partial = [r for r in clean if 0.0 < r["preserved_frac"] < 1.0]

    fracs = sorted(r["preserved_frac"] for r in clean)

    def pct(x, d):
        return f"{100*x/d:.1f}%" if d else "n/a"

    def quantile(sorted_vals, q):
        if not sorted_vals:
            return float("nan")
        i = min(len(sorted_vals) - 1, int(q * (len(sorted_vals) - 1)))
        return sorted_vals[i]

    print("\n" + "=" * 70)
    print("TRANSCRIPT-VERSION DRIFT — consecutive-version pairs")
    print("=" * 70)
    print(f"build                       : {args.build}")
    print(f"accessions analysed         : {len(accs):,}")
    print(f"consecutive pairs compared  : {n:,}")
    print(f"  pairs with alignment gaps : {len(with_gap):,} "
          f"(linear map inexact — excluded from %s below)")
    print(f"  clean pairs (exact map)   : {len(clean):,}")
    print()
    print("Over CLEAN pairs (genomic c.->g. comparison is exact):")
    print(f"  fully coordinate-preserving (every c. base same g.) : "
          f"{len(fully_preserved):,} ({pct(len(fully_preserved), len(clean))})")
    print(f"  partially preserving (some c. bases drift)          : "
          f"{len(partial):,} ({pct(len(partial), len(clean))})")
    print(f"  fully diverged (no c. base shares a g. coord)       : "
          f"{len(fully_diverged):,} ({pct(len(fully_diverged), len(clean))})")
    print()
    mean_frac = sum(fracs) / len(fracs)
    print(f"  mean   preserved fraction : {mean_frac:.4f}")
    print(f"  median preserved fraction : {quantile(fracs, 0.5):.4f}")
    print(f"  p05 / p25 preserved frac  : {quantile(fracs, 0.05):.4f} / {quantile(fracs, 0.25):.4f}")
    print()
    # Mechanism breakdown over the pairs that drift at all.
    drifting = [r for r in clean if r["diff"] > 0]
    print(f"Among {len(drifting):,} drifting clean pairs:")
    print(f"  contig changed            : {sum(r['contig_changed'] for r in drifting):,}")
    print(f"  strand changed            : {sum(r['strand_changed'] for r in drifting):,}")
    print(f"  start_codon (CDS) shifted : {sum(r['start_codon_changed'] for r in drifting):,}")
    print(f"  CDS length changed        : {sum(r['cds_len_changed'] for r in drifting):,}")
    print()
    # For the *non*-drifting pairs, did anything change at all?
    sc_shift_preserved = sum(r["start_codon_changed"] for r in fully_preserved)
    print(f"Among {len(fully_preserved):,} fully-preserving clean pairs:")
    print(f"  start_codon shifted but g. preserved (UTR-only edit) : {sc_shift_preserved:,}")
    print()
    # Curated (NM_) vs predicted (XM_) split — the safety question that matters
    # clinically is whether *curated* RefSeq bumps preserve coordinates.
    print("By accession class (clean pairs):")
    for prefix in ("NM_", "XM_", "ENST"):
        sub = [r for r in clean if r["acc"].startswith(prefix)]
        if not sub:
            continue
        pres = sum(1 for r in sub if r["diff"] == 0)
        div = sum(1 for r in sub if r["preserved_frac"] == 0.0)
        print(f"  {prefix:5s} pairs={len(sub):>6,}  "
              f"fully-preserving={pct(pres, len(sub))}  "
              f"fully-diverged={pct(div, len(sub))}")
    print()
    # Worst offenders for context (no private data — these are public accessions).
    worst = sorted(drifting, key=lambda r: r["preserved_frac"])[:10]
    print("Most-drifted example pairs (public accessions):")
    for r in worst:
        print(f"  {r['acc']}.{r['v_lo']}->{r['v_hi']}  "
              f"preserved={r['preserved_frac']:.3f}  "
              f"diff={r['diff']}/{r['compared']}  "
              f"first_diff=c.{r['first_diff_c']}  "
              f"{'CDSlen ' if r['cds_len_changed'] else ''}"
              f"{'contig ' if r['contig_changed'] else ''}")
    print("=" * 70)


if __name__ == "__main__":
    main()
