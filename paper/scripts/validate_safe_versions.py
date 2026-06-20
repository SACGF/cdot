#!/usr/bin/env python3
"""Empirical ClinVar validation of safe transcript-version substitution (R5b).

compute_version_stability.py (R5) proves *structurally* that a version bump is
coordinate-preserving iff the version's intrinsic CDS structure is unchanged. This
script is the empirical complement: take every real ClinVar variant where another
version of its transcript exists and the safety check (is_version_substitution_safe)
calls the V->W substitution safe, resolve the variant under W, and confirm the
genomic coordinate is byte-identical to resolving it under V. The structural Tier-1
claim then becomes an observed zero-error rate, the strongest reviewer defence of
the safe-version-fallback feature.

Consumes the per-variant table from resolve_clinvar_pass.py (paper infra A1): its
`converted_g` IS the under-V resolution (so V is not re-resolved here), and `tx` /
`version` give the candidate. For each safe W != V we resolve the same variant under
W and compare to that stored V coordinate.

Headline fact: N safe substitutions tested, K coordinate changes (expect K = 0).
Any K > 0 is a real bug in the safety check to chase. Written to
output/facts/version_safety_validation.csv. The benchmark's degraded `absent_version`
mode forces fallbacks but does NOT isolate the safe-identified subset, so it cannot
make this claim; this script can.

Use the SAME provider/build that built the A1 table, so `converted_g` (the V side)
and the W resolution are produced identically.

Per CLAUDE.md: never run against full datasets in development. Verify the logic
with a small table from tests/test_data/clinvar_hgvs/; the headline is a dedicated
run over the full-ClinVar A1 table.

Usage:
    python paper/scripts/validate_safe_versions.py \
        output/clinvar_pass/refseq.csv --rest \
        --out output/facts/version_safety_validation.csv
"""
import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
import hgvs.parser
from hgvs.assemblymapper import AssemblyMapper

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_resolution import build_provider, classify  # noqa: E402
from cdot.hgvs.gene_hgvs import _cited_position  # noqa: E402


def under_version(c_hgvs, tx, w):
    """Rewrite c_hgvs to cite version w of tx (keeps the c./n. body unchanged)."""
    rest = c_hgvs.split(":", 1)[1]
    return f"{tx}.{w}:{rest}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results", help="A1 per-variant table (resolve_clinvar_pass.py output)")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--rest", action="store_true", help="cdot REST (default)")
    g.add_argument("--rest-insecure", action="store_true")
    g.add_argument("--json", help="cdot JSON.gz file")
    g.add_argument("--uta", action="store_true")
    ap.add_argument("--build", default="GRCh38")
    ap.add_argument("--fasta", nargs="+", help="genome FASTA(s) for local sequence fetching")
    ap.add_argument("--out", default="output/facts/version_safety_validation.csv")
    ap.add_argument("--show-changes", action="store_true",
                    help="print every coordinate change (each is a safety-check bug)")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.debug else logging.WARNING)
    args.replace_reference = False  # build_provider reads this

    df = pd.read_csv(args.results, dtype=str).fillna("")
    # Only rows cdot resolved under the cited version V (converted_g is the V coord),
    # with a parseable integer version to substitute around.
    rows = df[df["bucket"].isin(["correct", "incorrect"])
              & (df["converted_g"] != "") & (df["version"] != "") & (df["tx"] != "")]
    print(f"== Validating safe substitutions over {len(rows)} resolved variants "
          f"(of {len(df)}) ==")

    provider, label = build_provider(args)
    if not hasattr(provider, "get_tx_versions"):
        sys.exit(f"provider {label} cannot enumerate transcript versions")
    am = AssemblyMapper(provider, assembly_name=args.build, alt_aln_method="splign",
                        replace_reference=False)
    hp = hgvs.parser.Parser()

    versions_cache = {}

    def tx_versions(tx):
        if tx not in versions_cache:
            try:
                versions_cache[tx] = provider.get_tx_versions(tx, args.build)
            except Exception as e:  # noqa: BLE001
                logging.debug("get_tx_versions(%s): %s", tx, e)
                versions_cache[tx] = []
        return versions_cache[tx]

    agg = {"variants_with_alt": 0, "safe_subs": 0, "unsafe_subs": 0,
           "resolved_under_w": 0, "coordinate_changes": 0, "unresolved_under_w": 0}
    changes = []
    for _, row in rows.iterrows():
        tx, v = row["tx"], int(row["version"])
        v_coord = row["converted_g"]
        others = [w for w in tx_versions(tx) if w != v]
        if not others:
            continue
        agg["variants_with_alt"] += 1
        cited_position = _cited_position(row["c_hgvs"])
        for w in others:
            safe, _reason = provider.is_version_substitution_safe(
                tx, v, w, args.build, cited_position=cited_position)
            if not safe:
                agg["unsafe_subs"] += 1
                continue
            agg["safe_subs"] += 1
            # V side reused from the table; resolve only the W side and compare.
            bucket, w_coord = classify(am, hp, v_coord, under_version(row["c_hgvs"], tx, w))
            if w_coord is None:
                agg["unresolved_under_w"] += 1
                continue
            agg["resolved_under_w"] += 1
            if bucket != "correct":  # classify scores against v_coord, so != correct means moved
                agg["coordinate_changes"] += 1
                changes.append((tx, v, w, v_coord, w_coord))
                if args.show_changes:
                    print(f"  CHANGE {tx} .{v}->.{w}: V={v_coord}  W={w_coord}")

    facts = {
        "n_variants_with_alt_version": [agg["variants_with_alt"]],
        "n_safe_substitutions":        [agg["safe_subs"]],
        "n_unsafe_substitutions":      [agg["unsafe_subs"]],
        "n_resolved_under_w":          [agg["resolved_under_w"]],
        "n_coordinate_changes":        [agg["coordinate_changes"]],
        "n_unresolved_under_w":        [agg["unresolved_under_w"]],
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(facts).to_csv(out, index=False)

    print(f"  provider              : {label}")
    print(f"  variants w/ alt version: {agg['variants_with_alt']}")
    print(f"  safe substitutions     : {agg['safe_subs']}  (tested)")
    print(f"  unsafe (skipped)       : {agg['unsafe_subs']}")
    print(f"  resolved under W       : {agg['resolved_under_w']}")
    print(f"  unresolved under W     : {agg['unresolved_under_w']}  (not a coord change)")
    print(f"  COORDINATE CHANGES (K) : {agg['coordinate_changes']}  (expect 0)")
    if changes:
        print("  !! safe substitutions that moved the coordinate (safety-check bugs):")
        for tx, v, w, vc, wc in changes[:20]:
            print(f"     {tx} .{v}->.{w}: {vc} -> {wc}")
    print(f"  Written: {out}")


if __name__ == "__main__":
    main()
