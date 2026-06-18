#!/usr/bin/env python3
"""Print a side-by-side table from benchmark_resolution.py --out JSON files.

Usage:
    python analysis/compare_providers.py output/facts/compare/*.json
"""
import json
import sys
from pathlib import Path


def pct(n, d):
    return f"{100 * n / d:.1f}%" if d else "n/a"


def main():
    paths = sys.argv[1:]
    rows = [json.loads(Path(p).read_text()) for p in paths]
    rows.sort(key=lambda r: r["provider"])

    hdr = f"{'provider':<34} {'n':>4} {'correct':>9} {'no_data':>9} {'error':>7} {'load_s':>7} {'pref_s':>7} {'HGVS/s':>8}"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        c = r["counts"]
        n = r["n"]
        print(f"{r['provider']:<34} {n:>4} "
              f"{pct(c['correct'], n):>9} {pct(c['no_data'], n):>9} {c['error']:>7} "
              f"{r.get('load_time_s', 0):>7.2f} "
              f"{(r.get('prefetch_s') or 0):>7.2f} "
              f"{r.get('throughput_hgvs_per_s', 0):>8.1f}")


if __name__ == "__main__":
    main()
