#!/usr/bin/env python3
"""Full-scale local-JSON throughput (``benchmark_fullscale.csv``), averaged over runs.

Runs ``resolve_clinvar_pass.py`` over the whole ClinVar pair set ``--runs + 1`` times,
each in a fresh process. The first pass warms the OS sequence cache and is reported as
``cold_tps``; ``hot_tps`` is the mean of the remaining passes, with their range, since a
single 2 h pass on a shared VM varies by several percent between days.

    HGVS_SEQREPO_DIR=... SEQREPO_FD_CACHE_MAXSIZE=128 \\
    python paper/scripts/measure_fullscale.py clinvar_pairs.GRCh38.tsv \\
        --json cdot-0.2.36.refseq.GRCh38.json.gz --runs 3 --out benchmark_fullscale.csv
"""
import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path
from statistics import mean

PASS_SCRIPT = Path(__file__).with_name("resolve_clinvar_pass.py")


def run_pass(pairs_file, json_file):
    out = subprocess.run([sys.executable, str(PASS_SCRIPT), pairs_file, "--json", json_file,
                          "--out", "/dev/null"], check=True, capture_output=True, text=True).stdout
    print(out, flush=True)
    tps, wall = re.search(r"throughput\s*: ([\d.]+) HGVS/s\s+\(wall ([\d.]+)s\)", out).groups()
    resolved_pct = re.search(r"resolved \(any g\.\)\s*:\s*\d+\s+\(([\d.]+)%\)", out).group(1)
    n_pairs = re.search(r"\((\d+) rows\)", out).group(1)
    return {"tps": float(tps), "wall": float(wall), "resolved_pct": resolved_pct, "n_pairs": n_pairs}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pairs_file")
    ap.add_argument("--json", required=True, help="cdot RefSeq GRCh38 JSON.gz")
    ap.add_argument("--runs", type=int, default=3, help="timed passes after the warm-up pass")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cold = run_pass(args.pairs_file, args.json)
    hot = [run_pass(args.pairs_file, args.json) for _ in range(args.runs)]
    hot_tps = [h["tps"] for h in hot]
    row = {
        "n_pairs": cold["n_pairs"],
        "resolved_pct": cold["resolved_pct"],
        "cold_tps": round(cold["tps"]),
        "hot_tps": round(mean(hot_tps)),
        "hot_tps_min": round(min(hot_tps)),
        "hot_tps_max": round(max(hot_tps)),
        "n_runs": args.runs,
        "wall_min": round(mean(h["wall"] for h in hot) / 60),
    }
    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["field", "value"])
        w.writerows(row.items())
    print(row)


if __name__ == "__main__":
    main()
