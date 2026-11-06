#!/usr/bin/env python3
"""Full-scale local-JSON throughput (``benchmark_fullscale.csv``), averaged over runs.

Runs ``resolve_clinvar_pass.py`` over the whole ClinVar pair set ``--runs + 1`` times,
each in a fresh process. The first pass warms the OS sequence cache and is reported as
``cold_tps``; ``hot_tps`` is the mean of the remaining passes, with their range, since a
single 2 h pass on a shared VM varies by several percent between days.

Also writes the values the R3 wording is generated from (see derived()): the percent
difference between the hot and cold passes and, from the public remote UTA throughput in
``--benchmark-csv`` (Table 1), how many days the same pass would take against that server.

    HGVS_SEQREPO_DIR=... SEQREPO_FD_CACHE_MAXSIZE=128 \\
    python paper/scripts/measure_fullscale.py clinvar_pairs.GRCh38.tsv \\
        --json cdot-0.2.36.refseq.GRCh38.json.gz --runs 3 --out benchmark_fullscale.csv

    # recompute only the derived rows of an existing CSV (no resolution passes)
    python paper/scripts/measure_fullscale.py --derive-only \\
        --out paper/empirical_results/benchmark_fullscale.csv \\
        --benchmark-csv paper/empirical_results/benchmark.csv
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


def derived(n_pairs, cold_tps, hot_tps, benchmark_csv):
    out = {"hot_vs_cold_diff_pct": 100 * abs(hot_tps - cold_tps) / cold_tps}
    if benchmark_csv and Path(benchmark_csv).exists():
        with open(benchmark_csv, newline="") as f:
            uta_remote_tps = float(next(csv.DictReader(f))["uta_remote_tps"])
        out["uta_remote_days"] = int(n_pairs) / uta_remote_tps / 86400
    return out


def write_rows(path, row):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["field", "value"])
        w.writerows(row.items())


def derive_only(path, benchmark_csv):
    with open(path, newline="") as f:
        row = {r["field"]: r["value"] for r in csv.DictReader(f)}
    row.update(derived(row["n_pairs"], float(row["cold_tps"]), float(row["hot_tps"]), benchmark_csv))
    write_rows(path, row)
    print(row)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pairs_file", nargs="?")
    ap.add_argument("--json", help="cdot RefSeq GRCh38 JSON.gz")
    ap.add_argument("--runs", type=int, default=3, help="timed passes after the warm-up pass")
    ap.add_argument("--out", required=True)
    ap.add_argument("--benchmark-csv", default="output/facts/benchmark.csv",
                    help="Table 1 facts, for the public remote UTA throughput")
    ap.add_argument("--derive-only", action="store_true",
                    help="only recompute the derived rows of an existing --out CSV")
    args = ap.parse_args()
    if args.derive_only:
        derive_only(args.out, args.benchmark_csv)
        return
    if not (args.pairs_file and args.json):
        ap.error("pairs_file and --json are required")

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
    row.update(derived(cold["n_pairs"], cold["tps"], mean(hot_tps), args.benchmark_csv))
    write_rows(args.out, row)
    print(row)


if __name__ == "__main__":
    main()
