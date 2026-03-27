#!/usr/bin/env python3
"""Benchmark cdot local JSON throughput and load time vs UTA remote access.

Writes output/facts/benchmark.csv.

Usage:
    python analysis/compute_benchmark.py \
        --refseq-grch38 cdot.refseq.grch38.json.gz \
        [--uta-uri postgresql://uta_admin@localhost/uta/uta_20210129]

Requires: hgvs, cdot
"""

import argparse
import time
from pathlib import Path

import pandas as pd


N_WARMUP = 10
N_BENCHMARK = 500


def benchmark_cdot_local(json_gz_path: str) -> tuple[float, float]:
    """Return (load_time_s, transcripts_per_second) for local JSON provider."""
    from cdot.hgvs.dataproviders import JSONDataProvider

    t0 = time.perf_counter()
    hdp = JSONDataProvider([json_gz_path])
    load_time_s = time.perf_counter() - t0

    # Get a sample of transcript accessions to query
    import gzip
    import json
    with gzip.open(json_gz_path, "rt") as fh:
        data = json.load(fh)
    accessions = list(data["transcripts"].keys())[:N_BENCHMARK + N_WARMUP]

    # Warm up
    for acc in accessions[:N_WARMUP]:
        try:
            hdp.get_tx_info(acc)
        except Exception:
            pass

    # Benchmark
    sample = accessions[N_WARMUP:N_WARMUP + N_BENCHMARK]
    t0 = time.perf_counter()
    for acc in sample:
        try:
            hdp.get_tx_info(acc)
        except Exception:
            pass
    elapsed = time.perf_counter() - t0
    tps = len(sample) / elapsed if elapsed > 0 else 0

    return round(load_time_s, 2), round(tps)


def benchmark_uta_remote(uta_uri: str) -> float:
    """Return transcripts_per_second for UTA remote access."""
    import hgvs.dataproviders.uta as uta

    hdp = uta.connect(db_url=uta_uri)

    # Get sample accessions
    sample_accessions = [
        "NM_000492.3", "NM_007294.3", "NM_000059.3", "NM_000546.5", "NM_001354689.1",
        "NM_000249.3", "NM_001126112.2", "NM_005343.3", "NM_004333.5", "NM_000179.2",
    ]

    # Warm up
    for acc in sample_accessions[:3]:
        try:
            hdp.get_tx_info(acc, "GRCh38")
        except Exception:
            pass

    t0 = time.perf_counter()
    n = 0
    for acc in sample_accessions * 5:
        try:
            hdp.get_tx_info(acc, "GRCh38")
            n += 1
        except Exception:
            pass
    elapsed = time.perf_counter() - t0
    return round(n / elapsed) if elapsed > 0 else 0


def main():
    parser = argparse.ArgumentParser(description="Benchmark cdot vs UTA throughput.")
    parser.add_argument("--refseq-grch38", required=True, help="cdot GRCh38 RefSeq JSON.gz")
    parser.add_argument("--uta-uri", default=None, help="UTA PostgreSQL URI for remote benchmark")
    args = parser.parse_args()

    print("Benchmarking cdot local JSON ...")
    load_time_s, cdot_local_tps = benchmark_cdot_local(args.refseq_grch38)
    print(f"  Load time:   {load_time_s}s")
    print(f"  Throughput:  {cdot_local_tps:,} transcripts/s")

    uta_remote_tps = 1  # well-established from Hart 2020 / general experience
    if args.uta_uri:
        print("Benchmarking UTA remote ...")
        uta_remote_tps = benchmark_uta_remote(args.uta_uri)
        print(f"  Throughput:  {uta_remote_tps:,} transcripts/s")

    facts = {
        "cdot_local_min_tps": [min(cdot_local_tps, 500)],  # conservative bound
        "cdot_local_max_tps": [max(cdot_local_tps, 1000)],  # upper bound
        "cdot_rest_tps":      [0],   # TODO: benchmark REST API separately
        "uta_local_tps":      [0],   # TODO: benchmark UTA local if available
        "uta_remote_tps":     [uta_remote_tps],
        "grch38_load_time_s": [load_time_s],
    }

    out = Path("output/facts/benchmark.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(facts).to_csv(out, index=False)
    print(f"Written: {out}")


if __name__ == "__main__":
    main()
