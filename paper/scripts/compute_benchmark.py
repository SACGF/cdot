#!/usr/bin/env python3
"""Throughput benchmark for Table 1: HGVS resolution speed by transcript backend.

Every configuration resolves the SAME fixed set of committed ClinVar (g.HGVS, c.HGVS)
pairs (tests/test_data/clinvar_hgvs/clinvar_hgvs_500.tsv) through the identical
biocommons/hgvs engine, with the sequence layer held constant: one shared local SeqRepo
SeqFetcher instance (file-descriptor caching enabled) is passed to every provider, so
the only thing that varies between configurations is the transcript-data layer.

Methodology (Results R3 / Methods):
  * N_REPEATS timed passes per configuration over the identical pair set.
  * Timed: the resolution loop only (parse + c_to_g per string, fresh AssemblyMapper
    per pass). Untimed: provider construction/load, REST prefetch(), and one warm-up
    pass per configuration (OS page cache, SeqRepo file handles).
  * The engine-side memo (the 100-entry LRU biocommons hgvs wraps around provider
    methods) is cleared before every timed pass - see clear_engine_caches().
  * Reported: median and IQR (Q1-Q3) of per-pass throughput.
  * Public remote UTA is ~0.1 HGVS/s, so a full-size pass is impractical; it is
    measured on the first UTA_REMOTE_N pairs of the same set (documented in Methods).

Writes output/facts/benchmark.csv, including the ratios and percent difference that the
R3 wording is generated from (see derived()).

Usage:
    python paper/scripts/compute_benchmark.py \
        --refseq-grch38 cdot.refseq.grch38.json.gz \
        [--uta-uri postgresql://postgres@127.0.0.1:5433/uta/uta_20241220] \
        [--skip-uta-remote]

    # recompute only the derived facts of an existing benchmark.csv (no benchmarking)
    python paper/scripts/compute_benchmark.py --derive-only paper/empirical_results/benchmark.csv

Requires: hgvs, cdot, a local SeqRepo (HGVS_SEQREPO_DIR).
"""

import argparse
import os
import statistics
import sys
import time
from pathlib import Path

# The sequence layer must be identical (and fast) for every configuration: enable
# SeqRepo file-descriptor caching before biocommons.seqrepo is imported. Without it
# every fetch re-opens a bgzf file (~1ms) and the sequence layer, not the transcript
# layer, dominates every row.
os.environ.setdefault("SEQREPO_FD_CACHE_MAXSIZE", "128")

from facts_io import read_facts, write_typed_facts  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PAIRS = REPO_ROOT / "tests/test_data/clinvar_hgvs/clinvar_hgvs_500.tsv"

N_REPEATS = 5
UTA_REMOTE_N = 25   # public remote UTA is ~0.1 HGVS/s; a 500-pair pass is impractical


def derived(tps):
    """The values the R3 wording is generated from, given median throughputs keyed like the
    facts (cdot_local_tps, ...): how many times faster local JSON is than local UTA, REST
    (one request per transcript) than the public remote UTA, and the percent difference
    between prefetched REST and local JSON."""
    out = {}
    if "cdot_local_tps" in tps and "uta_local_tps" in tps:
        out["cdot_local_vs_uta_local_ratio"] = tps["cdot_local_tps"] / tps["uta_local_tps"]
    if "cdot_rest_tps" in tps and "uta_remote_tps" in tps:
        out["cdot_rest_vs_uta_remote_ratio"] = tps["cdot_rest_tps"] / tps["uta_remote_tps"]
    if "cdot_rest_prefetch_tps" in tps and "cdot_local_tps" in tps:
        out["rest_prefetch_vs_local_diff_pct"] = (
            100 * abs(tps["cdot_rest_prefetch_tps"] - tps["cdot_local_tps"]) / tps["cdot_local_tps"])
    return out


def derive_only(csv_path):
    """Rewrite csv_path with its derived facts recomputed from the (rounded) medians
    it holds. A full run derives them from the unrounded medians instead."""
    facts = read_facts(csv_path)
    facts.update(derived({k: float(v) for k, v in facts.items() if k.endswith("_tps")}))
    write_typed_facts(csv_path, facts)
    print(f"Written: {csv_path}")


def clear_engine_caches(hdp):
    """biocommons hgvs wraps every data-provider method (get_tx_exons, get_seq, ...)
    in a per-instance 100-entry LRU at Interface construction. Clear that engine-side
    memo before each timed pass, so a pass measures the provider itself, including any
    internal caching the provider ships (cdot's in-memory transcript dict and exon
    cache, REST's prefetched transcript cache), rather than the engine memo. Without
    this, a working set smaller than 100 transcripts (the remote-UTA subset) is served
    entirely from the engine memo after the first pass and the backend is never hit."""
    for name in ("get_seq", "get_tx_exons", "get_tx_info", "get_tx_identity_info",
                 "get_tx_mapping_options", "get_pro_ac_for_tx_ac", "get_gene_info",
                 "get_acs_for_protein_seq", "get_tx_for_gene", "get_tx_for_region"):
        cache_clear = getattr(getattr(hdp, name, None), "cache_clear", None)
        if cache_clear is not None:
            cache_clear()


def resolve_pass(hdp, hp, pairs, build="GRCh38"):
    """One timed pass: parse + project every pair; fresh AssemblyMapper per pass.
    Returns (throughput_hgvs_per_s, counts)."""
    from hgvs.assemblymapper import AssemblyMapper
    from benchmark_resolution import classify

    am = AssemblyMapper(hdp, assembly_name=build, alt_aln_method="splign",
                        replace_reference=False)
    counts = {}
    t0 = time.perf_counter()
    for g_hgvs, c_hgvs in pairs:
        bucket, _ = classify(am, hp, g_hgvs, c_hgvs)
        counts[bucket] = counts.get(bucket, 0) + 1
    elapsed = time.perf_counter() - t0
    return len(pairs) / elapsed, counts


def run_config(name, make_provider, hp, pairs, repeats=N_REPEATS,
               fresh_provider_per_repeat=False, warmup=True):
    """Run one configuration: optional untimed warm-up pass, then `repeats` timed
    passes. `make_provider` is called once for persistent-provider configs, or once
    per repeat when the configuration is defined by a cold provider cache (REST)."""
    print(f"== {name} ==")
    provider = None
    if not fresh_provider_per_repeat:
        provider = make_provider()
    if warmup:
        warm_provider = provider if provider is not None else make_provider()
        tps, counts = resolve_pass(warm_provider, hp, pairs)
        print(f"  warm-up (untimed): {tps:7.1f} HGVS/s  {counts}")
    tps_values = []
    for i in range(repeats):
        hdp = provider if provider is not None else make_provider()
        clear_engine_caches(hdp)
        tps, counts = resolve_pass(hdp, hp, pairs)
        tps_values.append(tps)
        print(f"  repeat {i + 1}/{repeats}: {tps:7.1f} HGVS/s  {counts}")
    med = statistics.median(tps_values)
    q1, _q2, q3 = statistics.quantiles(tps_values, n=4, method="inclusive")
    print(f"  median {med:.1f}  IQR {q1:.1f}-{q3:.1f}")
    return {"median": med, "q1": q1, "q3": q3, "repeats": tps_values}


def main():
    parser = argparse.ArgumentParser(description="Benchmark HGVS resolution throughput by backend.")
    parser.add_argument("--refseq-grch38", help="cdot GRCh38 RefSeq JSON.gz")
    parser.add_argument("--uta-uri", default=None, help="local UTA PostgreSQL URI")
    parser.add_argument("--pairs", default=str(DEFAULT_PAIRS), help="fixed (g.HGVS, c.HGVS) pair set")
    parser.add_argument("--repeats", type=int, default=N_REPEATS)
    parser.add_argument("--skip-uta-remote", action="store_true",
                        help="skip the (slow, ~25 min) public remote UTA rows")
    parser.add_argument("--skip-rest", action="store_true", help="skip the REST configurations")
    parser.add_argument("--derive-only", metavar="CSV",
                        help="only recompute the derived facts of an existing benchmark.csv")
    args = parser.parse_args()
    if args.derive_only:
        derive_only(args.derive_only)
        return
    if not args.refseq_grch38:
        parser.error("--refseq-grch38 is required")

    sys.path.insert(0, str(Path(__file__).resolve().parent))

    import hgvs.parser
    from hgvs.dataproviders.seqfetcher import SeqFetcher
    from cdot.hgvs.dataproviders import JSONDataProvider, RESTDataProvider
    from benchmark_resolution import load_pairs

    if not os.environ.get("HGVS_SEQREPO_DIR"):
        sys.exit("Set HGVS_SEQREPO_DIR to a local SeqRepo: the sequence layer must be "
                 "local and identical for every configuration.")

    pairs = load_pairs(args.pairs)
    hp = hgvs.parser.Parser()
    # One shared SeqRepo seqfetcher for every provider: SeqRepo's fd cache is
    # per-instance, so per-provider seqfetchers would give configurations
    # differently-warmed sequence layers.
    shared_seqfetcher = SeqFetcher()

    facts = {"n_pairs": len(pairs), "n_repeats": args.repeats}
    medians = {}  # unrounded, for derived()

    # --- cdot local JSON -------------------------------------------------------
    t0 = time.perf_counter()
    json_provider = JSONDataProvider([args.refseq_grch38], seqfetcher=shared_seqfetcher)
    load_time_s = time.perf_counter() - t0
    print(f"cdot local JSON load: {load_time_s:.1f}s")
    r = run_config("cdot local JSON", lambda: json_provider, hp, pairs, args.repeats)
    medians["cdot_local_tps"] = r["median"]
    facts.update(cdot_local_tps=round(r["median"]), cdot_local_tps_q1=round(r["q1"]),
                 cdot_local_tps_q3=round(r["q3"]), grch38_load_time_s=round(load_time_s, 1))

    # --- cdot REST -------------------------------------------------------------
    if not args.skip_rest:
        def make_rest():
            return RESTDataProvider(seqfetcher=shared_seqfetcher)

        r = run_config("cdot REST (cold cache, one request per transcript)", make_rest,
                       hp, pairs, args.repeats, fresh_provider_per_repeat=True)
        medians["cdot_rest_tps"] = r["median"]
        facts.update(cdot_rest_tps=round(r["median"]), cdot_rest_tps_q1=round(r["q1"]),
                     cdot_rest_tps_q3=round(r["q3"]))

        tx_acs = sorted({c.split(":", 1)[0] for _g, c in pairs if ":" in c})

        def make_rest_prefetched():
            hdp = RESTDataProvider(seqfetcher=shared_seqfetcher)
            t0 = time.perf_counter()
            n = hdp.prefetch(tx_acs)
            print(f"  prefetched {n}/{len(tx_acs)} transcripts in {time.perf_counter() - t0:.2f}s (untimed)")
            return hdp

        r = run_config("cdot REST (after one batch prefetch())", make_rest_prefetched,
                       hp, pairs, args.repeats, fresh_provider_per_repeat=True)
        medians["cdot_rest_prefetch_tps"] = r["median"]
        facts.update(cdot_rest_prefetch_tps=round(r["median"]),
                     cdot_rest_prefetch_tps_q1=round(r["q1"]),
                     cdot_rest_prefetch_tps_q3=round(r["q3"]))

    # --- UTA local -------------------------------------------------------------
    import hgvs.dataproviders.uta as uta
    if args.uta_uri:
        def make_uta_local():
            hdp = uta.connect(db_url=args.uta_uri)
            hdp.seqfetcher = shared_seqfetcher  # identical sequence layer
            return hdp

        r = run_config("UTA local PostgreSQL", make_uta_local, hp, pairs, args.repeats)
        medians["uta_local_tps"] = r["median"]
        facts.update(uta_local_tps=round(r["median"]), uta_local_tps_q1=round(r["q1"]),
                     uta_local_tps_q3=round(r["q3"]))

    # --- UTA public remote (subset: ~0.1 HGVS/s) -------------------------------
    if not args.skip_uta_remote:
        remote_pairs = pairs[:UTA_REMOTE_N]

        def make_uta_remote():
            hdp = uta.connect()  # default public uta.biocommons.org
            hdp.seqfetcher = shared_seqfetcher
            return hdp

        print(f"(public remote UTA measured on the first {len(remote_pairs)} pairs of the same set)")
        r = run_config("UTA public remote", make_uta_remote, hp, remote_pairs,
                       args.repeats, warmup=False)
        medians["uta_remote_tps"] = r["median"]
        facts.update(uta_remote_tps=round(r["median"], 2),
                     uta_remote_tps_q1=round(r["q1"], 2),
                     uta_remote_tps_q3=round(r["q3"], 2),
                     uta_remote_n=len(remote_pairs))

    out = Path("output/facts/benchmark.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    facts.update(derived(medians))
    write_typed_facts(out, facts)
    print(f"Written: {out}")


if __name__ == "__main__":
    main()
