#!/usr/bin/env python3
"""Resolution / recovery / speed benchmark for cdot.

Answers the three paper questions over real ClinVar (g.HGVS, c.HGVS) pairs:

  1. RESOLUTION  - of real c.HGVS, how many convert to the correct genomic coord?
  2. RECOVERY    - how many failing strings does fix_hgvs() rescue, and via which mode?
  3. SPEED       - throughput (HGVS/s) and provider load time.

Ground truth: the g.HGVS column (from ClinVar CLNHGVS). A c->g conversion is "correct"
if it reproduces that g.HGVS.

Recovery (--recovery) routes every string through the single combined entry point
fix_hgvs() (clean_hgvs -> resolve_gene_hgvs -> version-bump fallback) and attributes
each rescue to a failure mode by reading the HGVSFix codes it returns. Two passes:
  * as-is    - run fix_hgvs() on the pairs unchanged; rescues + false-rescues + a
               regression guard. This is the production-mode measurement: on clean
               input nothing fires, so feed messy data to see real per-mode numbers.
  * degraded - degrade each clean pair per mode (inject a formatting error / request
               an absent transcript version), then fix_hgvs() -> resolve and score
               against ground truth, with a false-rescue count. Exercises the rescue
               path on the (already-clean) committed test data.
Note: the "more transcripts" (coverage) mode is NOT a string repair - it needs a
provider swap (run with --uta vs cdot and compare no_data), not fix_hgvs().

Provider is pluggable so we isolate the data-layer contribution (the biocommons engine
is held constant). Defaults to cdot REST (cdotlib.org) so it runs with no local data.

Usage:
    # quick: REST over the committed test pairs
    python analysis/benchmark_resolution.py tests/test_data/clinvar_hgvs/clinvar_hgvs_100.tsv --rest

    # local JSON (after downloading a release json.gz)
    python analysis/benchmark_resolution.py pairs.tsv --json cdot-X.refseq.grch38.json.gz

    # add the unified fix_hgvs() recovery passes (as-is + degraded)
    python analysis/benchmark_resolution.py pairs.tsv --rest --recovery

Input file: TSV with two columns  g.HGVS <TAB> c.HGVS  (one variant per line).
"""
import argparse
import json
import logging
import re
import time
from pathlib import Path

import hgvs.dataproviders.uta
import hgvs.parser
from hgvs.assemblymapper import AssemblyMapper
from hgvs.exceptions import HGVSDataNotAvailableError, HGVSInvalidVariantError

from cdot.hgvs.dataproviders import JSONDataProvider, RESTDataProvider, FastaSeqFetcher
from cdot.hgvs.clean import HGVSFixCode, HGVSFixSeverity, VersionStrategy
from cdot.hgvs.gene_hgvs import fix_hgvs


def build_provider(args):
    """Return (provider, label). replace_reference is handled by the caller."""
    # Local FASTA seqfetcher avoids remote NCBI sequence fetching (the real speed bottleneck).
    seqfetcher = FastaSeqFetcher(*args.fasta) if args.fasta else None
    seq_label = "fasta" if args.fasta else "ncbi"
    if args.uta:
        # NOTE: do NOT attach cdot's --fasta FastaSeqFetcher to UTA. That fetcher
        # reconstructs transcript (NM_) sequence from the genome using *cdot's*
        # exon data; bare UTA has none, so every transcript-sequence request
        # raises HGVSDataNotAvailableError and is mis-counted as no_data (it made
        # a clean UTA run look like 85% no_data — a pure artifact). Stock UTA must
        # use its own SeqRepo/remote seqfetcher. Set HGVS_SEQREPO_DIR for a fast,
        # offline, transcript-capable baseline; otherwise it falls back to the
        # (slow) remote bioutils seqfetcher.
        if args.fasta:
            print("  (--fasta ignored for --uta: cdot's genome fetcher can't serve "
                  "transcript sequence; UTA uses SeqRepo/remote)")
        provider = hgvs.dataproviders.uta.connect()
        import os
        seq_src = "seqrepo" if os.environ.get("HGVS_SEQREPO_DIR") else "remote"
        uta_url = os.environ.get("UTA_DB_URL", "")
        loc = "local" if ("localhost" in uta_url or "127.0.0.1" in uta_url) else "remote"
        return provider, f"uta_{loc}/{seq_src}"
    if args.json:
        return JSONDataProvider([args.json], seqfetcher=seqfetcher), f"cdot_json({Path(args.json).name})/{seq_label}"
    # default: REST
    return RESTDataProvider(secure=not args.rest_insecure, seqfetcher=seqfetcher), f"cdot_rest(cdotlib.org)/{seq_label}"


def load_pairs(path):
    pairs = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) >= 2:
                pairs.append((parts[0], parts[1]))
    return pairs


def classify(am, hp, g_hgvs, c_hgvs):
    """Return (bucket, converted_or_None). buckets: correct/incorrect/no_data/error."""
    try:
        var_c = hp.parse_hgvs_variant(c_hgvs)
        if ":c." in c_hgvs:
            converted = str(am.c_to_g(var_c))
        else:
            converted = str(am.n_to_g(var_c))
    except HGVSDataNotAvailableError:
        return "no_data", None
    except HGVSInvalidVariantError:
        return "error", None
    except Exception as e:  # noqa: BLE001 - benchmark wants to keep going
        logging.debug("error on %s: %s", c_hgvs, e)
        return "error", None
    return ("correct" if converted == g_hgvs else "incorrect"), converted


def run_resolution(am, hp, pairs, show_incorrect=False):
    counts = {"correct": 0, "incorrect": 0, "no_data": 0, "error": 0}
    times = []
    t0 = time.perf_counter()
    for g_hgvs, c_hgvs in pairs:
        s = time.perf_counter()
        bucket, converted = classify(am, hp, g_hgvs, c_hgvs)
        times.append(time.perf_counter() - s)
        counts[bucket] += 1
        if show_incorrect and bucket == "incorrect":
            print(f"  INCORRECT {c_hgvs}: expected {g_hgvs} got {converted}")
        if bucket == "no_data":
            logging.debug("no_data: %s", c_hgvs)
    wall = time.perf_counter() - t0
    return counts, times, wall


# ---- unified recovery via fix_hgvs ------------------------------------------
#
# fix_hgvs() is the single combined entry point: clean_hgvs -> resolve_gene_hgvs
# -> resolve_transcript_version (version bump). It returns the repaired string
# plus the list of HGVSFix applied, so we attribute each rescue to a failure
# mode just by reading the fix codes returned — no separate per-mode pipeline.

# WARNING-level fix codes grouped by the failure mode they address.
_VERSION_CODES = {HGVSFixCode.USED_ADJACENT_VERSION}
_GENE_CODES = {HGVSFixCode.RESOLVED_GENE_TO_TRANSCRIPT, HGVSFixCode.UPPERCASED_GENE_SYMBOL}


def _mode_of(code):
    if code in _VERSION_CODES:
        return "version_bump"
    if code in _GENE_CODES:
        return "gene_resolution"
    return "cleaning"  # every other WARNING code is a clean_hgvs() string fix


def _fix_and_classify(am, hp, provider, build, g_hgvs, c_input, version_fallback):
    """Run fix_hgvs() then resolve. Returns (bucket, fixed_string, modes_fired)."""
    try:
        fixed, fixes = fix_hgvs(c_input, provider, build, version_fallback=version_fallback)
    except Exception as e:  # noqa: BLE001 - benchmark keeps going
        logging.debug("fix_hgvs error on %s: %s", c_input, e)
        return "error", c_input, set()
    modes = {_mode_of(f.code) for f in fixes if f.severity == HGVSFixSeverity.WARNING}
    bucket, _converted = classify(am, hp, g_hgvs, fixed)
    return bucket, fixed, modes


def run_asis(am, hp, provider, build, pairs, version_fallback):
    """Production-mode pass: run fix_hgvs() on each pair as-is and read what the
    returned fixes rescued. Also guards against regressions (a string that
    resolved correctly before fix_hgvs but not after)."""
    agg = {"n": len(pairs), "baseline_correct": 0, "fixed_correct": 0,
           "rescued": 0, "regressed": 0, "false_rescue": 0,
           "by_mode": {}}
    for g_hgvs, c_hgvs in pairs:
        base_bucket, _ = classify(am, hp, g_hgvs, c_hgvs)
        fbucket, _fixed, modes = _fix_and_classify(
            am, hp, provider, build, g_hgvs, c_hgvs, version_fallback)
        if base_bucket == "correct":
            agg["baseline_correct"] += 1
        if fbucket == "correct":
            agg["fixed_correct"] += 1
        for m in modes:
            agg["by_mode"].setdefault(m, {"fired": 0, "rescued": 0, "false_rescue": 0})["fired"] += 1
        if base_bucket != "correct" and fbucket == "correct":
            agg["rescued"] += 1
            for m in modes:
                agg["by_mode"][m]["rescued"] += 1
        if base_bucket == "correct" and fbucket != "correct":
            agg["regressed"] += 1
        if fbucket == "incorrect":
            agg["false_rescue"] += 1
            for m in modes:
                agg["by_mode"][m]["false_rescue"] += 1
    return agg


def _degrade(c_hgvs):
    """Yield (mode, category, degraded_string) — only degradations that keep the
    SAME transcript, so the ground-truth g.HGVS still applies after repair."""
    m = re.match(r"([A-Za-z]+_?\d+)\.(\d+):(c\..+)", c_hgvs)
    if not m:
        return
    base, ver, rest = m.group(1), m.group(2), m.group(3)
    tx = f"{base}.{ver}"
    # cleaning — string-formatting damage clean_hgvs() should fully reverse
    yield "cleaning", "whitespace", f"{tx} {rest}"
    yield "cleaning", "lowercase_tx", f"{tx.lower()}:{rest}"
    yield "cleaning", "missing_kind", f"{tx}:{rest[2:]}"   # drop the "c."
    yield "cleaning", "double_colon", f"{tx}::{rest}"
    # version bump — request a version that is (almost certainly) absent, forcing
    # the adjacent-version fallback; correct iff the chosen version reproduces g.
    yield "version_bump", "absent_version", f"{base}.99:{rest}"


def run_degraded(am, hp, provider, build, pairs, version_fallback):
    """Exercise the fix_hgvs() rescue path on (clean) test data by degrading each
    pair per mode, then scoring the repaired string against the ground-truth g."""
    cats = {}
    for g_hgvs, c_hgvs in pairs:
        for mode, cat, broken in _degrade(c_hgvs):
            d = cats.setdefault((mode, cat),
                                {"n": 0, "resolved": 0, "correct": 0, "false_rescue": 0})
            d["n"] += 1
            bucket, _fixed, _modes = _fix_and_classify(
                am, hp, provider, build, g_hgvs, broken, version_fallback)
            if bucket in ("correct", "incorrect"):
                d["resolved"] += 1
            if bucket == "correct":
                d["correct"] += 1
            elif bucket == "incorrect":
                d["false_rescue"] += 1
    return cats


def pct(n, d):
    return f"{100 * n / d:.1f}%" if d else "n/a"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("pairs_file", help="TSV: g.HGVS <TAB> c.HGVS")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--rest", action="store_true", help="cdot REST (default)")
    g.add_argument("--rest-insecure", action="store_true")
    g.add_argument("--json", help="cdot JSON.gz file")
    g.add_argument("--uta", action="store_true")
    ap.add_argument("--build", default="GRCh38")
    ap.add_argument("--fasta", nargs="+", help="genome FASTA(s) for local sequence fetching (avoids remote NCBI)")
    ap.add_argument("--replace-reference", action="store_true",
                    help="validate ref base (needs a seqfetcher); off = pure coordinate projection")
    ap.add_argument("--prefetch", action="store_true", help="(REST) warm the transcript cache concurrently before resolving")
    ap.add_argument("--prefetch-workers", type=int, default=10)
    ap.add_argument("--recovery", action="store_true", help="also run cleaning + version-bump experiments")
    ap.add_argument("--show-incorrect", action="store_true")
    ap.add_argument("--out", help="write JSON results here")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.debug else logging.WARNING)

    pairs = load_pairs(args.pairs_file)

    t0 = time.perf_counter()
    provider, label = build_provider(args)
    load_time = time.perf_counter() - t0

    am = AssemblyMapper(provider, assembly_name=args.build, alt_aln_method="splign",
                        replace_reference=args.replace_reference)
    hp = hgvs.parser.Parser()

    prefetch_time = None
    if args.prefetch:
        if not hasattr(provider, "prefetch"):
            print(f"  (provider {label} has no prefetch(); ignoring --prefetch)")
        else:
            tx_acs = sorted({c.split(":", 1)[0] for _g, c in pairs if ":" in c})
            t0 = time.perf_counter()
            n = provider.prefetch(tx_acs, max_workers=args.prefetch_workers)
            prefetch_time = time.perf_counter() - t0
            print(f"  prefetched {n}/{len(tx_acs)} transcripts in {prefetch_time:.2f}s "
                  f"({args.prefetch_workers} workers)")

    print(f"== Resolution: {label}  |  {len(pairs)} pairs  |  {args.build}  "
          f"|  replace_reference={args.replace_reference} ==")
    counts, times, wall = run_resolution(am, hp, pairs, args.show_incorrect)
    resolved = counts["correct"] + counts["incorrect"]
    total = len(pairs)
    print(f"  load_time         : {load_time:.2f}s")
    print(f"  correct           : {counts['correct']:>5}  ({pct(counts['correct'], total)})")
    print(f"  incorrect         : {counts['incorrect']:>5}  ({pct(counts['incorrect'], total)})")
    print(f"  no_data (missing) : {counts['no_data']:>5}  ({pct(counts['no_data'], total)})")
    print(f"  error             : {counts['error']:>5}  ({pct(counts['error'], total)})")
    print(f"  resolved (any g.) : {resolved:>5}  ({pct(resolved, total)})")
    print(f"  throughput        : {total / wall:.1f} HGVS/s  (wall {wall:.2f}s)")

    results = {"provider": label, "build": args.build, "n": total,
               "replace_reference": args.replace_reference,
               "load_time_s": round(load_time, 2), "counts": counts,
               "prefetch_s": round(prefetch_time, 2) if prefetch_time is not None else None,
               "throughput_hgvs_per_s": round(total / wall, 1), "wall_s": round(wall, 2)}

    if args.recovery:
        vf = VersionStrategy.UP_THEN_DOWN

        print("\n== Recovery (as-is): fix_hgvs() on each pair, attributed by the fixes returned ==")
        asis = run_asis(am, hp, provider, args.build, pairs, vf)
        print(f"  baseline correct  : {asis['baseline_correct']:>5}  ({pct(asis['baseline_correct'], asis['n'])})")
        print(f"  after fix_hgvs    : {asis['fixed_correct']:>5}  ({pct(asis['fixed_correct'], asis['n'])})")
        print(f"  rescued by fix    : {asis['rescued']:>5}  ({pct(asis['rescued'], asis['n'])})")
        print(f"  false-rescue      : {asis['false_rescue']:>5}  (resolved, but to the wrong g.)")
        print(f"  regressed         : {asis['regressed']:>5}  (was correct, broken by fix_hgvs)")
        for mode, d in sorted(asis["by_mode"].items()):
            print(f"    mode {mode:<16} fired={d['fired']:<4} rescued={d['rescued']:<4} false_rescue={d['false_rescue']}")
        if not asis["by_mode"]:
            print("    (no fixes fired — test pairs are already clean; see the degraded pass below)")
        results["recovery_asis"] = asis

        print("\n== Recovery (degraded): degrade per mode -> fix_hgvs() -> resolve vs ground truth ==")
        degraded = run_degraded(am, hp, provider, args.build, pairs, vf)
        for (mode, cat), d in sorted(degraded.items()):
            print(f"  {mode:<14} {cat:<16} n={d['n']:<4} "
                  f"correct={pct(d['correct'], d['n']):>6} "
                  f"false_rescue={pct(d['false_rescue'], d['n']):>6}")
        results["recovery_degraded"] = {f"{m}/{c}": d for (m, c), d in degraded.items()}

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(results, indent=2))
        print(f"\nWritten: {args.out}")


if __name__ == "__main__":
    main()
