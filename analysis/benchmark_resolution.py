#!/usr/bin/env python3
"""Resolution / recovery / speed benchmark for cdot.

Answers the three paper questions over real ClinVar (g.HGVS, c.HGVS) pairs:

  1. RESOLUTION  - of real c.HGVS, how many convert to the correct genomic coord?
  2. RECOVERY    - of strings that fail, how many do clean_hgvs() / version-bumping rescue?
  3. SPEED       - throughput (HGVS/s) and provider load time.

Ground truth: the g.HGVS column (from ClinVar CLNHGVS). A c->g conversion is "correct"
if it reproduces that g.HGVS.

Provider is pluggable so we isolate the data-layer contribution (the biocommons engine
is held constant). Defaults to cdot REST (cdotlib.org) so it runs with no local data.

Usage:
    # quick: REST over the committed test pairs
    python analysis/benchmark_resolution.py tests/test_data/clinvar_hgvs/clinvar_hgvs_100.tsv --rest

    # local JSON (after downloading a release json.gz)
    python analysis/benchmark_resolution.py pairs.tsv --json cdot-X.refseq.grch38.json.gz

    # add the recovery experiments
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
from cdot.hgvs.clean import clean_hgvs, get_best_transcript_version, VersionStrategy


def build_provider(args):
    """Return (provider, label). replace_reference is handled by the caller."""
    # Local FASTA seqfetcher avoids remote NCBI sequence fetching (the real speed bottleneck).
    seqfetcher = FastaSeqFetcher(*args.fasta) if args.fasta else None
    seq_label = "fasta" if args.fasta else "ncbi"
    if args.uta:
        return hgvs.dataproviders.uta.connect(), f"uta_remote/{seq_label}"
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


# ---- recovery: cleaning -----------------------------------------------------

def _corrupt_variants(c_hgvs):
    """Yield (category, corrupted_string) for each clean.py fix category."""
    m = re.match(r"([A-Za-z]+_?\d+\.\d+):(c\..+)", c_hgvs)
    if not m:
        return
    tx, rest = m.groups()
    yield "stripped_whitespace", f"{tx} {rest}"          # space instead of colon
    yield "lowercase_transcript", f"{tx.lower()}:{rest}"
    yield "missing_kind", f"{tx}:{rest[2:]}"             # drop the "c."
    yield "double_colon", f"{tx}: :{rest}"
    yield "unbalanced_bracket", f"{tx}):{rest}"


def run_cleaning(am, hp, pairs):
    """For each variant, corrupt it per category, clean it, check it now resolves."""
    per_cat = {}
    for _g, c_hgvs in pairs:
        for cat, broken in _corrupt_variants(c_hgvs):
            d = per_cat.setdefault(cat, {"n": 0, "resolved_raw": 0, "cleaned_ok": 0, "resolved_clean": 0})
            d["n"] += 1
            # does the broken string parse as-is?
            try:
                hp.parse_hgvs_variant(broken)
                d["resolved_raw"] += 1
            except Exception:
                pass
            cleaned, _fixes = clean_hgvs(broken)
            if cleaned == c_hgvs:
                d["cleaned_ok"] += 1
            try:
                hp.parse_hgvs_variant(cleaned)
                d["resolved_clean"] += 1
            except Exception:
                pass
    return per_cat


# ---- recovery: version bumping ----------------------------------------------

def run_version_bump(provider, pairs):
    """Drop the exact requested version; does get_best_transcript_version find a neighbour?"""
    stats = {"n": 0, "had_alternatives": 0, "rescued": 0, "no_alternatives": 0}
    seen = set()
    for _g, c_hgvs in pairs:
        m = re.match(r"([A-Za-z]+_?\d+)\.(\d+):", c_hgvs)
        if not m:
            continue
        base, ver = m.group(1), int(m.group(2))
        if base in seen:
            continue
        seen.add(base)
        # discover available versions by probing the provider
        available = []
        for v in range(1, 12):
            try:
                provider._get_transcript(f"{base}.{v}")
                available.append(v)
            except Exception:
                pass
        if not available:
            continue
        stats["n"] += 1
        # simulate the requested version being absent
        remaining = [v for v in available if v != ver]
        if not remaining:
            stats["no_alternatives"] += 1
            continue
        stats["had_alternatives"] += 1
        try:
            best, fix = get_best_transcript_version(base, ver, remaining, VersionStrategy.UP_THEN_DOWN)
            if best in remaining:
                stats["rescued"] += 1
        except Exception:
            pass
    return stats


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
        print("\n== Recovery: cleaning (inject corruption -> clean_hgvs -> resolves?) ==")
        cleaning = run_cleaning(am, hp, pairs)
        for cat, d in sorted(cleaning.items()):
            print(f"  {cat:<22} n={d['n']:<4} parsed_broken={pct(d['resolved_raw'], d['n']):>6} "
                  f"cleaned_to_orig={pct(d['cleaned_ok'], d['n']):>6} "
                  f"parsed_after_clean={pct(d['resolved_clean'], d['n']):>6}")
        results["cleaning"] = cleaning

        print("\n== Recovery: version bumping (drop requested version -> find neighbour?) ==")
        vb = run_version_bump(provider, pairs)
        print(f"  transcripts tested      : {vb['n']}")
        print(f"  had alternative versions: {vb['had_alternatives']}  ({pct(vb['had_alternatives'], vb['n'])})")
        print(f"  rescued by bump         : {vb['rescued']}  ({pct(vb['rescued'], vb['n'])})")
        print(f"  no alternatives in data : {vb['no_alternatives']}")
        results["version_bump"] = vb

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(results, indent=2))
        print(f"\nWritten: {args.out}")


if __name__ == "__main__":
    main()
