#!/usr/bin/env python3
"""Single ClinVar resolution pass -> one per-variant results table (paper infra A1).

Several paper analyses each walk all of ClinVar and re-resolve the same (g.HGVS,
c.HGVS) pairs (benchmark_resolution.py, the planned R5b safe-version validation,
the R6 genomic-mismatch categorisation, the clinvar.csv aggregate). Resolving is
the slow part (one c_to_g per variant), so doing it once and writing the result
table that every downstream consumes is both faster and guarantees the sections
cannot drift out of sync.

This script does exactly one resolution pass and emits a per-variant table with:

    g_hgvs        ground-truth genomic HGVS (ClinVar CLNHGVS)
    c_hgvs        input transcript HGVS (as given)
    tx            transcript accession base (e.g. NM_000059), version stripped
    version       transcript version as an int (e.g. 4), or empty if versionless
    bucket        correct | incorrect | no_data | error
    converted_g   cdot's c_to_g output, or empty if it could not convert
    fix_codes     ';'-joined HGVSFix codes fix_hgvs() would apply (empty on clean
                  input); only populated with --with-fixes (off by default)

bucket is the BASELINE resolution of the c.HGVS exactly as given (no cleaning, no
version bump), scored against g_hgvs. That is what R5b (does a safe version bump
preserve converted_g?) and R6 (re-normalise the incorrect bucket) both build on.
fix_codes is recorded alongside so a recovery analysis can attribute rescues
without a second pass.

Provider is pluggable (shared with benchmark_resolution.py): defaults to cdot REST
so it runs with no local data; pass --json a cdot release for the offline run.

Per CLAUDE.md: never run against full datasets in development. Verify with the
committed pairs in tests/test_data/clinvar_hgvs/; the full pass is a dedicated run.

Usage:
    # quick: REST over the committed test pairs
    python paper/scripts/resolve_clinvar_pass.py \
        tests/test_data/clinvar_hgvs/clinvar_hgvs_100.tsv --rest \
        --out output/clinvar_pass/refseq_100.csv

    # offline over a downloaded release (fast, no network)
    python paper/scripts/resolve_clinvar_pass.py pairs.tsv \
        --json cdot.refseq.grch38.json.gz --fasta GRCh38.fa \
        --out output/clinvar_pass/refseq.csv

Input file: TSV with two columns  g.HGVS <TAB> c.HGVS  (one variant per line).
"""
import argparse
import logging
import re
import sys
import time
from pathlib import Path

import pandas as pd
import hgvs.parser
from hgvs.assemblymapper import AssemblyMapper

# Reuse the provider construction, pair loader and classifier from the sibling
# benchmark script so the two stay in lock-step (same buckets, same provider flags).
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_resolution import build_provider, load_pairs, classify  # noqa: E402

from cdot.hgvs.clean import HGVSFixSeverity, VersionStrategy  # noqa: E402
from cdot.hgvs.gene_hgvs import fix_hgvs  # noqa: E402


# Accession base + version from "NM_000059.4:c.123A>T" (version optional).
_TX_VERSION = re.compile(r"^([A-Za-z]+_?\d+)(?:\.(\d+))?:")


def parse_tx_version(c_hgvs):
    """Return (tx_base, version_int_or_None) for a transcript HGVS string."""
    m = _TX_VERSION.match(c_hgvs)
    if not m:
        return "", None
    return m.group(1), (int(m.group(2)) if m.group(2) else None)


def fix_codes_for(provider, build, c_hgvs, version_fallback):
    """Return the ';'-joined WARNING-level HGVSFix codes fix_hgvs() applies to
    c_hgvs (empty string when nothing fires, e.g. already-clean input)."""
    try:
        _fixed, fixes = fix_hgvs(c_hgvs, provider, build, version_fallback=version_fallback)
    except Exception as e:  # noqa: BLE001 - keep the pass going
        logging.debug("fix_hgvs error on %s: %s", c_hgvs, e)
        return ""
    return ";".join(f.code.name for f in fixes if f.severity == HGVSFixSeverity.WARNING)


def run_pass(am, hp, provider, build, pairs, with_fixes, version_fallback):
    """Resolve every pair once; return (rows, wall_seconds)."""
    rows = []
    t0 = time.perf_counter()
    for g_hgvs, c_hgvs in pairs:
        bucket, converted = classify(am, hp, g_hgvs, c_hgvs)
        tx, version = parse_tx_version(c_hgvs)
        rows.append({
            "g_hgvs": g_hgvs,
            "c_hgvs": c_hgvs,
            "tx": tx,
            "version": version if version is not None else "",
            "bucket": bucket,
            "converted_g": converted if converted is not None else "",
            "fix_codes": fix_codes_for(provider, build, c_hgvs, version_fallback) if with_fixes else "",
        })
    return rows, time.perf_counter() - t0


def check_out_format(out_path):
    """Fail fast (before the slow pass) if a parquet target has no engine installed."""
    if Path(out_path).suffix == ".parquet":
        try:
            import pyarrow  # noqa: F401
        except ImportError:
            try:
                import fastparquet  # noqa: F401
            except ImportError:
                sys.exit("--out is .parquet but neither pyarrow nor fastparquet is "
                         "installed; use a .csv path or `pip install pyarrow`.")


def write_table(rows, out_path):
    """Write the per-variant table as CSV or parquet (by extension)."""
    df = pd.DataFrame(rows, columns=[
        "g_hgvs", "c_hgvs", "tx", "version", "bucket", "converted_g", "fix_codes"])
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.suffix == ".parquet":
        df.to_parquet(out, index=False)
    else:
        df.to_csv(out, index=False)
    return df


def pct(n, d):
    return f"{100 * n / d:.1f}%" if d else "n/a"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pairs_file", help="TSV: g.HGVS <TAB> c.HGVS")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--rest", action="store_true", help="cdot REST (default)")
    g.add_argument("--rest-insecure", action="store_true")
    g.add_argument("--json", help="cdot JSON.gz file")
    g.add_argument("--uta", action="store_true")
    ap.add_argument("--build", default="GRCh38")
    ap.add_argument("--fasta", nargs="+", help="genome FASTA(s) for local sequence fetching")
    ap.add_argument("--replace-reference", action="store_true",
                    help="validate ref base (needs a seqfetcher); off = pure coordinate projection")
    ap.add_argument("--with-fixes", dest="with_fixes", action="store_true", default=False,
                    help="also run fix_hgvs() per variant and record its codes (off by "
                         "default: it doubles the resolution work and fires nothing on "
                         "already-clean ClinVar input)")
    ap.add_argument("--out", required=True, help="per-variant results table (.csv or .parquet)")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.debug else logging.WARNING)

    check_out_format(args.out)
    pairs = load_pairs(args.pairs_file)

    t0 = time.perf_counter()
    provider, label = build_provider(args)
    load_time = time.perf_counter() - t0

    am = AssemblyMapper(provider, assembly_name=args.build, alt_aln_method="splign",
                        replace_reference=args.replace_reference)
    hp = hgvs.parser.Parser()

    print(f"== ClinVar pass: {label}  |  {len(pairs)} pairs  |  {args.build}  "
          f"|  with_fixes={args.with_fixes} ==")
    rows, wall = run_pass(am, hp, provider, args.build, pairs, args.with_fixes,
                          VersionStrategy.UP_THEN_DOWN)
    df = write_table(rows, args.out)

    counts = df["bucket"].value_counts().to_dict()
    total = len(df)
    resolved = counts.get("correct", 0) + counts.get("incorrect", 0)
    print(f"  load_time         : {load_time:.2f}s")
    for b in ("correct", "incorrect", "no_data", "error"):
        print(f"  {b:<17} : {counts.get(b, 0):>5}  ({pct(counts.get(b, 0), total)})")
    print(f"  resolved (any g.) : {resolved:>5}  ({pct(resolved, total)})")
    print(f"  throughput        : {total / wall:.1f} HGVS/s  (wall {wall:.2f}s)")
    print(f"  Written: {args.out}  ({total} rows)")


if __name__ == "__main__":
    main()
