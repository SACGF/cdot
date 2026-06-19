#!/usr/bin/env python3
"""Per-source / failure-reason breakdown of the ClinVar pass (paper Supplementary S4).

Consumes the per-variant table from resolve_clinvar_pass.py (paper infra A1) and
aggregates it into the Supplementary Table S4 breakdown: for each transcript source
(RefSeq vs Ensembl) the resolution outcome (resolved / matched / no_data / error),
and the failure reasons behind the unresolved variants. This is the cdot-only view
at full ClinVar scale: unlike the R1 cdot-vs-UTA comparison (gated by UTA throughput
to a small sample), every resolved ClinVar variant can be summarised here because
only cdot is in the loop.

It is a pure aggregation of the A1 table (no resolution), so it is instant and needs
no provider. The one optional enrichment, splitting `no_data` into unknown-accession
vs unknown-version, does need a provider (it asks which versions exist for each
transcript); it is off by default.

Caveat for the paper: ClinVar submissions are dominated by a handful of large
(largely US) clinical laboratories, and the corpus cites mostly current RefSeq
versions. So the source mix here reflects that submitter population, not clinical
practice at large; it is a clean, public, reproducible sanity check that the
pipeline resolves real variants at scale, NOT an unbiased sample of the transcripts
clinical labs actually use. The unbiased real-world complement is the Shariant
historical corpus (Results R1, Tier 2), whose older-version, multi-lab traffic is
exactly what ClinVar's recency and submitter skew leave out. S4 should be read with
that caveat stated.

Per CLAUDE.md: never run against full datasets in development. Verify on a table
built from tests/test_data/clinvar_hgvs/; the real S4 is a dedicated full run.

Usage:
    # pure aggregation (no provider, instant)
    python paper/scripts/summarize_clinvar_pass.py output/clinvar_pass/refseq.csv \
        --out output/facts/clinvar_source_breakdown.csv

    # also split no_data into unknown-accession vs unknown-version
    python paper/scripts/summarize_clinvar_pass.py output/clinvar_pass/refseq.csv \
        --split-no-data --json cdot.refseq.grch38.json.gz \
        --out output/facts/clinvar_source_breakdown.csv
"""
import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

_REFSEQ_PREFIXES = ("NM_", "NR_", "XM_", "XR_")


def source_of(tx):
    """Map a transcript accession base to its annotation source."""
    if tx.startswith(_REFSEQ_PREFIXES):
        return "refseq"
    if tx.startswith("ENST"):
        return "ensembl"
    return "other"


def pct(n, d):
    return round(100 * n / d, 2) if d else 0.0


def split_no_data(df, provider, build):
    """Return (unknown_accession, unknown_version) counts over the no_data rows.

    A no_data variant is unknown-accession if the provider holds NO version of the
    transcript, else unknown-version (the accession exists, the cited version does
    not). Versions are cached per accession."""
    cache = {}

    def versions(tx):
        if tx not in cache:
            try:
                cache[tx] = set(provider.get_tx_versions(tx, build))
            except Exception as e:  # noqa: BLE001
                logging.debug("get_tx_versions(%s): %s", tx, e)
                cache[tx] = set()
        return cache[tx]

    unknown_acc = unknown_ver = 0
    for _, row in df[df["bucket"] == "no_data"].iterrows():
        if not versions(row["tx"]):
            unknown_acc += 1
        else:
            unknown_ver += 1
    return unknown_acc, unknown_ver


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results", help="A1 per-variant table (resolve_clinvar_pass.py output)")
    ap.add_argument("--out", default="output/facts/clinvar_source_breakdown.csv")
    ap.add_argument("--split-no-data", action="store_true",
                    help="split no_data into unknown-accession vs unknown-version (needs a provider)")
    # provider flags (only used with --split-no-data)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--rest", action="store_true", help="cdot REST (default for --split-no-data)")
    g.add_argument("--rest-insecure", action="store_true")
    g.add_argument("--json", help="cdot JSON.gz file")
    g.add_argument("--uta", action="store_true")
    ap.add_argument("--build", default="GRCh38")
    ap.add_argument("--fasta", nargs="+")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.debug else logging.WARNING)

    df = pd.read_csv(args.results, dtype=str).fillna("")
    df["source"] = df["tx"].map(source_of)

    rows = []
    for source in ["refseq", "ensembl", "other", "total"]:
        sub = df if source == "total" else df[df["source"] == source]
        n = len(sub)
        if n == 0 and source != "total":
            continue
        bc = sub["bucket"].value_counts().to_dict()
        correct, incorrect = bc.get("correct", 0), bc.get("incorrect", 0)
        no_data, error = bc.get("no_data", 0), bc.get("error", 0)
        resolved = correct + incorrect
        rows.append({
            "source": source, "n": n,
            "resolved": resolved, "resolved_pct": pct(resolved, n),
            "matched": correct, "matched_pct": pct(correct, n),
            "incorrect": incorrect, "no_data": no_data, "error": error,
        })

    breakdown = pd.DataFrame(rows)

    enrich = None
    if args.split_no_data:
        args.replace_reference = False
        from benchmark_resolution import build_provider  # noqa: E402
        provider, label = build_provider(args)
        if not hasattr(provider, "get_tx_versions"):
            sys.exit(f"provider {label} cannot enumerate versions; drop --split-no-data")
        unknown_acc, unknown_ver = split_no_data(df, provider, args.build)
        enrich = {"unknown_accession": unknown_acc, "unknown_version": unknown_ver,
                  "provider": label}

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    breakdown.to_csv(out, index=False)

    print("== ClinVar pass breakdown by source (Supplementary S4) ==")
    print("   NOTE: ClinVar is skewed to large (largely US) labs and current RefSeq")
    print("   versions; this is a scale sanity check, not an unbiased transcript sample.")
    print(breakdown.to_string(index=False))
    if enrich:
        n_no_data = int((df["bucket"] == "no_data").sum())
        print(f"\n  no_data failure reasons ({enrich['provider']}):")
        print(f"    unknown_accession : {enrich['unknown_accession']:>6}  "
              f"({pct(enrich['unknown_accession'], n_no_data)}% of no_data)")
        print(f"    unknown_version   : {enrich['unknown_version']:>6}  "
              f"({pct(enrich['unknown_version'], n_no_data)}% of no_data)")
        # append a second facts file so the table stays tidy
        enr_out = out.with_name(out.stem + "_failure_reasons.csv")
        pd.DataFrame([{
            "n_no_data": n_no_data,
            "unknown_accession": enrich["unknown_accession"],
            "unknown_version": enrich["unknown_version"],
        }]).to_csv(enr_out, index=False)
        print(f"  Written: {enr_out}")
    print(f"  Written: {out}")


if __name__ == "__main__":
    main()
