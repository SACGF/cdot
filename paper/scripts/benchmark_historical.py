#!/usr/bin/env python3
"""Historical clinical-lab resolution benchmark: cdot vs UTA.

Motivation
----------
cdot's collection of *historical* transcript versions was driven by real clinical
data: variant-curation labs accumulate HGVS over many years, citing whatever
transcript version was current when each variant was first classified. A single
annotation snapshot (one UTA release, one RefSeq release) does not necessarily
carry every version a lab's back-catalogue references.

This benchmark resolves a real corpus of such historical clinical-lab c.HGVS
descriptions (the alleles imported into the Shariant variant-sharing platform)
through the *identical* biocommons/hgvs engine, swapping only the transcript-data
layer (cdot local JSON vs a locally loaded UTA). It reports, per provider, the
fraction that produce a genomic coordinate (pure coordinate projection,
``replace_reference=False`` — so the sequence layer never differs and only the
transcript data matters), then cross-tabulates the two to isolate the strings one
backend resolves and the other does not, and classifies *why* the loser fails
(Ensembl transcript absent / RefSeq historical version not aligned / accession
absent entirely).

Unlike the ClinVar benchmark there is no ground-truth genomic coordinate for this
private corpus, so the metric is *resolution rate* (did the backend produce a
coordinate at all), not correctness. The corpus is not shareable, so the resulting
numbers are reported in the paper as frozen Tier-2 constants — this script is the
provenance, not a referee-runnable reproduction.

Usage
-----
    UTA_DB_URL=postgresql://uta:uta@localhost:5432/uta/uta_20241220 \
    HGVS_SEQREPO_DIR=/path/to/seqrepo \
    python paper/scripts/benchmark_historical.py CORPUS.txt \
        --json /path/cdot.refseq.grch38.json.gz \
        --uta-uri postgresql://uta:uta@localhost:5432/uta/uta_20241220 \
        --out output/facts/historical.json

Input: one HGVS string per line (c. or n.), optionally with a gene symbol in
parentheses, e.g. ``NM_000016.5(ACADM):c.127G>A`` — the biocommons parser keys on
the accession and ignores the parenthetical, so no cleaning is required.
"""
import argparse
import json
import logging
import random
import re
import time
from pathlib import Path

import hgvs.parser
from hgvs.assemblymapper import AssemblyMapper
from hgvs.exceptions import HGVSDataNotAvailableError, HGVSInvalidVariantError

from cdot.hgvs.dataproviders import JSONDataProvider, FastaSeqFetcher

# transcript accession at the start of a string, captured as (base, version)
_AC_RE = re.compile(r"^([A-Za-z]+_?\d+)\.(\d+)")
_REFSEQ_PREFIXES = ("NM_", "NR_", "XM_", "XR_")


def source_of(s):
    """RefSeq vs Ensembl vs other, from the leading transcript accession."""
    m = _AC_RE.match(s)
    if not m:
        return "other"
    base = m.group(1)
    if base.startswith("ENST"):
        return "ensembl"
    if base.startswith(_REFSEQ_PREFIXES):
        return "refseq"
    return "other"


def per_source_resolution(buckets, strings):
    """{source: {n, resolved, resolved_pct}} for one provider's per-string buckets."""
    res = {src: {"n": 0, "resolved": 0} for src in ("refseq", "ensembl", "other")}
    for s in strings:
        d = res[source_of(s)]
        d["n"] += 1
        if buckets[s] == "resolved":
            d["resolved"] += 1
    for d in res.values():
        d["resolved_pct"] = pct(d["resolved"], d["n"])
    return res


def load_strings(path, sample=None, seed=1):
    raw = [l.strip() for l in open(path) if l.strip()]
    if sample and sample < len(raw):
        random.seed(seed)
        raw = random.sample(raw, sample)
    return raw


def classify(am, hp, s):
    """Return one of: resolved / no_data / error / parsefail."""
    try:
        v = hp.parse_hgvs_variant(s)
    except Exception:  # noqa: BLE001 - corpus contains non-HGVS junk lines
        return "parsefail"
    try:
        if v.type == "c":
            am.c_to_g(v)
        elif v.type == "n":
            am.n_to_g(v)
        else:
            return "parsefail"
    except HGVSDataNotAvailableError:
        return "no_data"
    except HGVSInvalidVariantError:
        return "error"
    except Exception as e:  # noqa: BLE001 - benchmark keeps going
        logging.debug("error on %s: %s", s, e)
        return "error"
    return "resolved"


def run_provider(provider, strings, build, label):
    hp = hgvs.parser.Parser()
    am = AssemblyMapper(provider, assembly_name=build, alt_aln_method="splign",
                        replace_reference=False)
    buckets = {}
    t0 = time.perf_counter()
    for i, s in enumerate(strings):
        buckets[s] = classify(am, hp, s)
        if (i + 1) % 2000 == 0:
            logging.info("%s: %d/%d (%.0fs)", label, i + 1, len(strings),
                         time.perf_counter() - t0)
    wall = time.perf_counter() - t0
    counts = {k: 0 for k in ("resolved", "no_data", "error", "parsefail")}
    for v in buckets.values():
        counts[v] += 1
    return buckets, counts, wall


def uta_aligned_accessions(uri, build_contig_prefix="NC_0000"):
    """Set of transcript accessions UTA has a splign alignment for (any NC_ contig).
    Used only to explain *why* UTA misses a string, never for the headline rate."""
    import psycopg2
    # uri ends with /<db>/<schema>; psycopg2 wants /<db> and a separate search_path
    m = re.match(r"(.*/[^/]+)/([^/]+)$", uri)
    dsn, schema = (m.group(1), m.group(2)) if m else (uri, "uta_20241220")
    conn = psycopg2.connect(dsn)
    cur = conn.cursor()
    cur.execute(f"set search_path to {schema}")
    cur.execute("select distinct tx_ac from tx_exon_aln_v "
                "where alt_aln_method='splign' and alt_ac like %s",
                (build_contig_prefix + "%",))
    acc = {r[0] for r in cur.fetchall()}
    conn.close()
    return acc


def pct(n, d):
    return round(100 * n / d, 1) if d else 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("corpus", help="one HGVS string per line")
    ap.add_argument("--json", nargs="+", help="cdot JSON.gz file(s) for the target build "
                    "(e.g. RefSeq + Ensembl); loaded together as one provider")
    ap.add_argument("--uta-uri", help="UTA postgres URI .../<db>/<schema>")
    ap.add_argument("--build", default="GRCh38")
    ap.add_argument("--fasta", nargs="+", help="genome FASTA(s) for cdot (only needed "
                    "with --replace-reference; default projection needs none)")
    ap.add_argument("--sample", type=int, help="resolve a random N-line subset (testing)")
    ap.add_argument("--out", help="write JSON facts here")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.debug else logging.INFO,
                        format="%(message)s")
    # hgvs logs an INFO/WARNING per intronic n. variant ("Can't update reference
    # sequence ... without alt_ac"); harmless here (we never replace_reference), so mute it.
    logging.getLogger("hgvs").setLevel(logging.ERROR)

    strings = load_strings(args.corpus, args.sample)
    print(f"== Historical corpus: {len(strings)} lines  |  {args.build}  |  "
          f"pure projection (replace_reference=False) ==")

    # corpus source composition (RefSeq vs Ensembl): the mix of real clinical-lab
    # traffic, so the per-source resolution rates below can be read in context.
    comp = {src: 0 for src in ("refseq", "ensembl", "other")}
    for s in strings:
        comp[source_of(s)] += 1
    n = len(strings)
    print(f"\nCorpus source mix: "
          f"refseq {comp['refseq']} ({pct(comp['refseq'], n)}%), "
          f"ensembl {comp['ensembl']} ({pct(comp['ensembl'], n)}%), "
          f"other {comp['other']} ({pct(comp['other'], n)}%)")

    results = {"n_lines": len(strings), "build": args.build,
               "source_mix": {src: {"n": c, "pct": pct(c, n)} for src, c in comp.items()}}
    cdot_b = uta_b = None

    if args.json:
        sf = FastaSeqFetcher(*args.fasta) if args.fasta else None
        prov = JSONDataProvider(args.json, seqfetcher=sf)
        cdot_b, c, w = run_provider(prov, strings, args.build, "cdot")
        print(f"\ncdot ({', '.join(Path(p).name for p in args.json)}):")
        for k in ("resolved", "no_data", "error", "parsefail"):
            print(f"  {k:<10}: {c[k]:>6}  ({pct(c[k], len(strings))}%)")
        print(f"  throughput: {len(strings)/w:.1f} HGVS/s")
        cdot_by_source = per_source_resolution(cdot_b, strings)
        for src in ("refseq", "ensembl"):
            d = cdot_by_source[src]
            print(f"  {src:<8} resolved {d['resolved']}/{d['n']} ({d['resolved_pct']}%)")
        results["cdot"] = {"counts": c, "resolved_pct": pct(c["resolved"], len(strings)),
                           "by_source": cdot_by_source,
                           "throughput_hgvs_per_s": round(len(strings)/w, 1)}

    if args.uta_uri:
        import hgvs.dataproviders.uta
        prov = hgvs.dataproviders.uta.connect(args.uta_uri)
        uta_b, c, w = run_provider(prov, strings, args.build, "uta")
        print(f"\nUTA ({args.uta_uri.rsplit('/', 1)[-1]}):")
        for k in ("resolved", "no_data", "error", "parsefail"):
            print(f"  {k:<10}: {c[k]:>6}  ({pct(c[k], len(strings))}%)")
        print(f"  throughput: {len(strings)/w:.1f} HGVS/s")
        uta_by_source = per_source_resolution(uta_b, strings)
        for src in ("refseq", "ensembl"):
            d = uta_by_source[src]
            print(f"  {src:<8} resolved {d['resolved']}/{d['n']} ({d['resolved_pct']}%)")
        results["uta"] = {"counts": c, "resolved_pct": pct(c["resolved"], len(strings)),
                          "by_source": uta_by_source,
                          "throughput_hgvs_per_s": round(len(strings)/w, 1)}

    # ---- head-to-head: strings cdot resolves but UTA does not ----------------
    if cdot_b and uta_b:
        cdot_only = [s for s in strings
                     if cdot_b[s] == "resolved" and uta_b[s] != "resolved"]
        uta_only = [s for s in strings
                    if uta_b[s] == "resolved" and cdot_b[s] != "resolved"]
        both = sum(1 for s in strings if cdot_b[s] == "resolved" and uta_b[s] == "resolved")
        print(f"\n== Head-to-head ==")
        print(f"  both resolve     : {both}  ({pct(both, len(strings))}%)")
        print(f"  cdot only        : {len(cdot_only)}  ({pct(len(cdot_only), len(strings))}%)")
        print(f"  UTA only         : {len(uta_only)}  ({pct(len(uta_only), len(strings))}%)")

        # why does UTA miss the cdot-only strings? classify by transcript.
        cat = {"ensembl": 0, "refseq_historical": 0, "refseq_absent": 0, "other": 0}
        aligned = uta_aligned_accessions(args.uta_uri) if args.uta_uri else set()
        aligned_base = {a.rsplit(".", 1)[0] for a in aligned}
        for s in cdot_only:
            m = _AC_RE.match(s)
            if not m:
                cat["other"] += 1
                continue
            base, _ = m.group(1), m.group(2)
            ac = f"{m.group(1)}.{m.group(2)}"
            if base.startswith("ENST"):
                cat["ensembl"] += 1
            elif base in aligned_base:
                cat["refseq_historical"] += 1   # base aligned in UTA, this version not
            else:
                cat["refseq_absent"] += 1
        print("  cdot-only failures in UTA, by cause:")
        for k, n in cat.items():
            print(f"    {k:<18}: {n}  ({pct(n, len(cdot_only))}% of cdot-only)")
        results["head_to_head"] = {
            "both_resolved": both, "cdot_only": len(cdot_only), "uta_only": len(uta_only),
            "cdot_only_by_cause": cat,
        }

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(results, indent=2))
        print(f"\nWritten: {args.out}")


if __name__ == "__main__":
    main()
