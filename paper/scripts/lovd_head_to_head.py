#!/usr/bin/env python3
"""
Head-to-head: cdot ``clean_hgvs()`` vs the LOVD HGVS syntax checker (issue #112).

Runs the exact injection corpus from ``inject_and_clean.py`` (same seed, same
per-class caps, same injectors, so the case list is identical) through both
tools and scores them with the same rule.

LOVD HGVS syntax checker
------------------------
github.com/LOVDnl/HGVS-syntax-checker: a deliberately sequence-agnostic PHP
library that validates HGVS syntax and suggests corrections for invalid
descriptions, ranked by a confidence score. Run locally from the CLI::

    git clone https://github.com/LOVDnl/HGVS-syntax-checker
    (cd HGVS-syntax-checker && git checkout v1.2.2)
    php -f HGVS-syntax-checker/HGVS.php "NM_004006.3:c.157C>T"   # emits JSON

Version pinned for the paper: **v1.2.2** (commit 4cd074ba9cbd, 2026-06-15),
PHP 8.5 CLI, fully offline. Pass the path to ``HGVS.php`` via ``--checker``.

Scoring (decided up front, applied identically to both tools)
-------------------------------------------------------------
Each injected case has a known canonical target (``expected``). A tool
recovers a case iff its output string equals the target, compared
*gene-annotation-insensitively*: a parenthesised gene symbol between the
accession and the colon is stripped from both sides before comparing, because
the two ecosystems canonicalise that form in opposite directions (biocommons
keeps ``NM_x.y(GENE):c.``, LOVD corrects it to ``NM_x.y:c.``) and both denote
the same variant. Concretely:

- **cdot**: ``clean_hgvs(perturbed)`` == target.
- **LOVD top-1** (primary): the checker's highest-confidence entry in
  ``corrected_values`` == target.
- **LOVD any-rank** (secondary): target appears anywhere in the ranked
  ``corrected_values`` list.

False-correction check on the uncorrupted originals: a tool "falsely corrects"
a valid input if its output differs from the input (same insensitive
comparison). For LOVD this also counts inputs it flags invalid.

Stratification note: the fairness plan called for splitting string-repairable
vs data-requiring categories, but the injection corpus contains *only*
string-repairable categories by construction (``inject_and_clean.py`` already
excludes ops that need a data provider, eg accession-prefix restoration), so
both tools are eligible on every case and no split is needed. Per-category
results are still emitted so scope disputes stay visible.

Usage::

    python paper/scripts/lovd_head_to_head.py --checker /path/to/HGVS.php
"""

import argparse
import csv
import json
import random
import re
import subprocess
import sys
import time
from pathlib import Path

import hgvs.parser

from cdot.hgvs.clean import clean_hgvs

sys.path.insert(0, str(Path(__file__).resolve().parent))
import inject_and_clean as iac  # noqa: E402  (injectors, seed, sample, weights)

FACTS_DIR = iac.FACTS_DIR
LOVD_VERSION = "v1.2.2"

# Strip a parenthesised gene symbol between accession and colon:
# "NM_000059.4(BRCA2):c.316G>A" -> "NM_000059.4:c.316G>A"
_GENE_ANNOT = re.compile(r"^([A-Za-z]{2,4}_?\d+\.\d+)\(([^()]+)\):")


def norm(s):
    return _GENE_ANNOT.sub(r"\1:", s)


# ---------------------------------------------------------------------------
# Case generation: mirrors inject_and_clean.run() exactly (same seeds, same
# skip rules, same caps) so the attempted case list is identical.
# ---------------------------------------------------------------------------

def build_cases(parser, pool):
    cases = []  # (category, orig, perturbed, expected)
    for idx, (name, injector) in enumerate(iac.INJECTORS):
        order = list(pool)
        rng_c = random.Random(iac.SEED + idx)
        rng_c.shuffle(order)
        attempted = 0
        gi = 0
        for orig in order:
            res = injector(orig, iac.GENE_POOL[gi % len(iac.GENE_POOL)])
            if res is None:
                continue
            if attempted >= iac.TARGET_PER_CLASS:
                continue
            gi += 1
            perturbed, expected = res
            if not iac.parses(parser, expected):
                continue
            if perturbed == expected:
                continue
            attempted += 1
            cases.append((name, orig, perturbed, expected))
    return cases


# ---------------------------------------------------------------------------
# LOVD CLI driver
# ---------------------------------------------------------------------------

def run_lovd(checker, strings, chunk=500):
    """Run the PHP CLI over strings (batched argv), return {input: result}."""
    results = {}
    todo = [s for s in dict.fromkeys(strings)]  # dedupe, keep order
    for i in range(0, len(todo), chunk):
        batch = todo[i:i + chunk]
        proc = subprocess.run(
            ["php", "-f", str(checker), *batch],
            capture_output=True, text=True, check=True)
        for item in json.loads(proc.stdout):
            results[item["input"]] = item
    missing = [s for s in todo if s not in results]
    if missing:
        raise RuntimeError(f"LOVD checker returned no result for {len(missing)} "
                           f"inputs, eg {missing[0]!r}")
    return results


def lovd_ranked(result):
    """Ranked correction list (best first) from a checker result."""
    cv = result.get("corrected_values") or {}
    if not isinstance(cv, dict):  # empty PHP array serialises as []
        return []
    return sorted(cv, key=lambda k: -cv[k])


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def score(cases, lovd_results):
    per_class = {}
    for name, _ in iac.INJECTORS:
        per_class[name] = {
            "weight": iac.REAL_RESCUE_OP_COUNTS[name],
            "n_attempted": 0, "cdot_recovered": 0,
            "lovd_top1": 0, "lovd_anyrank": 0,
        }
    for name, orig, perturbed, expected in cases:
        c = per_class[name]
        c["n_attempted"] += 1
        target = norm(expected)
        cleaned, _ = clean_hgvs(perturbed)
        if norm(cleaned) == target:
            c["cdot_recovered"] += 1
        ranked = [norm(v) for v in lovd_ranked(lovd_results[perturbed])]
        if ranked and ranked[0] == target:
            c["lovd_top1"] += 1
        if target in ranked:
            c["lovd_anyrank"] += 1

    for c in per_class.values():
        n = c["n_attempted"]
        for k in ("cdot_recovered", "lovd_top1", "lovd_anyrank"):
            c[k + "_pct"] = round(100.0 * c[k] / n, 1) if n else None
    return per_class


def weighted(per_class, key):
    wsum = wnum = 0.0
    for c in per_class.values():
        if c["n_attempted"]:
            wsum += c["weight"]
            wnum += c["weight"] * c[key + "_pct"]
    return round(wnum / wsum, 1) if wsum else None


def false_corrections(pool, lovd_results):
    lovd_fc = cdot_fc = lovd_invalid = 0
    examples = []
    for orig in pool:
        target = norm(orig)
        cleaned, _ = clean_hgvs(orig)
        if norm(cleaned) != target:
            cdot_fc += 1
        r = lovd_results[orig]
        ranked = [norm(v) for v in lovd_ranked(r)]
        changed = not ranked or ranked[0] != target
        if not r.get("valid"):
            lovd_invalid += 1
        if changed:
            lovd_fc += 1
            if len(examples) < 10:
                examples.append({"input": orig, "valid": r.get("valid"),
                                 "top": ranked[0] if ranked else None})
    return lovd_fc, lovd_invalid, cdot_fc, examples


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--checker", required=True, type=Path,
                    help="path to the LOVD checker's HGVS.php (pinned v1.2.2)")
    ap.add_argument("--limit", type=int, default=None,
                    help="cap cases per category (smoke test only; not for facts)")
    args = ap.parse_args()
    if not args.checker.exists():
        sys.exit(f"checker not found: {args.checker}")

    parser = hgvs.parser.Parser()
    pool = iac.load_sample(parser, regenerate=False)
    cases = build_cases(parser, pool)
    if args.limit:
        by_cat = {}
        cases = [c for c in cases
                 if by_cat.setdefault(c[0], []).append(c) or len(by_cat[c[0]]) <= args.limit]

    t0 = time.time()
    lovd_results = run_lovd(args.checker, [c[2] for c in cases] + list(pool))
    lovd_seconds = round(time.time() - t0, 1)

    per_class = score(cases, lovd_results)
    n = sum(c["n_attempted"] for c in per_class.values())
    totals = {k: sum(c[k] for c in per_class.values())
              for k in ("cdot_recovered", "lovd_top1", "lovd_anyrank")}
    lovd_fc, lovd_invalid, cdot_fc, fc_examples = false_corrections(pool, lovd_results)

    facts = {
        "issue": "SACGF/cdot#112",
        "tier": 1,
        "description": (
            "Head-to-head on the reproducible injection corpus: cdot clean_hgvs() "
            "vs the LOVD HGVS syntax checker (local PHP CLI, pinned "
            f"{LOVD_VERSION}). Same cases, same gene-annotation-insensitive "
            "exact-match scoring; see paper/scripts/lovd_head_to_head.py."
        ),
        "lovd_version": LOVD_VERSION,
        "n_cases": n,
        "sample_size": len(pool),
        "seed": iac.SEED,
        "cdot_recovered": totals["cdot_recovered"],
        "cdot_pct": round(100.0 * totals["cdot_recovered"] / n, 1),
        "lovd_top1": totals["lovd_top1"],
        "lovd_top1_pct": round(100.0 * totals["lovd_top1"] / n, 1),
        "lovd_anyrank": totals["lovd_anyrank"],
        "lovd_anyrank_pct": round(100.0 * totals["lovd_anyrank"] / n, 1),
        "cdot_weighted_pct": weighted(per_class, "cdot_recovered"),
        "lovd_top1_weighted_pct": weighted(per_class, "lovd_top1"),
        "lovd_anyrank_weighted_pct": weighted(per_class, "lovd_anyrank"),
        "originals_n": len(pool),
        "lovd_false_corrections": lovd_fc,
        "lovd_flagged_invalid": lovd_invalid,
        "lovd_flagged_invalid_pct": round(100.0 * lovd_invalid / len(pool), 1),
        "cdot_false_corrections": cdot_fc,
        "lovd_seconds": lovd_seconds,
        "per_class": per_class,
        "false_correction_examples": fc_examples,
    }

    FACTS_DIR.mkdir(parents=True, exist_ok=True)
    json_path = FACTS_DIR / "lovd_comparison.json"
    json_path.write_text(json.dumps(facts, indent=2) + "\n")
    csv_path = FACTS_DIR / "lovd_comparison.csv"
    row = {k: facts[k] for k in (
        "lovd_version", "n_cases", "cdot_pct", "lovd_top1_pct",
        "lovd_anyrank_pct", "cdot_weighted_pct", "lovd_top1_weighted_pct",
        "lovd_anyrank_weighted_pct", "originals_n", "lovd_false_corrections",
        "lovd_flagged_invalid", "lovd_flagged_invalid_pct",
        "cdot_false_corrections")}
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row))
        w.writeheader()
        w.writerow(row)

    print(f"\nLOVD head-to-head ({LOVD_VERSION}, {n} injected cases, "
          f"LOVD wall time {lovd_seconds}s)")
    print(f"  {'category':<34} {'n':>4} {'cdot%':>7} {'lovd@1%':>8} {'lovd@any%':>9}")
    for name, c in per_class.items():
        if not c["n_attempted"]:
            continue
        print(f"  {name:<34} {c['n_attempted']:>4} {c['cdot_recovered_pct']:>7} "
              f"{c['lovd_top1_pct']:>8} {c['lovd_anyrank_pct']:>9}")
    print(f"  {'-'*66}")
    print(f"  overall : cdot {facts['cdot_pct']}%  lovd top-1 {facts['lovd_top1_pct']}%  "
          f"lovd any-rank {facts['lovd_anyrank_pct']}%")
    print(f"  weighted: cdot {facts['cdot_weighted_pct']}%  "
          f"lovd top-1 {facts['lovd_top1_weighted_pct']}%  "
          f"lovd any-rank {facts['lovd_anyrank_weighted_pct']}%")
    print(f"  originals ({len(pool)}): lovd false-corrections {lovd_fc} "
          f"(flagged invalid {lovd_invalid}), cdot false-corrections {cdot_fc}")
    if fc_examples:
        print("  lovd false-correction examples:")
        for e in fc_examples:
            print(f"    - {e['input']!r} -> {e['top']!r} (valid={e['valid']})")
    print(f"\nWrote: {json_path}")
    print(f"Wrote: {csv_path}")


if __name__ == "__main__":
    main()
