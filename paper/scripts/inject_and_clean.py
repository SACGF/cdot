#!/usr/bin/env python3
"""
Tier-1 reproducible cleaning benchmark (issue #112).

Takes public ClinVar c.HGVS strings (from ``tests/test_data/clinvar_hgvs/``),
injects one well-defined formatting error per ``cdot.hgvs.clean`` fix category
into each clean string, runs ``clean_hgvs()``, and records whether the perturbed
string is recovered.

Recovery metric (deterministic, parser-independent)
---------------------------------------------------
The biocommons HGVS parser is lenient (it tolerates lowercase bases, a missing
transcript underscore, etc.), so "does it parse again" is a weak signal. Instead
each injector declares the *expected canonical output* — the exact string a
correct cleaner should produce — and a case counts as **recovered** iff::

    clean_hgvs(perturbed) == expected      and      biocommons parses expected

For noise injections (whitespace, quotes, doubled punctuation, …) ``expected``
is the original clean string, so recovery == an exact round-trip. For injections
that add information (a gene symbol in the wrong slot) ``expected`` is the
canonical annotated form.

**Regression** = a baseline string that parsed before cleaning but no longer
parses (or is altered) after ``clean_hgvs()``. Expected to be 0.

Privacy / two-tier fact model (see ``claude/paper_plan.md`` §3)
--------------------------------------------------------------
This benchmark is Tier 1 — fully reproducible from public ClinVar data committed
to this repo. The per-category injection *weights* are frozen CONSTANTS derived
from the aggregate rescue-op distribution measured on a private production corpus
(``cdot_private/output/cleaning_analysis_20260617.txt``, N=32,752; issue #112).
Only the aggregate counts are cited here — no corpus string is copied into this
repo, per CLAUDE.md.

Usage::

    python paper/scripts/inject_and_clean.py                 # use committed sample
    python paper/scripts/inject_and_clean.py --regenerate     # rebuild sample list
"""

import argparse
import csv
import json
import random
import re
from pathlib import Path

import hgvs.parser

from cdot.hgvs.clean import clean_hgvs

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

REPO = Path(__file__).resolve().parent.parent.parent
SOURCE_TSVS = [
    REPO / "tests/test_data/clinvar_hgvs/clinvar_hgvs_500.tsv",          # RefSeq
    REPO / "tests/test_data/clinvar_hgvs/clinvar_hgvs_ensembl_500.tsv",  # Ensembl
]
SAMPLE_FILE = Path(__file__).resolve().parent / "clean_injection_sample.txt"
FACTS_DIR = REPO / "output/facts"

SEED = 112  # issue #112 — fixed for reproducibility
TARGET_PER_CLASS = 200  # cap per category (applicable pools may be smaller)

# Public gene symbols only (used to build gene-wrapper / gene-transcript-swap
# injections). They need not match the transcript — recovery is a string round
# trip, not a biological assertion.
GENE_POOL = ["BRCA2", "RUNX1", "BRCA1", "TP53", "MSH2"]

# ---------------------------------------------------------------------------
# Injection weights — FROZEN CONSTANTS (Tier 2, cited not shipped)
# Source: cdot_private/output/cleaning_analysis_20260617.txt
#   "Which clean ops do the rescuing?" (rescued rows, ops counted; N=32,752).
# Used only to weight the reproducible per-class recovery into a single
# headline number that reflects the real-world mix of errors. No corpus string
# is read or copied.
# ---------------------------------------------------------------------------

REAL_RESCUE_OP_COUNTS = {
    "STRIPPED_WHITESPACE": 940,
    "UPPERCASED_BASES": 553,
    "STRIPPED_PROTEIN_SUFFIX": 130,
    "FIXED_GENE_WRAPPER": 100,
    "SWAPPED_GENE_TRANSCRIPT": 99,
    "STRIPPED_SURROUNDING_PUNCTUATION": 84,
    "FIXED_MULTIPLE_COLON": 60,
    "STRIPPED_UNBALANCED_BRACKETS": 58,
    "FIXED_SEPARATOR_TYPO": 37,
    "FIXED_MULTIPLE_DOT": 35,
    "STRIPPED_LEADING_JUNK": 23,
    "FIXED_MULTIPLE_KIND": 13,
    "UPPERCASED_TRANSCRIPT": 10,
    "DROPPED_DEL_DUP_COUNT": 8,
    "ADDED_TRANSCRIPT_UNDERSCORE": 6,
    "FIXED_PREFIX_COLON": 3,
    "LOWERCASED_MUTATION_TYPE": 3,
    "ADDED_N_PREFIX": 2,
}

# ---------------------------------------------------------------------------
# Injectors: orig (clean, parseable) -> (perturbed, expected) or None.
# Each returns None when the category does not apply to that string.
# `expected` is the exact canonical string clean_hgvs() should produce.
# ---------------------------------------------------------------------------

_SUB_BASES = re.compile(r"([ACGT]+)>([ACGT]+)$")
_OP_BASES = re.compile(r"(del|ins|dup)([ACGT]+)$")
_KIND_DOT = re.compile(r":([cgnmpr])\.")
_RANGED_DEL_DUP = re.compile(r"\d+_\d[\d+\-_*]*(del|dup)$")
_TX_GENE_PREFIXES = ("NM_", "NR_", "XM_", "XR_", "ENST")


def _split(orig):
    """(prefix, allele) split on the first colon, or None if absent."""
    if ":" not in orig:
        return None
    return tuple(orig.split(":", 1))


def inj_whitespace(orig, gene):
    # Leading space + a space after the colon: "NM_x:c.1A>G" -> " NM_x: c.1A>G"
    return " " + orig.replace(":", ": ", 1), orig


def inj_uppercased_bases(orig, gene):
    m = _SUB_BASES.search(orig)
    if m:
        lowered = m.group(1).lower() + ">" + m.group(2).lower()
        return orig[: m.start()] + lowered, orig
    m = _OP_BASES.search(orig)
    if m:
        lowered = m.group(1) + m.group(2).lower()
        return orig[: m.start()] + lowered, orig
    return None


def inj_protein_suffix(orig, gene):
    # Append a synthetic trailing protein annotation (space + p.)
    return orig + " p.Arg100Ter", orig


def inj_gene_wrapper(orig, gene):
    s = _split(orig)
    if not s:
        return None
    tx, allele = s
    if not tx.upper().startswith(_TX_GENE_PREFIXES):
        return None
    # "NM_x.v:(GENE):c..."  ->  "NM_x.v(GENE):c..."
    return f"{tx}:({gene}):{allele}", f"{tx}({gene}):{allele}"


def inj_swapped_gene_transcript(orig, gene):
    s = _split(orig)
    if not s:
        return None
    tx, allele = s
    if not tx.upper().startswith(_TX_GENE_PREFIXES):
        return None
    # gene and transcript swapped: "GENE(NM_x.v):c..." -> "NM_x.v(GENE):c..."
    return f"{gene}({tx}):{allele}", f"{tx}({gene}):{allele}"


def inj_surrounding_punctuation(orig, gene):
    # Wrap in stray quotes (copy/paste artefact)
    return f'"{orig}"', orig


def inj_double_colon(orig, gene):
    return orig.replace(":", "::", 1), orig


def inj_unbalanced_brackets(orig, gene):
    if "(" in orig or ")" in orig:
        return None
    return "(" + orig, orig  # stray unbalanced open paren


def inj_separator_typo(orig, gene):
    m = _KIND_DOT.search(orig)
    if not m:
        return None
    # ":c." -> ":c," (comma in place of the kind dot)
    return orig[: m.start()] + f":{m.group(1)}," + orig[m.end():], orig


def inj_double_dot(orig, gene):
    if "." not in orig:
        return None
    i = orig.index(".")  # first dot = transcript version dot
    return orig[:i] + ".." + orig[i + 1:], orig


def inj_leading_junk(orig, gene):
    if not orig.upper().startswith(_TX_GENE_PREFIXES):
        return None
    return "GRCh38.p2 " + orig, orig


def inj_double_kind(orig, gene):
    m = _KIND_DOT.search(orig)
    if not m:
        return None
    kind = m.group(1)
    return orig[: m.start()] + f":{kind}.{kind}." + orig[m.end():], orig


def inj_uppercased_transcript(orig, gene):
    s = _split(orig)
    if not s:
        return None
    tx, allele = s
    if not tx.upper().startswith(_TX_GENE_PREFIXES):
        return None
    return f"{tx.lower()}:{allele}", orig


def inj_dropped_del_dup_count(orig, gene):
    if not _RANGED_DEL_DUP.search(orig):
        return None
    return orig + "23", orig  # redundant base count after a ranged del/dup


def inj_added_transcript_underscore(orig, gene):
    if orig.startswith("NM_"):
        return "NM" + orig[3:], orig
    if orig.startswith("NC_"):
        return "NC" + orig[3:], orig
    return None


def inj_fixed_prefix_colon(orig, gene):
    for pfx in ("NM_", "NR_", "XM_", "XR_"):
        if orig.startswith(pfx):
            return orig[:2] + ":_" + orig[3:], orig
    return None


def inj_lowercased_mutation_type(orig, gene):
    s = _split(orig)
    if not s:
        return None
    prefix, allele = s
    for op in ("delins", "del", "ins", "dup", "inv"):
        if op in allele:
            up = allele.replace(op, op.upper(), 1)
            return f"{prefix}:{up}", orig
    return None


def inj_added_n_prefix(orig, gene):
    if len(orig) > 2 and orig[0] == "N" and orig[1] in "MRC" and orig[2] == "_":
        return orig[1:], orig  # drop the leading N
    return None


# (category name == rescue-op key, injector). Ordered for stable reporting.
INJECTORS = [
    ("STRIPPED_WHITESPACE", inj_whitespace),
    ("UPPERCASED_BASES", inj_uppercased_bases),
    ("STRIPPED_PROTEIN_SUFFIX", inj_protein_suffix),
    ("FIXED_GENE_WRAPPER", inj_gene_wrapper),
    ("SWAPPED_GENE_TRANSCRIPT", inj_swapped_gene_transcript),
    ("STRIPPED_SURROUNDING_PUNCTUATION", inj_surrounding_punctuation),
    ("FIXED_MULTIPLE_COLON", inj_double_colon),
    ("STRIPPED_UNBALANCED_BRACKETS", inj_unbalanced_brackets),
    ("FIXED_SEPARATOR_TYPO", inj_separator_typo),
    ("FIXED_MULTIPLE_DOT", inj_double_dot),
    ("STRIPPED_LEADING_JUNK", inj_leading_junk),
    ("FIXED_MULTIPLE_KIND", inj_double_kind),
    ("UPPERCASED_TRANSCRIPT", inj_uppercased_transcript),
    ("DROPPED_DEL_DUP_COUNT", inj_dropped_del_dup_count),
    ("ADDED_TRANSCRIPT_UNDERSCORE", inj_added_transcript_underscore),
    ("FIXED_PREFIX_COLON", inj_fixed_prefix_colon),
    ("LOWERCASED_MUTATION_TYPE", inj_lowercased_mutation_type),
    ("ADDED_N_PREFIX", inj_added_n_prefix),
]


# ---------------------------------------------------------------------------
# Sample handling
# ---------------------------------------------------------------------------

def build_sample(parser):
    """Extract clean, parseable c.HGVS strings from the public TSV sources."""
    seen = set()
    for tsv in SOURCE_TSVS:
        for line in tsv.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) == 2:
                seen.add(parts[1].strip())
    pool = []
    for c in sorted(seen):
        try:
            parser.parse_hgvs_variant(c)
            pool.append(c)
        except Exception:
            pass  # only keep strings that parse cleanly as a baseline
    rng = random.Random(SEED)
    rng.shuffle(pool)
    return pool


def load_sample(parser, regenerate):
    if SAMPLE_FILE.exists() and not regenerate:
        return [ln.strip() for ln in SAMPLE_FILE.read_text().splitlines()
                if ln.strip() and not ln.startswith("#")]
    pool = build_sample(parser)
    SAMPLE_FILE.write_text(
        "# Seeded sample of clean public ClinVar c.HGVS strings for the\n"
        "# Tier-1 injection cleaning benchmark (paper/scripts/inject_and_clean.py).\n"
        f"# Source: {', '.join(p.name for p in SOURCE_TSVS)}; seed={SEED}.\n"
        "# Regenerate with: python paper/scripts/inject_and_clean.py --regenerate\n"
        + "\n".join(pool)
        + "\n"
    )
    print(f"Wrote seeded sample ({len(pool)} strings) -> {SAMPLE_FILE}")
    return pool


# ---------------------------------------------------------------------------
# Benchmark
# ---------------------------------------------------------------------------

def parses(parser, s):
    try:
        parser.parse_hgvs_variant(s)
        return True
    except Exception:
        return False


def run(parser, pool):
    rng = random.Random(SEED)
    per_class = {}
    notes = []

    for idx, (name, injector) in enumerate(INJECTORS):
        # Build the applicable pool for this category (seeded order per class)
        order = list(pool)
        rng_c = random.Random(SEED + idx)
        rng_c.shuffle(order)

        attempted = recovered = parses_after = ineffective = expected_unparseable = 0
        applicable = 0
        gi = 0
        for orig in order:
            res = injector(orig, GENE_POOL[gi % len(GENE_POOL)])
            if res is None:
                continue
            applicable += 1
            if attempted >= TARGET_PER_CLASS:
                continue
            gi += 1
            perturbed, expected = res
            # Sanity: the expected canonical target must itself be valid HGVS.
            if not parses(parser, expected):
                expected_unparseable += 1
                continue
            if perturbed == expected:
                ineffective += 1  # injection was a no-op for this string
                continue
            attempted += 1
            cleaned, _ = clean_hgvs(perturbed)
            if cleaned == expected:
                recovered += 1
            if parses(parser, cleaned):
                parses_after += 1

        capped = applicable > attempted + ineffective + expected_unparseable
        if capped:
            notes.append(
                f"{name}: capped at {TARGET_PER_CLASS} of {applicable} applicable strings"
            )
        if expected_unparseable:
            notes.append(
                f"{name}: skipped {expected_unparseable} strings whose canonical "
                f"target did not parse"
            )
        per_class[name] = {
            "weight": REAL_RESCUE_OP_COUNTS[name],
            "n_attempted": attempted,
            "n_applicable": applicable,
            "n_recovered": recovered,
            "n_parses_after_clean": parses_after,
            "recovery_pct": round(100.0 * recovered / attempted, 1) if attempted else None,
        }

    # Regression check: cleaning a clean baseline string must not break it.
    regressions = 0
    for orig in pool:
        cleaned, _ = clean_hgvs(orig)
        if not parses(parser, cleaned):
            regressions += 1
            notes.append(f"REGRESSION: {orig} -> {cleaned} no longer parses")

    # Aggregate recovery
    total_attempted = sum(c["n_attempted"] for c in per_class.values())
    total_recovered = sum(c["n_recovered"] for c in per_class.values())
    overall_pct = round(100.0 * total_recovered / total_attempted, 1) if total_attempted else None

    # Weighted by the real-world rescue-op distribution (Tier-2 constants)
    wsum = wnum = 0.0
    for c in per_class.values():
        if c["recovery_pct"] is not None:
            wsum += c["weight"]
            wnum += c["weight"] * c["recovery_pct"]
    weighted_pct = round(wnum / wsum, 1) if wsum else None

    return {
        "issue": "SACGF/cdot#112",
        "tier": 1,
        "description": (
            "Reproducible injection cleaning benchmark: perturb public ClinVar "
            "c.HGVS with each clean_hgvs() fix category, then measure exact "
            "recovery. Per-class weights are frozen constants from the aggregate "
            "production rescue-op distribution (Tier 2; "
            "cdot_private/output/cleaning_analysis_20260617.txt, N=32752)."
        ),
        "sample_size": len(pool),
        "seed": SEED,
        "target_per_class": TARGET_PER_CLASS,
        "n_categories": len(INJECTORS),
        "overall_recovery_pct": overall_pct,
        "weighted_recovery_pct": weighted_pct,
        "total_attempted": total_attempted,
        "total_recovered": total_recovered,
        "regressions": regressions,
        "per_class": per_class,
        "notes": notes,
    }


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_facts(facts):
    FACTS_DIR.mkdir(parents=True, exist_ok=True)
    json_path = FACTS_DIR / "cleaning.json"
    json_path.write_text(json.dumps(facts, indent=2) + "\n")

    # Single-row CSV for vibepaper {{ cleaning.* }} templating.
    # TIER 1 ONLY — every value here is reproducible from public data committed
    # to this repo. Tier-2 production-corpus numbers (91.5%->96.4%, the residual
    # taxonomy) are deliberately NOT emitted as facts; they appear in the
    # manuscript as flagged literal constants (plan §3).
    csv_path = FACTS_DIR / "cleaning.csv"
    row = {
        "inject_overall_pct": facts["overall_recovery_pct"],
        "inject_weighted_pct": facts["weighted_recovery_pct"],
        "inject_total_attempted": facts["total_attempted"],
        "inject_total_recovered": facts["total_recovered"],
        "inject_regressions": facts["regressions"],
        "inject_n_categories": facts["n_categories"],
        "inject_sample_size": facts["sample_size"],
    }
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row))
        w.writeheader()
        w.writerow(row)
    return json_path, csv_path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--regenerate", action="store_true",
                    help="rebuild the seeded sample list from the source TSVs")
    args = ap.parse_args()

    parser = hgvs.parser.Parser()
    pool = load_sample(parser, args.regenerate)
    facts = run(parser, pool)
    json_path, csv_path = write_facts(facts)

    # Report
    print(f"\nTier-1 injection cleaning benchmark (issue #112)")
    print(f"  sample: {facts['sample_size']} clean public ClinVar c.HGVS strings (seed={SEED})")
    print(f"  {'category':<34} {'n':>4} {'recov':>6}  recovery%")
    for name, c in facts["per_class"].items():
        pct = "-" if c["recovery_pct"] is None else f"{c['recovery_pct']:.1f}"
        print(f"  {name:<34} {c['n_attempted']:>4} {c['n_recovered']:>6}  {pct:>8}")
    print(f"  {'-'*34}")
    print(f"  overall recovery : {facts['overall_recovery_pct']}%  "
          f"({facts['total_recovered']}/{facts['total_attempted']})")
    print(f"  weighted recovery: {facts['weighted_recovery_pct']}%  "
          f"(by real rescue-op distribution)")
    print(f"  regressions      : {facts['regressions']}")
    if facts["notes"]:
        print("  notes:")
        for n in facts["notes"]:
            print(f"    - {n}")
    print(f"\nWrote: {json_path}")
    print(f"Wrote: {csv_path}")


if __name__ == "__main__":
    main()
