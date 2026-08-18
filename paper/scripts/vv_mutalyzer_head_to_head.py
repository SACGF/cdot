#!/usr/bin/env python3
"""
Head-to-head: cdot ``clean_hgvs()`` vs VariantValidator and Mutalyzer (issue #112).

Extends ``lovd_head_to_head.py`` to the two sequence-aware validation services,
run over their public REST APIs on the exact injection corpus from
``inject_and_clean.py`` (same seed, caps and cases; see that script). Because
these are remote services there is no version pin; the service date and the
version metadata each API reports are recorded in the facts instead
(VariantValidator returns ``variantvalidator_version``; measured 2026-08-18).

Services
--------
- VariantValidator (rest.variantvalidator.org): batch pipe-delimited GET,
  ``/VariantValidator/variantvalidator/GRCh38/{v1|v2|...}/all``. Ensembl
  transcripts need the ``variantvalidator_ensembl`` endpoint, so requests are
  routed by the accession family of the case's canonical target. Responses
  for unparseable input embed the LOVD syntax checker's suggestions
  (``lovd_corrections``, via api.lovd.nl), which is recorded per case.
- Mutalyzer (mutalyzer.nl): ``/api/normalize/{description}``, one call per
  string, a few polite concurrent workers.

Scoring (identical rule to lovd_head_to_head.py)
------------------------------------------------
A case is recovered iff the service returns a description equal to the known
canonical target, compared gene-annotation-insensitively (see
``lovd_head_to_head.norm``). For VariantValidator the compared value is the
validated transcript description; for Mutalyzer either
``corrected_description`` or ``normalized_description`` may match (the
normalizer may legitimately re-shift a representation; accepting either is the
service's best answer). False-correction check on the uncorrupted originals:
the service returns a description differing from the valid input.

Requests are cached to JSONL (``output/rest_cache/``) keyed by input string,
so an interrupted run resumes without re-querying; re-running only re-scores.
VariantValidator's public service blocks sustained automated use (it tarpits
the connection after a few hundred requests, which is how any scripted user
experiences it); after repeated dead batches the fetcher gives up, VV is
scored over the cases it answered, and the facts record that coverage
(``vv_scored_cases`` / ``vv_coverage_pct``). A later resume refetches the gap.

Usage::

    python paper/scripts/vv_mutalyzer_head_to_head.py            # full run
    python paper/scripts/vv_mutalyzer_head_to_head.py --limit 3  # smoke test
"""

import argparse
import concurrent.futures
import csv
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import hgvs.parser

from cdot.hgvs.clean import clean_hgvs

sys.path.insert(0, str(Path(__file__).resolve().parent))
import inject_and_clean as iac              # noqa: E402
from lovd_head_to_head import build_cases, norm  # noqa: E402

FACTS_DIR = iac.FACTS_DIR
CACHE_DIR = iac.REPO / "output/rest_cache"
SERVICE_DATE = "2026-08-18"
USER_AGENT = "cdot-paper-benchmark (github.com/SACGF/cdot; davmlaw@gmail.com)"

VV_BASE = "https://rest.variantvalidator.org/VariantValidator"
MUT_BASE = "https://mutalyzer.nl/api/normalize/"
VV_BATCH = 5


def _get(url, timeout=300):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read())


def _get_deadline(url, deadline):
    """_get with a hard wall-clock deadline. Socket timeouts only bound
    *inactivity*, so a tarpitting server that drip-feeds bytes can hold a
    request open indefinitely; this abandons the request outright instead."""
    ex = concurrent.futures.ThreadPoolExecutor(1)
    try:
        return ex.submit(_get, url, deadline).result(timeout=deadline)
    finally:
        ex.shutdown(wait=False, cancel_futures=True)


# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------

class Cache:
    def __init__(self, path):
        self.path = path
        self.data = {}
        if path.exists():
            for line in path.read_text().splitlines():
                if line.strip():
                    rec = json.loads(line)
                    if rec.get("error"):
                        continue  # transient failure: refetch on resume
                    self.data[rec["input"]] = rec
        path.parent.mkdir(parents=True, exist_ok=True)
        self.fh = open(path, "a")

    def put(self, rec):
        self.data[rec["input"]] = rec
        self.fh.write(json.dumps(rec) + "\n")
        self.fh.flush()


# ---------------------------------------------------------------------------
# VariantValidator driver
# ---------------------------------------------------------------------------

def vv_endpoint(target):
    fam = "variantvalidator_ensembl" if target.startswith("ENST") else "variantvalidator"
    return f"{VV_BASE}/{fam}/GRCh38/"


def vv_parse(response, meta_out):
    """Map submitted_variant -> record from a (possibly batch) VV response."""
    out = {}
    meta_out.update(response.get("metadata") or {})
    for key, obj in response.items():
        if key in ("flag", "metadata") or not isinstance(obj, dict):
            continue
        sub = obj.get("submitted_variant")
        if sub is None:
            continue
        rec = out.setdefault(sub, {"descriptions": [], "warnings": [], "lovd_suggestions": False})
        desc = obj.get("hgvs_transcript_variant") or ""
        if desc:
            rec["descriptions"].append(desc)
        rec["warnings"].extend(obj.get("validation_warnings") or [])
        if obj.get("lovd_corrections"):
            rec["lovd_suggestions"] = True
    return out


def _collapse(s):
    return "".join(s.split())


def vv_fetch_batch(cache, strings, endpoint, state):
    """Fetch a batch. Timeouts and connection failures are throttling signals:
    retry the same batch with growing backoff (never bisect, which multiplies
    request load). Only an HTTP error response (a content problem) bisects,
    down to singles, so one bad case cannot poison the batch. Three batches in
    a row exhausting their retries abandons the VV fetch for this run (the
    service is blocking sustained automated use; whatever is cached gets
    scored, with coverage reported, and a later resume refetches the rest)."""
    todo = [s for s in strings if s not in cache.data]
    if not todo or state["abandoned"]:
        return
    url = (endpoint + urllib.parse.quote("|".join(todo), safe="") +
           "/all?content-type=application%2Fjson")
    err = None
    for attempt in range(4):
        try:
            meta = {}
            parsed = vv_parse(_get_deadline(url, deadline=120), meta)
            # VV normalises whitespace before echoing submitted_variant, so
            # fall back to a whitespace-collapsed lookup (batches never contain
            # two strings sharing a collapsed form; vv_fetch_all guards that).
            collapsed = {_collapse(k): v for k, v in parsed.items()}
            for s in todo:
                rec = parsed.get(s) or collapsed.get(_collapse(s)) or \
                    {"descriptions": [], "warnings": ["NO_RESULT"],
                     "lovd_suggestions": False}
                rec = dict(rec)
                rec["input"] = s
                rec["vv_version"] = meta.get("variantvalidator_version", "")
                cache.put(rec)
            state["consec_fail"] = 0
            time.sleep(2)  # pacing between batch requests
            return
        except urllib.error.HTTPError as e:
            err = e
            break
        except Exception as e:  # timeout / tarpit / connection reset: back off
            err = e
            wait = min(600, 60 * 2 ** attempt)
            print(f"  VV throttled ({str(e)[:60] or type(e).__name__}); "
                  f"backing off {wait}s", flush=True)
            time.sleep(wait)
    else:
        # Deadline exhaustion is a service-level block, not a content problem:
        # bisecting would just multiply load. Leave the batch uncached and
        # count the failure toward giving up.
        state["consec_fail"] += 1
        if state["consec_fail"] >= 3:
            state["abandoned"] = True
            print("  VV: giving up for this run (service is blocking sustained "
                  "automated use); scoring the cases fetched so far", flush=True)
        return
    if len(todo) == 1:
        cache.put({"input": todo[0], "descriptions": [], "warnings": [],
                   "lovd_suggestions": False, "error": str(err)[:200]})
    else:
        mid = len(todo) // 2
        vv_fetch_batch(cache, todo[:mid], endpoint, state)
        vv_fetch_batch(cache, todo[mid:], endpoint, state)


def vv_fetch_all(cache, items, log_every=20):
    """items: list of (input_string, target). Routed and batched; a batch is
    flushed early if it would contain a duplicate input or predicted target
    (batch responses are keyed by corrected description, which must stay
    unambiguous within one request)."""
    groups = {}
    for s, target in items:
        groups.setdefault(vv_endpoint(target), []).append((s, target))
    n_req = 0
    state = {"consec_fail": 0, "abandoned": False}
    for endpoint, pairs in groups.items():
        batch, seen = [], set()
        pending = [p for p in pairs if p[0] not in cache.data]
        for s, target in pending:
            if state["abandoned"]:
                return
            key = norm(target)
            if len(batch) >= VV_BATCH or _collapse(s) in seen or key in seen:
                vv_fetch_batch(cache, batch, endpoint, state)
                n_req += 1
                if n_req % log_every == 0:
                    done = sum(1 for p in pairs if p[0] in cache.data)
                    print(f"  VV {endpoint.split('/')[-3]}: {done}/{len(pairs)}",
                          flush=True)
                batch, seen = [], set()
            batch.append(s)
            seen.add(_collapse(s))
            seen.add(key)
        if batch and not state["abandoned"]:
            vv_fetch_batch(cache, batch, endpoint, state)


# ---------------------------------------------------------------------------
# Mutalyzer driver
# ---------------------------------------------------------------------------

def mut_fetch_one(s):
    url = MUT_BASE + urllib.parse.quote(s, safe="")
    for attempt in range(3):
        try:
            d = _get(url, timeout=120)
            return {"input": s,
                    "normalized": d.get("normalized_description", ""),
                    "corrected": d.get("corrected_description", ""),
                    "errors": []}
        except urllib.error.HTTPError as e:
            if e.code < 500:
                try:
                    body = json.loads(e.read())
                    body = body.get("custom", body)
                    codes = [x.get("code") for x in (body.get("errors") or [])]
                except Exception:
                    codes = [f"HTTP{e.code}"]
                return {"input": s, "normalized": "", "corrected": "", "errors": codes}
            time.sleep(2 ** (attempt + 1))
        except Exception as e:
            if attempt == 2:
                return {"input": s, "normalized": "", "corrected": "",
                        "errors": [f"FAIL:{str(e)[:100]}"]}
            time.sleep(2 ** (attempt + 1))
    return {"input": s, "normalized": "", "corrected": "", "errors": ["RETRIES"]}


def mut_fetch_all(cache, strings, workers=5, log_every=200):
    todo = [s for s in dict.fromkeys(strings) if s not in cache.data]
    with concurrent.futures.ThreadPoolExecutor(workers) as pool:
        for i, rec in enumerate(pool.map(mut_fetch_one, todo)):
            cache.put(rec)
            if (i + 1) % log_every == 0:
                print(f"  Mutalyzer: {i + 1}/{len(todo)}", flush=True)


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def vv_hit(rec, target):
    return any(norm(d) == norm(target) for d in rec.get("descriptions", []))


def mut_hit(rec, target):
    t = norm(target)
    return norm(rec.get("corrected") or "") == t or norm(rec.get("normalized") or "") == t


def source_family(target):
    return "ensembl" if target.startswith("ENST") else "refseq"


_INTRONIC = __import__("re").compile(r"\d+[+-]\d+")


def is_intronic(target):
    """Intronic (or UTR-offset) position: Mutalyzer and LOVD both reject these
    on a transcript reference by design (EINTRONIC / EWRONGREFERENCE)."""
    return bool(_INTRONIC.search(target.split(":", 1)[-1]))


def score(cases, vv_cache, mut_cache):
    per_class = {}
    per_source = {"refseq": {}, "ensembl": {}}
    for name, _ in iac.INJECTORS:
        per_class[name] = {"weight": iac.REAL_RESCUE_OP_COUNTS[name],
                           "n_attempted": 0, "n_vv": 0, "cdot": 0, "vv": 0,
                           "mut": 0, "vv_lovd_suggestions": 0}
    for src in per_source.values():
        src.update(n=0, n_vv=0, cdot=0, vv=0, mut=0)
    nonintronic = {"n": 0, "n_vv": 0, "mut": 0, "vv": 0, "cdot": 0}

    for name, orig, perturbed, expected in cases:
        c = per_class[name]
        s = per_source[source_family(expected)]
        c["n_attempted"] += 1
        s["n"] += 1
        if not is_intronic(expected):
            nonintronic["n"] += 1
        cleaned, _ = clean_hgvs(perturbed)
        if norm(cleaned) == norm(expected):
            c["cdot"] += 1
            s["cdot"] += 1
            if not is_intronic(expected):
                nonintronic["cdot"] += 1
        # VV is scored only over the cases it answered before blocking
        # sustained automated use; coverage is reported alongside.
        vrec = vv_cache.data.get(perturbed)
        if vrec is not None:
            c["n_vv"] += 1
            s["n_vv"] += 1
            if not is_intronic(expected):
                nonintronic["n_vv"] += 1
            if vv_hit(vrec, expected):
                c["vv"] += 1
                s["vv"] += 1
                if not is_intronic(expected):
                    nonintronic["vv"] += 1
            if vrec.get("lovd_suggestions"):
                c["vv_lovd_suggestions"] += 1
        if mut_hit(mut_cache.data[perturbed], expected):
            c["mut"] += 1
            s["mut"] += 1
            if not is_intronic(expected):
                nonintronic["mut"] += 1

    def pct(num, den):
        return round(100.0 * num / den, 1) if den else None

    for c in per_class.values():
        c["cdot_pct"] = pct(c["cdot"], c["n_attempted"])
        c["mut_pct"] = pct(c["mut"], c["n_attempted"])
        c["vv_pct"] = pct(c["vv"], c["n_vv"])
    for s in per_source.values():
        s["cdot_pct"] = pct(s["cdot"], s["n"])
        s["mut_pct"] = pct(s["mut"], s["n"])
        s["vv_pct"] = pct(s["vv"], s["n_vv"])
    nonintronic["cdot_pct"] = pct(nonintronic["cdot"], nonintronic["n"])
    nonintronic["mut_pct"] = pct(nonintronic["mut"], nonintronic["n"])
    nonintronic["vv_pct"] = pct(nonintronic["vv"], nonintronic["n_vv"])
    return per_class, per_source, nonintronic


def weighted(per_class, key):
    wsum = wnum = 0.0
    for c in per_class.values():
        if c[key + "_pct"] is not None:
            wsum += c["weight"]
            wnum += c["weight"] * c[key + "_pct"]
    return round(wnum / wsum, 1) if wsum else None


def false_corrections(pool, vv_cache, mut_cache):
    """On valid inputs, distinguish *altered* (service returns a different
    description: a true false correction) from *rejected* (service returns
    nothing: a validity design position or coverage gap, eg Mutalyzer's
    EINTRONIC for intronic positions on a transcript reference, or a
    transcript missing from VV's database)."""
    counts = {"vv_scored": 0, "vv_altered": 0, "vv_rejected": 0,
              "mut_altered": 0, "mut_rejected": 0, "mut_repr_changes": 0}
    for orig in pool:
        t = norm(orig)
        vrec = vv_cache.data.get(orig)
        if vrec is not None:
            counts["vv_scored"] += 1
            if not vv_hit(vrec, orig):
                if vrec.get("descriptions"):
                    counts["vv_altered"] += 1
                else:
                    counts["vv_rejected"] += 1
        mrec = mut_cache.data[orig]
        if not mut_hit(mrec, orig):
            if mrec.get("normalized") or mrec.get("corrected"):
                counts["mut_altered"] += 1
            else:
                counts["mut_rejected"] += 1
        elif norm(mrec.get("normalized") or "") != t:
            counts["mut_repr_changes"] += 1  # accepted, representation re-shifted
    return counts


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit", type=int, default=None,
                    help="cap cases per category (smoke test only; not for facts)")
    ap.add_argument("--cache-dir", type=Path, default=CACHE_DIR)
    args = ap.parse_args()

    parser = hgvs.parser.Parser()
    pool = iac.load_sample(parser, regenerate=False)
    cases = build_cases(parser, pool)
    if args.limit:
        by_cat = {}
        cases = [c for c in cases
                 if by_cat.setdefault(c[0], []).append(c) or len(by_cat[c[0]]) <= args.limit]
        pool = pool[: args.limit * 10]

    vv_cache = Cache(args.cache_dir / "vv.jsonl")
    mut_cache = Cache(args.cache_dir / "mutalyzer.jsonl")

    items = [(perturbed, expected) for _, _, perturbed, expected in cases]
    items += [(orig, orig) for orig in pool]

    t0 = time.time()
    print(f"Mutalyzer: {sum(1 for s, _ in items if s not in mut_cache.data)} of "
          f"{len(items)} strings uncached")
    print(f"VariantValidator: {sum(1 for s, _ in items if s not in vv_cache.data)} of "
          f"{len(items)} strings uncached")
    # Mutalyzer first: it tolerates steady traffic, and the sequencing gives
    # VariantValidator's anti-abuse tarpit a quiet period before VV resumes.
    mut_fetch_all(mut_cache, [s for s, _ in items])
    vv_fetch_all(vv_cache, items)
    fetch_minutes = round((time.time() - t0) / 60, 1)

    per_class, per_source, nonintronic = score(cases, vv_cache, mut_cache)
    n = sum(c["n_attempted"] for c in per_class.values())
    n_vv = sum(c["n_vv"] for c in per_class.values())
    totals = {k: sum(c[k] for c in per_class.values()) for k in ("cdot", "vv", "mut")}
    fc = false_corrections(pool, vv_cache, mut_cache)
    vv_version = next((r.get("vv_version") for r in vv_cache.data.values()
                       if r.get("vv_version")), "")

    facts = {
        "issue": "SACGF/cdot#112",
        "tier": 1,
        "description": (
            "Head-to-head on the reproducible injection corpus: cdot clean_hgvs() "
            "vs VariantValidator and Mutalyzer public REST services, measured "
            f"{SERVICE_DATE}; same cases and scoring as the LOVD comparison "
            "(see paper/scripts/vv_mutalyzer_head_to_head.py)."
        ),
        "service_date": SERVICE_DATE,
        "vv_version": vv_version,
        "n_cases": n,
        "sample_size": len(pool),
        "vv_scored_cases": n_vv,
        "vv_coverage_pct": round(100.0 * n_vv / n, 1) if n else None,
        "cdot_pct": round(100.0 * totals["cdot"] / n, 1),
        "vv_pct": round(100.0 * totals["vv"] / n_vv, 1) if n_vv else None,
        "mut_pct": round(100.0 * totals["mut"] / n, 1),
        "cdot_weighted_pct": weighted(per_class, "cdot"),
        "vv_weighted_pct": weighted(per_class, "vv"),
        "mut_weighted_pct": weighted(per_class, "mut"),
        "vv_refseq_pct": per_source["refseq"]["vv_pct"],
        "vv_ensembl_pct": per_source["ensembl"]["vv_pct"],
        "mut_refseq_pct": per_source["refseq"]["mut_pct"],
        "mut_ensembl_pct": per_source["ensembl"]["mut_pct"],
        "n_nonintronic": nonintronic["n"],
        "cdot_nonintronic_pct": nonintronic["cdot_pct"],
        "vv_nonintronic_pct": nonintronic["vv_pct"],
        "mut_nonintronic_pct": nonintronic["mut_pct"],
        "originals_n": len(pool),
        "vv_originals_scored": fc["vv_scored"],
        "vv_false_corrections": fc["vv_altered"],
        "vv_rejected_valid": fc["vv_rejected"],
        "mut_false_corrections": fc["mut_altered"],
        "mut_rejected_valid": fc["mut_rejected"],
        "mut_representation_changes": fc["mut_repr_changes"],
        "vv_lovd_suggestion_cases": sum(c["vv_lovd_suggestions"] for c in per_class.values()),
        "fetch_minutes": fetch_minutes,
        "per_class": per_class,
        "per_source": per_source,
    }

    FACTS_DIR.mkdir(parents=True, exist_ok=True)
    json_path = FACTS_DIR / "vv_mutalyzer_comparison.json"
    json_path.write_text(json.dumps(facts, indent=2) + "\n")
    csv_path = FACTS_DIR / "vv_mutalyzer_comparison.csv"
    row = {k: facts[k] for k in (
        "service_date", "vv_version", "n_cases", "vv_scored_cases",
        "vv_coverage_pct", "cdot_pct", "vv_pct", "mut_pct",
        "cdot_weighted_pct", "vv_weighted_pct", "mut_weighted_pct",
        "vv_refseq_pct", "vv_ensembl_pct", "mut_refseq_pct", "mut_ensembl_pct",
        "n_nonintronic", "cdot_nonintronic_pct", "vv_nonintronic_pct",
        "mut_nonintronic_pct", "originals_n", "vv_originals_scored",
        "vv_false_corrections", "vv_rejected_valid", "mut_false_corrections",
        "mut_rejected_valid", "mut_representation_changes",
        "vv_lovd_suggestion_cases")}
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row))
        w.writeheader()
        w.writerow(row)

    def s(x):
        return "-" if x is None else str(x)

    print(f"\nVV/Mutalyzer head-to-head ({SERVICE_DATE}, {n} injected cases, "
          f"VV answered {n_vv} ({facts['vv_coverage_pct']}%), "
          f"fetch {fetch_minutes} min, VV {vv_version})")
    print(f"  {'category':<34} {'n':>4} {'cdot%':>7} {'n_vv':>5} {'vv%':>7} {'mut%':>7}")
    for name, c in per_class.items():
        if not c["n_attempted"]:
            continue
        print(f"  {name:<34} {c['n_attempted']:>4} {s(c['cdot_pct']):>7} "
              f"{c['n_vv']:>5} {s(c['vv_pct']):>7} {s(c['mut_pct']):>7}")
    print(f"  {'-'*68}")
    print(f"  overall : cdot {s(facts['cdot_pct'])}%  vv {s(facts['vv_pct'])}% "
          f"(of answered)  mut {s(facts['mut_pct'])}%")
    print(f"  weighted: cdot {s(facts['cdot_weighted_pct'])}%  "
          f"vv {s(facts['vv_weighted_pct'])}%  mut {s(facts['mut_weighted_pct'])}%")
    print(f"  by source: vv refseq {s(facts['vv_refseq_pct'])}% ensembl "
          f"{s(facts['vv_ensembl_pct'])}% | mut refseq {s(facts['mut_refseq_pct'])}% "
          f"ensembl {s(facts['mut_ensembl_pct'])}%")
    print(f"  non-intronic ({nonintronic['n']}): cdot "
          f"{s(facts['cdot_nonintronic_pct'])}%  vv {s(facts['vv_nonintronic_pct'])}%  "
          f"mut {s(facts['mut_nonintronic_pct'])}%")
    print(f"  originals ({len(pool)}, vv answered {fc['vv_scored']}): "
          f"vv altered {fc['vv_altered']} / rejected {fc['vv_rejected']}; "
          f"mut altered {fc['mut_altered']} / rejected {fc['mut_rejected']} "
          f"(+{fc['mut_repr_changes']} accepted with re-shifted representation)")
    print(f"  vv responses embedding LOVD suggestions: "
          f"{facts['vv_lovd_suggestion_cases']}/{n_vv}")
    print(f"\nWrote: {json_path}")
    print(f"Wrote: {csv_path}")


if __name__ == "__main__":
    main()
