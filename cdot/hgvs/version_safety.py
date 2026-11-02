"""
Transcript-version substitution safety (issue #28).

When the exact transcript version a variant cites is absent from the loaded data,
cdot's opt-in fallback (:func:`cdot.hgvs.gene_hgvs.resolve_transcript_version`)
substitutes the nearest available version.  That is only safe if the variant's
coding (c.) position still maps to the *same* genomic coordinate after the
substitution — otherwise the substitution silently changes what the variant means.

The study behind this (see ``docs/transcript_version_safety.md`` and
``analysis/transcript_version_drift.py`` / ``analysis/transcript_version_crossbuild.py``)
established a near-deterministic test for this: a version bump is
coordinate-preserving **iff** the version's *intrinsic CDS structure* — the CDS
length plus the coding-exon-segment lengths *in transcript coordinates* — is
unchanged.  Crucially that structure is **build-independent** (only the genomic
placement of a version changes between builds), so it can be read for the
requested version from *any* build that holds it (cdot's all-builds JSON, or the
REST ``/transcript/<ac>`` view which returns every version across builds).  This
catches the transient-revert versions (A→B→A) that a genomic-bracket check is
blind to, and needs no flanking pair.

This module provides the pure, build-independent structural helper
:func:`intrinsic_cds_structure` plus the private genomic-bracket helpers used only
for the residual case where the requested version exists in no build at all.  The
provider-level check that ties them together is
``AbstractJSONDataProvider.is_version_substitution_safe`` (in
``cdot/hgvs/dataproviders/json_data_provider.py``).

Method note (why this is exact and cheap, mirroring the analysis scripts): cdot
stores the genome *alignment*, not the sequence.  Each exon is
``[alt_start, alt_end, exon_id, cds_start, cds_end, gap]``, a piecewise-linear map
between 1-based spliced-transcript coordinates (``cds_start``/``cds_end``) and
genomic coordinates.  ``c.X`` anchors at the start codon.  The intrinsic structure
is derived purely from the transcript coordinates, which are identical across
builds.  Those transcript coordinates are unaffected by alignment gaps, but the
c.-to-genomic *mapping* is not: a transcript-vs-genome indel (an alignment gap)
shifts every coding base downstream of it, so two versions with identical
intrinsic structure can still place a coding base at different genomic
coordinates.  Identical intrinsic structure is therefore necessary but not
sufficient; the safety check (:meth:`is_version_substitution_safe`) also requires
the CDS alignment gaps to match (:func:`cds_alignment_gaps`).  Gaps are
build-specific, so that gap check is exact when both versions are placed in the
same build and conservative across builds.  Note the structure validates only
coding (c.) positions; a UTR/intronic position can still move if the surrounding
non-coding exon structure changes between versions.

Non-coding transcripts (``NR_``, n. positions) have no CDS, so they are compared on
their whole exon structure instead (:func:`intrinsic_exon_structure`: strand plus
per-exon transcript lengths), with the same alignment-gap and same-build placement
checks (:func:`exon_alignment_gaps`, :func:`exon_genomic_spans`).
"""

import re
from typing import Optional


# ---------------------------------------------------------------------------
# Build-independent intrinsic CDS structure (the primary safety signal)
# ---------------------------------------------------------------------------

def intrinsic_cds_structure(transcript, genome_build: Optional[str] = None):
    """Return a version's build-independent CDS structure, or ``None``.

    The structure is ``(strand, cds_len, seg_lengths)`` where ``seg_lengths`` is
    the tuple of coding-exon-segment lengths *in transcript coordinates* (the CDS
    partitioned at every coding-exon boundary).  Two transcript versions whose
    intrinsic structures are equal map every coding c. position to the same
    genomic coordinate, so substituting one for the other is coordinate-safe; a
    change in CDS length or in where the coding-exon boundaries fall is exactly
    what moves a coordinate.

    ``transcript`` is a cdot transcript record — either the plain JSON dict or a
    :class:`cdot.models.Transcript` (both support ``["key"]`` / ``.get(...)`` and
    positional exon access).  It must carry ``start_codon`` / ``stop_codon`` (a
    coding transcript) and at least one genome build.

    The structure is read from ``genome_build`` if that build is present;
    otherwise it is read from *any* available build, because the transcript
    coordinates are build-independent.  This is what lets a substitution be judged
    even when the requested version is absent from the target build but present in
    another (the whole point of the cross-build result).  Returns ``None`` for
    non-coding transcripts, a missing/zero-length CDS, a record with no usable
    build, or a malformed record (an exon missing its transcript coordinates).  The transcript coordinates this reads are exact regardless of alignment
    gaps, but gaps still move the genomic mapping, so an equal structure must be
    paired with :func:`cds_alignment_gaps` equality before a substitution is called
    coordinate-safe.
    """
    sc = transcript.get("start_codon")
    ec = transcript.get("stop_codon")
    if sc is None or ec is None:
        return None
    cds_len = ec - sc
    if cds_len <= 0:
        return None

    gb = _pick_build(transcript, genome_build)
    if gb is None:
        return None

    cds_lo = sc + 1          # first CDS transcript position (1-based)
    cds_hi = sc + cds_len    # last CDS transcript position

    # CDS-offset of the start of each coding-exon segment (0-based from the CDS
    # start). Derived only from the exons' transcript coordinates, so it is the
    # same in every build.
    exons = _exon_tx_coords(gb)
    if exons is None:
        return None
    seg_starts = []
    for tx_start, tx_end, _gap in exons:
        a = max(tx_start, cds_lo)
        b = min(tx_end, cds_hi)
        if a > b:
            continue
        seg_starts.append(a - cds_lo)

    if not seg_starts:
        return None
    seg_starts.sort()
    bounds = seg_starts + [cds_len]
    seg_lengths = tuple(hi - lo for lo, hi in zip(bounds, bounds[1:]))
    return (gb["strand"], cds_len, seg_lengths)


def cds_alignment_gaps(transcript, genome_build: Optional[str] = None):
    """Return a version's CDS alignment-gap signature in a build, or ``None``.

    :func:`intrinsic_cds_structure` compares *transcript* coordinates, which are
    blind to the transcript-to-genome alignment: two versions with identical CDS
    structure can still map a coding base to a *different* genomic coordinate when
    their alignment differs by a gap (a transcript-vs-genome indel shifts every
    coding base downstream of it by the gap length). So identical intrinsic
    structure is necessary but not sufficient for coordinate safety; the CDS
    alignment gaps must match too.

    The signature is a tuple of ``(coding_segment_index, gap_string)`` for every
    coding exon that carries an alignment gap, in CDS order (an empty tuple when
    the CDS aligns cleanly). Only exons that overlap the CDS are considered, since
    a gap entirely in a non-coding exon cannot move a coding coordinate.

    Unlike the intrinsic structure, alignment gaps are build-specific (a version
    can align cleanly in one build and with a gap in another), so this is read
    from ``genome_build`` when the record holds it, otherwise from any available
    build. Two versions both placed in ``genome_build`` are therefore compared
    exactly; a cross-build comparison (one side read from another build) is best
    effort and conservative. Returns ``None`` for the same non-coding / missing-CDS
    / no-build / malformed cases as :func:`intrinsic_cds_structure`.
    """
    sc = transcript.get("start_codon")
    ec = transcript.get("stop_codon")
    if sc is None or ec is None:
        return None
    cds_len = ec - sc
    if cds_len <= 0:
        return None

    gb = _pick_build(transcript, genome_build)
    if gb is None:
        return None

    cds_lo = sc + 1
    cds_hi = sc + cds_len

    exons = _exon_tx_coords(gb)
    if exons is None:
        return None
    gaps = []
    seg_idx = 0
    for tx_start, tx_end, gap in sorted(exons, key=lambda e: e[0]):
        if tx_start > cds_hi or tx_end < cds_lo:
            continue  # non-coding exon: its gap cannot move a coding coordinate
        if gap:
            gaps.append((seg_idx, gap))
        seg_idx += 1
    return tuple(gaps)


def utr_lengths(transcript, genome_build: Optional[str] = None):
    """Return ``(five_prime_len, three_prime_len)`` in transcript bases, or ``None``.

    The intrinsic CDS structure normalises the UTRs away, so it cannot tell whether
    a *non-coding* cited position (a 5'UTR ``c.-N`` or a 3'UTR ``c.*N``) still maps
    to the same genomic coordinate after a substitution. The 5'UTR length is the
    start-codon offset; the 3'UTR length is the total spliced cDNA length minus the
    stop-codon offset. Both are transcript coordinates, so (like
    :func:`intrinsic_cds_structure`) they are build-independent and gap-independent;
    a build only changes a version's genomic placement, not these lengths.

    Returns ``None`` for non-coding / missing-CDS / no-build / malformed records.
    """
    sc = transcript.get("start_codon")
    ec = transcript.get("stop_codon")
    if sc is None or ec is None:
        return None
    gb = _pick_build(transcript, genome_build)
    if gb is None:
        return None
    exons = _exon_tx_coords(gb)
    if not exons:
        return None
    # cds_start/cds_end are cumulative 1-based spliced transcript coordinates, so
    # the largest cds_end is the total spliced cDNA length.
    total_cdna = max(tx_end for _tx_start, tx_end, _gap in exons)
    return (sc, total_cdna - ec)


# A cited c. position is in the 5'UTR when it is negative (leading "-", or a range
# whose endpoint is negative) and in the 3'UTR when it carries a "*". Coding and
# CDS-intronic positions touch neither, so the CDS-structure check already covers
# them. We return the set of UTRs a position/range touches so a range spanning a
# UTR boundary is handled conservatively. A position can also start after "(" in an
# uncertain range ("(-30_-20)del", "(?_-30)_(5_?)del"), so "(" counts as a position
# start too. A CDS-intronic offset ("2604-1") follows a digit, so never matches.
_FIVE_PRIME = re.compile(r"(?:^|[_(])-\d")   # "-49", "-49+10", "-3_5", "(-30_-20)"...
_THREE_PRIME = re.compile(r"\*")          # any 3'UTR component


def utr_regions_touched(cited_position: str) -> set:
    """Return the UTRs (``{'five_prime', 'three_prime'}``) a cited c. position or
    range touches; empty for a purely coding/CDS-intronic position.

    ``cited_position`` is the text after ``c.``/``n.`` (the position plus edit, eg
    ``"-49="``, ``"*100A>G"``, ``"123+5G>A"``). Only the position syntax is
    inspected; ordinary edits (``>``, ``del``, ``dup``, ``ins``) carry no ``*`` and
    no leading ``-``.
    """
    regions = set()
    if not cited_position:
        return regions
    if _THREE_PRIME.search(cited_position):
        regions.add("three_prime")
    if _FIVE_PRIME.search(cited_position):
        regions.add("five_prime")
    return regions


def describe_structure_change(req_struct, sub_struct) -> str:
    """Return a short reason describing how two intrinsic structures differ.

    Assumes ``req_struct != sub_struct``.  Used to explain a not-coordinate-safe
    verdict to the caller (and ultimately in the HGVSFix message).
    """
    req_strand, req_len, _req_segs = req_struct
    sub_strand, sub_len, _sub_segs = sub_struct
    if req_strand != sub_strand:
        return f"strand changed {req_strand}→{sub_strand}"
    if req_len != sub_len:
        return f"CDS length changed {req_len}→{sub_len} bp"
    return "coding-exon structure changed (same CDS length)"


# ---------------------------------------------------------------------------
# Non-coding (n.) transcripts
# ---------------------------------------------------------------------------

def has_cds(transcript) -> bool:
    """True if the record carries a start and stop codon (a coding transcript)."""
    return transcript.get("start_codon") is not None and transcript.get("stop_codon") is not None


def intrinsic_exon_structure(transcript, genome_build: Optional[str] = None):
    """Return a version's build-independent exon structure, or ``None``.

    The n. analogue of :func:`intrinsic_cds_structure`, for non-coding transcripts.
    The structure is ``(strand, exon_lengths)`` where ``exon_lengths`` is the tuple
    of per-exon lengths *in transcript coordinates*, in transcript order. An n.
    position (exonic, intronic as an offset from an exon boundary, or ``n.-N`` /
    ``n.*N`` outside the transcript) is counted from the transcript start and the
    exon boundaries, so it depends only on these lengths plus the genomic placement
    and alignment gaps of the exons. Equal structure is therefore necessary for an
    n. substitution to be coordinate-safe, and is paired with
    :func:`exon_alignment_gaps` (and, when both versions are placed in the target
    build, :func:`exon_genomic_spans`) exactly like the coding check.

    Read from ``genome_build`` if present, else any build. Returns ``None`` for a
    record with no usable build or a malformed record.
    """
    gb = _pick_build(transcript, genome_build)
    if gb is None:
        return None
    exons = _exon_tx_coords(gb)
    if not exons:
        return None
    exons = sorted(exons, key=lambda e: e[0])
    return (gb["strand"], tuple(tx_end - tx_start + 1 for tx_start, tx_end, _gap in exons))


def exon_alignment_gaps(transcript, genome_build: Optional[str] = None):
    """Return ``(exon_index, gap_string)`` for every gapped exon, in transcript order.

    The n. analogue of :func:`cds_alignment_gaps`: every exon counts, since an n.
    position can sit in any exon. Build-specific, read from ``genome_build`` when
    held, else any build. Returns ``None`` for a no-build / malformed record.
    """
    gb = _pick_build(transcript, genome_build)
    if gb is None:
        return None
    exons = _exon_tx_coords(gb)
    if exons is None:
        return None
    exons = sorted(exons, key=lambda e: e[0])
    return tuple((i, gap) for i, (_s, _e, gap) in enumerate(exons) if gap)


def exon_genomic_spans(transcript, genome_build: str):
    """Return ``(contig, strand, spans)`` for a version placed in ``genome_build``.

    ``spans`` is the sorted tuple of ``(alt_start, alt_end)`` genomic exon spans.
    With equal :func:`intrinsic_exon_structure` and equal alignment gaps, equal
    spans mean every n. position maps to the same genomic coordinate. Returns
    ``None`` if the record is not placed in that build or is malformed.
    """
    gb = transcript["genome_builds"].get(genome_build)
    if gb is None or _exon_tx_coords(gb) is None:
        return None
    spans = []
    for exon in gb["exons"]:
        if exon[0] is None or exon[1] is None:
            return None
        spans.append((exon[0], exon[1]))
    return gb["contig"], gb["strand"], tuple(sorted(spans))


# ---------------------------------------------------------------------------
# Genomic-bracket fallback (only when the requested version is in no build)
# ---------------------------------------------------------------------------

def cds_genomic_map(transcript, genome_build: str):
    """Exact genomic signature of a version's CDS map in one build, or ``None``.

    Returns a dict with ``contig``, ``strand``, ``cds_len``, ``sign`` and ``segs``
    (a sorted list of ``(cds_offset, genomic_coordinate)`` breakpoints, one per
    coding-exon segment).  ``has_gap`` flags an alignment gap, which makes the
    unit-slope assumption inexact; callers treat a gapped map as non-comparable.
    Used only by the genomic-bracket fallback — the primary path is structural.
    Returns ``None`` for a non-coding or malformed record.
    """
    sc = transcript.get("start_codon")
    ec = transcript.get("stop_codon")
    if sc is None or ec is None:
        return None
    cds_len = ec - sc
    if cds_len <= 0:
        return None
    gb = transcript["genome_builds"].get(genome_build)
    if gb is None:
        return None

    strand = gb["strand"]
    sign = 1 if strand == "+" else -1
    cds_lo = sc + 1
    cds_hi = sc + cds_len

    if _exon_tx_coords(gb) is None:
        return None
    segs = []
    has_gap = False
    for exon in gb["exons"]:
        g_start, g_end, tx_start, tx_end = exon[0], exon[1], exon[3], exon[4]
        gap = exon[5] if len(exon) > 5 else None
        if g_start is None or g_end is None:
            return None
        if gap:
            has_gap = True
        a = max(tx_start, cds_lo)
        b = min(tx_end, cds_hi)
        if a > b:
            continue
        off_in_exon = a - tx_start
        g_at_a = (g_start + off_in_exon) if sign == 1 else (g_end - 1 - off_in_exon)
        segs.append((a - cds_lo, g_at_a))

    if not segs:
        return None
    segs.sort()
    return {
        "contig": gb["contig"], "strand": strand, "cds_len": cds_len,
        "sign": sign, "segs": segs, "has_gap": has_gap,
    }


def _genomic_at(sig, off: int) -> int:
    """Genomic coordinate of CDS offset ``off`` (0-based) under signature ``sig``."""
    segs = sig["segs"]
    lo, hi = 0, len(segs) - 1
    idx = 0
    while lo <= hi:
        mid = (lo + hi) // 2
        if segs[mid][0] <= off:
            idx = mid
            lo = mid + 1
        else:
            hi = mid - 1
    off0, g0 = segs[idx]
    return g0 + sig["sign"] * (off - off0)


def genomic_maps_fully_agree(a, b) -> bool:
    """True iff two gap-free genomic CDS maps agree at *every* shared coding base.

    Exact in O(exons): the maps are unit-slope piecewise-linear, so they agree on
    a sub-interval iff they agree at its left breakpoint.  Returns ``False`` if
    either map carries an alignment gap (non-comparable) or they sit on a
    different contig/strand.
    """
    if a is None or b is None or a["has_gap"] or b["has_gap"]:
        return False
    if a["contig"] != b["contig"] or a["strand"] != b["strand"]:
        return False
    L = min(a["cds_len"], b["cds_len"])
    if L <= 0:
        return False
    bounds = sorted({0, L}
                    | {o for o, _g in a["segs"] if 0 < o < L}
                    | {o for o, _g in b["segs"] if 0 < o < L})
    for x0 in bounds[:-1]:
        if _genomic_at(a, x0) != _genomic_at(b, x0):
            return False
    return True


# ---------------------------------------------------------------------------
# Internal
# ---------------------------------------------------------------------------

def _exon_tx_coords(gb):
    """Return ``[(tx_start, tx_end, gap), ...]`` for a build's exons, or ``None``.

    The single place exon records are unpacked, so every helper treats a malformed
    record the same way. An exon is ``[alt_start, alt_end, exon_id, cds_start,
    cds_end, gap]``; the trailing gap may be absent (treated as no gap). A record
    with an exon that is too short or whose transcript coordinates are ``None`` is
    malformed and yields ``None``, so callers report "cannot verify" rather than
    raise.
    """
    out = []
    for exon in gb.get("exons") or []:
        if len(exon) < 5:
            return None
        tx_start, tx_end = exon[3], exon[4]
        if tx_start is None or tx_end is None:
            return None
        out.append((tx_start, tx_end, exon[5] if len(exon) > 5 else None))
    return out


def _pick_build(transcript, genome_build: Optional[str]):
    """Return the requested build's coordinates, else any build's (tx coords are
    build-independent), else ``None``."""
    builds = transcript["genome_builds"]
    if genome_build is not None:
        gb = builds.get(genome_build)
        if gb is not None:
            return gb
    return next(iter(builds.values()), None)
