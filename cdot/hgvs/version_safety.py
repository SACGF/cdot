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

:func:`is_version_substitution_safe` is the check, on the two versions' records;
:func:`intrinsic_cds_structure` is the build-independent structural helper behind
it; :func:`genomic_bracket_safe` is the probabilistic fallback on the two
neighbouring records for the residual case where the requested version exists in
no build at all.  These take transcript records and are pure.
:func:`is_version_substitution_safe_by_accession` is the same check on two
versioned accessions: it fetches the records from a data provider and chooses
between the pure checks, and is what the version fallback in
:mod:`cdot.hgvs.gene_hgvs` calls.

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
sufficient; the safety check (:func:`is_version_substitution_safe`) also requires
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

from cdot.hgvs import _version_replacement_blocklist


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
# Genomic CDS maps (placement check, and the bracket fallback)
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
# The substitution safety check (pure: two transcript records)
# ---------------------------------------------------------------------------

_VERSIONED_ID = re.compile(r"^(.+)\.(\d+)$")


def is_known_replacement(transcript_id: str, version: int, substitute_version: int) -> bool:
    """Is this substitution on the shipped genomic re-placement blocklist?

    The blocklist (``cdot.hgvs._version_replacement_blocklist``) holds substitutions
    with the same intrinsic CDS structure but a different genomic placement, so
    they are refused even when the requested version is absent from the loaded data
    and the placement difference cannot be recomputed.
    """
    key = (transcript_id, version, substitute_version)
    return key in _version_replacement_blocklist.KNOWN_REPLACEMENT_SUBSTITUTIONS


def is_version_substitution_safe(transcript, substitute_transcript, genome_build: Optional[str] = None,
                                 cited_position: Optional[str] = None) -> tuple[bool, str]:
    """Is substituting ``substitute_transcript`` for ``transcript`` coordinate-safe?

    Returns ``(safe, reason)``.  "Safe" means the variant's position maps to the
    *same* genomic coordinate under the substitute, so the substitution does not
    silently change what the variant means.

    ``transcript`` and ``substitute_transcript`` are two versions of one transcript
    as cdot records (plain JSON dicts or :class:`cdot.models.Transcript`), eg
    ``NM_000059.3`` and ``NM_000059.4``.  Both must be present; when the requested
    version is held in no build at all there is nothing to compare, and the
    probabilistic :func:`genomic_bracket_safe` on its neighbours is the fallback
    (:func:`is_version_substitution_safe_by_accession` fetches the records from a
    data provider and chooses between the two).

    The check is structural and near-deterministic: compare the two versions'
    build-independent intrinsic CDS structure.  The requested version's structure
    is read from *any* build that holds it, so the check still decides when the
    requested version is absent from ``genome_build`` but present in another (JSON
    all-builds files and the REST ``/transcript/<ac>`` view both supply this).
    Identical structure (plus matching CDS alignment gaps and, when both versions
    are placed in ``genome_build``, matching genomic placement) is safe for coding
    positions; a CDS-length, coding-exon-structure, gap or placement change is not.
    A pair on the shipped re-placement blocklist is refused outright.

    ``cited_position`` is the text after ``c.``/``n.`` (eg ``"-49="``,
    ``"*100A>G"``).  The intrinsic structure normalises the UTRs away, so for a
    *non-coding* cited position it is not enough: a 5'UTR (``c.-N``) or 3'UTR
    (``c.*N``) coordinate can still move if the relevant UTR length changes
    between versions (very common; UTR annotations churn nearly every version),
    so when the position is in a UTR the matching UTR length must also be equal.
    Coding / CDS-intronic positions (and the default ``cited_position=None``) are
    covered by the CDS-structure check alone.

    Non-coding transcripts (both versions without a CDS, eg ``NR_``) are compared
    on their whole exon structure instead: the same strand and per-exon transcript
    lengths, matching alignment gaps, and (when both are placed in
    ``genome_build``) the same genomic exon spans.  An n. coordinate depends only
    on those, so this is exact.  A version pair where one is coding and the other
    is not is never safe.  A malformed record (an exon missing its transcript
    coordinates) is "cannot verify", never a crash.

    A not-safe verdict is surfaced by
    :func:`cdot.hgvs.gene_hgvs.resolve_transcript_version` (never silently
    applied) per its ``on_unsafe_version`` policy.
    """
    transcript_match = _VERSIONED_ID.match(transcript.get("id") or "")
    substitute_match = _VERSIONED_ID.match(substitute_transcript.get("id") or "")
    if transcript_match and substitute_match and transcript_match.group(1) == substitute_match.group(1):
        if is_known_replacement(transcript_match.group(1), int(transcript_match.group(2)), int(substitute_match.group(2))):
            return False, "known genomic re-placement (blocklisted)"
    # Version labels for the reasons: ".3" from the record id, else the role.
    transcript_label = f".{transcript_match.group(2)}" if transcript_match else "the requested version"
    substitute_label = f".{substitute_match.group(2)}" if substitute_match else "the substitute"

    transcript_coding, substitute_coding = has_cds(transcript), has_cds(substitute_transcript)
    if transcript_coding != substitute_coding:
        coding = transcript_label if transcript_coding else substitute_label
        return False, f"only {coding} is coding (c. and n. coordinates differ)"
    if not transcript_coding:
        # Non-coding (n.): an n. coordinate depends only on the exon lengths in
        # transcript order, the alignment gaps and the genomic exon placement.
        transcript_exons = intrinsic_exon_structure(transcript, genome_build)
        substitute_exons = intrinsic_exon_structure(substitute_transcript, genome_build)
        if transcript_exons is None or substitute_exons is None:
            bad = transcript_label if transcript_exons is None else substitute_label
            return False, f"version {bad} has a malformed exon record; cannot verify"
        if transcript_exons[0] != substitute_exons[0]:
            return False, f"strand changed {transcript_exons[0]}→{substitute_exons[0]}"
        if transcript_exons != substitute_exons:
            return False, "exon structure changed (n. coordinate would shift)"
        if exon_alignment_gaps(transcript, genome_build) != exon_alignment_gaps(substitute_transcript, genome_build):
            return False, "exon alignment gap differs (n. coordinate would shift)"
        if genome_build:
            transcript_spans = exon_genomic_spans(transcript, genome_build)
            substitute_spans = exon_genomic_spans(substitute_transcript, genome_build)
            if transcript_spans is not None and substitute_spans is not None and transcript_spans != substitute_spans:
                return False, "exon genomic placement differs (coordinate would shift)"
        return True, "identical exon structure (non-coding)"

    transcript_struct = intrinsic_cds_structure(transcript, genome_build)
    substitute_struct = intrinsic_cds_structure(substitute_transcript, genome_build)
    if substitute_struct is None:
        # The substitute claims a CDS but it is unusable (malformed exons /
        # zero-length CDS): nothing to compare against.
        return False, f"substitute version {substitute_label} has no usable CDS structure; cannot verify"
    if transcript_struct is None:
        # The requested version is present but its CDS is unusable. There is no
        # structure to compare; is_version_substitution_safe_by_accession brackets
        # it instead.
        return False, f"requested version {transcript_label} has no usable CDS structure; cannot verify"
    if transcript_struct != substitute_struct:
        return False, describe_structure_change(transcript_struct, substitute_struct)
    # Identical transcript structure is necessary but not sufficient: an
    # alignment gap (a transcript-vs-genome indel) present in one version
    # and not the other shifts every coding base downstream of it, so the
    # CDS alignment gaps must match too (gaps are build-specific, so this
    # is exact when both versions are placed in genome_build, conservative
    # across builds).
    if cds_alignment_gaps(transcript, genome_build) != cds_alignment_gaps(substitute_transcript, genome_build):
        return False, "CDS alignment gap differs (coding coordinate would shift)"
    # Intrinsic structure is build-independent and so is blind to genomic
    # *placement*: a version can re-align the same CDS to a different locus
    # (a rare re-placement, eg NM_018263). When both versions are placed in
    # genome_build we can see this directly by comparing their genomic CDS
    # maps - recomputed from the loaded annotation, so no precomputed list and
    # never stale. (Cross-build, where the requested version is absent from
    # genome_build, this cannot be checked and falls back to the build-
    # independent structure.)
    transcript_map = cds_genomic_map(transcript, genome_build) if genome_build else None
    substitute_map = cds_genomic_map(substitute_transcript, genome_build) if genome_build else None
    if transcript_map is not None and substitute_map is not None:
        if (transcript_map["contig"] != substitute_map["contig"]
                or transcript_map["segs"] != substitute_map["segs"]):
            return False, "CDS genomic placement differs (coordinate would shift)"
    # The CDS is coordinate-safe. A non-coding cited position additionally
    # needs the relevant UTR length unchanged, since the CDS check says
    # nothing about where a c.-N / c.*N position falls.
    regions = utr_regions_touched(cited_position)
    if regions:
        transcript_utr = utr_lengths(transcript, genome_build)
        substitute_utr = utr_lengths(substitute_transcript, genome_build)
        if transcript_utr is not None and substitute_utr is not None:
            if "five_prime" in regions and transcript_utr[0] != substitute_utr[0]:
                return False, "5'UTR length differs (UTR coordinate would shift)"
            if "three_prime" in regions and transcript_utr[1] != substitute_utr[1]:
                return False, "3'UTR length differs (UTR coordinate would shift)"
    return True, "identical CDS structure"


def is_version_substitution_safe_by_accession(data_provider, transcript_accession: str,
                                              substitute_transcript_accession: str,
                                              genome_build: Optional[str] = None,
                                              cited_position: Optional[str] = None) -> tuple[bool, str]:
    """:func:`is_version_substitution_safe` on two versioned accessions, fetched from a provider.

    ``transcript_accession`` and ``substitute_transcript_accession`` are two versions
    of one transcript (eg ``"NM_000059.3"`` and ``"NM_000059.4"``).  Their records
    come from ``data_provider._get_transcript``, and the pure check decides.  Only
    when there is no structure to compare (the requested version exists in no build
    the provider can see, or is present but its CDS is unusable) does this fall back
    to the probabilistic :func:`genomic_bracket_safe` on the nearest held versions
    either side, which needs ``genome_build`` and the provider's
    ``get_tx_versions``; without them the verdict is an honest "cannot verify".
    The shipped re-placement blocklist is applied first, so it holds even when the
    requested record is absent.

    This is the one function in the module that touches a data provider; it only
    fetches records and chooses between the pure checks.  Used by
    :func:`cdot.hgvs.gene_hgvs.resolve_transcript_version`.
    """
    transcript_match = _VERSIONED_ID.match(transcript_accession)
    substitute_match = _VERSIONED_ID.match(substitute_transcript_accession)
    if transcript_match is None or substitute_match is None:
        return False, "both accessions must be versioned (eg NM_000059.4); cannot verify"
    transcript_id, version = transcript_match.group(1), int(transcript_match.group(2))
    substitute_transcript_id, substitute_version = (substitute_match.group(1),
                                                    int(substitute_match.group(2)))
    if transcript_id != substitute_transcript_id:
        return False, (f"{transcript_accession} and {substitute_transcript_accession} are "
                       "versions of different transcripts")
    if is_known_replacement(transcript_id, version, substitute_version):
        # Applied here as well as in the pure check so it holds when the requested
        # record is absent from the loaded data.
        return False, "known genomic re-placement (blocklisted)"

    substitute_transcript = data_provider._get_transcript(substitute_transcript_accession)
    if substitute_transcript is None:
        return False, f"substitute version .{substitute_version} not found"

    transcript = data_provider._get_transcript(transcript_accession)
    if transcript is None:
        # Requested version absent from every build the provider can see.
        what = "requested version absent from all builds"
    elif (has_cds(transcript) and has_cds(substitute_transcript)
          and intrinsic_cds_structure(transcript, genome_build) is None
          and intrinsic_cds_structure(substitute_transcript, genome_build) is not None):
        # Requested version present but its CDS is unusable (a malformed record):
        # nothing to compare, so bracket it rather than guess. The reason says the
        # version is present rather than absent.
        what = "requested version has no usable CDS structure"
    else:
        return is_version_substitution_safe(transcript, substitute_transcript, genome_build,
                                            cited_position=cited_position)

    if genome_build is None:
        return False, f"{what}; cannot verify (no genome build for bracket fallback)"
    try:
        versions = data_provider.get_tx_versions(transcript_id)
    except NotImplementedError:
        return False, f"{what}; cannot verify structure"

    def _nearest(candidates):
        # Nearest held version (in candidate order) that is placed in genome_build.
        for v in candidates:
            candidate = data_provider._get_transcript(f"{transcript_id}.{v}")
            if candidate is not None and candidate["genome_builds"].get(genome_build) is not None:
                return candidate
        return None

    transcript_below = _nearest(sorted((v for v in versions if v < version), reverse=True))
    transcript_above = _nearest(sorted(v for v in versions if v > version))
    return genomic_bracket_safe(transcript_below, transcript_above, genome_build, what=what)


def genomic_bracket_safe(transcript_below, transcript_above, genome_build: str, *,
                         what: str = "requested version absent from all builds") -> tuple[bool, str]:
    """Genomic-bracket fallback when the requested version's structure is unknown.

    ``transcript_below`` and ``transcript_above`` are the records of the nearest held
    versions either side of the requested one (``None`` if there is none on that side).  If their
    genomic CDS maps in ``genome_build`` agree at every coding base, the gap
    between them is treated as coordinate-safe (probabilistic: a version can
    change and change back between its neighbours); otherwise not safe.  Either
    side missing, non-coding or not placed in ``genome_build`` is "cannot verify".

    ``what`` says why the structural check could not run, for the reason text:
    normally the requested version is present in no build; it is also used when
    the requested version is present but has no usable CDS structure (a malformed
    record).
    """
    lo_map = cds_genomic_map(transcript_below, genome_build) if transcript_below is not None else None
    hi_map = cds_genomic_map(transcript_above, genome_build) if transcript_above is not None else None
    if lo_map is None or hi_map is None:
        return False, f"{what} and cannot bracket it in {genome_build}; coordinate unverified"
    if genomic_maps_fully_agree(lo_map, hi_map):
        return True, "flanking versions agree (genomic bracket; probabilistic)"
    return False, "flanking versions disagree; coordinate may have moved"


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
