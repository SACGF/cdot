"""
HGVS fixes that need a data provider.

clean_hgvs() (cdot.hgvs.clean) works on the string alone; the fixes here look
things up in a data provider:

- resolve_gene_hgvs() handles HGVS strings that name a gene instead of a
  transcript (e.g. "BRCA2:c.36del"), substituting the best available transcript
  from the provider's tag data (MANE_Select > MANE_Plus_Clinical > RefSeq_Select >
  Ensembl_canonical > longest).
- resolve_missing_accession_prefix() restores a dropped RefSeq prefix when the
  bare number matches exactly one transcript.
- resolve_transcript_version() substitutes an adjacent transcript version when
  the cited one is not held, after cdot.hgvs.version_safety has judged the
  substitution coordinate-safe. It also gives a versionless transcript accession
  (e.g. "NM_000059:c.36del") a version when every held version agrees.

fix_hgvs() is the single combined entry point that chains clean_hgvs() and
these in one call, the recommended function for most callers.
"""

import inspect
import re
import warnings
from enum import Enum
from typing import Optional

from hgvs.exceptions import HGVSDataNotAvailableError

from cdot.hgvs import version_safety
from cdot.hgvs.clean import (
    HGVSFix,
    HGVSFixCode,
    HGVSFixSeverity,
    HGVSInputError,
    VersionStrategy,
    _looks_like_gene_symbol,
    _looks_like_transcript,
    _refseq_prefixes,
    clean_hgvs,
    get_best_transcript_version,
)

DEFAULT_TAG_PRIORITY: list[str] = [
    "MANE_Select",
    "MANE_Plus_Clinical",
    "RefSeq_Select",
    "Ensembl_canonical",
]


class Consortium(Enum):
    """Annotation consortium. Used as a hard filter when resolving a gene symbol
    to a transcript: setting a consortium guarantees the returned accession is
    from that consortium (RefSeq NM_/NR_/XM_/XR_ vs Ensembl ENST), never the
    other.
    """
    REFSEQ = "refseq"
    ENSEMBL = "ensembl"


# RefSeq is the default: clinical/diagnostic HGVS overwhelmingly uses NM_
# transcripts, and biocommons HGVS is typically driven with RefSeq accessions.
DEFAULT_CONSORTIUM: Consortium = Consortium.REFSEQ

# Marks the deprecated prefer_consortium alias as not passed.
_UNSET = object()


class UnsafeVersionPolicy(Enum):
    """What the adjacent-version fallback does when a substitution is NOT
    structure-verified coordinate-safe (the substitute's intrinsic CDS structure
    differs from the requested version's, or it could not be verified).

    REFUSE (default): do not substitute — leave the string unchanged and report an
        ERROR (``REFUSED_UNSAFE_VERSION``), which becomes an ``HGVSInputError`` under
        ``raise_on_errors=True``. Preserves exact-version semantics; a coordinate
        that may have moved is never applied automatically.
    SUBSTITUTE: apply the substitution anyway but report a WARNING
        (``USED_ADJACENT_VERSION_COORD_UNVERIFIED``) that the coordinate may have
        moved, leaving the caller to decide.

    A substitution that IS structure-verified safe is always applied (with a
    ``USED_ADJACENT_VERSION_COORD_SAFE`` WARNING) regardless of this policy.

    The policy does not apply to a versionless accession: when its held versions
    disagree at the cited position it is always refused.
    """
    REFUSE = "refuse"
    SUBSTITUTE = "substitute"


DEFAULT_UNSAFE_VERSION_POLICY: UnsafeVersionPolicy = UnsafeVersionPolicy.REFUSE


def consortium_of(tx_ac: str) -> Consortium:
    """Return the annotation consortium an accession belongs to.

    Ensembl transcripts start with 'ENST'; everything else (NM_/NR_/XM_/XR_/LRG)
    is treated as RefSeq for tie-breaking purposes. Exposed so consumers can
    label/group results by consortium and map cdot's Consortium onto their own enum.
    """
    return Consortium.ENSEMBL if tx_ac.startswith("ENST") else Consortium.REFSEQ


# Backwards-compatible alias for the previous private name.
_consortium_of = consortium_of


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _parse_gene_only_hgvs(hgvs_string: str) -> tuple[Optional[str], Optional[str]]:
    """
    Return (gene_symbol, allele) if hgvs_string is gene-only HGVS (no transcript),
    or (None, None) if it already contains a transcript accession.

    Gene-only HGVS: a c., n., r. or p. allele on a reference that looks like a gene
    symbol (see clean._looks_like_gene_symbol), e.g. "BRCA2:c.36del", "BRCA2:p.Arg10Gly".
    Genomic HGVS ("NC_000013.11:g.32316528del", "BRCA2:g.…") is not gene-only.
    """
    if ":" not in hgvs_string:
        return None, None
    prefix, allele = hgvs_string.split(":", 1)
    if not allele.startswith(("c.", "n.", "r.", "p.")) or not _looks_like_gene_symbol(prefix):
        return None, None
    return prefix, allele


def _normalize_tag(tag: str) -> str:
    """
    Normalise a transcript tag for comparison against tag_priority.

    RefSeq stores tags with spaces ("MANE Select", "RefSeq Select") while
    Ensembl uses underscores ("MANE_Select", "Ensembl_canonical").  Without
    this, RefSeq's MANE/canonical transcripts never match the underscore-form
    entries in DEFAULT_TAG_PRIORITY.
    """
    return tag.replace(" ", "_")


def _filter_by_consortium(
    tx_and_tags: list[tuple[str, list[str]]],
    consortium: Consortium,
) -> list[tuple[str, list[str]]]:
    """Keep only transcripts from the given consortium (hard filter)."""
    return [t for t in tx_and_tags if consortium_of(t[0]) == consortium]


def _annotation_consortium_arg(annotation_consortium, prefer_consortium,
                               stacklevel: int = 3) -> Optional[Consortium]:
    """
    The effective ``annotation_consortium`` once the deprecated ``prefer_consortium``
    alias is folded in: ``prefer_consortium`` maps onto ``annotation_consortium`` with a
    DeprecationWarning, and it is a TypeError to give both (detected when
    ``annotation_consortium`` is not its default).
    """
    if prefer_consortium is _UNSET:
        return annotation_consortium
    if annotation_consortium is not DEFAULT_CONSORTIUM:
        raise TypeError("annotation_consortium and prefer_consortium (its deprecated alias) were "
                        "both given; pass only annotation_consortium")
    warnings.warn("prefer_consortium is deprecated and will be removed in a future release; "
                  "use annotation_consortium (same values, it has always been a hard filter)",
                  DeprecationWarning, stacklevel=stacklevel)
    return prefer_consortium


def _tag_priority_arg(tag_priority: Optional[list[str]]) -> list[str]:
    """``tag_priority``, or a copy of DEFAULT_TAG_PRIORITY when None."""
    return list(DEFAULT_TAG_PRIORITY) if tag_priority is None else tag_priority


def _transcript_length(data_provider, tx_ac: str) -> Optional[int]:
    """
    Transcript (spliced cDNA) length in bases: the sum of the exon lengths from the
    provider's biocommons ``get_tx_identity_info``. None if the provider can't say
    (no such method, transcript not held, or no lengths in the record).
    """
    get_tx_identity_info = getattr(data_provider, "get_tx_identity_info", None)
    if get_tx_identity_info is None:
        return None
    try:
        tx_info = get_tx_identity_info(tx_ac)
        lengths = tx_info["lengths"] if tx_info else None
    except (HGVSDataNotAvailableError, KeyError, TypeError):
        return None
    if not lengths:
        return None
    return sum(lengths)


def _longest_transcript(data_provider, ranked: list[tuple[str, list[str], Optional[str]]]) -> tuple[str, bool]:
    """
    Pick the longest transcript of ``ranked`` (``[(tx_ac, tags, matched_tag), ...]``) by
    :func:`_transcript_length`; ties keep the earlier entry.

    Returns ``(tx_ac, length_verified)``. If any candidate's length can't be read the
    lengths can't be compared, so the first entry (the provider's order) is returned
    with ``length_verified=False``.
    """
    lengths = [(_transcript_length(data_provider, tx_ac), tx_ac) for tx_ac, _tags, _matched in ranked]
    if any(length is None for length, _tx_ac in lengths):
        return ranked[0][0], False
    return max(lengths, key=lambda item: item[0])[1], True


def _rank_transcripts_by_tags(
    tx_and_tags: list[tuple[str, list[str]]],
    tag_priority: list[str],
) -> list[tuple[str, list[str], Optional[str]]]:
    """
    Return [(tx_ac, tags, matched_tag_or_None), ...] sorted best-first.

    Transcripts whose tags include an earlier entry in tag_priority sort first.
    Transcripts with no matching priority tag sort last, in their input order
    (cdot's data providers return them longest-first from get_tx_ac_tags_for_gene,
    but that is not relied on: see :func:`_longest_transcript`).

    Tag comparison is done on normalised tags (see _normalize_tag), so RefSeq
    space-form tags match the underscore-form entries in tag_priority.  The
    original (un-normalised) tags are preserved in the returned tuples.
    """
    def matched_priority_tag(tags: list[str]) -> Optional[str]:
        norm = {_normalize_tag(t) for t in tags}
        return next((tag for tag in tag_priority if tag in norm), None)

    def sort_key(item: tuple[str, list[str]]) -> int:
        _tx_ac, tags = item
        matched = matched_priority_tag(tags)
        return tag_priority.index(matched) if matched else len(tag_priority)

    ranked = sorted(tx_and_tags, key=sort_key)
    result = []
    for tx_ac, tags in ranked:
        result.append((tx_ac, tags, matched_priority_tag(tags)))
    return result


# ---------------------------------------------------------------------------
# Public: rank_transcripts_for_gene
# ---------------------------------------------------------------------------

def rank_transcripts_for_gene(
    gene_symbol: str,
    data_provider,
    genome_build: str,
    tag_priority: Optional[list[str]] = None,
    annotation_consortium: Optional[Consortium] = DEFAULT_CONSORTIUM,
    prefer_consortium=_UNSET,
) -> tuple[list[tuple[str, list[str], Optional[str]]], list[HGVSFix]]:
    """
    Look up the transcripts of a gene symbol and return them ranked best-first.

    This is the gene→transcript *ranking* core, split out of the gene→transcript
    *resolution* done by :func:`resolve_gene_hgvs`. It performs the gene lookup
    (with an uppercase-gene retry after a case-sensitive miss), the consortium
    hard-filter, and the MANE/canonical tag ranking, but does NOT collapse to a
    single transcript or rewrite any HGVS string.

    Returns:
        (ranked, fixes)
        ``ranked`` is ``[(tx_ac, tags, matched_tag_or_None), ...]`` best-first
        (see :func:`_rank_transcripts_by_tags`). It is empty if no transcript
        could be found, in which case ``fixes`` carries an ERROR-level HGVSFix.
        ``fixes`` may also carry a WARNING UPPERCASED_GENE_SYMBOL fix if the
        gene symbol had to be uppercased to find a match.

    A consumer wanting the single best transcript takes ``ranked[0]``; a consumer
    running a multi-transcript search uses the whole ordered list. Either way it
    avoids re-implementing the lookup/uppercase/consortium/tag logic.

    tag_priority and annotation_consortium behave as documented on
    :func:`resolve_gene_hgvs`. prefer_consortium is a deprecated alias of
    annotation_consortium.
    """
    annotation_consortium = _annotation_consortium_arg(annotation_consortium, prefer_consortium)
    tag_priority = _tag_priority_arg(tag_priority)
    fixes: list[HGVSFix] = []

    tx_and_tags = data_provider.get_tx_ac_tags_for_gene(gene_symbol, genome_build)

    # Gene lookups are case-sensitive; gene symbols are conventionally uppercase
    # (some legitimately have lowercase parts, e.g. "C7orf26", so we only retry
    # uppercased after an exact-case miss rather than uppercasing unconditionally).
    if not tx_and_tags and gene_symbol != gene_symbol.upper():
        uc_gene_symbol = gene_symbol.upper()
        tx_and_tags = data_provider.get_tx_ac_tags_for_gene(uc_gene_symbol, genome_build)
        if tx_and_tags:
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=HGVSFixCode.UPPERCASED_GENE_SYMBOL,
                message=f"Uppercased gene symbol '{gene_symbol}' to '{uc_gene_symbol}'",
                original=gene_symbol,
                fixed=uc_gene_symbol,
            ))
            gene_symbol = uc_gene_symbol

    if not tx_and_tags:
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.ERROR,
            code=HGVSFixCode.NO_TRANSCRIPT_FOR_GENE,
            message=f"No transcripts found for gene '{gene_symbol}' in {genome_build}",
        ))
        return [], fixes

    # Consortium is a hard filter: with one set we never return a transcript
    # from the other consortium. If the chosen consortium has none (but the
    # other does), error with an actionable message rather than crossing over.
    if annotation_consortium is not None:
        filtered = _filter_by_consortium(tx_and_tags, annotation_consortium)
        if not filtered:
            other = (Consortium.ENSEMBL if annotation_consortium == Consortium.REFSEQ
                     else Consortium.REFSEQ)
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.ERROR,
                code=HGVSFixCode.NO_TRANSCRIPT_FOR_GENE,
                message=(
                    f"No {annotation_consortium.value} transcript found for gene "
                    f"'{gene_symbol}' in {genome_build} ({len(tx_and_tags)} "
                    f"{other.value} transcript(s) available). Pass "
                    f"annotation_consortium=Consortium.{other.name} or None to use them."
                ),
            ))
            return [], fixes
        tx_and_tags = filtered

    ranked = _rank_transcripts_by_tags(tx_and_tags, tag_priority)
    return ranked, fixes


# ---------------------------------------------------------------------------
# Public: resolve_gene_hgvs
# ---------------------------------------------------------------------------

def resolve_gene_hgvs(
    hgvs_string: str,
    data_provider,
    genome_build: str,
    tag_priority: Optional[list[str]] = None,
    fallback_to_longest: bool = False,
    annotation_consortium: Optional[Consortium] = DEFAULT_CONSORTIUM,
    prefer_consortium=_UNSET,
) -> tuple[str, list[HGVSFix]]:
    """
    If hgvs_string contains a gene symbol but no transcript (e.g. "BRCA2:c.36del"),
    look up the best available transcript from data_provider and substitute it. For
    a protein change ("BRCA2:p.Arg10Gly") the chosen transcript's protein is used.

    This is a thin wrapper over :func:`rank_transcripts_for_gene`: it takes the
    top-ranked transcript, applies the fallback_to_longest policy, and rewrites
    the HGVS string. Consumers that need the full candidate list (e.g. a
    multi-transcript search) should call rank_transcripts_for_gene directly.

    Returns:
        (resolved_string, fixes)
        fixes contains a WARNING-level HGVSFix describing the substitution,
        or an ERROR-level fix if no suitable transcript could be found.

    If hgvs_string already has a transcript accession, it is returned unchanged
    with an empty fixes list (this function is non-destructive).

    Transcript selection priority (configurable via tag_priority, which defaults
    to DEFAULT_TAG_PRIORITY):
        MANE_Select > MANE_Plus_Clinical > RefSeq_Select > Ensembl_canonical

    If no transcript carries a priority tag and fallback_to_longest=False (default),
    an ERROR-level fix is returned. Pass fallback_to_longest=True to use the longest
    available transcript instead, with a WARNING. Longest means the spliced
    transcript length (the sum of the exon lengths from the provider's
    ``get_tx_identity_info``), compared explicitly across the candidates. If the
    provider can't give a length (eg it has no ``get_tx_identity_info``), the first
    transcript the provider returned is used and the WARNING says so.

    annotation_consortium is a HARD filter on the returned transcript's consortium:
    with Consortium.REFSEQ (the default) you never get an Ensembl transcript back,
    and vice versa. If the chosen consortium has no transcript for the gene (but the
    other does), an ERROR fix is returned rather than crossing over. Pass None to
    consider both consortiums (ranked by tag, then input order).
    prefer_consortium is a deprecated alias of annotation_consortium
    (DeprecationWarning; giving both is a TypeError).

    Example::

        resolved, fixes = resolve_gene_hgvs("BRCA2:c.36del", data_provider, "GRCh38")
        # resolved = "NM_000059.4:c.36del"
        # fixes[0].message = "Resolved gene symbol 'BRCA2' to NM_000059.4 (MANE Select)"
    """
    annotation_consortium = _annotation_consortium_arg(annotation_consortium, prefer_consortium)
    gene_symbol, allele = _parse_gene_only_hgvs(hgvs_string)
    if gene_symbol is None:
        return hgvs_string, []  # already has transcript, nothing to do

    ranked, fixes = rank_transcripts_for_gene(
        gene_symbol, data_provider, genome_build,
        tag_priority=tag_priority,
        annotation_consortium=annotation_consortium,
    )
    if not ranked:
        return hgvs_string, fixes  # lookup failed, fixes carries the ERROR

    # The uppercase-gene retry may have changed the effective symbol; reflect it
    # in the resolution message.
    for f in fixes:
        if f.code == HGVSFixCode.UPPERCASED_GENE_SYMBOL:
            gene_symbol = f.fixed

    best_tx, _best_tags, matched_tag = ranked[0]

    if matched_tag is None and not fallback_to_longest:
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.ERROR,
            code=HGVSFixCode.NO_TRANSCRIPT_FOR_GENE,
            message=(
                f"No MANE or canonical transcript found for gene '{gene_symbol}' "
                f"in {genome_build}. Pass fallback_to_longest=True to use the "
                f"longest available transcript."
            ),
        ))
        return hgvs_string, fixes

    if matched_tag is not None:
        tag_desc = f" ({matched_tag})"
    else:
        # ranked is sorted tag-first, so if the top entry has no priority tag none
        # does: pick the longest of them by length rather than trusting the
        # provider's order.
        best_tx, length_verified = _longest_transcript(data_provider, ranked)
        tag_desc = (" (longest transcript, no canonical tag)" if length_verified
                    else " (first transcript returned by the data provider, no canonical tag)")
    if allele.startswith("p."):
        # A protein change goes on the chosen transcript's protein
        if not (protein := data_provider.get_pro_ac_for_tx_ac(best_tx)):
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.ERROR,
                code=HGVSFixCode.NO_TRANSCRIPT_FOR_GENE,
                message=(f"Gene symbol '{gene_symbol}' resolved to {best_tx}{tag_desc}, which has no "
                         f"protein, so a p. change can't be placed on it"),
            ))
            return hgvs_string, fixes
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.RESOLVED_GENE_TO_PROTEIN,
            message=f"Resolved gene symbol '{gene_symbol}' to {protein}, the protein of {best_tx}{tag_desc}",
            original=gene_symbol,
            fixed=protein,
        ))
        return f"{protein}:{allele}", fixes

    resolved = f"{best_tx}:{allele}"
    fixes.append(HGVSFix(
        severity=HGVSFixSeverity.WARNING,
        code=HGVSFixCode.RESOLVED_GENE_TO_TRANSCRIPT,
        message=f"Resolved gene symbol '{gene_symbol}' to {best_tx}{tag_desc}",
        original=gene_symbol,
        fixed=best_tx,
    ))
    return resolved, fixes


# ---------------------------------------------------------------------------
# Public: resolve_missing_accession_prefix — restore a dropped RefSeq prefix
# ---------------------------------------------------------------------------

# A bare-number accession: the user dropped the RefSeq prefix entirely (or left
# only its underscore), e.g. "000059.4:c.68del" / "_001128425.1:c.667A>G".
# RefSeq digit fields are 6 or 9 digits, but 7-8 are accepted as a 9-digit
# accession missing its leading zeros (see _digit_field_candidates).
_BARE_NUMBER_ACCESSION = re.compile(r"^_?(\d{6,9})((?:\.\d+)?):([cn])\.", re.IGNORECASE)

# The kind letter constrains which RefSeq prefixes could have been dropped:
# a coding "c." implies an mRNA accession, a non-coding "n." an RNA accession.
_PREFIXES_FOR_KIND = {
    "c": tuple(_refseq_prefixes("mRNA")),
    "n": tuple(_refseq_prefixes("RNA")),
}


def _digit_field_candidates(digits: str) -> list[str]:
    """Candidate digit fields for a bare-number accession.

    RefSeq accessions have 6 or 9 digits after the prefix. A 6 or 9 digit
    input is taken as written; a 7-8 digit input can only be a 9-digit
    accession whose leading zeros were also dropped, so it is zero-padded.
    """
    if len(digits) in (6, 9):
        return [digits]
    return [digits.zfill(9)]


def resolve_missing_accession_prefix(
    hgvs_string: str,
    data_provider,
    genome_build: Optional[str] = None,
) -> tuple[str, list[HGVSFix]]:
    """
    Restore a RefSeq prefix that was dropped from the transcript accession,
    e.g. "000059.4:c.68del" → "NM_000059.4:c.68del".

    The kind letter narrows the possibilities ("c." implies NM_/XM_, "n."
    implies NR_/XR_) but cannot fully disambiguate, so the candidates are
    checked against the data provider (via ``get_tx_versions``) and the fix is
    applied only when exactly one candidate accession exists there. With zero
    or several matches the string is returned unchanged - this function never
    guesses.

    Returns:
        (resolved_string, fixes)
        fixes contains a WARNING ADDED_ACCESSION_PREFIX HGVSFix when the prefix
        was restored, and is empty otherwise (no bare-number accession, no
        provider match, ambiguous match, or a provider that cannot enumerate
        versions).

    ``genome_build`` is optional and only used to restrict the existence check
    to accessions placeable in that build (for providers whose
    ``get_tx_versions`` accepts it).

    Example::

        resolved, fixes = resolve_missing_accession_prefix(
            "000059.4:c.68del", data_provider)
        # resolved = "NM_000059.4:c.68del"   (if only NM_000059 exists)
        # fixes[0].code = HGVSFixCode.ADDED_ACCESSION_PREFIX
    """
    m = _BARE_NUMBER_ACCESSION.match(hgvs_string)
    if m is None:
        return hgvs_string, []

    digits, version, kind = m.group(1), m.group(2), m.group(3).lower()
    candidates = []
    for prefix in _PREFIXES_FOR_KIND[kind]:
        for digit_field in _digit_field_candidates(digits):
            accession = f"{prefix}{digit_field}"
            try:
                versions = _get_tx_versions(data_provider, accession, genome_build)
            except NotImplementedError:
                return hgvs_string, []  # provider can't enumerate, don't guess
            if versions:
                candidates.append(accession)

    if len(candidates) != 1:
        return hgvs_string, []  # unknown or ambiguous, never guess

    accession = candidates[0]
    original_head = hgvs_string[:m.end(2)]
    fixed_head = f"{accession}{version}"
    resolved = fixed_head + hgvs_string[m.end(2):]
    return resolved, [HGVSFix(
        severity=HGVSFixSeverity.WARNING,
        code=HGVSFixCode.ADDED_ACCESSION_PREFIX,
        message=(
            f"Restored missing accession prefix: '{original_head}' → "
            f"'{fixed_head}' (unique match in data provider)"
        ),
        original=original_head,
        fixed=fixed_head,
    )]


# ---------------------------------------------------------------------------
# Public: resolve_transcript_version — adjacent-version fallback
# ---------------------------------------------------------------------------

def _parse_transcript_accession(hgvs_string: str) -> Optional[tuple[str, Optional[int]]]:
    """
    Return (versionless_accession, version) if hgvs_string starts with a transcript
    accession (e.g. "NM_000059.4:c.36del" → ("NM_000059", 4)), with version None when
    it has none ("NM_000059:c.36del" → ("NM_000059", None)), else None.

    Handles a gene-in-parens suffix ("NM_000059.4(BRCA2):c.36del"). Returns None for
    gene-only HGVS, genomic refs and LRG accessions (which have no versions).
    """
    if ":" not in hgvs_string:
        return None
    prefix, _allele = hgvs_string.split(":", 1)
    accession = prefix.split("(", 1)[0]  # drop "(GENE)" suffix if present
    if not _looks_like_transcript(accession) or accession.startswith("LRG_"):
        return None
    base, dot, version = accession.rpartition(".")
    if not dot:
        return accession, None
    if not version.isdigit():
        return None
    return base, int(version)


def _resolve_versionless_transcript(hgvs_string: str, accession: str, versions: list[int],
                                    data_provider, genome_build: Optional[str]) -> tuple[str, list[HGVSFix]]:
    """
    Give a versionless accession (eg "NM_000059:c.36del") its latest held version, if
    every held version gives the same coordinate at the cited position. Otherwise
    refuse: with no version given there is no nearest version to fall back to, so
    picking one would be a guess between different variants (``on_unsafe_version``
    does not apply).
    """
    if not versions:
        return hgvs_string, [HGVSFix(
            severity=HGVSFixSeverity.ERROR,
            code=HGVSFixCode.NO_TRANSCRIPT_VERSIONS,
            message=f"No versions of {accession} are available in the data provider",
        )]
    latest = f"{accession}.{max(versions)}"
    cited_position = _cited_position(hgvs_string)
    disagree = []
    for version in sorted(set(versions) - {max(versions)}):
        if hasattr(data_provider, "_get_transcript"):
            safe, reason = version_safety.is_version_substitution_safe_by_accession(
                data_provider, f"{accession}.{version}", latest, genome_build,
                cited_position=cited_position)
        else:
            safe, reason = False, "no transcript records to compare"
        if not safe:
            disagree.append(f".{version} ({reason})")

    if disagree:
        return hgvs_string, [HGVSFix(
            severity=HGVSFixSeverity.ERROR,
            code=HGVSFixCode.REFUSED_UNSAFE_VERSION,
            message=(f"Transcript {accession} has no version, and these versions may give a "
                     f"different coordinate to {latest}: {', '.join(disagree)}. "
                     f"Give the transcript version."),
            original=accession,
            fixed=latest,
        )]
    return latest + hgvs_string[len(accession):], [HGVSFix(
        severity=HGVSFixSeverity.WARNING,
        code=HGVSFixCode.RESOLVED_VERSIONLESS_TRANSCRIPT,
        message=(f"Transcript {accession} has no version; using {latest}, as every "
                 f"version held gives the same coordinate"),
        original=accession,
        fixed=latest,
    )]


def _check_kind_fits_transcript(hgvs_string: str, data_provider) -> list[HGVSFix]:
    """ An ERROR KIND_REFERENCE_MISMATCH if hgvs_string is c. on a non-coding transcript or n. on
        a coding one, going by the transcript record (the accession alone doesn't say for Ensembl,
        or for a transcript that gene resolution picked). n. counts from the transcript start and
        c. from the start codon, so the same number on the other one is a different base. """
    parsed = _parse_transcript_accession(hgvs_string)
    if parsed is None or parsed[1] is None or not hasattr(data_provider, "get_transcript"):
        return []
    tx_ac = f"{parsed[0]}.{parsed[1]}"
    kind = hgvs_string.split(":", 1)[1][:2]
    if kind not in ("c.", "n.") or (transcript := data_provider.get_transcript(tx_ac)) is None:
        return []
    coding = version_safety.has_cds(transcript)
    if coding and kind == "n.":
        message = (f"{tx_ac} is a coding transcript, so it can't take 'n.' (use c.): "
                   f"n. counts from the transcript start, not the start codon")
    elif not coding and kind == "c.":
        message = f"{tx_ac} is a non-coding transcript, so it can't take 'c.' (use n.): it has no start codon"
    else:
        return []
    return [HGVSFix(severity=HGVSFixSeverity.ERROR, code=HGVSFixCode.KIND_REFERENCE_MISMATCH,
                    message=message)]


def _cited_position(hgvs_string: str) -> Optional[str]:
    """Return the text after ``c.``/``n.`` (position plus edit) of a c./n. HGVS, else None.

    Used to make the version-substitution safety check position-aware: a 5'UTR
    (``c.-N``) or 3'UTR (``c.*N``) cited position needs the matching UTR length
    preserved, not just the CDS structure (see
    :func:`cdot.hgvs.version_safety.is_version_substitution_safe`).
    """
    if ":" not in hgvs_string:
        return None
    _prefix, allele = hgvs_string.split(":", 1)
    for kind in ("c.", "n."):
        if allele.startswith(kind):
            return allele[len(kind):]
    return None


def _get_tx_versions(data_provider, versionless: str, genome_build: Optional[str]) -> list[int]:
    """Call ``get_tx_versions``, passing ``genome_build`` only if the provider accepts it.

    The build kwarg was added later (so the fallback substitutes among versions
    placeable in the target build, not all builds); older / external providers
    implement ``get_tx_versions(accession)`` only, so we introspect rather than
    assume the wider signature.
    """
    if genome_build is not None:
        try:
            params = inspect.signature(data_provider.get_tx_versions).parameters
        except (ValueError, TypeError):
            params = {}
        if "genome_build" in params:
            return data_provider.get_tx_versions(versionless, genome_build=genome_build)
    return data_provider.get_tx_versions(versionless)


def resolve_transcript_version(
    hgvs_string: str,
    data_provider,
    strategy: VersionStrategy = VersionStrategy.UP_THEN_DOWN,
    genome_build: Optional[str] = None,
    on_unsafe_version: UnsafeVersionPolicy = DEFAULT_UNSAFE_VERSION_POLICY,
) -> tuple[str, list[HGVSFix]]:
    """
    If hgvs_string names a transcript version that the data provider doesn't have,
    substitute the best available adjacent version (per ``strategy``) and rewrite
    the string.

    Returns:
        (resolved_string, fixes)
        fixes is empty if the requested version exists (non-destructive) or if the
        string carries no transcript accession. When a substitution is made the fix
        depends on whether it is coordinate-safe (see below); an ERROR
        NO_TRANSCRIPT_VERSIONS fix (string unchanged) if the provider has no version
        of the accession at all.

    The data provider must implement ``get_tx_versions(accession)`` (JSONDataProvider
    and RESTDataProvider do). Providers that can't enumerate versions raise
    NotImplementedError, which is surfaced as an ERROR fix.

    Coordinate-safety (#28): if the provider can return transcript records
    (``_get_transcript``, as JSONDataProvider / RESTDataProvider do), the substitution
    is checked by :func:`cdot.hgvs.version_safety.is_version_substitution_safe_by_accession`
    (via the build-independent intrinsic CDS structure) before it is applied:

      * structure-verified safe → substitute, WARNING USED_ADJACENT_VERSION_COORD_SAFE;
      * not safe and on_unsafe_version=SUBSTITUTE → substitute anyway, WARNING
        USED_ADJACENT_VERSION_COORD_UNVERIFIED ("coordinate may have moved");
      * not safe and on_unsafe_version=REFUSE (default) → string unchanged, ERROR
        REFUSED_UNSAFE_VERSION.

    ``genome_build`` is only needed for the genomic-bracket fallback inside the safety
    check (the residual case where the requested version is absent from every build);
    the structural check itself is build-independent. Providers without the safety
    check fall back to a plain WARNING USED_ADJACENT_VERSION (unchanged behaviour).

    A transcript accession with no version (``"NM_000059:c.36del"``) gets its latest
    held version (WARNING RESOLVED_VERSIONLESS_TRANSCRIPT) only if every held version
    gives the same coordinate at the cited position; otherwise it is left unchanged
    with an ERROR REFUSED_UNSAFE_VERSION, whatever on_unsafe_version says.

    Example::

        # data provider has NM_000059 versions [3, 4] but not the requested .2
        resolved, fixes = resolve_transcript_version("NM_000059.2:c.36del", dp)
        # resolved = "NM_000059.4:c.36del"  (if .4 is structure-verified safe)
        # fixes[0].code = USED_ADJACENT_VERSION_COORD_SAFE
    """
    parsed = _parse_transcript_accession(hgvs_string)
    if parsed is None:
        return hgvs_string, []  # no transcript accession, nothing to do
    versionless, requested_version = parsed

    try:
        available_versions = _get_tx_versions(data_provider, versionless, genome_build)
    except NotImplementedError as e:
        return hgvs_string, [HGVSFix(
            severity=HGVSFixSeverity.ERROR,
            code=HGVSFixCode.NO_TRANSCRIPT_VERSIONS,
            message=str(e),
        )]

    if requested_version is None:
        return _resolve_versionless_transcript(hgvs_string, versionless, available_versions,
                                               data_provider, genome_build)

    try:
        best, fix = get_best_transcript_version(
            versionless, requested_version, available_versions, strategy=strategy,
        )
    except HGVSInputError as e:
        return hgvs_string, [HGVSFix(
            severity=HGVSFixSeverity.ERROR,
            code=HGVSFixCode.NO_TRANSCRIPT_VERSIONS,
            message=str(e),
        )]

    if fix is None:
        return hgvs_string, []  # requested version exists — non-destructive

    resolved = hgvs_string.replace(
        f"{versionless}.{requested_version}", f"{versionless}.{best}", 1
    )

    if not hasattr(data_provider, "_get_transcript"):
        # Provider can't supply transcript records, so coordinate-safety can't be
        # assessed; keep the plain substitution fix.
        return resolved, [fix]

    safe, reason = version_safety.is_version_substitution_safe_by_accession(
        data_provider, f"{versionless}.{requested_version}", f"{versionless}.{best}",
        genome_build, cited_position=_cited_position(hgvs_string))
    moved = (f"Transcript version {versionless}.{requested_version} not found; "
             f"using .{best} instead")
    if safe:
        return resolved, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.USED_ADJACENT_VERSION_COORD_SAFE,
            message=f"{moved} (coordinate-safe: {reason})",
            original=f"{versionless}.{requested_version}",
            fixed=f"{versionless}.{best}",
        )]
    if on_unsafe_version == UnsafeVersionPolicy.SUBSTITUTE:
        return resolved, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.USED_ADJACENT_VERSION_COORD_UNVERIFIED,
            message=f"{moved}; coordinate may have moved ({reason})",
            original=f"{versionless}.{requested_version}",
            fixed=f"{versionless}.{best}",
        )]
    # REFUSE (default): do not substitute; leave the string unchanged.
    return hgvs_string, [HGVSFix(
        severity=HGVSFixSeverity.ERROR,
        code=HGVSFixCode.REFUSED_UNSAFE_VERSION,
        message=(
            f"Transcript version {versionless}.{requested_version} not found and the "
            f"nearest available version .{best} is not coordinate-safe ({reason}); "
            f"not substituting. Pass on_unsafe_version=UnsafeVersionPolicy.SUBSTITUTE "
            f"to override."
        ),
        original=f"{versionless}.{requested_version}",
        fixed=f"{versionless}.{best}",
    )]


# ---------------------------------------------------------------------------
# Public: fix_hgvs — combined entry point
# ---------------------------------------------------------------------------

def fix_hgvs(
    hgvs_string: str,
    data_provider=None,
    genome_build: Optional[str] = None,
    tag_priority: Optional[list[str]] = None,
    fallback_to_longest: bool = False,
    raise_on_errors: bool = False,
    ops: Optional[set] = None,
    annotation_consortium: Optional[Consortium] = DEFAULT_CONSORTIUM,
    version_fallback: Optional[VersionStrategy] = None,
    on_unsafe_version: UnsafeVersionPolicy = DEFAULT_UNSAFE_VERSION_POLICY,
    prefer_consortium=_UNSET,
) -> tuple[str, list[HGVSFix]]:
    """
    Clean and resolve an HGVS string in one call.

    Runs clean_hgvs() first (always), then resolve_missing_accession_prefix()
    if a data_provider is supplied, then resolve_gene_hgvs() if a data_provider
    and genome_build are supplied.  All fixes from every step are returned in a
    single list.

    If data_provider/genome_build are omitted, only string cleaning is performed.
    This is useful when the caller knows the input already contains a transcript.
    With a data_provider but no genome_build, a gene-only input (eg "BRCA2:c.36del")
    can't be resolved to a transcript, and a WARNING SKIPPED_GENE_RESOLUTION fix
    says so; an input that already names a transcript gets no such fix.

    If raise_on_errors=True, raises HGVSInputError on the first ERROR-level fix
    instead of returning it in the list.

    ops selects the subset of clean_hgvs cleaning operations to apply (see
    clean_hgvs); None (default) runs them all.

    annotation_consortium hard-filters the resolved transcript to RefSeq (default)
    or Ensembl; pass None to allow both (see resolve_gene_hgvs). prefer_consortium
    is a deprecated alias of annotation_consortium.

    version_fallback (default None = off) opts in to adjacent transcript-version
    fallback: if set to a VersionStrategy and a data_provider is supplied, an
    HGVS string whose exact transcript version isn't in the data is rewritten to
    the best available version (see resolve_transcript_version). If the provider can
    assess coordinate-safety, a structure-verified-safe substitution is applied with
    a WARNING; otherwise on_unsafe_version decides (REFUSE by default → ERROR, string
    unchanged; SUBSTITUTE → apply with a "coordinate may have moved" WARNING). A
    genome_build, if supplied, is used only for the genomic-bracket fallback inside
    that check. A versionless accession (eg "NM_000059:c.36del") gets a version only
    when every held version agrees at the cited position.

    Example — gene-only input resolved via MANE::

        result, fixes = fix_hgvs("BRCA2:c.36DEL", data_provider, "GRCh38")
        # result  = "NM_000059.4:c.36del"
        # fixes   = [HGVSFix(WARNING, LOWERCASED_MUTATION_TYPE, ...),
        #            HGVSFix(WARNING, RESOLVED_GENE_TO_TRANSCRIPT, ...)]

    Example — cleaning only (no data provider)::

        result, fixes = fix_hgvs("NM_000059.4 c.316+5G>A")
        # result = "NM_000059.4:c.316+5G>A"

    Example — opt-in version fallback::

        result, fixes = fix_hgvs("NM_000059.2:c.36del", data_provider,
                                 version_fallback=VersionStrategy.UP_THEN_DOWN)
        # result = "NM_000059.4:c.36del"  (if .2 absent but .4 present)
    """
    annotation_consortium = _annotation_consortium_arg(annotation_consortium, prefer_consortium)
    result, fixes = clean_hgvs(hgvs_string, ops=ops)

    if data_provider is not None:
        # Restore a dropped RefSeq prefix (eg "000059.4:c.68del") before gene
        # resolution, so the bare number is not mistaken for a gene symbol.
        result, prefix_fixes = resolve_missing_accession_prefix(
            result, data_provider, genome_build=genome_build,
        )
        fixes.extend(prefix_fixes)

    if data_provider is not None and genome_build is not None:
        result, resolution_fixes = resolve_gene_hgvs(
            result, data_provider, genome_build,
            tag_priority=tag_priority,
            fallback_to_longest=fallback_to_longest,
            annotation_consortium=annotation_consortium,
        )
        fixes.extend(resolution_fixes)
    elif data_provider is not None:
        # Gene resolution needs a genome build; say so rather than silently
        # leaving a gene-only string as it is.
        gene_symbol, _allele = _parse_gene_only_hgvs(result)
        if gene_symbol is not None:
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=HGVSFixCode.SKIPPED_GENE_RESOLUTION,
                message=(
                    f"Gene symbol '{gene_symbol}' was not resolved to a transcript: "
                    f"fix_hgvs needs a genome_build for gene resolution"
                ),
                original=gene_symbol,
            ))

    if version_fallback is not None and data_provider is not None:
        result, version_fixes = resolve_transcript_version(
            result, data_provider, strategy=version_fallback,
            genome_build=genome_build, on_unsafe_version=on_unsafe_version,
        )
        fixes.extend(version_fixes)

    if data_provider is not None and HGVSFixCode.KIND_REFERENCE_MISMATCH not in {f.code for f in fixes}:
        fixes.extend(_check_kind_fits_transcript(result, data_provider))

    if raise_on_errors:
        errors = [f for f in fixes if f.severity == HGVSFixSeverity.ERROR]
        if errors:
            raise HGVSInputError(errors[0].message)

    return result, fixes
