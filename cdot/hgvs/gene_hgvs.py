"""
Gene-symbol HGVS resolution via MANE/canonical transcript tags.

resolve_gene_hgvs() handles the common case of HGVS strings that name a gene
instead of a transcript (e.g. "BRCA2:c.36del").  It looks up the best available
transcript from the data provider's tag data (MANE_Select > MANE_Plus_Clinical >
RefSeq_Select > Ensembl_canonical > longest) and substitutes it.

fix_hgvs() is the single combined entry point that chains clean_hgvs() and
resolve_gene_hgvs() in one call — the recommended function for most callers.
"""

from enum import Enum
from typing import Optional

from cdot.hgvs.clean import (
    HGVSFix,
    HGVSFixCode,
    HGVSFixSeverity,
    HGVSInputError,
    VersionStrategy,
    _looks_like_transcript,
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


def _consortium_of(tx_ac: str) -> Consortium:
    """Ensembl transcripts start with 'ENST'; everything else (NM_/NR_/XM_/XR_/LRG)
    is treated as RefSeq for tie-breaking purposes."""
    return Consortium.ENSEMBL if tx_ac.startswith("ENST") else Consortium.REFSEQ


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _parse_gene_only_hgvs(hgvs_string: str) -> tuple[Optional[str], Optional[str]]:
    """
    Return (gene_symbol, allele) if hgvs_string is gene-only HGVS (no transcript),
    or (None, None) if it already contains a transcript accession.

    Gene-only HGVS: the part before ':' does not look like a transcript accession
    and contains no parentheses (e.g. "BRCA2:c.36del").
    """
    if ":" not in hgvs_string:
        return None, None
    prefix, allele = hgvs_string.split(":", 1)
    # Parentheses mean a transcript(gene) pair — already has a transcript
    if "(" in prefix or ")" in prefix:
        return None, None
    if _looks_like_transcript(prefix):
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
    return [t for t in tx_and_tags if _consortium_of(t[0]) == consortium]


def _rank_transcripts_by_tags(
    tx_and_tags: list[tuple[str, list[str]]],
    tag_priority: list[str],
) -> list[tuple[str, list[str], Optional[str]]]:
    """
    Return [(tx_ac, tags, matched_tag_or_None), ...] sorted best-first.

    Transcripts whose tags include an earlier entry in tag_priority sort first.
    Transcripts with no matching priority tag sort last (preserving the
    input order, which is longest-first from get_tx_ac_tags_for_gene).

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
# Public: resolve_gene_hgvs
# ---------------------------------------------------------------------------

def resolve_gene_hgvs(
    hgvs_string: str,
    data_provider,
    genome_build: str,
    tag_priority: list[str] = DEFAULT_TAG_PRIORITY,
    fallback_to_longest: bool = False,
    prefer_consortium: Optional[Consortium] = DEFAULT_CONSORTIUM,
) -> tuple[str, list[HGVSFix]]:
    """
    If hgvs_string contains a gene symbol but no transcript (e.g. "BRCA2:c.36del"),
    look up the best available transcript from data_provider and substitute it.

    Returns:
        (resolved_string, fixes)
        fixes contains a WARNING-level HGVSFix describing the substitution,
        or an ERROR-level fix if no suitable transcript could be found.

    If hgvs_string already has a transcript accession, it is returned unchanged
    with an empty fixes list — this function is non-destructive.

    Transcript selection priority (configurable via tag_priority):
        MANE_Select > MANE_Plus_Clinical > RefSeq_Select > Ensembl_canonical

    If no transcript carries a priority tag and fallback_to_longest=False (default),
    an ERROR-level fix is returned. Pass fallback_to_longest=True to use the longest
    available transcript instead, with a WARNING.

    prefer_consortium is a HARD filter on the returned transcript's consortium:
    with Consortium.REFSEQ (the default) you never get an Ensembl transcript back,
    and vice versa. If the preferred consortium has no transcript for the gene
    (but the other does), an ERROR fix is returned rather than crossing over.
    Pass None to consider both consortiums (ranked by tag, then input order).

    Example::

        resolved, fixes = resolve_gene_hgvs("BRCA2:c.36del", data_provider, "GRCh38")
        # resolved = "NM_000059.4:c.36del"
        # fixes[0].message = "Resolved gene symbol 'BRCA2' to NM_000059.4 (MANE Select)"
    """
    fixes: list[HGVSFix] = []

    gene_symbol, allele = _parse_gene_only_hgvs(hgvs_string)
    if gene_symbol is None:
        return hgvs_string, fixes  # already has transcript — nothing to do

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
        return hgvs_string, fixes

    # Consortium is a hard filter: with a preference set we never return a
    # transcript from the other consortium. If the preferred consortium has
    # none (but the other does), error with an actionable message rather than
    # crossing over.
    if prefer_consortium is not None:
        filtered = _filter_by_consortium(tx_and_tags, prefer_consortium)
        if not filtered:
            other = (Consortium.ENSEMBL if prefer_consortium == Consortium.REFSEQ
                     else Consortium.REFSEQ)
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.ERROR,
                code=HGVSFixCode.NO_TRANSCRIPT_FOR_GENE,
                message=(
                    f"No {prefer_consortium.value} transcript found for gene "
                    f"'{gene_symbol}' in {genome_build} ({len(tx_and_tags)} "
                    f"{other.value} transcript(s) available). Pass "
                    f"prefer_consortium=Consortium.{other.name} or None to use them."
                ),
            ))
            return hgvs_string, fixes
        tx_and_tags = filtered

    ranked = _rank_transcripts_by_tags(tx_and_tags, tag_priority)
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

    tag_desc = f" ({matched_tag})" if matched_tag else " (longest transcript, no canonical tag)"
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
# Public: resolve_transcript_version — adjacent-version fallback
# ---------------------------------------------------------------------------

def _parse_versioned_transcript(hgvs_string: str) -> Optional[tuple[str, int]]:
    """
    Return (versionless_accession, version) if hgvs_string carries a versioned
    transcript accession (e.g. "NM_000059.4:c.36del" → ("NM_000059", 4)), else None.

    Handles a gene-in-parens suffix ("NM_000059.4(BRCA2):c.36del"). Returns None for
    gene-only HGVS, genomic refs, or transcripts with no version (nothing to fall back from).
    """
    if ":" not in hgvs_string:
        return None
    prefix, _allele = hgvs_string.split(":", 1)
    accession = prefix.split("(", 1)[0]  # drop "(GENE)" suffix if present
    if not _looks_like_transcript(accession):
        return None
    base, dot, version = accession.rpartition(".")
    if not dot or not version.isdigit():
        return None  # no version present — nothing to fall back from
    return base, int(version)


def resolve_transcript_version(
    hgvs_string: str,
    data_provider,
    strategy: VersionStrategy = VersionStrategy.UP_THEN_DOWN,
) -> tuple[str, list[HGVSFix]]:
    """
    If hgvs_string names a transcript version that the data provider doesn't have,
    substitute the best available adjacent version (per ``strategy``) and rewrite
    the string.

    Returns:
        (resolved_string, fixes)
        fixes is empty if the requested version exists (non-destructive) or if the
        string carries no versioned transcript. A WARNING USED_ADJACENT_VERSION fix
        is returned when a substitution is made; an ERROR NO_TRANSCRIPT_VERSIONS fix
        (string unchanged) if the provider has no version of the accession at all.

    The data provider must implement ``get_tx_versions(accession)`` (JSONDataProvider
    and RESTDataProvider do). Providers that can't enumerate versions raise
    NotImplementedError, which is surfaced as an ERROR fix.

    Example::

        # data provider has NM_000059 versions [3, 4] but not the requested .2
        resolved, fixes = resolve_transcript_version("NM_000059.2:c.36del", dp)
        # resolved = "NM_000059.4:c.36del"
        # fixes[0].code = USED_ADJACENT_VERSION
    """
    parsed = _parse_versioned_transcript(hgvs_string)
    if parsed is None:
        return hgvs_string, []  # no versioned transcript — nothing to do
    versionless, requested_version = parsed

    try:
        available_versions = data_provider.get_tx_versions(versionless)
    except NotImplementedError as e:
        return hgvs_string, [HGVSFix(
            severity=HGVSFixSeverity.ERROR,
            code=HGVSFixCode.NO_TRANSCRIPT_VERSIONS,
            message=str(e),
        )]

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
    return resolved, [fix]


# ---------------------------------------------------------------------------
# Public: fix_hgvs — combined entry point
# ---------------------------------------------------------------------------

def fix_hgvs(
    hgvs_string: str,
    data_provider=None,
    genome_build: Optional[str] = None,
    tag_priority: list[str] = DEFAULT_TAG_PRIORITY,
    fallback_to_longest: bool = False,
    raise_on_errors: bool = False,
    ops: Optional[set] = None,
    prefer_consortium: Optional[Consortium] = DEFAULT_CONSORTIUM,
    version_fallback: Optional[VersionStrategy] = None,
) -> tuple[str, list[HGVSFix]]:
    """
    Clean and resolve an HGVS string in one call.

    Runs clean_hgvs() first (always), then resolve_gene_hgvs() if a
    data_provider and genome_build are supplied.  All fixes from both steps
    are returned in a single list.

    If data_provider/genome_build are omitted, only string cleaning is performed.
    This is useful when the caller knows the input already contains a transcript.

    If raise_on_errors=True, raises HGVSInputError on the first ERROR-level fix
    instead of returning it in the list.

    ops selects the subset of clean_hgvs cleaning operations to apply (see
    clean_hgvs); None (default) runs them all.

    prefer_consortium hard-filters the resolved transcript to RefSeq (default)
    or Ensembl; pass None to allow both (see resolve_gene_hgvs).

    version_fallback (default None = off) opts in to adjacent transcript-version
    fallback: if set to a VersionStrategy and a data_provider is supplied, an
    HGVS string whose exact transcript version isn't in the data is rewritten to
    the best available version, with a WARNING (see resolve_transcript_version).
    Only needs a data_provider, not a genome_build.

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
    result, fixes = clean_hgvs(hgvs_string, ops=ops)

    if data_provider is not None and genome_build is not None:
        result, resolution_fixes = resolve_gene_hgvs(
            result, data_provider, genome_build,
            tag_priority=tag_priority,
            fallback_to_longest=fallback_to_longest,
            prefer_consortium=prefer_consortium,
        )
        fixes.extend(resolution_fixes)

    if version_fallback is not None and data_provider is not None:
        result, version_fixes = resolve_transcript_version(
            result, data_provider, strategy=version_fallback,
        )
        fixes.extend(version_fixes)

    if raise_on_errors:
        errors = [f for f in fixes if f.severity == HGVSFixSeverity.ERROR]
        if errors:
            raise HGVSInputError(errors[0].message)

    return result, fixes
