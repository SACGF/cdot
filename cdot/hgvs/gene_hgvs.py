"""
Gene-symbol HGVS resolution via MANE/canonical transcript tags.

resolve_gene_hgvs() handles the common case of HGVS strings that name a gene
instead of a transcript (e.g. "BRCA2:c.36del").  It looks up the best available
transcript from the data provider's tag data (MANE_Select > MANE_Plus_Clinical >
RefSeq_Select > Ensembl_canonical > longest) and substitutes it.

fix_hgvs() is the single combined entry point that chains clean_hgvs() and
resolve_gene_hgvs() in one call — the recommended function for most callers.
"""

from typing import Optional

from cdot.hgvs.clean import (
    HGVSFix,
    HGVSFixCode,
    HGVSFixSeverity,
    HGVSInputError,
    _looks_like_transcript,
    clean_hgvs,
)

DEFAULT_TAG_PRIORITY: list[str] = [
    "MANE_Select",
    "MANE_Plus_Clinical",
    "RefSeq_Select",
    "Ensembl_canonical",
]


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


def _rank_transcripts_by_tags(
    tx_and_tags: list[tuple[str, list[str]]],
    tag_priority: list[str],
) -> list[tuple[str, list[str], Optional[str]]]:
    """
    Return [(tx_ac, tags, matched_tag_or_None), ...] sorted best-first.

    Transcripts whose tags include an earlier entry in tag_priority sort first.
    Transcripts with no matching priority tag sort last (preserving the
    input order, which is longest-first from get_tx_ac_tags_for_gene).
    """
    def sort_key(item: tuple[str, list[str]]) -> int:
        _tx_ac, tags = item
        for i, tag in enumerate(tag_priority):
            if tag in tags:
                return i
        return len(tag_priority)

    ranked = sorted(tx_and_tags, key=sort_key)
    result = []
    for tx_ac, tags in ranked:
        matched = next((t for t in tag_priority if t in tags), None)
        result.append((tx_ac, tags, matched))
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

    if not tx_and_tags:
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.ERROR,
            code=HGVSFixCode.NO_TRANSCRIPT_FOR_GENE,
            message=f"No transcripts found for gene '{gene_symbol}' in {genome_build}",
        ))
        return hgvs_string, fixes

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
# Public: fix_hgvs — combined entry point
# ---------------------------------------------------------------------------

def fix_hgvs(
    hgvs_string: str,
    data_provider=None,
    genome_build: Optional[str] = None,
    tag_priority: list[str] = DEFAULT_TAG_PRIORITY,
    fallback_to_longest: bool = False,
    raise_on_errors: bool = False,
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

    Example — gene-only input resolved via MANE::

        result, fixes = fix_hgvs("BRCA2:c.36DEL", data_provider, "GRCh38")
        # result  = "NM_000059.4:c.36del"
        # fixes   = [HGVSFix(WARNING, LOWERCASED_MUTATION_TYPE, ...),
        #            HGVSFix(WARNING, RESOLVED_GENE_TO_TRANSCRIPT, ...)]

    Example — cleaning only (no data provider)::

        result, fixes = fix_hgvs("NM_000059.4 c.316+5G>A")
        # result = "NM_000059.4:c.316+5G>A"
    """
    result, fixes = clean_hgvs(hgvs_string)

    if data_provider is not None and genome_build is not None:
        result, resolution_fixes = resolve_gene_hgvs(
            result, data_provider, genome_build,
            tag_priority=tag_priority,
            fallback_to_longest=fallback_to_longest,
        )
        fixes.extend(resolution_fixes)

    if raise_on_errors:
        errors = [f for f in fixes if f.severity == HGVSFixSeverity.ERROR]
        if errors:
            raise HGVSInputError(errors[0].message)

    return result, fixes
