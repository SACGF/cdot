""" Biocommons hgvs integration: data providers, and tools for fixing real-world HGVS input.

    The data providers and sequence fetchers are in :mod:`cdot.hgvs.dataproviders`. This package
    re-exports the HGVS fixing API from its submodules:

    - :mod:`cdot.hgvs.clean`: :func:`clean_hgvs` (pure string clean-up, no data needed) and the
      :class:`HGVSFix` records it reports
    - :mod:`cdot.hgvs.gene_hgvs`: :func:`fix_hgvs` (clean-up plus gene-to-transcript and
      transcript-version resolution against a data provider) and its building blocks
    - :mod:`cdot.hgvs.version_safety`: whether substituting another transcript version keeps a
      coordinate's meaning
    - :mod:`cdot.hgvs.translation`: :func:`fix_c_to_p`, protein prediction with the transcript's
      translation exceptions
"""
from cdot.hgvs.clean import (
    ALL_CLEAN_OPS,
    clean_hgvs,
    error_messages,
    get_best_transcript_version,
    HGVSCleanOp,
    HGVSCleanOpWarning,
    HGVSFix,
    HGVSFixCode,
    HGVSFixSeverity,
    HGVSInputError,
    messages,
    rank_transcript_versions,
    VersionStrategy,
    warning_messages,
)
from cdot.hgvs.gene_hgvs import (
    Consortium,
    consortium_of,
    DEFAULT_CONSORTIUM,
    DEFAULT_TAG_PRIORITY,
    DEFAULT_UNSAFE_VERSION_POLICY,
    fix_hgvs,
    rank_transcripts_for_gene,
    resolve_gene_hgvs,
    resolve_missing_accession_prefix,
    resolve_transcript_version,
    UnsafeVersionPolicy,
)
from cdot.hgvs.version_safety import (
    intrinsic_cds_structure,
    is_version_substitution_safe,
    is_version_substitution_safe_by_accession,
)
from cdot.hgvs.translation import fix_c_to_p

__all__ = [
    # cdot.hgvs.clean
    "ALL_CLEAN_OPS",
    "clean_hgvs",
    "error_messages",
    "get_best_transcript_version",
    "HGVSCleanOp",
    "HGVSCleanOpWarning",
    "HGVSFix",
    "HGVSFixCode",
    "HGVSFixSeverity",
    "HGVSInputError",
    "messages",
    "rank_transcript_versions",
    "VersionStrategy",
    "warning_messages",
    # cdot.hgvs.gene_hgvs
    "Consortium",
    "consortium_of",
    "DEFAULT_CONSORTIUM",
    "DEFAULT_TAG_PRIORITY",
    "DEFAULT_UNSAFE_VERSION_POLICY",
    "fix_hgvs",
    "rank_transcripts_for_gene",
    "resolve_gene_hgvs",
    "resolve_missing_accession_prefix",
    "resolve_transcript_version",
    "UnsafeVersionPolicy",
    # cdot.hgvs.version_safety
    "intrinsic_cds_structure",
    "is_version_substitution_safe",
    "is_version_substitution_safe_by_accession",
    # cdot.hgvs.translation
    "fix_c_to_p",
]
