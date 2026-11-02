"""
HGVS string cleaning, normalisation, and transcript-version fallback.

clean_hgvs() is a pure string operation — no genome build, no HGVS parser,
no data provider required.  It fixes the most common real-world formatting
mistakes found in clinical reports, search boxes, and literature.

get_best_transcript_version() is an opt-in helper for callers that want to
fall back to an adjacent transcript version when the exact one is not found
in their data.  It is deliberately *not* called automatically so that callers
retain full control over HGVS compatibility.
"""

import re
import warnings
from dataclasses import dataclass, replace
from enum import Enum
from functools import lru_cache, wraps
from typing import Callable, Optional


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------

class HGVSInputError(ValueError):
    """Raised for malformed HGVS input that could not be automatically fixed."""
    pass


# ---------------------------------------------------------------------------
# Fix types
# ---------------------------------------------------------------------------

class HGVSFixSeverity(Enum):
    WARNING = "warning"  # Fixed something non-standard but unambiguously correctable
    ERROR   = "error"    # Could not fix; input remains invalid


class HGVSFixCode(Enum):
    # String-level fixes (all WARNING)
    STRIPPED_WHITESPACE          = "stripped_whitespace"
    STRIPPED_PROTEIN_SUFFIX      = "stripped_protein_suffix"
    FIXED_MULTIPLE_COLON         = "fixed_multiple_colon"
    FIXED_MULTIPLE_KIND          = "fixed_multiple_kind"
    ADDED_N_PREFIX               = "added_n_prefix"
    LOWERCASED_MUTATION_TYPE     = "lowercased_mutation_type"
    STRIPPED_UNBALANCED_BRACKETS = "stripped_unbalanced_brackets"
    STRIPPED_SURROUNDING_PUNCTUATION = "stripped_surrounding_punctuation"
    REPLACED_UNICODE_DASH        = "replaced_unicode_dash"
    FIXED_MULTIPLE_DOT           = "fixed_multiple_dot"
    FIXED_PREFIX_COLON           = "fixed_prefix_colon"
    STRIPPED_LEADING_JUNK        = "stripped_leading_junk"
    STRIPPED_CONTIG_NAME         = "stripped_contig_name"
    STRIPPED_BUILD_NAME          = "stripped_build_name"
    FIXED_GENE_WRAPPER           = "fixed_gene_wrapper"
    WRAPPED_GENE_SYMBOL          = "wrapped_gene_symbol"
    FIXED_SEPARATOR_TYPO         = "fixed_separator_typo"
    DROPPED_DEL_DUP_COUNT        = "dropped_del_dup_count"
    DROPPED_GENOMIC_REF_IN_PARENS = "dropped_genomic_ref_in_parens"
    ADDED_TRANSCRIPT_UNDERSCORE  = "added_transcript_underscore"
    FIXED_MULTIPLE_UNDERSCORE    = "fixed_multiple_underscore"
    DROPPED_EMPTY_VERSION        = "dropped_empty_version"
    UPPERCASED_BASES             = "uppercased_bases"
    RECONSTRUCTED_STRUCTURE      = "reconstructed_structure"
    ADDED_MISSING_KIND           = "added_missing_kind"
    # Gene/transcript fixes (all WARNING)
    SWAPPED_GENE_TRANSCRIPT      = "swapped_gene_transcript"
    # g./m. HGVS with a genomic and a transcript accession: the genomic one is kept
    DROPPED_TRANSCRIPT_FOR_GENOMIC_KIND = "dropped_transcript_for_genomic_kind"
    UPPERCASED_TRANSCRIPT        = "uppercased_transcript"
    UPPERCASED_HGVS_PREFIX       = "uppercased_hgvs_prefix"
    UPPERCASED_GENE_SYMBOL       = "uppercased_gene_symbol"
    # Transcript version fallback (WARNING on substitution, ERROR if none available)
    USED_ADJACENT_VERSION        = "used_adjacent_version"
    # Version fallback with a structural coordinate-safety verdict (#28). The plain
    # USED_ADJACENT_VERSION above is still emitted when no safety check could be run
    # (eg a provider that can't read the requested version's structure).
    USED_ADJACENT_VERSION_COORD_SAFE        = "used_adjacent_version_coord_safe"
    USED_ADJACENT_VERSION_COORD_UNVERIFIED  = "used_adjacent_version_coord_unverified"
    REFUSED_UNSAFE_VERSION                  = "refused_unsafe_version"
    NO_TRANSCRIPT_VERSIONS       = "no_transcript_versions"
    # Versionless accession (eg NM_000059:c.36del) given a version by the fallback,
    # because every held version agrees at the cited position
    RESOLVED_VERSIONLESS_TRANSCRIPT = "resolved_versionless_transcript"
    # Gene → transcript resolution (WARNING on success, ERROR on failure)
    RESOLVED_GENE_TO_TRANSCRIPT  = "resolved_gene_to_transcript"
    RESOLVED_GENE_TO_PROTEIN     = "resolved_gene_to_protein"  # gene-only p. HGVS
    NO_TRANSCRIPT_FOR_GENE       = "no_transcript_for_gene"
    SKIPPED_GENE_RESOLUTION      = "skipped_gene_resolution"  # WARNING: fix_hgvs got no genome_build
    # Accession prefix restoration for a bare-number accession, verified against
    # a data provider (see gene_hgvs.resolve_missing_accession_prefix; WARNING)
    ADDED_ACCESSION_PREFIX       = "added_accession_prefix"
    # c. to p. translation checks from cdot's translation data (see translation.fix_c_to_p)
    USED_TRANSLATION_TABLE       = "used_translation_table"         # WARNING
    NEW_UGA_READ_AS_STOP         = "new_uga_read_as_stop"           # WARNING
    UGA_MAY_BE_READ_AS_SEC       = "uga_may_be_read_as_sec"         # WARNING
    TRANSLATION_EXCEPTION_NOT_APPLIED = "translation_exception_not_applied"  # WARNING
    TRANSLATION_DATA_INCOMPLETE  = "translation_data_incomplete"    # WARNING
    REFERENCE_CODON_MISMATCH     = "reference_codon_mismatch"       # WARNING, ERROR if it makes a stop
    UNSUPPORTED_TRANSLATION_TABLE = "unsupported_translation_table"  # ERROR
    RIBOSOMAL_SLIPPAGE_UNSUPPORTED = "ribosomal_slippage_unsupported"  # ERROR
    # Problems a cleaning op found but could not fix without guessing (ERROR)
    DEL_DUP_COUNT_MISMATCH       = "del_dup_count_mismatch"
    UNBALANCED_BRACKETS          = "unbalanced_brackets"
    # Validation errors (ERROR)
    NO_COLON                     = "no_colon"
    MISSING_REFERENCE_SEQUENCE   = "missing_reference_sequence"
    KIND_REFERENCE_MISMATCH      = "kind_reference_mismatch"  # eg NM_ with g., NC_ with c.
    INS_WITH_INTEGER_LENGTH      = "ins_with_integer_length"
    INS_MISSING_SEQUENCE         = "ins_missing_sequence"


class HGVSCleanOp(Enum):
    """A single cleaning operation in the clean_hgvs() pipeline.

    Unlike HGVSFixCode (which describes *outcomes* and is not 1:1 with steps),
    each member here names exactly one step, so it is a stable key for selecting
    a subset of cleaning operations.  Members are listed in pipeline order, but
    selection only *filters* the pipeline — it never reorders it (several steps
    depend on earlier ones running first). Selecting a step without a step it
    relies on warns with HGVSCleanOpWarning.
    """
    STRIP_SURROUNDING_PUNCTUATION = "strip_surrounding_punctuation"
    REPLACE_UNICODE_DASH     = "replace_unicode_dash"
    STRIP_PROTEIN_SUFFIX     = "strip_protein_suffix"
    STRIP_CONTIG_NAME        = "strip_contig_name"
    STRIP_BUILD_NAME         = "strip_build_name"
    STRIP_LEADING_JUNK       = "strip_leading_junk"
    FIX_SYMBOL_BEFORE_ACCESSION = "fix_symbol_before_accession"
    FIX_SPACE_BEFORE_KIND    = "fix_space_before_kind"
    STRIP_WHITESPACE         = "strip_whitespace"
    ADD_N_PREFIX             = "add_n_prefix"
    ADD_TRANSCRIPT_UNDERSCORE = "add_transcript_underscore"
    FIX_MULTIPLE_COLON       = "fix_multiple_colon"
    FIX_GENE_WRAPPER         = "fix_gene_wrapper"
    FIX_MULTIPLE_DOT         = "fix_multiple_dot"
    FIX_MULTIPLE_UNDERSCORE  = "fix_multiple_underscore"
    FIX_PREFIX_COLON         = "fix_prefix_colon"
    FIX_MULTIPLE_KIND        = "fix_multiple_kind"
    FIX_SEPARATOR_TYPO       = "fix_separator_typo"
    DROP_EMPTY_VERSION       = "drop_empty_version"
    DROP_DEL_DUP_COUNT       = "drop_del_dup_count"
    STRIP_UNBALANCED_BRACKETS = "strip_unbalanced_brackets"
    RECONSTRUCT_STRUCTURE    = "reconstruct_structure"
    LOWERCASE_MUTATION_TYPE  = "lowercase_mutation_type"
    UPPERCASE_BASES          = "uppercase_bases"
    ADD_MISSING_KIND         = "add_missing_kind"
    FIX_GENE_TRANSCRIPT      = "fix_gene_transcript"
    DROP_GENOMIC_REF_IN_PARENS = "drop_genomic_ref_in_parens"


@dataclass(frozen=True)
class HGVSFix:
    severity: HGVSFixSeverity
    code: HGVSFixCode
    message: str
    original: Optional[str] = None  # before-value, if applicable
    fixed: Optional[str] = None     # after-value, if applicable

    def __str__(self) -> str:
        return self.message


def messages(fixes: list[HGVSFix], severity: Optional[HGVSFixSeverity] = None) -> list[str]:
    """Return the messages of ``fixes``, optionally filtered to one severity.

    Saves consumers repeating ``[f.message for f in fixes if f.severity == ...]``
    at every call site.
    """
    return [f.message for f in fixes if severity is None or f.severity == severity]


def warning_messages(fixes: list[HGVSFix]) -> list[str]:
    """Return the messages of all WARNING-level fixes."""
    return messages(fixes, HGVSFixSeverity.WARNING)


def error_messages(fixes: list[HGVSFix]) -> list[str]:
    """Return the messages of all ERROR-level fixes."""
    return messages(fixes, HGVSFixSeverity.ERROR)


# ---------------------------------------------------------------------------
# Transcript prefix helpers
# Ported from vg_code/transcripts_utils.py (itself from pyhgvs RefSeq table)
# ---------------------------------------------------------------------------

_REFSEQ_PREFIX_TYPE: dict[str, str] = {
    'AC_': 'genomic', 'NC_': 'genomic', 'NG_': 'genomic',
    'NT_': 'genomic', 'NW_': 'genomic', 'NS_': 'genomic', 'NZ_': 'genomic',
    'NM_': 'mRNA',    'NR_': 'RNA',
    'XM_': 'mRNA',    'XR_': 'RNA',
    'AP_': 'Protein', 'NP_': 'Protein', 'YP_': 'Protein',
    'XP_': 'Protein', 'ZP_': 'Protein',
}

# Non-RefSeq accession prefixes cdot handles. Ensembl transcripts have no
# underscore after the prefix; LRG references do, and may carry a "t1"/"p1"
# transcript/protein suffix.
_ENSEMBL_TRANSCRIPT_PREFIX = 'ENST'
_LRG_PREFIX = 'LRG_'

# Every HGVS kind letter ("c.", "g.", "m.", "n.", "p.", "r."), grouped by the kind of
# reference sequence its positions are on. This is the one place they are defined:
# every pattern below that matches a kind uses _KIND_LETTERS.
_GENOMIC_KINDS = frozenset("gm")
_TRANSCRIPT_KINDS = frozenset("cnr")
_PROTEIN_KINDS = frozenset("p")
_KIND_LETTERS = "".join(sorted(_GENOMIC_KINDS | _TRANSCRIPT_KINDS | _PROTEIN_KINDS))
_KIND = f"[{_KIND_LETTERS}]"

# The kind implied by each reference-sequence type, used to insert a missing one
_KIND_FOR_TYPE = {'mRNA': 'c', 'RNA': 'n', 'genomic': 'g'}


def _refseq_prefixes(*types: str) -> list[str]:
    """ RefSeq prefixes (eg 'NM_') of the given types, or all of them if none given """
    return sorted(p for p, t in _REFSEQ_PREFIX_TYPE.items() if not types or t in types)


def _alternation(prefixes) -> str:
    """ Non-capturing regex alternation of literal prefixes """
    return "(?:" + "|".join(re.escape(p) for p in prefixes) + ")"


# Any accession prefix: every RefSeq type, Ensembl transcript and LRG
_ANY_ACCESSION_PREFIX = _alternation(
    _refseq_prefixes() + [_ENSEMBL_TRANSCRIPT_PREFIX, _LRG_PREFIX])
# As above minus LRG, whose accessions legitimately contain a second underscore
_NON_LRG_ACCESSION_PREFIX = _alternation(_refseq_prefixes() + [_ENSEMBL_TRANSCRIPT_PREFIX])
_TRANSCRIPT_ACCESSION_PREFIX = _alternation(
    _refseq_prefixes('mRNA', 'RNA') + [_ENSEMBL_TRANSCRIPT_PREFIX, _LRG_PREFIX])
_GENOMIC_ACCESSION_PREFIX = _alternation(_refseq_prefixes('genomic'))


def _looks_like_transcript(s: str) -> bool:
    """True if s starts with a known mRNA/RNA RefSeq prefix, or is Ensembl/LRG.

    Case-sensitive — callers must uppercase s themselves if they want a
    case-insensitive check (matching VG behaviour).
    """
    if s.startswith(_ENSEMBL_TRANSCRIPT_PREFIX) or s.startswith(_LRG_PREFIX):
        return True
    # Keys in _REFSEQ_PREFIX_TYPE are uppercase 3-char prefixes e.g. 'NM_'
    return _REFSEQ_PREFIX_TYPE.get(s[:3]) in ('mRNA', 'RNA')


def _looks_like_hgvs_prefix(s: str) -> bool:
    """True if s starts with any known RefSeq prefix, or is Ensembl/LRG.

    Case-sensitive — callers must uppercase s themselves if needed.
    """
    if s.startswith(_ENSEMBL_TRANSCRIPT_PREFIX) or s.startswith(_LRG_PREFIX):
        return True
    return s[:3] in _REFSEQ_PREFIX_TYPE


# A known accession prefix (any case) at the start of a string, followed by the
# accession digits. Only the prefix is uppercased, so an LRG "t1" suffix stays.
_ACCESSION_PREFIX_AT_START = re.compile(rf"^{_ANY_ACCESSION_PREFIX}(?=\d)", re.IGNORECASE)


def _uppercase_accession_prefix(s: str) -> str:
    """ s with a leading known accession prefix uppercased (eg "ng_012337.3" ->
        "NG_012337.3"). Returned unchanged if it doesn't start with one. """
    if m := _ACCESSION_PREFIX_AT_START.match(s):
        return m.group(0).upper() + s[m.end():]
    return s


def _reference_type(s: str) -> Optional[str]:
    """ Reference-sequence type ('mRNA', 'RNA', 'genomic' or 'Protein') of the
        known accession prefix (any case) s starts with. An Ensembl transcript is
        'mRNA'. None for an LRG or anything without a known prefix. """
    if not (m := _ACCESSION_PREFIX_AT_START.match(s)):
        return None
    prefix = m.group(0).upper()
    if prefix == _ENSEMBL_TRANSCRIPT_PREFIX:
        return 'mRNA'
    return _REFSEQ_PREFIX_TYPE.get(prefix)


# ---------------------------------------------------------------------------
# Compiled regex patterns
# ---------------------------------------------------------------------------

# A trailing protein annotation after a nucleotide allele, eg
# "NM_000059.4:c.68del p.(Glu23fs)". Any whitespace separates them (a tab is as
# common as a space in pasted text), or a comma/semicolon, which is removed with
# the annotation ("NM_000059.4:c.68del, p.(Glu23fs)").
_P_HGVS_REMOVAL = re.compile(r"^(?P<main>.*?)(?:\s*[,;]\s*|\s+)(p\.|[(]p\.).*$")
# The part before the protein annotation must hold a nucleotide allele (a
# non-protein kind, or something other than "p." after the colon). Otherwise the
# "p." is the allele itself, eg "NP_000050.2 p.Arg1Gly", and must be kept.
_HAS_NUCLEOTIDE_ALLELE = re.compile(r":\s*(?!p\.)\S|[cgmnr]\.", re.IGNORECASE)

# Two-letter RefSeq accession prefixes (NM, NR, XM, XR, NC, NG, ...) taken from
# the type table, used to repair an accession missing its underscore. Requires
# 6+ digits (every RefSeq accession has at least six) so a real gene symbol that
# starts with these letters plus a stray digit (eg NR3C1, NS1, ZP3) is never
# mistaken for an accession and mangled. A stray bracket, dot or colon typed in
# place of the underscore ("NM)000059.4", "NM.000059.4", "NM:000059.4") is
# replaced by it.
_TWO_LETTER_PREFIXES = sorted({p[:2] for p in _REFSEQ_PREFIX_TYPE})
_TRANSCRIPT_NO_UNDERSCORE = re.compile(
    rf"^({'|'.join(_TWO_LETTER_PREFIXES)})[().:]?(\d{{6,}})", re.IGNORECASE)

# Second letters for which prepending "N" yields a real genomic/transcript prefix
# (NM, NC, NR, NG, NT, NW, NS, NZ, NP) — used to restore a dropped leading "N".
_N_PREFIX_SECOND_LETTERS = frozenset(
    p[1] for p in _REFSEQ_PREFIX_TYPE if p.startswith("N"))

# A coordinate may start with a digit, "-" (5' UTR), "*" (3' UTR), "(" (an
# uncertain position) or "[" (an allele).
_COORD_START = r"[\d*\-(\[]"

# Lenient pattern, matching HGVS strings whether or not they have a transcript,
# gene symbol, version, colon, or dot. _parse_loose() reads its groups into the
# fields, which the canonical form is serialised from.
#
# Kept tight on purpose (#150), as whatever it matches is rewritten:
#  * an accession needs a real prefix (with its underscore) and digits, so a
#    gene symbol starting with prefix letters (NR3C1, ACTB) is not read as one;
#  * the version is bounded, so its digits can't be split off into a coordinate;
#  * the gene symbol is letters, digits and "-" only (never '>', '+', ')');
#  * the kind letter must be followed by its '.' (or a stray ':'), or directly
#    by a coordinate. An uppercase kind letter is only accepted with the '.'/':',
#    so a base or the last letter of a gene symbol (eg "FANCG100A>G") is not
#    taken as the kind.
_LOOSE_HGVS_PATTERN = re.compile(
    # transcript_value allows an optional LRG transcript/protein suffix (eg the
    # "t1" in "LRG_199t1") so the suffix isn't mis-parsed as the gene symbol.
    rf"^((?P<transcript_prefix>(?i:{_ANY_ACCESSION_PREFIX}))(?P<transcript_value>\d+(?:[tp]\d+)?)"
    r"(\.(?P<transcript_version>\d+)(?!\d))?)?"
    r"([(]?(?P<gene_symbol>[A-Za-z][A-Za-z0-9-]*[A-Za-z0-9])[)]?)?"
    # A stray colon between the kind letter and the coordinate (eg "...2c:1188")
    # is accepted as the kind separator, so the same glued-kind input reconstructs
    # whether or not that colon is present (it is just misplaced, belonging before
    # the kind letter rather than after it).
    rf":?(?:(?P<letter>{_KIND})(?:[.:]|(?={_COORD_START}))|(?P<upper_letter>{_KIND.upper()})[.:])"
    r"(?P<nomen>[^:.]*)$",
)

# An allele with no kind letter: it starts with a coordinate, or a later colon on
# the same line is followed by one.
_ALLELE_WITHOUT_KIND = re.compile(rf"(?:.*:)?{_COORD_START}")

_REF_ALT_NUC     = re.compile(r"(?P<ref>[gatc]+)>(?P<alt>[gatc=]+)$", re.IGNORECASE)
_DEL_INS_DUP_NUC = re.compile(
    r"(del(?P<del>[gatc]+))?(?P<op>ins|dup|del)(?P<ins>[gatc]*)$",
    re.IGNORECASE,
)
# RNA ("r.") alleles are written in lowercase (eg r.76a>c), so their bases must
# not be uppercased
_RNA_ALLELE = re.compile(r":r\.")

_GENE_IN_PARENS = re.compile(r"(.+)\((.*)\)")

_INS_INTEGER = re.compile(r".*ins\d+$", re.IGNORECASE)
_INS_EMPTY   = re.compile(r".*ins$",    re.IGNORECASE)

# A string with no reference sequence: it starts with ':' or '.', or with a bare
# kind letter followed by its '.'/':' or directly by a coordinate. A lowercase
# accession that merely starts with a kind letter (eg "ng_012337.3") is not one.
_MISSING_REFERENCE = re.compile(rf"^(?:[:.]|{_KIND}(?:[.:]|{_COORD_START}))")

# Misplaced colon between a RefSeq prefix and its underscore, e.g. "NM:_000059"
_PREFIX_COLON_UNDERSCORE = re.compile(
    rf"^({'|'.join(_TWO_LETTER_PREFIXES)}):_", re.IGNORECASE)

# Leading junk before an accession, e.g. "#", ":", "-", "*", a bullet, a literal
# "\t", whitespace (": NM_..."), or a genome-build prefix like "GRCh38.p2 " /
# "GRCh37.p12#". Only stripped when an accession actually follows.
_LEADING_JUNK = re.compile(r"^(?:GRCh\d+(?:\.p\d+)?[#\s]*|\\t|[#:\-*\u2022\u00b7\s])+", re.IGNORECASE)

# Gene symbol wedged between extra colons, e.g. "tx:(GENE):c." / "tx:(GENE)c." /
# "tx:GENE:c." — should be "tx(GENE):c."
_GENE_WRAPPER = re.compile(
    rf"^({_ANY_ACCESSION_PREFIX}[\w.]+):\(?([A-Za-z][A-Za-z0-9-]+)\)?:?({_KIND}[.,])", re.IGNORECASE)
# Gene symbol repeated after its parens, e.g. "tx(GENE):GENE:c." -> "tx(GENE):c."
_GENE_REPEATED_AFTER_PARENS = re.compile(rf"\(([A-Za-z][A-Za-z0-9-]*)\):\1:(?={_KIND}\.)", re.IGNORECASE)
# Stray colon inside the gene parens, e.g. "(RAD51C:)" -> "(RAD51C)"
_GENE_PARENS_COLON = re.compile(r"\(([A-Za-z][A-Za-z0-9-]*):\)")
# Parens wrapping the whole accession, e.g. "(NM_000059.4):c." -> "NM_000059.4:c."
_PARENS_WRAPPED_ACCESSION = re.compile(rf"^\(({_ANY_ACCESSION_PREFIX}[\w.]+)\)(?=:)", re.IGNORECASE)
# Missing colon between a (gene) wrapper and the kind, e.g.
# "BRCA1(NM_000059.4)c.123A>G" -> "BRCA1(NM_000059.4):c.123A>G" (the
# gene/transcript swap is then done by _fix_gene_transcript).
_MISSING_COLON_BEFORE_KIND = re.compile(rf"(\([^)]+\))({_KIND}\.)", re.IGNORECASE)

# Comma used in place of the c./g. dot, e.g. ":c,1811" -> ":c.1811"
_KIND_COMMA = re.compile(rf":({_KIND}),", re.IGNORECASE)
# Colon used in place of the kind dot, e.g. ":c:1811" -> ":c.1811". The
# lookahead keeps this from touching a normal ":c." (the second char is a dot,
# not a coordinate).
_KIND_COLON = re.compile(rf":({_KIND}):(?=[\d*?(\-])", re.IGNORECASE)
# Period used in place of the substitution '>', e.g. "1030C.T" -> "1030C>T"
_SUB_PERIOD = re.compile(r"([ACGT])\.([ACGT])$", re.IGNORECASE)

# Semicolon or underscore used in place of the accession:allele colon, e.g.
# "NM_000059.4;c.68del" or "NM_000059.4_c.68del". Anchored to the start so a
# ';' inside a multi-variant allele is never touched. LRG is excluded (its
# accessions legitimately contain underscores).
_ACCESSION_KIND_SEPARATOR = re.compile(
    rf"^({_NON_LRG_ACCESSION_PREFIX}\d+(?:\.\d+)?)[;_](?={_KIND}\.)", re.IGNORECASE)

# Whitespace used in place of the reference:kind colon after a bare gene symbol,
# e.g. "BRCA2 c.68del" -> "BRCA2:c.68del". Restricted to a plain symbol (no
# '_'/'.'/':') so accessions are left to the whitespace/reconstruction steps,
# and runs before whitespace stripping (which would otherwise glue the symbol to
# the kind letter, mangling symbols that resemble an accession, eg NR2E3).
_SPACE_BEFORE_KIND = re.compile(rf"^([A-Za-z][A-Za-z0-9]{{1,9}})\s+(?={_KIND}\.)", re.IGNORECASE)

# Chromosome names of the human assemblies cdot supports, plus "M" (UCSC's chrM)
_CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y", "M", "MT"]
_CHROMOSOME_NAMES = frozenset(_CHROMOSOMES)

# A chromosome name, whitespace, then a chromosomal accession, e.g.
# "chr13 NC_000013.11:g.32316528del". Only the exact "chr" + chromosome forms,
# as many gene symbols start with CHR (CHRM1, CHRNA7, CHRD...); no HGNC symbol,
# alias or previous symbol matches this pattern.
_CHROM_NAME = rf"chr(?:{'|'.join(_CHROMOSOMES)})"
_CONTIG_NAME_BEFORE_ACCESSION = re.compile(
    rf"(?<![\w.])({_CHROM_NAME})\s+(?=(NC_\d+)(\.\d+)?(?![\w.]))", re.IGNORECASE)

# A gene symbol, whitespace (or a ':' or '.'), then an accession, e.g.
# "APC NM_001407446.1:c.-40G>A" or "APC:NM_001407446.1:c.-40G>A" ->
# "NM_001407446.1(APC):c.-40G>A". Runs before
# whitespace stripping, which would otherwise glue the symbol onto the accession
# ("APCNM_001407446.1"), leaving something that parses as neither.
_SYMBOL_BEFORE_ACCESSION = re.compile(
    rf"^(?!{_CHROM_NAME}[\s:.])([A-Za-z][A-Za-z0-9-]{{1,9}})(?:\s*[:.]\s*|\s+)({_ANY_ACCESSION_PREFIX}[\w.]+)", re.IGNORECASE)

# A genome build / assembly name at the start, e.g. "hg19 NM_000059.4:c.68del" or
# "GRCh38 BRCA2 NM_000059.4:c.68del". Without this the symbol-before-accession
# rule above would take the build name for a gene symbol. Only the exact build
# names (UCSC hg19/hg38, GRCh37/GRCh38 with an optional patch, T2T-CHM13), none
# of which is an HGNC gene symbol, followed by whitespace or a '#'/':' separator
# and then the rest of the string (never the whole input).
_BUILD_NAME = r"(?:GRCh3[78](?:\.p\d+)?|hg(?:19|37|38)|T2T(?:-CHM13(?:v[\d.]+)?)?|CHM13(?:v[\d.]+)?)"
_BUILD_NAME_BEFORE_REFERENCE = re.compile(rf"^{_BUILD_NAME}[\s#:]+(?=\S)", re.IGNORECASE)

# An accession with an empty version field: a dot and no digits after it, e.g.
# the "NM_001754." of "NM_001754.:c.749G>A".
_EMPTY_VERSION_ACCESSION = re.compile(
    rf"({_NON_LRG_ACCESSION_PREFIX}\d+)\.", re.IGNORECASE)

# del/dup with an explicit range followed by a redundant base count, e.g.
# "c.1315_1337dup23" -> "c.1315_1337dup" (the range already gives the length;
# normalise to plain del/dup). The range guard avoids changing the meaning of a
# single-position "c.123del5" (delete 5 bases). delins/ins are excluded — they
# genuinely need the inserted sequence, so they stay flagged as errors.
_DEL_DUP_RANGE_COUNT = re.compile(r"(\d+_\d[\d+\-_*]*(?:del|dup))\d+$", re.IGNORECASE)
# The subset of the above whose range is two plain positions (no 5'/3' UTR "-"
# or "*", no intronic offset), so its length is known without the transcript and
# the count can be checked against it before it is dropped.
_DEL_DUP_PLAIN_RANGE_COUNT = re.compile(
    r"(?<![\d+\-*])(?P<start>\d+)_(?P<end>\d+)(?P<op>del|dup)(?P<count>\d+)$", re.IGNORECASE)

# The kind letter and its dot, used to find where the allele starts when the
# reference:allele colon is missing, and to split the kind off an allele
_KIND_DOT = re.compile(rf"{_KIND}\.")

# A transcript accession in the reference, with a genomic accession in its
# gene-symbol parentheses (see _op_drop_genomic_ref_in_parens)
_TRANSCRIPT_ACCESSION = re.compile(rf"{_TRANSCRIPT_ACCESSION_PREFIX}[\w.]+", re.IGNORECASE)
_GENOMIC_ACCESSION = re.compile(rf"{_GENOMIC_ACCESSION_PREFIX}\d[\w.]*", re.IGNORECASE)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _fix_del_ins_bases(m: re.Match) -> str:
    parts = []
    if del_nucs := m.group('del'):
        parts.append(f"del{del_nucs.upper()}")
    parts.append(m.group('op').lower())
    parts.append(m.group('ins').upper())
    return "".join(parts)


# ---------------------------------------------------------------------------
# Parsed fields
#
# The structural ops work on an HGVS string split into its fields rather than
# on the string itself. _parse() splits a string with a colon into _HGVSParts,
# and str() serialises the fields back. Parsing is lossless, so a structural op
# that changes nothing returns the string it was given.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class _HGVSParts:
    """ The fields of "accession(gene_symbol):kind.nomen".

        accession is the reference outside the parentheses and gene_symbol what is
        inside them (None without parentheses). They hold whatever the input put
        there, eg a gene symbol before its accession until the gene/transcript
        fix swaps them. kind is the kind letter when the allele starts with one
        and its dot, otherwise None and nomen is the whole allele. """
    accession: str
    gene_symbol: Optional[str]
    kind: Optional[str]
    nomen: str

    @property
    def reference(self) -> str:
        if self.gene_symbol is None:
            return self.accession
        return f"{self.accession}({self.gene_symbol})"

    @property
    def allele(self) -> str:
        if self.kind:
            return f"{self.kind}.{self.nomen}"
        return self.nomen

    def __str__(self) -> str:
        return f"{self.reference}:{self.allele}"


def _parse(s: str) -> Optional[_HGVSParts]:
    """ s split into its fields at the first colon, or None if it has no colon.
        str() of the result is s. """
    reference, sep, allele = s.partition(":")
    if not sep:
        return None
    if m := _GENE_IN_PARENS.fullmatch(reference):
        accession, gene_symbol = m.groups()
    else:
        accession, gene_symbol = reference, None
    if _KIND_DOT.match(allele):
        kind, nomen = allele[0], allele[2:]
    else:
        kind, nomen = None, allele
    return _HGVSParts(accession, gene_symbol, kind, nomen)


def _parse_loose(s: str) -> Optional[tuple[_HGVSParts, bool]]:
    """ The canonical fields of a string matching _LOOSE_HGVS_PATTERN, which may
        lack the colon or the kind's dot, or have the kind letter or accession
        prefix in the wrong case, and whether the accession prefix had to be
        uppercased. None if it doesn't match, or has neither an accession nor a
        gene symbol. """
    if not (m := _LOOSE_HGVS_PATTERN.match(s)):
        return None
    kind = (m.group("letter") or m.group("upper_letter")).lower()
    gene_symbol = m.group("gene_symbol")
    if transcript_prefix := m.group("transcript_prefix"):
        accession = f"{transcript_prefix.upper()}{m.group('transcript_value')}"
        if transcript_version := m.group("transcript_version"):
            accession += f".{transcript_version}"
        parts = _HGVSParts(accession, gene_symbol, kind, m.group("nomen"))
        return parts, transcript_prefix != transcript_prefix.upper()
    if gene_symbol:
        return _HGVSParts(gene_symbol, None, kind, m.group("nomen")), False
    return None


_StringOp = Callable[[str], tuple[str, list[HGVSFix]]]
_StructuralOp = Callable[[_HGVSParts], tuple[_HGVSParts, list[HGVSFix]]]


def _on_fields(fn: _StructuralOp) -> _StringOp:
    """ A structural op (on _HGVSParts) as a string op: parse, run, serialise.
        A string with no colon has no fields and is returned unchanged. """
    @wraps(fn)
    def string_op(s: str) -> tuple[str, list[HGVSFix]]:
        if (parts := _parse(s)) is None:
            return s, []
        parts, fixes = fn(parts)
        return str(parts), fixes
    return string_op


# The kind letters each reference-sequence type can take, and how to name the type in a message.
# A coding transcript can't take n.: n. counts from the transcript start and c. from the start
# codon, so the same number is a different base.
_KINDS_FOR_TYPE = {'mRNA': (_TRANSCRIPT_KINDS - {"n"}, "coding transcript"),
                   'RNA': (_TRANSCRIPT_KINDS - {"c"}, "non-coding transcript"),
                   'genomic': (_GENOMIC_KINDS, "genomic"), 'Protein': (_PROTEIN_KINDS, "protein")}


def _kinds_for_reference(accession: str) -> Optional[tuple[frozenset[str], str]]:
    """ (kind letters that can follow accession, description of it), eg ({"n", "r"}, "non-coding
        transcript") for NR_, or None if it isn't a known RefSeq or Ensembl transcript accession.
        An Ensembl transcript may be coding or not, so takes either. """
    if accession.upper().startswith(_ENSEMBL_TRANSCRIPT_PREFIX):
        return _TRANSCRIPT_KINDS, "transcript"
    return _KINDS_FOR_TYPE.get(_reference_type(accession))


_N_ON_CODING = ": n. counts from the transcript start, not the start codon"


def _is_chromosome_name(reference: str) -> bool:
    """ True if reference is a chromosome name, with or without "chr" ("chr13", "13", "X", "chrM") """
    return reference.upper().removeprefix("CHR") in _CHROMOSOME_NAMES


# Gene names as HGNC and Ensembl write them ("BRCA2", "HLA-A", "7SK", "Y_RNA", "AC004556.1")
_GENE_SYMBOL = re.compile(r"[A-Za-z0-9][\w.-]*")


def _looks_like_gene_symbol(reference: str) -> bool:
    """ True if reference (the part of an HGVS string before ':') could be a gene symbol:
        no parentheses, not an accession in any case ("nm_000059.4"), and not a chromosome name. """
    return (bool(_GENE_SYMBOL.fullmatch(reference))
            and not _looks_like_hgvs_prefix(reference.upper())
            and not _is_chromosome_name(reference))


def _or_list(items: list[str]) -> str:
    """ "a", "a or b", "a, b or c" """
    return " or ".join([", ".join(items[:-1]), items[-1]] if len(items) > 2 else items)


def _validate(s: str) -> list[HGVSFix]:
    """Post-cleaning validation — returns ERROR-level fixes for unfixable problems."""
    fixes = []
    if _INS_INTEGER.match(s):
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.ERROR,
            code=HGVSFixCode.INS_WITH_INTEGER_LENGTH,
            message="Insertions require inserted sequence, not an integer length",
        ))
    elif _INS_EMPTY.match(s):
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.ERROR,
            code=HGVSFixCode.INS_MISSING_SEQUENCE,
            message="Insertions require inserted sequence",
        ))
    if ":" not in s:
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.ERROR,
            code=HGVSFixCode.NO_COLON,
            message="No colon ':' found in HGVS string",
        ))
    if (parts := _parse(s)) and parts.kind and (allowed := _kinds_for_reference(parts.accession)):
        kinds, description = allowed
        if parts.kind not in kinds:
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.ERROR,
                code=HGVSFixCode.KIND_REFERENCE_MISMATCH,
                message=(f"{parts.accession} is a {description} accession, so it can't take "
                         f"'{parts.kind}.' (use {_or_list([k + '.' for k in sorted(kinds)])})"
                         + (_N_ON_CODING if parts.kind == "n" and "c" in kinds else "")),
            ))
    elif parts and parts.kind in _GENOMIC_KINDS and _looks_like_gene_symbol(parts.reference):
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.ERROR,
            code=HGVSFixCode.KIND_REFERENCE_MISMATCH,
            message=(f"{parts.reference} looks like a gene symbol, so it can't take '{parts.kind}.': "
                     f"{parts.kind}. positions are on a genomic accession (eg NC_000013.11)"),
        ))
    elif parts and parts.kind in _TRANSCRIPT_KINDS | _PROTEIN_KINDS and _is_chromosome_name(parts.reference):
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.ERROR,
            code=HGVSFixCode.KIND_REFERENCE_MISMATCH,
            message=(f"{parts.reference} is a chromosome name, so it can't take '{parts.kind}.': "
                     + ("p. positions are on a protein accession (eg NP_000050.3)" if parts.kind in _PROTEIN_KINDS
                        else f"{parts.kind}. positions are on a transcript accession (eg NM_000059.4)")),
        ))
    if _MISSING_REFERENCE.match(s):
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.ERROR,
            code=HGVSFixCode.MISSING_REFERENCE_SEQUENCE,
            message=f"Missing reference sequence (string starts with '{s[0]}')",
        ))
    return fixes


# ---------------------------------------------------------------------------
# Cleaning operations
#
# Each op has the signature (s: str) -> tuple[str, list[HGVSFix]]: it returns
# the (possibly unchanged) string and a list of HGVSFix describing what it did
# (empty if it made no change), or an ERROR for a problem it found but could not
# fix without guessing.  Ops are pure and order-dependent: _PIPELINE
# fixes the canonical order; selection only filters which ops run. The ops here
# are lexical (they work on the raw string); the structural ops further down
# work on its parsed fields.
# ---------------------------------------------------------------------------

def _op_strip_protein_suffix(s: str) -> tuple[str, list[HGVSFix]]:
    # e.g. "NM_000059.4:c.316+5G>A p.Arg106*" → "NM_000059.4:c.316+5G>A"
    if (m := _P_HGVS_REMOVAL.match(s)) and _HAS_NUCLEOTIDE_ALLELE.search(m.group("main")):
        new_s = m.group("main")
        return new_s, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.STRIPPED_PROTEIN_SUFFIX,
            message="Removed trailing protein (p.) annotation",
            original=s, fixed=new_s,
        )]
    return s, []


def _op_strip_whitespace(s: str) -> tuple[str, list[HGVSFix]]:
    # Strip non-printable characters and all whitespace
    s2 = "".join(c for c in s if c.isprintable()).replace(" ", "")
    if s2 != s:
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.STRIPPED_WHITESPACE,
            message="Removed whitespace and non-printable characters",
            original=s, fixed=s2,
        )]
    return s, []


def _op_strip_surrounding_punctuation(s: str) -> tuple[str, list[HGVSFix]]:
    # Strip stray wrapping quotes/backticks and trailing separators left by
    # copy/paste or search boxes, e.g. "NM_000059.4:c.123A>G'" or
    # "NM_000059.4:c.123A>G;", "NM_000059.4:c.123A>G." or a footnote mark
    # ("NM_000059.4:c.123A>G†"). No valid HGVS string starts/ends with these.
    # Surrounding whitespace on its own is left for the whitespace step to strip
    # and report (this op runs first).
    s2 = s.strip().strip("'\"`").rstrip("/;,.\\\u2020\u2021").strip()
    if s2 != s.strip():
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.STRIPPED_SURROUNDING_PUNCTUATION,
            message="Stripped surrounding quotes/punctuation",
            original=s, fixed=s2,
        )]
    return s, []


def _op_fix_multiple_colon(s: str) -> tuple[str, list[HGVSFix]]:
    # Collapse a run of 2+ colons, e.g. "NM_000059.4::c.123" or ":::" → single colon
    if "::" in s:
        s2 = re.sub(r"::+", ":", s)
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.FIXED_MULTIPLE_COLON,
            message="Fixed repeated colon ':'", original=s, fixed=s2,
        )]
    return s, []


def _op_fix_multiple_dot(s: str) -> tuple[str, list[HGVSFix]]:
    # Collapse a run of 2+ dots, e.g. "NM_000059..4" or "c...123" → single dot
    if ".." in s:
        s2 = re.sub(r"\.\.+", ".", s)
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.FIXED_MULTIPLE_DOT,
            message="Fixed repeated dot '.'", original=s, fixed=s2,
        )]
    return s, []


def _op_fix_multiple_underscore(s: str) -> tuple[str, list[HGVSFix]]:
    # Collapse a run of 2+ underscores, e.g. "NM__000059" or "NM___000059" → "NM_"
    if "__" in s:
        s2 = re.sub(r"_{2,}", "_", s)
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.FIXED_MULTIPLE_UNDERSCORE,
            message="Fixed repeated underscore '_'", original=s, fixed=s2,
        )]
    return s, []


def _op_fix_prefix_colon(s: str) -> tuple[str, list[HGVSFix]]:
    # Misplaced colon in the transcript prefix, e.g. "NM:_000059.4" → "NM_000059.4"
    s2 = _PREFIX_COLON_UNDERSCORE.sub(lambda m: m.group(1) + "_", s)
    if s2 != s:
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.FIXED_PREFIX_COLON,
            message="Fixed misplaced colon in transcript prefix",
            original=s, fixed=s2,
        )]
    return s, []


def _op_fix_multiple_kind(s: str) -> tuple[str, list[HGVSFix]]:
    # Collapse a repeated kind prefix, e.g. "c.c." or "c.c.c." → "c."
    fixes = []
    for kind in _KIND_LETTERS:
        if f"{kind}.{kind}." in s.lower():
            # (?:c\.){2,} matches a run of 2+ of this kind prefix
            s2 = re.sub(rf"(?:{kind}\.){{2,}}", f"{kind}.", s, flags=re.IGNORECASE)
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=HGVSFixCode.FIXED_MULTIPLE_KIND,
                message=f"Fixed repeated kind notation '{kind}.'",
                original=s, fixed=s2,
            ))
            s = s2
        # The kind repeated on the second position of a range, eg
        # "c.67+3_c.67+6del" -> "c.67+3_67+6del"
        if f":{kind}." in s.lower():
            s2 = re.sub(rf"(?<=[\d)])_{kind}\.(?=[\d*\-(])", "_", s, flags=re.IGNORECASE)
            if s2 != s:
                fixes.append(HGVSFix(
                    severity=HGVSFixSeverity.WARNING,
                    code=HGVSFixCode.FIXED_MULTIPLE_KIND,
                    message=f"Removed kind notation '{kind}.' repeated inside a range",
                    original=s, fixed=s2,
                ))
                s = s2
    return s, fixes


_UNICODE_DASHES = re.compile("[\u2010-\u2015\u2212]")


def _op_replace_unicode_dash(s: str) -> tuple[str, list[HGVSFix]]:
    # A typographic dash (hyphen, en/em dash, minus sign) from a word processor in
    # place of "-", eg "c.67\u20102A>G" -> "c.67-2A>G"
    s2 = _UNICODE_DASHES.sub("-", s)
    if s2 != s:
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.REPLACED_UNICODE_DASH,
            message="Replaced typographic dash with '-'",
            original=s, fixed=s2,
        )]
    return s, []


def _op_add_n_prefix(s: str) -> tuple[str, list[HGVSFix]]:
    # Restore a dropped leading "N" on a RefSeq accession, eg "M_001234" →
    # "NM_001234", "G_012337" → "NG_012337", "R_002196" → "NR_002196". Only fires
    # when "N" + the first letter is a real genomic/transcript prefix and an
    # underscore follows, so gene symbols (no "_" in the second position) are
    # never touched.
    if len(s) >= 2 and s[1] == "_" and s[0].upper() in _N_PREFIX_SECOND_LETTERS:
        s2 = "N" + s
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.ADDED_N_PREFIX,
            message="Added missing 'N' prefix",
            original=s, fixed=s2,
        )]
    return s, []


def _split_reference_allele(s: str) -> tuple[str, str, str]:
    """ (reference, separator, allele): the allele is what follows the last ':',
        or, with no colon, starts at the kind letter and its dot ("c."). With
        neither the whole string is the reference. """
    reference, sep, allele = s.rpartition(":")
    if sep:
        return reference, sep, allele
    if m := _KIND_DOT.search(s):
        return s[:m.start()], "", s[m.start():]
    return s, "", ""


def _brackets_balanced(s: str) -> bool:
    return s.count("(") == s.count(")")


def _op_strip_unbalanced_brackets(s: str) -> tuple[str, list[HGVSFix]]:
    # Only act when brackets are genuinely unbalanced (a stray "(" or ")" from
    # copy/paste). Do NOT strip when they balance: multiple balanced pairs are
    # valid HGVS - a gene symbol in parens plus uncertain-range notation, e.g.
    # "NM_004006.2(DMD):c.(4071+1_4072-1)_(5154+1_5155-1)del".
    if _brackets_balanced(s):
        return s, []
    # In the reference a bracket is the gene-symbol wrapper or a stray, so
    # stripping is safe (a glued "NM_000059.4BRCA2" is re-wrapped by structure
    # reconstruction). In the allele a bracket marks an uncertain position
    # ("c.(100_200)del") or a predicted consequence ("p.(Glu23fs)"), and
    # stripping it would make an uncertain variant look certain. The only
    # unambiguous repair there is a bracket that closes/opens nothing: a
    # trailing ")" or "(", or a "(" before the kind letter. Anything else is left
    # alone and reported as an ERROR.
    reference, sep, allele = _split_reference_allele(s)
    fixes = []
    if not _brackets_balanced(reference):
        reference = reference.replace("(", "").replace(")", "")
    if not _brackets_balanced(allele):
        if allele.endswith((")", "(")) and _brackets_balanced(allele[:-1]):
            allele = allele[:-1]
        elif allele.startswith("(") and _brackets_balanced(allele[1:]):
            allele = allele[1:]
        else:
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.ERROR,
                code=HGVSFixCode.UNBALANCED_BRACKETS,
                message=f"Unbalanced brackets in '{allele}' not stripped: a bracket in the "
                        f"allele marks an uncertain position",
                original=s,
            ))
    s2 = f"{reference}{sep}{allele}"
    if s2 != s:
        fixes.insert(0, HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.STRIPPED_UNBALANCED_BRACKETS,
            message="Stripped unbalanced brackets",
            original=s, fixed=s2,
        ))
    return s2, fixes


def _op_add_transcript_underscore(s: str) -> tuple[str, list[HGVSFix]]:
    # Add a missing underscore to any RefSeq accession prefix, eg "NM000059" →
    # "NM_000059", "NR002196" → "NR_002196", "XM005260" → "XM_005260". Requires
    # 6+ digits so gene symbols (NR3C1, ZP3, ...) are left alone.
    s2 = _TRANSCRIPT_NO_UNDERSCORE.sub(r"\1_\2", s)
    if s2 != s:
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.ADDED_TRANSCRIPT_UNDERSCORE,
            message="Added missing underscore to transcript prefix",
            original=s, fixed=s2,
        )]
    return s, []


def _op_uppercase_bases(s: str) -> tuple[str, list[HGVSFix]]:
    fixes = []
    if _RNA_ALLELE.search(s):
        return s, fixes  # r. alleles use lowercase bases
    # Uppercase bases in ref>alt substitutions (e.g. "c.100a>g" → "c.100A>G")
    s2 = _REF_ALT_NUC.sub(
        lambda m: m.group("ref").upper() + ">" + m.group("alt").upper(), s
    )
    if s2 != s:
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.UPPERCASED_BASES,
            message="Uppercased bases in ref>alt",
            original=s, fixed=s2,
        ))
        s = s2

    # Uppercase bases in del/ins/dup operations (e.g. "c.100delg" → "c.100delG")
    s2 = _DEL_INS_DUP_NUC.sub(_fix_del_ins_bases, s)
    if s2 != s:
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.UPPERCASED_BASES,
            message="Uppercased bases in del/ins/dup",
            original=s, fixed=s2,
        ))
        s = s2
    return s, fixes


def _op_strip_leading_junk(s: str) -> tuple[str, list[HGVSFix]]:
    # Strip junk before the transcript accession, e.g. "#NM_...", ":NM_...",
    # "GRCh38.p2 NM_...". Only fires when a transcript prefix actually follows,
    # so a leading gene symbol (e.g. "BRCA1(NM_...)") is never touched.
    s2 = _LEADING_JUNK.sub("", s)
    if s2 != s and re.match(_ANY_ACCESSION_PREFIX, s2, re.IGNORECASE):
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.STRIPPED_LEADING_JUNK,
            message="Stripped leading junk before transcript", original=s, fixed=s2,
        )]
    return s, []


def _op_fix_gene_wrapper(s: str) -> tuple[str, list[HGVSFix]]:
    # Normalise a gene symbol wedged between extra colons / stray parens:
    #   "tx:(GENE):c."  "tx:(GENE)c."  "tx:GENE:c."  -> "tx(GENE):c."
    #   "tx(GENE:):c."                               -> "tx(GENE):c."
    #   "tx(GENE):GENE:c."                           -> "tx(GENE):c."
    #   "(NM_000059.4):c."                           -> "NM_000059.4:c."
    #   "BRCA1(NM_000059.4)c."                       -> "BRCA1(NM_000059.4):c."
    s2 = _GENE_PARENS_COLON.sub(r"(\1)", s)
    s2 = _GENE_REPEATED_AFTER_PARENS.sub(r"(\1):", s2)
    s2 = _GENE_WRAPPER.sub(r"\1(\2):\3", s2)
    s2 = _PARENS_WRAPPED_ACCESSION.sub(r"\1", s2)
    s2 = _MISSING_COLON_BEFORE_KIND.sub(r"\1:\2", s2)
    if s2 != s:
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.FIXED_GENE_WRAPPER,
            message="Fixed transcript(gene) wrapper punctuation",
            original=s, fixed=s2,
        )]
    return s, []


@lru_cache(maxsize=1)
def _chromosome_accessions() -> dict[str, str]:
    """ Chromosome name (eg "13", "X", "MT") for each NC_ accession of the human
        assemblies cdot supports, keyed by both versioned and unversioned accession """
    from bioutils.assemblies import make_ac_name_map

    ac_to_chrom = {}
    for assembly in ("GRCh37.p13", "GRCh38.p14", "T2T-CHM13v2.0"):
        for ac, name in make_ac_name_map(assembly).items():
            if ac.startswith("NC_"):
                ac_to_chrom[ac] = name.upper()
                ac_to_chrom[ac.split(".")[0]] = name.upper()
    return ac_to_chrom


def _op_strip_contig_name(s: str) -> tuple[str, list[HGVSFix]]:
    # A chromosome name before its own accession is redundant:
    #   "chr13 NC_000013.11:g.32316528del" -> "NC_000013.11:g.32316528del"
    # Only stripped when the name matches the accession's chromosome, so a
    # mismatch (eg "chr7 NC_000013.11"), which suggests a copy/paste mix-up, is
    # left alone rather than guessing which one was meant.
    def _strip_if_matching(m: re.Match) -> str:
        chrom = m.group(1)[3:].upper()
        if chrom == "M":
            chrom = "MT"
        accession = m.group(2) + (m.group(3) or "")
        if _chromosome_accessions().get(accession.upper()) == chrom:
            return ""
        return m.group(0)

    s2 = _CONTIG_NAME_BEFORE_ACCESSION.sub(_strip_if_matching, s)
    if s2 != s:
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.STRIPPED_CONTIG_NAME,
            message="Stripped chromosome name before its accession",
            original=s, fixed=s2,
        )]
    return s, []


def _op_strip_build_name(s: str) -> tuple[str, list[HGVSFix]]:
    # A genome build name before the reference is not part of the HGVS string:
    #   "hg19 NM_000059.4:c.68del"         -> "NM_000059.4:c.68del"
    #   "GRCh38 BRCA2 NM_000059.4:c.68del" -> "BRCA2 NM_000059.4:c.68del"
    # (the gene symbol is then wrapped by the symbol-before-accession step).
    s2 = _BUILD_NAME_BEFORE_REFERENCE.sub("", s)
    if s2 != s:
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.STRIPPED_BUILD_NAME,
            message="Stripped genome build name before the reference",
            original=s, fixed=s2,
        )]
    return s, []


def _op_fix_symbol_before_accession(s: str) -> tuple[str, list[HGVSFix]]:
    # Gene symbol written before the accession, separated by whitespace, ':' or '.':
    #   "APC NM_001407446.1:c.-40G>A"  -> "NM_001407446.1(APC):c.-40G>A"
    #   "APC:NM_001407446.1:c.-40G>A"  -> "NM_001407446.1(APC):c.-40G>A"
    #   "BRCA2 NM_000059.4 c.68del"    -> "NM_000059.4(BRCA2) c.68del"
    # (the remaining space is stripped and the colon restored by later steps).
    # An accession that already has its gene symbol keeps it: a matching symbol
    # before it is dropped ("BRCA2 NM_000059.4(BRCA2):c.68del"), a different one
    # is left alone rather than guessing which is right.
    def _move_symbol(m: re.Match) -> str:
        symbol, accession = m.groups()
        rest = s[m.end():]
        if not rest.startswith("("):
            return f"{accession}({symbol})"
        if rest[1:].upper().startswith(symbol.upper() + ")"):
            return accession
        return m.group(0)

    s2 = _SYMBOL_BEFORE_ACCESSION.sub(_move_symbol, s)
    if s2 != s:
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.WRAPPED_GENE_SYMBOL,
            message="Moved gene symbol before transcript into transcript(gene) form",
            original=s, fixed=s2,
        )]
    return s, []


def _op_fix_space_before_kind(s: str) -> tuple[str, list[HGVSFix]]:
    # Whitespace in place of the reference:kind colon after a bare gene symbol,
    # e.g. "BRCA2 c.68del" -> "BRCA2:c.68del". Must run before whitespace
    # stripping, which would glue the symbol to the kind letter and leave
    # accession-lookalike symbols (eg NR2E3) unrepairable.
    s2 = _SPACE_BEFORE_KIND.sub(r"\1:", s)
    if s2 != s:
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.FIXED_SEPARATOR_TYPO,
            message="Replaced space between reference and kind with ':'",
            original=s, fixed=s2,
        )]
    return s, []


def _op_fix_separator_typo(s: str) -> tuple[str, list[HGVSFix]]:
    # Comma ("c,1811") or colon ("c:1811") in place of the kind dot, period
    # in place of the substitution '>' ("1030C.T" -> "1030C>T"), and semicolon
    # or underscore in place of the accession:allele colon ("NM_000059.4;c.68del").
    s2 = _KIND_COMMA.sub(r":\1.", s)
    s2 = _KIND_COLON.sub(r":\1.", s2)
    s2 = _SUB_PERIOD.sub(r"\1>\2", s2)
    s2 = _ACCESSION_KIND_SEPARATOR.sub(r"\1:", s2)
    if s2 != s:
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.FIXED_SEPARATOR_TYPO,
            message="Fixed comma/colon/period separator typo", original=s, fixed=s2,
        )]
    return s, []


def _op_drop_del_dup_count(s: str) -> tuple[str, list[HGVSFix]]:
    # Drop a redundant base count after a ranged del/dup, normalising to plain
    # del/dup, e.g. "c.1315_1337dup23" -> "c.1315_1337dup". When both positions
    # are plain integers the range's length is known, so a count that disagrees
    # with it ("c.123_127del3", a 5 bp range) is a contradiction: the string is
    # left alone and reported as an ERROR rather than the count being dropped.
    # A range with an intronic offset or UTR position ("c.100+5_200-3del5") has
    # no length without the transcript, so its count is dropped unchecked.
    head, sep, allele = s.rpartition(":")
    target = allele if sep else s
    if _DEL_DUP_RANGE_COUNT.search(target):
        if m := _DEL_DUP_PLAIN_RANGE_COUNT.search(target):
            start, end, count = int(m["start"]), int(m["end"]), int(m["count"])
            length = end - start + 1
            if count != length:
                return s, [HGVSFix(
                    severity=HGVSFixSeverity.ERROR,
                    code=HGVSFixCode.DEL_DUP_COUNT_MISMATCH,
                    message=f"Base count {count} after {m['op'].lower()} contradicts the "
                            f"{length} bp range {start}_{end}",
                    original=s,
                )]
        new_allele = _DEL_DUP_RANGE_COUNT.sub(r"\1", target)
        s2 = f"{head}:{new_allele}" if sep else new_allele
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.DROPPED_DEL_DUP_COUNT,
            message="Dropped redundant base count after ranged del/dup",
            original=s, fixed=s2,
        )]
    return s, []


def _op_reconstruct_structure(s: str) -> tuple[str, list[HGVSFix]]:
    # Loose parse and serialise: uppercase the accession prefix, lowercase the
    # kind letter, insert the ':' and '.' where missing, and put the gene symbol
    # in its parentheses. The pattern's accession prefixes include their
    # underscore, so a gene symbol that merely starts with prefix letters plus a
    # digit (eg NR3C1, NC2) is read as a gene symbol, never as "NR3(C1)".
    if (parsed := _parse_loose(s)) is None:
        return s, []
    parts, uppercased_prefix = parsed
    s2 = str(parts)
    if s2 == s:
        return s, []
    if uppercased_prefix:
        code, message = HGVSFixCode.UPPERCASED_TRANSCRIPT, "Uppercased transcript prefix"
    else:
        code, message = (HGVSFixCode.RECONSTRUCTED_STRUCTURE,
                         "Reconstructed HGVS with correct casing and structure")
    return s2, [HGVSFix(
        severity=HGVSFixSeverity.WARNING,
        code=code, message=message,
        original=s, fixed=s2,
    )]


# ---------------------------------------------------------------------------
# Structural operations
#
# Each takes and returns _HGVSParts, and is wrapped by _on_fields into the
# string-op signature above, so it only runs on a string with a colon.
# ---------------------------------------------------------------------------

@_on_fields
def _op_drop_empty_version(parts: _HGVSParts) -> tuple[_HGVSParts, list[HGVSFix]]:
    # An accession with a dot but no version digits, e.g. "NM_001754.:c.749G>A"
    # -> "NM_001754:c.749G>A". The version was dropped; an unversioned accession
    # is valid HGVS, so drop the dangling dot rather than invent a version.
    if parts.gene_symbol is not None or not _EMPTY_VERSION_ACCESSION.fullmatch(parts.accession):
        return parts, []
    fixed = replace(parts, accession=parts.accession[:-1])
    return fixed, [HGVSFix(
        severity=HGVSFixSeverity.WARNING,
        code=HGVSFixCode.DROPPED_EMPTY_VERSION,
        message="Dropped dangling dot of an empty transcript version",
        original=str(parts), fixed=str(fixed),
    )]


@_on_fields
def _op_lowercase_mutation_type(parts: _HGVSParts) -> tuple[_HGVSParts, list[HGVSFix]]:
    # Lowercase mutation types (DUP→dup, DEL→del, INS→ins, INV→inv), which also
    # handles DELINS. Only the allele is changed: gene symbols such as INSR, INVS,
    # INSL3 or DELEC1 contain a mutation-type substring and must not be corrupted
    # (e.g. "NM_000208.4(INSR):c.215A>G" must not become "insR").
    fixes = []
    for mt in ["INS", "DEL", "DUP", "INV"]:
        if mt in parts.nomen:
            fixed = replace(parts, nomen=parts.nomen.replace(mt, mt.lower()))
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=HGVSFixCode.LOWERCASED_MUTATION_TYPE,
                message=f"Lowercased mutation type '{mt}'",
                original=str(parts), fixed=str(fixed),
            ))
            parts = fixed
    return parts, fixes


@_on_fields
def _op_add_missing_kind(parts: _HGVSParts) -> tuple[_HGVSParts, list[HGVSFix]]:
    # Insert a missing kind letter + dot after the colon, choosing the kind from
    # the reference-sequence type:
    #   mRNA / Ensembl transcript (NM_, XM_, ENST) → "c."  eg "NM_001754.5:557T>A"
    #   non-coding RNA           (NR_, XR_)        → "n."  eg "NR_002196.2:601A>G"
    #   genomic                  (NC_, NG_, ...)   → "g."  eg "NC_000007.13:1G>T"
    # A reference spanning lines is not one reference, so it is left alone.
    kind = _KIND_FOR_TYPE.get(_reference_type(parts.accession))
    if not kind or "\n" in parts.reference or not _ALLELE_WITHOUT_KIND.match(parts.allele):
        return parts, []
    fixed = replace(parts, kind=kind, nomen=parts.allele)
    return fixed, [HGVSFix(
        severity=HGVSFixSeverity.WARNING,
        code=HGVSFixCode.ADDED_MISSING_KIND,
        message=f"Added missing '{kind}.' kind prefix",
        original=str(parts), fixed=str(fixed),
    )]


@_on_fields
def _op_fix_gene_transcript(parts: _HGVSParts) -> tuple[_HGVSParts, list[HGVSFix]]:
    # Fix the common clinical mistake of swapping gene symbol and transcript, and
    # uppercase a lowercase accession prefix. Empty parentheses are dropped.
    # Ported from HGVSMatcher.fix_gene_transcript() in vg_code/hgvs/hgvs_matcher.py.
    fixes = []
    original = parts.reference
    accession, gene_symbol = parts.accession, parts.gene_symbol or None

    transcript_ok = _looks_like_transcript(accession)
    # A genomic accession is the right reference for g./m. HGVS, so it isn't swapped out
    genomic_ref_for_genomic_kind = parts.kind in ("g", "m") and _GENOMIC_ACCESSION.fullmatch(accession)

    if not transcript_ok and gene_symbol and not genomic_ref_for_genomic_kind:
        uc_gene = gene_symbol.upper()
        if _looks_like_transcript(uc_gene):
            accession, gene_symbol = gene_symbol, accession
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=HGVSFixCode.SWAPPED_GENE_TRANSCRIPT,
                message="Swapped gene/transcript",
                original=original,
                fixed=f"{accession}({gene_symbol})",
            ))
            transcript_ok = _looks_like_transcript(accession)
        elif _looks_like_hgvs_prefix(uc_gene):
            gene_symbol = uc_gene
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=HGVSFixCode.UPPERCASED_HGVS_PREFIX,
                message="Upper cased HGVS prefix",
                original=original,
                fixed=f"{accession}({gene_symbol})",
            ))

    if not transcript_ok:
        # Any known accession prefix, not only transcripts: a lowercase genomic or
        # protein accession (eg "ng_012337.3", "np_000050.2") is uppercased too.
        uc = _uppercase_accession_prefix(accession)
        if uc != accession:
            accession = uc
            if _looks_like_transcript(uc):
                code, message = HGVSFixCode.UPPERCASED_TRANSCRIPT, "Upper cased transcript"
            else:
                code, message = HGVSFixCode.UPPERCASED_HGVS_PREFIX, "Upper cased HGVS prefix"
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=code,
                message=message,
                original=original,
                fixed=accession,
            ))

    return replace(parts, accession=accession, gene_symbol=gene_symbol), fixes


@_on_fields
def _op_drop_genomic_ref_in_parens(parts: _HGVSParts) -> tuple[_HGVSParts, list[HGVSFix]]:
    # Drop a genomic accession (NC_/NG_/NW_, ...) sitting in the gene-symbol
    # parentheses of a transcript reference, e.g.
    # "NM_000059.4(NC_000013.11):c.68del" -> "NM_000059.4:c.68del". The
    # transcript on its own is a complete reference, so the genomic accession is
    # redundant and (being neither a gene symbol nor a valid selector here) only
    # blocks parsing. "NC_000013.11(NM_000059.4):c.68del" ends up the same way,
    # as the gene/transcript fix has already swapped it.
    #
    # For g./m. HGVS the genomic accession is the reference instead, so the
    # transcript is dropped: "NC_000013.11(NM_000059.4):g.32316528del" and
    # "NM_000059.4(NC_000013.11):g.32316528del" -> "NC_000013.11:g.32316528del".
    if parts.kind in ("g", "m") and parts.gene_symbol is not None:
        accessions = [parts.accession, parts.gene_symbol]
        genomic = [a for a in accessions if _GENOMIC_ACCESSION.fullmatch(a)]
        transcript = [a for a in accessions if _TRANSCRIPT_ACCESSION.fullmatch(a)]
        if len(genomic) == 1 and len(transcript) == 1:
            fixed = replace(parts, accession=genomic[0], gene_symbol=None)
            return fixed, [HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=HGVSFixCode.DROPPED_TRANSCRIPT_FOR_GENOMIC_KIND,
                message=f"Dropped transcript {transcript[0]}: a {parts.kind}. variant is on the genomic accession",
                original=str(parts), fixed=str(fixed),
            )]
    if (parts.gene_symbol is None
            or not _TRANSCRIPT_ACCESSION.fullmatch(parts.accession)
            or not _GENOMIC_ACCESSION.fullmatch(parts.gene_symbol)):
        return parts, []
    fixed = replace(parts, gene_symbol=None)
    return fixed, [HGVSFix(
        severity=HGVSFixSeverity.WARNING,
        code=HGVSFixCode.DROPPED_GENOMIC_REF_IN_PARENS,
        message="Dropped genomic accession from the gene-symbol parentheses",
        original=str(parts), fixed=str(fixed),
    )]


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

class HGVSCleanOpWarning(UserWarning):
    """ Warned by clean_hgvs() when the selected ops include one without an op
        it relies on, so it will miss inputs it would otherwise fix. """


@dataclass(frozen=True)
class _PipelineStep:
    """ One op of the pipeline. requires lists the earlier ops whose output this
        op relies on: selecting it without them still runs it, but with a
        HGVSCleanOpWarning. """
    op: HGVSCleanOp
    fn: _StringOp
    requires: frozenset = frozenset()


# Restore a dropped "N" or underscore in the accession prefix ("M_" -> "NM_",
# "NM000059" -> "NM_000059"), needed by the ops that look for a full prefix
_FULL_ACCESSION_PREFIX = frozenset({HGVSCleanOp.ADD_N_PREFIX, HGVSCleanOp.ADD_TRANSCRIPT_UNDERSCORE})

# Canonical pipeline, in the order the ops must run: lexical ops on the string,
# then the loose parse that reconstructs the structure, then structural ops on
# the fields. The allele-level lexical fixes (del/dup counts, brackets, base
# case) stay where they are, as moving them would change what the structural
# ops are given. clean_hgvs() filters this by the caller's `ops` set but never
# reorders it.
_PIPELINE: list[_PipelineStep] = [
    # Lexical ops
    # Wrapping quotes first: the ops below that are anchored to the start of the
    # string ('^') can't see an accession or symbol behind a leading quote.
    _PipelineStep(HGVSCleanOp.STRIP_SURROUNDING_PUNCTUATION, _op_strip_surrounding_punctuation),
    _PipelineStep(HGVSCleanOp.REPLACE_UNICODE_DASH, _op_replace_unicode_dash),
    _PipelineStep(HGVSCleanOp.STRIP_PROTEIN_SUFFIX, _op_strip_protein_suffix),
    _PipelineStep(HGVSCleanOp.STRIP_CONTIG_NAME, _op_strip_contig_name),
    _PipelineStep(HGVSCleanOp.STRIP_BUILD_NAME, _op_strip_build_name),
    _PipelineStep(HGVSCleanOp.STRIP_LEADING_JUNK, _op_strip_leading_junk),
    # Would otherwise take a build name for a gene symbol
    _PipelineStep(HGVSCleanOp.FIX_SYMBOL_BEFORE_ACCESSION, _op_fix_symbol_before_accession,
                  requires=frozenset({HGVSCleanOp.STRIP_BUILD_NAME})),
    # Before whitespace stripping, which would glue the symbol to the kind
    _PipelineStep(HGVSCleanOp.FIX_SPACE_BEFORE_KIND, _op_fix_space_before_kind),
    _PipelineStep(HGVSCleanOp.STRIP_WHITESPACE, _op_strip_whitespace),
    _PipelineStep(HGVSCleanOp.ADD_N_PREFIX, _op_add_n_prefix),
    _PipelineStep(HGVSCleanOp.ADD_TRANSCRIPT_UNDERSCORE, _op_add_transcript_underscore),
    _PipelineStep(HGVSCleanOp.FIX_MULTIPLE_COLON, _op_fix_multiple_colon),
    _PipelineStep(HGVSCleanOp.FIX_GENE_WRAPPER, _op_fix_gene_wrapper),
    _PipelineStep(HGVSCleanOp.FIX_MULTIPLE_DOT, _op_fix_multiple_dot),
    _PipelineStep(HGVSCleanOp.FIX_MULTIPLE_UNDERSCORE, _op_fix_multiple_underscore),
    _PipelineStep(HGVSCleanOp.FIX_PREFIX_COLON, _op_fix_prefix_colon),
    _PipelineStep(HGVSCleanOp.FIX_MULTIPLE_KIND, _op_fix_multiple_kind),
    _PipelineStep(HGVSCleanOp.FIX_SEPARATOR_TYPO, _op_fix_separator_typo,
                  requires=_FULL_ACCESSION_PREFIX),
    _PipelineStep(HGVSCleanOp.DROP_EMPTY_VERSION, _op_drop_empty_version,
                  requires=_FULL_ACCESSION_PREFIX),
    _PipelineStep(HGVSCleanOp.DROP_DEL_DUP_COUNT, _op_drop_del_dup_count),
    _PipelineStep(HGVSCleanOp.STRIP_UNBALANCED_BRACKETS, _op_strip_unbalanced_brackets),
    # Parse
    _PipelineStep(HGVSCleanOp.RECONSTRUCT_STRUCTURE, _op_reconstruct_structure,
                  requires=_FULL_ACCESSION_PREFIX),
    # Structural ops
    # Needs the ':' that reconstruction may have just inserted
    _PipelineStep(HGVSCleanOp.LOWERCASE_MUTATION_TYPE, _op_lowercase_mutation_type,
                  requires=frozenset({HGVSCleanOp.RECONSTRUCT_STRUCTURE})),
    _PipelineStep(HGVSCleanOp.UPPERCASE_BASES, _op_uppercase_bases),
    _PipelineStep(HGVSCleanOp.ADD_MISSING_KIND, _op_add_missing_kind),
    _PipelineStep(HGVSCleanOp.FIX_GENE_TRANSCRIPT, _op_fix_gene_transcript),
    # An "NC_…(NM_…)" reference is swapped to "NM_…(NC_…)" by the gene/transcript
    # fix, so dropping the genomic parenthetical after it collapses both the
    # broken and the swapped form to the same parseable transcript reference.
    _PipelineStep(HGVSCleanOp.DROP_GENOMIC_REF_IN_PARENS, _op_drop_genomic_ref_in_parens,
                  requires=frozenset({HGVSCleanOp.FIX_GENE_TRANSCRIPT})),
]

# All cleaning ops, the default for clean_hgvs(). Use set algebra for a
# blocklist, e.g. ALL_CLEAN_OPS - {HGVSCleanOp.STRIP_PROTEIN_SUFFIX}.
ALL_CLEAN_OPS: frozenset = frozenset(HGVSCleanOp)


def _warn_missing_prerequisites(ops) -> None:
    """ Warn (HGVSCleanOpWarning) about each selected op whose prerequisites are
        not all selected """
    problems = []
    for step in _PIPELINE:
        if step.op in ops and (missing := step.requires - ops):
            names = ", ".join(sorted(m.name for m in missing))
            problems.append(f"{step.op.name} (needs {names})")
    if problems:
        warnings.warn(
            "clean_hgvs ops selected without the ops they rely on, so some inputs they "
            f"would fix are left alone: {'; '.join(problems)}",
            HGVSCleanOpWarning, stacklevel=3)


# ---------------------------------------------------------------------------
# Public: clean_hgvs
# ---------------------------------------------------------------------------

def clean_hgvs(
    hgvs_string: str,
    ops: Optional[set] = None,
    raise_on_errors: bool = False,
    validate: bool = True,
) -> tuple[str, list[HGVSFix]]:
    """
    Attempt to fix common formatting problems in an HGVS string.

    Args:
        hgvs_string: the input to clean.
        ops: the subset of HGVSCleanOp to apply. None (default) runs them all
            (ALL_CLEAN_OPS). For an allowlist, pass just the ops you want, e.g.
            {HGVSCleanOp.STRIP_WHITESPACE}. For a blocklist, use set algebra,
            e.g. ALL_CLEAN_OPS - {HGVSCleanOp.STRIP_PROTEIN_SUFFIX}. Ops always
            run in the canonical pipeline order regardless of set order. Selecting
            an op without an op it relies on (eg LOWERCASE_MUTATION_TYPE without
            RECONSTRUCT_STRUCTURE) warns with HGVSCleanOpWarning.
        raise_on_errors: if True, raise HGVSInputError on the first ERROR-level
            fix instead of returning it in the list.
        validate: if True (default), run post-cleaning validation and include any
            ERROR-level problems in the returned fixes. Pass False to skip it
            (e.g. when you only want string normalisation).

    Returns:
        (cleaned_string, fixes)
        where fixes is a list of HGVSFix describing every change made and
        any validation errors found.

    The returned string is always the best attempt at a valid HGVS string.
    If ERROR-level fixes are present (and raise_on_errors=False), the
    returned string may still be invalid.

    Example::

        cleaned, fixes = clean_hgvs("nm_000059.4 c.316+5G>A")
        # cleaned == "NM_000059.4:c.316+5G>A"
        # fixes == [
        #   HGVSFix(WARNING, STRIPPED_WHITESPACE, ...),
        #   HGVSFix(WARNING, UPPERCASED_TRANSCRIPT, ...),
        # ]
    """
    if ops is None:
        ops = ALL_CLEAN_OPS
    else:
        _warn_missing_prerequisites(ops)

    fixes: list[HGVSFix] = []
    s = hgvs_string

    for step in _PIPELINE:
        if step.op in ops:
            s, op_fixes = step.fn(s)
            fixes.extend(op_fixes)

    # Validation — run after cleaning so both the fix and any remaining
    # problem are captured in one call
    if validate:
        fixes.extend(_validate(s))

    if raise_on_errors:
        errors = [f for f in fixes if f.severity == HGVSFixSeverity.ERROR]
        if errors:
            raise HGVSInputError(errors[0].message)

    return s, fixes


# ---------------------------------------------------------------------------
# Public: transcript version fallback
# ---------------------------------------------------------------------------

class VersionStrategy(Enum):
    UP_THEN_DOWN = "up_then_down"  # prefer higher versions first, then lower
    CLOSEST      = "closest"       # alternate ±1 from the requested version
    LATEST       = "latest"        # highest available, ignoring requested version


def rank_transcript_versions(
    requested_version: int,
    available_versions: list[int],
    strategy: VersionStrategy = VersionStrategy.UP_THEN_DOWN,
) -> list[int]:
    """
    Return ``available_versions`` ordered best-first for the strategy.

    This is the full ordering behind :func:`get_best_transcript_version` (which
    returns only the single best). Callers that try candidate versions in turn
    want the whole ranking, so it is exposed here.

    Ported from HGVSMatcher._get_sort_key_transcript_version_and_methods()
    in vg_code/hgvs/hgvs_matcher.py (minus the method/ClinGen sorting,
    which is VG-specific).
    """
    if strategy == VersionStrategy.LATEST:
        return sorted(available_versions, reverse=True)

    def sort_key(v: int) -> tuple:
        distance = abs(requested_version - v)
        prefer_later = v < requested_version  # False=later, True=earlier → sort asc = later first
        if strategy == VersionStrategy.CLOSEST:
            return (distance, prefer_later)
        else:  # UP_THEN_DOWN
            return (prefer_later, distance)

    return sorted(available_versions, key=sort_key)


# Backwards-compatible alias for the previous private name.
_rank_versions = rank_transcript_versions


def get_best_transcript_version(
    accession: str,
    requested_version: int,
    available_versions: list[int],
    strategy: VersionStrategy = VersionStrategy.UP_THEN_DOWN,
) -> tuple[int, Optional[HGVSFix]]:
    """
    Find the best available transcript version when the requested one is absent.

    Returns:
        (best_version, fix_or_None)
        fix is None if requested_version was found exactly.
        fix is HGVSFix(WARNING, USED_ADJACENT_VERSION, ...) otherwise.

    The caller is responsible for deciding whether to use the substituted
    version and for updating the HGVS string accordingly.  This function
    never silently changes anything — the caller must act on the returned fix.

    This picks a version from a bare ``list[int]`` and is deliberately
    structure-agnostic: it has no transcript records, so it cannot judge whether
    the substitution is *coordinate-safe*.  That check lives in the provider-aware
    layer — ``resolve_transcript_version`` calls
    ``cdot.hgvs.version_safety.is_version_substitution_safe_by_accession`` (which
    fetches the two versions' records and runs the pure check on them)
    and upgrades this plain USED_ADJACENT_VERSION
    fix to a coordinate-safe / unverified / refused outcome (#28).

    Example::

        version, fix = get_best_transcript_version("NM_000059", 4, [1, 2, 3, 5, 6])
        if fix:
            fixes.append(fix)
            hgvs_string = hgvs_string.replace(
                f"NM_000059.{4}", f"NM_000059.{version}"
            )
    """
    if requested_version in available_versions:
        return requested_version, None

    if not available_versions:
        raise HGVSInputError(
            f"No versions of {accession} are available in the data provider"
        )

    ranked = rank_transcript_versions(requested_version, available_versions, strategy)
    best = ranked[0]
    fix = HGVSFix(
        severity=HGVSFixSeverity.WARNING,
        code=HGVSFixCode.USED_ADJACENT_VERSION,
        message=(
            f"Transcript version {accession}.{requested_version} not found; "
            f"using .{best} instead"
        ),
        original=f"{accession}.{requested_version}",
        fixed=f"{accession}.{best}",
    )
    return best, fix
